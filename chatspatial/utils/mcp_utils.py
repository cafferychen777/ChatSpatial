"""
MCP utilities for ChatSpatial.

Tools for MCP server: error handling decorator and output suppression.

Error Handling Design:
======================
All tool errors are raised as exceptions, which MCPServer converts to
``CallToolResult(isError=True)`` protocol responses automatically.

Since MCP SDK 2.1.0, MCPServer forwards the text of an exception to the client
only when the exception is a ``ToolError``. Any other exception is reported as
``Error executing tool <name>`` and its text stays in the server log. The
``expose_anticipated_errors`` wrapper, applied to every tool registered on the
ChatSpatial server, therefore re-raises anticipated errors as ``ToolError``
with the original type name and message, chained to the original exception:

Client-visible errors (``CLIENT_VISIBLE_ERRORS``, type and message, no traceback):
- The ChatSpatialError hierarchy (ParameterError, DataError, DataNotFoundError,
  DataCompatibilityError, ProcessingError, DependencyError)
- ValueError (legacy), FileNotFoundError and PermissionError, which ChatSpatial
  still raises for invalid inputs, missing input files and unwritable outputs

Unexpected errors:
- Reach the client only as ``Error executing tool <name>`` (MCP SDK >= 2.1.0;
  SDK 2.0.x appends the exception message)
- Never carry a traceback unless explicitly enabled for local debugging

The ``mcp_tool_error_handler`` decorator records the complete traceback of
every non-user error in the server log before the exception reaches MCPServer.
"""

import contextvars
import inspect
import logging
import sys
import threading
import traceback
from collections.abc import AsyncIterator, Callable, Iterable, Iterator
from contextlib import asynccontextmanager, contextmanager
from functools import wraps
from typing import Any, TextIO

from mcp.server.mcpserver.exceptions import ToolError

from .exceptions import (
    ChatSpatialError,
    DataError,
    DependencyError,
    ParameterError,
)

logger = logging.getLogger(__name__)

# Exceptions that don't need traceback (message is self-explanatory)
# These are "user errors" - the error message is sufficient for understanding
USER_ERRORS = (
    ParameterError,
    DataError,
    DependencyError,
    ValueError,  # Legacy compatibility
)

# Exceptions whose message is written for the user and must reach the MCP
# client. ProcessingError belongs here as well: it is raised deliberately with a
# message that names the failed step, even though its traceback is logged.
CLIENT_VISIBLE_ERRORS = (
    ChatSpatialError,
    ValueError,  # Legacy compatibility
    FileNotFoundError,  # Missing input paths in load_data and reload_data
    PermissionError,  # Unwritable output directories
)


# =============================================================================
# Output Suppression
# =============================================================================
_OUTPUT_SUPPRESSED = contextvars.ContextVar(
    "chatspatial_output_suppressed",
    default=False,
)
_STREAM_STATE_LOCK = threading.RLock()
_ACTIVE_SUPPRESSORS = 0
_ORIGINAL_STREAMS: tuple[TextIO, TextIO] | None = None


class _ContextAwareTextStream:
    """Forward writes unless the current task or thread suppresses output."""

    def __init__(self, target: TextIO) -> None:
        self._target = target

    def write(self, data: str) -> int:
        if _OUTPUT_SUPPRESSED.get():
            return len(data)
        written = self._target.write(data)
        return len(data) if written is None else written

    def writelines(self, lines: Iterable[str]) -> None:
        if _OUTPUT_SUPPRESSED.get():
            for _line in lines:
                pass
            return
        self._target.writelines(lines)

    def flush(self) -> None:
        if not _OUTPUT_SUPPRESSED.get():
            self._target.flush()

    def __getattr__(self, name: str) -> Any:
        return getattr(self._target, name)


_STREAM_PROXIES: tuple[_ContextAwareTextStream, _ContextAwareTextStream] | None = None


def _install_context_aware_streams() -> None:
    global _ACTIVE_SUPPRESSORS, _ORIGINAL_STREAMS, _STREAM_PROXIES

    with _STREAM_STATE_LOCK:
        if _ACTIVE_SUPPRESSORS == 0:
            _ORIGINAL_STREAMS = (sys.stdout, sys.stderr)
            _STREAM_PROXIES = (
                _ContextAwareTextStream(sys.stdout),
                _ContextAwareTextStream(sys.stderr),
            )
            sys.stdout, sys.stderr = _STREAM_PROXIES
        _ACTIVE_SUPPRESSORS += 1


def _restore_context_aware_streams() -> None:
    global _ACTIVE_SUPPRESSORS, _ORIGINAL_STREAMS, _STREAM_PROXIES

    with _STREAM_STATE_LOCK:
        _ACTIVE_SUPPRESSORS -= 1
        if _ACTIVE_SUPPRESSORS != 0:
            return

        if _ORIGINAL_STREAMS is not None and _STREAM_PROXIES is not None:
            stdout_proxy, stderr_proxy = _STREAM_PROXIES
            original_stdout, original_stderr = _ORIGINAL_STREAMS
            if sys.stdout is stdout_proxy:
                sys.stdout = original_stdout
            if sys.stderr is stderr_proxy:
                sys.stderr = original_stderr
        _ORIGINAL_STREAMS = None
        _STREAM_PROXIES = None


@contextmanager
def suppress_output() -> Iterator[None]:
    """Suppress Python stdout and stderr for the current execution context.

    A process-global redirect is unsafe across concurrent MCP tasks. The stream
    proxy remains shared only while needed, while a ``ContextVar`` decides
    whether each task or thread is suppressed. ``asyncio.to_thread`` propagates
    this context automatically.

    Usage:
        with suppress_output():
            noisy_function()
    """
    token = _OUTPUT_SUPPRESSED.set(True)
    try:
        _install_context_aware_streams()
    except BaseException:
        _OUTPUT_SUPPRESSED.reset(token)
        raise

    try:
        yield
    finally:
        _OUTPUT_SUPPRESSED.reset(token)
        _restore_context_aware_streams()


@asynccontextmanager
async def suppress_output_async() -> AsyncIterator[None]:
    """Async form of :func:`suppress_output` for scopes containing ``await``."""
    with suppress_output():
        yield


# =============================================================================
# MCP Tool Error Handler
# =============================================================================
def mcp_tool_error_handler(include_traceback: bool = False):
    """
    Decorator for MCP tools that enriches error messages before re-raising.

    All exceptions are re-raised for MCPServer to convert into
    ``CallToolResult(isError=True)`` protocol responses. The decorator
    logs non-user errors with traceback details on the server. Tracebacks are
    excluded from client-visible messages by default to avoid disclosing local
    paths and implementation details.

    Args:
        include_traceback: Append traceback to non-user error messages. Enable
            only for trusted local debugging.
    """

    def decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            try:
                return await func(*args, **kwargs)
            except USER_ERRORS:
                # User errors already have clear messages — re-raise as-is
                raise
            except Exception as e:
                logger.exception("MCP tool %s failed", func.__name__)
                if include_traceback:
                    tb = traceback.format_exc()
                    # Enrich message in-place — preserves exception type
                    # and any custom attributes regardless of constructor
                    e.args = (f"{e}\n\nTraceback:\n{tb}",)
                raise

        return wrapper

    return decorator


def expose_anticipated_errors(func: Callable[..., Any]) -> Callable[..., Any]:
    """Re-raise client-visible errors as ``ToolError`` for MCP registration.

    ``ToolError`` is the MCP SDK's supported signal for an anticipated tool
    failure: MCPServer returns its message to the client, whereas the SDK
    replaces the message of any other exception with a generic text since
    release 2.1.0. The original exception remains available as ``__cause__``
    for server-side logging. Unexpected exceptions pass through unchanged, so
    the SDK keeps their details on the server.
    """

    def translate(exc: Exception) -> ToolError:
        return ToolError(f"{type(exc).__name__}: {exc}")

    if inspect.iscoroutinefunction(func):

        @wraps(func)
        async def async_wrapper(*args: Any, **kwargs: Any) -> Any:
            try:
                return await func(*args, **kwargs)
            except ToolError:
                raise
            except CLIENT_VISIBLE_ERRORS as exc:
                raise translate(exc) from exc

        return async_wrapper

    @wraps(func)
    def sync_wrapper(*args: Any, **kwargs: Any) -> Any:
        try:
            return func(*args, **kwargs)
        except ToolError:
            raise
        except CLIENT_VISIBLE_ERRORS as exc:
            raise translate(exc) from exc

    return sync_wrapper
