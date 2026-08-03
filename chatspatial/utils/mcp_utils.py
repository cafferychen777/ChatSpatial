"""
MCP utilities for ChatSpatial.

Tools for MCP server: error handling decorator and output suppression.

Error Handling Design:
======================
All tool errors are raised as exceptions, which MCPServer converts to
``CallToolResult(isError=True)`` protocol responses automatically.

The ``mcp_tool_error_handler`` decorator enriches error messages before
they reach MCPServer:

User-understandable errors (clean message, no traceback):
- ParameterError, DataError, DataNotFoundError, DataCompatibilityError
- DependencyError, ValueError (legacy)

Code/algorithm errors:
- Preserve the concise exception message for the client
- Record the complete traceback in server logs
- Expose a traceback only when explicitly enabled for local debugging
"""

import contextvars
import logging
import sys
import threading
import traceback
from collections.abc import AsyncIterator, Iterable, Iterator
from contextlib import asynccontextmanager, contextmanager
from functools import wraps
from typing import Any, TextIO

from .exceptions import (
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
