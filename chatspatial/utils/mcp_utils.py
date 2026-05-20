"""
MCP utilities for ChatSpatial.

Tools for MCP server: error handling decorator and output suppression.

Error Handling Design:
======================
All tool errors are raised as exceptions, which FastMCP converts to
``CallToolResult(isError=True)`` protocol responses automatically.

The ``mcp_tool_error_handler`` decorator enriches error messages before
they reach FastMCP:

User-understandable errors (clean message, no traceback):
- ParameterError, DataError, DataNotFoundError, DataCompatibilityError
- DependencyError, ValueError (legacy)

Code/algorithm errors (message + traceback for debugging):
- ProcessingError, all other exceptions
"""

import logging
import os
import traceback
import warnings
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from functools import wraps

from .exceptions import (
    DataCompatibilityError,
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
)

# Exceptions that don't need traceback (message is self-explanatory)
# These are "user errors" - the error message is sufficient for understanding
USER_ERRORS = (
    ParameterError,
    DataError,
    DataNotFoundError,
    DataCompatibilityError,
    DependencyError,
    ValueError,  # Legacy compatibility
)


# =============================================================================
# Output Suppression
# =============================================================================
@contextmanager
def suppress_output():
    """
    Context manager to suppress stdout, stderr, warnings, and logging.

    Usage:
        with suppress_output():
            noisy_function()
    """
    old_level = logging.getLogger().level
    logging.getLogger().setLevel(logging.ERROR)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        try:
            with open(os.devnull, "w") as devnull:
                with redirect_stdout(devnull), redirect_stderr(devnull):
                    yield
        finally:
            logging.getLogger().setLevel(old_level)


# =============================================================================
# MCP Tool Error Handler
# =============================================================================
def mcp_tool_error_handler(include_traceback: bool = True):
    """
    Decorator for MCP tools that enriches error messages before re-raising.

    All exceptions are re-raised for FastMCP to convert into
    ``CallToolResult(isError=True)`` protocol responses. The decorator
    only adds traceback detail to non-user errors for debugging.

    Args:
        include_traceback: Append traceback to non-user error messages.
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
                if include_traceback:
                    tb = traceback.format_exc()
                    # Enrich message in-place — preserves exception type
                    # and any custom attributes regardless of constructor
                    e.args = (f"{e}\n\nTraceback:\n{tb}",)
                raise

        return wrapper

    return decorator
