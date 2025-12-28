"""
MCP utilities for ChatSpatial.

Tools for MCP server: error handling decorator, output suppression.
"""

import io
import logging
import traceback
import warnings
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from functools import wraps
from typing import get_type_hints


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
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()

        try:
            with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
                yield
        finally:
            logging.getLogger().setLevel(old_level)


# =============================================================================
# MCP Tool Error Handler
# =============================================================================
def _get_return_type_category(func) -> str:
    """
    Determine what category of return type a function has.

    Returns one of: "image", "basemodel", "str", "simple", "unknown"
    """
    try:
        hints = get_type_hints(func)
        return_type = hints.get("return", None)

        if return_type is None:
            return "unknown"

        type_str = str(return_type)

        if "ImageContent" in type_str:
            return "image"

        # Pydantic model types from models/analysis.py
        basemodel_types = [
            "SpatialDataset",
            "PreprocessingResult",
            "AnnotationResult",
            "SpatialStatisticsResult",
            "DifferentialExpressionResult",
            "CNVResult",
            "DeconvolutionResult",
            "CellCommunicationResult",
            "EnrichmentResult",
            "RNAVelocityResult",
            "TrajectoryResult",
            "IntegrationResult",
            "SpatialDomainResult",
            "SpatialVariableGenesResult",
        ]

        for model_type in basemodel_types:
            if model_type in type_str:
                return "basemodel"

        if type_str in ["<class 'str'>", "str", "typing.Optional[str]"]:
            return "str"

        return "simple"

    except Exception:
        return "unknown"


def mcp_tool_error_handler(include_traceback: bool = True):
    """
    Decorator for MCP tools to handle errors gracefully.

    Errors are returned in the result object (not raised as exceptions),
    allowing LLMs to see and potentially handle them.
    """

    def decorator(func):
        return_type_category = _get_return_type_category(func)

        @wraps(func)
        async def wrapper(*args, **kwargs):
            try:
                result = await func(*args, **kwargs)

                # Return results directly - let FastMCP handle serialization
                from mcp.types import ImageContent
                from pydantic import BaseModel

                if isinstance(result, ImageContent):
                    return result
                if isinstance(result, tuple) and len(result) >= 1:
                    if isinstance(result[0], ImageContent):
                        return result
                if isinstance(result, BaseModel):
                    return result

                return result

            except Exception as e:
                error_msg = str(e)

                if return_type_category == "image":
                    # Return error as text (Union[ImageContent, str] allows this)
                    msg = f"Visualization Error: {error_msg}"
                    if include_traceback and not isinstance(e, (ValueError, KeyError)):
                        msg += f"\n\nDetails:\n{traceback.format_exc()}"
                    return msg

                elif return_type_category == "basemodel":
                    # Re-raise for FastMCP to handle
                    raise

                elif return_type_category == "str":
                    return f"Error: {error_msg}"

                else:
                    # Return error dict for simple types
                    content = [{"type": "text", "text": f"Error: {error_msg}"}]
                    if include_traceback and not isinstance(e, ValueError):
                        content.append({"type": "text", "text": f"Traceback:\n{traceback.format_exc()}"})
                    return {"content": content, "isError": True}

        return wrapper

    return decorator
