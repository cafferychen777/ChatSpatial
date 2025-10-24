"""
Enhanced error handling for MCP tools according to specification.

Tools should return errors in the result object, not as protocol-level errors.
This allows LLMs to see and potentially handle the error.

!!!!!!!!!! CRITICAL WARNING - IMAGE HANDLING !!!!!!!!!!
This module contains CRITICAL code for handling Image objects from FastMCP.
DO NOT modify the Image handling logic in mcp_tool_error_handler!
A bug in this code caused images to display as object strings for 2 WEEKS!
!!!!!!!!!! CRITICAL WARNING - IMAGE HANDLING !!!!!!!!!!
"""

import traceback
from dataclasses import dataclass
from functools import wraps
from typing import Any, Dict, List, get_type_hints


@dataclass
class ToolResult:
    """Standard tool result format according to MCP specification"""

    content: List[Dict[str, Any]]
    isError: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format expected by MCP"""
        return {"content": self.content, "isError": self.isError}


def create_error_result(error: Exception, include_traceback: bool = True) -> ToolResult:
    """Create a standardized error result for tools"""
    error_content = [{"type": "text", "text": f"Error: {str(error)}"}]

    if include_traceback and not isinstance(error, ValueError):
        # Include traceback for non-value errors
        tb = traceback.format_exc()
        error_content.append({"type": "text", "text": f"Traceback:\n{tb}"})

    return ToolResult(content=error_content, isError=True)


def _check_return_type_category(func) -> str:
    """Check what category of type a function returns

    Args:
        func: Function to check

    Returns:
        One of: "image", "basemodel", "str", "simple", "unknown"
        - "image": Returns ImageContent
        - "basemodel": Returns Pydantic BaseModel (but not ImageContent)
        - "str": Returns str (must return str even on error for MCP schema)
        - "simple": Returns simple types (int, dict, list)
        - "unknown": Cannot determine type
    """
    try:
        # Get type hints
        hints = get_type_hints(func)
        return_type = hints.get("return", None)

        if return_type is None:
            return "unknown"

        # Convert to string for easier checking
        type_str = str(return_type)

        # Check if ImageContent is in the return type
        if "ImageContent" in type_str:
            return "image"

        # Check if it's a Pydantic BaseModel (common result types)
        # These are defined in models/analysis.py and models/data.py
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

        # Check if it returns str specifically (important for MCP schema compliance)
        if type_str in ["<class 'str'>", "str", "typing.Optional[str]"]:
            return "str"

        # Otherwise it's a simple type
        return "simple"

    except Exception:
        return "unknown"


def _create_error_placeholder_image(error: Exception):
    """Create a placeholder image displaying error information

    Args:
        error: The exception that occurred

    Returns:
        ImageContent object with error message
    """
    try:
        from ..utils.image_utils import create_placeholder_image

        # Format error message
        error_msg = f"Error: {str(error)}"

        # Truncate if too long
        if len(error_msg) > 200:
            error_msg = error_msg[:197] + "..."

        # Create placeholder image
        return create_placeholder_image(message=error_msg, figsize=(8, 4))

    except Exception:
        # If we can't create placeholder image, return None and fall back to error dict
        return None


def mcp_tool_error_handler(include_traceback: bool = True):
    """
    Decorator for MCP tools to handle errors according to specification.

    This ensures that tool errors are returned in the result object
    rather than being raised as exceptions.
    """

    def decorator(func):
        # Check what category of type the function returns
        return_type_category = _check_return_type_category(func)

        @wraps(func)
        async def async_wrapper(*args, **kwargs):
            try:
                result = await func(*args, **kwargs)
                # If the function already returns a ToolResult, use it
                if isinstance(result, ToolResult):
                    return result.to_dict()
                # !!!!!!!!!! CRITICAL WARNING - DO NOT MODIFY !!!!!!!!!!
                # ImageContent objects MUST be returned directly without any wrapping!
                # FastMCP needs to see the raw ImageContent object to convert it properly.
                # If you wrap ImageContent objects in dictionaries or ToolResult,
                # Claude Desktop will show "<ImageContent object at 0x...>" instead of the actual image.
                # This bug took 2 WEEKS to find and fix. DO NOT CHANGE THIS!
                # !!!!!!!!!! CRITICAL WARNING - DO NOT MODIFY !!!!!!!!!!
                from mcp.types import ImageContent
                from pydantic import BaseModel

                if isinstance(result, ImageContent):
                    return result  # MUST return raw ImageContent object!
                # Support for Tuple[ImageContent, EmbeddedResource] returns
                if isinstance(result, tuple) and len(result) >= 1:
                    if isinstance(result[0], ImageContent):
                        return result  # Let FastMCP handle the tuple!
                # MCP 1.10+ can handle Pydantic models directly - don't wrap them!
                if isinstance(result, BaseModel):
                    return result  # Let FastMCP serialize Pydantic models
                # MCP 1.17+ can handle simple types (str, int, list, dict) directly
                # Return them without wrapping to avoid validation errors
                return result
            except Exception as e:
                # Handle errors based on return type category
                if return_type_category == "image":
                    # For ImageContent, return error as placeholder image
                    # This ensures type consistency - MCP expects ImageContent, not error dict
                    error_image = _create_error_placeholder_image(e)
                    if error_image is not None:
                        return error_image  # Return placeholder ImageContent
                    # Fall through to re-raise if placeholder creation fails

                elif return_type_category == "basemodel":
                    # For BaseModel types, re-raise the exception
                    # Let FastMCP handle it at a higher level
                    # This prevents MCP schema validation errors
                    raise

                elif return_type_category == "str":
                    # For str return type, must return str (not dict) for MCP schema
                    # Simply return error message as string
                    return f"Error: {str(e)}"

                # For simple types or unknown, return error in the result object
                return create_error_result(e, include_traceback).to_dict()

        # All MCP tools are async functions, return async wrapper
        return async_wrapper

    return decorator


# Convenience functions for common error scenarios
