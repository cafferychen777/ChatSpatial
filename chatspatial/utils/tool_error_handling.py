"""
Enhanced error handling for MCP tools according to specification.

Tools should return errors in the result object, not as protocol-level errors.
This allows LLMs to see and potentially handle the error.

IMAGE HANDLING NOTE (Updated 2024-12):
- Currently all visualizations save to disk and return file paths (DIRECT_EMBED_THRESHOLD = 0)
- ImageContent handling code is preserved for future use when MCP protocol improves
- The type-aware error handling ensures correct error format for different return types
"""

import traceback
from functools import wraps
from typing import get_type_hints


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
                # IMPORTANT: ImageContent handling
                # ImageContent objects must be returned directly without wrapping
                # FastMCP requires the raw ImageContent object for proper conversion
                # Currently not used (DIRECT_EMBED_THRESHOLD = 0) but preserved for future use
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
                    # Return error as readable text message (not placeholder image)
                    # The Union[ImageContent, str] return type allows str returns
                    # FastMCP auto-converts str to TextContent for display
                    error_msg = f"Visualization Error: {str(e)}"

                    # Include traceback for detailed debugging (except for common user errors)
                    if include_traceback and not isinstance(e, (ValueError, KeyError)):
                        tb = traceback.format_exc()
                        error_msg += f"\n\nTechnical Details:\n{tb}"

                    return error_msg  # Return error message as string

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
                # Inline the error result creation (was create_error_result)
                error_content = [{"type": "text", "text": f"Error: {str(e)}"}]
                if include_traceback and not isinstance(e, ValueError):
                    tb = traceback.format_exc()
                    error_content.append({"type": "text", "text": f"Traceback:\n{tb}"})
                return {"content": error_content, "isError": True}

        # All MCP tools are async functions, return async wrapper
        return async_wrapper

    return decorator
