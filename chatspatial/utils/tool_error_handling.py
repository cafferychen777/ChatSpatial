"""
Enhanced error handling for MCP tools according to specification.

Tools should return errors in the result object, not as protocol-level errors.
This allows LLMs to see and potentially handle the error.

!!!!!!!!!! CRITICAL WARNING - IMAGE HANDLING !!!!!!!!!!
This module contains CRITICAL code for handling Image objects from FastMCP.
DO NOT modify the Image handling logic in mcp_tool_error_handler without
reading /docs/CRITICAL_IMAGE_DISPLAY_BUG.md first!
A bug in this code caused images to display as object strings for 2 WEEKS!
!!!!!!!!!! CRITICAL WARNING - IMAGE HANDLING !!!!!!!!!!
"""

from typing import Any, Dict, List, Union, Optional
from dataclasses import dataclass
import traceback
from functools import wraps
import asyncio


@dataclass
class ToolResult:
    """Standard tool result format according to MCP specification"""
    content: List[Dict[str, Any]]
    isError: bool = False
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format expected by MCP"""
        return {
            "content": self.content,
            "isError": self.isError
        }


def create_error_result(error: Exception, include_traceback: bool = True) -> ToolResult:
    """Create a standardized error result for tools"""
    error_content = [
        {
            "type": "text",
            "text": f"Error: {str(error)}"
        }
    ]
    
    if include_traceback and not isinstance(error, ValueError):
        # Include traceback for non-value errors
        tb = traceback.format_exc()
        error_content.append({
            "type": "text", 
            "text": f"Traceback:\n{tb}"
        })
    
    return ToolResult(
        content=error_content,
        isError=True
    )


def create_success_result(content: Union[str, Dict[str, Any], Any]) -> ToolResult:
    """Create a standardized success result for tools"""
    if isinstance(content, str):
        return ToolResult(
            content=[{
                "type": "text",
                "text": content
            }],
            isError=False
        )
    elif isinstance(content, dict):
        # If it's already in the right format
        if "content" in content and isinstance(content["content"], list):
            return ToolResult(
                content=content["content"],
                isError=False
            )
        else:
            # Convert dict to text representation
            import json
            return ToolResult(
                content=[{
                    "type": "text",
                    "text": json.dumps(content, indent=2)
                }],
                isError=False
            )
    else:
        # For model objects, convert to dict first
        if hasattr(content, 'model_dump'):
            # Pydantic v2
            data = content.model_dump()
        elif hasattr(content, 'dict'):
            # Pydantic v1
            data = content.dict()
        else:
            # Fallback to string representation
            data = str(content)
        
        return create_success_result(data)


def mcp_tool_error_handler(include_traceback: bool = True):
    """
    Decorator for MCP tools to handle errors according to specification.
    
    This ensures that tool errors are returned in the result object
    rather than being raised as exceptions.
    """
    def decorator(func):
        @wraps(func)
        async def async_wrapper(*args, **kwargs):
            try:
                result = await func(*args, **kwargs)
                # If the function already returns a ToolResult, use it
                if isinstance(result, ToolResult):
                    return result.to_dict()
                # !!!!!!!!!! CRITICAL WARNING - DO NOT MODIFY !!!!!!!!!!
                # Image objects MUST be returned directly without any wrapping!
                # FastMCP needs to see the raw Image object to convert it properly.
                # If you wrap Image objects in dictionaries or ToolResult, 
                # Claude Desktop will show "<Image object at 0x...>" instead of the actual image.
                # This bug took 2 WEEKS to find and fix. DO NOT CHANGE THIS!
                # See /docs/CRITICAL_IMAGE_DISPLAY_BUG.md for full details.
                # !!!!!!!!!! CRITICAL WARNING - DO NOT MODIFY !!!!!!!!!!
                from mcp.server.fastmcp.utilities.types import Image
                if isinstance(result, Image):
                    return result  # MUST return raw Image object!
                # Otherwise, wrap the result
                return create_success_result(result).to_dict()
            except Exception as e:
                # Return error in the result object
                return create_error_result(e, include_traceback).to_dict()
        
        @wraps(func)
        def sync_wrapper(*args, **kwargs):
            try:
                result = func(*args, **kwargs)
                # If the function already returns a ToolResult, use it
                if isinstance(result, ToolResult):
                    return result.to_dict()
                # !!!!!!!!!! CRITICAL WARNING - DO NOT MODIFY !!!!!!!!!!
                # Image objects MUST be returned directly without any wrapping!
                # FastMCP needs to see the raw Image object to convert it properly.
                # If you wrap Image objects in dictionaries or ToolResult, 
                # Claude Desktop will show "<Image object at 0x...>" instead of the actual image.
                # This bug took 2 WEEKS to find and fix. DO NOT CHANGE THIS!
                # See /docs/CRITICAL_IMAGE_DISPLAY_BUG.md for full details.
                # !!!!!!!!!! CRITICAL WARNING - DO NOT MODIFY !!!!!!!!!!
                from mcp.server.fastmcp.utilities.types import Image
                if isinstance(result, Image):
                    return result  # MUST return raw Image object!
                # Otherwise, wrap the result
                return create_success_result(result).to_dict()
            except Exception as e:
                # Return error in the result object
                return create_error_result(e, include_traceback).to_dict()
        
        # Return appropriate wrapper based on function type
        if asyncio.iscoroutinefunction(func):
            return async_wrapper
        else:
            return sync_wrapper
    
    return decorator


# Convenience functions for common error scenarios
def dataset_not_found_error(data_id: str) -> ToolResult:
    """Create a dataset not found error result"""
    return ToolResult(
        content=[{
            "type": "text",
            "text": f"Error: Dataset '{data_id}' not found. Please load a dataset first using the 'load_data' tool."
        }],
        isError=True
    )


def invalid_parameter_error(param_name: str, expected: str, got: Any) -> ToolResult:
    """Create an invalid parameter error result"""
    return ToolResult(
        content=[{
            "type": "text",
            "text": f"Error: Invalid parameter '{param_name}'. Expected {expected}, got {type(got).__name__}: {got}"
        }],
        isError=True
    )


def analysis_failed_error(analysis_type: str, reason: str) -> ToolResult:
    """Create an analysis failed error result"""
    return ToolResult(
        content=[{
            "type": "text",
            "text": f"Error: {analysis_type} analysis failed. Reason: {reason}"
        }],
        isError=True
    )


def file_operation_error(operation: str, path: str, reason: str) -> ToolResult:
    """Create a file operation error result"""
    return ToolResult(
        content=[{
            "type": "text",
            "text": f"Error: Failed to {operation} file '{path}'. Reason: {reason}"
        }],
        isError=True
    )


