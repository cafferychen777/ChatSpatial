"""
MCP Error Handling System

This module provides error types and formatting functions for MCP protocol compliance.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Optional


class ErrorType(Enum):
    """MCP Error Types"""

    INVALID_REQUEST = -32600
    METHOD_NOT_FOUND = -32601
    INVALID_PARAMS = -32602
    INTERNAL_ERROR = -32603

    # Custom error codes (range -32000 to -32099)
    DATASET_NOT_FOUND = -32000
    INVALID_DATA_FORMAT = -32001
    ANALYSIS_FAILED = -32002
    VISUALIZATION_ERROR = -32003
    REFERENCE_DATA_ERROR = -32004

    # Additional error codes from error_handling.py usage
    DATA_NOT_FOUND = -32000  # Alias for DATASET_NOT_FOUND
    INVALID_DATASET = -32002  # Alias for compatibility
    PARAMETER_ERROR = -32004  # Alias for compatibility


@dataclass
class MCPError:
    """MCP Error Format"""

    code: int
    message: str
    data: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        result = {"code": self.code, "message": self.message}
        if self.data:
            result["data"] = self.data
        return result


def format_mcp_error(
    error_type: ErrorType, message: str, data: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """Format error according to MCP specification

    Args:
        error_type: The type of error from ErrorType enum
        message: Human-readable error message
        data: Optional additional error data

    Returns:
        Dictionary formatted according to MCP error specification
    """
    return MCPError(code=error_type.value, message=message, data=data).to_dict()
