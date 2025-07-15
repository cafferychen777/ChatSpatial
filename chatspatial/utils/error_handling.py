"""
Enhanced error handling for MCP server
"""

from typing import Any, Callable, Dict, Optional, Union, List
from functools import wraps
import traceback
import scanpy as sc
from ..mcp.errors import ErrorType, format_mcp_error


# Custom exception classes
class SpatialMCPError(Exception):
    """Base exception for spatial analysis errors"""
    pass


class DataNotFoundError(SpatialMCPError):
    """Exception raised when required data is not found"""
    pass


class InvalidParameterError(SpatialMCPError):
    """Exception raised for invalid parameters"""
    pass


class ProcessingError(SpatialMCPError):
    """Exception raised during processing"""
    pass


class DataCompatibilityError(SpatialMCPError):
    """Exception raised for data compatibility issues"""
    pass


def validate_adata(adata: sc.AnnData, required_keys: Dict[str, Union[str, List[str]]], 
                  context: Optional[Any] = None) -> None:
    """Validate AnnData object has required keys
    
    Args:
        adata: AnnData object to validate
        required_keys: Dict of required keys by category (obs, var, obsm, etc.)
        context: Optional context for logging
    
    Raises:
        DataNotFoundError: If required keys are missing
    """
    missing = []
    
    for category, keys in required_keys.items():
        if isinstance(keys, str):
            keys = [keys]
            
        attr = getattr(adata, category, None)
        if attr is None:
            missing.extend([f"{category}.{k}" for k in keys])
            continue
            
        for key in keys:
            if hasattr(attr, 'columns'):  # DataFrame
                if key not in attr.columns:
                    missing.append(f"{category}.{key}")
            elif hasattr(attr, 'keys'):  # Dict-like
                if key not in attr.keys():
                    missing.append(f"{category}.{key}")
            else:
                missing.append(f"{category}.{key}")
    
    if missing:
        raise DataNotFoundError(f"Missing required keys: {', '.join(missing)}")


def handle_error(error: Exception, context: Optional[Any] = None) -> None:
    """Handle errors with appropriate logging
    
    Args:
        error: The exception to handle
        context: Optional context for logging
    """
    if context:
        if isinstance(error, SpatialMCPError):
            context.warning(f"Analysis error: {str(error)}")
        else:
            context.error(f"Unexpected error: {str(error)}")
            context.debug(traceback.format_exc())


def try_except_with_feedback(func: Callable) -> Callable:
    """Decorator to add error handling with feedback
    
    Args:
        func: Function to wrap
        
    Returns:
        Wrapped function with error handling
    """
    @wraps(func)
    async def wrapper(*args, **kwargs):
        context = kwargs.get('context')
        try:
            return await func(*args, **kwargs)
        except Exception as e:
            handle_error(e, context)
            raise
    
    return wrapper


def mcp_error_handler(func: Callable) -> Callable:
    """Decorator to handle errors according to MCP specification"""
    @wraps(func)
    async def wrapper(*args, **kwargs):
        try:
            return await func(*args, **kwargs)
        except ValueError as e:
            # Check if it's already an MCP-formatted error
            error_str = str(e)
            if error_str.startswith('{"code":'):
                raise  # Re-raise as is
            
            # Convert to MCP error format
            if "not found" in error_str.lower():
                error = format_mcp_error(
                    ErrorType.DATASET_NOT_FOUND,
                    error_str
                )
            elif "invalid" in error_str.lower() or "format" in error_str.lower():
                error = format_mcp_error(
                    ErrorType.INVALID_DATA_FORMAT,
                    error_str
                )
            else:
                error = format_mcp_error(
                    ErrorType.INVALID_PARAMS,
                    error_str
                )
            raise ValueError(error)
            
        except FileNotFoundError as e:
            error = format_mcp_error(
                ErrorType.DATASET_NOT_FOUND,
                f"File not found: {str(e)}",
                {"path": str(e.filename) if hasattr(e, 'filename') else None}
            )
            raise ValueError(error)
            
        except MemoryError as e:
            error = format_mcp_error(
                ErrorType.INTERNAL_ERROR,
                "Out of memory while processing data",
                {"error_type": "MemoryError"}
            )
            raise ValueError(error)
            
        except Exception as e:
            # Log the full traceback for debugging
            tb = traceback.format_exc()
            
            # Determine appropriate error type
            error_msg = str(e)
            if "analysis" in error_msg.lower():
                error_type = ErrorType.ANALYSIS_FAILED
            elif "visual" in error_msg.lower() or "plot" in error_msg.lower():
                error_type = ErrorType.VISUALIZATION_ERROR
            elif "reference" in error_msg.lower():
                error_type = ErrorType.REFERENCE_DATA_ERROR
            else:
                error_type = ErrorType.INTERNAL_ERROR
            
            error = format_mcp_error(
                error_type,
                f"Operation failed: {error_msg}",
                {
                    "error_type": type(e).__name__,
                    "traceback": tb if len(tb) < 1000 else tb[-1000:]  # Limit traceback size
                }
            )
            raise ValueError(error)
    
    return wrapper


def sync_mcp_error_handler(func: Callable) -> Callable:
    """Synchronous version of the error handler"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ValueError as e:
            # Check if it's already an MCP-formatted error
            error_str = str(e)
            if error_str.startswith('{"code":'):
                raise  # Re-raise as is
            
            # Convert to MCP error format
            if "not found" in error_str.lower():
                error = format_mcp_error(
                    ErrorType.DATASET_NOT_FOUND,
                    error_str
                )
            else:
                error = format_mcp_error(
                    ErrorType.INVALID_PARAMS,
                    error_str
                )
            raise ValueError(error)
            
        except Exception as e:
            error = format_mcp_error(
                ErrorType.INTERNAL_ERROR,
                f"Operation failed: {str(e)}",
                {"error_type": type(e).__name__}
            )
            raise ValueError(error)
    
    return wrapper