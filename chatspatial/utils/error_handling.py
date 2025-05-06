"""
Error handling utilities for spatial transcriptomics MCP.

This module provides standardized error handling and user feedback functions.
"""

import traceback
import sys
import logging
from typing import Optional, Dict, Any, Callable, TypeVar, Awaitable, Union
from functools import wraps
import inspect
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from mcp.server.fastmcp import Context

# Type variables for generic functions
T = TypeVar('T')
R = TypeVar('R')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('spatial_transcriptomics_mcp')


class SpatialMCPError(Exception):
    """Base class for all spatial transcriptomics MCP errors."""
    def __init__(self, message: str, details: Optional[Dict[str, Any]] = None):
        self.message = message
        self.details = details or {}
        super().__init__(message)


class DataNotFoundError(SpatialMCPError):
    """Error raised when a dataset is not found."""
    pass


class InvalidParameterError(SpatialMCPError):
    """Error raised when an invalid parameter is provided."""
    pass


class ProcessingError(SpatialMCPError):
    """Error raised when processing fails."""
    pass


class DependencyError(SpatialMCPError):
    """Error raised when a required dependency is missing."""
    pass


class DataCompatibilityError(SpatialMCPError):
    """Error raised when data is not compatible with the requested operation."""
    pass


async def handle_error(
    error: Exception,
    context: Optional[Context] = None,
    fallback_value: Optional[T] = None,
    error_message: Optional[str] = None,
    log_traceback: bool = True
) -> Optional[T]:
    """Handle an error and provide user feedback.
    
    Args:
        error: The exception to handle
        context: MCP context for user feedback
        fallback_value: Value to return if error is handled
        error_message: Custom error message to display
        log_traceback: Whether to log the full traceback
        
    Returns:
        Fallback value if provided, otherwise re-raises the exception
    """
    # Get error details
    error_type = type(error).__name__
    error_msg = str(error)
    if not error_message:
        error_message = f"{error_type}: {error_msg}"
    
    # Log the error
    if log_traceback:
        logger.error(f"Error: {error_message}\n{traceback.format_exc()}")
    else:
        logger.error(f"Error: {error_message}")
    
    # Provide user feedback through context
    if context:
        if isinstance(error, DependencyError):
            await context.warning(f"Missing dependency: {error_msg}")
            await context.info("Please install the required dependency and try again.")
        elif isinstance(error, DataNotFoundError):
            await context.warning(f"Data not found: {error_msg}")
        elif isinstance(error, InvalidParameterError):
            await context.warning(f"Invalid parameter: {error_msg}")
        elif isinstance(error, DataCompatibilityError):
            await context.warning(f"Data compatibility issue: {error_msg}")
        else:
            await context.warning(f"Error: {error_message}")
    
    # Return fallback value or re-raise
    if fallback_value is not None:
        return fallback_value
    raise error


def check_dependency(package_name: str, import_name: Optional[str] = None) -> bool:
    """Check if a dependency is installed.
    
    Args:
        package_name: Name of the package to check
        import_name: Name to use for importing (if different from package_name)
        
    Returns:
        True if dependency is installed, False otherwise
    """
    import_name = import_name or package_name
    try:
        __import__(import_name)
        return True
    except ImportError:
        return False


def require_dependency(package_name: str, import_name: Optional[str] = None) -> None:
    """Require a dependency to be installed.
    
    Args:
        package_name: Name of the package to check
        import_name: Name to use for importing (if different from package_name)
        
    Raises:
        DependencyError: If the dependency is not installed
    """
    if not check_dependency(package_name, import_name):
        raise DependencyError(
            f"Required dependency '{package_name}' is not installed. "
            f"Please install it with 'pip install {package_name}'."
        )


async def try_except_with_feedback(
    func: Callable[..., Awaitable[R]],
    context: Optional[Context],
    error_message: str,
    fallback_func: Optional[Callable[..., Awaitable[R]]] = None,
    *args,
    **kwargs
) -> R:
    """Execute a function with error handling and user feedback.
    
    Args:
        func: Async function to execute
        context: MCP context for user feedback
        error_message: Error message to display if function fails
        fallback_func: Fallback function to execute if main function fails
        *args: Arguments to pass to the function
        **kwargs: Keyword arguments to pass to the function
        
    Returns:
        Result of the function or fallback function
    """
    try:
        return await func(*args, **kwargs)
    except Exception as e:
        if context:
            await context.warning(f"{error_message}: {str(e)}")
            if fallback_func:
                await context.info("Using fallback method...")
        
        if fallback_func:
            return await fallback_func(*args, **kwargs)
        raise


def validate_adata(
    adata: ad.AnnData,
    require_spatial: bool = False,
    require_var_names: bool = False,
    require_obs_names: bool = False,
    min_cells: int = 0,
    min_genes: int = 0
) -> None:
    """Validate an AnnData object.
    
    Args:
        adata: AnnData object to validate
        require_spatial: Whether spatial coordinates are required
        require_var_names: Whether var_names are required
        require_obs_names: Whether obs_names are required
        min_cells: Minimum number of cells required
        min_genes: Minimum number of genes required
        
    Raises:
        DataCompatibilityError: If validation fails
    """
    # Check if AnnData is valid
    if adata is None:
        raise DataCompatibilityError("AnnData object is None")
    
    # Check dimensions
    if min_cells > 0 and adata.n_obs < min_cells:
        raise DataCompatibilityError(f"AnnData has too few cells: {adata.n_obs} < {min_cells}")
    
    if min_genes > 0 and adata.n_vars < min_genes:
        raise DataCompatibilityError(f"AnnData has too few genes: {adata.n_vars} < {min_genes}")
    
    # Check spatial coordinates
    if require_spatial:
        has_spatial = False
        if 'spatial' in adata.obsm:
            has_spatial = True
        elif any('spatial' in key for key in adata.obsm.keys()):
            has_spatial = True
        elif 'spatial' in adata.uns:
            has_spatial = True
        
        if not has_spatial:
            raise DataCompatibilityError("AnnData does not contain spatial coordinates")
    
    # Check var_names
    if require_var_names and (adata.var_names is None or len(adata.var_names) == 0):
        raise DataCompatibilityError("AnnData does not have valid var_names")
    
    # Check obs_names
    if require_obs_names and (adata.obs_names is None or len(adata.obs_names) == 0):
        raise DataCompatibilityError("AnnData does not have valid obs_names")


def format_error_details(error: Exception) -> Dict[str, Any]:
    """Format error details for logging and user feedback.
    
    Args:
        error: Exception to format
        
    Returns:
        Dictionary with error details
    """
    return {
        'error_type': type(error).__name__,
        'error_message': str(error),
        'traceback': traceback.format_exc(),
        'module': getattr(error, '__module__', 'unknown'),
        'timestamp': pd.Timestamp.now().isoformat()
    }


def is_valid_parameter(value: Any, param_type: type, min_value: Optional[float] = None, max_value: Optional[float] = None) -> bool:
    """Check if a parameter value is valid.
    
    Args:
        value: Parameter value to check
        param_type: Expected parameter type
        min_value: Minimum allowed value (for numeric types)
        max_value: Maximum allowed value (for numeric types)
        
    Returns:
        True if parameter is valid, False otherwise
    """
    # Check type
    if not isinstance(value, param_type):
        return False
    
    # Check range for numeric types
    if isinstance(value, (int, float)):
        if min_value is not None and value < min_value:
            return False
        if max_value is not None and value > max_value:
            return False
    
    return True


def validate_parameters(params: Dict[str, Any], validation_rules: Dict[str, Dict[str, Any]]) -> None:
    """Validate parameters against a set of rules.
    
    Args:
        params: Dictionary of parameters to validate
        validation_rules: Dictionary of validation rules
        
    Raises:
        InvalidParameterError: If validation fails
    """
    for param_name, rules in validation_rules.items():
        if param_name not in params:
            if rules.get('required', False):
                raise InvalidParameterError(f"Required parameter '{param_name}' is missing")
            continue
        
        value = params[param_name]
        param_type = rules.get('type')
        min_value = rules.get('min_value')
        max_value = rules.get('max_value')
        allowed_values = rules.get('allowed_values')
        
        if param_type and not isinstance(value, param_type):
            raise InvalidParameterError(
                f"Parameter '{param_name}' has invalid type: expected {param_type.__name__}, got {type(value).__name__}"
            )
        
        if isinstance(value, (int, float)):
            if min_value is not None and value < min_value:
                raise InvalidParameterError(
                    f"Parameter '{param_name}' is too small: {value} < {min_value}"
                )
            if max_value is not None and value > max_value:
                raise InvalidParameterError(
                    f"Parameter '{param_name}' is too large: {value} > {max_value}"
                )
        
        if allowed_values is not None and value not in allowed_values:
            raise InvalidParameterError(
                f"Parameter '{param_name}' has invalid value: {value}. Allowed values: {allowed_values}"
            )
