"""
Enhanced error handling for MCP server
"""

from typing import Any, Callable, Dict, Optional, Union, List, TYPE_CHECKING
from functools import wraps
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import traceback
import io
import logging
import warnings
from ..mcp.errors import ErrorType, format_mcp_error

if TYPE_CHECKING:
    import scanpy as sc
    import numpy as np


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


def validate_adata(
    adata: 'sc.AnnData', 
    required_keys: Dict[str, Union[str, List[str]]], 
    context: Optional[Any] = None,
    # New parameters for enhanced validation (backward compatible)
    check_spatial: bool = False,
    check_velocity: bool = False,
    spatial_key: str = 'spatial'
) -> None:
    """Validate AnnData object has required keys and optional data integrity checks
    
    Args:
        adata: AnnData object to validate
        required_keys: Dict of required keys by category (obs, var, obsm, etc.)
        context: Optional context for logging
        check_spatial: Whether to validate spatial coordinates (default: False)
        check_velocity: Whether to validate velocity data layers (default: False) 
        spatial_key: Key for spatial coordinates in adata.obsm (default: 'spatial')
    
    Raises:
        DataNotFoundError: If required keys are missing or validation fails
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
    
    # Enhanced validation checks (new functionality, backward compatible)
    if check_spatial:
        _validate_spatial_data_internal(adata, spatial_key, missing)
    
    if check_velocity:
        _validate_velocity_data_internal(adata, missing)
    
    # Raise combined errors if any validation failed
    if missing:
        raise DataNotFoundError(f"Validation failed: {', '.join(missing)}")


def _validate_spatial_data_internal(adata: 'sc.AnnData', spatial_key: str, issues: List[str]) -> None:
    """Internal helper for spatial data validation"""
    import numpy as np
    
    if spatial_key not in adata.obsm:
        issues.append(f"Missing '{spatial_key}' coordinates in adata.obsm")
        return
    
    spatial_coords = adata.obsm[spatial_key]
    
    # Check dimensions
    if spatial_coords.shape[1] < 2:
        issues.append(f"Spatial coordinates should have at least 2 dimensions, found {spatial_coords.shape[1]}")
    
    # Check for NaN values
    if np.any(np.isnan(spatial_coords)):
        issues.append("Spatial coordinates contain NaN values")
    
    # Check for constant coordinates
    if np.std(spatial_coords[:, 0]) == 0 and np.std(spatial_coords[:, 1]) == 0:
        issues.append("All spatial coordinates are identical")


def _validate_velocity_data_internal(adata: 'sc.AnnData', issues: List[str]) -> None:
    """Internal helper for velocity data validation"""
    import numpy as np
    
    if 'spliced' not in adata.layers:
        issues.append("Missing 'spliced' layer required for RNA velocity")
    if 'unspliced' not in adata.layers:
        issues.append("Missing 'unspliced' layer required for RNA velocity")
    
    # Check if layers have data
    if 'spliced' in adata.layers and 'unspliced' in adata.layers:
        for layer_name in ['spliced', 'unspliced']:
            layer_data = adata.layers[layer_name]
            
            # Check for empty or all-zero layers  
            if hasattr(layer_data, 'nnz'):  # Sparse matrix
                if layer_data.nnz == 0:
                    issues.append(f"'{layer_name}' layer is empty (all zeros)")
            else:  # Dense matrix
                if np.all(layer_data == 0):
                    issues.append(f"'{layer_name}' layer is empty (all zeros)")
            
            # Check for NaN values
            if hasattr(layer_data, 'data'):  # Sparse matrix
                if np.any(np.isnan(layer_data.data)):
                    issues.append(f"'{layer_name}' layer contains NaN values")
            else:  # Dense matrix
                if np.any(np.isnan(layer_data)):
                    issues.append(f"'{layer_name}' layer contains NaN values"


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




@contextmanager
def suppress_output():
    """Context manager to suppress stdout, stderr, warnings, and logging during analysis.
    
    This combines the best features of both previous implementations:
    - Suppresses warnings (from previous utils implementation)
    - Controls logging levels (from previous local implementation)  
    - Captures stdout and stderr output
    
    Usage:
        with suppress_output():
            # Code that produces unwanted output
            pass
    """
    # Save original logging level
    old_level = logging.getLogger().level
    logging.getLogger().setLevel(logging.ERROR)
    
    # Suppress warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        # Capture stdout and stderr
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        
        try:
            with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
                yield
        finally:
            # Restore original logging level
            logging.getLogger().setLevel(old_level)
            # Note: Captured output is intentionally discarded for cleaner analysis output