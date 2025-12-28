"""
Utility functions for spatial transcriptomics data analysis.
"""

from .adata_utils import (
    # Constants
    ALTERNATIVE_BATCH_KEYS,
    ALTERNATIVE_CELL_TYPE_KEYS,
    ALTERNATIVE_CLUSTER_KEYS,
    ALTERNATIVE_SPATIAL_KEYS,
    BATCH_KEY,
    CELL_TYPE_KEY,
    CLUSTER_KEY,
    SPATIAL_KEY,
    # Field discovery
    get_batch_key,
    get_cell_type_key,
    get_cluster_key,
    get_spatial_key,
    # Data access
    get_spatial_coordinates,
    # Validation
    ensure_categorical,
    validate_adata,
    validate_adata_basics,
    validate_obs_column,
    validate_obs_columns,
    validate_var_column,
    validate_var_columns,
    # Ensure
    ensure_spatial_key,
    # Standardization
    standardize_adata,
)
from .dependency_manager import (
    DependencyCategory,
    DependencyInfo,
    DependencyManager,
    get,
    get_manager,
    is_available,
    require,
    validate_r_environment,
    validate_scvi_tools,
)
from .exceptions import (
    ChatSpatialError,
    DataCompatibilityError,
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)
from .mcp_utils import mcp_tool_error_handler, suppress_output

__all__ = [
    # Exceptions
    "ChatSpatialError",
    "DataError",
    "DataNotFoundError",
    "DataCompatibilityError",
    "ParameterError",
    "ProcessingError",
    "DependencyError",
    # MCP utilities
    "suppress_output",
    "mcp_tool_error_handler",
    # Constants
    "SPATIAL_KEY",
    "CELL_TYPE_KEY",
    "CLUSTER_KEY",
    "BATCH_KEY",
    "ALTERNATIVE_SPATIAL_KEYS",
    "ALTERNATIVE_CELL_TYPE_KEYS",
    "ALTERNATIVE_CLUSTER_KEYS",
    "ALTERNATIVE_BATCH_KEYS",
    # Field discovery
    "get_batch_key",
    "get_cell_type_key",
    "get_cluster_key",
    "get_spatial_key",
    # Data access
    "get_spatial_coordinates",
    # Validation
    "validate_adata",
    "validate_obs_column",
    "validate_var_column",
    "validate_obs_columns",
    "validate_var_columns",
    "validate_adata_basics",
    "ensure_categorical",
    # Ensure
    "ensure_spatial_key",
    # Standardization
    "standardize_adata",
    # Dependency management
    "DependencyManager",
    "DependencyInfo",
    "DependencyCategory",
    "get_manager",
    "require",
    "get",
    "is_available",
    "validate_r_environment",
    "validate_scvi_tools",
]
