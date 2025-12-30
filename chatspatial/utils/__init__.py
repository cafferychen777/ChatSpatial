"""
Utility functions for spatial transcriptomics data analysis.
"""

from .adata_utils import (  # Constants; Field discovery; Data access; Validation; Ensure; Standardization
    ALTERNATIVE_BATCH_KEYS,
    ALTERNATIVE_CELL_TYPE_KEYS,
    ALTERNATIVE_CLUSTER_KEYS,
    ALTERNATIVE_SPATIAL_KEYS,
    BATCH_KEY,
    CELL_TYPE_KEY,
    CLUSTER_KEY,
    SPATIAL_KEY,
    ensure_categorical,
    ensure_counts_layer,
    ensure_spatial_key,
    find_common_genes,
    get_batch_key,
    get_cell_type_key,
    get_cluster_key,
    get_gene_expression,
    get_genes_expression,
    get_spatial_coordinates,
    get_spatial_key,
    standardize_adata,
    to_dense,
    validate_adata,
    validate_adata_basics,
    validate_gene_overlap,
    validate_obs_column,
    validate_var_column,
)
from .dependency_manager import (
    DependencyCategory,
    DependencyInfo,
    DependencyManager,
    get,
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
    # Expression extraction
    "to_dense",
    "get_gene_expression",
    "get_genes_expression",
    # Validation
    "validate_adata",
    "validate_obs_column",
    "validate_var_column",
    "validate_adata_basics",
    "validate_gene_overlap",
    "ensure_categorical",
    # Gene overlap
    "find_common_genes",
    # Ensure
    "ensure_counts_layer",
    "ensure_spatial_key",
    # Standardization
    "standardize_adata",
    # Dependency management
    "DependencyManager",
    "DependencyInfo",
    "DependencyCategory",
    "require",
    "get",
    "is_available",
    "validate_r_environment",
    "validate_scvi_tools",
]
