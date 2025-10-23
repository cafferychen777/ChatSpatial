"""
Utility functions for spatial transcriptomics data analysis.
"""

from .data_adapter import (
    ensure_spatial_neighbors,
    get_cell_types,
    get_clusters,
    get_spatial_coordinates,
    standardize_adata,
)
from .data_validator import (
    raise_if_invalid,
    validate_for_cell_communication,
    validate_for_deconvolution,
    validate_for_spatial_statistics,
    validate_spatial_data,
)
from .error_handling import ProcessingError, suppress_output
from .tool_error_handling import (
    ToolResult,
    create_error_result,
    mcp_tool_error_handler,
)

__all__ = [
    "ProcessingError",
    "suppress_output",
    "ToolResult",
    "create_error_result",
    "mcp_tool_error_handler",
    "standardize_adata",
    "get_spatial_coordinates",
    "get_cell_types",
    "get_clusters",
    "ensure_spatial_neighbors",
    "validate_spatial_data",
    "validate_for_cell_communication",
    "validate_for_deconvolution",
    "validate_for_spatial_statistics",
    "raise_if_invalid",
]
