"""
Visualization module for spatial transcriptomics.

Refactored architecture with 11 unified plot_types:
- feature: Spatial/UMAP feature visualization (basis='spatial'|'umap')
- expression: Aggregated expression (subtype='heatmap'|'violin'|'dotplot'|'correlation')
- deconvolution: Cell type proportions (subtype='spatial_multi'|'pie'|'dominant'|'imputation')
- communication: Cell-cell communication patterns
- interaction: Spatial ligand-receptor pairs
- trajectory: Pseudotime and fate analysis
- velocity: RNA velocity visualization
- statistics: Spatial statistics (Moran's I, etc.)
- enrichment: Pathway/gene set enrichment
- cnv: Copy number variation (subtype='heatmap'|'spatial')
- integration: Batch integration quality

Usage:
    from chatspatial.tools.visualization import visualize_data, PLOT_HANDLERS
"""

from .cell_comm import create_cell_communication_visualization
from .cnv import create_cnv_visualization

# Core utilities and data classes
from .core import (
    FIGURE_DEFAULTS,
    CellCommunicationData,
    DeconvolutionData,
    add_colorbar,
    create_figure,
    create_figure_from_params,
    get_categorical_cmap,
    get_category_colors,
    get_colormap,
    get_diverging_colormap,
    get_validated_features,
    plot_spatial_feature,
    resolve_figure_size,
    setup_multi_panel_figure,
    validate_and_prepare_feature,
)
from .deconvolution import create_deconvolution_visualization
from .enrichment import create_pathway_enrichment_visualization
from .expression import create_expression_visualization

# Unified visualization handlers
from .feature import create_feature_visualization
from .integration import create_batch_integration_visualization

# Main entry point and handler registry
from .main import PLOT_HANDLERS, visualize_data
from .multi_gene import create_spatial_interaction_visualization
from .spatial_stats import create_spatial_statistics_visualization
from .trajectory import create_trajectory_visualization
from .velocity import create_rna_velocity_visualization

__all__ = [
    # Core utilities
    "FIGURE_DEFAULTS",
    "create_figure",
    "create_figure_from_params",
    "resolve_figure_size",
    "setup_multi_panel_figure",
    "add_colorbar",
    "get_colormap",
    "get_categorical_cmap",
    "get_category_colors",
    "get_diverging_colormap",
    "plot_spatial_feature",
    "get_validated_features",
    "validate_and_prepare_feature",
    # Data classes
    "DeconvolutionData",
    "CellCommunicationData",
    # Main entry point
    "visualize_data",
    # Handler registry
    "PLOT_HANDLERS",
    # Unified visualization handlers
    "create_feature_visualization",
    "create_expression_visualization",
    "create_deconvolution_visualization",
    "create_cell_communication_visualization",
    "create_spatial_interaction_visualization",
    "create_trajectory_visualization",
    "create_rna_velocity_visualization",
    "create_spatial_statistics_visualization",
    "create_pathway_enrichment_visualization",
    "create_cnv_visualization",
    "create_batch_integration_visualization",
]
