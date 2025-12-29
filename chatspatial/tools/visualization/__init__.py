"""
Visualization module for spatial transcriptomics.

This module provides visualization functions organized by analysis type:
- basic: Spatial plots, UMAP, heatmaps, violin plots, dotplots
- deconvolution: Cell type proportion visualizations
- cell_comm: Cell-cell communication visualizations
- velocity: RNA velocity visualizations
- trajectory: Trajectory and pseudotime visualizations
- spatial_stats: Spatial statistics visualizations
- enrichment: Pathway enrichment visualizations
- cnv: Copy number variation visualizations
- integration: Batch integration quality visualizations
- persistence: Visualization saving and export

Usage:
    from chatspatial.tools.visualization import (
        create_spatial_visualization,
        create_umap_visualization,
        create_deconvolution_visualization,
        # ... etc
    )
"""

# Basic visualizations
from .basic import (
    create_dotplot_visualization,
    create_heatmap_visualization,
    create_spatial_visualization,
    create_umap_visualization,
    create_violin_visualization,
)

# Cell communication visualizations
from .cell_comm import create_cell_communication_visualization

# CNV visualizations
from .cnv import create_cnv_heatmap_visualization, create_spatial_cnv_visualization

# Core utilities and data classes
from .core import (
    CellCommunicationData,
    DeconvolutionData,
    add_colorbar,
    create_figure,
    get_colormap,
    get_diverging_colormap,
    get_validated_features,
    plot_spatial_feature,
    setup_multi_panel_figure,
    validate_and_prepare_feature,
)

# CARD imputation (from deconvolution module)
# Deconvolution visualizations
from .deconvolution import (
    create_card_imputation_visualization,
    create_deconvolution_visualization,
)

# Enrichment visualizations
from .enrichment import (
    create_enrichment_visualization,
    create_pathway_enrichment_visualization,
)

# Batch integration visualizations
from .integration import create_batch_integration_visualization

# Main entry point and handler registry (from main.py to avoid circular imports)
from .main import PLOT_HANDLERS, visualize_data

# Multi-gene visualizations
from .multi_gene import (
    create_gene_correlation_visualization,
    create_lr_pairs_visualization,
    create_multi_gene_umap_visualization,
    create_multi_gene_visualization,
    create_spatial_interaction_visualization,
)

# Persistence functions
from .persistence import (
    clear_visualization_cache,
    export_all_visualizations,
    save_visualization,
)

# Spatial statistics visualizations
from .spatial_stats import create_spatial_statistics_visualization

# Trajectory visualizations
from .trajectory import create_trajectory_visualization

# RNA velocity visualizations
from .velocity import create_rna_velocity_visualization

__all__ = [
    # Core utilities
    "create_figure",
    "setup_multi_panel_figure",
    "add_colorbar",
    "get_colormap",
    "get_diverging_colormap",
    "plot_spatial_feature",
    "get_validated_features",
    "validate_and_prepare_feature",
    # Data classes
    "DeconvolutionData",
    "CellCommunicationData",
    # Basic visualizations
    "create_spatial_visualization",
    "create_umap_visualization",
    "create_heatmap_visualization",
    "create_violin_visualization",
    "create_dotplot_visualization",
    # Specialized visualizations
    "create_deconvolution_visualization",
    "create_cell_communication_visualization",
    "create_rna_velocity_visualization",
    "create_trajectory_visualization",
    "create_spatial_statistics_visualization",
    "create_enrichment_visualization",
    "create_pathway_enrichment_visualization",
    # CNV visualizations
    "create_card_imputation_visualization",
    "create_spatial_cnv_visualization",
    "create_cnv_heatmap_visualization",
    # Integration visualizations
    "create_batch_integration_visualization",
    # Multi-gene visualizations
    "create_multi_gene_visualization",
    "create_multi_gene_umap_visualization",
    "create_lr_pairs_visualization",
    "create_gene_correlation_visualization",
    "create_spatial_interaction_visualization",
    # Persistence functions
    "save_visualization",
    "export_all_visualizations",
    "clear_visualization_cache",
    # Main entry point
    "visualize_data",
    # Handler registry
    "PLOT_HANDLERS",
]
