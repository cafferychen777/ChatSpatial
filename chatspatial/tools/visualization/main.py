"""
Main visualization entry point.

This module contains the main visualize_data function that dispatches
to appropriate visualization handlers based on plot_type.

Refactored architecture (11 unified plot_types):
- feature: Spatial/UMAP feature visualization (basis='spatial'|'umap')
- expression: Aggregated expression (subtype='heatmap'|'violin'|'dotplot'|'correlation')
- deconvolution: Cell type proportions (subtype='spatial_multi'|'pie'|'dominant'|'diversity'|'umap'|'imputation')
- communication: Cell-cell communication (subtype='dotplot'|'tileplot'|'circle_plot')
- interaction: Spatial ligand-receptor pairs
- trajectory: Pseudotime and fate analysis (subtype='pseudotime'|'circular'|'fate_map'|'gene_trends'|'fate_heatmap'|'palantir')
- velocity: RNA velocity visualization (subtype='stream'|'phase'|'proportions'|'heatmap'|'paga')
- statistics: Spatial statistics (subtype='neighborhood'|'co_occurrence'|'ripley'|'moran'|'centrality'|'getis_ord')
- enrichment: Pathway enrichment (subtype='barplot'|'dotplot')
- cnv: Copy number variation (subtype='heatmap'|'spatial')
- integration: Batch integration quality (subtype='batch'|'cluster'|'highlight')
"""

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import scanpy as sc

from ...models.data import VisualizationParameters
from ...utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ...utils.image_utils import optimize_fig_to_image_with_cache

# Import unified visualization handlers
from .cell_comm import (
    create_cell_communication_visualization as create_communication_visualization,
)
from .cnv import create_cnv_visualization
from .deconvolution import create_deconvolution_visualization
from .enrichment import (
    create_pathway_enrichment_visualization as create_enrichment_visualization,
)
from .expression import create_expression_visualization
from .feature import create_feature_visualization
from .integration import (
    create_batch_integration_visualization as create_integration_visualization,
)
from .multi_gene import (
    create_spatial_interaction_visualization as create_interaction_visualization,
)
from .spatial_stats import (
    create_spatial_statistics_visualization as create_statistics_visualization,
)
from .trajectory import create_trajectory_visualization
from .velocity import create_rna_velocity_visualization as create_velocity_visualization

if TYPE_CHECKING:
    from ...spatial_mcp_adapter import ToolContext


# Handler registry for dispatch - 11 unified plot_types
PLOT_HANDLERS = {
    # Core feature visualization (replaces spatial, umap, multi_gene, lr_pairs)
    "feature": create_feature_visualization,
    # Aggregated expression (replaces heatmap, violin, dotplot, gene_correlation)
    "expression": create_expression_visualization,
    # Deconvolution (includes card_imputation as subtype)
    "deconvolution": create_deconvolution_visualization,
    # Cell-cell communication
    "communication": create_communication_visualization,
    # Spatial ligand-receptor interaction
    "interaction": create_interaction_visualization,
    # Trajectory/pseudotime
    "trajectory": create_trajectory_visualization,
    # RNA velocity
    "velocity": create_velocity_visualization,
    # Spatial statistics
    "statistics": create_statistics_visualization,
    # Pathway enrichment
    "enrichment": create_enrichment_visualization,
    # CNV (replaces cnv_heatmap, spatial_cnv)
    "cnv": create_cnv_visualization,
    # Batch integration quality
    "integration": create_integration_visualization,
}


async def visualize_data(
    data_id: str,
    ctx: "ToolContext",
    params: VisualizationParameters = VisualizationParameters(),
) -> str:
    """Visualize spatial transcriptomics data.

    Args:
        data_id: Dataset ID
        ctx: ToolContext for unified data access and logging
        params: Visualization parameters

    Returns:
        str: Path to saved visualization file with metadata

    Raises:
        DataNotFoundError: If the dataset is not found
        ParameterError: If parameters are invalid
        DataCompatibilityError: If data is not compatible with the visualization
        ProcessingError: If processing fails
    """
    # Validate parameters - use PLOT_HANDLERS as single source of truth
    if params.plot_type not in PLOT_HANDLERS:
        raise ParameterError(
            f"Invalid plot_type: {params.plot_type}. "
            f"Must be one of {list(PLOT_HANDLERS)}"
        )

    try:
        # Retrieve the AnnData object via ToolContext
        adata = await ctx.get_adata(data_id)

        # Validate AnnData object - basic validation
        if adata.n_obs < 5:
            raise DataNotFoundError("Dataset has too few cells (minimum 5 required)")
        if adata.n_vars < 5:
            raise DataNotFoundError("Dataset has too few genes (minimum 5 required)")

        # Set matplotlib style for better visualizations
        sc.settings.set_figure_params(dpi=params.dpi or 100, facecolor="white")

        # Dispatch to appropriate handler (pass ToolContext, not raw MCP Context)
        handler = PLOT_HANDLERS[params.plot_type]
        fig = await handler(adata, params, ctx)

        # Generate plot_type_key with subtype if applicable (for cache consistency)
        subtype = params.subtype
        plot_type_key = f"{params.plot_type}_{subtype}" if subtype else params.plot_type

        # Use the optimized conversion function
        return await optimize_fig_to_image_with_cache(
            fig,
            params,
            ctx,
            data_id=data_id,
            plot_type=plot_type_key,
        )

    except Exception as e:
        # Make sure to close any open figures in case of error
        plt.close("all")

        # Re-raise known error types directly
        if isinstance(e, (DataNotFoundError, ParameterError, DataCompatibilityError)):
            raise

        # Wrap unknown errors in ProcessingError
        raise ProcessingError(
            f"Failed to create {params.plot_type} visualization: {e}"
        ) from e
