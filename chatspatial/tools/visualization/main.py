"""
Main visualization entry point.

This module contains the main visualize_data function that dispatches
to appropriate visualization handlers based on plot_type.
"""

import traceback
from typing import TYPE_CHECKING, Tuple, Union

import matplotlib.pyplot as plt
import scanpy as sc
from mcp.server.fastmcp.utilities.types import ImageContent
from mcp.types import EmbeddedResource

from ...models.data import VisualizationParameters
from ...utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ...utils.image_utils import optimize_fig_to_image_with_cache

# Import all visualization handlers
from .basic import (
    create_dotplot_visualization,
    create_heatmap_visualization,
    create_spatial_visualization,
    create_umap_visualization,
    create_violin_visualization,
)
from .cell_comm import create_cell_communication_visualization
from .cnv import create_cnv_heatmap_visualization, create_spatial_cnv_visualization
from .deconvolution import (
    create_card_imputation_visualization,
    create_deconvolution_visualization,
)
from .enrichment import create_pathway_enrichment_visualization
from .integration import create_batch_integration_visualization
from .multi_gene import (
    create_gene_correlation_visualization,
    create_lr_pairs_visualization,
    create_multi_gene_visualization,
    create_spatial_interaction_visualization,
)
from .spatial_stats import create_spatial_statistics_visualization
from .trajectory import create_trajectory_visualization
from .velocity import create_rna_velocity_visualization

if TYPE_CHECKING:
    from ...spatial_mcp_adapter import ToolContext


# Handler registry for dispatch - defined here to avoid circular imports
PLOT_HANDLERS = {
    # Basic plots
    "spatial": create_spatial_visualization,
    "umap": create_umap_visualization,
    "heatmap": create_heatmap_visualization,
    "violin": create_violin_visualization,
    "dotplot": create_dotplot_visualization,
    # Analysis-specific plots
    "deconvolution": create_deconvolution_visualization,
    "cell_communication": create_cell_communication_visualization,
    "rna_velocity": create_rna_velocity_visualization,
    "trajectory": create_trajectory_visualization,
    "spatial_statistics": create_spatial_statistics_visualization,
    "pathway_enrichment": create_pathway_enrichment_visualization,
    # CNV plots
    "card_imputation": create_card_imputation_visualization,
    "spatial_cnv": create_spatial_cnv_visualization,
    "cnv_heatmap": create_cnv_heatmap_visualization,
    # Integration plots
    "batch_integration": create_batch_integration_visualization,
    # Multi-gene plots
    "multi_gene": create_multi_gene_visualization,
    "lr_pairs": create_lr_pairs_visualization,
    "gene_correlation": create_gene_correlation_visualization,
    "spatial_interaction": create_spatial_interaction_visualization,
}


async def visualize_data(
    data_id: str,
    ctx: "ToolContext",
    params: VisualizationParameters = VisualizationParameters(),  # type: ignore[call-arg]
) -> Union[ImageContent, Tuple[ImageContent, EmbeddedResource]]:
    """Visualize spatial transcriptomics data.

    Args:
        data_id: Dataset ID
        ctx: ToolContext for unified data access and logging
        params: Visualization parameters

    Returns:
        Union[ImageContent, Tuple[ImageContent, EmbeddedResource]]:
            - Small images (<100KB): ImageContent object
            - Large images (>=100KB): Tuple[Preview ImageContent, High-quality Resource]

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
            f"Must be one of {list(PLOT_HANDLERS.keys())}"
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

        # Dispatch to appropriate handler
        handler = PLOT_HANDLERS[params.plot_type]
        fig = await handler(adata, params, ctx._mcp_context)

        # Generate plot_type_key with subtype if applicable (for cache consistency)
        subtype = params.subtype
        plot_type_key = f"{params.plot_type}_{subtype}" if subtype else params.plot_type

        # Use the optimized conversion function
        return await optimize_fig_to_image_with_cache(
            fig,
            params,
            ctx._mcp_context,
            data_id=data_id,
            plot_type=plot_type_key,
            mode="auto",
        )

    except Exception as e:
        # Make sure to close any open figures in case of error
        plt.close("all")

        # For image conversion errors, return error message as string
        if "fig_to_image" in str(e) or "convert" in str(e).lower():
            error_details = traceback.format_exc()
            return (
                f"Error in {params.plot_type} visualization:\n\n"
                f"{str(e)}\n\n"
                f"Technical details:\n{error_details}"
            )

        # Wrap the error in a more informative exception
        if isinstance(e, (DataNotFoundError, ParameterError, DataCompatibilityError)):
            raise
        else:
            raise ProcessingError(
                f"Failed to create {params.plot_type} visualization: {str(e)}"
            ) from e
