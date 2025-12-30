"""
CNV (Copy Number Variation) visualization functions.

This module contains:
- CNV heatmap visualization
- Spatial CNV projection visualization
"""

from typing import TYPE_CHECKING, Any, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ...models.data import VisualizationParameters
from ...utils.adata_utils import require_spatial_coords, validate_obs_column
from ...utils.dependency_manager import require
from ...utils.exceptions import DataNotFoundError
from .core import plot_spatial_feature

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# Spatial CNV Visualization
# =============================================================================


async def create_spatial_cnv_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create spatial CNV projection visualization.

    Uses the unified plot_spatial_feature() helper for cleaner code
    and consistent parameter handling.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure object

    Raises:
        DataNotFoundError: If spatial coordinates or CNV features not found
    """
    if context:
        await context.info("Creating spatial CNV projection visualization")

    # Validate spatial coordinates
    require_spatial_coords(adata)

    # Determine feature to visualize
    feature_to_plot = params.feature

    # Auto-detect CNV-related features if none specified
    if not feature_to_plot:
        if "numbat_clone" in adata.obs:
            feature_to_plot = "numbat_clone"
            if context:
                await context.info(
                    "No feature specified, using 'numbat_clone' (Numbat clone assignment)"
                )
        elif "cnv_score" in adata.obs:
            feature_to_plot = "cnv_score"
            if context:
                await context.info(
                    "No feature specified, using 'cnv_score' (CNV score)"
                )
        elif "numbat_p_cnv" in adata.obs:
            feature_to_plot = "numbat_p_cnv"
            if context:
                await context.info(
                    "No feature specified, using 'numbat_p_cnv' (Numbat CNV probability)"
                )
        else:
            error_msg = "No CNV features found. Run analyze_cnv() first."
            if context:
                await context.warning(error_msg)
            raise DataNotFoundError(error_msg)

    # Validate feature exists
    validate_obs_column(adata, feature_to_plot, "CNV feature")

    if context:
        await context.info(f"Visualizing {feature_to_plot} on spatial coordinates")

    # Override colormap default for CNV data (RdBu_r is better for CNV scores)
    if not params.colormap:
        params.colormap = (
            "RdBu_r"
            if not pd.api.types.is_categorical_dtype(adata.obs[feature_to_plot])
            else "tab20"
        )

    figsize = params.figure_size if params.figure_size else (10, 8)
    fig, ax = plt.subplots(figsize=figsize)

    # Use the enhanced plot_spatial_feature helper
    plot_spatial_feature(adata, ax, feature=feature_to_plot, params=params)

    if context:
        await context.info(f"Spatial CNV projection created for {feature_to_plot}")

    return fig


# =============================================================================
# CNV Heatmap Visualization
# =============================================================================


async def create_cnv_heatmap_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CNV heatmap visualization.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure object

    Raises:
        DataNotFoundError: If CNV data not found
        DataCompatibilityError: If infercnvpy not installed
    """
    if context:
        await context.info("Creating CNV heatmap visualization")

    # Auto-detect CNV data source (infercnvpy or Numbat)
    cnv_method = None

    if "X_cnv" in adata.obsm:
        cnv_method = "infercnvpy"
    elif "X_cnv_numbat" in adata.obsm:
        cnv_method = "numbat"
    else:
        error_msg = "CNV data not found in obsm. Run analyze_cnv() first."
        if context:
            await context.warning(error_msg)
        raise DataNotFoundError(error_msg)

    if context:
        await context.info(f"Detected CNV data from {cnv_method} method")

    # Check if infercnvpy is available (needed for visualization)
    require("infercnvpy", feature="CNV heatmap visualization")

    # For Numbat data, temporarily copy to X_cnv for visualization
    if cnv_method == "numbat":
        if context:
            await context.info(
                "Converting Numbat CNV data to infercnvpy format for visualization"
            )
        adata.obsm["X_cnv"] = adata.obsm["X_cnv_numbat"]
        # Also ensure cnv metadata exists for infercnvpy plotting
        if "cnv" not in adata.uns:
            adata.uns["cnv"] = {
                "genomic_positions": False,
            }
            if context:
                await context.info(
                    "Note: Chromosome labels not available for Numbat heatmap. "
                    "Install R packages for full chromosome annotation."
                )

    # Check if CNV metadata exists
    if "cnv" not in adata.uns:
        error_msg = (
            "CNV metadata not found in adata.uns['cnv']. "
            "The CNV analysis may not have completed properly. "
            "Please re-run analyze_cnv()."
        )
        if context:
            await context.warning(error_msg)
        raise DataNotFoundError(error_msg)

    # Create CNV heatmap
    if context:
        await context.info("Generating CNV heatmap...")

    figsize = params.figure_size if params.figure_size else (12, 8)

    # For Numbat data without chromosome info, use aggregated heatmap by group
    if cnv_method == "numbat" and "chromosome" not in adata.var.columns:
        if context:
            await context.info(
                "Creating aggregated CNV heatmap by group (chromosome positions not available)"
            )

        # Get CNV matrix
        cnv_matrix = adata.obsm["X_cnv"]

        # Aggregate by feature (e.g., clone) for cleaner visualization
        if params.feature and params.feature in adata.obs.columns:
            # Group cells by feature and compute mean CNV per group
            feature_values = adata.obs[params.feature]
            unique_groups = sorted(feature_values.unique())

            # Compute mean CNV for each group
            aggregated_cnv_list: List[Any] = []
            group_labels: List[str] = []
            group_sizes: List[Any] = []

            for group in unique_groups:
                group_mask = feature_values == group
                group_cnv = cnv_matrix[group_mask, :].mean(axis=0)
                aggregated_cnv_list.append(group_cnv)
                group_labels.append(str(group))
                group_sizes.append(group_mask.sum())

            aggregated_cnv = np.array(aggregated_cnv_list)

            # Calculate appropriate figure width based on number of bins
            n_bins = aggregated_cnv.shape[1]
            fig_width = min(max(6, n_bins * 0.004), 12)
            fig_height = max(4, len(unique_groups) * 1.2)

            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

            # Plot aggregated heatmap with fixed aspect ratio
            im = ax.imshow(
                aggregated_cnv,
                cmap="RdBu_r",
                aspect="auto",
                vmin=-1,
                vmax=1,
                interpolation="nearest",
            )

            # Add colorbar
            plt.colorbar(im, ax=ax, label="Mean CNV state")

            # Set y-axis labels with group names and cell counts
            ax.set_yticks(range(len(group_labels)))
            ax.set_yticklabels(
                [
                    f"{label} (n={size})"
                    for label, size in zip(group_labels, group_sizes, strict=False)
                ]
            )
            feature_label = (
                params.feature
                if isinstance(params.feature, str)
                else ", ".join(params.feature) if params.feature else ""
            )
            ax.set_ylabel(feature_label, fontsize=12, fontweight="bold")

            # Set x-axis
            ax.set_xlabel("Genomic position (binned)", fontsize=12)
            ax.set_xticks([])  # Hide x-axis ticks for cleaner look

            # Add title
            ax.set_title(
                f"CNV Profile by {params.feature}\n(Numbat analysis, aggregated by group)",
                fontsize=14,
                fontweight="bold",
            )

            # Add gridlines between groups
            for i in range(len(group_labels) + 1):
                ax.axhline(i - 0.5, color="white", linewidth=2)

        else:
            # No grouping - show warning and plot all cells (not recommended)
            fig, ax = plt.subplots(figsize=figsize)

            sns.heatmap(
                cnv_matrix,
                cmap="RdBu_r",
                center=0,
                cbar_kws={"label": "CNV state"},
                yticklabels=False,
                xticklabels=False,
                ax=ax,
                vmin=-1,
                vmax=1,
            )

            ax.set_xlabel("Genomic position (binned)")
            ax.set_ylabel("Cells")
            ax.set_title("CNV Heatmap (Numbat)\nAll cells (ungrouped)")

        plt.tight_layout()

    else:
        # Use infercnvpy chromosome_heatmap for infercnvpy data or Numbat with chr info
        import infercnvpy as cnv

        if context:
            await context.info("Creating chromosome-organized CNV heatmap...")

        cnv.pl.chromosome_heatmap(
            adata,
            groupby=params.cluster_key,
            dendrogram=True,
            show=False,
            figsize=figsize,
        )
        # Get current figure
        fig = plt.gcf()

    if context:
        await context.info("CNV heatmap created successfully")

    return fig
