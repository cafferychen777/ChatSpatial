"""
Basic visualization functions for spatial transcriptomics.

This module contains:
- Spatial feature plots
- UMAP visualizations
- Heatmap visualizations
- Violin plots
- Dot plots
"""

from typing import TYPE_CHECKING, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...utils.adata_utils import (
    ensure_categorical,
    get_gene_expression,
    validate_obs_column,
)
from ...utils.compute import ensure_umap
from ...utils.exceptions import DataNotFoundError, ParameterError
from .core import (
    add_colorbar,
    create_figure,
    get_colormap,
    plot_spatial_feature,
    setup_multi_panel_figure,
)

# =============================================================================
# Spatial Visualization
# =============================================================================


async def create_spatial_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create spatial visualization for one or more features.

    Args:
        adata: AnnData object with spatial coordinates
        params: Visualization parameters
        context: Optional tool context for logging

    Returns:
        matplotlib Figure object
    """
    if params.feature is None:
        features: List[str] = []
    elif isinstance(params.feature, list):
        features = params.feature
    else:
        features = [params.feature]

    if not features:
        # Default to leiden clustering if available
        if "leiden" in adata.obs.columns:
            features = ["leiden"]
        else:
            raise ParameterError(
                "No features specified and no default clustering found"
            )

    n_features = len(features)
    fig, axes = setup_multi_panel_figure(n_features, params, "")

    for i, feature in enumerate(features):
        ax = axes[i]
        try:
            mappable = plot_spatial_feature(
                adata,
                ax,
                feature=feature,
                params=params,
                title=feature,
            )
            if mappable is not None and params.show_colorbar:
                add_colorbar(fig, ax, mappable, params)
        except Exception as e:
            ax.text(
                0.5,
                0.5,
                f"Error: {str(e)[:50]}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")

    plt.tight_layout()
    return fig


# =============================================================================
# UMAP Visualization
# =============================================================================


async def create_umap_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create UMAP visualization.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: Optional tool context for logging

    Returns:
        matplotlib Figure object
    """
    if context:
        feature_desc = params.feature if params.feature else "clusters"
        await context.info(f"Creating UMAP plot for {feature_desc}")

    # Ensure UMAP is computed (lazy computation)
    if ensure_umap(adata):
        if context:
            await context.info("Computed UMAP embedding")

    # Determine what to color by
    color_by = params.feature
    if color_by is None:
        # Default to leiden if available
        if "leiden" in adata.obs.columns:
            color_by = "leiden"
        elif "louvain" in adata.obs.columns:
            color_by = "louvain"

    # Create figure
    fig, ax = create_figure(params.figure_size or (10, 8))

    # Get UMAP coordinates
    umap_coords = adata.obsm["X_umap"]

    # Get color values
    spot_size = params.spot_size if params.spot_size is not None else 150.0
    if color_by is None:
        # No color - just plot points
        ax.scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            s=spot_size // 3,
            alpha=params.alpha,
            c="steelblue",
        )
    elif color_by in adata.var_names:
        # Gene expression - use unified utility
        values = get_gene_expression(adata, color_by)

        scatter = ax.scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c=values,
            cmap=params.colormap,
            s=spot_size // 3,
            alpha=params.alpha,
        )
        if params.show_colorbar:
            add_colorbar(fig, ax, scatter, params, label=color_by)
    elif color_by in adata.obs.columns:
        # Observation column
        values = adata.obs[color_by]
        is_categorical = (
            pd.api.types.is_categorical_dtype(values) or values.dtype == object
        )

        if is_categorical:
            ensure_categorical(adata, color_by)
            categories = adata.obs[color_by].cat.categories
            n_cats = len(categories)
            colors = get_colormap(params.colormap, n_colors=n_cats)

            for i, cat in enumerate(categories):
                mask = adata.obs[color_by] == cat
                ax.scatter(
                    umap_coords[mask, 0],
                    umap_coords[mask, 1],
                    c=[colors[i]],
                    s=spot_size // 3,
                    alpha=params.alpha,
                    label=cat,
                )

            if params.show_legend:
                ax.legend(
                    loc="center left",
                    bbox_to_anchor=(1, 0.5),
                    fontsize=8,
                    frameon=False,
                )
        else:
            scatter = ax.scatter(
                umap_coords[:, 0],
                umap_coords[:, 1],
                c=values,
                cmap=params.colormap,
                s=spot_size // 3,
                alpha=params.alpha,
            )
            if params.show_colorbar:
                add_colorbar(fig, ax, scatter, params, label=color_by)
    else:
        raise DataNotFoundError(f"Color feature '{color_by}' not found in genes or obs")

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(params.title or f"UMAP - {color_by}")

    if not params.show_axes:
        ax.axis("off")

    plt.tight_layout()
    return fig


# =============================================================================
# Heatmap Visualization
# =============================================================================


async def create_heatmap_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create heatmap visualization for gene expression.

    Args:
        adata: AnnData object
        params: Visualization parameters (requires cluster_key and feature list)
        context: Optional tool context for logging

    Returns:
        matplotlib Figure object
    """
    if not params.cluster_key:
        raise ParameterError("Heatmap requires cluster_key parameter")

    validate_obs_column(adata, params.cluster_key, "cluster_key")

    if params.feature is None:
        features_raw: List[str] = []
    elif isinstance(params.feature, list):
        features_raw = params.feature
    else:
        features_raw = [params.feature]
    features = [f for f in features_raw if f in adata.var_names]

    if not features:
        raise ParameterError("No valid gene features provided for heatmap")

    if context:
        await context.info(
            f"Creating heatmap for {len(features)} genes grouped by {params.cluster_key}"
        )

    # Use scanpy's heatmap function
    # Note: return_fig=True causes issues with newer matplotlib versions
    sc.pl.heatmap(
        adata,
        var_names=features,
        groupby=params.cluster_key,
        cmap=params.colormap,
        show=False,
        dendrogram=params.dotplot_dendrogram,
        swap_axes=params.dotplot_swap_axes,
        standard_scale=params.dotplot_standard_scale,
    )
    fig = plt.gcf()

    return fig


# =============================================================================
# Violin Plot Visualization
# =============================================================================


async def create_violin_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create violin plot visualization for gene expression.

    Args:
        adata: AnnData object
        params: Visualization parameters (requires cluster_key and feature)
        context: Optional tool context for logging

    Returns:
        matplotlib Figure object
    """
    if not params.cluster_key:
        raise ParameterError("Violin plot requires cluster_key parameter")

    validate_obs_column(adata, params.cluster_key, "cluster_key")

    if params.feature is None:
        features_raw: List[str] = []
    elif isinstance(params.feature, list):
        features_raw = params.feature
    else:
        features_raw = [params.feature]
    features = [f for f in features_raw if f in adata.var_names]

    if not features:
        raise ParameterError("No valid gene features provided for violin plot")

    if context:
        await context.info(
            f"Creating violin plot for {len(features)} genes grouped by {params.cluster_key}"
        )

    # Use scanpy's violin function
    # Note: return_fig=True causes issues with newer matplotlib/seaborn versions
    # Instead, we use plt.gcf() to get the current figure
    sc.pl.violin(
        adata,
        keys=features,
        groupby=params.cluster_key,
        show=False,
    )
    fig = plt.gcf()

    return fig


# =============================================================================
# Dot Plot Visualization
# =============================================================================


async def create_dotplot_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create dot plot visualization for gene expression.

    Args:
        adata: AnnData object
        params: Visualization parameters (requires cluster_key and feature list)
        context: Optional tool context for logging

    Returns:
        matplotlib Figure object
    """
    if not params.cluster_key:
        raise ParameterError("Dot plot requires cluster_key parameter")

    validate_obs_column(adata, params.cluster_key, "cluster_key")

    if params.feature is None:
        features_raw: List[str] = []
    elif isinstance(params.feature, list):
        features_raw = params.feature
    else:
        features_raw = [params.feature]
    features = [f for f in features_raw if f in adata.var_names]

    if not features:
        raise ParameterError("No valid gene features provided for dot plot")

    if context:
        await context.info(
            f"Creating dot plot for {len(features)} genes grouped by {params.cluster_key}"
        )

    # Build kwargs for dotplot
    # Note: return_fig=True returns DotPlot object, not Figure
    # Use plt.gcf() instead to get the figure
    dotplot_kwargs = {
        "adata": adata,
        "var_names": features,
        "groupby": params.cluster_key,
        "cmap": params.colormap,
        "show": False,
    }

    # Add optional parameters
    if params.dotplot_dendrogram:
        dotplot_kwargs["dendrogram"] = True
    if params.dotplot_swap_axes:
        dotplot_kwargs["swap_axes"] = True
    if params.dotplot_standard_scale:
        dotplot_kwargs["standard_scale"] = params.dotplot_standard_scale
    if params.dotplot_dot_min is not None:
        dotplot_kwargs["dot_min"] = params.dotplot_dot_min
    if params.dotplot_dot_max is not None:
        dotplot_kwargs["dot_max"] = params.dotplot_dot_max
    if params.dotplot_smallest_dot is not None:
        dotplot_kwargs["smallest_dot"] = params.dotplot_smallest_dot
    if params.dotplot_var_groups:
        dotplot_kwargs["var_group_positions"] = list(params.dotplot_var_groups.keys())
        dotplot_kwargs["var_group_labels"] = list(params.dotplot_var_groups.keys())

    sc.pl.dotplot(**dotplot_kwargs)
    fig = plt.gcf()

    return fig
