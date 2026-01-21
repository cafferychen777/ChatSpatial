"""
Unified expression visualization for spatial transcriptomics.

This module provides aggregated gene expression visualizations
grouped by cell clusters or cell types.

Replaces: heatmap, violin, dotplot, gene_correlation (from basic.py and multi_gene.py)
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from ...models.data import VisualizationParameters
from ...utils.adata_utils import get_genes_expression, validate_obs_column
from ...utils.exceptions import ParameterError
from .core import get_validated_features

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# Unified Expression Visualization
# =============================================================================


async def create_expression_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Unified expression visualization with multiple subtypes.

    Subtype selection via params.subtype:
        - "heatmap" (default): Clustered heatmap of gene expression
        - "violin": Violin plots grouped by cluster
        - "dotplot": Dot plot showing expression and percentage
        - "correlation": Gene-gene correlation matrix

    All subtypes except "correlation" require params.cluster_key.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: Optional tool context for logging

    Returns:
        matplotlib Figure

    Raises:
        ParameterError: If required parameters are missing
    """
    subtype = params.subtype or "heatmap"

    if context:
        await context.info(f"Creating expression visualization: {subtype}")

    if subtype == "heatmap":
        return await _create_heatmap(adata, params, context)
    elif subtype == "violin":
        return await _create_violin(adata, params, context)
    elif subtype == "dotplot":
        return await _create_dotplot(adata, params, context)
    elif subtype == "correlation":
        return await _create_correlation(adata, params, context)
    else:
        raise ParameterError(
            f"Invalid expression subtype: {subtype}. "
            "Use 'heatmap', 'violin', 'dotplot', or 'correlation'."
        )


# =============================================================================
# Heatmap Visualization
# =============================================================================


async def _create_heatmap(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"],
) -> plt.Figure:
    """Create heatmap visualization for gene expression.

    Args:
        adata: AnnData object
        params: Visualization parameters (requires cluster_key and feature list)
        context: Optional tool context for logging

    Returns:
        matplotlib Figure
    """
    if not params.cluster_key:
        raise ParameterError(
            "Heatmap requires cluster_key parameter. "
            "Example: plot_type='expression', subtype='heatmap', cluster_key='leiden'"
        )

    validate_obs_column(adata, params.cluster_key, "Cluster")

    features = await get_validated_features(adata, params, context, genes_only=True)
    if not features:
        raise ParameterError("No valid gene features provided for heatmap")

    if context:
        await context.info(
            f"Creating heatmap for {len(features)} genes grouped by {params.cluster_key}"
        )

    # Use scanpy's heatmap function
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


async def _create_violin(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"],
) -> plt.Figure:
    """Create violin plot visualization for gene expression.

    Args:
        adata: AnnData object
        params: Visualization parameters (requires cluster_key and feature)
        context: Optional tool context for logging

    Returns:
        matplotlib Figure
    """
    if not params.cluster_key:
        raise ParameterError(
            "Violin plot requires cluster_key parameter. "
            "Example: plot_type='expression', subtype='violin', cluster_key='leiden'"
        )

    validate_obs_column(adata, params.cluster_key, "Cluster")

    features = await get_validated_features(adata, params, context, genes_only=True)
    if not features:
        raise ParameterError("No valid gene features provided for violin plot")

    if context:
        await context.info(
            f"Creating violin plot for {len(features)} genes grouped by {params.cluster_key}"
        )

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


async def _create_dotplot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"],
) -> plt.Figure:
    """Create dot plot visualization for gene expression.

    Args:
        adata: AnnData object
        params: Visualization parameters (requires cluster_key and feature list)
        context: Optional tool context for logging

    Returns:
        matplotlib Figure
    """
    if not params.cluster_key:
        raise ParameterError(
            "Dot plot requires cluster_key parameter. "
            "Example: plot_type='expression', subtype='dotplot', cluster_key='leiden'"
        )

    validate_obs_column(adata, params.cluster_key, "Cluster")

    features = await get_validated_features(adata, params, context, genes_only=True)
    if not features:
        raise ParameterError("No valid gene features provided for dot plot")

    if context:
        await context.info(
            f"Creating dot plot for {len(features)} genes grouped by {params.cluster_key}"
        )

    # Build kwargs for dotplot
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


# =============================================================================
# Correlation Visualization
# =============================================================================


async def _create_correlation(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"],
) -> plt.Figure:
    """Create gene correlation visualization using seaborn clustermap.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: Optional tool context for logging

    Returns:
        matplotlib Figure with gene correlation clustermap
    """
    # Get validated genes
    available_genes = await get_validated_features(
        adata, params, max_features=10, context=context
    )

    if context:
        await context.info(
            f"Creating gene correlation visualization for {len(available_genes)} genes"
        )

    # Get expression matrix
    expr_matrix = get_genes_expression(adata, available_genes)

    # Apply color scaling
    if params.color_scale == "log":
        expr_matrix = np.log1p(expr_matrix)
    elif params.color_scale == "sqrt":
        expr_matrix = np.sqrt(expr_matrix)

    # Create DataFrame for correlation
    expr_df = pd.DataFrame(expr_matrix, columns=available_genes)

    # Calculate correlation
    corr_df = expr_df.corr(method=params.correlation_method)

    # Use seaborn clustermap
    figsize = params.figure_size or (
        max(8, len(available_genes)),
        max(8, len(available_genes)),
    )

    g = sns.clustermap(
        corr_df,
        cmap=params.colormap,
        center=0,
        annot=True,
        fmt=".2f",
        square=True,
        figsize=figsize,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
    )

    title = params.title or f"Gene Correlation ({params.correlation_method.title()})"
    g.fig.suptitle(title, y=1.02, fontsize=14)

    return g.fig
