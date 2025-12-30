"""
RNA velocity visualization functions for spatial transcriptomics.

This module contains:
- Velocity stream plots
- Phase plots (spliced vs unspliced)
- Proportions plots (pie charts)
- Velocity heatmaps
- PAGA with velocity arrows
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...utils.adata_utils import validate_obs_column
from ...utils.dependency_manager import require
from ...utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
)
from .core import get_categorical_columns, infer_basis

# =============================================================================
# Main Router
# =============================================================================


async def create_rna_velocity_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create RNA velocity visualization based on subtype.

    Dispatcher function that routes to appropriate scVelo visualization.

    Args:
        adata: AnnData object with computed RNA velocity
        params: Visualization parameters including subtype
        context: MCP context

    Returns:
        Matplotlib figure with RNA velocity visualization

    Subtypes:
        - stream (default): Velocity embedding stream plot
        - phase: Phase plot showing spliced vs unspliced
        - proportions: Pie chart of spliced/unspliced ratios
        - heatmap: Gene expression ordered by latent_time
        - paga: PAGA with velocity arrows
    """
    subtype = params.subtype or "stream"

    if context:
        await context.info(f"Creating RNA velocity visualization (subtype: {subtype})")

    if subtype == "stream":
        return await _create_velocity_stream_plot(adata, params, context)
    elif subtype == "phase":
        return await _create_velocity_phase_plot(adata, params, context)
    elif subtype == "proportions":
        return await _create_velocity_proportions_plot(adata, params, context)
    elif subtype == "heatmap":
        return await _create_velocity_heatmap(adata, params, context)
    elif subtype == "paga":
        return await _create_velocity_paga_plot(adata, params, context)
    else:
        raise ParameterError(
            f"Unsupported subtype for rna_velocity: '{subtype}'. "
            f"Available subtypes: stream, phase, proportions, heatmap, paga"
        )


# =============================================================================
# Visualization Functions
# =============================================================================


async def _create_velocity_stream_plot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create RNA velocity stream plot using scv.pl.velocity_embedding_stream.

    Data requirements:
        - adata.uns['velocity_graph']: Velocity transition graph
        - adata.obsm['X_umap'] or 'spatial': Embedding for visualization
    """
    require("scvelo", feature="RNA velocity visualization")
    import scvelo as scv

    if "velocity_graph" not in adata.uns:
        raise DataNotFoundError(
            "RNA velocity not computed. Run analyze_velocity_data first."
        )

    # Determine basis for plotting
    basis = infer_basis(adata, preferred=params.basis)
    if not basis:
        raise DataCompatibilityError(
            f"No valid embedding basis found. "
            f"Available keys: {list(adata.obsm.keys())}"
        )
    if context and basis != params.basis:
        await context.info(f"Using '{basis}' as basis")

    # Prepare feature for coloring
    feature = params.feature
    if not feature:
        categorical_cols = get_categorical_columns(adata)
        feature = categorical_cols[0] if categorical_cols else None
        if feature and context:
            await context.info(f"Using '{feature}' for coloring")

    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    scv.pl.velocity_embedding_stream(
        adata,
        basis=basis,
        color=feature,
        ax=ax,
        show=False,
        alpha=params.alpha,
        legend_loc="right margin" if feature and feature in adata.obs.columns else None,
        frameon=params.show_axes,
        title="",
    )

    title = params.title or f"RNA Velocity Stream on {basis.capitalize()}"
    ax.set_title(title, fontsize=14)

    if basis == "spatial":
        ax.invert_yaxis()

    plt.tight_layout()
    return fig


async def _create_velocity_phase_plot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create velocity phase plot using scv.pl.velocity.

    Shows spliced vs unspliced counts with fitted velocity model for specified genes.

    Data requirements:
        - adata.layers['velocity']: Velocity vectors
        - adata.layers['Ms']: Smoothed spliced counts
        - adata.layers['Mu']: Smoothed unspliced counts
    """
    require("scvelo", feature="velocity phase plots")
    import scvelo as scv

    required_layers = ["velocity", "Ms", "Mu"]
    missing_layers = [layer for layer in required_layers if layer not in adata.layers]
    if missing_layers:
        raise DataNotFoundError(
            f"Missing layers for phase plot: {missing_layers}. Run velocity analysis first."
        )

    if params.feature:
        if isinstance(params.feature, str):
            var_names = [params.feature]
        else:
            var_names = list(params.feature)
    else:
        if "velocity_genes" in adata.var.columns:
            velocity_genes = adata.var_names[adata.var["velocity_genes"]]
            var_names = list(velocity_genes[:4])
        else:
            var_names = list(adata.var_names[:4])

    valid_genes = [g for g in var_names if g in adata.var_names]
    if not valid_genes:
        raise DataNotFoundError(
            f"None of the specified genes found in data: {var_names}. "
            f"Available genes (first 10): {list(adata.var_names[:10])}"
        )

    if context:
        await context.info(f"Creating phase plot for genes: {valid_genes}")

    basis = infer_basis(adata, preferred=params.basis, priority=["umap", "spatial"])
    figsize = params.figure_size or (4 * len(valid_genes), 4)
    color = params.cluster_key if params.cluster_key else None

    scv.pl.velocity(
        adata,
        var_names=valid_genes,
        basis=basis,
        color=color,
        figsize=figsize,
        dpi=params.dpi,
        show=False,
        ncols=len(valid_genes),
    )

    fig = plt.gcf()
    title = params.title or "RNA Velocity Phase Plot"
    fig.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


async def _create_velocity_proportions_plot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create velocity proportions plot using scv.pl.proportions.

    Shows pie chart of spliced/unspliced RNA proportions per cluster.

    Data requirements:
        - adata.layers['spliced']: Spliced counts
        - adata.layers['unspliced']: Unspliced counts
        - adata.obs[cluster_key]: Cluster labels for grouping
    """
    require("scvelo", feature="proportions plot")
    import scvelo as scv

    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise DataNotFoundError(
            "Spliced and unspliced layers are required for proportions plot. "
            "Your data may not contain RNA velocity information."
        )

    cluster_key = params.cluster_key
    if not cluster_key:
        categorical_cols = get_categorical_columns(adata)
        if categorical_cols:
            cluster_key = categorical_cols[0]
            if context:
                await context.info(f"Using cluster_key: '{cluster_key}'")
        else:
            raise ParameterError(
                "cluster_key is required for proportions plot. "
                f"Available columns: {list(adata.obs.columns)[:10]}"
            )

    validate_obs_column(adata, cluster_key, f"Cluster key '{cluster_key}'")

    if context:
        await context.info(f"Creating proportions plot grouped by '{cluster_key}'")

    figsize = params.figure_size or (12, 4)

    scv.pl.proportions(
        adata,
        groupby=cluster_key,
        figsize=figsize,
        dpi=params.dpi,
        show=False,
    )

    fig = plt.gcf()
    title = params.title or f"Spliced/Unspliced Proportions by {cluster_key}"
    fig.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


async def _create_velocity_heatmap(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create velocity heatmap using scv.pl.heatmap.

    Shows gene expression patterns ordered by latent time.

    Data requirements:
        - adata.obs['latent_time']: Latent time from dynamical model
        - adata.var['velocity_genes']: Velocity genes (optional)
    """
    require("scvelo", feature="velocity heatmap")
    import scvelo as scv

    validate_obs_column(adata, "latent_time", "Latent time")

    if params.feature:
        if isinstance(params.feature, str):
            var_names = [params.feature]
        else:
            var_names = list(params.feature)
        valid_genes = [g for g in var_names if g in adata.var_names]
        if not valid_genes:
            raise DataNotFoundError(f"None of the specified genes found: {var_names}")
        var_names = valid_genes
    else:
        if "velocity_genes" in adata.var.columns:
            velocity_genes = adata.var_names[adata.var["velocity_genes"]]
            var_names = list(velocity_genes[:50])
        else:
            if "highly_variable" in adata.var.columns:
                hvg = adata.var_names[adata.var["highly_variable"]]
                var_names = list(hvg[:50])
            else:
                var_names = list(adata.var_names[:50])

    if context:
        await context.info(f"Creating velocity heatmap with {len(var_names)} genes")

    figsize = params.figure_size or (12, 8)

    scv.pl.heatmap(
        adata,
        var_names=var_names,
        sortby="latent_time",
        col_color=params.cluster_key,
        n_convolve=30,
        show=False,
        figsize=figsize,
    )

    fig = plt.gcf()
    fig.set_dpi(params.dpi)

    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


async def _create_velocity_paga_plot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create PAGA plot with velocity using scv.pl.paga.

    Shows partition-based graph abstraction with directed velocity arrows.

    Data requirements:
        - adata.uns['velocity_graph']: Velocity transition graph
        - adata.uns['paga']: PAGA results (computed by scv.tl.paga)
        - adata.obs[cluster_key]: Cluster labels used for PAGA
    """
    require("scvelo", feature="velocity PAGA plot")
    import scvelo as scv

    if "velocity_graph" not in adata.uns:
        raise DataNotFoundError("velocity_graph required. Run velocity analysis first.")

    cluster_key = params.cluster_key
    if not cluster_key:
        if "paga" in adata.uns and "groups" in adata.uns.get("paga", {}):
            cluster_key = adata.uns["paga"].get("groups")
        else:
            categorical_cols = get_categorical_columns(adata)
            if categorical_cols:
                cluster_key = categorical_cols[0]

    if not cluster_key or cluster_key not in adata.obs.columns:
        raise ParameterError(
            f"cluster_key is required for PAGA plot. "
            f"Available columns: {list(adata.obs.columns)[:10]}"
        )

    # Compute PAGA if not already done
    if "paga" not in adata.uns:
        if context:
            await context.info(f"Computing PAGA for cluster_key='{cluster_key}'")
        import scanpy as sc

        sc.tl.paga(adata, groups=cluster_key)
        scv.tl.paga(adata, groups=cluster_key)

    if context:
        await context.info(f"Creating velocity PAGA plot for '{cluster_key}'")

    basis = infer_basis(adata, preferred=params.basis, priority=["umap", "spatial"])
    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    scv.pl.paga(
        adata,
        basis=basis,
        color=cluster_key,
        ax=ax,
        show=False,
        frameon=params.show_axes,
    )

    if params.title:
        ax.set_title(params.title, fontsize=14)

    if basis == "spatial":
        ax.invert_yaxis()

    plt.tight_layout()
    return fig
