"""
Spatial statistics visualization functions for spatial transcriptomics.

This module contains:
- Neighborhood enrichment heatmaps
- Co-occurrence plots
- Ripley's function visualizations
- Moran's I scatter plots
- Centrality scores
- Getis-Ord Gi* hotspot maps
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...utils.adata_utils import require_spatial_coords
from ...utils.dependency_manager import require
from ...utils.exceptions import DataNotFoundError, ParameterError
from .core import get_categorical_columns, setup_multi_panel_figure

# =============================================================================
# Main Router
# =============================================================================


async def create_spatial_statistics_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create spatial statistics visualization based on subtype.

    Args:
        adata: AnnData object with spatial statistics results
        params: Visualization parameters including subtype
        context: MCP context

    Returns:
        Matplotlib figure with spatial statistics visualization

    Subtypes:
        - neighborhood: Neighborhood enrichment heatmap
        - co_occurrence: Co-occurrence analysis plot
        - ripley: Ripley's K/L function curves
        - moran: Moran's I scatter plot
        - centrality: Graph centrality scores
        - getis_ord: Getis-Ord Gi* hotspot/coldspot maps
    """
    subtype = params.subtype or "neighborhood"

    if context:
        await context.info(f"Creating {subtype} spatial statistics visualization")

    if subtype == "neighborhood":
        return await _create_neighborhood_enrichment_visualization(
            adata, params, context
        )
    elif subtype == "co_occurrence":
        return await _create_co_occurrence_visualization(adata, params, context)
    elif subtype == "ripley":
        return await _create_ripley_visualization(adata, params, context)
    elif subtype == "moran":
        return _create_moran_visualization(adata, params, context)
    elif subtype == "centrality":
        return await _create_centrality_visualization(adata, params, context)
    elif subtype == "getis_ord":
        return await _create_getis_ord_visualization(adata, params, context)
    else:
        raise ParameterError(
            f"Unsupported subtype for spatial_statistics: '{subtype}'. "
            f"Available subtypes: neighborhood, co_occurrence, ripley, moran, "
            f"centrality, getis_ord"
        )


# =============================================================================
# Visualization Functions
# =============================================================================


async def _create_neighborhood_enrichment_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create neighborhood enrichment visualization using squidpy.

    Data requirements:
        - adata.uns['{cluster_key}_nhood_enrichment']: Enrichment results
        - adata.obs[cluster_key]: Cluster labels
    """
    require("squidpy", feature="neighborhood enrichment visualization")
    import squidpy as sq

    # Infer cluster_key from params or existing results
    cluster_key = params.cluster_key
    if not cluster_key:
        enrichment_keys = [
            k for k in adata.uns.keys() if k.endswith("_nhood_enrichment")
        ]
        if enrichment_keys:
            cluster_key = enrichment_keys[0].replace("_nhood_enrichment", "")
            if context:
                await context.info(f"Inferred cluster_key: '{cluster_key}'")
        else:
            categorical_cols = get_categorical_columns(adata, limit=10)
            raise ParameterError(
                f"cluster_key required. Available: {', '.join(categorical_cols)}"
            )

    enrichment_key = f"{cluster_key}_nhood_enrichment"
    if enrichment_key not in adata.uns:
        raise DataNotFoundError(
            f"Neighborhood enrichment not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    sq.pl.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        cmap=params.colormap or "coolwarm",
        ax=ax,
        title=params.title or f"Neighborhood Enrichment ({cluster_key})",
    )

    plt.tight_layout()
    return fig


async def _create_co_occurrence_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create co-occurrence analysis visualization using squidpy.

    Data requirements:
        - adata.uns['{cluster_key}_co_occurrence']: Co-occurrence results
        - adata.obs[cluster_key]: Cluster labels
    """
    require("squidpy", feature="co-occurrence visualization")
    import squidpy as sq

    # Infer cluster_key from params or existing results
    cluster_key = params.cluster_key
    if not cluster_key:
        co_occurrence_keys = [
            k for k in adata.uns.keys() if k.endswith("_co_occurrence")
        ]
        if co_occurrence_keys:
            cluster_key = co_occurrence_keys[0].replace("_co_occurrence", "")
            if context:
                await context.info(f"Inferred cluster_key: '{cluster_key}'")
        else:
            categorical_cols = get_categorical_columns(adata, limit=10)
            raise ParameterError(
                f"cluster_key required. Available: {', '.join(categorical_cols)}"
            )

    co_occurrence_key = f"{cluster_key}_co_occurrence"
    if co_occurrence_key not in adata.uns:
        raise DataNotFoundError(
            f"Co-occurrence not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    categories = adata.obs[cluster_key].cat.categories.tolist()
    clusters_to_show = categories[: min(4, len(categories))]

    figsize = params.figure_size or (12, 10)

    sq.pl.co_occurrence(
        adata,
        cluster_key=cluster_key,
        clusters=clusters_to_show,
        figsize=figsize,
        dpi=params.dpi,
    )

    fig = plt.gcf()
    if params.title:
        fig.suptitle(params.title)

    plt.tight_layout()
    return fig


async def _create_ripley_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create Ripley's function visualization using squidpy.

    Data requirements:
        - adata.uns['{cluster_key}_ripley_L']: Ripley's L function results
        - adata.obs[cluster_key]: Cluster labels
    """
    require("squidpy", feature="Ripley visualization")
    import squidpy as sq

    # Infer cluster_key from params or existing results
    cluster_key = params.cluster_key
    if not cluster_key:
        ripley_keys = [k for k in adata.uns.keys() if k.endswith("_ripley_L")]
        if ripley_keys:
            cluster_key = ripley_keys[0].replace("_ripley_L", "")
            if context:
                await context.info(f"Inferred cluster_key: '{cluster_key}'")
        else:
            categorical_cols = get_categorical_columns(adata, limit=10)
            raise ParameterError(
                f"cluster_key required. Available: {', '.join(categorical_cols)}"
            )

    ripley_key = f"{cluster_key}_ripley_L"
    if ripley_key not in adata.uns:
        raise DataNotFoundError(
            f"Ripley results not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' and analysis_type='ripley' first."
        )

    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    sq.pl.ripley(adata, cluster_key=cluster_key, mode="L", plot_sims=True, ax=ax)

    if params.title:
        ax.set_title(params.title)

    plt.tight_layout()
    return fig


def _create_moran_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create Moran's I volcano-style visualization.

    Shows -log10(p-value) vs Moran's I for spatially variable genes.
    Color indicates Moran's I value (positive = clustered, negative = dispersed).

    Data requirements:
        - adata.uns['moranI']: DataFrame with I, pval_norm columns
    """
    if "moranI" not in adata.uns:
        raise DataNotFoundError("Moran's I results not found. Expected key: moranI")

    moran_data = adata.uns["moranI"]
    moran_i = moran_data["I"].values
    pvals = moran_data["pval_norm"].values

    # Handle zero/negative p-values for log transform
    pvals_safe = np.clip(pvals, 1e-300, 1.0)
    neg_log_pval = -np.log10(pvals_safe)

    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    # Color by Moran's I value (meaningful: positive=clustered, negative=dispersed)
    scatter = ax.scatter(
        neg_log_pval,
        moran_i,
        s=50,
        alpha=params.alpha,
        c=moran_i,
        cmap="RdBu_r",  # Diverging colormap centered at 0
        vmin=-max(abs(moran_i.min()), abs(moran_i.max())),
        vmax=max(abs(moran_i.min()), abs(moran_i.max())),
    )

    # Label top significant genes (high I and low p-value)
    gene_names = moran_data.index.tolist()
    sig_threshold = -np.log10(0.05)
    significant_mask = (neg_log_pval > sig_threshold) & (moran_i > 0)

    if np.any(significant_mask):
        # Sort by combined score (high I * high significance)
        scores = moran_i * neg_log_pval
        top_indices = np.argsort(scores)[::-1][:10]  # Top 10

        for idx in top_indices:
            if significant_mask[idx]:
                ax.annotate(
                    gene_names[idx],
                    (neg_log_pval[idx], moran_i[idx]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    alpha=0.8,
                )

    # Reference lines
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5, label="I=0")
    ax.axvline(x=sig_threshold, color="red", linestyle="--", alpha=0.5, label="p=0.05")

    title = params.title or "Moran's I Spatial Autocorrelation"
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("-log₁₀(p-value)", fontsize=12)
    ax.set_ylabel("Moran's I", fontsize=12)

    if params.show_colorbar:
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Moran's I (+ clustered, − dispersed)", fontsize=10)

    # Add legend for reference lines
    ax.legend(loc="upper left", fontsize=9)

    plt.tight_layout()
    return fig


async def _create_centrality_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create centrality scores visualization using squidpy.

    Data requirements:
        - adata.uns['{cluster_key}_centrality_scores']: Centrality scores
        - adata.obs[cluster_key]: Cluster labels
    """
    require("squidpy", feature="centrality visualization")
    import squidpy as sq

    # Infer cluster_key from params or existing results
    cluster_key = params.cluster_key
    if not cluster_key:
        centrality_keys = [
            k for k in adata.uns.keys() if k.endswith("_centrality_scores")
        ]
        if centrality_keys:
            cluster_key = centrality_keys[0].replace("_centrality_scores", "")
            if context:
                await context.info(f"Inferred cluster_key: '{cluster_key}'")
        else:
            categorical_cols = get_categorical_columns(adata, limit=10)
            raise ParameterError(
                f"cluster_key required. Available: {', '.join(categorical_cols)}"
            )

    centrality_key = f"{cluster_key}_centrality_scores"
    if centrality_key not in adata.uns:
        raise DataNotFoundError(
            f"Centrality scores not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    figsize = params.figure_size or (10, 8)

    sq.pl.centrality_scores(
        adata,
        cluster_key=cluster_key,
        figsize=figsize,
        dpi=params.dpi,
    )

    fig = plt.gcf()
    if params.title:
        fig.suptitle(params.title)

    plt.tight_layout()
    return fig


async def _create_getis_ord_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create Getis-Ord Gi* hotspot/coldspot visualization.

    Data requirements:
        - adata.obs['{gene}_getis_ord_z']: Z-scores for each gene
        - adata.obs['{gene}_getis_ord_p']: P-values for each gene
        - adata.obsm['spatial']: Spatial coordinates
    """
    # Find genes with Getis-Ord results
    getis_ord_genes = []
    for col in adata.obs.columns:
        if col.endswith("_getis_ord_z"):
            gene = col.replace("_getis_ord_z", "")
            if f"{gene}_getis_ord_p" in adata.obs.columns:
                getis_ord_genes.append(gene)

    if not getis_ord_genes:
        raise DataNotFoundError("No Getis-Ord results found in adata.obs")

    # Get genes to plot
    feature_list = (
        params.feature
        if isinstance(params.feature, list)
        else ([params.feature] if params.feature else [])
    )
    if feature_list:
        genes_to_plot = [g for g in feature_list if g in getis_ord_genes]
    else:
        genes_to_plot = getis_ord_genes[:6]

    if not genes_to_plot:
        raise DataNotFoundError(
            f"None of the specified genes have Getis-Ord results: {feature_list}"
        )

    if context:
        await context.info(
            f"Plotting Getis-Ord results for {len(genes_to_plot)} genes: {genes_to_plot}"
        )

    fig, axes = setup_multi_panel_figure(
        n_panels=len(genes_to_plot),
        params=params,
        default_title="Getis-Ord Gi* Hotspots/Coldspots",
    )

    coords = require_spatial_coords(adata)

    for i, gene in enumerate(genes_to_plot):
        if i < len(axes):
            ax = axes[i]
            z_key = f"{gene}_getis_ord_z"
            p_key = f"{gene}_getis_ord_p"

            if z_key not in adata.obs or p_key not in adata.obs:
                ax.text(
                    0.5,
                    0.5,
                    f"No Getis-Ord data for {gene}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{gene} (No Data)")
                continue

            z_scores = adata.obs[z_key].values
            p_vals = adata.obs[p_key].values

            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=z_scores,
                cmap="RdBu_r",
                s=params.spot_size or 20,
                alpha=params.alpha,
                vmin=-3,
                vmax=3,
            )

            if params.show_colorbar:
                plt.colorbar(scatter, ax=ax, label="Gi* Z-score")

            # Count significant hot and cold spots
            alpha = 0.05
            significant = p_vals < alpha
            hot_spots = np.sum((z_scores > 0) & significant)
            cold_spots = np.sum((z_scores < 0) & significant)

            if params.add_gene_labels:
                ax.set_title(f"{gene}\nHot: {hot_spots}, Cold: {cold_spots}")
            else:
                ax.set_title(f"{gene}")

            ax.set_xlabel("Spatial X")
            ax.set_ylabel("Spatial Y")
            ax.set_aspect("equal")
            ax.invert_yaxis()

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    return fig
