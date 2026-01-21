"""
Spatial statistics visualization functions for spatial transcriptomics.

This module contains:
- Neighborhood enrichment heatmaps
- Co-occurrence plots
- Ripley's function visualizations
- Moran's I barplots (standard spatial transcriptomics format)
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
from ...utils.adata_utils import (
    get_analysis_parameter,
    require_spatial_coords,
    validate_obs_column,
)
from ...utils.dependency_manager import require
from ...utils.exceptions import DataNotFoundError, ParameterError
from .core import (
    auto_spot_size,
    create_figure_from_params,
    get_categorical_columns,
    resolve_figure_size,
    setup_multi_panel_figure,
)


def _resolve_cluster_key(
    adata: "ad.AnnData",
    analysis_type: str,
    params_cluster_key: Optional[str],
) -> str:
    """Resolve cluster_key from params or stored metadata.

    Priority:
        1. User-provided cluster_key (params_cluster_key)
        2. cluster_key from analysis metadata

    Args:
        adata: AnnData object
        analysis_type: Analysis type (e.g., "neighborhood", "co_occurrence")
        params_cluster_key: User-provided cluster_key

    Returns:
        Resolved cluster_key

    Raises:
        ParameterError: If no cluster_key can be determined
    """
    cluster_key = params_cluster_key or get_analysis_parameter(
        adata, f"spatial_stats_{analysis_type}", "cluster_key"
    )
    if not cluster_key:
        categorical_cols = get_categorical_columns(adata, limit=10)
        raise ParameterError(
            f"cluster_key required for {analysis_type} visualization. "
            f"Available categorical columns: {', '.join(categorical_cols)}"
        )
    validate_obs_column(adata, cluster_key, "Cluster key")
    return cluster_key


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
        - moran: Moran's I barplot (top spatially variable genes)
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

    cluster_key = _resolve_cluster_key(adata, "neighborhood", params.cluster_key)

    enrichment_key = f"{cluster_key}_nhood_enrichment"
    if enrichment_key not in adata.uns:
        raise DataNotFoundError(
            f"Neighborhood enrichment not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    fig, axes = create_figure_from_params(params, "spatial")
    ax = axes[0]

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

    cluster_key = _resolve_cluster_key(adata, "co_occurrence", params.cluster_key)

    co_occurrence_key = f"{cluster_key}_co_occurrence"
    if co_occurrence_key not in adata.uns:
        raise DataNotFoundError(
            f"Co-occurrence not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    categories = adata.obs[cluster_key].cat.categories.tolist()
    clusters_to_show = categories[: min(4, len(categories))]

    # Calculate appropriate figsize based on number of clusters
    # squidpy default: (5 * n_clusters, 5) with constrained_layout=True
    # Only override if user explicitly provides figure_size
    if params.figure_size:
        figsize = params.figure_size
    else:
        # Let squidpy use its default sizing (5 inches per cluster, 5 height)
        figsize = None

    sq.pl.co_occurrence(
        adata,
        cluster_key=cluster_key,
        clusters=clusters_to_show,
        figsize=figsize,
        dpi=params.dpi,
    )

    fig = plt.gcf()
    if params.title:
        fig.suptitle(params.title, y=1.02)

    # Don't call tight_layout - squidpy uses constrained_layout=True internally
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

    cluster_key = _resolve_cluster_key(adata, "ripley", params.cluster_key)

    ripley_key = f"{cluster_key}_ripley_L"
    if ripley_key not in adata.uns:
        raise DataNotFoundError(
            f"Ripley results not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' and analysis_type='ripley' first."
        )

    fig, axes = create_figure_from_params(params, "spatial")
    ax = axes[0]

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
    """Create Moran's I barplot visualization (standard format).

    Shows top spatially variable genes ranked by Moran's I value.
    Color indicates statistical significance (-log10 p-value).
    This is the standard visualization format in spatial transcriptomics.

    Data requirements:
        - adata.uns['moranI']: DataFrame with I, pval_norm columns
    """

    if "moranI" not in adata.uns:
        raise DataNotFoundError("Moran's I results not found. Expected key: moranI")

    moran_data = adata.uns["moranI"].copy()

    # Prepare data for visualization
    moran_data["gene"] = moran_data.index
    pvals = moran_data["pval_norm"].values

    # Handle zero/negative p-values for log transform
    # Use data-driven minimum to avoid extreme values
    min_pval = max(1e-50, np.min(pvals[pvals > 0]) if np.any(pvals > 0) else 1e-50)
    pvals_safe = np.clip(pvals, min_pval, 1.0)
    moran_data["neg_log_pval"] = -np.log10(pvals_safe)

    # Mark significance
    moran_data["significant"] = pvals < 0.05

    # Sort by Moran's I (descending) and take top genes
    n_top = min(20, len(moran_data))  # Use actual data size if less than 20
    top_genes = moran_data.nlargest(n_top, "I")
    n_actual = len(top_genes)  # Actual number of genes to display

    # Create figure with appropriate size based on actual gene count
    # Width: 8 inches for gene names, Height: 0.4 per gene + margins
    if params.figure_size:
        figsize = params.figure_size
    else:
        # Minimum height of 3 for small gene counts, scale with actual genes
        figsize = (8, max(n_actual * 0.4 + 1.5, 3))

    fig, ax = plt.subplots(figsize=figsize)

    # Create horizontal barplot (easier to read gene names)
    # Color by -log10(p-value) to show significance
    norm = plt.Normalize(
        vmin=top_genes["neg_log_pval"].min(), vmax=top_genes["neg_log_pval"].max()
    )
    cmap = plt.colormaps.get_cmap(params.colormap or "viridis")
    colors = [cmap(norm(v)) for v in top_genes["neg_log_pval"].values]

    # Plot bars
    y_pos = np.arange(len(top_genes))
    ax.barh(
        y_pos,
        top_genes["I"].values,
        color=colors,
        alpha=params.alpha,
        edgecolor="black",
        linewidth=0.5,
    )

    # Set y-axis labels (gene names)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_genes["gene"].values)

    # Invert y-axis so highest Moran's I is at top
    ax.invert_yaxis()

    # Add significance markers
    for i, (idx, row) in enumerate(top_genes.iterrows()):
        if row["significant"]:
            ax.text(
                row["I"] + 0.01,
                i,
                "*",
                va="center",
                ha="left",
                fontsize=12,
                fontweight="bold",
            )

    # Labels and title
    title = params.title or "Top Spatially Variable Genes (Moran's I)"
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Moran's I (spatial autocorrelation)", fontsize=12)
    ax.set_ylabel("Gene", fontsize=12)

    # Add colorbar for significance
    if params.show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_label("-log10(p-value)", fontsize=10)

    # Add vertical line at I=0 for reference
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5, linewidth=1)

    # Add annotation for significance
    n_significant = top_genes["significant"].sum()
    ax.text(
        0.98,
        0.02,
        f"* p < 0.05 ({n_significant}/{n_actual} significant)",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        style="italic",
        color="gray",
    )

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

    cluster_key = _resolve_cluster_key(adata, "centrality", params.cluster_key)

    centrality_key = f"{cluster_key}_centrality_scores"
    if centrality_key not in adata.uns:
        raise DataNotFoundError(
            f"Centrality scores not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    figsize = resolve_figure_size(params, "spatial")

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

    # Only use figure suptitle for multi-panel plots
    # For single panel, axes title is sufficient
    fig, axes = setup_multi_panel_figure(
        n_panels=len(genes_to_plot),
        params=params,
        default_title="Getis-Ord Gi* Hotspots/Coldspots"
        if len(genes_to_plot) > 1
        else "",
    )

    coords = require_spatial_coords(adata)

    # Calculate spot size (auto or user-specified)
    spot_size = auto_spot_size(adata, params.spot_size, basis="spatial")

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
                s=spot_size,
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

            # Format title based on number of panels
            if len(genes_to_plot) == 1:
                # Single panel: include full description in title
                ax.set_title(
                    f"{gene}\nGetis-Ord Gi*: Hot: {hot_spots}, Cold: {cold_spots}"
                )
            elif params.add_gene_labels:
                ax.set_title(f"{gene}\nHot: {hot_spots}, Cold: {cold_spots}")
            else:
                ax.set_title(f"{gene}")

            ax.set_xlabel("Spatial X")
            ax.set_ylabel("Spatial Y")
            ax.set_aspect("equal")
            ax.invert_yaxis()

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    return fig
