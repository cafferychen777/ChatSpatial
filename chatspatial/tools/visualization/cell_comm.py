"""
Cell communication visualization functions for spatial transcriptomics.

This module contains:
- LIANA+ cluster-based visualizations (dotplot, tileplot, circle_plot)
- LIANA+ spatial bivariate visualizations
- CellPhoneDB visualizations (heatmap, dotplot, chord)
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...utils.adata_utils import get_spatial_coordinates, validate_obs_column
from ...utils.dependency_manager import require
from ...utils.exceptions import DataNotFoundError, ParameterError, ProcessingError
from .core import CellCommunicationData

# =============================================================================
# Data Retrieval
# =============================================================================


async def get_cell_communication_data(
    adata: "ad.AnnData",
    method: Optional[str] = None,
    context: Optional["ToolContext"] = None,
) -> CellCommunicationData:
    """
    Unified function to retrieve cell communication results from AnnData.

    This function consolidates all cell communication data retrieval logic into
    a single, consistent interface. It handles:
    - LIANA+ spatial bivariate analysis results
    - LIANA+ cluster-based analysis results
    - CellPhoneDB analysis results

    Args:
        adata: AnnData object with cell communication results
        method: Analysis method hint (optional)
        context: MCP context for logging

    Returns:
        CellCommunicationData object with results and metadata

    Raises:
        DataNotFoundError: No cell communication results found
    """
    # Check for LIANA+ spatial bivariate results (highest priority)
    if "liana_spatial_scores" in adata.obsm:
        spatial_scores = adata.obsm["liana_spatial_scores"]
        lr_pairs = adata.uns.get("liana_spatial_interactions", [])
        results_df = adata.uns.get("liana_spatial_res", pd.DataFrame())

        if not isinstance(results_df, pd.DataFrame):
            results_df = pd.DataFrame()

        if context:
            await context.info(
                f"Found LIANA+ spatial results: {len(lr_pairs)} LR pairs, "
                f"{spatial_scores.shape[0]} spots"
            )

        return CellCommunicationData(
            results=results_df,
            method="liana_spatial",
            analysis_type="spatial",
            lr_pairs=lr_pairs if lr_pairs else [],
            spatial_scores=spatial_scores,
            spatial_pvals=adata.obsm.get("liana_spatial_pvals"),
            results_key="liana_spatial_res",
        )

    # Check for LIANA+ cluster-based results
    if "liana_res" in adata.uns:
        results = adata.uns["liana_res"]
        if isinstance(results, pd.DataFrame) and len(results) > 0:
            if (
                "ligand_complex" in results.columns
                and "receptor_complex" in results.columns
            ):
                lr_pairs = (
                    (results["ligand_complex"] + "^" + results["receptor_complex"])
                    .unique()
                    .tolist()
                )
            else:
                lr_pairs = []

            source_labels = (
                results["source"].unique().tolist()
                if "source" in results.columns
                else None
            )
            target_labels = (
                results["target"].unique().tolist()
                if "target" in results.columns
                else None
            )

            if context:
                await context.info(
                    f"Found LIANA+ cluster results: {len(lr_pairs)} LR pairs"
                )

            return CellCommunicationData(
                results=results,
                method="liana_cluster",
                analysis_type="cluster",
                lr_pairs=lr_pairs,
                source_labels=source_labels,
                target_labels=target_labels,
                results_key="liana_res",
            )

    # Check for CellPhoneDB results
    if "cellphonedb_means" in adata.uns:
        means = adata.uns["cellphonedb_means"]
        if isinstance(means, pd.DataFrame):
            lr_pairs = means.index.tolist()

            if context:
                await context.info(
                    f"Found CellPhoneDB results: {len(lr_pairs)} LR pairs"
                )

            return CellCommunicationData(
                results=means,
                method="cellphonedb",
                analysis_type="cluster",
                lr_pairs=lr_pairs,
                results_key="cellphonedb_means",
            )

    # No results found
    raise DataNotFoundError(
        "No cell communication results found. "
        "Run analyze_cell_communication() first."
    )


# =============================================================================
# Main Router
# =============================================================================


async def create_cell_communication_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create cell communication visualization using unified data retrieval.

    Routes to appropriate visualization based on analysis type and subtype:
    - Spatial analysis: Multi-panel spatial plot
    - Cluster analysis: LIANA+ visualizations or CellPhoneDB

    Args:
        adata: AnnData object with cell communication results
        params: Visualization parameters (use params.subtype to select viz type)
        context: MCP context for logging

    Returns:
        matplotlib Figure object
    """
    if context:
        await context.info("Creating cell communication visualization")

    data = await get_cell_communication_data(adata, context=context)

    if context:
        await context.info(
            f"Using {data.method} results ({data.analysis_type} analysis, "
            f"{len(data.lr_pairs)} LR pairs)"
        )

    if data.analysis_type == "spatial":
        return _create_spatial_lr_visualization(adata, data, params, context)
    else:
        if data.method == "cellphonedb":
            subtype = params.subtype or "heatmap"
            if subtype == "dotplot":
                return _create_cellphonedb_dotplot(adata, data, params, context)
            elif subtype == "chord":
                return _create_cellphonedb_chord(adata, data, params, context)
            else:
                return _create_cellphonedb_heatmap(adata, data, params, context)
        else:
            subtype = params.subtype or "dotplot"
            if subtype == "tileplot":
                return await _create_liana_tileplot(adata, data, params, context)
            elif subtype == "circle_plot":
                return await _create_liana_circle_plot(adata, data, params, context)
            else:
                return await _create_cluster_lr_visualization(
                    adata, data, params, context
                )


# =============================================================================
# LIANA+ Visualizations
# =============================================================================


def _create_spatial_lr_visualization(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create spatial L-R visualization using scanpy (official LIANA+ approach)."""
    if data.spatial_scores is None or len(data.lr_pairs) == 0:
        raise DataNotFoundError(
            "No spatial communication scores found. Run spatial analysis first."
        )

    n_pairs = min(params.plot_top_pairs or 6, len(data.lr_pairs), 6)

    # Determine top pairs based on global metric
    if len(data.results) > 0:
        metric_col = None
        for col in ["morans", "lee", "global_score"]:
            if col in data.results.columns:
                metric_col = col
                break

        if metric_col:
            top_results = data.results.nlargest(n_pairs, metric_col)
            top_pairs = top_results.index.tolist()
        else:
            top_pairs = data.lr_pairs[:n_pairs]
    else:
        top_pairs = data.lr_pairs[:n_pairs]

    if not top_pairs:
        raise DataNotFoundError("No LR pairs found in spatial results.")

    # Get pair indices
    pair_indices = []
    valid_pairs = []
    for pair in top_pairs:
        if pair in data.lr_pairs:
            pair_indices.append(data.lr_pairs.index(pair))
            valid_pairs.append(pair)

    if not valid_pairs:
        valid_pairs = data.lr_pairs[:n_pairs]
        pair_indices = list(range(len(valid_pairs)))

    # Create figure
    n_panels = len(valid_pairs)
    n_cols = min(3, n_panels)
    n_rows = (n_panels + n_cols - 1) // n_cols

    figsize = params.figure_size or (5 * n_cols, 4 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)

    if n_panels == 1:
        axes = np.array([axes])
    axes = np.atleast_1d(axes).flatten()

    x_coords, y_coords = get_spatial_coordinates(adata)

    for i, (pair, pair_idx) in enumerate(zip(valid_pairs, pair_indices, strict=False)):
        ax = axes[i]

        if pair_idx < data.spatial_scores.shape[1]:
            scores = data.spatial_scores[:, pair_idx]
        else:
            scores = np.zeros(len(adata))

        scatter = ax.scatter(
            x_coords,
            y_coords,
            c=scores,
            cmap=params.colormap or "viridis",
            s=params.spot_size or 15,
            alpha=params.alpha or 0.8,
            edgecolors="none",
        )

        display_name = pair.replace("^", " → ").replace("_", " → ")

        if len(data.results) > 0 and pair in data.results.index:
            for metric in ["morans", "lee", "global_score"]:
                if metric in data.results.columns:
                    val = data.results.loc[pair, metric]
                    display_name += f"\n({metric}: {val:.3f})"
                    break

        ax.set_title(display_name, fontsize=10)
        ax.set_aspect("equal")
        ax.set_xlabel("")
        ax.set_ylabel("")
        plt.colorbar(scatter, ax=ax, shrink=0.7, label="Score")

    for i in range(n_panels, len(axes)):
        axes[i].set_visible(False)

    plt.suptitle("Spatial Cell Communication", fontsize=14, fontweight="bold")
    plt.tight_layout()

    return fig


async def _create_cluster_lr_visualization(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create cluster-based L-R visualization using LIANA+ dotplot."""
    require("liana", feature="LIANA+ plotting")
    require("plotnine", feature="LIANA+ plotting")
    import liana as li

    if context:
        await context.info("Using LIANA+ official dotplot")

    try:
        orderby_col = None
        for col in ["magnitude_rank", "specificity_rank", "lr_means"]:
            if col in data.results.columns:
                orderby_col = col
                break

        if orderby_col is None:
            raise DataNotFoundError("No valid orderby column found in LIANA results")

        p = li.pl.dotplot(
            adata=adata,
            uns_key=data.results_key,
            colour=(
                "magnitude_rank" if "magnitude_rank" in data.results.columns else None
            ),
            size=(
                "specificity_rank"
                if "specificity_rank" in data.results.columns
                else None
            ),
            orderby=orderby_col,
            orderby_ascending=True,
            top_n=params.plot_top_pairs or 20,
            inverse_colour=True,
            inverse_size=True,
            cmap=params.colormap or "viridis",
            figure_size=params.figure_size or (10, 8),
            return_fig=True,
        )

        fig = _plotnine_to_matplotlib(p, params)
        return fig

    except Exception as e:
        raise ProcessingError(
            f"LIANA+ dotplot failed: {e}\n\n"
            "Ensure cell communication analysis completed successfully."
        ) from e


async def _create_liana_tileplot(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create LIANA+ tileplot visualization."""
    try:
        import liana as li

        if context:
            await context.info("Creating LIANA+ tileplot")

        orderby_col = None
        for col in ["magnitude_rank", "specificity_rank", "lr_means"]:
            if col in data.results.columns:
                orderby_col = col
                break

        if orderby_col is None:
            raise DataNotFoundError("No valid orderby column found in LIANA results")

        fill_col = (
            "magnitude_rank"
            if "magnitude_rank" in data.results.columns
            else orderby_col
        )
        label_col = "lr_means" if "lr_means" in data.results.columns else fill_col

        p = li.pl.tileplot(
            adata=adata,
            uns_key=data.results_key,
            fill=fill_col,
            label=label_col,
            orderby=orderby_col,
            orderby_ascending=True,
            top_n=params.plot_top_pairs or 15,
            figure_size=params.figure_size or (14, 8),
            return_fig=True,
        )

        fig = _plotnine_to_matplotlib(p, params)
        return fig

    except Exception as e:
        raise ProcessingError(
            f"LIANA+ tileplot failed: {e}\n\n"
            "Ensure cell communication analysis completed successfully."
        ) from e


async def _create_liana_circle_plot(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create LIANA+ circle plot (network diagram) visualization."""
    try:
        import liana as li

        if context:
            await context.info("Creating LIANA+ circle plot")

        score_col = None
        for col in ["magnitude_rank", "specificity_rank", "lr_means"]:
            if col in data.results.columns:
                score_col = col
                break

        if score_col is None:
            raise DataNotFoundError("No valid score column found in LIANA results")

        groupby = params.cluster_key
        if groupby is None:
            if "source" in data.results.columns:
                groupby = (
                    data.results["source"].iloc[0] if len(data.results) > 0 else None
                )
            if groupby is None:
                raise ParameterError(
                    "cluster_key is required for circle_plot. "
                    "Specify the cell type column used in analysis."
                )

        fig_size = params.figure_size or (10, 10)
        fig, ax = plt.subplots(figsize=fig_size)

        li.pl.circle_plot(
            adata=adata,
            uns_key=data.results_key,
            groupby=groupby,
            score_key=score_col,
            inverse_score=True,
            top_n=params.plot_top_pairs * 3 if params.plot_top_pairs else 50,
            orderby=score_col,
            orderby_ascending=True,
            figure_size=fig_size,
        )

        fig = plt.gcf()
        return fig

    except Exception as e:
        raise ProcessingError(
            f"LIANA+ circle_plot failed: {e}\n\n"
            "Ensure cell communication analysis completed successfully."
        ) from e


# =============================================================================
# CellPhoneDB Visualizations
# =============================================================================


def _create_cellphonedb_heatmap(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CellPhoneDB heatmap visualization using ktplotspy."""
    import ktplotspy as kpy

    means = data.results

    if not isinstance(means, pd.DataFrame) or len(means) == 0:
        raise DataNotFoundError("CellPhoneDB results empty. Re-run analysis.")

    pvalues = adata.uns.get("cellphonedb_pvalues", None)

    if pvalues is None or not isinstance(pvalues, pd.DataFrame):
        raise DataNotFoundError("CellPhoneDB pvalues not found. Re-run analysis.")

    grid = kpy.plot_cpdb_heatmap(
        pvals=pvalues,
        title=params.title or "CellPhoneDB: Significant Interactions",
        alpha=0.05,
        symmetrical=True,
    )

    return grid.fig


def _create_cellphonedb_dotplot(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CellPhoneDB dotplot visualization using ktplotspy."""
    means = data.results

    if not isinstance(means, pd.DataFrame) or len(means) == 0:
        raise DataNotFoundError("CellPhoneDB results empty. Re-run analysis.")

    require("ktplotspy", feature="CellPhoneDB dotplot visualization")
    import ktplotspy as kpy

    try:
        pvalues = adata.uns.get("cellphonedb_pvalues", None)

        if pvalues is None or not isinstance(pvalues, pd.DataFrame):
            raise DataNotFoundError("Missing pvalues DataFrame for ktplotspy dotplot")

        cluster_key = params.cluster_key or "leiden"
        validate_obs_column(adata, cluster_key, "cluster_key")

        gg = kpy.plot_cpdb(
            adata=adata,
            cell_type1=".",
            cell_type2=".",
            means=means,
            pvals=pvalues,
            celltype_key=cluster_key,
            genes=None,
            figsize=params.figure_size or (12, 10),
            title="CellPhoneDB: L-R Interactions",
            max_size=10,
            alpha=0.05,
            keep_significant_only=True,
            standard_scale=True,
        )

        fig = gg.draw()
        return fig

    except Exception as e:
        raise ProcessingError(
            f"Failed to create CellPhoneDB dotplot: {str(e)}\n\n"
            "Try using subtype='heatmap' instead."
        ) from e


def _create_cellphonedb_chord(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CellPhoneDB chord/circos diagram using ktplotspy."""
    from matplotlib.lines import Line2D

    means = data.results

    if not isinstance(means, pd.DataFrame) or len(means) == 0:
        raise DataNotFoundError("CellPhoneDB results empty. Re-run analysis.")

    require("ktplotspy", feature="CellPhoneDB chord visualization")
    import ktplotspy as kpy
    import matplotlib.colors as mcolors

    try:
        pvalues = adata.uns.get("cellphonedb_pvalues", None)
        deconvoluted = adata.uns.get("cellphonedb_deconvoluted", None)

        if pvalues is None or not isinstance(pvalues, pd.DataFrame):
            raise DataNotFoundError(
                "Missing pvalues DataFrame for ktplotspy chord plot"
            )

        if deconvoluted is None or not isinstance(deconvoluted, pd.DataFrame):
            raise DataNotFoundError(
                "Missing deconvoluted DataFrame for chord plot. "
                "Re-run CellPhoneDB analysis."
            )

        cluster_key = params.cluster_key or "leiden"
        validate_obs_column(adata, cluster_key, "cluster_key")

        link_colors = None
        legend_items = []

        if "interacting_pair" in deconvoluted.columns:
            unique_pairs = deconvoluted["interacting_pair"].unique()
            n_pairs = min(params.plot_top_pairs or 50, len(unique_pairs))
            top_pairs = unique_pairs[:n_pairs]

            if n_pairs <= 10:
                cmap = plt.cm.get_cmap("tab10", 10)
            elif n_pairs <= 20:
                cmap = plt.cm.get_cmap("tab20", 20)
            else:
                cmap = plt.cm.get_cmap("nipy_spectral", n_pairs)

            link_colors = {}
            for i, pair in enumerate(top_pairs):
                color = mcolors.rgb2hex(cmap(i % cmap.N))
                link_colors[pair] = color
                legend_items.append((pair, color))

        circos = kpy.plot_cpdb_chord(
            adata=adata,
            means=means,
            pvals=pvalues,
            deconvoluted=deconvoluted,
            celltype_key=cluster_key,
            cell_type1=".",
            cell_type2=".",
            link_colors=link_colors,
        )

        fig = circos.ax.figure
        fig.set_size_inches(14, 10)

        if legend_items:
            line_handles = [
                Line2D([], [], color=color, label=label, linewidth=2)
                for label, color in legend_items
            ]

            legend = circos.ax.legend(
                handles=line_handles,
                loc="center left",
                bbox_to_anchor=(1.15, 0.5),
                fontsize=6,
                frameon=True,
                framealpha=0.9,
                title="L-R Pairs",
                title_fontsize=7,
            )

            fig._chatspatial_extra_artists = [legend]

        return fig

    except Exception as e:
        raise ProcessingError(
            f"Failed to create CellPhoneDB chord diagram: {str(e)}\n\n"
            "Try using subtype='heatmap' instead."
        ) from e


# =============================================================================
# Utilities
# =============================================================================


def _plotnine_to_matplotlib(p, params: VisualizationParameters) -> plt.Figure:
    """Convert plotnine ggplot object to matplotlib Figure.

    Uses plotnine's native draw() method which returns the underlying
    matplotlib Figure, avoiding rasterization through PNG buffer.
    """
    try:
        # plotnine's draw() returns the matplotlib Figure directly
        fig = p.draw()

        # Apply DPI setting if specified
        if params.dpi:
            fig.set_dpi(params.dpi)

        return fig

    except Exception as e:
        raise ProcessingError(f"Failed to convert plotnine figure: {e}") from e
