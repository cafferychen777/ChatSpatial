"""
Enrichment analysis visualization functions.

This module contains:
- Pathway enrichment barplots and dotplots
- GSEA enrichment score plots
- Spatial enrichment score visualization
- EnrichMap spatial autocorrelation plots
"""

from typing import TYPE_CHECKING, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...utils.adata_utils import validate_obs_column
from ...utils.exceptions import DataNotFoundError, ParameterError, ProcessingError
from .core import (
    create_figure,
    get_categorical_columns,
    plot_spatial_feature,
    setup_multi_panel_figure,
)

# =============================================================================
# Helper Functions
# =============================================================================


def _ensure_enrichmap_compatibility(adata: "ad.AnnData") -> None:
    """Ensure data has required metadata structure for EnrichMap visualization.

    EnrichMap and squidpy require:
    1. adata.obs['library_id'] - sample identifier column
    2. adata.uns['spatial'] - spatial metadata dictionary

    This function adds minimal metadata for single-sample data without these.
    """
    if "library_id" not in adata.obs.columns:
        adata.obs["library_id"] = "sample_1"

    if "spatial" not in adata.uns:
        library_ids = adata.obs["library_id"].unique()
        adata.uns["spatial"] = {}
        for lib_id in library_ids:
            adata.uns["spatial"][lib_id] = {
                "images": {},
                "scalefactors": {
                    "spot_diameter_fullres": 1.0,
                    "tissue_hires_scalef": 1.0,
                    "fiducial_diameter_fullres": 1.0,
                    "tissue_lowres_scalef": 1.0,
                },
            }


def _get_score_columns(adata: "ad.AnnData") -> List[str]:
    """Get all enrichment score columns from adata.obs."""
    return [col for col in adata.obs.columns if col.endswith("_score")]


def _resolve_score_column(
    adata: "ad.AnnData",
    feature: Optional[str],
    score_cols: List[str],
) -> str:
    """Resolve feature name to actual score column name."""
    if feature:
        if feature in adata.obs.columns:
            return feature
        if f"{feature}_score" in adata.obs.columns:
            return f"{feature}_score"
        raise DataNotFoundError(
            f"Score column '{feature}' not found. Available: {score_cols}"
        )
    if score_cols:
        return score_cols[0]
    raise DataNotFoundError("No enrichment scores found in adata.obs")


# =============================================================================
# Main Routers
# =============================================================================


async def create_enrichment_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create enrichment score visualization.

    Routes to appropriate visualization based on params:
    - violin: Enrichment scores violin plot by cluster
    - spatial_*: EnrichMap spatial visualizations
    - Default: Standard spatial scatter plot

    Args:
        adata: AnnData object with enrichment scores
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        Matplotlib figure
    """
    if context:
        await context.info("Creating enrichment visualization")

    score_cols = _get_score_columns(adata)
    if not score_cols:
        raise DataNotFoundError(
            "No enrichment scores found. Run 'analyze_enrichment' first."
        )

    # Route based on plot_type or subtype
    if params.plot_type == "violin":
        return _create_enrichment_violin(adata, params, score_cols, context)

    if params.subtype and params.subtype.startswith("spatial_"):
        return _create_enrichmap_spatial(adata, params, score_cols, context)

    # Default: spatial scatter plot
    return await _create_enrichment_spatial(adata, params, score_cols, context)


async def create_pathway_enrichment_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create pathway enrichment visualization (GSEA/ORA results).

    Supports multiple visualization types:
    Traditional:
    - barplot: Top enriched pathways barplot
    - dotplot: Multi-cluster enrichment dotplot
    - enrichment_plot: Classic GSEA running score plot

    Spatial EnrichMap:
    - spatial_score, spatial_correlogram, etc.

    Args:
        adata: AnnData object with enrichment results
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        Matplotlib figure
    """
    if context:
        await context.info("Creating pathway enrichment visualization")

    plot_type = params.subtype or "barplot"

    # Route spatial subtypes to enrichment visualization
    if plot_type.startswith("spatial_"):
        return await create_enrichment_visualization(adata, params, context)

    # Get GSEA/ORA results from adata.uns
    gsea_key = getattr(params, "gsea_results_key", "gsea_results")
    if gsea_key not in adata.uns:
        alt_keys = ["rank_genes_groups", "de_results", "pathway_enrichment"]
        for key in alt_keys:
            if key in adata.uns:
                gsea_key = key
                break
        else:
            raise DataNotFoundError(f"GSEA results not found. Expected key: {gsea_key}")

    gsea_results = adata.uns[gsea_key]

    if plot_type == "enrichment_plot":
        return _create_gsea_enrichment_plot(gsea_results, params)
    elif plot_type == "dotplot":
        return _create_gsea_dotplot(gsea_results, params)
    else:  # Default to barplot
        return _create_gsea_barplot(gsea_results, params)


# =============================================================================
# Enrichment Score Visualizations
# =============================================================================


def _create_enrichment_violin(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    score_cols: List[str],
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create violin plot of enrichment scores grouped by cluster."""
    if not params.cluster_key:
        categorical_cols = get_categorical_columns(adata, limit=15)
        raise ParameterError(
            "Enrichment violin plot requires 'cluster_key' parameter.\n"
            f"Available categorical columns: {', '.join(categorical_cols)}"
        )

    validate_obs_column(adata, params.cluster_key, "Grouping variable")

    # Determine scores to plot
    scores_to_plot = _resolve_feature_list(
        params.feature, adata.obs.columns, score_cols
    )
    if not scores_to_plot:
        scores_to_plot = score_cols[:3]

    n_scores = len(scores_to_plot)
    fig, axes = plt.subplots(1, n_scores, figsize=(5 * n_scores, 6))
    if n_scores == 1:
        axes = [axes]

    for i, score in enumerate(scores_to_plot):
        ax = axes[i]
        data = pd.DataFrame(
            {
                params.cluster_key: adata.obs[params.cluster_key],
                "Score": adata.obs[score],
            }
        )
        sns.violinplot(data=data, x=params.cluster_key, y="Score", ax=ax)

        sig_name = score.replace("_score", "")
        ax.set_title(f"{sig_name} by {params.cluster_key}")
        ax.set_xlabel(params.cluster_key)
        ax.set_ylabel("Enrichment Score")
        ax.tick_params(axis="x", rotation=45)
        for label in ax.get_xticklabels():
            label.set_horizontalalignment("right")

    plt.tight_layout()
    return fig


async def _create_enrichment_spatial(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    score_cols: List[str],
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create spatial scatter plot of enrichment scores."""
    feature_list = _resolve_feature_list(params.feature, adata.obs.columns, score_cols)

    if feature_list and len(feature_list) > 1:
        # Multi-score visualization
        scores_to_plot = []
        for feat in feature_list:
            if feat in adata.obs.columns:
                scores_to_plot.append(feat)
            elif f"{feat}_score" in adata.obs.columns:
                scores_to_plot.append(f"{feat}_score")

        if not scores_to_plot:
            raise DataNotFoundError(
                f"None of the specified scores found: {feature_list}"
            )

        fig, axes = setup_multi_panel_figure(
            n_panels=len(scores_to_plot),
            params=params,
            default_title="Enrichment Scores",
        )

        for i, score in enumerate(scores_to_plot):
            if i < len(axes):
                ax = axes[i]
                plot_spatial_feature(adata, feature=score, ax=ax, params=params)
                sig_name = score.replace("_score", "")
                ax.set_title(f"{sig_name} Enrichment")
    else:
        # Single score visualization
        score_col = _resolve_score_column(adata, params.feature, score_cols)
        if context:
            await context.info(f"Using score column: {score_col}")

        fig, ax = create_figure(figsize=(10, 8))
        plot_spatial_feature(adata, feature=score_col, ax=ax, params=params)

        sig_name = score_col.replace("_score", "")
        ax.set_title(f"{sig_name} Enrichment Score", fontsize=14)

        if params.show_colorbar and hasattr(ax, "collections") and ax.collections:
            cbar = plt.colorbar(ax.collections[0], ax=ax)
            cbar.set_label("Enrichment Score", fontsize=12)

    plt.tight_layout()
    return fig


def _create_enrichmap_spatial(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    score_cols: List[str],
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create EnrichMap spatial autocorrelation visualizations."""
    try:
        import enrichmap as em
    except ImportError as e:
        raise ProcessingError(
            f"Spatial enrichment visualization ('{params.subtype}') requires EnrichMap.\n"
            "Install with: pip install enrichmap"
        ) from e

    _ensure_enrichmap_compatibility(adata)
    library_id = adata.obs["library_id"].unique()[0]

    try:
        if params.subtype == "spatial_cross_correlation":
            return _create_enrichmap_cross_correlation(adata, params, library_id, em)
        else:
            return _create_enrichmap_single_score(
                adata, params, library_id, em, context
            )
    except DataNotFoundError:
        raise
    except Exception as e:
        plt.close("all")
        raise ProcessingError(
            f"EnrichMap {params.subtype} visualization failed: {e}\n\n"
            "Solutions:\n"
            "1. Verify the enrichment analysis completed successfully\n"
            "2. Check that spatial neighbors graph exists\n"
            "3. Ensure enrichment scores are stored in adata.obs"
        ) from e


def _create_enrichmap_cross_correlation(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    library_id: str,
    em,
) -> plt.Figure:
    """Create EnrichMap cross-correlation visualization."""
    if "enrichment_gene_sets" not in adata.uns:
        raise DataNotFoundError("enrichment_gene_sets not found in adata.uns")

    pathways = list(adata.uns["enrichment_gene_sets"].keys())
    if len(pathways) < 2:
        raise DataNotFoundError("Need at least 2 pathways for cross-correlation")

    score_x = f"{pathways[0]}_score"
    score_y = f"{pathways[1]}_score"

    em.pl.cross_moran_scatter(
        adata, score_x=score_x, score_y=score_y, library_id=library_id
    )

    fig = plt.gcf()
    if params.figure_size:
        fig.set_size_inches(params.figure_size)
    if params.dpi:
        fig.set_dpi(params.dpi)

    return fig


def _create_enrichmap_single_score(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    library_id: str,
    em,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create single-score EnrichMap visualization."""
    if not params.feature:
        raise DataNotFoundError(
            "Feature parameter required for spatial enrichment visualization"
        )

    score_col = f"{params.feature}_score"
    validate_obs_column(adata, score_col, "Score column")

    if params.subtype == "spatial_correlogram":
        em.pl.morans_correlogram(adata, score_key=score_col, library_id=library_id)
    elif params.subtype == "spatial_variogram":
        em.pl.variogram(adata, score_keys=[score_col])
    elif params.subtype == "spatial_score":
        spot_size = params.spot_size if params.spot_size is not None else 0.5
        em.pl.spatial_enrichmap(
            adata,
            score_key=score_col,
            library_id=library_id,
            cmap="seismic",
            vcenter=0,
            size=spot_size,
            img=False,
        )

    fig = plt.gcf()
    if params.figure_size:
        fig.set_size_inches(params.figure_size)
    if params.dpi:
        fig.set_dpi(params.dpi)

    return fig


# =============================================================================
# GSEA/ORA Pathway Visualizations
# =============================================================================


def _create_gsea_enrichment_plot(
    gsea_results,
    params: VisualizationParameters,
) -> plt.Figure:
    """Create classic GSEA running enrichment score plot.

    Requires full gseapy result object with RES and hits data.
    """
    pathway = params.feature if params.feature else None

    if isinstance(gsea_results, pd.DataFrame):
        raise DataNotFoundError(
            "Enrichment plot requires running enrichment scores (RES) data.\n"
            "The stored results contain only summary statistics.\n\n"
            "Solutions:\n"
            "1. Use subtype='barplot' or subtype='dotplot' instead\n"
            "2. Re-run GSEA analysis and store the full result object"
        )

    if isinstance(gsea_results, dict):
        if pathway and pathway in gsea_results:
            result = gsea_results[pathway]
        else:
            pathway = list(gsea_results.keys())[0]
            result = gsea_results[pathway]

        if not isinstance(result, dict) or "RES" not in result:
            raise DataNotFoundError(
                "Enrichment plot requires 'RES' (running enrichment scores) data.\n"
                "Use subtype='barplot' or subtype='dotplot' instead."
            )

        import gseapy as gp

        figsize = params.figure_size or (6, 5.5)
        fig = gp.gseaplot(
            term=pathway,
            hits=result.get("hits", result.get("hit_indices", [])),
            nes=result.get("NES", result.get("nes", 0)),
            pval=result.get("pval", result.get("NOM p-val", 0)),
            fdr=result.get("fdr", result.get("FDR q-val", 0)),
            RES=result["RES"],
            rank_metric=result.get("rank_metric"),
            figsize=figsize,
            ofname=None,
        )
        return fig

    raise ParameterError(f"Unsupported GSEA results format: {type(gsea_results)}")


def _create_gsea_barplot(
    gsea_results,
    params: VisualizationParameters,
) -> plt.Figure:
    """Create barplot of top enriched pathways."""
    import gseapy as gp

    n_top = getattr(params, "n_top_pathways", 10)
    df = _gsea_results_to_dataframe(gsea_results)

    if df.empty:
        raise DataNotFoundError("No enrichment results found")

    pval_col = _find_pvalue_column(df)
    _ensure_term_column(df)

    figsize = params.figure_size or (6, max(4, n_top * 0.4))
    color = params.colormap if params.colormap != "coolwarm" else "salmon"

    try:
        ax = gp.barplot(
            df=df,
            column=pval_col,
            title=params.title or "Top Enriched Pathways",
            cutoff=1.0,
            top_term=n_top,
            figsize=figsize,
            color=color,
            ofname=None,
        )
        fig = ax.get_figure()
        plt.tight_layout()
        return fig
    except Exception as e:
        raise ProcessingError(
            f"gseapy.barplot failed: {e}\n" f"Available columns: {list(df.columns)}"
        ) from e


def _create_gsea_dotplot(
    gsea_results,
    params: VisualizationParameters,
) -> plt.Figure:
    """Create dotplot of pathway enrichment."""
    import gseapy as gp

    n_top = getattr(params, "n_top_pathways", 10)

    # Handle nested dict (multi-condition)
    if isinstance(gsea_results, dict) and all(
        isinstance(v, dict) for v in gsea_results.values()
    ):
        df, x_col = _nested_dict_to_dataframe(gsea_results)
    else:
        df = _gsea_results_to_dataframe(gsea_results)
        x_col = None

    if df.empty:
        raise DataNotFoundError("No enrichment results found")

    _ensure_term_column(df)
    pval_col = _find_pvalue_column(df)

    figsize = params.figure_size or (6, 8)
    cmap = params.colormap if params.colormap != "coolwarm" else "viridis_r"

    try:
        ax = gp.dotplot(
            df=df,
            column=pval_col,
            x=x_col,
            y="Term",
            title=params.title or "Pathway Enrichment",
            cutoff=1.0,
            top_term=n_top,
            figsize=figsize,
            cmap=cmap,
            size=5,
            ofname=None,
        )
        fig = ax.get_figure()
        plt.tight_layout()
        return fig
    except Exception as e:
        raise ProcessingError(
            f"gseapy.dotplot failed: {e}\n" f"Available columns: {list(df.columns)}"
        ) from e


# =============================================================================
# Utility Functions
# =============================================================================


def _resolve_feature_list(
    feature,
    obs_columns: pd.Index,
    score_cols: List[str],
) -> List[str]:
    """Resolve feature parameter to list of valid score columns."""
    if feature is None:
        return []
    if isinstance(feature, list):
        return feature
    return [feature]


def _gsea_results_to_dataframe(gsea_results) -> pd.DataFrame:
    """Convert GSEA results to DataFrame."""
    if isinstance(gsea_results, pd.DataFrame):
        return gsea_results.copy()
    if isinstance(gsea_results, dict):
        rows = []
        for pathway, data in gsea_results.items():
            if isinstance(data, dict):
                row = {"Term": pathway}
                row.update(data)
                rows.append(row)
        return pd.DataFrame(rows)
    raise ParameterError("Unsupported GSEA results format")


def _nested_dict_to_dataframe(gsea_results: dict):
    """Convert nested dict (multi-condition) to DataFrame with Group column."""
    rows = []
    for condition, pathways in gsea_results.items():
        for pathway, data in pathways.items():
            if isinstance(data, dict):
                row = {"Term": pathway, "Group": condition}
                row.update(data)
                rows.append(row)
    return pd.DataFrame(rows), "Group"


def _find_pvalue_column(df: pd.DataFrame) -> str:
    """Find the p-value column in GSEA results DataFrame."""
    for col in ["Adjusted P-value", "FDR q-val", "fdr", "P-value", "NOM p-val", "pval"]:
        if col in df.columns:
            return col
    return "Adjusted P-value"


def _ensure_term_column(df: pd.DataFrame) -> None:
    """Ensure DataFrame has a 'Term' column."""
    if "Term" in df.columns:
        return
    if "pathway" in df.columns:
        df["Term"] = df["pathway"]
    elif df.index.name or not df.index.equals(pd.RangeIndex(len(df))):
        df["Term"] = df.index
    else:
        raise DataNotFoundError(
            "No pathway/term column found. Expected 'Term' or 'pathway' column."
        )
