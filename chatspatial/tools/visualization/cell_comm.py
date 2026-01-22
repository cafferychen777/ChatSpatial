"""
Cell communication visualization functions for spatial transcriptomics.

Supported methods:
- LIANA+: cluster-based (dotplot, tileplot, circle_plot) and spatial bivariate
- CellPhoneDB: heatmap, dotplot, chord (using ktplotspy)
- FastCCC: heatmap, dotplot (using seaborn - matrix format)
- CellChat R: heatmap, dotplot (using seaborn - matrix format)

Architecture notes:
- CellPhoneDB has rich metadata format with 'interacting_pair' column → ktplotspy
- FastCCC/CellChat R have simple matrix format (index=LR, columns=cell pairs) → seaborn
- Each method uses visualization tools that match its native data format
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...utils.adata_utils import (
    get_cluster_key,
    require_spatial_coords,
    validate_obs_column,
)
from ...utils.dependency_manager import require
from ...utils.exceptions import DataNotFoundError, ParameterError, ProcessingError
from .core import CellCommunicationData, auto_spot_size

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
    - FastCCC analysis results
    - CellChat R analysis results

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

    # Check for FastCCC results
    if "fastccc_interactions_strength" in adata.uns:
        strength = adata.uns["fastccc_interactions_strength"]
        if isinstance(strength, pd.DataFrame):
            # Extract LR pairs from index (format: "LIGAND_RECEPTOR")
            lr_pairs = []
            for pair_str in strength.index:
                parts = str(pair_str).split("_", 1)
                if len(parts) == 2:
                    lr_pairs.append(f"{parts[0]}^{parts[1]}")
                else:
                    lr_pairs.append(str(pair_str))

            if context:
                await context.info(f"Found FastCCC results: {len(lr_pairs)} LR pairs")

            return CellCommunicationData(
                results=strength,
                method="fastccc",
                analysis_type="cluster",
                lr_pairs=lr_pairs,
                results_key="fastccc_interactions_strength",
            )

    # Check for CellChat R results
    if "cellchat_r_lr_pairs" in adata.uns:
        lr_pairs_df = adata.uns["cellchat_r_lr_pairs"]
        if isinstance(lr_pairs_df, pd.DataFrame):
            # Extract LR pairs from DataFrame
            if "interaction_name" in lr_pairs_df.columns:
                lr_pairs = lr_pairs_df["interaction_name"].tolist()
            elif "ligand" in lr_pairs_df.columns and "receptor" in lr_pairs_df.columns:
                lr_pairs = [
                    f"{row['ligand']}^{row['receptor']}"
                    for _, row in lr_pairs_df.iterrows()
                ]
            else:
                lr_pairs = lr_pairs_df.index.tolist()

            if context:
                await context.info(
                    f"Found CellChat R results: {len(lr_pairs)} LR pairs"
                )

            return CellCommunicationData(
                results=lr_pairs_df,
                method="cellchat_r",
                analysis_type="cluster",
                lr_pairs=lr_pairs,
                results_key="cellchat_r_lr_pairs",
            )

    # No results found
    raise DataNotFoundError(
        "No cell communication results found. "
        "Run analyze_cell_communication() first with method='liana', "
        "'cellphonedb', 'fastccc', or 'cellchat_r'."
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

    Routes to appropriate visualization based on analysis method and subtype:
    - Spatial LIANA: Multi-panel spatial LR plot (only 'spatial' subtype supported)
    - CellPhoneDB/FastCCC/CellChat R: heatmap (default), chord, dotplot
    - LIANA cluster: dotplot (default), tileplot, circle_plot

    Note: LIANA spatial bivariate analysis measures LR pair co-localization
    (Moran's I), not cell type pair interactions. For heatmap/dotplot showing
    cell type interactions, use cluster-based analysis (perform_spatial_analysis=False)
    or other methods (cellphonedb, fastccc, cellchat_r).

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
        # LIANA spatial analysis supports specific visualization types
        subtype = params.subtype or "spatial"
        if subtype in ("spatial", None):
            return _create_spatial_lr_visualization(adata, data, params, context)
        else:
            raise ParameterError(
                f"LIANA spatial analysis does not support '{subtype}' visualization.\n\n"
                f"Available for spatial analysis:\n"
                f"  - 'spatial' (default): Multi-panel spatial LR plot showing "
                f"co-localization scores per spot\n\n"
                f"For heatmap/dotplot visualizations showing cell type pair "
                f"interactions:\n"
                f"  1. Re-run analyze_cell_communication with "
                f"perform_spatial_analysis=False\n"
                f"  2. Or use method='cellphonedb', 'fastccc', or 'cellchat_r'\n\n"
                f"Note: Spatial bivariate analysis (Moran's I) measures LR pair "
                f"co-localization patterns, which are best visualized spatially."
            )
    else:
        # CellPhoneDB uses ktplotspy (rich metadata format with 'interacting_pair' column)
        if data.method == "cellphonedb":
            subtype = params.subtype or "heatmap"
            if subtype == "dotplot":
                return _create_cellphonedb_dotplot(adata, data, params, context)
            elif subtype == "chord":
                return _create_cellphonedb_chord(adata, data, params, context)
            elif subtype == "heatmap":
                return _create_cellphonedb_heatmap(adata, data, params, context)
            else:
                raise ParameterError(
                    f"Unknown CellPhoneDB visualization type: {subtype}. "
                    f"Available: heatmap, chord, dotplot"
                )

        # FastCCC and CellChat R use seaborn (simple matrix format)
        if data.method in ("fastccc", "cellchat_r"):
            subtype = params.subtype or "heatmap"
            if subtype == "heatmap":
                return _create_matrix_heatmap(adata, data, params, context)
            elif subtype == "dotplot":
                return _create_matrix_dotplot(adata, data, params, context)
            elif subtype == "chord":
                raise ParameterError(
                    f"Chord diagram not supported for {data.method}. "
                    f"Chord requires CellPhoneDB deconvoluted format. "
                    f"Available for {data.method}: heatmap, dotplot"
                )
            else:
                raise ParameterError(
                    f"Unknown {data.method} visualization type: {subtype}. "
                    f"Available: heatmap, dotplot"
                )
        # LIANA cluster-based results use LIANA+ visualizations
        else:
            subtype = params.subtype or "dotplot"
            if subtype == "tileplot":
                return await _create_liana_tileplot(adata, data, params, context)
            elif subtype == "circle_plot":
                return await _create_liana_circle_plot(adata, data, params, context)
            elif subtype == "dotplot":
                return await _create_cluster_lr_visualization(
                    adata, data, params, context
                )
            else:
                raise ParameterError(
                    f"Unknown LIANA visualization type: {subtype}. "
                    f"Available: dotplot, tileplot, circle_plot"
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

    coords = require_spatial_coords(adata)
    x_coords, y_coords = coords[:, 0], coords[:, 1]

    # Calculate spot size (auto or user-specified)
    spot_size = auto_spot_size(adata, params.spot_size, basis="spatial")

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
            s=spot_size,
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
# Matrix-based Visualizations (FastCCC, CellChat R)
# =============================================================================


def _create_matrix_heatmap(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create heatmap for matrix-format cell communication results.

    Used for FastCCC and CellChat R which store results as simple matrices
    (index=LR pairs, columns=cell type pairs).

    Args:
        adata: AnnData object
        data: CellCommunicationData with method and results
        params: Visualization parameters
        context: Optional logging context

    Returns:
        matplotlib Figure
    """
    import seaborn as sns

    # Get p-values matrix based on method
    pvalues_key = {
        "fastccc": "fastccc_pvalues",
        "cellchat_r": "cellchat_r_pval",
    }.get(data.method)

    if pvalues_key is None:
        raise ParameterError(f"Unsupported method for matrix heatmap: {data.method}")

    pvalues = adata.uns.get(pvalues_key)

    if pvalues is None:
        raise DataNotFoundError(
            f"{data.method} pvalues not found in adata.uns['{pvalues_key}']. "
            f"Re-run {data.method} analysis."
        )

    # Convert to DataFrame if numpy array (CellChat R case)
    if isinstance(pvalues, np.ndarray):
        # Try to get labels from related data
        strength_key = {
            "fastccc": "fastccc_interactions_strength",
            "cellchat_r": "cellchat_r_prob",
        }.get(data.method)
        strength = adata.uns.get(strength_key)

        if isinstance(strength, pd.DataFrame):
            pvalues = pd.DataFrame(
                pvalues, index=strength.index, columns=strength.columns
            )
        else:
            # Use generic labels
            pvalues = pd.DataFrame(pvalues)

    if not isinstance(pvalues, pd.DataFrame) or pvalues.empty:
        raise DataNotFoundError(f"{data.method} pvalues DataFrame is empty or invalid.")

    # Select top LR pairs for visualization
    n_top_rows = min(params.plot_top_pairs or 30, len(pvalues))

    # Sort by minimum p-value (most significant)
    numeric_cols = pvalues.select_dtypes(include=[np.number]).columns

    # Also limit the number of columns (cell type pairs) to prevent huge figures
    max_cols = 50  # Reasonable maximum for visualization
    if len(numeric_cols) > max_cols:
        # Select columns with most significant interactions
        col_min_pvals = pvalues[numeric_cols].min(axis=0)
        numeric_cols = col_min_pvals.nsmallest(max_cols).index

    if len(numeric_cols) > 0:
        min_pvals = pvalues[numeric_cols].min(axis=1)
        top_idx = min_pvals.nsmallest(n_top_rows).index
        pvalues_top = pvalues.loc[top_idx, numeric_cols]
    else:
        pvalues_top = pvalues.iloc[:n_top_rows]

    # Transform to -log10(pvalue) for better visualization
    neg_log_pval = -np.log10(pvalues_top.astype(float) + 1e-10)

    # Create figure with reasonable size limits
    n_rows, n_cols = neg_log_pval.shape
    # Calculate ideal size but cap at reasonable limits
    fig_width = min(20, max(8, n_cols * 0.6 + 3))  # Max 20 inches wide
    fig_height = min(16, max(6, n_rows * 0.3 + 2))  # Max 16 inches tall
    figsize = params.figure_size or (fig_width, fig_height)

    # Ensure DPI is reasonable (cap at 300 for figures)
    effective_dpi = min(params.dpi or 300, 300)

    fig, ax = plt.subplots(figsize=figsize, dpi=effective_dpi)

    # Create heatmap
    cmap = params.colormap if params.colormap != "coolwarm" else "Reds"
    sns.heatmap(
        neg_log_pval,
        ax=ax,
        cmap=cmap,
        annot=n_rows <= 20 and n_cols <= 10,  # Show values if small enough
        fmt=".1f",
        linewidths=0.5,
        cbar_kws={"label": "-log10(p-value)"},
        vmin=params.vmin,
        vmax=params.vmax,
    )

    # Method-specific title
    method_name = {"fastccc": "FastCCC", "cellchat_r": "CellChat"}.get(
        data.method, data.method
    )
    ax.set_title(
        params.title
        or f"{method_name}: Interaction Significance\n(top {n_top_rows} pairs)",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_xlabel("Cell Type Pairs", fontsize=10)
    ax.set_ylabel("Ligand-Receptor Pairs", fontsize=10)

    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=8)

    plt.tight_layout()
    return fig


def _create_matrix_dotplot(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create dotplot for matrix-format cell communication results.

    Following best practices from sc-best-practices.org:
    - Dot SIZE represents SIGNIFICANCE (-log10 p-value, larger = more significant)
    - Dot COLOR represents STRENGTH/MAGNITUDE (interaction intensity)

    This matches the standard visualization approach used by LIANA, CellPhoneDB,
    and other cell communication tools.

    Args:
        adata: AnnData object
        data: CellCommunicationData with method and results
        params: Visualization parameters
        context: Optional logging context

    Returns:
        matplotlib Figure
    """
    # Get data based on method
    data_keys = {
        "fastccc": {
            "strength": "fastccc_interactions_strength",
            "pvalues": "fastccc_pvalues",
        },
        "cellchat_r": {
            "strength": "cellchat_r_prob",
            "pvalues": "cellchat_r_pval",
        },
    }.get(data.method)

    if data_keys is None:
        raise ParameterError(f"Unsupported method for matrix dotplot: {data.method}")

    strength = adata.uns.get(data_keys["strength"])
    pvalues = adata.uns.get(data_keys["pvalues"])

    if strength is None or pvalues is None:
        raise DataNotFoundError(
            f"{data.method} results not found. Re-run {data.method} analysis."
        )

    # Convert to DataFrame if numpy array
    if isinstance(strength, np.ndarray):
        strength = pd.DataFrame(strength)
    if isinstance(pvalues, np.ndarray):
        if isinstance(strength, pd.DataFrame):
            pvalues = pd.DataFrame(
                pvalues, index=strength.index, columns=strength.columns
            )
        else:
            pvalues = pd.DataFrame(pvalues)

    if not isinstance(strength, pd.DataFrame) or strength.empty:
        raise DataNotFoundError(f"{data.method} results are empty or invalid.")

    # Select top LR pairs (rows)
    n_top_rows = min(params.plot_top_pairs or 20, len(strength))
    numeric_cols = strength.select_dtypes(include=[np.number]).columns

    if len(numeric_cols) == 0:
        raise DataNotFoundError(f"{data.method} results have no numeric columns.")

    # Also limit the number of columns (cell type pairs) to prevent huge figures
    max_cols = 50  # Reasonable maximum for visualization
    if len(numeric_cols) > max_cols:
        # Select columns with highest mean interaction strength
        col_mean_strength = strength[numeric_cols].mean(axis=0)
        numeric_cols = col_mean_strength.nlargest(max_cols).index

    # Sort by mean strength
    mean_strength = strength[numeric_cols].mean(axis=1)
    top_idx = mean_strength.nlargest(n_top_rows).index

    strength_top = strength.loc[top_idx, numeric_cols]
    pvalues_top = pvalues.loc[top_idx, numeric_cols]

    # Prepare data for dotplot (long format)
    plot_data = []
    for lr_pair in strength_top.index:
        for cell_pair in numeric_cols:
            s = strength_top.loc[lr_pair, cell_pair]
            p = pvalues_top.loc[lr_pair, cell_pair]
            plot_data.append(
                {
                    "LR_pair": str(lr_pair),
                    "Cell_pair": str(cell_pair),
                    "strength": float(s) if pd.notna(s) else 0,
                    "neg_log_pval": -np.log10(float(p) + 1e-10) if pd.notna(p) else 0,
                }
            )

    plot_df = pd.DataFrame(plot_data)

    # Create figure with reasonable size limits
    n_rows = len(top_idx)
    n_cols_plot = len(numeric_cols)
    # Calculate ideal size but cap at reasonable limits
    fig_width = min(20, max(8, n_cols_plot * 0.8 + 3))  # Max 20 inches wide
    fig_height = min(16, max(6, n_rows * 0.4 + 2))  # Max 16 inches tall
    figsize = params.figure_size or (fig_width, fig_height)

    # Ensure DPI is reasonable (cap at 300 for figures)
    effective_dpi = min(params.dpi or 300, 300)

    fig, ax = plt.subplots(figsize=figsize, dpi=effective_dpi)

    # Best practice: SIZE = significance (-log10 p-value), COLOR = strength
    # Larger dots = more significant (smaller p-value)
    size_norm = plot_df["neg_log_pval"].max()
    if size_norm > 0:
        sizes = (plot_df["neg_log_pval"] / size_norm * 200) + 10
    else:
        sizes = 50

    scatter = ax.scatter(
        x=pd.Categorical(plot_df["Cell_pair"]).codes,
        y=pd.Categorical(plot_df["LR_pair"]).codes,
        s=sizes,
        c=plot_df["strength"],
        cmap="viridis",  # viridis for continuous strength values
        alpha=params.alpha,
        edgecolors="black",
        linewidths=0.5,
    )

    # Colorbar for strength
    plt.colorbar(scatter, ax=ax, label="Interaction Strength")

    # Labels
    unique_cell_pairs = plot_df["Cell_pair"].unique()
    unique_lr_pairs = plot_df["LR_pair"].unique()

    ax.set_xticks(range(len(unique_cell_pairs)))
    ax.set_xticklabels(unique_cell_pairs, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(unique_lr_pairs)))
    ax.set_yticklabels(unique_lr_pairs, fontsize=8)

    # Method-specific title
    method_name = {"fastccc": "FastCCC", "cellchat_r": "CellChat"}.get(
        data.method, data.method
    )
    ax.set_title(
        params.title
        or f"{method_name}: L-R Interactions\n(size=significance, color=strength)",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_xlabel("Cell Type Pairs", fontsize=10)
    ax.set_ylabel("Ligand-Receptor Pairs", fontsize=10)

    plt.tight_layout()
    return fig


# =============================================================================
# CellPhoneDB Visualizations (ktplotspy)
# =============================================================================


def _create_cellphonedb_heatmap(
    adata: "ad.AnnData",
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CellPhoneDB heatmap visualization using ktplotspy.

    Uses ktplotspy which requires CellPhoneDB-specific format with
    'interacting_pair' column and metadata.
    """
    import ktplotspy as kpy

    means = data.results

    if not isinstance(means, pd.DataFrame) or len(means) == 0:
        raise DataNotFoundError("CellPhoneDB results empty. Re-run analysis.")

    pvalues = adata.uns.get("cellphonedb_pvalues")

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
    """Create CellPhoneDB dotplot visualization using ktplotspy.

    Uses ktplotspy which requires CellPhoneDB-specific format with
    'interacting_pair' column and metadata.
    """
    means = data.results

    if not isinstance(means, pd.DataFrame) or len(means) == 0:
        raise DataNotFoundError("CellPhoneDB results empty. Re-run analysis.")

    require("ktplotspy", feature="CellPhoneDB dotplot visualization")
    import ktplotspy as kpy

    try:
        pvalues = adata.uns.get("cellphonedb_pvalues")

        if pvalues is None or not isinstance(pvalues, pd.DataFrame):
            raise DataNotFoundError("Missing pvalues DataFrame for CellPhoneDB dotplot")

        cluster_key = params.cluster_key or get_cluster_key(adata)
        if not cluster_key:
            raise ParameterError(
                "cluster_key required for CellPhoneDB dotplot. "
                "No default cluster key found in data."
            )
        validate_obs_column(adata, cluster_key, "Cluster")

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
            f"Failed to create CellPhoneDB dotplot: {e}\n\n"
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
        pvalues = adata.uns.get("cellphonedb_pvalues")
        deconvoluted = adata.uns.get("cellphonedb_deconvoluted")

        if pvalues is None or not isinstance(pvalues, pd.DataFrame):
            raise DataNotFoundError(
                "Missing pvalues DataFrame for ktplotspy chord plot"
            )

        if deconvoluted is None or not isinstance(deconvoluted, pd.DataFrame):
            raise DataNotFoundError(
                "Missing deconvoluted DataFrame for chord plot. "
                "Re-run CellPhoneDB analysis."
            )

        cluster_key = params.cluster_key or get_cluster_key(adata)
        if not cluster_key:
            raise ParameterError(
                "cluster_key required for CellPhoneDB chord plot. "
                "No default cluster key found in data."
            )
        validate_obs_column(adata, cluster_key, "Cluster")

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
            f"Failed to create CellPhoneDB chord diagram: {e}\n\n"
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
