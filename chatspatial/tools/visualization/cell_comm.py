"""
Cell communication visualization functions for spatial transcriptomics.

Architecture (First Principles):
    All CCC methods' results are converted to LIANA format for unified visualization.
    This enables using LIANA's native visualization functions (dotplot, tileplot)
    regardless of the analysis method used.

LIANA Format (canonical):
    DataFrame with columns:
    - source: Source cell type
    - target: Target cell type
    - ligand_complex: Ligand gene name
    - receptor_complex: Receptor gene name
    - lr_means: Expression score (higher = stronger)
    - magnitude_rank: Significance rank (lower = more significant)

Visualization routing:
    - dotplot: Use LIANA native li.pl.dotplot (magnitude + specificity)
    - tileplot: Use LIANA native li.pl.tileplot (ligand/receptor stats)
    - circle_plot: Custom network visualization (cell type connectivity)
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...tools.cell_communication import (
    CCC_SPATIAL_PVALS_KEY,
    CCC_SPATIAL_SCORES_KEY,
    get_ccc_results,
)
from ...utils.adata_utils import require_spatial_coords
from ...utils.dependency_manager import require
from ...utils.exceptions import DataNotFoundError, ParameterError, ProcessingError
from .core import CellCommunicationData, auto_spot_size

# =============================================================================
# Data Retrieval and Format Conversion
# =============================================================================


def _convert_to_liana_format(
    results: pd.DataFrame,
    pvalues: Optional[pd.DataFrame],
    method: str,
    method_data: Optional[dict] = None,
) -> pd.DataFrame:
    """Convert non-LIANA results to LIANA-compatible DataFrame format.

    This is the core function for unified visualization. All CCC methods'
    results are converted to LIANA format so we can use LIANA's native
    visualization functions (li.pl.dotplot, li.pl.tileplot).

    LIANA format required columns:
        - source: Source cell type
        - target: Target cell type
        - ligand_complex: Ligand gene name
        - receptor_complex: Receptor gene name
        - lr_means: Expression score (higher = stronger)
        - magnitude_rank: Significance rank (0-1, lower = more significant)

    Args:
        results: Original results DataFrame (format varies by method)
        pvalues: P-values DataFrame (optional, used for ranking)
        method: Analysis method name
        method_data: Method-specific data (required for CellChat R 3D matrices)

    Returns:
        LIANA-compatible DataFrame
    """
    if results is None or len(results) == 0:
        return pd.DataFrame()

    # LIANA results are already in correct format
    if method == "liana":
        return results

    # Check if already in LIANA format (has required columns)
    liana_cols = {"source", "target", "ligand_complex", "receptor_complex"}
    if liana_cols.issubset(results.columns):
        return results

    # CellChat R stores results in 3D matrices, requires special handling
    if method == "cellchat_r" and method_data is not None:
        return _cellchat_3d_to_liana_format(results, method_data)

    # Convert matrix-format results (CellPhoneDB, FastCCC)
    # Matrix format: rows = LR pairs, columns = cell type pairs (e.g., "CellA|CellB")
    return _matrix_to_liana_format(results, pvalues, method)


def _cellchat_3d_to_liana_format(
    results: pd.DataFrame,
    method_data: dict,
) -> pd.DataFrame:
    """Convert CellChat R 3D matrices to LIANA long format.

    CellChat R stores data as:
        - results: LR pairs metadata (interaction_name, ligand, receptor, etc.)
        - method_data['prob_matrix']: 3D array (n_sources × n_targets × n_lr_pairs)
        - method_data['pval_matrix']: 3D array (same shape)
        - method_data['cell_type_names']: Cell type labels for axes

    Args:
        results: LR pairs metadata DataFrame with 'interaction_name', 'ligand', 'receptor'
        method_data: Dict containing prob_matrix, pval_matrix, cell_type_names

    Returns:
        LIANA-format DataFrame
    """
    prob_matrix = method_data.get("prob_matrix")
    pval_matrix = method_data.get("pval_matrix")
    cell_type_names = method_data.get("cell_type_names")

    if prob_matrix is None or cell_type_names is None:
        return pd.DataFrame()

    # Ensure cell_type_names is a list of strings
    if hasattr(cell_type_names, "tolist"):
        cell_type_names = cell_type_names.tolist()
    cell_type_names = [str(x) for x in cell_type_names]

    n_sources, n_targets, n_lr_pairs = prob_matrix.shape

    # Get LR pair info from results DataFrame
    lr_info = {}
    if "interaction_name" in results.columns:
        for idx, row in results.iterrows():
            name = str(row["interaction_name"])
            ligand = str(row.get("ligand", name.split("_")[0]))
            receptor = str(row.get("receptor", "_".join(name.split("_")[1:])))
            lr_info[idx] = (name, ligand, receptor)

    rows = []
    for lr_idx in range(n_lr_pairs):
        # Get ligand/receptor from metadata
        if lr_idx in lr_info:
            lr_name, ligand, receptor = lr_info[lr_idx]
        else:
            lr_name = f"LR_{lr_idx}"
            ligand = lr_name
            receptor = lr_name

        for src_idx, source in enumerate(cell_type_names):
            for tgt_idx, target in enumerate(cell_type_names):
                prob = prob_matrix[src_idx, tgt_idx, lr_idx]
                if prob == 0 or np.isnan(prob):
                    continue  # Skip zero/nan interactions

                pval = 1.0
                if pval_matrix is not None:
                    pval = pval_matrix[src_idx, tgt_idx, lr_idx]
                    if np.isnan(pval):
                        pval = 1.0

                rows.append(
                    {
                        "source": source,
                        "target": target,
                        "ligand_complex": ligand,
                        "receptor_complex": receptor,
                        "lr_means": float(prob),
                        "magnitude_rank": float(pval),
                    }
                )

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    # Normalize magnitude_rank to 0-1 if needed (p-values already 0-1)
    if len(df) > 0 and df["magnitude_rank"].max() > 1:
        max_rank = df["magnitude_rank"].max()
        df["magnitude_rank"] = df["magnitude_rank"] / max_rank

    return df


def _matrix_to_liana_format(
    results: pd.DataFrame,
    pvalues: Optional[pd.DataFrame],
    method: str,
) -> pd.DataFrame:
    """Convert matrix-format CCC results to LIANA long format.

    Matrix format (CellPhoneDB, FastCCC, CellChat R):
        - Index: LR pair names (e.g., "LIGAND_RECEPTOR" or "LIGAND^RECEPTOR")
        - Columns: Cell type pairs (e.g., "CellA|CellB")
        - Values: Interaction strength/means

    Args:
        results: Matrix-format results DataFrame
        pvalues: P-values DataFrame (same shape as results)
        method: Analysis method name

    Returns:
        LIANA-format DataFrame
    """
    rows = []

    # Get numeric columns (cell type pairs)
    numeric_cols = results.select_dtypes(include=[np.number]).columns.tolist()

    # Also handle CellPhoneDB which has metadata columns
    # Filter to columns that look like cell type pairs (contain "|")
    cell_pair_cols = [c for c in numeric_cols if "|" in str(c)]
    if not cell_pair_cols:
        # Fallback: use all numeric columns
        cell_pair_cols = numeric_cols

    if not cell_pair_cols:
        return pd.DataFrame()

    # Get LR pair column or use index
    lr_col = None
    for col in ["interacting_pair", "interaction_name", "lr_pair"]:
        if col in results.columns:
            lr_col = col
            break

    for idx, row in results.iterrows():
        # Get LR pair name
        if lr_col:
            lr_pair = str(row[lr_col])
        else:
            lr_pair = str(idx)

        # Parse ligand and receptor from LR pair
        # Handle various separators: "_", "^", "-"
        ligand, receptor = _parse_lr_pair(lr_pair)

        for cell_pair in cell_pair_cols:
            # Parse source and target from cell pair
            # Format: "CellA|CellB" or "CellA_CellB"
            if "|" in str(cell_pair):
                parts = str(cell_pair).split("|")
            elif "_" in str(cell_pair):
                parts = str(cell_pair).split("_", 1)
            else:
                continue

            if len(parts) != 2:
                continue

            source, target = parts[0], parts[1]

            # Get expression value
            lr_means = row[cell_pair]
            if pd.isna(lr_means):
                continue

            # Get p-value for ranking
            pval = 1.0
            if pvalues is not None and cell_pair in pvalues.columns:
                try:
                    pval = pvalues.loc[idx, cell_pair]
                    if pd.isna(pval):
                        pval = 1.0
                except (KeyError, TypeError):
                    pval = 1.0

            rows.append(
                {
                    "source": source,
                    "target": target,
                    "ligand_complex": ligand,
                    "receptor_complex": receptor,
                    "lr_means": float(lr_means),
                    "magnitude_rank": float(
                        pval
                    ),  # Use p-value as rank (lower = better)
                }
            )

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    # Normalize magnitude_rank to 0-1 range if not already
    if len(df) > 0 and df["magnitude_rank"].max() > 1:
        max_rank = df["magnitude_rank"].max()
        df["magnitude_rank"] = df["magnitude_rank"] / max_rank

    return df


def _parse_lr_pair(lr_pair: str) -> tuple[str, str]:
    """Parse ligand and receptor from LR pair string.

    Handles various formats:
        - "LIGAND_RECEPTOR"
        - "LIGAND^RECEPTOR"
        - "LIGAND-RECEPTOR"
        - "complex:LIGAND_RECEPTOR" (CellPhoneDB)

    Returns:
        Tuple of (ligand, receptor)
    """
    # Remove prefix if present (e.g., "complex:")
    if ":" in lr_pair:
        lr_pair = lr_pair.split(":")[-1]

    # Try different separators
    for sep in ["^", "_", "-"]:
        if sep in lr_pair:
            parts = lr_pair.split(sep, 1)
            if len(parts) == 2:
                return parts[0].strip(), parts[1].strip()

    # Fallback: return as-is for both
    return lr_pair, lr_pair


async def get_cell_communication_data(
    adata: "ad.AnnData",
    method: Optional[str] = None,
    context: Optional["ToolContext"] = None,
) -> CellCommunicationData:
    """
    Unified function to retrieve cell communication results from AnnData.

    Reads from unified CCC storage at adata.uns["ccc"] (single source of truth).
    Converts all results to LIANA format for unified visualization.

    Args:
        adata: AnnData object with cell communication results
        method: Analysis method hint (optional, unused - kept for API compatibility)
        context: MCP context for logging

    Returns:
        CellCommunicationData with LIANA-format results for visualization

    Raises:
        DataNotFoundError: No cell communication results found
    """
    ccc = get_ccc_results(adata)

    if ccc is None:
        raise DataNotFoundError(
            "No cell communication results found. "
            "Run analyze_cell_communication() first with method='liana', "
            "'cellphonedb', 'fastccc', or 'cellchat_r'."
        )

    # Get spatial data from obsm if available
    spatial_scores = adata.obsm.get(CCC_SPATIAL_SCORES_KEY)
    spatial_pvals = adata.obsm.get(CCC_SPATIAL_PVALS_KEY)

    # Convert results to LIANA format for unified visualization
    # This is the key step for first-principles architecture:
    # all methods use the same visualization code path
    liana_results = _convert_to_liana_format(
        ccc.results,
        ccc.pvalues,
        ccc.method,
        ccc.method_data,  # Required for CellChat R 3D matrix conversion
    )

    # Extract source/target labels from converted results
    source_labels = None
    target_labels = None
    if len(liana_results) > 0:
        if "source" in liana_results.columns:
            source_labels = liana_results["source"].unique().tolist()
        if "target" in liana_results.columns:
            target_labels = liana_results["target"].unique().tolist()

    if context:
        await context.info(
            f"Found {ccc.method} results ({ccc.analysis_type} analysis, "
            f"{ccc.n_pairs} LR pairs) → converted to LIANA format"
        )

    return CellCommunicationData(
        results=liana_results,
        method=ccc.method,
        analysis_type=ccc.analysis_type,
        lr_pairs=ccc.lr_pairs,
        pvalues=ccc.pvalues,  # Keep original pvalues for methods that need them
        spatial_scores=spatial_scores,
        spatial_pvals=spatial_pvals,
        source_labels=source_labels,
        target_labels=target_labels,
        method_data=ccc.method_data,  # Keep original data for method-specific viz
    )


# =============================================================================
# Main Router
# =============================================================================


async def create_cell_communication_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create cell communication visualization from unified CCC storage.

    Architecture (First Principles):
        All CCC results are converted to LIANA format in get_cell_communication_data().
        This enables unified visualization regardless of the analysis method used.
        One data format → One set of visualization functions.

    Unified visualizations (all methods):
        - dotplot: LIANA native li.pl.dotplot (default for cluster analysis)
        - tileplot: LIANA native li.pl.tileplot
        - circle_plot: Cell type network diagram

    Spatial-only (LIANA li.mt.bivariate):
        - spatial: Per-spot LR scores on tissue

    Args:
        adata: AnnData object with cell communication results
        params: Visualization parameters
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

    # Determine subtype with method-aware defaults
    subtype = params.subtype
    if subtype is None:
        if data.analysis_type == "spatial":
            subtype = "spatial"
        else:
            subtype = "dotplot"

    # Route based on subtype (unified approach)
    # Data is already in LIANA format, so most visualizations work uniformly

    # Spatial visualization (LIANA spatial analysis only)
    if subtype == "spatial":
        if data.analysis_type != "spatial":
            raise ParameterError(
                f"Spatial visualization requires spatial analysis.\n"
                f"Current analysis type: {data.analysis_type}\n\n"
                f"Re-run with method='liana' and perform_spatial_analysis=True"
            )
        return _create_spatial_lr_visualization(adata, data, params, context)

    # Unified visualizations (all methods use same code path)
    # Data is already in LIANA format, enabling unified visualization
    if subtype == "dotplot":
        return _create_unified_dotplot(data, params, context)

    if subtype == "tileplot":
        return _create_unified_tileplot(data, params, context)

    if subtype == "circle_plot":
        return _create_unified_circle_plot(data, params, context)

    # Unknown subtype
    available = ["dotplot", "tileplot", "circle_plot"]
    if data.analysis_type == "spatial":
        available.insert(0, "spatial")

    raise ParameterError(
        f"Unknown visualization type: {subtype}. " f"Available: {', '.join(available)}"
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


def _create_unified_dotplot(
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create unified dotplot using LIANA's native visualization.

    Uses li.pl.dotplot() which works with any LIANA-format DataFrame.
    All CCC methods' results are converted to LIANA format in get_cell_communication_data().

    Best practice visualization:
    - Size = significance (smaller magnitude_rank = more significant = larger dot)
    - Color = expression (lr_means)
    """
    require("liana", feature="Cell communication dotplot")
    import liana as li

    df = data.results
    if df is None or len(df) == 0:
        raise DataNotFoundError(
            f"No {data.method} results found. Run analyze_cell_communication() first."
        )

    # Validate required columns
    required_cols = {"source", "target", "ligand_complex", "receptor_complex"}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise DataNotFoundError(f"Missing required columns: {missing}")

    # Determine columns for visualization
    size_col = "magnitude_rank" if "magnitude_rank" in df.columns else None
    color_col = "lr_means" if "lr_means" in df.columns else size_col

    if size_col is None and color_col is None:
        raise DataNotFoundError("No suitable columns for visualization")

    # Create dotplot using LIANA native function
    n_top = params.plot_top_pairs or 20
    figsize = params.figure_size or (12, 10)

    try:
        p = li.pl.dotplot(
            liana_res=df,
            colour=color_col,
            size=size_col,
            inverse_size=True,  # Smaller rank = larger dot
            top_n=n_top,
            orderby=size_col,
            orderby_ascending=True,
            cmap=params.colormap or "viridis",
            figure_size=figsize,
            return_fig=True,
        )

        # Convert plotnine to matplotlib
        return _plotnine_to_matplotlib(p, params)

    except Exception as e:
        # Fallback to custom implementation if LIANA fails
        return _create_fallback_dotplot(data, params, context, error=str(e))


def _create_fallback_dotplot(
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
    error: Optional[str] = None,
) -> plt.Figure:
    """Fallback dotplot using matplotlib when LIANA visualization fails."""
    df = data.results

    # Determine columns
    rank_col = "magnitude_rank" if "magnitude_rank" in df.columns else None
    expr_col = "lr_means" if "lr_means" in df.columns else rank_col

    if rank_col is None:
        raise DataNotFoundError("No ranking column found in results")

    # Select top interactions
    n_top = params.plot_top_pairs or 20
    top_df = df.nsmallest(n_top, rank_col).copy()

    # Create labels
    top_df["lr_label"] = top_df["ligand_complex"] + " → " + top_df["receptor_complex"]
    top_df["cell_pair"] = top_df["source"] + " → " + top_df["target"]

    # Calculate dot sizes
    max_rank = top_df[rank_col].max()
    if max_rank > 0:
        top_df["dot_size"] = (1 - top_df[rank_col] / max_rank) * 200 + 20
    else:
        top_df["dot_size"] = 100

    # Create figure
    figsize = params.figure_size or (12, 10)
    fig, ax = plt.subplots(figsize=figsize)

    lr_labels = top_df["lr_label"].unique()
    cell_pairs = top_df["cell_pair"].unique()

    lr_map = {label: i for i, label in enumerate(lr_labels)}
    cell_map = {pair: i for i, pair in enumerate(cell_pairs)}

    top_df["lr_pos"] = top_df["lr_label"].map(lr_map)
    top_df["cell_pos"] = top_df["cell_pair"].map(cell_map)

    scatter = ax.scatter(
        top_df["cell_pos"],
        top_df["lr_pos"],
        s=top_df["dot_size"],
        c=top_df[expr_col] if expr_col else "steelblue",
        cmap=params.colormap or "viridis",
        alpha=params.alpha or 0.8,
        edgecolors="black",
        linewidths=0.5,
    )

    if expr_col:
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(expr_col, fontsize=10)

    ax.set_xticks(range(len(cell_pairs)))
    ax.set_xticklabels(cell_pairs, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(lr_labels)))
    ax.set_yticklabels(lr_labels, fontsize=9)

    ax.set_xlabel("Cell Type Pairs", fontsize=11)
    ax.set_ylabel("Ligand-Receptor Pairs", fontsize=11)

    method_name = data.method.upper() if data.method else "CCC"
    ax.set_title(
        params.title or f"{method_name}: Top {n_top} Interactions",
        fontsize=12,
        fontweight="bold",
    )

    plt.tight_layout()
    return fig


def _create_unified_tileplot(
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create unified tileplot using LIANA's native visualization.

    Uses li.pl.tileplot() which works with any LIANA-format DataFrame.
    """
    require("liana", feature="Cell communication tileplot")
    import liana as li

    df = data.results
    if df is None or len(df) == 0:
        raise DataNotFoundError(
            f"No {data.method} results found. Run analyze_cell_communication() first."
        )

    # Determine columns for visualization
    fill_col = "magnitude_rank" if "magnitude_rank" in df.columns else "lr_means"
    label_col = "lr_means" if "lr_means" in df.columns else fill_col

    if fill_col not in df.columns:
        raise DataNotFoundError("No suitable columns for tileplot visualization")

    n_top = params.plot_top_pairs or 20
    figsize = params.figure_size or (10, 8)

    try:
        p = li.pl.tileplot(
            liana_res=df,
            fill=fill_col,
            label=label_col,
            top_n=n_top,
            orderby=fill_col,
            orderby_ascending=True,
            cmap=params.colormap or "viridis",
            figure_size=figsize,
            return_fig=True,
        )

        return _plotnine_to_matplotlib(p, params)

    except Exception as e:
        # Fallback to seaborn heatmap
        return _create_fallback_tileplot(data, params, context, error=str(e))


def _create_fallback_tileplot(
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
    error: Optional[str] = None,
) -> plt.Figure:
    """Fallback tileplot using seaborn when LIANA visualization fails."""
    import seaborn as sns

    df = data.results

    # Determine value column
    value_col = None
    for col in ["magnitude_rank", "lr_means"]:
        if col in df.columns:
            value_col = col
            break

    if value_col is None:
        raise DataNotFoundError("No suitable value column found")

    df = df.copy()
    df["lr_label"] = df["ligand_complex"] + "_" + df["receptor_complex"]
    df["cell_pair"] = df["source"] + " → " + df["target"]

    n_top = params.plot_top_pairs or 20
    top_df = df.nsmallest(n_top, value_col)

    pivot = top_df.pivot_table(
        index="lr_label", columns="cell_pair", values=value_col, aggfunc="mean"
    )

    figsize = params.figure_size or (14, 10)
    fig, ax = plt.subplots(figsize=figsize)

    cmap = params.colormap if params.colormap != "coolwarm" else "viridis_r"
    sns.heatmap(
        pivot,
        ax=ax,
        cmap=cmap,
        annot=len(pivot) <= 15 and len(pivot.columns) <= 10,
        fmt=".2f",
        linewidths=0.5,
        cbar_kws={"label": value_col},
    )

    method_name = data.method.upper() if data.method else "CCC"
    ax.set_title(
        params.title or f"{method_name}: Interaction Scores",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_xlabel("Cell Type Pairs", fontsize=10)
    ax.set_ylabel("Ligand-Receptor Pairs", fontsize=10)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)
    plt.setp(ax.get_yticklabels(), fontsize=9)

    plt.tight_layout()
    return fig


def _create_unified_circle_plot(
    data: CellCommunicationData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create unified circle plot (network diagram) for all CCC methods.

    Shows cell-cell communication as a network with cell types as nodes
    and interaction strength as edges. Works with all methods since data
    is already in LIANA format.

    Note: Uses custom implementation rather than LIANA's circle_plot because
    the latter requires specific adata setup that doesn't work with our
    unified storage architecture.
    """
    df = data.results
    if df is None or len(df) == 0:
        raise DataNotFoundError(
            f"No {data.method} results found. Run analyze_cell_communication() first."
        )

    if "source" not in df.columns or "target" not in df.columns:
        raise DataNotFoundError("Missing source/target columns in results")

    # Determine score column
    score_col = None
    for col in ["magnitude_rank", "specificity_rank", "lr_means"]:
        if col in df.columns:
            score_col = col
            break

    # Aggregate interactions by cell type pairs
    if score_col:
        # Use inverse rank as weight (lower rank = stronger interaction)
        df = df.copy()
        max_val = df[score_col].max()
        df["weight"] = max_val - df[score_col] + 1
        agg = df.groupby(["source", "target"])["weight"].sum().reset_index()
    else:
        # Count interactions
        agg = df.groupby(["source", "target"]).size().reset_index(name="weight")

    # Get unique cell types
    cell_types = list(set(agg["source"].tolist() + agg["target"].tolist()))
    n_types = len(cell_types)

    # Create figure
    figsize = params.figure_size or (10, 10)
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": "polar"})

    # Calculate positions on circle
    angles = np.linspace(0, 2 * np.pi, n_types, endpoint=False)
    type_angles = {ct: angles[i] for i, ct in enumerate(cell_types)}

    # Draw cell type nodes
    radius = 0.9
    for ct, angle in type_angles.items():
        ax.scatter(angle, radius, s=500, zorder=5, c="steelblue", edgecolors="white")
        ax.annotate(
            ct,
            xy=(angle, radius + 0.15),
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
        )

    # Draw edges (simplified chord representation using bezier-like curves)
    max_weight = agg["weight"].max()
    for _, row in agg.iterrows():
        src_angle = type_angles[row["source"]]
        tgt_angle = type_angles[row["target"]]
        weight = row["weight"]

        # Line width proportional to weight
        lw = (weight / max_weight) * 5 + 0.5

        # Draw arc connecting source and target
        if src_angle != tgt_angle:
            # Draw a curved line through center
            mid_angle = (src_angle + tgt_angle) / 2
            mid_radius = 0.3  # Curve through center region

            ax.plot(
                [src_angle, mid_angle, tgt_angle],
                [radius * 0.85, mid_radius, radius * 0.85],
                lw=lw,
                alpha=0.6,
                c="coral",
            )
        else:
            # Self-loop (autocrine)
            loop_angles = np.linspace(src_angle - 0.2, src_angle + 0.2, 20)
            loop_radii = radius * 0.85 + 0.1 * np.sin(np.linspace(0, np.pi, 20))
            ax.plot(loop_angles, loop_radii, lw=lw, alpha=0.6, c="coral")

    # Clean up polar plot
    ax.set_ylim(0, 1.2)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines["polar"].set_visible(False)

    method_name = data.method.upper() if data.method else "CCC"
    ax.set_title(
        params.title or f"{method_name}: Cell-Cell Communication Network",
        fontsize=12,
        fontweight="bold",
        y=1.08,
    )

    return fig


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
