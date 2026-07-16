"""
Trajectory visualization functions for spatial transcriptomics.

This module contains:
- Pseudotime visualizations
- CellRank circular projections
- Fate map visualizations
- Gene trends along lineages
- Fate heatmaps
- Palantir results visualization
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

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
from ...utils.image_utils import non_interactive_backend
from .core import (
    get_categorical_columns,
    infer_basis,
    resolve_figure_size,
)

# CellRank 2.x default forward fate key
_CELLRANK_FATE_KEY = "to_terminal_states"

# Candidate obsm keys for CellRank fate probabilities (priority order)
_FATE_KEY_CANDIDATES = ["lineages_fwd", "to_terminal_states", "fate_probabilities"]


def _resolve_cellrank_fate_key(adata: "ad.AnnData") -> str:
    """Find the CellRank fate-probabilities key in *adata.obsm* (read-only).

    The analysis pipeline now writes both ``"fate_probabilities"`` and the
    CellRank-standard ``"to_terminal_states"`` alias, so this function only
    needs to locate the first matching key without modifying *adata*.

    Returns:
        The obsm key that holds the fate probabilities.

    Raises:
        DataNotFoundError: If no fate probabilities are found.
    """
    for key in _FATE_KEY_CANDIDATES:
        if key in adata.obsm:
            return key

    raise DataNotFoundError(
        "CellRank fate probabilities not found. Run trajectory analysis first."
    )


# =============================================================================
# Main Router
# =============================================================================


async def create_trajectory_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create trajectory visualization based on subtype.

    Dispatcher function that routes to appropriate trajectory visualization.

    Args:
        adata: AnnData object with computed trajectory/pseudotime
        params: Visualization parameters including subtype
        context: MCP context

    Returns:
        Matplotlib figure with trajectory visualization

    Subtypes:
        - pseudotime (default): Pseudotime on embedding with optional velocity stream
        - circular: CellRank circular projection of fate probabilities
        - fate_map: CellRank aggregated fate probabilities (bar/paga/heatmap)
        - gene_trends: CellRank gene expression trends along lineages
        - fate_heatmap: CellRank smoothed expression heatmap by pseudotime
        - palantir: Palantir comprehensive results (pseudotime, entropy, fate probs)
    """
    subtype = params.subtype or "pseudotime"

    if context:
        await context.info(f"Creating trajectory visualization (subtype: {subtype})")

    if subtype == "pseudotime":
        return await _create_trajectory_pseudotime_plot(adata, params, context)
    elif subtype == "circular":
        return await _create_cellrank_circular_projection(adata, params, context)
    elif subtype == "fate_map":
        return await _create_cellrank_fate_map(adata, params, context)
    elif subtype == "gene_trends":
        return await _create_cellrank_gene_trends(adata, params, context)
    elif subtype == "fate_heatmap":
        return await _create_cellrank_fate_heatmap(adata, params, context)
    elif subtype == "palantir":
        return await _create_palantir_results(adata, params, context)
    else:
        raise ParameterError(
            f"Unsupported subtype for trajectory: '{subtype}'. "
            f"Available subtypes: pseudotime, circular, fate_map, gene_trends, "
            f"fate_heatmap, palantir"
        )


# =============================================================================
# Visualization Functions
# =============================================================================


async def _create_trajectory_pseudotime_plot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create trajectory pseudotime visualization.

    Shows pseudotime on embedding with optional velocity stream plot.

    Data requirements:
        - adata.obs['*pseudotime*']: Any pseudotime column
        - adata.obsm['X_umap'] or 'spatial': Embedding for visualization
        - adata.uns['velocity_graph']: Optional, for velocity stream panel
    """
    # Find pseudotime key (normalize list to single feature)
    pseudotime_key: str | None = None
    if params.feature is not None:
        pseudotime_key = (
            params.feature[0] if isinstance(params.feature, list) else params.feature
        )
    if not pseudotime_key:
        pseudotime_candidates = [
            k for k in adata.obs.columns if "pseudotime" in k.lower()
        ]
        if pseudotime_candidates:
            pseudotime_key = pseudotime_candidates[0]
            if context:
                await context.info(f"Found pseudotime column: {pseudotime_key}")
        else:
            raise DataNotFoundError(
                "No pseudotime found. Run trajectory analysis first."
            )

    if pseudotime_key is None:
        raise DataNotFoundError("No pseudotime found. Run trajectory analysis first.")
    validate_obs_column(adata, pseudotime_key, "Pseudotime")

    # Check if RNA velocity graph is a real transition matrix (not a bool placeholder)
    import scipy.sparse

    vg = adata.uns.get("velocity_graph")
    has_velocity = vg is not None and (
        scipy.sparse.issparse(vg) or isinstance(vg, np.ndarray)
    )

    # Determine basis for plotting
    basis = infer_basis(adata, preferred=params.basis)
    if not basis:
        raise DataCompatibilityError(
            f"No valid embedding basis found. Available keys: {list(adata.obsm.keys())}"
        )

    # Setup figure: 1 panel if no velocity, 2 panels if velocity exists
    n_panels = 2 if has_velocity else 1

    # Use explicit figure setup to avoid suptitle overlap with subplot titles
    figsize = resolve_figure_size(
        params, n_panels=n_panels, panel_width=6, panel_height=5
    )
    fig, axes = plt.subplots(1, n_panels, figsize=figsize, dpi=params.dpi)
    axes = [axes] if n_panels == 1 else list(axes)

    # Get colormap for consistent colorbars
    from matplotlib import colormaps
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    cmap = colormaps.get_cmap(params.colormap)

    # Panel 1: Pseudotime plot (disable auto colorbar for manual control)
    ax1 = axes[0]
    sc.pl.embedding(
        adata,
        basis=basis,
        color=pseudotime_key,
        cmap=params.colormap,
        ax=ax1,
        show=False,
        frameon=params.show_axes,
        alpha=params.alpha,
        colorbar_loc=None,  # Disable auto colorbar
        title=f"Pseudotime ({pseudotime_key})",
    )

    if basis == "spatial":
        ax1.invert_yaxis()

    # Add consistent colorbar for panel 1 (use actual data range)
    if params.show_colorbar:
        pt_vals = adata.obs[pseudotime_key].dropna()
        vmin, vmax = (
            (float(pt_vals.min()), float(pt_vals.max())) if len(pt_vals) else (0, 1)
        )
        divider1 = make_axes_locatable(ax1)
        cax1 = divider1.append_axes("right", size="4%", pad=0.05)
        sm1 = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))
        sm1.set_array([])
        fig.colorbar(sm1, cax=cax1)

    # Panel 2: Velocity stream plot (if available)
    if has_velocity and n_panels > 1:
        ax2 = axes[1]
        scv = require(
            "scvelo",
            context,
            feature="RNA velocity stream visualization",
        )

        # Note: scvelo uses 'color_map' not 'cmap', disable auto colorbar
        scv.pl.velocity_embedding_stream(
            adata,
            basis=basis,
            color=pseudotime_key,
            color_map=params.colormap,
            ax=ax2,
            show=False,
            alpha=params.alpha,
            frameon=params.show_axes,
            title="RNA Velocity Stream",
            colorbar=False,  # Disable scvelo's inset colorbar
        )

        if basis == "spatial":
            ax2.invert_yaxis()

        # Add consistent colorbar for panel 2 (same range as panel 1)
        if params.show_colorbar:
            divider2 = make_axes_locatable(ax2)
            cax2 = divider2.append_axes("right", size="4%", pad=0.05)
            sm2 = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))
            sm2.set_array([])
            fig.colorbar(sm2, cax=cax2)

    # Only add suptitle if explicitly provided by user
    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_cellrank_circular_projection(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CellRank circular projection using cr.pl.circular_projection.

    Shows fate probabilities in a circular layout.

    Data requirements:
        - adata.obs['terminal_states'] or 'term_states_fwd': Terminal state labels
        - adata.obsm['lineages_fwd'] or 'to_terminal_states': Fate probabilities
    """
    cr = require("cellrank", feature="circular projection")

    _resolve_cellrank_fate_key(adata)

    if context:
        await context.info("Creating CellRank circular projection")

    # Determine keys for coloring
    if params.cluster_key:
        validate_obs_column(adata, params.cluster_key, "Cluster")
    keys = [params.cluster_key] if params.cluster_key else None
    if not keys:
        categorical_cols = get_categorical_columns(adata, limit=3)
        keys = categorical_cols if categorical_cols else None

    # Use centralized figure size resolution
    figsize = resolve_figure_size(params, "trajectory")

    with non_interactive_backend():
        cr.pl.circular_projection(
            adata,
            keys=keys,
            figsize=figsize,
            dpi=params.dpi,
        )
        fig = plt.gcf()

    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_cellrank_fate_map(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CellRank aggregated fate probabilities.

    Shows fate probabilities aggregated by cluster using CellRank's default
    mode (paga_pie).

    Data requirements:
        - adata.obsm with CellRank fate probabilities (any standard key)
        - adata.obs[cluster_key]: Cluster labels for aggregation
    """
    cr = require("cellrank", feature="fate map")

    _resolve_cellrank_fate_key(adata)

    # Determine cluster key
    cluster_key = params.cluster_key
    if not cluster_key:
        categorical_cols = get_categorical_columns(adata)
        if categorical_cols:
            cluster_key = categorical_cols[0]
            if context:
                await context.info(f"Using cluster_key: '{cluster_key}'")
        else:
            raise ParameterError("cluster_key is required for fate map visualization.")
    validate_obs_column(adata, cluster_key, "Cluster")

    if context:
        await context.info(f"Creating CellRank fate map for '{cluster_key}'")

    # Use centralized figure size resolution
    figsize = resolve_figure_size(params, "violin")  # similar width to violin plots

    with non_interactive_backend():
        cr.pl.aggregate_fate_probabilities(
            adata,
            cluster_key=cluster_key,
            mode="heatmap",
            figsize=figsize,
            dpi=params.dpi,
        )
        fig = plt.gcf()

    # Only add suptitle if explicitly provided to avoid overlap with CellRank's title
    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_cellrank_gene_trends(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CellRank gene expression trends using cr.pl.gene_trends.

    Shows gene expression trends along lineages/pseudotime.

    Data requirements:
        - adata.obsm['lineages_fwd'] or 'to_terminal_states': Fate probabilities
        - adata.obs['latent_time'] or similar pseudotime
        - Gene expression in adata.X
    """
    cr = require("cellrank", feature="gene trends")

    # Import GAM model preparation from trajectory module
    from ..trajectory import prepare_gam_model_for_visualization

    fate_key = _resolve_cellrank_fate_key(adata)

    # Find time key
    time_key = None
    time_candidates = [
        "latent_time",
        "palantir_pseudotime",
        "dpt_pseudotime",
        "pseudotime",
    ]
    for key in time_candidates:
        if key in adata.obs.columns:
            time_key = key
            break

    if not time_key:
        raise DataNotFoundError("No pseudotime found. Run trajectory analysis first.")

    # Get genes to plot
    if params.feature:
        if isinstance(params.feature, str):
            genes = [params.feature]
        else:
            genes = list(params.feature)
        valid_genes = [g for g in genes if g in adata.var_names]
        if not valid_genes:
            raise DataNotFoundError(f"None of the specified genes found: {genes}")
        genes = valid_genes[:6]
    else:
        if "highly_variable" in adata.var.columns:
            hvg = adata.var_names[adata.var["highly_variable"]]
            genes = list(hvg[:6])
        else:
            genes = list(adata.var_names[:6])

    if context:
        await context.info(f"Creating gene trends for: {genes}")

    # Use centralized figure size resolution with dynamic panel height
    figsize = resolve_figure_size(
        params, n_panels=len(genes), panel_width=12, panel_height=3
    )

    model, lineage_names = prepare_gam_model_for_visualization(
        adata, genes, time_key=time_key, fate_key=fate_key
    )

    if context:
        await context.info(f"Lineages: {lineage_names}")

    with non_interactive_backend():
        cr.pl.gene_trends(
            adata,
            model=model,
            genes=genes,
            time_key=time_key,
            figsize=figsize,
            n_jobs=1,
            show_progress_bar=False,
        )
        fig = plt.gcf()
        fig.set_dpi(params.dpi)

    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_cellrank_fate_heatmap(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CellRank fate heatmap using cr.pl.heatmap.

    Shows smoothed gene expression ordered by pseudotime per lineage.

    Data requirements:
        - adata.obsm['lineages_fwd'] or 'to_terminal_states': Fate probabilities
        - adata.obs['latent_time'] or similar pseudotime
        - Gene expression in adata.X
    """
    cr = require("cellrank", feature="fate heatmap")

    # Import GAM model preparation from trajectory module
    from ..trajectory import prepare_gam_model_for_visualization

    fate_key = _resolve_cellrank_fate_key(adata)

    # Find time key
    time_key = None
    time_candidates = [
        "latent_time",
        "palantir_pseudotime",
        "dpt_pseudotime",
        "pseudotime",
    ]
    for key in time_candidates:
        if key in adata.obs.columns:
            time_key = key
            break

    if not time_key:
        raise DataNotFoundError("No pseudotime found for fate heatmap.")

    # Get genes
    if params.feature:
        if isinstance(params.feature, str):
            genes = [params.feature]
        else:
            genes = list(params.feature)
        valid_genes = [g for g in genes if g in adata.var_names]
        if not valid_genes:
            raise DataNotFoundError(f"None of the genes found: {genes}")
        genes = valid_genes[:50]
    else:
        if "highly_variable" in adata.var.columns:
            hvg = adata.var_names[adata.var["highly_variable"]]
            genes = list(hvg[:50])
        else:
            genes = list(adata.var_names[:50])

    if context:
        await context.info(f"Creating fate heatmap with {len(genes)} genes")

    # Use centralized figure size resolution
    figsize = resolve_figure_size(params, "heatmap")

    model, lineage_names = prepare_gam_model_for_visualization(
        adata, genes, time_key=time_key, fate_key=fate_key
    )

    if context:
        await context.info(f"Lineages: {lineage_names}")

    with non_interactive_backend():
        cr.pl.heatmap(
            adata,
            model=model,
            genes=genes,
            time_key=time_key,
            figsize=figsize,
            n_jobs=1,
            show_progress_bar=False,
        )
        fig = plt.gcf()
        fig.set_dpi(params.dpi)

    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_palantir_results(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create Palantir comprehensive results visualization.

    Shows pseudotime, entropy, and fate probabilities in a multi-panel figure.

    Data requirements:
        - adata.obs['palantir_pseudotime']: Pseudotime
        - adata.obs['palantir_entropy']: Differentiation entropy
        - adata.obsm['palantir_fate_probs'] or 'palantir_branch_probs': Fate probabilities
    """
    # Check for Palantir results
    has_pseudotime = "palantir_pseudotime" in adata.obs.columns
    has_entropy = "palantir_entropy" in adata.obs.columns

    # Check for fate probabilities and ensure they're not empty
    fate_key = None
    for key in ["palantir_fate_probs", "palantir_branch_probs"]:
        if key in adata.obsm:
            fate_data = adata.obsm[key]
            # Check if fate data has any columns (not empty)
            if hasattr(fate_data, "shape") and fate_data.shape[1] > 0:
                fate_key = key
                break

    if not has_pseudotime:
        raise DataNotFoundError(
            "Palantir results not found. Run trajectory analysis first."
        )

    if context:
        await context.info("Creating Palantir results visualization")

    # Determine basis
    basis = infer_basis(
        adata, preferred=params.basis, priority=["umap", "spatial", "pca"]
    )

    # Determine number of panels (only include fate panel if fate data is non-empty)
    n_panels = 1 + int(has_entropy) + (1 if fate_key else 0)

    # Create figure with centralized utility
    figsize = resolve_figure_size(
        params, n_panels=n_panels, panel_width=5, panel_height=5
    )
    fig, axes = plt.subplots(1, n_panels, figsize=figsize, dpi=params.dpi)
    if n_panels == 1:
        axes = [axes]

    panel_idx = 0

    # Panel 1: Pseudotime
    ax = axes[panel_idx]
    sc.pl.embedding(
        adata,
        basis=basis,
        color="palantir_pseudotime",
        cmap="viridis",
        ax=ax,
        show=False,
        frameon=params.show_axes,
        title="Palantir Pseudotime",
    )
    if basis == "spatial":
        ax.invert_yaxis()
    panel_idx += 1

    # Panel 2: Entropy (if available)
    if has_entropy and panel_idx < n_panels:
        ax = axes[panel_idx]
        sc.pl.embedding(
            adata,
            basis=basis,
            color="palantir_entropy",
            cmap="magma",
            ax=ax,
            show=False,
            frameon=params.show_axes,
            title="Differentiation Entropy",
        )
        if basis == "spatial":
            ax.invert_yaxis()
        panel_idx += 1

    # Panel 3: Fate probabilities summary (if available)
    if fate_key and panel_idx < n_panels:
        ax = axes[panel_idx]
        fate_probs = adata.obsm[fate_key]

        # Handle both DataFrame (Palantir) and ndarray (CellRank) formats
        import pandas as pd

        if isinstance(fate_probs, pd.DataFrame):
            # Palantir returns DataFrame - use idxmax
            # All-NaN rows produce NaN from idxmax; label as "unassigned"
            nan_rows = fate_probs.isna().all(axis=1)
            dominant_fate = fate_probs.idxmax(axis=1).astype(object)
            if nan_rows.any():
                dominant_fate[nan_rows] = "unassigned"
        else:
            # CellRank returns ndarray - use argmax
            dominant_fate = fate_probs.argmax(axis=1)

            # NaN rows (failed computation) should not be assigned a fate
            nan_rows = np.isnan(fate_probs).any(axis=1)
            if nan_rows.any():
                dominant_fate = dominant_fate.astype(object)
                dominant_fate[nan_rows] = "unassigned"

        _fate_temp_key = "_dominant_fate"
        adata.obs[_fate_temp_key] = dominant_fate.astype(str)
        try:
            sc.pl.embedding(
                adata,
                basis=basis,
                color=_fate_temp_key,
                ax=ax,
                show=False,
                frameon=params.show_axes,
                title="Dominant Fate",
            )
            if basis == "spatial":
                ax.invert_yaxis()
        finally:
            if _fate_temp_key in adata.obs.columns:
                del adata.obs[_fate_temp_key]

    title = params.title or "Palantir Trajectory Analysis"
    fig.suptitle(title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig
