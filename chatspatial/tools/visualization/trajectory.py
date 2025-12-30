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
from .core import get_categorical_columns, infer_basis, setup_multi_panel_figure

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
    # Find pseudotime key
    pseudotime_key = params.feature
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

    validate_obs_column(adata, pseudotime_key, "Pseudotime column")

    # Check if RNA velocity is available
    has_velocity = "velocity_graph" in adata.uns

    # Determine basis for plotting
    basis = infer_basis(adata, preferred=params.basis)
    if not basis:
        raise DataCompatibilityError(
            f"No valid embedding basis found. "
            f"Available keys: {list(adata.obsm.keys())}"
        )

    # Setup figure: 1 panel if no velocity, 2 panels if velocity exists
    n_panels = 2 if has_velocity else 1

    fig, axes = setup_multi_panel_figure(
        n_panels=n_panels,
        params=params,
        default_title=f"Trajectory Analysis - Pseudotime ({pseudotime_key})",
    )

    # Panel 1: Pseudotime plot
    ax1 = axes[0]
    try:
        sc.pl.embedding(
            adata,
            basis=basis,
            color=pseudotime_key,
            cmap=params.colormap,
            ax=ax1,
            show=False,
            frameon=params.show_axes,
            alpha=params.alpha,
            colorbar_loc="right" if params.show_colorbar else None,
        )

        if basis == "spatial":
            ax1.invert_yaxis()

    except Exception as e:
        ax1.text(
            0.5,
            0.5,
            f"Error plotting pseudotime:\n{str(e)}",
            ha="center",
            va="center",
            transform=ax1.transAxes,
        )
        ax1.set_title("Pseudotime (Error)", fontsize=12)

    # Panel 2: Velocity stream plot (if available)
    if has_velocity and n_panels > 1:
        ax2 = axes[1]
        try:
            import scvelo as scv

            scv.pl.velocity_embedding_stream(
                adata,
                basis=basis,
                color=pseudotime_key,
                cmap=params.colormap,
                ax=ax2,
                show=False,
                alpha=params.alpha,
                frameon=params.show_axes,
            )
            ax2.set_title("RNA Velocity Stream", fontsize=12)

            if basis == "spatial":
                ax2.invert_yaxis()

        except ImportError:
            ax2.text(
                0.5,
                0.5,
                "scvelo not installed",
                ha="center",
                va="center",
                transform=ax2.transAxes,
            )
        except Exception as e:
            ax2.text(
                0.5,
                0.5,
                f"Error: {str(e)[:50]}",
                ha="center",
                va="center",
                transform=ax2.transAxes,
            )

    plt.tight_layout(rect=(0, 0, 1, 0.95))
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
    require("cellrank", feature="circular projection")
    import cellrank as cr

    # Check for CellRank results
    fate_key_candidates = ["lineages_fwd", "to_terminal_states"]
    fate_key = None
    for key in fate_key_candidates:
        if key in adata.obsm:
            fate_key = key
            break

    if not fate_key:
        raise DataNotFoundError(
            "CellRank fate probabilities not found. Run trajectory analysis first."
        )

    if context:
        await context.info("Creating CellRank circular projection")

    # Determine keys for coloring
    keys = [params.cluster_key] if params.cluster_key else None
    if not keys:
        categorical_cols = get_categorical_columns(adata, limit=3)
        keys = categorical_cols if categorical_cols else None

    figsize = params.figure_size or (10, 10)

    import matplotlib

    backend = matplotlib.get_backend()
    matplotlib.use("Agg")

    cr.pl.circular_projection(
        adata,
        keys=keys,
        figsize=figsize,
        dpi=params.dpi,
    )
    fig = plt.gcf()
    matplotlib.use(backend)

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

    Shows fate probabilities aggregated by cluster as bar, paga, or heatmap.

    Data requirements:
        - adata.obsm['lineages_fwd'] or 'to_terminal_states': Fate probabilities
        - adata.obs[cluster_key]: Cluster labels for aggregation
    """
    require("cellrank", feature="fate map")
    import cellrank as cr

    # Check for CellRank results
    fate_key_candidates = ["lineages_fwd", "to_terminal_states"]
    fate_key = None
    for key in fate_key_candidates:
        if key in adata.obsm:
            fate_key = key
            break

    if not fate_key:
        raise DataNotFoundError(
            "CellRank fate probabilities not found. Run trajectory analysis first."
        )

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

    if context:
        await context.info(f"Creating CellRank fate map for '{cluster_key}'")

    figsize = params.figure_size or (12, 6)

    import matplotlib

    backend = matplotlib.get_backend()
    matplotlib.use("Agg")

    cr.pl.aggregate_fate_probabilities(
        adata,
        cluster_key=cluster_key,
        mode="bar",
        figsize=figsize,
        dpi=params.dpi,
    )
    fig = plt.gcf()
    matplotlib.use(backend)

    title = params.title or f"CellRank Fate Probabilities by {cluster_key}"
    fig.suptitle(title, fontsize=14, y=1.02)

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
    require("cellrank", feature="gene trends")
    import cellrank as cr

    # Import GAM model preparation from trajectory module
    from ..trajectory import prepare_gam_model_for_visualization

    # Check for fate probabilities
    fate_key_candidates = ["lineages_fwd", "to_terminal_states"]
    fate_key = None
    for key in fate_key_candidates:
        if key in adata.obsm:
            fate_key = key
            break

    if not fate_key:
        raise DataNotFoundError(
            "CellRank fate probabilities not found. Run trajectory analysis first."
        )

    # Find time key
    time_key = None
    time_candidates = ["latent_time", "palantir_pseudotime", "dpt_pseudotime"]
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

    figsize = params.figure_size or (12, 3 * len(genes))

    import matplotlib

    backend = matplotlib.get_backend()
    matplotlib.use("Agg")

    model, lineage_names = prepare_gam_model_for_visualization(
        adata, genes, time_key=time_key, fate_key=fate_key
    )

    if context:
        await context.info(f"Lineages: {lineage_names}")

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
    matplotlib.use(backend)

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
    require("cellrank", feature="fate heatmap")
    import cellrank as cr

    # Import GAM model preparation from trajectory module
    from ..trajectory import prepare_gam_model_for_visualization

    # Check for fate probabilities
    fate_key_candidates = ["lineages_fwd", "to_terminal_states"]
    fate_key = None
    for key in fate_key_candidates:
        if key in adata.obsm:
            fate_key = key
            break

    if not fate_key:
        raise DataNotFoundError(
            "CellRank fate probabilities not found. Run trajectory analysis first."
        )

    # Find time key
    time_key = None
    time_candidates = ["latent_time", "palantir_pseudotime", "dpt_pseudotime"]
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

    figsize = params.figure_size or (12, 10)

    import matplotlib

    backend = matplotlib.get_backend()
    matplotlib.use("Agg")

    model, lineage_names = prepare_gam_model_for_visualization(
        adata, genes, time_key=time_key, fate_key=fate_key
    )

    if context:
        await context.info(f"Lineages: {lineage_names}")

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
    matplotlib.use(backend)

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
    fate_key = None
    for key in ["palantir_fate_probs", "palantir_branch_probs"]:
        if key in adata.obsm:
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

    # Determine number of panels
    n_panels = 1 + int(has_entropy) + (1 if fate_key else 0)

    # Create figure
    figsize = params.figure_size or (5 * n_panels, 5)
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
        dominant_fate = fate_probs.argmax(axis=1)
        adata.obs["_dominant_fate"] = dominant_fate.astype(str)

        sc.pl.embedding(
            adata,
            basis=basis,
            color="_dominant_fate",
            ax=ax,
            show=False,
            frameon=params.show_axes,
            title="Dominant Fate",
        )
        if basis == "spatial":
            ax.invert_yaxis()

        # Clean up temporary column
        del adata.obs["_dominant_fate"]

    title = params.title or "Palantir Trajectory Analysis"
    fig.suptitle(title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig
