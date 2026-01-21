"""
Batch integration visualization functions.

Subtypes:
- batch: UMAP colored by batch (default, assess mixing)
- cluster: UMAP colored by cluster (assess biological conservation)
- highlight: Per-batch highlighted UMAP panels
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np

from ...models.data import VisualizationParameters
from ...utils.adata_utils import validate_obs_column
from ...utils.exceptions import DataNotFoundError, ParameterError
from .core import get_categorical_cmap

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# Main Entry Point
# =============================================================================


async def create_batch_integration_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create batch integration quality visualization.

    Subtypes:
        - batch: UMAP colored by batch (default) - assess mixing quality
        - cluster: UMAP colored by cluster - assess biological conservation
        - highlight: Per-batch highlighted panels - detailed batch distribution

    Args:
        adata: AnnData object with integrated samples
        params: Visualization parameters (batch_key required)
        context: MCP context for logging

    Returns:
        matplotlib Figure object
    """
    subtype = params.subtype or "batch"

    if subtype == "batch":
        return await _create_umap_by_batch(adata, params, context)
    elif subtype == "cluster":
        return await _create_umap_by_cluster(adata, params, context)
    elif subtype == "highlight":
        return await _create_batch_highlight(adata, params, context)
    else:
        raise ParameterError(
            f"Unknown integration subtype: {subtype}. "
            f"Available: batch, cluster, highlight"
        )


# =============================================================================
# Subtype: batch - UMAP colored by batch
# =============================================================================


async def _create_umap_by_batch(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create UMAP colored by batch to assess mixing quality.

    Good integration shows mixed colors (batches intermingled).
    Poor integration shows separated clusters by batch.
    """
    if context:
        await context.info("Creating UMAP colored by batch")

    # Validate requirements
    batch_key = params.batch_key
    validate_obs_column(adata, batch_key, "Batch")

    if "X_umap" not in adata.obsm:
        raise DataNotFoundError(
            "UMAP coordinates not found. Run compute_embeddings first."
        )

    # Setup figure
    figsize = params.figure_size or (8, 6)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    umap_coords = adata.obsm["X_umap"]
    batch_values = adata.obs[batch_key]
    unique_batches = (
        batch_values.cat.categories
        if hasattr(batch_values, "cat")
        else batch_values.unique()
    )

    # Get colors
    cmap_name = get_categorical_cmap(len(unique_batches))
    cmap = plt.cm.get_cmap(cmap_name)
    colors = [
        cmap(i / max(1, len(unique_batches) - 1)) for i in range(len(unique_batches))
    ]

    # Plot each batch
    for i, batch in enumerate(unique_batches):
        mask = batch_values == batch
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=[colors[i]],
            label=str(batch),
            s=params.spot_size or 10,
            alpha=params.alpha,
            rasterized=True,
        )

    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(
        params.title or "UMAP by Batch\n(Good integration = mixed colors)",
        fontsize=14,
    )

    # Legend
    if params.show_legend:
        ax.legend(
            title="Batch",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
        )

    ax.set_aspect("equal", adjustable="datalim")
    plt.tight_layout()
    return fig


# =============================================================================
# Subtype: cluster - UMAP colored by cluster
# =============================================================================


async def _create_umap_by_cluster(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create UMAP colored by cluster to assess biological conservation.

    Good integration preserves biological structure (clear cell type clusters).
    """
    if context:
        await context.info("Creating UMAP colored by cluster")

    if "X_umap" not in adata.obsm:
        raise DataNotFoundError(
            "UMAP coordinates not found. Run compute_embeddings first."
        )

    # Find cluster key
    cluster_key = params.cluster_key
    if not cluster_key:
        # Auto-detect cluster column
        for key in ["leiden", "louvain", "cluster", "cell_type", "celltype"]:
            if key in adata.obs.columns:
                cluster_key = key
                break

    if not cluster_key or cluster_key not in adata.obs.columns:
        raise DataNotFoundError(
            f"Cluster column not found. Available: {list(adata.obs.columns)[:10]}..."
        )

    # Setup figure
    figsize = params.figure_size or (8, 6)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    umap_coords = adata.obsm["X_umap"]
    cluster_values = adata.obs[cluster_key]
    unique_clusters = (
        cluster_values.cat.categories
        if hasattr(cluster_values, "cat")
        else cluster_values.unique()
    )

    # Get colors
    cmap_name = get_categorical_cmap(len(unique_clusters))
    cmap = plt.cm.get_cmap(cmap_name)
    colors = [
        cmap(i / max(1, len(unique_clusters) - 1)) for i in range(len(unique_clusters))
    ]

    # Plot each cluster
    for i, cluster in enumerate(unique_clusters):
        mask = cluster_values == cluster
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=[colors[i]],
            label=str(cluster),
            s=params.spot_size or 10,
            alpha=params.alpha,
            rasterized=True,
        )

    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title(
        params.title
        or f"UMAP by {cluster_key}\n(Good integration = preserved clusters)",
        fontsize=14,
    )

    # Legend
    if params.show_legend:
        ncol = 1 if len(unique_clusters) <= 10 else 2
        ax.legend(
            title=cluster_key,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            ncol=ncol,
        )

    ax.set_aspect("equal", adjustable="datalim")
    plt.tight_layout()
    return fig


# =============================================================================
# Subtype: highlight - Per-batch highlighted panels
# =============================================================================


async def _create_batch_highlight(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create multi-panel UMAP with each batch highlighted separately.

    Shows detailed distribution of each batch across the embedding.
    """
    if context:
        await context.info("Creating per-batch highlight visualization")

    # Validate requirements
    batch_key = params.batch_key
    validate_obs_column(adata, batch_key, "Batch")

    if "X_umap" not in adata.obsm:
        raise DataNotFoundError(
            "UMAP coordinates not found. Run compute_embeddings first."
        )

    umap_coords = adata.obsm["X_umap"]
    batch_values = adata.obs[batch_key]
    unique_batches = (
        batch_values.cat.categories
        if hasattr(batch_values, "cat")
        else batch_values.unique()
    )
    n_batches = len(unique_batches)

    # Calculate grid layout
    n_cols = min(4, n_batches)
    n_rows = (n_batches + n_cols - 1) // n_cols

    # Setup figure
    if params.figure_size:
        figsize = params.figure_size
    else:
        figsize = (4 * n_cols, 3.5 * n_rows)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, dpi=params.dpi)
    axes = np.atleast_2d(axes).flatten()

    # Get colors
    cmap_name = get_categorical_cmap(n_batches)
    cmap = plt.cm.get_cmap(cmap_name)
    colors = [cmap(i / max(1, n_batches - 1)) for i in range(n_batches)]

    # Plot each batch
    for i, batch in enumerate(unique_batches):
        ax = axes[i]
        mask = batch_values == batch

        # Background: all other cells in gray
        ax.scatter(
            umap_coords[~mask, 0],
            umap_coords[~mask, 1],
            c="lightgray",
            s=3,
            alpha=0.3,
            rasterized=True,
        )

        # Foreground: highlighted batch
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=[colors[i]],
            s=params.spot_size or 8,
            alpha=params.alpha,
            rasterized=True,
        )

        n_cells = mask.sum()
        ax.set_title(f"{batch}\n(n={n_cells:,})", fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal", adjustable="datalim")

    # Hide unused axes
    for i in range(n_batches, len(axes)):
        axes[i].axis("off")

    fig.suptitle(
        params.title or "Per-Batch Distribution",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout()
    return fig
