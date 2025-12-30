"""
Batch integration visualization functions.

This module contains:
- Batch integration quality assessment visualization
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import entropy

from ...models.data import VisualizationParameters
from ...utils.adata_utils import get_spatial_key, validate_obs_column

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# Batch Integration Visualization
# =============================================================================


async def create_batch_integration_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create multi-panel visualization to assess batch integration quality.

    This visualization is specifically for evaluating the quality of batch correction
    after integrating multiple samples. It requires proper batch information.

    Args:
        adata: AnnData object with integrated samples
        params: Visualization parameters (batch_key required)
        context: MCP context for logging

    Returns:
        matplotlib Figure object

    Raises:
        DataNotFoundError: If batch information not found
    """
    if context:
        await context.info("Creating batch integration quality visualization")

    # Validate batch key exists
    batch_key = params.batch_key
    validate_obs_column(adata, batch_key, "Batch key")

    # Create multi-panel figure (2x2 layout)
    figsize = params.figure_size if params.figure_size else (16, 12)
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    batch_values = adata.obs[batch_key]
    unique_batches = batch_values.unique()
    colors = plt.get_cmap("Set3")(np.linspace(0, 1, len(unique_batches)))

    # Panel 1: UMAP colored by batch (shows mixing)
    if "X_umap" in adata.obsm:
        umap_coords = adata.obsm["X_umap"]

        for i, batch in enumerate(unique_batches):
            mask = batch_values == batch
            axes[0, 0].scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                c=[colors[i]],
                label=f"{batch}",
                s=5,
                alpha=0.7,
            )

        axes[0, 0].set_title(
            "UMAP colored by batch\n(Good integration = mixed colors)", fontsize=12
        )
        axes[0, 0].set_xlabel("UMAP 1")
        axes[0, 0].set_ylabel("UMAP 2")
        axes[0, 0].legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        if context:
            await context.warning(
                "UMAP coordinates not available. "
                "Run preprocessing with UMAP computation for complete visualization."
            )
        axes[0, 0].text(
            0.5, 0.5, "UMAP coordinates not available", ha="center", va="center"
        )
        axes[0, 0].set_title("UMAP (Not Available)", fontsize=12)

    # Panel 2: Spatial plot colored by batch (if spatial data available)
    spatial_key = get_spatial_key(adata)
    if spatial_key:
        spatial_coords = adata.obsm[spatial_key]

        for i, batch in enumerate(unique_batches):
            mask = batch_values == batch
            axes[0, 1].scatter(
                spatial_coords[mask, 0],
                spatial_coords[mask, 1],
                c=[colors[i]],
                label=f"{batch}",
                s=10,
                alpha=0.7,
            )

        axes[0, 1].set_title("Spatial coordinates colored by batch", fontsize=12)
        axes[0, 1].set_xlabel("Spatial X")
        axes[0, 1].set_ylabel("Spatial Y")
        axes[0, 1].set_aspect("equal")
        axes[0, 1].legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        if context:
            await context.info(
                "Spatial coordinates not available. "
                "This is expected for non-spatial datasets."
            )
        axes[0, 1].text(
            0.5, 0.5, "Spatial coordinates not available", ha="center", va="center"
        )
        axes[0, 1].set_title("Spatial (Not Available)", fontsize=12)

    # Panel 3: Batch composition bar plot
    batch_counts = adata.obs[batch_key].value_counts()
    axes[1, 0].bar(
        range(len(batch_counts)),
        batch_counts.values,
        color=colors[: len(batch_counts)],
    )
    axes[1, 0].set_xticks(range(len(batch_counts)))
    axes[1, 0].set_xticklabels(batch_counts.index, rotation=45, ha="right")
    axes[1, 0].set_title("Cell counts per batch", fontsize=12)
    axes[1, 0].set_ylabel("Number of cells")

    # Panel 4: Integration quality metrics (if available)
    axes[1, 1].text(
        0.1,
        0.9,
        "Integration Quality Assessment:",
        fontsize=14,
        fontweight="bold",
        transform=axes[1, 1].transAxes,
    )

    metrics_text = f"Total cells: {adata.n_obs:,}\n"
    metrics_text += f"Total genes: {adata.n_vars:,}\n"
    metrics_text += (
        f"Batches: {len(unique_batches)} ({', '.join(map(str, unique_batches))})\n\n"
    )

    if params.integration_method:
        metrics_text += f"Integration method: {params.integration_method}\n"

    # Add basic mixing metrics
    if "X_umap" in adata.obsm:
        # Calculate simple mixing metric (entropy)
        umap_coords = adata.obsm["X_umap"]
        x_bins = np.linspace(umap_coords[:, 0].min(), umap_coords[:, 0].max(), 10)
        y_bins = np.linspace(umap_coords[:, 1].min(), umap_coords[:, 1].max(), 10)

        entropies = []
        for i in range(len(x_bins) - 1):
            for j in range(len(y_bins) - 1):
                mask = (
                    (umap_coords[:, 0] >= x_bins[i])
                    & (umap_coords[:, 0] < x_bins[i + 1])
                    & (umap_coords[:, 1] >= y_bins[j])
                    & (umap_coords[:, 1] < y_bins[j + 1])
                )
                if mask.sum() > 10:  # Only consider regions with enough cells
                    batch_props = adata.obs[batch_key][mask].value_counts(
                        normalize=True
                    )
                    entropies.append(entropy(batch_props))

        if entropies:
            avg_entropy = np.mean(entropies)
            max_entropy = np.log(len(unique_batches))  # Perfect mixing entropy
            mixing_score = avg_entropy / max_entropy if max_entropy > 0 else 0
            metrics_text += (
                f"Mixing score: {mixing_score:.3f} (0=segregated, 1=perfectly mixed)\n"
            )

    axes[1, 1].text(
        0.1,
        0.7,
        metrics_text,
        fontsize=10,
        transform=axes[1, 1].transAxes,
        verticalalignment="top",
        fontfamily="monospace",
    )
    axes[1, 1].set_xlim(0, 1)
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].set_xticks([])
    axes[1, 1].set_yticks([])
    axes[1, 1].set_title("Integration Metrics", fontsize=12)

    plt.tight_layout()
    return fig
