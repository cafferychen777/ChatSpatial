"""
Plotting utilities for spatial transcriptomics data.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Optional, Tuple, List, Dict, Any, Union
from mcp.server.fastmcp.utilities.types import Image

# Import standardized image utilities
from .image_utils import fig_to_image, fig_to_base64, create_placeholder_image


def create_spatial_plot(
    coordinates: np.ndarray,
    values: Optional[np.ndarray] = None,
    tissue_image: Optional[np.ndarray] = None,
    spot_size: float = 50.0,
    colormap: str = 'viridis',
    title: str = 'Spatial Plot',
    feature_name: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 8)
) -> plt.Figure:
    """Create a spatial plot for transcriptomics data

    Args:
        coordinates: Array of (x, y) coordinates with shape (n_spots, 2)
        values: Optional array of values to color the spots by
        tissue_image: Optional tissue image as background
        spot_size: Size of spots in the plot
        colormap: Colormap for the values
        title: Plot title
        feature_name: Name of the feature being plotted (for colorbar label)
        figsize: Figure size in inches

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    # If tissue image is available, show it as background
    if tissue_image is not None:
        ax.imshow(tissue_image, alpha=0.7)

    # Plot spots
    if values is not None:
        scatter = ax.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            c=values,
            s=spot_size,
            cmap=colormap,
            alpha=0.8
        )
        cbar = plt.colorbar(scatter, ax=ax)
        if feature_name:
            cbar.set_label(feature_name)
    else:
        ax.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            s=spot_size,
            alpha=0.8
        )

    ax.set_title(title)
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_aspect('equal')

    return fig


def create_cluster_plot(
    embedding: np.ndarray,
    clusters: np.ndarray,
    cluster_names: Optional[List[str]] = None,
    title: str = 'Cluster Plot',
    figsize: Tuple[int, int] = (10, 8)
) -> plt.Figure:
    """Create a cluster plot (e.g., UMAP, t-SNE) colored by clusters

    Args:
        embedding: Array of 2D embedding coordinates with shape (n_cells, 2)
        clusters: Array of cluster assignments
        cluster_names: Optional list of cluster names for legend
        title: Plot title
        figsize: Figure size in inches

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Convert clusters to numeric if they're not already
    if not np.issubdtype(clusters.dtype, np.number):
        unique_clusters = np.unique(clusters)
        cluster_map = {c: i for i, c in enumerate(unique_clusters)}
        numeric_clusters = np.array([cluster_map[c] for c in clusters])
    else:
        numeric_clusters = clusters

    # Plot clusters
    scatter = ax.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=numeric_clusters,
        cmap='tab20',
        alpha=0.7,
        s=10
    )

    # Add legend if cluster names are provided
    if cluster_names is not None:
        handles, labels = scatter.legend_elements()
        ax.legend(handles, cluster_names, title="Clusters", loc="upper right")

    ax.set_title(title)
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    return fig


def create_heatmap(
    data: np.ndarray,
    row_labels: Optional[List[str]] = None,
    col_labels: Optional[List[str]] = None,
    title: str = 'Heatmap',
    colormap: str = 'viridis',
    figsize: Tuple[int, int] = (12, 10)
) -> plt.Figure:
    """Create a heatmap for gene expression or other matrix data

    Args:
        data: 2D array of values
        row_labels: Optional list of row labels
        col_labels: Optional list of column labels
        title: Plot title
        colormap: Colormap for the heatmap
        figsize: Figure size in inches

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(data, cmap=colormap, aspect='auto')
    plt.colorbar(im, ax=ax)

    # Add row and column labels if provided
    if row_labels is not None:
        ax.set_yticks(np.arange(len(row_labels)))
        ax.set_yticklabels(row_labels)

    if col_labels is not None:
        ax.set_xticks(np.arange(len(col_labels)))
        ax.set_xticklabels(col_labels, rotation=45, ha='right')

    ax.set_title(title)

    # Adjust layout to make room for labels
    plt.tight_layout()

    return fig