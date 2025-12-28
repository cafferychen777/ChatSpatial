"""
Core visualization utilities and shared functions.

This module contains:
- Figure setup and utility functions
- Shared data classes
- Common visualization helpers
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING, List, Optional, Tuple

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ...models.data import VisualizationParameters
from ...utils.adata_utils import get_spatial_coordinates

plt.ioff()

if TYPE_CHECKING:
    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# Figure Creation Utilities
# =============================================================================


def create_figure(figsize: Tuple[int, int] = (10, 8)) -> Tuple[plt.Figure, plt.Axes]:
    """Create a matplotlib figure with the right size and style."""
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


def setup_multi_panel_figure(
    n_panels: int,
    params: VisualizationParameters,
    default_title: str,
    use_tight_layout: bool = False,
) -> Tuple[plt.Figure, np.ndarray]:
    """Sets up a multi-panel matplotlib figure.

    Args:
        n_panels: The total number of panels required.
        params: VisualizationParameters object with GridSpec spacing parameters.
        default_title: Default title for the figure if not provided in params.
        use_tight_layout: If True, skip gridspec_kw and use tight_layout.

    Returns:
        A tuple of (matplotlib.Figure, flattened numpy.ndarray of Axes).
    """
    if params.panel_layout:
        n_rows, n_cols = params.panel_layout
    else:
        n_cols = min(3, n_panels)
        n_rows = (n_panels + n_cols - 1) // n_cols

    if params.figure_size:
        figsize = params.figure_size
    else:
        figsize = (min(5 * n_cols, 15), min(4 * n_rows, 16))

    if not use_tight_layout:
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=figsize,
            dpi=params.dpi,
            squeeze=False,
            gridspec_kw={
                "wspace": params.subplot_wspace,
                "hspace": params.subplot_hspace,
            },
        )
    else:
        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=figsize, dpi=params.dpi, squeeze=False
        )

    axes = axes.flatten()

    title = params.title or default_title
    fig.suptitle(title, fontsize=16)

    for i in range(n_panels, len(axes)):
        axes[i].axis("off")

    return fig, axes


def add_colorbar(
    fig: plt.Figure,
    ax: plt.Axes,
    mappable,
    params: VisualizationParameters,
    label: str = "",
) -> None:
    """Add a colorbar to an axis with consistent styling.

    Args:
        fig: The figure object
        ax: The axes object to attach colorbar to
        mappable: The mappable object (from scatter, imshow, etc.)
        params: Visualization parameters for styling
        label: Colorbar label
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=params.colorbar_size, pad=params.colorbar_pad)
    cbar = fig.colorbar(mappable, cax=cax)
    if label:
        cbar.set_label(label, fontsize=10)


# =============================================================================
# Data Classes for Unified Data Access
# =============================================================================


@dataclass
class DeconvolutionData:
    """Unified representation of deconvolution results.

    Attributes:
        proportions: DataFrame with cell type proportions (n_spots x n_cell_types)
        method: Deconvolution method name (e.g., "cell2location", "rctd")
        cell_types: List of cell type names
        proportions_key: Key in adata.obsm where proportions are stored
        dominant_type_key: Key in adata.obs for dominant cell type (if exists)
    """

    proportions: pd.DataFrame
    method: str
    cell_types: List[str]
    proportions_key: str
    dominant_type_key: Optional[str] = None

    @property
    def n_cell_types(self) -> int:
        return len(self.cell_types)

    @property
    def n_spots(self) -> int:
        return len(self.proportions)


@dataclass
class CellCommunicationData:
    """Unified representation of cell communication analysis results.

    Attributes:
        results: Main results DataFrame (format varies by method)
        method: Analysis method name ("liana_cluster", "liana_spatial", "cellphonedb")
        analysis_type: Type of analysis ("cluster" or "spatial")
        lr_pairs: List of ligand-receptor pair names
        spatial_scores: Spatial communication scores array (n_spots x n_pairs)
        spatial_pvals: P-values for spatial scores (optional)
        source_labels: List of source cell type labels
        target_labels: List of target cell type labels
        results_key: Key in adata.uns where results are stored
    """

    results: pd.DataFrame
    method: str
    analysis_type: str  # "cluster" or "spatial"
    lr_pairs: List[str]
    spatial_scores: Optional[np.ndarray] = None
    spatial_pvals: Optional[np.ndarray] = None
    source_labels: Optional[List[str]] = None
    target_labels: Optional[List[str]] = None
    results_key: str = ""

    @property
    def n_interactions(self) -> int:
        return len(self.results)


# =============================================================================
# Feature Validation and Preparation
# =============================================================================


async def get_validated_features(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> List[str]:
    """Validate and return features for visualization.

    Args:
        adata: AnnData object
        params: Visualization parameters containing feature specification
        context: Optional tool context for logging

    Returns:
        List of validated feature names
    """
    features = params.feature if isinstance(params.feature, list) else [params.feature]
    validated = []

    for feat in features:
        if feat is None:
            continue

        # Check if feature is in var_names (genes)
        if feat in adata.var_names:
            validated.append(feat)
        # Check if feature is in obs columns
        elif feat in adata.obs.columns:
            validated.append(feat)
        # Check if feature is an obsm key
        elif feat in adata.obsm:
            validated.append(feat)
        else:
            if context:
                await context.warning(
                    f"Feature '{feat}' not found in genes, obs, or obsm"
                )

    return validated


async def validate_and_prepare_feature(
    adata: ad.AnnData,
    feature: str,
    context: Optional["ToolContext"] = None,
) -> Tuple[np.ndarray, str, bool]:
    """Validate a single feature and prepare its data for visualization.

    Args:
        adata: AnnData object
        feature: Feature name to validate
        context: Optional tool context for logging

    Returns:
        Tuple of (data array, display name, is_categorical)
    """
    # Gene expression
    if feature in adata.var_names:
        gene_idx = adata.var_names.get_loc(feature)
        if hasattr(adata.X, "toarray"):
            data = adata.X[:, gene_idx].toarray().flatten()
        else:
            data = adata.X[:, gene_idx].flatten()
        return data, feature, False

    # Observation column
    if feature in adata.obs.columns:
        data = adata.obs[feature]
        is_cat = pd.api.types.is_categorical_dtype(data) or data.dtype == object
        return data.values, feature, is_cat

    raise ValueError(f"Feature '{feature}' not found in data")


# =============================================================================
# Colormap Utilities
# =============================================================================


def get_colormap(name: str, n_colors: Optional[int] = None):
    """Get a matplotlib colormap by name.

    Args:
        name: Colormap name (supports matplotlib and seaborn palettes)
        n_colors: Number of discrete colors (for categorical data)

    Returns:
        Colormap object or list of colors
    """
    # Check if it's a seaborn palette
    if name in ["tab10", "tab20", "Set1", "Set2", "Set3", "Paired", "husl"]:
        if n_colors:
            return sns.color_palette(name, n_colors=n_colors)
        return sns.color_palette(name)

    # Otherwise use matplotlib
    return plt.get_cmap(name)


def get_diverging_colormap(center: float = 0.0) -> str:
    """Get an appropriate diverging colormap centered at a value."""
    return "RdBu_r"


# =============================================================================
# Spatial Plot Utilities
# =============================================================================


def plot_spatial_feature(
    adata: ad.AnnData,
    ax: plt.Axes,
    feature: Optional[str] = None,
    values: Optional[np.ndarray] = None,
    params: VisualizationParameters = None,
    spatial_key: str = "spatial",
    show_colorbar: bool = True,
    title: Optional[str] = None,
) -> Optional[plt.cm.ScalarMappable]:
    """Plot a feature on spatial coordinates.

    Args:
        adata: AnnData object with spatial coordinates
        ax: Matplotlib axes to plot on
        feature: Feature name (gene or obs column)
        values: Pre-computed values to plot (overrides feature)
        params: Visualization parameters
        spatial_key: Key for spatial coordinates in obsm
        show_colorbar: Whether to add a colorbar
        title: Plot title

    Returns:
        ScalarMappable for colorbar creation, or None for categorical data
    """
    if params is None:
        params = VisualizationParameters()

    # Get spatial coordinates
    coords = get_spatial_coordinates(adata, spatial_key)

    # Get values to plot
    if values is not None:
        plot_values = values
        is_categorical = pd.api.types.is_categorical_dtype(values)
    elif feature is not None:
        if feature in adata.var_names:
            gene_idx = adata.var_names.get_loc(feature)
            if hasattr(adata.X, "toarray"):
                plot_values = adata.X[:, gene_idx].toarray().flatten()
            else:
                plot_values = adata.X[:, gene_idx].flatten()
            is_categorical = False
        elif feature in adata.obs.columns:
            plot_values = adata.obs[feature].values
            is_categorical = pd.api.types.is_categorical_dtype(adata.obs[feature])
        else:
            raise ValueError(f"Feature '{feature}' not found")
    else:
        raise ValueError("Either feature or values must be provided")

    # Handle categorical data
    if is_categorical:
        categories = (
            plot_values.categories
            if hasattr(plot_values, "categories")
            else np.unique(plot_values)
        )
        n_cats = len(categories)
        colors = get_colormap(params.colormap, n_colors=n_cats)
        cat_to_idx = {cat: i for i, cat in enumerate(categories)}
        color_indices = [cat_to_idx[v] for v in plot_values]

        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=[colors[i] for i in color_indices],
            s=params.spot_size,
            alpha=params.alpha,
        )

        # Add legend for categorical
        if params.show_legend:
            handles = [
                plt.Line2D(
                    [0], [0], marker="o", color="w", markerfacecolor=colors[i], markersize=8
                )
                for i in range(n_cats)
            ]
            ax.legend(
                handles,
                categories,
                loc="center left",
                bbox_to_anchor=(1, 0.5),
                fontsize=8,
            )
        mappable = None
    else:
        # Continuous data
        cmap = get_colormap(params.colormap)
        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=plot_values,
            cmap=cmap,
            s=params.spot_size,
            alpha=params.alpha,
            vmin=params.vmin,
            vmax=params.vmax,
        )
        mappable = scatter

    ax.set_aspect("equal")
    ax.set_xlabel("")
    ax.set_ylabel("")

    if not params.show_axes:
        ax.axis("off")

    if title:
        ax.set_title(title, fontsize=12)

    return mappable
