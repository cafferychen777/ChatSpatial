"""
Core visualization utilities and shared functions.

This module contains:
- Figure setup and utility functions
- Shared data structures
- Common visualization helpers
"""

from typing import TYPE_CHECKING, Any, NamedTuple, Optional

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ...models.data import VisualizationParameters
from ...utils.adata_utils import get_gene_expression, require_spatial_coords
from ...utils.exceptions import DataNotFoundError, ParameterError

plt.ioff()

if TYPE_CHECKING:
    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# Figure Creation Utilities
# =============================================================================


# Default figure sizes by plot type for consistency
FIGURE_DEFAULTS = {
    "spatial": (10, 8),
    "umap": (10, 8),
    "heatmap": (12, 10),
    "violin": (12, 6),
    "dotplot": (10, 8),
    "trajectory": (10, 10),
    "gene_trends": (12, 6),
    "velocity": (10, 8),
    "deconvolution": (10, 8),
    "cell_communication": (10, 10),
    "enrichment": (6, 8),
    "cnv": (12, 8),
    "integration": (16, 12),
    "default": (10, 8),
}


def resolve_figure_size(
    params: VisualizationParameters,
    plot_type: str = "default",
    n_panels: Optional[int] = None,
    panel_width: float = 5.0,
    panel_height: float = 4.0,
) -> tuple[int, int]:
    """Resolve figure size from params with smart defaults.

    This centralizes figure size resolution logic to ensure consistency
    across all visualization modules.

    Args:
        params: VisualizationParameters with optional figure_size
        plot_type: Type of plot for default selection (e.g., "spatial", "heatmap")
        n_panels: Number of panels for multi-panel figures
        panel_width: Width per panel for multi-panel figures
        panel_height: Height per panel for multi-panel figures

    Returns:
        Tuple of (width, height) in inches

    Examples:
        >>> resolve_figure_size(params, "spatial")  # User override or (10, 8)
        >>> resolve_figure_size(params, n_panels=4)  # Compute from panel count
    """
    # User-specified size always takes precedence
    if params.figure_size:
        return params.figure_size

    # Multi-panel figure: compute from panel dimensions
    if n_panels is not None and n_panels > 1:
        n_cols = min(3, n_panels)
        n_rows = (n_panels + n_cols - 1) // n_cols
        width = min(panel_width * n_cols, 15)
        height = min(panel_height * n_rows, 16)
        return (int(width), int(height))

    # Use plot-type specific default
    return FIGURE_DEFAULTS.get(plot_type, FIGURE_DEFAULTS["default"])


def create_figure(figsize: tuple[int, int] = (10, 8)) -> tuple[plt.Figure, plt.Axes]:
    """Create a matplotlib figure with the right size and style."""
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


def create_figure_from_params(
    params: VisualizationParameters,
    plot_type: str = "default",
    n_panels: Optional[int] = None,
    n_rows: int = 1,
    n_cols: int = 1,
    squeeze: bool = True,
) -> tuple[plt.Figure, np.ndarray]:
    """Create a figure with axes from visualization parameters.

    This is the preferred way to create figures in visualization modules.
    It centralizes figure size resolution and applies consistent settings.

    Args:
        params: VisualizationParameters
        plot_type: Type of plot for default size selection
        n_panels: Number of panels (for auto-layout calculation)
        n_rows: Number of subplot rows
        n_cols: Number of subplot columns
        squeeze: Whether to squeeze single-element arrays

    Returns:
        Tuple of (Figure, array of Axes)

    Examples:
        >>> fig, axes = create_figure_from_params(params, "spatial")
        >>> fig, axes = create_figure_from_params(params, n_rows=2, n_cols=3)
    """
    figsize = resolve_figure_size(params, plot_type, n_panels)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=figsize,
        dpi=params.dpi,
        squeeze=squeeze,
    )

    # Ensure axes is always an array for consistent handling
    if squeeze and n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    elif squeeze and (n_rows == 1 or n_cols == 1):
        axes = np.atleast_1d(axes)

    return fig, axes


def setup_multi_panel_figure(
    n_panels: int,
    params: VisualizationParameters,
    default_title: str,
    use_tight_layout: bool = False,
) -> tuple[plt.Figure, np.ndarray]:
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

    # Only set suptitle if title is explicitly provided and non-empty
    # y=1.02 places title above figure to avoid overlap with subplot titles
    title = params.title or default_title
    if title:
        fig.suptitle(title, fontsize=16, y=1.02)

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
    cax = divider.append_axes(
        "right", size=params.colorbar_size, pad=params.colorbar_pad
    )
    cbar = fig.colorbar(mappable, cax=cax)
    if label:
        cbar.set_label(label, fontsize=10)


# =============================================================================
# Data Structures for Unified Data Access
# =============================================================================


class DeconvolutionData(NamedTuple):
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
    cell_types: list[str]
    proportions_key: str
    dominant_type_key: Optional[str] = None


class CellCommunicationData(NamedTuple):
    """Unified representation of cell communication analysis results.

    All visualization data comes from CCCStorage (single source of truth).
    Visualization functions should only read from this object, never from adata.uns directly.

    Attributes:
        results: Main results DataFrame (format varies by method)
        method: Analysis method ("liana", "cellphonedb", "fastccc", "cellchat_r")
        analysis_type: Type of analysis ("cluster" or "spatial")
        lr_pairs: List of ligand-receptor pair names (standardized: LIGAND_RECEPTOR)
        pvalues: P-values DataFrame/array (method-specific format)
        spatial_scores: Spatial communication scores array (n_spots x n_pairs)
        spatial_pvals: P-values for spatial scores (LIANA spatial only)
        source_labels: List of source cell type labels
        target_labels: List of target cell type labels
        method_data: Method-specific additional data (e.g., deconvoluted for CellPhoneDB)
    """

    results: pd.DataFrame
    method: str
    analysis_type: str  # "cluster" or "spatial"
    lr_pairs: list[str]
    pvalues: Optional[pd.DataFrame] = None
    spatial_scores: Optional[np.ndarray] = None
    spatial_pvals: Optional[np.ndarray] = None
    source_labels: Optional[list[str]] = None
    target_labels: Optional[list[str]] = None
    method_data: Optional[dict[str, Any]] = None


# =============================================================================
# Feature Validation and Preparation
# =============================================================================


async def get_validated_features(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
    max_features: Optional[int] = None,
    genes_only: bool = False,
) -> list[str]:
    """Validate and return features for visualization.

    Args:
        adata: AnnData object
        params: Visualization parameters containing feature specification
        context: Optional tool context for logging
        max_features: Maximum number of features to return (truncates if exceeded)
        genes_only: If True, only validate against var_names (genes).
                   If False, also check obs columns and obsm keys.

    Returns:
        List of validated feature names
    """
    if params.feature is None:
        features: list[str] = []
    elif isinstance(params.feature, list):
        features = params.feature
    else:
        features = [params.feature]
    validated: list[str] = []

    for feat in features:
        # Check if feature is in var_names (genes)
        if feat in adata.var_names:
            validated.append(feat)
        elif not genes_only:
            # Also check obs columns and obsm keys
            if feat in adata.obs.columns:
                validated.append(feat)
            elif feat in adata.obsm:
                validated.append(feat)
            else:
                if context:
                    await context.warning(
                        f"Feature '{feat}' not found in genes, obs, or obsm"
                    )
        else:
            if context:
                await context.warning(f"Gene '{feat}' not found in var_names")

    # Truncate if max_features specified
    if max_features is not None and len(validated) > max_features:
        if context:
            await context.warning(
                f"Too many features ({len(validated)}), limiting to {max_features}"
            )
        validated = validated[:max_features]

    return validated


def validate_and_prepare_feature(
    adata: ad.AnnData,
    feature: str,
    context: Optional["ToolContext"] = None,
) -> tuple[np.ndarray, str, bool]:
    """Validate a single feature and prepare its data for visualization.

    Args:
        adata: AnnData object
        feature: Feature name to validate
        context: Optional tool context for logging

    Returns:
        Tuple of (data array, display name, is_categorical)
    """
    # Gene expression - use unified utility
    if feature in adata.var_names:
        data = get_gene_expression(adata, feature)
        return data, feature, False

    # Observation column
    if feature in adata.obs.columns:
        data = adata.obs[feature]
        is_cat = pd.api.types.is_categorical_dtype(data) or data.dtype == object
        return data.values, feature, is_cat

    raise DataNotFoundError(f"Feature '{feature}' not found in data")


# =============================================================================
# Colormap Utilities
# =============================================================================

# Categorical colormaps by size threshold
_CATEGORICAL_CMAPS = {
    10: "tab10",  # Best for <= 10 categories
    20: "tab20",  # Best for 11-20 categories
    40: "tab20b",  # Extended palette for more categories
}


def get_categorical_cmap(n_categories: int, user_cmap: Optional[str] = None) -> str:
    """Select the best categorical colormap based on number of categories.

    This centralizes the categorical colormap selection logic that was
    previously scattered across visualization modules.

    Args:
        n_categories: Number of distinct categories to color
        user_cmap: User-specified colormap (takes precedence if provided
                  and is a known categorical palette)

    Returns:
        Colormap name suitable for categorical data

    Examples:
        >>> get_categorical_cmap(5)   # Returns "tab10"
        >>> get_categorical_cmap(15)  # Returns "tab20"
        >>> get_categorical_cmap(8, user_cmap="Set2")  # Returns "Set2"
    """
    # Known categorical palettes that user might specify
    categorical_palettes = {
        "tab10",
        "tab20",
        "tab20b",
        "tab20c",
        "Set1",
        "Set2",
        "Set3",
        "Paired",
        "Accent",
        "Dark2",
        "Pastel1",
        "Pastel2",
    }

    # User preference takes precedence if it's a categorical palette
    if user_cmap and user_cmap in categorical_palettes:
        return user_cmap

    # Auto-select based on category count
    for threshold, cmap in sorted(_CATEGORICAL_CMAPS.items()):
        if n_categories <= threshold:
            return cmap

    # Fallback for very large category counts
    return "tab20"


def get_category_colors(
    n_categories: int,
    cmap_name: Optional[str] = None,
) -> list:
    """Get a list of colors for categorical data.

    This is the primary function for obtaining colors for categorical
    visualizations. It handles colormap selection and color extraction.

    Args:
        n_categories: Number of categories to color
        cmap_name: Colormap name (auto-selected if None)

    Returns:
        List of colors (can be used with matplotlib scatter, legend, etc.)

    Examples:
        >>> colors = get_category_colors(5)  # 5 distinct colors
        >>> colors = get_category_colors(15, "tab20")  # 15 colors from tab20
    """
    # Select appropriate colormap
    if cmap_name is None:
        cmap_name = get_categorical_cmap(n_categories)

    # Seaborn palettes
    if cmap_name in ["tab10", "tab20", "Set1", "Set2", "Set3", "Paired", "husl"]:
        return sns.color_palette(cmap_name, n_colors=n_categories)

    # Matplotlib colormaps
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i / max(n_categories - 1, 1)) for i in range(n_categories)]


def get_colormap(name: str, n_colors: Optional[int] = None):
    """Get a matplotlib colormap by name.

    For categorical data, prefer using get_category_colors() instead.
    This function is for backward compatibility and continuous colormaps.

    Args:
        name: Colormap name (supports matplotlib and seaborn palettes)
        n_colors: Number of discrete colors (for categorical data)

    Returns:
        If n_colors is specified: List of colors (always indexable)
        Otherwise: Colormap object (for continuous data)
    """
    # For categorical with n_colors, delegate to specialized function
    if n_colors:
        return get_category_colors(n_colors, name)

    # Check if it's a seaborn palette (return as palette for consistency)
    if name in ["tab10", "tab20", "Set1", "Set2", "Set3", "Paired", "husl"]:
        return sns.color_palette(name)

    # For matplotlib colormaps, return the colormap object
    return plt.get_cmap(name)


def get_diverging_colormap() -> str:
    """Get an appropriate diverging colormap for symmetric data."""
    return "RdBu_r"


# =============================================================================
# Spatial Plot Utilities
# =============================================================================


def auto_spot_size(
    adata: ad.AnnData,
    user_spot_size: Optional[float] = None,
    basis: str = "spatial",
) -> float:
    """Calculate optimal spot size for visualization.

    Follows scanpy/squidpy conventions with a priority-based approach:
    1. User-specified value takes highest priority
    2. For spatial data with metadata: use spot_diameter_fullres * scale_factor
    3. Fallback: adaptive formula 120000 / n_cells (clamped to 5-200)

    Args:
        adata: AnnData object
        user_spot_size: User-specified spot size (takes priority if provided)
        basis: Embedding basis ("spatial", "umap", etc.)

    Returns:
        Calculated spot size for matplotlib scatter's s parameter

    Examples:
        >>> auto_spot_size(adata)  # Auto-calculate
        44.3
        >>> auto_spot_size(adata, user_spot_size=100)  # User override
        100.0
        >>> auto_spot_size(adata, basis="umap")  # UMAP uses formula only
        41.0
    """
    # Priority 1: User-specified value
    if user_spot_size is not None:
        return user_spot_size

    # Priority 2: For spatial basis, try to get from metadata
    if basis == "spatial" and "spatial" in adata.uns:
        spatial_data = adata.uns["spatial"]
        if spatial_data and isinstance(spatial_data, dict):
            # Get first library_id
            library_ids = list(spatial_data.keys())
            if library_ids:
                lib_data = spatial_data[library_ids[0]]

                if isinstance(lib_data, dict) and "scalefactors" in lib_data:
                    scalefactors = lib_data["scalefactors"]

                    # Get spot diameter
                    spot_diameter = scalefactors.get("spot_diameter_fullres")
                    if spot_diameter and spot_diameter > 0:
                        # Get scale factor (prefer hires, fallback to lowres)
                        scale_factor = scalefactors.get(
                            "tissue_hires_scalef",
                            scalefactors.get("tissue_lowres_scalef", 1.0),
                        )

                        # Calculate: scatter s parameter is area (diameter^2 based)
                        # Apply 0.5 adjustment factor to match typical visual expectations
                        calculated_size = (spot_diameter * scale_factor * 0.5) ** 2
                        return max(calculated_size, 5.0)  # minimum size of 5

    # Priority 3: Adaptive formula based on cell count
    n_cells = adata.n_obs
    adaptive_size = 120000 / n_cells

    # Clamp to reasonable range [5, 200]
    return max(min(adaptive_size, 200.0), 5.0)


def plot_spatial_feature(
    adata: ad.AnnData,
    ax: plt.Axes,
    feature: Optional[str] = None,
    values: Optional[np.ndarray] = None,
    params: Optional[VisualizationParameters] = None,
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

    # Calculate spot size (auto or user-specified)
    spot_size = auto_spot_size(adata, params.spot_size, basis=spatial_key)

    # Get spatial coordinates
    coords = require_spatial_coords(adata, spatial_key=spatial_key)

    # Get values to plot
    if values is not None:
        plot_values = values
        is_categorical = pd.api.types.is_categorical_dtype(values)
    elif feature is not None:
        if feature in adata.var_names:
            # Use unified utility for gene expression extraction
            plot_values = get_gene_expression(adata, feature)
            is_categorical = False
        elif feature in adata.obs.columns:
            plot_values = adata.obs[feature].values
            is_categorical = pd.api.types.is_categorical_dtype(adata.obs[feature])
        else:
            raise DataNotFoundError(f"Feature '{feature}' not found")
    else:
        raise ParameterError("Either feature or values must be provided")

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
            s=spot_size,
            alpha=params.alpha,
        )

        # Add legend for categorical
        if params.show_legend:
            handles = [
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=colors[i],
                    markersize=8,
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
            s=spot_size,
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


# =============================================================================
# Data Inference Utilities
# =============================================================================


def get_categorical_columns(
    adata: ad.AnnData,
    limit: Optional[int] = None,
) -> list[str]:
    """Get categorical column names from adata.obs.

    Args:
        adata: AnnData object
        limit: Maximum number of columns to return (None for all)

    Returns:
        List of categorical column names
    """
    categorical_cols = [
        col
        for col in adata.obs.columns
        if adata.obs[col].dtype.name in ["object", "category"]
    ]
    if limit is not None:
        return categorical_cols[:limit]
    return categorical_cols


def infer_basis(
    adata: ad.AnnData,
    preferred: Optional[str] = None,
    priority: Optional[list[str]] = None,
) -> Optional[str]:
    """Infer the best embedding basis from available options.

    Args:
        adata: AnnData object
        preferred: User-specified preferred basis (returned if valid)
        priority: Priority order for basis selection.
                  Default: ["spatial", "umap", "pca"]

    Returns:
        Best available basis name (without X_ prefix), or None if none found

    Examples:
        >>> infer_basis(adata)  # Auto-detect: spatial > umap > pca
        'umap'
        >>> infer_basis(adata, preferred='tsne')  # Use if valid
        'tsne'
        >>> infer_basis(adata, priority=['umap', 'spatial'])  # Custom order
        'umap'
    """
    if priority is None:
        priority = ["spatial", "umap", "pca"]

    # Check preferred basis first
    if preferred:
        key = preferred if preferred == "spatial" else f"X_{preferred}"
        if key in adata.obsm:
            return preferred

    # Check priority list
    for basis in priority:
        key = basis if basis == "spatial" else f"X_{basis}"
        if key in adata.obsm:
            return basis

    # Fallback: return first available X_* key
    for key in adata.obsm.keys():
        if key.startswith("X_"):
            return key[2:]  # Strip X_ prefix

    return None
