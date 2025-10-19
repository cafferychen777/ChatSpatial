"""
Visualization tools for spatial transcriptomics data.
"""

import os
import traceback
from datetime import datetime
from functools import wraps
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import anndata as ad
# Set non-interactive backend for matplotlib to prevent GUI popups on macOS
import matplotlib

matplotlib.use("Agg")  # Non-interactive backend - prevents Dock popup
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.ioff()  # Turn off interactive mode
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import seaborn as sns  # noqa: E402
from mcp.server.fastmcp import Context  # noqa: E402
from mcp.types import EmbeddedResource, ImageContent  # noqa: E402
from scipy.stats import pearsonr, spearmanr  # noqa: E402

# Optional CNV analysis visualization
try:
    import infercnvpy as cnv

    INFERCNVPY_AVAILABLE = True
except ImportError:
    INFERCNVPY_AVAILABLE = False

from ..models.data import VisualizationParameters  # noqa: E402
# Import spatial coordinates helper from data adapter
from ..utils.data_adapter import get_spatial_coordinates  # noqa: E402
# Import error handling utilities
from ..utils.error_handling import (DataCompatibilityError, DataNotFoundError,  # noqa: E402
                                    InvalidParameterError, ProcessingError)
# Import standardized image utilities
from ..utils.image_utils import (create_placeholder_image,  # noqa: E402
                                 optimize_fig_to_image_with_cache)
# Import path utilities for safe file operations
from ..utils.path_utils import (get_output_dir_from_config,  # noqa: E402
                                get_safe_output_path, is_safe_output_path)
# Import color utilities for categorical data
from ._color_utils import _ensure_categorical_colors  # noqa: E402

# Import publication export utilities


# Helper function to create a figure with the right size
def create_figure(figsize=(10, 8)):
    """Create a matplotlib figure with the right size and style"""
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


# New helper functions to reduce code duplication
def setup_multi_panel_figure(
    n_panels: int,
    params: VisualizationParameters,
    default_title: str,
    use_tight_layout: bool = False,
) -> Tuple[plt.Figure, np.ndarray]:
    """Sets up a multi-panel matplotlib figure.

    Args:
        n_panels: The total number of panels required.
        params: VisualizationParameters object.
        default_title: Default title for the figure if not provided in params.
        use_tight_layout: If True, skip subplots_adjust and use tight_layout later (for plots with colorbars).

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
        # Adjust size to prevent images from becoming too large
        figsize = (min(5 * n_cols, 15), min(4 * n_rows, 16))

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=figsize, dpi=params.dpi, squeeze=False
    )
    axes = axes.flatten()

    title = params.title or default_title
    fig.suptitle(title, fontsize=16)

    # Only use subplots_adjust if not using tight_layout
    # tight_layout is better for plots with colorbars (L-R pairs, multi-gene spatial)
    if not use_tight_layout:
        # Adjust subplot to make room for suptitle
        # Reduced spacing for tighter layout: hspace=0.3 (was 0.4), wspace=0.2 (was 0.3)
        plt.subplots_adjust(top=0.92, hspace=0.3, wspace=0.2)

    # Hide empty axes from the start
    for i in range(n_panels, len(axes)):
        axes[i].axis("off")

    return fig, axes


async def get_validated_features(
    adata: ad.AnnData,
    params: VisualizationParameters,
    max_features: int = 12,
    context: Optional[Context] = None,
) -> List[str]:
    """Gets a validated list of features (genes) for visualization."""

    # Ensure unique var_names before proceeding
    if not adata.var_names.is_unique:
        if context:
            await context.info("Making gene names unique to avoid indexing errors.")
        adata.var_names_make_unique()

    # Check if user specified any features
    if not params.feature:
        raise DataNotFoundError(
            "❌ No genes specified for visualization.\n\n"
            "🔧 SOLUTIONS:\n"
            "1. Specify genes: visualize_data(data_id, params={'feature': ['CD3D', 'CD8A']})\n"
            "2. Find markers: find_markers(data_id, group_key='cell_type')\n"
            "3. Run preprocessing: preprocess_data(data_id, params={'n_top_genes': 2000})\n\n"
            "📋 CONTEXT: Explicit gene selection required for scientific accuracy."
        )

    # Parse user features
    features = params.feature if isinstance(params.feature, list) else [params.feature]

    # Check which features are available
    available_features = [f for f in features if f in adata.var_names]

    if not available_features:
        missing = [f for f in features if f not in adata.var_names]
        examples = list(adata.var_names[:10])
        raise DataNotFoundError(
            f"❌ Genes not found: {missing}\n\n"
            f"🔧 SOLUTIONS:\n"
            f"1. Check names (available: {examples})\n"
            f"2. Search for similar gene names\n"
            f"3. Use gene discovery tools\n\n"
            f"📋 CONTEXT: Dataset has {adata.n_vars:,} genes."
        )

    if len(available_features) > max_features:
        raise ValueError(
            f"❌ Too many genes ({len(available_features)} > {max_features}).\n\n"
            f"🔧 SOLUTIONS:\n"
            f"1. Select {max_features} genes from: {available_features}\n"
            f"2. Create multiple visualizations\n"
            f"3. Use heatmap for larger sets"
        )

    return available_features


def plot_spatial_feature(
    adata: ad.AnnData,
    feature: Optional[str],
    ax: plt.Axes,
    params: VisualizationParameters,
    force_no_colorbar: bool = False,
):
    """Plots a feature on spatial coordinates, handling background images.

    Note: This function sets the title to an empty string. Callers should
    manually set ax.set_title() after calling this function if a specific
    title is desired.

    Args:
        adata: AnnData object
        feature: Feature to plot
        ax: Matplotlib axes
        params: Visualization parameters
        force_no_colorbar: If True, disable colorbar even if params.show_colorbar is True.
                          Used when caller wants to manually add colorbar with make_axes_locatable.
    """
    has_image = "spatial" in adata.uns and "images" in adata.uns["spatial"]

    # Base kwargs for both functions
    plot_kwargs = {
        "color": feature,
        "ax": ax,
        "show": False,
        "cmap": params.colormap,
        "alpha": params.alpha,
        "frameon": params.show_axes,
        "colorbar_loc": (
            None if force_no_colorbar else ("right" if params.show_colorbar else None)
        ),
        "title": "",  # We will set the title manually
    }

    # For spatial plots with tissue image, use 'spot_size'; for embedding use 'size'
    if has_image:
        # Add alpha_img for tissue background dimming (only supported by sc.pl.spatial)
        plot_kwargs["alpha_img"] = params.alpha_img
        # Only add spot_size if explicitly set (None = auto-determined, recommended)
        if params.spot_size is not None:
            plot_kwargs["spot_size"] = params.spot_size
        sc.pl.spatial(adata, img_key="hires", **plot_kwargs)
    else:
        # For embeddings, use 'size' parameter if spot_size is set
        if params.spot_size is not None:
            plot_kwargs["size"] = params.spot_size
        sc.pl.embedding(adata, basis="spatial", **plot_kwargs)
        ax.set_aspect("equal")

    # Note: We don't set the title here to allow callers to set their own titles
    # The old code: if params.add_gene_labels and feature: ax.set_title(feature, fontsize=12)
    # is removed so callers can manually set appropriate titles


def handle_visualization_errors(plot_title: str):
    """
    Improved decorator to catch errors in visualization functions.

    Strategy:
    - Scientific/data errors (DataNotFoundError, etc.): Propagate to MCP layer
    - Technical/rendering errors (matplotlib, memory): Return error figure
    - Unknown errors: Propagate (conservative approach)
    """

    def decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            try:
                # The wrapped function is async, so we need to await it
                return await func(*args, **kwargs)

            except (
                DataNotFoundError,
                InvalidParameterError,
                DataCompatibilityError,
            ):
                # Scientific/user errors must propagate for proper error handling
                # These indicate issues with data or parameters that the user must fix
                raise

            except (ValueError, MemoryError, RuntimeError, ImportError) as e:
                # Technical/rendering errors - return a friendly error figure
                # These are often transient or due to system limitations
                fig, ax = plt.subplots(figsize=(8, 6))

                # Create clear error message with visual emphasis
                error_text = (
                    "⚠️ TECHNICAL ERROR - NOT A DATA ISSUE ⚠️\n\n"
                    f"{type(e).__name__}: {str(e)}\n\n"
                    "This is a rendering/technical issue, not a problem with your data.\n"
                    "The analysis may have succeeded but visualization failed."
                )

                ax.text(
                    0.5,
                    0.5,
                    error_text,
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    wrap=True,
                    fontsize=10,
                    color="darkred",
                    weight="bold",
                )

                ax.set_title(
                    f"{plot_title} - Technical Error",
                    color="red",
                    fontsize=14,
                    weight="bold",
                )
                ax.axis("off")

                # Add red border to emphasize this is an error
                rect = plt.Rectangle(
                    (0.02, 0.02),
                    0.96,
                    0.96,
                    transform=ax.transAxes,
                    linewidth=3,
                    edgecolor="red",
                    facecolor="none",
                    linestyle="--",
                )
                ax.add_patch(rect)

                return fig

            except Exception:
                # Unknown exceptions are propagated (conservative strategy)
                # This ensures we don't hide unexpected errors
                raise

        return wrapper

    return decorator


# Helper function to validate and prepare feature for visualization
async def validate_and_prepare_feature(
    adata,
    feature: Optional[str],
    context: Optional[Context] = None,
    default_feature: Optional[str] = None,
) -> Optional[str]:
    """Validate and prepare a feature for visualization

    This function handles various feature formats including:
    - Regular gene names or observation annotations
    - Deconvolution results with cell type using colon format (deconvolution_key:cell_type)
    - Deconvolution results with cell type using underscore format (deconvolution_key_cell_type)
    - Deconvolution keys in obsm

    Args:
        adata: AnnData object
        feature: Feature name or None
        context: MCP context for logging
        default_feature: Default feature to use if the provided feature is not found

    Returns:
        Validated feature name or default_feature or None
    """
    if not feature:
        return default_feature

    # Check if it's a deconvolution result with cell type using colon format
    if ":" in feature:
        # Format: "deconvolution_method:cell_type"
        parts = feature.split(":")
        if len(parts) == 2:
            deconv_key, cell_type = parts

            # First check if we already have this in obs (from previous deconvolution)
            obs_key = f"{deconv_key}_{cell_type}"
            if obs_key in adata.obs.columns:
                if context:
                    await context.info(
                        f"Found cell type {cell_type} in obs as {obs_key}"
                    )
                return obs_key
            else:
                # Get deconvolution results as DataFrame
                deconv_df = get_deconvolution_dataframe(adata, deconv_key)
                if deconv_df is not None and cell_type in deconv_df.columns:
                    # Add the cell type proportion to obs for visualization
                    adata.obs[feature] = deconv_df[cell_type]
                    return feature
                else:
                    if context:
                        await context.warning(
                            f"Deconvolution result {feature} not found. Using default feature."
                        )
                    return default_feature
        else:
            if context:
                await context.warning(
                    f"Invalid feature format {feature}. Use 'deconvolution_key:cell_type'. Using default feature."
                )
            return default_feature

    # Check if it's a deconvolution result with cell type using underscore format
    elif "_" in feature and feature.startswith("deconvolution_"):
        # Format: "deconvolution_method_cell_type"
        if feature in adata.obs.columns:
            if context:
                await context.info(f"Found cell type proportion in obs as {feature}")
            return feature
        else:
            if context:
                await context.warning(
                    f"Feature {feature} not found in dataset. Using default feature."
                )
            return default_feature

    # Check if it's a deconvolution key
    elif feature in adata.obsm:
        # Try to get cell types from uns
        cell_types_key = f"{feature}_cell_types"
        if cell_types_key in adata.uns and len(adata.uns[cell_types_key]) > 0:
            # Get top cell type by mean proportion
            deconv_df = get_deconvolution_dataframe(adata, feature)
            if deconv_df is not None:
                top_cell_type = deconv_df.mean().sort_values(ascending=False).index[0]
                obs_key = f"{feature}_{top_cell_type}"

                # Add to obs if not already there
                if obs_key not in adata.obs.columns:
                    adata.obs[obs_key] = deconv_df[top_cell_type].values

                if context:
                    await context.info(
                        f"Found deconvolution result {feature}. Showing top cell type: {top_cell_type}"
                    )
                return obs_key
            else:
                if context:
                    await context.info(
                        f"Found deconvolution result {feature}, but could not determine cell types. Please specify a cell type using 'deconvolution_key:cell_type' format."
                    )
                return default_feature
        else:
            if context:
                await context.info(
                    f"Found deconvolution result {feature}. Please specify a cell type using 'deconvolution_key:cell_type' format."
                )
            return default_feature

    # Check if it's a regular feature
    elif feature not in adata.var_names and feature not in adata.obs.columns:
        if context:
            await context.warning(
                f"Feature {feature} not found in dataset. Using default feature."
            )
        return default_feature
    else:
        return feature


# Helper function to get deconvolution results as DataFrame
def get_deconvolution_dataframe(adata, deconv_key):
    """Get deconvolution results as DataFrame

    Args:
        adata: AnnData object
        deconv_key: Key in adata.obsm for deconvolution results

    Returns:
        DataFrame with deconvolution results or None if not found
    """
    if deconv_key not in adata.obsm:
        return None

    # Get deconvolution results
    deconv_results = adata.obsm[deconv_key]

    # Convert to DataFrame if it's a numpy array
    if isinstance(deconv_results, np.ndarray):
        # Try to get cell types from uns
        if f"{deconv_key}_cell_types" in adata.uns:
            cell_types = adata.uns[f"{deconv_key}_cell_types"]
            return pd.DataFrame(
                deconv_results, index=adata.obs_names, columns=cell_types
            )
        else:
            # Use generic column names
            return pd.DataFrame(
                deconv_results,
                index=adata.obs_names,
                columns=[f"Cell_Type_{i}" for i in range(deconv_results.shape[1])],
            )

    # If it's already a DataFrame, return it
    elif isinstance(deconv_results, pd.DataFrame):
        return deconv_results

    # Otherwise, try to convert it
    else:
        try:
            return pd.DataFrame(deconv_results, index=adata.obs_names)
        except Exception:
            return None


async def visualize_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: VisualizationParameters = VisualizationParameters(),
    context: Optional[Context] = None,
) -> Union[ImageContent, Tuple[ImageContent, EmbeddedResource]]:
    """Visualize spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Visualization parameters
        context: MCP context

    Returns:
        Union[ImageContent, Tuple[ImageContent, EmbeddedResource]]:
            - Small images (<100KB): ImageContent object
            - Large images (>=100KB): Tuple[Preview ImageContent, High-quality Resource]

    Raises:
        DataNotFoundError: If the dataset is not found
        InvalidParameterError: If parameters are invalid
        DataCompatibilityError: If data is not compatible with the visualization
        ProcessingError: If processing fails
    """
    # Validate parameters
    valid_plot_types = [
        "spatial",
        "umap",
        "heatmap",
        "violin",
        "deconvolution",
        "spatial_domains",
        "cell_communication",
        "trajectory",
        "rna_velocity",
        "spatial_statistics",
        "multi_gene",
        "lr_pairs",
        "gene_correlation",
        "pathway_enrichment",
        "spatial_enrichment",
        "spatial_interaction",
        "batch_integration",  # Batch integration quality assessment
        "cnv_heatmap",  # Copy number variation heatmap
        "spatial_cnv",  # Spatial CNV projection
        "card_imputation",  # CARD imputation high-resolution results
    ]
    if params.plot_type not in valid_plot_types:
        error_msg = (
            f"Invalid plot_type: {params.plot_type}. Must be one of {valid_plot_types}"
        )
        if context:
            await context.warning(error_msg)
        raise InvalidParameterError(error_msg)

    if context:
        await context.info(f"Visualizing {params.plot_type} plot for dataset {data_id}")
        await context.info(
            f"Parameters: feature={params.feature}, colormap={params.colormap}"
        )

    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        error_msg = f"Dataset {data_id} not found in data store"
        if context:
            await context.warning(error_msg)
        raise DataNotFoundError(error_msg)

    try:
        adata = data_store[data_id]["adata"]

        # Validate AnnData object - basic validation
        if adata.n_obs < 5:
            raise DataNotFoundError("Dataset has too few cells (minimum 5 required)")
        if adata.n_vars < 5:
            raise DataNotFoundError("Dataset has too few genes (minimum 5 required)")

        # Set matplotlib style for better visualizations
        sc.settings.set_figure_params(dpi=120, facecolor="white")

        # Create figure based on plot type
        if params.plot_type == "spatial":
            # Check if this should be a multi-gene visualization
            feature_list = (
                params.feature
                if isinstance(params.feature, list)
                else ([params.feature] if params.feature else [])
            )
            if feature_list and params.multi_panel and len(feature_list) > 1:
                if context:
                    await context.info(
                        f"Creating multi-gene spatial plot for {len(feature_list)} genes"
                    )
                # Create multi-gene visualization
                fig = await create_multi_gene_visualization(adata, params, context)
            else:
                single_feature = feature_list[0] if feature_list else None
                if context:
                    await context.info(
                        f"Creating spatial plot for {single_feature if single_feature else 'tissue'}"
                    )

                # Validate and prepare feature for visualization
                feature = await validate_and_prepare_feature(
                    adata, single_feature, context, default_feature=None
                )

                # If user explicitly requested a feature but it wasn't found, raise error
                if single_feature and not feature:
                    # Provide helpful error message
                    examples = list(adata.var_names[:10])
                    obs_examples = list(adata.obs.columns[:10])
                    raise DataNotFoundError(
                        f"❌ Feature '{single_feature}' not found in dataset.\n\n"
                        f"🔧 SOLUTIONS:\n"
                        f"1. Check spelling and case (gene names are case-sensitive)\n"
                        f"   - Mouse genes: First letter uppercase (e.g., 'Cd5l', 'Gbp2b')\n"
                        f"   - Human genes: All uppercase (e.g., 'CD5L', 'GBP2B')\n"
                        f"2. Available genes (first 10): {examples}\n"
                        f"3. Available annotations (first 10): {obs_examples}\n\n"
                        f"📋 CONTEXT: Dataset has {adata.n_vars:,} genes and "
                        f"{len(adata.obs.columns)} annotations."
                    )

                # Create spatial plot
                # For 10x Visium data - check if any sample has tissue images
                has_tissue_image = False
                if "spatial" in adata.uns:
                    for key in adata.uns["spatial"].keys():
                        if "images" in adata.uns["spatial"][key]:
                            has_tissue_image = True
                            break

                if has_tissue_image:

                    # Get sample key for library_id
                    sample_key = (
                        list(adata.uns["spatial"].keys())[0]
                        if "spatial" in adata.uns
                        else None
                    )

                    # Build common kwargs for sc.pl.spatial with beautification parameters
                    # NOTE: Do NOT use return_fig=True as it breaks spot rendering!
                    spatial_kwargs = {
                        "img_key": "hires",
                        "show": False,
                        # NO return_fig parameter - it breaks spot rendering!
                        "alpha": params.alpha,  # Spot transparency
                        "alpha_img": params.alpha_img,  # Dim background for better visibility
                        "frameon": params.show_axes,
                    }

                    # Only add spot_size if explicitly set (None = auto-determined, recommended)
                    if params.spot_size is not None:
                        spatial_kwargs["spot_size"] = params.spot_size

                    # Add library_id for multi-sample reliability
                    if sample_key:
                        spatial_kwargs["library_id"] = sample_key

                    # Add outline parameters if requested
                    if params.add_outline:
                        spatial_kwargs["na_in_legend"] = False

                    # With tissue image background
                    if feature:
                        if feature in adata.var_names:
                            # Gene expression (continuous data)
                            spatial_kwargs["cmap"] = params.colormap or "viridis"
                            sc.pl.spatial(adata, color=feature, **spatial_kwargs)
                        elif feature in adata.obs.columns:
                            # Observation annotation (categorical data)
                            _ensure_categorical_colors(adata, feature)
                            # Use palette for categorical data
                            if pd.api.types.is_categorical_dtype(adata.obs[feature]):
                                spatial_kwargs["palette"] = params.colormap or "Set2"
                            else:
                                spatial_kwargs["cmap"] = params.colormap or "viridis"
                            sc.pl.spatial(adata, color=feature, **spatial_kwargs)
                    else:
                        # Just tissue image with spots
                        sc.pl.spatial(adata, **spatial_kwargs)

                    # Get figure using plt.gcf() since we didn't use return_fig
                    fig = plt.gcf()

                    # Add outline/contour overlay if requested
                    if (
                        params.add_outline
                        and params.outline_cluster_key
                        and params.outline_cluster_key in adata.obs.columns
                    ):
                        if context:
                            await context.info(
                                f"Adding cluster outline overlay for {params.outline_cluster_key}"
                            )
                        try:
                            # Get the first (and usually only) axis from the figure
                            ax = fig.get_axes()[0] if fig.get_axes() else None
                            if ax:
                                # Use scanpy's add_outline functionality
                                sc.pl.spatial(
                                    adata,
                                    img_key="hires",
                                    color=params.outline_cluster_key,
                                    add_outline=True,
                                    outline_color=params.outline_color,
                                    outline_width=params.outline_width,
                                    show=False,
                                    ax=ax,
                                    legend_loc=None,
                                    colorbar=False,
                                )
                        except Exception as e:
                            if context:
                                await context.warning(
                                    f"Failed to add outline overlay: {str(e)}"
                                )

                else:
                    # For other spatial data without tissue image
                    fig, ax = create_figure()
                    if feature:
                        if feature in adata.var_names:
                            # Gene expression
                            sc.pl.embedding(
                                adata,
                                basis="spatial",
                                color=feature,
                                cmap=params.colormap,
                                show=False,
                                ax=ax,
                            )
                        elif feature in adata.obs.columns:
                            # Observation annotation (like clusters or cell types)
                            # Ensure categorical features have proper colors
                            _ensure_categorical_colors(adata, feature)
                            sc.pl.embedding(
                                adata, basis="spatial", color=feature, show=False, ax=ax
                            )
                    else:
                        # Just spatial coordinates
                        sc.pl.embedding(adata, basis="spatial", show=False, ax=ax)
                        ax.set_aspect("equal")
                        ax.set_title("Spatial coordinates")

                    # Add outline/contour overlay for non-tissue data if requested
                    if (
                        params.add_outline
                        and params.outline_cluster_key
                        and params.outline_cluster_key in adata.obs.columns
                    ):
                        if context:
                            await context.info(
                                f"Adding cluster outline for {params.outline_cluster_key}"
                            )
                        try:
                            # For non-tissue data, we can use scanpy's embedding plot with groups to highlight boundaries

                            # Get spatial coordinates
                            spatial_coords = adata.obsm["spatial"]
                            cluster_labels = adata.obs[
                                params.outline_cluster_key
                            ].values

                            # For each unique cluster, create a boundary
                            for cluster in np.unique(cluster_labels):
                                if pd.isna(cluster):
                                    continue
                                cluster_mask = cluster_labels == cluster
                                cluster_coords = spatial_coords[cluster_mask]

                                if (
                                    len(cluster_coords) > 2
                                ):  # Need at least 3 points for a boundary
                                    try:
                                        # Create convex hull for boundary
                                        from scipy.spatial import ConvexHull

                                        hull = ConvexHull(cluster_coords)
                                        boundary_coords = cluster_coords[hull.vertices]

                                        # Plot boundary
                                        boundary_coords = np.vstack(
                                            [boundary_coords, boundary_coords[0]]
                                        )  # Close the polygon
                                        ax.plot(
                                            boundary_coords[:, 0],
                                            boundary_coords[:, 1],
                                            color=params.outline_color,
                                            linewidth=params.outline_width,
                                            alpha=0.8,
                                        )
                                    except Exception:
                                        # Fallback: just plot the cluster points with a different marker
                                        ax.scatter(
                                            cluster_coords[:, 0],
                                            cluster_coords[:, 1],
                                            s=20,
                                            facecolors="none",
                                            edgecolors=params.outline_color,
                                            linewidth=params.outline_width,
                                            alpha=0.6,
                                        )
                        except Exception as e:
                            if context:
                                await context.warning(
                                    f"Failed to add outline overlay: {str(e)}"
                                )

        elif params.plot_type == "umap":
            if context:
                await context.info(
                    f"Creating UMAP plot for {params.feature if params.feature else 'clusters'}"
                )

            # Check if UMAP has been computed
            if "X_umap" not in adata.obsm:
                # STRICT: NO FALLBACK - UMAP and PCA are fundamentally different algorithms

                # Check prerequisites for UMAP
                if "neighbors" not in adata.uns:
                    error_msg = (
                        "❌ UMAP visualization requires neighborhood graph.\n\n"
                        "🔧 SOLUTION:\n"
                        "Run preprocessing first:\n"
                        "1. sc.pp.neighbors(adata)\n"
                        "2. sc.tl.umap(adata)\n\n"
                        "📋 SCIENTIFIC INTEGRITY: UMAP requires a k-nearest neighbor graph "
                        "to preserve local manifold structure. Cannot proceed without it."
                    )
                    if context:
                        await context.error(error_msg)
                    raise ValueError(error_msg)

                # Try to compute UMAP with proper error handling
                if context:
                    await context.info(
                        "UMAP not found. Attempting to compute UMAP coordinates..."
                    )

                try:
                    # Clean data before computing UMAP to handle extreme values
                    if hasattr(adata.X, "data"):
                        # Sparse matrix
                        adata.X.data = np.nan_to_num(
                            adata.X.data, nan=0.0, posinf=0.0, neginf=0.0
                        )
                    else:
                        # Dense matrix
                        adata.X = np.nan_to_num(
                            adata.X, nan=0.0, posinf=0.0, neginf=0.0
                        )

                    # Compute UMAP (scanpy already imported globally)
                    sc.tl.umap(adata)

                    if context:
                        await context.info("Successfully computed UMAP coordinates")

                except Exception as e:
                    # NO FALLBACK: Honest error reporting with alternatives
                    error_msg = (
                        "❌ Failed to compute UMAP for visualization.\n\n"
                        f"Error: {str(e)}\n\n"
                        "🔧 ALTERNATIVES:\n"
                        "1. Use PCA visualization instead (LINEAR method):\n"
                        "   visualize_data(data_id, params={'plot_type': 'pca'})\n\n"
                        "2. Use t-SNE visualization (NON-LINEAR, different from UMAP):\n"
                        "   visualize_data(data_id, params={'plot_type': 'tsne'})\n\n"
                        "3. Fix UMAP computation in preprocessing:\n"
                        "   - Ensure data is properly normalized\n"
                        "   - Check for extreme values or NaN\n"
                        "   - Verify neighbors graph is computed correctly\n"
                        "   - Try different UMAP parameters\n\n"
                        "📋 SCIENTIFIC INTEGRITY: UMAP (Uniform Manifold Approximation) and PCA "
                        "(Principal Component Analysis) are fundamentally different algorithms:\n"
                        "• UMAP: Non-linear, preserves local structure, reveals clusters\n"
                        "• PCA: Linear, preserves global variance, shows major axes of variation\n"
                        "We cannot substitute one for another without explicit user consent."
                    )
                    if context:
                        await context.error(error_msg)
                    raise RuntimeError(error_msg)

            # Check if we should create multi-panel plot for multiple features
            feature_list = (
                params.feature
                if isinstance(params.feature, list)
                else ([params.feature] if params.feature else [])
            )
            if feature_list and params.multi_panel and len(feature_list) > 1:
                # Use the new dedicated function for multi-gene UMAP
                fig = await create_multi_gene_umap_visualization(adata, params, context)
            else:
                # REQUIRE explicit feature specification - no automatic defaults
                if not feature_list:
                    categorical_cols = [
                        col
                        for col in adata.obs.columns
                        if adata.obs[col].dtype.name in ["object", "category"]
                    ]
                    raise ValueError(
                        "❌ UMAP visualization requires 'feature' parameter.\n\n"
                        f"Available categorical columns ({len(categorical_cols)} total):\n"
                        f"  {', '.join(categorical_cols[:15])}\n\n"
                        "🔧 SOLUTION: Specify feature explicitly:\n"
                        "  visualize_data(data_id, params={'plot_type': 'umap', 'feature': 'your_column_name'})"
                    )

                single_feature = feature_list[0]
                # Directly call the helper, it will handle all cases including deconvolution keys
                feature = await validate_and_prepare_feature(
                    adata, single_feature, context, default_feature=None
                )

                # Create UMAP plot with potential dual encoding (color + size)
                plot_kwargs = {"show": False, "return_fig": True}

                if feature:
                    plot_kwargs["color"] = feature
                    if feature in adata.var_names:
                        plot_kwargs["cmap"] = params.colormap
                    elif feature in adata.obs.columns:
                        # Ensure categorical features have proper colors
                        _ensure_categorical_colors(adata, feature)

                # Add size encoding if requested
                if params.size_by:
                    if (
                        params.size_by in adata.var_names
                        or params.size_by in adata.obs.columns
                    ):
                        # For scanpy compatibility, we need to pass size values correctly
                        if params.size_by in adata.var_names:
                            # Gene expression - get the values
                            size_values = adata[:, params.size_by].X
                            if hasattr(size_values, "toarray"):
                                size_values = size_values.toarray().flatten()
                            else:
                                size_values = size_values.flatten()
                            # Normalize to reasonable size range (10-100)
                            size_values = 10 + 90 * (
                                size_values - size_values.min()
                            ) / (size_values.max() - size_values.min() + 1e-8)
                            plot_kwargs["size"] = size_values
                        else:
                            # Observation column
                            if adata.obs[params.size_by].dtype.name in [
                                "category",
                                "object",
                            ]:
                                # Categorical - use different sizes for different categories
                                unique_vals = adata.obs[params.size_by].unique()
                                size_map = {
                                    val: 20 + i * 15
                                    for i, val in enumerate(unique_vals)
                                }
                                plot_kwargs["size"] = (
                                    adata.obs[params.size_by].map(size_map).values
                                )
                            else:
                                # Numeric - normalize to size range
                                size_values = adata.obs[params.size_by].values
                                size_values = 10 + 90 * (
                                    size_values - size_values.min()
                                ) / (size_values.max() - size_values.min() + 1e-8)
                                plot_kwargs["size"] = size_values

                        if context:
                            await context.info(
                                f"Using dual encoding: color={feature}, size={params.size_by}"
                            )
                    elif context:
                        await context.warning(
                            f"Size feature '{params.size_by}' not found in data"
                        )

                # Create the plot
                fig = sc.pl.umap(adata, **plot_kwargs)

                # Add velocity or trajectory overlays if requested
                if params.show_velocity or params.show_trajectory:
                    try:
                        # Get the axis from the figure
                        ax = fig.get_axes()[0] if fig.get_axes() else None
                        if ax:
                            # Velocity overlay
                            if params.show_velocity:
                                if "velocity_umap" in adata.obsm:
                                    try:
                                        # Try to use scvelo if available
                                        import scvelo as scv

                                        scv.pl.velocity_embedding(
                                            adata,
                                            basis="umap",
                                            ax=ax,
                                            show=False,
                                            arrow_length=params.velocity_scale,
                                            arrow_size=params.velocity_scale,
                                        )
                                        if context:
                                            await context.info(
                                                "Added RNA velocity vectors to UMAP"
                                            )
                                    except ImportError:
                                        if context:
                                            await context.warning(
                                                "scvelo not available for velocity overlay"
                                            )
                                    except Exception as e:
                                        if context:
                                            await context.warning(
                                                f"Failed to add velocity overlay: {str(e)}"
                                            )
                                elif context:
                                    await context.warning(
                                        "Velocity data (velocity_umap) not found in adata.obsm"
                                    )

                            # Trajectory overlay (PAGA connections)
                            if params.show_trajectory:
                                if "paga" in adata.uns:
                                    try:
                                        # Add PAGA trajectory connections

                                        # Get PAGA connectivity matrix
                                        paga_adj = adata.uns["paga"][
                                            "connectivities"
                                        ].toarray()

                                        # Get cluster centroids in UMAP space
                                        if feature and feature in adata.obs.columns:
                                            umap_coords = adata.obsm["X_umap"]
                                            clusters = adata.obs[feature].astype(str)
                                            unique_clusters = clusters.unique()

                                            # Calculate centroids
                                            centroids = {}
                                            for cluster in unique_clusters:
                                                mask = clusters == cluster
                                                if mask.sum() > 0:
                                                    centroids[cluster] = umap_coords[
                                                        mask
                                                    ].mean(axis=0)

                                            # Draw connections based on PAGA
                                            threshold = (
                                                np.percentile(
                                                    paga_adj[paga_adj > 0], 75
                                                )
                                                if np.any(paga_adj > 0)
                                                else 0
                                            )
                                            for i, cluster_i in enumerate(
                                                unique_clusters
                                            ):
                                                for j, cluster_j in enumerate(
                                                    unique_clusters
                                                ):
                                                    if (
                                                        i < j
                                                        and i < paga_adj.shape[0]
                                                        and j < paga_adj.shape[1]
                                                    ):
                                                        if paga_adj[i, j] > threshold:
                                                            if (
                                                                cluster_i in centroids
                                                                and cluster_j
                                                                in centroids
                                                            ):
                                                                x_coords = [
                                                                    centroids[
                                                                        cluster_i
                                                                    ][0],
                                                                    centroids[
                                                                        cluster_j
                                                                    ][0],
                                                                ]
                                                                y_coords = [
                                                                    centroids[
                                                                        cluster_i
                                                                    ][1],
                                                                    centroids[
                                                                        cluster_j
                                                                    ][1],
                                                                ]
                                                                ax.plot(
                                                                    x_coords,
                                                                    y_coords,
                                                                    "k-",
                                                                    alpha=0.6,
                                                                    linewidth=2,
                                                                )

                                            if context:
                                                await context.info(
                                                    "Added PAGA trajectory connections to UMAP"
                                                )
                                    except Exception as e:
                                        if context:
                                            await context.warning(
                                                f"Failed to add trajectory overlay: {str(e)}"
                                            )
                                elif context:
                                    await context.warning(
                                        "PAGA trajectory data not found in adata.uns"
                                    )
                    except Exception as e:
                        if context:
                            await context.warning(f"Failed to add overlays: {str(e)}")

        elif params.plot_type == "heatmap":
            if context:
                await context.info("Creating heatmap plot")

            # For heatmap, we need highly variable genes for meaningful clustering
            hvg_exists = "highly_variable" in adata.var
            hvg_any = adata.var["highly_variable"].any() if hvg_exists else False

            if not hvg_exists or not hvg_any:
                # Provide detailed diagnostic information
                diagnostic_info = "\n📋 DIAGNOSTIC INFO:\n"
                diagnostic_info += (
                    f"  - 'highly_variable' column exists: {hvg_exists}\n"
                )
                if hvg_exists:
                    diagnostic_info += f"  - Number of HVGs marked: {adata.var['highly_variable'].sum()}\n"
                    diagnostic_info += f"  - Total genes in dataset: {adata.n_vars}\n"
                    diagnostic_info += (
                        f"  - HVG column dtype: {adata.var['highly_variable'].dtype}\n"
                    )
                diagnostic_info += (
                    f"  - Available var columns: {list(adata.var.columns)[:10]}\n"
                )

                raise ValueError(
                    "❌ Heatmap requires highly variable genes but none found.\n\n"
                    "🔧 SOLUTIONS:\n"
                    "1. Run preprocessing: preprocess_data(data_id, params={'n_top_genes': 2000})\n"
                    "2. Specify genes: visualize_data(data_id, params={'feature': ['CD3D']})\n"
                    "3. Use spatial plot: visualize_data(data_id, params={'plot_type': 'spatial'})\n"
                    + diagnostic_info
                )

            # Create heatmap of top genes across groups
            # REQUIRE explicit cluster_key specification for grouping
            if not params.cluster_key:
                categorical_cols = [
                    col
                    for col in adata.obs.columns
                    if adata.obs[col].dtype.name in ["object", "category"]
                ]
                raise ValueError(
                    "❌ Heatmap visualization requires 'cluster_key' parameter.\n\n"
                    f"Available categorical columns ({len(categorical_cols)} total):\n"
                    f"  {', '.join(categorical_cols[:15])}\n\n"
                    "🔧 SOLUTION: Specify cluster_key explicitly:\n"
                    "  visualize_data(data_id, params={'plot_type': 'heatmap', 'cluster_key': 'your_column_name'})\n\n"
                    "💡 NOTE: ChatSpatial uses 'cluster_key' (not 'groupby' as in Scanpy).\n"
                    "   This maintains consistency with Squidpy spatial analysis functions."
                )

            if params.cluster_key not in adata.obs.columns:
                raise ValueError(
                    f"❌ Cluster key '{params.cluster_key}' not found in data.\n\n"
                    f"Available columns: {', '.join(list(adata.obs.columns)[:20])}"
                )

            groupby = params.cluster_key
            # Limit the number of groups to avoid oversized responses
            n_groups = len(adata.obs[groupby].cat.categories)
            if n_groups > 10:
                if context:
                    await context.warning(
                        f"Too many groups ({n_groups}). Limiting to 10 groups."
                    )
                # Get the 10 largest groups
                group_counts = adata.obs[groupby].value_counts().nlargest(10).index
                # Subset the data to include only these groups
                adata = adata[adata.obs[groupby].isin(group_counts)].copy()

            # Limit the number of genes to avoid oversized responses
            # Use a smaller number (20) instead of 50 to reduce response size
            n_genes = 20
            if context:
                await context.info(
                    f"Using top {n_genes} highly variable genes for heatmap"
                )

            # Get genes for heatmap
            feature_list = (
                params.feature
                if isinstance(params.feature, list)
                else ([params.feature] if params.feature else [])
            )
            if feature_list and len(feature_list) > 0:
                # Use user-specified genes
                available_genes = [
                    gene for gene in feature_list if gene in adata.var_names
                ]
                if not available_genes:
                    if context:
                        await context.warning(
                            f"None of specified genes found: {feature_list}. Using highly variable genes."
                        )
                    # Fall back to highly variable genes
                    available_genes = adata.var_names[adata.var.highly_variable][
                        :n_genes
                    ].tolist()
                else:
                    # Limit to reasonable number for visualization
                    available_genes = available_genes[:n_genes]
            else:
                # Use highly variable genes
                available_genes = adata.var_names[adata.var.highly_variable][
                    :n_genes
                ].tolist()

            if context:
                await context.info(
                    f"Creating heatmap with {len(available_genes)} genes: {available_genes[:5]}..."
                )

            # Create heatmap with improved settings
            plt.figure(
                figsize=(
                    max(8, len(available_genes) * 0.5),
                    (
                        max(6, len(adata.obs[groupby].cat.categories) * 0.3)
                        if groupby
                        else 6
                    ),
                )
            )

            try:
                # Use scanpy's heatmap with better parameters and enhanced annotations
                heatmap_kwargs = {
                    "var_names": available_genes,
                    "groupby": groupby,
                    "cmap": params.colormap,
                    "show": False,
                    "dendrogram": False,
                    "standard_scale": "var",  # Standardize genes (rows)
                    "figsize": None,  # Let matplotlib handle the size
                    "swap_axes": False,  # Keep genes as rows
                    "vmin": None,  # Let scanpy determine range
                    "vmax": None,
                }

                # Add obs annotations if specified
                if params.obs_annotation:
                    available_obs_annotations = [
                        col for col in params.obs_annotation if col in adata.obs.columns
                    ]
                    if available_obs_annotations:
                        # For scanpy heatmap, we can use the var_group_labels parameter
                        # But scanpy heatmap has limited annotation support, so we'll enhance post-creation
                        if context:
                            await context.info(
                                f"Will add obs annotations: {available_obs_annotations}"
                            )

                sc.pl.heatmap(adata, **heatmap_kwargs)
                fig = plt.gcf()

                # Add custom annotations if requested
                if params.obs_annotation or params.var_annotation:
                    try:
                        # Get the main heatmap axis
                        axes = fig.get_axes()
                        if axes:
                            main_ax = axes[
                                0
                            ]  # Usually the main heatmap is the first axis

                            # Add obs annotations (column annotations)
                            if params.obs_annotation:
                                available_obs = [
                                    col
                                    for col in params.obs_annotation
                                    if col in adata.obs.columns
                                ]
                                if available_obs:
                                    # Create a simple annotation bar above the heatmap
                                    # Get the position of the main heatmap
                                    pos = main_ax.get_position()

                                    # Create annotation axis above the heatmap
                                    ann_height = 0.02 * len(available_obs)
                                    ann_ax = fig.add_axes(
                                        [pos.x0, pos.y1 + 0.01, pos.width, ann_height]
                                    )

                                    # Create annotation data
                                    ann_data = adata.obs[available_obs].copy()

                                    # Convert categorical to numeric for visualization
                                    for col in available_obs:
                                        if (
                                            ann_data[col].dtype == "category"
                                            or ann_data[col].dtype == "object"
                                        ):
                                            unique_vals = ann_data[col].unique()
                                            color_map = {
                                                val: i
                                                for i, val in enumerate(unique_vals)
                                            }
                                            ann_data[col] = ann_data[col].map(color_map)

                                    # Create annotation heatmap
                                    im = ann_ax.imshow(
                                        ann_data.T,
                                        aspect="auto",
                                        cmap="Set3",
                                        interpolation="nearest",
                                    )
                                    ann_ax.set_xlim(main_ax.get_xlim())
                                    ann_ax.set_xticks([])
                                    ann_ax.set_yticks(range(len(available_obs)))
                                    ann_ax.set_yticklabels(available_obs, fontsize=8)
                                    ann_ax.tick_params(left=False)

                            # Add var annotations (row annotations) - similar approach
                            if params.var_annotation:
                                available_var = [
                                    col
                                    for col in params.var_annotation
                                    if col in adata.var.columns
                                ]
                                if available_var and context:
                                    await context.info(
                                        f"Adding var annotations: {available_var}"
                                    )
                                    # This would require similar logic for row annotations

                    except Exception as e:
                        if context:
                            await context.warning(
                                f"Failed to add annotations: {str(e)}"
                            )

                # Improve the layout
                plt.tight_layout()

            except Exception as e:
                if context:
                    await context.warning(
                        f"Scanpy heatmap failed: {e}. Creating custom heatmap..."
                    )

                # Fallback: create custom heatmap using seaborn
                # Note: pandas and seaborn already imported at module level
                import seaborn as sns

                # Get expression data
                expr_data = adata[:, available_genes].X
                if hasattr(expr_data, "toarray"):
                    expr_data = expr_data.toarray()

                # Create DataFrame
                expr_df = pd.DataFrame(
                    expr_data.T,  # Transpose so genes are rows
                    index=available_genes,
                    columns=adata.obs.index,
                )

                # Add group information if available
                if groupby:
                    # Sort by groups
                    group_order = adata.obs[groupby].cat.categories
                    sorted_cells = []
                    for group in group_order:
                        group_cells = adata.obs[adata.obs[groupby] == group].index
                        sorted_cells.extend(group_cells)

                    expr_df = expr_df[sorted_cells]

                # Create heatmap
                fig, ax = plt.subplots(
                    figsize=(
                        max(8, len(available_genes) * 0.5),
                        max(6, len(available_genes) * 0.3),
                    )
                )

                # Standardize data (z-score normalization across genes)
                from scipy.stats import zscore

                expr_df_norm = expr_df.apply(zscore, axis=1)

                # Create heatmap with seaborn
                sns.heatmap(
                    expr_df_norm,
                    cmap=params.colormap,
                    center=0,
                    robust=True,
                    ax=ax,
                    cbar_kws={"shrink": 0.8},
                    xticklabels=False,  # Too many cells to show labels
                    yticklabels=True,
                )

                ax.set_title(f"Gene Expression Heatmap ({len(available_genes)} genes)")
                ax.set_xlabel("Cells")
                ax.set_ylabel("Genes")

                plt.tight_layout()

        elif params.plot_type == "violin":
            if context:
                await context.info(
                    f"Creating violin plot for {params.feature if params.feature else 'top genes'}"
                )

            # Normalize feature to list format (same logic as heatmap)
            feature_list = (
                params.feature
                if isinstance(params.feature, list)
                else ([params.feature] if params.feature else [])
            )

            # Use specified genes or default to first 3 genes
            if feature_list:
                # Check which genes exist in dataset
                available_genes = [
                    gene for gene in feature_list if gene in adata.var_names
                ]
                missing_genes = [
                    gene for gene in feature_list if gene not in adata.var_names
                ]

                if not available_genes:
                    # None of the specified genes found
                    examples = list(adata.var_names[:10])
                    raise DataNotFoundError(
                        f"❌ Genes not found in dataset: {missing_genes}\n\n"
                        f"🔧 SOLUTIONS:\n"
                        f"1. Choose from available genes (examples: {examples})\n"
                        f"2. Check gene name format (e.g., 'CD3D' vs 'Cd3d')\n"
                        f"3. Use gene search functionality\n\n"
                        f"📋 CONTEXT: Dataset contains {adata.n_vars:,} genes total."
                    )
                elif missing_genes and context:
                    # Some genes missing - warn but continue
                    await context.warning(
                        f"Genes not found (skipping): {missing_genes}"
                    )

                genes = available_genes
            else:
                # No genes specified, use first 3 genes
                genes = list(adata.var_names[:3])

            # REQUIRE explicit cluster_key specification for grouping
            if not params.cluster_key:
                categorical_cols = [
                    col
                    for col in adata.obs.columns
                    if adata.obs[col].dtype.name in ["object", "category"]
                ]
                raise ValueError(
                    "❌ Violin plot requires 'cluster_key' parameter.\n\n"
                    f"Available categorical columns ({len(categorical_cols)} total):\n"
                    f"  {', '.join(categorical_cols[:15])}\n\n"
                    "🔧 SOLUTION: Specify cluster_key explicitly:\n"
                    "  visualize_data(data_id, params={'plot_type': 'violin', 'cluster_key': 'your_column_name'})\n\n"
                    "💡 NOTE: ChatSpatial uses 'cluster_key' (not 'groupby' as in Scanpy).\n"
                    "   This maintains consistency with Squidpy spatial analysis functions."
                )

            if params.cluster_key not in adata.obs.columns:
                raise ValueError(
                    f"❌ Cluster key '{params.cluster_key}' not found in data.\n\n"
                    f"Available columns: {', '.join(list(adata.obs.columns)[:20])}"
                )

            groupby = params.cluster_key
            # Limit the number of groups to avoid oversized responses
            n_groups = len(adata.obs[groupby].cat.categories)
            if n_groups > 8:
                if context:
                    await context.warning(
                        f"Too many groups ({n_groups}). Limiting to 8 groups."
                    )
                # Get the 8 largest groups
                group_counts = adata.obs[groupby].value_counts().nlargest(8).index
                # Subset the data to include only these groups
                adata = adata[adata.obs[groupby].isin(group_counts)].copy()

            if False:  # Disabled the no-grouping path - cluster_key now required
                groupby = None

            # Create violin plot with smaller figure size
            # Reduce figure size from (12, 10) to (8, 6) to decrease image size
            plt.figure(figsize=(8, 6))
            ax = sc.pl.violin(
                adata,
                genes,
                groupby=groupby,
                show=False,
                jitter=0.2,  # Reduce jitter amount
                scale="width",
            )  # Scale violins to same width
            fig = plt.gcf()  # Get the current figure

            # Fix overlapping x-axis labels
            for ax_item in fig.get_axes():
                # Rotate x-axis labels 45 degrees for better readability
                ax_item.tick_params(axis="x", labelrotation=45)
                # Align labels to the right for better appearance
                for label in ax_item.get_xticklabels():
                    label.set_horizontalalignment("right")

            # Adjust layout to prevent label cutoff
            plt.tight_layout()

        elif params.plot_type == "deconvolution":
            if context:
                await context.info("Creating deconvolution visualization")
            fig = await create_deconvolution_visualization(adata, params, context)

        elif params.plot_type == "spatial_domains":
            if context:
                await context.info("Creating spatial domains visualization")
            fig = await create_spatial_domains_visualization(adata, params, context)

        elif params.plot_type == "cell_communication":
            if context:
                await context.info("Creating cell communication visualization")
            fig = await create_cell_communication_visualization(adata, params, context)

        elif params.plot_type == "multi_gene":
            if context:
                await context.info("Creating multi-gene visualization")
            fig = await create_multi_gene_visualization(adata, params, context)

        elif params.plot_type == "lr_pairs":
            if context:
                await context.info("Creating ligand-receptor pairs visualization")
            fig = await create_lr_pairs_visualization(adata, params, context)

        elif params.plot_type == "gene_correlation":
            if context:
                await context.info("Creating gene correlation visualization")
            fig = await create_gene_correlation_visualization(adata, params, context)

        elif params.plot_type == "rna_velocity":
            if context:
                await context.info("Creating RNA velocity visualization")
            fig = await create_rna_velocity_visualization(adata, params, context)

        elif params.plot_type == "trajectory":
            if context:
                await context.info("Creating trajectory visualization")
            fig = await create_trajectory_visualization(adata, params, context)

        elif params.plot_type == "spatial_statistics":
            if context:
                await context.info("Creating spatial statistics visualization")
            fig = await create_spatial_statistics_visualization(adata, params, context)

        elif params.plot_type == "spatial_enrichment":
            if context:
                await context.info("Creating spatial enrichment visualization")
            fig = await create_enrichment_visualization(adata, params, context)

        elif params.plot_type == "pathway_enrichment":
            if context:
                await context.info("Creating pathway enrichment visualization")
            fig = await create_gsea_visualization(adata, params, context)

        elif params.plot_type == "spatial_interaction":
            if context:
                await context.info("Creating spatial interaction visualization")
            fig = await create_spatial_interaction_visualization(adata, params, context)

        elif params.plot_type == "batch_integration":
            if context:
                await context.info("Creating batch integration quality visualization")
            fig = await create_batch_integration_visualization(adata, params, context)

        elif params.plot_type == "cnv_heatmap":
            if context:
                await context.info("Creating CNV heatmap visualization")

            # Auto-detect CNV data source (infercnvpy or Numbat)
            cnv_method = None

            if "X_cnv" in adata.obsm:
                cnv_method = "infercnvpy"
            elif "X_cnv_numbat" in adata.obsm:
                cnv_method = "numbat"
            else:
                error_msg = (
                    "CNV data not found. Expected 'X_cnv' (infercnvpy) or "
                    "'X_cnv_numbat' (Numbat) in adata.obsm. "
                    "Please run CNV analysis first using analyze_cnv()."
                )
                if context:
                    await context.warning(error_msg)
                raise DataNotFoundError(error_msg)

            if context:
                await context.info(f"Detected CNV data from {cnv_method} method")

            # Check if infercnvpy is available (needed for visualization)
            if not INFERCNVPY_AVAILABLE:
                error_msg = (
                    "infercnvpy is not installed. Please install it with:\n"
                    "  pip install 'chatspatial[cnv]'\n"
                    "or:\n"
                    "  pip install infercnvpy>=0.4.0"
                )
                if context:
                    await context.warning(error_msg)
                raise DataCompatibilityError(error_msg)

            # For Numbat data, temporarily copy to X_cnv for visualization
            if cnv_method == "numbat":
                if context:
                    await context.info(
                        "Converting Numbat CNV data to infercnvpy format for visualization"
                    )
                adata.obsm["X_cnv"] = adata.obsm["X_cnv_numbat"]
                # Also ensure cnv metadata exists for infercnvpy plotting
                # If chromosome info is not available, infercnvpy will still plot but without chromosome labels
                if "cnv" not in adata.uns:
                    # Create minimal metadata structure for infercnvpy
                    adata.uns["cnv"] = {
                        "genomic_positions": False,  # No genomic position info available
                    }
                    if context:
                        await context.info(
                            "Note: Chromosome labels not available for Numbat heatmap. "
                            "Install R packages for full chromosome annotation."
                        )

            # Check if CNV metadata exists
            if "cnv" not in adata.uns:
                error_msg = (
                    "CNV metadata not found in adata.uns['cnv']. "
                    "The CNV analysis may not have completed properly. "
                    "Please re-run analyze_cnv()."
                )
                if context:
                    await context.warning(error_msg)
                raise DataNotFoundError(error_msg)

            # Create CNV heatmap
            if context:
                await context.info("Generating CNV heatmap...")

            figsize = params.figure_size if params.figure_size else (12, 8)

            # For Numbat data without chromosome info, use aggregated heatmap by group
            if cnv_method == "numbat" and "chromosome" not in adata.var.columns:
                if context:
                    await context.info(
                        "Creating aggregated CNV heatmap by group (chromosome positions not available)"
                    )

                # Create aggregated heatmap of CNV matrix
                import seaborn as sns

                # Get CNV matrix
                cnv_matrix = adata.obsm["X_cnv"]

                # Aggregate by feature (e.g., clone) for cleaner visualization
                if params.feature and params.feature in adata.obs.columns:
                    # Group cells by feature and compute mean CNV per group
                    feature_values = adata.obs[params.feature]
                    unique_groups = sorted(feature_values.unique())

                    # Compute mean CNV for each group
                    aggregated_cnv = []
                    group_labels = []
                    group_sizes = []

                    for group in unique_groups:
                        group_mask = feature_values == group
                        group_cnv = cnv_matrix[group_mask, :].mean(axis=0)
                        aggregated_cnv.append(group_cnv)
                        group_labels.append(str(group))
                        group_sizes.append(group_mask.sum())

                    aggregated_cnv = np.array(aggregated_cnv)

                    # Calculate appropriate figure width based on number of bins
                    # Use narrower bins for better visualization (0.004 inches per bin)
                    n_bins = aggregated_cnv.shape[1]
                    fig_width = min(max(6, n_bins * 0.004), 12)  # Between 6-12 inches
                    fig_height = max(4, len(unique_groups) * 1.2)

                    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                    # Plot aggregated heatmap with fixed aspect ratio
                    im = ax.imshow(
                        aggregated_cnv,
                        cmap="RdBu_r",
                        aspect="auto",
                        vmin=-1,
                        vmax=1,
                        interpolation="nearest",
                    )

                    # Add colorbar
                    cbar = plt.colorbar(im, ax=ax, label="Mean CNV state")

                    # Set y-axis labels with group names and cell counts
                    ax.set_yticks(range(len(group_labels)))
                    ax.set_yticklabels(
                        [
                            f"{label} (n={size})"
                            for label, size in zip(group_labels, group_sizes)
                        ]
                    )
                    ax.set_ylabel(params.feature, fontsize=12, fontweight="bold")

                    # Set x-axis
                    ax.set_xlabel("Genomic position (binned)", fontsize=12)
                    ax.set_xticks([])  # Hide x-axis ticks for cleaner look

                    # Add title
                    ax.set_title(
                        f"CNV Profile by {params.feature}\n(Numbat analysis, aggregated by group)",
                        fontsize=14,
                        fontweight="bold",
                    )

                    # Add gridlines between groups
                    for i in range(len(group_labels) + 1):
                        ax.axhline(i - 0.5, color="white", linewidth=2)

                else:
                    # No grouping - show warning and plot all cells (not recommended)
                    fig, ax = plt.subplots(figsize=figsize)

                    sns.heatmap(
                        cnv_matrix,
                        cmap="RdBu_r",
                        center=0,
                        cbar_kws={"label": "CNV state"},
                        yticklabels=False,
                        xticklabels=False,
                        ax=ax,
                        vmin=-1,
                        vmax=1,
                    )

                    ax.set_xlabel("Genomic position (binned)")
                    ax.set_ylabel("Cells")
                    ax.set_title("CNV Heatmap (Numbat)\nAll cells (ungrouped)")

                plt.tight_layout()

            else:
                # Use infercnvpy chromosome_heatmap for infercnvpy data or Numbat with chr info
                if context:
                    await context.info("Creating chromosome-organized CNV heatmap...")

                cnv.pl.chromosome_heatmap(
                    adata,
                    groupby=params.feature if params.feature else None,
                    dendrogram=True,
                    show=False,
                    figsize=figsize,
                )
                # Get current figure
                fig = plt.gcf()

            if context:
                await context.info("CNV heatmap created successfully")

        elif params.plot_type == "spatial_cnv":
            if context:
                await context.info("Creating spatial CNV projection visualization")

            # Check if spatial coordinates exist
            if "spatial" not in adata.obsm:
                error_msg = (
                    "Spatial coordinates not found in adata.obsm['spatial']. "
                    "This plot type requires spatial transcriptomics data."
                )
                if context:
                    await context.warning(error_msg)
                raise DataNotFoundError(error_msg)

            # Determine feature to visualize
            feature_to_plot = params.feature

            # Auto-detect CNV-related features if none specified
            if not feature_to_plot:
                if "numbat_clone" in adata.obs:
                    feature_to_plot = "numbat_clone"
                    if context:
                        await context.info(
                            "No feature specified, using 'numbat_clone' (Numbat clone assignment)"
                        )
                elif "cnv_score" in adata.obs:
                    feature_to_plot = "cnv_score"
                    if context:
                        await context.info(
                            "No feature specified, using 'cnv_score' (CNV score)"
                        )
                elif "numbat_p_cnv" in adata.obs:
                    feature_to_plot = "numbat_p_cnv"
                    if context:
                        await context.info(
                            "No feature specified, using 'numbat_p_cnv' (Numbat CNV probability)"
                        )
                else:
                    error_msg = (
                        "No CNV-related features found in adata.obs. "
                        "Expected one of: 'numbat_clone', 'cnv_score', 'numbat_p_cnv'. "
                        "Please run CNV analysis first using analyze_cnv()."
                    )
                    if context:
                        await context.warning(error_msg)
                    raise DataNotFoundError(error_msg)

            # Validate feature exists
            if feature_to_plot not in adata.obs.columns:
                error_msg = (
                    f"Feature '{feature_to_plot}' not found in adata.obs. "
                    f"Available CNV features: {[col for col in adata.obs.columns if 'cnv' in col.lower() or 'numbat' in col.lower()]}"
                )
                if context:
                    await context.warning(error_msg)
                raise DataNotFoundError(error_msg)

            if context:
                await context.info(
                    f"Visualizing {feature_to_plot} on spatial coordinates"
                )

            # Reuse spatial plot logic
            # Check for tissue images
            has_tissue_image = False
            if "spatial" in adata.uns:
                for key in adata.uns["spatial"].keys():
                    if "images" in adata.uns["spatial"][key]:
                        has_tissue_image = True
                        break

            if has_tissue_image:
                # With tissue image background
                sample_key = (
                    list(adata.uns["spatial"].keys())[0]
                    if "spatial" in adata.uns
                    else None
                )

                spatial_kwargs = {
                    "img_key": "hires",
                    "show": False,
                    "alpha": params.alpha,
                    "alpha_img": params.alpha_img,
                    "frameon": params.show_axes,
                }

                if params.spot_size is not None:
                    spatial_kwargs["spot_size"] = params.spot_size

                if sample_key:
                    spatial_kwargs["library_id"] = sample_key

                # Determine if categorical or continuous
                if (
                    pd.api.types.is_categorical_dtype(adata.obs[feature_to_plot])
                    or adata.obs[feature_to_plot].dtype == "object"
                ):
                    # Categorical data (e.g., clone assignments)
                    _ensure_categorical_colors(adata, feature_to_plot)
                    spatial_kwargs["palette"] = params.colormap or "tab20"
                else:
                    # Continuous data (e.g., CNV scores, probabilities)
                    spatial_kwargs["cmap"] = params.colormap or "RdBu_r"

                sc.pl.spatial(adata, color=feature_to_plot, **spatial_kwargs)
                fig = plt.gcf()

            else:
                # Without tissue image - use embedding plot
                figsize = params.figure_size if params.figure_size else (10, 8)
                fig, ax = plt.subplots(figsize=figsize)

                if (
                    pd.api.types.is_categorical_dtype(adata.obs[feature_to_plot])
                    or adata.obs[feature_to_plot].dtype == "object"
                ):
                    # Categorical
                    sc.pl.embedding(
                        adata,
                        basis="spatial",
                        color=feature_to_plot,
                        palette=params.colormap or "tab20",
                        show=False,
                        ax=ax,
                    )
                else:
                    # Continuous
                    sc.pl.embedding(
                        adata,
                        basis="spatial",
                        color=feature_to_plot,
                        cmap=params.colormap or "RdBu_r",
                        show=False,
                        ax=ax,
                    )

            if context:
                await context.info(
                    f"Spatial CNV projection created for {feature_to_plot}"
                )

        elif params.plot_type == "card_imputation":
            if context:
                await context.info("Creating CARD imputation visualization")

            # Check if CARD imputation data exists
            if "card_imputation" not in adata.uns:
                error_msg = (
                    "CARD imputation data not found in adata.uns['card_imputation']. "
                    "Please run CARD deconvolution with card_imputation=True first."
                )
                if context:
                    await context.warning(error_msg)
                raise DataNotFoundError(error_msg)

            # Extract imputation data
            impute_data = adata.uns["card_imputation"]
            imputed_proportions = impute_data["proportions"]
            imputed_coords = impute_data["coordinates"]

            # Determine what to visualize
            feature = params.feature
            if not feature:
                # Default: show dominant cell types
                feature = "dominant"

            figsize = params.figure_size if params.figure_size else (12, 10)
            fig, ax = plt.subplots(figsize=figsize)

            if feature == "dominant":
                # Show dominant cell types
                dominant_types = imputed_proportions.idxmax(axis=1)
                unique_types = dominant_types.unique()

                # Create color map
                import seaborn as sns

                colors = sns.color_palette("tab20", n_colors=len(unique_types))
                color_map = {ct: colors[i] for i, ct in enumerate(unique_types)}
                point_colors = [color_map[ct] for ct in dominant_types]

                scatter = ax.scatter(
                    imputed_coords["x"],
                    imputed_coords["y"],
                    c=point_colors,
                    s=25,
                    edgecolors="none",
                    alpha=0.7,
                )

                ax.set_title(
                    f"CARD Imputation: Dominant Cell Types\n"
                    f"({len(imputed_coords)} locations, "
                    f"{impute_data['resolution_increase']:.1f}x resolution)",
                    fontsize=14,
                    fontweight="bold",
                )

                # Create legend
                from matplotlib.patches import Patch

                legend_elements = [
                    Patch(facecolor=color_map[ct], label=ct)
                    for ct in sorted(unique_types)
                ]
                ax.legend(
                    handles=legend_elements,
                    bbox_to_anchor=(1.05, 1),
                    loc="upper left",
                    fontsize=9,
                )

            elif feature in imputed_proportions.columns:
                # Show specific cell type proportion
                scatter = ax.scatter(
                    imputed_coords["x"],
                    imputed_coords["y"],
                    c=imputed_proportions[feature],
                    s=30,
                    cmap=params.colormap or "viridis",
                    vmin=0,
                    vmax=imputed_proportions[feature].quantile(0.95),
                    edgecolors="none",
                    alpha=0.8,
                )

                ax.set_title(
                    f"CARD Imputation: {feature}\n"
                    f"(Mean: {imputed_proportions[feature].mean():.3f}, "
                    f"{len(imputed_coords)} locations)",
                    fontsize=14,
                    fontweight="bold",
                )

                # Add colorbar
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label("Proportion", fontsize=12)

            else:
                error_msg = (
                    f"Feature '{feature}' not found in imputation results. "
                    f"Available cell types: {list(imputed_proportions.columns)}\n"
                    f"Use feature='dominant' to show dominant cell types."
                )
                if context:
                    await context.warning(error_msg)
                raise DataNotFoundError(error_msg)

            ax.set_xlabel("X coordinate", fontsize=12)
            ax.set_ylabel("Y coordinate", fontsize=12)
            ax.set_aspect("equal")
            plt.tight_layout()

            if context:
                await context.info("CARD imputation visualization created successfully")

        else:
            # This should never happen due to parameter validation at the beginning
            error_msg = f"Unsupported plot type: {params.plot_type}"
            if context:
                await context.warning(error_msg)
            raise InvalidParameterError(error_msg)

        # Convert figure with optimization (preview + resource for large images)
        if context:
            await context.info(
                f"Converting {params.plot_type} figure with token optimization..."
            )
        # Use the optimized conversion function
        return await optimize_fig_to_image_with_cache(
            fig,
            params,
            context,
            data_id=data_id,
            plot_type=params.plot_type,
            mode="auto",
        )

    except Exception as e:
        # Make sure to close any open figures in case of error
        plt.close("all")

        # Log the error
        error_msg = f"Error in {params.plot_type} visualization: {str(e)}"
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {traceback.format_exc()}")

        # For image conversion errors, return a placeholder image
        if "fig_to_image" in str(e) or "convert" in str(e).lower():
            return create_placeholder_image(
                f"Error in {params.plot_type} visualization: {str(e)}"
            )

        # Wrap the error in a more informative exception
        if isinstance(
            e, (DataNotFoundError, InvalidParameterError, DataCompatibilityError)
        ):
            # Re-raise specific errors
            raise
        else:
            # Wrap generic errors
            raise ProcessingError(
                f"Failed to create {params.plot_type} visualization: {str(e)}"
            ) from e


@handle_visualization_errors("Deconvolution")
async def create_deconvolution_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create deconvolution results visualization

    Args:
        adata: AnnData object with deconvolution results
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with deconvolution visualization
    """
    # Find deconvolution results in obsm
    deconv_keys = [
        key
        for key in adata.obsm.keys()
        if "deconvolution" in key.lower() or "proportions" in key.lower()
    ]

    # Import pandas at the top to ensure it's always available
    import pandas as pd

    if not deconv_keys:
        # Look for individual cell type proportions in obs
        # Support all deconvolution methods
        deconv_methods = [
            "cell2location",
            "rctd",
            "destvi",
            "stereoscope",
            "spotlight",
            "tangram",
            "card",
        ]
        cell_type_cols = [
            col
            for col in adata.obs.columns
            if any(method in col.lower() for method in deconv_methods)
        ]

        if not cell_type_cols:
            raise DataNotFoundError(
                f"No deconvolution results found in adata.obsm or adata.obs.\n"
                f"Looking for columns containing: {', '.join(deconv_methods)}\n"
                f"Available columns: {list(adata.obs.columns[:10])}"
            )

        # Create proportions DataFrame from obs columns
        # Extract base key (e.g., 'deconvolution_cell2location' from 'deconvolution_cell2location_T_cell')
        base_keys = set()
        for col in cell_type_cols:
            parts = col.split("_")
            if len(parts) >= 3:
                base_key = "_".join(parts[:-1])  # Remove last part (cell type name)
                base_keys.add(base_key)

        if not base_keys:
            raise DataNotFoundError("Could not identify deconvolution result structure")

        # Use the first base key found
        base_key = list(base_keys)[0]
        method_name = base_key.split("_")[0].upper()

        # Get all columns for this base key
        relevant_cols = [col for col in cell_type_cols if col.startswith(base_key)]
        cell_types = [col.replace(f"{base_key}_", "") for col in relevant_cols]

        # Create proportions DataFrame
        proportions = pd.DataFrame(
            {
                cell_type: adata.obs[f"{base_key}_{cell_type}"].values
                for cell_type in cell_types
            },
            index=adata.obs.index,
        )
    else:
        # Use the first deconvolution key found
        deconv_key = deconv_keys[0]
        proportions = pd.DataFrame(
            adata.obsm[deconv_key],
            index=adata.obs.index,
            columns=adata.uns.get(
                f"{deconv_key}_cell_types",
                [f"CellType_{i}" for i in range(adata.obsm[deconv_key].shape[1])],
            ),
        )
        method_name = deconv_key.split("_")[0].upper()

    # Get top cell types by mean proportion
    n_cell_types = min(params.n_cell_types, proportions.shape[1])
    top_cell_types = (
        proportions.mean().sort_values(ascending=False).index[:n_cell_types]
    )

    # USE THE NEW HELPER
    fig, axes = setup_multi_panel_figure(
        n_panels=len(top_cell_types),
        params=params,
        default_title=f"{method_name} Cell Type Proportions",
    )

    # Plot each cell type
    # Use a more unique temporary column name to avoid conflicts
    temp_feature_key = "deconv_prop_temp_viz_99_unique"

    for i, cell_type in enumerate(top_cell_types):
        if i < len(axes):
            ax = axes[i]
            try:
                # Get cell type proportions
                proportions_values = proportions[cell_type].values

                # Check for NaN values
                if pd.isna(proportions_values).any():
                    proportions_values = pd.Series(proportions_values).fillna(0).values

                # Add temporary column to original adata (more efficient than copying)
                adata.obs[temp_feature_key] = proportions_values

                # Plot spatial distribution
                if "spatial" in adata.obsm:
                    # USE THE NEW SPATIAL PLOT HELPER (but we'll override the title)
                    plot_spatial_feature(
                        adata, feature=temp_feature_key, ax=ax, params=params
                    )
                    # Manually set title to show actual cell type name
                    ax.set_title(cell_type)
                    ax.invert_yaxis()
                else:
                    # Fallback: bar plot
                    sorted_props = proportions[cell_type].sort_values(ascending=False)
                    ax.bar(
                        range(len(sorted_props)),
                        sorted_props.values,
                        alpha=params.alpha,
                    )
                    ax.set_title(cell_type)
                    ax.set_xlabel("Spots (sorted)")
                    ax.set_ylabel("Proportion")

            except Exception as e:
                # Handle individual cell type plotting errors
                ax.text(
                    0.5,
                    0.5,
                    f"Error plotting {cell_type}:\n{str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{cell_type} (Error)")

    # Clean up the temporary column
    if temp_feature_key in adata.obs:
        del adata.obs[temp_feature_key]

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("Spatial Domains")
async def create_spatial_domains_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create spatial domains visualization"""
    # Look for spatial domain results in adata.obs
    domain_keys = [
        col
        for col in adata.obs.columns
        if "spatial_domains" in col.lower() or "domain" in col.lower()
    ]

    # Also check for any categorical clustering results that might represent domains
    if not domain_keys:
        # Find categorical columns (potential clustering/domain results)
        categorical_cols = [
            col
            for col in adata.obs.columns
            if adata.obs[col].dtype.name in ["object", "category"]
        ]
        domain_keys = categorical_cols[
            :5
        ]  # Take first 5 categorical columns as potential domains

    if not domain_keys:
        # No spatial domains found, suggest running domain identification first
        fig, ax = plt.subplots(figsize=params.figure_size or (10, 8))
        ax.text(
            0.5,
            0.5,
            "No spatial domains found in dataset.\n\n"
            "Please run spatial domain identification first:\n"
            'identify_spatial_domains(data_id="data_1", params={"method": "spagcn"})',
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        ax.set_title("Spatial Domains - Not Available")
        ax.axis("off")
        return fig

    # Use the first available domain key
    domain_key = domain_keys[0]
    if context:
        await context.info(f"Visualizing spatial domains using column: {domain_key}")

    # Get spatial coordinates - STRICT: NO FALLBACK to non-spatial data
    try:
        x_coords, y_coords = get_spatial_coordinates(adata)
        coord_type = "spatial"
    except Exception as e:
        # NO FALLBACK: Spatial domains require REAL spatial coordinates
        error_msg = (
            "❌ Spatial domain visualization requires actual spatial coordinates.\n\n"
            f"Error details: {str(e)}\n\n"
            "🔧 SOLUTIONS:\n"
            "1. Ensure your data contains spatial coordinates:\n"
            "   - Standard location: adata.obsm['spatial']\n"
            "   - Alternative: adata.obs['x'] and adata.obs['y']\n\n"
            "2. For Visium data, ensure proper loading:\n"
            "   sc.read_visium(path_to_data)\n\n"
            "3. For other spatial platforms, add coordinates manually:\n"
            "   adata.obsm['spatial'] = np.column_stack([x_coords, y_coords])\n\n"
            "📋 SCIENTIFIC INTEGRITY: Spatial domains are tissue regions defined by\n"
            "PHYSICAL LOCATION. PCA coordinates represent gene expression similarity,\n"
            "NOT spatial position. We cannot substitute one for the other.\n\n"
            "• Spatial coords: Real (x,y) positions on tissue\n"
            "• PCA coords: Abstract dimensions of gene expression variance\n"
            "Using PCA as spatial would completely invalidate spatial analysis."
        )
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)

    # Create figure
    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    # Get domain labels
    domains = adata.obs[domain_key].astype(str)
    unique_domains = sorted(domains.unique())
    n_domains = len(unique_domains)

    if context:
        await context.info(f"Found {n_domains} spatial domains: {unique_domains}")

    # Use a colormap with distinct colors
    if n_domains <= 10:
        colors = plt.cm.tab10(np.linspace(0, 1, n_domains))
    elif n_domains <= 20:
        colors = plt.cm.tab20(np.linspace(0, 1, n_domains))
    else:
        colors = plt.cm.viridis(np.linspace(0, 1, n_domains))

    # Create scatter plot
    for i, domain in enumerate(unique_domains):
        mask = domains == domain
        n_spots = mask.sum()
        ax.scatter(
            x_coords[mask],
            y_coords[mask],
            c=[colors[i]],
            label=f"Domain {domain} (n={n_spots})",
            s=params.spot_size or 50,
            alpha=params.alpha,
            edgecolors="none",
        )

    # Set labels and title
    if coord_type == "spatial":
        ax.set_xlabel("Spatial X")
        ax.set_ylabel("Spatial Y")
        # Invert y-axis for proper spatial orientation
        ax.invert_yaxis()
    elif coord_type == "X_spatial":
        ax.set_xlabel("Spatial X")
        ax.set_ylabel("Spatial Y")
        ax.invert_yaxis()
    else:
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    # Set title
    title = params.title or f"Spatial Domains ({domain_key})"
    ax.set_title(title, fontsize=14)

    # Add legend
    if (
        params.show_legend and n_domains <= 15
    ):  # Only show legend if not too many domains
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10)
    elif n_domains > 15:
        # Add text indicating number of domains
        ax.text(
            0.02,
            0.98,
            f"{n_domains} domains",
            transform=ax.transAxes,
            verticalalignment="top",
            fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

    ax.set_aspect("equal")

    # Adjust layout
    plt.tight_layout()

    return fig


@handle_visualization_errors("Cell Communication")
async def create_cell_communication_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create cell communication visualization"""
    if context:
        await context.info("Creating cell communication visualization")

    # Look for LIANA+ results in adata.uns and adata.obsm
    liana_keys = [key for key in adata.uns.keys() if "liana" in key.lower()]
    spatial_score_keys = [
        key for key in adata.obsm.keys() if "liana_spatial_scores" in key.lower()
    ]

    # Try to find communication results
    communication_results = None
    spatial_scores = None
    lr_pairs = None
    analysis_type = "none"

    # Check for LIANA+ spatial results
    if spatial_score_keys or "liana_spatial_scores" in adata.obsm:
        spatial_scores_key = "liana_spatial_scores"
        if spatial_scores_key in adata.obsm:
            spatial_scores = adata.obsm[spatial_scores_key]
            # Get LR pair names from uns
            if "liana_spatial_interactions" in adata.uns:
                lr_pairs = adata.uns["liana_spatial_interactions"]
            analysis_type = "spatial"
            if context:
                await context.info(
                    f"Found LIANA+ spatial scores in obsm: {spatial_scores_key}"
                )
                if lr_pairs:
                    await context.info(f"Found {len(lr_pairs)} LR pairs")

    # Check for general LIANA+ results
    elif liana_keys:
        liana_key = liana_keys[0]
        if liana_key in adata.uns:
            communication_results = adata.uns[liana_key]
            analysis_type = "cluster"
            if context:
                await context.info(f"Found LIANA+ results: {liana_key}")

    # Check for ligand-receptor pair results in obs
    lr_columns = [
        col
        for col in adata.obs.columns
        if any(
            pattern in col.lower()
            for pattern in ["ligand", "receptor", "lr_", "communication"]
        )
    ]

    if analysis_type == "spatial" and spatial_scores is not None:
        return _create_spatial_communication_plot(
            adata, spatial_scores, lr_pairs, params, context
        )
    elif analysis_type == "cluster" and communication_results is not None:
        return _create_cluster_communication_plot(
            adata, communication_results, params, context
        )
    elif lr_columns:
        return _create_lr_expression_plot(adata, lr_columns, params, context)
    else:
        # No communication data found, create instructional plot
        fig, ax = plt.subplots(figsize=params.figure_size or (10, 8))
        ax.text(
            0.5,
            0.5,
            "No cell communication results found in dataset.\n\n"
            "To analyze cell communication, first run:\n"
            'analyze_cell_communication(data_id="data_1", params={"method": "liana"})\n\n'
            "This will generate LIANA+ results for visualization.",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        ax.set_title("Cell Communication - Not Available")
        ax.axis("off")
        return fig


def _create_spatial_communication_plot(
    adata, spatial_scores, lr_pairs, params, context
):
    """Create spatial communication visualization using LIANA+ spatial scores"""
    try:
        # Get spatial coordinates
        x_coords, y_coords = get_spatial_coordinates(adata)

        # Get LR results with significance info
        liana_res = None
        global_metric_col = None
        if "liana_spatial_res" in adata.uns:
            liana_res = adata.uns["liana_spatial_res"]

        # Get top communication pairs
        if lr_pairs and liana_res is not None:
            # Get top pairs based on global metric (morans or lee)
            import pandas as pd

            if isinstance(liana_res, pd.DataFrame):
                # Find global metric column
                global_metric_col = None
                for col in ["morans", "lee", "global_score"]:
                    if col in liana_res.columns:
                        global_metric_col = col
                        break

                if global_metric_col:
                    # Sort by global metric and get top pairs
                    n_top = getattr(params, "plot_top_pairs", 6)
                    top_results = liana_res.nlargest(n_top, global_metric_col)
                    top_pairs = top_results.index.tolist()
                    # Get the indices of these pairs in the original list
                    pair_indices = [
                        lr_pairs.index(pair) if pair in lr_pairs else -1
                        for pair in top_pairs
                    ]
                    pair_indices = [idx for idx in pair_indices if idx >= 0]
                else:
                    # Use first N pairs
                    n_top = getattr(params, "plot_top_pairs", 6)
                    top_pairs = lr_pairs[:n_top]
                    pair_indices = list(range(len(top_pairs)))
            else:
                # Use first N pairs
                n_top = getattr(params, "plot_top_pairs", 6)
                top_pairs = lr_pairs[:n_top]
                pair_indices = list(range(len(top_pairs)))
        elif lr_pairs:
            # Use provided LR pairs
            n_top = getattr(params, "plot_top_pairs", 6)
            top_pairs = lr_pairs[:n_top]
            pair_indices = list(range(len(top_pairs)))
        else:
            top_pairs = []
            pair_indices = []

        if not top_pairs:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(
                0.5,
                0.5,
                "No communication pairs found in spatial scores",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title("Cell Communication - Spatial")
            ax.axis("off")
            return fig

        # Create subplot layout
        n_pairs = min(len(top_pairs), 6)  # Limit to 6 for display
        n_cols = min(3, n_pairs)
        n_rows = (n_pairs + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
        if n_pairs == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()

        # Ensure pair_indices matches top_pairs
        if len(pair_indices) < len(top_pairs):
            pair_indices.extend(list(range(len(pair_indices), len(top_pairs))))

        for i, (pair, pair_idx) in enumerate(
            zip(top_pairs[:n_pairs], pair_indices[:n_pairs])
        ):
            ax = axes[i]

            try:
                # Get spatial scores for this pair
                # spatial_scores is a numpy array where columns are LR pairs
                if pair_idx < spatial_scores.shape[1]:
                    scores = spatial_scores[:, pair_idx]
                else:
                    scores = np.zeros(len(adata))

                # Create scatter plot
                scatter = ax.scatter(
                    x_coords,
                    y_coords,
                    c=scores,
                    cmap=params.colormap or "viridis",
                    s=15,
                    alpha=0.8,
                )

                # Format pair name for display
                # LR pairs are typically formatted as "Ligand^Receptor"
                if "^" in pair:
                    display_name = pair.replace("^", " → ")
                elif "_" in pair:
                    display_name = pair.replace("_", " → ")
                else:
                    display_name = pair

                # Add global metric value if available
                if (
                    liana_res is not None
                    and isinstance(liana_res, pd.DataFrame)
                    and pair in liana_res.index
                ):
                    if global_metric_col and global_metric_col in liana_res.columns:
                        metric_value = liana_res.loc[pair, global_metric_col]
                        display_name += f"\n({global_metric_col}: {metric_value:.3f})"

                ax.set_title(display_name, fontsize=10)
                ax.set_xlabel("X coordinate")
                ax.set_ylabel("Y coordinate")

                # Add colorbar
                plt.colorbar(scatter, ax=ax, shrink=0.7, label="Communication Score")

            except Exception as e:
                ax.text(
                    0.5,
                    0.5,
                    f"Error: {str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(pair, fontsize=10)

        # Hide unused subplots
        for i in range(n_pairs, len(axes)):
            axes[i].set_visible(False)

        plt.suptitle("Cell Communication - Spatial Distribution", fontsize=14)
        plt.tight_layout()
        return fig

    except Exception as e:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(
            0.5,
            0.5,
            f"Error creating spatial communication plot:\n{str(e)}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title("Cell Communication - Error")
        ax.axis("off")
        return fig


def _create_cluster_communication_plot(adata, communication_results, params, context):
    """Create cluster-based communication visualization"""
    try:
        # Try to extract top communication pairs from results
        if isinstance(communication_results, pd.DataFrame):
            # Results are in DataFrame format
            df = communication_results

            # Look for common LIANA+ columns
            score_cols = [
                col
                for col in df.columns
                if any(
                    term in col.lower() for term in ["score", "pvalue", "significant"]
                )
            ]

            if score_cols and len(df) > 0:
                # Sort by first score column and get top pairs
                top_results = df.nlargest(8, score_cols[0])

                # Create ligand-receptor pair names
                if "ligand" in df.columns and "receptor" in df.columns:
                    pair_names = [
                        f"{row['ligand']} → {row['receptor']}"
                        for _, row in top_results.iterrows()
                    ]
                    scores = top_results[score_cols[0]].values
                elif "lr_pair" in df.columns:
                    pair_names = top_results["lr_pair"].tolist()
                    scores = top_results[score_cols[0]].values
                else:
                    pair_names = [f"Pair {i+1}" for i in range(len(top_results))]
                    scores = top_results[score_cols[0]].values
            else:
                pair_names = ["No significant pairs"]
                scores = [0]
        else:
            # Results in other format
            pair_names = ["Communication analysis completed"]
            scores = [1.0]

        # Create horizontal bar plot
        fig, ax = plt.subplots(figsize=params.figure_size or (12, 8))

        y_pos = np.arange(len(pair_names))
        bars = ax.barh(y_pos, scores, color="steelblue", alpha=0.7)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(pair_names, fontsize=10)
        ax.set_xlabel("Communication Strength")
        ax.set_title("Top Cell Communication Pairs")
        ax.invert_yaxis()

        # Add value labels on bars
        for i, (bar, score) in enumerate(zip(bars, scores)):
            if score > 0:
                ax.text(
                    score + max(scores) * 0.01,
                    bar.get_y() + bar.get_height() / 2,
                    f"{score:.3f}",
                    va="center",
                    fontsize=9,
                )

        plt.tight_layout()
        return fig

    except Exception as e:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(
            0.5,
            0.5,
            f"Error creating cluster communication plot:\n{str(e)}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title("Cell Communication - Error")
        ax.axis("off")
        return fig


def _create_lr_expression_plot(adata, lr_columns, params, context):
    """Create ligand-receptor expression visualization"""
    try:
        # Get spatial coordinates
        x_coords, y_coords = get_spatial_coordinates(adata)

        # Select top LR columns (limit to 6)
        selected_cols = lr_columns[:6]

        # Create subplot layout
        n_cols = min(3, len(selected_cols))
        n_rows = (len(selected_cols) + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
        if len(selected_cols) == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()

        for i, col in enumerate(selected_cols):
            ax = axes[i]

            try:
                values = adata.obs[col].values

                # Handle categorical or numeric data
                if pd.api.types.is_numeric_dtype(values):
                    scatter = ax.scatter(
                        x_coords,
                        y_coords,
                        c=values,
                        cmap=params.colormap or "viridis",
                        s=15,
                        alpha=0.8,
                    )
                    plt.colorbar(scatter, ax=ax, shrink=0.7)
                else:
                    # Categorical data
                    unique_vals = pd.unique(values)
                    colors = plt.cm.Set1(np.linspace(0, 1, len(unique_vals)))

                    for j, val in enumerate(unique_vals):
                        mask = values == val
                        ax.scatter(
                            x_coords[mask],
                            y_coords[mask],
                            c=[colors[j]],
                            label=str(val),
                            s=15,
                            alpha=0.8,
                        )

                    if len(unique_vals) <= 10:
                        ax.legend(
                            bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8
                        )

                # Clean up column name for display
                display_name = col.replace("_", " ").title()
                ax.set_title(display_name, fontsize=10)
                ax.set_xlabel("X coordinate")
                ax.set_ylabel("Y coordinate")

            except Exception as e:
                ax.text(
                    0.5,
                    0.5,
                    f"Error: {str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(col, fontsize=10)

        # Hide unused subplots
        for i in range(len(selected_cols), len(axes)):
            axes[i].set_visible(False)

        plt.suptitle("Ligand-Receptor Expression Patterns", fontsize=14)
        plt.tight_layout()
        return fig

    except Exception as e:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(
            0.5,
            0.5,
            f"Error creating LR expression plot:\n{str(e)}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_title("Cell Communication - Error")
        ax.axis("off")
        return fig


@handle_visualization_errors("Multi-Gene UMAP")
async def create_multi_gene_umap_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create multi-gene UMAP visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with multi-gene UMAP visualization
    """
    # USE THE NEW FEATURE HELPER
    available_genes = await get_validated_features(
        adata, params, max_features=12, context=context
    )

    if context:
        await context.info(
            f"Creating multi-panel UMAP plot for features: {available_genes}"
        )

    # USE THE NEW PANEL HELPER
    fig, axes = setup_multi_panel_figure(
        n_panels=len(available_genes),
        params=params,
        default_title=f"Multi-Gene Expression UMAP ({len(available_genes)} genes)",
    )

    # Plot each gene
    # Use a more unique temporary column name to avoid conflicts
    temp_feature_key = "umap_gene_temp_viz_99_unique"

    for i, gene in enumerate(available_genes):
        if i < len(axes):
            ax = axes[i]
            try:
                # Get gene expression to determine color scale
                gene_expr = adata[:, gene].X
                if hasattr(gene_expr, "toarray"):
                    gene_expr = gene_expr.toarray().flatten()
                elif hasattr(gene_expr, "A1"):
                    gene_expr = gene_expr.A1
                else:
                    gene_expr = np.array(gene_expr).flatten()

                # Apply color scaling if specified
                if params.color_scale == "log":
                    gene_expr = np.log1p(gene_expr)
                elif params.color_scale == "sqrt":
                    gene_expr = np.sqrt(gene_expr)

                # Add temporary column
                adata.obs[temp_feature_key] = gene_expr

                # Set appropriate color scale for sparse data
                vmin = 0  # Always start from 0 for gene expression
                vmax = max(gene_expr.max(), 0.1)  # Ensure we have some range

                # Use percentile-based scaling for better visualization
                if np.sum(gene_expr > 0) > 10:  # If we have enough expressing cells
                    vmax = np.percentile(gene_expr[gene_expr > 0], 95)

                sc.pl.umap(
                    adata,
                    color=temp_feature_key,
                    cmap=params.colormap,
                    ax=ax,
                    show=False,
                    frameon=False,
                    vmin=vmin,
                    vmax=vmax,
                    colorbar_loc="right",
                )
                # Manually set title to show actual gene name
                ax.set_title(gene)

            except Exception as e:
                # Handle individual gene plotting errors
                ax.text(
                    0.5,
                    0.5,
                    f"Error plotting {gene}:\n{str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{gene} (Error)")

    # Clean up the temporary column
    if temp_feature_key in adata.obs:
        del adata.obs[temp_feature_key]

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("Multi-Gene")
async def create_multi_gene_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create multi-gene visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with multi-gene visualization
    """
    # USE THE NEW FEATURE HELPER
    available_genes = await get_validated_features(
        adata, params, max_features=12, context=context
    )

    if context:
        await context.info(
            f"Visualizing {len(available_genes)} genes: {available_genes}"
        )

    # USE THE NEW PANEL HELPER
    # Use tight_layout for better handling of colorbars in spatial plots
    fig, axes = setup_multi_panel_figure(
        n_panels=len(available_genes),
        params=params,
        default_title=f"Multi-Gene Expression ({len(available_genes)} genes)",
        use_tight_layout=True,
    )

    # Plot each gene
    # Use a more unique temporary column name to avoid conflicts
    temp_feature_key = "multi_gene_expr_temp_viz_99_unique"

    for i, gene in enumerate(available_genes):
        if i < len(axes):
            ax = axes[i]
            try:
                # Get gene expression
                gene_expr = adata[:, gene].X
                if hasattr(gene_expr, "toarray"):
                    gene_expr = gene_expr.toarray().flatten()
                elif hasattr(gene_expr, "A1"):
                    gene_expr = gene_expr.A1
                else:
                    gene_expr = np.array(gene_expr).flatten()

                # Apply color scaling
                if params.color_scale == "log":
                    gene_expr = np.log1p(gene_expr)
                elif params.color_scale == "sqrt":
                    gene_expr = np.sqrt(gene_expr)

                # Add temporary column to original adata (more efficient than copying)
                adata.obs[temp_feature_key] = gene_expr

                # Create spatial plot
                if "spatial" in adata.obsm:
                    # USE THE NEW SPATIAL PLOT HELPER (disable colorbar for manual control)
                    plot_spatial_feature(
                        adata,
                        feature=temp_feature_key,
                        ax=ax,
                        params=params,
                        force_no_colorbar=True,
                    )

                    # Set color limits manually for better visualization
                    vmin = (
                        params.vmin
                        if params.vmin is not None
                        else np.percentile(gene_expr, 1)
                    )
                    vmax = (
                        params.vmax
                        if params.vmax is not None
                        else np.percentile(gene_expr, 99)
                    )

                    # Update the colorbar limits
                    scatter = ax.collections[0] if ax.collections else None
                    if scatter:
                        scatter.set_clim(vmin, vmax)

                        # Add colorbar with exact height matching using make_axes_locatable
                        if params.show_colorbar:
                            divider = make_axes_locatable(ax)
                            cax = divider.append_axes("right", size="5%", pad=0.05)
                            plt.colorbar(scatter, cax=cax)

                    ax.invert_yaxis()
                    # Manually set title to show actual gene name
                    if params.add_gene_labels:
                        ax.set_title(gene, fontsize=12)
                else:
                    # Fallback: histogram
                    ax.hist(gene_expr, bins=30, alpha=params.alpha, color="steelblue")
                    ax.set_xlabel("Expression")
                    ax.set_ylabel("Frequency")
                    if params.add_gene_labels:
                        ax.set_title(gene, fontsize=12)

            except Exception as e:
                # Handle individual gene plotting errors
                ax.text(
                    0.5,
                    0.5,
                    f"Error plotting {gene}:\n{str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                if params.add_gene_labels:
                    ax.set_title(f"{gene} (Error)", fontsize=12)

    # Clean up the temporary column
    if temp_feature_key in adata.obs:
        del adata.obs[temp_feature_key]

    # Adjust spacing for compact layout with proper colorbar spacing
    # Using subplots_adjust directly (no tight_layout to avoid conflicts)
    fig.subplots_adjust(top=0.92, wspace=0.2, hspace=0.3, right=0.98)
    return fig


@handle_visualization_errors("Ligand-Receptor Pairs")
async def create_lr_pairs_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create ligand-receptor pairs visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with LR pairs visualization
    """
    # Ensure unique gene names to avoid indexing errors
    if not adata.var_names.is_unique:
        if context:
            await context.info("Making gene names unique to avoid indexing errors")
        adata.var_names_make_unique()

    # Parse LR pairs from various sources
    lr_pairs = []

    # 1. Check for explicit lr_pairs parameter
    if params.lr_pairs:
        lr_pairs = params.lr_pairs
    else:
        # 2. Try to parse from feature parameter
        feature_list = (
            params.feature
            if isinstance(params.feature, list)
            else ([params.feature] if params.feature else [])
        )

        if feature_list:
            # First check if any items have special format
            has_special_format = any(
                "^" in str(f) or ("_" in str(f) and not str(f).startswith("_"))
                for f in feature_list
            )

            if has_special_format:
                # Parse different L-R pair formats
                for item in feature_list:
                    # Handle "Ligand^Receptor" format from LIANA
                    if "^" in str(item):
                        ligand, receptor = str(item).split("^", 1)
                        lr_pairs.append((ligand, receptor))
                    # Handle "Ligand_Receptor" format
                    elif "_" in str(item) and not str(item).startswith("_"):
                        parts = str(item).split("_")
                        if len(parts) == 2:
                            lr_pairs.append((parts[0], parts[1]))

        # 3. Try to get from stored analysis results if no pairs found yet
        if not lr_pairs and hasattr(adata, "uns"):
            # Check for standardized storage
            if "detected_lr_pairs" in adata.uns:
                lr_pairs = adata.uns["detected_lr_pairs"]
            # Try to parse from cell communication results
            elif "cell_communication_results" in adata.uns:
                comm_results = adata.uns["cell_communication_results"]
                if "top_lr_pairs" in comm_results:
                    for pair_str in comm_results["top_lr_pairs"]:
                        if "^" in pair_str:
                            ligand, receptor = pair_str.split("^", 1)
                            lr_pairs.append((ligand, receptor))
                        elif "_" in pair_str:
                            parts = pair_str.split("_")
                            if len(parts) == 2:
                                lr_pairs.append((parts[0], parts[1]))

    # CRITICAL: No hardcoded defaults! Scientific integrity requires real data only
    if not lr_pairs:
        raise DataNotFoundError(
            "No ligand-receptor pairs to visualize.\n\n"
            "This tool requires actual L-R pairs from your analysis. Options:\n"
            "1. Run cell communication analysis first (recommended)\n"
            "2. Specify lr_pairs parameter: lr_pairs=[('Ligand', 'Receptor')]\n"
            "3. Use LIANA format: feature=['Ligand^Receptor']\n"
            "4. Use underscore format: feature=['Ligand_Receptor']\n\n"
            "Note: This tool does NOT provide demo/default data.\n"
            "Only real experimental results should be visualized in scientific analysis."
        )

    # Filter pairs where both genes exist
    available_pairs = []
    for ligand, receptor in lr_pairs:
        if ligand in adata.var_names and receptor in adata.var_names:
            available_pairs.append((ligand, receptor))

    if not available_pairs:
        raise DataNotFoundError(
            f"None of the specified LR pairs found in data: {lr_pairs}"
        )

    # Limit to avoid overly large plots
    max_pairs = 4
    if len(available_pairs) > max_pairs:
        if context:
            await context.warning(
                f"Too many LR pairs ({len(available_pairs)}). Limiting to first {max_pairs}."
            )
        available_pairs = available_pairs[:max_pairs]

    if context:
        await context.info(
            f"Visualizing {len(available_pairs)} LR pairs: {available_pairs}"
        )

    # Each pair gets 3 panels: ligand, receptor, correlation
    n_panels = len(available_pairs) * 3

    # USE THE NEW PANEL HELPER
    # Use tight_layout for better handling of colorbars in spatial plots
    fig, axes = setup_multi_panel_figure(
        n_panels=n_panels,
        params=params,
        default_title=f"Ligand-Receptor Pairs ({len(available_pairs)} pairs)",
        use_tight_layout=True,
    )

    # Plot each LR pair
    ax_idx = 0
    # Use a more unique temporary column name to avoid conflicts
    temp_feature_key = "lr_expr_temp_viz_99_unique"

    for pair_idx, (ligand, receptor) in enumerate(available_pairs):
        try:
            # Get expression data
            ligand_expr = adata[:, ligand].X
            receptor_expr = adata[:, receptor].X

            if hasattr(ligand_expr, "toarray"):
                ligand_expr = ligand_expr.toarray().flatten()
                receptor_expr = receptor_expr.toarray().flatten()
            elif hasattr(ligand_expr, "A1"):
                ligand_expr = ligand_expr.A1
                receptor_expr = receptor_expr.A1
            else:
                ligand_expr = np.array(ligand_expr).flatten()
                receptor_expr = np.array(receptor_expr).flatten()

            # Apply color scaling
            if params.color_scale == "log":
                ligand_expr = np.log1p(ligand_expr)
                receptor_expr = np.log1p(receptor_expr)
            elif params.color_scale == "sqrt":
                ligand_expr = np.sqrt(ligand_expr)
                receptor_expr = np.sqrt(receptor_expr)

            # Plot ligand
            if ax_idx < len(axes) and "spatial" in adata.obsm:
                ax = axes[ax_idx]
                # Add temporary column for ligand expression
                adata.obs[temp_feature_key] = ligand_expr
                # USE THE NEW SPATIAL PLOT HELPER (disable colorbar for manual control)
                plot_spatial_feature(
                    adata,
                    feature=temp_feature_key,
                    ax=ax,
                    params=params,
                    force_no_colorbar=True,
                )

                # Add colorbar with exact height matching using make_axes_locatable
                if params.show_colorbar and ax.collections:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    plt.colorbar(ax.collections[-1], cax=cax)

                ax.invert_yaxis()
                # Manually set the title to show actual ligand name
                if params.add_gene_labels:
                    ax.set_title(f"{ligand} (Ligand)", fontsize=10)
                ax_idx += 1

            # Plot receptor
            if ax_idx < len(axes) and "spatial" in adata.obsm:
                ax = axes[ax_idx]
                # Add temporary column for receptor expression
                adata.obs[temp_feature_key] = receptor_expr
                # USE THE NEW SPATIAL PLOT HELPER (disable colorbar for manual control)
                plot_spatial_feature(
                    adata,
                    feature=temp_feature_key,
                    ax=ax,
                    params=params,
                    force_no_colorbar=True,
                )

                # Add colorbar with exact height matching using make_axes_locatable
                if params.show_colorbar and ax.collections:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    plt.colorbar(ax.collections[-1], cax=cax)

                ax.invert_yaxis()
                # Manually set the title to show actual receptor name
                if params.add_gene_labels:
                    ax.set_title(f"{receptor} (Receptor)", fontsize=10)
                ax_idx += 1

            # Plot correlation
            if ax_idx < len(axes):
                ax = axes[ax_idx]

                # Calculate correlation
                from scipy.stats import kendalltau, pearsonr, spearmanr

                if params.correlation_method == "pearson":
                    corr, p_value = pearsonr(ligand_expr, receptor_expr)
                elif params.correlation_method == "spearman":
                    corr, p_value = spearmanr(ligand_expr, receptor_expr)
                else:  # kendall
                    corr, p_value = kendalltau(ligand_expr, receptor_expr)

                # Create scatter plot
                ax.scatter(ligand_expr, receptor_expr, alpha=params.alpha, s=20)
                ax.set_xlabel(f"{ligand} Expression")
                ax.set_ylabel(f"{receptor} Expression")

                if params.show_correlation_stats:
                    ax.set_title(
                        f"Correlation: {corr:.3f}\np-value: {p_value:.2e}", fontsize=10
                    )
                else:
                    ax.set_title(f"{ligand} vs {receptor}", fontsize=10)

                # Add trend line
                z = np.polyfit(ligand_expr, receptor_expr, 1)
                p = np.poly1d(z)
                ax.plot(ligand_expr, p(ligand_expr), "r--", alpha=0.8)

                ax_idx += 1

        except Exception as e:
            # Handle individual pair plotting errors
            if ax_idx < len(axes):
                ax = axes[ax_idx]
                ax.text(
                    0.5,
                    0.5,
                    f"Error plotting {ligand}-{receptor}:\n{str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{ligand}-{receptor} (Error)", fontsize=10)
                ax_idx += 1

    # Clean up the temporary column
    if temp_feature_key in adata.obs:
        del adata.obs[temp_feature_key]

    # Adjust spacing for compact layout with proper colorbar spacing
    # Using subplots_adjust directly (no tight_layout to avoid conflicts)
    fig.subplots_adjust(top=0.92, wspace=0.2, hspace=0.3, right=0.98)
    return fig


@handle_visualization_errors("RNA Velocity")
async def create_rna_velocity_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create RNA velocity stream plot

    Args:
        adata: AnnData object with computed RNA velocity
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with RNA velocity stream plot
    """
    try:
        import scvelo as scv
    except ImportError:
        raise DataCompatibilityError(
            "scvelo package is required for RNA velocity visualization. Install it with: pip install scvelo>=0.2.5"
        )

    # Check if RNA velocity has been computed
    if "velocity_graph" not in adata.uns:
        raise DataNotFoundError(
            "RNA velocity has not been computed. Please run 'analyze_rna_velocity' first."
        )

    if context:
        await context.info("Creating RNA velocity stream plot")

    # Determine basis for plotting
    basis = params.basis or "spatial"
    basis_key = f"X_{basis}" if basis != "spatial" else "spatial"

    if basis_key not in adata.obsm:
        # Try to find an alternative basis
        if "spatial" in adata.obsm:
            basis = "spatial"
            if context:
                await context.info(
                    "Using 'spatial' as basis for RNA velocity visualization"
                )
        elif "X_umap" in adata.obsm:
            basis = "umap"
            if context:
                await context.info(
                    "Using 'umap' as basis for RNA velocity visualization"
                )
        elif "X_pca" in adata.obsm:
            basis = "pca"
            if context:
                await context.info(
                    "Using 'pca' as basis for RNA velocity visualization"
                )
        else:
            available_bases = [
                k.replace("X_", "") for k in adata.obsm.keys() if k.startswith("X_")
            ]
            available_bases.extend([k for k in adata.obsm.keys() if k == "spatial"])
            raise DataCompatibilityError(
                f"Basis '{params.basis or 'spatial'}' not found. Available bases: {available_bases}"
            )

    # Prepare feature for coloring
    # Find first categorical column as default (no hardcoded assumptions)
    default_feature = None
    if not params.feature:
        categorical_cols = [
            col
            for col in adata.obs.columns
            if adata.obs[col].dtype.name in ["object", "category"]
        ]
        default_feature = categorical_cols[0] if categorical_cols else None

    feature = await validate_and_prepare_feature(
        adata,
        params.feature,
        context,
        default_feature=default_feature,
    )

    # Create figure
    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    # Create RNA velocity stream plot
    scv.pl.velocity_embedding_stream(
        adata,
        basis=basis,
        color=feature,
        ax=ax,
        show=False,
        alpha=params.alpha,
        legend_loc="right margin" if feature and feature in adata.obs.columns else None,
        frameon=params.show_axes,
        title="",  # We'll set our own title
    )

    # Set title
    title = params.title or f"RNA Velocity Stream on {basis.capitalize()}"
    ax.set_title(title, fontsize=14)

    # Handle spatial coordinates orientation
    if basis == "spatial":
        ax.invert_yaxis()

    plt.tight_layout()
    return fig


@handle_visualization_errors("Spatial Statistics")
async def create_spatial_statistics_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create spatial statistics visualization based on analysis_type

    Args:
        adata: AnnData object with spatial statistics results
        params: Visualization parameters including analysis_type
        context: MCP context

    Returns:
        Matplotlib figure with spatial statistics visualization
    """
    if context:
        await context.info(
            f"Creating {params.analysis_type} spatial statistics visualization"
        )

    if params.analysis_type == "neighborhood":
        return await create_neighborhood_enrichment_visualization(
            adata, params, context
        )
    elif params.analysis_type == "co_occurrence":
        return await create_co_occurrence_visualization(adata, params, context)
    elif params.analysis_type == "ripley":
        return await create_ripley_visualization(adata, params, context)
    elif params.analysis_type == "moran":
        return await create_moran_visualization(adata, params, context)
    elif params.analysis_type == "centrality":
        return await create_centrality_visualization(adata, params, context)
    elif params.analysis_type == "getis_ord":
        return await create_getis_ord_visualization(adata, params, context)
    else:
        raise InvalidParameterError(
            f"Unsupported analysis type: {params.analysis_type}"
        )


@handle_visualization_errors("Neighborhood Enrichment")
async def create_neighborhood_enrichment_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create neighborhood enrichment visualization with optional network view"""
    # Infer cluster_key from params or from available results in adata.uns
    cluster_key = params.cluster_key
    if not cluster_key:
        # Try to infer from adata.uns keys
        enrichment_keys = [
            k for k in adata.uns.keys() if k.endswith("_nhood_enrichment")
        ]
        if enrichment_keys:
            cluster_key = enrichment_keys[0].replace("_nhood_enrichment", "")
            if context:
                await context.info(
                    f"Inferred cluster_key: '{cluster_key}' from existing results"
                )
        else:
            raise ValueError(
                "cluster_key parameter is required but not provided.\n\n"
                "Available categorical columns in data:\n  "
                + ", ".join(
                    [
                        col
                        for col in adata.obs.columns
                        if adata.obs[col].dtype.name in ["object", "category"]
                    ][:10]
                )
                + "\n\nPlease specify: params={'cluster_key': 'your_column_name'}"
            )

    enrichment_key = f"{cluster_key}_nhood_enrichment"
    if enrichment_key not in adata.uns:
        raise DataNotFoundError(
            f"Neighborhood enrichment results not found. Expected key: {enrichment_key}"
        )

    enrichment_matrix = adata.uns[enrichment_key]["zscore"]
    categories = adata.obs[cluster_key].cat.categories

    # Check if network visualization is requested
    if params.show_network:
        if context:
            await context.info(
                "Creating network-style neighborhood enrichment visualization"
            )
        return await create_neighborhood_network_visualization(
            enrichment_matrix, categories, params, context
        )
    else:
        # Standard heatmap visualization
        figsize = params.figure_size or (10, 8)
        fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

        im = ax.imshow(enrichment_matrix, cmap=params.colormap)

        if params.show_colorbar:
            plt.colorbar(im, ax=ax, label="Z-score")

        ax.set_xticks(np.arange(len(categories)))
        ax.set_yticks(np.arange(len(categories)))
        ax.set_xticklabels(categories, rotation=45, ha="right")
        ax.set_yticklabels(categories)

        title = params.title or f"Neighborhood Enrichment ({cluster_key})"
        ax.set_title(title)
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Cluster")

        plt.tight_layout()
        return fig


async def create_neighborhood_network_visualization(
    enrichment_matrix: np.ndarray,
    categories,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """Create network-style visualization of neighborhood enrichment"""
    try:
        import networkx as nx
    except ImportError:
        if context:
            await context.warning(
                "NetworkX not available. Falling back to heatmap visualization."
            )
        # Fallback to standard heatmap
        figsize = params.figure_size or (10, 8)
        fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)
        ax.imshow(enrichment_matrix, cmap=params.colormap)
        ax.set_title("Neighborhood Enrichment (NetworkX not available)")
        return fig

    # Create network graph
    G = nx.Graph()

    # Add nodes (cell types/clusters)
    for i, category in enumerate(categories):
        G.add_node(category)

    # Add edges based on enrichment scores above threshold
    threshold = params.network_threshold
    for i in range(len(categories)):
        for j in range(i + 1, len(categories)):
            zscore = enrichment_matrix[i, j]
            if abs(zscore) > threshold:
                G.add_edge(
                    categories[i], categories[j], weight=abs(zscore), zscore=zscore
                )

    # Create visualization
    figsize = params.figure_size or (12, 10)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    # Choose layout
    if params.network_layout == "spring":
        pos = nx.spring_layout(G, k=2, iterations=50)
    elif params.network_layout == "circular":
        pos = nx.circular_layout(G)
    elif params.network_layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)
    elif params.network_layout == "fruchterman_reingold":
        pos = nx.fruchterman_reingold_layout(G)
    else:
        pos = nx.spring_layout(G)

    # Draw nodes
    node_sizes = [
        500 + 100 * len(categories) for _ in G.nodes()
    ]  # Size based on number of categories
    nx.draw_networkx_nodes(
        G, pos, node_color="lightblue", node_size=node_sizes, alpha=0.8, ax=ax
    )

    # Draw edges with colors based on enrichment/depletion
    edges = G.edges()
    edge_colors = []
    edge_widths = []

    for edge in edges:
        zscore = G[edge[0]][edge[1]]["zscore"]
        weight = G[edge[0]][edge[1]]["weight"]

        # Color: red for enrichment (positive), blue for depletion (negative)
        edge_colors.append("red" if zscore > 0 else "blue")
        # Width based on absolute z-score
        edge_widths.append(min(5, max(0.5, weight / 2)))

    if edges:
        nx.draw_networkx_edges(
            G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.6, ax=ax
        )

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight="bold", ax=ax)

    # Set title and clean up
    title = params.title or f"Neighborhood Enrichment Network (threshold: {threshold})"
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")

    # Add legend
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], color="red", lw=2, label="Enrichment (Z > 0)"),
        Line2D([0], [0], color="blue", lw=2, label="Depletion (Z < 0)"),
        Line2D([0], [0], color="gray", lw=1, label=f"|Z-score| > {threshold}"),
    ]
    ax.legend(handles=legend_elements, loc="upper right")

    # Add network statistics as text
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    stats_text = f"Nodes: {n_nodes}\nEdges: {n_edges}\nThreshold: {threshold}"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    return fig


@handle_visualization_errors("Co-occurrence")
async def create_co_occurrence_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create co-occurrence analysis visualization"""
    # Infer cluster_key from params or from available results in adata.uns
    cluster_key = params.cluster_key
    if not cluster_key:
        # Try to infer from adata.uns keys
        co_occurrence_keys = [
            k for k in adata.uns.keys() if k.endswith("_co_occurrence")
        ]
        if co_occurrence_keys:
            cluster_key = co_occurrence_keys[0].replace("_co_occurrence", "")
            if context:
                await context.info(
                    f"Inferred cluster_key: '{cluster_key}' from existing results"
                )
        else:
            raise ValueError(
                "cluster_key parameter is required but not provided.\n\n"
                "Available categorical columns in data:\n  "
                + ", ".join(
                    [
                        col
                        for col in adata.obs.columns
                        if adata.obs[col].dtype.name in ["object", "category"]
                    ][:10]
                )
                + "\n\nPlease specify: params={'cluster_key': 'your_column_name'}"
            )

    co_occurrence_key = f"{cluster_key}_co_occurrence"
    if co_occurrence_key not in adata.uns:
        raise DataNotFoundError(
            f"Co-occurrence results not found. Expected key: {co_occurrence_key}"
        )

    co_occurrence_matrix = adata.uns[co_occurrence_key]["occ"]
    interval = adata.uns[co_occurrence_key]["interval"]
    categories = adata.obs[cluster_key].cat.categories

    # The co-occurrence matrix is 3D (clusters x clusters x distances)
    n_clusters = co_occurrence_matrix.shape[0]
    n_distances = co_occurrence_matrix.shape[2]

    # Limit to at most 4 distance panels to keep the figure manageable
    max_panels = min(4, n_distances)
    selected_distances = np.linspace(0, n_distances - 1, max_panels, dtype=int)

    fig, axes = setup_multi_panel_figure(
        n_panels=max_panels,
        params=params,
        default_title=f"Co-occurrence at different distances ({cluster_key})",
    )

    for i, dist_idx in enumerate(selected_distances):
        if i < len(axes):
            ax = axes[i]
            dist_value = interval[dist_idx]

            im = ax.imshow(
                co_occurrence_matrix[:, :, dist_idx],
                cmap=params.colormap,
                vmin=0,
                vmax=np.max(co_occurrence_matrix),
            )

            # Set ticks and labels
            if n_clusters <= 10:
                ax.set_xticks(np.arange(len(categories)))
                ax.set_yticks(np.arange(len(categories)))
                ax.set_xticklabels(categories, rotation=45, ha="right")
                ax.set_yticklabels(categories)
            else:
                step = max(1, n_clusters // 5)
                ax.set_xticks(np.arange(0, n_clusters, step))
                ax.set_yticks(np.arange(0, n_clusters, step))
                ax.set_xticklabels(categories[::step], rotation=45, ha="right")
                ax.set_yticklabels(categories[::step])

            ax.set_title(f"Distance: {dist_value:.2f}")

            if i == 0:
                ax.set_ylabel("Cluster")
            if i == max_panels // 2:
                ax.set_xlabel("Cluster")

    # Add colorbar to the last subplot
    if params.show_colorbar:
        plt.colorbar(im, ax=axes[-1], label="Co-occurrence probability")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("Ripley Function")
async def create_ripley_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create Ripley's function visualization"""
    # Infer cluster_key from params or from available results in adata.uns
    cluster_key = params.cluster_key
    if not cluster_key:
        # Try to infer from adata.uns keys
        ripley_keys = [k for k in adata.uns.keys() if k.endswith("_ripley_L")]
        if ripley_keys:
            cluster_key = ripley_keys[0].replace("_ripley_L", "")
            if context:
                await context.info(
                    f"Inferred cluster_key: '{cluster_key}' from existing results"
                )
        else:
            raise ValueError(
                "cluster_key parameter is required but not provided.\n\n"
                "Available categorical columns in data:\n  "
                + ", ".join(
                    [
                        col
                        for col in adata.obs.columns
                        if adata.obs[col].dtype.name in ["object", "category"]
                    ][:10]
                )
                + "\n\nPlease specify: params={'cluster_key': 'your_column_name'}"
            )

    ripley_key = f"{cluster_key}_ripley_L"
    if ripley_key not in adata.uns:
        raise DataNotFoundError(
            f"Ripley analysis results not found. Expected key: {ripley_key}"
        )

    try:
        import squidpy as sq

        figsize = params.figure_size or (10, 8)
        fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

        sq.pl.ripley(adata, cluster_key=cluster_key, mode="L", plot_sims=True, ax=ax)

        title = params.title or f"Ripley's L Function ({cluster_key})"
        ax.set_title(title)

    except Exception as e:
        # Ripley's L function plotting failed - fail honestly instead of generating fake data
        error_msg = (
            f"Ripley's L function visualization failed: {e}. "
            f"This requires squidpy for proper spatial statistics calculation. "
            f"Please install squidpy (pip install squidpy) and ensure Ripley analysis "
            f"has been performed first using analyze_spatial_statistics with analysis_type='ripley'."
        )
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg)

    plt.tight_layout()
    return fig


@handle_visualization_errors("Moran's I")
async def create_moran_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create Moran's I visualization"""
    if "moranI" not in adata.uns:
        raise DataNotFoundError("Moran's I results not found. Expected key: moranI")

    moran_data = adata.uns["moranI"]

    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    # Plot -log10(p-value) vs Moran's I
    scatter = ax.scatter(
        -np.log10(moran_data["pval_norm"]),
        moran_data["I"],
        s=50,
        alpha=params.alpha,
        c=range(len(moran_data["I"])),
        cmap=params.colormap,
    )

    # Add labels for top significant genes
    gene_names = (
        moran_data.index
        if hasattr(moran_data, "index")
        else range(len(moran_data["I"]))
    )
    for i, gene in enumerate(gene_names[:5]):  # Label top 5 genes
        if moran_data["pval_norm"][i] < 0.05:
            ax.annotate(
                str(gene),
                (-np.log10(moran_data["pval_norm"][i]), moran_data["I"][i]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
            )

    # Add reference lines
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(x=-np.log10(0.05), color="red", linestyle="--", alpha=0.5)

    title = params.title or "Moran's I Spatial Autocorrelation"
    ax.set_title(title)
    ax.set_xlabel("-log10(p-value)")
    ax.set_ylabel("Moran's I")

    if params.show_colorbar:
        plt.colorbar(scatter, ax=ax, label="Gene index")

    plt.tight_layout()
    return fig


@handle_visualization_errors("Centrality Scores")
async def create_centrality_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create centrality scores visualization"""
    # Infer cluster_key from params or from available results in adata.uns
    cluster_key = params.cluster_key
    if not cluster_key:
        # Try to infer from adata.uns keys
        centrality_keys = [
            k for k in adata.uns.keys() if k.endswith("_centrality_scores")
        ]
        if centrality_keys:
            cluster_key = centrality_keys[0].replace("_centrality_scores", "")
            if context:
                await context.info(
                    f"Inferred cluster_key: '{cluster_key}' from existing results"
                )
        else:
            raise ValueError(
                "cluster_key parameter is required but not provided.\n\n"
                "Available categorical columns in data:\n  "
                + ", ".join(
                    [
                        col
                        for col in adata.obs.columns
                        if adata.obs[col].dtype.name in ["object", "category"]
                    ][:10]
                )
                + "\n\nPlease specify: params={'cluster_key': 'your_column_name'}"
            )

    centrality_key = f"{cluster_key}_centrality_scores"
    if centrality_key not in adata.uns:
        raise DataNotFoundError(
            f"Centrality scores not found. Expected key: {centrality_key}"
        )

    centrality_data = adata.uns[centrality_key]
    categories = adata.obs[cluster_key].cat.categories

    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    centrality_types = list(centrality_data.keys())
    x = np.arange(len(categories))
    width = 0.8 / len(centrality_types)

    for i, ctype in enumerate(centrality_types):
        values = centrality_data[ctype]
        offset = width * i - width * len(centrality_types) / 2 + width / 2
        ax.bar(x + offset, values, width, label=ctype, alpha=params.alpha)

    ax.set_xlabel("Cluster")
    ax.set_ylabel("Centrality Score")
    title = params.title or f"Centrality Scores by Cluster ({cluster_key})"
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, rotation=45, ha="right")

    if params.show_legend:
        ax.legend()

    plt.tight_layout()
    return fig


@handle_visualization_errors("Getis-Ord Gi*")
async def create_getis_ord_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create Getis-Ord Gi* visualization"""
    # Find genes with Getis-Ord results
    getis_ord_genes = []
    for col in adata.obs.columns:
        if col.endswith("_getis_ord_z"):
            gene = col.replace("_getis_ord_z", "")
            if f"{gene}_getis_ord_p" in adata.obs.columns:
                getis_ord_genes.append(gene)

    if not getis_ord_genes:
        raise DataNotFoundError("No Getis-Ord results found in adata.obs")

    # Get genes to plot
    feature_list = (
        params.feature
        if isinstance(params.feature, list)
        else ([params.feature] if params.feature else [])
    )
    if feature_list:
        genes_to_plot = [g for g in feature_list if g in getis_ord_genes]
    else:
        genes_to_plot = getis_ord_genes[:6]  # Default to first 6 genes

    if not genes_to_plot:
        raise DataNotFoundError(
            f"None of the specified genes have Getis-Ord results: {feature_list}"
        )

    if context:
        await context.info(
            f"Plotting Getis-Ord results for {len(genes_to_plot)} genes: {genes_to_plot}"
        )

    fig, axes = setup_multi_panel_figure(
        n_panels=len(genes_to_plot),
        params=params,
        default_title="Getis-Ord Gi* Hotspots/Coldspots",
    )

    if "spatial" not in adata.obsm:
        raise DataCompatibilityError(
            "Spatial coordinates not found in adata.obsm['spatial']"
        )

    coords = adata.obsm["spatial"]

    for i, gene in enumerate(genes_to_plot):
        if i < len(axes):
            ax = axes[i]
            z_key = f"{gene}_getis_ord_z"
            p_key = f"{gene}_getis_ord_p"

            if z_key not in adata.obs or p_key not in adata.obs:
                ax.text(
                    0.5,
                    0.5,
                    f"No Getis-Ord data for {gene}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{gene} (No Data)")
                continue

            z_scores = adata.obs[z_key].values
            p_vals = adata.obs[p_key].values

            # Create scatter plot with Z-scores as colors
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=z_scores,
                cmap="RdBu_r",  # Red for hot spots, blue for cold spots
                s=params.spot_size or 20,
                alpha=params.alpha,
                vmin=-3,
                vmax=3,  # Standard Z-score range
            )

            if params.show_colorbar:
                plt.colorbar(scatter, ax=ax, label="Gi* Z-score")

            # Count significant hot and cold spots
            alpha = 0.05  # Default alpha for significance
            significant = p_vals < alpha
            hot_spots = np.sum((z_scores > 0) & significant)
            cold_spots = np.sum((z_scores < 0) & significant)

            if params.add_gene_labels:
                ax.set_title(f"{gene}\nHot: {hot_spots}, Cold: {cold_spots}")
            else:
                ax.set_title(f"{gene}")

            ax.set_xlabel("Spatial X")
            ax.set_ylabel("Spatial Y")
            ax.set_aspect("equal")
            ax.invert_yaxis()  # Proper spatial orientation

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("Trajectory")
async def create_trajectory_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create trajectory pseudotime visualization

    Args:
        adata: AnnData object with computed trajectory/pseudotime
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with trajectory visualization
    """
    if context:
        await context.info("Creating trajectory pseudotime visualization")

    # Find pseudotime key
    pseudotime_key = params.feature
    if not pseudotime_key:
        # Look for pseudotime columns
        pseudotime_candidates = [
            k for k in adata.obs.columns if "pseudotime" in k.lower()
        ]
        if pseudotime_candidates:
            pseudotime_key = pseudotime_candidates[0]
            if context:
                await context.info(f"Found pseudotime column: {pseudotime_key}")
        else:
            raise DataNotFoundError(
                "No pseudotime information found. Please run trajectory analysis first or specify a pseudotime column in the 'feature' parameter."
            )

    if pseudotime_key not in adata.obs.columns:
        raise DataNotFoundError(
            f"Pseudotime column '{pseudotime_key}' not found. Please run trajectory analysis first."
        )

    # Check if RNA velocity is available
    has_velocity = "velocity_graph" in adata.uns

    # Determine basis for plotting
    basis = params.basis or "spatial"
    basis_key = f"X_{basis}" if basis != "spatial" else "spatial"

    if basis_key not in adata.obsm:
        # Try to find an alternative basis
        if "spatial" in adata.obsm:
            basis = "spatial"
            if context:
                await context.info(
                    "Using 'spatial' as basis for trajectory visualization"
                )
        elif "X_umap" in adata.obsm:
            basis = "umap"
            if context:
                await context.info("Using 'umap' as basis for trajectory visualization")
        elif "X_pca" in adata.obsm:
            basis = "pca"
            if context:
                await context.info("Using 'pca' as basis for trajectory visualization")
        else:
            available_bases = [
                k.replace("X_", "") for k in adata.obsm.keys() if k.startswith("X_")
            ]
            available_bases.extend([k for k in adata.obsm.keys() if k == "spatial"])
            raise DataCompatibilityError(
                f"Basis '{params.basis or 'spatial'}' not found. Available bases: {available_bases}"
            )

    # Setup figure: 1 panel if no velocity, 2 panels if velocity exists
    n_panels = 2 if has_velocity else 1

    # USE THE NEW PANEL HELPER
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
        # Don't set individual panel title since we have a comprehensive main title
        # ax1.set_title(f"Pseudotime ({pseudotime_key})", fontsize=12)

        # Handle spatial coordinates orientation
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

            # Handle spatial coordinates orientation
            if basis == "spatial":
                ax2.invert_yaxis()

        except ImportError:
            ax2.text(
                0.5,
                0.5,
                "scvelo not installed\nfor velocity visualization",
                ha="center",
                va="center",
                transform=ax2.transAxes,
            )
            ax2.set_title("Velocity Stream - Error", fontsize=12)
        except Exception as e:
            ax2.text(
                0.5,
                0.5,
                f"Error plotting velocity:\n{str(e)}",
                ha="center",
                va="center",
                transform=ax2.transAxes,
            )
            ax2.set_title("Velocity Stream - Error", fontsize=12)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


@handle_visualization_errors("Gene Correlation")
async def create_gene_correlation_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create gene correlation visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with gene correlation visualization
    """
    # USE THE NEW FEATURE HELPER
    available_genes = await get_validated_features(
        adata, params, max_features=10, context=context
    )

    if context:
        await context.info(
            f"Computing correlations for {len(available_genes)} genes: {available_genes}"
        )

    # Get expression matrix
    expr_matrix = adata[:, available_genes].X
    if hasattr(expr_matrix, "toarray"):
        expr_matrix = expr_matrix.toarray()

    # Apply color scaling
    if params.color_scale == "log":
        expr_matrix = np.log1p(expr_matrix)
    elif params.color_scale == "sqrt":
        expr_matrix = np.sqrt(expr_matrix)

    # Calculate correlation matrix

    n_genes = len(available_genes)
    corr_matrix = np.zeros((n_genes, n_genes))
    p_value_matrix = np.zeros((n_genes, n_genes))

    for i in range(n_genes):
        for j in range(n_genes):
            if i == j:
                corr_matrix[i, j] = 1.0
                p_value_matrix[i, j] = 0.0
            else:
                if params.correlation_method == "pearson":
                    corr, p_val = pearsonr(expr_matrix[:, i], expr_matrix[:, j])
                else:  # spearman
                    corr, p_val = spearmanr(expr_matrix[:, i], expr_matrix[:, j])

                corr_matrix[i, j] = corr
                p_value_matrix[i, j] = p_val

    # Determine figure size
    if params.figure_size:
        figsize = params.figure_size
    else:
        figsize = (max(8, n_genes), max(6, n_genes))

    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, dpi=params.dpi)

    # Set figure title
    title = (
        params.title
        or f"Gene Correlation Analysis ({params.correlation_method.title()})"
    )
    fig.suptitle(title, fontsize=16)

    # Plot correlation heatmap
    import seaborn as sns

    # Correlation heatmap
    sns.heatmap(
        corr_matrix,
        annot=True,
        cmap=params.colormap,
        center=0,
        square=True,
        xticklabels=available_genes,
        yticklabels=available_genes,
        ax=ax1,
        cbar_kws={"shrink": 0.8},
    )
    ax1.set_title(f"{params.correlation_method.title()} Correlation")
    ax1.tick_params(axis="x", rotation=45)
    ax1.tick_params(axis="y", rotation=0)

    # P-value heatmap
    if params.show_correlation_stats:
        # Create significance mask
        sig_mask = p_value_matrix < 0.05

        sns.heatmap(
            -np.log10(p_value_matrix + 1e-10),  # -log10 p-values
            annot=True,
            cmap="Reds",
            square=True,
            xticklabels=available_genes,
            yticklabels=available_genes,
            ax=ax2,
            cbar_kws={"shrink": 0.8, "label": "-log10(p-value)"},
        )
        ax2.set_title("Statistical Significance")
        ax2.tick_params(axis="x", rotation=45)
        ax2.tick_params(axis="y", rotation=0)

        # Add significance markers
        for i in range(n_genes):
            for j in range(n_genes):
                if sig_mask[i, j] and i != j:
                    ax2.text(
                        j + 0.5,
                        i + 0.5,
                        "*",
                        ha="center",
                        va="center",
                        color="white",
                        fontsize=16,
                        fontweight="bold",
                    )
    else:
        # Show clustered correlation
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import squareform

        # Convert correlation to distance
        distance_matrix = 1 - np.abs(corr_matrix)
        condensed_distances = squareform(distance_matrix, checks=False)

        # Perform hierarchical clustering
        linkage_matrix = linkage(condensed_distances, method="average")

        # Plot dendrogram
        dendrogram(linkage_matrix, labels=available_genes, ax=ax2, orientation="top")
        ax2.set_title("Gene Clustering")
        ax2.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    return fig


@handle_visualization_errors("EnrichMap")
async def create_enrichment_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create EnrichMap enrichment visualization

    Supports multiple visualization types:
    - Default: Spatial enrichment map
    - violin: Enrichment scores by cluster
    - heatmap: Gene contributions heatmap
    - correlation: Signature correlation heatmap

    Args:
        adata: AnnData object with enrichment scores
        params: Visualization parameters
            - feature: Score column name or signature name
            - feature: Score column name, signature name, or list of scores for multi-panel plot
            - cluster_key: For violin plots, the grouping variable (default: leiden)
            - show_gene_contributions: Show gene contribution heatmap
        context: MCP context

    Returns:
        Matplotlib figure with enrichment visualization
    """
    if context:
        await context.info("Creating EnrichMap enrichment visualization")

    # Import EnrichMap for specialized visualizations
    try:
        # Import EnrichMap directly (should be installed via pip)
        import enrichmap as em
    except ImportError:
        # Fallback to basic visualization
        if context:
            await context.info(
                "EnrichMap not available. Install with: pip install enrichmap"
            )
        em = None

    # Find enrichment score columns
    score_cols = [col for col in adata.obs.columns if col.endswith("_score")]

    if not score_cols:
        raise DataNotFoundError(
            "No enrichment scores found. Please run 'analyze_enrichment' first."
        )

    # Check if user wants gene contribution visualization
    if hasattr(params, "show_gene_contributions") and params.show_gene_contributions:
        if "gene_contributions" not in adata.uns:
            raise DataNotFoundError(
                "Gene contributions not found. Please run 'analyze_enrichment' first."
            )

        # Create gene contribution heatmap
        gene_contribs = adata.uns["gene_contributions"]

        if not gene_contribs:
            raise DataNotFoundError("No gene contributions available.")

        # Convert to DataFrame for visualization
        import pandas as pd

        contrib_data = {}
        for sig, genes in gene_contribs.items():
            contrib_data[sig] = genes

        if not contrib_data:
            raise DataNotFoundError("No gene contribution data to visualize.")

        df = pd.DataFrame(contrib_data).fillna(0)

        # Create heatmap
        fig, ax = create_figure(
            figsize=(max(8, len(df.columns) * 1.5), max(6, len(df) * 0.3))
        )

        sns.heatmap(
            df,
            annot=True,
            fmt=".3f",
            cmap="RdBu_r",
            center=0,
            ax=ax,
            cbar_kws={"label": "Gene Weight"},
        )
        ax.set_title("Gene Contributions to Enrichment Scores", fontsize=14)
        ax.set_xlabel("Signatures", fontsize=12)
        ax.set_ylabel("Genes", fontsize=12)

        plt.tight_layout()
        return fig

    # Check if user wants violin plot by cluster
    if params.plot_type == "violin" or (
        hasattr(params, "show_violin") and params.show_violin
    ):
        # Determine grouping variable - require explicit specification
        if hasattr(params, "cluster_key") and params.cluster_key:
            group_by = params.cluster_key
        else:
            # No cluster_key provided - show error with available options
            categorical_cols = [
                col
                for col in adata.obs.columns
                if adata.obs[col].dtype.name in ["object", "category"]
            ]
            raise ValueError(
                "❌ Enrichment violin plot requires 'cluster_key' parameter.\n\n"
                f"Available categorical columns ({len(categorical_cols)} total):\n"
                f"  {', '.join(categorical_cols[:15])}\n\n"
                "🔧 SOLUTION: Specify cluster_key explicitly:\n"
                "  params={'cluster_key': 'your_column_name'}"
            )

        if group_by not in adata.obs.columns:
            raise DataNotFoundError(
                f"Grouping variable '{group_by}' not found in adata.obs"
            )

        # Determine which scores to plot
        feature_list = (
            params.feature
            if isinstance(params.feature, list)
            else ([params.feature] if params.feature else [])
        )
        if feature_list and len(feature_list) > 1:
            scores_to_plot = []
            for feat in feature_list:
                if feat in adata.obs.columns:
                    scores_to_plot.append(feat)
                elif f"{feat}_score" in adata.obs.columns:
                    scores_to_plot.append(f"{feat}_score")
        elif feature_list:
            if params.feature in adata.obs.columns:
                scores_to_plot = [params.feature]
            elif f"{params.feature}_score" in adata.obs.columns:
                scores_to_plot = [f"{params.feature}_score"]
            else:
                scores_to_plot = [score_cols[0]]
        else:
            scores_to_plot = score_cols[:3]  # Limit to 3 scores for clarity

        # Create violin plot
        n_scores = len(scores_to_plot)
        fig, axes = plt.subplots(1, n_scores, figsize=(5 * n_scores, 6))

        if n_scores == 1:
            axes = [axes]

        for i, score in enumerate(scores_to_plot):
            ax = axes[i]

            # Create violin plot using seaborn
            data_for_plot = pd.DataFrame(
                {group_by: adata.obs[group_by], "Score": adata.obs[score]}
            )

            sns.violinplot(data=data_for_plot, x=group_by, y="Score", ax=ax)

            sig_name = score.replace("_score", "")
            ax.set_title(f"{sig_name} by {group_by}", fontsize=12)
            ax.set_xlabel(group_by, fontsize=10)
            ax.set_ylabel("Enrichment Score", fontsize=10)
            ax.tick_params(axis="x", rotation=45)
            # Align rotated labels to the right for better readability
            for label in ax.get_xticklabels():
                label.set_horizontalalignment("right")

        plt.tight_layout()
        return fig

    # Check if EnrichMap is available for advanced visualizations
    if em is not None and hasattr(params, "enrichmap_plot_type"):
        plot_type = params.enrichmap_plot_type

        # Determine score to visualize
        if params.feature:
            if params.feature in adata.obs.columns:
                score_col = params.feature
            elif f"{params.feature}_score" in adata.obs.columns:
                score_col = f"{params.feature}_score"
            else:
                score_col = score_cols[0]
        else:
            score_col = score_cols[0]

        fig = plt.figure(figsize=(10, 8))

        try:
            if plot_type == "correlogram":
                # Moran's I correlogram
                em.pl.morans_correlogram(adata, score_key=score_col)
            elif plot_type == "variogram":
                # Variogram
                em.pl.variogram(adata, score_key=score_col)
            elif plot_type == "cross_correlation":
                # Cross-correlation between scores
                if len(score_cols) >= 2:
                    em.pl.cross_moran_scatter(adata, score_keys=score_cols[:2])
                else:
                    raise DataNotFoundError(
                        "Need at least 2 enrichment scores for cross-correlation"
                    )
            else:
                # Default spatial enrichment
                em.pl.spatial_enrichmap(adata, score_key=score_col)
        except Exception as e:
            if context:
                await context.warning(
                    f"EnrichMap visualization failed: {e}. Using fallback."
                )
            plt.close(fig)
            # Fall back to standard visualization
            fig, ax = create_figure(figsize=(10, 8))
            plot_spatial_feature(adata, feature=score_col, ax=ax, params=params)

        return fig

    # Standard spatial visualization (original code)
    # Determine which score to visualize
    if params.feature:
        # User specified a score
        if params.feature in adata.obs.columns:
            score_col = params.feature
        elif f"{params.feature}_score" in adata.obs.columns:
            score_col = f"{params.feature}_score"
        else:
            raise DataNotFoundError(
                f"Score column '{params.feature}' not found. Available scores: {score_cols}"
            )
    else:
        # Use the first available score
        score_col = score_cols[0]
        if context:
            await context.info(f"Using score column: {score_col}")

    # Check if we should create multi-panel plot for multiple scores
    feature_list = (
        params.feature
        if isinstance(params.feature, list)
        else ([params.feature] if params.feature else [])
    )
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

        # Create multi-panel figure
        fig, axes = setup_multi_panel_figure(
            n_panels=len(scores_to_plot),
            params=params,
            default_title="EnrichMap Enrichment Scores",
        )

        # Plot each score
        for i, score in enumerate(scores_to_plot):
            if i < len(axes):
                ax = axes[i]
                plot_spatial_feature(adata, feature=score, ax=ax, params=params)

                # Extract signature name from score column
                sig_name = score.replace("_score", "")
                ax.set_title(f"{sig_name} Enrichment")
    else:
        # Single score visualization
        fig, ax = create_figure(figsize=(10, 8))
        plot_spatial_feature(adata, feature=score_col, ax=ax, params=params)

        # Extract signature name from score column
        sig_name = score_col.replace("_score", "")
        ax.set_title(f"{sig_name} Enrichment Score", fontsize=14)

        # Add colorbar label if shown
        if params.show_colorbar:
            # Find the colorbar and update its label
            if hasattr(ax, "collections") and ax.collections:
                cbar = plt.colorbar(ax.collections[0], ax=ax)
                cbar.set_label("Enrichment Score", fontsize=12)

    plt.tight_layout()
    return fig


@handle_visualization_errors("GSEA")
async def create_gsea_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create GSEA (Gene Set Enrichment Analysis) visualization

    Supports multiple visualization types:
    - enrichment_plot: Classic GSEA enrichment score plot
    - barplot: Top enriched pathways barplot
    - dotplot: Multi-cluster enrichment dotplot
    - spatial: Spatial distribution of enrichment scores

    Args:
        adata: AnnData object with GSEA results
        params: Visualization parameters
            - gsea_results_key: Key in adata.uns for GSEA results (default: 'gsea_results')
            - gsea_plot_type: Type of plot ('enrichment_plot', 'barplot', 'dotplot', 'spatial')
            - feature: Specific pathway/gene set to visualize
            - n_top_pathways: Number of top pathways to show (default: 10)
        context: MCP context

    Returns:
        Matplotlib figure with GSEA visualization
    """
    if context:
        await context.info("Creating GSEA visualization")

    # Get GSEA results
    gsea_key = getattr(params, "gsea_results_key", "gsea_results")
    if gsea_key not in adata.uns:
        # Try common alternative keys
        alt_keys = ["rank_genes_groups", "de_results", "pathway_enrichment"]
        for key in alt_keys:
            if key in adata.uns:
                gsea_key = key
                break
        else:
            raise DataNotFoundError(f"GSEA results not found. Expected key: {gsea_key}")

    gsea_results = adata.uns[gsea_key]
    plot_type = getattr(params, "gsea_plot_type", "barplot")

    if plot_type == "enrichment_plot":
        # Classic GSEA enrichment score plot
        return _create_gsea_enrichment_plot(adata, gsea_results, params, context)
    elif plot_type == "barplot":
        # Top pathways barplot
        return _create_gsea_barplot(adata, gsea_results, params, context)
    elif plot_type == "dotplot":
        # Multi-cluster dotplot
        return _create_gsea_dotplot(adata, gsea_results, params, context)
    elif plot_type == "spatial":
        # Spatial distribution of pathway scores
        return _create_gsea_spatial_plot(adata, gsea_results, params, context)
    else:
        # Default to barplot
        return _create_gsea_barplot(adata, gsea_results, params, context)


def _create_gsea_enrichment_plot(adata, gsea_results, params, context):
    """Create classic GSEA enrichment score plot"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), height_ratios=[3, 1])

    # Get specific pathway if specified
    pathway = params.feature if params.feature else None

    if isinstance(gsea_results, dict):
        if pathway and pathway in gsea_results:
            result = gsea_results[pathway]
        else:
            # Use first available pathway
            pathway = list(gsea_results.keys())[0]
            result = gsea_results[pathway]
    else:
        # Assume it's a single result
        result = gsea_results
        pathway = pathway or "Gene Set"

    # Extract enrichment data - fail if data is missing
    if not isinstance(result, dict):
        raise ValueError("Enrichment result is not a valid dictionary structure")

    # Require actual enrichment data - no fake generation
    if "es" not in result:
        raise ValueError(
            "Enrichment score (es) not found in result. Cannot visualize without real enrichment data."
        )

    es = result["es"]
    nes = result.get("nes", result.get("NES"))
    pval = result.get("pval", result.get("pvalue"))
    positions = result.get("positions")

    if nes is None or pval is None:
        raise ValueError(
            "Missing required enrichment statistics (NES or p-value) for visualization."
        )

    # Normalize ES if needed
    if isinstance(es, (list, np.ndarray)):
        es = np.array(es)
        if es.max() > 1 or es.min() < -1:
            es = es / np.abs(es).max()

    # Plot enrichment score
    x = np.arange(len(es))
    ax1.plot(x, es, "b-", linewidth=2)
    ax1.axhline(y=0, color="gray", linestyle="-", alpha=0.5)

    # Mark maximum/minimum ES
    if es.max() > abs(es.min()):
        max_idx = np.argmax(es)
        ax1.axhline(y=es[max_idx], color="red", linestyle="--", alpha=0.8)
        ax1.text(
            max_idx,
            es[max_idx],
            f"ES = {es[max_idx]:.3f}",
            ha="left",
            va="bottom",
            fontsize=10,
        )
    else:
        min_idx = np.argmin(es)
        ax1.axhline(y=es[min_idx], color="red", linestyle="--", alpha=0.8)
        ax1.text(
            min_idx,
            es[min_idx],
            f"ES = {es[min_idx]:.3f}",
            ha="left",
            va="top",
            fontsize=10,
        )

    # Fill under curve
    ax1.fill_between(x, 0, es, alpha=0.3, color="blue")

    ax1.set_xlim(0, len(es) - 1)
    ax1.set_ylabel("Enrichment Score (ES)", fontsize=12)
    ax1.set_title(
        f"{pathway}\nNES = {nes:.2f}, p-value = {pval:.3e}",
        fontsize=14,
        fontweight="bold",
    )

    # Gene hit positions
    ax2.eventplot(positions, colors="black", linewidths=0.5)
    ax2.set_xlim(0, len(es) - 1)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("Rank in Ordered Gene List", fontsize=12)
    ax2.set_ylabel("Hits", fontsize=10)
    ax2.set_yticks([])

    plt.tight_layout()
    return fig


def _create_gsea_barplot(adata, gsea_results, params, context):
    """Create barplot of top enriched pathways"""
    n_top = getattr(params, "n_top_pathways", 10)

    # Convert results to DataFrame for easier handling
    import pandas as pd

    if isinstance(gsea_results, pd.DataFrame):
        df = gsea_results
    elif isinstance(gsea_results, dict):
        # Convert dict to DataFrame
        rows = []
        for pathway, data in gsea_results.items():
            if isinstance(data, dict):
                row = {"pathway": pathway}
                row.update(data)
                rows.append(row)
        df = pd.DataFrame(rows)
    else:
        raise ValueError("Unsupported GSEA results format")

    # Determine score column
    score_cols = ["NES", "nes", "enrichment_score", "score"]
    score_col = None
    for col in score_cols:
        if col in df.columns:
            score_col = col
            break

    if not score_col:
        raise DataNotFoundError("No enrichment score column found in results")

    # Convert score column to numeric type
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    # Remove rows with NaN scores
    df = df.dropna(subset=[score_col])

    if df.empty:
        raise DataNotFoundError("No valid enrichment scores found in results")

    # Sort by score and get top pathways
    df_sorted = df.nlargest(min(n_top, len(df)), score_col)

    # Create barplot
    fig, ax = plt.subplots(figsize=(10, max(6, n_top * 0.4)))

    y_pos = np.arange(len(df_sorted))
    scores = df_sorted[score_col].values
    # Handle different column names for pathways
    if "pathway" in df_sorted.columns:
        pathways = df_sorted["pathway"].values
    elif "Term" in df_sorted.columns:
        pathways = df_sorted["Term"].values
    else:
        pathways = df_sorted.index

    # Color based on score
    colors = [
        "darkred" if s > 2 else "red" if s > 1.5 else "orange" if s > 1 else "gray"
        for s in scores
    ]

    ax.barh(y_pos, scores, color=colors, alpha=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(pathways, fontsize=10)
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=12)
    ax.set_title("Top Enriched Pathways", fontsize=14, fontweight="bold")
    ax.axvline(x=0, color="black", linewidth=0.8)

    # Add significance markers
    if "pval" in df_sorted.columns or "pvalue" in df_sorted.columns:
        pval_col = "pval" if "pval" in df_sorted.columns else "pvalue"
        pvals = df_sorted[pval_col].values

        for i, (score, pval) in enumerate(zip(scores, pvals)):
            if pval < 0.001:
                sig_marker = "***"
            elif pval < 0.01:
                sig_marker = "**"
            elif pval < 0.05:
                sig_marker = "*"
            else:
                sig_marker = ""

            if sig_marker:
                ax.text(
                    score + max(scores) * 0.02, i, sig_marker, va="center", fontsize=10
                )

    # Add FDR info if available
    if "FDR" in df_sorted.columns or "fdr" in df_sorted.columns:
        ax.text(
            0.95,
            0.05,
            "FDR < 0.05",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

    plt.tight_layout()
    return fig


def _create_gsea_dotplot(adata, gsea_results, params, context):
    """Create dotplot showing enrichment across multiple conditions"""
    fig, ax = plt.subplots(figsize=(12, 8))

    # This assumes results are organized by cluster/condition
    if not isinstance(gsea_results, dict) or not all(
        isinstance(v, dict) for v in gsea_results.values()
    ):
        raise ValueError(
            "Dotplot requires nested dict structure: {condition: {pathway: scores}}"
        )

    # Extract data for plotting
    conditions = list(gsea_results.keys())
    all_pathways = set()
    for cond_results in gsea_results.values():
        all_pathways.update(cond_results.keys())
    pathways = sorted(all_pathways)[:20]  # Limit to top 20 pathways

    # Create matrices for plotting
    nes_matrix = np.zeros((len(conditions), len(pathways)))
    pval_matrix = np.ones((len(conditions), len(pathways)))

    for i, condition in enumerate(conditions):
        for j, pathway in enumerate(pathways):
            if pathway in gsea_results[condition]:
                result = gsea_results[condition][pathway]
                nes_matrix[i, j] = result.get("NES", result.get("nes", 0))
                pval_matrix[i, j] = result.get("pval", result.get("pvalue", 1))

    # Create dotplot
    for i, condition in enumerate(conditions):
        for j, pathway in enumerate(pathways):
            nes = nes_matrix[i, j]
            pval = pval_matrix[i, j]

            if pval < 0.05:  # Only show significant results
                size = -np.log10(pval) * 100  # Size based on significance
                color = nes  # Color based on NES

                scatter = ax.scatter(
                    j,
                    i,
                    s=size,
                    c=color,
                    cmap="RdBu_r",
                    vmin=-3,
                    vmax=3,
                    alpha=0.8,
                    edgecolors="black",
                    linewidth=0.5,
                )

    # Formatting
    ax.set_xticks(range(len(pathways)))
    ax.set_xticklabels(pathways, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(conditions)))
    ax.set_yticklabels(conditions, fontsize=10)
    ax.set_xlim(-0.5, len(pathways) - 0.5)
    ax.set_ylim(-0.5, len(conditions) - 0.5)
    ax.set_title("Pathway Enrichment Across Conditions", fontsize=14, fontweight="bold")

    # Add colorbar
    plt.colorbar(scatter, ax=ax, label="NES")

    # Add size legend
    for size, label in [(100, "p < 0.1"), (200, "p < 0.01"), (300, "p < 0.001")]:
        ax.scatter(
            [],
            [],
            s=size,
            c="gray",
            alpha=0.8,
            edgecolors="black",
            linewidth=0.5,
            label=label,
        )
    ax.legend(loc="upper left", bbox_to_anchor=(1.2, 1), title="Significance")

    plt.tight_layout()
    return fig


def _create_gsea_spatial_plot(adata, gsea_results, params, context):
    """Create spatial visualization of pathway enrichment scores"""
    # Get spatial coordinates
    x_coords, y_coords = get_spatial_coordinates(adata)

    # Get pathway scores (these should be pre-computed and stored in adata.obs)
    pathway = params.feature if params.feature else None

    # Look for pathway scores in obs
    pathway_score_cols = [
        col
        for col in adata.obs.columns
        if "pathway" in col.lower() or "gsea" in col.lower()
    ]

    if pathway and f"{pathway}_score" in adata.obs.columns:
        score_col = f"{pathway}_score"
    elif pathway and pathway in adata.obs.columns:
        score_col = pathway
    elif pathway_score_cols:
        score_col = pathway_score_cols[0]
    else:
        # No pathway scores available - fail honestly
        raise ValueError(
            "No pathway scores found in the dataset. "
            "Available columns: " + ", ".join(adata.obs.columns) + ". "
            "Cannot visualize pathway scores without real enrichment analysis results."
        )

    scores = adata.obs[score_col].values

    # Create spatial plot
    fig, ax = plt.subplots(figsize=(10, 8))

    scatter = ax.scatter(
        x_coords,
        y_coords,
        c=scores,
        cmap="RdBu_r",
        s=50,
        alpha=0.8,
        vmin=np.percentile(scores, 5),
        vmax=np.percentile(scores, 95),
    )

    ax.set_xlabel("Spatial X")
    ax.set_ylabel("Spatial Y")
    ax.set_title(
        f'Spatial Distribution of {score_col.replace("_", " ").title()}',
        fontsize=14,
        fontweight="bold",
    )
    ax.set_aspect("equal")

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Enrichment Score", fontsize=12)

    # Add statistics
    stats_text = f"Mean: {scores.mean():.3f}\nStd: {scores.std():.3f}"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        verticalalignment="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    return fig


async def create_spatial_interaction_visualization(
    adata, params, context: Optional[Context] = None
):
    """Create spatial visualization showing ligand-receptor interactions"""

    if context:
        await context.info("Creating spatial interaction plot")

    try:
        # Get spatial coordinates
        spatial_coords = adata.obsm["spatial"]

        # Validate that lr_pairs are provided
        if not params.lr_pairs or len(params.lr_pairs) == 0:
            raise ValueError(
                "No ligand-receptor pairs provided for spatial interaction visualization"
            )

        # Create figure
        fig, ax = create_figure(figsize=(12, 10))

        # Plot all cells as background
        ax.scatter(
            spatial_coords[:, 0],
            spatial_coords[:, 1],
            c="lightgray",
            s=10,
            alpha=0.5,
            label="All cells",
        )

        # Color mapping for different LR pairs
        colors = plt.cm.Set3(np.linspace(0, 1, len(params.lr_pairs)))

        interaction_count = 0
        for i, (ligand, receptor) in enumerate(params.lr_pairs):
            color = colors[i]

            # Check if ligand and receptor are in gene expression data
            if ligand in adata.var_names and receptor in adata.var_names:
                # Get cells expressing ligand and receptor
                ligand_expr = (
                    adata[:, ligand].X.toarray().flatten()
                    if hasattr(adata.X, "toarray")
                    else adata[:, ligand].X.flatten()
                )
                receptor_expr = (
                    adata[:, receptor].X.toarray().flatten()
                    if hasattr(adata.X, "toarray")
                    else adata[:, receptor].X.flatten()
                )

                # Define expression threshold (e.g., > median expression)
                ligand_threshold = (
                    np.median(ligand_expr[ligand_expr > 0])
                    if np.any(ligand_expr > 0)
                    else 0
                )
                receptor_threshold = (
                    np.median(receptor_expr[receptor_expr > 0])
                    if np.any(receptor_expr > 0)
                    else 0
                )

                # Find expressing cells
                ligand_cells = ligand_expr > ligand_threshold
                receptor_cells = receptor_expr > receptor_threshold

                if np.any(ligand_cells) and np.any(receptor_cells):
                    # Plot ligand-expressing cells
                    ligand_coords = spatial_coords[ligand_cells]
                    ax.scatter(
                        ligand_coords[:, 0],
                        ligand_coords[:, 1],
                        c=color,
                        s=50,
                        alpha=0.7,
                        marker="o",
                        label=f"{ligand}+ (Ligand)",
                    )

                    # Plot receptor-expressing cells
                    receptor_coords = spatial_coords[receptor_cells]
                    ax.scatter(
                        receptor_coords[:, 0],
                        receptor_coords[:, 1],
                        c=color,
                        s=50,
                        alpha=0.7,
                        marker="^",
                        label=f"{receptor}+ (Receptor)",
                    )

                    # Draw lines between nearby ligand and receptor cells (within distance threshold)
                    if len(ligand_coords) > 0 and len(receptor_coords) > 0:
                        # Calculate pairwise distances
                        from scipy.spatial.distance import cdist

                        distances = cdist(ligand_coords, receptor_coords)

                        # Set distance threshold (e.g., 10th percentile of all distances)
                        distance_threshold = np.percentile(distances, 10)

                        # Draw connections
                        ligand_indices, receptor_indices = np.where(
                            distances <= distance_threshold
                        )
                        for li, ri in zip(
                            ligand_indices[:50], receptor_indices[:50]
                        ):  # Limit connections to avoid clutter
                            ax.plot(
                                [ligand_coords[li, 0], receptor_coords[ri, 0]],
                                [ligand_coords[li, 1], receptor_coords[ri, 1]],
                                color=color,
                                alpha=0.3,
                                linewidth=0.5,
                            )
                            interaction_count += 1

            elif context:
                await context.warning(
                    f"Genes {ligand} or {receptor} not found in expression data"
                )

        ax.set_xlabel("Spatial X")
        ax.set_ylabel("Spatial Y")
        ax.set_title(
            f"Spatial Ligand-Receptor Interactions\n({interaction_count} connections shown)",
            fontsize=14,
            fontweight="bold",
        )
        ax.set_aspect("equal")
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()
        return fig

    except Exception as e:
        if context:
            await context.warning(
                f"Failed to create spatial interaction plot: {str(e)}"
            )
        # Fallback to basic spatial plot
        fig, ax = create_figure()
        ax.text(
            0.5,
            0.5,
            f"Spatial interaction visualization failed:\n{str(e)}",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
        )
        ax.set_title("Spatial Interaction (Error)", fontsize=14)
        return fig


@handle_visualization_errors("Batch Integration")
async def create_batch_integration_visualization(
    adata, params, context: Optional[Context] = None
):
    """Create multi-panel visualization to assess batch integration quality

    This visualization is specifically for evaluating the quality of batch correction
    after integrating multiple samples. It requires proper batch information.
    """

    if context:
        await context.info("Creating batch integration quality visualization")

    # Check if batch information exists - STRICT validation, no fallback
    batch_key = params.batch_key
    if batch_key not in adata.obs.columns:
        raise DataNotFoundError(
            f"Batch key '{batch_key}' not found in data. "
            f"This visualization requires proper batch information from sample integration. "
            f"Available columns: {', '.join(adata.obs.columns[:10])}{'...' if len(adata.obs.columns) > 10 else ''}. "
            f"Please ensure you have run 'integrate_samples' first or specify the correct batch_key."
        )

    # Create multi-panel figure (2x2 layout)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Panel 1: UMAP colored by batch (shows mixing)
    if "X_umap" in adata.obsm:
        umap_coords = adata.obsm["X_umap"]
        batch_values = adata.obs[batch_key]
        unique_batches = batch_values.unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_batches)))

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
        axes[0, 0].text(
            0.5, 0.5, "UMAP coordinates not available", ha="center", va="center"
        )
        axes[0, 0].set_title("UMAP (Not Available)", fontsize=12)

    # Panel 2: Spatial plot colored by batch (if spatial data available)
    if "spatial" in adata.obsm:
        spatial_coords = adata.obsm["spatial"]
        batch_values = adata.obs[batch_key]
        unique_batches = batch_values.unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_batches)))

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
        from scipy.stats import entropy

        # Divide UMAP space into grid and calculate batch entropy per grid cell
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


# ============================================================================
# FILE PERSISTENCE FUNCTIONS
# ============================================================================


async def save_visualization(
    data_id: str,
    plot_type: str,
    output_dir: str = "./outputs",
    filename: Optional[str] = None,
    format: str = "png",
    dpi: Optional[int] = None,
    visualization_cache: Optional[Dict[str, Any]] = None,
    data_store: Optional[Dict[str, Any]] = None,
    regenerate_high_quality: bool = False,
    context: Optional[Context] = None,
) -> str:
    """Save a visualization from cache to disk

    Args:
        data_id: Dataset ID
        plot_type: Type of plot to save
        output_dir: Directory to save the file (default: ./outputs)
        filename: Custom filename (optional, auto-generated if not provided)
        format: Image format (png, jpg, pdf, svg)
        dpi: DPI for saved image (default: 100 for png/jpg, 300 for pdf/svg)
              For publication quality, use 300+ DPI
        visualization_cache: Cache dictionary containing visualizations
        data_store: Data store for regenerating visualization at higher quality
        regenerate_high_quality: If True, regenerate the plot at specified DPI
        context: MCP context for logging

    Returns:
        Path to the saved file

    Raises:
        DataNotFoundError: If visualization not found in cache
        ProcessingError: If saving fails

    Examples:
        # Save for publication (high DPI)
        save_visualization("data1", "spatial", dpi=300, format="pdf")

        # Save with custom filename
        save_visualization("data1", "umap", filename="figure_1a")
    """
    try:
        # Use environment variable for output_dir if default value was passed
        if output_dir == "./outputs":
            output_dir = get_output_dir_from_config(default="./outputs")
            if context and output_dir != "./outputs":
                await context.info(
                    f"Using output directory from configuration: {output_dir}"
                )

        # Validate format
        valid_formats = ["png", "jpg", "jpeg", "pdf", "svg"]
        if format.lower() not in valid_formats:
            raise InvalidParameterError(
                f"Invalid format: {format}. Must be one of {valid_formats}"
            )

        # Check cache
        if visualization_cache is None:
            raise ProcessingError("Visualization cache not provided")

        cache_key = f"{data_id}_{plot_type}"
        if cache_key not in visualization_cache:
            available_keys = [
                k for k in visualization_cache.keys() if k.startswith(data_id)
            ]
            if not available_keys:
                raise DataNotFoundError(
                    f"No visualizations found for dataset '{data_id}'"
                )
            else:
                raise DataNotFoundError(
                    f"Visualization '{plot_type}' not found for dataset '{data_id}'. "
                    f"Available: {', '.join([k.replace(data_id + '_', '') for k in available_keys])}"
                )

        # Set default DPI based on format
        if dpi is None:
            if format.lower() in ["pdf", "svg"]:
                dpi = 300  # High quality for vector formats
            else:
                dpi = 100  # Standard quality for raster formats

        # For publication quality, recommend at least 300 DPI
        if dpi >= 300 and context:
            await context.info(f"Using publication-quality DPI: {dpi}")

        # Create output directory using safe path handling
        # This resolves relative paths against project root (not cwd)
        # and automatically falls back to /tmp if the directory is not writable
        try:
            output_path = get_safe_output_path(
                output_dir, fallback_to_tmp=True, create_if_missing=True
            )

            # Inform user if fallback was used
            if context and str(output_path) != str(Path(output_dir).resolve()):
                await context.info(
                    f"Using output directory: {output_path} "
                    f"(fallback from requested path)"
                )

        except PermissionError as e:
            raise ProcessingError(
                f"Cannot save visualization to {output_dir}: {str(e)}. "
                f"Please check directory permissions or specify a different output_dir."
            ) from e

        # Generate filename if not provided
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            if dpi != 100:
                filename = f"{data_id}_{plot_type}_{dpi}dpi_{timestamp}.{format}"
            else:
                filename = f"{data_id}_{plot_type}_{timestamp}.{format}"
        else:
            # Ensure filename has correct extension
            if not filename.endswith(f".{format}"):
                filename = f"{filename}.{format}"

        # Full path for the file
        file_path = output_path / filename

        # Try to get cached figure object for high-quality export

        from ..utils.image_utils import get_cached_figure, load_figure_pickle

        cached_fig = None

        # 1. Try to get figure from in-memory cache
        if dpi > 100 or format.lower() in ["pdf", "svg"]:
            cached_fig = get_cached_figure(cache_key)

            # 2. If not in memory, try pickle file
            if cached_fig is None:
                pickle_path = f"/tmp/chatspatial/figures/{cache_key}.pkl"
                if os.path.exists(pickle_path):
                    try:
                        cached_fig = load_figure_pickle(pickle_path)
                        if context:
                            await context.info(
                                "Loaded figure from pickle cache for high-quality export"
                            )
                    except Exception as e:
                        if context:
                            await context.warning(
                                f"Failed to load figure from pickle: {str(e)}"
                            )

        # If we have the original figure, regenerate at high quality
        if cached_fig is not None:
            if context:
                await context.info(
                    f"Regenerating visualization at {dpi} DPI for {format} format from cached figure..."
                )

            try:
                cached_fig.savefig(
                    str(file_path),
                    format=format,
                    dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white",
                )

                # Get file size
                file_size_kb = os.path.getsize(file_path) / 1024

                if context:
                    await context.info(
                        f"✅ Saved high-quality visualization to {file_path} ({file_size_kb:.1f} KB at {dpi} DPI)"
                    )

                return str(file_path)

            except Exception as e:
                if context:
                    await context.warning(
                        f"Failed to regenerate from figure: {str(e)}. Falling back to cached version."
                    )

        # Fallback: Save cached compressed version
        if context:
            await context.warning(
                "⚠️  Original figure not found in cache. Saving compressed version from MCP transmission."
            )
            if dpi > 100:
                await context.warning(
                    f"⚠️  Requested {dpi} DPI but saving at lower quality (~80 DPI). "
                    f"To save at {dpi} DPI, regenerate the visualization immediately before saving."
                )

        viz_data = visualization_cache[cache_key]
        if not viz_data:
            raise ProcessingError("Visualization data is empty")

        with open(file_path, "wb") as f:
            f.write(viz_data)

        if context:
            file_size_kb = len(viz_data) / 1024
            await context.info(
                f"Saved visualization to {file_path} ({file_size_kb:.1f} KB)"
            )
            await context.info(
                "💡 Tip: For publication-quality exports, save immediately after generating the visualization."
            )

        return str(file_path)

    except (DataNotFoundError, InvalidParameterError):
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to save visualization: {str(e)}")


async def list_saved_visualizations(
    output_dir: str = "./outputs",
    pattern: Optional[str] = None,
    context: Optional[Context] = None,
) -> List[Dict[str, Any]]:
    """List all saved visualizations in the output directory

    Args:
        output_dir: Directory to search for saved files
        pattern: Optional glob pattern to filter files (e.g., "*spatial*")
        context: MCP context for logging

    Returns:
        List of dictionaries with file information
    """
    try:
        # Resolve path against current working directory (respects user config)
        user_path = Path(output_dir)
        if user_path.is_absolute():
            output_path = user_path
        else:
            output_path = Path.cwd() / user_path

        # Check if directory exists (don't create, just list)
        if not output_path.exists():
            if context:
                await context.info(f"Output directory does not exist: {output_path}")
            return []

        # Security check using new helper
        if not is_safe_output_path(output_path):
            raise InvalidParameterError(
                f"Output directory must be within project or /tmp/chatspatial: {output_path}"
            )

        # Get all image files
        image_extensions = ["*.png", "*.jpg", "*.jpeg", "*.pdf", "*.svg"]
        files = []

        for ext in image_extensions:
            search_pattern = pattern if pattern else ext
            for file_path in output_path.glob(search_pattern if pattern else ext):
                if file_path.is_file():
                    stat = file_path.stat()
                    files.append(
                        {
                            "filename": file_path.name,
                            "path": str(file_path),
                            "size_kb": stat.st_size / 1024,
                            "modified": datetime.fromtimestamp(
                                stat.st_mtime
                            ).isoformat(),
                            "format": file_path.suffix[1:],  # Remove the dot
                        }
                    )

        # Sort by modification time (newest first)
        files.sort(key=lambda x: x["modified"], reverse=True)

        if context:
            await context.info(
                f"Found {len(files)} saved visualization(s) in {output_dir}"
            )

        return files

    except InvalidParameterError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to list visualizations: {str(e)}")


async def export_all_visualizations(
    data_id: str,
    output_dir: str = "./exports",
    visualization_cache: Optional[Dict[str, Any]] = None,
    context: Optional[Context] = None,
) -> List[str]:
    """Export all cached visualizations for a dataset to disk

    Args:
        data_id: Dataset ID to export visualizations for
        output_dir: Directory to save files
        visualization_cache: Cache dictionary containing visualizations
        context: MCP context for logging

    Returns:
        List of paths to saved files
    """
    try:
        if visualization_cache is None:
            raise ProcessingError("Visualization cache not provided")

        # Find all visualizations for this dataset
        relevant_keys = [
            k for k in visualization_cache.keys() if k.startswith(f"{data_id}_")
        ]

        if not relevant_keys:
            if context:
                await context.warning(
                    f"No visualizations found for dataset '{data_id}'"
                )
            return []

        saved_files = []

        for cache_key in relevant_keys:
            # Extract plot type from cache key
            plot_type = cache_key.replace(f"{data_id}_", "")

            try:
                saved_path = await save_visualization(
                    data_id=data_id,
                    plot_type=plot_type,
                    output_dir=output_dir,
                    visualization_cache=visualization_cache,
                    context=context,
                )
                saved_files.append(saved_path)
            except Exception as e:
                if context:
                    await context.warning(f"Failed to export {cache_key}: {str(e)}")

        if context:
            await context.info(
                f"Exported {len(saved_files)} visualization(s) for dataset '{data_id}' to {output_dir}"
            )

        return saved_files

    except ProcessingError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to export visualizations: {str(e)}")


async def clear_visualization_cache(
    data_id: Optional[str] = None,
    visualization_cache: Optional[Dict[str, Any]] = None,
    context: Optional[Context] = None,
) -> int:
    """Clear visualization cache to free memory

    Args:
        data_id: Optional dataset ID to clear specific visualizations
        visualization_cache: Cache dictionary to clear
        context: MCP context for logging

    Returns:
        Number of visualizations cleared
    """
    try:
        if visualization_cache is None:
            return 0

        if data_id:
            # Clear specific dataset visualizations
            keys_to_remove = [
                k for k in visualization_cache.keys() if k.startswith(f"{data_id}_")
            ]
            for key in keys_to_remove:
                del visualization_cache[key]
            cleared_count = len(keys_to_remove)

            if context:
                await context.info(
                    f"Cleared {cleared_count} visualization(s) for dataset '{data_id}'"
                )
        else:
            # Clear all visualizations
            cleared_count = len(visualization_cache)
            visualization_cache.clear()

            if context:
                await context.info(
                    f"Cleared all {cleared_count} visualization(s) from cache"
                )

        return cleared_count

    except Exception as e:
        raise ProcessingError(f"Failed to clear cache: {str(e)}")
