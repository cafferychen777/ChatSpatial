"""
Visualization tools for spatial transcriptomics data.
"""

import os
import traceback
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Tuple, Union

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

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

from ..models.data import VisualizationParameters  # noqa: E402
# Import adata utilities (spatial coordinates, categorical handling)
from ..utils.adata_utils import (ensure_categorical,  # noqa: E402
                                 get_spatial_coordinates)
# Use centralized dependency manager
from ..utils.dependency_manager import require  # noqa: E402
# Import error handling utilities
from ..utils.exceptions import (DataCompatibilityError,  # noqa: E402
                                DataNotFoundError, ParameterError,
                                ProcessingError)
# Import standardized image utilities
from ..utils.image_utils import optimize_fig_to_image_with_cache  # noqa: E402
# Import path utilities for safe file operations
from ..utils.path_utils import (get_output_dir_from_config,  # noqa: E402
                                get_safe_output_path)

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
        params: VisualizationParameters object with GridSpec spacing parameters:
                - subplot_wspace: horizontal spacing (default 0.0)
                - subplot_hspace: vertical spacing (default 0.3)
        default_title: Default title for the figure if not provided in params.
        use_tight_layout: If True, skip gridspec_kw and use tight_layout (legacy mode).
                         For spatial plots with equal aspect ratio, gridspec_kw is recommended.

    Returns:
        A tuple of (matplotlib.Figure, flattened numpy.ndarray of Axes).

    Notes:
        - For spatial plots with ax.set_aspect('equal'), gridspec_kw is the ONLY way
          to control spacing. subplots_adjust() is ineffective due to aspect constraints.
        - Default spacing (wspace=0.0, hspace=0.3) optimized for spatial plots with colorbars.
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

    # Use GridSpec for proper spacing control with equal aspect ratio
    # This is critical for spatial plots where ax.set_aspect('equal') is used
    if not use_tight_layout:
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=figsize,
            dpi=params.dpi,
            squeeze=False,
            gridspec_kw={
                "wspace": params.subplot_wspace,  # Horizontal spacing
                "hspace": params.subplot_hspace,  # Vertical spacing
            },
        )
    else:
        # Legacy mode: use tight_layout (not recommended for spatial plots)
        fig, axes = plt.subplots(
            n_rows, n_cols, figsize=figsize, dpi=params.dpi, squeeze=False
        )

    axes = axes.flatten()

    title = params.title or default_title
    fig.suptitle(title, fontsize=16)

    # Hide empty axes from the start
    for i in range(n_panels, len(axes)):
        axes[i].axis("off")

    return fig, axes


# =============================================================================
# Unified Deconvolution Data Retrieval
# =============================================================================


@dataclass
class DeconvolutionData:
    """Unified representation of deconvolution results.

    This dataclass provides a single, consistent interface for accessing
    deconvolution data regardless of the method used.

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


async def get_deconvolution_data(
    adata: ad.AnnData,
    method: Optional[str] = None,
    context: Optional[Context] = None,
) -> DeconvolutionData:
    """
    Unified function to retrieve deconvolution results from AnnData.

    This function consolidates all deconvolution data retrieval logic into
    a single, consistent interface. It handles:
    - Auto-detection when only one result exists
    - Explicit method specification
    - Clear error messages with solutions

    Args:
        adata: AnnData object with deconvolution results
        method: Deconvolution method name (e.g., "cell2location", "rctd").
                If None and only one result exists, auto-selects it.
                If None and multiple results exist, raises ValueError.
        context: MCP context for logging

    Returns:
        DeconvolutionData object with proportions and metadata

    Raises:
        DataNotFoundError: No deconvolution results found
        ValueError: Multiple results found but method not specified

    Example:
        >>> data = await get_deconvolution_data(adata, method="cell2location")
        >>> print(data.proportions.head())
        >>> print(data.cell_types)
    """
    # 1. Find all deconvolution results in obsm
    deconv_keys = [key for key in adata.obsm.keys() if key.startswith("deconvolution_")]

    # 2. Handle method specification
    if method is not None:
        # User specified method explicitly
        target_key = f"deconvolution_{method}"
        if target_key not in adata.obsm:
            available = [k.replace("deconvolution_", "") for k in deconv_keys]
            raise DataNotFoundError(
                f"Deconvolution results for method '{method}' not found.\n\n"
                f"Available methods: {available if available else 'None'}\n\n"
                "SOLUTION: Run deconvolution first or specify an available method."
            )
        proportions_key = target_key
    else:
        # Auto-detect
        if not deconv_keys:
            raise DataNotFoundError(
                "No deconvolution results found in adata.obsm.\n\n"
                "SOLUTION: Run deconvolution first:\n"
                '  deconvolve_data(data_id="your_data", '
                'params={"method": "cell2location", "cell_type_key": "..."})\n\n'
                "Available methods: cell2location, rctd, destvi, stereoscope, "
                "spotlight, card, tangram"
            )

        if len(deconv_keys) > 1:
            available = [k.replace("deconvolution_", "") for k in deconv_keys]
            raise ValueError(
                f"Multiple deconvolution results found: {available}\n\n"
                f"SOLUTION: Specify which method to visualize:\n"
                f"  params={{'deconv_method': '{available[0]}'}}\n\n"
                f"NOTE: Different methods have different assumptions.\n"
                f"      Always explicitly specify which result to visualize."
            )

        # Single result - auto-select
        proportions_key = deconv_keys[0]
        method = proportions_key.replace("deconvolution_", "")

        if context:
            await context.info(f"Auto-selected deconvolution method: {method}")

    # 3. Get cell type names
    cell_types_key = f"{proportions_key}_cell_types"
    if cell_types_key in adata.uns:
        cell_types = list(adata.uns[cell_types_key])
    else:
        # Fallback: generate generic names from shape
        n_cell_types = adata.obsm[proportions_key].shape[1]
        cell_types = [f"CellType_{i}" for i in range(n_cell_types)]
        if context:
            await context.warning(
                f"Cell type names not found in adata.uns['{cell_types_key}']. "
                f"Using generic names."
            )

    # 4. Create DataFrame
    proportions = pd.DataFrame(
        adata.obsm[proportions_key], index=adata.obs_names, columns=cell_types
    )

    # 5. Check if dominant type annotation exists
    dominant_type_key = f"dominant_celltype_{method}"
    if dominant_type_key not in adata.obs.columns:
        dominant_type_key = None

    return DeconvolutionData(
        proportions=proportions,
        method=method,
        cell_types=cell_types,
        proportions_key=proportions_key,
        dominant_type_key=dominant_type_key,
    )


# =============================================================================
# Cell Communication Data Retrieval (Unified Interface)
# =============================================================================


@dataclass
class CellCommunicationData:
    """Unified representation of cell communication analysis results.

    This dataclass provides a single, consistent interface for accessing
    cell communication data regardless of the analysis method used.

    Attributes:
        results: Main results DataFrame (format varies by method)
        method: Analysis method name ("liana_cluster", "liana_spatial", "cellphonedb")
        analysis_type: Type of analysis ("cluster" or "spatial")
        lr_pairs: List of ligand-receptor pair names

        # For spatial analysis
        spatial_scores: Spatial communication scores array (n_spots x n_pairs)
        spatial_pvals: P-values for spatial scores (optional)

        # For cluster analysis
        source_labels: List of source cell type labels
        target_labels: List of target cell type labels

        # Storage reference
        results_key: Key in adata.uns where results are stored
    """

    results: pd.DataFrame
    method: str
    analysis_type: str  # "cluster" or "spatial"
    lr_pairs: List[str]

    # Spatial analysis specific
    spatial_scores: Optional[np.ndarray] = None
    spatial_pvals: Optional[np.ndarray] = None

    # Cluster analysis specific
    source_labels: Optional[List[str]] = None
    target_labels: Optional[List[str]] = None

    # Storage key reference
    results_key: str = ""


async def get_cell_communication_data(
    adata: ad.AnnData,
    method: Optional[str] = None,
    context: Optional[Context] = None,
) -> CellCommunicationData:
    """
    Unified function to retrieve cell communication results from AnnData.

    This function consolidates all cell communication data retrieval logic into
    a single, consistent interface. It handles:
    - LIANA+ spatial bivariate analysis results
    - LIANA+ cluster-based analysis results
    - CellPhoneDB analysis results
    - Auto-detection when results exist
    - Clear error messages with solutions

    Args:
        adata: AnnData object with cell communication results
        method: Analysis method hint (optional, for future multi-result scenarios)
        context: MCP context for logging

    Returns:
        CellCommunicationData object with results and metadata

    Raises:
        DataNotFoundError: No cell communication results found

    Example:
        >>> data = await get_cell_communication_data(adata)
        >>> print(f"Found {data.method} results with {len(data.lr_pairs)} LR pairs")
        >>> if data.analysis_type == "spatial":
        ...     print(f"Spatial scores shape: {data.spatial_scores.shape}")
    """
    # 1. Check for LIANA+ spatial bivariate results (highest priority)
    if "liana_spatial_scores" in adata.obsm:
        spatial_scores = adata.obsm["liana_spatial_scores"]
        lr_pairs = adata.uns.get("liana_spatial_interactions", [])
        results_df = adata.uns.get("liana_spatial_res", pd.DataFrame())

        # Ensure results is a DataFrame
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

    # 2. Check for LIANA+ cluster-based results
    if "liana_res" in adata.uns:
        results = adata.uns["liana_res"]
        if isinstance(results, pd.DataFrame) and len(results) > 0:
            # Extract unique LR pairs
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

            # Extract source/target labels
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

    # 3. Check for CellPhoneDB results
    if "cellphonedb_means" in adata.uns:
        means = adata.uns["cellphonedb_means"]
        if isinstance(means, pd.DataFrame):
            # CellPhoneDB format: rows are LR pairs, columns are cell type pairs
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

    # 4. No results found - provide helpful error
    raise DataNotFoundError(
        "No cell communication results found in dataset.\n\n"
        "SOLUTIONS:\n"
        "1. Run cell communication analysis first:\n"
        "   analyze_cell_communication(\n"
        '       data_id="your_data_id",\n'
        "       params={\n"
        '           "species": "human",  # or "mouse"\n'
        '           "cell_type_key": "leiden",  # your cluster column\n'
        "       }\n"
        "   )\n\n"
        "2. Ensure analysis completed successfully\n\n"
        "Available methods: liana, cellphonedb, cellchat_r"
    )


def _ensure_enrichmap_compatibility(adata: ad.AnnData) -> None:
    """
    Ensure data has required metadata structure for EnrichMap visualization.

    EnrichMap and squidpy require:
    1. adata.obs['library_id'] - sample identifier column
    2. adata.uns['spatial'] - spatial metadata dictionary

    This function adds minimal metadata for single-sample data without these structures.
    """
    # Add library_id column if missing (single-sample case)
    if "library_id" not in adata.obs.columns:
        adata.obs["library_id"] = "sample_1"

    # Build spatial metadata structure if missing
    if "spatial" not in adata.uns:
        library_ids = adata.obs["library_id"].unique()
        adata.uns["spatial"] = {}
        for lib_id in library_ids:
            adata.uns["spatial"][lib_id] = {
                "images": {},  # No tissue image
                "scalefactors": {
                    "spot_diameter_fullres": 1.0,
                    "tissue_hires_scalef": 1.0,
                    "fiducial_diameter_fullres": 1.0,
                    "tissue_lowres_scalef": 1.0,
                },
            }


async def get_validated_features(
    adata: ad.AnnData,
    params: VisualizationParameters,
    max_features: int = 12,
    allow_empty: bool = False,
    default_genes: Optional[List[str]] = None,
    context: Optional[Context] = None,
) -> List[str]:
    """Gets a validated list of features (genes) for visualization.

    Args:
        adata: AnnData object
        params: Visualization parameters
        max_features: Maximum number of features allowed
        allow_empty: If True, allow params.feature to be None/empty and use default_genes
        default_genes: Default genes to use if params.feature is empty (requires allow_empty=True)
        context: MCP context for logging

    Returns:
        List of validated gene names

    Raises:
        DataNotFoundError: If genes not found or no genes specified when required
        ValueError: If too many genes specified or invalid configuration
    """
    # Ensure unique var_names before proceeding
    if not adata.var_names.is_unique:
        if context:
            await context.info("Making gene names unique to avoid indexing errors.")
        adata.var_names_make_unique()

    # Check if user specified any features
    if not params.feature:
        if allow_empty and default_genes:
            # Use default genes
            if context:
                await context.info(
                    f"No genes specified, using {len(default_genes)} default genes: "
                    f"{default_genes[:3]}{'...' if len(default_genes) > 3 else ''}"
                )
            return default_genes[:max_features]
        elif allow_empty and not default_genes:
            raise ValueError(
                "No default genes available and no features specified.\n\n"
                "INTERNAL ERROR: allow_empty=True but default_genes=None.\n"
                "This is a configuration error in the visualization function."
            )
        else:
            # Require explicit feature specification
            raise DataNotFoundError(
                "No genes specified for visualization.\n\n"
                "SOLUTIONS:\n"
                "1. Specify genes: visualize_data(data_id, params={'feature': ['CD3D', 'CD8A']})\n"
                "2. Find markers: find_markers(data_id, group_key='cell_type')\n"
                "3. Run preprocessing: preprocess_data(data_id, params={'n_top_genes': 2000})\n\n"
                "CONTEXT: Explicit gene selection required for scientific accuracy."
            )

    # Parse user features
    features = params.feature if isinstance(params.feature, list) else [params.feature]

    # Check which features are available
    available_features = [f for f in features if f in adata.var_names]
    missing_features = [f for f in features if f not in adata.var_names]

    if not available_features:
        # All features missing
        examples = list(adata.var_names[:10])
        raise DataNotFoundError(
            f"Genes not found: {missing_features}\n\n"
            f"SOLUTIONS:\n"
            f"1. Check names (available: {examples})\n"
            f"2. Search for similar gene names\n"
            f"3. Use gene discovery tools\n\n"
            f"CONTEXT: Dataset has {adata.n_vars:,} genes."
        )

    # Warn about missing features
    if missing_features and context:
        await context.warning(
            f"WARNING: {len(missing_features)} gene(s) not found and will be skipped: {missing_features}\n"
            f"Proceeding with {len(available_features)} available gene(s): {available_features}"
        )

    # Limit to max_features
    if len(available_features) > max_features:
        if context:
            await context.warning(
                f"Too many genes ({len(available_features)} > {max_features}), "
                f"limiting to first {max_features} genes"
            )
        available_features = available_features[:max_features]

    return available_features


def plot_spatial_feature(
    adata: ad.AnnData,
    feature: Optional[str],
    ax: plt.Axes,
    params: VisualizationParameters,
    force_no_colorbar: bool = False,
):
    """Plots a feature on spatial coordinates, handling background images.

    Enhanced version that supports:
    - Automatic detection of categorical vs continuous data
    - Color palette for categorical data
    - vmin/vmax for color range control
    - Tissue image handling

    Note: This function sets the title to an empty string. Callers should
    manually set ax.set_title() after calling this function if a specific
    title is desired.

    Args:
        adata: AnnData object
        feature: Feature to plot (gene name or obs column)
        ax: Matplotlib axes
        params: Visualization parameters
        force_no_colorbar: If True, disable colorbar even if params.show_colorbar is True.
                          Used when caller wants to manually add colorbar with make_axes_locatable.

    Raises:
        DataNotFoundError: If spatial coordinates are not available in adata.obsm
    """
    # Check for spatial coordinates first
    if "spatial" not in adata.obsm:
        raise DataNotFoundError(
            "Spatial coordinates not found in adata.obsm['spatial']. "
            "Spatial visualization requires 2D coordinates.\n\n"
            "Solutions:\n"
            "1. For single-cell data without spatial info, use plot_type='umap' or 'heatmap'\n"
            "2. Ensure your data has spatial coordinates in obsm['spatial']\n"
            "3. For Visium data, load with the spatial folder containing tissue_positions_list.csv"
        )

    # CRITICAL FIX: Check for tissue images correctly
    # Structure is: adata.uns["spatial"][library_id]["images"]
    # NOT: adata.uns["spatial"]["images"]
    # Must also verify that the images dict is not empty and contains actual image data
    has_image = False
    if "spatial" in adata.uns and isinstance(adata.uns["spatial"], dict):
        # Check if any library has images
        for lib_id in adata.uns["spatial"].keys():
            if isinstance(adata.uns["spatial"][lib_id], dict):
                images_dict = adata.uns["spatial"][lib_id].get("images", {})
                # Check that images dict is not empty and has hires or lowres
                if images_dict and ("hires" in images_dict or "lowres" in images_dict):
                    has_image = True
                    break

    # Detect if feature is categorical or continuous
    is_categorical = False
    if feature and feature in adata.obs.columns:
        # Ensure categorical colors are set for obs columns
        ensure_categorical(adata, feature)
        is_categorical = pd.api.types.is_categorical_dtype(adata.obs[feature])

    # Base kwargs for both functions
    plot_kwargs = {
        "color": feature,
        "ax": ax,
        "show": False,
        "alpha": params.alpha,
        "frameon": params.show_axes,
        "colorbar_loc": (
            None if force_no_colorbar else ("right" if params.show_colorbar else None)
        ),
        "title": "",  # We will set the title manually
    }

    # Add color mapping (cmap for continuous, palette for categorical)
    if is_categorical:
        # Use palette for categorical data
        plot_kwargs["palette"] = params.colormap or "Set2"
    else:
        # Use cmap for continuous data (genes, scores, etc.)
        plot_kwargs["cmap"] = params.colormap or "viridis"

    # Add vmin/vmax for continuous data color range control
    if not is_categorical:
        if params.vmin is not None:
            plot_kwargs["vmin"] = params.vmin
        if params.vmax is not None:
            plot_kwargs["vmax"] = params.vmax

    # For spatial plots with tissue image, use 'spot_size'; for embedding use 'size'
    if has_image and params.show_tissue_image:
        # Show tissue image in background
        # Add alpha_img for tissue background dimming (only supported by sc.pl.spatial)
        plot_kwargs["alpha_img"] = params.alpha_img
        # Only add spot_size if explicitly set (None = auto-determined, recommended)
        if params.spot_size is not None:
            plot_kwargs["spot_size"] = params.spot_size

        # CRITICAL FIX: Add library_id for proper spot rendering
        # Get sample key from adata.uns["spatial"]
        sample_key = None
        if "spatial" in adata.uns and isinstance(adata.uns["spatial"], dict):
            keys = list(adata.uns["spatial"].keys())
            if keys:
                sample_key = keys[0]

        # Add library_id for multi-sample reliability and proper spot rendering
        if sample_key:
            plot_kwargs["library_id"] = sample_key

        # Determine which image key to use (prefer hires, fallback to lowres)
        img_key = "hires"
        if sample_key and sample_key in adata.uns.get("spatial", {}):
            images_dict = adata.uns["spatial"][sample_key].get("images", {})
            if "hires" not in images_dict and "lowres" in images_dict:
                img_key = "lowres"

        sc.pl.spatial(adata, img_key=img_key, **plot_kwargs)
    else:
        # No image available, or user chose not to show tissue image
        # For embeddings, use 'size' parameter if spot_size is set
        if params.spot_size is not None:
            plot_kwargs["size"] = params.spot_size
        sc.pl.embedding(adata, basis="spatial", **plot_kwargs)
        ax.set_aspect("equal")

    # Note: We don't set the title here to allow callers to set their own titles
    # The old code: if params.add_gene_labels and feature: ax.set_title(feature, fontsize=12)
    # is removed so callers can manually set appropriate titles


# Helper function to validate and prepare feature for visualization
async def validate_and_prepare_feature(
    adata,
    feature: Optional[str],
    context: Optional[Context] = None,
    default_feature: Optional[str] = None,
) -> Optional[str]:
    """Validate and prepare a feature for visualization.

    This function validates that a feature exists in the dataset as either:
    - A gene name in adata.var_names
    - An observation column in adata.obs

    For deconvolution visualization, use plot_type='deconvolution' with
    the deconv_method parameter instead of passing deconvolution features here.

    Args:
        adata: AnnData object
        feature: Feature name (gene or obs column) or None
        context: MCP context for logging
        default_feature: Default feature to use if feature is None

    Returns:
        Validated feature name or default_feature or None

    Raises:
        DataNotFoundError: If feature is not found in var_names or obs columns
    """
    if not feature:
        return default_feature

    # Check if it's a valid gene
    if feature in adata.var_names:
        return feature

    # Check if it's a valid obs column
    if feature in adata.obs.columns:
        return feature

    # Feature not found - provide helpful error message
    available_genes = list(adata.var_names[:10])
    available_obs = [
        col
        for col in adata.obs.columns
        if adata.obs[col].dtype.name in ["object", "category"]
    ][:10]

    raise DataNotFoundError(
        f"Feature '{feature}' not found in dataset.\n\n"
        f"AVAILABLE DATA:\n"
        f"  - Genes (example): {available_genes}\n"
        f"  - Annotations (example): {available_obs}\n\n"
        "SOLUTIONS:\n"
        "1. Check spelling and capitalization of feature name\n"
        "2. Use an available gene or annotation column\n"
        "3. For deconvolution visualization, use:\n"
        "   plot_type='deconvolution', deconv_method='cell2location'"
    )


# ============================================================
# Private Visualization Functions (Extracted from visualize_data)
# ============================================================
# These functions handle specific plot types and follow a unified signature:
#   async def _create_<plot_type>_visualization(adata, params, context) -> plt.Figure
# ============================================================


async def _create_violin_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None,
) -> plt.Figure:
    """
    Create violin plot visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure object

    Raises:
        DataNotFoundError: If genes not found
        ValueError: If cluster_key not specified or not found
    """
    if context:
        await context.info(
            f"Creating violin plot for {params.feature if params.feature else 'top genes'}"
        )

    # Get genes for violin plot using unified validation function
    # Default: use first 3 genes if not specified
    default_genes = list(adata.var_names[:3])

    genes = await get_validated_features(
        adata,
        params,
        max_features=20,  # Reasonable limit for violin plots
        allow_empty=True,
        default_genes=default_genes,
        context=context,
    )

    # REQUIRE explicit cluster_key specification for grouping
    if not params.cluster_key:
        categorical_cols = [
            col
            for col in adata.obs.columns
            if adata.obs[col].dtype.name in ["object", "category"]
        ]
        raise ValueError(
            "Violin plot requires 'cluster_key' parameter.\n\n"
            f"Available categorical columns ({len(categorical_cols)} total):\n"
            f"  {', '.join(categorical_cols[:15])}\n\n"
            "SOLUTION: Specify cluster_key explicitly:\n"
            "  visualize_data(data_id, params={'plot_type': 'violin', 'cluster_key': 'your_column_name'})\n\n"
            "NOTE: ChatSpatial uses 'cluster_key' (not 'groupby' as in Scanpy).\n"
            "   This maintains consistency with Squidpy spatial analysis functions."
        )

    if params.cluster_key not in adata.obs.columns:
        raise ValueError(
            f"Cluster key '{params.cluster_key}' not found in data.\n\n"
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

    # Create violin plot
    # Use user's figure size if specified, otherwise default to (8, 6)
    figsize = params.figure_size if params.figure_size else (8, 6)
    plt.figure(figsize=figsize, dpi=params.dpi)
    sc.pl.violin(
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

    return fig


async def _create_dotplot_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None,
) -> plt.Figure:
    """
    Create dotplot visualization for marker gene expression.

    Dotplot displays:
    - Dot size: Fraction of cells expressing the gene (percent expressed)
    - Dot color: Mean expression level in expressing cells

    This is one of the most common visualizations for showing marker genes
    across different cell types or clusters in single-cell/spatial analysis.

    Args:
        adata: AnnData object with gene expression data
        params: Visualization parameters including:
            - feature: Gene(s) to visualize (required, or uses HVGs)
            - cluster_key: Column for grouping cells (required)
            - colormap: Color scheme (default: 'Reds')
            - dotplot_dendrogram: Show gene clustering dendrogram
            - dotplot_swap_axes: Swap genes and groups axes
            - dotplot_standard_scale: Standardize by 'var' or 'group'
            - dotplot_var_groups: Dict for grouping genes by category
        context: MCP context for logging

    Returns:
        matplotlib Figure object

    Raises:
        DataNotFoundError: If specified genes not found
        ValueError: If cluster_key not specified or not found

    Examples:
        # Basic dotplot
        params = {"plot_type": "dotplot", "feature": ["CD3D", "CD4", "CD8A"],
                  "cluster_key": "cell_type"}

        # With gene grouping
        params = {"plot_type": "dotplot", "cluster_key": "leiden",
                  "dotplot_var_groups": {"T cells": ["CD3D", "CD4"],
                                         "B cells": ["CD19", "MS4A1"]}}
    """
    if context:
        await context.info(
            f"Creating dotplot for "
            f"{params.feature if params.feature else 'marker genes'}"
        )

    # Determine genes to plot
    # If dotplot_var_groups is provided, use that structure
    if params.dotplot_var_groups:
        var_names = params.dotplot_var_groups
        # Validate all genes exist
        all_genes = []
        for group_genes in var_names.values():
            all_genes.extend(group_genes)
        missing_genes = [g for g in all_genes if g not in adata.var_names]
        if missing_genes:
            if context:
                await context.warning(
                    f"Some genes not found and will be skipped: "
                    f"{missing_genes[:10]}"
                )
            # Filter out missing genes from each group
            var_names = {
                group: [g for g in genes if g in adata.var_names]
                for group, genes in var_names.items()
            }
            # Remove empty groups
            var_names = {k: v for k, v in var_names.items() if v}
            if not var_names:
                raise DataNotFoundError(
                    "No valid genes found in dotplot_var_groups. "
                    "Please check gene names."
                )
    else:
        # Use feature parameter or default to HVGs
        default_genes = None
        if "highly_variable" in adata.var.columns:
            hvg_genes = adata.var_names[adata.var["highly_variable"]].tolist()
            default_genes = hvg_genes[:20] if hvg_genes else None

        if default_genes is None:
            default_genes = list(adata.var_names[:20])

        var_names = await get_validated_features(
            adata,
            params,
            max_features=50,  # Reasonable limit for dotplot
            allow_empty=True,
            default_genes=default_genes,
            context=context,
        )

    # REQUIRE explicit cluster_key specification for grouping
    if not params.cluster_key:
        categorical_cols = [
            col
            for col in adata.obs.columns
            if adata.obs[col].dtype.name in ["object", "category"]
        ]
        raise ValueError(
            "Dotplot requires 'cluster_key' parameter for grouping.\n\n"
            f"Available categorical columns ({len(categorical_cols)} total):\n"
            f"  {', '.join(categorical_cols[:15])}\n\n"
            "SOLUTION: Specify cluster_key explicitly:\n"
            "  visualize_data(data_id, params={'plot_type': 'dotplot', "
            "'feature': ['gene1', 'gene2'], 'cluster_key': 'your_column'})\n\n"
            "Example: cluster_key='leiden' or cluster_key='cell_type'"
        )

    if params.cluster_key not in adata.obs.columns:
        raise ValueError(
            f"Cluster key '{params.cluster_key}' not found in data.\n\n"
            f"Available columns: {', '.join(list(adata.obs.columns)[:20])}"
        )

    groupby = params.cluster_key

    # Ensure groupby column is categorical
    if adata.obs[groupby].dtype.name != "category":
        adata.obs[groupby] = adata.obs[groupby].astype("category")

    # Limit number of groups to avoid oversized plots
    n_groups = len(adata.obs[groupby].cat.categories)
    max_groups = 15
    categories_order = params.dotplot_categories_order

    if n_groups > max_groups and categories_order is None:
        if context:
            await context.warning(
                f"Too many groups ({n_groups}). Limiting to {max_groups} "
                f"largest groups."
            )
        # Get the largest groups
        group_counts = adata.obs[groupby].value_counts().nlargest(max_groups)
        categories_order = group_counts.index.tolist()

    # Determine colormap - default to 'Reds' for dotplot
    cmap = params.colormap if params.colormap != "coolwarm" else "Reds"

    # Calculate figure size based on number of genes and groups
    if params.figure_size:
        figsize = params.figure_size
    else:
        n_genes = (
            len(var_names)
            if isinstance(var_names, list)
            else sum(len(v) for v in var_names.values())
        )
        n_display_groups = (
            len(categories_order) if categories_order else min(n_groups, max_groups)
        )
        # Width based on groups, height based on genes
        width = max(4, min(12, n_display_groups * 0.8 + 2))
        height = max(4, min(16, n_genes * 0.35 + 2))
        figsize = (width, height)

    # Create dotplot
    plt.figure(figsize=figsize, dpi=params.dpi)

    # Build kwargs for sc.pl.dotplot
    dotplot_kwargs = {
        "var_names": var_names,
        "groupby": groupby,
        "show": False,
        "cmap": cmap,
        "dendrogram": params.dotplot_dendrogram,
        "swap_axes": params.dotplot_swap_axes,
    }

    # Add optional parameters if specified
    if params.dotplot_standard_scale:
        dotplot_kwargs["standard_scale"] = params.dotplot_standard_scale

    if params.dotplot_dot_max is not None:
        dotplot_kwargs["dot_max"] = params.dotplot_dot_max

    if params.dotplot_dot_min is not None:
        dotplot_kwargs["dot_min"] = params.dotplot_dot_min

    if params.dotplot_smallest_dot > 0:
        dotplot_kwargs["smallest_dot"] = params.dotplot_smallest_dot

    if categories_order:
        dotplot_kwargs["categories_order"] = categories_order

    if params.vmin is not None:
        dotplot_kwargs["vmin"] = params.vmin

    if params.vmax is not None:
        dotplot_kwargs["vmax"] = params.vmax

    # Set title
    if params.title:
        dotplot_kwargs["title"] = params.title

    # Create the dotplot
    sc.pl.dotplot(adata, **dotplot_kwargs)

    # Get the figure
    fig = plt.gcf()

    # Adjust layout - use subplots_adjust instead of tight_layout
    # because scanpy dotplot creates complex axes that aren't compatible
    # with tight_layout
    try:
        plt.tight_layout()
    except Exception:
        # Fallback for complex layouts (e.g., with dendrogram)
        fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.15)

    return fig


async def _create_card_imputation_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None,
) -> plt.Figure:
    """
    Create CARD imputation visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure object

    Raises:
        DataNotFoundError: If CARD imputation data not found or feature not found
    """
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
            Patch(facecolor=color_map[ct], label=ct) for ct in sorted(unique_types)
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

    return fig


async def _create_spatial_cnv_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None,
) -> plt.Figure:
    """
    Create spatial CNV projection visualization (REFACTORED to use plot_spatial_feature helper)

    This function now uses the unified plot_spatial_feature() helper for cleaner
    code and consistent parameter handling.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure object

    Raises:
        DataNotFoundError: If spatial coordinates or CNV features not found
    """
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
        await context.info(f"Visualizing {feature_to_plot} on spatial coordinates")

    # REFACTORED: Use unified plot_spatial_feature() helper
    # Override colormap default for CNV data (RdBu_r is better for CNV scores)
    if not params.colormap:
        # Only override if user hasn't specified a colormap
        params.colormap = (
            "RdBu_r"
            if not pd.api.types.is_categorical_dtype(adata.obs[feature_to_plot])
            else "tab20"
        )

    figsize = params.figure_size if params.figure_size else (10, 8)
    fig, ax = plt.subplots(figsize=figsize)

    # Use the enhanced plot_spatial_feature helper
    # It automatically handles: tissue images, categorical/continuous data, palette/cmap, vmin/vmax, spot_size
    plot_spatial_feature(adata, feature_to_plot, ax, params)

    if context:
        await context.info(f"Spatial CNV projection created for {feature_to_plot}")

    return fig


async def _create_cnv_heatmap_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None,
) -> plt.Figure:
    """
    Create CNV heatmap visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure object

    Raises:
        DataNotFoundError: If CNV data not found
        DataCompatibilityError: If infercnvpy not installed
    """
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
    require("infercnvpy", feature="CNV heatmap visualization")

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
            plt.colorbar(im, ax=ax, label="Mean CNV state")

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
        import infercnvpy as cnv

        if context:
            await context.info("Creating chromosome-organized CNV heatmap...")

        cnv.pl.chromosome_heatmap(
            adata,
            groupby=params.cluster_key,
            dendrogram=True,
            show=False,
            figsize=figsize,
        )
        # Get current figure
        fig = plt.gcf()

    if context:
        await context.info("CNV heatmap created successfully")

    return fig


# ============================================================
# Plot Type Handlers - Dispatch Dictionary
# ============================================================
# Maps plot_type string to corresponding visualization function
# All handlers follow the same signature:
#   async def handler(adata, params, context) -> plt.Figure
#
# NOTE: This is a partial implementation. Some plot types (umap, spatial, heatmap)
#       still use inline logic in visualize_data and will be migrated in future.
# ============================================================


async def _create_umap_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None,
) -> plt.Figure:
    """
    Create UMAP visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure object
    """
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
                "UMAP visualization requires neighborhood graph.\n\n"
                "SOLUTION:\n"
                "Run preprocessing first:\n"
                "1. sc.pp.neighbors(adata)\n"
                "2. sc.tl.umap(adata)\n\n"
                "SCIENTIFIC INTEGRITY: UMAP requires a k-nearest neighbor graph "
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
                adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

            # Compute UMAP (scanpy already imported globally)
            sc.tl.umap(adata)

            if context:
                await context.info("Successfully computed UMAP coordinates")

        except Exception as e:
            # NO FALLBACK: Honest error reporting with alternatives
            error_msg = (
                "Failed to compute UMAP for visualization.\n\n"
                f"Error: {str(e)}\n\n"
                "ALTERNATIVES:\n"
                "1. Use PCA visualization instead (LINEAR method):\n"
                "   visualize_data(data_id, params={'plot_type': 'pca'})\n\n"
                "2. Use t-SNE visualization (NON-LINEAR, different from UMAP):\n"
                "   visualize_data(data_id, params={'plot_type': 'tsne'})\n\n"
                "3. Fix UMAP computation in preprocessing:\n"
                "   - Ensure data is properly normalized\n"
                "   - Check for extreme values or NaN\n"
                "   - Verify neighbors graph is computed correctly\n"
                "   - Try different UMAP parameters\n\n"
                "SCIENTIFIC INTEGRITY: UMAP (Uniform Manifold Approximation) and PCA "
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
        fig = await _create_multi_gene_umap_visualization(adata, params, context)
    else:
        # REQUIRE explicit feature specification - no automatic defaults
        if not feature_list:
            categorical_cols = [
                col
                for col in adata.obs.columns
                if adata.obs[col].dtype.name in ["object", "category"]
            ]
            raise ValueError(
                "UMAP visualization requires 'feature' parameter.\n\n"
                f"Available categorical columns ({len(categorical_cols)} total):\n"
                f"  {', '.join(categorical_cols[:15])}\n\n"
                "SOLUTION: Specify feature explicitly:\n"
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
                ensure_categorical(adata, feature)

        # Add size encoding if requested
        if params.size_by:
            if params.size_by in adata.var_names or params.size_by in adata.obs.columns:
                # For scanpy compatibility, we need to pass size values correctly
                if params.size_by in adata.var_names:
                    # Gene expression - get the values
                    size_values = adata[:, params.size_by].X
                    if hasattr(size_values, "toarray"):
                        size_values = size_values.toarray().flatten()
                    else:
                        size_values = size_values.flatten()
                    # Normalize to reasonable size range (10-100)
                    size_values = 10 + 90 * (size_values - size_values.min()) / (
                        size_values.max() - size_values.min() + 1e-8
                    )
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
                            val: 20 + i * 15 for i, val in enumerate(unique_vals)
                        }
                        plot_kwargs["size"] = (
                            adata.obs[params.size_by].map(size_map).values
                        )
                    else:
                        # Numeric - normalize to size range
                        size_values = adata.obs[params.size_by].values
                        size_values = 10 + 90 * (size_values - size_values.min()) / (
                            size_values.max() - size_values.min() + 1e-8
                        )
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

        # Add velocity overlay if requested
        if params.show_velocity:
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
            except Exception as e:
                if context:
                    await context.warning(f"Failed to add overlays: {str(e)}")

    return fig


async def _create_spatial_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None,
) -> plt.Figure:
    """
    Create spatial visualization (REFACTORED to use plot_spatial_feature helper)

    This function now uses the unified plot_spatial_feature() helper for cleaner
    code and consistent parameter handling across all spatial visualizations.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure object
    """
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
        fig = await _create_multi_gene_visualization(adata, params, context)
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
                f"Feature '{single_feature}' not found in dataset.\n\n"
                f"SOLUTIONS:\n"
                f"1. Check spelling and case (gene names are case-sensitive)\n"
                f"   - Mouse genes: First letter uppercase (e.g., 'Cd5l', 'Gbp2b')\n"
                f"   - Human genes: All uppercase (e.g., 'CD5L', 'GBP2B')\n"
                f"2. Available genes (first 10): {examples}\n"
                f"3. Available annotations (first 10): {obs_examples}\n\n"
                f"CONTEXT: Dataset has {adata.n_vars:,} genes and "
                f"{len(adata.obs.columns)} annotations."
            )

        # Create spatial plot using unified helper function
        # REFACTORED: Now uses plot_spatial_feature() for all plotting
        fig, ax = create_figure()

        # Use the enhanced plot_spatial_feature helper
        # It handles: tissue images, categorical/continuous data, palette/cmap, vmin/vmax, spot_size
        plot_spatial_feature(adata, feature, ax, params)

        # Set title if no feature specified
        if not feature:
            ax.set_title("Spatial coordinates")

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
                # Check if we have tissue images (must have actual hires or lowres)
                has_tissue_image = False
                img_key_to_use = "hires"
                if "spatial" in adata.uns:
                    for key in adata.uns["spatial"].keys():
                        images_dict = adata.uns["spatial"][key].get("images", {})
                        if images_dict and (
                            "hires" in images_dict or "lowres" in images_dict
                        ):
                            has_tissue_image = True
                            # Prefer hires, fallback to lowres
                            img_key_to_use = (
                                "hires" if "hires" in images_dict else "lowres"
                            )
                            break

                if has_tissue_image and params.show_tissue_image:
                    # Use scanpy's add_outline functionality for tissue data with image
                    sc.pl.spatial(
                        adata,
                        img_key=img_key_to_use,
                        color=params.outline_cluster_key,
                        add_outline=True,
                        outline_color=params.outline_color,
                        outline_width=params.outline_width,
                        show=False,
                        ax=ax,
                        legend_loc=None,
                        colorbar=False,
                    )
                else:
                    # For non-tissue data, create manual outlines using ConvexHull
                    spatial_coords = adata.obsm["spatial"]
                    cluster_labels = adata.obs[params.outline_cluster_key].values

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
                    await context.warning(f"Failed to add outline overlay: {str(e)}")

    return fig


async def _create_heatmap_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context: Optional[Context] = None,
) -> plt.Figure:
    """Create heatmap visualization using scanpy."""
    # Check for highly variable genes
    hvg_exists = "highly_variable" in adata.var
    hvg_any = adata.var["highly_variable"].any() if hvg_exists else False

    if not hvg_exists or not hvg_any:
        raise ValueError(
            "Heatmap requires highly variable genes but none found.\n"
            "Run preprocess_data first or specify genes via 'feature' parameter."
        )

    # Require cluster_key
    if not params.cluster_key:
        categorical_cols = [
            col
            for col in adata.obs.columns
            if adata.obs[col].dtype.name in ["object", "category"]
        ][:15]
        raise ValueError(
            f"cluster_key required. Available: {', '.join(categorical_cols)}"
        )

    if params.cluster_key not in adata.obs.columns:
        raise ValueError(f"'{params.cluster_key}' not found in data.")

    groupby = params.cluster_key

    # Limit groups to avoid oversized plots
    n_groups = len(adata.obs[groupby].cat.categories)
    if n_groups > 10:
        if context:
            await context.warning(f"Limiting to 10 largest groups (from {n_groups})")
        group_counts = adata.obs[groupby].value_counts().nlargest(10).index
        adata = adata[adata.obs[groupby].isin(group_counts)].copy()

    # Get genes for heatmap
    n_genes = 20
    default_genes = adata.var_names[adata.var.highly_variable][:n_genes].tolist()

    available_genes = await get_validated_features(
        adata,
        params,
        max_features=n_genes,
        allow_empty=True,
        default_genes=default_genes,
        context=context,
    )

    # Use scanpy's heatmap directly
    figsize = params.figure_size or (
        max(8, len(available_genes) * 0.5),
        max(6, len(adata.obs[groupby].cat.categories) * 0.3),
    )

    sc.pl.heatmap(
        adata,
        var_names=available_genes,
        groupby=groupby,
        cmap=params.colormap,
        show=False,
        dendrogram=False,
        standard_scale="var",
        figsize=figsize,
    )

    fig = plt.gcf()
    plt.tight_layout()
    return fig


async def visualize_data(
    data_id: str,
    ctx: "ToolContext",
    params: VisualizationParameters = VisualizationParameters(),
) -> Union[ImageContent, Tuple[ImageContent, EmbeddedResource]]:
    """Visualize spatial transcriptomics data

    Args:
        data_id: Dataset ID
        ctx: ToolContext for unified data access and logging
        params: Visualization parameters

    Returns:
        Union[ImageContent, Tuple[ImageContent, EmbeddedResource]]:
            - Small images (<100KB): ImageContent object
            - Large images (>=100KB): Tuple[Preview ImageContent, High-quality Resource]

    Raises:
        DataNotFoundError: If the dataset is not found
        ParameterError: If parameters are invalid
        DataCompatibilityError: If data is not compatible with the visualization
        ProcessingError: If processing fails
    """
    # Validate parameters - use PLOT_HANDLERS as single source of truth
    if params.plot_type not in PLOT_HANDLERS:
        error_msg = (
            f"Invalid plot_type: {params.plot_type}. "
            f"Must be one of {list(PLOT_HANDLERS.keys())}"
        )
        await ctx.warning(error_msg)
        raise ParameterError(error_msg)

    await ctx.info(f"Visualizing {params.plot_type} plot for dataset {data_id}")
    await ctx.info(f"Parameters: feature={params.feature}, colormap={params.colormap}")

    try:
        # Retrieve the AnnData object via ToolContext
        adata = await ctx.get_adata(data_id)

        # Validate AnnData object - basic validation
        if adata.n_obs < 5:
            raise DataNotFoundError("Dataset has too few cells (minimum 5 required)")
        if adata.n_vars < 5:
            raise DataNotFoundError("Dataset has too few genes (minimum 5 required)")

        # Set matplotlib style for better visualizations
        # Use user's DPI setting if provided, otherwise default to 100
        sc.settings.set_figure_params(dpi=params.dpi or 100, facecolor="white")

        # ============================================================
        # Dispatch Logic - Use PLOT_HANDLERS for all plot types
        # ============================================================
        if params.plot_type in PLOT_HANDLERS:
            # Use dispatch dictionary for all visualization functions
            handler = PLOT_HANDLERS[params.plot_type]
            fig = await handler(adata, params, ctx._mcp_context)
        else:
            raise ValueError(f"Unknown plot type: {params.plot_type}")

        # Convert figure with optimization (preview + resource for large images)
        await ctx.info(f"Converting {params.plot_type} figure...")

        # Generate plot_type_key with subtype if applicable (for cache consistency)
        subtype = params.subtype
        plot_type_key = f"{params.plot_type}_{subtype}" if subtype else params.plot_type

        # Use the optimized conversion function
        return await optimize_fig_to_image_with_cache(
            fig,
            params,
            ctx._mcp_context,
            data_id=data_id,
            plot_type=plot_type_key,
            mode="auto",
        )

    except Exception as e:
        # Make sure to close any open figures in case of error
        plt.close("all")

        # Log the error
        error_msg = f"Error in {params.plot_type} visualization: {str(e)}"
        await ctx.warning(error_msg)
        await ctx.info(f"Error details: {traceback.format_exc()}")

        # For image conversion errors, return error message as string
        if "fig_to_image" in str(e) or "convert" in str(e).lower():
            error_details = traceback.format_exc()
            return (
                f"Error in {params.plot_type} visualization:\n\n"
                f"{str(e)}\n\n"
                f"Technical details:\n{error_details}"
            )

        # Wrap the error in a more informative exception
        if isinstance(e, (DataNotFoundError, ParameterError, DataCompatibilityError)):
            # Re-raise specific errors
            raise
        else:
            # Wrap generic errors
            raise ProcessingError(
                f"Failed to create {params.plot_type} visualization: {str(e)}"
            ) from e


async def _create_dominant_celltype_map(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """
    Create dominant cell type map (CARD-style)

    Shows the dominant cell type at each spatial location, optionally
    marking "pure" vs "mixed" spots based on proportion threshold.

    Args:
        adata: AnnData with deconvolution results
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure

    Raises:
        DataNotFoundError: If deconvolution results not found
    """
    # Get deconvolution data using unified function
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Get dominant cell type
    dominant_idx = data.proportions.values.argmax(axis=1)
    dominant_types = data.proportions.columns[dominant_idx].values
    dominant_proportions = data.proportions.values.max(axis=1)

    # Mark pure vs mixed spots
    if params.show_mixed_spots:
        spot_categories = np.where(
            dominant_proportions >= params.min_proportion_threshold,
            dominant_types,
            "Mixed",
        )
    else:
        spot_categories = dominant_types

    # Get spatial coordinates
    if "spatial" not in adata.obsm:
        raise DataNotFoundError(
            "Spatial coordinates not found in adata.obsm['spatial']"
        )
    spatial_coords = adata.obsm["spatial"]

    # Create figure
    figsize = params.figure_size if params.figure_size else (10, 8)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Get unique categories
    unique_categories = np.unique(spot_categories)
    n_categories = len(unique_categories)

    # Create colormap
    if params.show_mixed_spots and "Mixed" in unique_categories:
        # Separate colors for Mixed
        cell_type_categories = [c for c in unique_categories if c != "Mixed"]
        n_cell_types = len(cell_type_categories)

        # Use tab20 for cell types
        if n_cell_types <= 20:
            cell_type_cmap = plt.cm.get_cmap("tab20", n_cell_types)
            cell_type_colors = {
                ct: cell_type_cmap(i) for i, ct in enumerate(cell_type_categories)
            }
        else:
            cell_type_colors = {
                ct: plt.cm.get_cmap(params.colormap or "tab20", n_cell_types)(i)
                for i, ct in enumerate(cell_type_categories)
            }

        # Gray for mixed
        cell_type_colors["Mixed"] = (0.7, 0.7, 0.7, 1.0)

        # Plot
        for category in unique_categories:
            mask = spot_categories == category
            ax.scatter(
                spatial_coords[mask, 0],
                spatial_coords[mask, 1],
                c=[cell_type_colors[category]],
                s=params.spot_size or 10,
                alpha=0.8 if category == "Mixed" else 1.0,
                label=category,
                edgecolors="none",
            )
    else:
        # No mixed spots - standard categorical plot

        if n_categories <= 20:
            cmap = plt.cm.get_cmap("tab20", n_categories)
        else:
            cmap = plt.cm.get_cmap(params.colormap or "tab20", n_categories)

        colors = {cat: cmap(i) for i, cat in enumerate(unique_categories)}

        for category in unique_categories:
            mask = spot_categories == category
            ax.scatter(
                spatial_coords[mask, 0],
                spatial_coords[mask, 1],
                c=[colors[category]],
                s=params.spot_size or 10,
                alpha=1.0,
                label=category,
                edgecolors="none",
            )

    # Formatting
    ax.set_xlabel("Spatial X")
    ax.set_ylabel("Spatial Y")
    ax.set_title(
        f"Dominant Cell Type Map ({data.method})\n"
        f"Threshold: {params.min_proportion_threshold:.2f}"
        if params.show_mixed_spots
        else f"Dominant Cell Type Map ({data.method})"
    )
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        ncol=1 if n_categories <= 15 else 2,
        fontsize=8,
        markerscale=0.5,
    )
    ax.set_aspect("equal")

    plt.tight_layout()

    if context:
        await context.info(
            f"Created dominant cell type map with {n_categories} categories "
            f"({len([c for c in unique_categories if c != 'Mixed'])} cell types"
            + (", 1 mixed category)" if "Mixed" in unique_categories else ")")
        )

    return fig


async def _create_diversity_map(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """
    Create Shannon entropy diversity map

    Shows cell type diversity at each spatial location using Shannon entropy.
    Higher entropy = more diverse/mixed cell types.
    Lower entropy = more homogeneous/dominated by single type.

    Args:
        adata: AnnData with deconvolution results
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure

    Raises:
        DataNotFoundError: If deconvolution results not found
    """
    from scipy.stats import entropy

    # Get deconvolution data using unified function
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Calculate Shannon entropy for each spot
    # Add small epsilon to avoid log(0)
    epsilon = 1e-10
    proportions_safe = data.proportions.values + epsilon

    # Shannon entropy: -sum(p * log(p))
    spot_entropy = entropy(proportions_safe.T, base=2)  # Base 2 for bits

    # Normalize to [0, 1] range
    # Max entropy = log2(n_cell_types)
    max_entropy = np.log2(data.proportions.shape[1])
    normalized_entropy = spot_entropy / max_entropy

    # Get spatial coordinates
    if "spatial" not in adata.obsm:
        raise DataNotFoundError(
            "Spatial coordinates not found in adata.obsm['spatial']"
        )
    spatial_coords = adata.obsm["spatial"]

    # Create figure
    figsize = params.figure_size if params.figure_size else (10, 8)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Plot entropy
    scatter = ax.scatter(
        spatial_coords[:, 0],
        spatial_coords[:, 1],
        c=normalized_entropy,
        cmap=params.colormap or "viridis",
        s=params.spot_size or 10,
        alpha=1.0,
        edgecolors="none",
    )

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Cell Type Diversity (Shannon Entropy)", rotation=270, labelpad=20)

    # Formatting
    ax.set_xlabel("Spatial X")
    ax.set_ylabel("Spatial Y")
    ax.set_title(
        f"Cell Type Diversity Map ({data.method})\n"
        f"Shannon Entropy (0=homogeneous, 1=maximally diverse)"
    )
    ax.set_aspect("equal")

    plt.tight_layout()

    if context:
        # Calculate statistics
        mean_entropy = normalized_entropy.mean()
        std_entropy = normalized_entropy.std()
        high_diversity_pct = (
            (normalized_entropy > 0.7).sum() / len(normalized_entropy) * 100
        )
        low_diversity_pct = (
            (normalized_entropy < 0.3).sum() / len(normalized_entropy) * 100
        )

        await context.info(
            f"Created diversity map:\n"
            f"  Mean entropy: {mean_entropy:.3f} ± {std_entropy:.3f}\n"
            f"  High diversity (>0.7): {high_diversity_pct:.1f}% of spots\n"
            f"  Low diversity (<0.3): {low_diversity_pct:.1f}% of spots"
        )

    return fig


async def _create_stacked_barplot(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """
    Create stacked barplot of cell type proportions

    Shows cell type proportions for each spot as stacked bars.
    Spots can be sorted by dominant cell type, spatial order, or cluster.

    Args:
        adata: AnnData with deconvolution results
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure

    Raises:
        DataNotFoundError: If deconvolution results not found
    """
    # Get deconvolution data using unified function
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Limit number of spots for readability
    n_spots = len(data.proportions)
    if n_spots > params.max_spots:
        # Sample spots
        sample_indices = np.random.choice(n_spots, size=params.max_spots, replace=False)
        proportions_plot = data.proportions.iloc[sample_indices]
        if context:
            await context.warning(
                f"Sampled {params.max_spots} spots out of {n_spots} for readability. "
                f"Adjust max_spots parameter to show more."
            )
    else:
        proportions_plot = data.proportions

    # Sort spots based on sort_by parameter
    if params.sort_by == "dominant_type":
        # Sort by dominant cell type
        dominant_idx = proportions_plot.values.argmax(axis=1)
        dominant_types = proportions_plot.columns[dominant_idx]
        sort_order = np.argsort(dominant_types)
    elif params.sort_by == "spatial":
        # Sort by spatial distance (hierarchical clustering on coordinates)
        if "spatial" in adata.obsm:
            from scipy.cluster.hierarchy import dendrogram, linkage

            spatial_coords = adata.obsm["spatial"][proportions_plot.index]
            linkage_matrix = linkage(spatial_coords, method="ward")
            dend = dendrogram(linkage_matrix, no_plot=True)
            sort_order = dend["leaves"]
        else:
            sort_order = np.arange(len(proportions_plot))
    elif params.sort_by == "cluster":
        # Sort by cluster if available
        cluster_key = params.cluster_key or "leiden"
        if cluster_key in adata.obs.columns:
            cluster_values = adata.obs.loc[proportions_plot.index, cluster_key]
            sort_order = np.argsort(cluster_values.astype(str))
        else:
            if context:
                await context.warning(
                    f"Cluster key '{cluster_key}' not found. Using default order."
                )
            sort_order = np.arange(len(proportions_plot))
    else:
        sort_order = np.arange(len(proportions_plot))

    proportions_sorted = proportions_plot.iloc[sort_order]

    # Create figure
    figsize = params.figure_size if params.figure_size else (12, 6)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Get cell types
    cell_types = proportions_sorted.columns.tolist()
    n_cell_types = len(cell_types)

    # Create colormap
    if n_cell_types <= 20:
        cmap = plt.cm.get_cmap("tab20", n_cell_types)
    else:
        cmap = plt.cm.get_cmap(params.colormap or "tab20", n_cell_types)
    colors = [cmap(i) for i in range(n_cell_types)]

    # Plot stacked bars
    x_positions = np.arange(len(proportions_sorted))
    bottom = np.zeros(len(proportions_sorted))

    for i, cell_type in enumerate(cell_types):
        values = proportions_sorted[cell_type].values
        ax.bar(
            x_positions,
            values,
            bottom=bottom,
            color=colors[i],
            label=cell_type,
            width=1.0,
            edgecolor="none",
        )
        bottom += values

    # Formatting
    ax.set_xlabel(
        "Spot Index"
        if params.sort_by == "spatial"
        else params.sort_by.replace("_", " ").title()
    )
    ax.set_ylabel("Cell Type Proportion")
    ax.set_title(
        f"Cell Type Proportions ({data.method})\n"
        f"Sorted by: {params.sort_by.replace('_', ' ').title()}"
    )
    ax.set_ylim([0, 1])
    ax.set_xlim([0, len(proportions_sorted)])

    # Legend
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        ncol=1 if n_cell_types <= 15 else 2,
        fontsize=8,
    )

    # Remove x-tick labels for readability (too many spots)
    ax.set_xticks([])

    plt.tight_layout()

    if context:
        await context.info(
            f"Created stacked barplot with {len(proportions_sorted)} spots and {n_cell_types} cell types"
        )

    return fig


async def _create_scatterpie_plot(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """
    Create spatial scatterpie plot (SPOTlight-style)

    Shows cell type proportions as pie charts at each spatial location.
    Each spot is represented as a pie chart showing the composition of
    different cell types.

    Args:
        adata: AnnData with deconvolution results
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure

    Raises:
        DataNotFoundError: If deconvolution or spatial data not found
    """
    from matplotlib.patches import Wedge

    # Get deconvolution data using unified function
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Get spatial coordinates
    if "spatial" not in adata.obsm:
        raise DataNotFoundError(
            "Spatial coordinates not found in adata.obsm['spatial']"
        )
    spatial_coords = adata.obsm["spatial"]

    # Use all spots (no sampling)
    proportions_plot = data.proportions
    coords_plot = spatial_coords

    # Get cell types
    cell_types = proportions_plot.columns.tolist()
    n_cell_types = len(cell_types)

    # Create colormap
    if n_cell_types <= 20:
        cmap = plt.cm.get_cmap("tab20", n_cell_types)
    else:
        cmap = plt.cm.get_cmap(params.colormap or "tab20", n_cell_types)
    colors = {cell_type: cmap(i) for i, cell_type in enumerate(cell_types)}

    # Create figure
    figsize = params.figure_size if params.figure_size else (12, 10)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Calculate pie radius based on spatial scale
    # Use params.pie_scale to adjust size
    coord_range = np.ptp(coords_plot, axis=0).max()
    # Use 2% (1/50th) to match R scatterpie package behavior
    # This ensures consistency with SPOTlight, STdeconvolve, and standard R workflows
    base_radius = coord_range * 0.02  # 2% of coordinate range (matches R scatterpie)
    pie_radius = base_radius * params.pie_scale

    # Draw pie charts at each location
    for idx in range(len(proportions_plot)):
        x, y = coords_plot[idx]
        prop_values = proportions_plot.iloc[idx].values

        # Skip if all proportions are zero
        if prop_values.sum() == 0:
            continue

        # Normalize proportions
        prop_normalized = prop_values / prop_values.sum()

        # Draw wedges
        start_angle = 0
        for cell_type, proportion in zip(cell_types, prop_normalized):
            if proportion > 0.01:  # Only draw if >1%
                angle = proportion * 360
                wedge = Wedge(
                    center=(x, y),
                    r=pie_radius,
                    theta1=start_angle,
                    theta2=start_angle + angle,
                    facecolor=colors[cell_type],
                    edgecolor="white",
                    linewidth=0.5,
                    alpha=params.scatterpie_alpha,
                )
                ax.add_patch(wedge)
                start_angle += angle

    # Set axis limits with padding
    x_min, x_max = coords_plot[:, 0].min(), coords_plot[:, 0].max()
    y_min, y_max = coords_plot[:, 1].min(), coords_plot[:, 1].max()
    padding = pie_radius * 2
    ax.set_xlim([x_min - padding, x_max + padding])
    ax.set_ylim([y_min - padding, y_max + padding])

    # Formatting
    ax.set_xlabel("Spatial X")
    ax.set_ylabel("Spatial Y")
    ax.set_title(
        f"Spatial Scatterpie Plot ({data.method})\n"
        f"Cell Type Composition (pie scale: {params.pie_scale:.2f})"
    )
    ax.set_aspect("equal")

    # Create legend
    from matplotlib.patches import Patch

    legend_elements = [Patch(facecolor=colors[ct], label=ct) for ct in cell_types]
    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        ncol=1 if n_cell_types <= 15 else 2,
        fontsize=8,
    )

    plt.tight_layout()

    if context:
        await context.info(
            f"Created scatterpie plot with {len(proportions_plot)} spots and {n_cell_types} cell types"
        )

    return fig


async def _create_umap_proportions(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """
    Create UMAP colored by cell type proportions

    Shows UMAP embeddings in multi-panel format, with each panel showing
    the proportion of a specific cell type.

    Args:
        adata: AnnData with deconvolution results and UMAP
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure

    Raises:
        DataNotFoundError: If deconvolution or UMAP not found
    """
    # Get deconvolution data using unified function
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Check for UMAP coordinates
    if "X_umap" not in adata.obsm:
        raise DataNotFoundError(
            "UMAP coordinates not found in adata.obsm['X_umap']. "
            "Run UMAP dimensionality reduction first."
        )
    umap_coords = adata.obsm["X_umap"]

    # Get cell types (limit to n_cell_types)
    cell_types = data.proportions.columns.tolist()
    n_cell_types_total = len(cell_types)

    # Select top cell types by mean proportion
    mean_proportions = data.proportions.mean(axis=0).sort_values(ascending=False)
    top_cell_types = mean_proportions.head(params.n_cell_types).index.tolist()

    # Create multi-panel figure
    n_panels = len(top_cell_types)
    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(ncols * 4, nrows * 3.5), squeeze=False
    )
    axes = axes.flatten()

    # Plot each cell type
    for idx, cell_type in enumerate(top_cell_types):
        ax = axes[idx]

        # Get proportions for this cell type
        prop_values = data.proportions[cell_type].values

        # Create scatter plot
        scatter = ax.scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c=prop_values,
            cmap=params.colormap or "viridis",
            s=params.spot_size or 5,
            alpha=0.8,
            vmin=0,
            vmax=1,
            edgecolors="none",
        )

        # Formatting
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"{cell_type}\n(mean: {mean_proportions[cell_type]:.3f})")
        ax.set_aspect("equal")

        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Proportion", rotation=270, labelpad=15, fontsize=8)

    # Hide extra axes
    for idx in range(n_panels, len(axes)):
        axes[idx].axis("off")

    # Overall title
    fig.suptitle(
        f"UMAP Cell Type Proportions ({data.method})\n"
        f"Top {n_panels} cell types (out of {n_cell_types_total})",
        fontsize=12,
        y=0.995,
    )

    plt.tight_layout()

    if context:
        await context.info(
            f"Created UMAP proportions plot with {n_panels} cell types "
            f"(showing top {n_panels}/{n_cell_types_total} by mean proportion)"
        )

    return fig


async def _create_deconvolution_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create deconvolution results visualization

    Routes to appropriate visualization based on params.subtype:
    - spatial_multi: Multi-panel spatial maps (default)
    - dominant_type: Dominant cell type map (CARD-style)
    - diversity: Shannon entropy diversity map
    - stacked_bar: Stacked barplot
    - scatterpie: Spatial scatterpie (SPOTlight-style)
    - umap: UMAP colored by proportions

    Args:
        adata: AnnData object with deconvolution results
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with deconvolution visualization
    """
    # Route to appropriate visualization
    viz_type = params.subtype

    if viz_type == "dominant_type":
        return await _create_dominant_celltype_map(adata, params, context)
    elif viz_type == "diversity":
        return await _create_diversity_map(adata, params, context)
    elif viz_type == "stacked_bar":
        return await _create_stacked_barplot(adata, params, context)
    elif viz_type == "scatterpie":
        return await _create_scatterpie_plot(adata, params, context)
    elif viz_type == "umap":
        return await _create_umap_proportions(adata, params, context)
    elif viz_type == "spatial_multi":
        # Original multi-panel spatial implementation
        return await _create_spatial_multi_deconvolution(adata, params, context)
    else:
        raise ValueError(
            f"Unknown deconvolution visualization type: {viz_type}. "
            f"Available: spatial_multi, dominant_type, diversity, stacked_bar, scatterpie, umap"
        )


async def _create_spatial_multi_deconvolution(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Multi-panel spatial deconvolution visualization.

    Shows top N cell types as separate spatial plots.

    Args:
        adata: AnnData object with deconvolution results
        params: Visualization parameters (uses params.deconv_method)
        context: MCP context

    Returns:
        Matplotlib figure with multi-panel spatial visualization
    """
    # Use unified data retrieval function
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Get top cell types by mean proportion
    n_cell_types = min(params.n_cell_types, len(data.cell_types))
    top_cell_types = (
        data.proportions.mean().sort_values(ascending=False).index[:n_cell_types]
    )

    # Setup multi-panel figure
    fig, axes = setup_multi_panel_figure(
        n_panels=len(top_cell_types),
        params=params,
        default_title=f"{data.method.upper()} Cell Type Proportions",
    )

    # Plot each cell type using temporary column
    temp_feature_key = "_deconv_viz_temp"

    for i, cell_type in enumerate(top_cell_types):
        if i < len(axes):
            ax = axes[i]
            try:
                # Get cell type proportions
                proportions_values = data.proportions[cell_type].values

                # Handle NaN values
                if pd.isna(proportions_values).any():
                    proportions_values = pd.Series(proportions_values).fillna(0).values

                # Add temporary column for plotting
                adata.obs[temp_feature_key] = proportions_values

                # Plot spatial distribution
                if "spatial" in adata.obsm:
                    plot_spatial_feature(
                        adata, feature=temp_feature_key, ax=ax, params=params
                    )
                    ax.set_title(cell_type)
                    ax.invert_yaxis()
                else:
                    # Fallback: bar plot for non-spatial data
                    sorted_props = data.proportions[cell_type].sort_values(
                        ascending=False
                    )
                    ax.bar(
                        range(len(sorted_props)),
                        sorted_props.values,
                        alpha=params.alpha,
                    )
                    ax.set_title(cell_type)
                    ax.set_xlabel("Spots (sorted)")
                    ax.set_ylabel("Proportion")

            except Exception as e:
                ax.text(
                    0.5,
                    0.5,
                    f"Error plotting {cell_type}:\n{str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{cell_type} (Error)")

    # Clean up temporary column
    if temp_feature_key in adata.obs.columns:
        del adata.obs[temp_feature_key]

    fig.subplots_adjust(top=0.92, wspace=0.1, hspace=0.3, right=0.98)
    return fig


async def _create_cell_communication_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create cell communication visualization using unified data retrieval.

    This function routes to appropriate visualization based on analysis type and subtype:
    - Spatial analysis: Multi-panel spatial plot using scanpy (LIANA+ official approach)
    - Cluster analysis: LIANA+ visualizations (dotplot, tileplot, circle_plot) or CellPhoneDB

    Supported subtypes for LIANA+ cluster analysis:
    - 'dotplot' (default): Matrix of L-R pairs with size/color encoding
    - 'tileplot': Heatmap tiles showing Source/Target cell type interactions
    - 'circle_plot': Network diagram showing cell-cell communication patterns

    Args:
        adata: AnnData object with cell communication results
        params: Visualization parameters (use params.subtype to select viz type)
        context: MCP context for logging

    Returns:
        matplotlib Figure object
    """
    if context:
        await context.info("Creating cell communication visualization")

    # Use unified data retrieval function
    data = await get_cell_communication_data(adata, context=context)

    if context:
        await context.info(
            f"Using {data.method} results ({data.analysis_type} analysis, "
            f"{len(data.lr_pairs)} LR pairs)"
        )

    # Route to appropriate visualization based on analysis type
    if data.analysis_type == "spatial":
        return await _create_spatial_lr_visualization(adata, data, params, context)
    else:
        # Cluster-based analysis
        if data.method == "cellphonedb":
            # CellPhoneDB supports multiple visualization types via ktplotspy
            subtype = params.subtype or "heatmap"
            if subtype == "dotplot":
                return _create_cellphonedb_dotplot(adata, data, params, context)
            elif subtype == "chord":
                return _create_cellphonedb_chord(adata, data, params, context)
            else:
                # Default: heatmap
                return _create_cellphonedb_heatmap(adata, data, params, context)
        else:
            # LIANA+ cluster results - route based on subtype
            subtype = params.subtype or "dotplot"
            if subtype == "tileplot":
                return await _create_liana_tileplot(adata, data, params, context)
            elif subtype == "circle_plot":
                return await _create_liana_circle_plot(adata, data, params, context)
            else:
                # Default: dotplot
                return await _create_cluster_lr_visualization(
                    adata, data, params, context
                )


async def _create_spatial_lr_visualization(
    adata: ad.AnnData,
    data: CellCommunicationData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """Create spatial L-R visualization using scanpy (official LIANA+ approach).

    Following LIANA+ documentation, spatial bivariate analysis results are best
    visualized using scanpy's sc.pl.spatial() function.

    Args:
        adata: Original AnnData object with spatial coordinates
        data: CellCommunicationData from get_cell_communication_data()
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure with multi-panel spatial LR visualization
    """
    if data.spatial_scores is None or len(data.lr_pairs) == 0:
        raise DataNotFoundError(
            "No spatial communication scores found.\n\n"
            "SOLUTION: Run spatial cell communication analysis first:\n"
            '  analyze_cell_communication(data_id="...", params={"species": "...", ...})'
        )

    # Select top LR pairs to visualize
    n_pairs = min(params.plot_top_pairs or 6, len(data.lr_pairs), 6)

    # Determine top pairs based on global metric
    if len(data.results) > 0:
        # Find global metric column
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
        raise DataNotFoundError(
            "No LR pairs found in spatial results.\n\n"
            "POSSIBLE CAUSES:\n"
            "1. Analysis generated empty results\n"
            "2. Parameters too stringent\n\n"
            "SOLUTION: Re-run analysis with adjusted parameters"
        )

    # Get pair indices in spatial_scores array
    pair_indices = []
    valid_pairs = []
    for pair in top_pairs:
        if pair in data.lr_pairs:
            pair_indices.append(data.lr_pairs.index(pair))
            valid_pairs.append(pair)

    if not valid_pairs:
        # Fallback to first N pairs
        valid_pairs = data.lr_pairs[:n_pairs]
        pair_indices = list(range(len(valid_pairs)))

    # Create figure with subplots
    n_panels = len(valid_pairs)
    n_cols = min(3, n_panels)
    n_rows = (n_panels + n_cols - 1) // n_cols

    figsize = params.figure_size or (5 * n_cols, 4 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)

    # Flatten axes for consistent indexing
    if n_panels == 1:
        axes = np.array([axes])
    axes = np.atleast_1d(axes).flatten()

    # Get spatial coordinates
    x_coords, y_coords = get_spatial_coordinates(adata)

    # Plot each LR pair
    for i, (pair, pair_idx) in enumerate(zip(valid_pairs, pair_indices)):
        ax = axes[i]

        # Get scores for this pair
        if pair_idx < data.spatial_scores.shape[1]:
            scores = data.spatial_scores[:, pair_idx]
        else:
            scores = np.zeros(len(adata))

        # Create scatter plot
        scatter = ax.scatter(
            x_coords,
            y_coords,
            c=scores,
            cmap=params.colormap or "viridis",
            s=params.spot_size or 15,
            alpha=params.alpha or 0.8,
            edgecolors="none",
        )

        # Format pair name for display
        display_name = pair.replace("^", " → ").replace("_", " → ")

        # Add global metric value if available
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

        # Add colorbar
        plt.colorbar(scatter, ax=ax, shrink=0.7, label="Score")

    # Hide unused subplots
    for i in range(n_panels, len(axes)):
        axes[i].set_visible(False)

    plt.suptitle("Spatial Cell Communication", fontsize=14, fontweight="bold")
    plt.tight_layout()

    if context:
        await context.info(f"Created spatial visualization for {n_panels} LR pairs")

    return fig


async def _create_cluster_lr_visualization(
    adata: ad.AnnData,
    data: CellCommunicationData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """Create cluster-based L-R visualization.

    Uses LIANA+ official dotplot for visualization.

    Args:
        adata: AnnData object
        data: CellCommunicationData from get_cell_communication_data()
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure
    """
    require("liana", feature="LIANA+ plotting")
    require("plotnine", feature="LIANA+ plotting")
    import liana as li

    if context:
        await context.info("Using LIANA+ official dotplot")

    try:

        # Determine orderby column (required by LIANA+ dotplot)
        orderby_col = None
        for col in ["magnitude_rank", "specificity_rank", "lr_means"]:
            if col in data.results.columns:
                orderby_col = col
                break

        if orderby_col is None:
            raise ValueError("No valid orderby column found in LIANA results")

        # LIANA+ dotplot returns plotnine.ggplot object
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
            orderby_ascending=True,  # Lower rank = better
            top_n=params.plot_top_pairs or 20,
            inverse_colour=True,
            inverse_size=True,
            cmap=params.colormap or "viridis",
            figure_size=params.figure_size or (10, 8),
            return_fig=True,
        )

        # Convert plotnine to matplotlib Figure
        fig = _plotnine_to_matplotlib(p, params)
        return fig

    except Exception as e:
        raise ProcessingError(
            f"LIANA+ dotplot failed: {e}\n\n"
            "This may be due to incompatible data format or missing columns.\n"
            "Ensure cell communication analysis completed successfully."
        ) from e


async def _create_liana_tileplot(
    adata: ad.AnnData,
    data: CellCommunicationData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """Create LIANA+ tileplot visualization.

    Tileplot shows L-R interactions as heatmap tiles, with separate panels
    for Source and Target cell types.

    Args:
        adata: AnnData object
        data: CellCommunicationData from get_cell_communication_data()
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure
    """
    try:
        import liana as li

        if context:
            await context.info("Creating LIANA+ tileplot")

        # Determine orderby column
        orderby_col = None
        for col in ["magnitude_rank", "specificity_rank", "lr_means"]:
            if col in data.results.columns:
                orderby_col = col
                break

        if orderby_col is None:
            raise ValueError("No valid orderby column found in LIANA results")

        # Determine fill and label columns
        fill_col = (
            "magnitude_rank"
            if "magnitude_rank" in data.results.columns
            else orderby_col
        )
        label_col = "lr_means" if "lr_means" in data.results.columns else fill_col

        # Create tileplot
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

        # Convert plotnine to matplotlib Figure
        fig = _plotnine_to_matplotlib(p, params)
        return fig

    except Exception as e:
        raise ProcessingError(
            f"LIANA+ tileplot failed: {e}\n\n"
            "This may be due to incompatible data format or missing columns.\n"
            "Ensure cell communication analysis completed successfully."
        ) from e


async def _create_liana_circle_plot(
    adata: ad.AnnData,
    data: CellCommunicationData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """Create LIANA+ circle plot (network diagram) visualization.

    Circle plot shows cell-cell communication patterns as a network,
    with nodes representing cell types and edges representing interactions.

    Args:
        adata: AnnData object
        data: CellCommunicationData from get_cell_communication_data()
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure
    """
    try:
        import liana as li

        if context:
            await context.info("Creating LIANA+ circle plot")

        # Determine score and orderby columns
        score_col = None
        for col in ["magnitude_rank", "specificity_rank", "lr_means"]:
            if col in data.results.columns:
                score_col = col
                break

        if score_col is None:
            raise ValueError("No valid score column found in LIANA results")

        # Get groupby key from results
        groupby = params.cluster_key
        if groupby is None:
            # Try to infer from source column
            if "source" in data.results.columns:
                groupby = (
                    data.results["source"].iloc[0] if len(data.results) > 0 else None
                )
            if groupby is None:
                raise ValueError(
                    "cluster_key is required for circle_plot. "
                    "Specify the cell type column used in analysis."
                )

        # Circle plot returns matplotlib axes directly (not plotnine)
        fig_size = params.figure_size or (10, 10)
        fig, ax = plt.subplots(figsize=fig_size)

        li.pl.circle_plot(
            adata=adata,
            uns_key=data.results_key,
            groupby=groupby,
            score_key=score_col,
            inverse_score=True,  # Lower rank = better
            top_n=params.plot_top_pairs * 3 if params.plot_top_pairs else 50,
            orderby=score_col,
            orderby_ascending=True,
            figure_size=fig_size,
        )

        # Get current figure (circle_plot modifies current figure)
        fig = plt.gcf()
        return fig

    except Exception as e:
        raise ProcessingError(
            f"LIANA+ circle_plot failed: {e}\n\n"
            "This may be due to incompatible data format or missing columns.\n"
            "Ensure cell communication analysis completed successfully."
        ) from e


def _plotnine_to_matplotlib(p, params: VisualizationParameters) -> plt.Figure:
    """Convert plotnine ggplot object to matplotlib Figure.

    Args:
        p: plotnine.ggplot object
        params: Visualization parameters for DPI

    Returns:
        matplotlib Figure
    """
    import io

    try:
        from PIL import Image

        # Save plotnine figure to buffer
        buf = io.BytesIO()
        dpi = params.dpi or 300
        p.save(buf, format="png", dpi=dpi, verbose=False)
        buf.seek(0)

        # Load as PIL Image and create matplotlib figure
        img = Image.open(buf)
        fig, ax = plt.subplots(figsize=(img.width / dpi, img.height / dpi), dpi=dpi)
        ax.imshow(img)
        ax.axis("off")
        plt.tight_layout(pad=0)

        buf.close()
        return fig

    except Exception as e:
        raise ProcessingError(
            f"Failed to convert plotnine figure: {e}\n\n"
            "SOLUTION: Ensure Pillow is installed: pip install Pillow"
        ) from e


def _create_cellphonedb_heatmap(
    adata: ad.AnnData,
    data: CellCommunicationData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """Create CellPhoneDB heatmap visualization using ktplotspy.

    Uses ktplotspy.plot_cpdb_heatmap() to create a clustered heatmap showing
    the count of significant interactions between cell type pairs.

    Args:
        adata: AnnData object with cell type annotations
        data: CellCommunicationData with CellPhoneDB results
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure

    References:
        - ktplotspy: https://ktplotspy.readthedocs.io/
        - GitHub: https://github.com/zktuong/ktplotspy
    """
    import ktplotspy as kpy

    means = data.results

    if not isinstance(means, pd.DataFrame) or len(means) == 0:
        raise DataNotFoundError(
            "CellPhoneDB results are empty or invalid format.\n\n"
            "SOLUTION: Re-run CellPhoneDB analysis"
        )

    # Get pvalues (required for ktplotspy heatmap)
    pvalues = adata.uns.get("cellphonedb_pvalues", None)

    if pvalues is None or not isinstance(pvalues, pd.DataFrame):
        raise DataNotFoundError(
            "CellPhoneDB pvalues not found in adata.uns['cellphonedb_pvalues'].\n\n"
            "SOLUTION: Re-run CellPhoneDB analysis to generate pvalues"
        )

    # Use ktplotspy heatmap (shows count of significant interactions)
    grid = kpy.plot_cpdb_heatmap(
        pvals=pvalues,
        title=params.title or "CellPhoneDB: Significant Interactions",
        alpha=0.05,
        symmetrical=True,
    )

    # ClusterGrid has a .fig attribute
    return grid.fig


def _create_cellphonedb_dotplot(
    adata: ad.AnnData,
    data: CellCommunicationData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """Create CellPhoneDB dotplot visualization using ktplotspy.

    Shows L-R interactions as dots where:
    - Size encodes -log10(p-value) or proportion
    - Color encodes mean expression

    Args:
        adata: AnnData object with cell type annotations
        data: CellCommunicationData with CellPhoneDB results
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure
    """
    means = data.results

    if not isinstance(means, pd.DataFrame) or len(means) == 0:
        raise DataNotFoundError(
            "CellPhoneDB results are empty or invalid format.\n\n"
            "SOLUTION: Re-run CellPhoneDB analysis"
        )

    require("ktplotspy", feature="CellPhoneDB dotplot visualization")
    import ktplotspy as kpy

    try:
        # Get pvalues
        pvalues = adata.uns.get("cellphonedb_pvalues", None)

        if pvalues is None or not isinstance(pvalues, pd.DataFrame):
            raise ValueError("Missing pvalues DataFrame for ktplotspy dotplot")

        # Get cell types from cluster_key
        cluster_key = params.cluster_key or "leiden"
        if cluster_key not in adata.obs.columns:
            # Try to find a suitable column
            for key in ["cell_type", "celltype", "leiden", "louvain"]:
                if key in adata.obs.columns:
                    cluster_key = key
                    break

        # ktplotspy.plot_cpdb returns a plotnine ggplot object, not matplotlib Figure
        # We need to convert it to matplotlib Figure
        gg = kpy.plot_cpdb(
            adata=adata,
            cell_type1=".",  # All cell types
            cell_type2=".",  # All cell types
            means=means,
            pvals=pvalues,
            celltype_key=cluster_key,
            genes=None,  # Show all genes
            figsize=params.figure_size or (12, 10),
            title="CellPhoneDB: L-R Interactions",
            max_size=10,
            alpha=0.05,
            keep_significant_only=True,
            standard_scale=True,
        )

        # Convert plotnine ggplot to matplotlib Figure
        fig = gg.draw()
        return fig

    except Exception as e:
        raise ProcessingError(
            f"Failed to create CellPhoneDB dotplot: {str(e)}\n\n"
            "Try using subtype='heatmap' instead."
        )


def _create_cellphonedb_chord(
    adata: ad.AnnData,
    data: CellCommunicationData,
    params: VisualizationParameters,
    context=None,
) -> plt.Figure:
    """Create CellPhoneDB chord/circos diagram using ktplotspy.

    Shows cell-cell communication as a circular network diagram with:
    - Nodes representing cell types
    - Chords/arcs representing interactions between cell types
    - Width/color encoding interaction strength
    - Legend showing L-R pair names on the right side

    Args:
        adata: AnnData object with cell type annotations
        data: CellCommunicationData with CellPhoneDB results
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure
    """
    from matplotlib.lines import Line2D

    means = data.results

    if not isinstance(means, pd.DataFrame) or len(means) == 0:
        raise DataNotFoundError(
            "CellPhoneDB results are empty or invalid format.\n\n"
            "SOLUTION: Re-run CellPhoneDB analysis"
        )

    require("ktplotspy", feature="CellPhoneDB chord visualization")
    import ktplotspy as kpy
    import matplotlib.colors as mcolors

    try:
        # Get pvalues and deconvoluted (required for chord plot)
        pvalues = adata.uns.get("cellphonedb_pvalues", None)
        deconvoluted = adata.uns.get("cellphonedb_deconvoluted", None)

        if pvalues is None or not isinstance(pvalues, pd.DataFrame):
            raise ValueError("Missing pvalues DataFrame for ktplotspy chord plot")

        if deconvoluted is None or not isinstance(deconvoluted, pd.DataFrame):
            raise ValueError(
                "Missing deconvoluted DataFrame for chord plot. "
                "Re-run CellPhoneDB analysis."
            )

        # Get cell types from cluster_key
        cluster_key = params.cluster_key or "leiden"
        if cluster_key not in adata.obs.columns:
            for key in ["cell_type", "celltype", "leiden", "louvain"]:
                if key in adata.obs.columns:
                    cluster_key = key
                    break

        # Generate link_colors dict for coloring links
        # We will create our own legend AFTER ktplotspy returns the circos object
        # This allows us to properly include the legend in bbox_extra_artists
        link_colors = None
        legend_items = []  # Store (label, color) for our custom legend

        if "interacting_pair" in deconvoluted.columns:
            unique_pairs = deconvoluted["interacting_pair"].unique()
            # Limit to top N pairs for readability (configurable via plot_top_pairs)
            n_pairs = min(params.plot_top_pairs or 50, len(unique_pairs))
            top_pairs = unique_pairs[:n_pairs]

            # Use a colormap with enough distinct colors
            if n_pairs <= 10:
                cmap = plt.cm.get_cmap("tab10", 10)
            elif n_pairs <= 20:
                cmap = plt.cm.get_cmap("tab20", 20)
            else:
                cmap = plt.cm.get_cmap("nipy_spectral", n_pairs)

            link_colors = {}
            for i, pair in enumerate(top_pairs):
                # ktplotspy uses 'converted_pair' format internally (hyphen-separated)
                # but the exact format depends on internal processing
                # We provide colors with underscore format as that's what CellPhoneDB uses
                color = mcolors.rgb2hex(cmap(i % cmap.N))
                link_colors[pair] = color
                # Store for our custom legend (display with underscore as in CellPhoneDB)
                legend_items.append((pair, color))

        # ktplotspy.plot_cpdb_chord creates a chord diagram
        # It returns a pycirclize Circos object, not matplotlib Figure
        # NOTE: We don't pass legend_kwargs here - we'll create our own legend
        circos = kpy.plot_cpdb_chord(
            adata=adata,
            means=means,
            pvals=pvalues,
            deconvoluted=deconvoluted,
            celltype_key=cluster_key,
            cell_type1=".",  # All cell types
            cell_type2=".",  # All cell types
            link_colors=link_colors,  # Colors for links (ktplotspy may or may not match)
        )

        # Extract matplotlib Figure from Circos object
        fig = circos.ax.figure

        # Set figure size - make it wider to accommodate legend on right
        fig.set_size_inches(14, 10)

        # Create our own legend with proper positioning
        # This approach gives us control over the legend and ensures it's included
        # when saving with bbox_inches='tight' via bbox_extra_artists
        if legend_items:
            line_handles = [
                Line2D([], [], color=color, label=label, linewidth=2)
                for label, color in legend_items
            ]

            # Create legend and position it further right to avoid overlap with chord
            legend = circos.ax.legend(
                handles=line_handles,
                loc="center left",
                bbox_to_anchor=(1.15, 0.5),  # Move further right to avoid overlap
                fontsize=6,
                frameon=True,
                framealpha=0.9,
                title="L-R Pairs",
                title_fontsize=7,
            )

            # Store legend as figure attribute for bbox_extra_artists in image_utils
            # This allows the save function to include the legend in the bounding box
            fig._chatspatial_extra_artists = [legend]

        return fig

    except Exception as e:
        raise ProcessingError(
            f"Failed to create CellPhoneDB chord diagram: {str(e)}\n\n"
            "Try using subtype='heatmap' instead."
        )


async def _create_multi_gene_umap_visualization(
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
        default_title="",  # No default title for cleaner visualization
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

    # No need for subplots_adjust - spacing is controlled by GridSpec in setup_multi_panel_figure
    # GridSpec wspace/hspace parameters are the ONLY effective way to control spacing
    # when using ax.set_aspect('equal') for spatial plots
    return fig


async def _create_multi_gene_visualization(
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
    # Use explicit subplots_adjust for better spacing control (like LR pairs)
    fig, axes = setup_multi_panel_figure(
        n_panels=len(available_genes),
        params=params,
        default_title="",  # No default title for cleaner visualization
        use_tight_layout=False,  # Changed from True to use explicit spacing control
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
                        # Use configurable colorbar size and padding for optimal spacing
                        if params.show_colorbar:
                            divider = make_axes_locatable(ax)
                            cax = divider.append_axes(
                                "right",
                                size=params.colorbar_size,
                                pad=params.colorbar_pad,
                            )
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
    # Reduced wspace from 0.2 to 0.1 for tighter left-right spacing
    fig.subplots_adjust(top=0.92, wspace=0.1, hspace=0.3, right=0.98)
    return fig


async def _create_lr_pairs_visualization(
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
                # Use configurable colorbar size and padding for optimal spacing
                if params.show_colorbar and ax.collections:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right",
                        size=params.colorbar_size,
                        pad=params.colorbar_pad,
                    )
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
                # Use configurable colorbar size and padding for optimal spacing
                if params.show_colorbar and ax.collections:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right",
                        size=params.colorbar_size,
                        pad=params.colorbar_pad,
                    )
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
    # Reduced wspace from 0.2 to 0.1 for tighter left-right spacing
    fig.subplots_adjust(top=0.92, wspace=0.1, hspace=0.3, right=0.98)
    return fig


async def _create_rna_velocity_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create RNA velocity visualization based on subtype

    Dispatcher function that routes to appropriate scVelo visualization based on subtype.

    Args:
        adata: AnnData object with computed RNA velocity
        params: Visualization parameters including subtype
        context: MCP context

    Returns:
        Matplotlib figure with RNA velocity visualization

    Subtypes:
        - stream (default): Velocity embedding stream plot (scv.pl.velocity_embedding_stream)
        - phase: Phase plot showing spliced vs unspliced (scv.pl.velocity)
        - proportions: Pie chart of spliced/unspliced ratios (scv.pl.proportions)
        - heatmap: Gene expression ordered by latent_time (scv.pl.heatmap)
        - paga: PAGA with velocity arrows (scv.pl.paga)
    """
    # Default to stream if no subtype specified
    subtype = params.subtype or "stream"

    if context:
        await context.info(f"Creating RNA velocity visualization (subtype: {subtype})")

    if subtype == "stream":
        return await _create_velocity_stream_plot(adata, params, context)
    elif subtype == "phase":
        return await _create_velocity_phase_plot(adata, params, context)
    elif subtype == "proportions":
        return await _create_velocity_proportions_plot(adata, params, context)
    elif subtype == "heatmap":
        return await _create_velocity_heatmap(adata, params, context)
    elif subtype == "paga":
        return await _create_velocity_paga_plot(adata, params, context)
    else:
        raise ParameterError(
            f"Unsupported subtype for rna_velocity: '{subtype}'. "
            f"Available subtypes: stream, phase, proportions, heatmap, paga"
        )


async def _create_velocity_stream_plot(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create RNA velocity stream plot using scv.pl.velocity_embedding_stream

    Args:
        adata: AnnData object with computed RNA velocity
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with RNA velocity stream plot

    Data requirements:
        - adata.uns['velocity_graph']: Velocity transition graph
        - adata.obsm['X_umap'] or 'spatial': Embedding for visualization
    """
    require("scvelo", feature="RNA velocity visualization")
    import scvelo as scv

    # Check if RNA velocity has been computed
    if "velocity_graph" not in adata.uns:
        raise DataNotFoundError(
            "RNA velocity has not been computed. Please run 'analyze_velocity_data' first."
        )

    # Determine basis for plotting
    basis = params.basis or "spatial"
    basis_key = f"X_{basis}" if basis != "spatial" else "spatial"

    if basis_key not in adata.obsm:
        # Try to find an alternative basis
        if "spatial" in adata.obsm:
            basis = "spatial"
            if context:
                await context.info("Using 'spatial' as basis")
        elif "X_umap" in adata.obsm:
            basis = "umap"
            if context:
                await context.info("Using 'umap' as basis")
        elif "X_pca" in adata.obsm:
            basis = "pca"
            if context:
                await context.info("Using 'pca' as basis")
        else:
            available_bases = [
                k.replace("X_", "") for k in adata.obsm.keys() if k.startswith("X_")
            ]
            available_bases.extend([k for k in adata.obsm.keys() if k == "spatial"])
            raise DataCompatibilityError(
                f"Basis '{params.basis or 'spatial'}' not found. "
                f"Available bases: {available_bases}"
            )

    # Prepare feature for coloring
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

    # Create RNA velocity stream plot using official scVelo function
    scv.pl.velocity_embedding_stream(
        adata,
        basis=basis,
        color=feature,
        ax=ax,
        show=False,
        alpha=params.alpha,
        legend_loc="right margin" if feature and feature in adata.obs.columns else None,
        frameon=params.show_axes,
        title="",
    )

    # Set title
    title = params.title or f"RNA Velocity Stream on {basis.capitalize()}"
    ax.set_title(title, fontsize=14)

    # Handle spatial coordinates orientation
    if basis == "spatial":
        ax.invert_yaxis()

    plt.tight_layout()
    return fig


async def _create_velocity_phase_plot(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create velocity phase plot using scv.pl.velocity

    Shows spliced vs unspliced counts with fitted velocity model for specified genes.

    Args:
        adata: AnnData object with computed RNA velocity
        params: Visualization parameters (feature should be gene name(s))
        context: MCP context

    Returns:
        Matplotlib figure with phase plot

    Data requirements:
        - adata.layers['velocity']: Velocity vectors
        - adata.layers['Ms']: Smoothed spliced counts
        - adata.layers['Mu']: Smoothed unspliced counts
        - adata.var_names: Gene names for feature parameter
    """
    require("scvelo", feature="velocity phase plots")
    import scvelo as scv

    # Check required layers
    required_layers = ["velocity", "Ms", "Mu"]
    missing_layers = [l for l in required_layers if l not in adata.layers]
    if missing_layers:
        raise DataNotFoundError(
            f"Missing required layers for phase plot: {missing_layers}. "
            f"Please run 'analyze_velocity_data' first."
        )

    # Get genes to plot
    if params.feature:
        if isinstance(params.feature, str):
            var_names = [params.feature]
        else:
            var_names = list(params.feature)
    else:
        # Select top velocity genes if no feature specified
        if "velocity_genes" in adata.var.columns:
            velocity_genes = adata.var_names[adata.var["velocity_genes"]]
            var_names = list(velocity_genes[:4])  # Top 4 velocity genes
        else:
            var_names = list(adata.var_names[:4])

    # Validate genes exist
    valid_genes = [g for g in var_names if g in adata.var_names]
    if not valid_genes:
        raise DataNotFoundError(
            f"None of the specified genes found in data: {var_names}. "
            f"Available genes (first 10): {list(adata.var_names[:10])}"
        )

    if context:
        await context.info(f"Creating phase plot for genes: {valid_genes}")

    # Determine basis for background embedding
    basis = params.basis or "umap"
    if f"X_{basis}" not in adata.obsm and basis != "spatial":
        if "X_umap" in adata.obsm:
            basis = "umap"
        elif "spatial" in adata.obsm:
            basis = "spatial"

    # Create figure using scVelo's velocity plot (phase plot)
    # scv.pl.velocity creates its own figure when multiple genes specified
    figsize = params.figure_size or (4 * len(valid_genes), 4)

    # Determine color for plots
    color = params.cluster_key if params.cluster_key else None

    # scv.pl.velocity handles multi-gene layouts internally
    scv.pl.velocity(
        adata,
        var_names=valid_genes,
        basis=basis,
        color=color,
        figsize=figsize,
        dpi=params.dpi,
        show=False,
        ncols=len(valid_genes),
    )

    fig = plt.gcf()
    title = params.title or "RNA Velocity Phase Plot"
    fig.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


async def _create_velocity_proportions_plot(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create velocity proportions plot using scv.pl.proportions

    Shows pie chart of spliced/unspliced RNA proportions per cluster.

    Args:
        adata: AnnData object with spliced/unspliced layers
        params: Visualization parameters (cluster_key for grouping)
        context: MCP context

    Returns:
        Matplotlib figure with proportions pie chart

    Data requirements:
        - adata.layers['spliced']: Spliced counts
        - adata.layers['unspliced']: Unspliced counts
        - adata.obs[cluster_key]: Cluster labels for grouping
    """
    require("scvelo", feature="proportions plot")
    import scvelo as scv

    # Check required layers
    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise DataNotFoundError(
            "Spliced and unspliced layers are required for proportions plot. "
            "Your data may not contain RNA velocity information."
        )

    # Determine cluster key
    cluster_key = params.cluster_key
    if not cluster_key:
        # Try to find a categorical column
        categorical_cols = [
            col
            for col in adata.obs.columns
            if adata.obs[col].dtype.name in ["object", "category"]
        ]
        if categorical_cols:
            cluster_key = categorical_cols[0]
            if context:
                await context.info(f"Using cluster_key: '{cluster_key}'")
        else:
            raise ParameterError(
                "cluster_key is required for proportions plot. "
                f"Available columns: {list(adata.obs.columns)[:10]}"
            )

    if cluster_key not in adata.obs.columns:
        raise DataNotFoundError(
            f"Cluster key '{cluster_key}' not found in data. "
            f"Available columns: {list(adata.obs.columns)[:10]}"
        )

    if context:
        await context.info(f"Creating proportions plot grouped by '{cluster_key}'")

    # Create figure using scVelo's proportions plot
    # Note: scv.pl.proportions does NOT support ax parameter - it creates its own figure
    figsize = params.figure_size or (12, 4)

    scv.pl.proportions(
        adata,
        groupby=cluster_key,
        figsize=figsize,
        dpi=params.dpi,
        show=False,
    )

    fig = plt.gcf()
    title = params.title or f"Spliced/Unspliced Proportions by {cluster_key}"
    fig.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


async def _create_velocity_heatmap(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create velocity heatmap using scv.pl.heatmap

    Shows gene expression patterns ordered by latent time.

    Args:
        adata: AnnData object with latent_time computed
        params: Visualization parameters (feature for gene selection)
        context: MCP context

    Returns:
        Matplotlib figure with velocity heatmap

    Data requirements:
        - adata.obs['latent_time']: Latent time from dynamical model
        - adata.var['velocity_genes']: Velocity genes (optional, for gene selection)
    """
    require("scvelo", feature="velocity heatmap")
    import scvelo as scv

    # Check for latent_time
    if "latent_time" not in adata.obs.columns:
        raise DataNotFoundError(
            "latent_time is required for velocity heatmap. "
            "Please run 'analyze_velocity_data' with dynamical mode first, "
            "or use analyze_trajectory_data for pseudotime calculation."
        )

    # Get genes to plot
    if params.feature:
        if isinstance(params.feature, str):
            var_names = [params.feature]
        else:
            var_names = list(params.feature)
        # Validate genes
        valid_genes = [g for g in var_names if g in adata.var_names]
        if not valid_genes:
            raise DataNotFoundError(f"None of the specified genes found: {var_names}")
        var_names = valid_genes
    else:
        # Select top velocity genes
        if "velocity_genes" in adata.var.columns:
            velocity_genes = adata.var_names[adata.var["velocity_genes"]]
            var_names = list(velocity_genes[:50])  # Top 50 velocity genes
        else:
            # Select highly variable genes
            if "highly_variable" in adata.var.columns:
                hvg = adata.var_names[adata.var["highly_variable"]]
                var_names = list(hvg[:50])
            else:
                var_names = list(adata.var_names[:50])

    if context:
        await context.info(f"Creating velocity heatmap with {len(var_names)} genes")

    # Create figure using scVelo's heatmap
    figsize = params.figure_size or (12, 8)

    # scvelo.pl.heatmap creates its own figure, capture it
    scv.pl.heatmap(
        adata,
        var_names=var_names,
        sortby="latent_time",
        col_color=params.cluster_key,
        n_convolve=30,
        show=False,
        figsize=figsize,
    )

    fig = plt.gcf()
    fig.set_dpi(params.dpi)

    # Only add title if explicitly provided by user
    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


async def _create_velocity_paga_plot(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create PAGA plot with velocity using scv.pl.paga

    Shows partition-based graph abstraction with directed velocity arrows.

    Args:
        adata: AnnData object with velocity and PAGA computed
        params: Visualization parameters (cluster_key for PAGA groups)
        context: MCP context

    Returns:
        Matplotlib figure with velocity PAGA plot

    Data requirements:
        - adata.uns['velocity_graph']: Velocity transition graph
        - adata.uns['paga']: PAGA results (computed by scv.tl.paga)
        - adata.obs[cluster_key]: Cluster labels used for PAGA
    """
    require("scvelo", feature="velocity PAGA plot")
    import scvelo as scv

    # Check for velocity graph
    if "velocity_graph" not in adata.uns:
        raise DataNotFoundError(
            "velocity_graph is required for PAGA visualization. "
            "Please run 'analyze_velocity_data' first."
        )

    # Determine cluster key
    cluster_key = params.cluster_key
    if not cluster_key:
        # Try to find from PAGA results
        if "paga" in adata.uns and "groups" in adata.uns.get("paga", {}):
            cluster_key = adata.uns["paga"].get("groups")
        else:
            categorical_cols = [
                col
                for col in adata.obs.columns
                if adata.obs[col].dtype.name in ["object", "category"]
            ]
            if categorical_cols:
                cluster_key = categorical_cols[0]

    if not cluster_key or cluster_key not in adata.obs.columns:
        raise ParameterError(
            f"cluster_key is required for PAGA plot. "
            f"Available columns: {list(adata.obs.columns)[:10]}"
        )

    # Compute PAGA if not already done
    if "paga" not in adata.uns:
        if context:
            await context.info(f"Computing PAGA for cluster_key='{cluster_key}'")
        import scanpy as sc

        sc.tl.paga(adata, groups=cluster_key)
        scv.tl.paga(adata, groups=cluster_key)

    if context:
        await context.info(f"Creating velocity PAGA plot for '{cluster_key}'")

    # Determine basis
    basis = params.basis or "umap"
    if f"X_{basis}" not in adata.obsm:
        if "X_umap" in adata.obsm:
            basis = "umap"
        elif "spatial" in adata.obsm:
            basis = "spatial"

    # Create figure
    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    # Create PAGA plot using scVelo
    scv.pl.paga(
        adata,
        basis=basis,
        color=cluster_key,
        ax=ax,
        show=False,
        frameon=params.show_axes,
    )

    # Only add title if explicitly provided by user
    if params.title:
        ax.set_title(params.title, fontsize=14)

    # Handle spatial coordinates orientation
    if basis == "spatial":
        ax.invert_yaxis()

    plt.tight_layout()
    return fig


async def _create_spatial_statistics_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create spatial statistics visualization based on subtype

    Args:
        adata: AnnData object with spatial statistics results
        params: Visualization parameters including subtype
        context: MCP context

    Returns:
        Matplotlib figure with spatial statistics visualization
    """
    if context:
        await context.info(
            f"Creating {params.subtype} spatial statistics visualization"
        )

    if params.subtype == "neighborhood":
        return await _create_neighborhood_enrichment_visualization(
            adata, params, context
        )
    elif params.subtype == "co_occurrence":
        return await _create_co_occurrence_visualization(adata, params, context)
    elif params.subtype == "ripley":
        return await _create_ripley_visualization(adata, params, context)
    elif params.subtype == "moran":
        return await _create_moran_visualization(adata, params, context)
    elif params.subtype == "centrality":
        return await _create_centrality_visualization(adata, params, context)
    elif params.subtype == "getis_ord":
        return await _create_getis_ord_visualization(adata, params, context)
    else:
        raise ParameterError(
            f"Unsupported subtype for spatial_statistics: {params.subtype}"
        )


async def _create_neighborhood_enrichment_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create neighborhood enrichment visualization using squidpy."""
    require("squidpy", feature="neighborhood enrichment visualization")
    import squidpy as sq

    # Infer cluster_key from params or existing results
    cluster_key = params.cluster_key
    if not cluster_key:
        enrichment_keys = [
            k for k in adata.uns.keys() if k.endswith("_nhood_enrichment")
        ]
        if enrichment_keys:
            cluster_key = enrichment_keys[0].replace("_nhood_enrichment", "")
            if context:
                await context.info(f"Inferred cluster_key: '{cluster_key}'")
        else:
            categorical_cols = [
                col
                for col in adata.obs.columns
                if adata.obs[col].dtype.name in ["object", "category"]
            ][:10]
            raise ValueError(
                f"cluster_key required. Available: {', '.join(categorical_cols)}"
            )

    enrichment_key = f"{cluster_key}_nhood_enrichment"
    if enrichment_key not in adata.uns:
        raise DataNotFoundError(
            f"Neighborhood enrichment not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    # Use squidpy's heatmap visualization
    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    sq.pl.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        cmap=params.colormap or "coolwarm",
        ax=ax,
        title=params.title or f"Neighborhood Enrichment ({cluster_key})",
    )

    plt.tight_layout()
    return fig


async def _create_co_occurrence_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create co-occurrence analysis visualization using squidpy."""
    require("squidpy", feature="co-occurrence visualization")
    import squidpy as sq

    # Infer cluster_key from params or existing results
    cluster_key = params.cluster_key
    if not cluster_key:
        co_occurrence_keys = [
            k for k in adata.uns.keys() if k.endswith("_co_occurrence")
        ]
        if co_occurrence_keys:
            cluster_key = co_occurrence_keys[0].replace("_co_occurrence", "")
            if context:
                await context.info(f"Inferred cluster_key: '{cluster_key}'")
        else:
            categorical_cols = [
                col
                for col in adata.obs.columns
                if adata.obs[col].dtype.name in ["object", "category"]
            ][:10]
            raise ValueError(
                f"cluster_key required. Available: {', '.join(categorical_cols)}"
            )

    co_occurrence_key = f"{cluster_key}_co_occurrence"
    if co_occurrence_key not in adata.uns:
        raise DataNotFoundError(
            f"Co-occurrence not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    # Use squidpy's co_occurrence plot
    categories = adata.obs[cluster_key].cat.categories.tolist()
    clusters_to_show = categories[: min(4, len(categories))]

    figsize = params.figure_size or (12, 10)

    # sq.pl.co_occurrence manages its own figure, doesn't accept ax parameter
    sq.pl.co_occurrence(
        adata,
        cluster_key=cluster_key,
        clusters=clusters_to_show,
        figsize=figsize,
        dpi=params.dpi,
    )

    fig = plt.gcf()
    if params.title:
        fig.suptitle(params.title)

    plt.tight_layout()
    return fig


async def _create_ripley_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create Ripley's function visualization using squidpy."""
    require("squidpy", feature="Ripley visualization")
    import squidpy as sq

    # Infer cluster_key from params or existing results
    cluster_key = params.cluster_key
    if not cluster_key:
        ripley_keys = [k for k in adata.uns.keys() if k.endswith("_ripley_L")]
        if ripley_keys:
            cluster_key = ripley_keys[0].replace("_ripley_L", "")
            if context:
                await context.info(f"Inferred cluster_key: '{cluster_key}'")
        else:
            categorical_cols = [
                col
                for col in adata.obs.columns
                if adata.obs[col].dtype.name in ["object", "category"]
            ][:10]
            raise ValueError(
                f"cluster_key required. Available: {', '.join(categorical_cols)}"
            )

    ripley_key = f"{cluster_key}_ripley_L"
    if ripley_key not in adata.uns:
        raise DataNotFoundError(
            f"Ripley results not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' and analysis_type='ripley' first."
        )

    # Use squidpy's ripley plot
    figsize = params.figure_size or (10, 8)
    fig, ax = plt.subplots(figsize=figsize, dpi=params.dpi)

    sq.pl.ripley(adata, cluster_key=cluster_key, mode="L", plot_sims=True, ax=ax)

    if params.title:
        ax.set_title(params.title)

    plt.tight_layout()
    return fig


async def _create_moran_visualization(
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


async def _create_centrality_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create centrality scores visualization using squidpy."""
    require("squidpy", feature="centrality visualization")
    import squidpy as sq

    # Infer cluster_key from params or existing results
    cluster_key = params.cluster_key
    if not cluster_key:
        centrality_keys = [
            k for k in adata.uns.keys() if k.endswith("_centrality_scores")
        ]
        if centrality_keys:
            cluster_key = centrality_keys[0].replace("_centrality_scores", "")
            if context:
                await context.info(f"Inferred cluster_key: '{cluster_key}'")
        else:
            categorical_cols = [
                col
                for col in adata.obs.columns
                if adata.obs[col].dtype.name in ["object", "category"]
            ][:10]
            raise ValueError(
                f"cluster_key required. Available: {', '.join(categorical_cols)}"
            )

    centrality_key = f"{cluster_key}_centrality_scores"
    if centrality_key not in adata.uns:
        raise DataNotFoundError(
            f"Centrality scores not found. Run analyze_spatial_statistics "
            f"with cluster_key='{cluster_key}' first."
        )

    # Use squidpy's centrality_scores plot
    figsize = params.figure_size or (10, 8)

    # sq.pl.centrality_scores manages its own figure, doesn't accept ax parameter
    sq.pl.centrality_scores(
        adata,
        cluster_key=cluster_key,
        figsize=figsize,
        dpi=params.dpi,
    )

    fig = plt.gcf()
    if params.title:
        fig.suptitle(params.title)

    plt.tight_layout()
    return fig


async def _create_getis_ord_visualization(
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


async def _create_trajectory_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create trajectory visualization based on subtype

    Dispatcher function that routes to appropriate trajectory visualization based on subtype.

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
    # Default to pseudotime if no subtype specified
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


async def _create_trajectory_pseudotime_plot(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create trajectory pseudotime visualization

    Shows pseudotime on embedding with optional velocity stream plot.

    Args:
        adata: AnnData object with computed trajectory/pseudotime
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with pseudotime visualization

    Data requirements:
        - adata.obs['*pseudotime*']: Any pseudotime column
        - adata.obsm['X_umap'] or 'spatial': Embedding for visualization
        - adata.uns['velocity_graph']: Optional, for velocity stream panel
    """
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
                "No pseudotime information found. Please run trajectory analysis "
                "first or specify a pseudotime column in the 'feature' parameter."
            )

    if pseudotime_key not in adata.obs.columns:
        raise DataNotFoundError(
            f"Pseudotime column '{pseudotime_key}' not found. "
            "Please run trajectory analysis first."
        )

    # Check if RNA velocity is available
    has_velocity = "velocity_graph" in adata.uns

    # Determine basis for plotting
    basis = params.basis or "spatial"
    basis_key = f"X_{basis}" if basis != "spatial" else "spatial"

    if basis_key not in adata.obsm:
        if "spatial" in adata.obsm:
            basis = "spatial"
        elif "X_umap" in adata.obsm:
            basis = "umap"
        elif "X_pca" in adata.obsm:
            basis = "pca"
        else:
            available_bases = [
                k.replace("X_", "") for k in adata.obsm.keys() if k.startswith("X_")
            ]
            available_bases.extend([k for k in adata.obsm.keys() if k == "spatial"])
            raise DataCompatibilityError(
                f"Basis not found. Available bases: {available_bases}"
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

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


async def _create_cellrank_circular_projection(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create CellRank circular projection using cr.pl.circular_projection

    Shows fate probabilities in a circular layout.

    Args:
        adata: AnnData object with CellRank analysis
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with circular projection

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
            "CellRank fate probabilities not found. "
            "Please run 'analyze_trajectory_data' with method='cellrank' first."
        )

    if context:
        await context.info("Creating CellRank circular projection")

    # Determine keys for coloring
    keys = [params.cluster_key] if params.cluster_key else None
    if not keys:
        # Try to find categorical columns
        categorical_cols = [
            col
            for col in adata.obs.columns
            if adata.obs[col].dtype.name in ["object", "category"]
        ][:3]
        keys = categorical_cols if categorical_cols else None

    # Create figure - CellRank creates its own figure, no ax parameter supported
    # Note: CellRank has compatibility issues with show parameter in some versions
    figsize = params.figure_size or (10, 10)

    import matplotlib

    backend = matplotlib.get_backend()
    matplotlib.use("Agg")  # Force non-interactive backend

    cr.pl.circular_projection(
        adata,
        keys=keys,
        figsize=figsize,
        dpi=params.dpi,
    )
    fig = plt.gcf()
    matplotlib.use(backend)  # Restore backend

    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_cellrank_fate_map(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create CellRank aggregated fate probabilities using cr.pl.aggregate_fate_probabilities

    Shows fate probabilities aggregated by cluster as bar, paga, or heatmap.

    Args:
        adata: AnnData object with CellRank analysis
        params: Visualization parameters (cluster_key for grouping)
        context: MCP context

    Returns:
        Matplotlib figure with aggregated fate map

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
            "CellRank fate probabilities not found. "
            "Please run 'analyze_trajectory_data' with method='cellrank' first."
        )

    # Determine cluster key
    cluster_key = params.cluster_key
    if not cluster_key:
        categorical_cols = [
            col
            for col in adata.obs.columns
            if adata.obs[col].dtype.name in ["object", "category"]
        ]
        if categorical_cols:
            cluster_key = categorical_cols[0]
            if context:
                await context.info(f"Using cluster_key: '{cluster_key}'")
        else:
            raise ParameterError("cluster_key is required for fate map visualization.")

    if context:
        await context.info(f"Creating CellRank fate map for '{cluster_key}'")

    # Create figure - CellRank creates its own figure, no ax parameter supported
    # Note: CellRank has compatibility issues with show parameter in some versions
    figsize = params.figure_size or (12, 6)

    import matplotlib

    backend = matplotlib.get_backend()
    matplotlib.use("Agg")  # Force non-interactive backend

    # Default to 'bar' mode - most informative
    cr.pl.aggregate_fate_probabilities(
        adata,
        cluster_key=cluster_key,
        mode="bar",
        figsize=figsize,
        dpi=params.dpi,
    )
    fig = plt.gcf()
    matplotlib.use(backend)  # Restore backend

    title = params.title or f"CellRank Fate Probabilities by {cluster_key}"
    fig.suptitle(title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_cellrank_gene_trends(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create CellRank gene expression trends using cr.pl.gene_trends

    Shows gene expression trends along lineages/pseudotime.

    Args:
        adata: AnnData object with CellRank analysis
        params: Visualization parameters (feature for gene selection)
        context: MCP context

    Returns:
        Matplotlib figure with gene trends

    Data requirements:
        - adata.obsm['lineages_fwd'] or 'to_terminal_states': Fate probabilities
        - adata.obs['latent_time'] or similar pseudotime
        - Gene expression in adata.X
    """
    require("cellrank", feature="gene trends")
    import cellrank as cr

    # Import GAM model preparation from trajectory module (separation of concerns)
    from .trajectory import prepare_gam_model_for_visualization

    # Check for fate probabilities
    fate_key_candidates = ["lineages_fwd", "to_terminal_states"]
    fate_key = None
    for key in fate_key_candidates:
        if key in adata.obsm:
            fate_key = key
            break

    if not fate_key:
        raise DataNotFoundError(
            "CellRank fate probabilities not found. "
            "Please run 'analyze_trajectory_data' with method='cellrank' first."
        )

    # Find time key
    time_key = None
    time_candidates = ["latent_time", "palantir_pseudotime", "dpt_pseudotime"]
    for key in time_candidates:
        if key in adata.obs.columns:
            time_key = key
            break

    if not time_key:
        raise DataNotFoundError(
            "No pseudotime found for gene trends. "
            "Please ensure trajectory analysis has been run."
        )

    # Get genes to plot
    if params.feature:
        if isinstance(params.feature, str):
            genes = [params.feature]
        else:
            genes = list(params.feature)
        # Validate genes
        valid_genes = [g for g in genes if g in adata.var_names]
        if not valid_genes:
            raise DataNotFoundError(f"None of the specified genes found: {genes}")
        genes = valid_genes[:6]  # Limit to 6 genes
    else:
        # Select highly variable genes
        if "highly_variable" in adata.var.columns:
            hvg = adata.var_names[adata.var["highly_variable"]]
            genes = list(hvg[:6])
        else:
            genes = list(adata.var_names[:6])

    if context:
        await context.info(f"Creating gene trends for: {genes}")

    # Create figure - CellRank creates its own
    figsize = params.figure_size or (12, 3 * len(genes))

    # Force Agg backend to avoid display issues
    import matplotlib

    backend = matplotlib.get_backend()
    matplotlib.use("Agg")

    # Use trajectory module for GAM model preparation (computation logic)
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
        n_jobs=1,  # Avoid multiprocessing issues in MCP context
        show_progress_bar=False,
    )

    fig = plt.gcf()
    fig.set_dpi(params.dpi)
    matplotlib.use(backend)  # Restore backend

    # Only add title if explicitly provided by user
    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_cellrank_fate_heatmap(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create CellRank fate heatmap using cr.pl.heatmap

    Shows smoothed gene expression ordered by pseudotime per lineage.

    Args:
        adata: AnnData object with CellRank analysis
        params: Visualization parameters (feature for gene selection)
        context: MCP context

    Returns:
        Matplotlib figure with fate heatmap

    Data requirements:
        - adata.obsm['lineages_fwd'] or 'to_terminal_states': Fate probabilities
        - adata.obs['latent_time'] or similar pseudotime
        - Gene expression in adata.X
    """
    require("cellrank", feature="fate heatmap")
    import cellrank as cr

    # Import GAM model preparation from trajectory module (separation of concerns)
    from .trajectory import prepare_gam_model_for_visualization

    # Check for fate probabilities
    fate_key_candidates = ["lineages_fwd", "to_terminal_states"]
    fate_key = None
    for key in fate_key_candidates:
        if key in adata.obsm:
            fate_key = key
            break

    if not fate_key:
        raise DataNotFoundError(
            "CellRank fate probabilities not found. "
            "Please run 'analyze_trajectory_data' with method='cellrank' first."
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

    # Force Agg backend to avoid display issues
    import matplotlib

    backend = matplotlib.get_backend()
    matplotlib.use("Agg")

    # Use trajectory module for GAM model preparation (computation logic)
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
        n_jobs=1,  # Avoid multiprocessing issues in MCP context
        show_progress_bar=False,
    )

    fig = plt.gcf()
    fig.set_dpi(params.dpi)
    matplotlib.use(backend)  # Restore backend

    # Only add title if explicitly provided by user
    if params.title:
        fig.suptitle(params.title, fontsize=14, y=1.02)

    plt.tight_layout()
    return fig


async def _create_palantir_results(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create Palantir comprehensive results visualization

    Shows pseudotime, entropy, and fate probabilities in a multi-panel figure.

    Args:
        adata: AnnData object with Palantir analysis
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with Palantir results

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
            "Palantir results not found. "
            "Please run 'analyze_trajectory_data' with method='palantir' first."
        )

    if context:
        await context.info("Creating Palantir results visualization")

    # Determine basis
    basis = params.basis or "umap"
    if f"X_{basis}" not in adata.obsm:
        if "X_umap" in adata.obsm:
            basis = "umap"
        elif "spatial" in adata.obsm:
            basis = "spatial"
        else:
            basis = "pca"

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
        # Show dominant fate
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


async def _create_gene_correlation_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create gene correlation visualization using seaborn clustermap."""
    import pandas as pd
    import seaborn as sns

    # Get validated genes
    available_genes = await get_validated_features(
        adata, params, max_features=10, context=context
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

    # Create DataFrame for correlation
    expr_df = pd.DataFrame(expr_matrix, columns=available_genes)

    # Calculate correlation using pandas (handles method selection)
    corr_df = expr_df.corr(method=params.correlation_method)

    # Use seaborn clustermap for visualization with clustering
    figsize = params.figure_size or (
        max(8, len(available_genes)),
        max(8, len(available_genes)),
    )

    g = sns.clustermap(
        corr_df,
        cmap=params.colormap,
        center=0,
        annot=True,
        fmt=".2f",
        square=True,
        figsize=figsize,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
    )

    title = params.title or f"Gene Correlation ({params.correlation_method.title()})"
    g.fig.suptitle(title, y=1.02, fontsize=14)

    return g.fig


async def _create_enrichment_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create spatial enrichment visualization (unified for EnrichMap and standard spatial)

    Routes to appropriate visualization based on params:
    - If plot_type="violin": Enrichment scores violin plot by cluster
    - If subtype starts with "spatial_": EnrichMap spatial visualizations
      (spatial_score, spatial_correlogram, spatial_variogram, spatial_cross_correlation)
    - Default: Standard spatial scatter plot of enrichment scores

    Args:
        adata: AnnData object with enrichment scores
        params: Visualization parameters
            - subtype: Spatial EnrichMap type ('spatial_score', etc.) or None for standard
            - feature: Score column name or signature name, or list for multi-panel
            - cluster_key: For violin plots, the grouping variable
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

    # Check if user wants violin plot by cluster
    if params.plot_type == "violin":
        # Determine grouping variable - require explicit specification
        if params.cluster_key:
            group_by = params.cluster_key
        else:
            # No cluster_key provided - show error with available options
            categorical_cols = [
                col
                for col in adata.obs.columns
                if adata.obs[col].dtype.name in ["object", "category"]
            ]
            raise ValueError(
                "Enrichment violin plot requires 'cluster_key' parameter.\n\n"
                f"Available categorical columns ({len(categorical_cols)} total):\n"
                f"  {', '.join(categorical_cols[:15])}\n\n"
                "SOLUTION: Specify cluster_key explicitly:\n"
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

    # Handle spatial EnrichMap visualizations through unified subtype
    if params.subtype and params.subtype.startswith("spatial_"):
        # spatial_ requests require EnrichMap - do not fallback silently
        if em is None:
            raise ProcessingError(
                f"Spatial enrichment visualization ('{params.subtype}') requires EnrichMap.\n\n"
                "INSTALLATION:\n"
                "  pip install enrichmap\n\n"
                "EnrichMap provides spatial autocorrelation analysis (Moran's I, variogram, etc.) "
                "for enrichment scores. Without it, these statistical visualizations cannot be generated."
            )

        # Ensure EnrichMap compatibility: add required metadata
        _ensure_enrichmap_compatibility(adata)

        # Get library_id for single-sample case
        library_id = adata.obs["library_id"].unique()[0]

        # EnrichMap functions create their own figures, so we don't pre-create
        # We'll capture the figure after EnrichMap creates it
        try:
            # spatial_cross_correlation doesn't need feature - uses first 2 pathways
            if params.subtype == "spatial_cross_correlation":
                # Need at least 2 pathways
                if "enrichment_gene_sets" not in adata.uns:
                    raise DataNotFoundError(
                        "enrichment_gene_sets not found in adata.uns"
                    )
                pathways = list(adata.uns["enrichment_gene_sets"].keys())
                if len(pathways) < 2:
                    raise DataNotFoundError(
                        "Need at least 2 pathways for cross-correlation"
                    )
                score_x = f"{pathways[0]}_score"
                score_y = f"{pathways[1]}_score"
                em.pl.cross_moran_scatter(
                    adata, score_x=score_x, score_y=score_y, library_id=library_id
                )
            else:
                # Other spatial enrichment types need feature parameter
                if params.feature:
                    score_col = f"{params.feature}_score"
                    if score_col not in adata.obs.columns:
                        raise DataNotFoundError(f"Score column '{score_col}' not found")
                else:
                    raise DataNotFoundError(
                        "Feature parameter is required for spatial enrichment visualization"
                    )

                if params.subtype == "spatial_correlogram":
                    em.pl.morans_correlogram(
                        adata, score_key=score_col, library_id=library_id
                    )
                elif params.subtype == "spatial_variogram":
                    em.pl.variogram(
                        adata,
                        score_keys=[score_col],  # Note: score_keys (plural) and as list
                    )
                elif params.subtype == "spatial_score":
                    # Use EnrichMap native spatial visualization
                    # Default size=0.5 for better visibility (EnrichMap default is 2)
                    spot_size = (
                        params.spot_size if params.spot_size is not None else 0.5
                    )
                    em.pl.spatial_enrichmap(
                        adata,
                        score_key=score_col,
                        library_id=library_id,
                        cmap="seismic",
                        vcenter=0,
                        size=spot_size,
                        img=False,  # Skip tissue image loading
                    )

            # Get the figure that EnrichMap created
            fig = plt.gcf()

            # Apply user-requested figure size and DPI if specified
            if params.figure_size:
                fig.set_size_inches(params.figure_size)
            if params.dpi:
                fig.set_dpi(params.dpi)

        except DataNotFoundError:
            # Re-raise data validation errors without modification
            raise
        except Exception as e:
            # Close any figure that might have been created
            plt.close("all")
            # Provide clear error message explaining why we cannot fallback
            score_info = ""
            if params.subtype == "spatial_cross_correlation":
                score_info = "Using first 2 pathways from enrichment results"
            elif params.feature:
                score_info = f"Score column: {params.feature}_score"

            raise ProcessingError(
                f"EnrichMap {params.subtype} visualization failed: {str(e)}\n\n"
                f"CONTEXT:\n"
                f"You requested a specific spatial enrichment visualization ('{params.subtype}').\n"
                f"{score_info}\n\n"
                f"SOLUTIONS:\n"
                f"1. Verify the enrichment analysis completed successfully\n"
                f"2. Check that spatial neighbors graph exists: adata.obsp['spatial_connectivities']\n"
                f"3. Ensure enrichment scores are properly stored in adata.obs\n"
                f"4. Try a different subtype: 'spatial_score', 'spatial_correlogram', 'spatial_variogram'\n\n"
                f"SCIENTIFIC INTEGRITY: Statistical visualizations (correlogram, variogram) "
                f"convey specific spatial patterns. We refuse to silently substitute them with "
                f"standard plots as this would misrepresent the analysis type."
            ) from e

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


async def _create_gsea_visualization(
    adata: ad.AnnData, params: VisualizationParameters, context=None
) -> plt.Figure:
    """Create pathway enrichment visualization (unified for ORA/GSEA and spatial EnrichMap)

    Supports multiple visualization types:
    Traditional (ORA/GSEA):
    - barplot: Top enriched pathways barplot
    - dotplot: Multi-cluster enrichment dotplot
    - heatmap: Pathway enrichment heatmap
    - enrichment_plot: Classic GSEA enrichment score plot

    Spatial EnrichMap:
    - spatial_score: Spatial distribution of enrichment scores
    - spatial_correlogram: Moran's I correlogram
    - spatial_variogram: Variogram analysis
    - spatial_cross_correlation: Cross-correlation between pathways

    Args:
        adata: AnnData object with enrichment results
        params: Visualization parameters
            - subtype: Type of plot (see above)
            - feature: Specific pathway/gene set to visualize
            - n_top_pathways: Number of top pathways to show (default: 10)
        context: MCP context

    Returns:
        Matplotlib figure with enrichment visualization
    """
    if context:
        await context.info("Creating pathway enrichment visualization")

    plot_type = params.subtype if params.subtype else "barplot"

    # Route to appropriate visualization based on subtype
    if plot_type.startswith("spatial_"):
        # Spatial EnrichMap visualizations (use per-cell scores in adata.obs)
        return await _create_enrichment_visualization(adata, params, context)
    else:
        # Traditional ORA/GSEA visualizations (use pathway statistics in adata.uns)
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
                raise DataNotFoundError(
                    f"GSEA results not found. Expected key: {gsea_key}"
                )

        gsea_results = adata.uns[gsea_key]

        if plot_type == "enrichment_plot":
            # Classic GSEA enrichment score plot
            return _create_gsea_enrichment_plot(gsea_results, params)
        elif plot_type == "barplot":
            # Top pathways barplot
            return _create_gsea_barplot(gsea_results, params)
        elif plot_type == "dotplot":
            # Multi-cluster dotplot
            return _create_gsea_dotplot(gsea_results)
        else:
            # Default to barplot
            return _create_gsea_barplot(gsea_results, params)


def _create_gsea_enrichment_plot(gsea_results, params):
    """Create classic GSEA enrichment score plot.

    Note: This plot requires running enrichment scores (RES) and hit positions
    which are typically only available from gseapy.prerank() or gseapy.gsea()
    result objects. For standard barplot/dotplot visualizations, use those
    subtypes instead.

    For this visualization to work, results must include 'RES' and 'hits' data,
    which requires storing the full gseapy result object (not just res2d DataFrame).
    """
    # Get specific pathway if specified
    pathway = params.feature if params.feature else None

    # Handle DataFrame format (standard stored format)
    if isinstance(gsea_results, pd.DataFrame):
        # DataFrame doesn't contain RES/hits data needed for enrichment plot
        # Provide helpful error message
        raise DataNotFoundError(
            "Enrichment plot requires running enrichment scores (RES) data.\n\n"
            "The stored results contain only summary statistics (DataFrame format).\n\n"
            "SOLUTIONS:\n"
            "1. Use subtype='barplot' or subtype='dotplot' instead (recommended)\n"
            "2. Re-run GSEA analysis and store the full result object\n\n"
            "Example: params={'subtype': 'barplot', 'n_top_pathways': 15}"
        )

    # Handle dict format with RES data
    if isinstance(gsea_results, dict):
        if pathway and pathway in gsea_results:
            result = gsea_results[pathway]
        else:
            # Use first available pathway
            pathway = list(gsea_results.keys())[0]
            result = gsea_results[pathway]

        # Check for required data
        if not isinstance(result, dict) or "RES" not in result:
            raise DataNotFoundError(
                f"Enrichment plot requires 'RES' (running enrichment scores) data.\n\n"
                f"Available keys: {list(result.keys()) if isinstance(result, dict) else 'N/A'}\n\n"
                "SOLUTION: Use subtype='barplot' or subtype='dotplot' instead."
            )

        # Use gseapy.gseaplot for professional visualization
        try:
            import gseapy as gp

            term = pathway
            hits = result.get("hits", result.get("hit_indices", []))
            nes = result.get("NES", result.get("nes", 0))
            pval = result.get("pval", result.get("NOM p-val", 0))
            fdr = result.get("fdr", result.get("FDR q-val", 0))
            RES = result["RES"]
            rank_metric = result.get("rank_metric")

            figsize = params.figure_size if params.figure_size else (6, 5.5)

            fig = gp.gseaplot(
                term=term,
                hits=hits,
                nes=nes,
                pval=pval,
                fdr=fdr,
                RES=RES,
                rank_metric=rank_metric,
                figsize=figsize,
                ofname=None,  # Don't save to file
            )
            return fig

        except Exception as e:
            raise ProcessingError(
                f"Failed to create GSEA enrichment plot: {e}\n\n"
                "SOLUTION: Use subtype='barplot' or subtype='dotplot' instead."
            ) from e

    raise ValueError(
        f"Unsupported GSEA results format: {type(gsea_results)}\n\n"
        "SOLUTION: Use subtype='barplot' or subtype='dotplot' instead."
    )


def _create_gsea_barplot(gsea_results, params):
    """Create barplot of top enriched pathways using gseapy.barplot.

    Uses the official gseapy visualization for professional, publication-ready output.
    """
    import gseapy as gp

    n_top = getattr(params, "n_top_pathways", 10)

    # Convert results to DataFrame if needed
    if isinstance(gsea_results, pd.DataFrame):
        df = gsea_results.copy()
    elif isinstance(gsea_results, dict):
        # Convert dict to DataFrame
        rows = []
        for pathway, data in gsea_results.items():
            if isinstance(data, dict):
                row = {"Term": pathway}
                row.update(data)
                rows.append(row)
        df = pd.DataFrame(rows)
    else:
        raise ValueError("Unsupported GSEA results format")

    if df.empty:
        raise DataNotFoundError("No enrichment results found")

    # Determine the p-value column for gseapy.barplot
    # gseapy expects: 'Adjusted P-value', 'P-value', 'NOM p-val', or 'FDR q-val'
    pval_col = "Adjusted P-value"  # default
    for col in ["Adjusted P-value", "FDR q-val", "fdr", "P-value", "NOM p-val", "pval"]:
        if col in df.columns:
            pval_col = col
            break

    # Ensure we have a 'Term' column (gseapy expects this)
    if "Term" not in df.columns:
        if "pathway" in df.columns:
            df["Term"] = df["pathway"]
        elif df.index.name or not df.index.equals(pd.RangeIndex(len(df))):
            df["Term"] = df.index
        else:
            raise DataNotFoundError(
                "No pathway/term column found in results. "
                "Expected 'Term' or 'pathway' column."
            )

    # Use gseapy.barplot for professional visualization
    figsize = params.figure_size if params.figure_size else (6, max(4, n_top * 0.4))

    try:
        ax = gp.barplot(
            df=df,
            column=pval_col,
            title=params.title or "Top Enriched Pathways",
            cutoff=1.0,  # Show all, we control via top_term
            top_term=n_top,
            figsize=figsize,
            color=params.colormap if params.colormap != "coolwarm" else "salmon",
            ofname=None,  # Don't save to file, return axes
        )

        # Get the figure from axes
        fig = ax.get_figure()
        plt.tight_layout()
        return fig

    except Exception as e:
        # Fallback to simple barplot if gseapy fails
        raise ProcessingError(
            f"gseapy.barplot failed: {e}\n\n"
            "This may be due to incompatible column names in the results DataFrame.\n"
            f"Available columns: {list(df.columns)}"
        ) from e


def _create_gsea_dotplot(gsea_results, params=None):
    """Create dotplot showing enrichment using gseapy.dotplot.

    Uses the official gseapy visualization for professional, publication-ready output.
    Supports both single-condition and multi-condition (grouped) results.
    """
    import gseapy as gp

    n_top = getattr(params, "n_top_pathways", 10) if params else 10

    # Handle nested dict structure (multi-condition)
    if isinstance(gsea_results, dict) and all(
        isinstance(v, dict) for v in gsea_results.values()
    ):
        # Convert nested dict to DataFrame with group column
        rows = []
        for condition, pathways in gsea_results.items():
            for pathway, data in pathways.items():
                if isinstance(data, dict):
                    row = {"Term": pathway, "Group": condition}
                    row.update(data)
                    rows.append(row)
        df = pd.DataFrame(rows)

        if df.empty:
            raise DataNotFoundError("No enrichment results found in nested structure")

        # Use x parameter for grouped dotplot
        x_col = "Group"
    elif isinstance(gsea_results, pd.DataFrame):
        df = gsea_results.copy()
        x_col = None  # Single condition
    elif isinstance(gsea_results, dict):
        # Single-level dict
        rows = []
        for pathway, data in gsea_results.items():
            if isinstance(data, dict):
                row = {"Term": pathway}
                row.update(data)
                rows.append(row)
        df = pd.DataFrame(rows)
        x_col = None
    else:
        raise ValueError("Unsupported GSEA results format")

    if df.empty:
        raise DataNotFoundError("No enrichment results found")

    # Ensure we have a 'Term' column
    if "Term" not in df.columns:
        if "pathway" in df.columns:
            df["Term"] = df["pathway"]
        elif df.index.name or not df.index.equals(pd.RangeIndex(len(df))):
            df["Term"] = df.index

    # Determine the p-value column
    pval_col = "Adjusted P-value"
    for col in ["Adjusted P-value", "FDR q-val", "fdr", "P-value", "NOM p-val", "pval"]:
        if col in df.columns:
            pval_col = col
            break

    figsize = params.figure_size if params and params.figure_size else (6, 8)
    cmap = params.colormap if params and params.colormap != "coolwarm" else "viridis_r"

    try:
        ax = gp.dotplot(
            df=df,
            column=pval_col,
            x=x_col,
            y="Term",
            title=params.title if params and params.title else "Pathway Enrichment",
            cutoff=1.0,  # Show all, control via top_term
            top_term=n_top,
            figsize=figsize,
            cmap=cmap,
            size=5,
            ofname=None,  # Don't save to file
        )

        fig = ax.get_figure()
        plt.tight_layout()
        return fig

    except Exception as e:
        raise ProcessingError(
            f"gseapy.dotplot failed: {e}\n\n"
            "This may be due to incompatible column names in the results DataFrame.\n"
            f"Available columns: {list(df.columns)}"
        ) from e


def _create_gsea_spatial_plot(adata, params):
    """Create spatial visualization of pathway enrichment scores"""
    # Get spatial coordinates
    x_coords, y_coords = get_spatial_coordinates(adata)

    # Get pathway to visualize
    pathway = params.feature if params.feature else None

    if not pathway:
        raise ValueError(
            "No pathway specified. Please provide a 'feature' parameter "
            "with the pathway name to visualize spatially."
        )

    # First try to get scores from adata.obs (pre-computed scores)
    if f"{pathway}_score" in adata.obs.columns:
        scores = adata.obs[f"{pathway}_score"].values
        score_col = f"{pathway}_score"
    elif pathway in adata.obs.columns:
        scores = adata.obs[pathway].values
        score_col = pathway
    # If not in obs, try to extract from GSEA results in uns
    elif "gsea_results" in adata.uns:
        gsea_results = adata.uns["gsea_results"]
        if isinstance(gsea_results, pd.DataFrame):
            # GSEA results are stored as DataFrame
            if pathway in gsea_results.index:
                # For spatial visualization, we need per-cell scores
                # GSEA doesn't provide per-cell scores, only pathway-level statistics
                raise ValueError(
                    f"Pathway '{pathway}' found in GSEA results, but GSEA analysis "
                    "does not provide per-cell enrichment scores for spatial visualization. "
                    "GSEA provides pathway-level statistics (enrichment scores, p-values) "
                    "but not cell-level scores needed for spatial mapping. "
                    "Use 'barplot', 'dotplot', or 'enrichment_plot' for GSEA results visualization."
                )
            else:
                available_pathways = list(gsea_results.index[:10])
                raise ValueError(
                    f"Pathway '{pathway}' not found in GSEA results. "
                    f"Available pathways: {', '.join(available_pathways)}..."
                )
    else:
        # No pathway scores available at all
        raise ValueError(
            "No pathway scores found in the dataset. "
            "GSEA analysis results are present but do not include per-cell spatial scores. "
            "Use 'barplot', 'dotplot', or 'enrichment_plot' for pathway enrichment visualization."
        )

    score_col = pathway

    # Create spatial plot
    figsize = params.figure_size if params.figure_size else (10, 8)
    fig, ax = plt.subplots(figsize=figsize)

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


async def _create_spatial_interaction_visualization(
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

    except ValueError:
        # Re-raise validation errors (e.g., missing lr_pairs)
        raise
    except Exception as e:
        # Provide clear error for any other failures
        raise ProcessingError(
            f"Spatial ligand-receptor interaction visualization failed: {str(e)}\n\n"
            f"CONTEXT:\n"
            f"This visualization maps spatial locations of cells expressing ligands and receptors.\n\n"
            f"COMMON CAUSES:\n"
            f"1. Ligand/receptor genes not found in dataset\n"
            f"2. No cells expressing the specified ligands/receptors above threshold\n"
            f"3. Missing spatial coordinates in adata.obsm['spatial']\n"
            f"4. Invalid lr_pairs format (should be list of (ligand, receptor) tuples)\n\n"
            f"SOLUTIONS:\n"
            f"1. Verify gene names match exactly (case-sensitive)\n"
            f"2. Check gene expression: adata[:, 'GENE_NAME'].X\n"
            f"3. Ensure spatial coordinates exist\n"
            f"4. Use valid lr_pairs: [('Ligand1', 'Receptor1'), ('Ligand2', 'Receptor2')]"
        ) from e


async def _create_batch_integration_visualization(
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
    figsize = params.figure_size if params.figure_size else (16, 12)
    fig, axes = plt.subplots(2, 2, figsize=figsize)

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
        if context:
            await context.warning(
                "UMAP coordinates not available. Run preprocessing with UMAP computation for complete visualization."
            )
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
        if context:
            await context.info(
                "Spatial coordinates not available. This is expected for non-spatial datasets."
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
# PLOT TYPE DISPATCH DICTIONARY
# ============================================================================
# This dictionary maps plot_type strings to their handler functions.
# Serves as the single source of truth for all supported visualization types.
# Adding a new plot type: 1) Create handler function, 2) Add to this dictionary
PLOT_HANDLERS = {
    # ============================================================
    # All visualization handlers are now private functions
    # Main entry point: visualize_data()
    # ============================================================
    "violin": _create_violin_visualization,
    "dotplot": _create_dotplot_visualization,
    "card_imputation": _create_card_imputation_visualization,
    "spatial_cnv": _create_spatial_cnv_visualization,
    "cnv_heatmap": _create_cnv_heatmap_visualization,
    "umap": _create_umap_visualization,
    "spatial": _create_spatial_visualization,
    "heatmap": _create_heatmap_visualization,
    "deconvolution": _create_deconvolution_visualization,
    "cell_communication": _create_cell_communication_visualization,
    "multi_gene": _create_multi_gene_visualization,
    "lr_pairs": _create_lr_pairs_visualization,
    "gene_correlation": _create_gene_correlation_visualization,
    "rna_velocity": _create_rna_velocity_visualization,
    "trajectory": _create_trajectory_visualization,
    "spatial_statistics": _create_spatial_statistics_visualization,
    "pathway_enrichment": _create_gsea_visualization,
    "spatial_interaction": _create_spatial_interaction_visualization,
    "batch_integration": _create_batch_integration_visualization,
}


# ============================================================================
# FILE PERSISTENCE FUNCTIONS
# ============================================================================


async def _regenerate_figure_for_export(
    adata: ad.AnnData,
    params: "VisualizationParameters",
    context: Optional[Context] = None,
) -> plt.Figure:
    """Regenerate a matplotlib figure from saved parameters for high-quality export.

    This is an internal helper function used by save_visualization to recreate
    figures from JSON metadata. It directly returns the matplotlib Figure object
    (instead of ImageContent) so it can be exported at arbitrary DPI/format.

    Args:
        adata: AnnData object containing the data
        params: VisualizationParameters reconstructed from saved metadata
        context: MCP context for logging

    Returns:
        Matplotlib Figure object ready for export

    Raises:
        ValueError: If plot_type is unknown
    """
    plot_type = params.plot_type

    if plot_type not in PLOT_HANDLERS:
        raise ValueError(f"Unknown plot type: {plot_type}")

    # Get the appropriate visualization function
    viz_func = PLOT_HANDLERS[plot_type]

    # Call the visualization function to get the figure
    # These functions return matplotlib Figure objects
    fig = viz_func(adata, params, context)

    return fig


async def save_visualization(
    data_id: str,
    ctx: "ToolContext",
    plot_type: str,
    subtype: Optional[str] = None,
    output_dir: str = "./outputs",
    filename: Optional[str] = None,
    format: str = "png",
    dpi: Optional[int] = None,
) -> str:
    """Save a visualization to disk at publication quality by regenerating from metadata.

    This function regenerates visualizations from stored metadata (JSON) and the original
    data, then exports at the requested quality. This approach is more secure than
    loading serialized figure objects (pickle) because:
    1. JSON metadata cannot contain executable code
    2. Regeneration uses the trusted visualization codebase
    3. All parameters are human-readable and auditable

    Supports multiple formats including vector (PDF, SVG, EPS) and raster (PNG, JPEG, TIFF)
    with publication-ready metadata.

    Args:
        data_id: Dataset ID
        ctx: ToolContext for unified data access and logging
        plot_type: Type of plot to save (e.g., 'spatial', 'deconvolution', 'spatial_statistics')
        subtype: Optional subtype for plot types with variants (e.g., 'neighborhood', 'scatterpie')
                 - For deconvolution: 'spatial_multi', 'dominant_type', 'diversity', 'stacked_bar', 'scatterpie', 'umap'
                 - For spatial_statistics: 'neighborhood', 'co_occurrence', 'ripley', 'moran', 'centrality', 'getis_ord'
        output_dir: Directory to save the file (default: ./outputs)
        filename: Custom filename (optional, auto-generated if not provided)
        format: Image format (png, jpg, jpeg, pdf, svg, eps, ps, tiff)
        dpi: DPI for raster formats (default: 300 for publication quality)
              Vector formats (PDF, SVG, EPS, PS) ignore DPI

    Returns:
        Path to the saved file

    Raises:
        DataNotFoundError: If visualization metadata not found or dataset not available
        ProcessingError: If regeneration or saving fails

    Examples:
        # Simple visualization
        save_visualization("data1", ctx, "spatial", format="pdf", dpi=300)

        # Visualization with subtype
        save_visualization("data1", ctx, "spatial_statistics", subtype="neighborhood", format="png")

        # Deconvolution scatterpie
        save_visualization("data1", ctx, "deconvolution", subtype="scatterpie", format="svg")

        # For high-res raster (PowerPoint, posters)
        save_visualization("data1", ctx, "umap", format="png", dpi=600)
    """
    try:
        # Use environment variable for output_dir if default value was passed
        if output_dir == "./outputs":
            output_dir = get_output_dir_from_config(default="./outputs")
            if output_dir != "./outputs":
                await ctx.info(
                    f"Using output directory from configuration: {output_dir}"
                )

        # Validate format
        valid_formats = ["png", "jpg", "jpeg", "pdf", "svg", "eps", "ps", "tiff"]
        if format.lower() not in valid_formats:
            raise ParameterError(
                f"Invalid format: {format}. Must be one of {valid_formats}"
            )

        # Generate cache key with subtype if provided
        if subtype:
            cache_key = f"{data_id}_{plot_type}_{subtype}"
        else:
            cache_key = f"{data_id}_{plot_type}"

        # Check if visualization exists in cache
        visualization_cache = ctx.get_visualization_cache()
        cache_key_exists = cache_key in visualization_cache

        # Set default DPI based on format
        if dpi is None:
            dpi = 300  # High quality for all formats (publication-ready)

        # For publication quality, recommend at least 300 DPI
        if dpi >= 300:
            await ctx.info(f"Using publication-quality DPI: {dpi}")

        # Create output directory using safe path handling
        # This resolves relative paths against project root (not cwd)
        # and automatically falls back to /tmp if the directory is not writable
        try:
            output_path = get_safe_output_path(
                output_dir, fallback_to_tmp=True, create_if_missing=True
            )

            # Inform user if fallback was used
            if str(output_path) != str(Path(output_dir).resolve()):
                await ctx.info(
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
            # Include subtype in filename if provided
            if subtype:
                plot_name = f"{plot_type}_{subtype}"
            else:
                plot_name = plot_type

            if dpi != 100:
                filename = f"{data_id}_{plot_name}_{dpi}dpi_{timestamp}.{format}"
            else:
                filename = f"{data_id}_{plot_name}_{timestamp}.{format}"
        else:
            # Ensure filename has correct extension
            if not filename.endswith(f".{format}"):
                filename = f"{filename}.{format}"

        # Full path for the file
        file_path = output_path / filename

        # Get cached figure object for export
        from ..utils.image_utils import (get_cached_figure,
                                         load_visualization_metadata)

        # 1. Try to get figure from in-memory cache (fast path, within same session)
        cached_fig = get_cached_figure(cache_key) if cache_key_exists else None

        # 2. If not in memory, regenerate from JSON metadata (secure approach)
        if cached_fig is None:
            metadata_path = f"/tmp/chatspatial/figures/{cache_key}.json"
            if os.path.exists(metadata_path):
                try:
                    # Load metadata
                    metadata = load_visualization_metadata(metadata_path)

                    await ctx.info(
                        "Regenerating figure from saved metadata (secure approach)"
                    )

                    # Reconstruct VisualizationParameters from saved metadata
                    from ..models.data import VisualizationParameters

                    saved_params = metadata.get("params", {})

                    # Override DPI with user's request for export
                    if dpi is not None:
                        saved_params["dpi"] = dpi

                    viz_params = VisualizationParameters(**saved_params)

                    # Regenerate the figure (this creates a new matplotlib figure)
                    # We need direct access to the figure for export, not the ImageContent
                    adata = await ctx.get_adata(data_id)
                    fig = await _regenerate_figure_for_export(
                        adata, viz_params, ctx._mcp_context
                    )
                    cached_fig = fig

                    await ctx.info(
                        f"Successfully regenerated {plot_type} visualization"
                    )

                except Exception as e:
                    await ctx.warning(
                        f"Failed to regenerate figure from metadata: {str(e)}"
                    )

        # Figure must exist (either in memory or regenerated)
        if cached_fig is None:
            # Provide helpful error message
            if not cache_key_exists:
                raise DataNotFoundError(
                    f"Visualization '{plot_type}' not found for dataset '{data_id}'. "
                    f"Please create the visualization first using visualize_data tool."
                )
            else:
                raise ProcessingError(
                    f"Cannot regenerate visualization for {cache_key}. "
                    f"Please ensure the dataset is loaded and regenerate using visualize_data tool."
                )

        # Export from cached figure with format-specific parameters
        format_info = format.upper()
        if format.lower() in ["pdf", "svg", "eps", "ps"]:
            format_info += " (vector)"
        else:
            format_info += f" ({dpi} DPI)"
        await ctx.info(f"Exporting visualization as {format_info}...")

        try:
            # Prepare save parameters
            save_params = {
                "bbox_inches": "tight",
                "facecolor": "white",
                "edgecolor": "none",
                "transparent": False,
                "pad_inches": 0.1,
            }

            # Format-specific settings
            if format.lower() == "pdf":
                import matplotlib

                # PDF metadata for publication
                save_params["dpi"] = dpi
                save_params["format"] = "pdf"
                save_params["metadata"] = {
                    "Title": f"{plot_type} visualization of {data_id}",
                    "Author": "ChatSpatial MCP",
                    "Subject": "Spatial Transcriptomics Analysis",
                    "Keywords": f"{plot_type}, {data_id}, spatial transcriptomics",
                    "Creator": "ChatSpatial with matplotlib",
                    "Producer": f"matplotlib {matplotlib.__version__}",
                }
            elif format.lower() == "svg":
                # SVG is vector, doesn't need DPI
                save_params["format"] = "svg"
            elif format.lower() in ["eps", "ps"]:
                # PostScript formats (vector)
                save_params["format"] = format.lower()
            elif format.lower() in ["png", "jpg", "jpeg", "tiff"]:
                # Raster formats need DPI
                save_params["dpi"] = dpi
                save_params["format"] = format.lower()
                if format.lower() in ["jpg", "jpeg"]:
                    save_params["pil_kwargs"] = {"quality": 95}  # High quality JPEG

            # Save the figure
            cached_fig.savefig(str(file_path), **save_params)

            # Get file size
            file_size_kb = os.path.getsize(file_path) / 1024

            size_str = (
                f"{file_size_kb:.1f} KB"
                if file_size_kb < 1024
                else f"{file_size_kb/1024:.1f} MB"
            )
            await ctx.info(f"Saved visualization to {file_path} ({size_str})")

            # Format-specific advice
            if format.lower() == "pdf":
                await ctx.info(
                    "PDF ready for journal submission (Nature/Science/Cell compatible)"
                )
            elif format.lower() == "svg":
                await ctx.info("SVG ready for web or editing in Illustrator/Inkscape")
            elif format.lower() in ["eps", "ps"]:
                await ctx.info("PostScript ready for LaTeX or professional printing")
            elif dpi >= 300:
                await ctx.info(f"High-resolution ({dpi} DPI) suitable for publication")

            return str(file_path)

        except Exception as e:
            raise ProcessingError(f"Failed to export visualization: {str(e)}") from e

    except (DataNotFoundError, ParameterError):
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to save visualization: {str(e)}")


async def export_all_visualizations(
    data_id: str,
    ctx: "ToolContext",
    output_dir: str = "./exports",
    format: str = "png",
    dpi: Optional[int] = None,
) -> List[str]:
    """Export all cached visualizations for a dataset to disk

    Args:
        data_id: Dataset ID to export visualizations for
        ctx: ToolContext for unified data access and logging
        output_dir: Directory to save files
        format: Image format (png, jpg, pdf, svg)
        dpi: DPI for saved images (default: 300 for publication quality)

    Returns:
        List of paths to saved files
    """
    try:
        visualization_cache = ctx.get_visualization_cache()

        # Strategy 1: Check session cache first (fast path)
        relevant_keys = [
            k for k in visualization_cache.keys() if k.startswith(f"{data_id}_")
        ]

        # Strategy 2: Fallback to JSON metadata files (persistent storage)
        if not relevant_keys:
            await ctx.info(
                f"Session cache empty, scanning metadata files for dataset '{data_id}'..."
            )

            # Scan metadata directory for this dataset
            import glob

            metadata_dir = "/tmp/chatspatial/figures"
            metadata_pattern = f"{metadata_dir}/{data_id}_*.json"
            metadata_files = glob.glob(metadata_pattern)

            if metadata_files:
                # Extract cache keys from metadata filenames
                relevant_keys = [
                    os.path.basename(f).replace(".json", "") for f in metadata_files
                ]

                await ctx.info(
                    f"Found {len(relevant_keys)} visualization(s) from metadata files"
                )
            else:
                await ctx.warning(
                    f"No visualizations found for dataset '{data_id}' "
                    f"(neither in session cache nor metadata files)"
                )
                return []

        saved_files = []

        for cache_key in relevant_keys:
            # Extract plot_type and subtype from cache key
            # Cache key format: {data_id}_{plot_type} or {data_id}_{plot_type}_{subtype}
            remainder = cache_key.replace(f"{data_id}_", "")

            # Known plot types that support subtypes
            known_plot_types_with_subtype = ["deconvolution", "spatial_statistics"]

            plot_type = None
            subtype = None

            # Try to match known plot types with subtypes
            for known_type in known_plot_types_with_subtype:
                if remainder.startswith(f"{known_type}_"):
                    plot_type = known_type
                    subtype = remainder[len(known_type) + 1 :]  # +1 for the underscore
                    break

            # If no match, treat the entire remainder as plot_type
            if plot_type is None:
                plot_type = remainder
                subtype = None

            try:
                saved_path = await save_visualization(
                    data_id=data_id,
                    ctx=ctx,
                    plot_type=plot_type,
                    subtype=subtype,
                    output_dir=output_dir,
                    format=format,
                    dpi=dpi,
                )
                saved_files.append(saved_path)
            except Exception as e:
                await ctx.warning(f"Failed to export {cache_key}: {str(e)}")

        await ctx.info(
            f"Exported {len(saved_files)} visualization(s) for dataset '{data_id}' to {output_dir}"
        )

        return saved_files

    except ProcessingError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to export visualizations: {str(e)}")


async def clear_visualization_cache(
    ctx: "ToolContext",
    data_id: Optional[str] = None,
) -> int:
    """Clear visualization cache to free memory

    Args:
        ctx: ToolContext for unified data access and logging
        data_id: Optional dataset ID to clear specific visualizations

    Returns:
        Number of visualizations cleared
    """
    try:
        if data_id:
            # Clear specific dataset visualizations using prefix
            cleared_count = ctx.clear_visualizations(prefix=f"{data_id}_")
            await ctx.info(
                f"Cleared {cleared_count} visualization(s) for dataset '{data_id}'"
            )
        else:
            # Clear all visualizations
            cleared_count = ctx.clear_visualizations()
            await ctx.info(f"Cleared all {cleared_count} visualization(s) from cache")

        return cleared_count

    except Exception as e:
        raise ProcessingError(f"Failed to clear cache: {str(e)}")
