"""
Deconvolution visualization functions for spatial transcriptomics.

This module contains:
- Cell type proportion spatial maps
- Dominant cell type visualization
- Diversity/entropy maps
- Stacked barplots
- Scatterpie plots (SPOTlight-style)
- UMAP proportion plots
- CARD imputation visualization
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Wedge
from scipy.stats import entropy

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...utils.adata_utils import (
    get_analysis_metadata_field,
    get_spatial_key,
    require_spatial_coords,
)
from ...utils.exceptions import DataNotFoundError, ParameterError
from .core import (
    DeconvolutionData,
    auto_spot_size,
    create_figure_from_params,
    get_category_colors,
    plot_spatial_feature,
    resolve_figure_size,
    setup_multi_panel_figure,
)

# =============================================================================
# Data Retrieval
# =============================================================================


def _get_available_runs(
    adata: "ad.AnnData",
) -> list[tuple[str, str]]:
    """Get available deconvolution runs as (method, analysis_key) pairs.

    Reads the ``method`` field from stored metadata to avoid brittle
    key-name parsing (analysis keys may contain reference_data_id).

    Returns:
        List of (method_name, analysis_key) tuples.
    """
    runs: list[tuple[str, str]] = []
    seen_keys: set[str] = set()

    # Primary: read from stored metadata
    for key in adata.uns:
        if key.startswith("deconvolution_") and key.endswith("_metadata"):
            analysis_key = key.removesuffix("_metadata")
            meta = adata.uns[key]
            # "method" is always stored by store_analysis_metadata;
            # skip entries without it (corrupt/foreign metadata).
            if not isinstance(meta, dict):
                continue
            method = meta.get("method")
            if method is None:
                continue
            if analysis_key not in seen_keys:
                runs.append((method, analysis_key))
                seen_keys.add(analysis_key)

    # Fallback: search obsm keys (legacy data without metadata)
    if not runs:
        for key in adata.obsm:
            if key.startswith("deconvolution_") and key not in seen_keys:
                method = key.removeprefix("deconvolution_")
                runs.append((method, key))
                seen_keys.add(key)

    return runs


async def get_deconvolution_data(
    adata: "ad.AnnData",
    method: Optional[str] = None,
    context: Optional["ToolContext"] = None,
) -> DeconvolutionData:
    """
    Unified function to retrieve deconvolution results from AnnData.

    This function consolidates all deconvolution data retrieval logic into
    a single, consistent interface. It handles:
    - Auto-detection when only one result exists
    - Explicit method specification
    - Clear error messages with solutions

    Priority for reading data:
        1. Read from stored metadata (most reliable)
        2. Fall back to key name inference (for legacy data)

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
    """
    runs = _get_available_runs(adata)
    method_names = [m for m, _k in runs]

    # Handle method specification
    if method is not None:
        # Find runs matching the requested method name
        matching = [(m, k) for m, k in runs if m == method]
        if not matching:
            raise DataNotFoundError(
                f"Deconvolution '{method}' not found. "
                f"Available: {method_names if method_names else 'None'}. "
                f"Run deconvolve_data() first."
            )
        # Use the last matching run (most recent)
        _method, analysis_name = matching[-1]
    else:
        # Auto-detect
        if not runs:
            raise DataNotFoundError(
                "No deconvolution results found. Run deconvolve_data() first."
            )

        if len(runs) > 1:
            raise ParameterError(
                f"Multiple deconvolution results: {method_names}. "
                f"Specify deconv_method parameter."
            )

        # Single result - auto-select
        method, analysis_name = runs[0]
        if context:
            await context.info(f"Auto-selected deconvolution method: {method}")

    # Try to get from stored metadata (these live under "statistics", not "parameters")
    stats = get_analysis_metadata_field(adata, analysis_name, "statistics") or {}
    proportions_key = stats.get("proportions_key")
    cell_types = stats.get("cell_types")
    dominant_type_key = stats.get("dominant_type_key")

    # Fall back to convention-based keys
    if not proportions_key:
        proportions_key = f"deconvolution_{method}"

    if proportions_key not in adata.obsm:
        raise DataNotFoundError(
            f"Proportions data '{proportions_key}' not found in adata.obsm"
        )

    # Get cell type names
    if not cell_types:
        cell_types_key = f"{proportions_key}_cell_types"
        if cell_types_key in adata.uns:
            cell_types = list(adata.uns[cell_types_key])
        else:
            # Fallback: generate generic names from shape
            n_cell_types = adata.obsm[proportions_key].shape[1]
            cell_types = [f"CellType_{i}" for i in range(n_cell_types)]
            if context:
                await context.warning("Cell type names not found. Using generic names.")

    # Check dominant type key
    if not dominant_type_key:
        dominant_type_key = f"dominant_celltype_{method}"

    if dominant_type_key not in adata.obs.columns:
        dominant_type_key = None

    # Create DataFrame
    proportions = pd.DataFrame(
        adata.obsm[proportions_key], index=adata.obs_names, columns=cell_types
    )

    return DeconvolutionData(
        proportions=proportions,
        method=method,
        cell_types=cell_types,
        proportions_key=proportions_key,
        dominant_type_key=dominant_type_key,
    )


# =============================================================================
# Visualization Functions
# =============================================================================


async def create_deconvolution_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create deconvolution results visualization.

    Routes to appropriate visualization based on params.subtype:
    - spatial_multi: Multi-panel spatial maps (default)
    - dominant_type: Dominant cell type map (CARD-style)
    - diversity: Shannon entropy diversity map
    - scatterpie: Spatial scatterpie (SPOTlight-style)
    - pie: Alias for scatterpie
    - umap: UMAP colored by proportions
    - imputation: CARD high-resolution imputation results

    Args:
        adata: AnnData object with deconvolution results
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with deconvolution visualization
    """
    viz_type = params.subtype or "spatial_multi"

    if viz_type == "dominant_type" or viz_type == "dominant":
        return await _create_dominant_celltype_map(adata, params, context)
    elif viz_type == "diversity":
        return await _create_diversity_map(adata, params, context)
    elif viz_type in ("scatterpie", "pie"):
        return await _create_scatterpie_plot(adata, params, context)
    elif viz_type == "umap":
        return await _create_umap_proportions(adata, params, context)
    elif viz_type == "spatial_multi":
        return await _create_spatial_multi_deconvolution(adata, params, context)
    elif viz_type == "imputation":
        return await _create_card_imputation(adata, params, context)
    else:
        raise ParameterError(
            f"Unknown deconvolution visualization type: {viz_type}. "
            f"Available: spatial_multi, pie, dominant, diversity, umap, imputation"
        )


async def _create_dominant_celltype_map(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create dominant cell type map (CARD-style).

    Shows the dominant cell type at each spatial location, optionally
    marking "pure" vs "mixed" spots based on proportion threshold.
    """
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Get dominant cell type (all-zero rows → "unassigned")
    prop_values = data.proportions.values
    row_sums = prop_values.sum(axis=1)
    zero_mask = row_sums == 0
    dominant_idx = prop_values.argmax(axis=1)
    dominant_types = data.proportions.columns[dominant_idx].values.copy()
    dominant_types[zero_mask] = "unassigned"
    dominant_proportions = prop_values.max(axis=1)

    # Mark pure vs mixed spots
    if params.show_mixed_spots:
        spot_categories = np.where(
            dominant_proportions >= params.min_proportion_threshold,
            dominant_types,
            "Mixed",
        )
        # Zero rows should stay "unassigned", not "Mixed"
        spot_categories[zero_mask] = "unassigned"
    else:
        spot_categories = dominant_types

    # Get spatial coordinates
    spatial_coords = require_spatial_coords(adata)

    # Calculate spot size (auto or user-specified)
    spot_size = auto_spot_size(adata, params.spot_size, basis="spatial")

    # Create figure
    fig, axes = create_figure_from_params(params, "deconvolution")
    ax = axes[0]

    # Get unique categories
    unique_categories = np.unique(spot_categories)
    n_categories = len(unique_categories)

    # Create colormap using centralized utility
    if params.show_mixed_spots and "Mixed" in unique_categories:
        cell_type_categories = [c for c in unique_categories if c != "Mixed"]
        n_cell_types = len(cell_type_categories)

        colors = get_category_colors(n_cell_types, params.colormap)
        cell_type_colors = {ct: colors[i] for i, ct in enumerate(cell_type_categories)}
        cell_type_colors["Mixed"] = (0.7, 0.7, 0.7, 1.0)

        for category in unique_categories:
            mask = spot_categories == category
            ax.scatter(
                spatial_coords[mask, 0],
                spatial_coords[mask, 1],
                c=[cell_type_colors[category]],
                s=spot_size,
                alpha=0.8 if category == "Mixed" else 1.0,
                label=category,
                edgecolors="none",
            )
    else:
        colors = get_category_colors(n_categories, params.colormap)
        color_map = {cat: colors[i] for i, cat in enumerate(unique_categories)}

        for category in unique_categories:
            mask = spot_categories == category
            ax.scatter(
                spatial_coords[mask, 0],
                spatial_coords[mask, 1],
                c=[color_map[category]],
                s=spot_size,
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
    return fig


async def _create_diversity_map(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create Shannon entropy diversity map.

    Shows cell type diversity at each spatial location using Shannon entropy.
    Higher entropy = more diverse/mixed cell types.
    Lower entropy = more homogeneous/dominated by single type.
    """
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Calculate Shannon entropy for each spot
    # All-zero rows (no deconvolution signal) get entropy = 0, not max
    prop_vals = data.proportions.values
    row_sums = prop_vals.sum(axis=1)
    zero_mask = row_sums == 0

    epsilon = 1e-10
    proportions_safe = prop_vals + epsilon
    spot_entropy = entropy(proportions_safe.T, base=2)

    # Zero-proportion rows: no information → entropy = 0
    spot_entropy[zero_mask] = 0.0

    # Normalize to [0, 1] range
    n_cell_types = data.proportions.shape[1]
    max_entropy = np.log2(n_cell_types) if n_cell_types > 1 else 1.0
    normalized_entropy = spot_entropy / max_entropy

    # Get spatial coordinates
    spatial_coords = require_spatial_coords(adata)

    # Calculate spot size (auto or user-specified)
    spot_size = auto_spot_size(adata, params.spot_size, basis="spatial")

    # Create figure
    fig, axes = create_figure_from_params(params, "deconvolution")
    ax = axes[0]

    scatter = ax.scatter(
        spatial_coords[:, 0],
        spatial_coords[:, 1],
        c=normalized_entropy,
        cmap=params.colormap or "viridis",
        s=spot_size,
        alpha=1.0,
        edgecolors="none",
    )

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Cell Type Diversity (Shannon Entropy)", rotation=270, labelpad=20)

    ax.set_xlabel("Spatial X")
    ax.set_ylabel("Spatial Y")
    ax.set_title(
        f"Cell Type Diversity Map ({data.method})\n"
        f"Shannon Entropy (0=homogeneous, 1=maximally diverse)"
    )
    ax.set_aspect("equal")

    plt.tight_layout()

    if context:
        mean_entropy = normalized_entropy.mean()
        std_entropy = normalized_entropy.std()
        high_div_pct = (normalized_entropy > 0.7).sum() / len(normalized_entropy) * 100
        low_div_pct = (normalized_entropy < 0.3).sum() / len(normalized_entropy) * 100
        await context.info(
            f"Created diversity map:\n"
            f"  Mean entropy: {mean_entropy:.3f} ± {std_entropy:.3f}\n"
            f"  High diversity (>0.7): {high_div_pct:.1f}% of spots\n"
            f"  Low diversity (<0.3): {low_div_pct:.1f}% of spots"
        )

    return fig


async def _create_scatterpie_plot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create spatial scatterpie plot (SPOTlight-style).

    Shows cell type proportions as pie charts at each spatial location.
    """
    data = await get_deconvolution_data(adata, params.deconv_method, context)
    spatial_coords = require_spatial_coords(adata)

    proportions_plot = data.proportions
    coords_plot = spatial_coords

    cell_types = proportions_plot.columns.tolist()
    n_cell_types = len(cell_types)

    # Use centralized colormap utility
    color_list = get_category_colors(n_cell_types, params.colormap)
    colors = {cell_type: color_list[i] for i, cell_type in enumerate(cell_types)}

    # Create figure
    fig, axes = create_figure_from_params(params, "deconvolution")
    ax = axes[0]

    # Calculate pie radius based on spatial scale
    coord_range = np.ptp(coords_plot, axis=0).max()
    base_radius = coord_range * 0.02
    pie_radius = base_radius * params.pie_scale

    for (x, y), (_, prop_row) in zip(
        coords_plot, proportions_plot.iterrows(), strict=True
    ):
        prop_values = prop_row.values

        if prop_values.sum() == 0:
            continue

        prop_normalized = prop_values / prop_values.sum()

        start_angle = 0
        for cell_type, proportion in zip(cell_types, prop_normalized, strict=True):
            if proportion > 0.01:
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

    x_min, x_max = coords_plot[:, 0].min(), coords_plot[:, 0].max()
    y_min, y_max = coords_plot[:, 1].min(), coords_plot[:, 1].max()
    padding = pie_radius * 2
    ax.set_xlim((x_min - padding, x_max + padding))
    ax.set_ylim((y_min - padding, y_max + padding))

    ax.set_xlabel("Spatial X")
    ax.set_ylabel("Spatial Y")
    ax.set_title(
        f"Spatial Scatterpie Plot ({data.method})\n"
        f"Cell Type Composition (pie scale: {params.pie_scale:.2f})"
    )
    ax.set_aspect("equal")

    legend_elements = [Patch(facecolor=colors[ct], label=ct) for ct in cell_types]
    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        ncol=1 if n_cell_types <= 15 else 2,
        fontsize=8,
    )

    plt.tight_layout()
    return fig


async def _create_umap_proportions(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create UMAP colored by cell type proportions.

    Shows UMAP embeddings in multi-panel format, with each panel showing
    the proportion of a specific cell type.
    """
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    if "X_umap" not in adata.obsm:
        raise DataNotFoundError(
            "UMAP coordinates not found in adata.obsm['X_umap']. "
            "Run UMAP dimensionality reduction first."
        )
    umap_coords = adata.obsm["X_umap"]

    # Calculate spot size (auto or user-specified, for UMAP basis)
    spot_size = auto_spot_size(adata, params.spot_size, basis="umap")

    # Select top cell types by mean proportion
    mean_proportions = data.proportions.mean(axis=0).sort_values(ascending=False)
    top_cell_types = mean_proportions.head(params.n_cell_types).index.tolist()

    n_panels = len(top_cell_types)
    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    # Use centralized figure size resolution
    figsize = resolve_figure_size(
        params, n_panels=n_panels, panel_width=4, panel_height=3.5
    )
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    for idx, cell_type in enumerate(top_cell_types):
        ax = axes[idx]
        prop_values = data.proportions[cell_type].values

        scatter = ax.scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c=prop_values,
            cmap=params.colormap or "viridis",
            s=spot_size // 3,  # Smaller for UMAP (similar to basic.py convention)
            alpha=0.8,
            vmin=0,
            vmax=1,
            edgecolors="none",
        )

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"{cell_type}\n(mean: {mean_proportions[cell_type]:.3f})")
        ax.set_aspect("equal")

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Proportion", rotation=270, labelpad=15, fontsize=8)

    for idx in range(n_panels, len(axes)):
        axes[idx].axis("off")

    fig.suptitle(
        f"UMAP Cell Type Proportions ({data.method})\n"
        f"Top {n_panels} cell types (out of {len(data.cell_types)})",
        fontsize=12,
        y=0.995,
    )

    plt.tight_layout()
    return fig


async def _create_spatial_multi_deconvolution(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Multi-panel spatial deconvolution visualization.

    Shows top N cell types as separate spatial plots.
    """
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    n_cell_types = min(params.n_cell_types, len(data.cell_types))
    top_cell_types = (
        data.proportions.mean().sort_values(ascending=False).index[:n_cell_types]
    )

    fig, axes = setup_multi_panel_figure(
        n_panels=len(top_cell_types),
        params=params,
        default_title=f"{data.method.upper()} Cell Type Proportions",
    )

    for i, cell_type in enumerate(top_cell_types):
        if i < len(axes):
            ax = axes[i]
            # Let errors propagate - don't silently create placeholder images
            proportions_values = data.proportions[cell_type].values

            if pd.isna(proportions_values).any():
                proportions_values = pd.Series(proportions_values).fillna(0).values

            if get_spatial_key(adata) is not None:
                plot_spatial_feature(
                    adata, values=proportions_values, ax=ax, params=params
                )
                ax.set_title(cell_type)
                ax.invert_yaxis()
            else:
                sorted_props = data.proportions[cell_type].sort_values(ascending=False)
                ax.bar(
                    range(len(sorted_props)),
                    sorted_props.values,
                    alpha=params.alpha,
                )
                ax.set_title(cell_type)
                ax.set_xlabel("Spots (sorted)")
                ax.set_ylabel("Proportion")

    fig.subplots_adjust(top=0.92, wspace=0.1, hspace=0.3, right=0.98)
    return fig


# =============================================================================
# CARD Imputation Visualization (subtype="imputation")
# =============================================================================


async def _create_card_imputation(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CARD imputation visualization.

    CARD's unique CAR model allows imputation at unmeasured locations,
    creating enhanced high-resolution spatial maps.

    Note: This is now called via deconvolution with subtype="imputation".

    Args:
        adata: AnnData object with CARD imputation results
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
        raise DataNotFoundError(
            "CARD imputation data not found. Run CARD with card_imputation=True."
        )

    # Extract imputation data
    impute_data = adata.uns["card_imputation"]
    imputed_proportions = impute_data["proportions"]
    imputed_coords = impute_data["coordinates"]

    # Determine what to visualize (normalize to string)
    feature = params.feature
    if isinstance(feature, list):
        feature = feature[0] if feature else "dominant"
    if not feature:
        feature = "dominant"

    # Create figure using centralized utility
    fig, axes = create_figure_from_params(params, "deconvolution")
    ax = axes[0]

    if feature == "dominant":
        # Show dominant cell types (all-zero rows → "unassigned")
        row_sums = imputed_proportions.sum(axis=1)
        zero_mask = row_sums == 0
        dominant_types = imputed_proportions.idxmax(axis=1)
        if zero_mask.any():
            dominant_types = dominant_types.copy()
            dominant_types[zero_mask] = "unassigned"
        unique_types = dominant_types.unique()

        # Use centralized colormap utility
        colors = get_category_colors(len(unique_types), params.colormap)
        color_map = {ct: colors[i] for i, ct in enumerate(unique_types)}
        if "unassigned" in color_map:
            color_map["unassigned"] = "#CCCCCC"
        point_colors = [color_map[ct] for ct in dominant_types]

        ax.scatter(
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

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Proportion", fontsize=12)

    else:
        raise DataNotFoundError(
            f"Feature '{feature}' not found. "
            f"Available: {list(imputed_proportions.columns)[:5]}..."
        )

    ax.set_xlabel("X coordinate", fontsize=12)
    ax.set_ylabel("Y coordinate", fontsize=12)
    ax.set_aspect("equal")
    plt.tight_layout()

    if context:
        await context.info("CARD imputation visualization created successfully")

    return fig
