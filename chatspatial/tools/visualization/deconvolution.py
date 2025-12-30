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
import seaborn as sns
from matplotlib.patches import Patch, Wedge
from scipy.stats import entropy

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.data import VisualizationParameters
from ...utils.adata_utils import get_spatial_key, require_spatial_coords
from ...utils.exceptions import DataNotFoundError, ParameterError
from .core import DeconvolutionData, plot_spatial_feature, setup_multi_panel_figure

# =============================================================================
# Data Retrieval
# =============================================================================


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
    # Find all deconvolution results in obsm
    deconv_keys = [key for key in adata.obsm.keys() if key.startswith("deconvolution_")]

    # Handle method specification
    if method is not None:
        target_key = f"deconvolution_{method}"
        if target_key not in adata.obsm:
            available = [k.replace("deconvolution_", "") for k in deconv_keys]
            raise DataNotFoundError(
                f"Deconvolution '{method}' not found. "
                f"Available: {available if available else 'None'}. "
                f"Run deconvolve_data() first."
            )
        proportions_key = target_key
    else:
        # Auto-detect
        if not deconv_keys:
            raise DataNotFoundError(
                "No deconvolution results found. Run deconvolve_data() first."
            )

        if len(deconv_keys) > 1:
            available = [k.replace("deconvolution_", "") for k in deconv_keys]
            raise ParameterError(
                f"Multiple deconvolution results: {available}. "
                f"Specify deconv_method parameter."
            )

        # Single result - auto-select
        proportions_key = deconv_keys[0]
        method = proportions_key.replace("deconvolution_", "")

        if context:
            await context.info(f"Auto-selected deconvolution method: {method}")

    # Get cell type names
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

    # Create DataFrame
    proportions = pd.DataFrame(
        adata.obsm[proportions_key], index=adata.obs_names, columns=cell_types
    )

    # Check if dominant type annotation exists
    dominant_type_key: Optional[str] = f"dominant_celltype_{method}"
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
    viz_type = params.subtype or "spatial_multi"

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
        return await _create_spatial_multi_deconvolution(adata, params, context)
    else:
        raise ParameterError(
            f"Unknown deconvolution visualization type: {viz_type}. "
            f"Available: spatial_multi, dominant_type, diversity, stacked_bar, "
            f"scatterpie, umap"
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
    spatial_coords = require_spatial_coords(adata)

    # Create figure
    figsize = params.figure_size if params.figure_size else (10, 8)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Get unique categories
    unique_categories = np.unique(spot_categories)
    n_categories = len(unique_categories)

    # Create colormap
    if params.show_mixed_spots and "Mixed" in unique_categories:
        cell_type_categories = [c for c in unique_categories if c != "Mixed"]
        n_cell_types = len(cell_type_categories)

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

        cell_type_colors["Mixed"] = (0.7, 0.7, 0.7, 1.0)

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
    epsilon = 1e-10
    proportions_safe = data.proportions.values + epsilon
    spot_entropy = entropy(proportions_safe.T, base=2)

    # Normalize to [0, 1] range
    max_entropy = np.log2(data.proportions.shape[1])
    normalized_entropy = spot_entropy / max_entropy

    # Get spatial coordinates
    spatial_coords = require_spatial_coords(adata)

    # Create figure
    figsize = params.figure_size if params.figure_size else (10, 8)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    scatter = ax.scatter(
        spatial_coords[:, 0],
        spatial_coords[:, 1],
        c=normalized_entropy,
        cmap=params.colormap or "viridis",
        s=params.spot_size or 10,
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


async def _create_stacked_barplot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create stacked barplot of cell type proportions.

    Shows cell type proportions for each spot as stacked bars.
    Spots can be sorted by dominant cell type, spatial order, or cluster.
    """
    data = await get_deconvolution_data(adata, params.deconv_method, context)

    # Limit number of spots for readability
    n_spots = len(data.proportions)
    if n_spots > params.max_spots:
        sample_indices = np.random.choice(n_spots, size=params.max_spots, replace=False)
        proportions_plot = data.proportions.iloc[sample_indices]
        if context:
            await context.warning(
                f"Sampled {params.max_spots} spots out of {n_spots} for readability."
            )
    else:
        proportions_plot = data.proportions

    # Sort spots based on sort_by parameter
    if params.sort_by == "dominant_type":
        dominant_idx = proportions_plot.values.argmax(axis=1)
        dominant_types = proportions_plot.columns[dominant_idx]
        sort_order = np.argsort(dominant_types)
    elif params.sort_by == "spatial":
        spatial_key = get_spatial_key(adata)
        if spatial_key:
            from scipy.cluster.hierarchy import dendrogram, linkage

            spatial_coords = adata.obsm[spatial_key][proportions_plot.index]
            linkage_matrix = linkage(spatial_coords, method="ward")
            dend = dendrogram(linkage_matrix, no_plot=True)
            sort_order = dend["leaves"]
        else:
            sort_order = np.arange(len(proportions_plot))
    elif params.sort_by == "cluster":
        cluster_key = params.cluster_key or "leiden"
        if cluster_key in adata.obs.columns:
            cluster_values = adata.obs.loc[proportions_plot.index, cluster_key]
            sort_order = np.argsort(cluster_values.astype(str))
        else:
            sort_order = np.arange(len(proportions_plot))
    else:
        sort_order = np.arange(len(proportions_plot))

    proportions_sorted = proportions_plot.iloc[sort_order]

    # Create figure
    figsize = params.figure_size if params.figure_size else (12, 6)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    cell_types = proportions_sorted.columns.tolist()
    n_cell_types = len(cell_types)

    if n_cell_types <= 20:
        cmap = plt.cm.get_cmap("tab20", n_cell_types)
    else:
        cmap = plt.cm.get_cmap(params.colormap or "tab20", n_cell_types)
    colors = [cmap(i) for i in range(n_cell_types)]

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

    ax.set_xlabel(params.sort_by.replace("_", " ").title())
    ax.set_ylabel("Cell Type Proportion")
    ax.set_title(
        f"Cell Type Proportions ({data.method})\n"
        f"Sorted by: {params.sort_by.replace('_', ' ').title()}"
    )
    ax.set_ylim((0, 1))
    ax.set_xlim((0, len(proportions_sorted)))
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        ncol=1 if n_cell_types <= 15 else 2,
        fontsize=8,
    )
    ax.set_xticks([])

    plt.tight_layout()
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

    if n_cell_types <= 20:
        cmap = plt.cm.get_cmap("tab20", n_cell_types)
    else:
        cmap = plt.cm.get_cmap(params.colormap or "tab20", n_cell_types)
    colors = {cell_type: cmap(i) for i, cell_type in enumerate(cell_types)}

    figsize = params.figure_size if params.figure_size else (12, 10)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Calculate pie radius based on spatial scale
    coord_range = np.ptp(coords_plot, axis=0).max()
    base_radius = coord_range * 0.02
    pie_radius = base_radius * params.pie_scale

    for idx in range(len(proportions_plot)):
        x, y = coords_plot[idx]
        prop_values = proportions_plot.iloc[idx].values

        if prop_values.sum() == 0:
            continue

        prop_normalized = prop_values / prop_values.sum()

        start_angle = 0
        for cell_type, proportion in zip(cell_types, prop_normalized, strict=False):
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

    # Select top cell types by mean proportion
    mean_proportions = data.proportions.mean(axis=0).sort_values(ascending=False)
    top_cell_types = mean_proportions.head(params.n_cell_types).index.tolist()

    n_panels = len(top_cell_types)
    ncols = min(3, n_panels)
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(ncols * 4, nrows * 3.5), squeeze=False
    )
    axes = axes.flatten()

    for idx, cell_type in enumerate(top_cell_types):
        ax = axes[idx]
        prop_values = data.proportions[cell_type].values

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

    temp_feature_key = "_deconv_viz_temp"

    for i, cell_type in enumerate(top_cell_types):
        if i < len(axes):
            ax = axes[i]
            try:
                proportions_values = data.proportions[cell_type].values

                if pd.isna(proportions_values).any():
                    proportions_values = pd.Series(proportions_values).fillna(0).values

                adata.obs[temp_feature_key] = proportions_values

                if "spatial" in adata.obsm:
                    plot_spatial_feature(
                        adata, feature=temp_feature_key, ax=ax, params=params
                    )
                    ax.set_title(cell_type)
                    ax.invert_yaxis()
                else:
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

    if temp_feature_key in adata.obs.columns:
        del adata.obs[temp_feature_key]

    fig.subplots_adjust(top=0.92, wspace=0.1, hspace=0.3, right=0.98)
    return fig


# =============================================================================
# CARD Imputation Visualization
# =============================================================================


async def create_card_imputation_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create CARD imputation visualization.

    CARD's unique CAR model allows imputation at unmeasured locations,
    creating enhanced high-resolution spatial maps.

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

    # Determine what to visualize
    feature = params.feature
    if not feature:
        feature = "dominant"

    figsize = params.figure_size if params.figure_size else (12, 10)
    fig, ax = plt.subplots(figsize=figsize)

    if feature == "dominant":
        # Show dominant cell types
        dominant_types = imputed_proportions.idxmax(axis=1)
        unique_types = dominant_types.unique()

        colors = sns.color_palette("tab20", n_colors=len(unique_types))
        color_map = {ct: colors[i] for i, ct in enumerate(unique_types)}
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
