"""
Unified feature visualization for spatial transcriptomics.

This module provides a single entry point for visualizing features
(genes, obs columns, or LR pairs) on different coordinate systems
(spatial, UMAP, PCA).

Replaces: spatial, umap, multi_gene (from basic.py and multi_gene.py)
"""

from typing import TYPE_CHECKING, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import kendalltau, pearsonr, spearmanr

from ...models.data import VisualizationParameters
from ...utils.adata_utils import (
    ensure_categorical,
    get_cluster_key,
    get_gene_expression,
)
from ...utils.compute import ensure_umap
from ...utils.exceptions import DataNotFoundError, ParameterError
from .core import (
    add_colorbar,
    auto_spot_size,
    create_figure,
    get_colormap,
    get_validated_features,
    plot_spatial_feature,
    setup_multi_panel_figure,
)

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# LR Pairs Parsing Helper
# =============================================================================


def _parse_lr_pairs_from_features(
    features: list[str],
) -> tuple[list[str], list[tuple[str, str]]]:
    """Parse LR pairs from feature list.

    Detects and extracts LR pair format ("Ligand^Receptor" or "Ligand_Receptor")
    from feature list.

    Args:
        features: List of feature names

    Returns:
        Tuple of (regular_features, lr_pairs)
        - regular_features: Features that are not LR pairs
        - lr_pairs: List of (ligand, receptor) tuples
    """
    regular_features = []
    lr_pairs = []

    for feature in features:
        if "^" in feature:
            # LIANA format: "Ligand^Receptor"
            ligand, receptor = feature.split("^", 1)
            lr_pairs.append((ligand, receptor))
        elif "_" in feature and not feature.startswith("_"):
            # Try underscore format, but only if it's clearly a pair
            parts = feature.split("_")
            if len(parts) == 2 and all(p[0].isupper() for p in parts if p):
                # Likely a gene pair (both start with uppercase)
                lr_pairs.append((parts[0], parts[1]))
            else:
                regular_features.append(feature)
        else:
            regular_features.append(feature)

    return regular_features, lr_pairs


# =============================================================================
# Unified Feature Visualization
# =============================================================================


async def create_feature_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Unified feature visualization on spatial or embedding coordinates.

    This function provides a single entry point for visualizing features
    (genes, obs columns, or LR pairs) on different coordinate systems.

    Coordinate selection via params.basis:
        - "spatial" (default): Spatial coordinates from adata.obsm['spatial']
        - "umap": UMAP embedding from adata.obsm['X_umap']
        - "pca": PCA embedding from adata.obsm['X_pca']

    Feature types:
        - Single gene: params.feature = "CD8A"
        - Multiple genes: params.feature = ["CD8A", "CD4", "FOXP3"]
        - Obs column: params.feature = "leiden" (cluster labels)
        - LR pairs: params.feature = ["CCL5^CCR5"] (auto-expands to panels)

    Args:
        adata: AnnData object with coordinates
        params: Visualization parameters
        context: Optional tool context for logging

    Returns:
        matplotlib Figure with feature visualization

    Raises:
        DataNotFoundError: If coordinates or features not found
        ParameterError: If parameters are invalid
    """
    # Determine coordinate basis (default: spatial)
    basis = params.basis or "spatial"

    # Validate basis and get coordinates
    if basis == "spatial":
        if "spatial" not in adata.obsm:
            raise DataNotFoundError(
                "Spatial coordinates not found in adata.obsm['spatial']."
            )
        coords = adata.obsm["spatial"]
    elif basis == "umap":
        if "X_umap" not in adata.obsm:
            # Try to compute UMAP
            if ensure_umap(adata) and context:
                await context.info("Computed UMAP embedding")
            else:
                raise DataNotFoundError(
                    "UMAP embedding not found. Run preprocessing with compute_umap=True."
                )
        coords = adata.obsm["X_umap"]
    elif basis == "pca":
        if "X_pca" not in adata.obsm:
            raise DataNotFoundError("PCA embedding not found. Run preprocessing first.")
        coords = adata.obsm["X_pca"][:, :2]  # First 2 PCs
    else:
        raise ParameterError(f"Invalid basis: {basis}. Use 'spatial', 'umap', or 'pca'")

    # Parse features
    if params.feature is None:
        features: list[str] = []
    elif isinstance(params.feature, list):
        features = params.feature
    else:
        features = [params.feature]

    # Default to cluster key if no features specified
    if not features:
        default_cluster = get_cluster_key(adata)
        if default_cluster:
            features = [default_cluster]
        else:
            raise ParameterError(
                "No features specified and no default clustering found"
            )

    # Parse LR pairs from features
    regular_features, lr_pairs = _parse_lr_pairs_from_features(features)

    # If we have LR pairs, use specialized LR visualization
    if lr_pairs:
        return await _create_lr_pairs_visualization(
            adata, params, context, lr_pairs, basis, coords
        )

    # Validate regular features
    validated_features = await get_validated_features(
        adata, params, max_features=12, context=context
    )

    if context:
        await context.info(
            f"Visualizing {len(validated_features)} features on {basis}: {validated_features}"
        )

    # Single feature: simple plot
    if len(validated_features) == 1:
        return await _create_single_feature_plot(
            adata, params, validated_features[0], basis, coords
        )

    # Multiple features: multi-panel plot
    return await _create_multi_feature_plot(
        adata, params, context, validated_features, basis, coords
    )


# =============================================================================
# Single Feature Plot
# =============================================================================


async def _create_single_feature_plot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    feature: str,
    basis: str,
    coords: np.ndarray,
) -> plt.Figure:
    """Create a single feature plot.

    Args:
        adata: AnnData object
        params: Visualization parameters
        feature: Feature to plot
        basis: Coordinate basis
        coords: Coordinate array

    Returns:
        matplotlib Figure
    """
    fig, ax = create_figure(params.figure_size or (10, 8))

    # Calculate spot size
    spot_size = auto_spot_size(adata, params.spot_size, basis=basis)

    # Determine if feature is a gene or obs column
    if feature in adata.var_names:
        # Gene expression
        values = get_gene_expression(adata, feature)

        # Apply color scaling
        if params.color_scale == "log":
            values = np.log1p(values)
        elif params.color_scale == "sqrt":
            values = np.sqrt(values)

        scatter = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=values,
            cmap=params.colormap,
            s=spot_size if basis == "spatial" else spot_size // 3,
            alpha=params.alpha,
            vmin=params.vmin,
            vmax=params.vmax,
        )

        if params.show_colorbar:
            add_colorbar(fig, ax, scatter, params, label=feature)

    elif feature in adata.obs.columns:
        # Observation column
        values = adata.obs[feature]
        is_categorical = (
            pd.api.types.is_categorical_dtype(values) or values.dtype == object
        )

        if is_categorical:
            ensure_categorical(adata, feature)
            categories = adata.obs[feature].cat.categories
            n_cats = len(categories)
            colors = get_colormap(params.colormap, n_colors=n_cats)

            for i, cat in enumerate(categories):
                mask = adata.obs[feature] == cat
                ax.scatter(
                    coords[mask, 0],
                    coords[mask, 1],
                    c=[colors[i]],
                    s=spot_size if basis == "spatial" else spot_size // 3,
                    alpha=params.alpha,
                    label=cat,
                )

            if params.show_legend:
                ax.legend(
                    loc="center left",
                    bbox_to_anchor=(1, 0.5),
                    fontsize=8,
                    frameon=False,
                )
        else:
            # Numeric obs column
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                cmap=params.colormap,
                s=spot_size if basis == "spatial" else spot_size // 3,
                alpha=params.alpha,
                vmin=params.vmin,
                vmax=params.vmax,
            )

            if params.show_colorbar:
                add_colorbar(fig, ax, scatter, params, label=feature)
    else:
        raise DataNotFoundError(f"Feature '{feature}' not found in genes or obs")

    # Set axis labels based on basis
    if basis == "spatial":
        ax.set_xlabel("Spatial X")
        ax.set_ylabel("Spatial Y")
        ax.invert_yaxis()  # Spatial coordinates typically need inversion
    elif basis == "umap":
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    elif basis == "pca":
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    ax.set_title(params.title or f"{feature} ({basis})")

    if not params.show_axes:
        ax.axis("off")

    plt.tight_layout()
    return fig


# =============================================================================
# Multi-Feature Plot
# =============================================================================


async def _create_multi_feature_plot(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"],
    features: list[str],
    basis: str,
    coords: np.ndarray,
) -> plt.Figure:
    """Create multi-panel feature plot.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: Tool context
        features: List of features to plot
        basis: Coordinate basis
        coords: Coordinate array

    Returns:
        matplotlib Figure with multi-panel layout
    """
    # Setup multi-panel figure
    fig, axes = setup_multi_panel_figure(
        n_panels=len(features),
        params=params,
        default_title="",
        use_tight_layout=False,
    )

    # Temporary column for expression values
    temp_key = "_feature_viz_temp_99"

    for i, feature in enumerate(features):
        if i >= len(axes):
            break

        ax = axes[i]

        if feature in adata.var_names:
            # Gene expression
            values = get_gene_expression(adata, feature)

            # Apply color scaling
            if params.color_scale == "log":
                values = np.log1p(values)
            elif params.color_scale == "sqrt":
                values = np.sqrt(values)

            # Color limits
            vmin = params.vmin if params.vmin is not None else np.percentile(values, 1)
            vmax = params.vmax if params.vmax is not None else np.percentile(values, 99)

            # UMAP-specific scaling for sparse data
            if basis == "umap" and np.sum(values > 0) > 10:
                vmax = np.percentile(values[values > 0], 95)

            if basis == "spatial":
                # Use plot_spatial_feature for spatial plots
                adata.obs[temp_key] = values
                plot_spatial_feature(
                    adata,
                    ax=ax,
                    feature=temp_key,
                    params=params,
                    show_colorbar=False,
                )

                # Update color limits
                if ax.collections:
                    scatter = ax.collections[0]
                    scatter.set_clim(vmin, vmax)

                    if params.show_colorbar:
                        divider = make_axes_locatable(ax)
                        cax = divider.append_axes(
                            "right",
                            size=params.colorbar_size,
                            pad=params.colorbar_pad,
                        )
                        plt.colorbar(scatter, cax=cax)

                ax.invert_yaxis()
            else:
                # UMAP/PCA plot
                scatter = ax.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    c=values,
                    cmap=params.colormap,
                    s=20,
                    alpha=params.alpha,
                    vmin=vmin,
                    vmax=vmax,
                )

                if params.show_colorbar:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right",
                        size=params.colorbar_size,
                        pad=params.colorbar_pad,
                    )
                    plt.colorbar(scatter, cax=cax)

        elif feature in adata.obs.columns:
            # Categorical obs column
            values = adata.obs[feature]
            is_categorical = (
                pd.api.types.is_categorical_dtype(values) or values.dtype == object
            )

            if is_categorical:
                ensure_categorical(adata, feature)
                categories = adata.obs[feature].cat.categories
                n_cats = len(categories)
                colors = get_colormap(params.colormap, n_colors=n_cats)

                for j, cat in enumerate(categories):
                    mask = adata.obs[feature] == cat
                    ax.scatter(
                        coords[mask, 0],
                        coords[mask, 1],
                        c=[colors[j]],
                        s=20,
                        alpha=params.alpha,
                        label=cat,
                    )
            else:
                scatter = ax.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    c=values,
                    cmap=params.colormap,
                    s=20,
                    alpha=params.alpha,
                )
                if params.show_colorbar:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right",
                        size=params.colorbar_size,
                        pad=params.colorbar_pad,
                    )
                    plt.colorbar(scatter, cax=cax)

            if basis == "spatial":
                ax.invert_yaxis()

        if params.add_gene_labels:
            ax.set_title(feature, fontsize=12)

        if not params.show_axes:
            ax.axis("off")

    # Clean up temporary column
    if temp_key in adata.obs:
        del adata.obs[temp_key]

    # Adjust spacing
    fig.subplots_adjust(
        top=0.92,
        wspace=params.subplot_wspace + 0.1,
        hspace=params.subplot_hspace,
        right=0.98,
    )
    return fig


# =============================================================================
# LR Pairs Visualization
# =============================================================================


async def _create_lr_pairs_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"],
    lr_pairs: list[tuple[str, str]],
    basis: str,
    coords: np.ndarray,
) -> plt.Figure:
    """Create LR pairs visualization with ligand, receptor, and correlation panels.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: Tool context
        lr_pairs: List of (ligand, receptor) tuples
        basis: Coordinate basis
        coords: Coordinate array

    Returns:
        matplotlib Figure with LR pairs visualization
    """
    # Filter pairs where both genes exist
    available_pairs = [
        (ligand, receptor)
        for ligand, receptor in lr_pairs
        if ligand in adata.var_names and receptor in adata.var_names
    ]

    if not available_pairs:
        raise DataNotFoundError(
            f"None of the specified LR pairs found in data: {lr_pairs}"
        )

    # Limit to avoid overly large plots
    max_pairs = 4
    if len(available_pairs) > max_pairs:
        if context:
            await context.warning(
                f"Too many LR pairs ({len(available_pairs)}). Limiting to {max_pairs}."
            )
        available_pairs = available_pairs[:max_pairs]

    if context:
        await context.info(
            f"Visualizing {len(available_pairs)} LR pairs on {basis}: {available_pairs}"
        )

    # Each pair gets 3 panels: ligand, receptor, correlation
    n_panels = len(available_pairs) * 3

    fig, axes = setup_multi_panel_figure(
        n_panels=n_panels,
        params=params,
        default_title=f"Ligand-Receptor Pairs ({len(available_pairs)} pairs)",
        use_tight_layout=True,
    )

    temp_key = "_lr_viz_temp_99"
    ax_idx = 0

    for ligand, receptor in available_pairs:
        # Get expression data
        ligand_expr = get_gene_expression(adata, ligand)
        receptor_expr = get_gene_expression(adata, receptor)

        # Apply color scaling
        if params.color_scale == "log":
            ligand_expr = np.log1p(ligand_expr)
            receptor_expr = np.log1p(receptor_expr)
        elif params.color_scale == "sqrt":
            ligand_expr = np.sqrt(ligand_expr)
            receptor_expr = np.sqrt(receptor_expr)

        # Plot ligand
        if ax_idx < len(axes):
            ax = axes[ax_idx]
            if basis == "spatial":
                adata.obs[temp_key] = ligand_expr
                plot_spatial_feature(
                    adata, ax=ax, feature=temp_key, params=params, show_colorbar=False
                )
                if params.show_colorbar and ax.collections:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right", size=params.colorbar_size, pad=params.colorbar_pad
                    )
                    plt.colorbar(ax.collections[-1], cax=cax)
                ax.invert_yaxis()
            else:
                scatter = ax.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    c=ligand_expr,
                    cmap=params.colormap,
                    s=20,
                    alpha=params.alpha,
                )
                if params.show_colorbar:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right", size=params.colorbar_size, pad=params.colorbar_pad
                    )
                    plt.colorbar(scatter, cax=cax)

            if params.add_gene_labels:
                ax.set_title(f"{ligand} (Ligand)", fontsize=10)
            ax_idx += 1

        # Plot receptor
        if ax_idx < len(axes):
            ax = axes[ax_idx]
            if basis == "spatial":
                adata.obs[temp_key] = receptor_expr
                plot_spatial_feature(
                    adata, ax=ax, feature=temp_key, params=params, show_colorbar=False
                )
                if params.show_colorbar and ax.collections:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right", size=params.colorbar_size, pad=params.colorbar_pad
                    )
                    plt.colorbar(ax.collections[-1], cax=cax)
                ax.invert_yaxis()
            else:
                scatter = ax.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    c=receptor_expr,
                    cmap=params.colormap,
                    s=20,
                    alpha=params.alpha,
                )
                if params.show_colorbar:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right", size=params.colorbar_size, pad=params.colorbar_pad
                    )
                    plt.colorbar(scatter, cax=cax)

            if params.add_gene_labels:
                ax.set_title(f"{receptor} (Receptor)", fontsize=10)
            ax_idx += 1

        # Plot correlation
        if ax_idx < len(axes):
            ax = axes[ax_idx]

            # Calculate correlation
            if params.correlation_method == "pearson":
                corr, p_value = pearsonr(ligand_expr, receptor_expr)
            elif params.correlation_method == "spearman":
                corr, p_value = spearmanr(ligand_expr, receptor_expr)
            else:  # kendall
                corr, p_value = kendalltau(ligand_expr, receptor_expr)

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

    # Clean up
    if temp_key in adata.obs:
        del adata.obs[temp_key]

    fig.subplots_adjust(top=0.92, wspace=0.1, hspace=0.3, right=0.98)
    return fig
