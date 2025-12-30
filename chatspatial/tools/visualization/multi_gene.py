"""
Multi-gene and ligand-receptor visualization functions.

This module contains:
- Multi-gene spatial visualization
- Multi-gene UMAP visualization
- Ligand-receptor pairs visualization
- Gene correlation visualization
- Spatial interaction visualization
"""

from typing import TYPE_CHECKING, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ...models.data import VisualizationParameters
from ...utils.adata_utils import (
    ensure_unique_var_names,
    get_gene_expression,
    get_genes_expression,
    require_spatial_coords,
)
from ...utils.exceptions import DataNotFoundError, ParameterError, ProcessingError
from .core import (
    create_figure,
    get_validated_features,
    plot_spatial_feature,
    setup_multi_panel_figure,
)

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# Multi-Gene Visualization
# =============================================================================


async def create_multi_gene_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create multi-gene spatial visualization.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure with multi-gene spatial visualization
    """
    # Get validated features
    available_genes = await get_validated_features(
        adata, params, max_features=12, context=context
    )

    if context:
        await context.info(
            f"Visualizing {len(available_genes)} genes: {available_genes}"
        )

    # Setup multi-panel figure
    fig, axes = setup_multi_panel_figure(
        n_panels=len(available_genes),
        params=params,
        default_title="",
        use_tight_layout=False,
    )

    # Use unique temporary column name to avoid conflicts
    temp_feature_key = "multi_gene_expr_temp_viz_99_unique"

    for i, gene in enumerate(available_genes):
        if i < len(axes):
            ax = axes[i]
            try:
                # Get gene expression using unified utility
                gene_expr = get_gene_expression(adata, gene)

                # Apply color scaling
                if params.color_scale == "log":
                    gene_expr = np.log1p(gene_expr)
                elif params.color_scale == "sqrt":
                    gene_expr = np.sqrt(gene_expr)

                # Add temporary column
                adata.obs[temp_feature_key] = gene_expr

                # Create spatial plot
                if "spatial" in adata.obsm:
                    plot_spatial_feature(
                        adata,
                        ax=ax,
                        feature=temp_feature_key,
                        params=params,
                        show_colorbar=False,
                    )

                    # Set color limits
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

                    # Update colorbar limits
                    scatter = ax.collections[0] if ax.collections else None
                    if scatter:
                        scatter.set_clim(vmin, vmax)

                        # Add colorbar
                        if params.show_colorbar:
                            divider = make_axes_locatable(ax)
                            cax = divider.append_axes(
                                "right",
                                size=params.colorbar_size,
                                pad=params.colorbar_pad,
                            )
                            plt.colorbar(scatter, cax=cax)

                    ax.invert_yaxis()
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

    # Clean up temporary column
    if temp_feature_key in adata.obs:
        del adata.obs[temp_feature_key]

    # Adjust spacing
    fig.subplots_adjust(top=0.92, wspace=0.1, hspace=0.3, right=0.98)
    return fig


async def create_multi_gene_umap_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create multi-gene UMAP visualization.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure with multi-gene UMAP visualization
    """
    # Get validated features
    available_genes = await get_validated_features(
        adata, params, max_features=12, context=context
    )

    if context:
        await context.info(
            f"Creating multi-panel UMAP plot for features: {available_genes}"
        )

    # Setup multi-panel figure
    fig, axes = setup_multi_panel_figure(
        n_panels=len(available_genes),
        params=params,
        default_title="",
    )

    # Use unique temporary column name
    temp_feature_key = "umap_gene_temp_viz_99_unique"

    for i, gene in enumerate(available_genes):
        if i < len(axes):
            ax = axes[i]
            try:
                # Get gene expression using unified utility
                gene_expr = get_gene_expression(adata, gene)

                # Apply color scaling
                if params.color_scale == "log":
                    gene_expr = np.log1p(gene_expr)
                elif params.color_scale == "sqrt":
                    gene_expr = np.sqrt(gene_expr)

                # Add temporary column
                adata.obs[temp_feature_key] = gene_expr

                # Set color scale
                vmin = 0
                vmax = max(gene_expr.max(), 0.1)

                # Use percentile-based scaling for sparse data
                if np.sum(gene_expr > 0) > 10:
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
                ax.set_title(gene)

            except Exception as e:
                ax.text(
                    0.5,
                    0.5,
                    f"Error plotting {gene}:\n{str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )
                ax.set_title(f"{gene} (Error)")

    # Clean up temporary column
    if temp_feature_key in adata.obs:
        del adata.obs[temp_feature_key]

    return fig


# =============================================================================
# Ligand-Receptor Pairs Visualization
# =============================================================================


def _parse_lr_pairs(
    adata: "ad.AnnData",
    params: VisualizationParameters,
) -> List[Tuple[str, str]]:
    """Parse ligand-receptor pairs from various sources.

    Args:
        adata: AnnData object
        params: Visualization parameters

    Returns:
        List of (ligand, receptor) tuples

    Raises:
        DataNotFoundError: If no LR pairs found
    """
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
            has_special_format = any(
                "^" in str(f) or ("_" in str(f) and not str(f).startswith("_"))
                for f in feature_list
            )

            if has_special_format:
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

        # 3. Try to get from stored analysis results
        if not lr_pairs and hasattr(adata, "uns"):
            if "detected_lr_pairs" in adata.uns:
                lr_pairs = adata.uns["detected_lr_pairs"]
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

    # No hardcoded defaults - scientific integrity
    if not lr_pairs:
        raise DataNotFoundError(
            "No ligand-receptor pairs to visualize.\n\n"
            "Options:\n"
            "1. Run cell communication analysis first\n"
            "2. Specify lr_pairs parameter: lr_pairs=[('Ligand', 'Receptor')]\n"
            "3. Use LIANA format: feature=['Ligand^Receptor']\n"
            "4. Use underscore format: feature=['Ligand_Receptor']"
        )

    return lr_pairs


async def create_lr_pairs_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create ligand-receptor pairs visualization.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure with LR pairs visualization
    """
    from scipy.stats import kendalltau, pearsonr, spearmanr

    # Ensure unique gene names
    ensure_unique_var_names(adata)

    # Parse LR pairs
    lr_pairs = _parse_lr_pairs(adata, params)

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

    fig, axes = setup_multi_panel_figure(
        n_panels=n_panels,
        params=params,
        default_title=f"Ligand-Receptor Pairs ({len(available_pairs)} pairs)",
        use_tight_layout=True,
    )

    # Use unique temporary column name
    temp_feature_key = "lr_expr_temp_viz_99_unique"
    ax_idx = 0

    for _pair_idx, (ligand, receptor) in enumerate(available_pairs):
        try:
            # Get expression data using unified utility
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
            if ax_idx < len(axes) and "spatial" in adata.obsm:
                ax = axes[ax_idx]
                adata.obs[temp_feature_key] = ligand_expr
                plot_spatial_feature(
                    adata,
                    ax=ax,
                    feature=temp_feature_key,
                    params=params,
                    show_colorbar=False,
                )

                if params.show_colorbar and ax.collections:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right",
                        size=params.colorbar_size,
                        pad=params.colorbar_pad,
                    )
                    plt.colorbar(ax.collections[-1], cax=cax)

                ax.invert_yaxis()
                if params.add_gene_labels:
                    ax.set_title(f"{ligand} (Ligand)", fontsize=10)
                ax_idx += 1

            # Plot receptor
            if ax_idx < len(axes) and "spatial" in adata.obsm:
                ax = axes[ax_idx]
                adata.obs[temp_feature_key] = receptor_expr
                plot_spatial_feature(
                    adata,
                    ax=ax,
                    feature=temp_feature_key,
                    params=params,
                    show_colorbar=False,
                )

                if params.show_colorbar and ax.collections:
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes(
                        "right",
                        size=params.colorbar_size,
                        pad=params.colorbar_pad,
                    )
                    plt.colorbar(ax.collections[-1], cax=cax)

                ax.invert_yaxis()
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

                # Create scatter plot
                ax.scatter(ligand_expr, receptor_expr, alpha=params.alpha, s=20)
                ax.set_xlabel(f"{ligand} Expression")
                ax.set_ylabel(f"{receptor} Expression")

                if params.show_correlation_stats:
                    ax.set_title(
                        f"Correlation: {corr:.3f}\np-value: {p_value:.2e}",
                        fontsize=10,
                    )
                else:
                    ax.set_title(f"{ligand} vs {receptor}", fontsize=10)

                # Add trend line
                z = np.polyfit(ligand_expr, receptor_expr, 1)
                p = np.poly1d(z)
                ax.plot(ligand_expr, p(ligand_expr), "r--", alpha=0.8)

                ax_idx += 1

        except Exception as e:
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

    # Clean up temporary column
    if temp_feature_key in adata.obs:
        del adata.obs[temp_feature_key]

    # Adjust spacing
    fig.subplots_adjust(top=0.92, wspace=0.1, hspace=0.3, right=0.98)
    return fig


# =============================================================================
# Gene Correlation Visualization
# =============================================================================


async def create_gene_correlation_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create gene correlation visualization using seaborn clustermap.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure with gene correlation clustermap
    """
    import pandas as pd
    import seaborn as sns

    # Get validated genes
    available_genes = await get_validated_features(
        adata, params, max_features=10, context=context
    )

    if context:
        await context.info(
            f"Creating gene correlation visualization for {len(available_genes)} genes"
        )

    # Get expression matrix using unified utility
    expr_matrix = get_genes_expression(adata, available_genes)

    # Apply color scaling
    if params.color_scale == "log":
        expr_matrix = np.log1p(expr_matrix)
    elif params.color_scale == "sqrt":
        expr_matrix = np.sqrt(expr_matrix)

    # Create DataFrame for correlation
    expr_df = pd.DataFrame(expr_matrix, columns=available_genes)

    # Calculate correlation
    corr_df = expr_df.corr(method=params.correlation_method)

    # Use seaborn clustermap
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


# =============================================================================
# Spatial Interaction Visualization
# =============================================================================


async def create_spatial_interaction_visualization(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Create spatial visualization showing ligand-receptor interactions.

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context for logging

    Returns:
        matplotlib Figure with spatial interaction visualization
    """
    from scipy.spatial.distance import cdist

    if context:
        await context.info("Creating spatial interaction plot")

    try:
        # Get spatial coordinates
        spatial_coords = require_spatial_coords(adata)

        # Validate lr_pairs
        if not params.lr_pairs or len(params.lr_pairs) == 0:
            raise ParameterError(
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
        colors = plt.get_cmap("Set3")(np.linspace(0, 1, len(params.lr_pairs)))

        interaction_count = 0
        for i, (ligand, receptor) in enumerate(params.lr_pairs):
            color = colors[i]

            # Check if genes exist
            if ligand in adata.var_names and receptor in adata.var_names:
                # Get expression using unified utility
                ligand_expr = get_gene_expression(adata, ligand)
                receptor_expr = get_gene_expression(adata, receptor)

                # Define expression threshold
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
                        c=[color],
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
                        c=[color],
                        s=50,
                        alpha=0.7,
                        marker="^",
                        label=f"{receptor}+ (Receptor)",
                    )

                    # Draw connections
                    if len(ligand_coords) > 0 and len(receptor_coords) > 0:
                        distances = cdist(ligand_coords, receptor_coords)
                        distance_threshold = np.percentile(distances, 10)

                        ligand_indices, receptor_indices = np.where(
                            distances <= distance_threshold
                        )
                        for li, ri in zip(
                            ligand_indices[:50], receptor_indices[:50], strict=False
                        ):
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
        raise
    except Exception as e:
        raise ProcessingError(
            f"Spatial ligand-receptor interaction visualization failed: {str(e)}\n\n"
            f"Check gene names exist and spatial coordinates are available."
        ) from e
