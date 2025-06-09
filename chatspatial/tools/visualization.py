"""
Visualization tools for spatial transcriptomics data.
"""

from typing import Dict, Optional, Any, List, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import squidpy as sq
import anndata as ad
import traceback
import pandas as pd
import warnings
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, kendalltau
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from ..models.data import VisualizationParameters

# Import standardized image utilities
from ..utils.image_utils import fig_to_image, fig_to_base64, create_placeholder_image

# Import error handling utilities
from ..utils.error_handling import (
    SpatialMCPError, DataNotFoundError, InvalidParameterError,
    ProcessingError, DataCompatibilityError, validate_adata,
    handle_error, try_except_with_feedback
)


# Helper function to create a figure with the right size
def create_figure(figsize=(10, 8)):
    """Create a matplotlib figure with the right size and style"""
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


# Helper function to validate and prepare feature for visualization
async def validate_and_prepare_feature(
    adata,
    feature: Optional[str],
    context: Optional[Context] = None,
    default_feature: Optional[str] = None
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
                    await context.info(f"Found cell type {cell_type} in obs as {obs_key}")
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
                        await context.warning(f"Deconvolution result {feature} not found. Using default feature.")
                    return default_feature
        else:
            if context:
                await context.warning(f"Invalid feature format {feature}. Use 'deconvolution_key:cell_type'. Using default feature.")
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
                await context.warning(f"Feature {feature} not found in dataset. Using default feature.")
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
                    await context.info(f"Found deconvolution result {feature}. Showing top cell type: {top_cell_type}")
                return obs_key
            else:
                if context:
                    await context.info(f"Found deconvolution result {feature}, but could not determine cell types. Please specify a cell type using 'deconvolution_key:cell_type' format.")
                return default_feature
        else:
            if context:
                await context.info(f"Found deconvolution result {feature}. Please specify a cell type using 'deconvolution_key:cell_type' format.")
            return default_feature

    # Check if it's a regular feature
    elif feature not in adata.var_names and feature not in adata.obs.columns:
        if context:
            await context.warning(f"Feature {feature} not found in dataset. Using default feature.")
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
            return pd.DataFrame(deconv_results, index=adata.obs_names, columns=cell_types)
        else:
            # Use generic column names
            return pd.DataFrame(
                deconv_results,
                index=adata.obs_names,
                columns=[f"Cell_Type_{i}" for i in range(deconv_results.shape[1])]
            )

    # If it's already a DataFrame, return it
    elif isinstance(deconv_results, pd.DataFrame):
        return deconv_results

    # Otherwise, try to convert it
    else:
        try:
            return pd.DataFrame(deconv_results, index=adata.obs_names)
        except:
            return None


# Helper function to visualize deconvolution results
async def visualize_top_cell_types(adata, deconv_key, n_cell_types=4, context=None):
    """Visualize top cell types from deconvolution results

    Args:
        adata: AnnData object
        deconv_key: Key in adata.obsm for deconvolution results
        n_cell_types: Number of top cell types to visualize (will be limited to 6 max)
        context: MCP context for logging

    Returns:
        List of feature names for visualization
    """
    # Limit the number of cell types to avoid oversized responses
    n_cell_types = min(n_cell_types, 6)  # Limit to maximum 6 cell types

    if context:
        await context.info(f"Limiting to top {n_cell_types} cell types to avoid oversized responses")

    # Get deconvolution results as DataFrame
    deconv_df = get_deconvolution_dataframe(adata, deconv_key)
    if deconv_df is None:
        if context:
            await context.warning(f"Deconvolution result {deconv_key} not found")
        return []

    # Get top cell types by mean proportion
    top_cell_types = deconv_df.mean().sort_values(ascending=False).index[:n_cell_types]

    # Add cell type proportions to obs for visualization
    features = []
    for cell_type in top_cell_types:
        obs_key = f"{deconv_key}_{cell_type}"

        # Add to obs if not already there
        if obs_key not in adata.obs.columns:
            adata.obs[obs_key] = deconv_df[cell_type].values

        features.append(obs_key)

    if context:
        await context.info(f"Added top {len(features)} cell types to obs: {', '.join(features)}")

    return features


async def visualize_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: VisualizationParameters = VisualizationParameters(),
    context: Optional[Context] = None
) -> Image:
    """Visualize spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Visualization parameters
        context: MCP context

    Returns:
        Visualization image

    Raises:
        DataNotFoundError: If the dataset is not found
        InvalidParameterError: If parameters are invalid
        DataCompatibilityError: If data is not compatible with the visualization
        ProcessingError: If processing fails
    """
    # Validate parameters
    valid_plot_types = [
        "spatial", "umap", "heatmap", "violin",
        "deconvolution", "spatial_domains", "cell_communication",
        "trajectory", "gaston_isodepth", "spatial_analysis",
        "multi_gene", "lr_pairs", "gene_correlation"
    ]
    if params.plot_type not in valid_plot_types:
        error_msg = f"Invalid plot_type: {params.plot_type}. Must be one of {valid_plot_types}"
        if context:
            await context.warning(error_msg)
        raise InvalidParameterError(error_msg)

    if context:
        await context.info(f"Visualizing {params.plot_type} plot for dataset {data_id}")
        await context.info(f"Parameters: feature={params.feature}, colormap={params.colormap}")

        # Log deconvolution parameters if enabled
        if params.show_deconvolution:
            await context.info(f"Showing top {params.n_cell_types} cell types from deconvolution results")

    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        error_msg = f"Dataset {data_id} not found in data store"
        if context:
            await context.warning(error_msg)
        raise DataNotFoundError(error_msg)

    try:
        adata = data_store[data_id]["adata"]

        # Validate AnnData object
        validate_adata(adata, min_cells=5, min_genes=5)

        # Set matplotlib style for better visualizations
        sc.settings.set_figure_params(dpi=120, facecolor='white')

        # Check if we should show deconvolution results
        if params.show_deconvolution:
            # Find deconvolution results in obsm
            deconv_keys = [key for key in adata.obsm.keys() if key.startswith("deconvolution_")]
            if deconv_keys:
                # Use the first deconvolution result
                deconv_key = deconv_keys[0]
                if context:
                    await context.info(f"Found deconvolution result: {deconv_key}")

                # Get top cell types
                features = await visualize_top_cell_types(adata, deconv_key, params.n_cell_types, context)

                if features:
                    # Create a multi-panel figure with smaller size per panel
                    n_cols = min(3, len(features))  # Use up to 3 columns instead of 2
                    n_rows = (len(features) + n_cols - 1) // n_cols

                    # Reduce figure size from 6x6 per panel to 4x4 to decrease overall image size
                    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
                    axes = axes.flatten() if n_rows * n_cols > 1 else [axes]

                    # Plot each cell type
                    for i, feature in enumerate(features):
                        if i < len(axes):
                            ax = axes[i]
                            if 'spatial' in adata.uns and 'images' in adata.uns['spatial']:
                                # With tissue image background - use lower spot size
                                sc.pl.spatial(adata, img_key="hires", color=feature, ax=ax,
                                             show=False, spot_size=None)  # Let scanpy determine optimal spot size
                            else:
                                # Without tissue image - use smaller point size
                                sc.pl.embedding(adata, basis="spatial", color=feature, ax=ax,
                                              show=False, size=30)  # Smaller point size
                                ax.set_aspect('equal')

                            # Get cell type name from feature
                            cell_type = feature.split('_')[-1] if '_' in feature else feature
                            ax.set_title(f"{cell_type}", fontsize=10)  # Smaller font size

                            # Simplify axis labels
                            ax.tick_params(labelsize=8)  # Smaller tick labels

                    # Hide empty axes
                    for i in range(len(features), len(axes)):
                        axes[i].axis('off')

                    plt.tight_layout()

                    # Convert figure to image with lower DPI and PNG format (Claude doesn't support JPG)
                    if context:
                        await context.info(f"Created multi-panel figure with {len(features)} cell types")
                    return fig_to_image(fig, dpi=80, format="png", max_size_kb=300)

        # Create figure based on plot type
        if params.plot_type == "spatial":
            # Check if this should be a multi-gene visualization
            if params.features and params.multi_panel and len(params.features) > 1:
                if context:
                    await context.info(f"Creating multi-gene spatial plot for {len(params.features)} genes")
                # Create multi-gene visualization
                fig = await create_multi_gene_visualization(adata, params, context)
            else:
                if context:
                    await context.info(f"Creating spatial plot for {params.feature if params.feature else 'tissue'}")

                # Validate and prepare feature for visualization
                feature = await validate_and_prepare_feature(adata, params.feature, context, default_feature=None)

                # Create spatial plot
                # For 10x Visium data
                if 'spatial' in adata.uns and 'images' in adata.uns['spatial']:
                    # With tissue image background
                    if feature:
                        if feature in adata.var_names:
                            # Gene expression
                            fig = sc.pl.spatial(adata, img_key="hires", color=feature, cmap=params.colormap,
                                                show=False, return_fig=True)
                        elif feature in adata.obs.columns:
                            # Observation annotation (like clusters)
                            fig = sc.pl.spatial(adata, img_key="hires", color=feature,
                                                show=False, return_fig=True)
                    else:
                        # Just tissue image with spots
                        fig = sc.pl.spatial(adata, img_key="hires", show=False, return_fig=True)
                else:
                    # For other spatial data without tissue image
                    fig, ax = create_figure()
                    if feature:
                        if feature in adata.var_names:
                            # Gene expression
                            sc.pl.embedding(adata, basis="spatial", color=feature, cmap=params.colormap,
                                            show=False, ax=ax)
                        elif feature in adata.obs.columns:
                            # Observation annotation
                            sc.pl.embedding(adata, basis="spatial", color=feature,
                                            show=False, ax=ax)
                    else:
                        # Just spatial coordinates
                        sc.pl.embedding(adata, basis="spatial", show=False, ax=ax)
                        ax.set_aspect('equal')
                        ax.set_title("Spatial coordinates")

        elif params.plot_type == "umap":
            if context:
                await context.info(f"Creating UMAP plot for {params.feature if params.feature else 'clusters'}")

            # Check if UMAP has been computed
            if 'X_umap' not in adata.obsm:
                if context:
                    await context.warning("UMAP not found in dataset. Computing UMAP...")
                try:
                    sc.pp.neighbors(adata)
                    sc.tl.umap(adata)
                except Exception as e:
                    if context:
                        await context.warning(f"Failed to compute UMAP: {str(e)}")
                        await context.info("Creating fallback UMAP using PCA...")

                    # Create a fallback UMAP using PCA
                    from sklearn.decomposition import PCA

                    # Compute PCA if not already done
                    if 'X_pca' not in adata.obsm:
                        pca = PCA(n_components=min(50, adata.n_vars - 1))
                        X_pca = pca.fit_transform(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X)
                        adata.obsm['X_pca'] = X_pca

                    # Use PCA as UMAP for visualization
                    adata.obsm['X_umap'] = adata.obsm['X_pca'][:, :2]

            # Check if we should create multi-panel plot for multiple features
            if params.features and params.multi_panel:
                if context:
                    await context.info(f"Creating multi-panel UMAP plot for features: {params.features}")

                # Filter genes that exist in the data
                available_genes = [gene for gene in params.features if gene in adata.var_names]

                if not available_genes:
                    if context:
                        await context.warning(f"None of the specified genes found in data: {params.features}")
                    # Default to leiden clusters if available
                    default_feature = 'leiden' if 'leiden' in adata.obs.columns else None
                    feature = await validate_and_prepare_feature(adata, params.feature, context, default_feature=default_feature)

                    # Create single UMAP plot
                    if feature:
                        fig = sc.pl.umap(adata, color=feature, cmap=params.colormap if feature in adata.var_names else None,
                                        show=False, return_fig=True)
                    else:
                        fig = sc.pl.umap(adata, show=False, return_fig=True)
                else:
                    # Create multi-panel UMAP plot
                    n_features = len(available_genes)
                    if params.panel_layout:
                        n_rows, n_cols = params.panel_layout
                    else:
                        n_cols = min(3, n_features)
                        n_rows = (n_features + n_cols - 1) // n_cols

                    # Determine figure size
                    if params.figure_size:
                        figsize = params.figure_size
                    else:
                        figsize = (4*n_cols, 4*n_rows)

                    # Create figure
                    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, dpi=params.dpi)
                    axes = axes.flatten() if n_features > 1 else [axes]

                    # Set figure title
                    title = params.title or f"Multi-Gene Expression UMAP ({len(available_genes)} genes)"
                    fig.suptitle(title, fontsize=16)
                    plt.subplots_adjust(top=0.9)

                    # Plot each gene
                    for i, gene in enumerate(available_genes):
                        if i < len(axes):
                            ax = axes[i]
                            try:
                                # Get gene expression to determine color scale
                                gene_idx = adata.var_names.get_loc(gene)
                                gene_expr = adata.X[:, gene_idx]

                                if hasattr(gene_expr, 'toarray'):
                                    gene_expr = gene_expr.toarray().flatten()
                                elif hasattr(gene_expr, 'A1'):
                                    gene_expr = gene_expr.A1
                                else:
                                    gene_expr = np.array(gene_expr).flatten()

                                # Set appropriate color scale for sparse data
                                vmin = 0  # Always start from 0 for gene expression
                                vmax = max(gene_expr.max(), 0.1)  # Ensure we have some range

                                # Use percentile-based scaling for better visualization
                                if np.sum(gene_expr > 0) > 10:  # If we have enough expressing cells
                                    vmax = np.percentile(gene_expr[gene_expr > 0], 95)

                                sc.pl.umap(adata, color=gene, cmap=params.colormap,
                                          ax=ax, show=False, frameon=False,
                                          vmin=vmin, vmax=vmax,
                                          colorbar_loc='right')
                                ax.set_title(gene)
                            except Exception as e:
                                # Handle individual gene plotting errors
                                ax.text(0.5, 0.5, f'Error plotting {gene}:\n{str(e)}',
                                       ha='center', va='center', transform=ax.transAxes)
                                ax.set_title(f'{gene} (Error)')

                    # Hide empty axes
                    for i in range(n_features, len(axes)):
                        axes[i].axis('off')

                    plt.tight_layout()
            else:
                # Validate and prepare feature for visualization
                # Check if we should use cell_type from deconvolution
                if params.feature and params.feature.startswith('deconvolution_'):
                    # Check if we already have cell_type from deconvolution
                    if 'cell_type' in adata.obs:
                        if context:
                            await context.info("Using cell_type from deconvolution for UMAP visualization")
                        feature = 'cell_type'
                    else:
                        # Try to create cell_type from deconvolution
                        deconv_key = params.feature
                        if deconv_key in adata.obsm:
                            if context:
                                await context.info(f"Creating cell_type from {deconv_key} for UMAP visualization")

                            # Get cell types from uns
                            cell_types_key = f"{deconv_key}_cell_types"
                            if cell_types_key in adata.uns:
                                # Get deconvolution results as DataFrame
                                deconv_df = get_deconvolution_dataframe(adata, deconv_key)

                                if deconv_df is not None:
                                    # Determine the dominant cell type for each spot
                                    dominant_cell_types = []
                                    for i in range(deconv_df.shape[0]):
                                        row = deconv_df.iloc[i]
                                        max_idx = row.argmax()
                                        dominant_cell_types.append(deconv_df.columns[max_idx])

                                    # Add to adata.obs
                                    adata.obs['cell_type'] = dominant_cell_types

                                    # Make it categorical
                                    adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

                                    if context:
                                        await context.info(f"Created cell_type annotation with {len(adata.uns[cell_types_key])} cell types")

                                    feature = 'cell_type'
                                else:
                                    if context:
                                        await context.warning(f"Could not get deconvolution dataframe for {deconv_key}")
                                    # Default to leiden clusters if available
                                    default_feature = 'leiden' if 'leiden' in adata.obs.columns else None
                                    feature = await validate_and_prepare_feature(adata, params.feature, context, default_feature=default_feature)
                            else:
                                if context:
                                    await context.warning(f"Cell types not found for {deconv_key}")
                                # Default to leiden clusters if available
                                default_feature = 'leiden' if 'leiden' in adata.obs.columns else None
                                feature = await validate_and_prepare_feature(adata, params.feature, context, default_feature=default_feature)
                        else:
                            # Default to leiden clusters if available
                            default_feature = 'leiden' if 'leiden' in adata.obs.columns else None
                            feature = await validate_and_prepare_feature(adata, params.feature, context, default_feature=default_feature)
                else:
                    # Default to leiden clusters if available
                    default_feature = 'leiden' if 'leiden' in adata.obs.columns else None
                    feature = await validate_and_prepare_feature(adata, params.feature, context, default_feature=default_feature)

                # Create UMAP plot
                if feature:
                    fig = sc.pl.umap(adata, color=feature, cmap=params.colormap if feature in adata.var_names else None,
                                    show=False, return_fig=True)
                else:
                    fig = sc.pl.umap(adata, show=False, return_fig=True)

        elif params.plot_type == "heatmap":
            if context:
                await context.info("Creating heatmap plot")

            # For heatmap, we need to select top genes or use highly variable genes
            highly_variable = adata.var.get('highly_variable', None)
            if highly_variable is None or not highly_variable.any():
                if context:
                    await context.warning("No highly variable genes found. Computing...")
                try:
                    # Clean data before computing highly variable genes
                    # Replace infinite values with NaN, then fill with 0
                    if hasattr(adata.X, 'data'):
                        # Sparse matrix
                        adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
                    else:
                        # Dense matrix
                        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

                    sc.pp.highly_variable_genes(adata, n_top_genes=50)
                except Exception as e:
                    if context:
                        await context.warning(f"Failed to compute highly variable genes: {e}. Using top expressed genes instead.")
                    # Fallback: use top expressed genes
                    gene_means = np.array(adata.X.mean(axis=0)).flatten()
                    top_gene_indices = np.argsort(gene_means)[-50:]
                    adata.var['highly_variable'] = False
                    adata.var.iloc[top_gene_indices, adata.var.columns.get_loc('highly_variable')] = True

            # Create heatmap of top genes across groups
            if 'leiden' in adata.obs.columns:
                # Use leiden clusters for grouping
                groupby = 'leiden'

                # Limit the number of groups to avoid oversized responses
                n_groups = len(adata.obs[groupby].cat.categories)
                if n_groups > 10:
                    if context:
                        await context.warning(f"Too many groups ({n_groups}). Limiting to 10 groups.")
                    # Get the 10 largest groups
                    group_counts = adata.obs[groupby].value_counts().nlargest(10).index
                    # Subset the data to include only these groups
                    adata = adata[adata.obs[groupby].isin(group_counts)].copy()
            else:
                # No grouping
                groupby = None

            # Limit the number of genes to avoid oversized responses
            # Use a smaller number (20) instead of 50 to reduce response size
            n_genes = 20
            if context:
                await context.info(f"Using top {n_genes} highly variable genes for heatmap")

            # Get genes for heatmap
            if params.features and len(params.features) > 0:
                # Use user-specified genes
                available_genes = [gene for gene in params.features if gene in adata.var_names]
                if not available_genes:
                    if context:
                        await context.warning(f"None of specified genes found: {params.features}. Using highly variable genes.")
                    # Fall back to highly variable genes
                    available_genes = adata.var_names[adata.var.highly_variable][:n_genes].tolist()
                else:
                    # Limit to reasonable number for visualization
                    available_genes = available_genes[:n_genes]
            else:
                # Use highly variable genes
                available_genes = adata.var_names[adata.var.highly_variable][:n_genes].tolist()

            if context:
                await context.info(f"Creating heatmap with {len(available_genes)} genes: {available_genes[:5]}...")

            # Create heatmap with improved settings
            plt.figure(figsize=(max(8, len(available_genes) * 0.5), max(6, len(adata.obs[groupby].cat.categories) * 0.3) if groupby else 6))

            try:
                # Use scanpy's heatmap with better parameters
                ax_dict = sc.pl.heatmap(
                    adata,
                    var_names=available_genes,
                    groupby=groupby,
                    cmap=params.colormap,
                    show=False,
                    dendrogram=False,
                    standard_scale='var',  # Standardize genes (rows)
                    figsize=None,  # Let matplotlib handle the size
                    swap_axes=False,  # Keep genes as rows
                    vmin=None,  # Let scanpy determine range
                    vmax=None
                )
                fig = plt.gcf()

                # Improve the layout
                plt.tight_layout()

            except Exception as e:
                if context:
                    await context.warning(f"Scanpy heatmap failed: {e}. Creating custom heatmap...")

                # Fallback: create custom heatmap using seaborn
                import seaborn as sns
                import pandas as pd

                # Get expression data
                expr_data = adata[:, available_genes].X
                if hasattr(expr_data, 'toarray'):
                    expr_data = expr_data.toarray()

                # Create DataFrame
                expr_df = pd.DataFrame(
                    expr_data.T,  # Transpose so genes are rows
                    index=available_genes,
                    columns=adata.obs.index
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
                fig, ax = plt.subplots(figsize=(max(8, len(available_genes) * 0.5), max(6, len(available_genes) * 0.3)))

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
                    cbar_kws={'shrink': 0.8},
                    xticklabels=False,  # Too many cells to show labels
                    yticklabels=True
                )

                ax.set_title(f'Gene Expression Heatmap ({len(available_genes)} genes)')
                ax.set_xlabel('Cells')
                ax.set_ylabel('Genes')

                plt.tight_layout()

        elif params.plot_type == "violin":
            if context:
                await context.info(f"Creating violin plot for {params.feature if params.feature else 'top genes'}")

            # Check if feature exists for violin plot
            if params.feature and params.feature not in adata.var_names:
                if context:
                    await context.warning(f"Gene {params.feature} not found in dataset. Showing top genes instead.")
                # Use top highly variable genes
                highly_variable = adata.var.get('highly_variable', None)
                if highly_variable is None or not highly_variable.any():
                    try:
                        # Clean data before computing highly variable genes
                        if hasattr(adata.X, 'data'):
                            adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
                        else:
                            adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

                        sc.pp.highly_variable_genes(adata, n_top_genes=50)
                    except Exception:
                        # Fallback: use top expressed genes
                        gene_means = np.array(adata.X.mean(axis=0)).flatten()
                        top_gene_indices = np.argsort(gene_means)[-50:]
                        adata.var['highly_variable'] = False
                        adata.var.iloc[top_gene_indices, adata.var.columns.get_loc('highly_variable')] = True
                # Limit to 3 genes to reduce image size
                genes = adata.var_names[adata.var.highly_variable][:3]
            else:
                genes = [params.feature] if params.feature else adata.var_names[:3]  # Limit to 3 genes

            # Check if we have clusters for grouping
            if 'leiden' in adata.obs.columns:
                groupby = 'leiden'

                # Limit the number of groups to avoid oversized responses
                n_groups = len(adata.obs[groupby].cat.categories)
                if n_groups > 8:
                    if context:
                        await context.warning(f"Too many groups ({n_groups}). Limiting to 8 groups.")
                    # Get the 8 largest groups
                    group_counts = adata.obs[groupby].value_counts().nlargest(8).index
                    # Subset the data to include only these groups
                    adata = adata[adata.obs[groupby].isin(group_counts)].copy()
            else:
                groupby = None

            # Create violin plot with smaller figure size
            # Reduce figure size from (12, 10) to (8, 6) to decrease image size
            plt.figure(figsize=(8, 6))
            ax = sc.pl.violin(adata, genes, groupby=groupby, show=False,
                             jitter=0.2,  # Reduce jitter amount
                             scale='width')  # Scale violins to same width
            fig = plt.gcf()  # Get the current figure

        elif params.plot_type == "deconvolution":
            if context:
                await context.info("Creating deconvolution visualization")

            # Create deconvolution visualization
            fig = await create_deconvolution_visualization(adata, params, context)

        elif params.plot_type == "spatial_domains":
            if context:
                await context.info("Creating spatial domains visualization")

            # Create spatial domains visualization
            fig = await create_spatial_domains_visualization(adata, params, context)

        elif params.plot_type == "cell_communication":
            if context:
                await context.info("Creating cell communication visualization")

            # Create cell communication visualization
            fig = await create_cell_communication_visualization(adata, params, context)

        elif params.plot_type == "multi_gene":
            if context:
                await context.info("Creating multi-gene visualization")

            # Create multi-gene visualization
            fig = await create_multi_gene_visualization(adata, params, context)

        elif params.plot_type == "lr_pairs":
            if context:
                await context.info("Creating ligand-receptor pairs visualization")

            # Create LR pairs visualization
            fig = await create_lr_pairs_visualization(adata, params, context)

        elif params.plot_type == "gene_correlation":
            if context:
                await context.info("Creating gene correlation visualization")

            # Create gene correlation visualization
            fig = await create_gene_correlation_visualization(adata, params, context)

        else:
            # This should never happen due to parameter validation at the beginning
            error_msg = f"Unsupported plot type: {params.plot_type}"
            if context:
                await context.warning(error_msg)
            raise InvalidParameterError(error_msg)

        # Convert figure directly to Image object
        if context:
            await context.info(f"Converting {params.plot_type} figure to image...")
        return fig_to_image(fig, format="png", dpi=params.dpi)

    except Exception as e:
        # Make sure to close any open figures in case of error
        plt.close('all')

        # Log the error
        error_msg = f"Error in {params.plot_type} visualization: {str(e)}"
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {traceback.format_exc()}")

        # For image conversion errors, return a placeholder image
        if "fig_to_image" in str(e) or "convert" in str(e).lower():
            return create_placeholder_image(f"Error in {params.plot_type} visualization: {str(e)}")

        # Wrap the error in a more informative exception
        if isinstance(e, (DataNotFoundError, InvalidParameterError, DataCompatibilityError)):
            # Re-raise specific errors
            raise
        else:
            # Wrap generic errors
            raise ProcessingError(
                f"Failed to create {params.plot_type} visualization: {str(e)}",
                details={
                    "plot_type": params.plot_type,
                    "feature": params.feature,
                    "data_id": data_id,
                    "error_type": type(e).__name__,
                    "traceback": traceback.format_exc()
                }
            ) from e


async def create_deconvolution_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context = None
) -> plt.Figure:
    """Create deconvolution results visualization

    Args:
        adata: AnnData object with deconvolution results
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with deconvolution visualization
    """
    try:
        # Find deconvolution results in obsm
        deconv_keys = [key for key in adata.obsm.keys() if 'deconvolution' in key.lower() or 'proportions' in key.lower()]

        if not deconv_keys:
            # Look for individual cell type proportions in obs
            cell_type_cols = [col for col in adata.obs.columns if any(method in col.lower() for method in ['cell2location', 'spotiphy'])]

            if not cell_type_cols:
                raise DataNotFoundError("No deconvolution results found in adata.obsm or adata.obs")

            # Create proportions DataFrame from obs columns
            import pandas as pd
            # Extract base key (e.g., 'deconvolution_cell2location' from 'deconvolution_cell2location_T_cell')
            base_keys = set()
            for col in cell_type_cols:
                parts = col.split('_')
                if len(parts) >= 3:
                    base_key = '_'.join(parts[:-1])  # Remove last part (cell type name)
                    base_keys.add(base_key)

            if not base_keys:
                raise DataNotFoundError("Could not identify deconvolution result structure")

            # Use the first base key found
            base_key = list(base_keys)[0]
            method_name = base_key.split('_')[0].upper()

            # Get all columns for this base key
            relevant_cols = [col for col in cell_type_cols if col.startswith(base_key)]
            cell_types = [col.replace(f"{base_key}_", "") for col in relevant_cols]

            # Create proportions DataFrame
            proportions = pd.DataFrame(
                {cell_type: adata.obs[f"{base_key}_{cell_type}"].values for cell_type in cell_types},
                index=adata.obs.index
            )
        else:
            # Use the first deconvolution key found
            deconv_key = deconv_keys[0]
            proportions = pd.DataFrame(
                adata.obsm[deconv_key],
                index=adata.obs.index,
                columns=adata.uns.get(f"{deconv_key}_cell_types", [f"CellType_{i}" for i in range(adata.obsm[deconv_key].shape[1])])
            )
            method_name = deconv_key.split('_')[0].upper()

        # Get top cell types by mean proportion
        n_cell_types = min(params.n_cell_types, proportions.shape[1])
        top_cell_types = proportions.mean().sort_values(ascending=False).index[:n_cell_types]

        # Determine figure layout
        if params.panel_layout:
            n_rows, n_cols = params.panel_layout
        else:
            n_cols = min(3, n_cell_types)
            n_rows = (n_cell_types + n_cols - 1) // n_cols

        # Determine figure size
        if params.figure_size:
            figsize = params.figure_size
        else:
            figsize = (4*n_cols, 4*n_rows)

        # Create figure
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, dpi=params.dpi)
        axes = axes.flatten() if n_cell_types > 1 else [axes]

        # Set figure title
        title = params.title or f"{method_name} Cell Type Proportions"
        fig.suptitle(title, fontsize=16)
        plt.subplots_adjust(top=0.9)

        # Plot each cell type
        for i, cell_type in enumerate(top_cell_types):
            if i < len(axes):
                ax = axes[i]
                try:
                    # Create temporary AnnData with cell type proportions
                    temp_adata = adata.copy()
                    temp_adata.obs['proportion'] = proportions[cell_type].values

                    # Check for NaN values
                    if temp_adata.obs['proportion'].isna().any():
                        temp_adata.obs['proportion'] = temp_adata.obs['proportion'].fillna(0)

                    # Plot spatial distribution
                    if 'spatial' in temp_adata.obsm:
                        sc.pl.embedding(
                            temp_adata,
                            basis='spatial',
                            color='proportion',
                            cmap=params.colormap,
                            alpha=params.alpha,
                            size=params.spot_size,
                            ax=ax,
                            show=False,
                            frameon=params.show_axes,
                            colorbar_loc='right' if params.show_colorbar else None
                        )
                        ax.set_title(cell_type)
                        ax.invert_yaxis()
                    else:
                        # Fallback: bar plot
                        sorted_props = proportions[cell_type].sort_values(ascending=False)
                        ax.bar(range(len(sorted_props)), sorted_props.values, alpha=params.alpha)
                        ax.set_title(cell_type)
                        ax.set_xlabel('Spots (sorted)')
                        ax.set_ylabel('Proportion')

                except Exception as e:
                    # Handle individual cell type plotting errors
                    ax.text(0.5, 0.5, f'Error plotting {cell_type}:\n{str(e)}',
                           ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f'{cell_type} (Error)')

        # Hide empty axes
        for i in range(n_cell_types, len(axes)):
            axes[i].axis('off')

        plt.tight_layout()
        return fig

    except Exception as e:
        # Create error figure
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f'Error creating deconvolution visualization:\n{str(e)}',
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Deconvolution Visualization Error')
        ax.axis('off')
        return fig


async def create_spatial_domains_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context = None
) -> plt.Figure:
    """Create spatial domains visualization (placeholder)"""
    # TODO: Implement spatial domains visualization
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(0.5, 0.5, 'Spatial domains visualization\n(To be implemented)',
           ha='center', va='center', transform=ax.transAxes)
    ax.set_title('Spatial Domains')
    ax.axis('off')
    return fig


async def create_cell_communication_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context = None
) -> plt.Figure:
    """Create cell communication visualization (placeholder)"""
    # TODO: Implement cell communication visualization
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(0.5, 0.5, 'Cell communication visualization\n(To be implemented)',
           ha='center', va='center', transform=ax.transAxes)
    ax.set_title('Cell Communication')
    ax.axis('off')
    return fig


async def create_multi_gene_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context = None
) -> plt.Figure:
    """Create multi-gene visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with multi-gene visualization
    """
    try:
        # Ensure unique gene names to avoid indexing errors
        if not adata.var_names.is_unique:
            if context:
                await context.info("Making gene names unique to avoid indexing errors")
            adata.var_names_make_unique()

        # Get genes to visualize
        if params.features:
            genes = params.features
        elif params.feature:
            genes = [params.feature]
        else:
            # Use top highly variable genes as default
            if 'highly_variable' in adata.var.columns:
                genes = adata.var_names[adata.var.highly_variable][:6].tolist()
            else:
                genes = adata.var_names[:6].tolist()

        # Filter genes that exist in the data
        available_genes = [gene for gene in genes if gene in adata.var_names]

        if not available_genes:
            raise DataNotFoundError(f"None of the specified genes found in data: {genes}")

        if context:
            await context.info(f"Visualizing {len(available_genes)} genes: {available_genes}")

        # Determine figure layout
        n_genes = len(available_genes)
        if params.panel_layout:
            n_rows, n_cols = params.panel_layout
        else:
            n_cols = min(3, n_genes)
            n_rows = (n_genes + n_cols - 1) // n_cols

        # Determine figure size
        if params.figure_size:
            figsize = params.figure_size
        else:
            figsize = (4*n_cols, 4*n_rows)

        # Create figure
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, dpi=params.dpi)
        axes = axes.flatten() if n_genes > 1 else [axes]

        # Set figure title
        title = params.title or f"Multi-Gene Expression ({len(available_genes)} genes)"
        fig.suptitle(title, fontsize=16)
        plt.subplots_adjust(top=0.9)

        # Plot each gene
        for i, gene in enumerate(available_genes):
            if i < len(axes):
                ax = axes[i]
                try:
                    # Get gene expression
                    gene_expr = adata[:, gene].X
                    if hasattr(gene_expr, 'toarray'):
                        gene_expr = gene_expr.toarray().flatten()
                    elif hasattr(gene_expr, 'A1'):
                        gene_expr = gene_expr.A1
                    else:
                        gene_expr = np.array(gene_expr).flatten()

                    # Apply color scaling
                    if params.color_scale == "log":
                        gene_expr = np.log1p(gene_expr)
                    elif params.color_scale == "sqrt":
                        gene_expr = np.sqrt(gene_expr)

                    # Set color limits
                    vmin = params.vmin if params.vmin is not None else np.percentile(gene_expr, 1)
                    vmax = params.vmax if params.vmax is not None else np.percentile(gene_expr, 99)

                    # Create spatial plot
                    if 'spatial' in adata.obsm:
                        coords = adata.obsm['spatial']
                        scatter = ax.scatter(
                            coords[:, 0],
                            coords[:, 1],
                            c=gene_expr,
                            cmap=params.colormap,
                            alpha=params.alpha,
                            s=params.spot_size or 20,
                            vmin=vmin,
                            vmax=vmax
                        )

                        # Add colorbar if requested
                        if params.show_colorbar:
                            plt.colorbar(scatter, ax=ax, shrink=0.8)

                        ax.set_aspect('equal')
                        ax.invert_yaxis()

                        if params.add_gene_labels:
                            ax.set_title(gene, fontsize=12)

                        if not params.show_axes:
                            ax.set_xticks([])
                            ax.set_yticks([])
                    else:
                        # Fallback: histogram
                        ax.hist(gene_expr, bins=30, alpha=params.alpha, color='steelblue')
                        ax.set_xlabel('Expression')
                        ax.set_ylabel('Frequency')
                        if params.add_gene_labels:
                            ax.set_title(gene, fontsize=12)

                except Exception as e:
                    # Handle individual gene plotting errors
                    ax.text(0.5, 0.5, f'Error plotting {gene}:\n{str(e)}',
                           ha='center', va='center', transform=ax.transAxes)
                    if params.add_gene_labels:
                        ax.set_title(f'{gene} (Error)', fontsize=12)

        # Hide empty axes
        for i in range(n_genes, len(axes)):
            axes[i].axis('off')

        plt.tight_layout()
        return fig

    except Exception as e:
        # Create error figure
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f'Error creating multi-gene visualization:\n{str(e)}',
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Multi-Gene Visualization Error')
        ax.axis('off')
        return fig


async def create_lr_pairs_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context = None
) -> plt.Figure:
    """Create ligand-receptor pairs visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with LR pairs visualization
    """
    try:
        # Ensure unique gene names to avoid indexing errors
        if not adata.var_names.is_unique:
            if context:
                await context.info("Making gene names unique to avoid indexing errors")
            adata.var_names_make_unique()

        # Get LR pairs to visualize
        if params.lr_pairs:
            lr_pairs = params.lr_pairs
        elif params.features and len(params.features) >= 2:
            # Assume pairs of genes in features list
            lr_pairs = [(params.features[i], params.features[i+1])
                       for i in range(0, len(params.features)-1, 2)]
        else:
            # Default LR pairs for demonstration
            lr_pairs = [
                ("CCL21", "CCR7"),
                ("CCL19", "CCR7"),
                ("CXCL13", "CXCR5"),
                ("FN1", "CD79A")
            ]

        # Filter pairs where both genes exist
        available_pairs = []
        for ligand, receptor in lr_pairs:
            if ligand in adata.var_names and receptor in adata.var_names:
                available_pairs.append((ligand, receptor))

        if not available_pairs:
            raise DataNotFoundError(f"None of the specified LR pairs found in data: {lr_pairs}")

        if context:
            await context.info(f"Visualizing {len(available_pairs)} LR pairs: {available_pairs}")

        # Determine figure layout
        n_pairs = len(available_pairs)
        if params.panel_layout:
            n_rows, n_cols = params.panel_layout
        else:
            # Each pair gets 3 panels: ligand, receptor, correlation
            n_panels = n_pairs * 3
            n_cols = min(3, n_panels)
            n_rows = (n_panels + n_cols - 1) // n_cols

        # Determine figure size
        if params.figure_size:
            figsize = params.figure_size
        else:
            figsize = (4*n_cols, 4*n_rows)

        # Create figure
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, dpi=params.dpi)
        axes = axes.flatten() if n_panels > 1 else [axes]

        # Set figure title
        title = params.title or f"Ligand-Receptor Pairs ({len(available_pairs)} pairs)"
        fig.suptitle(title, fontsize=16)
        plt.subplots_adjust(top=0.9)

        # Plot each LR pair
        ax_idx = 0
        for pair_idx, (ligand, receptor) in enumerate(available_pairs):
            try:
                # Get expression data
                ligand_expr = adata[:, ligand].X
                receptor_expr = adata[:, receptor].X

                if hasattr(ligand_expr, 'toarray'):
                    ligand_expr = ligand_expr.toarray().flatten()
                    receptor_expr = receptor_expr.toarray().flatten()
                elif hasattr(ligand_expr, 'A1'):
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
                if ax_idx < len(axes) and 'spatial' in adata.obsm:
                    coords = adata.obsm['spatial']
                    ax = axes[ax_idx]

                    vmin = params.vmin if params.vmin is not None else np.percentile(ligand_expr, 1)
                    vmax = params.vmax if params.vmax is not None else np.percentile(ligand_expr, 99)

                    scatter = ax.scatter(
                        coords[:, 0], coords[:, 1],
                        c=ligand_expr,
                        cmap=params.colormap,
                        alpha=params.alpha,
                        s=params.spot_size or 20,
                        vmin=vmin, vmax=vmax
                    )

                    if params.show_colorbar:
                        plt.colorbar(scatter, ax=ax, shrink=0.8)

                    ax.set_aspect('equal')
                    ax.invert_yaxis()
                    if params.add_gene_labels:
                        ax.set_title(f"{ligand} (Ligand)", fontsize=10)
                    if not params.show_axes:
                        ax.set_xticks([])
                        ax.set_yticks([])

                    ax_idx += 1

                # Plot receptor
                if ax_idx < len(axes) and 'spatial' in adata.obsm:
                    ax = axes[ax_idx]

                    vmin = params.vmin if params.vmin is not None else np.percentile(receptor_expr, 1)
                    vmax = params.vmax if params.vmax is not None else np.percentile(receptor_expr, 99)

                    scatter = ax.scatter(
                        coords[:, 0], coords[:, 1],
                        c=receptor_expr,
                        cmap=params.colormap,
                        alpha=params.alpha,
                        s=params.spot_size or 20,
                        vmin=vmin, vmax=vmax
                    )

                    if params.show_colorbar:
                        plt.colorbar(scatter, ax=ax, shrink=0.8)

                    ax.set_aspect('equal')
                    ax.invert_yaxis()
                    if params.add_gene_labels:
                        ax.set_title(f"{receptor} (Receptor)", fontsize=10)
                    if not params.show_axes:
                        ax.set_xticks([])
                        ax.set_yticks([])

                    ax_idx += 1

                # Plot correlation
                if ax_idx < len(axes):
                    ax = axes[ax_idx]

                    # Calculate correlation
                    from scipy.stats import pearsonr, spearmanr, kendalltau

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
                        ax.set_title(f"Correlation: {corr:.3f}\np-value: {p_value:.2e}", fontsize=10)
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
                    ax.text(0.5, 0.5, f'Error plotting {ligand}-{receptor}:\n{str(e)}',
                           ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(f'{ligand}-{receptor} (Error)', fontsize=10)
                    ax_idx += 1

        # Hide empty axes
        for i in range(ax_idx, len(axes)):
            axes[i].axis('off')

        plt.tight_layout()
        return fig

    except Exception as e:
        # Create error figure
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f'Error creating LR pairs visualization:\n{str(e)}',
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('LR Pairs Visualization Error')
        ax.axis('off')
        return fig


async def create_gene_correlation_visualization(
    adata: ad.AnnData,
    params: VisualizationParameters,
    context = None
) -> plt.Figure:
    """Create gene correlation visualization

    Args:
        adata: AnnData object
        params: Visualization parameters
        context: MCP context

    Returns:
        Matplotlib figure with gene correlation visualization
    """
    try:
        # Ensure unique gene names to avoid indexing errors
        if not adata.var_names.is_unique:
            if context:
                await context.info("Making gene names unique to avoid indexing errors")
            adata.var_names_make_unique()

        # Get genes to analyze
        if params.features:
            genes = params.features
        elif params.feature:
            genes = [params.feature]
        else:
            # Use top highly variable genes as default
            if 'highly_variable' in adata.var.columns:
                genes = adata.var_names[adata.var.highly_variable][:10].tolist()
            else:
                genes = adata.var_names[:10].tolist()

        # Filter genes that exist in the data
        available_genes = [gene for gene in genes if gene in adata.var_names]

        if len(available_genes) < 2:
            raise DataNotFoundError(f"Need at least 2 genes for correlation analysis. Found: {available_genes}")

        if context:
            await context.info(f"Computing correlations for {len(available_genes)} genes: {available_genes}")

        # Get expression matrix
        expr_matrix = adata[:, available_genes].X
        if hasattr(expr_matrix, 'toarray'):
            expr_matrix = expr_matrix.toarray()

        # Apply color scaling
        if params.color_scale == "log":
            expr_matrix = np.log1p(expr_matrix)
        elif params.color_scale == "sqrt":
            expr_matrix = np.sqrt(expr_matrix)

        # Calculate correlation matrix
        from scipy.stats import pearsonr, spearmanr

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
        title = params.title or f"Gene Correlation Analysis ({params.correlation_method.title()})"
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
            cbar_kws={'shrink': 0.8}
        )
        ax1.set_title(f'{params.correlation_method.title()} Correlation')
        ax1.tick_params(axis='x', rotation=45)
        ax1.tick_params(axis='y', rotation=0)

        # P-value heatmap
        if params.show_correlation_stats:
            # Create significance mask
            sig_mask = p_value_matrix < 0.05

            sns.heatmap(
                -np.log10(p_value_matrix + 1e-10),  # -log10 p-values
                annot=True,
                cmap='Reds',
                square=True,
                xticklabels=available_genes,
                yticklabels=available_genes,
                ax=ax2,
                cbar_kws={'shrink': 0.8, 'label': '-log10(p-value)'}
            )
            ax2.set_title('Statistical Significance')
            ax2.tick_params(axis='x', rotation=45)
            ax2.tick_params(axis='y', rotation=0)

            # Add significance markers
            for i in range(n_genes):
                for j in range(n_genes):
                    if sig_mask[i, j] and i != j:
                        ax2.text(j + 0.5, i + 0.5, '*',
                               ha='center', va='center',
                               color='white', fontsize=16, fontweight='bold')
        else:
            # Show clustered correlation
            from scipy.cluster.hierarchy import linkage, dendrogram
            from scipy.spatial.distance import squareform

            # Convert correlation to distance
            distance_matrix = 1 - np.abs(corr_matrix)
            condensed_distances = squareform(distance_matrix, checks=False)

            # Perform hierarchical clustering
            linkage_matrix = linkage(condensed_distances, method='average')

            # Plot dendrogram
            dendrogram(linkage_matrix, labels=available_genes, ax=ax2, orientation='top')
            ax2.set_title('Gene Clustering')
            ax2.tick_params(axis='x', rotation=45)

        plt.tight_layout()
        return fig

    except Exception as e:
        # Create error figure
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, f'Error creating gene correlation visualization:\n{str(e)}',
               ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Gene Correlation Visualization Error')
        ax.axis('off')
        return fig
