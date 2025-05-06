"""
Visualization tools for spatial transcriptomics data.
"""

from typing import Dict, Optional, Any
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import squidpy as sq
import traceback
import pandas as pd
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
    valid_plot_types = ["spatial", "umap", "heatmap", "violin"]
    if params.plot_type not in valid_plot_types:
        error_msg = f"Invalid plot_type: {params.plot_type}. Must be one of {valid_plot_types}"
        if context:
            await context.warning(error_msg)
        raise InvalidParameterError(error_msg)

    if context:
        await context.info(f"Visualizing {params.plot_type} plot for dataset {data_id}")
        await context.info(f"Parameters: feature={params.feature}, colormap={params.colormap}")

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

        # Create figure based on plot type
        if params.plot_type == "spatial":
            if context:
                await context.info(f"Creating spatial plot for {params.feature if params.feature else 'tissue'}")

            # Check if feature is provided and exists in the data
            if params.feature and params.feature not in adata.var_names and params.feature not in adata.obs.columns:
                if context:
                    await context.warning(f"Feature {params.feature} not found in dataset. Showing tissue plot instead.")
                feature = None
            else:
                feature = params.feature

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

            # Check if feature exists
            if params.feature and params.feature not in adata.var_names and params.feature not in adata.obs.columns:
                if context:
                    await context.warning(f"Feature {params.feature} not found in dataset. Showing clusters instead.")
                # Default to leiden clusters if available
                if 'leiden' in adata.obs.columns:
                    feature = 'leiden'
                else:
                    feature = None
            else:
                feature = params.feature

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
            if not adata.var.get('highly_variable', None).any():
                if context:
                    await context.warning("No highly variable genes found. Computing...")
                sc.pp.highly_variable_genes(adata, n_top_genes=50)

            # Create heatmap of top genes across groups
            if 'leiden' in adata.obs.columns:
                # Use leiden clusters for grouping
                groupby = 'leiden'
            else:
                # No grouping
                groupby = None

            # Create heatmap - note that heatmap doesn't support return_fig or ax parameter
            # We'll create a figure and let scanpy handle the axes
            plt.figure(figsize=(12, 10))
            ax_dict = sc.pl.heatmap(adata, var_names=adata.var_names[adata.var.highly_variable][:50],
                         groupby=groupby, cmap=params.colormap, show=False)
            fig = plt.gcf()  # Get the current figure

        elif params.plot_type == "violin":
            if context:
                await context.info(f"Creating violin plot for {params.feature if params.feature else 'top genes'}")

            # Check if feature exists for violin plot
            if params.feature and params.feature not in adata.var_names:
                if context:
                    await context.warning(f"Gene {params.feature} not found in dataset. Showing top genes instead.")
                # Use top highly variable genes
                if not adata.var.get('highly_variable', None).any():
                    sc.pp.highly_variable_genes(adata, n_top_genes=50)
                genes = adata.var_names[adata.var.highly_variable][:5]
            else:
                genes = [params.feature] if params.feature else adata.var_names[:5]

            # Check if we have clusters for grouping
            if 'leiden' in adata.obs.columns:
                groupby = 'leiden'
            else:
                groupby = None

            # Create violin plot - violin plot doesn't support return_fig
            # Create a figure first and let scanpy handle the axes
            plt.figure(figsize=(12, 10))
            ax = sc.pl.violin(adata, genes, groupby=groupby, show=False)
            fig = plt.gcf()  # Get the current figure

        else:
            # This should never happen due to parameter validation at the beginning
            error_msg = f"Unsupported plot type: {params.plot_type}"
            if context:
                await context.warning(error_msg)
            raise InvalidParameterError(error_msg)

        # Convert figure directly to Image object
        if context:
            await context.info(f"Converting {params.plot_type} figure to image...")
        return fig_to_image(fig, format="png")

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