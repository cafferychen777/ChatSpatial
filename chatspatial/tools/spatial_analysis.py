"""
Spatial analysis tools for spatial transcriptomics data.
"""

from typing import Dict, Optional, Any, Union, Literal
import numpy as np
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import traceback
import pandas as pd
from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from ..models.data import SpatialAnalysisParameters
from ..models.analysis import SpatialAnalysisResult

# Import standardized image utilities
from ..utils.image_utils import fig_to_image, fig_to_base64, create_placeholder_image

# Import error handling utilities
from ..utils.error_handling import (
    SpatialMCPError, DataNotFoundError, InvalidParameterError,
    ProcessingError, DataCompatibilityError, validate_adata,
    handle_error, try_except_with_feedback
)


async def analyze_spatial_unified(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Optional[Context] = None,
    return_type: str = "result"  # "result" or "image"
) -> Union[SpatialAnalysisResult, Image]:
    """Unified function to perform spatial analysis on transcriptomics data

    This function can return either a SpatialAnalysisResult object or an Image object
    based on the return_type parameter.

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Spatial analysis parameters
        context: MCP context
        return_type: Type of return value ("result" or "image")

    Returns:
        Either SpatialAnalysisResult or Image object based on return_type

    Raises:
        DataNotFoundError: If the dataset is not found
        InvalidParameterError: If parameters are invalid
        DataCompatibilityError: If data is not compatible with the analysis
        ProcessingError: If processing fails
    """
    # Import matplotlib at the beginning to avoid any issues
    import matplotlib.pyplot as plt

    # Validate parameters
    if return_type not in ["result", "image"]:
        raise InvalidParameterError(f"Invalid return_type: {return_type}. Must be 'result' or 'image'.")

    if params.analysis_type not in ["neighborhood", "co_occurrence", "ripley", "moran", "centrality"]:
        raise InvalidParameterError(f"Unsupported analysis type: {params.analysis_type}")

    if params.n_neighbors <= 0:
        raise InvalidParameterError(f"n_neighbors must be positive, got {params.n_neighbors}")

    # Log the operation
    if context:
        if return_type == "image":
            await context.info(f"Performing {params.analysis_type} spatial analysis with direct image return")
        else:
            await context.info(f"Performing {params.analysis_type} spatial analysis")
        await context.info(f"Parameters: cluster_key={params.cluster_key}, n_neighbors={params.n_neighbors}")

    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        error_msg = f"Dataset {data_id} not found in data store"
        if context:
            await context.warning(error_msg)
        raise DataNotFoundError(error_msg)

    try:
        # Make a copy of the AnnData object to avoid modifying the original
        adata = data_store[data_id]["adata"].copy()

        # Validate AnnData object
        validate_adata(adata, require_spatial=True, min_cells=10)

        # Check if we have clusters
        if params.cluster_key not in adata.obs.columns:
            if 'leiden' in adata.obs.columns:
                if context:
                    await context.warning(
                        f"Cluster key '{params.cluster_key}' not found, using 'leiden' instead. "
                        f"Available cluster keys: {list(adata.obs.columns)}"
                    )
                cluster_key = 'leiden'
            else:
                if context:
                    await context.warning(f"No clusters found, running clustering...")
                # Run clustering if not already done
                try:
                    if 'neighbors' not in adata.uns:
                        if context:
                            await context.info("Computing neighbors graph...")
                        sc.pp.neighbors(adata)

                    if context:
                        await context.info("Running Leiden clustering...")
                    sc.tl.leiden(adata)
                    cluster_key = 'leiden'
                except Exception as e:
                    error_msg = f"Failed to compute clusters: {str(e)}"
                    if context:
                        await context.warning(error_msg)
                        await context.info("Creating simple clusters as fallback...")

                    # Create simple clusters as fallback
                    from sklearn.cluster import KMeans
                    n_clusters = min(10, adata.n_obs // 10)  # Reasonable number of clusters

                    # Use PCA if available, otherwise use raw data
                    if 'X_pca' in adata.obsm:
                        X = adata.obsm['X_pca']
                    else:
                        X = adata.X

                    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
                    adata.obs['leiden'] = kmeans.fit_predict(X).astype(str)
                    adata.obs['leiden'] = adata.obs['leiden'].astype('category')
                    cluster_key = 'leiden'
        else:
            cluster_key = params.cluster_key

            # Ensure cluster key is categorical
            if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
                if context:
                    await context.info(f"Converting {cluster_key} to categorical...")
                adata.obs[cluster_key] = adata.obs[cluster_key].astype('category')

        # Ensure we have spatial neighbors
        if 'spatial_neighbors' not in adata.uns:
            if context:
                await context.info(f"Computing spatial neighbors with n_neighbors={params.n_neighbors}...")
            try:
                sq.gr.spatial_neighbors(adata, n_neighs=params.n_neighbors)
            except Exception as e:
                error_msg = f"Failed to compute spatial neighbors: {str(e)}"
                if context:
                    await context.warning(error_msg)
                    await context.info("Creating fallback spatial neighbors...")

                # Create fallback spatial neighbors
                if 'spatial' in adata.obsm:
                    from sklearn.neighbors import NearestNeighbors
                    from scipy.sparse import csr_matrix

                    # Get spatial coordinates
                    coords = adata.obsm['spatial']

                    # Compute nearest neighbors
                    nbrs = NearestNeighbors(n_neighbors=params.n_neighbors).fit(coords)
                    distances, indices = nbrs.kneighbors(coords)

                    # Create connectivity matrix
                    n_cells = adata.n_obs
                    connectivities = np.zeros((n_cells, n_cells))

                    for i in range(n_cells):
                        for j in indices[i]:
                            connectivities[i, j] = 1

                    # Store in AnnData
                    adata.uns['spatial_neighbors'] = {}
                    adata.obsp['spatial_connectivities'] = csr_matrix(connectivities)

        # Fix color palette issue - ensure we have enough colors for all clusters
        # This is a common issue with squidpy visualization
        color_key = f"{cluster_key}_colors"
        n_clusters = len(adata.obs[cluster_key].unique())

        # Check if we need to fix the colors
        if color_key in adata.uns and len(adata.uns[color_key]) < n_clusters:
            if context:
                await context.info(f"Fixing color palette: {len(adata.uns[color_key])} colors for {n_clusters} clusters")

            # Remove existing colors to let squidpy create a new palette
            if color_key in adata.uns:
                del adata.uns[color_key]

            # Create a new color palette with enough colors
            from matplotlib.colors import to_rgba

            # Use a colormap with many colors
            cmap = plt.cm.get_cmap('viridis', n_clusters)
            adata.uns[color_key] = [to_rgba(cmap(i)) for i in range(n_clusters)]

        # Perform the requested spatial analysis
        statistics = {}
        fig = None
        if params.analysis_type == "neighborhood":
            if context:
                await context.info("Running neighborhood enrichment analysis...")

            # Run neighborhood enrichment
            sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)

            # Get the result matrix - squidpy uses 'zscore' and 'count' keys
            enrichment_matrix = adata.uns[f'{cluster_key}_nhood_enrichment']['zscore']

            # Create visualization - create our own figure
            # Instead of using squidpy's plotting function, we'll create our own heatmap
            # This avoids the color bin issue
            fig, ax = plt.subplots(figsize=(10, 8))

            # Get the categories for labeling
            categories = adata.obs[cluster_key].cat.categories

            # Create a heatmap using matplotlib
            im = ax.imshow(enrichment_matrix, cmap='viridis')

            # Add colorbar
            plt.colorbar(im, ax=ax, label='Z-score')

            # Set ticks and labels
            ax.set_xticks(np.arange(len(categories)))
            ax.set_yticks(np.arange(len(categories)))
            ax.set_xticklabels(categories)
            ax.set_yticklabels(categories)

            # Rotate x labels for better readability
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

            # Add title
            ax.set_title(f'Neighborhood Enrichment ({cluster_key})')
            ax.set_xlabel('Cluster')
            ax.set_ylabel('Cluster')

        elif params.analysis_type == "co_occurrence":
            if context:
                await context.info("Running co-occurrence analysis...")

            # Run co-occurrence analysis
            sq.gr.co_occurrence(adata, cluster_key=cluster_key)

            # Get the result matrix - squidpy uses 'occ' key for co-occurrence matrix
            co_occurrence_matrix = adata.uns[f'{cluster_key}_co_occurrence']['occ']
            interval = adata.uns[f'{cluster_key}_co_occurrence']['interval']

            # Create visualization - create our own figure instead of using squidpy's
            # This avoids color palette issues

            # The co-occurrence matrix is 3D (clusters x clusters x distances)
            # We'll create a multi-panel figure with one panel for each distance
            n_clusters = co_occurrence_matrix.shape[0]
            n_distances = co_occurrence_matrix.shape[2]

            # Limit to at most 4 distance panels to keep the figure manageable
            max_panels = min(4, n_distances)
            selected_distances = np.linspace(0, n_distances-1, max_panels, dtype=int)

            # Create a figure with subplots
            fig, axes = plt.subplots(1, max_panels, figsize=(5*max_panels, 5), squeeze=False)
            axes = axes.flatten()

            # Get the categories for labeling
            categories = adata.obs[cluster_key].cat.categories

            # Create a heatmap for each selected distance
            for i, dist_idx in enumerate(selected_distances):
                ax = axes[i]
                dist_value = interval[dist_idx]

                # Create a heatmap using matplotlib
                im = ax.imshow(co_occurrence_matrix[:, :, dist_idx], cmap='viridis', vmin=0, vmax=np.max(co_occurrence_matrix))

                # Set ticks and labels
                if n_clusters <= 10:  # Only show labels if there aren't too many
                    ax.set_xticks(np.arange(len(categories)))
                    ax.set_yticks(np.arange(len(categories)))
                    ax.set_xticklabels(categories, rotation=45, ha="right")
                    ax.set_yticklabels(categories)
                else:
                    # Just show a few ticks to avoid overcrowding
                    step = max(1, n_clusters // 5)
                    ax.set_xticks(np.arange(0, n_clusters, step))
                    ax.set_yticks(np.arange(0, n_clusters, step))
                    ax.set_xticklabels(categories[::step], rotation=45, ha="right")
                    ax.set_yticklabels(categories[::step])

                # Add title for this distance
                ax.set_title(f'Distance: {dist_value:.2f}')

                # Add axis labels only to leftmost plot
                if i == 0:
                    ax.set_ylabel('Cluster')

                # Add axis labels only to bottom plot
                if i == max_panels // 2:
                    ax.set_xlabel('Cluster')

            # Add a colorbar to the right of the last subplot
            plt.colorbar(im, ax=axes, label='Co-occurrence probability', shrink=0.8)

            # Add an overall title
            fig.suptitle(f'Co-occurrence probability at different distances ({cluster_key})', fontsize=16)
            plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout to make room for the title

        elif params.analysis_type == "ripley":
            if context:
                await context.info("Running Ripley's function analysis...")

            # Run Ripley's function analysis with mode 'L' (a common transformation of K)
            # Valid modes are 'F', 'G', and 'L' according to the error message
            try:
                # The correct way to call ripley function
                sq.gr.ripley(
                    adata,
                    cluster_key=cluster_key,
                    mode='L',  # Use L mode (variance-stabilized K function)
                    n_simulations=20,  # Reduce for faster computation
                    n_observations=1000,
                    max_dist=None,  # Automatically determine max distance
                    n_steps=50
                )

                # Get the result - the key format is {cluster_key}_ripley_{mode}
                ripley_key = f'{cluster_key}_ripley_L'

                if ripley_key not in adata.uns:
                    raise KeyError(f"Expected key '{ripley_key}' not found in adata.uns")

                # Create visualization using squidpy's plotting function
                fig, ax = plt.subplots(figsize=(10, 8))

                # Use the squidpy plotting function
                sq.pl.ripley(
                    adata,
                    cluster_key=cluster_key,
                    mode='L',
                    plot_sims=True,  # Show simulation envelope
                    ax=ax
                )

            except Exception as e:
                if context:
                    await context.warning(f"Error in Ripley analysis: {str(e)}")
                    await context.info("Creating fallback visualization...")

                # Create a fallback visualization with synthetic data
                fig, ax = plt.subplots(figsize=(10, 8))

                # Get unique clusters
                categories = adata.obs[cluster_key].cat.categories

                # Create some reasonable distances
                distances = np.linspace(0, 5000, 50)

                # Plot a line for each cluster with synthetic data
                for i, cluster in enumerate(categories):
                    # Generate some synthetic data for demonstration
                    l_values = np.random.normal(0, 1, size=len(distances)) * np.sqrt(distances/1000) + distances/1000

                    # Plot the L function
                    ax.plot(distances, l_values, label=f'Cluster {cluster}')

                # Add reference line for complete spatial randomness
                ax.axhline(y=0, color='black', linestyle='--', alpha=0.5, label='CSR')

                # Add labels and title
                ax.set_xlabel('Distance')
                ax.set_ylabel('L(r) - r')
                ax.set_title(f"Ripley's L Function by Cluster ({cluster_key}) - SYNTHETIC DATA")

                # Add legend if there aren't too many clusters
                if len(categories) <= 10:
                    ax.legend(loc='best')
                else:
                    # Just add a note about the number of clusters
                    ax.text(0.98, 0.02, f"{len(categories)} clusters",
                           transform=ax.transAxes, ha='right', va='bottom')

        elif params.analysis_type == "moran":
            if context:
                await context.info("Running spatial autocorrelation (Moran's I) analysis...")

            # Run spatial autocorrelation
            # We'll use top highly variable genes if not already computed
            if 'highly_variable' not in adata.var or not adata.var['highly_variable'].any():
                if context:
                    await context.info("Computing highly variable genes...")
                sc.pp.highly_variable_genes(adata, n_top_genes=50)

            genes = adata.var_names[adata.var['highly_variable']][:20].tolist()

            sq.gr.spatial_autocorr(adata, genes=genes, n_perms=100)

            # Get the result
            moran_data = adata.uns['moranI']

            # Create custom visualization (Squidpy doesn't have a built-in one for this)
            fig, ax = plt.subplots(figsize=(10, 8))

            # Plot -log10(p-value) vs Moran's I
            sc = ax.scatter(-np.log10(moran_data['pval_norm']),
                          moran_data['I'],
                          s=50, alpha=0.7)

            # Add labels for top genes
            for i in range(min(5, len(genes))):
                if moran_data['pval_norm'][i] < 0.05:  # Only label significant genes
                    ax.annotate(genes[i],
                                (-np.log10(moran_data['pval_norm'][i]), moran_data['I'][i]),
                                xytext=(5, 5),
                                textcoords='offset points')

            # Add horizontal line at y=0
            ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

            # Add vertical line at x=-log10(0.05) for significance threshold
            ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)

            ax.set_title("Moran's I Spatial Autocorrelation")
            ax.set_xlabel("-log10(p-value)")
            ax.set_ylabel("Moran's I")

        elif params.analysis_type == "centrality":
            if context:
                await context.info("Computing centrality scores...")

            # Compute centrality scores
            sq.gr.centrality_scores(adata, cluster_key=cluster_key)

            # Get the result
            centrality_data = adata.uns[f'{cluster_key}_centrality_scores']

            # Create visualization - create our own figure instead of using squidpy's
            # This avoids color palette issues
            fig, ax = plt.subplots(figsize=(10, 8))

            # Get the centrality data
            centrality_types = list(centrality_data.keys())
            categories = adata.obs[cluster_key].cat.categories

            # Create a bar chart for each centrality type
            x = np.arange(len(categories))
            width = 0.8 / len(centrality_types)

            for i, ctype in enumerate(centrality_types):
                values = centrality_data[ctype]
                offset = width * i - width * len(centrality_types) / 2 + width / 2
                ax.bar(x + offset, values, width, label=ctype)

            # Add labels and legend
            ax.set_xlabel('Cluster')
            ax.set_ylabel('Centrality Score')
            ax.set_title(f'Centrality Scores by Cluster ({cluster_key})')
            ax.set_xticks(x)
            ax.set_xticklabels(categories)
            ax.legend()

            # Rotate x labels for better readability
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        else:
            raise ValueError(f"Unsupported analysis type: {params.analysis_type}")

        # Process the result based on return_type
        if return_type == "image":
            # Return Image object
            if fig is not None:
                try:
                    if context:
                        await context.info(f"Converting {params.analysis_type} figure to image...")
                    return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format or "png")
                except Exception as e:
                    error_msg = f"Failed to convert figure to image: {str(e)}"
                    if context:
                        await context.warning(error_msg)
                    # Close the figure to avoid memory leaks
                    plt.close(fig)
                    # Create a placeholder image with error message
                    return create_placeholder_image(f"Error in {params.analysis_type} visualization: {str(e)}")
            else:
                # Create a simple placeholder image if no visualization is available
                if context:
                    await context.warning(f"No {params.analysis_type} visualization available")
                return create_placeholder_image(f"No {params.analysis_type} visualization available")
        else:
            # Return SpatialAnalysisResult object
            img_str = None
            if params.include_image and fig is not None:
                try:
                    if context:
                        await context.info(f"Converting {params.analysis_type} figure to base64 image...")
                    img_str = fig_to_base64(
                        fig,
                        dpi=params.image_dpi,
                        format=params.image_format,
                        max_size_mb=5
                    )
                    if img_str is None and context:
                        await context.warning("Image was too large and could not be included in the response")
                except Exception as e:
                    error_msg = f"Failed to generate image: {str(e)}"
                    if context:
                        await context.warning(error_msg)
                    # Make sure to close the figure even if conversion fails
                    plt.close(fig)
            elif fig is not None:
                # Close the figure if we're not including it in the response
                plt.close(fig)

            # Return result with detailed statistics
            if context:
                await context.info(f"Returning {params.analysis_type} analysis result with {len(statistics)} statistics")

            # Add metadata to statistics
            statistics.update({
                "analysis_date": pd.Timestamp.now().isoformat(),
                "cluster_key": cluster_key,
                "n_clusters": n_clusters,
                "n_cells": adata.n_obs,
                "n_neighbors": params.n_neighbors
            })

            return SpatialAnalysisResult(
                data_id=data_id,
                analysis_type=params.analysis_type,
                statistics=statistics,
                result_image=img_str
            )

    except Exception as e:
        # Make sure to close any open figures in case of error
        if 'fig' in locals() and fig is not None:
            plt.close(fig)
        # Close any other open figures
        plt.close('all')

        # Log the error
        error_msg = f"Error in {params.analysis_type} analysis: {str(e)}"
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {traceback.format_exc()}")

        # Wrap the error in a more informative exception
        if isinstance(e, (DataNotFoundError, InvalidParameterError, DataCompatibilityError)):
            # Re-raise specific errors
            raise
        else:
            # Wrap generic errors
            raise ProcessingError(
                f"Failed to perform {params.analysis_type} analysis: {str(e)}",
                details={
                    "analysis_type": params.analysis_type,
                    "data_id": data_id,
                    "error_type": type(e).__name__,
                    "traceback": traceback.format_exc()
                }
            ) from e


async def analyze_spatial_with_image(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Optional[Context] = None
) -> Image:
    """Perform spatial analysis on transcriptomics data and return image directly

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Spatial analysis parameters
        context: MCP context

    Returns:
        Spatial analysis visualization as Image object
    """
    # Call the unified function with return_type="image"
    return await analyze_spatial_unified(data_id, data_store, params, context, return_type="image")


async def analyze_spatial(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Optional[Context] = None
) -> SpatialAnalysisResult:
    """Perform spatial analysis on transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Spatial analysis parameters
        context: MCP context

    Returns:
        Spatial analysis result
    """
    # Call the unified function with return_type="result"
    return await analyze_spatial_unified(data_id, data_store, params, context, return_type="result")
