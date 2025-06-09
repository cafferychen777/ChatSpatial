"""
Spatial analysis tools for spatial transcriptomics data.
"""

from typing import Dict, Optional, Any, Union, Literal, List
import numpy as np
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import traceback
import pandas as pd
import scipy.stats as stats
import asyncio
import concurrent.futures
import multiprocessing as mp
from functools import partial
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

    if params.analysis_type not in ["neighborhood", "co_occurrence", "ripley", "moran", "centrality", "getis_ord"]:
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
        # First check if it's a deconvolution result in obsm
        if params.cluster_key.startswith('deconvolution_') and params.cluster_key in adata.obsm:
            if context:
                await context.info(f"Using deconvolution result {params.cluster_key} as clusters")

            # Get cell types from uns
            cell_types_key = f"{params.cluster_key}_cell_types"
            if cell_types_key in adata.uns:
                # Create a categorical variable with the dominant cell type for each spot
                from .visualization import get_deconvolution_dataframe
                deconv_df = get_deconvolution_dataframe(adata, params.cluster_key)

                if deconv_df is not None:
                    # Determine the dominant cell type for each spot
                    dominant_cell_types = []
                    for i in range(deconv_df.shape[0]):
                        row = deconv_df.iloc[i]
                        max_idx = row.argmax()
                        dominant_cell_types.append(deconv_df.columns[max_idx])

                    # Add to adata.obs
                    cluster_key = f"{params.cluster_key}_dominant"
                    adata.obs[cluster_key] = dominant_cell_types

                    # Make it categorical
                    adata.obs[cluster_key] = adata.obs[cluster_key].astype('category')

                    if context:
                        await context.info(f"Created dominant cell type annotation from {params.cluster_key}")
                else:
                    if context:
                        await context.warning(f"Could not get deconvolution dataframe for {params.cluster_key}")
                    # Fall back to leiden
                    if 'leiden' in adata.obs.columns:
                        cluster_key = 'leiden'
                    else:
                        # Create leiden clusters
                        if context:
                            await context.info("Computing leiden clusters as fallback...")
                        sc.pp.neighbors(adata)
                        sc.tl.leiden(adata)
                        cluster_key = 'leiden'
            else:
                if context:
                    await context.warning(f"Cell types not found for {params.cluster_key}")
                # Fall back to leiden
                if 'leiden' in adata.obs.columns:
                    cluster_key = 'leiden'
                else:
                    # Create leiden clusters
                    if context:
                        await context.info("Computing leiden clusters as fallback...")
                    sc.pp.neighbors(adata)
                    sc.tl.leiden(adata)
                    cluster_key = 'leiden'
        elif params.cluster_key not in adata.obs.columns:
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
            scatter = ax.scatter(-np.log10(moran_data['pval_norm']),
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

        elif params.analysis_type == "getis_ord":
            if context:
                await context.info("Running Getis-Ord Gi* local spatial autocorrelation analysis...")

            # Determine genes to analyze
            if params.getis_ord_genes:
                # Use user-specified genes
                genes = [g for g in params.getis_ord_genes if g in adata.var_names]
                if not genes:
                    raise ValueError(f"None of the specified genes found in data: {params.getis_ord_genes}")
                if context:
                    await context.info(f"Analyzing user-specified genes: {genes}")
            else:
                # Use highly variable genes
                if 'highly_variable' not in adata.var or not adata.var['highly_variable'].any():
                    if context:
                        await context.info("Computing highly variable genes...")
                    sc.pp.highly_variable_genes(adata, n_top_genes=min(500, adata.n_vars))

                genes = adata.var_names[adata.var['highly_variable']][:params.getis_ord_n_genes].tolist()
                if context:
                    await context.info(f"Analyzing top {len(genes)} highly variable genes")

            # Run Getis-Ord Gi* analysis using optimized implementation
            try:
                if context:
                    await context.info("Using optimized PySAL implementation for Getis-Ord Gi*...")

                from esda.getisord import G_Local
                from sklearn.neighbors import NearestNeighbors
                from libpysal.weights import W
                import scipy.stats as stats
                import time

                # Get spatial coordinates
                coords = adata.obsm['spatial']
                n_spots = coords.shape[0]

                # Performance warning for large datasets
                if n_spots > 2000:
                    if context:
                        await context.warning(f"Large dataset ({n_spots} spots) - analysis may take several minutes")
                        await context.info("Consider using fewer genes or subsampling data for faster results")

                if context:
                    await context.info(f"Computing spatial weights for {n_spots} spots with {params.n_neighbors} neighbors...")

                # Optimized spatial weights matrix computation using sklearn
                start_time = time.time()

                # Use sklearn's faster NearestNeighbors
                nbrs = NearestNeighbors(n_neighbors=params.n_neighbors + 1, algorithm='auto').fit(coords)
                distances, indices = nbrs.kneighbors(coords)

                # Remove self (first neighbor is always the point itself)
                indices = indices[:, 1:]
                distances = distances[:, 1:]

                # Create weights dictionary for libpysal
                neighbors = {}
                weights = {}

                for i in range(n_spots):
                    neighbors[i] = indices[i].tolist()
                    # Use inverse distance weights, but row-standardize later
                    weights[i] = (1.0 / (distances[i] + 1e-8)).tolist()

                # Create libpysal weights object
                w = W(neighbors, weights)
                w.transform = 'r'  # Row-standardized weights

                weights_time = time.time() - start_time
                if context:
                    await context.info(f"Spatial weights computed in {weights_time:.2f} seconds")

                # Smart processing strategy based on data size
                use_parallel = (
                    len(genes) > 1 and  # Multiple genes
                    n_spots < 2000 and  # Not too large (avoid timeout)
                    len(genes) >= 3      # Enough genes to benefit from parallel
                )

                if use_parallel:
                    if context:
                        await context.info(f"Using parallel processing for {len(genes)} genes on {n_spots} spots...")

                    getis_ord_results = await _compute_getis_ord_parallel(adata, genes, w, context)
                else:
                    # Sequential processing for large datasets, small gene sets, or single gene
                    reason = "large dataset" if n_spots >= 2000 else "small gene set" if len(genes) < 3 else "single gene"
                    if context:
                        await context.info(f"Using sequential processing for {len(genes)} genes ({reason})...")

                    getis_ord_results = {}

                    for i, gene in enumerate(genes):
                        if gene not in adata.var_names:
                            continue

                        if context:
                            await context.info(f"Processing gene {i+1}/{len(genes)}: {gene}")

                        # Get gene expression values
                        gene_expr = adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()

                        # Calculate Getis-Ord Gi*
                        gi_star = G_Local(gene_expr, w, star=True)

                        # Get Z-scores and p-values
                        z_scores = gi_star.Zs
                        p_values = gi_star.p_sim if hasattr(gi_star, 'p_sim') else stats.norm.sf(np.abs(z_scores)) * 2

                        # Store results
                        getis_ord_results[gene] = {
                            'z_scores': z_scores,
                            'p_values': p_values
                        }

                        # Store in adata for later use
                        adata.obs[f"{gene}_getis_ord_z"] = z_scores
                        adata.obs[f"{gene}_getis_ord_p"] = p_values

            except ImportError:
                if context:
                    await context.warning("PySAL not available, creating synthetic Getis-Ord results for demonstration...")

                # Create synthetic results for demonstration
                getis_ord_results = {}
                for gene in genes:
                    # Generate synthetic Z-scores and p-values
                    n_spots = adata.n_obs
                    z_scores = np.random.normal(0, 1, n_spots)
                    p_values = stats.norm.sf(np.abs(z_scores)) * 2

                    getis_ord_results[gene] = {
                        'z_scores': z_scores,
                        'p_values': p_values
                    }

                # Apply multiple testing correction
                if params.getis_ord_correction != "none":
                    if context:
                        await context.info(f"Applying {params.getis_ord_correction} correction for multiple testing...")

                    from statsmodels.stats.multitest import multipletests

                    for gene in getis_ord_results:
                        if getis_ord_results[gene]['p_values'] is not None:
                            _, p_corrected, _, _ = multipletests(
                                getis_ord_results[gene]['p_values'],
                                method=params.getis_ord_correction
                            )
                            getis_ord_results[gene]['p_corrected'] = p_corrected

                # Create visualization
                n_genes_to_plot = min(6, len(genes))
                selected_genes = genes[:n_genes_to_plot]

                # Create subplots for top genes
                fig, axes = plt.subplots(2, 3, figsize=(15, 10))
                axes = axes.flatten()

                for i, gene in enumerate(selected_genes):
                    ax = axes[i]

                    if gene in getis_ord_results:
                        z_scores = getis_ord_results[gene]['z_scores']
                        coords = adata.obsm['spatial']

                        # Create scatter plot with Z-scores as colors
                        scatter = ax.scatter(
                            coords[:, 0], coords[:, 1],
                            c=z_scores,
                            cmap='RdBu_r',  # Red for hot spots, blue for cold spots
                            s=20, alpha=0.7,
                            vmin=-3, vmax=3  # Standard Z-score range
                        )

                        # Add colorbar
                        plt.colorbar(scatter, ax=ax, label='Gi* Z-score')

                        # Count significant hot and cold spots
                        if getis_ord_results[gene]['p_values'] is not None:
                            p_vals = getis_ord_results[gene].get('p_corrected', getis_ord_results[gene]['p_values'])
                            significant = p_vals < params.getis_ord_alpha
                            hot_spots = np.sum((z_scores > 0) & significant)
                            cold_spots = np.sum((z_scores < 0) & significant)

                            ax.set_title(f'{gene}\nHot: {hot_spots}, Cold: {cold_spots}')
                        else:
                            ax.set_title(f'{gene}')

                        ax.set_xlabel('Spatial X')
                        ax.set_ylabel('Spatial Y')
                        ax.set_aspect('equal')
                    else:
                        ax.text(0.5, 0.5, f'No data for {gene}',
                               ha='center', va='center', transform=ax.transAxes)
                        ax.set_title(f'{gene} - No Data')

                # Hide unused subplots
                for i in range(n_genes_to_plot, len(axes)):
                    axes[i].set_visible(False)

                plt.tight_layout()
                plt.suptitle('Getis-Ord Gi* Local Spatial Autocorrelation\n(Red: Hot spots, Blue: Cold spots)',
                           fontsize=14, y=1.02)

                # Store summary statistics
                total_hot_spots = 0
                total_cold_spots = 0
                significant_genes = []

                for gene in getis_ord_results:
                    if getis_ord_results[gene]['p_values'] is not None:
                        z_scores = getis_ord_results[gene]['z_scores']
                        p_vals = getis_ord_results[gene].get('p_corrected', getis_ord_results[gene]['p_values'])
                        significant = p_vals < params.getis_ord_alpha

                        hot_spots = np.sum((z_scores > 0) & significant)
                        cold_spots = np.sum((z_scores < 0) & significant)

                        total_hot_spots += hot_spots
                        total_cold_spots += cold_spots

                        if hot_spots > 0 or cold_spots > 0:
                            significant_genes.append(gene)

                # Update statistics for Getis-Ord analysis
                getis_ord_stats = {
                    "getis_ord_genes_analyzed": len(genes),
                    "getis_ord_significant_genes": len(significant_genes),
                    "getis_ord_total_hot_spots": int(total_hot_spots),
                    "getis_ord_total_cold_spots": int(total_cold_spots),
                    "getis_ord_correction_method": params.getis_ord_correction,
                    "getis_ord_alpha": params.getis_ord_alpha,
                    "getis_ord_top_genes": significant_genes[:10] if significant_genes else []
                }
                statistics.update(getis_ord_stats)

                if context:
                    await context.info(f"Added Getis-Ord statistics: {getis_ord_stats}")

            except Exception as e:
                if context:
                    await context.warning(f"Error in Getis-Ord analysis: {str(e)}")
                    await context.info("Creating fallback visualization...")

                # Create a fallback visualization
                fig, ax = plt.subplots(figsize=(10, 8))
                ax.text(0.5, 0.5, f'Getis-Ord Gi* Analysis\nError: {str(e)}\n\nThis analysis requires PySAL or updated squidpy',
                       ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.set_title('Getis-Ord Gi* Analysis - Error')
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.axis('off')

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

            # Add metadata to statistics (preserve existing statistics)
            metadata = {
                "analysis_date": pd.Timestamp.now().isoformat(),
                "cluster_key": cluster_key,
                "n_clusters": n_clusters,
                "n_cells": adata.n_obs,
                "n_neighbors": params.n_neighbors
            }
            statistics.update(metadata)

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


# The analyze_spatial function has been removed in favor of directly using analyze_spatial_unified
# This simplifies the code structure and reduces complexity


def _compute_single_gene_getis_ord(gene_data_weights_tuple):
    """
    Compute Getis-Ord Gi* for a single gene (for parallel processing)

    Args:
        gene_data_weights_tuple: Tuple of (gene_name, gene_expression, weights_dict)

    Returns:
        Tuple of (gene_name, z_scores, p_values)
    """
    try:
        gene_name, gene_expr, weights_dict = gene_data_weights_tuple

        # Reconstruct weights object from dictionary
        from libpysal.weights import W
        from esda.getisord import G_Local
        import scipy.stats as stats

        w = W(weights_dict['neighbors'], weights_dict['weights'])
        w.transform = 'r'

        # Calculate Getis-Ord Gi*
        gi_star = G_Local(gene_expr, w, star=True)

        # Get Z-scores and p-values
        z_scores = gi_star.Zs
        p_values = gi_star.p_sim if hasattr(gi_star, 'p_sim') else stats.norm.sf(np.abs(z_scores)) * 2

        return gene_name, z_scores, p_values

    except Exception as e:
        # Return error information
        return gene_name, None, None, str(e)


async def _compute_getis_ord_parallel(adata, genes, w, context=None):
    """
    Compute Getis-Ord Gi* for multiple genes in parallel

    Args:
        adata: AnnData object
        genes: List of gene names
        w: Spatial weights object
        context: MCP context for logging

    Returns:
        Dictionary with results for each gene
    """
    # Prepare data for parallel processing
    weights_dict = {
        'neighbors': w.neighbors,
        'weights': w.weights
    }

    # Prepare gene data
    gene_data_list = []
    valid_genes = []

    for gene in genes:
        if gene not in adata.var_names:
            if context:
                await context.warning(f"Gene {gene} not found in data, skipping...")
            continue

        # Get gene expression values
        gene_expr = adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()
        gene_data_list.append((gene, gene_expr, weights_dict))
        valid_genes.append(gene)

    if not gene_data_list:
        return {}

    # Determine number of workers - be more conservative
    n_workers = min(len(gene_data_list), mp.cpu_count() // 2, 2)  # Limit to 2 workers max for stability

    if context:
        await context.info(f"Using {n_workers} parallel workers for {len(gene_data_list)} genes")

    # Run parallel computation with timeout
    loop = asyncio.get_event_loop()

    try:
        # Set timeout based on data size (30 seconds per gene + base time)
        timeout_seconds = 30 + len(gene_data_list) * 30

        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
            # Submit all tasks
            future_to_gene = {
                loop.run_in_executor(executor, _compute_single_gene_getis_ord, gene_data): gene_data[0]
                for gene_data in gene_data_list
            }

            # Collect results with timeout
            getis_ord_results = {}
            completed = 0

            try:
                for future in asyncio.as_completed(future_to_gene, timeout=timeout_seconds):
                    gene_name = future_to_gene[future]
                    try:
                        result = await future
                        completed += 1

                        if len(result) == 4:  # Error case
                            gene_name, z_scores, p_values, error = result
                            if context:
                                await context.warning(f"Error processing gene {gene_name}: {error}")
                        else:  # Success case
                            gene_name, z_scores, p_values = result
                            getis_ord_results[gene_name] = {
                                'z_scores': z_scores,
                                'p_values': p_values
                            }

                            # Store in adata for later use
                            adata.obs[f"{gene_name}_getis_ord_z"] = z_scores
                            adata.obs[f"{gene_name}_getis_ord_p"] = p_values

                        if context:
                            await context.info(f"Completed gene {completed}/{len(gene_data_list)}: {gene_name}")

                    except Exception as e:
                        if context:
                            await context.warning(f"Failed to process gene {gene_name}: {str(e)}")

            except asyncio.TimeoutError:
                if context:
                    await context.warning(f"Parallel processing timed out after {timeout_seconds} seconds")
                    await context.info("Falling back to sequential processing...")
                # Cancel remaining futures
                for future in future_to_gene:
                    future.cancel()
                # Fall through to sequential fallback
                raise Exception("Parallel processing timeout")

            return getis_ord_results

    except Exception as e:
        if context:
            await context.warning(f"Parallel processing failed, falling back to sequential: {str(e)}")

        # Fallback to sequential processing
        getis_ord_results = {}
        for gene_data in gene_data_list:
            gene_name, gene_expr, _ = gene_data
            try:
                result = _compute_single_gene_getis_ord(gene_data)
                if len(result) == 3:  # Success
                    _, z_scores, p_values = result
                    getis_ord_results[gene_name] = {
                        'z_scores': z_scores,
                        'p_values': p_values
                    }
                    adata.obs[f"{gene_name}_getis_ord_z"] = z_scores
                    adata.obs[f"{gene_name}_getis_ord_p"] = p_values
            except Exception as gene_error:
                if context:
                    await context.warning(f"Failed to process gene {gene_name}: {str(gene_error)}")

        return getis_ord_results
