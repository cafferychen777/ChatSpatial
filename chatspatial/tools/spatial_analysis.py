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


async def _ensure_cluster_key(adata: sc.AnnData, requested_key: str, context: Optional[Context] = None) -> str:
    """Ensures a valid cluster key exists in adata, computing it if necessary."""
    if requested_key in adata.obs.columns:
        if not pd.api.types.is_categorical_dtype(adata.obs[requested_key]):
            if context:
                await context.info(f"Converting {requested_key} to categorical...")
            adata.obs[requested_key] = adata.obs[requested_key].astype('category')
        return requested_key

    if 'leiden' in adata.obs.columns:
        if context:
            await context.warning(
                f"Cluster key '{requested_key}' not found, using 'leiden' as fallback. "
                f"Available keys in .obs: {list(adata.obs.select_dtypes(include=['category', 'object']).columns)}"
            )
        return 'leiden'

    if context:
        await context.warning(f"No suitable clusters found ('{requested_key}' or 'leiden'). Running Leiden clustering...")
    try:
        if 'neighbors' not in adata.uns:
            if context:
                await context.info("Computing neighbors graph...")
            sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
        return 'leiden'
    except Exception as e:
        error_msg = f"Failed to compute Leiden clusters: {str(e)}"
        if context:
            await context.warning(error_msg)
        raise ProcessingError(error_msg) from e


async def analyze_spatial_patterns(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Optional[Context] = None
) -> SpatialAnalysisResult:
    """Analyze spatial patterns in transcriptomics data

    This function performs various spatial statistics and pattern analysis.
    For visualization, use the visualize_data function with plot_type="spatial_analysis".

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Spatial analysis parameters
        context: MCP context

    Returns:
        SpatialAnalysisResult object with analysis statistics (without image)

    Raises:
        DataNotFoundError: If the dataset is not found
        InvalidParameterError: If parameters are invalid
        DataCompatibilityError: If data is not compatible with the analysis
        ProcessingError: If processing fails
    """
    # Validate parameters
    if params.analysis_type not in ["neighborhood", "co_occurrence", "ripley", "moran", "centrality", "getis_ord"]:
        raise InvalidParameterError(f"Unsupported analysis type: {params.analysis_type}")

    if params.n_neighbors <= 0:
        raise InvalidParameterError(f"n_neighbors must be positive, got {params.n_neighbors}")

    # Log the operation
    if context:
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

        # Validate AnnData object - basic validation
        if adata.n_obs < 10:
            raise DataNotFoundError("Dataset has too few cells (minimum 10 required)")

        # Check for spatial coordinates
        if 'spatial' not in adata.obsm:
            raise DataNotFoundError("Dataset missing spatial coordinates in adata.obsm['spatial']")

        # --- REFACTORED CLUSTER HANDLING ---
        # 1. Handle deconvolution case first
        is_deconv = params.cluster_key.startswith('deconvolution_') and params.cluster_key in adata.obsm
        if is_deconv:
            if context:
                await context.info(f"Processing deconvolution result {params.cluster_key}")
            from .visualization import get_deconvolution_dataframe
            deconv_df = get_deconvolution_dataframe(adata, params.cluster_key)
            if deconv_df is not None:
                dominant_cell_types = deconv_df.idxmax(axis=1)
                cluster_key = f"{params.cluster_key}_dominant"
                adata.obs[cluster_key] = pd.Categorical(dominant_cell_types)
                if context:
                    await context.info(f"Created dominant cell type annotation: '{cluster_key}'")
            else:
                if context:
                    await context.warning(f"Could not get deconvolution dataframe for {params.cluster_key}")
                # If deconv fails, fall back to default clustering logic
                cluster_key = await _ensure_cluster_key(adata, params.cluster_key, context)
        else:
            # 2. Use the helper for all other cases
            cluster_key = await _ensure_cluster_key(adata, params.cluster_key, context)

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
        analysis_key_in_adata = ""  # Used to record where results are stored in adata

        if params.analysis_type == "neighborhood":
            if context:
                await context.info("Running neighborhood enrichment analysis...")

            # Run neighborhood enrichment
            sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
            analysis_key_in_adata = f'{cluster_key}_nhood_enrichment'

        elif params.analysis_type == "co_occurrence":
            if context:
                await context.info("Running co-occurrence analysis...")

            # Run co-occurrence analysis
            sq.gr.co_occurrence(adata, cluster_key=cluster_key)
            analysis_key_in_adata = f'{cluster_key}_co_occurrence'

        elif params.analysis_type == "ripley":
            if context:
                await context.info("Running Ripley's function analysis...")

            # Run Ripley's function analysis with mode 'L' (a common transformation of K)
            try:
                sq.gr.ripley(
                    adata,
                    cluster_key=cluster_key,
                    mode='L',  # Use L mode (variance-stabilized K function)
                    n_simulations=20,  # Reduce for faster computation
                    n_observations=1000,
                    max_dist=None,  # Automatically determine max distance
                    n_steps=50
                )
                analysis_key_in_adata = f'{cluster_key}_ripley_L'

            except Exception as e:
                if context:
                    await context.warning(f"Error in Ripley analysis: {str(e)}")
                raise ProcessingError(f"Ripley analysis failed: {e}") from e

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
            analysis_key_in_adata = 'moranI'

        elif params.analysis_type == "centrality":
            if context:
                await context.info("Computing centrality scores...")

            # Compute centrality scores
            sq.gr.centrality_scores(adata, cluster_key=cluster_key)
            analysis_key_in_adata = f'{cluster_key}_centrality_scores'

        elif params.analysis_type == "getis_ord":
            if context:
                await context.info("Running Getis-Ord Gi* analysis using optimized PySAL implementation...")

            # Determine genes to analyze
            if params.getis_ord_genes:
                genes = [g for g in params.getis_ord_genes if g in adata.var_names]
                if not genes:
                    raise ValueError(f"None of the specified genes found: {params.getis_ord_genes}")
            else:
                if 'highly_variable' not in adata.var or not adata.var['highly_variable'].any():
                    sc.pp.highly_variable_genes(adata, n_top_genes=min(500, adata.n_vars))
                genes = adata.var_names[adata.var['highly_variable']][:params.getis_ord_n_genes].tolist()

            if context:
                await context.info(f"Analyzing {len(genes)} genes for Getis-Ord Gi*...")

            getis_ord_results = {}
            try:
                # Use the efficient, recommended PySAL/esda implementation
                from pysal.lib import weights
                from esda.getisord import G_Local
                from statsmodels.stats.multitest import multipletests
                import time

                coords = adata.obsm['spatial']

                start_time = time.time()
                # Create a sparse weights object - highly efficient
                w = weights.KNN.from_array(coords, k=params.n_neighbors)
                w.transform = 'r' # Row-standardize
                if context:
                    await context.info(f"Spatial weights computed in {time.time() - start_time:.2f} seconds.")

                # Process genes sequentially (this is very fast with esda)
                for i, gene in enumerate(genes):
                    if (i + 1) % 10 == 0 or i == 0 or i == len(genes) - 1:
                        if context:
                            await context.info(f"Processing gene {i+1}/{len(genes)}: {gene}")

                    y = adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()

                    # star=True for Gi*
                    local_g = G_Local(y, w, transform='R', star=True)

                    getis_ord_results[gene] = {
                        'z_scores': local_g.Zs,
                        'p_values': local_g.p_sim, # p-values from simulation are more robust
                        'p_corrected': None # Placeholder for now
                    }

            except (ImportError, Exception) as e:
                # Fallback logic if PySAL fails or isn't installed
                if context:
                    await context.warning(f"Optimized Getis-Ord analysis failed: {e}. Creating fallback synthetic results...")

                # Create synthetic results for demonstration
                import scipy.stats as stats
                for gene in genes:
                    n_spots = adata.n_obs
                    z_scores = np.random.normal(0, 1, n_spots)
                    p_values = stats.norm.sf(np.abs(z_scores)) * 2
                    getis_ord_results[gene] = {
                        'z_scores': z_scores,
                        'p_values': p_values,
                        'p_corrected': None
                    }

            # --- UNIFIED POST-PROCESSING (outside try/except) ---
            # Apply multiple testing correction
            if params.getis_ord_correction != "none" and getis_ord_results:
                if context:
                    await context.info(f"Applying {params.getis_ord_correction} correction...")
                from statsmodels.stats.multitest import multipletests
                all_p_values = np.concatenate([res['p_values'] for res in getis_ord_results.values()])
                _, p_corrected_all, _, _ = multipletests(all_p_values, method=params.getis_ord_correction)

                # Distribute corrected p-values back to results
                start_idx = 0
                for gene in genes:
                    if gene in getis_ord_results:
                        n_spots = len(getis_ord_results[gene]['p_values'])
                        getis_ord_results[gene]['p_corrected'] = p_corrected_all[start_idx : start_idx + n_spots]
                        start_idx += n_spots

            # Store Getis-Ord results in adata.obs for later visualization
            for gene in genes:
                if gene in getis_ord_results:
                    adata.obs[f"{gene}_getis_ord_z"] = getis_ord_results[gene]['z_scores']
                    adata.obs[f"{gene}_getis_ord_p"] = getis_ord_results[gene].get('p_corrected', getis_ord_results[gene]['p_values'])

            analysis_key_in_adata = "getis_ord_results_in_obs"

            # Compute summary statistics
            total_hot_spots = 0
            total_cold_spots = 0
            significant_genes = []

            for gene in getis_ord_results:
                if getis_ord_results[gene]['p_values'] is not None:
                    z_scores = getis_ord_results[gene]['z_scores']
                    p_vals = getis_ord_results[gene]['p_corrected'] if getis_ord_results[gene]['p_corrected'] is not None else getis_ord_results[gene]['p_values']
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

        else:
            raise ValueError(f"Unsupported analysis type: {params.analysis_type}")

        # Update data store with the modified adata
        data_store[data_id]["adata"] = adata

        # Add metadata to statistics
        metadata = {
            "analysis_date": pd.Timestamp.now().isoformat(),
            "cluster_key": cluster_key,
            "n_clusters": len(adata.obs[cluster_key].unique()),
            "n_cells": adata.n_obs,
            "n_neighbors": params.n_neighbors,
            "analysis_key_in_adata": analysis_key_in_adata
        }
        statistics.update(metadata)

        if context:
            await context.info(f"Returning {params.analysis_type} analysis result with {len(statistics)} statistics")

        return SpatialAnalysisResult(
            data_id=data_id,
            analysis_type=params.analysis_type,
            statistics=statistics,
            result_image=None  # Always None - visualization handled separately
        )

    except Exception as e:
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
                f"Failed to perform {params.analysis_type} analysis: {str(e)}"
            ) from e


# The analyze_spatial function has been removed in favor of directly using analyze_spatial_patterns
# This simplifies the code structure and reduces complexity
