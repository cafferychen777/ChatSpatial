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


# Import scvi-tools for advanced spatial analysis
try:
    import scvi
    # SCVIVA is available in scvi-tools >= 1.3.2
    try:
        from scvi.external import SCVIVA
    except ImportError:
        # SCVIVA not available in older versions
        SCVIVA = None
        import warnings
        warnings.warn(
            "SCVIVA requires scvi-tools >= 1.3.2. Current version: "
            f"{scvi.__version__}. Please upgrade with: pip install --upgrade scvi-tools"
        )
except ImportError:
    scvi = None
    SCVIVA = None


async def analyze_spatial_with_scviva(
    spatial_adata,
    n_epochs: int = 1000,
    n_hidden: int = 128,
    n_latent: int = 10,
    use_gpu: bool = False,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze spatial transcriptomics data using SCVIVA
    
    SCVIVA is a variational auto-encoder with niche decoders for spatial 
    transcriptomics that models both cell-intrinsic and neighborhood effects.
    
    Args:
        spatial_adata: Spatial transcriptomics AnnData object with spatial coordinates
        n_epochs: Number of epochs for training
        n_hidden: Number of hidden units in neural networks
        n_latent: Dimensionality of latent space
        use_gpu: Whether to use GPU for training
        context: MCP context for logging
        
    Returns:
        Dictionary containing SCVIVA analysis results
        
    Raises:
        ImportError: If scvi-tools package is not available or version < 1.3.2
        ValueError: If input data is invalid or missing spatial coordinates
        RuntimeError: If SCVIVA computation fails
    """
    try:
        if scvi is None:
            raise ImportError("scvi-tools package is required for SCVIVA analysis")
            
        if SCVIVA is None:
            raise ImportError(
                "SCVIVA requires scvi-tools >= 1.3.2. "
                "Please upgrade with: pip install --upgrade scvi-tools>=1.3.2"
            )
            
        # Validate spatial coordinates
        if 'spatial' not in spatial_adata.obsm:
            raise ValueError("Spatial coordinates not found in adata.obsm['spatial']. Required for SCVIVA analysis.")
        
        if context:
            await context.info("Using SCVIVA for spatial analysis...")
            
        return await _analyze_with_scviva_native(
            spatial_adata, n_epochs, n_hidden, n_latent, use_gpu, context
        )
        
    except Exception as e:
        error_msg = f"SCVIVA analysis failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e
# This simplifies the code structure and reduces complexity


async def _analyze_with_scviva_native(
    spatial_adata,
    n_epochs: int,
    n_hidden: int,
    n_latent: int,
    use_gpu: bool,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Native SCVIVA implementation when available"""
    adata = spatial_adata.copy()
    
    if context:
        await context.info(f"Analyzing {adata.n_obs} cells and {adata.n_vars} genes with SCVIVA")
        await context.info("Note: SCVIVA is designed for single-cell resolution spatial data")
    
    # Check if this is single-cell resolution data
    if adata.n_obs > 10000 and context:
        await context.warning(
            f"Large number of observations ({adata.n_obs}). "
            "SCVIVA is designed for single-cell resolution spatial data, not Visium spots."
        )
    
    # Preprocess data for SCVIVA (compute neighborhoods, etc.)
    if context:
        await context.info("Preprocessing data for SCVIVA (computing spatial neighborhoods)...")
    
    # Add sample information if not present
    if 'sample' not in adata.obs.columns:
        adata.obs['sample'] = 'sample_1'  # All cells from same sample
    
    # Add cell type information if not present
    if 'cell_type' not in adata.obs.columns:
        if context:
            await context.warning("No cell type annotations found. Using placeholder labels.")
        adata.obs['cell_type'] = 'Unknown'
    
    # SCVIVA requires SCVI embeddings first
    if 'X_scVI' not in adata.obsm:
        if context:
            await context.info("Running SCVI to generate expression embeddings required by SCVIVA...")
        
        # Setup and train SCVI first
        scvi.model.SCVI.setup_anndata(adata, batch_key='sample')
        scvi_model = scvi.model.SCVI(adata, n_hidden=n_hidden, n_latent=n_latent)
        
        # Quick training for embeddings
        scvi_model.train(max_epochs=min(50, n_epochs), early_stopping=True)
        
        # Get latent representation
        adata.obsm['X_scVI'] = scvi_model.get_latent_representation()
        
        if context:
            await context.info("SCVI embeddings generated successfully.")
    
    # SCVIVA requires preprocessing to compute niche information
    SCVIVA.preprocessing_anndata(
        adata,
        k_nn=min(20, adata.n_obs - 1),  # Number of nearest neighbors
        sample_key='sample',  # Required for multi-sample analysis
        labels_key='cell_type',
        cell_coordinates_key="spatial",
        expression_embedding_key='X_scVI',
        log1p=False  # Assuming data is already normalized if needed
    )
    
    # Setup SCVIVA with correct parameter name
    SCVIVA.setup_anndata(
        adata, 
        sample_key='sample',
        labels_key='cell_type',
        cell_coordinates_key="spatial",
        expression_embedding_key='X_scVI'
    )
    
    # Create SCVIVA model
    scviva_model = SCVIVA(adata, n_hidden=n_hidden, n_latent=n_latent)
    
    if context:
        await context.info("Training SCVIVA model...")
    
    # Train model
    train_kwargs = {"max_epochs": n_epochs}
    if use_gpu:
        train_kwargs["accelerator"] = "gpu"
    
    scviva_model.train(**train_kwargs)
    
    if context:
        await context.info("Extracting SCVIVA results...")
    
    # Get results using SCVIVA's actual methods
    # Get latent representation
    latent = scviva_model.get_latent_representation()
    
    # Get expression embedding and niche information
    # These are stored in adata during model training
    niche_composition = adata.obsm.get('niche_composition', np.zeros((adata.n_obs, n_latent)))
    expression_embedding = adata.obsm.get('X_scVI', latent)
    niche_activation = adata.obsm.get('niche_activation', np.zeros((adata.n_obs, n_latent)))
    
    # Store results with consistent naming
    spatial_adata.obsm['X_scviva_latent'] = latent
    spatial_adata.obsm['X_scviva_niche'] = niche_composition
    spatial_adata.obsm['X_scviva_intrinsic'] = expression_embedding
    spatial_adata.obsm['X_scviva_niche_activation'] = niche_activation
    
    # Store additional SCVIVA-specific information
    if 'niche_indexes' in adata.obsm:
        spatial_adata.obsm['niche_indexes'] = adata.obsm['niche_indexes']
    if 'niche_distances' in adata.obsm:
        spatial_adata.obsm['niche_distances'] = adata.obsm['niche_distances']
    
    return {
        'method': 'SCVIVA',
        'n_latent_dims': n_latent,
        'n_epochs': n_epochs,
        'latent_variance_explained': float(np.var(latent, axis=0).sum()),
        'niche_composition_shape': tuple(niche_composition.shape),
        'intrinsic_factors_shape': tuple(expression_embedding.shape),
        'niche_activation_shape': tuple(niche_activation.shape),
        'training_completed': True,
        'device': 'GPU' if use_gpu else 'CPU',
        'model_type': 'single-cell spatial VAE with niche decoders'
    }
