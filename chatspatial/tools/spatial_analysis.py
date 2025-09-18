"""
Unified spatial analysis module for spatial transcriptomics data.

This module provides comprehensive spatial analysis capabilities including:
- Spatial autocorrelation (Moran's I, Geary's C)
- Local spatial statistics (LISA, Getis-Ord)
- Spatial pattern analysis (Ripley's K, neighborhood enrichment)
- Advanced analysis (Bivariate Moran's I, Join Count, Network properties)
- Deep learning analysis (SCVIVA)

All functions are accessible through the main entry point: analyze_spatial_patterns()
"""

from typing import Dict, Optional, Any, Union, List, Tuple
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
import matplotlib.pyplot as plt
import traceback
import scipy.stats as stats
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
import logging
import datetime

from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from ..models.data import SpatialAnalysisParameters
from ..models.analysis import SpatialAnalysisResult

# Import standardized utilities
from ..utils.image_utils import fig_to_image, fig_to_base64, create_placeholder_image
from ..utils.error_handling import (
    SpatialMCPError, DataNotFoundError, InvalidParameterError,
    ProcessingError, DataCompatibilityError, validate_adata,
    handle_error, try_except_with_feedback
)

logger = logging.getLogger(__name__)


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

async def analyze_spatial_patterns(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Optional[Context] = None
) -> SpatialAnalysisResult:
    """
    Unified entry point for all spatial analysis methods.
    
    This function performs various spatial statistics and pattern analysis.
    For visualization, use the visualize_data function with plot_type="spatial_analysis".
    
    Parameters
    ----------
    data_id : str
        Dataset identifier
    data_store : Dict[str, Any]
        Dictionary storing loaded datasets
    params : SpatialAnalysisParameters
        Analysis parameters including analysis_type
    context : Optional[Context]
        MCP context for logging
    
    Returns
    -------
    SpatialAnalysisResult
        Analysis result with statistics and metadata
        
    Raises
    ------
    DataNotFoundError
        If dataset not found
    InvalidParameterError
        If parameters are invalid
    ProcessingError
        If analysis fails
    """
    # Validate parameters
    supported_types = [
        "neighborhood", "co_occurrence", "ripley", "moran", "geary",
        "centrality", "getis_ord", "bivariate_moran", "join_count",
        "network_properties", "spatial_centrality", "scviva"
    ]
    
    if params.analysis_type not in supported_types:
        raise InvalidParameterError(f"Unsupported analysis type: {params.analysis_type}")
    
    if params.n_neighbors <= 0:
        raise InvalidParameterError(f"n_neighbors must be positive, got {params.n_neighbors}")
    
    # Log operation
    if context:
        await context.info(f"Performing {params.analysis_type} spatial analysis")
    
    # Retrieve dataset
    if data_id not in data_store:
        raise DataNotFoundError(f"Dataset {data_id} not found in data store")
    
    try:
        # Get AnnData object (make a copy to avoid modifying original)
        adata = data_store[data_id]["adata"].copy()
        
        # Basic validation
        _validate_spatial_data(adata)
        
        # Ensure cluster key if needed
        cluster_key = await _ensure_cluster_key(adata, params.cluster_key, context)
        
        # Ensure spatial neighbors
        await _ensure_spatial_neighbors(adata, params.n_neighbors, context)
        
        # Route to appropriate analysis function
        if params.analysis_type == "moran":
            result = await _analyze_morans_i(adata, params, context)
        elif params.analysis_type == "geary":
            result = await _analyze_gearys_c(adata, params, context)
        elif params.analysis_type == "neighborhood":
            result = await _analyze_neighborhood_enrichment(adata, cluster_key, params, context)
        elif params.analysis_type == "co_occurrence":
            result = await _analyze_co_occurrence(adata, cluster_key, params, context)
        elif params.analysis_type == "ripley":
            result = await _analyze_ripleys_k(adata, cluster_key, params, context)
        elif params.analysis_type == "getis_ord":
            result = await _analyze_getis_ord(adata, params, context)
        elif params.analysis_type == "centrality":
            result = await _analyze_centrality(adata, cluster_key, params, context)
        elif params.analysis_type == "bivariate_moran":
            result = await _analyze_bivariate_moran(adata, params, context)
        elif params.analysis_type == "join_count":
            result = await _analyze_join_count(adata, cluster_key, params, context)
        elif params.analysis_type == "network_properties":
            result = await _analyze_network_properties(adata, cluster_key, params, context)
        elif params.analysis_type == "spatial_centrality":
            result = await _analyze_spatial_centrality(adata, cluster_key, params, context)
        elif params.analysis_type == "scviva":
            result = await _analyze_with_scviva(adata, params, context)
        else:
            raise ValueError(f"Analysis type {params.analysis_type} not implemented")
        
        # Update data store with modified adata
        data_store[data_id]["adata"] = adata
        
        # Ensure result is a dictionary
        if not isinstance(result, dict):
            result = result.dict() if hasattr(result, 'dict') else {"error": "Invalid result format"}
        
        # Add metadata
        result.update({
            "analysis_date": pd.Timestamp.now().isoformat(),
            "n_cells": adata.n_obs,
            "n_neighbors": params.n_neighbors
        })
        
        if context:
            await context.info(f"Analysis completed: {params.analysis_type}")
        
        return SpatialAnalysisResult(
            data_id=data_id,
            analysis_type=params.analysis_type,
            statistics=result,
            result_image=None  # Visualization handled separately
        )
        
    except Exception as e:
        error_msg = f"Error in {params.analysis_type} analysis: {str(e)}"
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {traceback.format_exc()}")
        
        if isinstance(e, (DataNotFoundError, InvalidParameterError, DataCompatibilityError)):
            raise
        else:
            raise ProcessingError(error_msg) from e


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _validate_spatial_data(adata: ad.AnnData) -> None:
    """Validate AnnData object for spatial analysis."""
    if adata.n_obs < 10:
        raise DataNotFoundError("Dataset has too few cells (minimum 10 required)")
    
    if 'spatial' not in adata.obsm:
        raise DataNotFoundError("Dataset missing spatial coordinates in adata.obsm['spatial']")
    
    coords = adata.obsm['spatial']
    if np.any(np.isnan(coords)) or np.any(np.isinf(coords)):
        raise DataCompatibilityError("Spatial coordinates contain NaN or infinite values")


async def _ensure_cluster_key(
    adata: ad.AnnData,
    requested_key: str,
    context: Optional[Context] = None
) -> str:
    """Ensure a valid cluster key exists in adata."""
    if requested_key in adata.obs.columns:
        if not pd.api.types.is_categorical_dtype(adata.obs[requested_key]):
            if context:
                await context.info(f"Converting {requested_key} to categorical...")
            adata.obs[requested_key] = adata.obs[requested_key].astype('category')
        return requested_key
    
    if 'leiden' in adata.obs.columns:
        if context:
            await context.warning(f"Using 'leiden' as fallback for missing '{requested_key}'")
        return 'leiden'
    
    raise ValueError(f"Cluster key '{requested_key}' not found and no 'leiden' clustering available")


async def _ensure_spatial_neighbors(
    adata: ad.AnnData,
    n_neighbors: int,
    context: Optional[Context] = None
) -> None:
    """Ensure spatial neighbors are computed."""
    if 'spatial_neighbors' not in adata.uns:
        if context:
            await context.info(f"Computing spatial neighbors with n_neighbors={n_neighbors}...")
        
        try:
            sq.gr.spatial_neighbors(adata, n_neighs=n_neighbors)
        except Exception as e:
            if context:
                await context.warning(f"Failed to compute spatial neighbors: {e}")
                await context.info("Creating fallback spatial neighbors...")
            
            # Fallback implementation
            coords = adata.obsm['spatial']
            nbrs = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
            distances, indices = nbrs.kneighbors(coords)
            
            n_cells = adata.n_obs
            connectivities = np.zeros((n_cells, n_cells))
            for i in range(n_cells):
                for j in indices[i]:
                    connectivities[i, j] = 1
            
            adata.uns['spatial_neighbors'] = {}
            adata.obsp['spatial_connectivities'] = csr_matrix(connectivities)


def _get_optimal_n_jobs(n_obs: int, requested_n_jobs: Optional[int] = None) -> int:
    """Determine optimal number of parallel jobs based on data size."""
    import os
    
    if requested_n_jobs is not None:
        if requested_n_jobs == -1:
            return os.cpu_count() or 1
        return requested_n_jobs
    
    # Smart defaults based on data size
    if n_obs < 1000:
        return 1  # Single thread for small data
    elif n_obs < 5000:
        return min(2, os.cpu_count() or 1)
    else:
        return min(4, os.cpu_count() or 1)


# ============================================================================
# CORE ANALYSIS FUNCTIONS
# ============================================================================

async def _analyze_morans_i(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Compute Moran's I spatial autocorrelation using squidpy.
    
    Unified implementation using only squidpy for consistency and performance.
    """
    if context:
        await context.info("Running Moran's I spatial autocorrelation analysis...")
    
    # Determine genes to analyze
    if params.moran_genes:
        genes = [g for g in params.moran_genes if g in adata.var_names]
        if not genes:
            raise ValueError(f"None of the specified genes found: {params.moran_genes}")
    else:
        # Use highly variable genes
        if 'highly_variable' not in adata.var or not adata.var['highly_variable'].any():
            raise ValueError("Highly variable genes not found. Please run preprocessing first.")
        genes = adata.var_names[adata.var['highly_variable']][:params.moran_n_genes].tolist()
    
    if context:
        await context.info(f"Analyzing {len(genes)} genes for Moran's I...")
    
    # Optimize parallelization
    n_jobs = _get_optimal_n_jobs(adata.n_obs, params.n_jobs)
    
    # Run spatial autocorrelation
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        genes=genes,
        n_perms=params.moran_n_perms,
        two_tailed=params.moran_two_tailed,
        n_jobs=n_jobs,
        backend=params.backend,
        show_progress_bar=False
    )
    
    # Extract results
    moran_key = 'moranI'
    if moran_key in adata.uns:
        results_df = adata.uns[moran_key]
        
        # Get top significant genes
        significant_genes = results_df[results_df['pval_norm'] < 0.05].index.tolist()
        
        return {
            "n_genes_analyzed": len(genes),
            "n_significant": len(significant_genes),
            "top_positive": results_df.nlargest(10, 'I').index.tolist(),
            "top_negative": results_df.nsmallest(10, 'I').index.tolist(),
            "mean_morans_i": float(results_df['I'].mean()),
            "analysis_key": moran_key
        }
    
    return {"error": "Moran's I computation did not produce results"}


async def _analyze_gearys_c(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Compute Geary's C spatial autocorrelation."""
    if context:
        await context.info("Running Geary's C spatial autocorrelation analysis...")
    
    # Similar to Moran's I but with mode="geary"
    genes = params.moran_genes if params.moran_genes else adata.var_names[:params.moran_n_genes]
    
    sq.gr.spatial_autocorr(
        adata,
        mode="geary",
        genes=genes,
        n_perms=params.moran_n_perms,
        n_jobs=_get_optimal_n_jobs(adata.n_obs, params.n_jobs),
        show_progress_bar=False
    )
    
    # Extract results
    geary_key = 'gearyC'
    if geary_key in adata.uns:
        results_df = adata.uns[geary_key]
        
        return {
            "n_genes_analyzed": len(genes),
            "mean_gearys_c": float(results_df['C'].mean()),
            "analysis_key": geary_key
        }
    
    return {"error": "Geary's C computation did not produce results"}


async def _analyze_neighborhood_enrichment(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Compute neighborhood enrichment analysis."""
    if context:
        await context.info("Running neighborhood enrichment analysis...")
    
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    
    analysis_key = f'{cluster_key}_nhood_enrichment'
    if analysis_key in adata.uns:
        z_scores = adata.uns[analysis_key]['zscore']
        
        return {
            "n_clusters": len(z_scores),
            "max_enrichment": float(np.max(z_scores)),
            "min_enrichment": float(np.min(z_scores)),
            "analysis_key": analysis_key
        }
    
    return {"error": "Neighborhood enrichment did not produce results"}


async def _analyze_co_occurrence(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Compute co-occurrence analysis."""
    if context:
        await context.info("Running co-occurrence analysis...")
    
    sq.gr.co_occurrence(adata, cluster_key=cluster_key)
    
    analysis_key = f'{cluster_key}_co_occurrence'
    if analysis_key in adata.uns:
        co_occurrence = adata.uns[analysis_key]['occ']
        
        return {
            "n_clusters": len(co_occurrence),
            "analysis_key": analysis_key
        }
    
    return {"error": "Co-occurrence analysis did not produce results"}


async def _analyze_ripleys_k(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Compute Ripley's K function."""
    if context:
        await context.info("Running Ripley's K function analysis...")
    
    try:
        sq.gr.ripley(
            adata,
            cluster_key=cluster_key,
            mode='L',  # L-function (variance-stabilized)
            n_simulations=20,
            n_observations=min(1000, adata.n_obs),
            max_dist=None,
            n_steps=50
        )
        
        analysis_key = f'{cluster_key}_ripley_L'
        return {
            "analysis_completed": True,
            "analysis_key": analysis_key
        }
    except Exception as e:
        if context:
            await context.warning(f"Ripley's K analysis failed: {e}")
        return {"error": str(e)}


async def _analyze_getis_ord(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Compute Getis-Ord Gi* hot spot analysis."""
    if context:
        await context.info("Running Getis-Ord Gi* analysis...")
    
    # Determine genes
    if params.getis_ord_genes:
        genes = [g for g in params.getis_ord_genes if g in adata.var_names]
    else:
        if 'highly_variable' not in adata.var:
            raise ValueError("Highly variable genes required for Getis-Ord analysis")
        genes = adata.var_names[adata.var['highly_variable']][:params.getis_ord_n_genes].tolist()
    
    getis_ord_results = {}
    
    try:
        from pysal.lib import weights
        from esda.getisord import G_Local
        
        coords = adata.obsm['spatial']
        w = weights.KNN.from_array(coords, k=params.n_neighbors)
        w.transform = 'r'
        
        for gene in genes:
            if context:
                await context.info(f"Processing gene: {gene}")
            
            y = adata[:, gene].X
            if hasattr(y, 'toarray'):
                y = y.toarray().flatten()
            else:
                y = y.flatten()
            
            local_g = G_Local(y, w, transform='R', star=True)
            
            # Store results in adata.obs
            adata.obs[f"{gene}_getis_ord_z"] = local_g.Zs
            adata.obs[f"{gene}_getis_ord_p"] = local_g.p_sim
            
            getis_ord_results[gene] = {
                'mean_z': float(np.mean(local_g.Zs)),
                'n_hot_spots': int(np.sum(local_g.Zs > 1.96)),
                'n_cold_spots': int(np.sum(local_g.Zs < -1.96))
            }
    
    except ImportError:
        # PySAL dependency missing - fail honestly
        raise ImportError(
            "PySAL (Python Spatial Analysis Library) is required for Getis-Ord analysis but is not installed. "
            "Please install it with: pip install 'libpysal' 'esda'. "
            "Cannot perform spatial autocorrelation analysis without proper statistical methods."
        )
    
    return {
        "n_genes_analyzed": len(getis_ord_results),
        "genes_analyzed": list(getis_ord_results.keys()),
        "results": getis_ord_results
    }


async def _analyze_centrality(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Compute centrality scores."""
    if context:
        await context.info("Computing centrality scores...")
    
    sq.gr.centrality_scores(adata, cluster_key=cluster_key)
    
    analysis_key = f'{cluster_key}_centrality_scores'
    if analysis_key in adata.uns:
        scores = adata.uns[analysis_key]
        
        return {
            "analysis_completed": True,
            "analysis_key": analysis_key,
            "n_clusters": len(scores) if isinstance(scores, dict) else "computed"
        }
    
    return {"error": "Centrality analysis did not produce results"}


# ============================================================================
# ADVANCED ANALYSIS FUNCTIONS (from spatial_statistics.py)
# ============================================================================

async def _analyze_bivariate_moran(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Compute Bivariate Moran's I for gene pairs.
    
    Migrated from spatial_statistics.py
    """
    if context:
        await context.info("Running Bivariate Moran's I analysis...")
    
    # Get gene pairs from parameters
    if not hasattr(params, 'gene_pairs') or not params.gene_pairs:
        # Use top highly variable genes to create pairs
        if 'highly_variable' not in adata.var:
            raise ValueError("Highly variable genes required for bivariate analysis")
        
        hvg = adata.var_names[adata.var['highly_variable']][:10]
        gene_pairs = [(hvg[i], hvg[i+1]) for i in range(0, len(hvg)-1, 2)]
    else:
        gene_pairs = params.gene_pairs
    
    results = {}
    
    try:
        from libpysal.weights import KNN
        
        coords = adata.obsm['spatial']
        w = KNN.from_array(coords, k=params.n_neighbors)
        w.transform = 'R'
        
        for gene1, gene2 in gene_pairs:
            if gene1 in adata.var_names and gene2 in adata.var_names:
                x = adata[:, gene1].X
                y = adata[:, gene2].X
                
                if hasattr(x, 'toarray'):
                    x = x.toarray().flatten()
                    y = y.toarray().flatten()
                else:
                    x = x.flatten()
                    y = y.flatten()
                
                # Compute bivariate Moran's I
                n = len(x)
                x_mean = np.mean(x)
                y_mean = np.mean(y)
                
                numerator = 0
                for i in range(n):
                    for j in range(n):
                        if w.sparse[i, j] > 0:
                            numerator += w.sparse[i, j] * (x[i] - x_mean) * (y[j] - y_mean)
                
                denominator = np.sum((x - x_mean) ** 2)
                
                if denominator > 0:
                    I = (n / w.sparse.sum()) * (numerator / denominator)
                else:
                    I = 0
                
                results[f"{gene1}_vs_{gene2}"] = float(I)
    
    except Exception as e:
        if context:
            await context.warning(f"Bivariate Moran's I failed: {e}")
        return {"error": str(e)}
    
    return {
        "n_pairs_analyzed": len(results),
        "bivariate_morans_i": results,
        "mean_bivariate_i": float(np.mean(list(results.values()))) if results else 0
    }


async def _analyze_join_count(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Compute Join Count statistics for categorical spatial data.
    
    Migrated from spatial_statistics.py
    """
    if context:
        await context.info("Running Join Count analysis...")
    
    try:
        from esda.join_counts import Join_Counts
        from libpysal.weights import KNN
        
        coords = adata.obsm['spatial']
        w = KNN.from_array(coords, k=params.n_neighbors)
        
        # Get categorical data
        y = adata.obs[cluster_key].cat.codes.values
        
        # Compute join counts
        jc = Join_Counts(y, w)
        
        return {
            "bb": float(jc.bb),  # Black-Black joins
            "ww": float(jc.ww),  # White-White joins
            "bw": float(jc.bw),  # Black-White joins
            "J": float(jc.J),    # Total joins
            "p_value": float(jc.p_sim) if hasattr(jc, 'p_sim') else None
        }
    
    except ImportError:
        if context:
            await context.warning("Join Count requires esda package")
        return {"error": "esda package not installed"}
    except Exception as e:
        return {"error": str(e)}


async def _analyze_network_properties(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Analyze network properties of spatial graph.
    
    Migrated from spatial_statistics.py
    """
    if context:
        await context.info("Analyzing network properties...")
    
    try:
        import networkx as nx
        from scipy.sparse import csr_matrix
        
        # Get or create spatial connectivity
        if 'spatial_connectivities' in adata.obsp:
            conn_matrix = adata.obsp['spatial_connectivities']
        else:
            # Create connectivity matrix
            coords = adata.obsm['spatial']
            from sklearn.neighbors import kneighbors_graph
            conn_matrix = kneighbors_graph(coords, n_neighbors=params.n_neighbors, mode='connectivity')
        
        # Convert to networkx graph
        G = nx.from_scipy_sparse_array(conn_matrix)
        
        # Compute properties
        properties = {
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "density": float(nx.density(G)),
            "is_connected": nx.is_connected(G),
            "n_components": nx.number_connected_components(G)
        }
        
        # Additional metrics for connected graphs
        if properties["is_connected"]:
            properties["diameter"] = nx.diameter(G)
            properties["radius"] = nx.radius(G)
        else:
            # Analyze largest component
            largest_cc = max(nx.connected_components(G), key=len)
            G_cc = G.subgraph(largest_cc)
            properties["largest_component_size"] = len(largest_cc)
            properties["largest_component_fraction"] = len(largest_cc) / G.number_of_nodes()
        
        # Clustering coefficient
        try:
            properties["avg_clustering"] = float(nx.average_clustering(G))
        except:
            properties["avg_clustering"] = None
        
        # Degree statistics
        degrees = dict(G.degree())
        degree_values = list(degrees.values())
        properties["degree_mean"] = float(np.mean(degree_values))
        properties["degree_std"] = float(np.std(degree_values))
        
        return properties
    
    except ImportError:
        if context:
            await context.warning("NetworkX not installed")
        return {"error": "networkx package required"}
    except Exception as e:
        return {"error": str(e)}


async def _analyze_spatial_centrality(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Compute various centrality measures for spatial network.
    
    Migrated from spatial_statistics.py
    """
    if context:
        await context.info("Computing spatial centrality measures...")
    
    try:
        import networkx as nx
        
        # Get connectivity matrix
        if 'spatial_connectivities' in adata.obsp:
            conn_matrix = adata.obsp['spatial_connectivities']
        else:
            coords = adata.obsm['spatial']
            from sklearn.neighbors import kneighbors_graph
            conn_matrix = kneighbors_graph(coords, n_neighbors=params.n_neighbors, mode='connectivity')
        
        # Convert to networkx
        G = nx.from_scipy_sparse_array(conn_matrix)
        
        # Compute centrality measures
        degree_centrality = nx.degree_centrality(G)
        closeness_centrality = nx.closeness_centrality(G)
        betweenness_centrality = nx.betweenness_centrality(G)
        
        # Store in adata.obs
        adata.obs['degree_centrality'] = pd.Series(degree_centrality)
        adata.obs['closeness_centrality'] = pd.Series(closeness_centrality)
        adata.obs['betweenness_centrality'] = pd.Series(betweenness_centrality)
        
        # Compute statistics by cluster
        centrality_stats = {}
        for cluster in adata.obs[cluster_key].unique():
            mask = adata.obs[cluster_key] == cluster
            centrality_stats[str(cluster)] = {
                "mean_degree": float(adata.obs.loc[mask, 'degree_centrality'].mean()),
                "mean_closeness": float(adata.obs.loc[mask, 'closeness_centrality'].mean()),
                "mean_betweenness": float(adata.obs.loc[mask, 'betweenness_centrality'].mean())
            }
        
        return {
            "centrality_computed": True,
            "cluster_centrality": centrality_stats,
            "global_stats": {
                "mean_degree": float(np.mean(list(degree_centrality.values()))),
                "mean_closeness": float(np.mean(list(closeness_centrality.values()))),
                "mean_betweenness": float(np.mean(list(betweenness_centrality.values())))
            }
        }
    
    except ImportError:
        return {"error": "NetworkX required for centrality analysis"}
    except Exception as e:
        return {"error": str(e)}


# ============================================================================
# DEEP LEARNING ANALYSIS (preserved from original)
# ============================================================================

# Import scvi-tools for advanced spatial analysis
try:
    import scvi
    try:
        from scvi.external import MRVI  # Multi-resolution Variational Inference for spatial data
        SCVIVA = MRVI  # Use MRVI as the spatial analysis method
    except ImportError:
        SCVIVA = None
except ImportError:
    scvi = None
    SCVIVA = None


async def _analyze_with_scviva(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Analyze spatial data using SCVIVA deep learning model.
    
    ⚠️ TODO: UNDER DEVELOPMENT - DO NOT USE IN PRODUCTION
    Known Issues:
    - Produces NaN values during training with preprocessed data
    - Numerical instability in the initial embedding generation
    - Needs proper data normalization handling for scVI
    
    Status: Disabled for user access until issues are resolved
    """
    # TODO: Temporarily disabled until NaN issue is resolved
    return {
        "error": "SCVIVA analysis is currently under development and disabled. "
                "Please use other spatial analysis methods."
    }
    
    if scvi is None:
        return {"error": "scvi-tools package required for SCVIVA analysis"}
    
    if SCVIVA is None:
        return {"error": "MRVI spatial analysis requires scvi-tools >= 1.3.0"}
    
    if context:
        await context.info("Running SCVIVA deep learning analysis...")
    
    # Parameters from params object
    n_epochs = getattr(params, 'scviva_n_epochs', 1000)
    n_hidden = getattr(params, 'scviva_n_hidden', 128)
    n_latent = getattr(params, 'scviva_n_latent', 10)
    use_gpu = getattr(params, 'scviva_use_gpu', False)
    
    try:
        adata_copy = adata.copy()
        
        # Add required columns if missing
        if 'sample' not in adata_copy.obs.columns:
            adata_copy.obs['sample'] = 'sample_1'
        
        if 'cell_type' not in adata_copy.obs.columns:
            adata_copy.obs['cell_type'] = 'Unknown'
        
        # Generate SCVI embeddings if needed
        if 'X_scVI' not in adata_copy.obsm:
            if context:
                await context.info("Generating SCVI embeddings...")
            
            scvi.model.SCVI.setup_anndata(adata_copy, batch_key='sample')
            scvi_model = scvi.model.SCVI(adata_copy, n_hidden=n_hidden, n_latent=n_latent)
            scvi_model.train(max_epochs=min(50, n_epochs), early_stopping=True)
            adata_copy.obsm['X_scVI'] = scvi_model.get_latent_representation()
        
        # For now, use standard scVI with spatial-aware setup
        # MRVI requires specific setup that might differ from SCVIVA
        if context:
            await context.info("Using scVI for spatial analysis (MRVI/SCVIVA not available)")
        
        # Setup scVI with spatial awareness
        scvi.model.SCVI.setup_anndata(
            adata_copy,
            batch_key='sample',
            categorical_covariate_keys=['cell_type'] if 'cell_type' in adata_copy.obs.columns else None
        )
        
        # Create and train model
        model = scvi.model.SCVI(
            adata_copy, 
            n_hidden=n_hidden, 
            n_latent=n_latent,
            gene_likelihood="zinb"  # Zero-inflated negative binomial
        )
        
        # Train the model (use_gpu parameter removed in newer scvi versions)
        train_kwargs = {
            "max_epochs": n_epochs,
            "early_stopping": True
        }
        if use_gpu:
            train_kwargs["accelerator"] = "gpu"
        
        model.train(**train_kwargs)
        
        # Get latent representation
        latent = model.get_latent_representation()
        
        # Store results back
        adata.obsm['X_scviva_latent'] = latent
        
        return {
            "method": "scVI (spatial-aware)",
            "n_latent_dims": n_latent,
            "n_epochs": n_epochs,
            "latent_variance_explained": float(np.var(latent, axis=0).sum()),
            "training_completed": True,
            "device": "GPU" if use_gpu else "CPU",
            "note": "Using standard scVI with spatial context"
        }
    
    except Exception as e:
        if context:
            await context.error(f"SCVIVA analysis failed: {e}")
        return {"error": str(e)}


# ============================================================================
# BACKWARD COMPATIBILITY
# ============================================================================

async def analyze_spatial_with_scviva(
    spatial_adata,
    n_epochs: int = 1000,
    n_hidden: int = 128,
    n_latent: int = 10,
    use_gpu: bool = False,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Backward compatibility wrapper for SCVIVA analysis."""
    params = SpatialAnalysisParameters(analysis_type="scviva")
    params.scviva_n_epochs = n_epochs
    params.scviva_n_hidden = n_hidden
    params.scviva_n_latent = n_latent
    params.scviva_use_gpu = use_gpu
    
    return await _analyze_with_scviva(spatial_adata, params, context)


async def calculate_spatial_stats(
    data_id: str,
    data_store: Dict[str, Any],
    feature: str,
    statistic: str = "gearys_c",
    n_neighbors: int = 6,
    context = None
) -> Dict[str, Any]:
    """
    Calculate specialized spatial statistics for a single feature/gene.
    
    This function provides advanced spatial statistics not available in the main
    spatial analysis tool. For Moran's I analysis, use analyze_spatial_patterns()
    with analysis_type="moran".
    
    Parameters
    ----------
    data_id : str
        Dataset identifier
    data_store : Dict[str, Any]
        Data storage dictionary
    feature : str
        Gene or feature to analyze
    statistic : str
        Type of statistic to calculate (gearys_c, local_morans)
    n_neighbors : int
        Number of neighbors for spatial graph
    context : Optional
        MCP context
    
    Returns
    -------
    Dict[str, Any]
        Statistics results
    """
    # Get data
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")
    
    adata = data_store[data_id]["adata"]
    
    # Ensure spatial graph exists
    _ensure_spatial_graph(adata, n_neighbors)
    
    # Calculate statistics based on type
    if statistic == "gearys_c":
        # Calculate Geary's C using squidpy
        import squidpy as sq
        
        # Ensure feature exists
        if feature not in adata.var_names:
            raise ValueError(f"Feature {feature} not found in dataset")
        
        # Calculate Geary's C
        sq.gr.spatial_autocorr(
            adata,
            mode="geary",
            genes=[feature],
            n_jobs=1
        )
        
        # Extract results
        if "gearyC" in adata.uns and feature in adata.uns["gearyC"]:
            geary_c = float(adata.uns["gearyC"][feature])
            pval = float(adata.uns["gearyC_pval"][feature]) if "gearyC_pval" in adata.uns else None
            
            return {
                "statistic": "gearys_c",
                "feature": feature,
                "value": geary_c,
                "pvalue": pval,
                "interpretation": "Values < 1 indicate positive spatial autocorrelation",
                "n_neighbors": n_neighbors
            }
    
    elif statistic == "local_morans":
        # Calculate Local Moran's I
        import numpy as np
        from scipy.stats import zscore
        from scipy.sparse import issparse
        
        # Ensure feature exists
        if feature not in adata.var_names:
            raise ValueError(f"Feature {feature} not found in dataset")
        
        # Get expression values
        if issparse(adata.X):
            expr = adata[:, feature].X.toarray().flatten()
        else:
            expr = adata[:, feature].X.flatten()
        
        # Standardize expression
        z_expr = zscore(expr)
        
        # Get spatial weights
        W = adata.obsp["spatial_connectivities"]
        
        # Calculate local Moran's I
        n = len(expr)
        local_i = np.zeros(n)
        
        for i in range(n):
            neighbors = W[i].nonzero()[1]
            if len(neighbors) > 0:
                local_i[i] = z_expr[i] * np.mean(z_expr[neighbors])
        
        # Store in adata
        adata.obs[f"{feature}_local_morans"] = local_i
        
        return {
            "statistic": "local_morans",
            "feature": feature,
            "values": local_i.tolist(),
            "mean": float(np.mean(local_i)),
            "std": float(np.std(local_i)),
            "min": float(np.min(local_i)),
            "max": float(np.max(local_i)),
            "n_neighbors": n_neighbors
        }
    
    else:
        raise ValueError(f"Unsupported statistic: {statistic}")
    
    return {
        "error": f"Failed to calculate {statistic} for {feature}"
    }