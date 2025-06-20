"""
Enrichment analysis tool for spatial transcriptomics data using EnrichMap.

This module provides spatially-aware gene set enrichment analysis capabilities
using the EnrichMap package.
"""

from typing import Dict, List, Optional, Union, Any, Tuple
import logging
import numpy as np
from pathlib import Path
import sys
import os

# Add EnrichMap to Python path
current_dir = Path(__file__).parent
project_root = current_dir.parent.parent
enrichmap_path = os.path.join(project_root, "third_party", "EnrichMap")
sys.path.insert(0, enrichmap_path)

from mcp.server.fastmcp import Context
from ..utils.error_handling import ProcessingError

logger = logging.getLogger(__name__)


def is_enrichmap_available() -> Tuple[bool, str]:
    """Check if EnrichMap is available and all dependencies are met."""
    try:
        import enrichmap as em
        
        # Check for required dependencies
        # Map package names to their import names
        module_mapping = {
            'scanpy': 'scanpy',
            'squidpy': 'squidpy', 
            'scipy': 'scipy',
            'scikit-learn': 'sklearn',
            'statsmodels': 'statsmodels',
            'pygam': 'pygam',
            'scikit-gstat': 'skgstat',
            'adjustText': 'adjustText',
            'splot': 'splot'
        }
        
        missing = []
        for package, module in module_mapping.items():
            try:
                __import__(module)
            except ImportError:
                missing.append(package)
        
        if missing:
            return False, f"Missing EnrichMap dependencies: {', '.join(missing)}"
        
        return True, ""
    except ImportError:
        return False, "EnrichMap not found. Please check third_party/EnrichMap installation"


async def perform_enrichment_analysis(
    data_id: str,
    data_store: Dict[str, Any],
    gene_sets: Union[List[str], Dict[str, List[str]]],
    score_keys: Optional[Union[str, List[str]]] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = True,
    batch_key: Optional[str] = None,
    gene_weights: Optional[Dict[str, Dict[str, float]]] = None,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Perform spatially-aware gene set enrichment analysis using EnrichMap.
    
    Parameters
    ----------
    data_id : str
        Identifier for the spatial data in the data store
    data_store : Dict[str, Any]
        Dictionary containing the data
    gene_sets : Union[List[str], Dict[str, List[str]]]
        Either a single gene list or a dictionary of gene sets where keys are 
        signature names and values are lists of genes
    score_keys : Optional[Union[str, List[str]]]
        Names for the gene signatures if gene_sets is a list. Ignored if gene_sets 
        is already a dictionary
    spatial_key : str
        Key in adata.obsm containing spatial coordinates (default: "spatial")
    n_neighbors : int
        Number of nearest spatial neighbors for smoothing (default: 6)
    smoothing : bool
        Whether to perform spatial smoothing (default: True)
    correct_spatial_covariates : bool
        Whether to correct for spatial covariates using GAM (default: True)
    batch_key : Optional[str]
        Column in adata.obs for batch-wise normalization
    gene_weights : Optional[Dict[str, Dict[str, float]]]
        Pre-computed gene weights for each signature
    context : Optional[Context]
        Execution context
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - data_id: ID of the data with enrichment scores
        - signatures: List of computed signatures
        - score_columns: List of column names containing scores
        - gene_contributions: Dictionary of gene contributions per signature
        - summary_stats: Summary statistics for each signature
    """
    # Check if EnrichMap is available
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")
    
    # Import EnrichMap
    import enrichmap as em
    
    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")
    
    adata = data_store[data_id]["adata"]
    
    # Validate spatial coordinates
    if spatial_key not in adata.obsm:
        raise ProcessingError(f"Spatial coordinates '{spatial_key}' not found in adata.obsm")
    
    # Convert single gene list to dictionary format
    if isinstance(gene_sets, list):
        if score_keys is None:
            score_keys = "enrichmap_signature"
        gene_sets = {score_keys: gene_sets}
    
    # Validate gene sets
    available_genes = set(adata.var_names)
    validated_gene_sets = {}
    
    for sig_name, genes in gene_sets.items():
        common_genes = list(set(genes).intersection(available_genes))
        if len(common_genes) < 2:
            logger.warning(f"Signature '{sig_name}' has {len(common_genes)} genes in the dataset. Skipping.")
            continue
        validated_gene_sets[sig_name] = common_genes
        logger.info(f"Signature '{sig_name}': {len(common_genes)}/{len(genes)} genes found")
    
    if not validated_gene_sets:
        raise ProcessingError("No valid gene signatures found with at least 2 genes")
    
    # Run EnrichMap scoring
    try:
        em.tl.score(
            adata=adata,
            gene_set=validated_gene_sets,
            gene_weights=gene_weights,
            score_key=None,  # Not needed as we pass dict
            spatial_key=spatial_key,
            n_neighbors=n_neighbors,
            smoothing=smoothing,
            correct_spatial_covariates=correct_spatial_covariates,
            batch_key=batch_key
        )
    except Exception as e:
        raise ProcessingError(f"EnrichMap scoring failed: {str(e)}")
    
    # Collect results
    score_columns = [f"{sig}_score" for sig in validated_gene_sets.keys()]
    
    # Calculate summary statistics
    summary_stats = {}
    for sig_name in validated_gene_sets.keys():
        score_col = f"{sig_name}_score"
        scores = adata.obs[score_col]
        
        summary_stats[sig_name] = {
            "mean": float(scores.mean()),
            "std": float(scores.std()),
            "min": float(scores.min()),
            "max": float(scores.max()),
            "median": float(scores.median()),
            "q25": float(scores.quantile(0.25)),
            "q75": float(scores.quantile(0.75)),
            "n_genes": len(validated_gene_sets[sig_name])
        }
    
    # Get gene contributions
    gene_contributions = {}
    if "gene_contributions" in adata.uns:
        gene_contributions = {
            sig: {gene: float(contrib.mean()) for gene, contrib in contribs.items()}
            for sig, contribs in adata.uns["gene_contributions"].items()
        }
    
    return {
        "data_id": data_id,
        "signatures": list(validated_gene_sets.keys()),
        "score_columns": score_columns,
        "gene_contributions": gene_contributions,
        "summary_stats": summary_stats,
        "parameters": {
            "n_neighbors": n_neighbors,
            "smoothing": smoothing,
            "correct_spatial_covariates": correct_spatial_covariates,
            "batch_key": batch_key
        }
    }


async def compute_spatial_metrics(
    data_id: str,
    data_store: Dict[str, Any],
    score_key: str,
    metrics: Optional[List[str]] = None,
    n_neighbors: int = 6,
    n_perms: int = 999,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Compute spatial metrics for enrichment scores.
    
    Parameters
    ----------
    data_id : str
        Identifier for the spatial data
    data_store : Dict[str, Any]
        Dictionary containing the data
    score_key : str
        Column name containing the enrichment scores
    metrics : Optional[List[str]]
        List of metrics to compute. Options: ['morans_i', 'getis_ord', 'variance']
        If None, computes all metrics
    n_neighbors : int
        Number of spatial neighbors (default: 6)
    n_perms : int
        Number of permutations for significance testing (default: 999)
    context : Optional[Context]
        Execution context
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing computed spatial metrics
    """
    # Check if EnrichMap is available
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")
    
    import enrichmap as em
    
    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")
    
    adata = data_store[data_id]["adata"]
    
    # Validate score column
    if score_key not in adata.obs.columns:
        raise ProcessingError(f"Score column '{score_key}' not found in adata.obs")
    
    # Default metrics
    if metrics is None:
        metrics = ['morans_i', 'getis_ord', 'variance']
    
    try:
        # Compute spatial metrics
        result = em.tl.compute_spatial_metrics(
            adata=adata,
            score_keys=[score_key],
            metrics=metrics,
            n_neighs=n_neighbors,
            n_perms=n_perms
        )
        
        # Extract results for the single score key
        metric_results = {}
        for metric in metrics:
            if metric in result:
                metric_results[metric] = {
                    "value": float(result[metric][score_key]),
                    "p_value": float(result.get(f"{metric}_pval", {}).get(score_key, np.nan))
                }
        
        return {
            "data_id": data_id,
            "score_key": score_key,
            "metrics": metric_results,
            "parameters": {
                "n_neighbors": n_neighbors,
                "n_perms": n_perms
            }
        }
        
    except Exception as e:
        raise ProcessingError(f"Spatial metrics computation failed: {str(e)}")


async def cluster_gene_correlation(
    data_id: str,
    data_store: Dict[str, Any],
    signature_name: str,
    cluster_key: str = "leiden",
    correlation_method: str = "pearson",
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """
    Compute correlation between gene expression and enrichment scores per cluster.
    
    Parameters
    ----------
    data_id : str
        Identifier for the spatial data
    data_store : Dict[str, Any]
        Dictionary containing the data
    signature_name : str
        Name of the signature to analyze
    cluster_key : str
        Column in adata.obs containing cluster labels (default: "leiden")
    correlation_method : str
        Correlation method: 'pearson' or 'spearman' (default: "pearson")
    context : Optional[Context]
        Execution context
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing correlation results per cluster
    """
    # Check if EnrichMap is available
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")
    
    import enrichmap as em
    
    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")
    
    adata = data_store[data_id]["adata"]
    
    # Validate inputs
    score_key = f"{signature_name}_score"
    if score_key not in adata.obs.columns:
        raise ProcessingError(f"Score column '{score_key}' not found. Run enrichment analysis first.")
    
    if cluster_key not in adata.obs.columns:
        raise ProcessingError(f"Cluster column '{cluster_key}' not found in adata.obs")
    
    try:
        # Compute cluster-gene correlations
        result = em.tl.cluster_gene_correlation(
            adata=adata,
            signature_name=signature_name,
            cluster_key=cluster_key,
            correlation_method=correlation_method
        )
        
        # Convert results to serializable format
        correlation_results = {}
        for cluster, corr_df in result.items():
            correlation_results[str(cluster)] = {
                "genes": corr_df.index.tolist(),
                "correlations": corr_df["correlation"].tolist(),
                "top_positive": corr_df.nlargest(10, "correlation")[["correlation"]].to_dict(),
                "top_negative": corr_df.nsmallest(10, "correlation")[["correlation"]].to_dict()
            }
        
        return {
            "data_id": data_id,
            "signature_name": signature_name,
            "cluster_key": cluster_key,
            "correlation_method": correlation_method,
            "cluster_correlations": correlation_results
        }
        
    except Exception as e:
        raise ProcessingError(f"Cluster gene correlation failed: {str(e)}")