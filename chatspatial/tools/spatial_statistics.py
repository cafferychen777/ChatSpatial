"""
A module for quantitative spatial analysis of spatial transcriptomics data.

This module provides a collection of functions to compute various spatial
statistics. It includes methods for assessing global and local spatial
autocorrelation, analyzing neighborhood compositions, and evaluating spatial
patterns of cell clusters.

Key functionalities include:
- Global spatial autocorrelation (Moran's I, Geary's C).
- Local spatial autocorrelation (Local Moran's I / LISA for cluster detection).
- Local spatial statistics for hotspot detection (Getis-Ord Gi*).
- Cluster-based analysis (Neighborhood Enrichment, Co-occurrence, Ripley's K).
- Spatial network analysis (Centrality Scores, Network Properties).
- Bivariate spatial correlation analysis (Bivariate Moran's I).
- Categorical spatial analysis (Join Count statistics).
- Spatial centrality measures for tissue architecture.

The primary entry point is the `analyze_spatial_patterns` function, which
dispatches tasks to the appropriate analysis function based on user parameters.
All 12 analysis types are accessible through this unified interface with a
new unified 'genes' parameter for consistent gene selection across methods.
"""

import logging
import traceback
from typing import Any, Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq
from mcp.server.fastmcp import Context
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors

from ..models.analysis import SpatialAnalysisResult
from ..models.data import SpatialAnalysisParameters
from ..utils.error_handling import (DataCompatibilityError, DataNotFoundError,
                                    InvalidParameterError, ProcessingError)
# Import standardized utilities

logger = logging.getLogger(__name__)


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================


async def analyze_spatial_patterns(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Optional[Context] = None,
) -> SpatialAnalysisResult:
    """
    Serves as the central dispatcher for executing various spatial analysis methods.

    This function validates the input data, computes a spatial neighbor graph if one
    does not exist, and routes the analysis to the appropriate specialized function
    based on the `analysis_type` parameter. Results from the analysis are added to
    the `AnnData` object within the data store. Note that visualization is handled
    by a separate function.

    Parameters
    ----------
    data_id : str
        The identifier for the dataset.
    data_store : Dict[str, Any]
        A dictionary that stores the loaded datasets.
    params : SpatialAnalysisParameters
        An object containing the parameters for the analysis, including the
        specific `analysis_type` to perform.
    context : Optional[Context]
        The MCP context for logging and communication.

    Returns
    -------
    SpatialAnalysisResult
        An object containing the statistical results and metadata from the analysis.

    Raises
    ------
    DataNotFoundError
        If the specified dataset is not found in the data store.
    InvalidParameterError
        If the provided parameters are not valid for the requested analysis.
    ProcessingError
        If an error occurs during the execution of the analysis.
    """
    # Validate parameters
    supported_types = [
        "neighborhood",
        "co_occurrence",
        "ripley",
        "moran",
        "local_moran",  # Added Local Moran's I
        "geary",
        "centrality",
        "getis_ord",
        "bivariate_moran",
        "join_count",
        "network_properties",
        "spatial_centrality",
    ]

    if params.analysis_type not in supported_types:
        raise InvalidParameterError(
            f"Unsupported analysis type: {params.analysis_type}"
        )

    if params.n_neighbors <= 0:
        raise InvalidParameterError(
            f"n_neighbors must be positive, got {params.n_neighbors}"
        )

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
        elif params.analysis_type == "local_moran":
            result = await _analyze_local_moran(adata, params, context)
        elif params.analysis_type == "geary":
            result = await _analyze_gearys_c(adata, params, context)
        elif params.analysis_type == "neighborhood":
            result = await _analyze_neighborhood_enrichment(
                adata, cluster_key, params, context
            )
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
            result = await _analyze_network_properties(
                adata, cluster_key, params, context
            )
        elif params.analysis_type == "spatial_centrality":
            result = await _analyze_spatial_centrality(
                adata, cluster_key, params, context
            )
        else:
            raise ValueError(f"Analysis type {params.analysis_type} not implemented")

        # Update data store with modified adata
        data_store[data_id]["adata"] = adata

        # Ensure result is a dictionary
        if not isinstance(result, dict):
            result = (
                result.dict()
                if hasattr(result, "dict")
                else {"error": "Invalid result format"}
            )

        # Add metadata
        result.update(
            {
                "analysis_date": pd.Timestamp.now().isoformat(),
                "n_cells": adata.n_obs,
                "n_neighbors": params.n_neighbors,
            }
        )

        if context:
            await context.info(f"Analysis completed: {params.analysis_type}")

        return SpatialAnalysisResult(
            data_id=data_id,
            analysis_type=params.analysis_type,
            statistics=result,
            result_image=None,  # Visualization handled separately
        )

    except Exception as e:
        error_msg = f"Error in {params.analysis_type} analysis: {str(e)}"
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {traceback.format_exc()}")

        if isinstance(
            e, (DataNotFoundError, InvalidParameterError, DataCompatibilityError)
        ):
            raise
        else:
            raise ProcessingError(error_msg) from e


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def _validate_spatial_data(adata: ad.AnnData) -> None:
    """
    Checks if the AnnData object meets minimum requirements for spatial analysis.

    This function performs the following checks:
    1. Ensures the dataset contains a minimum number of cells (10).
    2. Verifies the presence of spatial coordinates in `adata.obsm['spatial']`.
    3. Confirms that spatial coordinates do not contain NaN or infinite values.
    """
    if adata.n_obs < 10:
        raise DataNotFoundError("Dataset has too few cells (minimum 10 required)")

    if "spatial" not in adata.obsm:
        raise DataNotFoundError(
            "Dataset missing spatial coordinates in adata.obsm['spatial']"
        )

    coords = adata.obsm["spatial"]
    if np.any(np.isnan(coords)) or np.any(np.isinf(coords)):
        raise DataCompatibilityError(
            "Spatial coordinates contain NaN or infinite values"
        )


async def _ensure_cluster_key(
    adata: ad.AnnData, requested_key: str, context: Optional[Context] = None
) -> str:
    """Ensure a valid cluster key exists in adata."""
    if requested_key in adata.obs.columns:
        if not pd.api.types.is_categorical_dtype(adata.obs[requested_key]):
            if context:
                await context.info(f"Converting {requested_key} to categorical...")
            adata.obs[requested_key] = adata.obs[requested_key].astype("category")
        return requested_key

    if "leiden" in adata.obs.columns:
        if context:
            await context.warning(
                f"Using 'leiden' as fallback for missing '{requested_key}'"
            )
        return "leiden"

    raise ValueError(
        f"Cluster key '{requested_key}' not found and no 'leiden' clustering available"
    )


async def _ensure_spatial_neighbors(
    adata: ad.AnnData, n_neighbors: int, context: Optional[Context] = None
) -> None:
    """
    Ensure spatial neighbors are computed using squidpy's scientifically validated methods.
    
    This function strictly uses squidpy.gr.spatial_neighbors, which is specifically
    designed for spatial transcriptomics data and supports advanced spatial topologies
    like Delaunay triangulation, adaptive neighborhood sizes, and biologically-relevant
    spatial relationships.
    
    IMPORTANT: This function will NOT fallback to sklearn.NearestNeighbors as the two
    methods are scientifically incompatible and could lead to incorrect spatial analysis
    results. If squidpy fails, the analysis should be terminated with proper error guidance.
    """
    if "spatial_neighbors" not in adata.uns:
        if context:
            await context.info(
                f"Computing spatial neighbors with squidpy (n_neighbors={n_neighbors})..."
            )

        try:
            # Use squidpy's spatial-transcriptomics optimized method
            sq.gr.spatial_neighbors(
                adata, 
                n_neighs=n_neighbors,
                coord_type="generic",  # Works for all spatial data types
                set_diag=False,        # Exclude self-neighbors
                key_added="spatial"    # Standard key for spatial neighbors
            )
            
            if context:
                await context.info("✅ Spatial neighbors computed successfully with squidpy")
                
        except Exception as e:
            error_msg = (
                f"❌ CRITICAL: Spatial neighbor computation failed: {e}\n\n"
                "🔬 SCIENTIFIC INTEGRITY NOTICE:\n"
                "Spatial neighbor graphs are fundamental to all spatial transcriptomics analyses.\n"
                "Using alternative methods (like sklearn) would compromise result validity.\n\n"
                "💡 SOLUTIONS:\n"
                "1. Check squidpy installation: pip install 'squidpy>=1.3.0'\n"
                "2. Verify spatial coordinates in adata.obsm['spatial']\n"
                "3. Ensure coordinates don't contain NaN/infinite values\n"
                "4. For Visium data, try coord_type='grid'\n"
                "5. Reduce n_neighbors if dataset is very small\n\n"
                "🚫 Analysis cannot proceed without proper spatial neighbors."
            )
            
            if context:
                await context.error(error_msg)
            
            raise ProcessingError(
                f"Spatial neighbor computation failed. Cannot proceed with spatial analysis. "
                f"Original error: {e}"
            ) from e


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
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """
    Calculates Moran's I to measure global spatial autocorrelation for genes.

    Moran's I is a statistic that indicates whether the expression of a gene is
    spatially clustered, dispersed, or randomly distributed.
    - A value near +1.0 indicates strong clustering of similar expression values.
    - A value near -1.0 indicates dispersion (a checkerboard-like pattern).
    - A value near 0 indicates a random spatial distribution.

    The analysis is performed on highly variable genes by default, but a
    specific gene list can be provided.
    """
    if context:
        await context.info("Running Moran's I spatial autocorrelation analysis...")

    # Determine genes to analyze (unified gene selection)
    if params.genes:
        # Use specific genes
        genes = [g for g in params.genes if g in adata.var_names]
        if not genes:
            raise ValueError(f"None of the specified genes found: {params.genes}")
    else:
        # Use highly variable genes
        if "highly_variable" not in adata.var or not adata.var["highly_variable"].any():
            raise ValueError(
                "Highly variable genes not found. Please run preprocessing first."
            )
        genes = adata.var_names[adata.var["highly_variable"]][
            : params.moran_n_genes
        ].tolist()

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
        show_progress_bar=False,
    )

    # Extract results
    moran_key = "moranI"
    if moran_key in adata.uns:
        results_df = adata.uns[moran_key]

        # Get top significant genes
        significant_genes = results_df[results_df["pval_norm"] < 0.05].index.tolist()

        return {
            "n_genes_analyzed": len(genes),
            "n_significant": len(significant_genes),
            "top_positive": results_df.nlargest(10, "I").index.tolist(),
            "top_negative": results_df.nsmallest(10, "I").index.tolist(),
            "mean_morans_i": float(results_df["I"].mean()),
            "analysis_key": moran_key,
        }

    return {"error": "Moran's I computation did not produce results"}


async def _analyze_gearys_c(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """Compute Geary's C spatial autocorrelation."""
    if context:
        await context.info("Running Geary's C spatial autocorrelation analysis...")

    # Determine genes to analyze (unified gene selection)
    if params.genes:
        # Use specific genes
        genes = [g for g in params.genes if g in adata.var_names]
        if not genes:
            raise ValueError(f"None of the specified genes found: {params.genes}")
    else:
        # Use highly variable genes
        if "highly_variable" in adata.var and adata.var["highly_variable"].any():
            genes = adata.var_names[adata.var["highly_variable"]][
                : params.moran_n_genes
            ].tolist()
        else:
            # Fallback to top genes
            genes = adata.var_names[: params.moran_n_genes].tolist()

    sq.gr.spatial_autocorr(
        adata,
        mode="geary",
        genes=genes,
        n_perms=params.moran_n_perms,
        n_jobs=_get_optimal_n_jobs(adata.n_obs, params.n_jobs),
        show_progress_bar=False,
    )

    # Extract results (squidpy returns DataFrame, not dict)
    geary_key = "gearyC"
    if geary_key in adata.uns:
        import pandas as pd
        results_df = adata.uns[geary_key]
        if isinstance(results_df, pd.DataFrame):
            return {
                "n_genes_analyzed": len(genes),
                "mean_gearys_c": float(results_df["C"].mean()),
                "analysis_key": geary_key,
            }

    return {"error": "Geary's C computation did not produce results"}


async def _analyze_neighborhood_enrichment(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """Compute neighborhood enrichment analysis."""
    if context:
        await context.info("Running neighborhood enrichment analysis...")

    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)

    analysis_key = f"{cluster_key}_nhood_enrichment"
    if analysis_key in adata.uns:
        z_scores = adata.uns[analysis_key]["zscore"]

        return {
            "n_clusters": len(z_scores),
            "max_enrichment": float(np.max(z_scores)),
            "min_enrichment": float(np.min(z_scores)),
            "analysis_key": analysis_key,
        }

    return {"error": "Neighborhood enrichment did not produce results"}


async def _analyze_co_occurrence(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """Compute co-occurrence analysis."""
    if context:
        await context.info("Running co-occurrence analysis...")

    sq.gr.co_occurrence(adata, cluster_key=cluster_key)

    analysis_key = f"{cluster_key}_co_occurrence"
    if analysis_key in adata.uns:
        co_occurrence = adata.uns[analysis_key]["occ"]

        return {"n_clusters": len(co_occurrence), "analysis_key": analysis_key}

    return {"error": "Co-occurrence analysis did not produce results"}


async def _analyze_ripleys_k(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """Compute Ripley's K function."""
    if context:
        await context.info("Running Ripley's K function analysis...")

    try:
        sq.gr.ripley(
            adata,
            cluster_key=cluster_key,
            mode="L",  # L-function (variance-stabilized)
            n_simulations=20,
            n_observations=min(1000, adata.n_obs),
            max_dist=None,
            n_steps=50,
        )

        analysis_key = f"{cluster_key}_ripley_L"
        return {"analysis_completed": True, "analysis_key": analysis_key}
    except Exception as e:
        if context:
            await context.warning(f"Ripley's K analysis failed: {e}")
        return {"error": str(e)}


async def _analyze_getis_ord(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """
    Performs Getis-Ord Gi* analysis to identify local spatial clusters.

    This method identifies statistically significant hot spots (clusters of high
    gene expression) and cold spots (clusters of low gene expression). It computes
    a Z-score for each spot, where high positive Z-scores indicate hot spots and
    low negative Z-scores indicate cold spots. This analysis requires the `esda`
    and `libpysal` libraries.
    """
    if context:
        await context.info("Running Getis-Ord Gi* analysis...")

    # Determine genes to analyze (unified gene selection)
    if params.genes:
        # Use specific genes
        genes = [g for g in params.genes if g in adata.var_names]
        if not genes:
            raise ValueError(f"None of the specified genes found: {params.genes}")
    else:
        # Use highly variable genes
        if "highly_variable" not in adata.var:
            raise ValueError("Highly variable genes required for Getis-Ord analysis")
        genes = adata.var_names[adata.var["highly_variable"]][
            : params.getis_ord_n_genes
        ].tolist()

    getis_ord_results = {}

    try:
        from esda.getisord import G_Local
        from pysal.lib import weights

        coords = adata.obsm["spatial"]
        w = weights.KNN.from_array(coords, k=params.n_neighbors)
        w.transform = "r"

        for gene in genes:
            if context:
                await context.info(f"Processing gene: {gene}")

            y = adata[:, gene].X
            if hasattr(y, "toarray"):
                y = y.toarray().flatten()
            else:
                y = y.flatten()
            
            # Ensure y is float64 to match weights matrix dtype
            y = y.astype(np.float64)

            local_g = G_Local(y, w, transform="R", star=True)

            # Store results in adata.obs
            adata.obs[f"{gene}_getis_ord_z"] = local_g.Zs
            adata.obs[f"{gene}_getis_ord_p"] = local_g.p_sim

            getis_ord_results[gene] = {
                "mean_z": float(np.mean(local_g.Zs)),
                "n_hot_spots": int(np.sum(local_g.Zs > 1.96)),
                "n_cold_spots": int(np.sum(local_g.Zs < -1.96)),
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
        "results": getis_ord_results,
    }


async def _analyze_centrality(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """Compute centrality scores."""
    if context:
        await context.info("Computing centrality scores...")

    sq.gr.centrality_scores(adata, cluster_key=cluster_key)

    analysis_key = f"{cluster_key}_centrality_scores"
    if analysis_key in adata.uns:
        scores = adata.uns[analysis_key]

        return {
            "analysis_completed": True,
            "analysis_key": analysis_key,
            "n_clusters": len(scores) if isinstance(scores, dict) else "computed",
        }

    return {"error": "Centrality analysis did not produce results"}


# ============================================================================
# ADVANCED ANALYSIS FUNCTIONS (from spatial_statistics.py)
# ============================================================================


async def _analyze_bivariate_moran(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """
    Calculates Bivariate Moran's I to assess spatial correlation between two genes.

    This statistic measures how the expression of one gene in a specific location
    relates to the expression of a second gene in neighboring locations. It is useful
    for identifying pairs of genes that are co-localized or spatially exclusive.
    A positive value suggests that high expression of gene A is surrounded by high
    expression of gene B.
    """
    if context:
        await context.info("Running Bivariate Moran's I analysis...")

    # Get gene pairs from parameters - NO ARBITRARY DEFAULTS
    if not hasattr(params, "gene_pairs") or not params.gene_pairs:
        raise ValueError(
            "Bivariate Moran's I analysis requires explicit gene pairs.\n\n"
            "This tool analyzes spatial correlation between TWO specific genes. "
            "Please provide gene_pairs parameter with scientifically meaningful pairs.\n\n"
            "Example: gene_pairs=[('Gene1', 'Gene2'), ('Gene3', 'Gene4')]\n\n"
            "Note: This tool does NOT create arbitrary gene pairs from HVG list.\n"
            "Only scientifically justified gene pairs should be analyzed."
        )
    else:
        gene_pairs = params.gene_pairs

    results = {}

    try:
        from libpysal.weights import KNN

        coords = adata.obsm["spatial"]
        w = KNN.from_array(coords, k=params.n_neighbors)
        w.transform = "R"

        for gene1, gene2 in gene_pairs:
            if gene1 in adata.var_names and gene2 in adata.var_names:
                x = adata[:, gene1].X
                y = adata[:, gene2].X

                if hasattr(x, "toarray"):
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
                            numerator += (
                                w.sparse[i, j] * (x[i] - x_mean) * (y[j] - y_mean)
                            )

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
        "mean_bivariate_i": float(np.mean(list(results.values()))) if results else 0,
    }


async def _analyze_join_count(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
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

        coords = adata.obsm["spatial"]
        w = KNN.from_array(coords, k=params.n_neighbors)

        # Get categorical data
        y = adata.obs[cluster_key].cat.codes.values

        # Compute join counts
        jc = Join_Counts(y, w)

        return {
            "bb": float(jc.bb),  # Black-Black joins
            "ww": float(jc.ww),  # White-White joins
            "bw": float(jc.bw),  # Black-White joins
            "J": float(jc.J),  # Total joins
            "p_value": float(jc.p_sim) if hasattr(jc, "p_sim") else None,
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
    context: Optional[Context] = None,
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
        if "spatial_connectivities" in adata.obsp:
            conn_matrix = adata.obsp["spatial_connectivities"]
        else:
            # Create connectivity matrix
            coords = adata.obsm["spatial"]
            from sklearn.neighbors import kneighbors_graph

            conn_matrix = kneighbors_graph(
                coords, n_neighbors=params.n_neighbors, mode="connectivity"
            )

        # Convert to networkx graph
        G = nx.from_scipy_sparse_array(conn_matrix)

        # Compute properties
        properties = {
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "density": float(nx.density(G)),
            "is_connected": nx.is_connected(G),
            "n_components": nx.number_connected_components(G),
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
            properties["largest_component_fraction"] = (
                len(largest_cc) / G.number_of_nodes()
            )

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
    context: Optional[Context] = None,
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
        if "spatial_connectivities" in adata.obsp:
            conn_matrix = adata.obsp["spatial_connectivities"]
        else:
            coords = adata.obsm["spatial"]
            from sklearn.neighbors import kneighbors_graph

            conn_matrix = kneighbors_graph(
                coords, n_neighbors=params.n_neighbors, mode="connectivity"
            )

        # Convert to networkx
        G = nx.from_scipy_sparse_array(conn_matrix)

        # Compute centrality measures
        degree_centrality = nx.degree_centrality(G)
        closeness_centrality = nx.closeness_centrality(G)
        betweenness_centrality = nx.betweenness_centrality(G)

        # Store in adata.obs
        adata.obs["degree_centrality"] = pd.Series(degree_centrality)
        adata.obs["closeness_centrality"] = pd.Series(closeness_centrality)
        adata.obs["betweenness_centrality"] = pd.Series(betweenness_centrality)

        # Compute statistics by cluster
        centrality_stats = {}
        for cluster in adata.obs[cluster_key].unique():
            mask = adata.obs[cluster_key] == cluster
            centrality_stats[str(cluster)] = {
                "mean_degree": float(adata.obs.loc[mask, "degree_centrality"].mean()),
                "mean_closeness": float(
                    adata.obs.loc[mask, "closeness_centrality"].mean()
                ),
                "mean_betweenness": float(
                    adata.obs.loc[mask, "betweenness_centrality"].mean()
                ),
            }

        return {
            "centrality_computed": True,
            "cluster_centrality": centrality_stats,
            "global_stats": {
                "mean_degree": float(np.mean(list(degree_centrality.values()))),
                "mean_closeness": float(np.mean(list(closeness_centrality.values()))),
                "mean_betweenness": float(
                    np.mean(list(betweenness_centrality.values()))
                ),
            },
        }

    except ImportError:
        return {"error": "NetworkX required for centrality analysis"}
    except Exception as e:
        return {"error": str(e)}


async def _analyze_local_moran(
    adata: ad.AnnData,
    params: SpatialAnalysisParameters,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """
    Calculate Local Moran's I (LISA) for spatial clustering detection.
    
    Local Moran's I identifies spatial clusters and outliers by measuring
    the local spatial autocorrelation for each observation.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data object
    params : SpatialAnalysisParameters
        Analysis parameters including genes to analyze
    context : Optional[Context]
        MCP context for logging
    
    Returns
    -------
    Dict[str, Any]
        Results including Local Moran's I values and statistics for each gene
    """
    import numpy as np
    from scipy.sparse import issparse
    from scipy.stats import zscore

    try:
        # Ensure spatial neighbors exist
        await _ensure_spatial_neighbors(
            adata, params.n_neighbors, context
        )
        
        # Determine genes to analyze
        if params.genes:
            genes = params.genes
        else:
            # Use top highly variable genes if available
            if "highly_variable" in adata.var.columns:
                hvg_genes = adata.var_names[
                    adata.var["highly_variable"]
                ].tolist()[:5]
                genes = hvg_genes if hvg_genes else adata.var_names[:5].tolist()
            else:
                genes = adata.var_names[:5].tolist()
        
        # Validate genes exist
        valid_genes = []
        for gene in genes:
            if gene in adata.var_names:
                valid_genes.append(gene)
            elif context:
                await context.warning(f"Gene {gene} not found in dataset")
        
        if not valid_genes:
            return {"error": "No valid genes found for analysis"}
        
        # Get spatial weights
        W = adata.obsp["spatial_connectivities"]
        
        results = {}
        for gene in valid_genes:
            # Get expression values
            if issparse(adata.X):
                expr = adata[:, gene].X.toarray().flatten()
            else:
                expr = adata[:, gene].X.flatten()
            
            # Standardize expression
            z_expr = zscore(expr)
            
            # Calculate local Moran's I
            n = len(expr)
            local_i = np.zeros(n)
            
            for i in range(n):
                neighbors = W[i].nonzero()[1]
                if len(neighbors) > 0:
                    local_i[i] = z_expr[i] * np.mean(z_expr[neighbors])
            
            # Store in adata.obs
            adata.obs[f"{gene}_local_morans"] = local_i
            
            # Identify clusters (significant positive values) and outliers (significant negative values)
            threshold = np.percentile(np.abs(local_i), 95)
            hotspots = np.where(local_i > threshold)[0].tolist()
            coldspots = np.where(local_i < -threshold)[0].tolist()
            
            results[gene] = {
                "mean": float(np.mean(local_i)),
                "std": float(np.std(local_i)),
                "min": float(np.min(local_i)),
                "max": float(np.max(local_i)),
                "n_hotspots": len(hotspots),
                "n_coldspots": len(coldspots),
            }
        
        # Store summary in uns
        adata.uns["local_moran"] = {
            "genes_analyzed": valid_genes,
            "n_neighbors": params.n_neighbors,
            "results": results,
        }
        
        if context:
            await context.info(
                f"Calculated Local Moran's I for {len(valid_genes)} genes"
            )
        
        return {
            "analysis_type": "local_moran",
            "genes_analyzed": valid_genes,
            "results": results,
            "interpretation": (
                "Positive values indicate spatial clustering (similar values nearby), "
                "negative values indicate spatial outliers (different values nearby)"
            ),
        }
    
    except Exception as e:
        if context:
            await context.warning(f"Local Moran's I analysis failed: {str(e)}")
        return {"error": f"Local Moran's I analysis failed: {str(e)}"}

