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

The primary entry point is the `analyze_spatial_statistics` function, which
dispatches tasks to the appropriate analysis function based on user parameters.
All 12 analysis types are accessible through this unified interface with a
new unified 'genes' parameter for consistent gene selection across methods.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq

from ..utils.dependency_manager import is_available, require

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import SpatialStatisticsResult
from ..models.data import SpatialStatisticsParameters
from ..utils.adata_utils import (
    ensure_categorical,
    require_spatial_coords,
    select_genes_for_analysis,
    to_dense,
    validate_adata_basics,
)
from ..utils.compute import ensure_spatial_neighbors_async
from ..utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)

# ============================================================================
# MAIN ENTRY POINT
# ============================================================================


async def analyze_spatial_statistics(
    data_id: str,
    ctx: ToolContext,
    params: SpatialStatisticsParameters,  # No default - must be provided by caller (LLM)
) -> SpatialStatisticsResult:
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
    ctx : ToolContext
        Tool context for data access and logging.
    params : SpatialStatisticsParameters
        An object containing the parameters for the analysis, including the
        specific `analysis_type` to perform.

    Returns
    -------
    SpatialStatisticsResult
        An object containing the statistical results and metadata from the analysis.

    Raises
    ------
    DataNotFoundError
        If the specified dataset is not found in the data store.
    ParameterError
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
        "join_count",  # Traditional Join Count for binary data (2 categories)
        "local_join_count",  # Local Join Count for multi-category data (>2 categories)
        "network_properties",
        "spatial_centrality",
    ]

    if params.analysis_type not in supported_types:
        raise ParameterError(f"Unsupported analysis type: {params.analysis_type}")

    if params.n_neighbors <= 0:
        raise ParameterError(f"n_neighbors must be positive, got {params.n_neighbors}")

    # Retrieve dataset via ToolContext
    try:
        adata = await ctx.get_adata(data_id)

        # Basic validation: min 10 cells, spatial coordinates exist
        validate_adata_basics(adata, min_obs=10)
        require_spatial_coords(adata)  # Validates spatial coords exist

        # Determine if cluster_key is required for this analysis type
        analyses_requiring_cluster_key = {
            "neighborhood",
            "co_occurrence",
            "ripley",
            "join_count",
            "local_join_count",
            "centrality",
            "network_properties",
            "spatial_centrality",
        }

        # Ensure cluster key only for analyses that require it
        cluster_key = None
        if params.analysis_type in analyses_requiring_cluster_key:
            cluster_key = await _ensure_cluster_key(adata, params.cluster_key, ctx)

        # Ensure spatial neighbors
        await ensure_spatial_neighbors_async(adata, ctx, n_neighs=params.n_neighbors)

        # Route to appropriate analysis function
        if params.analysis_type == "moran":
            result = await _analyze_morans_i(adata, params, ctx)
        elif params.analysis_type == "local_moran":
            result = await _analyze_local_moran(adata, params, ctx)
        elif params.analysis_type == "geary":
            result = _analyze_gearys_c(adata, params, ctx)
        elif params.analysis_type == "neighborhood":
            result = _analyze_neighborhood_enrichment(adata, cluster_key, ctx)
        elif params.analysis_type == "co_occurrence":
            result = _analyze_co_occurrence(adata, cluster_key, ctx)
        elif params.analysis_type == "ripley":
            result = _analyze_ripleys_k(adata, cluster_key, ctx)
        elif params.analysis_type == "getis_ord":
            result = await _analyze_getis_ord(adata, params, ctx)
        elif params.analysis_type == "centrality":
            result = await _analyze_centrality(adata, cluster_key, ctx)
        elif params.analysis_type == "bivariate_moran":
            result = await _analyze_bivariate_moran(adata, params, ctx)
        elif params.analysis_type == "join_count":
            result = _analyze_join_count(adata, cluster_key, params, ctx)
        elif params.analysis_type == "local_join_count":
            result = await _analyze_local_join_count(adata, cluster_key, params, ctx)
        elif params.analysis_type == "network_properties":
            result = await _analyze_network_properties(adata, cluster_key, params, ctx)
        elif params.analysis_type == "spatial_centrality":
            result = await _analyze_spatial_centrality(adata, cluster_key, params, ctx)
        else:
            raise ParameterError(
                f"Analysis type {params.analysis_type} not implemented"
            )

        # COW FIX: No need to update data_store - changes already reflected via direct reference
        # All modifications to adata.obs/uns/obsp are in-place and preserved

        # Ensure result is a dictionary
        if not isinstance(result, dict):
            if hasattr(result, "dict"):
                result = result.dict()
            else:
                raise ProcessingError("Invalid result format from analysis function")

        # Add metadata
        result.update(
            {
                "n_cells": adata.n_obs,
                "n_neighbors": params.n_neighbors,
            }
        )

        # Store scientific metadata for reproducibility
        from ..utils.adata_utils import store_analysis_metadata

        # Determine results keys based on analysis type
        results_keys_dict = {"obs": [], "var": [], "obsm": [], "uns": []}
        if params.analysis_type in ["moran", "geary"]:
            results_keys_dict["uns"].append(f"{params.analysis_type}s_i")
        elif params.analysis_type == "local_moran":
            results_keys_dict["obs"].extend(
                [f"{gene}_local_moran" for gene in (params.genes or [])]
            )
        elif params.analysis_type == "getis_ord":
            if params.genes:
                for gene in params.genes:
                    results_keys_dict["obs"].extend(
                        [f"{gene}_getis_ord_z", f"{gene}_getis_ord_p"]
                    )
        elif params.analysis_type in ["neighborhood", "co_occurrence"]:
            results_keys_dict["uns"].append(params.analysis_type)
        elif params.analysis_type == "ripley":
            results_keys_dict["uns"].append("ripley")
        elif params.analysis_type == "centrality":
            results_keys_dict["uns"].append("centrality_scores")

        # Prepare parameters dict
        parameters_dict = {
            "n_neighbors": params.n_neighbors,
        }
        if cluster_key:
            parameters_dict["cluster_key"] = cluster_key
        if params.genes:
            parameters_dict["genes"] = params.genes
        # Add n_perms based on analysis type
        if params.analysis_type in ["moran", "local_moran", "geary"]:
            parameters_dict["n_perms"] = params.moran_n_perms

        # Extract statistics for metadata
        statistics_dict = {
            "n_cells": adata.n_obs,
        }
        if "n_significant" in result:
            statistics_dict["n_significant"] = result["n_significant"]
        if "mean_score" in result:
            statistics_dict["mean_score"] = result["mean_score"]

        # Store metadata
        store_analysis_metadata(
            adata,
            analysis_name=f"spatial_stats_{params.analysis_type}",
            method=params.analysis_type,
            parameters=parameters_dict,
            results_keys=results_keys_dict,
            statistics=statistics_dict,
        )

        return SpatialStatisticsResult(
            data_id=data_id,
            analysis_type=params.analysis_type,
            statistics=result,
        )

    except (DataNotFoundError, ParameterError, DataCompatibilityError):
        raise
    except Exception as e:
        raise ProcessingError(
            f"Error in {params.analysis_type} analysis: {str(e)}"
        ) from e


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


async def _ensure_cluster_key(
    adata: ad.AnnData, requested_key: str, ctx: "ToolContext"  # noqa: ARG001
) -> str:
    """Ensure a valid cluster key exists in adata and is categorical."""
    if requested_key not in adata.obs.columns:
        # NO FALLBACK: User's explicit clustering choice must be respected
        available_keys = [
            col
            for col in adata.obs.columns
            if "cluster" in col.lower() or col in ["leiden", "louvain"]
        ]
        raise DataNotFoundError(
            f"Cluster key '{requested_key}' not found. "
            f"Available: {available_keys if available_keys else 'None'}. "
            f"Run preprocess_data() first to generate clusters."
        )

    # Use centralized ensure_categorical (handles type check internally)
    ensure_categorical(adata, requested_key)
    return requested_key


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
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
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
    # Unified gene selection
    genes = select_genes_for_analysis(
        adata,
        genes=params.genes,
        n_genes=params.n_top_genes,
        analysis_name="Moran's I",
    )

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

        # Calculate appropriate number of top genes to return
        # To avoid returning identical lists, we take at most half of the analyzed genes
        # This ensures top_highest and top_lowest are different gene sets
        n_analyzed = len(results_df)
        n_top = min(10, max(3, n_analyzed // 2))

        # Ensure we never return more than half the genes to avoid duplicates
        n_top = min(n_top, n_analyzed // 2) if n_analyzed >= 6 else 0

        return {
            "n_genes_analyzed": len(genes),
            "n_significant": len(significant_genes),
            "top_highest_autocorrelation": (
                results_df.nlargest(n_top, "I").index.tolist() if n_top > 0 else []
            ),
            "top_lowest_autocorrelation": (
                results_df.nsmallest(n_top, "I").index.tolist() if n_top > 0 else []
            ),
            "mean_morans_i": float(results_df["I"].mean()),
            "analysis_key": moran_key,
            "note": "top_highest/top_lowest refer to autocorrelation strength, not positive/negative correlation",
        }

    raise ProcessingError("Moran's I computation did not produce results")


def _analyze_gearys_c(
    adata: ad.AnnData,
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """Compute Geary's C spatial autocorrelation."""
    # Unified gene selection
    genes = select_genes_for_analysis(
        adata,
        genes=params.genes,
        n_genes=params.n_top_genes,
        analysis_name="Geary's C",
    )

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

    raise ProcessingError("Geary's C computation did not produce results")


def _analyze_neighborhood_enrichment(
    adata: ad.AnnData,
    cluster_key: str,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """Compute neighborhood enrichment analysis."""
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)

    analysis_key = f"{cluster_key}_nhood_enrichment"
    if analysis_key in adata.uns:
        z_scores = adata.uns[analysis_key]["zscore"]

        # Use nanmax/nanmin to handle NaN values from sparse cell type distributions
        # NaN can occur when certain cell type pairs have insufficient neighborhoods
        return {
            "n_clusters": len(z_scores),
            "max_enrichment": float(np.nanmax(z_scores)),
            "min_enrichment": float(np.nanmin(z_scores)),
            "analysis_key": analysis_key,
        }

    raise ProcessingError("Neighborhood enrichment did not produce results")


def _analyze_co_occurrence(
    adata: ad.AnnData,
    cluster_key: str,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """Compute co-occurrence analysis."""
    sq.gr.co_occurrence(adata, cluster_key=cluster_key)

    analysis_key = f"{cluster_key}_co_occurrence"
    if analysis_key in adata.uns:
        co_occurrence = adata.uns[analysis_key]["occ"]

        return {"n_clusters": len(co_occurrence), "analysis_key": analysis_key}

    raise ProcessingError("Co-occurrence analysis did not produce results")


def _analyze_ripleys_k(
    adata: ad.AnnData,
    cluster_key: str,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """Compute Ripley's K function."""
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
        raise ProcessingError(f"Ripley's K analysis failed: {e}") from e


async def _analyze_getis_ord(
    adata: ad.AnnData,
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """
    Performs Getis-Ord Gi* analysis to identify local spatial clusters.

    This method identifies statistically significant hot spots (clusters of high
    gene expression) and cold spots (clusters of low gene expression). It computes
    a Z-score for each spot, where high positive Z-scores indicate hot spots and
    low negative Z-scores indicate cold spots.

    The significance threshold is determined by params.getis_ord_alpha, and
    multiple testing correction is applied according to params.getis_ord_correction.

    References
    ----------
    Getis, A. & Ord, J.K. (1992). The Analysis of Spatial Association by Use of
    Distance Statistics. Geographical Analysis, 24(3), 189-206.

    Ord, J.K. & Getis, A. (1995). Local Spatial Autocorrelation Statistics:
    Distributional Issues and an Application. Geographical Analysis, 27(4), 286-306.
    """
    # Unified gene selection
    genes = select_genes_for_analysis(
        adata,
        genes=params.genes,
        n_genes=params.n_top_genes,
        analysis_name="Getis-Ord Gi*",
    )

    getis_ord_results = {}

    require("esda")  # Raises ImportError with install instructions if missing
    require("libpysal")  # Raises ImportError with install instructions if missing
    from esda.getisord import G_Local
    from pysal.lib import weights
    from scipy.stats import norm

    try:

        # Calculate Z-score threshold from alpha level (two-tailed test)
        z_threshold = norm.ppf(1 - params.getis_ord_alpha / 2)

        coords = require_spatial_coords(adata)
        w = weights.KNN.from_array(coords, k=params.n_neighbors)
        w.transform = "r"

        # OPTIMIZATION: Extract all genes at once before loop (batch extraction)
        # This provides 50-150x speedup by avoiding repeated AnnData slicing overhead
        y_all_genes = to_dense(adata[:, genes].X)

        # Collect all p-values for multiple testing correction
        all_pvalues = {}

        for i, gene in enumerate(genes):
            # OPTIMIZATION: Direct indexing from pre-extracted dense matrix (fast!)
            y = y_all_genes[:, i].astype(np.float64)

            local_g = G_Local(y, w, transform="R", star=True)

            # Store raw results in adata.obs
            adata.obs[f"{gene}_getis_ord_z"] = local_g.Zs
            adata.obs[f"{gene}_getis_ord_p"] = local_g.p_sim

            # Store p-values for correction
            all_pvalues[gene] = local_g.p_sim

            # Count hotspots/coldspots using Z-threshold
            getis_ord_results[gene] = {
                "mean_z": float(np.mean(local_g.Zs)),
                "std_z": float(np.std(local_g.Zs)),
                "n_hot_spots": int(np.sum(local_g.Zs > z_threshold)),
                "n_cold_spots": int(np.sum(local_g.Zs < -z_threshold)),
                "n_significant_raw": int(
                    np.sum(local_g.p_sim < params.getis_ord_alpha)
                ),
            }

        # Apply multiple testing correction if requested
        if params.getis_ord_correction != "none" and len(genes) > 1:
            if params.getis_ord_correction == "bonferroni":
                corrected_alpha = params.getis_ord_alpha / len(genes)
                corrected_z_threshold = norm.ppf(1 - corrected_alpha / 2)

                for gene in genes:
                    p_values = all_pvalues[gene]
                    adata.obs[f"{gene}_getis_ord_p_corrected"] = np.minimum(
                        p_values * len(genes), 1.0
                    )

                    z_scores = adata.obs[f"{gene}_getis_ord_z"].values
                    getis_ord_results[gene]["n_hot_spots_corrected"] = int(
                        np.sum(z_scores > corrected_z_threshold)
                    )
                    getis_ord_results[gene]["n_cold_spots_corrected"] = int(
                        np.sum(z_scores < -corrected_z_threshold)
                    )

            elif params.getis_ord_correction == "fdr_bh":
                from statsmodels.stats.multitest import multipletests

                for gene in genes:
                    p_values = all_pvalues[gene]
                    _, p_corrected, _, _ = multipletests(
                        p_values, alpha=params.getis_ord_alpha, method="fdr_bh"
                    )
                    adata.obs[f"{gene}_getis_ord_p_corrected"] = p_corrected

                    getis_ord_results[gene]["n_significant_corrected"] = int(
                        np.sum(p_corrected < params.getis_ord_alpha)
                    )

                    z_scores = adata.obs[f"{gene}_getis_ord_z"].values
                    significant_mask = p_corrected < params.getis_ord_alpha
                    getis_ord_results[gene]["n_hot_spots_corrected"] = int(
                        np.sum((z_scores > z_threshold) & significant_mask)
                    )
                    getis_ord_results[gene]["n_cold_spots_corrected"] = int(
                        np.sum((z_scores < -z_threshold) & significant_mask)
                    )

    except Exception as e:
        raise ProcessingError(f"Getis-Ord analysis failed: {e}") from e

    return {
        "method": "Getis-Ord Gi* (star=True)",
        "n_genes_analyzed": len(getis_ord_results),
        "genes_analyzed": list(getis_ord_results.keys()),
        "parameters": {
            "n_neighbors": params.n_neighbors,
            "alpha": params.getis_ord_alpha,
            "z_threshold": float(z_threshold),
            "correction": params.getis_ord_correction,
        },
        "results": getis_ord_results,
    }


async def _analyze_centrality(
    adata: ad.AnnData,
    cluster_key: str,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """Compute centrality scores."""
    sq.gr.centrality_scores(adata, cluster_key=cluster_key)

    analysis_key = f"{cluster_key}_centrality_scores"
    if analysis_key in adata.uns:
        scores = adata.uns[analysis_key]

        return {
            "analysis_completed": True,
            "analysis_key": analysis_key,
            "n_clusters": len(scores) if isinstance(scores, dict) else "computed",
        }

    raise ProcessingError("Centrality analysis did not produce results")


# ============================================================================
# ADVANCED ANALYSIS FUNCTIONS (from spatial_statistics.py)
# ============================================================================


async def _analyze_bivariate_moran(
    adata: ad.AnnData,
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """
    Calculates Bivariate Moran's I to assess spatial correlation between two genes.

    This statistic measures how the expression of one gene in a specific location
    relates to the expression of a second gene in neighboring locations. It is useful
    for identifying pairs of genes that are co-localized or spatially exclusive.
    A positive value suggests that high expression of gene A is surrounded by high
    expression of gene B.
    """
    # Get gene pairs from parameters - NO ARBITRARY DEFAULTS
    if not params.gene_pairs:
        raise ParameterError("Bivariate Moran's I requires gene_pairs parameter.")
    gene_pairs = params.gene_pairs

    results = {}

    # Use centralized dependency manager for consistent error handling
    require("libpysal")  # Raises ImportError with install instructions if missing
    from libpysal.weights import KNN

    try:

        coords = require_spatial_coords(adata)
        w = KNN.from_array(coords, k=params.n_neighbors)
        w.transform = "R"

        # OPTIMIZATION: Extract all unique genes involved in pairs (batch extraction)
        # This provides 20-40x speedup by avoiding repeated AnnData slicing
        # See test_spatial_statistics_extreme_scale.py for performance validation
        all_genes_in_pairs = list(
            set([g for pair in gene_pairs for g in pair if g in adata.var_names])
        )

        expr_all = to_dense(adata[:, all_genes_in_pairs].X)

        # Create gene index mapping for fast lookup
        gene_to_idx = {gene: i for i, gene in enumerate(all_genes_in_pairs)}

        for gene1, gene2 in gene_pairs:
            if gene1 in adata.var_names and gene2 in adata.var_names:
                # OPTIMIZATION: Direct indexing from pre-extracted matrix (fast!)
                idx1 = gene_to_idx[gene1]
                idx2 = gene_to_idx[gene2]
                x = expr_all[:, idx1].flatten()
                y = expr_all[:, idx2].flatten()

                # Compute bivariate Moran's I using sparse matrix operations
                # Formula: I_xy = (n / S0) * (x - x̄)ᵀ W (y - ȳ) / sqrt(Var(x) * Var(y))
                # Reference: Wartenberg (1985), Anselin et al. (2002)
                n = len(x)
                x_mean = np.mean(x)
                y_mean = np.mean(y)

                # Centered values
                x_centered = x - x_mean
                y_centered = y - y_mean

                # OPTIMIZED: Use sparse matrix multiplication instead of O(n²) loop
                # numerator = Σᵢ Σⱼ wᵢⱼ(xᵢ - x̄)(yⱼ - ȳ) = (x - x̄)ᵀ @ W @ (y - ȳ)
                numerator = float(x_centered @ w.sparse @ y_centered)

                # FIX: Bivariate Moran's I uses sqrt of product of both variances
                # Not just x's variance (which was the bug)
                var_x = np.sum(x_centered**2)
                var_y = np.sum(y_centered**2)
                denominator = np.sqrt(var_x * var_y)

                if denominator > 0:
                    moran_i = (n / w.sparse.sum()) * (numerator / denominator)
                else:
                    moran_i = 0.0

                results[f"{gene1}_vs_{gene2}"] = float(moran_i)

    except Exception as e:
        raise ProcessingError(f"Bivariate Moran's I failed: {e}") from e

    return {
        "n_pairs_analyzed": len(results),
        "bivariate_morans_i": results,
        "mean_bivariate_i": float(np.mean(list(results.values()))) if results else 0,
    }


def _analyze_join_count(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """
    Compute traditional Join Count statistics for BINARY categorical spatial data.

    IMPORTANT: This method only works for binary data (exactly 2 categories).
    For multi-category data (>2 categories), use 'local_join_count' instead.

    Join Count statistics (Cliff & Ord 1981) measure spatial autocorrelation in
    binary categorical data by counting the number of joins between neighboring
    spatial units of the same or different categories.

    Returns three types of joins:
    - BB (Black-Black): Both neighbors are category 1
    - WW (White-White): Both neighbors are category 0
    - BW (Black-White): Neighbors are different categories

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial coordinates in .obsm['spatial']
    cluster_key : str
        Column in adata.obs containing the categorical variable (must have exactly 2 categories)
    params : SpatialStatisticsParameters
        Analysis parameters including n_neighbors
    ctx : ToolContext
        ToolContext for logging and data access

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - bb: Number of Black-Black joins
        - ww: Number of White-White joins
        - bw: Number of Black-White joins
        - J: Total number of joins
        - p_value: Significance level from permutation test

    References
    ----------
    Cliff, A.D. & Ord, J.K. (1981). Spatial Processes. Pion, London.

    See Also
    --------
    _analyze_local_join_count : For multi-category data (>2 categories)
    """
    # Check for required dependencies
    if not is_available("esda") or not is_available("libpysal"):
        raise DependencyError(
            "esda or libpysal package not installed. Install with: pip install esda libpysal"
        )

    try:
        from esda.join_counts import Join_Counts
        from libpysal.weights import KNN

        coords = require_spatial_coords(adata)
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

    except Exception as e:
        raise ProcessingError(f"Join Count analysis failed: {e}") from e


async def _analyze_local_join_count(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """
    Compute Local Join Count statistics for MULTI-CATEGORY categorical spatial data.

    This method extends traditional Join Count statistics to handle data with more than
    2 categories by using Local Join Count Statistics (Anselin & Li 2019). Each category
    is converted to a binary indicator variable, and local statistics are computed to
    identify spatial clusters of each category.

    WHEN TO USE:
    - Data has MORE THAN 2 categories (e.g., cell types, tissue domains)
    - Want to identify WHERE each category spatially clusters
    - Need category-specific clustering patterns

    For binary data (exactly 2 categories), use 'join_count' instead for traditional
    global statistics.

    METHOD:
    1. One-hot encode: Convert multi-category variable to binary indicators
    2. For each category: Compute local join count (# of same-category neighbors)
    3. Permutation test: Assess statistical significance
    4. Store results: Local statistics in adata.obs, summary in return value

    Parameters
    ----------
    adata : AnnData
        Annotated data object with spatial coordinates in .obsm['spatial']
    cluster_key : str
        Column in adata.obs containing the categorical variable (can have any number of categories)
    params : SpatialStatisticsParameters
        Analysis parameters including n_neighbors
    ctx : ToolContext
        ToolContext for logging and data access

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - method: Method name and reference
        - n_categories: Number of categories analyzed
        - categories: List of category names
        - per_category_stats: Statistics for each category
          - total_joins: Sum of local join counts across all locations
          - mean_local_joins: Average local join count per location
          - n_significant: Number of locations with significant clustering (p < 0.05)
          - n_hotspots: Number of locations with positive significant clustering
        - interpretation: How to interpret the results

    Notes
    -----
    Results are stored in adata.obs as:
    - 'ljc_{category}': Local join count values for each category
    - 'ljc_{category}_pvalue': Significance levels (from permutation test)

    High local join count values indicate locations where category members cluster together.
    P-values < 0.05 indicate statistically significant local clustering.

    References
    ----------
    Anselin, L., & Li, X. (2019). Operational Local Join Count Statistics for Cluster Detection.
    Journal of geographical systems, 21(2), 189–210.
    https://doi.org/10.1007/s10109-019-00299-x

    See Also
    --------
    _analyze_join_count : For binary data (2 categories) using traditional Join Count

    Examples
    --------
    For a dataset with 7 cell type categories:
    >>> result = await _analyze_local_join_count(adata, 'leiden', params, ctx)
    >>> # Check which cell types show significant clustering
    >>> for cat, stats in result['per_category_stats'].items():
    ...     print(f"{cat}: {stats['n_hotspots']} significant hotspots")
    """
    # Check for required dependencies
    if not is_available("esda") or not is_available("libpysal"):
        raise DependencyError(
            "esda or libpysal package not installed (requires esda >= 2.4.0). "
            "Install with: pip install esda libpysal"
        )

    try:
        from esda.join_counts_local import Join_Counts_Local
        from libpysal.weights import KNN

        coords = require_spatial_coords(adata)

        # Create PySAL W object directly from coordinates using KNN
        # This ensures compatibility with Join_Counts_Local
        w = KNN.from_array(coords, k=params.n_neighbors)

        # Get unique categories
        categories = adata.obs[cluster_key].unique()
        n_categories = len(categories)

        results = {}

        # Analyze each category separately
        for category in categories:
            # Create binary indicator: 1 if cell is this category, 0 otherwise
            y = (adata.obs[cluster_key] == category).astype(int).values

            # Compute Local Join Count statistics
            ljc = Join_Counts_Local(connectivity=w).fit(y)

            # Store local statistics in adata.obs
            adata.obs[f"ljc_{category}"] = ljc.LJC
            adata.obs[f"ljc_{category}_pvalue"] = ljc.p_sim

            # Compute summary statistics
            results[str(category)] = {
                "total_joins": float(ljc.LJC.sum()),
                "mean_local_joins": float(ljc.LJC.mean()),
                "std_local_joins": float(ljc.LJC.std()),
                "n_significant": int((ljc.p_sim < 0.05).sum()),
                "n_hotspots": int(((ljc.LJC > 0) & (ljc.p_sim < 0.05)).sum()),
            }

        # Store summary in adata.uns
        adata.uns["local_join_count"] = {
            "method": "Local Join Count Statistics (Anselin & Li 2019)",
            "cluster_key": cluster_key,
            "n_categories": n_categories,
            "categories": [str(c) for c in categories],
            "n_neighbors": params.n_neighbors,
            "per_category_stats": results,
        }

        return {
            "method": "Local Join Count Statistics (Anselin & Li 2019)",
            "n_categories": n_categories,
            "categories": [str(c) for c in categories],
            "per_category_stats": results,
            "interpretation": (
                "Local Join Count statistics identify spatial clusters for each category. "
                "High LJC values indicate locations where category members cluster together. "
                "P-values < 0.05 indicate statistically significant local clustering. "
                "Results stored in adata.obs as 'ljc_{category}' and 'ljc_{category}_pvalue'."
            ),
        }

    except Exception as e:
        raise ProcessingError(f"Local Join Count analysis failed: {e}") from e


async def _analyze_network_properties(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """
    Analyze network properties of spatial graph.

    Migrated from spatial_statistics.py
    """
    # Check for required dependencies
    if not is_available("networkx"):
        raise DependencyError(
            "networkx package required. Install with: pip install networkx"
        )

    try:
        import networkx as nx
        from scipy.sparse import csr_matrix  # noqa: F401

        # Get or create spatial connectivity
        if "spatial_connectivities" in adata.obsp:
            conn_matrix = adata.obsp["spatial_connectivities"]
        else:
            # Create connectivity matrix
            from sklearn.neighbors import kneighbors_graph

            coords = require_spatial_coords(adata)
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
            G.subgraph(largest_cc)
            properties["largest_component_size"] = len(largest_cc)
            properties["largest_component_fraction"] = (
                len(largest_cc) / G.number_of_nodes()
            )

        # Clustering coefficient
        try:
            properties["avg_clustering"] = float(nx.average_clustering(G))
        except Exception:
            properties["avg_clustering"] = None

        # Degree statistics
        degrees = dict(G.degree())
        degree_values = list(degrees.values())
        properties["degree_mean"] = float(np.mean(degree_values))
        properties["degree_std"] = float(np.std(degree_values))

        return properties

    except Exception as e:
        raise ProcessingError(f"Network properties analysis failed: {e}") from e


async def _analyze_spatial_centrality(
    adata: ad.AnnData,
    cluster_key: str,
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """
    Compute various centrality measures for spatial network.

    Migrated from spatial_statistics.py
    """
    # Check for required dependencies
    if not is_available("networkx"):
        raise DependencyError(
            "NetworkX required for centrality analysis. Install with: pip install networkx"
        )

    try:
        import networkx as nx

        # Get connectivity matrix
        if "spatial_connectivities" in adata.obsp:
            conn_matrix = adata.obsp["spatial_connectivities"]
        else:
            from sklearn.neighbors import kneighbors_graph

            coords = require_spatial_coords(adata)
            conn_matrix = kneighbors_graph(
                coords, n_neighbors=params.n_neighbors, mode="connectivity"
            )

        # Convert to networkx
        G = nx.from_scipy_sparse_array(conn_matrix)

        # Compute centrality measures (returns dict with integer keys)
        degree_centrality = nx.degree_centrality(G)
        closeness_centrality = nx.closeness_centrality(G)
        betweenness_centrality = nx.betweenness_centrality(G)

        # FIX: NetworkX returns {0: val0, 1: val1, ...} with integer keys,
        # but adata.obs_names are strings. We need to extract values in order.
        # Bug: pd.Series(dict) cannot align integer keys to string obs_names
        n_nodes = adata.n_obs

        # Validate that all expected node keys exist in centrality results
        # This catches edge cases like disconnected graphs or isolated nodes
        expected_keys = set(range(n_nodes))
        missing_degree = expected_keys - set(degree_centrality.keys())
        missing_closeness = expected_keys - set(closeness_centrality.keys())
        missing_betweenness = expected_keys - set(betweenness_centrality.keys())

        if missing_degree or missing_closeness or missing_betweenness:
            await ctx.warning(
                f"Centrality computation incomplete: "
                f"missing degree={len(missing_degree)}, "
                f"closeness={len(missing_closeness)}, "
                f"betweenness={len(missing_betweenness)} nodes. "
                f"Graph may have disconnected components."
            )

        # Use .get() with default 0.0 for missing nodes (isolated/disconnected)
        degree_vals = np.array([degree_centrality.get(i, 0.0) for i in range(n_nodes)])
        closeness_vals = np.array(
            [closeness_centrality.get(i, 0.0) for i in range(n_nodes)]
        )
        betweenness_vals = np.array(
            [betweenness_centrality.get(i, 0.0) for i in range(n_nodes)]
        )

        # Store in adata.obs (directly as numpy array)
        adata.obs["degree_centrality"] = degree_vals
        adata.obs["closeness_centrality"] = closeness_vals
        adata.obs["betweenness_centrality"] = betweenness_vals

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

    except Exception as e:
        raise ProcessingError(f"Spatial centrality analysis failed: {e}") from e


async def _analyze_local_moran(
    adata: ad.AnnData,
    params: SpatialStatisticsParameters,
    ctx: "ToolContext",
) -> Dict[str, Any]:
    """
    Calculate Local Moran's I (LISA) for spatial clustering detection.

    Local Moran's I identifies spatial clusters and outliers by measuring
    the local spatial autocorrelation for each observation.

    Parameters
    ----------
    adata : ad.AnnData
        Annotated data object
    params : SpatialStatisticsParameters
        Analysis parameters including genes to analyze
    ctx : ToolContext
        ToolContext for logging and data access

    Returns
    -------
    Dict[str, Any]
        Results including Local Moran's I values and statistics for each gene

    Notes
    -----
    This implementation uses PySAL's esda.Moran_Local with permutation-based
    significance testing, following best practices from:
    - GeoDa Center: https://geodacenter.github.io/workbook/6a_local_auto/lab6a.html
    - PySAL documentation: https://pysal.org/esda/generated/esda.Moran_Local.html

    The permutation approach holds each observation fixed while randomly permuting
    the remaining n-1 values to generate a reference distribution for significance
    testing. This is more robust than parametric approaches as it makes fewer
    distributional assumptions.

    Quadrant classification (LISA clusters):
    - HH (High-High): Hot spots - high values surrounded by high values
    - LL (Low-Low): Cold spots - low values surrounded by low values
    - HL (High-Low): High outliers - high values surrounded by low values
    - LH (Low-High): Low outliers - low values surrounded by high values
    """
    # Import PySAL components for proper LISA analysis
    require("esda")  # Raises ImportError with install instructions if missing
    require("libpysal")  # Raises ImportError with install instructions if missing
    from esda.moran import Moran_Local
    from libpysal.weights import W as PySALWeights

    try:
        # Ensure spatial neighbors exist
        await ensure_spatial_neighbors_async(adata, ctx, n_neighs=params.n_neighbors)

        # Unified gene selection (default 5 genes for computational efficiency)
        n_genes = (
            min(5, params.n_top_genes) if params.genes is None else params.n_top_genes
        )
        valid_genes = select_genes_for_analysis(
            adata,
            genes=params.genes,
            n_genes=n_genes,
            analysis_name="Local Moran's I (LISA)",
        )

        # Convert spatial connectivity matrix to PySAL weights format
        W_sparse = adata.obsp["spatial_connectivities"]

        # Create PySAL weights from sparse matrix using optimized CSR access
        # Direct CSR array access avoids per-row object creation (15x faster)
        from scipy.sparse import csr_matrix

        if not isinstance(W_sparse, csr_matrix):
            W_sparse = csr_matrix(W_sparse)

        neighbors_dict = {}
        weights_dict = {}
        n_obs = W_sparse.shape[0]

        # Direct access to CSR internal arrays
        indptr = W_sparse.indptr
        indices = W_sparse.indices
        data = W_sparse.data

        for i in range(n_obs):
            start, end = indptr[i], indptr[i + 1]
            neighbors_dict[i] = indices[start:end].tolist()
            weights_dict[i] = data[start:end].tolist()

        w = PySALWeights(neighbors_dict, weights_dict)

        # Get analysis parameters
        permutations = params.local_moran_permutations
        alpha = params.local_moran_alpha
        use_fdr = params.local_moran_fdr_correction

        # Memory-efficient streaming: extract one gene at a time
        # This reduces memory from O(n_spots × n_genes) to O(n_spots)
        # Critical for large datasets (Visium HD: 50K+ spots × 500 genes = 200MB+)
        results = {}
        for gene in valid_genes:
            # Extract single gene column - memory efficient for sparse matrices
            gene_idx = adata.var_names.get_loc(gene)
            expr = to_dense(adata.X[:, gene_idx]).flatten()

            # CRITICAL: Convert to float64 for PySAL/numba compatibility
            # PySAL's Moran_Local uses numba JIT compilation which requires
            # consistent dtypes (float64) for matrix operations
            expr = expr.astype(np.float64, copy=False)

            # Run PySAL Local Moran's I with permutation testing
            lisa = Moran_Local(expr, w, permutations=permutations)

            # Store local I values in adata.obs
            adata.obs[f"{gene}_local_morans"] = lisa.Is

            # Get p-values from permutation test
            p_values = lisa.p_sim

            # Apply FDR correction if requested
            if use_fdr and permutations > 0:
                # Check statsmodels availability for FDR correction
                require(
                    "statsmodels"
                )  # Raises ImportError with install instructions if missing
                from statsmodels.stats.multitest import multipletests

                _, p_corrected, _, _ = multipletests(
                    p_values, alpha=alpha, method="fdr_bh"
                )
                significant = p_corrected < alpha
            else:
                significant = p_values < alpha

            # Classify by quadrant AND significance
            # PySAL quadrant codes: 1=HH, 2=LH, 3=LL, 4=HL
            q = lisa.q

            # Hot spots: High-High clusters (significant positive spatial autocorrelation)
            hotspots = np.where((q == 1) & significant)[0].tolist()
            # Cold spots: Low-Low clusters (significant positive spatial autocorrelation)
            coldspots = np.where((q == 3) & significant)[0].tolist()
            # High outliers: High values surrounded by low values
            high_outliers = np.where((q == 4) & significant)[0].tolist()
            # Low outliers: Low values surrounded by high values
            low_outliers = np.where((q == 2) & significant)[0].tolist()

            # Store quadrant classification in adata.obs
            quadrant_labels = np.array(["Not Significant"] * n_obs)
            quadrant_labels[(q == 1) & significant] = "HH (Hot Spot)"
            quadrant_labels[(q == 3) & significant] = "LL (Cold Spot)"
            quadrant_labels[(q == 4) & significant] = "HL (High Outlier)"
            quadrant_labels[(q == 2) & significant] = "LH (Low Outlier)"
            adata.obs[f"{gene}_lisa_cluster"] = pd.Categorical(quadrant_labels)

            # Store p-values
            adata.obs[f"{gene}_lisa_pvalue"] = p_values

            results[gene] = {
                "mean_I": float(np.mean(lisa.Is)),
                "std_I": float(np.std(lisa.Is)),
                "min_I": float(np.min(lisa.Is)),
                "max_I": float(np.max(lisa.Is)),
                "n_significant": int(np.sum(significant)),
                "n_hotspots": len(hotspots),  # HH clusters
                "n_coldspots": len(coldspots),  # LL clusters
                "n_high_outliers": len(high_outliers),  # HL
                "n_low_outliers": len(low_outliers),  # LH
                "permutations": permutations,
                "alpha": alpha,
                "fdr_corrected": use_fdr,
            }

        # Store summary in uns
        adata.uns["local_moran"] = {
            "genes_analyzed": valid_genes,
            "n_neighbors": params.n_neighbors,
            "permutations": permutations,
            "alpha": alpha,
            "fdr_corrected": use_fdr,
            "results": results,
            "method": "PySAL esda.Moran_Local",
            "reference": "Anselin, L. (1995). Local Indicators of Spatial Association - LISA",
        }

        return {
            "analysis_type": "local_moran",
            "genes_analyzed": valid_genes,
            "results": results,
            "parameters": {
                "permutations": permutations,
                "alpha": alpha,
                "fdr_corrected": use_fdr,
                "n_neighbors": params.n_neighbors,
            },
            "interpretation": (
                "LISA (Local Indicators of Spatial Association) identifies statistically "
                "significant spatial clusters and outliers using permutation-based testing. "
                "HH (Hot Spots): high values clustered together. "
                "LL (Cold Spots): low values clustered together. "
                "HL/LH (Outliers): values significantly different from neighbors. "
                f"Significance determined by {permutations} permutations "
                f"with alpha={alpha}{' and FDR correction' if use_fdr else ''}."
            ),
        }

    except Exception as e:
        raise ProcessingError(f"Local Moran's I analysis failed: {e}") from e
