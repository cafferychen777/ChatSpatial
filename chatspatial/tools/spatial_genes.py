"""
Spatial Variable Genes (SVG) identification for ChatSpatial MCP.

This module provides implementations for SVG detection methods including SpatialDE and SPARK-X,
enabling comprehensive spatial transcriptomics analysis. Each method offers distinct advantages
for identifying genes with spatial expression patterns.

Methods Overview:
    - SPARK-X (default): Non-parametric statistical method, best accuracy, requires R
    - SpatialDE: Gaussian process-based kernel method, statistically rigorous

The module integrates these tools into the ChatSpatial MCP framework, handling data preparation,
execution, result formatting, and error management across different computational backends.
"""

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from collections import Counter

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..models.analysis import SpatialVariableGenesResult  # noqa: E402
from ..models.data import SpatialVariableGenesParameters  # noqa: E402
from ..utils import validate_var_column  # noqa: E402
from ..utils.adata_utils import (  # noqa: E402
    get_raw_data_source,
    require_spatial_coords,
    to_dense,
)
from ..utils.dependency_manager import require  # noqa: E402
from ..utils.exceptions import DataNotFoundError  # noqa: E402
from ..utils.exceptions import DataError, ParameterError, ProcessingError
from ..utils.mcp_utils import suppress_output  # noqa: E402

# =============================================================================
# Shared Utilities for Spatial Variable Gene Detection
# =============================================================================

# Default limit for spatial_genes list returned to LLM
# Full results stored in adata.var for complete access
DEFAULT_TOP_GENES_LIMIT = 500


def _ensure_unique_gene_names(gene_names: list[str]) -> list[str]:
    """Ensure gene names are unique by adding suffixes to duplicates.

    Required for R-based methods (SPARK-X) that use gene names as rownames.

    Args:
        gene_names: List of gene names (may contain duplicates)

    Returns:
        List of unique gene names with suffixes added to duplicates
    """
    if len(gene_names) == len(set(gene_names)):
        return gene_names

    gene_counts = Counter(gene_names)
    unique_names = []
    seen_counts: dict[str, int] = {}

    for gene in gene_names:
        if gene_counts[gene] > 1:
            if gene not in seen_counts:
                seen_counts[gene] = 0
                unique_names.append(gene)
            else:
                seen_counts[gene] += 1
                unique_names.append(f"{gene}_{seen_counts[gene]}")
        else:
            unique_names.append(gene)

    return unique_names


def _calculate_sparse_gene_stats(X) -> tuple[np.ndarray, np.ndarray]:
    """Calculate gene statistics on sparse or dense matrix.

    Efficiently computes gene totals and expression counts without densifying
    the entire matrix.

    Args:
        X: Gene expression matrix (cells × genes), sparse or dense

    Returns:
        Tuple of (gene_totals, n_expressed_per_gene) as 1D arrays
    """
    is_sparse = sp.issparse(X)

    if is_sparse:
        gene_totals = np.array(X.sum(axis=0)).flatten()
        n_expressed = np.array((X > 0).sum(axis=0)).flatten()
    else:
        gene_totals = np.asarray(X.sum(axis=0)).flatten()
        n_expressed = np.asarray((X > 0).sum(axis=0)).flatten()

    return gene_totals, n_expressed


async def identify_spatial_genes(
    data_id: str,
    ctx: "ToolContext",
    params: SpatialVariableGenesParameters,
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using statistical methods.

    This is the main entry point for spatial gene detection, routing to the appropriate
    method based on params.method. Each method has different strengths:

    Method Selection Guide:
        - SPARK-X (default): Best for accuracy, handles large datasets efficiently
        - SpatialDE: Best for statistical rigor in publication-ready analyses

    Data Requirements:
        - SPARK-X: Works with raw counts or normalized data
        - SpatialDE: Works with raw count data

    Args:
        data_id: Dataset identifier in data store
        ctx: ToolContext for data access and logging
        params: Method-specific parameters (see SpatialVariableGenesParameters)

    Returns:
        SpatialVariableGenesResult containing:
            - List of significant spatial genes
            - Statistical metrics (p-values, q-values)
            - Method-specific results

    Raises:
        ValueError: If dataset not found or spatial coordinates missing
        ImportError: If required method dependencies not installed

    Performance Notes:
        - SPARK-X: ~2-5 min for 3000 spots × 20000 genes
        - SpatialDE: ~15-30 min (scales with spot count squared)
    """
    # Get data via ToolContext
    adata = await ctx.get_adata(data_id)

    # Validate spatial coordinates exist
    require_spatial_coords(adata, spatial_key=params.spatial_key)

    # Route to appropriate method
    if params.method == "spatialde":
        return await _identify_spatial_genes_spatialde(data_id, adata, params, ctx)
    elif params.method == "sparkx":
        return await _identify_spatial_genes_sparkx(data_id, adata, params, ctx)
    else:
        raise ParameterError(
            f"Unsupported method: {params.method}. Available methods: spatialde, sparkx"
        )


async def _identify_spatial_genes_spatialde(
    data_id: str,
    adata: Any,
    params: SpatialVariableGenesParameters,
    ctx: "ToolContext",
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using the SpatialDE statistical framework.

    SpatialDE employs Gaussian process regression with spatial kernels to decompose
    gene expression variance into spatial and non-spatial components. It provides
    rigorous statistical testing for spatial expression patterns with multiple
    testing correction.

    Official Preprocessing Workflow (Implemented):
        This implementation follows the official SpatialDE best practices:
        1. Filter low-expression genes (total_counts >= 3)
        2. Variance stabilization (NaiveDE.stabilize)
        3. Regress out library size effects (NaiveDE.regress_out)
        4. Run SpatialDE spatial covariance test
        5. Apply FDR correction (Storey q-value)

    Method Details:
        - Models spatial correlation using squared exponential kernel
        - Tests significance via likelihood ratio test
        - Applies FDR correction for multiple testing
        - Returns both raw and adjusted p-values

    Key Parameters:
        - n_top_genes: Limit analysis to top N genes (for performance)
            * If provided, preferentially uses HVGs if available
            * Recommended: 1000-3000 for quick analysis
            * None (default): Test all genes (may take 15-30 min for large datasets)

    Performance Notes:
        - ~10 minutes for 14,000 genes (official benchmark)
        - Scales approximately linearly with gene count
        - Performance warning issued when n_genes > 5000
        - Tip: Use n_top_genes parameter to reduce runtime

    Data Requirements:
        - Raw count data (from adata.raw or adata.X)
        - 2D spatial coordinates in adata.obsm['spatial']
        - Data will be automatically preprocessed using official workflow

    Returns:
        Results including:
            - List of significant spatial genes (q-value < 0.05)
            - Log-likelihood ratios as test statistics
            - Raw p-values and FDR-corrected q-values
            - Spatial correlation length scale per gene

    Requirements:
        - SpatialDE package with NaiveDE module
        - 2D spatial coordinates
        - Raw count data (not normalized)

    References:
        Svensson et al. (2018) "SpatialDE: identification of spatially variable genes"
        Nature Methods, DOI: 10.1038/nmeth.4636
        Official tutorial: https://github.com/Teichlab/SpatialDE
    """
    # Use centralized dependency manager for consistent error handling
    require("spatialde")  # Raises ImportError with install instructions if missing

    # Apply scipy compatibility patch for SpatialDE (scipy >= 1.14 removed scipy.misc.derivative)
    from ..utils.compat import ensure_spatialde_compat

    ensure_spatialde_compat()

    import NaiveDE
    import SpatialDE
    from SpatialDE.util import qvalue

    # Prepare spatial coordinates
    coords = pd.DataFrame(
        adata.obsm[params.spatial_key][:, :2],  # Ensure 2D coordinates
        columns=["x", "y"],
        index=adata.obs_names,
    )

    # Get raw count data for SpatialDE preprocessing
    # Use get_raw_data_source (single source of truth) for complete gene coverage
    # OPTIMIZATION: Filter genes on SPARSE matrix first, then convert only selected genes to dense
    raw_result = get_raw_data_source(adata, prefer_complete_genes=True)
    raw_data = raw_result.X
    var_names = raw_result.var_names
    var_df = adata.var  # For HVG lookup (HVG info is in adata.var)

    # Step 1: Filter low-expression genes ON SPARSE MATRIX (Official recommendation)
    # SpatialDE README: "Filter practically unobserved genes" with total_counts >= 3
    gene_totals, _ = _calculate_sparse_gene_stats(raw_data)

    keep_genes_mask = gene_totals >= 3
    selected_var_names = var_names[keep_genes_mask]
    # Step 2: Select top N HVGs ON SPARSE MATRIX (if requested)
    # This further reduces genes BEFORE densification
    final_genes = selected_var_names

    if params.n_top_genes is not None and params.n_top_genes < len(selected_var_names):
        if "highly_variable" in var_df.columns:
            # Prioritize HVGs if available
            hvg_mask = var_df.loc[selected_var_names, "highly_variable"]
            hvg_genes = selected_var_names[hvg_mask]

            if len(hvg_genes) >= params.n_top_genes:
                # Use HVGs
                final_genes = hvg_genes[: params.n_top_genes]
            else:
                # Not enough HVGs, select by expression
                gene_totals_filtered = gene_totals[keep_genes_mask]
                top_indices = np.argsort(gene_totals_filtered)[-params.n_top_genes :][
                    ::-1
                ]
                final_genes = selected_var_names[top_indices]
        else:
            # Select by expression
            gene_totals_filtered = gene_totals[keep_genes_mask]
            top_indices = np.argsort(gene_totals_filtered)[-params.n_top_genes :][::-1]
            final_genes = selected_var_names[top_indices]

    # Step 3: Slice sparse matrix to final genes, THEN convert to dense
    # This is where the memory optimization happens: only convert selected genes
    # Directly use raw_result (single source of truth) - no need for double access
    gene_indices = [raw_result.var_names.get_loc(g) for g in final_genes]

    # Now create DataFrame from the SUBSET (much smaller memory footprint)
    counts = pd.DataFrame(
        to_dense(raw_result.X[:, gene_indices]),
        columns=final_genes,
        index=adata.obs_names,
    )

    # Performance warning for large gene sets
    n_genes = counts.shape[1]
    n_spots = counts.shape[0]
    if n_genes > 5000:
        estimated_time = int(n_genes / 14000 * 10)  # Based on 14k genes = 10 min
        await ctx.warning(
            f"WARNING:Running SpatialDE on {n_genes} genes × {n_spots} spots may take {estimated_time}-{estimated_time*2} minutes.\n"
            f"   • Official benchmark: ~10 min for 14,000 genes\n"
            f"   • Tip: Use n_top_genes=1000-3000 to test fewer genes\n"
            f"   • Or use method='sparkx' for faster analysis (2-5 min)"
        )

    # Calculate total counts per spot for regress_out
    total_counts = pd.DataFrame(
        {"total_counts": counts.sum(axis=1)}, index=counts.index
    )

    # Apply official SpatialDE preprocessing workflow
    # Step 1: Variance stabilization
    norm_expr = NaiveDE.stabilize(counts.T).T

    # Step 2: Regress out library size effects
    resid_expr = NaiveDE.regress_out(
        total_counts, norm_expr.T, "np.log(total_counts)"
    ).T

    # Step 3: Run SpatialDE
    results = SpatialDE.run(coords.values, resid_expr)

    # Multiple testing correction using Storey q-value method
    if params.spatialde_pi0 is not None:
        # User-specified pi0 value
        results["qval"] = qvalue(results["pval"].values, pi0=params.spatialde_pi0)
    else:
        # Adaptive pi0 estimation (SpatialDE default, recommended)
        results["qval"] = qvalue(results["pval"].values)

    # Sort by q-value
    results = results.sort_values("qval")

    # Filter significant genes
    significant_genes_all = results[results["qval"] < 0.05]["g"].tolist()

    # Limit for MCP response (full results stored in adata.var)
    limit = params.n_top_genes or DEFAULT_TOP_GENES_LIMIT
    significant_genes = significant_genes_all[:limit]

    # Store results in adata
    results_key = f"spatialde_results_{data_id}"
    adata.var["spatialde_pval"] = results.set_index("g")["pval"]
    adata.var["spatialde_qval"] = results.set_index("g")["qval"]
    adata.var["spatialde_l"] = results.set_index("g")["l"]

    # Store scientific metadata for reproducibility
    from ..utils.adata_utils import store_analysis_metadata
    from ..utils.results_export import export_analysis_result

    store_analysis_metadata(
        adata,
        analysis_name="spatial_genes_spatialde",
        method="spatialde_official_workflow",
        parameters={
            "kernel": params.spatialde_kernel,
            "preprocessing": "NaiveDE.stabilize + NaiveDE.regress_out",
            "gene_filter_threshold": 3,
            "n_genes_tested": n_genes,
            "n_spots": n_spots,
            "pi0": (
                params.spatialde_pi0 if params.spatialde_pi0 is not None else "adaptive"
            ),
        },
        results_keys={
            "var": ["spatialde_pval", "spatialde_qval", "spatialde_l"],
            "obs": [],
            "obsm": [],
            "uns": [],
        },
        statistics={
            "n_genes_analyzed": len(results),
            "n_significant_genes": len(
                results[results["qval"] < 0.05]  # FDR standard threshold
            ),
        },
    )

    # Export results to CSV for reproducibility
    export_analysis_result(adata, data_id, "spatial_genes_spatialde")

    # Note: Detailed statistics (gene_statistics, p_values, q_values) are excluded
    # from MCP response via Field(exclude=True) in SpatialVariableGenesResult.
    # Full results are accessible via adata.var['spatialde_pval', 'spatialde_qval'].

    result = SpatialVariableGenesResult(
        data_id=data_id,
        method="spatialde",
        n_genes_analyzed=len(results),
        n_significant_genes=len(significant_genes_all),
        spatial_genes=significant_genes,
        results_key=results_key,
    )

    return result


async def _identify_spatial_genes_sparkx(
    data_id: str,
    adata: Any,
    params: SpatialVariableGenesParameters,
    ctx: "ToolContext",
) -> SpatialVariableGenesResult:
    """
    Identify spatial variable genes using the SPARK-X non-parametric method.

    SPARK-X is an efficient non-parametric method for detecting spatially variable
    genes without assuming specific distribution models. It uses spatial covariance
    testing and is particularly effective for large-scale datasets. The method is
    implemented in R and accessed via rpy2.

    Method Advantages:
        - Non-parametric: No distributional assumptions required
        - Computationally efficient: Scales well with gene count
        - Robust: Handles various spatial patterns effectively
        - Flexible: Works with both single and mixture spatial kernels

    Gene Filtering Pipeline (based on SPARK-X paper + 2024 best practices):
        TIER 1 - Standard Filtering (SPARK-X paper):
            - filter_mt_genes: Remove mitochondrial genes (MT-*, mt-*) [default: True]
            - filter_ribo_genes: Remove ribosomal genes (RPS*, RPL*) [default: False]
            - Expression filtering: Min percentage + total counts

        TIER 2 - Advanced Options (2024 best practice from PMC11537352):
            - test_only_hvg: Test only highly variable genes [default: False]
              * Reduces housekeeping gene dominance
              * Requires prior HVG computation in preprocessing

        TIER 3 - Quality Warnings:
            - warn_housekeeping: Warn if >30% top genes are housekeeping [default: True]
              * Alerts about potential biological interpretation issues

    Key Parameters:
        - sparkx_option: 'single' or 'mixture' kernel (default: 'mixture')
        - sparkx_percentage: Min percentage of cells expressing gene (default: 0.1)
        - sparkx_min_total_counts: Min total counts per gene (default: 10)
        - sparkx_n_cores: Number of CPU cores for parallel processing
        - filter_mt_genes: Filter mitochondrial genes (default: True)
        - filter_ribo_genes: Filter ribosomal genes (default: False)
        - test_only_hvg: Test only HVGs (default: False)
        - warn_housekeeping: Warn about housekeeping dominance (default: True)

    Data Processing:
        - Automatically filters low-expression genes based on parameters
        - Uses raw counts when available (adata.raw), otherwise current matrix
        - Handles duplicate gene names by adding suffixes

    Returns:
        Results including:
            - List of significant spatial genes (adjusted p-value < 0.05)
            - Raw p-values from spatial covariance test
            - Bonferroni-adjusted p-values
            - Results dataframe with all tested genes
            - Quality warnings if housekeeping genes dominate

    Requirements:
        - R installation with SPARK package
        - rpy2 Python package for R integration
        - Raw count data preferred (will use adata.raw if available)

    Performance:
        - Fastest among the three methods
        - ~2-5 minutes for typical datasets (3000 spots × 20000 genes)
        - Memory efficient through gene filtering

    References:
        - SPARK-X paper: Sun et al. (2021) Genome Biology
        - HVG+SVG best practice: PMC11537352 (2024)
    """
    # Use centralized dependency manager for consistent error handling
    require("rpy2")  # Raises ImportError with install instructions if missing
    from rpy2 import robjects as ro
    from rpy2.rinterface_lib import openrlib  # For thread safety
    from rpy2.robjects import conversion, default_converter
    from rpy2.robjects.packages import importr

    # Prepare spatial coordinates - SPARK needs data.frame format
    coords_array = adata.obsm[params.spatial_key][:, :2].astype(float)
    n_spots, n_genes = adata.shape

    # ==================== OPTIMIZED: Filter on sparse matrix, then convert ====================
    # Strategy: Keep data sparse throughout filtering, only convert final filtered result
    # Benefit: For 30k cells × 20k genes → 3k genes: save ~15GB memory

    # Get sparse count matrix using get_raw_data_source (single source of truth)
    raw_result = get_raw_data_source(adata, prefer_complete_genes=True)
    sparse_counts = raw_result.X  # Keep sparse!
    gene_names = [str(name) for name in raw_result.var_names]
    n_genes = len(gene_names)

    # Ensure gene names are unique (required for SPARK-X R rownames)
    gene_names = _ensure_unique_gene_names(gene_names)

    # ==================== Gene Filtering Pipeline (ON SPARSE MATRIX) ====================
    # Following SPARK-X paper best practices + 2024 literature recommendations
    # All filtering done on sparse matrix to minimize memory usage

    # Initialize gene mask (all True = keep all genes initially)
    gene_mask = np.ones(len(gene_names), dtype=bool)

    # TIER 1: Mitochondrial gene filtering (SPARK-X paper standard practice)
    # Use pattern-based detection on gene names (works regardless of data source)
    if params.filter_mt_genes:
        mt_mask = np.array([gene.startswith(("MT-", "mt-")) for gene in gene_names])
        n_mt_genes = mt_mask.sum()
        if n_mt_genes > 0:
            gene_mask &= ~mt_mask  # Exclude MT genes

    # TIER 1: Ribosomal gene filtering (optional)
    # Use pattern-based detection on gene names
    if params.filter_ribo_genes:
        ribo_mask = np.array(
            [gene.startswith(("RPS", "RPL", "Rps", "Rpl")) for gene in gene_names]
        )

        n_ribo_genes = ribo_mask.sum()
        if n_ribo_genes > 0:
            gene_mask &= ~ribo_mask  # Exclude ribosomal genes

    # TIER 2: HVG-only testing (2024 best practice from PMC11537352)
    if params.test_only_hvg:
        # Check if HVGs are available in adata.var (the preprocessed data)
        validate_var_column(
            adata,
            "highly_variable",
            "Highly variable genes marker (test_only_hvg=True requires this)",
        )

        # Get HVG list from preprocessed data (adata.var)
        hvg_genes_set = set(adata.var_names[adata.var["highly_variable"]])

        if len(hvg_genes_set) == 0:
            raise DataNotFoundError("No HVGs found. Run preprocessing first.")

        # Filter gene_names to only include HVGs
        hvg_mask = np.array([gene in hvg_genes_set for gene in gene_names])
        n_hvg = hvg_mask.sum()

        if n_hvg == 0:
            # No overlap between current gene list and HVGs
            raise DataError(
                f"test_only_hvg=True but no overlap found between current gene list ({len(gene_names)} genes) "
                f"and HVGs ({len(hvg_genes_set)} genes). "
                "This may occur if adata.raw contains different genes than the preprocessed data. "
                "Try setting test_only_hvg=False or ensure adata.raw is None."
            )

        gene_mask &= hvg_mask  # Keep only HVGs

    # TIER 1: Apply SPARK-X standard filtering (expression-based) - ON SPARSE MATRIX
    percentage = params.sparkx_percentage
    min_total_counts = params.sparkx_min_total_counts

    # Calculate gene statistics on sparse matrix (efficient!)
    gene_totals, n_expressed = _calculate_sparse_gene_stats(sparse_counts)

    # Filter genes: must be expressed in at least percentage of cells AND have min total counts
    min_cells = int(np.ceil(n_spots * percentage))
    expr_mask = (n_expressed >= min_cells) & (gene_totals >= min_total_counts)

    gene_mask &= expr_mask  # Combine with previous filters

    # Apply combined filter mask to sparse matrix (still sparse!)
    if gene_mask.sum() < len(gene_names):
        filtered_sparse = sparse_counts[:, gene_mask]
        gene_names = [
            gene for gene, keep in zip(gene_names, gene_mask, strict=False) if keep
        ]
    else:
        filtered_sparse = sparse_counts

    # NOW convert filtered sparse matrix to dense (much smaller!)
    # copy=True ensures we don't modify original for dense input
    counts_matrix = to_dense(filtered_sparse, copy=True)

    # Ensure counts are non-negative integers
    counts_matrix = np.maximum(counts_matrix, 0).astype(int)

    # Update gene count after filtering
    n_genes = len(gene_names)

    # Transpose for SPARK format (genes × spots)
    counts_transposed = counts_matrix.T

    # Create spot names
    spot_names = [str(name) for name in adata.obs_names]

    # Wrap ALL R operations in thread lock and localconverter for proper contextvars handling
    # This prevents "Conversion rules missing" errors in multithreaded/async environments
    with openrlib.rlock:  # Thread safety lock
        with conversion.localconverter(default_converter):  # Conversion context
            # Import SPARK package inside context (FIX for contextvars issue)
            try:
                spark = importr("SPARK")
            except Exception as e:
                raise ImportError(
                    f"SPARK not installed in R. Install with: install.packages('SPARK'). Error: {e}"
                ) from e

            # Convert to R format (already in context)
            # Count matrix: genes × spots
            r_counts = ro.r.matrix(
                ro.IntVector(counts_transposed.flatten()),
                nrow=n_genes,
                ncol=n_spots,
                byrow=True,
            )
            r_counts.rownames = ro.StrVector(gene_names)
            r_counts.colnames = ro.StrVector(spot_names)

            # Coordinates as data.frame (SPARK requirement)
            coords_df = pd.DataFrame(coords_array, columns=["x", "y"], index=spot_names)
            r_coords = ro.r["data.frame"](
                x=ro.FloatVector(coords_df["x"]),
                y=ro.FloatVector(coords_df["y"]),
                row_names=ro.StrVector(coords_df.index),
            )

            try:
                # Execute SPARK-X analysis inside context (FIX for contextvars issue)
                # Keep suppress_output for MCP communication compatibility
                with suppress_output():
                    results = spark.sparkx(
                        count_in=r_counts,
                        locus_in=r_coords,
                        X_in=ro.NULL,  # No additional covariates (could be extended in future)
                        numCores=params.sparkx_n_cores,
                        option=params.sparkx_option,
                        verbose=False,  # Ensure verbose is off for cleaner MCP communication
                    )

                # Extract p-values from results (inside context for proper conversion)
                # SPARK-X returns res_mtest as a data.frame with columns:
                # - combinedPval: combined p-values across spatial kernels
                # - adjustedPval: BY-adjusted p-values (Benjamini-Yekutieli FDR correction)
                # Reference: SPARK R package documentation
                try:
                    pvals = results.rx2("res_mtest")
                    if pvals is None:
                        raise ProcessingError(
                            "SPARK-X returned None for res_mtest. "
                            "This may indicate the analysis failed silently."
                        )

                    # Verify expected data.frame format
                    is_dataframe = ro.r["is.data.frame"](pvals)[0]
                    if not is_dataframe:
                        raise ProcessingError(
                            "SPARK-X output format error. Requires SPARK >= 1.1.0."
                        )

                    # Extract combinedPval (raw p-values combined across kernels)
                    combined_pvals = ro.r["$"](pvals, "combinedPval")
                    if combined_pvals is None:
                        raise ProcessingError(
                            "SPARK-X res_mtest missing 'combinedPval' column. "
                            "This is required for spatial gene identification."
                        )
                    pval_list = [float(p) for p in combined_pvals]

                    # Extract adjustedPval (BY-corrected p-values from SPARK-X)
                    adjusted_pvals = ro.r["$"](pvals, "adjustedPval")
                    if adjusted_pvals is None:
                        raise ProcessingError(
                            "SPARK-X res_mtest missing 'adjustedPval' column. "
                            "This column contains BY-corrected p-values for multiple testing."
                        )
                    adjusted_pval_list = [float(p) for p in adjusted_pvals]

                    # Create results dataframe
                    results_df = pd.DataFrame(
                        {
                            "gene": gene_names[: len(pval_list)],
                            "pvalue": pval_list,
                            "adjusted_pvalue": adjusted_pval_list,  # BY-corrected by SPARK-X
                        }
                    )

                    # Warn if returned genes much fewer than input genes
                    if len(results_df) < n_genes * 0.5:
                        await ctx.warning(
                            f"SPARK-X returned results for only {len(results_df)}/{n_genes} genes. "
                            f"This may indicate a problem with the R environment, SPARK package, or input data. "
                            f"Consider checking R logs or trying SpatialDE as an alternative method."
                        )

                except Exception as e:
                    # P-value extraction failed - provide clear error message
                    raise ProcessingError(
                        f"SPARK-X p-value extraction failed: {e}\n\n"
                        f"Expected SPARK-X output format:\n"
                        f"SPARK-X output invalid. Requires SPARK >= 1.1.0."
                    ) from e

            except Exception as e:
                raise ProcessingError(f"SPARK-X analysis failed: {e}") from e

    # Sort by adjusted p-value
    results_df = results_df.sort_values("adjusted_pvalue")

    # Filter significant genes
    significant_genes_all = results_df[results_df["adjusted_pvalue"] < 0.05][
        "gene"
    ].tolist()

    # Limit for MCP response (full results stored in adata.var)
    limit = params.n_top_genes or DEFAULT_TOP_GENES_LIMIT
    significant_genes = significant_genes_all[:limit]

    # TIER 3: Housekeeping gene warnings (post-processing quality check)
    if params.warn_housekeeping and len(results_df) > 0:
        # Define housekeeping gene patterns (based on literature)
        housekeeping_patterns = [
            "RPS",  # Ribosomal protein small subunit
            "RPL",  # Ribosomal protein large subunit
            "Rps",  # Mouse ribosomal small
            "Rpl",  # Mouse ribosomal large
            "MT-",  # Mitochondrial (human)
            "mt-",  # Mitochondrial (mouse)
            "ACTB",  # Beta-actin
            "GAPDH",  # Glyceraldehyde-3-phosphate dehydrogenase
            "EEF1A1",  # Eukaryotic translation elongation factor 1 alpha 1
            "TUBA1B",  # Tubulin alpha 1b
            "B2M",  # Beta-2-microglobulin
        ]

        # Check top significant genes (up to 50)
        top_genes_to_check = results_df.head(50)["gene"].tolist()

        # Mark housekeeping genes
        housekeeping_genes = [
            gene
            for gene in top_genes_to_check
            if any(
                gene.startswith(pattern) or gene == pattern
                for pattern in housekeeping_patterns
            )
        ]

        n_housekeeping = len(housekeeping_genes)
        n_top = len(top_genes_to_check)
        housekeeping_ratio = n_housekeeping / n_top if n_top > 0 else 0

        # Warn if >30% are housekeeping genes
        if housekeeping_ratio > 0.3:
            await ctx.warning(
                f"WARNING:Housekeeping gene dominance detected: {n_housekeeping}/{n_top} ({housekeeping_ratio*100:.1f}%) of top genes are housekeeping genes.\n"
                f"   • Housekeeping genes found: {', '.join(housekeeping_genes[:10])}{'...' if len(housekeeping_genes) > 10 else ''}\n"
                f"   • These genes may not represent true spatial patterns\n"
                f"   • Recommendations:\n"
                f"     1. Use test_only_hvg=True to reduce housekeeping dominance (2024 best practice)\n"
                f"     2. Use filter_ribo_genes=True to filter ribosomal genes\n"
                f"     3. Focus on genes with clear biological relevance\n"
                f"   • Note: This is a quality warning, not an error"
            )

    # Store results in adata
    results_key = f"sparkx_results_{data_id}"
    adata.var["sparkx_pval"] = pd.Series(
        dict(zip(results_df["gene"], results_df["pvalue"], strict=False)),
        name="sparkx_pval",
    ).reindex(adata.var_names, fill_value=1.0)

    adata.var["sparkx_qval"] = pd.Series(
        dict(zip(results_df["gene"], results_df["adjusted_pvalue"], strict=False)),
        name="sparkx_qval",
    ).reindex(adata.var_names, fill_value=1.0)

    # Store scientific metadata for reproducibility
    from ..utils.adata_utils import store_analysis_metadata
    from ..utils.results_export import export_analysis_result

    store_analysis_metadata(
        adata,
        analysis_name="spatial_genes_sparkx",
        method="sparkx",
        parameters={
            "num_core": params.sparkx_n_cores,
            "percentage": params.sparkx_percentage,
            "min_total_counts": params.sparkx_min_total_counts,
            "option": params.sparkx_option,
            "filter_mt_genes": params.filter_mt_genes,
            "filter_ribo_genes": params.filter_ribo_genes,
            "test_only_hvg": params.test_only_hvg,
            "warn_housekeeping": params.warn_housekeeping,
        },
        results_keys={
            "var": ["sparkx_pval", "sparkx_qval"],
            "obs": [],
            "obsm": [],
            "uns": [],
        },
        statistics={
            "n_genes_analyzed": len(results_df),
            "n_significant_genes": len(significant_genes_all),
        },
    )

    # Export results to CSV for reproducibility
    export_analysis_result(adata, data_id, "spatial_genes_sparkx")

    # Note: Detailed statistics (gene_statistics, p_values, q_values) are excluded
    # from MCP response via Field(exclude=True) in SpatialVariableGenesResult.
    # Full results are accessible via adata.var['sparkx_pval', 'sparkx_qval'].

    result = SpatialVariableGenesResult(
        data_id=data_id,
        method="sparkx",
        n_genes_analyzed=len(results_df),
        n_significant_genes=len(significant_genes_all),
        spatial_genes=significant_genes,
        results_key=results_key,
    )

    return result
