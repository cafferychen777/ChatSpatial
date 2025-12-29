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

from ..models.analysis import SpatialVariableGenesResult  # noqa: E402
from ..models.data import SpatialVariableGenesParameters  # noqa: E402
from ..utils import validate_var_column  # noqa: E402
from ..utils.adata_utils import require_spatial_coords  # noqa: E402
from ..utils.dependency_manager import require  # noqa: E402
from ..utils.exceptions import (  # noqa: E402
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ..utils.mcp_utils import suppress_output  # noqa: E402


def _extract_counts_matrix(X):
    """
    Extract count matrix from AnnData X, handling both sparse and dense matrices efficiently.

    For sparse matrices, toarray() already returns a new numpy array (no copy needed).
    For dense matrices, we need to copy to avoid modifying the original data.

    Args:
        X: Count matrix (sparse or dense)

    Returns:
        numpy.ndarray: Count matrix as a numpy array
    """
    if hasattr(X, "toarray"):
        # Sparse matrix: toarray() already returns a new array
        return X.toarray()
    else:
        # Dense matrix: need to copy to avoid modifying original
        return X.copy()


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
    await ctx.info(
        f"Starting spatial variable genes identification using {params.method}"
    )

    # Get data via ToolContext
    adata = await ctx.get_adata(data_id)

    # Validate spatial coordinates exist
    require_spatial_coords(adata, spatial_key=params.spatial_key)

    # Log data information
    await ctx.info(f"Processing data: {adata.n_obs} spots, {adata.n_vars} genes")

    # Route to appropriate method (pass ctx for logging)
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
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd

    # Use centralized dependency manager for consistent error handling
    require("spatialde")  # Raises ImportError with install instructions if missing
    import NaiveDE
    import SpatialDE
    from SpatialDE.util import qvalue

    # adata is now passed directly from main function (via ctx.get_adata())

    await ctx.info("Running SpatialDE analysis")

    # Prepare data
    coords = pd.DataFrame(
        adata.obsm[params.spatial_key][:, :2],  # Ensure 2D coordinates
        columns=["x", "y"],
        index=adata.obs_names,
    )

    # Get raw count data for SpatialDE preprocessing
    # OPTIMIZATION: Filter genes on SPARSE matrix first, then convert only selected genes to dense
    # This provides 95%+ memory savings for large datasets (tested with 5000×20000 genes)
    await ctx.info("Preparing count data for SpatialDE (memory-optimized)")

    # Get raw data source (keep as sparse matrix)
    if adata.raw is not None:
        raw_data = adata.raw.X
        var_names = adata.raw.var_names
        var_df = adata.var  # For HVG lookup
        await ctx.info("Using raw count matrix from adata.raw")
    else:
        # Check if current data appears to be raw counts
        data_max = adata.X.max() if hasattr(adata.X, "max") else np.max(adata.X)
        if data_max <= 10:  # Likely already normalized
            raise DataError(
                "SpatialDE requires raw counts. Data appears normalized (max<=10)."
            )

        raw_data = adata.X
        var_names = adata.var_names
        var_df = adata.var
        await ctx.info("Using current count matrix from adata.X")

    # Step 1: Filter low-expression genes ON SPARSE MATRIX (Official recommendation)
    # SpatialDE README: "Filter practically unobserved genes" with total_counts >= 3
    # This is done on sparse matrix to avoid memory overhead
    import scipy.sparse as sp

    is_sparse = sp.issparse(raw_data)

    if is_sparse:
        gene_totals = np.array(raw_data.sum(axis=0)).flatten()
    else:
        gene_totals = raw_data.sum(axis=0)

    keep_genes_mask = gene_totals >= 3
    selected_var_names = var_names[keep_genes_mask]
    n_filtered = keep_genes_mask.sum()
    n_total = len(keep_genes_mask)

    if n_filtered < n_total:
        await ctx.info(
            f"Filtered to {n_filtered}/{n_total} genes (total counts ≥ 3, official threshold)"
        )

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
                await ctx.info(
                    f"Selected {params.n_top_genes} highly variable genes (for performance)"
                )
            else:
                # Not enough HVGs, select by expression
                gene_totals_filtered = gene_totals[keep_genes_mask]
                top_indices = np.argsort(gene_totals_filtered)[-params.n_top_genes :][
                    ::-1
                ]
                final_genes = selected_var_names[top_indices]
                await ctx.info(
                    f"Selected top {params.n_top_genes} expressed genes (performance optimization)"
                )
        else:
            # Select by expression
            gene_totals_filtered = gene_totals[keep_genes_mask]
            top_indices = np.argsort(gene_totals_filtered)[-params.n_top_genes :][::-1]
            final_genes = selected_var_names[top_indices]
            await ctx.info(
                f"Selected top {params.n_top_genes} expressed genes (performance optimization)"
            )

    # Step 3: Slice sparse matrix to final genes, THEN convert to dense
    # This is where the memory optimization happens: only convert selected genes
    if adata.raw is not None:
        final_adata_subset = adata.raw[:, final_genes]
    else:
        final_adata_subset = adata[:, final_genes]

    # Now create DataFrame from the SUBSET (much smaller memory footprint)
    counts = pd.DataFrame(
        (
            final_adata_subset.X.toarray()
            if hasattr(final_adata_subset.X, "toarray")
            else final_adata_subset.X
        ),
        columns=final_adata_subset.var_names,
        index=final_adata_subset.obs_names,
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

    await ctx.info(
        "Applying SpatialDE official preprocessing workflow (variance stabilization + regress_out)"
    )

    # Step 1: Variance stabilization (Official: NaiveDE.stabilize)
    # This transforms count data to approximately normal distribution
    await ctx.info("Step 1/3: Variance stabilization (NaiveDE.stabilize)")

    norm_expr = NaiveDE.stabilize(counts.T).T

    # Step 2: Regress out library size effects (Official: NaiveDE.regress_out)
    # This removes technical variation from sequencing depth differences
    await ctx.info(
        "Step 2/3: Regressing out library size effects (NaiveDE.regress_out)"
    )

    resid_expr = NaiveDE.regress_out(
        total_counts, norm_expr.T, "np.log(total_counts)"
    ).T

    # Step 3: Run SpatialDE with preprocessed data
    await ctx.info(f"Step 3/3: Running SpatialDE on {n_genes} genes × {n_spots} spots")

    results = SpatialDE.run(coords.values, resid_expr)

    # Multiple testing correction using Storey q-value method
    # pi0 = proportion of genes under null hypothesis (no spatial pattern)
    if params.spatialde_pi0 is not None:
        # User-specified pi0 value
        results["qval"] = qvalue(results["pval"].values, pi0=params.spatialde_pi0)
        await ctx.info(
            f"Using user-specified pi0={params.spatialde_pi0} for q-value estimation"
        )
    else:
        # Adaptive pi0 estimation (SpatialDE default, recommended)
        # This estimates pi0 from the p-value distribution
        results["qval"] = qvalue(results["pval"].values)
        await ctx.info(
            "Using adaptive pi0 estimation for q-value calculation (recommended)"
        )

    # Sort by q-value
    results = results.sort_values("qval")

    # Filter significant genes
    significant_genes_all = results[results["qval"] < 0.05]["g"].tolist()

    # IMPORTANT: Limit returned gene list to avoid MCP token overflow
    # Return top 500 significant genes by default (full list stored in adata.var)
    MAX_GENES_TO_RETURN = 500
    significant_genes = significant_genes_all[:MAX_GENES_TO_RETURN]

    # Get top genes if requested
    # IMPORTANT: n_top_genes only limits the number of _significant_ genes returned
    # If no significant genes exist, return empty list (not non-significant genes)
    if params.n_top_genes is not None and len(significant_genes) > 0:
        significant_genes = significant_genes[: params.n_top_genes]

    # Store results in adata
    results_key = f"spatialde_results_{data_id}"
    adata.var["spatialde_pval"] = results.set_index("g")["pval"]
    adata.var["spatialde_qval"] = results.set_index("g")["qval"]
    adata.var["spatialde_l"] = results.set_index("g")["l"]

    # Store scientific metadata for reproducibility
    from ..utils.adata_utils import store_analysis_metadata

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

    # Create gene statistics dictionaries
    # IMPORTANT: Only include top genes to avoid exceeding MCP token limits
    # Return detailed statistics for top 100 genes only (full results in adata.var)
    MAX_STATS_TO_RETURN = 100
    top_stats_genes = significant_genes[:MAX_STATS_TO_RETURN]
    top_stats_results = results[results["g"].isin(top_stats_genes)]
    gene_statistics = dict(
        zip(top_stats_results["g"], top_stats_results["LLR"])
    )  # Log-likelihood ratio
    p_values = dict(zip(top_stats_results["g"], top_stats_results["pval"]))
    q_values = dict(zip(top_stats_results["g"], top_stats_results["qval"]))

    # Create SpatialDE-specific results
    # Only return summary statistics (top 10 genes) to avoid exceeding MCP token limit
    top_results = results.head(10)
    spatialde_results = {
        "top_genes_summary": {
            "genes": top_results["g"].tolist(),
            "pvalues": top_results["pval"].tolist(),
            "qvalues": top_results["qval"].tolist(),
            "log_likelihood_ratios": top_results["LLR"].tolist(),
        },
        "kernel": params.spatialde_kernel,
        "preprocessing_workflow": "Official: NaiveDE.stabilize + NaiveDE.regress_out",
        "n_genes_tested": n_genes,
        "n_spots": n_spots,
        "note": "Full results stored in adata.var['spatialde_pval', 'spatialde_qval', 'spatialde_l']. Top 500 significant genes returned to avoid MCP token limits.",
    }

    result = SpatialVariableGenesResult(
        data_id=data_id,
        method="spatialde",
        n_genes_analyzed=len(results),
        n_significant_genes=len(significant_genes_all),  # Total significant genes found
        n_returned_genes=len(significant_genes),  # Actually returned (may be truncated)
        spatial_genes=significant_genes,
        gene_statistics=gene_statistics,
        p_values=p_values,
        q_values=q_values,
        results_key=results_key,
        spatialde_results=spatialde_results,
    )

    await ctx.info("SpatialDE analysis completed")
    await ctx.info(f"Found {len(significant_genes_all)} significant spatial genes")
    if len(significant_genes_all) > MAX_GENES_TO_RETURN:
        await ctx.info(
            f"Returning top {MAX_GENES_TO_RETURN} genes (full results in adata.var)"
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
        - sparkx_num_core: Number of CPU cores for parallel processing
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
    # Import dependencies at runtime
    import numpy as np
    import pandas as pd

    # Use centralized dependency manager for consistent error handling
    require("rpy2")  # Raises ImportError with install instructions if missing
    from rpy2 import robjects as ro
    from rpy2.rinterface_lib import openrlib  # For thread safety
    from rpy2.robjects import conversion, default_converter
    from rpy2.robjects.packages import importr

    await ctx.info("Running SPARK-X non-parametric analysis")

    # adata is now passed directly from main function (via ctx.get_adata())

    # Prepare spatial coordinates - SPARK needs data.frame format
    coords_array = adata.obsm[params.spatial_key][:, :2].astype(float)
    n_spots, n_genes = adata.shape

    await ctx.info(f"Preparing data: {n_spots} spots × {n_genes} genes")

    # ==================== OPTIMIZED: Filter on sparse matrix, then convert ====================
    # Strategy: Keep data sparse throughout filtering, only convert final filtered result
    # Benefit: For 30k cells × 20k genes → 3k genes: save ~15GB memory

    # Get sparse count matrix - DO NOT convert to dense yet!
    if adata.raw is not None:
        sparse_counts = adata.raw.X  # Keep sparse!
        gene_names = [str(name) for name in adata.raw.var_names]
        n_genes = len(gene_names)
    else:
        sparse_counts = adata.X  # Keep sparse!
        gene_names = [str(name) for name in adata.var_names]
        n_genes = len(gene_names)

    # Ensure gene names are unique (required for SPARK-X R rownames)
    if len(gene_names) != len(set(gene_names)):
        from collections import Counter

        gene_counts = Counter(gene_names)
        unique_names = []
        seen_counts = {}
        for gene in gene_names:
            if gene_counts[gene] > 1:
                # Add suffix for duplicates
                if gene not in seen_counts:
                    seen_counts[gene] = 0
                    unique_names.append(gene)
                else:
                    seen_counts[gene] += 1
                    unique_names.append(f"{gene}_{seen_counts[gene]}")
            else:
                unique_names.append(gene)
        gene_names = unique_names
        await ctx.info(
            f"Made duplicate gene names unique (found {sum(1 for c in gene_counts.values() if c > 1)} duplicates)"
        )

    await ctx.info(
        f"Using {'raw' if adata.raw is not None else 'current'} count matrix (keeping sparse for efficient filtering)"
    )

    # ==================== Gene Filtering Pipeline (ON SPARSE MATRIX) ====================
    # Following SPARK-X paper best practices + 2024 literature recommendations
    # All filtering done on sparse matrix to minimize memory usage

    # Initialize gene mask (all True = keep all genes initially)
    gene_mask = np.ones(len(gene_names), dtype=bool)

    # Get var annotation source (prefer raw for complete gene annotations)
    var_source = adata.raw if adata.raw is not None else adata

    # TIER 1: Mitochondrial gene filtering (SPARK-X paper standard practice)
    # Reuse preprocessing annotations when available for consistency
    if params.filter_mt_genes:
        mt_mask = None

        # Try to reuse preprocessing annotations (elegant consistency)
        if "mt" in var_source.var.columns:
            mt_mask = var_source.var["mt"].values
            annotation_source = "preprocessing"
        else:
            # Fallback to pattern-based detection
            mt_mask = np.array([gene.startswith(("MT-", "mt-")) for gene in gene_names])
            annotation_source = "pattern detection"

        n_mt_genes = mt_mask.sum()
        if n_mt_genes > 0:
            gene_mask &= ~mt_mask  # Exclude MT genes
            await ctx.info(
                f"Marked {n_mt_genes} mitochondrial genes for filtering "
                f"(via {annotation_source}, SPARK-X paper standard)"
            )

    # TIER 1: Ribosomal gene filtering (optional)
    # Reuse preprocessing annotations when available for consistency
    if params.filter_ribo_genes:
        ribo_mask = None

        # Try to reuse preprocessing annotations (elegant consistency)
        if "ribo" in var_source.var.columns:
            ribo_mask = var_source.var["ribo"].values
            annotation_source = "preprocessing"
        else:
            # Fallback to pattern-based detection
            ribo_mask = np.array(
                [gene.startswith(("RPS", "RPL", "Rps", "Rpl")) for gene in gene_names]
            )
            annotation_source = "pattern detection"

        n_ribo_genes = ribo_mask.sum()
        if n_ribo_genes > 0:
            gene_mask &= ~ribo_mask  # Exclude ribosomal genes
            await ctx.info(
                f"Marked {n_ribo_genes} ribosomal genes for filtering "
                f"(via {annotation_source}, reduces housekeeping dominance)"
            )

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
        await ctx.info(
            f"Marked for HVG-only testing: {n_hvg} highly variable genes (2024 best practice - PMC11537352)"
        )

        # Smart detection: Check if preprocessing already excluded mt/ribo from HVGs
        # This provides transparency about the filtering pipeline
        preprocessing_info = adata.uns.get("preprocessing", {})
        gene_annot = preprocessing_info.get("gene_annotations", {})

        excluded_from_hvg = []
        if gene_annot.get("mt_column") and gene_annot.get("n_mt_genes", 0) > 0:
            # Check if any mt genes are in HVGs (if not, preprocessing excluded them)
            if "mt" in var_source.var.columns:
                mt_in_hvg = (
                    var_source.var["mt"] & adata.var.get("highly_variable", False)
                ).sum()
                if mt_in_hvg == 0:
                    excluded_from_hvg.append(
                        f"mitochondrial ({gene_annot.get('n_mt_genes', 0)} genes)"
                    )

        if gene_annot.get("ribo_column") and gene_annot.get("n_ribo_genes", 0) > 0:
            if "ribo" in var_source.var.columns:
                ribo_in_hvg = (
                    var_source.var["ribo"] & adata.var.get("highly_variable", False)
                ).sum()
                if ribo_in_hvg == 0:
                    excluded_from_hvg.append(
                        f"ribosomal ({gene_annot.get('n_ribo_genes', 0)} genes)"
                    )

        if excluded_from_hvg:
            await ctx.info(
                f"Note: {' and '.join(excluded_from_hvg)} already excluded from HVGs "
                "during preprocessing (unified filtering pipeline)"
            )

    # TIER 1: Apply SPARK-X standard filtering (expression-based) - ON SPARSE MATRIX
    percentage = params.sparkx_percentage
    min_total_counts = params.sparkx_min_total_counts

    # Calculate on sparse matrix (efficient!)
    # sum(axis=0) on sparse matrix returns matrix, need to convert to 1D array
    gene_totals = np.array(sparse_counts.sum(axis=0)).flatten()
    # Count non-zeros per gene (efficient on sparse)
    n_expressed = np.array((sparse_counts > 0).sum(axis=0)).flatten()

    # Filter genes: must be expressed in at least percentage of cells AND have min total counts
    min_cells = int(np.ceil(n_spots * percentage))
    expr_mask = (n_expressed >= min_cells) & (gene_totals >= min_total_counts)

    gene_mask &= expr_mask  # Combine with previous filters

    n_filtered_out = len(gene_names) - gene_mask.sum()
    if n_filtered_out > 0:
        await ctx.info(
            f"Marked {n_filtered_out} genes for filtering (expression threshold: >{percentage*100:.0f}% cells, >{min_total_counts} counts)"
        )

    # Apply combined filter mask to sparse matrix (still sparse!)
    if gene_mask.sum() < len(gene_names):
        filtered_sparse = sparse_counts[:, gene_mask]
        gene_names = [gene for gene, keep in zip(gene_names, gene_mask) if keep]
        n_genes_after_filter = len(gene_names)

        await ctx.info(
            f"Filtered {len(gene_mask)} → {n_genes_after_filter} genes (saved {len(gene_mask) - n_genes_after_filter} genes = {(1 - n_genes_after_filter/len(gene_mask))*100:.1f}% reduction)"
        )
    else:
        filtered_sparse = sparse_counts
        n_genes_after_filter = len(gene_names)

    # NOW convert filtered sparse matrix to dense (much smaller!)
    counts_matrix = _extract_counts_matrix(filtered_sparse)

    # Ensure counts are non-negative integers
    counts_matrix = np.maximum(counts_matrix, 0).astype(int)

    await ctx.info(
        f"Converted filtered sparse matrix to dense: {counts_matrix.shape} (saves ~{(len(gene_mask) - n_genes_after_filter) * n_spots * 4 / (1024**3):.1f}GB vs converting before filtering)"
    )

    # Update gene count after filtering
    n_genes = len(gene_names)

    # Transpose for SPARK format (genes × spots)
    counts_transposed = counts_matrix.T

    await ctx.info(f"Count matrix shape: {counts_transposed.shape} (genes × spots)")
    await ctx.info(f"Passing {n_genes} genes to SPARK-X for analysis")

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
                )

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

            await ctx.info(
                "Running SPARK-X analysis using sparkx function (output suppressed for MCP compatibility)"
            )

            try:
                # Execute SPARK-X analysis inside context (FIX for contextvars issue)
                # Keep suppress_output for MCP communication compatibility
                with suppress_output():
                    results = spark.sparkx(
                        count_in=r_counts,
                        locus_in=r_coords,
                        X_in=ro.NULL,  # No additional covariates (could be extended in future)
                        numCores=params.sparkx_num_core,
                        option=params.sparkx_option,
                        verbose=False,  # Ensure verbose is off for cleaner MCP communication
                    )

                await ctx.info("SPARK-X analysis completed successfully")

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

                    await ctx.info(
                        f"Extracted results for {len(results_df)} genes "
                        f"(p-values BY-corrected by SPARK-X)"
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
                    error_msg = (
                        f"SPARK-X p-value extraction failed: {e}\n\n"
                        f"Expected SPARK-X output format:\n"
                        f"SPARK-X output invalid. Requires SPARK >= 1.1.0."
                    )
                    await ctx.error(error_msg)
                    raise ProcessingError(error_msg)

            except Exception as e:
                await ctx.info(f"SPARK-X analysis failed: {e}")
                raise ProcessingError(f"SPARK-X analysis failed: {e}")

    # Sort by adjusted p-value
    results_df = results_df.sort_values("adjusted_pvalue")

    # Filter significant genes
    significant_genes_all = results_df[results_df["adjusted_pvalue"] < 0.05][
        "gene"
    ].tolist()

    # IMPORTANT: Limit returned gene list to avoid MCP token overflow
    # Return top 500 significant genes by default (full list stored in adata.var)
    MAX_GENES_TO_RETURN = 500
    significant_genes = significant_genes_all[:MAX_GENES_TO_RETURN]

    # Get top genes if requested
    # IMPORTANT: n_top_genes only limits the number of _significant_ genes returned
    # If no significant genes exist, return empty list (not non-significant genes)
    if params.n_top_genes is not None and len(significant_genes) > 0:
        significant_genes = significant_genes[: params.n_top_genes]

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
        dict(zip(results_df["gene"], results_df["pvalue"])), name="sparkx_pval"
    ).reindex(adata.var_names, fill_value=1.0)

    adata.var["sparkx_qval"] = pd.Series(
        dict(zip(results_df["gene"], results_df["adjusted_pvalue"])), name="sparkx_qval"
    ).reindex(adata.var_names, fill_value=1.0)

    # Store scientific metadata for reproducibility
    from ..utils.adata_utils import store_analysis_metadata

    store_analysis_metadata(
        adata,
        analysis_name="spatial_genes_sparkx",
        method="sparkx",
        parameters={
            "num_core": params.sparkx_num_core,
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

    # Create gene statistics dictionaries
    # IMPORTANT: Only include top genes to avoid exceeding MCP token limits
    # Return detailed statistics for top 100 genes only (full results in adata.var)
    MAX_STATS_TO_RETURN = 100
    top_stats_genes = significant_genes[:MAX_STATS_TO_RETURN]
    top_stats_results = results_df[results_df["gene"].isin(top_stats_genes)]
    gene_statistics = dict(zip(top_stats_results["gene"], top_stats_results["pvalue"]))
    p_values = dict(zip(top_stats_results["gene"], top_stats_results["pvalue"]))
    q_values = dict(
        zip(top_stats_results["gene"], top_stats_results["adjusted_pvalue"])
    )

    # Create SPARK-X specific results
    # Only return summary statistics (top 10 genes) to avoid exceeding MCP token limit
    top_results = results_df.head(10)
    sparkx_results = {
        "top_genes_summary": {
            "genes": top_results["gene"].tolist(),
            "pvalues": top_results["pvalue"].tolist(),
            "adjusted_pvalues": top_results["adjusted_pvalue"].tolist(),
            "combined_pvalues": (
                top_results["combined_pvalue"].tolist()
                if "combined_pvalue" in top_results.columns
                else None
            ),
        },
        "method": "sparkx",
        "num_core": params.sparkx_num_core,
        "option": params.sparkx_option,
        "data_format": "genes_x_spots",
        "note": "Full results stored in adata.var['sparkx_pval', 'sparkx_qval']. Top 500 significant genes returned to avoid MCP token limits.",
    }

    result = SpatialVariableGenesResult(
        data_id=data_id,
        method="sparkx",
        n_genes_analyzed=len(results_df),
        n_significant_genes=len(significant_genes_all),  # Total significant genes found
        n_returned_genes=len(significant_genes),  # Actually returned (may be truncated)
        spatial_genes=significant_genes,
        gene_statistics=gene_statistics,
        p_values=p_values,
        q_values=q_values,
        results_key=results_key,
        sparkx_results=sparkx_results,
    )

    await ctx.info("SPARK-X analysis completed successfully")
    await ctx.info(f"Analyzed {len(results_df)} genes")
    await ctx.info(
        f"Found {len(significant_genes_all)} significant spatial genes (q < 0.05)"
    )
    if len(significant_genes_all) > MAX_GENES_TO_RETURN:
        await ctx.info(
            f"Returning top {MAX_GENES_TO_RETURN} genes (full results in adata.var)"
        )

    return result
