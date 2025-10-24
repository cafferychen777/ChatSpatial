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

import logging
from typing import TYPE_CHECKING, Any, Dict

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)

from ..models.analysis import SpatialVariableGenesResult  # noqa: E402
from ..models.data import SpatialVariableGenesParameters  # noqa: E402
from ..utils import validate_var_column  # noqa: E402
from ..utils.error_handling import suppress_output  # noqa: E402


async def identify_spatial_genes(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None,
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
        data_store: Dictionary containing loaded datasets
        params: Method-specific parameters (see SpatialVariableGenesParameters)
        context: MCP context for logging and status updates

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
    if context:
        await context.info(
            f"Starting spatial variable genes identification using {params.method}"
        )

    # Get data
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Work with the original adata object, not a copy
    adata = data_store[data_id]["adata"]

    # Validate spatial coordinates
    if params.spatial_key not in adata.obsm:
        raise ValueError(
            f"Spatial coordinates not found in adata.obsm['{params.spatial_key}']"
        )

    # Extract spatial coordinates
    spatial_coords = adata.obsm[params.spatial_key]
    if spatial_coords.shape[1] != 2:
        raise ValueError("Spatial coordinates must be 2D (x, y)")

    # Log data information
    if context:
        await context.info(
            f"Processing data: {adata.n_obs} spots, {adata.n_vars} genes"
        )

    # Route to appropriate method
    if params.method == "spatialde":
        return await _identify_spatial_genes_spatialde(
            data_id, data_store, params, context
        )
    elif params.method == "sparkx":
        return await _identify_spatial_genes_sparkx(
            data_id, data_store, params, context
        )
    else:
        raise ValueError(
            f"Unsupported method: {params.method}. Available methods: spatialde, sparkx"
        )


async def _identify_spatial_genes_spatialde(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None,
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

    try:
        import NaiveDE
        import SpatialDE
        from SpatialDE.util import qvalue
    except ImportError:
        raise ImportError(
            "SpatialDE not installed. Install with: pip install spatialde"
        )

    adata = data_store[data_id]["adata"]

    if context:
        await context.info("Running SpatialDE analysis")

    # Prepare data
    coords = pd.DataFrame(
        adata.obsm[params.spatial_key][:, :2],  # Ensure 2D coordinates
        columns=["x", "y"],
        index=adata.obs_names,
    )

    # Get raw count data for SpatialDE preprocessing
    if context:
        await context.info("Preparing count data for SpatialDE")

    if adata.raw is not None:
        counts = pd.DataFrame(
            (adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X),
            columns=adata.raw.var_names,
            index=adata.obs_names,
        )
        if context:
            await context.info("Using raw count matrix from adata.raw")
    else:
        # Check if current data appears to be raw counts
        data_max = adata.X.max() if hasattr(adata.X, "max") else np.max(adata.X)
        if data_max <= 10:  # Likely already normalized
            raise ValueError(
                "SpatialDE requires raw count data but normalized data detected in adata.X. "
                "Please ensure adata.raw contains raw counts, or reload data before preprocessing."
            )

        counts = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
            columns=adata.var_names,
            index=adata.obs_names,
        )
        if context:
            await context.info("Using current count matrix from adata.X")

    # Step 0: Filter low-expression genes (Official recommendation)
    # SpatialDE README: "Filter practically unobserved genes" with total_counts >= 3
    gene_totals = counts.sum(axis=0)
    keep_genes = gene_totals >= 3
    n_filtered = keep_genes.sum()
    n_total = len(keep_genes)

    if n_filtered < n_total:
        counts = counts.loc[:, keep_genes]
        if context:
            await context.info(
                f"Filtered to {n_filtered}/{n_total} genes (total counts ≥ 3, official threshold)"
            )

    # Optional: Test only HVGs for performance (if requested via n_top_genes)
    # This is done BEFORE expensive preprocessing to save time
    if params.n_top_genes is not None and params.n_top_genes < counts.shape[1]:
        if "highly_variable" in adata.var.columns:
            # Prioritize HVGs if available
            hvg_mask = adata.var.loc[counts.columns, "highly_variable"]
            hvg_genes = counts.columns[hvg_mask]

            if len(hvg_genes) >= params.n_top_genes:
                # Use HVGs
                selected_genes = hvg_genes[: params.n_top_genes]
                counts = counts[selected_genes]
                if context:
                    await context.info(
                        f"Testing {params.n_top_genes} highly variable genes only (for performance)"
                    )
            else:
                # Not enough HVGs, select by expression
                top_genes = gene_totals.nlargest(params.n_top_genes).index
                counts = counts[top_genes]
                if context:
                    await context.info(
                        f"Testing top {params.n_top_genes} expressed genes (performance optimization)"
                    )
        else:
            # Select by expression
            top_genes = gene_totals.nlargest(params.n_top_genes).index
            counts = counts[top_genes]
            if context:
                await context.info(
                    f"Testing top {params.n_top_genes} expressed genes (performance optimization)"
                )

    # Performance warning for large gene sets
    n_genes = counts.shape[1]
    n_spots = counts.shape[0]
    if n_genes > 5000 and context:
        estimated_time = int(n_genes / 14000 * 10)  # Based on 14k genes = 10 min
        await context.warning(
            f"WARNING:Running SpatialDE on {n_genes} genes × {n_spots} spots may take {estimated_time}-{estimated_time*2} minutes.\n"
            f"   • Official benchmark: ~10 min for 14,000 genes\n"
            f"   • Tip: Use n_top_genes=1000-3000 to test fewer genes\n"
            f"   • Or use method='sparkx' for faster analysis (2-5 min)"
        )

    # Calculate total counts per spot for regress_out
    total_counts = pd.DataFrame(
        {"total_counts": counts.sum(axis=1)}, index=counts.index
    )

    if context:
        await context.info(
            "Applying SpatialDE official preprocessing workflow (variance stabilization + regress_out)"
        )

    # Step 1: Variance stabilization (Official: NaiveDE.stabilize)
    # This transforms count data to approximately normal distribution
    if context:
        await context.info("Step 1/3: Variance stabilization (NaiveDE.stabilize)")

    norm_expr = NaiveDE.stabilize(counts.T).T

    # Step 2: Regress out library size effects (Official: NaiveDE.regress_out)
    # This removes technical variation from sequencing depth differences
    if context:
        await context.info(
            "Step 2/3: Regressing out library size effects (NaiveDE.regress_out)"
        )

    resid_expr = NaiveDE.regress_out(
        total_counts, norm_expr.T, "np.log(total_counts)"
    ).T

    # Step 3: Run SpatialDE with preprocessed data
    if context:
        await context.info(
            f"Step 3/3: Running SpatialDE on {n_genes} genes × {n_spots} spots"
        )

    results = SpatialDE.run(coords.values, resid_expr)

    # Multiple testing correction
    results["qval"] = qvalue(results["pval"].values, pi0=0.1)

    # Sort by q-value
    results = results.sort_values("qval")

    # Filter significant genes
    significant_genes_all = results[results["qval"] < 0.05]["g"].tolist()

    # IMPORTANT: Limit returned gene list to avoid MCP token overflow
    # Return top 500 significant genes by default (full list stored in adata.var)
    MAX_GENES_TO_RETURN = 500
    significant_genes = significant_genes_all[:MAX_GENES_TO_RETURN]

    # Get top genes if requested
    if params.n_top_genes is not None:
        results = results.head(params.n_top_genes)
        significant_genes = results["g"].tolist()

    # Store results in adata
    results_key = f"spatialde_results_{data_id}"
    adata.var["spatialde_pval"] = results.set_index("g")["pval"]
    adata.var["spatialde_qval"] = results.set_index("g")["qval"]
    adata.var["spatialde_l"] = results.set_index("g")["l"]

    # Store scientific metadata for reproducibility
    from ..utils.metadata_storage import store_analysis_metadata

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
                results[results["qval"] < params.pvalue_threshold]
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
        "n_significant_genes_total": len(significant_genes_all),
        "n_significant_genes_returned": len(significant_genes),
        "note": "Full results stored in adata.var['spatialde_pval', 'spatialde_qval', 'spatialde_l']. Top 500 significant genes returned to avoid MCP token limits.",
    }

    result = SpatialVariableGenesResult(
        data_id=data_id,
        method="spatialde",
        n_genes_analyzed=len(results),
        n_significant_genes=len(significant_genes),
        spatial_genes=significant_genes,
        gene_statistics=gene_statistics,
        p_values=p_values,
        q_values=q_values,
        results_key=results_key,
        spatialde_results=spatialde_results,
    )

    if context:
        await context.info("SpatialDE analysis completed")
        await context.info(
            f"Found {len(significant_genes_all)} significant spatial genes"
        )
        if len(significant_genes_all) > MAX_GENES_TO_RETURN:
            await context.info(
                f"Returning top {MAX_GENES_TO_RETURN} genes (full results in adata.var)"
            )

    return result


async def _identify_spatial_genes_sparkx(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialVariableGenesParameters,
    context=None,
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

    try:
        from rpy2 import robjects as ro
        from rpy2.rinterface_lib import openrlib  # For thread safety
        from rpy2.robjects import conversion, default_converter
        from rpy2.robjects.packages import importr
    except ImportError:
        raise ImportError("rpy2 not installed. Install with: pip install rpy2")

    if context:
        await context.info("Running SPARK-X non-parametric analysis")

    adata = data_store[data_id]["adata"]

    # Prepare spatial coordinates - SPARK needs data.frame format
    coords_array = adata.obsm[params.spatial_key][:, :2].astype(float)
    n_spots, n_genes = adata.shape

    if context:
        await context.info(f"Preparing data: {n_spots} spots × {n_genes} genes")

    # Get count matrix - use raw counts if available, otherwise current matrix
    if adata.raw is not None:
        # Use raw counts for SPARK
        if hasattr(adata.raw.X, "toarray"):
            counts_matrix = adata.raw.X.toarray()
        else:
            counts_matrix = adata.raw.X.copy()
        # Use raw gene names
        gene_names = [str(name) for name in adata.raw.var_names]
        n_genes = len(gene_names)
    else:
        # Fallback to current matrix
        if hasattr(adata.X, "toarray"):
            counts_matrix = adata.X.toarray()
        else:
            counts_matrix = adata.X.copy()
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
        if context:
            await context.info(
                f"Made duplicate gene names unique (found {sum(1 for c in gene_counts.values() if c > 1)} duplicates)"
            )

    # Ensure counts are non-negative integers
    counts_matrix = np.maximum(counts_matrix, 0).astype(int)

    if context:
        await context.info(
            f"Using {'raw' if adata.raw is not None else 'current'} count matrix"
        )

    # ==================== Gene Filtering Pipeline ====================
    # Following SPARK-X paper best practices + 2024 literature recommendations

    # TIER 1: Mitochondrial gene filtering (SPARK-X paper standard practice)
    if params.filter_mt_genes:
        mt_mask = np.array([gene.startswith(("MT-", "mt-")) for gene in gene_names])
        n_mt_genes = mt_mask.sum()
        if n_mt_genes > 0:
            counts_matrix = counts_matrix[:, ~mt_mask]
            gene_names = [gene for gene, is_mt in zip(gene_names, mt_mask) if not is_mt]
            if context:
                await context.info(
                    f"Filtered {n_mt_genes} mitochondrial genes (SPARK-X paper standard)"
                )

    # TIER 1: Ribosomal gene filtering (optional)
    if params.filter_ribo_genes:
        ribo_mask = np.array(
            [gene.startswith(("RPS", "RPL", "Rps", "Rpl")) for gene in gene_names]
        )
        n_ribo_genes = ribo_mask.sum()
        if n_ribo_genes > 0:
            counts_matrix = counts_matrix[:, ~ribo_mask]
            gene_names = [
                gene for gene, is_ribo in zip(gene_names, ribo_mask) if not is_ribo
            ]
            if context:
                await context.info(
                    f"Filtered {n_ribo_genes} ribosomal genes (reduces housekeeping dominance)"
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
            raise ValueError(
                "test_only_hvg=True but no highly variable genes found. "
                "Please run preprocessing with n_top_genes parameter first."
            )

        # Filter gene_names to only include HVGs
        # gene_names comes from adata.raw if it exists, otherwise from adata.var
        hvg_mask = np.array([gene in hvg_genes_set for gene in gene_names])
        n_hvg = hvg_mask.sum()

        if n_hvg == 0:
            # No overlap between current gene list and HVGs
            # This can happen if adata.raw has different genes than adata.var
            raise ValueError(
                f"test_only_hvg=True but no overlap found between current gene list ({len(gene_names)} genes) "
                f"and HVGs ({len(hvg_genes_set)} genes). "
                "This may occur if adata.raw contains different genes than the preprocessed data. "
                "Try setting test_only_hvg=False or ensure adata.raw is None."
            )

        counts_matrix = counts_matrix[:, hvg_mask]
        gene_names = [gene for gene, is_hvg in zip(gene_names, hvg_mask) if is_hvg]

        if context:
            await context.info(
                f"Testing only {n_hvg} highly variable genes (2024 best practice - PMC11537352)"
            )

    # Update gene count after filtering
    n_genes = len(gene_names)

    # TIER 1: Apply SPARK-X standard filtering (expression-based)
    # Apply gene filtering based on SPARK-X parameters (like CreateSPARKObject in R)
    percentage = params.sparkx_percentage
    min_total_counts = params.sparkx_min_total_counts

    # Calculate total counts per gene
    gene_totals = counts_matrix.sum(axis=0)
    n_expressed = (counts_matrix > 0).sum(axis=0)

    # Filter genes: must be expressed in at least percentage of cells AND have min total counts
    min_cells = int(np.ceil(n_spots * percentage))
    keep_genes = (n_expressed >= min_cells) & (gene_totals >= min_total_counts)

    if keep_genes.sum() < len(gene_names):
        # Apply filtering
        counts_matrix = counts_matrix[:, keep_genes]
        gene_names = [gene for gene, keep in zip(gene_names, keep_genes) if keep]
        n_filtered = keep_genes.sum()

        if context:
            await context.info(
                f"Filtered to {n_filtered}/{len(keep_genes)} genes (>{percentage*100:.0f}% cells, >{min_total_counts} counts)"
            )

    # Update gene count after filtering
    n_genes = len(gene_names)

    # Transpose for SPARK format (genes × spots)
    counts_transposed = counts_matrix.T

    if context:
        await context.info(
            f"Count matrix shape: {counts_transposed.shape} (genes × spots)"
        )
        await context.info(f"Passing {n_genes} genes to SPARK-X for analysis")

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

            if context:
                await context.info(
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

                if context:
                    await context.info("SPARK-X analysis completed successfully")

                # Extract p-values from results (inside context for proper conversion)
                try:
                    pvals = results.rx2("res_mtest")
                    if pvals:
                        # SPARK-X returns res_mtest as a data.frame with columns:
                        # - combinedPval: combined p-values across kernels
                        # - adjustedPval: adjusted p-values
                        # We need to extract the combinedPval column

                        # Check if it's a data.frame (which it should be for SPARK-X)
                        is_dataframe = ro.r["is.data.frame"](pvals)[0]

                        if is_dataframe:
                            # Extract combinedPval column
                            combined_pvals = ro.r["$"](pvals, "combinedPval")
                            pval_list = [float(p) for p in combined_pvals]

                            # Also extract adjustedPval if available
                            adjusted_pvals = ro.r["$"](pvals, "adjustedPval")
                            adjusted_pval_list = [float(p) for p in adjusted_pvals]

                            # Create results dataframe
                            results_df = pd.DataFrame(
                                {
                                    "gene": gene_names[: len(pval_list)],
                                    "pvalue": pval_list,
                                    "adjusted_pvalue": adjusted_pval_list,
                                }
                            )
                        else:
                            # Fallback for older format (numeric vector)
                            pval_list = []
                            try:
                                pvals_numeric = ro.r["as.numeric"](pvals)
                            except Exception:
                                pvals_numeric = pvals

                            for i in range(len(pvals_numeric)):
                                val = pvals_numeric[i]
                                if hasattr(val, "__len__") and hasattr(
                                    val, "__getitem__"
                                ):
                                    try:
                                        pval_list.append(float(val[0]))
                                    except Exception:
                                        pval_list.append(float(val))
                                else:
                                    pval_list.append(float(val))

                            # Create results dataframe
                            results_df = pd.DataFrame(
                                {
                                    "gene": gene_names[: len(pval_list)],
                                    "pvalue": pval_list,
                                }
                            )

                            # Add adjusted p-values (Bonferroni correction)
                            n_tests = len(pval_list)
                            results_df["adjusted_pvalue"] = results_df["pvalue"].apply(
                                lambda p: min(p * n_tests, 1.0)
                            )

                        if context:
                            await context.info(
                                f"Extracted results for {len(results_df)} genes"
                            )
                            # Warn if returned genes much fewer than input genes
                            if len(results_df) < n_genes * 0.5:
                                await context.warning(
                                    f"SPARK-X returned results for only {len(results_df)}/{n_genes} genes. "
                                    f"This may indicate a problem with the R environment, SPARK package, or input data. "
                                    f"Consider checking R logs or trying SpatialDE as an alternative method."
                                )
                    else:
                        # SPARK-X results format not recognized - fail honestly instead of fake results
                        error_msg = (
                            "SPARK-X results format not recognized. Expected 'res_mtest' component with p-values. "
                            "This may indicate an issue with the SPARK-X R package, rpy2 integration, or input data. "
                            "Please check the R environment and SPARK-X installation."
                        )
                        if context:
                            await context.error(error_msg)
                        raise RuntimeError(error_msg)

                except Exception as e:
                    # P-value extraction failed - fail honestly instead of creating fake results
                    error_msg = (
                        f"SPARK-X p-value extraction failed: {e}. "
                        f"This indicates an issue with R-Python communication or SPARK-X result format. "
                        f"Please check rpy2 installation and R package versions."
                    )
                    if context:
                        await context.error(error_msg)
                    raise RuntimeError(error_msg)

            except Exception as e:
                if context:
                    await context.info(f"SPARK-X analysis failed: {e}")
                raise RuntimeError(f"SPARK-X analysis failed: {e}")

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
    if params.n_top_genes is not None:
        results_df = results_df.head(params.n_top_genes)
        significant_genes = results_df["gene"].tolist()

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
        if housekeeping_ratio > 0.3 and context:
            await context.warning(
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
    from ..utils.metadata_storage import store_analysis_metadata

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
        "n_significant_genes_total": len(significant_genes_all),
        "n_significant_genes_returned": len(significant_genes),
        "data_format": "genes_x_spots",
        "note": "Full results stored in adata.var['sparkx_pval', 'sparkx_qval']. Top 500 significant genes returned to avoid MCP token limits.",
    }

    result = SpatialVariableGenesResult(
        data_id=data_id,
        method="sparkx",
        n_genes_analyzed=len(results_df),
        n_significant_genes=len(significant_genes),
        spatial_genes=significant_genes,
        gene_statistics=gene_statistics,
        p_values=p_values,
        q_values=q_values,
        results_key=results_key,
        sparkx_results=sparkx_results,
    )

    if context:
        await context.info("SPARK-X analysis completed successfully")
        await context.info(f"Analyzed {len(results_df)} genes")
        await context.info(
            f"Found {len(significant_genes_all)} significant spatial genes (q < 0.05)"
        )
        if len(significant_genes_all) > MAX_GENES_TO_RETURN:
            await context.info(
                f"Returning top {MAX_GENES_TO_RETURN} genes (full results in adata.var)"
            )

    return result
