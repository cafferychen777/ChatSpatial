"""
Differential expression analysis tools for spatial transcriptomics data.
"""

import numpy as np
import pandas as pd
import scanpy as sc

from ..models.analysis import DifferentialExpressionResult
from ..models.data import DifferentialExpressionParameters
from ..spatial_mcp_adapter import ToolContext
from ..utils import validate_obs_column
from ..utils.adata_utils import (
    check_is_integer_counts,
    get_raw_data_source,
    shallow_copy_adata,
    store_analysis_metadata,
    to_dense,
)
from ..utils.dependency_manager import require
from ..utils.exceptions import (
    DataError,
    ParameterError,
    ProcessingError,
)
from ..utils.results_export import export_analysis_result


def _build_de_key(
    method: str,
    group1: str | None,
    group2: str | None,
) -> str:
    """Build a parametric analysis key for differential expression.

    Encodes method + comparison so that multiple DE runs on the same
    dataset coexist in metadata, export, and server cache.
    """
    if group1 is None:
        return f"de_{method}_all_groups"
    g2 = group2 or "rest"
    return f"de_{method}_{group1}_vs_{g2}"


async def differential_expression(
    data_id: str,
    ctx: ToolContext,
    params: DifferentialExpressionParameters,
) -> DifferentialExpressionResult:
    """Perform differential expression analysis.

    Args:
        data_id: Dataset ID
        ctx: Tool context for data access and logging
        params: Differential expression parameters

    Returns:
        Differential expression analysis result
    """
    # Extract parameters from params object
    group_key = params.group_key
    group1 = params.group1
    group2 = params.group2
    method = params.method
    n_top_genes = params.n_top_genes
    pseudocount = params.pseudocount
    min_cells = params.min_cells

    # Dispatch to pydeseq2 if requested
    if method == "pydeseq2":
        return await _run_pydeseq2(data_id, ctx, params)

    # Get AnnData directly via ToolContext (no redundant dict wrapping)
    adata = await ctx.get_adata(data_id)

    # Check if the group_key exists in adata.obs
    validate_obs_column(adata, group_key, "Group")

    # Check if dtype conversion is needed (numba doesn't support float16)
    # Defer conversion to after subsetting for memory efficiency
    needs_dtype_fix = hasattr(adata.X, "dtype") and adata.X.dtype == np.float16

    # If group1 is None, find markers for all groups
    if group1 is None:

        # Filter out groups with too few cells (user-configurable threshold)
        group_sizes = adata.obs[group_key].value_counts()
        # min_cells is now a parameter (default=3, minimum for Wilcoxon test)
        valid_groups = group_sizes[group_sizes >= min_cells]
        skipped_groups = group_sizes[group_sizes < min_cells]

        # Warn about skipped groups
        if len(skipped_groups) > 0:
            skipped_list = "\n".join(
                [f"  - {g}: {n} cell(s)" for g, n in skipped_groups.items()]
            )
            await ctx.warning(
                f"Skipped {len(skipped_groups)} group(s) with <{min_cells} cells:\n{skipped_list}"
            )

        # Check if any valid groups remain
        if len(valid_groups) == 0:
            all_sizes = "\n".join(
                [f"  • {g}: {n} cell(s)" for g, n in group_sizes.items()]
            )
            raise DataError(
                f"All groups have <{min_cells} cells. Cannot perform {method} test.\n\n"
                f"Group sizes:\n{all_sizes}\n\n"
                f"Try: find_markers(group_key='leiden') or merge small groups"
            )

        # Filter data to only include valid groups
        adata_filtered = shallow_copy_adata(
            adata[adata.obs[group_key].isin(valid_groups.index)]
        )

        # Convert dtype after subsetting (4x more memory efficient than copying first)
        if needs_dtype_fix:
            adata_filtered.X = adata_filtered.X.astype(np.float32)

        # Run rank_genes_groups on filtered data
        sc.tl.rank_genes_groups(
            adata_filtered,
            groupby=group_key,
            method=method,
            n_genes=n_top_genes,
            reference="rest",
        )

        # Get all groups (from filtered data)
        groups = adata_filtered.obs[group_key].unique()

        # Collect top genes from all groups
        all_top_genes = []
        if (
            "rank_genes_groups" in adata_filtered.uns
            and "names" in adata_filtered.uns["rank_genes_groups"]
        ):
            gene_names = adata_filtered.uns["rank_genes_groups"]["names"]
            for group in groups:
                if str(group) in gene_names.dtype.names:
                    genes = list(gene_names[str(group)][:n_top_genes])
                    all_top_genes.extend(genes)

        # Remove duplicates while preserving order
        seen = set()
        top_genes = []
        for gene in all_top_genes:
            if gene not in seen:
                seen.add(gene)
                top_genes.append(gene)

        # Limit to n_top_genes
        top_genes = top_genes[:n_top_genes]

        # Copy results back to original adata for persistence
        adata.uns["rank_genes_groups"] = adata_filtered.uns["rank_genes_groups"]

        # Per-run copy for provenance (shared key kept for viz)
        analysis_key = _build_de_key(method, None, None)
        per_run_key = f"rank_genes_groups_{analysis_key.removeprefix('de_')}"
        adata.uns[per_run_key] = adata_filtered.uns["rank_genes_groups"]

        # Store metadata for scientific provenance tracking
        store_analysis_metadata(
            adata,
            analysis_name=analysis_key,
            method=method,
            parameters={
                "group_key": group_key,
                "comparison_type": "all_groups",
                "n_top_genes": n_top_genes,
            },
            results_keys={"uns": [per_run_key]},
            statistics={
                "method": method,
                "n_groups": len(groups),
                "groups": list(map(str, groups)),
                "n_cells_analyzed": adata_filtered.n_obs,
                "n_genes_analyzed": adata_filtered.n_vars,
            },
        )

        # Export results to CSV for reproducibility
        export_analysis_result(adata, data_id, analysis_key)

        return DifferentialExpressionResult(
            data_id=data_id,
            comparison=f"All groups in {group_key}",
            n_genes=len(top_genes),
            top_genes=top_genes,
            statistics={
                "method": method,
                "n_groups": len(groups),
                "groups": list(map(str, groups)),
            },
        )

    # Original logic for specific group comparison
    # Check if the groups exist in the group_key
    if group1 not in adata.obs[group_key].values:
        raise ParameterError(f"Group '{group1}' not found in adata.obs['{group_key}']")

    # Special case for 'rest' as group2 or if group2 is None
    use_rest_as_reference = False
    if group2 is None or group2 == "rest":
        use_rest_as_reference = True
        group2 = "rest"  # Set it explicitly for display purposes
    elif group2 not in adata.obs[group_key].values:
        raise ParameterError(f"Group '{group2}' not found in adata.obs['{group_key}']")

    # Perform differential expression analysis

    # Prepare the AnnData object for analysis
    if use_rest_as_reference:
        # Use the full AnnData object when comparing with 'rest'
        temp_adata = shallow_copy_adata(adata)
    else:
        # Create a temporary copy of the AnnData object with only the two groups
        temp_adata = shallow_copy_adata(
            adata[adata.obs[group_key].isin([group1, group2])]
        )

    # Convert dtype after subsetting (4x more memory efficient than copying first)
    if needs_dtype_fix:
        temp_adata.X = temp_adata.X.astype(np.float32)

    # Run rank_genes_groups
    sc.tl.rank_genes_groups(
        temp_adata,
        groupby=group_key,
        groups=[group1],
        reference="rest" if use_rest_as_reference else group2,
        method=method,
        n_genes=n_top_genes,
    )

    # Extract results

    # Get the top genes
    top_genes = []
    if (
        hasattr(temp_adata, "uns")
        and "rank_genes_groups" in temp_adata.uns
        and "names" in temp_adata.uns["rank_genes_groups"]
    ):
        # Get the top genes for the first group (should be group1)
        gene_names = temp_adata.uns["rank_genes_groups"]["names"]
        if group1 in gene_names.dtype.names:
            top_genes = list(gene_names[group1][:n_top_genes])
        else:
            # If group1 is not in the names, use the first column
            top_genes = list(gene_names[gene_names.dtype.names[0]][:n_top_genes])

    # If no genes were found, fail honestly
    if not top_genes:
        raise ProcessingError(
            f"No DE genes found between {group1} and {group2}. "
            f"Check sample sizes and expression differences."
        )

    # Get statistics
    n_cells_group1 = np.sum(adata.obs[group_key] == group1)
    if use_rest_as_reference:
        n_cells_group2 = adata.n_obs - n_cells_group1  # All cells except group1
    else:
        n_cells_group2 = np.sum(adata.obs[group_key] == group2)

    # Get p-values from scanpy results
    pvals = []
    if (
        hasattr(temp_adata, "uns")
        and "rank_genes_groups" in temp_adata.uns
        and "pvals_adj" in temp_adata.uns["rank_genes_groups"]
        and group1 in temp_adata.uns["rank_genes_groups"]["pvals_adj"].dtype.names
    ):
        pvals = list(
            temp_adata.uns["rank_genes_groups"]["pvals_adj"][group1][:n_top_genes]
        )

    # Calculate fold change from expression data
    # Prefer raw counts for unbiased FC; fall back to normalized data with warning
    raw_result = get_raw_data_source(adata, prefer_complete_genes=True)
    if not getattr(raw_result, "is_integer_counts", True):
        await ctx.warning(
            f"Using {raw_result.source} (not raw counts) for fold change "
            "calculation. FC values reflect normalized expression ratios, "
            "not count-based fold changes."
        )
    raw_X = raw_result.X
    raw_var_names = raw_result.var_names
    log2fc_values = []

    # Create masks for the two groups
    if use_rest_as_reference:
        group1_mask = adata.obs[group_key] == group1
        group2_mask = ~group1_mask
    else:
        group1_mask = adata.obs[group_key] == group1
        group2_mask = adata.obs[group_key] == group2

    # CRITICAL: Normalize by library size to avoid composition bias
    lib_sizes = np.asarray(raw_X.sum(axis=1)).ravel()

    # Guard against zero library sizes (empty cells/spots) which cause inf
    zero_lib_mask = lib_sizes == 0
    if zero_lib_mask.any():
        n_zero = int(zero_lib_mask.sum())
        await ctx.warning(
            f"{n_zero} cells have zero total counts and will be excluded "
            "from fold change calculation"
        )
        lib_sizes_nonzero = lib_sizes[~zero_lib_mask]
        if len(lib_sizes_nonzero) == 0:
            median_lib_size = 0.0
        else:
            median_lib_size = float(np.median(lib_sizes_nonzero))
    else:
        median_lib_size = float(np.median(lib_sizes))

    # Calculate fold change for all top genes (batch + vectorized)
    # Pre-filter genes that exist in raw data using set for O(1) lookup
    var_names_set = set(raw_var_names)
    genes_in_raw = [g for g in top_genes if g in var_names_set]
    genes_missing = [g for g in top_genes if g not in var_names_set]

    # Warn about missing genes once (not per-gene in loop)
    if genes_missing:
        await ctx.warning(
            f"{len(genes_missing)} genes not found in raw data, "
            "skipping fold change calculation for them"
        )

    if genes_in_raw and median_lib_size > 0:
        # Batch extract expression matrix for all genes at once
        gene_indices = [list(raw_var_names).index(g) for g in genes_in_raw]
        gene_expr_matrix = to_dense(raw_X[:, gene_indices])  # cells × genes

        # Vectorized normalization by library size (exclude zero-lib cells)
        valid_mask = ~zero_lib_mask
        lib_size_factors = np.where(valid_mask, median_lib_size / lib_sizes, 0.0)
        gene_norm_matrix = gene_expr_matrix * lib_size_factors[:, np.newaxis]

        # Vectorized mean calculation for each group (only valid cells)
        g1_valid = group1_mask & valid_mask
        g2_valid = group2_mask & valid_mask
        mean_group1 = (
            gene_norm_matrix[g1_valid].mean(axis=0)
            if g1_valid.any()
            else np.zeros(len(genes_in_raw))
        )
        mean_group2 = (
            gene_norm_matrix[g2_valid].mean(axis=0)
            if g2_valid.any()
            else np.zeros(len(genes_in_raw))
        )

        # Vectorized log2 fold change calculation
        log2fc_array = np.log2(
            (mean_group1 + pseudocount) / (mean_group2 + pseudocount)
        )

        # Build result dict mapping gene -> log2fc
        gene_to_log2fc = dict(zip(genes_in_raw, log2fc_array.astype(float)))
        # Preserve order of top_genes, None for missing genes
        log2fc_values = [gene_to_log2fc.get(g) for g in top_genes]
    else:
        log2fc_values = [None] * len(top_genes)

    # Calculate mean log2fc (filtering out None values)
    valid_log2fc = [fc for fc in log2fc_values if fc is not None]
    mean_log2fc = np.mean(valid_log2fc) if valid_log2fc else None
    median_pvalue = np.median(pvals) if pvals else None

    # Warn if fold change values are suspiciously high (indicating calculation errors)
    if mean_log2fc is not None and abs(mean_log2fc) > 10:
        await ctx.warning(
            f"Extreme fold change: mean log2FC = {mean_log2fc:.2f} (>1024x). "
            f"May indicate sparse expression or low cell counts."
        )

    # Create statistics dictionary
    statistics = {
        "method": method,
        "n_cells_group1": int(n_cells_group1),
        "n_cells_group2": int(n_cells_group2),
        "mean_log2fc": float(mean_log2fc) if mean_log2fc is not None else None,
        "median_pvalue": float(median_pvalue) if median_pvalue is not None else None,
    }

    # Create comparison string
    comparison = f"{group1} vs {group2}"

    # Copy results back to original adata for persistence
    adata.uns["rank_genes_groups"] = temp_adata.uns["rank_genes_groups"]

    # Per-run copy for provenance (shared key kept for viz)
    analysis_key = _build_de_key(method, group1, group2)
    per_run_key = f"rank_genes_groups_{analysis_key.removeprefix('de_')}"
    adata.uns[per_run_key] = temp_adata.uns["rank_genes_groups"]

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name=analysis_key,
        method=method,
        parameters={
            "group_key": group_key,
            "group1": group1,
            "group2": group2,
            "comparison_type": "specific_groups",
            "n_top_genes": n_top_genes,
            "pseudocount": pseudocount,  # Track for reproducibility
        },
        results_keys={"uns": [per_run_key]},
        statistics={
            "method": method,
            "group1": group1,
            "group2": group2,
            "n_cells_group1": int(n_cells_group1),
            "n_cells_group2": int(n_cells_group2),
            "n_genes_analyzed": temp_adata.n_vars,
            "mean_log2fc": float(mean_log2fc) if mean_log2fc is not None else None,
            "median_pvalue": (
                float(median_pvalue) if median_pvalue is not None else None
            ),
            "pseudocount_used": pseudocount,  # Document in statistics
        },
    )

    # Export results to CSV for reproducibility
    export_analysis_result(adata, data_id, analysis_key)

    return DifferentialExpressionResult(
        data_id=data_id,
        comparison=comparison,
        n_genes=len(top_genes),
        top_genes=top_genes,
        statistics=statistics,
    )


async def _run_pydeseq2(
    data_id: str,
    ctx: ToolContext,
    params: DifferentialExpressionParameters,
) -> DifferentialExpressionResult:
    """Run PyDESeq2 pseudobulk differential expression analysis.

    This function performs pseudobulk aggregation by summing raw counts within
    each sample/group combination, then uses PyDESeq2 for DE analysis.

    Args:
        data_id: Dataset ID
        ctx: Tool context for data access and logging
        params: Differential expression parameters

    Returns:
        Differential expression analysis result

    Raises:
        ParameterError: If sample_key is not provided
        ImportError: If pydeseq2 is not installed
    """
    # Validate sample_key is provided
    if params.sample_key is None:
        raise ParameterError(
            "sample_key is required for pydeseq2 method.\n"
            "Provide a column in adata.obs that identifies biological replicates "
            "(e.g., 'sample', 'patient_id', 'batch').\n"
            "Example: find_markers(group_key='cell_type', method='pydeseq2', "
            "sample_key='sample')"
        )

    # Get data
    adata = await ctx.get_adata(data_id)

    # Validate columns
    validate_obs_column(adata, params.group_key, "Group")
    validate_obs_column(adata, params.sample_key, "Sample")

    # Get raw counts (required for DESeq2)
    # Use get_raw_data_source (single source of truth) with require_integer_counts
    # Keep sparse to avoid memory explosion - downstream handles both formats
    raw_result = get_raw_data_source(
        adata, prefer_complete_genes=False, require_integer_counts=True
    )
    raw_X = raw_result.X
    var_names = raw_result.var_names

    # Validate counts are integers (DESeq2 requirement)
    is_int, _, _ = check_is_integer_counts(raw_X)
    if not is_int:
        await ctx.warning(
            "Data appears to be normalized. DESeq2 requires raw integer counts. "
            "Results may be inaccurate."
        )

    # Determine comparison groups
    group_key = params.group_key
    sample_key = params.sample_key
    group1 = params.group1
    group2 = params.group2

    # If group1 is None, find first two groups for pairwise comparison
    unique_groups = sorted(adata.obs[group_key].unique().astype(str))
    if group1 is None:
        if len(unique_groups) < 2:
            raise DataError(
                f"Need at least 2 groups for DE analysis, "
                f"found {len(unique_groups)}"
            )
        group1 = unique_groups[0]
        group2 = unique_groups[1]
        await ctx.info(
            f"No group specified, comparing first two groups "
            f"(alphabetical): {group1} vs {group2}"
        )
    elif group2 is None or group2 == "rest":
        # Compare group1 vs all others combined as "rest"
        group2 = "rest"

    # Validate that specified groups actually exist
    if group1 not in unique_groups:
        raise ParameterError(
            f"Group '{group1}' not found in adata.obs['{group_key}']. "
            f"Available groups: {unique_groups}"
        )
    if group2 != "rest" and group2 not in unique_groups:
        raise ParameterError(
            f"Group '{group2}' not found in adata.obs['{group_key}']. "
            f"Available groups: {unique_groups}"
        )

    # Create pseudobulk aggregation

    # Build aggregation key
    # Cast obs values to str for comparison (group1/group2 are already str
    # from the .astype(str) on unique_groups above)
    obs_group_str = adata.obs[group_key].astype(str)
    if group2 == "rest":
        # Binary comparison: group1 vs rest
        condition = np.where(obs_group_str.to_numpy() == group1, group1, "rest")
        sample_values = adata.obs[sample_key].astype(str).to_numpy()
        raw_X_work = raw_X
    else:
        # Pairwise comparison: filter to only group1 and group2
        mask = obs_group_str.isin([group1, group2]).to_numpy()
        condition = obs_group_str[mask].to_numpy()
        sample_values = adata.obs.loc[mask, sample_key].astype(str).to_numpy()
        raw_X_work = raw_X[mask]

    # Build pseudobulk labels without mutating adata.obs
    condition_values = np.asarray(condition, dtype=str)
    pseudobulk_obs = pd.DataFrame(
        {
            "condition": condition_values,
            "sample": sample_values,
        }
    )

    # Aggregate counts using (sample, condition) tuple to avoid
    # string-concatenation collisions (e.g. sample="A_B" + cond="C"
    # vs sample="A" + cond="B_C" both producing "A_B_C").
    pseudobulk_groups = pseudobulk_obs.groupby(["sample", "condition"]).indices
    pseudobulk_tuples = list(pseudobulk_groups.keys())
    # Create collision-free string IDs using unit separator (\x1f)
    pseudobulk_ids = [f"{s}\x1f{c}" for s, c in pseudobulk_tuples]
    n_samples = len(pseudobulk_ids)

    if n_samples < 4:
        raise DataError(
            f"DESeq2 requires at least 2 samples per group. "
            f"Found only {n_samples} total pseudobulk samples. "
            f"Add more biological replicates or use a different method (wilcoxon)."
        )

    # Build pseudobulk count matrix
    pseudobulk_counts = np.zeros((n_samples, raw_X_work.shape[1]), dtype=np.int64)
    pseudobulk_metadata = []

    for i, pb_tuple in enumerate(pseudobulk_tuples):
        # Get integer positional indices for this pseudobulk group
        int_idx = pseudobulk_groups[pb_tuple]
        # Sum counts (handles both sparse and dense matrices)
        pseudobulk_counts[i] = (
            np.asarray(raw_X_work[int_idx].sum(axis=0)).flatten().astype(np.int64)
        )
        # Get condition from first cell in this group
        first_int_idx = int(int_idx[0])
        pseudobulk_metadata.append(
            {
                "sample_id": pseudobulk_ids[i],
                "condition": condition_values[first_int_idx],
                "sample": sample_values[first_int_idx],
            }
        )

    # Create metadata DataFrame
    metadata_df = pd.DataFrame(pseudobulk_metadata)
    metadata_df = metadata_df.set_index("sample_id")

    # Create count DataFrame
    counts_df = pd.DataFrame(pseudobulk_counts, index=pseudobulk_ids, columns=var_names)

    # Check sample counts per condition
    condition_counts = metadata_df["condition"].value_counts()

    if any(condition_counts < 2):
        raise DataError(
            f"DESeq2 requires at least 2 samples per condition. "
            f"Current counts: {condition_counts.to_dict()}"
        )

    # Import pydeseq2 only after precondition checks so users get actionable
    # validation errors even when optional dependency is unavailable.
    require("pydeseq2", ctx, feature="DESeq2 differential expression")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Run PyDESeq2
    try:
        # Design: ~condition (standard for pseudobulk DE).
        # Samples are biological replicates — including them as
        # covariates would consume degrees of freedom without benefit.
        # For paired designs, users should encode pairing in condition.
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata_df,
            design_factors="condition",
        )

        # Run DESeq2 pipeline
        dds.deseq2()

        # Get results
        stat_res = DeseqStats(dds, contrast=["condition", group1, group2])
        stat_res.summary()

        # Get results DataFrame
        results_df = stat_res.results_df

    except Exception as e:
        raise ProcessingError(
            f"PyDESeq2 analysis failed: {e}\n"
            "This may be due to low sample counts or data issues."
        ) from e

    # Extract top DE genes
    # Sort by adjusted p-value, filter significant genes
    results_df = results_df.dropna(subset=["padj"])
    results_df = results_df.sort_values("padj")

    top_genes = results_df.head(params.n_top_genes).index.tolist()

    if not top_genes:
        raise ProcessingError(
            f"No DE genes found between {group1} and {group2}. "
            "Check sample sizes and expression differences."
        )

    # Get statistics
    n_sig_genes = (results_df["padj"] < 0.05).sum()
    mean_log2fc = results_df.head(params.n_top_genes)["log2FoldChange"].mean()
    median_pvalue = results_df.head(params.n_top_genes)["padj"].median()

    # Store results as DataFrame directly — anndata supports DataFrames in uns
    # and _extract_from_uns() handles them natively (line 257–258).
    # Use parametric uns key so multiple pydeseq2 comparisons coexist.
    analysis_key = _build_de_key("pydeseq2", group1, group2)
    uns_results_key = f"pydeseq2_results_{group1}_vs_{group2}"
    adata.uns[uns_results_key] = results_df.copy()

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name=analysis_key,
        method="pydeseq2",
        parameters={
            "group_key": group_key,
            "sample_key": sample_key,
            "group1": group1,
            "group2": group2,
            "comparison_type": "pseudobulk",
            "n_top_genes": params.n_top_genes,
        },
        results_keys={"uns": [uns_results_key]},
        statistics={
            "method": "pydeseq2",
            "group1": group1,
            "group2": group2,
            "n_pseudobulk_samples": n_samples,
            "n_significant_genes": int(n_sig_genes),
            "mean_log2fc": float(mean_log2fc) if not np.isnan(mean_log2fc) else None,
            "median_padj": (
                float(median_pvalue) if not np.isnan(median_pvalue) else None
            ),
        },
    )

    # Export results to CSV for reproducibility
    export_analysis_result(adata, data_id, analysis_key)

    return DifferentialExpressionResult(
        data_id=data_id,
        comparison=f"{group1} vs {group2}",
        n_genes=len(top_genes),
        top_genes=top_genes,
        statistics={
            "method": "pydeseq2",
            "n_pseudobulk_samples": n_samples,
            "n_significant_genes": int(n_sig_genes),
            "mean_log2fc": float(mean_log2fc) if not np.isnan(mean_log2fc) else None,
            "median_padj": (
                float(median_pvalue) if not np.isnan(median_pvalue) else None
            ),
        },
    )
