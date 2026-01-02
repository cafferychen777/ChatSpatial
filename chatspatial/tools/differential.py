"""
Differential expression analysis tools for spatial transcriptomics data.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from ..models.analysis import DifferentialExpressionResult
from ..models.data import DifferentialExpressionParameters
from ..spatial_mcp_adapter import ToolContext
from ..utils import validate_obs_column
from ..utils.adata_utils import store_analysis_metadata, to_dense
from ..utils.dependency_manager import require
from ..utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)


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
    validate_obs_column(adata, group_key, "Group key")

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
        adata_filtered = adata[adata.obs[group_key].isin(valid_groups.index)].copy()

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

        # Store metadata for scientific provenance tracking
        store_analysis_metadata(
            adata,
            analysis_name="differential_expression",
            method=method,
            parameters={
                "group_key": group_key,
                "comparison_type": "all_groups",
                "n_top_genes": n_top_genes,
            },
            results_keys={"uns": ["rank_genes_groups"]},
            statistics={
                "method": method,
                "n_groups": len(groups),
                "groups": list(map(str, groups)),
                "n_cells_analyzed": adata_filtered.n_obs,
                "n_genes_analyzed": adata_filtered.n_vars,
            },
        )

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
        temp_adata = adata.copy()
    else:
        # Create a temporary copy of the AnnData object with only the two groups
        temp_adata = adata[adata.obs[group_key].isin([group1, group2])].copy()

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
    if hasattr(temp_adata, "uns") and "rank_genes_groups" in temp_adata.uns:
        if "names" in temp_adata.uns["rank_genes_groups"]:
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
    if hasattr(temp_adata, "uns") and "rank_genes_groups" in temp_adata.uns:
        if (
            "pvals_adj" in temp_adata.uns["rank_genes_groups"]
            and group1 in temp_adata.uns["rank_genes_groups"]["pvals_adj"].dtype.names
        ):
            pvals = list(
                temp_adata.uns["rank_genes_groups"]["pvals_adj"][group1][:n_top_genes]
            )

    # Calculate TRUE fold change from raw counts (Bug #3 Fix)
    # Issue: scanpy's logfoldchanges uses mean(log(counts)) which is mathematically incorrect
    # Solution: Calculate log(mean(counts1) / mean(counts2)) from raw data

    # Check if raw count data is available
    if adata.raw is None:
        raise DataNotFoundError(
            "Raw count data (adata.raw) required for fold change calculation. "
            "Run preprocess_data() first to preserve raw counts."
        )

    # Get raw count data
    raw_adata = adata.raw
    log2fc_values = []

    # Create masks for the two groups
    if use_rest_as_reference:
        group1_mask = adata.obs[group_key] == group1
        group2_mask = ~group1_mask
    else:
        group1_mask = adata.obs[group_key] == group1
        group2_mask = adata.obs[group_key] == group2

    # CRITICAL: Normalize by library size to avoid composition bias
    # Library size = total UMI counts per spot
    if hasattr(raw_adata.X, "toarray"):
        lib_sizes = np.array(raw_adata.X.sum(axis=1)).flatten()
    else:
        lib_sizes = raw_adata.X.sum(axis=1).flatten()

    median_lib_size = float(np.median(lib_sizes))

    # Calculate fold change for each top gene
    for gene in top_genes:
        if gene in raw_adata.var_names:
            gene_idx = raw_adata.var_names.get_loc(gene)

            # Get raw counts for this gene
            gene_raw_counts = to_dense(raw_adata.X[:, gene_idx]).flatten()

            # Normalize by library size (CPM-like normalization)
            # normalized_counts = raw_counts * (median_lib_size / spot_lib_size)
            gene_norm_counts = gene_raw_counts * (median_lib_size / lib_sizes)

            # Calculate mean normalized counts for each group
            mean_group1 = float(gene_norm_counts[group1_mask].mean())
            mean_group2 = float(gene_norm_counts[group2_mask].mean())

            # Calculate true log2 fold change from normalized counts
            # Add user-configurable pseudocount to avoid log(0)
            true_log2fc = np.log2(
                (mean_group1 + pseudocount) / (mean_group2 + pseudocount)
            )
            log2fc_values.append(float(true_log2fc))
        else:
            # Gene not in raw data (should not happen, but handle gracefully)
            await ctx.warning(
                f"Gene {gene} not found in raw data, skipping fold change calculation"
            )
            log2fc_values.append(None)

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

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name="differential_expression",
        method=method,
        parameters={
            "group_key": group_key,
            "group1": group1,
            "group2": group2,
            "comparison_type": "specific_groups",
            "n_top_genes": n_top_genes,
            "pseudocount": pseudocount,  # Track for reproducibility
        },
        results_keys={"uns": ["rank_genes_groups"]},
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

    # Import pydeseq2 (require() raises ImportError if not available)
    require("pydeseq2", ctx, feature="DESeq2 differential expression")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Get data
    adata = await ctx.get_adata(data_id)

    # Validate columns
    validate_obs_column(adata, params.group_key, "Group key")
    validate_obs_column(adata, params.sample_key, "Sample key")

    # Get raw counts (required for DESeq2)
    if adata.raw is not None:
        raw_X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        raw_X = adata.X
        var_names = adata.var_names

    # Convert to dense if sparse
    if sparse.issparse(raw_X):
        raw_X = raw_X.toarray()

    # Validate counts are integers (DESeq2 requirement)
    if not np.allclose(raw_X, raw_X.astype(int)):
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
    unique_groups = adata.obs[group_key].unique()
    if group1 is None:
        if len(unique_groups) < 2:
            raise DataError(
                f"Need at least 2 groups for DE analysis, found {len(unique_groups)}"
            )
        group1 = str(unique_groups[0])
        group2 = str(unique_groups[1])
        await ctx.info(
            f"No group specified, comparing first two groups: {group1} vs {group2}"
        )
    elif group2 is None or group2 == "rest":
        # Compare group1 vs all others combined as "rest"
        group2 = "rest"

    # Create pseudobulk aggregation
    await ctx.info(f"Creating pseudobulk samples by {sample_key} and {group_key}...")

    # Build aggregation key
    if group2 == "rest":
        # Binary comparison: group1 vs rest
        condition = adata.obs[group_key].apply(
            lambda x: group1 if x == group1 else "rest"
        )
    else:
        # Pairwise comparison: filter to only group1 and group2
        mask = adata.obs[group_key].isin([group1, group2])
        adata = adata[mask].copy()
        raw_X = raw_X[mask.values]
        condition = adata.obs[group_key].astype(str)

    # Create pseudobulk by aggregating (summing) counts per sample+condition
    adata.obs["_de_condition"] = condition.values
    adata.obs["_pseudobulk_id"] = (
        adata.obs[sample_key].astype(str) + "_" + adata.obs["_de_condition"].astype(str)
    )

    # Aggregate counts
    pseudobulk_groups = adata.obs.groupby("_pseudobulk_id")
    pseudobulk_ids = list(pseudobulk_groups.groups.keys())
    n_samples = len(pseudobulk_ids)

    await ctx.info(f"Aggregated into {n_samples} pseudobulk samples")

    if n_samples < 4:
        raise DataError(
            f"DESeq2 requires at least 2 samples per group. "
            f"Found only {n_samples} total pseudobulk samples. "
            f"Add more biological replicates or use a different method (wilcoxon)."
        )

    # Build pseudobulk count matrix
    pseudobulk_counts = np.zeros((n_samples, raw_X.shape[1]), dtype=np.int64)
    pseudobulk_metadata = []

    for i, pb_id in enumerate(pseudobulk_ids):
        # Get indices for this pseudobulk group
        group_labels = pseudobulk_groups.groups[pb_id]
        # Convert pandas Index to integer positional indices for numpy array indexing
        int_idx = adata.obs.index.get_indexer(group_labels)
        pseudobulk_counts[i] = raw_X[int_idx].sum(axis=0).astype(np.int64)
        # Get condition from first cell in this group
        first_int_idx = int_idx[0]
        pseudobulk_metadata.append(
            {
                "sample_id": pb_id,
                "condition": adata.obs.iloc[first_int_idx]["_de_condition"],
                "sample": adata.obs.iloc[first_int_idx][sample_key],
            }
        )

    # Create metadata DataFrame
    metadata_df = pd.DataFrame(pseudobulk_metadata)
    metadata_df = metadata_df.set_index("sample_id")

    # Create count DataFrame
    counts_df = pd.DataFrame(pseudobulk_counts, index=pseudobulk_ids, columns=var_names)

    # Check sample counts per condition
    condition_counts = metadata_df["condition"].value_counts()
    await ctx.info(f"Samples per condition: {condition_counts.to_dict()}")

    if any(condition_counts < 2):
        raise DataError(
            f"DESeq2 requires at least 2 samples per condition. "
            f"Current counts: {condition_counts.to_dict()}"
        )

    # Run PyDESeq2
    await ctx.info("Running PyDESeq2 differential expression analysis...")

    try:
        # Create DESeq2 dataset
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
            f"PyDESeq2 analysis failed: {str(e)}\n"
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

    # Store results in adata.uns for persistence
    adata.uns["pydeseq2_results"] = {
        "results_df": results_df.to_dict(),
        "comparison": f"{group1} vs {group2}",
        "n_samples": n_samples,
    }

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name="differential_expression",
        method="pydeseq2",
        parameters={
            "group_key": group_key,
            "sample_key": sample_key,
            "group1": group1,
            "group2": group2,
            "comparison_type": "pseudobulk",
            "n_top_genes": params.n_top_genes,
        },
        results_keys={"uns": ["pydeseq2_results"]},
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
