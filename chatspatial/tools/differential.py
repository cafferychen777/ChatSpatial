"""
Differential expression analysis tools for spatial transcriptomics data.
"""

from typing import Optional

import numpy as np
import scanpy as sc

from ..models.analysis import DifferentialExpressionResult
from ..spatial_mcp_adapter import ToolContext
from ..utils import validate_obs_column
from ..utils.adata_utils import store_analysis_metadata


async def differential_expression(
    data_id: str,
    ctx: ToolContext,
    group_key: str,
    group1: Optional[str] = None,
    group2: Optional[str] = None,
    method: str = "wilcoxon",
    n_top_genes: int = 50,
    pseudocount: float = 1.0,
    min_cells: int = 3,
) -> DifferentialExpressionResult:
    """Perform differential expression analysis

    Args:
        data_id: Dataset ID
        ctx: Tool context for data access and logging
        group_key: Key in adata.obs for grouping cells
        group1: First group for comparison (if None, find markers for all groups)
        group2: Second group for comparison (if None, compare against rest)
        n_top_genes: Number of top differentially expressed genes to return
        method: Statistical method for DE analysis
        pseudocount: Pseudocount added before log2 fold change calculation.
                    Default: 1.0 (standard practice, prevents log(0)).
                    Lower values (0.1-0.5): More sensitive to low-expression changes.
                    Higher values (1-10): More stable for sparse/noisy data.
        min_cells: Minimum number of cells per group for statistical testing.
                  Default: 3 (minimum required for Wilcoxon test).
                  Increase to 10-30 for more robust results.
                  Groups with fewer cells are skipped with a warning.

    Returns:
        Differential expression analysis result
    """
    # Get AnnData directly via ToolContext (no redundant dict wrapping)
    adata = await ctx.get_adata(data_id)

    # Check if the group_key exists in adata.obs
    validate_obs_column(adata, group_key, "Group key")

    # IMPORTANT: Handle float16 data type (numba doesn't support float16)
    # Convert to float32 if needed for differential expression analysis
    if hasattr(adata.X, "dtype") and adata.X.dtype == np.float16:
        await ctx.info(
            "Converting data from float16 to float32 for compatibility with rank_genes_groups"
        )
        # Create a copy to avoid modifying the original data
        adata = adata.copy()
        adata.X = adata.X.astype(np.float32)

    # If group1 is None, find markers for all groups
    if group1 is None:
        await ctx.info(f"Finding marker genes for all groups in '{group_key}'")

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
            raise ValueError(
                f"All groups have <{min_cells} cells. Cannot perform {method} test.\n\n"
                f"Group sizes:\n{all_sizes}\n\n"
                f"Try: find_markers(group_key='leiden') or merge small groups"
            )

        # Filter data to only include valid groups
        adata_filtered = adata[adata.obs[group_key].isin(valid_groups.index)].copy()

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
    await ctx.info(
        f"Performing differential expression analysis between {group1} and {group2}"
    )

    # Check if the groups exist in the group_key
    if group1 not in adata.obs[group_key].values:
        raise ValueError(f"Group '{group1}' not found in adata.obs['{group_key}']")

    # Special case for 'rest' as group2 or if group2 is None
    use_rest_as_reference = False
    if group2 is None or group2 == "rest":
        use_rest_as_reference = True
        group2 = "rest"  # Set it explicitly for display purposes
    elif group2 not in adata.obs[group_key].values:
        raise ValueError(f"Group '{group2}' not found in adata.obs['{group_key}']")

    # Perform differential expression analysis
    await ctx.info(f"Running rank_genes_groups with method '{method}'")

    # Prepare the AnnData object for analysis
    if use_rest_as_reference:
        # Use the full AnnData object when comparing with 'rest'
        temp_adata = adata.copy()
    else:
        # Create a temporary copy of the AnnData object with only the two groups
        temp_adata = adata[adata.obs[group_key].isin([group1, group2])].copy()

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
    await ctx.info("Extracting differential expression results")

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
        raise RuntimeError(
            f"Failed to identify differentially expressed genes for comparison between {group1} and {group2}. "
            f"This could be due to: 1) Insufficient statistical power, 2) No meaningful expression differences, "
            f"or 3) Issues with the analysis method. Please check data quality, sample sizes, or try different parameters."
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
        raise ValueError(
            "CRITICAL: Raw count data required for accurate fold change calculation\n\n"
            "Differential expression analysis requires raw counts (adata.raw) to calculate "
            "scientifically accurate fold change values. Your data lacks this.\n\n"
            "WHY THIS MATTERS:\n"
            "• Fold change from log-normalized data produces mathematically incorrect results\n"
            "• Can lead to extreme overestimation (e.g., 1000× instead of 3×)\n"
            "• Results would not be scientifically valid for publication\n\n"
            "SOLUTION:\n"
            "1. Reload your data and run preprocessing again:\n"
            "   preprocess_data(data_id='your_data', params={'normalization': 'log'})\n"
            "2. Preprocessing automatically saves raw counts to adata.raw\n"
            "3. Then run differential expression analysis\n\n"
            "WARNING: We refuse to return incorrect fold change values. "
            "Scientific integrity requires accurate calculations."
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

    await ctx.info(
        f"Calculating fold changes with library size normalization "
        f"(median library size: {median_lib_size:.0f} UMIs)"
    )
    await ctx.info(
        f"Group 1: {np.sum(group1_mask)} cells, Group 2: {np.sum(group2_mask)} cells"
    )

    # Calculate fold change for each top gene
    for gene in top_genes:
        if gene in raw_adata.var_names:
            gene_idx = raw_adata.var_names.get_loc(gene)

            # Get raw counts for this gene
            if hasattr(raw_adata.X, "toarray"):
                gene_raw_counts = raw_adata.X[:, gene_idx].toarray().flatten()
            else:
                gene_raw_counts = raw_adata.X[:, gene_idx].flatten()

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
            f"EXTREME FOLD CHANGE DETECTED: mean log2FC = {mean_log2fc:.2f}\n"
            f"   - This indicates fold change > 2^10 = 1024x\n"
            f"   - Biologically plausible range: 2-32x (log2FC = 1-5)\n"
            f"   - Possible causes:\n"
            f"     - Very sparse gene expression in one group\n"
            f"     - Low cell counts in comparison groups\n"
            f"     - Data quality issues\n"
            f"   - Consider filtering genes with low expression before DE analysis"
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
