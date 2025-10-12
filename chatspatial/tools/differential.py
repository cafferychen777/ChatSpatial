"""
Differential expression analysis tools for spatial transcriptomics data.
"""

from typing import Any, Dict, Optional

import numpy as np
import scanpy as sc
from mcp.server.fastmcp import Context

from ..models.analysis import DifferentialExpressionResult


async def differential_expression(
    data_id: str,
    data_store: Dict[str, Any],
    group_key: str,
    group1: Optional[str] = None,
    group2: Optional[str] = None,
    method: str = "wilcoxon",
    n_top_genes: int = 50,
    context: Optional[Context] = None,
) -> DifferentialExpressionResult:
    """Perform differential expression analysis

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        group_key: Key in adata.obs for grouping cells
        group1: First group for comparison (if None, find markers for all groups)
        group2: Second group for comparison (if None, compare against rest)
        n_top_genes: Number of top differentially expressed genes to return
        method: Statistical method for DE analysis
        context: MCP context

    Returns:
        Differential expression analysis result
    """
    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    adata = data_store[data_id]["adata"]

    # Check if the group_key exists in adata.obs
    if group_key not in adata.obs.columns:
        raise ValueError(f"Group key '{group_key}' not found in adata.obs")

    # IMPORTANT: Handle float16 data type (numba doesn't support float16)
    # Convert to float32 if needed for differential expression analysis
    if hasattr(adata.X, 'dtype') and adata.X.dtype == np.float16:
        if context:
            await context.info(
                "⚙️ Converting data from float16 to float32 for compatibility with rank_genes_groups"
            )
        # Create a copy to avoid modifying the original data
        adata = adata.copy()
        adata.X = adata.X.astype(np.float32)

    # If group1 is None, find markers for all groups
    if group1 is None:
        if context:
            await context.info(f"Finding marker genes for all groups in '{group_key}'")

        # Filter out groups with too few cells
        group_sizes = adata.obs[group_key].value_counts()
        min_cells = 3  # Minimum for Wilcoxon test
        valid_groups = group_sizes[group_sizes >= min_cells]
        skipped_groups = group_sizes[group_sizes < min_cells]

        # Warn about skipped groups
        if len(skipped_groups) > 0:
            if context:
                skipped_list = "\n".join(
                    [f"  • {g}: {n} cell(s)" for g, n in skipped_groups.items()]
                )
                await context.warning(
                    f"⚠️ Skipped {len(skipped_groups)} group(s) with <{min_cells} cells:\n{skipped_list}"
                )

        # Check if any valid groups remain
        if len(valid_groups) == 0:
            all_sizes = "\n".join(
                [f"  • {g}: {n} cell(s)" for g, n in group_sizes.items()]
            )
            raise ValueError(
                f"❌ All groups have <{min_cells} cells. Cannot perform {method} test.\n\n"
                f"📊 Group sizes:\n{all_sizes}\n\n"
                f"💡 Try: find_markers(group_key='leiden') or merge small groups"
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
    if context:
        await context.info(
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
    if context:
        await context.info(f"Running rank_genes_groups with method '{method}'")

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
    if context:
        await context.info("Extracting differential expression results")

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

    # Get log fold changes and p-values if available
    log2fc_values = []
    pvals = []
    if hasattr(temp_adata, "uns") and "rank_genes_groups" in temp_adata.uns:
        if (
            "logfoldchanges" in temp_adata.uns["rank_genes_groups"]
            and group1
            in temp_adata.uns["rank_genes_groups"]["logfoldchanges"].dtype.names
        ):
            log2fc_values = list(
                temp_adata.uns["rank_genes_groups"]["logfoldchanges"][group1][
                    :n_top_genes
                ]
            )

        if (
            "pvals_adj" in temp_adata.uns["rank_genes_groups"]
            and group1 in temp_adata.uns["rank_genes_groups"]["pvals_adj"].dtype.names
        ):
            pvals = list(
                temp_adata.uns["rank_genes_groups"]["pvals_adj"][group1][:n_top_genes]
            )

    # Calculate mean log2fc and median p-value
    mean_log2fc = np.mean(log2fc_values) if log2fc_values else None
    median_pvalue = np.median(pvals) if pvals else None

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

    return DifferentialExpressionResult(
        data_id=data_id,
        comparison=comparison,
        n_genes=len(top_genes),
        top_genes=top_genes,
        statistics=statistics,
    )
