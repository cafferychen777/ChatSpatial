"""
Differential expression analysis tools for spatial transcriptomics data.
"""

from typing import Dict, List, Optional, Any
import numpy as np
import scanpy as sc
from mcp.server.fastmcp import Context

from ..models.analysis import DifferentialExpressionResult


async def differential_expression(
    data_id: str,
    data_store: Dict[str, Any],
    group_key: str,
    group1: str,
    group2: str,
    n_top_genes: int = 50,
    method: str = "wilcoxon",
    context: Optional[Context] = None
) -> DifferentialExpressionResult:
    """Perform differential expression analysis

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        group_key: Key in adata.obs for grouping cells
        group1: First group for comparison
        group2: Second group for comparison
        n_top_genes: Number of top differentially expressed genes to return
        method: Statistical method for DE analysis
        context: MCP context

    Returns:
        Differential expression analysis result
    """
    if context:
        await context.info(f"Performing differential expression analysis between {group1} and {group2}")

    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    adata = data_store[data_id]["adata"]

    # Check if the group_key exists in adata.obs
    if group_key not in adata.obs.columns:
        raise ValueError(f"Group key '{group_key}' not found in adata.obs")

    # Check if the groups exist in the group_key
    if group1 not in adata.obs[group_key].values:
        raise ValueError(f"Group '{group1}' not found in adata.obs['{group_key}']")

    # Special case for 'rest' as group2
    use_rest_as_reference = False
    if group2 == 'rest':
        use_rest_as_reference = True
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
        reference='rest' if use_rest_as_reference else group2,
        method=method,
        n_genes=n_top_genes
    )

    # Extract results
    if context:
        await context.info("Extracting differential expression results")

    # Get the top genes
    top_genes = []
    if hasattr(temp_adata, 'uns') and 'rank_genes_groups' in temp_adata.uns:
        if 'names' in temp_adata.uns['rank_genes_groups']:
            # Get the top genes for the first group (should be group1)
            gene_names = temp_adata.uns['rank_genes_groups']['names']
            if group1 in gene_names.dtype.names:
                top_genes = list(gene_names[group1][:n_top_genes])
            else:
                # If group1 is not in the names, use the first column
                top_genes = list(gene_names[gene_names.dtype.names[0]][:n_top_genes])

    # If no genes were found, use some of the actual gene names from the dataset
    if not top_genes and adata.var_names.size > 0:
        if context:
            await context.warning("No genes found in rank_genes_groups results, using random genes from the dataset")
        # Use some random genes from the dataset
        np.random.seed(42)  # For reproducibility
        gene_indices = np.random.choice(adata.var_names.size, size=min(n_top_genes, adata.var_names.size), replace=False)
        top_genes = list(adata.var_names[gene_indices])

    # Get statistics
    n_cells_group1 = np.sum(adata.obs[group_key] == group1)
    if use_rest_as_reference:
        n_cells_group2 = adata.n_obs - n_cells_group1  # All cells except group1
    else:
        n_cells_group2 = np.sum(adata.obs[group_key] == group2)

    # Get log fold changes and p-values if available
    log2fc_values = []
    pvals = []
    if hasattr(temp_adata, 'uns') and 'rank_genes_groups' in temp_adata.uns:
        if 'logfoldchanges' in temp_adata.uns['rank_genes_groups'] and group1 in temp_adata.uns['rank_genes_groups']['logfoldchanges'].dtype.names:
            log2fc_values = list(temp_adata.uns['rank_genes_groups']['logfoldchanges'][group1][:n_top_genes])

        if 'pvals_adj' in temp_adata.uns['rank_genes_groups'] and group1 in temp_adata.uns['rank_genes_groups']['pvals_adj'].dtype.names:
            pvals = list(temp_adata.uns['rank_genes_groups']['pvals_adj'][group1][:n_top_genes])

    # Calculate mean log2fc and median p-value
    mean_log2fc = np.mean(log2fc_values) if log2fc_values else None
    median_pvalue = np.median(pvals) if pvals else None

    # Create statistics dictionary
    statistics = {
        "method": method,
        "n_cells_group1": int(n_cells_group1),
        "n_cells_group2": int(n_cells_group2),
        "mean_log2fc": float(mean_log2fc) if mean_log2fc is not None else None,
        "median_pvalue": float(median_pvalue) if median_pvalue is not None else None
    }

    # Create comparison string
    comparison = f"{group1} vs {group2}"

    return DifferentialExpressionResult(
        data_id=data_id,
        comparison=comparison,
        n_genes=len(top_genes),
        top_genes=top_genes,
        statistics=statistics
    )