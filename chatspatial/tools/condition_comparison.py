"""
Multi-sample condition comparison analysis for spatial transcriptomics data.

This module implements pseudobulk differential expression analysis for comparing
experimental conditions (e.g., Treatment vs Control) across biological samples.
"""

from typing import Optional, Union

import numpy as np
import pandas as pd
from scipy import sparse

from ..models.analysis import (
    CellTypeComparisonResult,
    ConditionComparisonResult,
    DEGene,
)
from ..models.data import ConditionComparisonParameters
from ..spatial_mcp_adapter import ToolContext
from ..utils import validate_obs_column
from ..utils.adata_utils import check_is_integer_counts, store_analysis_metadata
from ..utils.dependency_manager import require
from ..utils.exceptions import DataError, ParameterError, ProcessingError


async def compare_conditions(
    data_id: str,
    ctx: ToolContext,
    params: ConditionComparisonParameters,
) -> ConditionComparisonResult:
    """Compare experimental conditions across multiple biological samples.

    This function performs pseudobulk differential expression analysis using DESeq2.
    It aggregates cells by sample, then compares conditions (e.g., Treatment vs Control).

    Optionally, analysis can be stratified by cell type to identify cell type-specific
    condition effects.

    Args:
        data_id: Dataset ID
        ctx: Tool context for data access and logging
        params: Condition comparison parameters

    Returns:
        ConditionComparisonResult with differential expression results

    Example:
        # Global comparison (all cells)
        compare_conditions(
            data_id="data1",
            condition_key="treatment",
            condition1="Drug",
            condition2="Control",
            sample_key="patient_id"
        )

        # Cell type stratified comparison
        compare_conditions(
            data_id="data1",
            condition_key="treatment",
            condition1="Drug",
            condition2="Control",
            sample_key="patient_id",
            cell_type_key="cell_type"
        )
    """
    # Check pydeseq2 availability early (required for pseudobulk analysis)
    require("pydeseq2", ctx, feature="Condition comparison with DESeq2")

    # Get data
    adata = await ctx.get_adata(data_id)

    # Validate required columns
    validate_obs_column(adata, params.condition_key, "Condition")
    validate_obs_column(adata, params.sample_key, "Sample")
    if params.cell_type_key is not None:
        validate_obs_column(adata, params.cell_type_key, "Cell type")

    # Validate conditions exist
    unique_conditions = adata.obs[params.condition_key].unique()
    if params.condition1 not in unique_conditions:
        raise ParameterError(
            f"Condition '{params.condition1}' not found in '{params.condition_key}'.\n"
            f"Available conditions: {list(unique_conditions)}"
        )
    if params.condition2 not in unique_conditions:
        raise ParameterError(
            f"Condition '{params.condition2}' not found in '{params.condition_key}'.\n"
            f"Available conditions: {list(unique_conditions)}"
        )

    # Filter to only the two conditions of interest
    mask = adata.obs[params.condition_key].isin([params.condition1, params.condition2])
    adata_filtered = adata[mask].copy()

    await ctx.info(
        f"Comparing {params.condition1} vs {params.condition2}: "
        f"{adata_filtered.n_obs} cells from {adata_filtered.obs[params.sample_key].nunique()} samples"
    )

    # Get raw counts (required for DESeq2)
    raw_X, var_names = _get_raw_counts(adata_filtered)

    # Validate counts are integers (handles sparse matrices)
    is_int, _, _ = check_is_integer_counts(raw_X)
    if not is_int:
        await ctx.warning(
            "Data appears to be normalized. DESeq2 requires raw integer counts. "
            "Results may be inaccurate. Consider using adata.raw."
        )

    # Count samples per condition
    sample_condition_map = adata_filtered.obs.groupby(params.sample_key)[
        params.condition_key
    ].first()
    n_samples_cond1 = (sample_condition_map == params.condition1).sum()
    n_samples_cond2 = (sample_condition_map == params.condition2).sum()

    await ctx.info(
        f"Sample distribution: {params.condition1}={n_samples_cond1}, "
        f"{params.condition2}={n_samples_cond2}"
    )

    # Check minimum samples requirement
    if n_samples_cond1 < params.min_samples_per_condition:
        raise DataError(
            f"Insufficient samples for {params.condition1}: {n_samples_cond1} "
            f"(minimum: {params.min_samples_per_condition})"
        )
    if n_samples_cond2 < params.min_samples_per_condition:
        raise DataError(
            f"Insufficient samples for {params.condition2}: {n_samples_cond2} "
            f"(minimum: {params.min_samples_per_condition})"
        )

    # Determine analysis mode
    if params.cell_type_key is None:
        # Global analysis (all cells together)
        result = await _run_global_comparison(
            adata_filtered, raw_X, var_names, ctx, params
        )
    else:
        # Cell type stratified analysis
        result = await _run_stratified_comparison(
            adata_filtered, raw_X, var_names, ctx, params
        )

    # Update result with common fields
    result.data_id = data_id
    result.n_samples_condition1 = int(n_samples_cond1)
    result.n_samples_condition2 = int(n_samples_cond2)

    # Store results in adata
    results_key = f"condition_comparison_{params.condition1}_vs_{params.condition2}"
    adata.uns[results_key] = {
        "comparison": result.comparison,
        "method": result.method,
        "statistics": result.statistics,
    }

    # Store metadata for provenance
    store_analysis_metadata(
        adata,
        analysis_name="condition_comparison",
        method="pseudobulk_deseq2",
        parameters={
            "condition_key": params.condition_key,
            "condition1": params.condition1,
            "condition2": params.condition2,
            "sample_key": params.sample_key,
            "cell_type_key": params.cell_type_key,
        },
        results_keys={"uns": [results_key]},
        statistics=result.statistics,
    )

    result.results_key = results_key
    return result


def _get_raw_counts(
    adata,
) -> tuple[Union[np.ndarray, sparse.spmatrix], pd.Index]:
    """Extract raw count matrix from AnnData.

    Returns sparse matrix directly to avoid memory explosion.
    Downstream code handles both sparse and dense via np.asarray().

    Args:
        adata: AnnData object

    Returns:
        Tuple of (count_matrix, var_names) - matrix may be sparse or dense
    """
    if adata.raw is not None:
        return adata.raw.X, adata.raw.var_names
    return adata.X, adata.var_names


def _create_pseudobulk(
    adata,
    raw_X: Union[np.ndarray, sparse.spmatrix],
    var_names: pd.Index,
    sample_key: str,
    condition_key: str,
    cell_type: Optional[str] = None,
    cell_type_key: Optional[str] = None,
    min_cells_per_sample: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Create pseudobulk count matrix by aggregating cells per sample.

    Args:
        adata: AnnData object
        raw_X: Raw count matrix
        var_names: Gene names
        sample_key: Column for sample identification
        condition_key: Column for condition
        cell_type: Specific cell type to filter (optional)
        cell_type_key: Column for cell type (required if cell_type is provided)
        min_cells_per_sample: Minimum cells required per sample

    Returns:
        Tuple of (counts_df, metadata_df, cell_counts)
    """
    # Filter to specific cell type if provided
    if cell_type is not None and cell_type_key is not None:
        mask = adata.obs[cell_type_key] == cell_type
        adata = adata[mask].copy()
        raw_X = raw_X[mask.values]

    # Group by sample
    sample_groups = adata.obs.groupby(sample_key)

    pseudobulk_data = []
    metadata_list = []
    cell_counts = {}

    for sample_id, group in sample_groups:
        n_cells = len(group)
        if n_cells < min_cells_per_sample:
            continue

        # Get integer indices for this sample
        int_idx = adata.obs.index.get_indexer(group.index)

        # Sum counts (handles both sparse and dense matrices)
        sample_counts = (
            np.asarray(raw_X[int_idx].sum(axis=0)).flatten().astype(np.int64)
        )

        # Get condition for this sample
        condition = group[condition_key].iloc[0]

        pseudobulk_data.append(sample_counts)
        metadata_list.append(
            {
                "sample_id": sample_id,
                "condition": condition,
            }
        )
        cell_counts[str(sample_id)] = n_cells

    if len(pseudobulk_data) == 0:
        raise DataError(
            f"No samples have >= {min_cells_per_sample} cells. "
            "Try lowering min_cells_per_sample."
        )

    # Create DataFrames
    sample_ids = [m["sample_id"] for m in metadata_list]
    counts_df = pd.DataFrame(
        np.array(pseudobulk_data),
        index=sample_ids,
        columns=var_names,
    )
    metadata_df = pd.DataFrame(metadata_list).set_index("sample_id")

    return counts_df, metadata_df, cell_counts


def _run_deseq2(
    counts_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    condition1: str,
    condition2: str,
    n_top_genes: int = 50,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 0.0,
) -> tuple[list[DEGene], list[DEGene], int, pd.DataFrame]:
    """Run DESeq2 analysis on pseudobulk data.

    Args:
        counts_df: Pseudobulk count matrix
        metadata_df: Sample metadata with condition column
        condition1: First condition (experimental)
        condition2: Second condition (reference/control)
        n_top_genes: Number of top genes to return
        padj_threshold: Adjusted p-value threshold for significance
        log2fc_threshold: Log2 fold change threshold

    Returns:
        Tuple of (top_upregulated, top_downregulated, n_significant, results_df)
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Create DESeq2 dataset
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        design_factors="condition",
    )

    # Run DESeq2 pipeline
    dds.deseq2()

    # Get results (condition1 vs condition2)
    stat_res = DeseqStats(dds, contrast=["condition", condition1, condition2])
    stat_res.summary()

    results_df = stat_res.results_df.dropna(subset=["padj"])

    # Filter by thresholds
    sig_mask = (results_df["padj"] < padj_threshold) & (
        np.abs(results_df["log2FoldChange"]) > log2fc_threshold
    )
    n_significant = sig_mask.sum()

    # Separate upregulated and downregulated
    upregulated = results_df[
        (results_df["padj"] < padj_threshold)
        & (results_df["log2FoldChange"] > log2fc_threshold)
    ].sort_values("padj")

    downregulated = results_df[
        (results_df["padj"] < padj_threshold)
        & (results_df["log2FoldChange"] < -log2fc_threshold)
    ].sort_values("padj")

    # Convert to DEGene objects
    def df_to_degenes(df: pd.DataFrame, n: int) -> list[DEGene]:
        genes = []
        for gene_name, row in df.head(n).iterrows():
            genes.append(
                DEGene(
                    gene=str(gene_name),
                    log2fc=float(row["log2FoldChange"]),
                    pvalue=float(row["pvalue"]),
                    padj=float(row["padj"]),
                )
            )
        return genes

    top_up = df_to_degenes(upregulated, n_top_genes)
    top_down = df_to_degenes(downregulated, n_top_genes)

    return top_up, top_down, int(n_significant), results_df


async def _run_global_comparison(
    adata,
    raw_X: Union[np.ndarray, sparse.spmatrix],
    var_names: pd.Index,
    ctx: ToolContext,
    params: ConditionComparisonParameters,
) -> ConditionComparisonResult:
    """Run global comparison (all cells, no cell type stratification).

    Args:
        adata: Filtered AnnData object
        raw_X: Raw count matrix
        var_names: Gene names
        ctx: Tool context
        params: Comparison parameters

    Returns:
        ConditionComparisonResult
    """
    await ctx.info("Running global pseudobulk analysis (all cells)...")

    # Create pseudobulk
    counts_df, metadata_df, cell_counts = _create_pseudobulk(
        adata,
        raw_X,
        var_names,
        sample_key=params.sample_key,
        condition_key=params.condition_key,
        min_cells_per_sample=params.min_cells_per_sample,
    )

    # Check sample distribution
    cond_counts = metadata_df["condition"].value_counts()
    n_cond1 = cond_counts.get(params.condition1, 0)
    n_cond2 = cond_counts.get(params.condition2, 0)

    if n_cond1 < 2 or n_cond2 < 2:
        raise DataError(
            f"DESeq2 requires at least 2 samples per condition. "
            f"Found: {params.condition1}={n_cond1}, {params.condition2}={n_cond2}"
        )

    await ctx.info(
        f"Created {len(counts_df)} pseudobulk samples "
        f"({params.condition1}={n_cond1}, {params.condition2}={n_cond2})"
    )

    # Run DESeq2
    try:
        top_up, top_down, n_significant, results_df = _run_deseq2(
            counts_df,
            metadata_df,
            condition1=params.condition1,
            condition2=params.condition2,
            n_top_genes=params.n_top_genes,
            padj_threshold=params.padj_threshold,
            log2fc_threshold=params.log2fc_threshold,
        )
    except Exception as e:
        raise ProcessingError(f"DESeq2 analysis failed: {e}") from e

    await ctx.info(f"Found {n_significant} significant DE genes")

    # Build result
    comparison = f"{params.condition1} vs {params.condition2}"

    return ConditionComparisonResult(
        data_id="",  # Will be filled by caller
        method="pseudobulk",
        comparison=comparison,
        condition_key=params.condition_key,
        condition1=params.condition1,
        condition2=params.condition2,
        sample_key=params.sample_key,
        cell_type_key=None,
        n_samples_condition1=0,  # Will be filled by caller
        n_samples_condition2=0,  # Will be filled by caller
        global_n_significant=n_significant,
        global_top_upregulated=top_up,
        global_top_downregulated=top_down,
        cell_type_results=None,
        results_key="",  # Will be filled by caller
        statistics={
            "analysis_type": "global",
            "n_pseudobulk_samples": len(counts_df),
            "n_significant_genes": n_significant,
            "n_upregulated": len([g for g in top_up if g.padj < params.padj_threshold]),
            "n_downregulated": len(
                [g for g in top_down if g.padj < params.padj_threshold]
            ),
        },
    )


async def _run_stratified_comparison(
    adata,
    raw_X: Union[np.ndarray, sparse.spmatrix],
    var_names: pd.Index,
    ctx: ToolContext,
    params: ConditionComparisonParameters,
) -> ConditionComparisonResult:
    """Run cell type stratified comparison.

    Args:
        adata: Filtered AnnData object
        raw_X: Raw count matrix
        var_names: Gene names
        ctx: Tool context
        params: Comparison parameters

    Returns:
        ConditionComparisonResult with cell type stratified results
    """
    await ctx.info(f"Running stratified analysis by {params.cell_type_key}...")

    cell_types = adata.obs[params.cell_type_key].unique()
    await ctx.info(f"Found {len(cell_types)} cell types")

    cell_type_results: list[CellTypeComparisonResult] = []
    total_significant = 0

    for ct in cell_types:
        ct_mask = adata.obs[params.cell_type_key] == ct
        n_cells_ct = ct_mask.sum()

        if n_cells_ct < params.min_cells_per_sample * 2:
            await ctx.warning(
                f"Skipping {ct}: only {n_cells_ct} cells "
                f"(need {params.min_cells_per_sample * 2})"
            )
            continue

        try:
            # Create pseudobulk for this cell type
            counts_df, metadata_df, cell_counts = _create_pseudobulk(
                adata,
                raw_X,
                var_names,
                sample_key=params.sample_key,
                condition_key=params.condition_key,
                cell_type=ct,
                cell_type_key=params.cell_type_key,
                min_cells_per_sample=params.min_cells_per_sample,
            )

            # Check sample distribution
            cond_counts = metadata_df["condition"].value_counts()
            n_cond1 = cond_counts.get(params.condition1, 0)
            n_cond2 = cond_counts.get(params.condition2, 0)

            if n_cond1 < 2 or n_cond2 < 2:
                await ctx.warning(
                    f"Skipping {ct}: insufficient samples "
                    f"({params.condition1}={n_cond1}, {params.condition2}={n_cond2})"
                )
                continue

            # Run DESeq2
            top_up, top_down, n_significant, results_df = _run_deseq2(
                counts_df,
                metadata_df,
                condition1=params.condition1,
                condition2=params.condition2,
                n_top_genes=params.n_top_genes,
                padj_threshold=params.padj_threshold,
                log2fc_threshold=params.log2fc_threshold,
            )

            total_significant += n_significant

            # Count cells per condition for this cell type
            ct_adata = adata[ct_mask]
            n_cells_cond1 = (
                ct_adata.obs[params.condition_key] == params.condition1
            ).sum()
            n_cells_cond2 = (
                ct_adata.obs[params.condition_key] == params.condition2
            ).sum()

            cell_type_results.append(
                CellTypeComparisonResult(
                    cell_type=str(ct),
                    n_cells_condition1=int(n_cells_cond1),
                    n_cells_condition2=int(n_cells_cond2),
                    n_samples_condition1=int(n_cond1),
                    n_samples_condition2=int(n_cond2),
                    n_significant_genes=n_significant,
                    top_upregulated=top_up,
                    top_downregulated=top_down,
                )
            )

            await ctx.info(
                f"{ct}: {n_significant} significant genes "
                f"({len(top_up)} up, {len(top_down)} down)"
            )

        except Exception as e:
            await ctx.warning(f"Analysis failed for {ct}: {e}")
            continue

    if not cell_type_results:
        raise ProcessingError(
            "No cell types had sufficient samples for DESeq2 analysis. "
            "Try lowering min_cells_per_sample or min_samples_per_condition."
        )

    comparison = f"{params.condition1} vs {params.condition2}"

    return ConditionComparisonResult(
        data_id="",  # Will be filled by caller
        method="pseudobulk",
        comparison=comparison,
        condition_key=params.condition_key,
        condition1=params.condition1,
        condition2=params.condition2,
        sample_key=params.sample_key,
        cell_type_key=params.cell_type_key,
        n_samples_condition1=0,  # Will be filled by caller
        n_samples_condition2=0,  # Will be filled by caller
        global_n_significant=None,
        global_top_upregulated=None,
        global_top_downregulated=None,
        cell_type_results=cell_type_results,
        results_key="",  # Will be filled by caller
        statistics={
            "analysis_type": "cell_type_stratified",
            "n_cell_types_analyzed": len(cell_type_results),
            "total_significant_genes": total_significant,
            "cell_types_with_de_genes": len(
                [r for r in cell_type_results if r.n_significant_genes > 0]
            ),
        },
    )
