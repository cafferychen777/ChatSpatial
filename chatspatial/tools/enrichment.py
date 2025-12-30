"""
Enrichment analysis tools for spatial transcriptomics data.

This module provides both standard and spatially-aware enrichment analysis methods:
- Standard methods: GSEA, ORA, ssGSEA, Enrichr (via gseapy)
- Spatial methods: EnrichMap-based spatial enrichment analysis
"""

import logging
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import stats

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext
from statsmodels.stats.multitest import multipletests

from ..models.analysis import EnrichmentResult
from ..utils.adata_utils import store_analysis_metadata, to_dense
from ..utils.dependency_manager import require
from ..utils.exceptions import ParameterError, ProcessingError

logger = logging.getLogger(__name__)


# ============================================================================
# MCP RESPONSE OPTIMIZATION
# ============================================================================


def _filter_significant_statistics(
    gene_set_statistics: dict,
    enrichment_scores: dict,
    pvalues: dict,
    adjusted_pvalues: dict,
    method: str,
    fdr_threshold: Optional[float] = None,
) -> tuple:
    """
    Filter all enrichment result dictionaries to only include significant pathways.

    This dramatically reduces MCP response size for large gene set databases
    (e.g., KEGG 311 pathways, GO 10,000 terms) while preserving all important
    information for users.

    Args:
        gene_set_statistics: Full statistics for all gene sets
        enrichment_scores: Enrichment scores for all gene sets
        pvalues: P-values for all gene sets
        adjusted_pvalues: FDR-corrected p-values for all gene sets
        method: Enrichment method used ('gsea', 'ora', 'enrichr', 'ssgsea')
        fdr_threshold: FDR threshold for significance (default: None for method-based auto)
                       Method-based defaults (based on statistical best practices):
                       - GSEA: FDR < 0.25 (official recommendation from Subramanian et al. 2005)
                       - ORA/Enrichr: FDR < 0.05 (standard statistical threshold)
                       - ssGSEA: No filtering (no p-values produced)

    Returns:
        Tuple of (filtered_statistics, filtered_scores, filtered_pvals, filtered_adj_pvals)

    Example:
        Before: 311 pathways × 4 dicts × 100 chars = 124KB (KEGG)
        After: ~15 significant pathways × 4 dicts × 100 chars = 6KB (95% reduction)

    References:
        - GSEA: Subramanian et al. (2005) PNAS 102(43):15545-15550
          "We recommend an FDR cutoff of 25% when dealing with a single database"
        - ORA: Standard multiple testing correction threshold (Benjamini & Hochberg 1995)
    """
    if not adjusted_pvalues:
        # No p-values available (e.g., ssGSEA), return all results without filtering
        return gene_set_statistics, enrichment_scores, pvalues, adjusted_pvalues

    # Auto-determine threshold based on ANALYSIS METHOD if not specified
    # This is statistically principled: different methods have different FDR standards
    if fdr_threshold is None:
        method_lower = method.lower()
        if method_lower == "gsea":
            # GSEA official recommendation: FDR < 0.25
            # From Subramanian et al. 2005: "An FDR of 25% indicates that the result
            # is likely to be valid 3 out of 4 times"
            fdr_threshold = 0.25
        elif method_lower in ("ora", "enrichr", "pathway_ora", "pathway_enrichr"):
            # ORA and Enrichr: standard statistical threshold
            # Based on Benjamini-Hochberg FDR control at 5%
            fdr_threshold = 0.05
        else:
            # Default fallback for unknown methods
            fdr_threshold = 0.05
        logger.info(
            f"Auto-selected FDR threshold {fdr_threshold} for method '{method}' "
            f"(method-based default, not database-size-based)"
        )

    # Find significant pathways
    significant = {
        name
        for name, fdr in adjusted_pvalues.items()
        if fdr is not None and fdr < fdr_threshold
    }

    # Filter all dictionaries
    filtered_stats = {
        name: stats
        for name, stats in gene_set_statistics.items()
        if name in significant
    }

    filtered_scores = {
        name: score for name, score in enrichment_scores.items() if name in significant
    }

    filtered_pvals = {
        name: pval for name, pval in pvalues.items() if name in significant
    }

    filtered_adj_pvals = {
        name: adj_pval
        for name, adj_pval in adjusted_pvalues.items()
        if name in significant
    }

    logger.info(
        f"Filtered enrichment results: {len(adjusted_pvalues)} total → "
        f"{len(filtered_adj_pvals)} significant (FDR < {fdr_threshold})"
    )

    return filtered_stats, filtered_scores, filtered_pvals, filtered_adj_pvals


# ============================================================================
# GENE SET UTILITIES
# ============================================================================


def _filter_gene_sets_by_size(
    gene_sets: Dict[str, List[str]], min_size: int, max_size: int
) -> Dict[str, List[str]]:
    """
    Filter gene sets by size constraints.

    Parameters
    ----------
    gene_sets : Dict[str, List[str]]
        Dictionary mapping gene set names to gene lists
    min_size : int
        Minimum number of genes required
    max_size : int
        Maximum number of genes allowed

    Returns
    -------
    Dict[str, List[str]]
        Filtered gene sets within size constraints
    """
    return {
        name: genes
        for name, genes in gene_sets.items()
        if min_size <= len(genes) <= max_size
    }


# ============================================================================
# SPARSE MATRIX UTILITIES
# ============================================================================


def _compute_std_sparse_compatible(X, axis=0, ddof=1):
    """
    Compute standard deviation compatible with both dense and sparse matrices.

    For sparse matrices, uses the formula: std = sqrt(E[X^2] - E[X]^2) with Bessel correction.
    For dense matrices, uses numpy's built-in std method.

    Args:
        X: Input matrix (can be sparse or dense)
        axis: Axis along which to compute std (0 for columns, 1 for rows)
        ddof: Delta Degrees of Freedom for Bessel correction (default: 1)

    Returns:
        1D numpy array of standard deviations
    """
    import scipy.sparse as sp

    if sp.issparse(X):
        # Sparse matrix: use mathematical formula
        n = X.shape[axis]
        mean = np.array(X.mean(axis=axis)).flatten()
        mean_of_squares = np.array(X.power(2).mean(axis=axis)).flatten()

        # Compute variance with Bessel correction: n/(n-ddof)
        variance = mean_of_squares - np.power(mean, 2)
        variance = np.maximum(variance, 0)  # Avoid numerical errors
        if ddof > 0:
            variance = variance * n / (n - ddof)  # Bessel correction

        return np.sqrt(variance)
    else:
        # Dense matrix: use numpy's built-in method
        return np.array(X.std(axis=axis, ddof=ddof)).flatten()


# ============================================================================
# GENE FORMAT CONVERSION UTILITIES
# ============================================================================


def _convert_gene_format_for_matching(
    pathway_genes: List[str], dataset_genes: set, species: str
) -> Tuple[List[str], Dict[str, str]]:
    """
    Rule-based gene format conversion to match dataset format.

    Handles common gene format variations between pathway databases and datasets:
    - Uppercase (GENE) vs Title case (Gene) vs lowercase (gene)
    - Species-specific formatting rules
    - Special prefixes like Gm/GM/gm for mouse genes

    Args:
        pathway_genes: Gene names from pathway database (usually uppercase from gseapy)
        dataset_genes: Available gene names in dataset
        species: Species specified by user ("mouse" or "human")

    Returns:
        (dataset_format_genes, conversion_map)
        dataset_format_genes: Gene names in dataset format that can be found
        conversion_map: Maps dataset_format -> original_pathway_format
    """
    dataset_format_genes = []
    conversion_map = {}

    for gene in pathway_genes:
        # Try direct match first
        if gene in dataset_genes:
            dataset_format_genes.append(gene)
            conversion_map[gene] = gene
            continue

        # Apply multiple format conversion rules
        format_variations = []

        if species == "mouse":
            # Mouse-specific format rules (order matters for efficiency)
            # Rule 1: Title case (most common): Cd5l, Gbp2b
            if len(gene) > 1:
                format_variations.append(gene[0].upper() + gene[1:].lower())
            # Rule 2: All lowercase: cd5l, gbp2b
            format_variations.append(gene.lower())
            # Rule 3: All uppercase: CD5L, GBP2B
            format_variations.append(gene.upper())
            # Rule 4: Capitalize first letter only
            format_variations.append(gene.capitalize())

            # Special rule for Gm-prefixed genes (common in mouse)
            if gene.upper().startswith("GM"):
                format_variations.extend(
                    [
                        "gm" + gene[2:].lower(),  # gm42418
                        "Gm" + gene[2:].lower(),  # Gm42418
                        "GM" + gene[2:].upper(),  # GM42418
                    ]
                )

        elif species == "human":
            # Human-specific format rules
            # Rule 1: All uppercase (most common): HES1, FABP4
            format_variations.append(gene.upper())
            # Rule 2: All lowercase: hes1, fabp4
            format_variations.append(gene.lower())
            # Rule 3: Capitalize first letter
            format_variations.append(gene.capitalize())

        # Remove duplicates while preserving order
        seen = set()
        unique_variations = []
        for variation in format_variations:
            if variation not in seen and variation != gene:  # Skip if same as original
                seen.add(variation)
                unique_variations.append(variation)

        # Try each format variation against dataset
        for variant in unique_variations:
            if variant in dataset_genes:
                dataset_format_genes.append(variant)  # Use dataset's actual format
                conversion_map[variant] = gene
                break  # Stop after first match

    return dataset_format_genes, conversion_map


# ============================================================================
# ENRICHR DATABASE MAPPING
# ============================================================================


def map_gene_set_database_to_enrichr_library(database_name: str, species: str) -> str:
    """Map user-friendly database names to actual Enrichr library names.

    Args:
        database_name: User-friendly database name from MCP interface
        species: Species ('human', 'mouse', or 'zebrafish')

    Returns:
        Actual Enrichr library name

    Raises:
        ValueError: If database_name is not supported
    """
    mapping = {
        "GO_Biological_Process": "GO_Biological_Process_2025",
        "GO_Molecular_Function": "GO_Molecular_Function_2025",
        "GO_Cellular_Component": "GO_Cellular_Component_2025",
        "KEGG_Pathways": (
            "KEGG_2021_Human" if species.lower() == "human" else "KEGG_2019_Mouse"
        ),
        "Reactome_Pathways": "Reactome_Pathways_2024",
        "MSigDB_Hallmark": "MSigDB_Hallmark_2020",
        "Cell_Type_Markers": "CellMarker_Augmented_2021",
    }

    if database_name not in mapping:
        available_options = list(mapping.keys())
        raise ParameterError(
            f"Unknown gene set database: {database_name}. "
            f"Available options: {available_options}"
        )

    return mapping[database_name]


# ============================================================================
# ENRICHMENT ANALYSIS FUNCTIONS
# ============================================================================


def perform_gsea(
    adata,
    gene_sets: Dict[str, List[str]],
    ranking_key: Optional[str] = None,
    method: str = "signal_to_noise",
    permutation_num: int = 1000,
    min_size: int = 10,
    max_size: int = 500,
    species: Optional[str] = None,
    database: Optional[str] = None,
    ctx: "ToolContext" = None,
) -> "EnrichmentResult":
    """
    Perform Gene Set Enrichment Analysis (GSEA).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : Dict[str, List[str]]
        Gene sets to test
    ranking_key : Optional[str]
        Key in adata.var for pre-computed ranking. If None, compute from expression
    method : str
        Method for ranking genes if ranking_key is None
    permutation_num : int
        Number of permutations
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size
    species : Optional[str]
        Species for the analysis (e.g., 'mouse', 'human')
    database : Optional[str]
        Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process')
    ctx : ToolContext
        MCP tool context for logging

    Returns
    -------
    Dict containing enrichment results
    """
    require("gseapy", ctx, feature="GSEA analysis")
    import gseapy as gp

    # Prepare ranking
    if ranking_key and ranking_key in adata.var:
        # Use pre-computed ranking
        ranking = adata.var[ranking_key].to_dict()
    else:
        # Compute ranking from expression data
        # Use raw data when available and not log-normalized (full gene set for GSEA)
        # IMPORTANT: Keep X and var_names consistent to avoid dimension mismatch
        if "log1p" not in adata.uns and adata.raw is not None:
            X = adata.raw.X
            var_names = adata.raw.var_names
        else:
            X = adata.X
            var_names = adata.var_names

        # Compute gene ranking metric
        # IMPORTANT: GSEA requires biologically meaningful ranking, not just variance
        # Reference: Subramanian et al. (2005) PNAS, GSEA-MSIGDB documentation

        if "condition" in adata.obs or "group" in adata.obs:
            group_key = "condition" if "condition" in adata.obs else "group"
            groups = adata.obs[group_key].unique()

            if len(groups) == 2:
                # Binary comparison: Use Signal-to-Noise Ratio (GSEA default)
                # S2N = (μ1 - μ2) / (σ1 + σ2)
                # This captures both differential expression AND expression stability
                group1_mask = adata.obs[group_key] == groups[0]
                group2_mask = adata.obs[group_key] == groups[1]

                # Compute means
                mean1 = np.array(X[group1_mask, :].mean(axis=0)).flatten()
                mean2 = np.array(X[group2_mask, :].mean(axis=0)).flatten()

                # Compute standard deviations (sparse-compatible)
                std1 = _compute_std_sparse_compatible(X[group1_mask, :], axis=0, ddof=1)
                std2 = _compute_std_sparse_compatible(X[group2_mask, :], axis=0, ddof=1)

                # Apply minimum std threshold (GSEA standard: 0.2 * |mean|)
                # This prevents division by zero and reduces noise from low-variance genes
                min_std_factor = 0.2
                std1 = np.maximum(std1, min_std_factor * np.abs(mean1))
                std2 = np.maximum(std2, min_std_factor * np.abs(mean2))

                # Compute Signal-to-Noise Ratio
                s2n = (mean1 - mean2) / (std1 + std2)
                ranking = dict(zip(var_names, s2n, strict=True))

            else:
                # Multi-group: Use Coefficient of Variation (normalized variance)
                # CV = σ / μ - accounts for mean-variance relationship
                # This is more appropriate than raw variance for genes with different expression levels
                mean = np.array(X.mean(axis=0)).flatten()
                std = _compute_std_sparse_compatible(X, axis=0, ddof=1)

                # Compute CV (avoid division by zero)
                cv = np.zeros_like(mean)
                nonzero_mask = np.abs(mean) > 1e-10
                cv[nonzero_mask] = std[nonzero_mask] / np.abs(mean[nonzero_mask])

                ranking = dict(zip(var_names, cv, strict=False))
        else:
            # No group information: Use best available ranking method
            if "highly_variable_rank" in adata.var:
                # Prefer pre-computed HVG ranking (most robust)
                ranking = adata.var["highly_variable_rank"].to_dict()
            elif "dispersions_norm" in adata.var:
                # Use Seurat-style normalized dispersion
                ranking = adata.var["dispersions_norm"].to_dict()
            else:
                # Fallback: Coefficient of Variation (better than raw variance)
                # Use sparse-compatible std calculation
                mean = np.array(X.mean(axis=0)).flatten()
                std = _compute_std_sparse_compatible(X, axis=0, ddof=1)

                cv = np.zeros_like(mean)
                nonzero_mask = np.abs(mean) > 1e-10
                cv[nonzero_mask] = std[nonzero_mask] / np.abs(mean[nonzero_mask])

                ranking = dict(zip(var_names, cv, strict=False))

    # Run GSEA preranked
    try:
        # Convert ranking dict to DataFrame for gseapy
        ranking_df = pd.DataFrame.from_dict(ranking, orient="index", columns=["score"])
        ranking_df.index.name = "gene"
        ranking_df = ranking_df.sort_values("score", ascending=False)

        res = gp.prerank(
            rnk=ranking_df,  # Pass DataFrame instead of dict
            gene_sets=gene_sets,
            processes=1,
            permutation_num=permutation_num,
            min_size=min_size,
            max_size=max_size,
            seed=42,
            verbose=False,
            no_plot=True,
            outdir=None,
        )

        # Extract results
        results_df = res.res2d

        # Prepare output
        enrichment_scores = {}
        pvalues = {}
        adjusted_pvalues = {}
        gene_set_statistics = {}

        for _idx, row in results_df.iterrows():
            term = row["Term"]
            enrichment_scores[term] = row["ES"]
            pvalues[term] = row["NOM p-val"]
            adjusted_pvalues[term] = row["FDR q-val"]
            gene_set_statistics[term] = {
                "es": row["ES"],
                "nes": row["NES"],
                "pval": row["NOM p-val"],
                "fdr": row["FDR q-val"],
                "size": row.get(
                    "Matched_size", row.get("Gene %", 0)
                ),  # Different versions use different column names
                "lead_genes": (
                    row.get("Lead_genes", "").split(";")[:10]
                    if "Lead_genes" in row
                    else []
                ),
            }

        # Get top enriched and depleted
        results_df_sorted = results_df.sort_values("NES", ascending=False)
        top_enriched = (
            results_df_sorted[results_df_sorted["NES"] > 0].head(10)["Term"].tolist()
        )
        top_depleted = (
            results_df_sorted[results_df_sorted["NES"] < 0].head(10)["Term"].tolist()
        )

        # Save results to adata.uns for visualization
        # Store full results DataFrame for visualization
        adata.uns["gsea_results"] = results_df

        # Store gene set membership for validation
        adata.uns["enrichment_gene_sets"] = gene_sets

        # Store metadata for scientific provenance tracking
        store_analysis_metadata(
            adata,
            analysis_name="enrichment_gsea",
            method="gsea",
            parameters={
                "permutation_num": permutation_num,
                "ranking_method": method,
                "min_size": min_size,
                "max_size": max_size,
                "ranking_key": ranking_key,
            },
            results_keys={"uns": ["gsea_results", "enrichment_gene_sets"]},
            statistics={
                "n_gene_sets": len(gene_sets),
                "n_significant": len(results_df[results_df["FDR q-val"] < 0.05]),
            },
            species=species,
            database=database,
        )

        # Filter all result dictionaries to only significant pathways (reduces MCP response size)
        # Uses method-based FDR threshold: GSEA = 0.25 (Subramanian et al. 2005)
        (
            filtered_statistics,
            filtered_scores,
            filtered_pvals,
            filtered_adj_pvals,
        ) = _filter_significant_statistics(
            gene_set_statistics,
            enrichment_scores,
            pvalues,
            adjusted_pvalues,
            method="gsea",  # Method-based FDR: 0.25 for GSEA
        )

        return EnrichmentResult(
            method="gsea",
            n_gene_sets=len(gene_sets),
            n_significant=len(results_df[results_df["FDR q-val"] < 0.05]),
            enrichment_scores=filtered_scores,
            pvalues=filtered_pvals,
            adjusted_pvalues=filtered_adj_pvals,
            gene_set_statistics=filtered_statistics,
            top_gene_sets=top_enriched,
            top_depleted_sets=top_depleted,
        )

    except Exception as e:
        logger.error(f"GSEA failed: {e}")
        raise


async def perform_ora(
    adata,
    gene_sets: Dict[str, List[str]],
    gene_list: Optional[List[str]] = None,
    pvalue_threshold: float = 0.05,
    min_size: int = 10,
    max_size: int = 500,
    species: Optional[str] = None,
    database: Optional[str] = None,
    ctx: "ToolContext" = None,
) -> "EnrichmentResult":
    """
    Perform Over-Representation Analysis (ORA).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : Dict[str, List[str]]
        Gene sets to test
    gene_list : Optional[List[str]]
        List of genes to test. If None, use DEGs from rank_genes_groups
    pvalue_threshold : float
        P-value threshold for selecting DEGs (only used if rank_genes_groups exists)
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size
    species : Optional[str]
        Species for the analysis (e.g., 'mouse', 'human')
    database : Optional[str]
        Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process')
    ctx : ToolContext
        MCP tool context for logging

    Returns
    -------
    Dict containing enrichment results

    Notes
    -----
    LogFC filtering removed: ORA should use genes pre-filtered by find_markers.
    Different statistical methods (Wilcoxon, t-test) produce different logFC scales,
    making a fixed threshold inappropriate. Gene filtering is the responsibility of
    differential expression analysis, not enrichment analysis.
    """
    # Get gene list if not provided
    if gene_list is None:
        # Try to get DEGs from adata
        if "rank_genes_groups" in adata.uns:
            # Get DEGs
            result = adata.uns["rank_genes_groups"]
            names = result["names"]

            # Check if pvals exist (not all rank_genes_groups have pvals)
            pvals = None
            if "pvals_adj" in result:
                pvals = result["pvals_adj"]
            elif "pvals" in result:
                pvals = result["pvals"]

            # Get DEGs from all groups and merge
            # IMPORTANT: names is a numpy recarray with shape (n_genes,)
            # and dtype.names contains group names as fields
            # Access genes by group name: names[group_name][i]
            degs = []

            # Iterate over all groups
            for group_name in names.dtype.names:
                for i in range(len(names)):
                    # If pvals exist, filter by p-value threshold only
                    # LogFC filtering removed - see docstring Notes
                    if pvals is not None:
                        if pvals[group_name][i] < pvalue_threshold:
                            gene = names[group_name][i]
                            if gene not in degs:  # Avoid duplicates
                                degs.append(gene)
                    else:
                        # If no pvals, use top N genes per group
                        if i < 100:
                            gene = names[group_name][i]
                            if gene not in degs:
                                degs.append(gene)

            gene_list = degs
        else:
            # Use highly variable genes
            if "highly_variable" in adata.var:
                gene_list = adata.var_names[adata.var["highly_variable"]].tolist()
            else:
                # Use top variable genes (based on Coefficient of Variation)
                # CV = σ/μ is more appropriate than raw variance
                mean = np.array(adata.X.mean(axis=0)).flatten()
                std = _compute_std_sparse_compatible(adata.X, axis=0, ddof=1)

                # Compute CV (avoid division by zero)
                cv = np.zeros_like(mean)
                nonzero_mask = np.abs(mean) > 1e-10
                cv[nonzero_mask] = std[nonzero_mask] / np.abs(mean[nonzero_mask])

                top_indices = np.argsort(cv)[-500:]
                gene_list = adata.var_names[top_indices].tolist()

    # Background genes
    # IMPORTANT: Use adata.raw if available, as rank_genes_groups may have been
    # run on raw data with different gene name casing (e.g., uppercase)
    # while filtered adata.var_names may be lowercase
    if adata.raw is not None:
        background_genes = set(adata.raw.var_names)
    else:
        background_genes = set(adata.var_names)

    # Case-insensitive matching as fallback for gene name format differences
    # (e.g., MT.CO1 vs MT-CO1, uppercase vs lowercase)
    query_genes = set(gene_list) & background_genes

    # If no direct matches, try case-insensitive matching
    if len(query_genes) == 0 and len(gene_list) > 0:
        # Create case-insensitive lookup
        gene_name_map = {g.upper(): g for g in background_genes}
        query_genes = set()
        for gene in gene_list:
            if gene.upper() in gene_name_map:
                query_genes.add(gene_name_map[gene.upper()])

    # Perform hypergeometric test for each gene set
    enrichment_scores = {}
    pvalues = {}
    gene_set_statistics = {}

    for gs_name, gs_genes in gene_sets.items():
        gs_genes_set = set(gs_genes) & background_genes

        if len(gs_genes_set) < min_size or len(gs_genes_set) > max_size:
            continue

        # Hypergeometric test
        # a: genes in both query and gene set
        # b: genes in query but not in gene set
        # c: genes in gene set but not in query
        # d: genes in neither

        a = len(query_genes & gs_genes_set)
        b = len(query_genes - gs_genes_set)
        c = len(gs_genes_set - query_genes)
        d = len(background_genes - query_genes - gs_genes_set)

        # Fisher's exact test
        odds_ratio, p_value = stats.fisher_exact(
            [[a, b], [c, d]], alternative="greater"
        )

        enrichment_scores[gs_name] = odds_ratio
        pvalues[gs_name] = p_value

        gene_set_statistics[gs_name] = {
            "odds_ratio": odds_ratio,
            "pval": p_value,
            "overlap": a,
            "query_size": len(query_genes),
            "gs_size": len(gs_genes_set),
            "overlapping_genes": list(query_genes & gs_genes_set)[:20],  # Top 20
        }

    # Multiple testing correction
    if pvalues:
        pval_array = np.array(list(pvalues.values()))
        _, adjusted_pvals, _, _ = multipletests(pval_array, method="fdr_bh")
        adjusted_pvalues = dict(zip(pvalues.keys(), adjusted_pvals, strict=False))
    else:
        adjusted_pvalues = {}

    # Get top results
    sorted_by_pval = sorted(pvalues.items(), key=lambda x: x[1])
    top_gene_sets = [x[0] for x in sorted_by_pval[:10]]

    # Save results to adata.uns for visualization
    # Create DataFrame for visualization compatibility
    ora_df = pd.DataFrame(
        {
            "pathway": list(enrichment_scores.keys()),
            "odds_ratio": list(enrichment_scores.values()),
            "pvalue": [pvalues.get(k, 1.0) for k in enrichment_scores.keys()],
            "adjusted_pvalue": [
                adjusted_pvalues.get(k, 1.0) for k in enrichment_scores.keys()
            ],
        }
    )
    ora_df["NES"] = ora_df["odds_ratio"]  # Use odds_ratio as score for visualization
    ora_df = ora_df.sort_values("pvalue")

    adata.uns["ora_results"] = ora_df
    adata.uns["gsea_results"] = (
        ora_df  # Also save as gsea_results for visualization compatibility
    )

    # Store gene set membership for validation
    adata.uns["enrichment_gene_sets"] = gene_sets

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name="enrichment_ora",
        method="ora",
        parameters={
            "pvalue_threshold": pvalue_threshold,
            "min_size": min_size,
            "max_size": max_size,
            "n_query_genes": len(query_genes),
        },
        results_keys={"uns": ["ora_results", "enrichment_gene_sets"]},
        statistics={
            "n_gene_sets": len(gene_sets),
            "n_significant": sum(
                1 for p in adjusted_pvalues.values() if p is not None and p < 0.05
            ),
            "n_query_genes": len(query_genes),
        },
        species=species,
        database=database,
    )

    # Filter all result dictionaries to only significant pathways (reduces MCP response size)
    # Uses method-based FDR threshold: ORA = 0.05 (standard statistical threshold)
    (
        filtered_statistics,
        filtered_scores,
        filtered_pvals,
        filtered_adj_pvals,
    ) = _filter_significant_statistics(
        gene_set_statistics,
        enrichment_scores,
        pvalues,
        adjusted_pvalues,
        method="ora",  # Method-based FDR: 0.05 for ORA
    )

    return EnrichmentResult(
        method="ora",
        n_gene_sets=len(gene_sets),
        n_significant=sum(
            1 for p in adjusted_pvalues.values() if p is not None and p < 0.05
        ),
        enrichment_scores=filtered_scores,
        pvalues=filtered_pvals,
        adjusted_pvalues=filtered_adj_pvals,
        gene_set_statistics=filtered_statistics,
        top_gene_sets=top_gene_sets,
        top_depleted_sets=[],  # ORA does not produce depleted gene sets
    )


def perform_ssgsea(
    adata,
    gene_sets: Dict[str, List[str]],
    min_size: int = 10,
    max_size: int = 500,
    species: Optional[str] = None,
    database: Optional[str] = None,
    ctx: "ToolContext" = None,
) -> "EnrichmentResult":
    """
    Perform single-sample Gene Set Enrichment Analysis (ssGSEA).

    This calculates enrichment scores for each sample independently.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : Dict[str, List[str]]
        Gene sets to test
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size
    species : Optional[str]
        Species for the analysis (e.g., 'mouse', 'human')
    database : Optional[str]
        Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process')
    ctx : ToolContext
        MCP tool context for logging

    Returns
    -------
    Dict containing enrichment results
    """
    require("gseapy", ctx, feature="ssGSEA analysis")
    import gseapy as gp

    # Memory-efficient batch processing for large datasets
    # Threshold: process in batches if > 1000 samples to avoid OOM
    BATCH_SIZE = 500
    n_samples = adata.n_obs

    # Run ssGSEA (with batch processing for large datasets)
    try:
        if n_samples <= BATCH_SIZE:
            # Small dataset: process all at once (original behavior)
            expr_df = pd.DataFrame(
                to_dense(adata.X).T, index=adata.var_names, columns=adata.obs_names
            )
            res = gp.ssgsea(
                data=expr_df,
                gene_sets=gene_sets,
                min_size=min_size,
                max_size=max_size,
                permutation_num=0,
                no_plot=True,
                threads=1,
                seed=42,
            )
        else:
            # Large dataset: batch processing to reduce peak memory
            # Memory reduction: O(n_genes × n_samples) -> O(n_genes × batch_size)
            all_batch_results = []

            for batch_start in range(0, n_samples, BATCH_SIZE):
                batch_end = min(batch_start + BATCH_SIZE, n_samples)
                batch_indices = list(range(batch_start, batch_end))

                # Extract batch - only convert this batch to dense
                batch_X = to_dense(adata.X[batch_indices, :])
                batch_df = pd.DataFrame(
                    batch_X.T,
                    index=adata.var_names,
                    columns=adata.obs_names[batch_indices],
                )

                batch_res = gp.ssgsea(
                    data=batch_df,
                    gene_sets=gene_sets,
                    min_size=min_size,
                    max_size=max_size,
                    permutation_num=0,
                    no_plot=True,
                    threads=1,
                    seed=42,
                )

                if hasattr(batch_res, "results"):
                    all_batch_results.append(batch_res.results)

                # Free batch memory
                del batch_X, batch_df

            # Merge batch results into unified format
            # Create a mock result object with combined results
            class CombinedResult:
                def __init__(self, results_list):
                    self.results = {}
                    for batch_results in results_list:
                        if isinstance(batch_results, dict):
                            self.results.update(batch_results)

            res = CombinedResult(all_batch_results)

        # Extract results - ssGSEA stores enrichment scores in res.results
        if hasattr(res, "results") and isinstance(res.results, dict):
            # res.results is a dict where keys are sample names and values are DataFrames
            # We need to reorganize this into gene sets x samples format
            all_samples = list(res.results.keys())
            all_gene_sets = set()

            # Get all gene sets
            for sample_df in res.results.values():
                if isinstance(sample_df, pd.DataFrame) and "Term" in sample_df.columns:
                    all_gene_sets.update(sample_df["Term"].values)

            all_gene_sets = list(all_gene_sets)

            # Create scores matrix
            scores_matrix = pd.DataFrame(
                index=all_gene_sets, columns=all_samples, dtype=float
            )

            # Fill in scores
            for sample, df in res.results.items():
                if (
                    isinstance(df, pd.DataFrame)
                    and "Term" in df.columns
                    and "ES" in df.columns
                ):
                    for _, row in df.iterrows():
                        if row["Term"] in scores_matrix.index:
                            scores_matrix.loc[row["Term"], sample] = row["ES"]

            scores_df = scores_matrix.fillna(0)  # Fill missing values with 0
        else:
            error_msg = "ssGSEA results format not recognized."
            logger.error(error_msg)
            raise ProcessingError(error_msg)

        # Calculate statistics across samples
        enrichment_scores = {}
        gene_set_statistics = {}

        if not scores_df.empty:
            for gs_name in scores_df.index:
                scores = scores_df.loc[gs_name].values
                enrichment_scores[gs_name] = float(np.mean(scores))

                gene_set_statistics[gs_name] = {
                    "mean_score": float(np.mean(scores)),
                    "std_score": float(np.std(scores)),
                    "min_score": float(np.min(scores)),
                    "max_score": float(np.max(scores)),
                    "size": len(gene_sets.get(gs_name, [])),
                }

            # Add scores to adata
            for gs_name in scores_df.index:
                adata.obs[f"ssgsea_{gs_name}"] = scores_df.loc[gs_name].values

            # Store gene set membership for validation
            adata.uns["enrichment_gene_sets"] = gene_sets

            # Store metadata for scientific provenance tracking
            obs_keys = [f"ssgsea_{gs_name}" for gs_name in scores_df.index]
            store_analysis_metadata(
                adata,
                analysis_name="enrichment_ssgsea",
                method="ssgsea",
                parameters={
                    "min_size": min_size,
                    "max_size": max_size,
                },
                results_keys={"obs": obs_keys, "uns": ["enrichment_gene_sets"]},
                statistics={
                    "n_gene_sets": len(gene_sets),
                    "n_samples": adata.n_obs,
                },
                species=species,
                database=database,
            )

        # Get top gene sets by mean enrichment
        sorted_by_mean = sorted(
            enrichment_scores.items(), key=lambda x: x[1], reverse=True
        )
        top_gene_sets = [x[0] for x in sorted_by_mean[:10]]

        # ssGSEA doesn't provide p-values, so return empty gene_set_statistics
        # to reduce MCP response size (no significance filtering possible)
        pvalues = None
        adjusted_pvalues = None

        return EnrichmentResult(
            method="ssgsea",
            n_gene_sets=len(gene_sets),
            # IMPORTANT: ssGSEA does NOT perform significance testing
            # Setting n_significant=0 is honest: no pathways are "statistically significant"
            # All gene sets receive enrichment scores, but these are sample-level metrics
            # without associated p-values. Use GSEA or ORA for significance testing.
            n_significant=0,  # ssGSEA doesn't test significance - no p-values produced
            enrichment_scores=enrichment_scores,  # Mean scores per gene set
            pvalues=pvalues,
            adjusted_pvalues=adjusted_pvalues,
            gene_set_statistics={},  # Empty to reduce response size (no p-values available)
            top_gene_sets=top_gene_sets,
            top_depleted_sets=[],  # ssGSEA doesn't produce depleted sets
        )

    except Exception as e:
        logger.error(f"ssGSEA failed: {e}")
        raise


def perform_enrichr(
    gene_list: List[str],
    gene_sets: Optional[str] = None,
    organism: str = "human",
    ctx: "ToolContext" = None,
) -> "EnrichmentResult":
    """
    Perform enrichment analysis using Enrichr web service.

    Parameters
    ----------
    gene_list : List[str]
        List of genes to analyze
    gene_sets : Optional[str]
        Enrichr library name. If None, use default libraries
    organism : str
        Organism ('human' or 'mouse')
    ctx : ToolContext
        MCP tool context for logging

    Returns
    -------
    Dict containing enrichment results
    """
    require("gseapy", ctx, feature="Enrichr analysis")
    import gseapy as gp

    # Default gene set libraries
    if gene_sets is None:
        gene_sets = [
            "GO_Biological_Process_2023",
            "GO_Molecular_Function_2023",
            "GO_Cellular_Component_2023",
            "KEGG_2021_Human" if organism == "human" else "KEGG_2019_Mouse",
            "Reactome_2022",
            "MSigDB_Hallmark_2020",
        ]
    elif isinstance(gene_sets, str):
        # Map user-friendly database name to actual Enrichr library name
        enrichr_library = map_gene_set_database_to_enrichr_library(gene_sets, organism)
        gene_sets = [enrichr_library]

    # Run Enrichr
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism.capitalize(),
            outdir=None,
            cutoff=0.05,
        )

        # Get results - enr.results is already a DataFrame
        all_results = enr.results

        # Prepare output
        enrichment_scores = {}
        pvalues = {}
        adjusted_pvalues = {}
        gene_set_statistics = {}

        # Process all results in a single pass (optimized: 3 loops -> 1)
        genes_found_in_results = []
        for _idx, row in all_results.iterrows():
            term = row["Term"]
            enrichment_scores[term] = row["Combined Score"]
            pvalues[term] = row["P-value"]
            adjusted_pvalues[term] = row["Adjusted P-value"]

            genes_str = row["Genes"]
            genes_list = genes_str.split(";") if isinstance(genes_str, str) else []
            genes_found_in_results.extend(genes_list)

            gene_set_statistics[term] = {
                "combined_score": row["Combined Score"],
                "pval": row["P-value"],
                "adjusted_pval": row["Adjusted P-value"],
                "z_score": row.get("Z-score", np.nan),
                "overlap": row["Overlap"],
                "genes": genes_list,
                "odds_ratio": row.get("Odds Ratio", 1.0),
            }

        # Get top results
        all_results_sorted = all_results.sort_values("Combined Score", ascending=False)
        top_gene_sets = all_results_sorted.head(10)["Term"].tolist()

        # Filter all result dictionaries to only significant pathways (reduces MCP response size)
        # Uses method-based FDR threshold: Enrichr = 0.05 (same as ORA, hypergeometric-based)
        (
            filtered_statistics,
            filtered_scores,
            filtered_pvals,
            filtered_adj_pvals,
        ) = _filter_significant_statistics(
            gene_set_statistics,
            enrichment_scores,
            pvalues,
            adjusted_pvalues,
            method="enrichr",  # Method-based FDR: 0.05 for Enrichr
        )

        return EnrichmentResult(
            method="enrichr",
            n_gene_sets=len(all_results),
            n_significant=len(all_results[all_results["Adjusted P-value"] < 0.05]),
            enrichment_scores=filtered_scores,
            pvalues=filtered_pvals,
            adjusted_pvalues=filtered_adj_pvals,
            gene_set_statistics=filtered_statistics,
            top_gene_sets=top_gene_sets,
            top_depleted_sets=[],  # Enrichr doesn't produce depleted sets
        )

    except Exception as e:
        logger.error(f"Enrichr failed: {e}")
        raise


# ============================================================================
# Spatial Enrichment Analysis Functions (EnrichMap-based)
# ============================================================================


async def perform_spatial_enrichment(
    data_id: str,
    ctx: "ToolContext",
    gene_sets: Union[List[str], Dict[str, List[str]]],
    score_keys: Optional[Union[str, List[str]]] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = True,
    batch_key: Optional[str] = None,
    species: str = "unknown",
    database: Optional[str] = None,
) -> "EnrichmentResult":
    """
    Perform spatially-aware gene set enrichment analysis using EnrichMap.

    Parameters
    ----------
    data_id : str
        Identifier for the spatial data in the data store
    ctx : ToolContext
        MCP tool context for data access and logging
    gene_sets : Union[List[str], Dict[str, List[str]]]
        Either a single gene list or a dictionary of gene sets where keys are
        signature names and values are lists of genes
    score_keys : Optional[Union[str, List[str]]]
        Names for the gene signatures if gene_sets is a list. Ignored if gene_sets
        is already a dictionary
    spatial_key : str
        Key in adata.obsm containing spatial coordinates (default: "spatial")
    n_neighbors : int
        Number of nearest spatial neighbors for smoothing (default: 6)
    smoothing : bool
        Whether to perform spatial smoothing (default: True)
    correct_spatial_covariates : bool
        Whether to correct for spatial covariates using GAM (default: True)
    batch_key : Optional[str]
        Column in adata.obs for batch-wise normalization
    species : str
        Species for the analysis (e.g., 'mouse', 'human')
    database : Optional[str]
        Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process')

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - data_id: ID of the data with enrichment scores
        - signatures: List of computed signatures
        - score_columns: List of column names containing scores
        - gene_contributions: Dictionary of gene contributions per signature
        - summary_stats: Summary statistics for each signature
    """
    # Check if EnrichMap is available
    require("enrichmap", ctx, feature="spatial enrichment analysis")

    # Import EnrichMap
    import enrichmap as em

    # Get data using standard ctx pattern
    adata = await ctx.get_adata(data_id)

    # Validate spatial coordinates
    if spatial_key not in adata.obsm:
        raise ProcessingError(
            f"Spatial coordinates '{spatial_key}' not found in adata.obsm"
        )

    # Convert single gene list to dictionary format
    if isinstance(gene_sets, list):
        if score_keys is None:
            score_keys = "enrichmap_signature"
        gene_sets = {score_keys: gene_sets}

    # Validate gene sets with format conversion
    available_genes = set(adata.var_names)
    validated_gene_sets = {}

    for sig_name, genes in gene_sets.items():
        # Try direct matching first
        common_genes = [gene for gene in genes if gene in available_genes]

        # If few matches and we know the species, try format conversion
        if len(common_genes) < len(genes) * 0.5 and species != "unknown":
            dataset_format_genes, _ = _convert_gene_format_for_matching(
                genes, available_genes, species
            )

            if len(dataset_format_genes) > len(common_genes):
                # Format conversion helped, use dataset format genes for EnrichMap
                common_genes = dataset_format_genes

        logger.info(
            f"Checking signature '{sig_name}': requested {genes[:3]}... found {len(common_genes)}/{len(genes)}"
        )
        if len(common_genes) < 2:
            logger.warning(
                f"Signature '{sig_name}' has {len(common_genes)} genes in the dataset. Skipping."
            )
            continue
        validated_gene_sets[sig_name] = common_genes
        logger.info(
            f"Signature '{sig_name}': {len(common_genes)}/{len(genes)} genes found"
        )

    if not validated_gene_sets:
        raise ProcessingError(
            f"No valid gene signatures found (≥2 genes). "
            f"Dataset: {len(available_genes)} genes, requested: {len(gene_sets)} signatures. "
            f"Check species (human/mouse) and gene name format."
        )

    # Run EnrichMap scoring - process each gene set individually
    failed_signatures = []
    successful_signatures = []

    for sig_name, genes in validated_gene_sets.items():
        try:
            em.tl.score(
                adata=adata,
                gene_set=genes,  # Fixed: use gene_set (correct API parameter name)
                score_key=sig_name,  # Fixed: provide explicit score_key
                spatial_key=spatial_key,
                n_neighbors=n_neighbors,
                smoothing=smoothing,
                correct_spatial_covariates=correct_spatial_covariates,
                batch_key=batch_key,
            )
            successful_signatures.append(sig_name)

        except Exception as e:
            logger.error(f"EnrichMap failed for '{sig_name}': {e}")
            failed_signatures.append((sig_name, str(e)))

    # Check if any signatures were processed successfully
    if not successful_signatures:
        error_details = "; ".join(
            [f"{name}: {error}" for name, error in failed_signatures]
        )
        raise ProcessingError(
            f"All EnrichMap scoring failed. This may indicate:\n"
            f"1. EnrichMap package installation issues\n"
            f"2. Incompatible gene names or data format\n"
            f"3. Insufficient spatial information\n"
            f"Details: {error_details}"
        )

    # Update validated_gene_sets to only include successful ones
    validated_gene_sets = {
        sig: validated_gene_sets[sig] for sig in successful_signatures
    }

    if ctx and failed_signatures:
        await ctx.warning(
            f"Failed to process {len(failed_signatures)} gene sets: {[name for name, _ in failed_signatures]}"
        )

    # Collect results
    score_columns = [f"{sig}_score" for sig in validated_gene_sets.keys()]

    # Calculate summary statistics
    summary_stats = {}
    for sig_name in validated_gene_sets.keys():
        score_col = f"{sig_name}_score"
        scores = adata.obs[score_col]

        summary_stats[sig_name] = {
            "mean": float(scores.mean()),
            "std": float(scores.std()),
            "min": float(scores.min()),
            "max": float(scores.max()),
            "median": float(scores.median()),
            "q25": float(scores.quantile(0.25)),
            "q75": float(scores.quantile(0.75)),
            "n_genes": len(validated_gene_sets[sig_name]),
        }

    # Store gene set membership for validation
    adata.uns["enrichment_gene_sets"] = validated_gene_sets

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name="enrichment_spatial",
        method="spatial_enrichmap",
        parameters={
            "spatial_key": spatial_key,
            "n_neighbors": n_neighbors,
            "smoothing": smoothing,
            "correct_spatial_covariates": correct_spatial_covariates,
            "batch_key": batch_key,
        },
        results_keys={
            "obs": score_columns,
            "uns": ["gene_contributions", "enrichment_gene_sets"],
        },
        statistics={
            "n_gene_sets": len(validated_gene_sets),
            "n_successful_signatures": len(successful_signatures),
            "n_failed_signatures": len(failed_signatures),
        },
        species=species,
        database=database,
    )

    # Create enrichment scores (use max score per gene set)
    enrichment_scores = {
        sig_name: float(stats["max"]) for sig_name, stats in summary_stats.items()
    }

    # Sort by enrichment score to get top gene sets
    sorted_sigs = sorted(enrichment_scores.items(), key=lambda x: x[1], reverse=True)
    top_gene_sets = [sig_name for sig_name, _ in sorted_sigs[:10]]

    # Spatial enrichment doesn't provide p-values, so return empty gene_set_statistics
    # to reduce MCP response size (no significance filtering possible)
    pvalues = None
    adjusted_pvalues = None

    return EnrichmentResult(
        method="spatial_enrichmap",
        n_gene_sets=len(validated_gene_sets),
        n_significant=len(successful_signatures),
        enrichment_scores=enrichment_scores,
        pvalues=pvalues,
        adjusted_pvalues=adjusted_pvalues,
        gene_set_statistics={},  # Empty to reduce response size (no p-values available)
        spatial_scores_key=None,  # Scores are in obs columns, not obsm
        top_gene_sets=top_gene_sets,
        top_depleted_sets=[],  # Spatial enrichment doesn't produce depleted sets
    )


# ============================================================================
# Gene Set Loading Functions
# ============================================================================
# Simplified from GeneSetLoader class - no need for class overhead when
# functions are only called once from load_gene_sets()


def _get_organism_name(species: str) -> str:
    """Get organism name for gseapy from species code."""
    return "Homo sapiens" if species.lower() == "human" else "Mus musculus"


def load_msigdb_gene_sets(
    species: str,
    collection: str = "H",
    subcollection: Optional[str] = None,
    min_size: int = 10,
    max_size: int = 500,
) -> Dict[str, List[str]]:
    """
    Load gene sets from MSigDB using gseapy.

    Parameters
    ----------
    species : str
        Species for gene sets ('human' or 'mouse')
    collection : str
        MSigDB collection name:
        - H: hallmark gene sets
        - C1: positional gene sets
        - C2: curated gene sets (e.g., CGP, CP:KEGG, CP:REACTOME)
        - C3: motif gene sets
        - C4: computational gene sets
        - C5: GO gene sets (CC, BP, MF)
        - C6: oncogenic signatures
        - C7: immunologic signatures
        - C8: cell type signatures
    subcollection : Optional[str]
        Subcollection for specific databases (e.g., 'CP:KEGG', 'GO:BP')
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of gene sets
    """
    try:
        import gseapy as gp

        organism = _get_organism_name(species)
        gene_sets_dict = {}

        if collection == "H":
            # Hallmark gene sets
            gene_sets = gp.get_library_name(organism=organism)
            if "MSigDB_Hallmark_2020" in gene_sets:
                gene_sets_dict = gp.get_library(
                    "MSigDB_Hallmark_2020", organism=organism
                )

        elif collection == "C2" and subcollection == "CP:KEGG":
            # KEGG pathways
            if species.lower() == "human":
                gene_sets_dict = gp.get_library("KEGG_2021_Human", organism=organism)
            else:
                gene_sets_dict = gp.get_library("KEGG_2019_Mouse", organism=organism)

        elif collection == "C2" and subcollection == "CP:REACTOME":
            # Reactome pathways
            gene_sets_dict = gp.get_library("Reactome_2022", organism=organism)

        elif collection == "C5":
            # GO gene sets
            if subcollection == "GO:BP" or subcollection is None:
                gene_sets_dict.update(
                    gp.get_library("GO_Biological_Process_2023", organism=organism)
                )
            if subcollection == "GO:MF" or subcollection is None:
                gene_sets_dict.update(
                    gp.get_library("GO_Molecular_Function_2023", organism=organism)
                )
            if subcollection == "GO:CC" or subcollection is None:
                gene_sets_dict.update(
                    gp.get_library("GO_Cellular_Component_2023", organism=organism)
                )

        elif collection == "C8":
            # Cell type signatures
            gene_sets_dict = gp.get_library(
                "CellMarker_Augmented_2021", organism=organism
            )

        # Filter by size
        filtered_sets = _filter_gene_sets_by_size(gene_sets_dict, min_size, max_size)

        logger.info(f"Loaded {len(filtered_sets)} gene sets from MSigDB {collection}")
        return filtered_sets

    except Exception as e:
        raise ProcessingError(f"Failed to load MSigDB gene sets: {e}") from e


def load_go_gene_sets(
    species: str,
    aspect: str = "BP",
    min_size: int = 10,
    max_size: int = 500,
) -> Dict[str, List[str]]:
    """
    Load GO terms using gseapy.

    Parameters
    ----------
    species : str
        Species for gene sets ('human' or 'mouse')
    aspect : str
        GO aspect: 'BP' (biological process), 'MF' (molecular function),
        'CC' (cellular component)
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of GO gene sets
    """
    aspect_map = {
        "BP": "GO_Biological_Process_2023",
        "MF": "GO_Molecular_Function_2023",
        "CC": "GO_Cellular_Component_2023",
    }

    if aspect not in aspect_map:
        raise ParameterError(f"Invalid GO aspect: {aspect}")

    try:
        import gseapy as gp

        organism = _get_organism_name(species)
        gene_sets = gp.get_library(aspect_map[aspect], organism=organism)

        # Filter by size
        filtered_sets = _filter_gene_sets_by_size(gene_sets, min_size, max_size)

        logger.info(f"Loaded {len(filtered_sets)} GO {aspect} gene sets")
        return filtered_sets

    except Exception as e:
        raise ProcessingError(f"Failed to load GO gene sets: {e}") from e


def load_kegg_gene_sets(
    species: str, min_size: int = 10, max_size: int = 500
) -> Dict[str, List[str]]:
    """
    Load KEGG pathways using gseapy.

    Parameters
    ----------
    species : str
        Species for gene sets ('human' or 'mouse')
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of KEGG pathway gene sets
    """
    try:
        import gseapy as gp

        organism = _get_organism_name(species)

        if species.lower() == "human":
            gene_sets = gp.get_library("KEGG_2021_Human", organism=organism)
        else:
            gene_sets = gp.get_library("KEGG_2019_Mouse", organism=organism)

        # Filter by size
        filtered_sets = _filter_gene_sets_by_size(gene_sets, min_size, max_size)

        logger.info(f"Loaded {len(filtered_sets)} KEGG pathways")
        return filtered_sets

    except Exception as e:
        raise ProcessingError(f"Failed to load KEGG pathways: {e}") from e


def load_reactome_gene_sets(
    species: str, min_size: int = 10, max_size: int = 500
) -> Dict[str, List[str]]:
    """
    Load Reactome pathways using gseapy.

    Parameters
    ----------
    species : str
        Species for gene sets ('human' or 'mouse')
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of Reactome pathway gene sets
    """
    try:
        import gseapy as gp

        organism = _get_organism_name(species)
        gene_sets = gp.get_library("Reactome_2022", organism=organism)

        # Filter by size (use shared utility for consistency)
        filtered_sets = _filter_gene_sets_by_size(gene_sets, min_size, max_size)

        logger.info(f"Loaded {len(filtered_sets)} Reactome pathways")
        return filtered_sets

    except Exception as e:
        raise ProcessingError(f"Failed to load Reactome pathways: {e}") from e


def load_cell_marker_gene_sets(
    species: str, min_size: int = 5, max_size: int = 200
) -> Dict[str, List[str]]:
    """
    Load cell type marker gene sets using gseapy.

    Parameters
    ----------
    species : str
        Species for gene sets ('human' or 'mouse')
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of cell type marker gene sets
    """
    try:
        import gseapy as gp

        organism = _get_organism_name(species)
        gene_sets = gp.get_library("CellMarker_Augmented_2021", organism=organism)

        # Filter by size
        filtered_sets = _filter_gene_sets_by_size(gene_sets, min_size, max_size)

        logger.info(f"Loaded {len(filtered_sets)} cell type marker sets")
        return filtered_sets

    except Exception as e:
        raise ProcessingError(f"Failed to load cell markers: {e}") from e


def load_gene_sets(
    database: str,
    species: str = "human",
    min_genes: int = 10,
    max_genes: int = 500,
    ctx: "ToolContext" = None,
) -> Dict[str, List[str]]:
    """
    Load gene sets from specified database.

    Parameters
    ----------
    database : str
        Database name:
        - GO_Biological_Process, GO_Molecular_Function, GO_Cellular_Component
        - KEGG_Pathways
        - Reactome_Pathways
        - MSigDB_Hallmark
        - Cell_Type_Markers
    species : str
        Species ('human' or 'mouse')
    min_genes : int
        Minimum gene set size
    max_genes : int
        Maximum gene set size
    ctx : ToolContext
        MCP tool context for logging

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of gene sets
    """
    # Direct function calls - no class overhead
    database_map = {
        "GO_Biological_Process": lambda: load_go_gene_sets(
            species, "BP", min_genes, max_genes
        ),
        "GO_Molecular_Function": lambda: load_go_gene_sets(
            species, "MF", min_genes, max_genes
        ),
        "GO_Cellular_Component": lambda: load_go_gene_sets(
            species, "CC", min_genes, max_genes
        ),
        "KEGG_Pathways": lambda: load_kegg_gene_sets(species, min_genes, max_genes),
        "Reactome_Pathways": lambda: load_reactome_gene_sets(
            species, min_genes, max_genes
        ),
        "MSigDB_Hallmark": lambda: load_msigdb_gene_sets(
            species, "H", None, min_genes, max_genes
        ),
        "Cell_Type_Markers": lambda: load_cell_marker_gene_sets(
            species, min_genes, max_genes
        ),
    }

    if database not in database_map:
        raise ParameterError(
            f"Unknown database: {database}. Available: {list(database_map.keys())}"
        )

    gene_sets = database_map[database]()
    return gene_sets
