"""
Enrichment analysis tools for spatial transcriptomics data.

This module provides both standard and spatially-aware enrichment analysis methods:
- Standard methods: GSEA, ORA, ssGSEA, Enrichr (via gseapy)
- Spatial methods: EnrichMap-based spatial enrichment analysis
"""

import logging
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from mcp.server.fastmcp import Context
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..utils.error_handling import ProcessingError

logger = logging.getLogger(__name__)


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
        "KEGG_Pathways": "KEGG_2021_Human" if species.lower() == "human" else "KEGG_2019_Mouse",
        "Reactome_Pathways": "Reactome_Pathways_2024",
        "MSigDB_Hallmark": "MSigDB_Hallmark_2020",
        "Cell_Type_Markers": "CellMarker_Augmented_2021"
    }
    
    if database_name not in mapping:
        available_options = list(mapping.keys())
        raise ValueError(
            f"Unknown gene set database: {database_name}. "
            f"Available options: {available_options}"
        )
    
    return mapping[database_name]


# ============================================================================
# Standard Enrichment Analysis Functions (Non-spatial)
# ============================================================================


def is_gseapy_available() -> Tuple[bool, str]:
    """Check if gseapy is available."""
    try:
        import gseapy

        return True, ""
    except ImportError:
        return False, "gseapy not installed. Install with: pip install gseapy"


async def perform_gsea(
    adata,
    gene_sets: Dict[str, List[str]],
    ranking_key: Optional[str] = None,
    method: str = "signal_to_noise",
    permutation_num: int = 1000,
    min_size: int = 10,
    max_size: int = 500,
    context=None,
) -> Dict[str, Any]:
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
    context : Optional
        MCP context

    Returns
    -------
    Dict containing enrichment results
    """
    is_available, error_msg = is_gseapy_available()
    if not is_available:
        raise ImportError(error_msg)

    import gseapy as gp

    if context:
        await context.info("Running GSEA analysis...")

    # Prepare ranking
    if ranking_key and ranking_key in adata.var:
        # Use pre-computed ranking
        ranking = adata.var[ranking_key].to_dict()
    else:
        # Compute ranking from expression data
        if "log1p" in adata.uns:
            X = adata.X
        else:
            X = adata.raw.X if adata.raw else adata.X

        # Simple fold change if we have groups
        if "condition" in adata.obs or "group" in adata.obs:
            group_key = "condition" if "condition" in adata.obs else "group"
            groups = adata.obs[group_key].unique()
            if len(groups) == 2:
                # Binary comparison
                group1_mask = adata.obs[group_key] == groups[0]
                group2_mask = adata.obs[group_key] == groups[1]

                mean1 = np.array(X[group1_mask, :].mean(axis=0)).flatten()
                mean2 = np.array(X[group2_mask, :].mean(axis=0)).flatten()

                # Compute fold change
                fc = np.log2((mean2 + 1) / (mean1 + 1))
                ranking = dict(zip(adata.var_names, fc))
            else:
                # Use variance as ranking for multi-group
                # Handle sparse matrices
                if hasattr(X, "todense"):
                    var_scores = np.array(X.todense().var(axis=0)).flatten()
                else:
                    var_scores = np.array(X.var(axis=0)).flatten()
                ranking = dict(zip(adata.var_names, var_scores))
        else:
            # Use variance as default ranking
            # Handle sparse matrices
            if hasattr(X, "todense"):
                var_scores = np.array(X.todense().var(axis=0)).flatten()
            else:
                var_scores = np.array(X.var(axis=0)).flatten()
            ranking = dict(zip(adata.var_names, var_scores))

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

        for idx, row in results_df.iterrows():
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

        # Also store as dict format for backward compatibility
        adata.uns["gsea_results_dict"] = {
            "results_df": results_df,
            "enrichment_scores": enrichment_scores,
            "pvalues": pvalues,
            "adjusted_pvalues": adjusted_pvalues,
            "top_enriched": top_enriched,
            "top_depleted": top_depleted,
            "method": "gsea",
        }

        # Inform user about visualization options
        if context:
            await context.info(
                "GSEA analysis complete. Use create_visualization tool with plot_type='pathway_enrichment' to visualize results"
            )

        return {
            "method": "gsea",
            "n_gene_sets": len(gene_sets),
            "n_significant": len(results_df[results_df["FDR q-val"] < 0.05]),
            "enrichment_scores": enrichment_scores,
            "pvalues": pvalues,
            "adjusted_pvalues": adjusted_pvalues,
            "gene_set_statistics": gene_set_statistics,
            "gene_sets_used": {
                k: len(v) for k, v in gene_sets.items()
            },  # Only return gene set sizes
            "top_gene_sets": top_enriched,
            "top_depleted_sets": top_depleted,
            # Don't return the full DataFrame - it's too large
        }

    except Exception as e:
        logger.error(f"GSEA failed: {e}")
        raise


async def perform_ora(
    adata,
    gene_sets: Dict[str, List[str]],
    gene_list: Optional[List[str]] = None,
    pvalue_threshold: float = 0.05,
    logfc_threshold: float = 1.0,
    min_size: int = 10,
    max_size: int = 500,
    context=None,
) -> Dict[str, Any]:
    """
    Perform Over-Representation Analysis (ORA).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : Dict[str, List[str]]
        Gene sets to test
    gene_list : Optional[List[str]]
        List of genes to test. If None, use DEGs
    pvalue_threshold : float
        P-value threshold for selecting DEGs
    logfc_threshold : float
        Log fold change threshold for selecting DEGs
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size
    context : Optional
        MCP context

    Returns
    -------
    Dict containing enrichment results
    """
    if context:
        await context.info("Running Over-Representation Analysis...")

    # Get gene list if not provided
    if gene_list is None:
        # Try to get DEGs from adata
        if "rank_genes_groups" in adata.uns:
            # Get DEGs
            result = adata.uns["rank_genes_groups"]
            names = result["names"]
            pvals = result["pvals_adj"] if "pvals_adj" in result else result["pvals"]
            logfcs = result["logfoldchanges"]

            # Get first group's DEGs
            degs = []
            for i in range(len(names[0])):
                if (
                    pvals[0][i] < pvalue_threshold
                    and abs(logfcs[0][i]) > logfc_threshold
                ):
                    degs.append(names[0][i])

            gene_list = degs

            if context:
                await context.info(f"Using {len(gene_list)} DEGs for ORA")
        else:
            # Use highly variable genes
            if "highly_variable" in adata.var:
                gene_list = adata.var_names[adata.var["highly_variable"]].tolist()
            else:
                # Use top variable genes
                # Handle sparse matrices
                if hasattr(adata.X, "todense"):
                    var_scores = np.array(adata.X.todense().var(axis=0)).flatten()
                else:
                    var_scores = np.array(adata.X.var(axis=0)).flatten()
                top_indices = np.argsort(var_scores)[-500:]
                gene_list = adata.var_names[top_indices].tolist()

    # Background genes
    background_genes = set(adata.var_names)
    query_genes = set(gene_list) & background_genes

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
        adjusted_pvalues = dict(zip(pvalues.keys(), adjusted_pvals))
    else:
        adjusted_pvalues = {}

    # Get top results
    sorted_by_pval = sorted(pvalues.items(), key=lambda x: x[1])
    top_gene_sets = [x[0] for x in sorted_by_pval[:10]]

    # Save results to adata.uns for visualization
    # Create DataFrame for visualization compatibility
    import pandas as pd

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

    # Inform user about visualization options
    if context:
        await context.info(
            "ORA analysis complete. Use create_visualization tool with plot_type='pathway_enrichment' to visualize results"
        )

    return {
        "method": "ora",
        "n_gene_sets": len(gene_sets),
        "n_significant": sum(1 for p in adjusted_pvalues.values() if p is not None and p < 0.05),
        "enrichment_scores": enrichment_scores,
        "pvalues": pvalues,
        "adjusted_pvalues": adjusted_pvalues,
        "gene_set_statistics": gene_set_statistics,
        "gene_sets_used": {
            k: len(v) for k, v in gene_sets.items()
        },  # Only return gene set sizes
        "query_genes": list(query_genes),
        "top_gene_sets": top_gene_sets,
        "top_depleted_sets": [],  # ORA doesn't have depleted sets
    }


async def perform_ssgsea(
    adata,
    gene_sets: Dict[str, List[str]],
    min_size: int = 10,
    max_size: int = 500,
    context=None,
) -> Dict[str, Any]:
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
    context : Optional
        MCP context

    Returns
    -------
    Dict containing enrichment results
    """
    is_available, error_msg = is_gseapy_available()
    if not is_available:
        raise ImportError(error_msg)

    import gseapy as gp

    if context:
        await context.info("Running ssGSEA analysis...")

    # Prepare expression data
    if hasattr(adata.X, "todense"):
        expr_df = pd.DataFrame(
            adata.X.todense().T, index=adata.var_names, columns=adata.obs_names
        )
    else:
        expr_df = pd.DataFrame(
            adata.X.T, index=adata.var_names, columns=adata.obs_names
        )

    # Run ssGSEA
    try:
        res = gp.ssgsea(
            data=expr_df,
            gene_sets=gene_sets,
            min_size=min_size,
            max_size=max_size,
            permutation_num=0,  # No permutation for ssGSEA
            no_plot=True,
            processes=1,
            seed=42,
        )

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
            # ssGSEA results format not recognized - fail honestly instead of returning empty results
            error_msg = (
                "ssGSEA results format not recognized. Expected 'res' attribute with enrichment scores matrix. "
                "This may indicate an issue with the gseapy installation, gene set database, or input data format. "
                "Please check your gene sets and ensure gseapy is properly configured."
            )
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

        # Get top gene sets by mean enrichment
        sorted_by_mean = sorted(
            enrichment_scores.items(), key=lambda x: x[1], reverse=True
        )
        top_gene_sets = [x[0] for x in sorted_by_mean[:10]]

        return {
            "method": "ssgsea",
            "n_gene_sets": len(gene_sets),
            "n_significant": len(gene_sets),  # All gene sets get scores in ssGSEA
            "enrichment_scores": enrichment_scores,
            "pvalues": {},  # ssGSEA doesn't provide p-values
            "adjusted_pvalues": {},
            "gene_set_statistics": gene_set_statistics,
            "gene_sets_used": {
                k: len(v) for k, v in gene_sets.items()
            },  # Only return gene set sizes
            "top_gene_sets": top_gene_sets,
            "top_depleted_sets": [],
            # Don't return the full DataFrame - it's too large
            "scores_added_to_obs": True,
        }

    except Exception as e:
        logger.error(f"ssGSEA failed: {e}")
        raise


async def perform_enrichr(
    gene_list: List[str],
    gene_sets: Optional[str] = None,
    organism: str = "human",
    context=None,
) -> Dict[str, Any]:
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
    context : Optional
        MCP context

    Returns
    -------
    Dict containing enrichment results
    """
    is_available, error_msg = is_gseapy_available()
    if not is_available:
        raise ImportError(error_msg)

    import gseapy as gp

    if context:
        await context.info("Running Enrichr analysis...")

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

        for idx, row in all_results.iterrows():
            term = row["Term"]
            enrichment_scores[term] = row["Combined Score"]
            pvalues[term] = row["P-value"]
            adjusted_pvalues[term] = row["Adjusted P-value"]

            gene_set_statistics[term] = {
                "combined_score": row["Combined Score"],
                "pval": row["P-value"],
                "adjusted_pval": row["Adjusted P-value"],
                "z_score": row.get("Z-score", np.nan),
                "overlap": row["Overlap"],
                "genes": (
                    row["Genes"].split(";") if isinstance(row["Genes"], str) else []
                ),
            }

        # Get top results
        all_results_sorted = all_results.sort_values("Combined Score", ascending=False)
        top_gene_sets = all_results_sorted.head(10)["Term"].tolist()

        return {
            "method": "enrichr",
            "n_gene_sets": len(all_results),
            "n_significant": len(all_results[all_results["Adjusted P-value"] < 0.05]),
            "enrichment_scores": enrichment_scores,
            "pvalues": pvalues,
            "adjusted_pvalues": adjusted_pvalues,
            "gene_set_statistics": gene_set_statistics,
            "query_genes": gene_list,
            "top_gene_sets": top_gene_sets,
            "top_depleted_sets": [],
            # Don't return the full DataFrame - it's too large
        }

    except Exception as e:
        logger.error(f"Enrichr failed: {e}")
        raise


# ============================================================================
# Spatial Enrichment Analysis Functions (EnrichMap-based)
# ============================================================================


def is_enrichmap_available() -> Tuple[bool, str]:
    """Check if EnrichMap is available and all dependencies are met."""
    try:
        import enrichmap as em

        # Check for required dependencies
        # Map package names to their import names
        module_mapping = {
            "scanpy": "scanpy",
            "squidpy": "squidpy",
            "scipy": "scipy",
            "scikit-learn": "sklearn",
            "statsmodels": "statsmodels",
            "pygam": "pygam",
            "scikit-gstat": "skgstat",
            "adjustText": "adjustText",
            "splot": "splot",
        }

        missing = []
        for package, module in module_mapping.items():
            try:
                __import__(module)
            except ImportError:
                missing.append(package)

        if missing:
            return False, f"Missing EnrichMap dependencies: {', '.join(missing)}"

        return True, ""
    except ImportError:
        return False, "EnrichMap not installed. Install with: pip install enrichmap"


async def perform_spatial_enrichment(
    data_id: str,
    data_store: Dict[str, Any],
    gene_sets: Union[List[str], Dict[str, List[str]]],
    score_keys: Optional[Union[str, List[str]]] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = True,
    batch_key: Optional[str] = None,
    gene_weights: Optional[Dict[str, Dict[str, float]]] = None,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """
    Perform spatially-aware gene set enrichment analysis using EnrichMap.

    Parameters
    ----------
    data_id : str
        Identifier for the spatial data in the data store
    data_store : Dict[str, Any]
        Dictionary containing the data
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
    gene_weights : Optional[Dict[str, Dict[str, float]]]
        Pre-computed gene weights for each signature
    context : Optional[Context]
        Execution context

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
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")

    # Import EnrichMap
    import enrichmap as em

    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

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

    # Validate gene sets
    available_genes = set(adata.var_names)
    validated_gene_sets = {}

    # Debug: log total available genes
    logger.info(f"Total available genes in dataset: {len(available_genes)}")
    logger.info(f"First 10 genes: {list(available_genes)[:10]}")

    for sig_name, genes in gene_sets.items():
        # Preserve original gene order while filtering for available genes
        common_genes = [gene for gene in genes if gene in available_genes]
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
        # Collect diagnostic information for better error reporting
        sample_dataset_genes = list(available_genes)[:10]
        sample_requested_genes = []
        for sig_name, genes in gene_sets.items():
            sample_requested_genes.extend(genes[:3])
        sample_requested_genes = list(set(sample_requested_genes))[:10]
        
        error_msg = (
            f"No valid gene signatures found with at least 2 genes.\n\n"
            f"Diagnostic information:\n"
            f"• Dataset contains {len(available_genes)} genes\n"
            f"• Sample dataset genes: {sample_dataset_genes}\n"
            f"• Sample requested genes: {sample_requested_genes}\n"
            f"• Gene signatures provided: {list(gene_sets.keys())}\n\n"
            f"Common causes and solutions:\n"
            f"1. Species mismatch (human vs mouse gene naming)\n"
            f"   - Check if your data species matches the gene set database\n"
            f"2. Gene name format differences (e.g., 'CD3D' vs 'Cd3d')\n"
            f"   - Ensure gene names match your dataset's format\n"
            f"3. Outdated or incompatible gene sets\n"
            f"   - Try a different gene_set_database option\n"
            f"4. Custom gene sets with incorrect gene symbols\n"
            f"   - Verify gene symbols are current and correct"
        )
        raise ProcessingError(error_msg)

    # Run EnrichMap scoring - process each gene set individually
    failed_signatures = []
    successful_signatures = []
    
    for sig_name, genes in validated_gene_sets.items():
        try:
            if context:
                await context.info(f"Processing gene set '{sig_name}' with {len(genes)} genes")
            
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
        error_details = "; ".join([f"{name}: {error}" for name, error in failed_signatures])
        raise ProcessingError(
            f"All EnrichMap scoring failed. This may indicate:\n"
            f"1. EnrichMap package installation issues\n"
            f"2. Incompatible gene names or data format\n"
            f"3. Insufficient spatial information\n"
            f"Details: {error_details}"
        )
    
    # Update validated_gene_sets to only include successful ones
    validated_gene_sets = {sig: validated_gene_sets[sig] for sig in successful_signatures}
    
    if context and failed_signatures:
        await context.warning(f"Failed to process {len(failed_signatures)} gene sets: {[name for name, _ in failed_signatures]}")

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

    # Get gene contributions
    gene_contributions = {}
    if "gene_contributions" in adata.uns:
        gene_contributions = {
            sig: {gene: float(contrib.mean()) for gene, contrib in contribs.items()}
            for sig, contribs in adata.uns["gene_contributions"].items()
        }

    # Inform user about visualization options
    if context:
        await context.info(
            "Spatial enrichment analysis complete. Use create_visualization tool with plot_type='enrichment' to visualize results"
        )

    # Convert to the expected format for the server
    # Create mock p-values and enrichment scores based on the spatial scores
    enrichment_scores = {}
    pvalues = {}
    adjusted_pvalues = {}
    gene_set_statistics = {}

    for sig_name, stats in summary_stats.items():
        # Use the max score as enrichment score (normalized)
        enrichment_scores[sig_name] = (
            stats["max"] / (stats["max"] - stats["min"])
            if (stats["max"] - stats["min"]) > 0
            else 0
        )

        # Spatial enrichment analysis does not provide statistical p-values
        # Real p-values require proper null hypothesis testing with background distributions
        # Users requiring statistical significance should use dedicated enrichment tools
        pvalues[sig_name] = (
            None  # Explicitly set to None to indicate no statistical testing
        )
        adjusted_pvalues[sig_name] = (
            None  # No p-values available for multiple testing correction
        )

        gene_set_statistics[sig_name] = {
            "mean_score": stats["mean"],
            "std_score": stats["std"],
            "min_score": stats["min"],
            "max_score": stats["max"],
            "median_score": stats["median"],
            "q25_score": stats["q25"],
            "q75_score": stats["q75"],
            "n_genes": stats["n_genes"],
            "genes": validated_gene_sets.get(sig_name, []),
        }

    # Sort signatures by their max scores
    sorted_sigs = sorted(
        summary_stats.keys(), key=lambda x: summary_stats[x]["max"], reverse=True
    )

    return {
        "method": "spatial_enrichmap",
        "n_gene_sets": len(validated_gene_sets),
        "n_significant": len([p for p in pvalues.values() if p is not None and p < 0.05]),
        "enrichment_scores": enrichment_scores,
        "pvalues": pvalues,
        "adjusted_pvalues": adjusted_pvalues,
        "gene_set_statistics": gene_set_statistics,
        "gene_sets_used": validated_gene_sets,
        "genes_found": validated_gene_sets,  # For spatial enrichment, all validated genes are used
        "top_gene_sets": sorted_sigs[:10] if sorted_sigs else [],
        "top_depleted_sets": [],  # Spatial enrichment doesn't have depleted sets
        "spatial_metrics": summary_stats,
        "spatial_scores_key": ",".join(
            [f"{sig}_score" for sig in validated_gene_sets.keys()]
        ),
        # Additional spatial-specific information
        "signatures": list(validated_gene_sets.keys()),
        "score_columns": score_columns,
        "gene_contributions": gene_contributions,
        "summary_stats": summary_stats,
        "parameters": {
            "n_neighbors": n_neighbors,
            "smoothing": smoothing,
            "correct_spatial_covariates": correct_spatial_covariates,
            "batch_key": batch_key,
        },
    }


async def compute_spatial_metrics(
    data_id: str,
    data_store: Dict[str, Any],
    score_key: str,
    metrics: Optional[List[str]] = None,
    n_neighbors: int = 6,
    n_perms: int = 999,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """
    Compute spatial metrics for enrichment scores.

    Parameters
    ----------
    data_id : str
        Identifier for the spatial data
    data_store : Dict[str, Any]
        Dictionary containing the data
    score_key : str
        Column name containing the enrichment scores
    metrics : Optional[List[str]]
        List of metrics to compute. Options: ['morans_i', 'getis_ord', 'variance']
        If None, computes all metrics
    n_neighbors : int
        Number of spatial neighbors (default: 6)
    n_perms : int
        Number of permutations for significance testing (default: 999)
    context : Optional[Context]
        Execution context

    Returns
    -------
    Dict[str, Any]
        Dictionary containing computed spatial metrics
    """
    # Check if EnrichMap is available
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")

    import enrichmap as em

    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

    # Validate score column
    if score_key not in adata.obs.columns:
        raise ProcessingError(f"Score column '{score_key}' not found in adata.obs")

    # Default metrics
    if metrics is None:
        metrics = ["morans_i", "getis_ord", "variance"]

    try:
        # Compute spatial metrics
        result = em.tl.compute_spatial_metrics(
            adata=adata,
            score_keys=[score_key],
            metrics=metrics,
            n_neighs=n_neighbors,
            n_perms=n_perms,
        )

        # Extract results for the single score key
        metric_results = {}
        for metric in metrics:
            if metric in result:
                metric_results[metric] = {
                    "value": float(result[metric][score_key]),
                    "p_value": float(
                        result.get(f"{metric}_pval", {}).get(score_key, np.nan)
                    ),
                }

        return {
            "data_id": data_id,
            "score_key": score_key,
            "metrics": metric_results,
            "parameters": {"n_neighbors": n_neighbors, "n_perms": n_perms},
        }

    except Exception as e:
        raise ProcessingError(f"Spatial metrics computation failed: {str(e)}")


async def cluster_gene_correlation(
    data_id: str,
    data_store: Dict[str, Any],
    signature_name: str,
    cluster_key: str = "leiden",
    correlation_method: str = "pearson",
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """
    Compute correlation between gene expression and enrichment scores per cluster.

    Parameters
    ----------
    data_id : str
        Identifier for the spatial data
    data_store : Dict[str, Any]
        Dictionary containing the data
    signature_name : str
        Name of the signature to analyze
    cluster_key : str
        Column in adata.obs containing cluster labels (default: "leiden")
    correlation_method : str
        Correlation method: 'pearson' or 'spearman' (default: "pearson")
    context : Optional[Context]
        Execution context

    Returns
    -------
    Dict[str, Any]
        Dictionary containing correlation results per cluster
    """
    # Check if EnrichMap is available
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")

    import enrichmap as em

    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

    # Validate inputs
    score_key = f"{signature_name}_score"
    if score_key not in adata.obs.columns:
        raise ProcessingError(
            f"Score column '{score_key}' not found. Run enrichment analysis first."
        )

    if cluster_key not in adata.obs.columns:
        raise ProcessingError(f"Cluster column '{cluster_key}' not found in adata.obs")

    try:
        # Compute cluster-gene correlations
        result = em.tl.cluster_gene_correlation(
            adata=adata,
            signature_name=signature_name,
            cluster_key=cluster_key,
            correlation_method=correlation_method,
        )

        # Convert results to serializable format
        correlation_results = {}
        for cluster, corr_df in result.items():
            correlation_results[str(cluster)] = {
                "genes": corr_df.index.tolist(),
                "correlations": corr_df["correlation"].tolist(),
                "top_positive": corr_df.nlargest(10, "correlation")[
                    ["correlation"]
                ].to_dict(),
                "top_negative": corr_df.nsmallest(10, "correlation")[
                    ["correlation"]
                ].to_dict(),
            }

        return {
            "data_id": data_id,
            "signature_name": signature_name,
            "cluster_key": cluster_key,
            "correlation_method": correlation_method,
            "cluster_correlations": correlation_results,
        }

    except Exception as e:
        raise ProcessingError(f"Cluster gene correlation failed: {str(e)}")


# ============================================================================
# Gene Set Loading Functions
# ============================================================================


class GeneSetLoader:
    """Load gene sets from various databases and sources."""

    def __init__(self, species: str = "human"):
        """
        Initialize gene set loader.

        Parameters
        ----------
        species : str
            Species for gene sets ('human' or 'mouse')
        """
        self.species = species.lower()
        self.organism = "Homo sapiens" if species == "human" else "Mus musculus"

    def load_msigdb(
        self,
        collection: str = "H",
        subcollection: Optional[str] = None,
        min_size: int = 10,
        max_size: int = 500,
    ) -> Dict[str, List[str]]:
        """
        Load gene sets from MSigDB using gseapy.

        Parameters
        ----------
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

            # Get available gene sets
            gene_sets_dict = {}

            if collection == "H":
                # Hallmark gene sets
                gene_sets = gp.get_library_name(organism=self.organism)
                if "MSigDB_Hallmark_2020" in gene_sets:
                    gene_sets_dict = gp.get_library(
                        "MSigDB_Hallmark_2020", organism=self.organism
                    )

            elif collection == "C2" and subcollection == "CP:KEGG":
                # KEGG pathways
                if self.species == "human":
                    gene_sets_dict = gp.get_library(
                        "KEGG_2021_Human", organism=self.organism
                    )
                else:
                    gene_sets_dict = gp.get_library(
                        "KEGG_2019_Mouse", organism=self.organism
                    )

            elif collection == "C2" and subcollection == "CP:REACTOME":
                # Reactome pathways
                gene_sets_dict = gp.get_library("Reactome_2022", organism=self.organism)

            elif collection == "C5":
                # GO gene sets
                if subcollection == "GO:BP" or subcollection is None:
                    gene_sets_dict.update(
                        gp.get_library(
                            "GO_Biological_Process_2023", organism=self.organism
                        )
                    )
                if subcollection == "GO:MF" or subcollection is None:
                    gene_sets_dict.update(
                        gp.get_library(
                            "GO_Molecular_Function_2023", organism=self.organism
                        )
                    )
                if subcollection == "GO:CC" or subcollection is None:
                    gene_sets_dict.update(
                        gp.get_library(
                            "GO_Cellular_Component_2023", organism=self.organism
                        )
                    )

            elif collection == "C8":
                # Cell type signatures
                gene_sets_dict = gp.get_library(
                    "CellMarker_Augmented_2021", organism=self.organism
                )

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets_dict.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(
                f"Loaded {len(filtered_sets)} gene sets from MSigDB {collection}"
            )
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load MSigDB gene sets: {e}")
            return {}

    def load_go_terms(
        self, aspect: str = "BP", min_size: int = 10, max_size: int = 500
    ) -> Dict[str, List[str]]:
        """
        Load GO terms using gseapy.

        Parameters
        ----------
        aspect : str
            GO aspect: 'BP' (biological process), 'MF' (molecular function), 'CC' (cellular component)
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
            raise ValueError(f"Invalid GO aspect: {aspect}")

        try:
            import gseapy as gp

            gene_sets = gp.get_library(aspect_map[aspect], organism=self.organism)

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(f"Loaded {len(filtered_sets)} GO {aspect} gene sets")
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load GO gene sets: {e}")
            return {}

    def load_kegg_pathways(
        self, min_size: int = 10, max_size: int = 500
    ) -> Dict[str, List[str]]:
        """Load KEGG pathways."""
        try:
            import gseapy as gp

            if self.species == "human":
                gene_sets = gp.get_library("KEGG_2021_Human", organism=self.organism)
            else:
                gene_sets = gp.get_library("KEGG_2019_Mouse", organism=self.organism)

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(f"Loaded {len(filtered_sets)} KEGG pathways")
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load KEGG pathways: {e}")
            return {}

    def load_reactome_pathways(
        self, min_size: int = 10, max_size: int = 500
    ) -> Dict[str, List[str]]:
        """Load Reactome pathways."""
        try:
            import gseapy as gp

            gene_sets = gp.get_library("Reactome_2022", organism=self.organism)

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(f"Loaded {len(filtered_sets)} Reactome pathways")
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load Reactome pathways: {e}")
            return {}

    def load_cell_markers(
        self, min_size: int = 5, max_size: int = 200
    ) -> Dict[str, List[str]]:
        """Load cell type marker gene sets."""
        try:
            import gseapy as gp

            gene_sets = gp.get_library(
                "CellMarker_Augmented_2021", organism=self.organism
            )

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(f"Loaded {len(filtered_sets)} cell type marker sets")
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load cell markers: {e}")
            return {}


async def load_gene_sets(
    database: str,
    species: str = "human",
    min_genes: int = 10,
    max_genes: int = 500,
    context=None,
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
    context : Optional
        MCP context for logging

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of gene sets
    """
    loader = GeneSetLoader(species=species)

    database_map = {
        "GO_Biological_Process": lambda: loader.load_go_terms(
            "BP", min_genes, max_genes
        ),
        "GO_Molecular_Function": lambda: loader.load_go_terms(
            "MF", min_genes, max_genes
        ),
        "GO_Cellular_Component": lambda: loader.load_go_terms(
            "CC", min_genes, max_genes
        ),
        "KEGG_Pathways": lambda: loader.load_kegg_pathways(min_genes, max_genes),
        "Reactome_Pathways": lambda: loader.load_reactome_pathways(
            min_genes, max_genes
        ),
        "MSigDB_Hallmark": lambda: loader.load_msigdb("H", None, min_genes, max_genes),
        "Cell_Type_Markers": lambda: loader.load_cell_markers(min_genes, max_genes),
    }

    if database not in database_map:
        raise ValueError(
            f"Unknown database: {database}. Available: {list(database_map.keys())}"
        )

    if context:
        await context.info(f"Loading gene sets from {database} for {species}")

    gene_sets = database_map[database]()

    if context:
        await context.info(f"Loaded {len(gene_sets)} gene sets from {database}")

    return gene_sets
