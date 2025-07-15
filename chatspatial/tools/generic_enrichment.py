"""
Generic enrichment analysis tools for non-spatial data.

This module provides standard enrichment analysis methods including GSEA, ORA, and ssGSEA
that don't require spatial coordinates.
"""

import logging
from typing import Dict, List, Optional, Union, Any, Tuple
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


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
    method: str = 'signal_to_noise',
    permutation_num: int = 1000,
    min_size: int = 10,
    max_size: int = 500,
    context = None
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
        if 'log1p' in adata.uns:
            X = adata.X
        else:
            X = adata.raw.X if adata.raw else adata.X
        
        # Simple fold change if we have groups
        if 'condition' in adata.obs or 'group' in adata.obs:
            group_key = 'condition' if 'condition' in adata.obs else 'group'
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
                if hasattr(X, 'todense'):
                    var_scores = np.array(X.todense().var(axis=0)).flatten()
                else:
                    var_scores = np.array(X.var(axis=0)).flatten()
                ranking = dict(zip(adata.var_names, var_scores))
        else:
            # Use variance as default ranking
            # Handle sparse matrices
            if hasattr(X, 'todense'):
                var_scores = np.array(X.todense().var(axis=0)).flatten()
            else:
                var_scores = np.array(X.var(axis=0)).flatten()
            ranking = dict(zip(adata.var_names, var_scores))
    
    # Run GSEA preranked
    try:
        # Convert ranking dict to DataFrame for gseapy
        ranking_df = pd.DataFrame.from_dict(ranking, orient='index', columns=['score'])
        ranking_df.index.name = 'gene'
        ranking_df = ranking_df.sort_values('score', ascending=False)
        
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
            outdir=None
        )
        
        # Extract results
        results_df = res.res2d
        
        # Prepare output
        enrichment_scores = {}
        pvalues = {}
        adjusted_pvalues = {}
        gene_set_statistics = {}
        
        for idx, row in results_df.iterrows():
            term = row['Term']
            enrichment_scores[term] = row['ES']
            pvalues[term] = row['NOM p-val']
            adjusted_pvalues[term] = row['FDR q-val']
            gene_set_statistics[term] = {
                'es': row['ES'],
                'nes': row['NES'],
                'pval': row['NOM p-val'],
                'fdr': row['FDR q-val'],
                'size': row.get('Matched_size', row.get('Gene %', 0)),  # Different versions use different column names
                'lead_genes': row.get('Lead_genes', '').split(';')[:10] if 'Lead_genes' in row else []
            }
        
        # Get top enriched and depleted
        results_df_sorted = results_df.sort_values('NES', ascending=False)
        top_enriched = results_df_sorted[results_df_sorted['NES'] > 0].head(10)['Term'].tolist()
        top_depleted = results_df_sorted[results_df_sorted['NES'] < 0].head(10)['Term'].tolist()
        
        return {
            'method': 'gsea',
            'n_gene_sets': len(gene_sets),
            'n_significant': len(results_df[results_df['FDR q-val'] < 0.05]),
            'enrichment_scores': enrichment_scores,
            'pvalues': pvalues,
            'adjusted_pvalues': adjusted_pvalues,
            'gene_set_statistics': gene_set_statistics,
            'gene_sets_used': gene_sets,
            'top_gene_sets': top_enriched,
            'top_depleted_sets': top_depleted,
            'results_df': results_df
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
    context = None
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
        if 'rank_genes_groups' in adata.uns:
            # Get DEGs
            result = adata.uns['rank_genes_groups']
            names = result['names']
            pvals = result['pvals_adj'] if 'pvals_adj' in result else result['pvals']
            logfcs = result['logfoldchanges']
            
            # Get first group's DEGs
            degs = []
            for i in range(len(names[0])):
                if pvals[0][i] < pvalue_threshold and abs(logfcs[0][i]) > logfc_threshold:
                    degs.append(names[0][i])
            
            gene_list = degs
            
            if context:
                await context.info(f"Using {len(gene_list)} DEGs for ORA")
        else:
            # Use highly variable genes
            if 'highly_variable' in adata.var:
                gene_list = adata.var_names[adata.var['highly_variable']].tolist()
            else:
                # Use top variable genes
                # Handle sparse matrices
                if hasattr(adata.X, 'todense'):
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
        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        
        enrichment_scores[gs_name] = odds_ratio
        pvalues[gs_name] = p_value
        
        gene_set_statistics[gs_name] = {
            'odds_ratio': odds_ratio,
            'pval': p_value,
            'overlap': a,
            'query_size': len(query_genes),
            'gs_size': len(gs_genes_set),
            'overlapping_genes': list(query_genes & gs_genes_set)[:20]  # Top 20
        }
    
    # Multiple testing correction
    if pvalues:
        pval_array = np.array(list(pvalues.values()))
        _, adjusted_pvals, _, _ = multipletests(pval_array, method='fdr_bh')
        adjusted_pvalues = dict(zip(pvalues.keys(), adjusted_pvals))
    else:
        adjusted_pvalues = {}
    
    # Get top results
    sorted_by_pval = sorted(pvalues.items(), key=lambda x: x[1])
    top_gene_sets = [x[0] for x in sorted_by_pval[:10]]
    
    return {
        'method': 'ora',
        'n_gene_sets': len(gene_sets),
        'n_significant': sum(1 for p in adjusted_pvalues.values() if p < 0.05),
        'enrichment_scores': enrichment_scores,
        'pvalues': pvalues,
        'adjusted_pvalues': adjusted_pvalues,
        'gene_set_statistics': gene_set_statistics,
        'gene_sets_used': gene_sets,
        'query_genes': list(query_genes),
        'top_gene_sets': top_gene_sets,
        'top_depleted_sets': []  # ORA doesn't have depleted sets
    }


async def perform_ssgsea(
    adata,
    gene_sets: Dict[str, List[str]],
    min_size: int = 10,
    max_size: int = 500,
    context = None
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
    if hasattr(adata.X, 'todense'):
        expr_df = pd.DataFrame(
            adata.X.todense().T,
            index=adata.var_names,
            columns=adata.obs_names
        )
    else:
        expr_df = pd.DataFrame(
            adata.X.T,
            index=adata.var_names,
            columns=adata.obs_names
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
            seed=42
        )
        
        # Extract results - ssGSEA stores enrichment scores in res.results
        if hasattr(res, 'results') and isinstance(res.results, dict):
            # res.results is a dict where keys are sample names and values are DataFrames
            # We need to reorganize this into gene sets x samples format
            all_samples = list(res.results.keys())
            all_gene_sets = set()
            
            # Get all gene sets
            for sample_df in res.results.values():
                if isinstance(sample_df, pd.DataFrame) and 'Term' in sample_df.columns:
                    all_gene_sets.update(sample_df['Term'].values)
            
            all_gene_sets = list(all_gene_sets)
            
            # Create scores matrix
            scores_matrix = pd.DataFrame(
                index=all_gene_sets,
                columns=all_samples,
                dtype=float
            )
            
            # Fill in scores
            for sample, df in res.results.items():
                if isinstance(df, pd.DataFrame) and 'Term' in df.columns and 'ES' in df.columns:
                    for _, row in df.iterrows():
                        if row['Term'] in scores_matrix.index:
                            scores_matrix.loc[row['Term'], sample] = row['ES']
            
            scores_df = scores_matrix.fillna(0)  # Fill missing values with 0
        else:
            # Fallback: try to extract from res2d
            scores_df = pd.DataFrame()
            if hasattr(res, 'res2d') and isinstance(res.res2d, pd.DataFrame):
                # res2d contains top results per sample, not full matrix
                # We'll return empty scores for now
                logger.warning("ssGSEA results format not as expected, returning summary only")
        
        # Calculate statistics across samples
        enrichment_scores = {}
        gene_set_statistics = {}
        
        if not scores_df.empty:
            for gs_name in scores_df.index:
                scores = scores_df.loc[gs_name].values
                enrichment_scores[gs_name] = float(np.mean(scores))
                
                gene_set_statistics[gs_name] = {
                'mean_score': float(np.mean(scores)),
                'std_score': float(np.std(scores)),
                'min_score': float(np.min(scores)),
                'max_score': float(np.max(scores)),
                    'size': len(gene_sets.get(gs_name, []))
                }
            
            # Add scores to adata
            for gs_name in scores_df.index:
                adata.obs[f'ssgsea_{gs_name}'] = scores_df.loc[gs_name].values
        
        # Get top gene sets by mean enrichment
        sorted_by_mean = sorted(enrichment_scores.items(), key=lambda x: x[1], reverse=True)
        top_gene_sets = [x[0] for x in sorted_by_mean[:10]]
        
        return {
            'method': 'ssgsea',
            'n_gene_sets': len(gene_sets),
            'n_significant': len(gene_sets),  # All gene sets get scores in ssGSEA
            'enrichment_scores': enrichment_scores,
            'pvalues': {},  # ssGSEA doesn't provide p-values
            'adjusted_pvalues': {},
            'gene_set_statistics': gene_set_statistics,
            'gene_sets_used': gene_sets,
            'top_gene_sets': top_gene_sets,
            'top_depleted_sets': [],
            'scores_df': scores_df,
            'scores_added_to_obs': True
        }
        
    except Exception as e:
        logger.error(f"ssGSEA failed: {e}")
        raise


async def perform_enrichr(
    gene_list: List[str],
    gene_sets: Optional[str] = None,
    organism: str = 'human',
    context = None
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
            'GO_Biological_Process_2023',
            'GO_Molecular_Function_2023',
            'GO_Cellular_Component_2023',
            'KEGG_2021_Human' if organism == 'human' else 'KEGG_2019_Mouse',
            'Reactome_2022',
            'MSigDB_Hallmark_2020'
        ]
    elif isinstance(gene_sets, str):
        gene_sets = [gene_sets]
    
    # Run Enrichr
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism.capitalize(),
            outdir=None,
            cutoff=0.05
        )
        
        # Combine results from all libraries
        all_results = pd.concat(enr.results, ignore_index=True)
        
        # Prepare output
        enrichment_scores = {}
        pvalues = {}
        adjusted_pvalues = {}
        gene_set_statistics = {}
        
        for idx, row in all_results.iterrows():
            term = row['Term']
            enrichment_scores[term] = row['Combined Score']
            pvalues[term] = row['P-value']
            adjusted_pvalues[term] = row['Adjusted P-value']
            
            gene_set_statistics[term] = {
                'combined_score': row['Combined Score'],
                'pval': row['P-value'],
                'adjusted_pval': row['Adjusted P-value'],
                'z_score': row.get('Z-score', np.nan),
                'overlap': row['Overlap'],
                'genes': row['Genes'].split(';') if isinstance(row['Genes'], str) else []
            }
        
        # Get top results
        all_results_sorted = all_results.sort_values('Combined Score', ascending=False)
        top_gene_sets = all_results_sorted.head(10)['Term'].tolist()
        
        return {
            'method': 'enrichr',
            'n_gene_sets': len(all_results),
            'n_significant': len(all_results[all_results['Adjusted P-value'] < 0.05]),
            'enrichment_scores': enrichment_scores,
            'pvalues': pvalues,
            'adjusted_pvalues': adjusted_pvalues,
            'gene_set_statistics': gene_set_statistics,
            'query_genes': gene_list,
            'top_gene_sets': top_gene_sets,
            'top_depleted_sets': [],
            'results_df': all_results
        }
        
    except Exception as e:
        logger.error(f"Enrichr failed: {e}")
        raise