"""
Spatial Statistics Tool

This module provides spatial statistical methods for analyzing spatial transcriptomics data.
Focuses on spatial autocorrelation, local indicators of spatial association (LISA),
and spatial point pattern analysis.

Note: Spatial variable gene identification has been moved to spatial_genes.py
"""

import logging
from typing import List, Optional, Dict, Any, Tuple, Union
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.neighbors import NearestNeighbors
import datetime
from dataclasses import asdict
import json

logger = logging.getLogger(__name__)

class SpatialStatistics:
    """
    A class for performing spatial statistical analysis.

    This class provides methods for:
    - Spatial autocorrelation analysis (Moran's I, Geary's C)
    - Local indicators of spatial association (LISA)
    - Neighborhood enrichment analysis (Squidpy)
    - Bivariate Moran's I analysis (global and local)
    - Spatial point pattern analysis (Ripley's K)

    All local statistics results are automatically written back to adata.obs 
    for downstream visualization and analysis.

    Note: Spatial variable gene identification methods have been moved to spatial_genes.py
    """

    def __init__(self):
        """Initialize the SpatialStatistics tool."""
        pass
    
    def validate_adata_for_spatial_analysis(self, adata: ad.AnnData, spatial_key: str = 'spatial') -> None:
        """
        Validate AnnData object for spatial analysis.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
            
        Raises
        ------
        ValueError
            If adata is not suitable for spatial analysis
        """
        if not isinstance(adata, ad.AnnData):
            raise ValueError(f"Expected AnnData object, got {type(adata)}")
        
        if adata.n_obs == 0:
            raise ValueError("AnnData object contains no observations")
        
        if adata.n_vars == 0:
            raise ValueError("AnnData object contains no variables")
        
        if spatial_key not in adata.obsm:
            raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")
        
        spatial_coords = adata.obsm[spatial_key]
        if spatial_coords.shape[0] != adata.n_obs:
            raise ValueError(f"Spatial coordinates shape {spatial_coords.shape} does not match n_obs {adata.n_obs}")
        
        if spatial_coords.shape[1] < 2:
            raise ValueError(f"Spatial coordinates must be at least 2D, got {spatial_coords.shape[1]}D")
        
        # Check for NaN or infinite values in spatial coordinates
        if np.any(np.isnan(spatial_coords)) or np.any(np.isinf(spatial_coords)):
            raise ValueError("Spatial coordinates contain NaN or infinite values")
        
        logger.info(f"Validation passed: {adata.n_obs} observations, {adata.n_vars} variables, spatial shape: {spatial_coords.shape}")
    
    def validate_genes_in_adata(self, adata: ad.AnnData, genes: List[str]) -> List[str]:
        """
        Validate and filter genes that exist in adata.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        genes : List[str]
            List of genes to validate
            
        Returns
        -------
        List[str]
            List of valid genes that exist in adata
        """
        valid_genes = []
        missing_genes = []
        
        for gene in genes:
            if gene in adata.var_names:
                valid_genes.append(gene)
            else:
                missing_genes.append(gene)
        
        if missing_genes:
            logger.warning(f"Missing genes: {missing_genes[:10]}{'...' if len(missing_genes) > 10 else ''} ({len(missing_genes)} total)")
        
        if len(valid_genes) == 0:
            raise ValueError(f"None of the provided genes found in adata. Check gene names.")
        
        logger.info(f"Found {len(valid_genes)} valid genes out of {len(genes)} requested")
        return valid_genes
    
    def neighborhood_enrichment(
        self,
        adata: ad.AnnData,
        cluster_key: str = 'leiden',
        spatial_key: str = 'spatial',
        n_neighbors: int = 6,
        seed: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Compute neighborhood enrichment analysis using Squidpy.
        
        This method analyzes whether specific cell types are more likely
        to be neighbors of each other than expected by chance.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix with cluster annotations
        cluster_key : str
            Key in adata.obs containing cluster/cell type information
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
        n_neighbors : int
            Number of spatial neighbors to consider
        seed : Optional[int]
            Random seed for reproducibility
            
        Returns
        -------
        pd.DataFrame
            Neighborhood enrichment z-scores matrix
        """
        try:
            import squidpy as sq
        except ImportError:
            logger.error("squidpy not installed. Install with: pip install squidpy")
            raise ImportError("Please install squidpy: pip install squidpy")
        
        if cluster_key not in adata.obs.columns:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
        
        # Create a copy to avoid modifying original data
        adata_copy = adata.copy()
        
        # Build spatial graph
        sq.gr.spatial_neighbors(
            adata_copy,
            coord_type='generic',
            spatial_key=spatial_key,
            n_neighs=n_neighbors,
            set_diag=False
        )
        
        # Compute neighborhood enrichment
        sq.gr.nhood_enrichment(
            adata_copy,
            cluster_key=cluster_key,
            seed=seed
        )
        
        # Extract enrichment z-scores
        enrichment_key = f"{cluster_key}_nhood_enrichment"
        if enrichment_key not in adata_copy.uns:
            raise RuntimeError("Neighborhood enrichment analysis failed")
            
        z_scores = adata_copy.uns[enrichment_key]['zscore']
        
        # Store results in original adata for downstream use
        adata.uns[enrichment_key] = adata_copy.uns[enrichment_key]
        
        logger.info(f"Neighborhood enrichment analysis completed for {len(z_scores)} clusters")
        
        return z_scores
    
    def bivariate_moran(
        self,
        adata: ad.AnnData,
        gene_pairs: List[Tuple[str, str]],
        spatial_key: str = 'spatial',
        n_neighbors: int = 30,
        analysis_type: str = 'both'
    ) -> Dict[str, Any]:
        """
        Compute bivariate Moran's I analysis for gene pairs.
        
        This method analyzes spatial correlation between two different genes,
        useful for studying co-localization patterns.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        gene_pairs : List[Tuple[str, str]]
            List of gene pairs to analyze as (gene1, gene2) tuples
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
        n_neighbors : int
            Number of spatial neighbors to consider
        analysis_type : str
            Type of analysis: 'global', 'local', or 'both'
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing global and/or local bivariate Moran's I results
        """
        try:
            from esda.moran import Moran_BV, Moran_Local_BV
            from libpysal.weights import KNN
        except ImportError:
            logger.error("pysal not installed. Install with: pip install pysal")
            raise ImportError("Please install pysal: pip install pysal")
        
        # Get spatial coordinates
        coords = adata.obsm[spatial_key]
        
        # Create spatial weights with row standardization
        w = KNN.from_array(coords, k=n_neighbors)
        w.transform = 'R'  # Row standardization
        
        results = {}
        
        for gene1, gene2 in gene_pairs:
            if gene1 not in adata.var_names or gene2 not in adata.var_names:
                logger.warning(f"Gene pair ({gene1}, {gene2}) not found in data")
                continue
                
            # Get expression values
            expr1 = adata[:, gene1].X
            expr2 = adata[:, gene2].X
            
            # Convert to arrays if sparse
            if hasattr(expr1, 'toarray'):
                expr1 = expr1.toarray().flatten()
            else:
                expr1 = expr1.flatten()
                
            if hasattr(expr2, 'toarray'):
                expr2 = expr2.toarray().flatten()
            else:
                expr2 = expr2.flatten()
                
            # Clean expression data
            expr1 = np.array(expr1, dtype=float).flatten()
            expr2 = np.array(expr2, dtype=float).flatten()
            
            # Handle infinite values and NaNs
            expr1 = np.where(np.isinf(expr1), np.nan, expr1)
            expr1 = np.nan_to_num(expr1, nan=0.0)
            expr2 = np.where(np.isinf(expr2), np.nan, expr2)
            expr2 = np.nan_to_num(expr2, nan=0.0)
            
            # Skip if either gene has no expression variation
            if np.all(expr1 == 0) or np.all(expr2 == 0):
                logger.warning(f"Gene pair ({gene1}, {gene2}) has no expression variation. Skipping.")
                continue
                
            pair_key = f"{gene1}_{gene2}"
            results[pair_key] = {}
            
            try:
                # Global bivariate Moran's I
                if analysis_type in ['global', 'both']:
                    moran_bv = Moran_BV(expr1, expr2, w)
                    
                    # Handle version compatibility
                    variance = getattr(moran_bv, 'VI_norm', getattr(moran_bv, 'VI_rand', getattr(moran_bv, 'VI', None)))
                    z_score = getattr(moran_bv, 'z_norm', getattr(moran_bv, 'ZI_norm', None))
                    p_value = getattr(moran_bv, 'p_norm', getattr(moran_bv, 'p_sim', None))
                    
                    results[pair_key]['global'] = {
                        'I_bivariate': moran_bv.I,
                        'expected_I': moran_bv.EI_sim if hasattr(moran_bv, 'EI_sim') else moran_bv.EI,
                        'variance': variance,
                        'z_score': z_score,
                        'p_value': p_value
                    }
                
                # Local bivariate Moran's I
                if analysis_type in ['local', 'both']:
                    moran_local_bv = Moran_Local_BV(expr1, expr2, w)
                    
                    results[pair_key]['local'] = {
                        'Is': moran_local_bv.Is,
                        'p_values': moran_local_bv.p_sim,
                        'z_scores': moran_local_bv.z_sim if hasattr(moran_local_bv, 'z_sim') else None
                    }
                    
                    # Write local results back to adata.obs
                    adata.obs[f'{pair_key}_bivariate_moran_I'] = moran_local_bv.Is
                    adata.obs[f'{pair_key}_bivariate_moran_pval'] = moran_local_bv.p_sim
                    if hasattr(moran_local_bv, 'z_sim'):
                        adata.obs[f'{pair_key}_bivariate_moran_zscore'] = moran_local_bv.z_sim
                    
            except Exception as e:
                logger.warning(f"Bivariate Moran's I analysis failed for gene pair ({gene1}, {gene2}): {str(e)}")
                continue
        
        logger.info(f"Bivariate Moran's I analysis completed for {len(results)} gene pairs")
        
        return results
    
    def compute_spatial_autocorrelation(
        self,
        adata: ad.AnnData,
        genes: Optional[List[str]] = None,
        spatial_key: str = 'spatial',
        method: str = 'moran',
        n_neighbors: int = 30
    ) -> pd.DataFrame:
        """
        Compute spatial autocorrelation for genes.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        genes : Optional[List[str]]
            Genes to analyze (if None, use all)
        spatial_key : str
            Spatial coordinate key
        method : str
            Method to use ('moran', 'geary')
        n_neighbors : int
            Number of neighbors for spatial weights
            
        Returns
        -------
        pd.DataFrame
            Autocorrelation statistics
        """
        try:
            from esda.moran import Moran
            from esda.geary import Geary
            from libpysal.weights import KNN
        except ImportError:
            logger.error("pysal not installed. Install with: pip install pysal")
            raise ImportError("Please install pysal: pip install pysal")
        
        # Get spatial coordinates
        coords = adata.obsm[spatial_key]
        
        # Create spatial weights with row standardization for robustness
        w = KNN.from_array(coords, k=n_neighbors)
        w.transform = 'R'  # Row standardization ensures each row sums to 1
        
        # Select genes
        if genes is None:
            # Use highly variable genes if available, otherwise first 100 genes
            if 'highly_variable' in adata.var and adata.var['highly_variable'].any():
                hvg_mask = adata.var['highly_variable']
                genes = adata.var_names[hvg_mask][:100]
                logger.info(f"Computing autocorrelation for top {len(genes)} highly variable genes")
            else:
                genes = adata.var_names[:100]
                logger.info(f"Computing autocorrelation for first 100 genes (HVGs not available)")
        
        results = []
        
        for gene in genes:
            if gene not in adata.var_names:
                continue
            
            # Get expression values and ensure clean numeric data
            expr = adata[:, gene].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()
            
            # Clean expression data: ensure float type and handle NaN/Inf
            expr = np.array(expr, dtype=float).flatten()
            
            # Handle infinite values by replacing with NaN, then NaN with 0
            expr = np.where(np.isinf(expr), np.nan, expr)
            expr = np.nan_to_num(expr, nan=0.0)
            
            # Skip if all values are zero (no expression variation)
            if np.all(expr == 0):
                logger.warning(f"Gene {gene} has no expression variation after cleaning. Skipping.")
                continue
            
            if method == 'moran':
                mi = Moran(expr, w)
                
                # esda version compatibility: handle different attribute names
                variance = getattr(mi, 'VI_norm', getattr(mi, 'VI_rand', getattr(mi, 'VI', None)))
                z_score = getattr(mi, 'z_norm', getattr(mi, 'ZI_norm', None))
                p_value = getattr(mi, 'p_norm', getattr(mi, 'p_sim', None))
                
                if variance is None or z_score is None or p_value is None:
                    logger.warning(f"Missing attributes for gene {gene} in esda Moran. Skipping.")
                    continue
                
                results.append({
                    'gene': gene,
                    'moran_I': mi.I,
                    'expected_I': mi.EI,
                    'variance': variance,
                    'z_score': z_score,
                    'p_value': p_value
                })
            elif method == 'geary':
                c = Geary(expr, w)
                
                # esda version compatibility: handle different attribute names
                variance = getattr(c, 'VC_norm', getattr(c, 'VC_rand', getattr(c, 'VC', None)))
                z_score = getattr(c, 'z_norm', None)
                p_value = getattr(c, 'p_norm', getattr(c, 'p_sim', None))
                
                if variance is None or z_score is None or p_value is None:
                    logger.warning(f"Missing attributes for gene {gene} in esda Geary. Skipping.")
                    continue
                    
                results.append({
                    'gene': gene,
                    'geary_C': c.C,
                    'expected_C': c.EC,
                    'variance': variance,
                    'z_score': z_score,
                    'p_value': p_value
                })
        
        results_df = pd.DataFrame(results)
        
        # Check if results are empty
        if results_df.empty:
            logger.warning("No valid genes found for spatial autocorrelation analysis")
            return results_df  # Return empty DataFrame
        
        # Multiple testing correction
        from statsmodels.stats.multitest import multipletests
        _, results_df['q_value'], _, _ = multipletests(
            results_df['p_value'], method='fdr_bh'
        )
        
        return results_df.sort_values('p_value')
    
    def local_spatial_statistics(
        self,
        adata: ad.AnnData,
        genes: List[str],
        spatial_key: str = 'spatial',
        method: str = 'local_moran',
        n_neighbors: int = 30
    ) -> Dict[str, np.ndarray]:
        """
        Compute local indicators of spatial association (LISA).
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        genes : List[str]
            Genes to analyze
        spatial_key : str
            Spatial coordinate key
        method : str
            Method ('local_moran', 'getis_ord')
        n_neighbors : int
            Number of neighbors
            
        Returns
        -------
        Dict[str, np.ndarray]
            Local statistics for each gene
        """
        try:
            from esda.moran import Moran_Local
            from esda.getisord import G_Local
            from libpysal.weights import KNN
        except ImportError:
            logger.error("pysal not installed")
            raise ImportError("Please install pysal: pip install pysal")
        
        # Get spatial coordinates
        coords = adata.obsm[spatial_key]
        
        # Create spatial weights with row standardization for robustness
        w = KNN.from_array(coords, k=n_neighbors)
        w.transform = 'R'  # Row standardization ensures each row sums to 1
        
        results = {}
        
        for gene in genes:
            if gene not in adata.var_names:
                logger.warning(f"Gene {gene} not found")
                continue
            
            # Get expression values and ensure clean numeric data
            expr = adata[:, gene].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()
            
            # Clean expression data: ensure float type and handle NaN/Inf
            expr = np.array(expr, dtype=float).flatten()
            
            # Handle infinite values by replacing with NaN, then NaN with 0
            expr = np.where(np.isinf(expr), np.nan, expr)
            expr = np.nan_to_num(expr, nan=0.0)
            
            # Skip if all values are zero (no expression variation)
            if np.all(expr == 0):
                logger.warning(f"Gene {gene} has no expression variation after cleaning. Skipping.")
                continue
            
            if method == 'local_moran':
                lisa = Moran_Local(expr, w)
                results[gene] = {
                    'Is': lisa.Is,
                    'p_values': lisa.p_sim,
                    'clusters': lisa.q
                }
            elif method == 'getis_ord':
                g = G_Local(expr, w)
                results[gene] = {
                    'Gi': g.Gs,
                    'p_values': g.p_sim,
                    'z_scores': g.Zs
                }
        
        # Write results back to adata.obs with standardized naming
        for gene in results:
            if method == 'local_moran':
                # Local Moran's I results
                adata.obs[f'{gene}_local_moran_I'] = results[gene]['Is']
                adata.obs[f'{gene}_local_moran_pval'] = results[gene]['p_values']
                adata.obs[f'{gene}_local_moran_cluster'] = results[gene]['clusters']
            elif method == 'getis_ord':
                # Getis-Ord G* results
                adata.obs[f'{gene}_getis_ord_Gi'] = results[gene]['Gi']
                adata.obs[f'{gene}_getis_ord_pval'] = results[gene]['p_values']
                adata.obs[f'{gene}_getis_ord_zscore'] = results[gene]['z_scores']
        
        logger.info(f"Local statistics results written to adata.obs for {len(results)} genes")
        
        return results
    
    def spatial_point_patterns(
        self,
        adata: ad.AnnData,
        cell_type_key: str,
        spatial_key: str = 'spatial',
        method: str = 'ripley'
    ) -> Dict[str, Any]:
        """
        Analyze spatial point patterns for cell types.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        cell_type_key : str
            Key in obs containing cell types
        spatial_key : str
            Spatial coordinate key
        method : str
            Method to use ('ripley', 'nearest_neighbor')
            
        Returns
        -------
        Dict[str, Any]
            Point pattern statistics
        """
        try:
            from pointpats import k, l
            use_pointpats = True
        except ImportError:
            logger.warning("pointpats not installed, using basic implementation")
            use_pointpats = False
        
        coords = adata.obsm[spatial_key]
        cell_types = adata.obs[cell_type_key]
        
        results = {}
        
        for cell_type in cell_types.unique():
            # Get coordinates for this cell type
            mask = cell_types == cell_type
            ct_coords = coords[mask]
            
            if len(ct_coords) < 3:
                continue
                
            try:
                if use_pointpats and method == 'ripley':
                    # Use numpy arrays directly to avoid PointPattern bugs
                    try:
                        # pointpats k and l functions work with raw numpy arrays
                        k_result = k(ct_coords)
                        l_result = l(ct_coords)
                        
                        # k and l return tuples: (support, statistic)
                        k_support, k_statistic = k_result
                        l_support, l_statistic = l_result
                        
                        # Both k and l should have the same support (distance values)
                        # Use k_support as the reference
                        # Calculate enhanced L-function: L(r) = sqrt(K(r)/π) - r
                        l_enhanced = np.sqrt(k_statistic / np.pi) - k_support
                        
                        # Calculate Clark-Evans nearest neighbor index
                        clark_evans_idx = self._calculate_clark_evans_index(ct_coords, coords)
                        
                        # Calculate Hopkins-Skellam index for clustering assessment
                        hopkins_skellam_idx = self._calculate_hopkins_skellam_index(ct_coords, coords)
                        
                        results[cell_type] = {
                            'support': k_support,
                            'K': k_statistic,
                            'L': l_statistic,
                            'L_enhanced': l_enhanced,  # L(r) = sqrt(K(r)/π) - r
                            'clark_evans_index': clark_evans_idx,
                            'hopkins_skellam_index': hopkins_skellam_idx,
                            'l_support': l_support  # Store both for validation if needed
                        }
                    except Exception as e:
                        logger.warning(f"pointpats analysis failed for {cell_type}: {str(e)}")
                        # Fall back to basic implementation for this cell type
                        basic_result = self._basic_point_patterns_single(ct_coords, coords)
                        results[cell_type] = basic_result
                        
                else:
                    # Use basic implementation
                    basic_result = self._basic_point_patterns_single(ct_coords, coords)
                    results[cell_type] = basic_result
                    
            except Exception as e:
                logger.warning(f"Point pattern analysis failed for {cell_type}: {str(e)}")
                continue
        
        return results
    
    def cross_type_point_patterns(
        self,
        adata: ad.AnnData,
        cell_type_key: str,
        type_pairs: List[Tuple[str, str]],
        spatial_key: str = 'spatial'
    ) -> Dict[str, Any]:
        """
        Analyze spatial relationships between different cell types using Ripley's Cross-K function.
        
        This method computes cross-type point pattern statistics to assess whether
        two different cell types show spatial attraction, repulsion, or independence.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix containing spatial coordinates and cell type annotations
        cell_type_key : str
            Key in adata.obs containing cell type information
        type_pairs : List[Tuple[str, str]]
            List of cell type pairs to analyze as (type1, type2) tuples
        spatial_key : str, default 'spatial'
            Key in adata.obsm containing spatial coordinates
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing cross-K statistics for each cell type pair.
            Each entry contains 'support' (distance values) and 'K_cross' (cross-K values).
        """
        try:
            from pointpats import k_cross
            use_pointpats = True
        except ImportError:
            logger.warning("pointpats not installed, falling back to basic cross-correlation analysis")
            use_pointpats = False
        
        if cell_type_key not in adata.obs.columns:
            raise ValueError(f"Cell type key '{cell_type_key}' not found in adata.obs")
        
        coords = adata.obsm[spatial_key]
        cell_types = adata.obs[cell_type_key]
        
        results = {}
        
        for type1, type2 in type_pairs:
            # Check if both cell types exist
            if type1 not in cell_types.values or type2 not in cell_types.values:
                logger.warning(f"Cell type pair ({type1}, {type2}) not found in data")
                continue
            
            # Get coordinates for each cell type
            mask1 = cell_types == type1
            mask2 = cell_types == type2
            coords1 = coords[mask1]
            coords2 = coords[mask2]
            
            # Need at least 3 points of each type for meaningful analysis
            if len(coords1) < 3 or len(coords2) < 3:
                logger.warning(f"Insufficient points for cell type pair ({type1}, {type2})")
                continue
            
            pair_key = f"{type1}_{type2}"
            
            try:
                if use_pointpats:
                    # Use pointpats k_cross function
                    try:
                        cross_k_result = k_cross(coords1, coords2)
                        # k_cross returns (support, statistic)
                        support, k_cross_values = cross_k_result
                        
                        results[pair_key] = {
                            'support': support,
                            'K_cross': k_cross_values,
                            'n_type1': len(coords1),
                            'n_type2': len(coords2)
                        }
                        
                    except Exception as e:
                        logger.warning(f"pointpats k_cross failed for pair ({type1}, {type2}): {str(e)}")
                        # Fall back to basic cross-correlation
                        basic_result = self._basic_cross_correlation(coords1, coords2)
                        results[pair_key] = basic_result
                else:
                    # Use basic cross-correlation implementation
                    basic_result = self._basic_cross_correlation(coords1, coords2)
                    results[pair_key] = basic_result
                    
            except Exception as e:
                logger.warning(f"Cross-type analysis failed for pair ({type1}, {type2}): {str(e)}")
                continue
        
        logger.info(f"Cross-type point pattern analysis completed for {len(results)} cell type pairs")
        return results
    
    def _basic_cross_correlation(
        self, 
        coords1: np.ndarray, 
        coords2: np.ndarray
    ) -> Dict[str, Any]:
        """
        Basic cross-correlation analysis when pointpats is not available.
        
        Parameters
        ----------
        coords1 : np.ndarray
            Coordinates for first cell type
        coords2 : np.ndarray
            Coordinates for second cell type
            
        Returns
        -------
        Dict[str, Any]
            Basic cross-correlation metrics
        """
        # Compute pairwise distances between the two point sets
        from sklearn.metrics.pairwise import euclidean_distances
        
        cross_distances = euclidean_distances(coords1, coords2)
        
        # Get minimum distance from each point in type1 to nearest point in type2
        min_cross_distances = cross_distances.min(axis=1)
        
        return {
            'n_type1': len(coords1),
            'n_type2': len(coords2),
            'mean_cross_distance': float(min_cross_distances.mean()),
            'std_cross_distance': float(min_cross_distances.std()),
            'median_cross_distance': float(np.median(min_cross_distances))
        }
    
    def colocalization_analysis(
        self,
        adata: ad.AnnData,
        cluster_key: str,
        spatial_key: str = 'spatial',
        n_rings: int = 5
    ) -> pd.DataFrame:
        """
        Compute co-occurrence/colocalization metrics using squidpy.
        
        This method analyzes spatial co-occurrence patterns between different 
        cell types using distance-based ring analysis. Higher co-occurrence
        values indicate spatial attraction between cell types.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix containing spatial coordinates and cluster annotations
        cluster_key : str
            Key in adata.obs containing cluster/cell type information
        spatial_key : str, default 'spatial'
            Key in adata.obsm containing spatial coordinates
        n_rings : int, default 5
            Number of distance rings to use for co-occurrence analysis
            
        Returns
        -------
        pd.DataFrame
            Co-occurrence probability matrix between all cluster pairs.
            Values > 1 indicate attraction, < 1 indicate repulsion.
        """
        try:
            import squidpy as sq
        except ImportError:
            logger.warning("squidpy not installed, falling back to basic co-occurrence analysis")
            return self._basic_cooccurrence_analysis(adata, cluster_key, spatial_key, n_rings)
        
        if cluster_key not in adata.obs.columns:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
        
        # Create a copy to avoid modifying original data
        adata_copy = adata.copy()
        
        try:
            # Build spatial graph first
            sq.gr.spatial_neighbors(
                adata_copy,
                coord_type='generic',
                spatial_key=spatial_key,
                n_rings=n_rings,
                set_diag=False
            )
            
            # Compute co-occurrence analysis
            sq.gr.co_occurrence(
                adata_copy,
                cluster_key=cluster_key
            )
            
            # Extract co-occurrence results
            co_occurrence_key = f"{cluster_key}_co_occurrence"
            if co_occurrence_key not in adata_copy.uns:
                raise RuntimeError("Co-occurrence analysis failed - no results found")
            
            # Get the co-occurrence probability matrix
            cooccurrence_data = adata_copy.uns[co_occurrence_key]
            
            # Extract the probability matrix (occ key contains the main results)
            if 'occ' in cooccurrence_data:
                prob_matrix = cooccurrence_data['occ']
            else:
                # Fallback to the first available data
                prob_matrix = list(cooccurrence_data.values())[0]
            
            # Store results in original adata for downstream use
            adata.uns[co_occurrence_key] = adata_copy.uns[co_occurrence_key]
            
            # Also add co-occurrence probabilities to obs if local statistics are available
            if hasattr(prob_matrix, 'values'):  # If it's a pandas DataFrame
                # Add aggregate co-occurrence scores to adata.obs
                cluster_types = adata.obs[cluster_key].unique()
                for cluster_type in cluster_types:
                    if cluster_type in prob_matrix.index:
                        # Mean co-occurrence with all other cell types
                        mean_cooccurrence = prob_matrix.loc[cluster_type].mean()
                        adata.obs[f'{cluster_type}_mean_cooccurrence'] = (
                            adata.obs[cluster_key] == cluster_type
                        ).astype(float) * mean_cooccurrence
            
            logger.info(f"Co-occurrence analysis completed for {len(prob_matrix)} clusters")
            
            return prob_matrix
            
        except Exception as e:
            logger.warning(f"squidpy co-occurrence analysis failed: {str(e)}")
            # Fallback to basic distance-based co-occurrence
            logger.info("Falling back to basic distance-based co-occurrence analysis")
            return self._basic_cooccurrence_analysis(adata, cluster_key, spatial_key, n_rings)
    
    def _basic_cooccurrence_analysis(
        self,
        adata: ad.AnnData,
        cluster_key: str,
        spatial_key: str,
        n_rings: int
    ) -> pd.DataFrame:
        """
        Basic co-occurrence analysis when squidpy fails or produces errors.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        cluster_key : str
            Key in adata.obs containing cluster information
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
        n_rings : int
            Number of distance rings for analysis
            
        Returns
        -------
        pd.DataFrame
            Basic co-occurrence probability matrix
        """
        from sklearn.neighbors import NearestNeighbors
        
        coords = adata.obsm[spatial_key]
        clusters = adata.obs[cluster_key]
        unique_clusters = clusters.unique()
        
        # Create distance rings based on nearest neighbor distances
        nbrs = NearestNeighbors(n_neighbors=min(20, len(coords)))
        nbrs.fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        
        # Determine ring boundaries
        max_distance = np.percentile(distances[:, -1], 90)  # Use 90th percentile as max
        ring_boundaries = np.linspace(0, max_distance, n_rings + 1)
        
        # Initialize co-occurrence matrix
        cooccurrence_matrix = pd.DataFrame(
            0.0, 
            index=unique_clusters, 
            columns=unique_clusters
        )
        
        # Count co-occurrences within each distance ring
        for i, coord in enumerate(coords):
            center_cluster = clusters.iloc[i]
            
            # Find neighbors within each ring
            for ring_idx in range(n_rings):
                min_dist = ring_boundaries[ring_idx]
                max_dist = ring_boundaries[ring_idx + 1]
                
                # Find neighbors in this ring
                ring_distances = distances[i]
                ring_indices = indices[i]
                
                in_ring = (ring_distances >= min_dist) & (ring_distances < max_dist)
                neighbor_clusters = clusters.iloc[ring_indices[in_ring]]
                
                # Count co-occurrences
                for neighbor_cluster in neighbor_clusters:
                    cooccurrence_matrix.loc[center_cluster, neighbor_cluster] += 1
        
        # Normalize to probabilities
        row_sums = cooccurrence_matrix.sum(axis=1)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        cooccurrence_probs = cooccurrence_matrix.div(row_sums, axis=0)
        
        # Write basic results back to adata.obs
        for cluster_type in unique_clusters:
            mean_cooccurrence = cooccurrence_probs.loc[cluster_type].mean()
            adata.obs[f'{cluster_type}_basic_cooccurrence'] = (
                clusters == cluster_type
            ).astype(float) * mean_cooccurrence
        
        logger.info(f"Basic co-occurrence analysis completed for {len(unique_clusters)} clusters")
        
        return cooccurrence_probs
    
    def _calculate_clark_evans_index(
        self, 
        ct_coords: np.ndarray, 
        all_coords: np.ndarray
    ) -> float:
        """
        Calculate Clark-Evans nearest neighbor index.
        
        The Clark-Evans index measures the departure from a random distribution:
        - R = 1: random distribution
        - R < 1: clustered distribution  
        - R > 1: regular/dispersed distribution
        
        Parameters
        ----------
        ct_coords : np.ndarray
            Coordinates for specific cell type
        all_coords : np.ndarray
            All coordinates (for computing study area)
            
        Returns
        -------
        float
            Clark-Evans nearest neighbor index
        """
        if len(ct_coords) < 2:
            return np.nan
        
        # Calculate observed mean nearest neighbor distance
        nbrs = NearestNeighbors(n_neighbors=2)
        nbrs.fit(ct_coords)
        distances, _ = nbrs.kneighbors(ct_coords)
        
        # Exclude self-distance (distance to itself = 0)
        nn_distances = distances[:, 1]
        observed_mean_distance = np.mean(nn_distances)
        
        # Calculate expected mean distance under random distribution
        # Expected distance = 1 / (2 * sqrt(density))
        area = self._compute_area(all_coords)
        density = len(ct_coords) / area
        expected_mean_distance = 1 / (2 * np.sqrt(density))
        
        # Clark-Evans index R
        clark_evans_index = observed_mean_distance / expected_mean_distance
        
        return float(clark_evans_index)
    
    def _calculate_hopkins_skellam_index(
        self, 
        ct_coords: np.ndarray, 
        all_coords: np.ndarray,
        n_sample: int = None
    ) -> float:
        """
        Calculate Hopkins-Skellam index for clustering assessment.
        
        The Hopkins-Skellam index measures spatial clustering:
        - H ≈ 0.5: random distribution
        - H < 0.5: regular/dispersed distribution
        - H > 0.5: clustered distribution
        
        Parameters
        ----------
        ct_coords : np.ndarray
            Coordinates for specific cell type
        all_coords : np.ndarray
            All coordinates (for defining study area bounds)
        n_sample : int, optional
            Number of random points to sample (default: min(30, n_points//2))
            
        Returns
        -------
        float
            Hopkins-Skellam clustering index
        """
        if len(ct_coords) < 3:
            return np.nan
        
        # Determine sample size
        if n_sample is None:
            n_sample = min(30, len(ct_coords) // 2)
        
        if n_sample < 1:
            return np.nan
        
        # Get bounds of the study area
        min_coords = all_coords.min(axis=0)
        max_coords = all_coords.max(axis=0)
        
        # Sample random points within the study area bounds
        np.random.seed(42)  # For reproducibility
        random_points = np.random.uniform(
            low=min_coords,
            high=max_coords,
            size=(n_sample, ct_coords.shape[1])
        )
        
        # Sample n_sample points from the actual point pattern
        sample_indices = np.random.choice(len(ct_coords), n_sample, replace=False)
        sampled_points = ct_coords[sample_indices]
        
        # Calculate distances to nearest neighbors
        nbrs_actual = NearestNeighbors(n_neighbors=1)
        nbrs_actual.fit(ct_coords)
        
        # Distance from random points to nearest actual point
        u_distances, _ = nbrs_actual.kneighbors(random_points)
        u_distances = u_distances.flatten()
        
        # Distance from sampled actual points to nearest other actual point
        nbrs_others = NearestNeighbors(n_neighbors=2)  # 2 to exclude self
        nbrs_others.fit(ct_coords)
        w_distances, _ = nbrs_others.kneighbors(sampled_points)
        w_distances = w_distances[:, 1]  # Take second nearest (first is self)
        
        # Hopkins-Skellam statistic
        u_sum = np.sum(u_distances)
        w_sum = np.sum(w_distances)
        
        hopkins_index = u_sum / (u_sum + w_sum)
        
        return float(hopkins_index)
    
    def _basic_point_patterns_single(
        self,
        ct_coords: np.ndarray,
        all_coords: np.ndarray
    ) -> Dict[str, float]:
        """
        Basic point pattern analysis for a single cell type.
        """
        # Basic nearest neighbor analysis
        nbrs = NearestNeighbors(n_neighbors=min(2, len(ct_coords)))
        nbrs.fit(ct_coords)
        distances, _ = nbrs.kneighbors(ct_coords)
        
        # Exclude self-distance
        nn_distances = distances[:, 1] if distances.shape[1] > 1 else distances[:, 0]
        
        # Calculate additional indices
        clark_evans_idx = self._calculate_clark_evans_index(ct_coords, all_coords)
        hopkins_skellam_idx = self._calculate_hopkins_skellam_index(ct_coords, all_coords)
        
        return {
            'n_cells': len(ct_coords),
            'mean_nn_distance': float(nn_distances.mean()),
            'std_nn_distance': float(nn_distances.std()),
            'density': len(ct_coords) / self._compute_area(all_coords),
            'clark_evans_index': clark_evans_idx,
            'hopkins_skellam_index': hopkins_skellam_idx
        }
    
    def _basic_point_patterns(
        self,
        adata: ad.AnnData,
        cell_type_key: str,
        spatial_key: str
    ) -> Dict[str, Any]:
        """
        Basic point pattern analysis without pointpats.
        """
        coords = adata.obsm[spatial_key]
        cell_types = adata.obs[cell_type_key]
        
        results = {}
        
        for cell_type in cell_types.unique():
            mask = cell_types == cell_type
            ct_coords = coords[mask]
            
            if len(ct_coords) < 3:
                continue
            
            # Basic nearest neighbor analysis
            nbrs = NearestNeighbors(n_neighbors=min(2, len(ct_coords)))
            nbrs.fit(ct_coords)
            distances, _ = nbrs.kneighbors(ct_coords)
            
            # Exclude self-distance
            nn_distances = distances[:, 1] if distances.shape[1] > 1 else distances[:, 0]
            
            # Calculate additional indices
            clark_evans_idx = self._calculate_clark_evans_index(ct_coords, coords)
            hopkins_skellam_idx = self._calculate_hopkins_skellam_index(ct_coords, coords)
            
            results[cell_type] = {
                'n_cells': len(ct_coords),
                'mean_nn_distance': float(nn_distances.mean()),
                'std_nn_distance': float(nn_distances.std()),
                'density': len(ct_coords) / self._compute_area(coords),
                'clark_evans_index': clark_evans_idx,
                'hopkins_skellam_index': hopkins_skellam_idx
            }
        
        return results
    
    def _compute_area(self, coords: np.ndarray) -> float:
        """Compute bounding box area of coordinates."""
        min_coords = coords.min(axis=0)
        max_coords = coords.max(axis=0)
        return np.prod(max_coords - min_coords)
    
    def integrate_analysis_results(
        self,
        adata: ad.AnnData,
        analysis_name: str = "spatial_analysis",
        timestamp: Optional[str] = None
    ) -> str:
        """
        Systematically organize all spatial statistics results in adata.
        
        This method consolidates all spatial analysis results into a standardized format
        with proper naming conventions, metadata, and result organization for easy
        visualization and downstream analysis.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix with spatial statistics results
        analysis_name : str
            Name for this analysis session (used as prefix for result keys)
        timestamp : Optional[str]
            Timestamp string, if None current time will be used
            
        Returns
        -------
        str
            Key prefix used for storing integrated results
        """
        if timestamp is None:
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        
        result_key = f"{analysis_name}_{timestamp}"
        
        logger.info(f"Integrating spatial analysis results with key: {result_key}")
        
        # Initialize integration metadata
        integration_metadata = {
            'analysis_name': analysis_name,
            'timestamp': timestamp,
            'n_obs': adata.n_obs,
            'n_vars': adata.n_vars,
            'spatial_key': 'spatial',  # Standard spatial key
            'methods_used': [],
            'local_statistics': {},
            'global_statistics': {},
            'result_keys': {}
        }
        
        # Scan adata.obs for existing spatial statistics results
        obs_columns = adata.obs.columns.tolist()
        
        # Identify local Moran's I results
        moran_i_cols = [col for col in obs_columns if col.endswith('_local_moran_I')]
        moran_pval_cols = [col for col in obs_columns if col.endswith('_local_moran_pval')]
        moran_cluster_cols = [col for col in obs_columns if col.endswith('_local_moran_cluster')]
        
        # Identify Getis-Ord results  
        getis_gi_cols = [col for col in obs_columns if col.endswith('_getis_ord_Gi')]
        getis_pval_cols = [col for col in obs_columns if col.endswith('_getis_ord_pval')]
        getis_zscore_cols = [col for col in obs_columns if col.endswith('_getis_ord_zscore')]
        
        # Identify bivariate Moran results
        biv_moran_i_cols = [col for col in obs_columns if '_bivariate_moran_I' in col]
        biv_moran_pval_cols = [col for col in obs_columns if '_bivariate_moran_pval' in col]
        
        # Organize local Moran's I results
        if moran_i_cols:
            integration_metadata['methods_used'].append('local_moran')
            moran_genes = []
            for col in moran_i_cols:
                gene = col.replace('_local_moran_I', '')
                moran_genes.append(gene)
                
                # Store in standardized format
                adata.obs[f"{result_key}_local_moran_{gene}_I"] = adata.obs[col]
                if f"{gene}_local_moran_pval" in obs_columns:
                    adata.obs[f"{result_key}_local_moran_{gene}_pval"] = adata.obs[f"{gene}_local_moran_pval"]
                if f"{gene}_local_moran_cluster" in obs_columns:
                    adata.obs[f"{result_key}_local_moran_{gene}_cluster"] = adata.obs[f"{gene}_local_moran_cluster"]
            
            integration_metadata['local_statistics']['moran_genes'] = moran_genes
            integration_metadata['result_keys']['local_moran'] = f"{result_key}_local_moran"
        
        # Organize Getis-Ord results
        if getis_gi_cols:
            integration_metadata['methods_used'].append('getis_ord')
            getis_genes = []
            for col in getis_gi_cols:
                gene = col.replace('_getis_ord_Gi', '')
                getis_genes.append(gene)
                
                # Store in standardized format
                adata.obs[f"{result_key}_getis_ord_{gene}_Gi"] = adata.obs[col]
                if f"{gene}_getis_ord_pval" in obs_columns:
                    adata.obs[f"{result_key}_getis_ord_{gene}_pval"] = adata.obs[f"{gene}_getis_ord_pval"]
                if f"{gene}_getis_ord_zscore" in obs_columns:
                    adata.obs[f"{result_key}_getis_ord_{gene}_zscore"] = adata.obs[f"{gene}_getis_ord_zscore"]
            
            integration_metadata['local_statistics']['getis_genes'] = getis_genes  
            integration_metadata['result_keys']['getis_ord'] = f"{result_key}_getis_ord"
        
        # Organize bivariate Moran results
        if biv_moran_i_cols:
            integration_metadata['methods_used'].append('bivariate_moran')
            bivariate_pairs = []
            for col in biv_moran_i_cols:
                pair_key = col.replace('_bivariate_moran_I', '')
                bivariate_pairs.append(pair_key)
                
                # Store in standardized format
                adata.obs[f"{result_key}_bivariate_moran_{pair_key}_I"] = adata.obs[col]
                if f"{pair_key}_bivariate_moran_pval" in obs_columns:
                    adata.obs[f"{result_key}_bivariate_moran_{pair_key}_pval"] = adata.obs[f"{pair_key}_bivariate_moran_pval"]
            
            integration_metadata['local_statistics']['bivariate_pairs'] = bivariate_pairs
            integration_metadata['result_keys']['bivariate_moran'] = f"{result_key}_bivariate_moran"
        
        # Check for global statistics in adata.uns
        if hasattr(adata, 'uns') and adata.uns:
            # Look for neighborhood enrichment results
            enrichment_keys = [key for key in adata.uns.keys() if 'nhood_enrichment' in key]
            if enrichment_keys:
                integration_metadata['methods_used'].append('neighborhood_enrichment')
                integration_metadata['result_keys']['neighborhood_enrichment'] = enrichment_keys[0]
                integration_metadata['global_statistics']['neighborhood_enrichment'] = {
                    'enrichment_matrix_shape': adata.uns[enrichment_keys[0]]['zscore'].shape
                }
        
        # Store integration metadata
        adata.uns[f"{result_key}_integration_metadata"] = integration_metadata
        
        # Create visualization-ready feature lists
        all_features = []
        if moran_i_cols:
            all_features.extend([col.replace('_local_moran_I', '') for col in moran_i_cols])
        if getis_gi_cols:
            all_features.extend([col.replace('_getis_ord_Gi', '') for col in getis_gi_cols])
        
        adata.uns[f"{result_key}_visualization_features"] = {
            'all_spatial_features': list(set(all_features)),
            'moran_features': [col.replace('_local_moran_I', '') for col in moran_i_cols],
            'getis_features': [col.replace('_getis_ord_Gi', '') for col in getis_gi_cols],
            'bivariate_features': [col.replace('_bivariate_moran_I', '') for col in biv_moran_i_cols]
        }
        
        logger.info(f"Integration completed. Found {len(integration_metadata['methods_used'])} spatial analysis methods")
        logger.info(f"Total spatial features: {len(all_features)}")
        
        return result_key
    
    def batch_spatial_analysis(
        self,
        adata: ad.AnnData,
        analysis_types: List[str],
        feature_lists: Optional[Dict[str, List[str]]] = None,
        spatial_key: str = 'spatial',
        n_neighbors: int = 30,
        optimize_performance: bool = True,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Run multiple spatial statistics in one optimized call.
        
        This method efficiently computes multiple spatial statistics by reusing
        spatial graph construction and optimizing data access patterns.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        analysis_types : List[str]
            List of analysis types to run. Options: 
            - 'autocorrelation_moran': Global Moran's I
            - 'autocorrelation_geary': Global Geary's C  
            - 'local_moran': Local Moran's I (LISA)
            - 'getis_ord': Getis-Ord G*
            - 'bivariate_moran': Bivariate Moran's I
            - 'neighborhood_enrichment': Cell type neighborhood enrichment
            - 'point_patterns': Spatial point pattern analysis
        feature_lists : Optional[Dict[str, List[str]]]
            Dictionary mapping analysis types to specific gene lists
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
        n_neighbors : int
            Number of spatial neighbors to consider
        optimize_performance : bool
            Whether to optimize for performance (reuse spatial weights)
        **kwargs
            Additional parameters passed to individual analysis methods
            
        Returns
        -------
        Dict[str, Any]
            Comprehensive results dictionary with all analysis results
        """
        logger.info(f"Starting batch spatial analysis with {len(analysis_types)} analysis types")
        
        # Validate adata for spatial analysis
        try:
            self.validate_adata_for_spatial_analysis(adata, spatial_key)
        except ValueError as e:
            logger.error(f"AnnData validation failed: {e}")
            raise
        
        # Get spatial coordinates
        coords = adata.obsm[spatial_key]
        
        # Initialize results
        batch_results = {
            'analysis_types': analysis_types,
            'n_neighbors': n_neighbors,
            'n_obs': adata.n_obs,
            'n_vars': adata.n_vars,
            'spatial_coordinates_shape': coords.shape,
            'results': {},
            'performance_metrics': {
                'start_time': datetime.datetime.now().isoformat(),
                'analysis_times': {}
            }
        }
        
        # Create spatial weights once if optimizing performance
        spatial_weights = None
        if optimize_performance and any(analysis in analysis_types for analysis in 
                                      ['autocorrelation_moran', 'autocorrelation_geary', 
                                       'local_moran', 'getis_ord', 'bivariate_moran']):
            try:
                from libpysal.weights import KNN
                spatial_weights = KNN.from_array(coords, k=n_neighbors)
                spatial_weights.transform = 'R'  # Row standardization
                logger.info("Created shared spatial weights matrix for optimized performance")
            except ImportError:
                logger.warning("libpysal not available, spatial weights optimization disabled")
                optimize_performance = False
        
        # Process each analysis type
        for analysis_type in analysis_types:
            start_time = datetime.datetime.now()
            logger.info(f"Running {analysis_type} analysis")
            
            try:
                # Get feature list for this analysis
                if feature_lists and analysis_type in feature_lists:
                    features = feature_lists[analysis_type]
                    # Validate genes exist in data
                    features = self.validate_genes_in_adata(adata, features)
                else:
                    # Use default feature selection
                    if hasattr(adata.var, 'highly_variable') and adata.var['highly_variable'].any():
                        features = adata.var_names[adata.var['highly_variable']][:100].tolist()
                    else:
                        features = adata.var_names[:100].tolist()
                
                if analysis_type == 'autocorrelation_moran':
                    result = self.compute_spatial_autocorrelation(
                        adata, genes=features, spatial_key=spatial_key,
                        method='moran', n_neighbors=n_neighbors
                    )
                    batch_results['results']['autocorrelation_moran'] = result
                    
                elif analysis_type == 'autocorrelation_geary':
                    result = self.compute_spatial_autocorrelation(
                        adata, genes=features, spatial_key=spatial_key,
                        method='geary', n_neighbors=n_neighbors
                    )
                    batch_results['results']['autocorrelation_geary'] = result
                    
                elif analysis_type == 'local_moran':
                    result = self.local_spatial_statistics(
                        adata, genes=features, spatial_key=spatial_key,
                        method='local_moran', n_neighbors=n_neighbors
                    )
                    batch_results['results']['local_moran'] = result
                    
                elif analysis_type == 'getis_ord':
                    result = self.local_spatial_statistics(
                        adata, genes=features, spatial_key=spatial_key,
                        method='getis_ord', n_neighbors=n_neighbors
                    )
                    batch_results['results']['getis_ord'] = result
                    
                elif analysis_type == 'bivariate_moran':
                    # Create gene pairs from features
                    if len(features) >= 2:
                        gene_pairs = [(features[i], features[j]) 
                                    for i in range(len(features)//2) 
                                    for j in range(i+1, min(i+3, len(features)))][:10]  # Limit pairs
                        
                        # Filter kwargs to only include valid parameters for bivariate_moran
                        valid_params = {'analysis_type'}
                        filtered_kwargs = {k: v for k, v in kwargs.items() if k in valid_params}
                        
                        result = self.bivariate_moran(
                            adata, gene_pairs=gene_pairs, spatial_key=spatial_key,
                            n_neighbors=n_neighbors, **filtered_kwargs
                        )
                        batch_results['results']['bivariate_moran'] = result
                    else:
                        logger.warning("Need at least 2 features for bivariate analysis")
                        
                elif analysis_type == 'neighborhood_enrichment':
                    # Requires cluster information
                    cluster_key = kwargs.get('cluster_key', 'leiden')
                    if cluster_key in adata.obs.columns:
                        result = self.neighborhood_enrichment(
                            adata, cluster_key=cluster_key, spatial_key=spatial_key,
                            n_neighbors=n_neighbors, **kwargs
                        )
                        batch_results['results']['neighborhood_enrichment'] = result
                    else:
                        logger.warning(f"Cluster key '{cluster_key}' not found in adata.obs")
                        
                elif analysis_type == 'point_patterns':
                    # Requires cell type information
                    cell_type_key = kwargs.get('cell_type_key', 'leiden')
                    if cell_type_key in adata.obs.columns:
                        result = self.spatial_point_patterns(
                            adata, cell_type_key=cell_type_key, spatial_key=spatial_key,
                            **kwargs
                        )
                        batch_results['results']['point_patterns'] = result
                    else:
                        logger.warning(f"Cell type key '{cell_type_key}' not found in adata.obs")
                
                # Record timing
                end_time = datetime.datetime.now()
                duration = (end_time - start_time).total_seconds()
                batch_results['performance_metrics']['analysis_times'][analysis_type] = duration
                logger.info(f"Completed {analysis_type} analysis in {duration:.2f}s")
                
            except Exception as e:
                logger.error(f"Failed to run {analysis_type} analysis: {str(e)}")
                batch_results['results'][f'{analysis_type}_error'] = str(e)
        
        # Record total time
        batch_results['performance_metrics']['end_time'] = datetime.datetime.now().isoformat()
        total_time = sum(batch_results['performance_metrics']['analysis_times'].values())
        batch_results['performance_metrics']['total_analysis_time'] = total_time
        
        # Integrate results into adata
        integration_key = self.integrate_analysis_results(adata, analysis_name="batch_analysis")
        batch_results['integration_key'] = integration_key
        
        logger.info(f"Batch spatial analysis completed in {total_time:.2f}s")
        logger.info(f"Successfully completed {len(batch_results['results'])} analyses")
        
        return batch_results
    
    def prepare_visualization_data(
        self,
        adata: ad.AnnData,
        result_key: Optional[str] = None,
        feature_types: List[str] = None,
        output_format: str = "dict"
    ) -> Dict[str, Any]:
        """
        Format spatial statistics results for direct visualization.
        
        This method creates visualization-ready data structures from spatial
        analysis results, handling different data types and providing
        standardized formats for plotting libraries.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix with spatial statistics results
        result_key : Optional[str]
            Specific integration key to format (if None, searches for latest)
        feature_types : List[str]
            Types of features to include: ['moran', 'getis', 'bivariate']
        output_format : str
            Output format: 'dict', 'dataframe', or 'plotly'
            
        Returns
        -------
        Dict[str, Any]
            Visualization-ready data structure
        """
        # Find integration metadata
        if result_key is None:
            # Find most recent integration
            integration_keys = [key for key in adata.uns.keys() if key.endswith('_integration_metadata')]
            if not integration_keys:
                raise ValueError("No spatial analysis results found. Run spatial analysis first.")
            result_key = integration_keys[-1].replace('_integration_metadata', '')
        
        metadata_key = f"{result_key}_integration_metadata"
        if metadata_key not in adata.uns:
            raise ValueError(f"Integration metadata not found for key: {result_key}")
        
        metadata = adata.uns[metadata_key]
        
        # Initialize visualization data
        viz_data = {
            'metadata': metadata,
            'spatial_coordinates': adata.obsm['spatial'],
            'features': {},
            'color_maps': {},
            'data_types': {},
            'plot_configs': {}
        }
        
        # Set default feature types if not specified
        if feature_types is None:
            feature_types = metadata.get('methods_used', [])
        
        # Process each feature type
        for method in feature_types:
            if method == 'local_moran':
                viz_data['features']['local_moran'] = self._prepare_local_moran_viz(adata, result_key, metadata)
                viz_data['color_maps']['local_moran'] = 'RdBu_r'  # Diverging colormap for Moran's I
                viz_data['data_types']['local_moran'] = 'continuous'
                
            elif method == 'getis_ord':
                viz_data['features']['getis_ord'] = self._prepare_getis_ord_viz(adata, result_key, metadata)
                viz_data['color_maps']['getis_ord'] = 'RdYlBu_r'  # Hot-cold colormap for G*
                viz_data['data_types']['getis_ord'] = 'continuous'
                
            elif method == 'bivariate_moran':
                viz_data['features']['bivariate_moran'] = self._prepare_bivariate_moran_viz(adata, result_key, metadata)
                viz_data['color_maps']['bivariate_moran'] = 'RdBu_r'
                viz_data['data_types']['bivariate_moran'] = 'continuous'
                
            elif method == 'neighborhood_enrichment':
                viz_data['features']['neighborhood_enrichment'] = self._prepare_enrichment_viz(adata, metadata)
                viz_data['color_maps']['neighborhood_enrichment'] = 'RdBu_r'
                viz_data['data_types']['neighborhood_enrichment'] = 'matrix'
        
        # Add plot configuration suggestions
        viz_data['plot_configs'] = self._generate_plot_configs(viz_data)
        
        # Format output based on requested format
        if output_format == "dataframe":
            return self._format_as_dataframes(viz_data)
        elif output_format == "plotly":
            return self._format_for_plotly(viz_data)
        else:
            return viz_data
    
    def _prepare_local_moran_viz(self, adata: ad.AnnData, result_key: str, metadata: Dict) -> Dict[str, Any]:
        """Prepare local Moran's I results for visualization."""
        viz_features = {}
        
        if 'moran_genes' in metadata.get('local_statistics', {}):
            for gene in metadata['local_statistics']['moran_genes']:
                i_col = f"{result_key}_local_moran_{gene}_I"
                pval_col = f"{result_key}_local_moran_{gene}_pval"
                cluster_col = f"{result_key}_local_moran_{gene}_cluster"
                
                feature_data = {}
                if i_col in adata.obs.columns:
                    feature_data['values'] = adata.obs[i_col].values
                    feature_data['value_type'] = 'moran_I'
                    feature_data['range'] = [float(adata.obs[i_col].min()), float(adata.obs[i_col].max())]
                
                if pval_col in adata.obs.columns:
                    feature_data['p_values'] = adata.obs[pval_col].values
                    feature_data['significant_mask'] = adata.obs[pval_col] < 0.05
                
                if cluster_col in adata.obs.columns:
                    feature_data['clusters'] = adata.obs[cluster_col].values
                    feature_data['cluster_labels'] = {
                        0: 'Non-significant',
                        1: 'High-High',
                        2: 'Low-Low', 
                        3: 'Low-High',
                        4: 'High-Low'
                    }
                
                viz_features[gene] = feature_data
                
        return viz_features
    
    def _prepare_getis_ord_viz(self, adata: ad.AnnData, result_key: str, metadata: Dict) -> Dict[str, Any]:
        """Prepare Getis-Ord G* results for visualization."""
        viz_features = {}
        
        if 'getis_genes' in metadata.get('local_statistics', {}):
            for gene in metadata['local_statistics']['getis_genes']:
                gi_col = f"{result_key}_getis_ord_{gene}_Gi"
                pval_col = f"{result_key}_getis_ord_{gene}_pval"
                zscore_col = f"{result_key}_getis_ord_{gene}_zscore"
                
                feature_data = {}
                if gi_col in adata.obs.columns:
                    feature_data['values'] = adata.obs[gi_col].values
                    feature_data['value_type'] = 'getis_ord_Gi'
                    feature_data['range'] = [float(adata.obs[gi_col].min()), float(adata.obs[gi_col].max())]
                
                if pval_col in adata.obs.columns:
                    feature_data['p_values'] = adata.obs[pval_col].values
                    feature_data['significant_mask'] = adata.obs[pval_col] < 0.05
                
                if zscore_col in adata.obs.columns:
                    feature_data['z_scores'] = adata.obs[zscore_col].values
                    # Create hotspot/coldspot labels
                    z_vals = adata.obs[zscore_col].values
                    feature_data['hotspot_labels'] = np.select(
                        [z_vals > 1.96, z_vals < -1.96], 
                        ['Hot Spot', 'Cold Spot'], 
                        default='Not Significant'
                    )
                
                viz_features[gene] = feature_data
                
        return viz_features
    
    def _prepare_bivariate_moran_viz(self, adata: ad.AnnData, result_key: str, metadata: Dict) -> Dict[str, Any]:
        """Prepare bivariate Moran's I results for visualization."""
        viz_features = {}
        
        if 'bivariate_pairs' in metadata.get('local_statistics', {}):
            for pair in metadata['local_statistics']['bivariate_pairs']:
                i_col = f"{result_key}_bivariate_moran_{pair}_I"
                pval_col = f"{result_key}_bivariate_moran_{pair}_pval"
                
                feature_data = {}
                if i_col in adata.obs.columns:
                    feature_data['values'] = adata.obs[i_col].values
                    feature_data['value_type'] = 'bivariate_moran_I'
                    feature_data['range'] = [float(adata.obs[i_col].min()), float(adata.obs[i_col].max())]
                    feature_data['pair_name'] = pair
                
                if pval_col in adata.obs.columns:
                    feature_data['p_values'] = adata.obs[pval_col].values
                    feature_data['significant_mask'] = adata.obs[pval_col] < 0.05
                
                viz_features[pair] = feature_data
                
        return viz_features
    
    def _prepare_enrichment_viz(self, adata: ad.AnnData, metadata: Dict) -> Dict[str, Any]:
        """Prepare neighborhood enrichment results for visualization."""
        viz_data = {}
        
        if 'neighborhood_enrichment' in metadata.get('result_keys', {}):
            enrichment_key = metadata['result_keys']['neighborhood_enrichment']
            if enrichment_key in adata.uns:
                enrichment_result = adata.uns[enrichment_key]
                viz_data['zscore_matrix'] = enrichment_result['zscore']
                viz_data['count_matrix'] = enrichment_result.get('count', None)
                viz_data['cluster_names'] = list(enrichment_result['zscore'].index)
                viz_data['data_type'] = 'heatmap'
                
        return viz_data
    
    def _generate_plot_configs(self, viz_data: Dict) -> Dict[str, Any]:
        """Generate plot configuration suggestions."""
        configs = {}
        
        for feature_type, features in viz_data['features'].items():
            if feature_type in ['local_moran', 'getis_ord', 'bivariate_moran']:
                configs[feature_type] = {
                    'plot_type': 'spatial_scatter',
                    'color_by': 'values',
                    'size': 6,
                    'alpha': 0.8,
                    'colormap': viz_data['color_maps'][feature_type],
                    'show_colorbar': True,
                    'title_template': f'{feature_type.title()} Statistics for {{feature}}'
                }
            elif feature_type == 'neighborhood_enrichment':
                configs[feature_type] = {
                    'plot_type': 'heatmap',
                    'colormap': viz_data['color_maps'][feature_type],
                    'show_colorbar': True,
                    'title': 'Neighborhood Enrichment Z-scores'
                }
        
        return configs
    
    def _format_as_dataframes(self, viz_data: Dict) -> Dict[str, Any]:
        """Format visualization data as pandas DataFrames."""
        formatted_data = {
            'metadata': viz_data['metadata'],
            'coordinates': pd.DataFrame(
                viz_data['spatial_coordinates'], 
                columns=['x', 'y']
            ),
            'features': {}
        }
        
        for feature_type, features in viz_data['features'].items():
            if feature_type in ['local_moran', 'getis_ord', 'bivariate_moran']:
                feature_dfs = {}
                for feature_name, feature_data in features.items():
                    df_data = {'x': viz_data['spatial_coordinates'][:, 0],
                              'y': viz_data['spatial_coordinates'][:, 1]}
                    df_data.update(feature_data)
                    feature_dfs[feature_name] = pd.DataFrame(df_data)
                formatted_data['features'][feature_type] = feature_dfs
            else:
                formatted_data['features'][feature_type] = features
                
        return formatted_data
    
    def _format_for_plotly(self, viz_data: Dict) -> Dict[str, Any]:
        """Format visualization data for Plotly visualization."""
        plotly_data = {
            'metadata': viz_data['metadata'],
            'figures': {}
        }
        
        coords = viz_data['spatial_coordinates']
        
        for feature_type, features in viz_data['features'].items():
            if feature_type in ['local_moran', 'getis_ord', 'bivariate_moran']:
                plotly_data['figures'][feature_type] = []
                for feature_name, feature_data in features.items():
                    fig_data = {
                        'type': 'scatter',
                        'mode': 'markers',
                        'x': coords[:, 0].tolist(),
                        'y': coords[:, 1].tolist(),
                        'marker': {
                            'color': feature_data['values'].tolist(),
                            'colorscale': viz_data['color_maps'][feature_type],
                            'showscale': True,
                            'colorbar': {'title': feature_data['value_type']}
                        },
                        'name': feature_name,
                        'text': [f"{feature_name}: {val:.3f}" for val in feature_data['values']]
                    }
                    plotly_data['figures'][feature_type].append(fig_data)
                    
        return plotly_data
    
    def generate_analysis_summary(
        self,
        adata: ad.AnnData,
        result_key: Optional[str] = None,
        output_format: str = "dict",
        include_interpretation: bool = True
    ) -> Union[Dict[str, Any], str]:
        """
        Generate comprehensive summary of spatial analysis results.
        
        This method produces human-readable reports including statistical
        significance, effect sizes, and interpretation hints for spatial
        analysis results.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix with spatial statistics results
        result_key : Optional[str]
            Specific integration key to summarize (if None, uses latest)
        output_format : str
            Output format: 'dict', 'markdown', 'html', or 'json'
        include_interpretation : bool
            Whether to include interpretation guidelines
            
        Returns
        -------
        Union[Dict[str, Any], str]
            Analysis summary in requested format
        """
        # Find integration metadata
        if result_key is None:
            integration_keys = [key for key in adata.uns.keys() if key.endswith('_integration_metadata')]
            if not integration_keys:
                raise ValueError("No spatial analysis results found")
            result_key = integration_keys[-1].replace('_integration_metadata', '')
        
        metadata_key = f"{result_key}_integration_metadata"
        metadata = adata.uns[metadata_key]
        
        # Generate comprehensive summary
        summary = {
            'analysis_overview': {
                'analysis_name': metadata['analysis_name'],
                'timestamp': metadata['timestamp'],
                'dataset_info': {
                    'n_observations': metadata['n_obs'],
                    'n_variables': metadata['n_vars']
                },
                'methods_used': metadata['methods_used'],
                'total_features_analyzed': 0
            },
            'method_summaries': {},
            'significant_findings': {},
            'interpretation_guide': {} if include_interpretation else None
        }
        
        # Summarize each method
        for method in metadata['methods_used']:
            method_summary = self._summarize_method_results(adata, result_key, method, metadata)
            summary['method_summaries'][method] = method_summary
            summary['analysis_overview']['total_features_analyzed'] += method_summary.get('n_features', 0)
            
            # Extract significant findings
            if method_summary.get('significant_features'):
                summary['significant_findings'][method] = {
                    'n_significant': len(method_summary['significant_features']),
                    'top_features': method_summary['significant_features'][:5],
                    'effect_sizes': method_summary.get('effect_size_summary', {})
                }
        
        # Add interpretation guide
        if include_interpretation:
            summary['interpretation_guide'] = self._generate_interpretation_guide(metadata['methods_used'])
        
        # Format output
        if output_format == "markdown":
            return self._format_as_markdown(summary)
        elif output_format == "html":
            return self._format_as_html(summary)
        elif output_format == "json":
            return json.dumps(summary, indent=2, default=str)
        else:
            return summary
    
    def _summarize_method_results(self, adata: ad.AnnData, result_key: str, method: str, metadata: Dict) -> Dict[str, Any]:
        """Summarize results for a specific method."""
        method_summary = {'method': method, 'n_features': 0, 'significant_features': []}
        
        if method == 'local_moran':
            if 'moran_genes' in metadata.get('local_statistics', {}):
                genes = metadata['local_statistics']['moran_genes']
                method_summary['n_features'] = len(genes)
                
                significant_genes = []
                effect_sizes = {}
                
                for gene in genes:
                    i_col = f"{result_key}_local_moran_{gene}_I"
                    pval_col = f"{result_key}_local_moran_{gene}_pval"
                    
                    if i_col in adata.obs.columns:
                        i_values = adata.obs[i_col].values
                        effect_sizes[gene] = {
                            'mean_I': float(np.mean(i_values)),
                            'std_I': float(np.std(i_values)),
                            'max_I': float(np.max(i_values)),
                            'min_I': float(np.min(i_values))
                        }
                        
                        if pval_col in adata.obs.columns:
                            n_significant = np.sum(adata.obs[pval_col] < 0.05)
                            effect_sizes[gene]['n_significant_spots'] = int(n_significant)
                            effect_sizes[gene]['proportion_significant'] = float(n_significant / len(i_values))
                            
                            if n_significant > len(i_values) * 0.05:  # More than expected by chance
                                significant_genes.append(gene)
                
                method_summary['significant_features'] = sorted(significant_genes, 
                    key=lambda g: effect_sizes[g]['proportion_significant'], reverse=True)
                method_summary['effect_size_summary'] = effect_sizes
                
        elif method == 'getis_ord':
            if 'getis_genes' in metadata.get('local_statistics', {}):
                genes = metadata['local_statistics']['getis_genes']
                method_summary['n_features'] = len(genes)
                
                significant_genes = []
                effect_sizes = {}
                
                for gene in genes:
                    gi_col = f"{result_key}_getis_ord_{gene}_Gi"
                    pval_col = f"{result_key}_getis_ord_{gene}_pval"
                    
                    if gi_col in adata.obs.columns:
                        gi_values = adata.obs[gi_col].values
                        effect_sizes[gene] = {
                            'mean_Gi': float(np.mean(gi_values)),
                            'std_Gi': float(np.std(gi_values)),
                            'max_Gi': float(np.max(gi_values)),
                            'min_Gi': float(np.min(gi_values))
                        }
                        
                        if pval_col in adata.obs.columns:
                            n_significant = np.sum(adata.obs[pval_col] < 0.05)
                            effect_sizes[gene]['n_significant_spots'] = int(n_significant)
                            effect_sizes[gene]['proportion_significant'] = float(n_significant / len(gi_values))
                            
                            if n_significant > len(gi_values) * 0.05:
                                significant_genes.append(gene)
                
                method_summary['significant_features'] = sorted(significant_genes,
                    key=lambda g: effect_sizes[g]['proportion_significant'], reverse=True)
                method_summary['effect_size_summary'] = effect_sizes
        
        elif method == 'neighborhood_enrichment':
            if 'neighborhood_enrichment' in metadata.get('result_keys', {}):
                enrichment_key = metadata['result_keys']['neighborhood_enrichment']
                if enrichment_key in adata.uns:
                    zscore_matrix = adata.uns[enrichment_key]['zscore']
                    method_summary['n_features'] = zscore_matrix.shape[0] * zscore_matrix.shape[1]
                    
                    # Find significant enrichments (|z| > 2)
                    significant_pairs = []
                    for i, cluster1 in enumerate(zscore_matrix.index):
                        for j, cluster2 in enumerate(zscore_matrix.columns):
                            z_score = zscore_matrix.iloc[i, j]
                            if abs(z_score) > 2:
                                significant_pairs.append(f"{cluster1}-{cluster2}")
                    
                    method_summary['significant_features'] = significant_pairs
                    method_summary['enrichment_matrix_shape'] = zscore_matrix.shape
        
        return method_summary
    
    def _generate_interpretation_guide(self, methods_used: List[str]) -> Dict[str, str]:
        """Generate interpretation guide for spatial analysis methods."""
        guide = {}
        
        if 'local_moran' in methods_used:
            guide['local_moran'] = (
                "Local Moran's I identifies spatial clusters and outliers. "
                "Positive values indicate spatial clustering (similar values nearby), "
                "negative values indicate spatial dispersion (dissimilar values nearby). "
                "Clusters: HH (high-high), LL (low-low), HL (high-low outlier), LH (low-high outlier)."
            )
        
        if 'getis_ord' in methods_used:
            guide['getis_ord'] = (
                "Getis-Ord G* identifies hot spots (high values surrounded by high values) "
                "and cold spots (low values surrounded by low values). "
                "Positive G* values indicate hot spots, negative values indicate cold spots. "
                "Significance determined by p-values < 0.05."
            )
        
        if 'bivariate_moran' in methods_used:
            guide['bivariate_moran'] = (
                "Bivariate Moran's I measures spatial correlation between two different variables. "
                "Positive values indicate that high values of one variable tend to be near "
                "high values of another variable spatially."
            )
        
        if 'neighborhood_enrichment' in methods_used:
            guide['neighborhood_enrichment'] = (
                "Neighborhood enrichment analysis identifies whether cell types are more "
                "or less likely to be neighbors than expected by chance. "
                "Z-scores > 2 indicate significant enrichment, Z-scores < -2 indicate depletion."
            )
        
        return guide
    
    def _format_as_markdown(self, summary: Dict[str, Any]) -> str:
        """Format summary as markdown string."""
        md = f"""# Spatial Analysis Summary

## Analysis Overview
- **Analysis Name**: {summary['analysis_overview']['analysis_name']}
- **Timestamp**: {summary['analysis_overview']['timestamp']}
- **Dataset**: {summary['analysis_overview']['dataset_info']['n_observations']} observations, {summary['analysis_overview']['dataset_info']['n_variables']} variables
- **Methods Used**: {', '.join(summary['analysis_overview']['methods_used'])}
- **Total Features Analyzed**: {summary['analysis_overview']['total_features_analyzed']}

"""

        # Add method summaries
        md += "## Method Summaries\n\n"
        for method, method_summary in summary['method_summaries'].items():
            md += f"### {method.replace('_', ' ').title()}\n"
            md += f"- Features analyzed: {method_summary['n_features']}\n"
            if method_summary['significant_features']:
                md += f"- Significant features: {len(method_summary['significant_features'])}\n"
                md += f"- Top significant features: {', '.join(method_summary['significant_features'][:5])}\n"
            md += "\n"
        
        # Add significant findings
        if summary['significant_findings']:
            md += "## Significant Findings\n\n"
            for method, findings in summary['significant_findings'].items():
                md += f"### {method.replace('_', ' ').title()}\n"
                md += f"- Number of significant features: {findings['n_significant']}\n"
                md += f"- Top features: {', '.join(findings['top_features'])}\n\n"
        
        # Add interpretation guide
        if summary['interpretation_guide']:
            md += "## Interpretation Guide\n\n"
            for method, interpretation in summary['interpretation_guide'].items():
                md += f"### {method.replace('_', ' ').title()}\n"
                md += f"{interpretation}\n\n"
        
        return md
    
    def _format_as_html(self, summary: Dict[str, Any]) -> str:
        """Format summary as HTML string."""
        html = f"""
        <html>
        <head><title>Spatial Analysis Summary</title></head>
        <body>
        <h1>Spatial Analysis Summary</h1>
        
        <h2>Analysis Overview</h2>
        <ul>
            <li><strong>Analysis Name:</strong> {summary['analysis_overview']['analysis_name']}</li>
            <li><strong>Timestamp:</strong> {summary['analysis_overview']['timestamp']}</li>
            <li><strong>Dataset:</strong> {summary['analysis_overview']['dataset_info']['n_observations']} observations, {summary['analysis_overview']['dataset_info']['n_variables']} variables</li>
            <li><strong>Methods Used:</strong> {', '.join(summary['analysis_overview']['methods_used'])}</li>
            <li><strong>Total Features Analyzed:</strong> {summary['analysis_overview']['total_features_analyzed']}</li>
        </ul>
        """
        
        # Add method summaries
        html += "<h2>Method Summaries</h2>"
        for method, method_summary in summary['method_summaries'].items():
            html += f"<h3>{method.replace('_', ' ').title()}</h3><ul>"
            html += f"<li>Features analyzed: {method_summary['n_features']}</li>"
            if method_summary['significant_features']:
                html += f"<li>Significant features: {len(method_summary['significant_features'])}</li>"
                html += f"<li>Top significant features: {', '.join(method_summary['significant_features'][:5])}</li>"
            html += "</ul>"
        
        # Add interpretation guide if present
        if summary['interpretation_guide']:
            html += "<h2>Interpretation Guide</h2>"
            for method, interpretation in summary['interpretation_guide'].items():
                html += f"<h3>{method.replace('_', ' ').title()}</h3><p>{interpretation}</p>"
        
        html += "</body></html>"
        return html
    
    def join_count_analysis(
        self,
        adata: ad.AnnData,
        category_key: str,
        spatial_key: str = 'spatial',
        n_neighbors: int = 6,
        case_values: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Perform join count analysis for categorical data.
        
        This method analyzes spatial autocorrelation of categorical variables
        such as cell types or clusters using join count statistics.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix with categorical annotations
        category_key : str
            Key in adata.obs containing categorical data (cell types/clusters)
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
        n_neighbors : int
            Number of spatial neighbors to consider
        case_values : Optional[List[str]]
            Specific category values to analyze. If None, analyzes all categories
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing join count statistics for each category
        """
        try:
            from esda.join_counts import Join_Counts
            from libpysal.weights import KNN
        except ImportError as e:
            logger.error(f"Required packages not installed: {e}")
            logger.error("Install with: pip install esda libpysal")
            raise ImportError("Please install esda and libpysal: pip install esda libpysal")
        
        if category_key not in adata.obs.columns:
            raise ValueError(f"Category key '{category_key}' not found in adata.obs")
        
        # Get spatial coordinates
        coords = adata.obsm[spatial_key]
        
        # Create spatial weights
        w = KNN.from_array(coords, k=n_neighbors)
        w.transform = 'B'  # Binary weights for join counts
        
        # Get categorical data
        categories = adata.obs[category_key].astype('category')
        
        # Determine which categories to analyze
        if case_values is None:
            case_values = categories.cat.categories.tolist()
        
        results = {}
        
        for case_value in case_values:
            if case_value not in categories.cat.categories:
                logger.warning(f"Category '{case_value}' not found in {category_key}")
                continue
            
            # Create binary indicator (1 for case, 0 for non-case)
            y = (categories == case_value).astype(int).values
            
            # Skip if all values are the same (no variation)
            if len(np.unique(y)) < 2:
                logger.warning(f"Category '{case_value}' has no variation. Skipping.")
                continue
            
            try:
                # Perform join count analysis
                jc = Join_Counts(y, w)
                
                # Store results - use lowercase attributes for current esda API
                results[case_value] = {
                    'BB': jc.bb,  # Black-Black joins (case-case)
                    'BW': jc.bw,  # Black-White joins (case-non-case)
                    'WW': jc.ww,  # White-White joins (non-case-non-case)
                    'J': jc.J,    # Total joins
                    'mean_BB': jc.mean_bb if hasattr(jc, 'mean_bb') else jc.mean_BB if hasattr(jc, 'mean_BB') else None,  # Expected BB under null
                    'var_BB': jc.var_bb if hasattr(jc, 'var_bb') else jc.var_BB if hasattr(jc, 'var_BB') else None,    # Variance of BB
                    'z_BB': jc.z_sim_bb if hasattr(jc, 'z_sim_bb') else jc.z_sim_BB if hasattr(jc, 'z_sim_BB') else None,  # Z-score for BB
                    'p_BB': jc.p_sim_bb if hasattr(jc, 'p_sim_bb') else jc.p_sim_BB if hasattr(jc, 'p_sim_BB') else None,  # P-value for BB
                    'n_case': int(y.sum()),
                    'n_total': len(y),
                    'case_proportion': float(y.mean())
                }
                
            except Exception as e:
                logger.warning(f"Join count analysis failed for category '{case_value}': {str(e)}")
                continue
        
        # Store results in adata.uns for downstream use
        adata.uns[f'{category_key}_join_counts'] = results
        
        logger.info(f"Join count analysis completed for {len(results)} categories")
        
        return results
    
    def spatial_centrality_analysis(
        self,
        adata: ad.AnnData,
        cluster_key: str,
        spatial_key: str = 'spatial',
        n_neighbors: int = 6,
        metrics: List[str] = ['degree', 'closeness', 'betweenness']
    ) -> Dict[str, Any]:
        """
        Calculate graph centrality measures for spatial networks.
        
        This method computes various centrality measures for the spatial graph,
        providing insights into the structural importance of different locations.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        cluster_key : str
            Key in adata.obs containing cluster/cell type information
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
        n_neighbors : int
            Number of spatial neighbors to consider
        metrics : List[str]
            Centrality metrics to compute: 'degree', 'closeness', 'betweenness'
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing centrality scores and cluster-level aggregations
        """
        try:
            import squidpy as sq
            logger.debug(f"Using squidpy version: {sq.__version__}")
        except ImportError as e:
            logger.error(f"squidpy not installed: {e}")
            logger.error("Install with: pip install squidpy")
            raise ImportError("Please install squidpy: pip install squidpy")
        
        if cluster_key not in adata.obs.columns:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
        
        # Create a copy to avoid modifying original data
        adata_copy = adata.copy()
        
        # Ensure cluster_key is categorical (required by squidpy)
        if not adata_copy.obs[cluster_key].dtype.name == 'category':
            adata_copy.obs[cluster_key] = adata_copy.obs[cluster_key].astype('category')
        
        # Build spatial graph
        sq.gr.spatial_neighbors(
            adata_copy,
            coord_type='generic',
            spatial_key=spatial_key,
            n_neighs=n_neighbors,
            set_diag=False
        )
        
        # Compute centrality scores - use 'score' parameter instead of 'metrics'
        try:
            sq.gr.centrality_scores(
                adata_copy,
                cluster_key=cluster_key,
                score=metrics
            )
        except Exception as e:
            logger.error(f"Squidpy centrality analysis failed: {str(e)}")
            # Fall back to manual centrality calculation
            return self._manual_centrality_analysis(adata_copy, cluster_key, spatial_key, n_neighbors, metrics)
        
        # Extract results
        results = {
            'node_centrality': {},
            'cluster_centrality': {}
        }
        
        # Node-level centrality scores
        for metric in metrics:
            centrality_key = f'{cluster_key}_{metric}_centrality'
            if centrality_key in adata_copy.obs.columns:
                results['node_centrality'][metric] = adata_copy.obs[centrality_key].values
                # Write back to original adata
                adata.obs[centrality_key] = adata_copy.obs[centrality_key]
            else:
                logger.warning(f"Centrality metric '{metric}' not found in results")
        
        # Cluster-level aggregations
        clusters = adata.obs[cluster_key].unique()
        for metric in metrics:
            if metric in results['node_centrality']:
                cluster_means = {}
                cluster_stds = {}
                
                for cluster in clusters:
                    mask = adata.obs[cluster_key] == cluster
                    values = results['node_centrality'][metric][mask]
                    cluster_means[str(cluster)] = float(np.mean(values))
                    cluster_stds[str(cluster)] = float(np.std(values))
                
                results['cluster_centrality'][f'{metric}_mean'] = cluster_means
                results['cluster_centrality'][f'{metric}_std'] = cluster_stds
        
        # Store results in adata.uns
        adata.uns[f'{cluster_key}_centrality_analysis'] = results
        
        logger.info(f"Centrality analysis completed for {len(metrics)} metrics across {len(clusters)} clusters")
        
        return results
    
    def _manual_centrality_analysis(
        self,
        adata: ad.AnnData,
        cluster_key: str,
        spatial_key: str,
        n_neighbors: int,
        metrics: List[str]
    ) -> Dict[str, Any]:
        """
        Manual centrality calculation fallback using networkx.
        """
        try:
            import networkx as nx
            from libpysal.weights import KNN
        except ImportError:
            logger.error("networkx not installed. Install with: pip install networkx")
            raise ImportError("Please install networkx: pip install networkx")
        
        # Get spatial coordinates
        coords = adata.obsm[spatial_key]
        
        # Create spatial weights
        w = KNN.from_array(coords, k=n_neighbors)
        
        # Convert to networkx graph
        G = nx.from_scipy_sparse_array(w.sparse)
        
        results = {
            'node_centrality': {},
            'cluster_centrality': {}
        }
        
        # Calculate centrality metrics
        for metric in metrics:
            if metric == 'degree':
                centrality = dict(G.degree())
                centrality_values = np.array([centrality[i] for i in range(len(adata))])
            elif metric == 'closeness':
                centrality = nx.closeness_centrality(G)
                centrality_values = np.array([centrality[i] for i in range(len(adata))])
            elif metric == 'betweenness':
                centrality = nx.betweenness_centrality(G)
                centrality_values = np.array([centrality[i] for i in range(len(adata))])
            else:
                logger.warning(f"Unknown centrality metric: {metric}")
                continue
            
            # Store node-level results
            results['node_centrality'][metric] = centrality_values
            
            # Write back to adata
            adata.obs[f'{cluster_key}_{metric}_centrality'] = centrality_values
        
        # Cluster-level aggregations
        clusters = adata.obs[cluster_key].unique()
        for metric in metrics:
            if metric in results['node_centrality']:
                cluster_means = {}
                cluster_stds = {}
                
                for cluster in clusters:
                    mask = adata.obs[cluster_key] == cluster
                    values = results['node_centrality'][metric][mask]
                    cluster_means[str(cluster)] = float(np.mean(values))
                    cluster_stds[str(cluster)] = float(np.std(values))
                
                results['cluster_centrality'][f'{metric}_mean'] = cluster_means
                results['cluster_centrality'][f'{metric}_std'] = cluster_stds
        
        # Store results in adata.uns
        adata.uns[f'{cluster_key}_centrality_analysis'] = results
        
        return results
    
    def network_properties_analysis(
        self,
        adata: ad.AnnData,
        spatial_key: str = 'spatial',
        cluster_key: Optional[str] = None,
        n_neighbors: int = 6
    ) -> Dict[str, Any]:
        """
        Calculate global network properties of the spatial graph.
        
        This method computes various graph-level metrics that characterize
        the overall structure and connectivity of the spatial network.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
        cluster_key : Optional[str]
            Key in adata.obs containing cluster information (for modularity)
        n_neighbors : int
            Number of spatial neighbors to consider
            
        Returns
        -------
        Dict[str, Any]
            Dictionary containing global network properties
        """
        try:
            import networkx as nx
            from libpysal.weights import KNN
        except ImportError:
            logger.error("networkx not installed. Install with: pip install networkx")
            raise ImportError("Please install networkx: pip install networkx")
        
        # Get spatial coordinates
        coords = adata.obsm[spatial_key]
        
        # Create spatial weights
        w = KNN.from_array(coords, k=n_neighbors)
        
        # Convert to networkx graph
        G = nx.from_scipy_sparse_array(w.sparse)
        
        results = {
            'n_nodes': G.number_of_nodes(),
            'n_edges': G.number_of_edges(),
            'density': nx.density(G),
            'is_connected': nx.is_connected(G)
        }
        
        # Clustering coefficient
        try:
            results['avg_clustering'] = nx.average_clustering(G)
            results['clustering_coefficients'] = dict(nx.clustering(G))
        except Exception as e:
            logger.warning(f"Clustering coefficient calculation failed: {str(e)}")
            results['avg_clustering'] = None
        
        # Assortativity (degree correlation)
        try:
            results['degree_assortativity'] = nx.degree_assortativity_coefficient(G)
        except Exception as e:
            logger.warning(f"Degree assortativity calculation failed: {str(e)}")
            results['degree_assortativity'] = None
        
        # Global efficiency
        try:
            results['global_efficiency'] = nx.global_efficiency(G)
        except Exception as e:
            logger.warning(f"Global efficiency calculation failed: {str(e)}")
            results['global_efficiency'] = None
        
        # Average path length (only for connected graphs)
        if results['is_connected']:
            try:
                results['avg_path_length'] = nx.average_shortest_path_length(G)
                results['diameter'] = nx.diameter(G)
                results['radius'] = nx.radius(G)
            except Exception as e:
                logger.warning(f"Path length calculations failed: {str(e)}")
                results['avg_path_length'] = None
                results['diameter'] = None
                results['radius'] = None
        else:
            # For disconnected graphs, compute for largest component
            largest_cc = max(nx.connected_components(G), key=len)
            G_cc = G.subgraph(largest_cc)
            results['n_components'] = nx.number_connected_components(G)
            results['largest_component_size'] = len(largest_cc)
            results['largest_component_fraction'] = len(largest_cc) / G.number_of_nodes()
            
            try:
                results['avg_path_length_largest_cc'] = nx.average_shortest_path_length(G_cc)
                results['diameter_largest_cc'] = nx.diameter(G_cc)
                results['radius_largest_cc'] = nx.radius(G_cc)
            except Exception as e:
                logger.warning(f"Largest component path calculations failed: {str(e)}")
        
        # Modularity (if cluster information is available)
        if cluster_key and cluster_key in adata.obs.columns:
            try:
                # Create partition from cluster labels
                clusters = adata.obs[cluster_key].astype('category')
                partition = {}
                for i, cluster_value in enumerate(clusters):
                    partition[i] = clusters.cat.codes[i]
                
                # Calculate modularity
                results['modularity'] = nx.community.modularity(G, 
                    [set([i for i, c in partition.items() if c == cluster_id]) 
                     for cluster_id in np.unique(list(partition.values()))])
                
            except Exception as e:
                logger.warning(f"Modularity calculation failed: {str(e)}")
                results['modularity'] = None
        
        # Transitivity (alternative to clustering coefficient)
        try:
            results['transitivity'] = nx.transitivity(G)
        except Exception as e:
            logger.warning(f"Transitivity calculation failed: {str(e)}")
            results['transitivity'] = None
        
        # Degree distribution statistics
        degrees = dict(G.degree())
        degree_values = list(degrees.values())
        results['degree_stats'] = {
            'mean': float(np.mean(degree_values)),
            'std': float(np.std(degree_values)),
            'min': int(np.min(degree_values)),
            'max': int(np.max(degree_values)),
            'median': float(np.median(degree_values))
        }
        
        # Store results in adata.uns
        storage_key = f'network_properties_{cluster_key}' if cluster_key else 'network_properties'
        adata.uns[storage_key] = results
        
        logger.info(f"Network properties analysis completed. Graph has {results['n_nodes']} nodes and {results['n_edges']} edges")
        
        return results

# Note: Spatial variable gene identification functions have been moved to spatial_genes.py
# Use the find_spatial_genes function from the main server interface instead

def compute_spatial_autocorrelation(
    adata: ad.AnnData,
    genes: Optional[List[str]] = None,
    **kwargs
) -> pd.DataFrame:
    """
    Compute spatial autocorrelation.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix
    genes : Optional[List[str]]
        Genes to analyze
    **kwargs
        Additional parameters
        
    Returns
    -------
    pd.DataFrame
        Autocorrelation statistics
    """
    stats_tool = SpatialStatistics()
    return stats_tool.compute_spatial_autocorrelation(adata, genes=genes, **kwargs)


async def calculate_spatial_stats(
    data_id: str,
    data_store: Dict[str, Any],
    feature: str,
    statistic: str = "morans_i",
    n_neighbors: int = 6,
    context = None
) -> Dict[str, Any]:
    """
    Calculate spatial statistics for a single feature/gene.
    
    Parameters
    ----------
    data_id : str
        Dataset identifier
    data_store : Dict[str, Any]
        Data storage dictionary
    feature : str
        Feature/gene name to analyze
    statistic : str
        Type of statistic to compute
    n_neighbors : int
        Number of spatial neighbors
    context : optional
        MCP context for logging
        
    Returns
    -------
    Dict[str, Any]
        Statistics result
    """
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")
    
    adata = data_store[data_id]["adata"]
    
    if feature not in adata.var_names:
        raise ValueError(f"Feature '{feature}' not found in dataset")
    
    if context:
        await context.info(f"Computing {statistic} for {feature}")
    
    # Map statistic names to methods
    statistic_map = {
        "morans_i": "moran",
        "gearys_c": "geary",
        "local_morans": "local_moran"
    }
    
    method = statistic_map.get(statistic, "moran")
    
    try:
        stats_tool = SpatialStatistics()
        
        if statistic == "local_morans":
            # Local statistics
            result = stats_tool.local_spatial_statistics(
                adata=adata,
                genes=[feature],
                spatial_key='spatial',
                method='local_moran',
                n_neighbors=n_neighbors
            )
            
            return {
                "feature": feature,
                "statistic": statistic,
                "local_values": result[feature]['Is'].tolist(),
                "p_values": result[feature]['p_values'].tolist(),
                "clusters": result[feature]['clusters'].tolist() if 'clusters' in result[feature] else None,
                "n_neighbors": n_neighbors
            }
        else:
            # Global statistics
            result_df = stats_tool.compute_spatial_autocorrelation(
                adata=adata,
                genes=[feature],
                spatial_key='spatial',
                method=method,
                n_neighbors=n_neighbors
            )
            
            if len(result_df) > 0:
                row = result_df.iloc[0]
                
                if method == "moran":
                    return {
                        "feature": feature,
                        "statistic": statistic,
                        "value": float(row['moran_I']),
                        "expected": float(row['expected_I']),
                        "variance": float(row['variance']),
                        "z_score": float(row['z_score']),
                        "p_value": float(row['p_value']),
                        "q_value": float(row['q_value']) if 'q_value' in row else None,
                        "n_neighbors": n_neighbors
                    }
                else:  # geary
                    return {
                        "feature": feature,
                        "statistic": statistic,
                        "value": float(row['geary_C']),
                        "expected": float(row['expected_C']),
                        "variance": float(row['variance']),
                        "z_score": float(row['z_score']),
                        "p_value": float(row['p_value']),
                        "q_value": float(row['q_value']) if 'q_value' in row else None,
                        "n_neighbors": n_neighbors
                    }
            else:
                raise ValueError(f"No statistics computed for {feature}")
                
    except ImportError as e:
        if context:
            await context.info(f"Required package not installed: {e}")
        raise ImportError("Please install pysal: pip install pysal")
    except Exception as e:
        if context:
            await context.info(f"Error computing statistics: {e}")
        raise


def neighborhood_enrichment_analysis(
    adata: ad.AnnData,
    cluster_key: str = 'leiden',
    **kwargs
) -> pd.DataFrame:
    """
    Convenience function for neighborhood enrichment analysis.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix
    cluster_key : str
        Key in adata.obs containing cluster/cell type information
    **kwargs
        Additional parameters passed to SpatialStatistics.neighborhood_enrichment()
        
    Returns
    -------
    pd.DataFrame
        Neighborhood enrichment z-scores matrix
    """
    stats_tool = SpatialStatistics()
    return stats_tool.neighborhood_enrichment(adata, cluster_key=cluster_key, **kwargs)


def bivariate_moran_analysis(
    adata: ad.AnnData,
    gene_pairs: List[Tuple[str, str]],
    **kwargs
) -> Dict[str, Any]:
    """
    Convenience function for bivariate Moran's I analysis.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix
    gene_pairs : List[Tuple[str, str]]
        List of gene pairs to analyze as (gene1, gene2) tuples
    **kwargs
        Additional parameters passed to SpatialStatistics.bivariate_moran()
        
    Returns
    -------
    Dict[str, Any]
        Dictionary containing global and/or local bivariate Moran's I results
    """
    stats_tool = SpatialStatistics()
    return stats_tool.bivariate_moran(adata, gene_pairs=gene_pairs, **kwargs)


def integrate_spatial_analysis_results(
    adata: ad.AnnData,
    analysis_name: str = "spatial_analysis",
    timestamp: Optional[str] = None
) -> str:
    """
    Convenience function for integrating spatial analysis results.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix with spatial statistics results
    analysis_name : str
        Name for this analysis session
    timestamp : Optional[str]
        Timestamp string, if None current time will be used
        
    Returns
    -------
    str
        Key prefix used for storing integrated results
    """
    stats_tool = SpatialStatistics()
    return stats_tool.integrate_analysis_results(adata, analysis_name, timestamp)


def batch_spatial_analysis(
    adata: ad.AnnData,
    analysis_types: List[str],
    feature_lists: Optional[Dict[str, List[str]]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Convenience function for running batch spatial analysis.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix
    analysis_types : List[str]
        List of analysis types to run
    feature_lists : Optional[Dict[str, List[str]]]
        Dictionary mapping analysis types to specific gene lists
    **kwargs
        Additional parameters
        
    Returns
    -------
    Dict[str, Any]
        Comprehensive results dictionary with all analysis results
    """
    stats_tool = SpatialStatistics()
    return stats_tool.batch_spatial_analysis(adata, analysis_types, feature_lists, **kwargs)


def prepare_spatial_visualization_data(
    adata: ad.AnnData,
    result_key: Optional[str] = None,
    feature_types: List[str] = None,
    output_format: str = "dict"
) -> Dict[str, Any]:
    """
    Convenience function for preparing visualization data.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix with spatial statistics results
    result_key : Optional[str]
        Specific integration key to format (if None, searches for latest)
    feature_types : List[str]
        Types of features to include: ['moran', 'getis', 'bivariate']
    output_format : str
        Output format: 'dict', 'dataframe', or 'plotly'
        
    Returns
    -------
    Dict[str, Any]
        Visualization-ready data structure
    """
    stats_tool = SpatialStatistics()
    return stats_tool.prepare_visualization_data(adata, result_key, feature_types, output_format)


def generate_spatial_analysis_summary(
    adata: ad.AnnData,
    result_key: Optional[str] = None,
    output_format: str = "dict",
    include_interpretation: bool = True
) -> Union[Dict[str, Any], str]:
    """
    Convenience function for generating analysis summary.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix with spatial statistics results
    result_key : Optional[str]
        Specific integration key to summarize (if None, uses latest)
    output_format : str
        Output format: 'dict', 'markdown', 'html', or 'json'
    include_interpretation : bool
        Whether to include interpretation guidelines
        
    Returns
    -------
    Union[Dict[str, Any], str]
        Analysis summary in requested format
    """
    stats_tool = SpatialStatistics()
    return stats_tool.generate_analysis_summary(adata, result_key, output_format, include_interpretation)