"""
Spatial Statistics Tool

This module provides spatial statistical methods for analyzing spatial transcriptomics data.
Focuses on spatial autocorrelation, local indicators of spatial association (LISA),
and spatial point pattern analysis.

Note: Spatial variable gene identification has been moved to spatial_genes.py
"""

import logging
from typing import List, Optional, Dict, Any, Union, Tuple
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import stats
from sklearn.neighbors import NearestNeighbors
import warnings

logger = logging.getLogger(__name__)

class SpatialStatistics:
    """
    A class for performing spatial statistical analysis.

    This class provides methods for:
    - Spatial autocorrelation analysis (Moran's I, Geary's C)
    - Local indicators of spatial association (LISA)
    - Spatial point pattern analysis (Ripley's K)

    Note: Spatial variable gene identification methods have been moved to spatial_genes.py
    """

    def __init__(self):
        """Initialize the SpatialStatistics tool."""
        pass
    


    
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
        
        # Create spatial weights
        w = KNN.from_array(coords, k=n_neighbors)
        
        # Select genes
        if genes is None:
            genes = adata.var_names[:100]  # Default to top 100 HVGs
            logger.info(f"Computing autocorrelation for top 100 genes")
        
        results = []
        
        for gene in genes:
            if gene not in adata.var_names:
                continue
            
            # Get expression values
            expr = adata[:, gene].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()
            
            if method == 'moran':
                mi = Moran(expr, w)
                results.append({
                    'gene': gene,
                    'moran_I': mi.I,
                    'expected_I': mi.EI,
                    'variance': mi.VI,
                    'z_score': mi.z_norm,
                    'p_value': mi.p_norm
                })
            elif method == 'geary':
                c = Geary(expr, w)
                results.append({
                    'gene': gene,
                    'geary_C': c.C,
                    'expected_C': c.EC,
                    'variance': c.VC,
                    'z_score': c.z_norm,
                    'p_value': c.p_norm
                })
        
        results_df = pd.DataFrame(results)
        
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
        
        # Create spatial weights
        w = KNN.from_array(coords, k=n_neighbors)
        
        results = {}
        
        for gene in genes:
            if gene not in adata.var_names:
                logger.warning(f"Gene {gene} not found")
                continue
            
            # Get expression values
            expr = adata[:, gene].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()
            else:
                expr = expr.flatten()
            
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
                        
                        results[cell_type] = {
                            'support': k_support,
                            'K': k_statistic,
                            'L': l_statistic
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
        
        return {
            'n_cells': len(ct_coords),
            'mean_nn_distance': float(nn_distances.mean()),
            'std_nn_distance': float(nn_distances.std()),
            'density': len(ct_coords) / self._compute_area(all_coords)
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
            
            results[cell_type] = {
                'n_cells': len(ct_coords),
                'mean_nn_distance': float(nn_distances.mean()),
                'std_nn_distance': float(nn_distances.std()),
                'density': len(ct_coords) / self._compute_area(coords)
            }
        
        return results
    
    def _compute_area(self, coords: np.ndarray) -> float:
        """Compute bounding box area of coordinates."""
        min_coords = coords.min(axis=0)
        max_coords = coords.max(axis=0)
        return np.prod(max_coords - min_coords)

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