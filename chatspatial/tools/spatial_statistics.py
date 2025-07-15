"""
Spatial Statistics Tool

This module provides advanced spatial statistical methods for analyzing spatial transcriptomics data.
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
    A class for performing advanced spatial statistical analysis.
    
    This class provides methods for:
    - Spatial variable gene detection (SpatialDE, SPARK)
    - Spatial autocorrelation analysis
    - Spatial point pattern analysis
    - Local indicators of spatial association
    """
    
    def __init__(self):
        """Initialize the SpatialStatistics tool."""
        self.method_info = {
            'spatialDE': {
                'name': 'SpatialDE',
                'description': 'Gaussian process-based spatial gene expression analysis',
                'reference': 'Svensson et al., Nature Methods 2018',
                'type': 'python'
            },
            'spark': {
                'name': 'SPARK',
                'description': 'Spatial Pattern Recognition via Kernels',
                'reference': 'Sun et al., Nature Methods 2020',
                'type': 'R'
            },
            'somde': {
                'name': 'SOMDE',
                'description': 'Spatial Omnibus test for Differential Expression',
                'reference': 'Hao et al., Genome Biology 2021',
                'type': 'python'
            }
        }
    
    def find_spatial_genes(
        self,
        adata: ad.AnnData,
        method: str = 'spatialDE',
        n_genes: Optional[int] = None,
        spatial_key: str = 'spatial',
        **kwargs
    ) -> pd.DataFrame:
        """
        Find spatially variable genes.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        method : str
            Method to use ('spatialDE', 'spark', 'somde')
        n_genes : Optional[int]
            Number of top genes to return
        spatial_key : str
            Key in obsm containing spatial coordinates
        **kwargs
            Method-specific parameters
            
        Returns
        -------
        pd.DataFrame
            DataFrame with spatial statistics for each gene
        """
        if method not in self.method_info:
            raise ValueError(f"Unknown method: {method}. Available: {list(self.method_info.keys())}")
        
        logger.info(f"Finding spatial genes using {method}")
        
        if method == 'spatialDE':
            return self._run_spatialDE(adata, n_genes, spatial_key, **kwargs)
        elif method == 'spark':
            return self._run_spark(adata, n_genes, spatial_key, **kwargs)
        elif method == 'somde':
            return self._run_somde(adata, n_genes, spatial_key, **kwargs)
        else:
            raise NotImplementedError(f"Method {method} not implemented")
    
    def _run_spatialDE(
        self,
        adata: ad.AnnData,
        n_genes: Optional[int],
        spatial_key: str,
        normalized: bool = True,
        **kwargs
    ) -> pd.DataFrame:
        """
        Run SpatialDE for spatial gene detection.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        n_genes : Optional[int]
            Number of top genes to return
        spatial_key : str
            Key in obsm containing spatial coordinates
        normalized : bool
            Whether data is already normalized
        **kwargs
            Additional SpatialDE parameters
            
        Returns
        -------
        pd.DataFrame
            SpatialDE results
        """
        try:
            import SpatialDE
        except ImportError:
            logger.error("SpatialDE not installed. Install with: pip install spatialde")
            raise ImportError("Please install SpatialDE: pip install spatialde")
        
        # Prepare data
        coords = pd.DataFrame(
            adata.obsm[spatial_key][:, :2],  # Ensure 2D coordinates
            columns=['x', 'y'],
            index=adata.obs_names
        )
        
        # Get expression data
        if normalized:
            # Use normalized data
            if 'log1p' in adata.uns_keys():
                counts = pd.DataFrame(
                    adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
                    columns=adata.var_names,
                    index=adata.obs_names
                )
            else:
                logger.warning("Data may not be log-normalized")
                counts = pd.DataFrame(
                    adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
                    columns=adata.var_names,
                    index=adata.obs_names
                )
        else:
            # Use raw counts and normalize
            if adata.raw is not None:
                raw_counts = pd.DataFrame(
                    adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X,
                    columns=adata.raw.var_names,
                    index=adata.obs_names
                )
            else:
                raw_counts = pd.DataFrame(
                    adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
                    columns=adata.var_names,
                    index=adata.obs_names
                )
            
            # Normalize using standard approach (similar to scanpy)
            # Total count normalization and log1p transformation
            total_counts = raw_counts.sum(axis=1)
            norm_counts = raw_counts.div(total_counts, axis=0) * np.median(total_counts)
            counts = np.log1p(norm_counts)
        
        # Run SpatialDE
        logger.info("Running SpatialDE analysis")
        # SpatialDE expects numpy arrays for coordinates
        results = SpatialDE.run(coords.values, counts, **kwargs)
        
        # Multiple testing correction
        from SpatialDE.util import qvalue
        results['qval'] = qvalue(results['pval'].values, pi0=0.1)
        
        # Sort by q-value
        results = results.sort_values('qval')
        
        # Get top genes if requested
        if n_genes is not None:
            results = results.head(n_genes)
        
        # Add to adata
        adata.var['spatial_pval'] = results.set_index('g')['pval']
        adata.var['spatial_qval'] = results.set_index('g')['qval']
        adata.var['spatial_l'] = results.set_index('g')['l']
        
        return results
    
    def _run_spark(
        self,
        adata: ad.AnnData,
        n_genes: Optional[int],
        spatial_key: str,
        **kwargs
    ) -> pd.DataFrame:
        """
        Run SPARK using rpy2.
        
        Parameters
        ----------
        adata : ad.AnnData
            Annotated data matrix
        n_genes : Optional[int]
            Number of top genes
        spatial_key : str
            Spatial coordinate key
        **kwargs
            SPARK parameters
            
        Returns
        -------
        pd.DataFrame
            SPARK results
        """
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.packages import importr
            pandas2ri.activate()
        except ImportError:
            logger.error("rpy2 not installed. Install with: pip install rpy2")
            raise ImportError("Please install rpy2 for R integration")
        
        logger.info("Running SPARK through R")
        
        # Check if SPARK is installed in R
        try:
            spark = importr('SPARK')
        except:
            logger.error("SPARK not installed in R. Install with: install.packages('SPARK')")
            raise ImportError("Please install SPARK in R")
        
        # Prepare data for R - SPARK has specific requirements
        # 1. Coordinates must be a data.frame with row names matching column names of counts
        coords = pd.DataFrame(
            adata.obsm[spatial_key], 
            columns=['x', 'y'],
            index=adata.obs_names
        )
        
        # 2. Counts must be genes x cells matrix
        counts = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
            columns=adata.var_names,
            index=adata.obs_names
        ).T  # SPARK expects genes x cells
        
        # 3. Calculate library sizes (total counts per cell)
        lib_sizes = counts.sum(axis=0).values  # Sum for each cell
        
        # Convert to R objects
        r_coords = pandas2ri.py2rpy(coords)
        r_counts = pandas2ri.py2rpy(counts)
        r_lib_sizes = ro.FloatVector(lib_sizes)
        
        # Run SPARK with correct parameters
        ro.r('''
        run_spark <- function(counts, coords, lib_sizes) {
            library(SPARK)
            
            # Create SPARK object
            spark_obj <- CreateSPARKObject(counts=counts, 
                                         location=coords,
                                         percentage=0.1,
                                         min_total_counts=10)
            
            # Run SPARK analysis with proper lib_size vector
            spark_obj <- spark.vc(spark_obj, 
                                covariates=NULL, 
                                lib_size=lib_sizes,  # Use actual library sizes
                                num_core=1,
                                verbose=FALSE)
                                
            spark_obj <- spark.test(spark_obj, 
                                   check_positive=TRUE, 
                                   verbose=FALSE)
                                   
            return(spark_obj@res_mtest)
        }
        ''')
        
        run_spark = ro.r['run_spark']
        results = run_spark(r_counts, r_coords, r_lib_sizes)
        
        # Convert back to pandas
        results_df = pandas2ri.rpy2py(results)
        results_df = results_df.reset_index().rename(columns={'index': 'gene'})
        
        # Sort by adjusted p-value
        results_df = results_df.sort_values('adjusted_pvalue')
        
        if n_genes is not None:
            results_df = results_df.head(n_genes)
        
        # Add to adata
        for _, row in results_df.iterrows():
            gene = row['gene']
            if gene in adata.var_names:
                adata.var.loc[gene, 'spark_pval'] = row['combined_pvalue']
                adata.var.loc[gene, 'spark_qval'] = row['adjusted_pvalue']
        
        return results_df
    
    def _run_somde(
        self,
        adata: ad.AnnData,
        n_genes: Optional[int],
        spatial_key: str,
        **kwargs
    ) -> pd.DataFrame:
        """
        Run SOMDE for spatial gene detection.
        
        Note: This is a placeholder for SOMDE integration.
        """
        logger.warning("SOMDE integration not yet implemented")
        raise NotImplementedError("SOMDE integration coming soon")
    
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

# Convenience functions
def find_spatial_variable_genes(
    adata: ad.AnnData,
    method: str = 'spatialDE',
    **kwargs
) -> pd.DataFrame:
    """
    Find spatially variable genes.
    
    Parameters
    ----------
    adata : ad.AnnData
        Annotated data matrix
    method : str
        Method to use
    **kwargs
        Method-specific parameters
        
    Returns
    -------
    pd.DataFrame
        Spatial gene statistics
    """
    stats_tool = SpatialStatistics()
    return stats_tool.find_spatial_genes(adata, method=method, **kwargs)

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