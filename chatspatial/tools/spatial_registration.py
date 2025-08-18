"""
Spatial Registration Tool

This module provides functionality for aligning and registering multiple spatial transcriptomics slices.
"""

import logging
from typing import List, Optional, Dict, Any, Union, Tuple
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from pathlib import Path

logger = logging.getLogger(__name__)

class SpatialRegistration:
    """
    A class for performing spatial registration of multiple tissue slices.
    
    This class provides methods for:
    - Pairwise alignment of spatial transcriptomics slices
    - Multi-slice integration into a common coordinate system
    - Spatial coordinate transformation
    """
    
    def __init__(self):
        """Initialize the SpatialRegistration tool."""
        self.method_info = {
            'paste': {
                'name': 'PASTE',
                'description': 'Probabilistic Alignment of Spatial Transcriptomics Experiments',
                'reference': 'Zeira et al., Nature Methods 2022',
                'type': 'python'
            },
            'stalign': {
                'name': 'STalign',
                'description': 'Spatial transcriptomics alignment using image registration',
                'reference': 'Clifton et al., bioRxiv 2023',
                'type': 'python'
            }
        }
    
    def register_slices(
        self,
        adata_list: List[ad.AnnData],
        method: str = 'paste',
        reference_idx: Optional[int] = None,
        **kwargs
    ) -> List[ad.AnnData]:
        """
        Register multiple spatial transcriptomics slices.
        
        Parameters
        ----------
        adata_list : List[ad.AnnData]
            List of AnnData objects to register
        method : str
            Registration method to use ('paste', 'stalign')
        reference_idx : Optional[int]
            Index of reference slice (if None, use first slice)
        **kwargs
            Additional method-specific parameters
            
        Returns
        -------
        List[ad.AnnData]
            List of registered AnnData objects with aligned coordinates
        """
        if method not in self.method_info:
            raise ValueError(f"Unknown method: {method}. Available methods: {list(self.method_info.keys())}")
        
        logger.info(f"Performing spatial registration using {method}")
        
        if method == 'paste':
            return self._register_with_paste(adata_list, reference_idx, **kwargs)
        elif method == 'stalign':
            return self._register_with_stalign(adata_list, reference_idx, **kwargs)
        else:
            raise NotImplementedError(f"Method {method} not yet implemented")
    
    def _register_with_paste(
        self,
        adata_list: List[ad.AnnData],
        reference_idx: Optional[int] = None,
        alpha: float = 0.1,
        n_components: int = 30,
        use_gpu: bool = False,
        **kwargs
    ) -> List[ad.AnnData]:
        """
        Register slices using PASTE algorithm.
        
        Parameters
        ----------
        adata_list : List[ad.AnnData]
            List of AnnData objects to register
        reference_idx : Optional[int]
            Index of reference slice
        alpha : float
            Spatial regularization parameter
        n_components : int
            Number of components for dimensionality reduction
        use_gpu : bool
            Whether to use GPU acceleration
        **kwargs
            Additional PASTE parameters
            
        Returns
        -------
        List[ad.AnnData]
            Registered AnnData objects
        """
        try:
            import paste as pst
        except ImportError:
            logger.error("PASTE not installed. Install with: pip install paste-bio")
            raise ImportError("Please install paste-bio: pip install paste-bio")
        
        if reference_idx is None:
            reference_idx = 0
        
        # Prepare data
        logger.info("Preparing data for PASTE registration")
        registered_list = [adata.copy() for adata in adata_list]
        
        # Ensure all slices have spatial coordinates
        for i, adata in enumerate(registered_list):
            if 'spatial' not in adata.obsm:
                raise ValueError(f"Slice {i} missing spatial coordinates in obsm['spatial']")
        
        # Perform pairwise alignment
        if len(registered_list) == 2:
            logger.info("Performing pairwise alignment")
            
            # Get common genes
            common_genes = list(set(registered_list[0].var_names) & set(registered_list[1].var_names))
            logger.info(f"Using {len(common_genes)} common genes")
            
            slice1 = registered_list[0][:, common_genes]
            slice2 = registered_list[1][:, common_genes]
            
            # Create backend object
            import ot
            if use_gpu:
                try:
                    import torch
                    if torch.cuda.is_available():
                        backend = ot.backend.TorchBackend()
                    else:
                        backend = ot.backend.NumpyBackend()
                except ImportError:
                    backend = ot.backend.NumpyBackend()
            else:
                backend = ot.backend.NumpyBackend()
            
            # Filter out parameters that will be passed explicitly
            filtered_kwargs = {k: v for k, v in kwargs.items() 
                             if k not in ['alpha', 'backend', 'use_gpu', 'verbose', 'gpu_verbose']}
            
            # Run PASTE
            pi = pst.pairwise_align(
                slice1, slice2,
                alpha=alpha,
                backend=backend,
                use_gpu=use_gpu,
                verbose=False,
                gpu_verbose=False,
                **filtered_kwargs
            )
            
            # Apply transformation
            aligned_slices = pst.stack_slices_pairwise([slice1, slice2], [pi])
            
            # Extract the aligned spatial coordinates
            registered_list[0].obsm['spatial_registered'] = aligned_slices[0].obsm['spatial']
            registered_list[1].obsm['spatial_registered'] = aligned_slices[1].obsm['spatial']
            
        else:
            # Multi-slice registration using center alignment
            logger.info(f"Performing multi-slice registration with {len(registered_list)} slices")
            
            # Get common genes
            common_genes = set(registered_list[0].var_names)
            for adata in registered_list[1:]:
                common_genes = common_genes & set(adata.var_names)
            common_genes = list(common_genes)
            logger.info(f"Using {len(common_genes)} common genes")
            
            # Subset to common genes
            slices = [adata[:, common_genes] for adata in registered_list]
            
            # Initial pairwise alignments to reference
            pis = []
            
            # Create backend object
            import ot
            if use_gpu:
                try:
                    import torch
                    if torch.cuda.is_available():
                        backend = ot.backend.TorchBackend()
                    else:
                        backend = ot.backend.NumpyBackend()
                except ImportError:
                    backend = ot.backend.NumpyBackend()
            else:
                backend = ot.backend.NumpyBackend()
            
            for i, slice_data in enumerate(slices):
                if i == reference_idx:
                    # Identity mapping for reference slice
                    n = slices[i].shape[0]
                    pi = np.eye(n)
                else:
                    # Filter out parameters that will be passed explicitly
                    filtered_kwargs = {k: v for k, v in kwargs.items() 
                                     if k not in ['alpha', 'backend', 'use_gpu', 'verbose', 'gpu_verbose']}
                    pi = pst.pairwise_align(
                        slices[reference_idx], slice_data,
                        alpha=alpha,
                        backend=backend,
                        use_gpu=use_gpu,
                        verbose=False,
                        gpu_verbose=False,
                        **filtered_kwargs
                    )
                pis.append(pi)
            
            # Center alignment
            center_slice, pis_new = pst.center_align(
                slices[reference_idx], slices,
                pis_init=pis,
                alpha=alpha,
                backend=backend,
                use_gpu=use_gpu,
                n_components=n_components,
                verbose=False,
                gpu_verbose=False
            )
            
            # Apply transformations
            for i, (adata, pi) in enumerate(zip(registered_list, pis_new)):
                if i == reference_idx:
                    adata.obsm['spatial_registered'] = adata.obsm['spatial'].copy()
                else:
                    # Transform coordinates based on alignment
                    coords = adata.obsm['spatial']
                    # Apply transformation from pi
                    adata.obsm['spatial_registered'] = self._transform_coordinates(
                        coords, pi, slices[reference_idx].obsm['spatial']
                    )
        
        logger.info("PASTE registration completed")
        return registered_list
    
    def _register_with_stalign(
        self,
        adata_list: List[ad.AnnData],
        reference_idx: Optional[int] = None,
        **kwargs
    ) -> List[ad.AnnData]:
        """
        Register slices using STalign algorithm.
        
        Note: This is a placeholder for STalign integration.
        """
        logger.warning("STalign integration not yet implemented")
        raise NotImplementedError("STalign integration coming soon")
    
    def _transform_coordinates(
        self,
        coords: np.ndarray,
        transport_matrix: np.ndarray,
        reference_coords: np.ndarray
    ) -> np.ndarray:
        """
        Transform coordinates based on optimal transport matrix.
        
        Parameters
        ----------
        coords : np.ndarray
            Original coordinates
        transport_matrix : np.ndarray
            Optimal transport matrix from PASTE
        reference_coords : np.ndarray
            Reference slice coordinates
            
        Returns
        -------
        np.ndarray
            Transformed coordinates
        """
        # Normalize transport matrix
        transport_matrix = transport_matrix / transport_matrix.sum(axis=1, keepdims=True)
        
        # Compute weighted average of reference coordinates
        transformed = transport_matrix @ reference_coords
        
        return transformed
    
    def compute_registration_quality(
        self,
        adata_list: List[ad.AnnData],
        use_registered: bool = True
    ) -> Dict[str, float]:
        """
        Compute quality metrics for registration.
        
        Parameters
        ----------
        adata_list : List[ad.AnnData]
            List of registered AnnData objects
        use_registered : bool
            Whether to use registered coordinates
            
        Returns
        -------
        Dict[str, float]
            Dictionary of quality metrics
        """
        coord_key = 'spatial_registered' if use_registered else 'spatial'
        
        metrics = {}
        
        # Compute pairwise distances between slices
        for i in range(len(adata_list) - 1):
            coords1 = adata_list[i].obsm[coord_key]
            coords2 = adata_list[i + 1].obsm[coord_key]
            
            # Compute mean distance to nearest neighbor
            from sklearn.neighbors import NearestNeighbors
            nbrs = NearestNeighbors(n_neighbors=1).fit(coords2)
            distances, _ = nbrs.kneighbors(coords1)
            
            metrics[f'mean_nn_distance_{i}_{i+1}'] = float(distances.mean())
        
        # Compute overlap scores if cell types are available
        if all('cell_type' in adata.obs for adata in adata_list):
            for i in range(len(adata_list) - 1):
                ct1 = adata_list[i].obs['cell_type'].value_counts(normalize=True)
                ct2 = adata_list[i + 1].obs['cell_type'].value_counts(normalize=True)
                
                # Compute overlap coefficient
                common_types = set(ct1.index) & set(ct2.index)
                overlap = sum(min(ct1[ct], ct2[ct]) for ct in common_types)
                metrics[f'cell_type_overlap_{i}_{i+1}'] = float(overlap)
        
        return metrics
    
    def visualize_registration(
        self,
        adata_list: List[ad.AnnData],
        use_registered: bool = True,
        save: Optional[str] = None,
        **kwargs
    ):
        """
        Visualize registration results.
        
        Parameters
        ----------
        adata_list : List[ad.AnnData]
            List of AnnData objects
        use_registered : bool
            Whether to use registered coordinates
        save : Optional[str]
            Path to save figure
        **kwargs
            Additional plotting parameters
        """
        import matplotlib.pyplot as plt
        
        coord_key = 'spatial_registered' if use_registered else 'spatial'
        n_slices = len(adata_list)
        
        fig, axes = plt.subplots(1, n_slices, figsize=(5*n_slices, 5))
        if n_slices == 1:
            axes = [axes]
        
        for i, (adata, ax) in enumerate(zip(adata_list, axes)):
            coords = adata.obsm[coord_key]
            
            # Plot with colors if available
            if 'cell_type' in adata.obs:
                import seaborn as sns
                palette = sns.color_palette('tab20', n_colors=adata.obs['cell_type'].nunique())
                color_map = dict(zip(adata.obs['cell_type'].unique(), palette))
                colors = [color_map[ct] for ct in adata.obs['cell_type']]
                ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=10, alpha=0.6)
            else:
                ax.scatter(coords[:, 0], coords[:, 1], s=10, alpha=0.6)
            
            ax.set_title(f'Slice {i+1}')
            ax.set_aspect('equal')
            ax.axis('off')
        
        plt.suptitle('Registered' if use_registered else 'Original')
        plt.tight_layout()
        
        if save:
            plt.savefig(save, dpi=300, bbox_inches='tight')
        plt.show()

# Create convenience functions for direct use
def register_spatial_slices(
    adata_list: List[ad.AnnData],
    method: str = 'paste',
    **kwargs
) -> List[ad.AnnData]:
    """
    Register multiple spatial transcriptomics slices.
    
    Parameters
    ----------
    adata_list : List[ad.AnnData]
        List of AnnData objects to register
    method : str
        Registration method ('paste', 'stalign')
    **kwargs
        Method-specific parameters
        
    Returns
    -------
    List[ad.AnnData]
        Registered AnnData objects
    """
    registration_tool = SpatialRegistration()
    return registration_tool.register_slices(adata_list, method=method, **kwargs)

def evaluate_registration(
    adata_list: List[ad.AnnData],
    use_registered: bool = True
) -> Dict[str, float]:
    """
    Evaluate registration quality.
    
    Parameters
    ----------
    adata_list : List[ad.AnnData]
        List of registered AnnData objects
    use_registered : bool
        Whether to evaluate registered coordinates
        
    Returns
    -------
    Dict[str, float]
        Registration quality metrics
    """
    registration_tool = SpatialRegistration()
    return registration_tool.compute_registration_quality(adata_list, use_registered)


# Standard architecture-compliant wrapper function for MCP server
async def register_spatial_slices_mcp(
    source_id: str,
    target_id: str,
    data_store: dict,
    method: str = "paste",
    landmarks: Optional[dict] = None,
    context = None
) -> dict:
    """
    Register spatial slices following the standard MCP tool architecture.
    
    This function follows the standard pattern used by other tools:
    - Accepts data IDs and data_store dictionary
    - Extracts AnnData objects internally
    - Returns results in standard format
    
    Args:
        source_id: Source dataset ID
        target_id: Target dataset ID to align to
        data_store: Dictionary containing dataset information
        method: Registration method (paste, stalign)
        landmarks: Additional parameters (currently unused)
        context: MCP context for logging
        
    Returns:
        Dictionary with registration results
    """
    if context:
        await context.info(f"Registering {source_id} to {target_id} using method: {method}")
    
    # Validate data exists
    if source_id not in data_store:
        raise ValueError(f"Source dataset {source_id} not found in data store")
    if target_id not in data_store:
        raise ValueError(f"Target dataset {target_id} not found in data store")
    
    # Extract AnnData objects
    source_adata = data_store[source_id]["adata"].copy()
    target_adata = data_store[target_id]["adata"].copy()
    adata_list = [source_adata, target_adata]
    
    try:
        # Call the core registration function
        registered_adata_list = register_spatial_slices(
            adata_list, 
            method=method
        )
        
        # Update data store with registered data
        data_store[source_id]["adata"] = registered_adata_list[0]  # source (registered)
        data_store[target_id]["adata"] = registered_adata_list[1]  # target (reference)
        
        # Create result dictionary
        result = {
            "method": method,
            "source_id": source_id,
            "target_id": target_id,
            "n_source_spots": registered_adata_list[0].n_obs,
            "n_target_spots": registered_adata_list[1].n_obs,
            "registration_completed": True,
            "spatial_key_registered": "spatial_registered"
        }
        
        if context:
            await context.info(f"Registration completed. Registered coordinates stored in 'spatial_registered' key.")
        
        return result
        
    except Exception as e:
        error_msg = f"Registration failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e