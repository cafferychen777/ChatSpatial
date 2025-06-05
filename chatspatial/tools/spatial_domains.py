"""
Spatial domain identification tools for spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
import traceback

from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from ..models.data import SpatialDomainParameters
from ..models.analysis import SpatialDomainResult
from ..utils.image_utils import fig_to_image, create_placeholder_image


# Import optional dependencies with fallbacks
try:
    import STAGATE as stg
    STAGATE_AVAILABLE = True
except ImportError:
    STAGATE_AVAILABLE = False

try:
    import SpaGCN as spg
    SPAGCN_AVAILABLE = True
except ImportError:
    SPAGCN_AVAILABLE = False


async def identify_spatial_domains(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialDomainParameters = SpatialDomainParameters(),
    context: Optional[Context] = None
) -> SpatialDomainResult:
    """Identify spatial domains in spatial transcriptomics data
    
    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Spatial domain identification parameters
        context: MCP context
        
    Returns:
        Spatial domain identification result
    """
    if context:
        await context.info(f"Identifying spatial domains using {params.method} method")
    
    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")
    
    adata = data_store[data_id]["adata"].copy()
    
    try:
        # Check if spatial coordinates exist
        if 'spatial' not in adata.obsm and not any('spatial' in key for key in adata.obsm.keys()):
            raise ValueError("No spatial coordinates found in the dataset")
        
        # Prepare data for domain identification
        if context:
            await context.info("Preparing data for spatial domain identification...")
        
        # Use highly variable genes if requested and available
        if params.use_highly_variable and 'highly_variable' in adata.var.columns:
            adata_subset = adata[:, adata.var['highly_variable']].copy()
            if context:
                await context.info(f"Using {sum(adata.var['highly_variable'])} highly variable genes")
        else:
            adata_subset = adata.copy()
            if context:
                await context.info(f"Using all {adata.n_vars} genes")
        
        # Ensure data is properly normalized and log-transformed
        if 'log1p' not in adata_subset.uns:
            if context:
                await context.info("Normalizing and log-transforming data...")
            sc.pp.normalize_total(adata_subset, target_sum=1e4)
            sc.pp.log1p(adata_subset)
        
        # Identify domains based on method
        if params.method == "stagate":
            domain_labels, embeddings_key, statistics = await _identify_domains_stagate(
                adata_subset, params, context
            )
        elif params.method == "spagcn":
            domain_labels, embeddings_key, statistics = await _identify_domains_spagcn(
                adata_subset, params, context
            )
        elif params.method in ["leiden", "louvain"]:
            domain_labels, embeddings_key, statistics = await _identify_domains_clustering(
                adata_subset, params, context
            )
        else:
            raise ValueError(f"Unsupported method: {params.method}")
        
        # Store domain labels in original adata
        domain_key = f"spatial_domains_{params.method}"
        adata.obs[domain_key] = domain_labels
        adata.obs[domain_key] = adata.obs[domain_key].astype('category')
        
        # Store embeddings if available
        if embeddings_key and embeddings_key in adata_subset.obsm:
            adata.obsm[embeddings_key] = adata_subset.obsm[embeddings_key]
        
        # Refine domains if requested
        refined_domain_key = None
        if params.refine_domains:
            if context:
                await context.info("Refining spatial domains using spatial smoothing...")
            refined_domain_key = f"{domain_key}_refined"
            refined_labels = _refine_spatial_domains(adata, domain_key, refined_domain_key)
            adata.obs[refined_domain_key] = refined_labels
            adata.obs[refined_domain_key] = adata.obs[refined_domain_key].astype('category')
        
        # Get domain counts
        domain_counts = adata.obs[domain_key].value_counts().to_dict()
        domain_counts = {str(k): int(v) for k, v in domain_counts.items()}
        
        # Create visualization if requested
        visualization = None
        if params.include_image:
            if context:
                await context.info("Creating spatial domain visualization...")
            visualization = _create_domain_visualization(
                adata, 
                domain_key if not refined_domain_key else refined_domain_key,
                params
            )
        
        # Update data store
        data_store[data_id]["adata"] = adata
        
        # Create result
        result = SpatialDomainResult(
            data_id=data_id,
            method=params.method,
            n_domains=len(domain_counts),
            domain_key=domain_key,
            domain_counts=domain_counts,
            refined_domain_key=refined_domain_key,
            visualization=visualization,
            statistics=statistics,
            embeddings_key=embeddings_key
        )
        
        if context:
            await context.info(f"Successfully identified {len(domain_counts)} spatial domains")
        
        return result
        
    except Exception as e:
        error_msg = f"Error in spatial domain identification: {str(e)}"
        if context:
            await context.warning(error_msg)
        raise RuntimeError(error_msg)


async def _identify_domains_stagate(
    adata: Any,
    params: SpatialDomainParameters,
    context: Optional[Context] = None
) -> tuple:
    """Identify spatial domains using STAGATE"""
    if not STAGATE_AVAILABLE:
        raise ImportError("STAGATE is not installed. Please install it with: pip install STAGATE")
    
    if context:
        await context.info("Running STAGATE for spatial domain identification...")
    
    try:
        # Calculate spatial network
        if context:
            await context.info("Computing spatial neighbor network...")
        stg.Cal_Spatial_Net(adata, rad_cutoff=150)  # Default radius cutoff
        
        # Train STAGATE
        if context:
            await context.info(f"Training STAGATE with {params.stagate_n_epochs} epochs...")
        
        adata = stg.train_STAGATE(
            adata,
            alpha=params.stagate_alpha,
            hidden_dims=params.stagate_hidden_dims,
            n_epochs=params.stagate_n_epochs,
            lr=params.stagate_lr,
            random_seed=params.stagate_random_seed,
            verbose=False
        )
        
        # Perform clustering on STAGATE embeddings
        if context:
            await context.info("Clustering STAGATE embeddings...")
        
        sc.pp.neighbors(adata, use_rep='STAGATE')
        sc.tl.leiden(adata, resolution=params.resolution, key_added='leiden_stagate')
        
        # Use mclust for better clustering if available
        try:
            adata = stg.mclust_R(adata, used_obsm='STAGATE', num_cluster=params.n_domains)
            domain_labels = adata.obs['mclust'].astype(str)
        except:
            # Fallback to leiden clustering
            domain_labels = adata.obs['leiden_stagate'].astype(str)
        
        statistics = {
            "method": "stagate",
            "alpha": params.stagate_alpha,
            "n_epochs": params.stagate_n_epochs,
            "hidden_dims": params.stagate_hidden_dims,
            "embedding_dim": adata.obsm['STAGATE'].shape[1]
        }
        
        return domain_labels, 'STAGATE', statistics
        
    except Exception as e:
        raise RuntimeError(f"STAGATE failed: {str(e)}")


async def _identify_domains_spagcn(
    adata: Any,
    params: SpatialDomainParameters,
    context: Optional[Context] = None
) -> tuple:
    """Identify spatial domains using SpaGCN"""
    if not SPAGCN_AVAILABLE:
        raise ImportError("SpaGCN is not installed. Please install it with: pip install SpaGCN")
    
    if context:
        await context.info("Running SpaGCN for spatial domain identification...")
    
    try:
        # Get spatial coordinates
        if 'spatial' in adata.obsm:
            x_array = adata.obsm['spatial'][:, 0].tolist()
            y_array = adata.obsm['spatial'][:, 1].tolist()
        else:
            raise ValueError("Spatial coordinates not found in adata.obsm['spatial']")
        
        # For SpaGCN, we need pixel coordinates for histology
        # If not available, use array coordinates
        x_pixel = x_array.copy()
        y_pixel = y_array.copy()
        
        # Create a dummy histology image if not available
        img = None
        if params.spagcn_use_histology:
            # Try to get histology image from adata.uns
            if 'spatial' in adata.uns and 'images' in adata.uns['spatial']:
                # This is for 10x Visium data
                img_key = list(adata.uns['spatial'].keys())[0]
                if 'images' in adata.uns['spatial'][img_key]:
                    img_dict = adata.uns['spatial'][img_key]['images']
                    if 'hires' in img_dict:
                        img = img_dict['hires']
                    elif 'lowres' in img_dict:
                        img = img_dict['lowres']
        
        if img is None:
            # Create dummy image or disable histology
            if context:
                await context.info("No histology image found, running SpaGCN without histology")
            params.spagcn_use_histology = False
            img = np.ones((100, 100, 3), dtype=np.uint8) * 255  # White dummy image
        
        # Run SpaGCN easy mode
        if context:
            await context.info("Running SpaGCN domain detection...")
        
        domain_labels = spg.detect_spatial_domains_ez_mode(
            adata,
            img,
            x_array,
            y_array,
            x_pixel,
            y_pixel,
            n_clusters=params.n_domains,
            histology=params.spagcn_use_histology,
            s=params.spagcn_s,
            b=params.spagcn_b,
            p=params.spagcn_p,
            r_seed=params.spagcn_random_seed
        )
        
        domain_labels = pd.Series(domain_labels, index=adata.obs.index).astype(str)
        
        statistics = {
            "method": "spagcn",
            "n_clusters": params.n_domains,
            "s_parameter": params.spagcn_s,
            "b_parameter": params.spagcn_b,
            "p_parameter": params.spagcn_p,
            "use_histology": params.spagcn_use_histology
        }
        
        return domain_labels, None, statistics
        
    except Exception as e:
        raise RuntimeError(f"SpaGCN failed: {str(e)}")


async def _identify_domains_clustering(
    adata: Any,
    params: SpatialDomainParameters,
    context: Optional[Context] = None
) -> tuple:
    """Identify spatial domains using standard clustering methods"""
    if context:
        await context.info(f"Running {params.method} clustering for spatial domain identification...")
    
    try:
        # Compute PCA if not already done
        if 'X_pca' not in adata.obsm:
            if context:
                await context.info("Computing PCA...")
            sc.tl.pca(adata, svd_solver='arpack')
        
        # Compute neighborhood graph
        if context:
            await context.info("Computing neighborhood graph...")
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
        
        # Add spatial information to the neighborhood graph
        if 'spatial' in adata.obsm:
            if context:
                await context.info("Adding spatial constraints to neighborhood graph...")
            sq.gr.spatial_neighbors(adata, coord_type='generic')
            
            # Combine expression and spatial graphs
            spatial_weight = 0.3  # Weight for spatial information
            expr_weight = 1 - spatial_weight
            
            if 'spatial_connectivities' in adata.obsp:
                combined_conn = (expr_weight * adata.obsp['connectivities'] + 
                               spatial_weight * adata.obsp['spatial_connectivities'])
                adata.obsp['connectivities'] = combined_conn
        
        # Perform clustering
        if context:
            await context.info(f"Performing {params.method} clustering...")
        
        if params.method == "leiden":
            sc.tl.leiden(adata, resolution=params.resolution, key_added='spatial_leiden')
            domain_labels = adata.obs['spatial_leiden'].astype(str)
        else:  # louvain
            sc.tl.louvain(adata, resolution=params.resolution, key_added='spatial_louvain')
            domain_labels = adata.obs['spatial_louvain'].astype(str)
        
        statistics = {
            "method": params.method,
            "resolution": params.resolution,
            "n_neighbors": 15,
            "spatial_weight": 0.3 if 'spatial' in adata.obsm else 0.0
        }
        
        return domain_labels, 'X_pca', statistics
        
    except Exception as e:
        raise RuntimeError(f"{params.method} clustering failed: {str(e)}")


def _refine_spatial_domains(adata: Any, domain_key: str, refined_key: str) -> pd.Series:
    """Refine spatial domains using spatial smoothing"""
    try:
        # Get spatial coordinates
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
        else:
            # Use first two principal components as proxy
            coords = adata.obsm['X_pca'][:, :2]
        
        # Get domain labels
        labels = adata.obs[domain_key].astype(str)
        
        # Simple spatial smoothing: assign each spot to the most common domain in its neighborhood
        from sklearn.neighbors import NearestNeighbors
        
        # Find k nearest neighbors
        k = min(10, len(labels) - 1)
        nbrs = NearestNeighbors(n_neighbors=k).fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        
        refined_labels = []
        for i, neighbors in enumerate(indices):
            neighbor_labels = labels.iloc[neighbors]
            # Get most common label among neighbors
            most_common = neighbor_labels.mode()
            if len(most_common) > 0:
                refined_labels.append(most_common.iloc[0])
            else:
                refined_labels.append(labels.iloc[i])
        
        return pd.Series(refined_labels, index=labels.index)
        
    except Exception as e:
        # If refinement fails, return original labels
        return adata.obs[domain_key].astype(str)


def _create_domain_visualization(
    adata: Any,
    domain_key: str,
    params: SpatialDomainParameters
) -> Image:
    """Create visualization of spatial domains"""
    try:
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Get spatial coordinates
        if 'spatial' in adata.obsm:
            x_coords = adata.obsm['spatial'][:, 0]
            y_coords = adata.obsm['spatial'][:, 1]
        else:
            # Fallback to first two PCs
            x_coords = adata.obsm['X_pca'][:, 0]
            y_coords = adata.obsm['X_pca'][:, 1]
        
        # Create scatter plot
        domains = adata.obs[domain_key]
        unique_domains = domains.unique()
        
        # Use a colormap with distinct colors
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_domains)))
        
        for i, domain in enumerate(unique_domains):
            mask = domains == domain
            ax.scatter(
                x_coords[mask], 
                y_coords[mask], 
                c=[colors[i]], 
                label=f'Domain {domain}',
                s=50,
                alpha=0.8
            )
        
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.set_title(f'Spatial Domains ({params.method.upper()})')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Invert y-axis for proper spatial orientation
        ax.invert_yaxis()
        ax.set_aspect('equal')
        
        plt.tight_layout()
        
        # Convert to Image object
        return fig_to_image(fig, dpi=params.image_dpi, format=params.image_format)
        
    except Exception as e:
        # Return placeholder image if visualization fails
        return create_placeholder_image(f"Visualization failed: {str(e)}")
