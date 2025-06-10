"""
Spatial domain identification tools for spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import traceback

from mcp.server.fastmcp import Context

from ..models.data import SpatialDomainParameters
from ..models.analysis import SpatialDomainResult


# Import optional dependencies with fallbacks
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

        # Ensure data is float type for SpaGCN compatibility
        if adata_subset.X.dtype != np.float32 and adata_subset.X.dtype != np.float64:
            if context:
                await context.info("Converting data to float32 for SpaGCN compatibility...")
            adata_subset.X = adata_subset.X.astype(np.float32)
        
        # Identify domains based on method
        if params.method == "spagcn":
            domain_labels, embeddings_key, statistics = await _identify_domains_spagcn(
                adata_subset, params, context
            )
        elif params.method in ["leiden", "louvain"]:
            domain_labels, embeddings_key, statistics = await _identify_domains_clustering(
                adata_subset, params, context
            )
        else:
            raise ValueError(f"Unsupported method: {params.method}. Available methods: spagcn, leiden, louvain")
        
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
            try:
                refined_domain_key = f"{domain_key}_refined"
                refined_labels = _refine_spatial_domains(adata, domain_key, refined_domain_key)
                adata.obs[refined_domain_key] = refined_labels
                adata.obs[refined_domain_key] = adata.obs[refined_domain_key].astype('category')
            except Exception as e:
                if context:
                    await context.warning(f"Domain refinement failed: {e}. Proceeding with unrefined domains.")
                refined_domain_key = None  # Reset key if refinement failed
        
        # Get domain counts
        domain_counts = adata.obs[domain_key].value_counts().to_dict()
        domain_counts = {str(k): int(v) for k, v in domain_counts.items()}
        
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

        # Import and call SpaGCN function directly
        from SpaGCN.ez_mode import detect_spatial_domains_ez_mode

        # Add detailed logging for debugging
        if context:
            await context.info(f"SpaGCN parameters: n_clusters={params.n_domains}, s={params.spagcn_s}, b={params.spagcn_b}, p={params.spagcn_p}")
            await context.info(f"Data shape: {adata.shape}, spatial coords: {len(x_array)} spots")

        # Call SpaGCN with error handling
        try:
            domain_labels = detect_spatial_domains_ez_mode(
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
        except Exception as spagcn_error:
            # Capture and re-raise with more details
            error_msg = f"SpaGCN detect_spatial_domains_ez_mode failed: {str(spagcn_error)}"
            if context:
                await context.warning(error_msg)
            raise RuntimeError(error_msg) from spagcn_error

        if context:
            await context.info(f"SpaGCN completed, got {len(set(domain_labels))} domains")

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
        # Enhanced error reporting
        error_msg = f"SpaGCN execution failed: {str(e)}"
        if context:
            await context.warning(f"Full error details: {traceback.format_exc()}")
        raise RuntimeError(error_msg) from e


async def _identify_domains_clustering(
    adata: Any,
    params: SpatialDomainParameters,
    context: Optional[Context] = None
) -> tuple:
    """Identify spatial domains using standard clustering methods"""
    if context:
        await context.info(f"Running {params.method} clustering for spatial domain identification...")
    
    try:
        # Get parameters from params, use defaults if not provided
        n_neighbors = params.cluster_n_neighbors or 15
        spatial_weight = params.cluster_spatial_weight or 0.3
        
        # Compute PCA if not already done
        if 'X_pca' not in adata.obsm:
            if context:
                await context.info("Computing PCA...")
            sc.tl.pca(adata, svd_solver='arpack')
        
        # Compute neighborhood graph
        if context:
            await context.info("Computing neighborhood graph...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca')
        
        # Add spatial information to the neighborhood graph
        if 'spatial' in adata.obsm:
            if context:
                await context.info("Adding spatial constraints to neighborhood graph...")
            sq.gr.spatial_neighbors(adata, coord_type='generic')
            
            # Combine expression and spatial graphs
            expr_weight = 1 - spatial_weight
            
            if 'spatial_connectivities' in adata.obsp:
                combined_conn = (expr_weight * adata.obsp['connectivities'] + 
                               spatial_weight * adata.obsp['spatial_connectivities'])
                adata.obsp['connectivities'] = combined_conn
        
        # Perform clustering
        if context:
            await context.info(f"Performing {params.method} clustering...")
        
        # Use a variable to store key_added to ensure consistency
        key_added = f"spatial_{params.method}"  # e.g., "spatial_leiden" or "spatial_louvain"
        
        if params.method == "leiden":
            sc.tl.leiden(adata, resolution=params.resolution, key_added=key_added)
        else:  # louvain
            sc.tl.louvain(adata, resolution=params.resolution, key_added=key_added)
        
        domain_labels = adata.obs[key_added].astype(str)
        
        statistics = {
            "method": params.method,
            "resolution": params.resolution,
            "n_neighbors": n_neighbors,
            "spatial_weight": spatial_weight if 'spatial' in adata.obsm else 0.0
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
        # Raise error instead of silently failing
        raise RuntimeError(f"Failed to refine spatial domains: {str(e)}") from e
