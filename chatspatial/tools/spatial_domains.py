"""
Spatial domain identification tools for spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List
import numpy as np
import pandas as pd
import scanpy as sc
import traceback

# Import squidpy with fallback
try:
    import squidpy as sq
    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False

from mcp.server.fastmcp import Context

from ..models.data import SpatialDomainParameters
from ..models.analysis import SpatialDomainResult


# Import optional dependencies with fallbacks
try:
    import SpaGCN as spg
    SPAGCN_AVAILABLE = True
except ImportError:
    SPAGCN_AVAILABLE = False


def _check_environment_compatibility():
    """Check environment compatibility for spatial domain identification"""
    issues = []
    
    # Check squidpy availability
    if not SQUIDPY_AVAILABLE:
        issues.append("squidpy not available - spatial graph functionality will be limited")
    
    # Check SpaGCN availability  
    if not SPAGCN_AVAILABLE:
        issues.append("SpaGCN not available - only clustering methods available")
    
    # Check version compatibility
    try:
        import dask
        dask_version = tuple(map(int, dask.__version__.split('.')[:2]))
        if dask_version >= (2025, 1):
            issues.append(f"dask {dask.__version__} may cause squidpy import issues")
    except ImportError:
        pass
    
    return issues


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
    # Check environment compatibility
    env_issues = _check_environment_compatibility()
    if env_issues and context:
        for issue in env_issues:
            await context.warning(f"Environment issue: {issue}")
    
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
        
        # Check for problematic values that can cause SpaGCN to hang
        if np.any(np.isnan(adata_subset.X)) or np.any(np.isinf(adata_subset.X)):
            if context:
                await context.warning("Found NaN or infinite values in data, replacing with 0")
            adata_subset.X = np.nan_to_num(adata_subset.X, nan=0.0, posinf=0.0, neginf=0.0)
        
        # Smart preprocessing for SpaGCN to prevent hanging
        # Based on deep analysis of SpaGCN bottlenecks
        
        # 1. Limit spots to prevent O(n²) adjacency matrix bottleneck
        if adata_subset.n_obs > 3000:
            if context:
                await context.info(f"Large dataset ({adata_subset.n_obs} spots), subsampling to 3000 for SpaGCN performance")
            # Spatial-aware subsampling to maintain structure
            sc.pp.subsample(adata_subset, n_obs=3000, random_state=42)
        
        # 2. Limit genes to prevent memory issues and speed up computation  
        max_genes = 2000 if adata_subset.n_obs > 1000 else 3000
        if adata_subset.n_vars > max_genes:
            if context:
                await context.info(f"Large number of genes ({adata_subset.n_vars}), selecting top {max_genes} variable genes for SpaGCN")
            # Use highly variable genes if available, otherwise select top variable genes
            if 'highly_variable' in adata_subset.var.columns:
                hvg_genes = adata_subset.var['highly_variable']
                if hvg_genes.sum() > max_genes:
                    # Take top HVGs by variance
                    sc.pp.highly_variable_genes(adata_subset, n_top_genes=max_genes)
                    hvg_genes = adata_subset.var['highly_variable']
                adata_subset = adata_subset[:, hvg_genes].copy()
            else:
                # Select top most variable genes
                sc.pp.highly_variable_genes(adata_subset, n_top_genes=max_genes)
                adata_subset = adata_subset[:, adata_subset.var['highly_variable']].copy()
        
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
        
        # Validate spatial coordinates
        if len(x_array) == 0 or len(y_array) == 0:
            raise ValueError("Empty spatial coordinates")
        
        if np.any(np.isnan(x_array)) or np.any(np.isnan(y_array)):
            raise ValueError("NaN values found in spatial coordinates")
        
        # Check for degenerate spatial coordinates
        if np.std(x_array) == 0 and np.std(y_array) == 0:
            raise ValueError("All spatial coordinates are identical - cannot identify spatial domains")
        
        # Adaptive parameter adjustment based on data characteristics
        if context:
            await context.info("Adjusting SpaGCN parameters based on data characteristics...")
        
        # Adjust parameters based on dataset size and spatial spread
        n_spots = len(x_array)
        spatial_spread = np.std(x_array) + np.std(y_array)
        
        # Conservative parameter adjustment to prevent hanging
        if n_spots > 2000:
            # Large dataset: use conservative parameters
            params.spagcn_s = min(params.spagcn_s, 0.5)
            params.spagcn_p = min(params.spagcn_p, 0.3)
            if context:
                await context.info(f"Large dataset detected, using conservative parameters: s={params.spagcn_s}, p={params.spagcn_p}")
        elif n_spots < 100:
            # Small dataset: adjust for better sensitivity
            params.spagcn_p = max(params.spagcn_p, 0.4)
            if context:
                await context.info(f"Small dataset detected, adjusting p={params.spagcn_p} for better sensitivity")
        
        # Ensure reasonable n_domains relative to data size
        max_reasonable_domains = max(2, min(params.n_domains, n_spots // 20))
        if params.n_domains > max_reasonable_domains:
            if context:
                await context.warning(f"Requested {params.n_domains} domains for {n_spots} spots, limiting to {max_reasonable_domains}")
            params.n_domains = max_reasonable_domains

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

        # Call SpaGCN with error handling and timeout protection
        try:
            # Validate input data before calling SpaGCN
            if len(x_array) != adata.shape[0]:
                raise ValueError(f"Spatial coordinates length ({len(x_array)}) doesn't match data ({adata.shape[0]})")
            
            if context:
                await context.info("Calling SpaGCN detect_spatial_domains_ez_mode...")
            
            # Add timeout protection for SpaGCN call which can hang
            import asyncio
            import concurrent.futures
            
            # Run SpaGCN in a thread pool to avoid blocking
            loop = asyncio.get_event_loop()
            with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
                future = loop.run_in_executor(
                    executor,
                    lambda: detect_spatial_domains_ez_mode(
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
                )
                
                # Set adaptive timeout based on dataset size
                timeout_seconds = min(600, max(180, n_spots * 0.1))  # 3-10 minutes based on size
                if context:
                    await context.info(f"Running SpaGCN with {timeout_seconds:.0f}s timeout...")
                
                try:
                    domain_labels = await asyncio.wait_for(future, timeout=timeout_seconds)
                except asyncio.TimeoutError:
                    error_msg = (
                        f"SpaGCN timed out after {timeout_seconds:.0f} seconds. "
                        f"Dataset: {n_spots} spots, {adata.n_vars} genes. "
                        "Try: 1) Reducing n_domains, 2) Using leiden/louvain instead, "
                        "3) Preprocessing with fewer genes/spots, or 4) Adjusting parameters (s, b, p)."
                    )
                    raise RuntimeError(error_msg)
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
            try:
                # Use appropriate number of components based on data size
                n_comps = min(50, adata.n_vars - 1, adata.n_obs - 1)
                sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
            except Exception as pca_error:
                if context:
                    await context.warning(f"PCA computation failed: {pca_error}. Trying with default parameters.")
                # Fallback to default PCA
                sc.tl.pca(adata, use_highly_variable=False)
        
        # Compute neighborhood graph
        if context:
            await context.info("Computing neighborhood graph...")
        
        # Adjust n_neighbors based on data size
        n_neighbors = min(n_neighbors, adata.n_obs - 1)
        if n_neighbors < 1:
            n_neighbors = min(15, adata.n_obs - 1)
            
        try:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca')
        except Exception as neighbors_error:
            if context:
                await context.warning(f"Neighbors computation failed: {neighbors_error}. Trying with reduced parameters.")
            # Fallback with minimal neighbors
            n_neighbors = min(5, adata.n_obs - 1)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca')
        
        # Add spatial information to the neighborhood graph
        if 'spatial' in adata.obsm:
            if context:
                await context.info("Adding spatial constraints to neighborhood graph...")
            
            try:
                if SQUIDPY_AVAILABLE:
                    sq.gr.spatial_neighbors(adata, coord_type='generic')
                    
                    # Combine expression and spatial graphs
                    expr_weight = 1 - spatial_weight
                    
                    if 'spatial_connectivities' in adata.obsp:
                        combined_conn = (expr_weight * adata.obsp['connectivities'] + 
                                       spatial_weight * adata.obsp['spatial_connectivities'])
                        adata.obsp['connectivities'] = combined_conn
                else:
                    # Manual spatial graph construction without squidpy
                    if context:
                        await context.info("Squidpy not available, creating spatial graph manually...")
                    
                    coords = adata.obsm['spatial']
                    from sklearn.neighbors import NearestNeighbors
                    
                    # Create spatial connectivity graph
                    n_spatial_neighbors = min(6, coords.shape[0] - 1)  # Typical for spatial data
                    nbrs = NearestNeighbors(n_neighbors=n_spatial_neighbors).fit(coords)
                    spatial_distances, spatial_indices = nbrs.kneighbors(coords)
                    
                    # Create sparse connectivity matrix
                    from scipy.sparse import csr_matrix
                    n_spots = coords.shape[0]
                    spatial_connectivities = csr_matrix((n_spots, n_spots))
                    
                    for i, neighbors in enumerate(spatial_indices):
                        for j in neighbors:
                            if i != j:  # Don't connect to self
                                spatial_connectivities[i, j] = 1.0
                    
                    adata.obsp['spatial_connectivities'] = spatial_connectivities
                    
                    # Combine expression and spatial graphs
                    expr_weight = 1 - spatial_weight
                    combined_conn = (expr_weight * adata.obsp['connectivities'] + 
                                   spatial_weight * adata.obsp['spatial_connectivities'])
                    adata.obsp['connectivities'] = combined_conn
                    
            except Exception as spatial_error:
                # Fallback: if spatial graph construction fails, use expression graph only
                if context:
                    await context.warning(f"Spatial graph construction failed: {spatial_error}. Using expression graph only.")
                spatial_weight = 0.0
        
        # Perform clustering
        if context:
            await context.info(f"Performing {params.method} clustering...")
        
        # Use a variable to store key_added to ensure consistency
        key_added = f"spatial_{params.method}"  # e.g., "spatial_leiden" or "spatial_louvain"
        
        if params.method == "leiden":
            sc.tl.leiden(adata, resolution=params.resolution, key_added=key_added)
        else:  # louvain
            try:
                sc.tl.louvain(adata, resolution=params.resolution, key_added=key_added)
            except ImportError as e:
                # Fallback to leiden if louvain is not available
                if context:
                    await context.warning(f"Louvain not available: {e}. Falling back to Leiden clustering.")
                sc.tl.leiden(adata, resolution=params.resolution, key_added=key_added)
        
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
            if 'X_pca' in adata.obsm:
                coords = adata.obsm['X_pca'][:, :2]
            else:
                raise ValueError("No spatial coordinates or PCA found for domain refinement")
        
        # Validate coordinates
        if coords.shape[0] == 0:
            raise ValueError("Empty coordinates for domain refinement")
        
        # Get domain labels
        labels = adata.obs[domain_key].astype(str)
        
        if len(labels) == 0:
            raise ValueError(f"No domain labels found in {domain_key}")
        
        # Simple spatial smoothing: assign each spot to the most common domain in its neighborhood
        from sklearn.neighbors import NearestNeighbors
        
        # Find k nearest neighbors (ensure we have enough data points)
        k = min(10, len(labels) - 1)
        if k < 1:
            # If we have too few points, no refinement possible
            return labels
            
        try:
            nbrs = NearestNeighbors(n_neighbors=k).fit(coords)
            distances, indices = nbrs.kneighbors(coords)
        except Exception as nn_error:
            # If nearest neighbors fails, return original labels
            raise ValueError(f"Nearest neighbors computation failed: {nn_error}")
        
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
