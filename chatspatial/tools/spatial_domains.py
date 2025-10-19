"""
A module for identifying spatial domains in spatial transcriptomics data.

This module provides an interface to several algorithms designed to partition
spatial data into distinct domains based on gene expression and spatial proximity.
It includes graph-based clustering methods (SpaGCN, STAGATE) and standard clustering
algorithms (Leiden, Louvain) adapted for spatial data. The primary entry point is the `identify_spatial_domains`
function, which handles data preparation and dispatches to the selected method.
"""

import traceback
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
import scanpy as sc

# Import squidpy with fallback
try:
    import squidpy as sq

    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False

from mcp.server.fastmcp import Context

from ..models.analysis import SpatialDomainResult
from ..models.data import SpatialDomainParameters

# Import optional dependencies with fallbacks
try:
    import SpaGCN as spg

    SPAGCN_AVAILABLE = True
except ImportError:
    SPAGCN_AVAILABLE = False

try:
    import STAGATE_pyG

    STAGATE_AVAILABLE = True
except ImportError:
    STAGATE_AVAILABLE = False

try:
    from GraphST.GraphST import GraphST
    from GraphST.utils import clustering as graphst_clustering

    GRAPHST_AVAILABLE = True
except ImportError:
    GRAPHST_AVAILABLE = False

# BANKSY support has been completely removed
# Use alternative methods: spagcn, leiden, louvain, stagate, or graphst


def _check_environment_compatibility():
    """Check environment compatibility for spatial domain identification"""
    issues = []

    # Check squidpy availability
    if not SQUIDPY_AVAILABLE:
        issues.append(
            "squidpy not available - spatial graph functionality will be limited"
        )

    # Check SpaGCN availability
    if not SPAGCN_AVAILABLE:
        issues.append("SpaGCN not available - only clustering methods available")

    # Check STAGATE availability
    if not STAGATE_AVAILABLE:
        issues.append("STAGATE not available - graph attention method unavailable")

    # BANKSY has been removed - use alternative methods

    # Check version compatibility
    try:
        import dask

        dask_version = tuple(map(int, dask.__version__.split(".")[:2]))
        if dask_version >= (2025, 1):
            issues.append(f"dask {dask.__version__} may cause squidpy import issues")
    except ImportError:
        pass

    return issues


async def identify_spatial_domains(
    data_id: str,
    data_store: Dict[str, Any],
    params: SpatialDomainParameters = SpatialDomainParameters(),
    context: Optional[Context] = None,
) -> SpatialDomainResult:
    """
    Identifies spatial domains by clustering spots based on gene expression and location.

    This function serves as the main entry point for various spatial domain
    identification methods. It performs initial data validation and preparation,
    including checks for required preprocessing steps like normalization and
    highly variable gene selection. It then calls the specific algorithm
    requested by the user. The resulting domain labels are stored back in the
    AnnData object.

    Args:
        data_id: The identifier for the dataset.
        data_store: A dictionary that stores the loaded datasets.
        params: An object containing parameters for the analysis, including the
                method to use and its specific settings.
        context: The MCP context for logging and communication.

    Returns:
        A SpatialDomainResult object containing the identified domains and
        associated metadata.
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

    # ✅ COW FIX: Direct reference instead of copy
    # Only add metadata to adata.obs/obsm/obsp, never overwrite entire adata
    adata = data_store[data_id]["adata"]

    try:
        # Check if spatial coordinates exist
        if "spatial" not in adata.obsm and not any(
            "spatial" in key for key in adata.obsm.keys()
        ):
            raise ValueError("No spatial coordinates found in the dataset")

        # Prepare data for domain identification
        if context:
            await context.info("Preparing data for spatial domain identification...")

        # Use highly variable genes if requested and available
        if params.use_highly_variable and "highly_variable" in adata.var.columns:
            adata_subset = adata[:, adata.var["highly_variable"]].copy()
            if context:
                await context.info(
                    f"Using {sum(adata.var['highly_variable'])} highly variable genes"
                )
        else:
            adata_subset = adata.copy()
            if context:
                await context.info(f"Using all {adata.n_vars} genes")

        # Check if data has been scaled (z-score normalized)
        # Scaled data typically has negative values and is centered around 0
        from scipy.sparse import issparse

        # Validate data preprocessing state
        data_min = (
            adata_subset.X.min()
            if not issparse(adata_subset.X)
            else adata_subset.X.data.min()
        )
        data_max = (
            adata_subset.X.max()
            if not issparse(adata_subset.X)
            else adata_subset.X.data.max()
        )

        # Report data statistics for LLM awareness
        has_negatives = data_min < 0
        has_large_values = data_max > 100

        if context:
            await context.info(
                f"Data statistics: min={data_min:.2f}, max={data_max:.2f}, "
                f"has_negatives={has_negatives}, has_large_values={has_large_values}"
            )

        # Provide informative warnings without enforcing
        if has_negatives:
            if context:
                await context.warning(
                    f"Data contains negative values (min={data_min:.2f}). "
                    "This might indicate scaled/z-scored data. "
                    "SpaGCN typically works best with normalized, log-transformed data."
                )

            # Check if raw data is available
            if hasattr(adata, "raw") and adata.raw is not None:
                if context:
                    await context.info(
                        "Raw data is available. Using it may provide better results."
                    )
                gene_mask = adata.raw.var_names.isin(adata_subset.var_names)
                adata_subset = adata.raw[:, gene_mask].to_adata()

        elif has_large_values:
            if context:
                await context.warning(
                    f"Data contains large values (max={data_max:.2f}). "
                    "This might indicate raw count data. "
                    "Consider normalizing and log-transforming for better results."
                )

        # Ensure data is float type for SpaGCN compatibility
        if adata_subset.X.dtype != np.float32 and adata_subset.X.dtype != np.float64:
            if context:
                await context.info(
                    "Converting data to float32 for SpaGCN compatibility..."
                )
            adata_subset.X = adata_subset.X.astype(np.float32)

        # Check for problematic values that can cause SpaGCN to hang
        # Handle both dense and sparse matrices
        from scipy.sparse import issparse

        if issparse(adata_subset.X):
            # For sparse matrices, check the data attribute
            if np.any(np.isnan(adata_subset.X.data)) or np.any(
                np.isinf(adata_subset.X.data)
            ):
                if context:
                    await context.warning(
                        "Found NaN or infinite values in sparse data, replacing with 0"
                    )
                adata_subset.X.data = np.nan_to_num(
                    adata_subset.X.data, nan=0.0, posinf=0.0, neginf=0.0
                )
        else:
            # For dense matrices
            if np.any(np.isnan(adata_subset.X)) or np.any(np.isinf(adata_subset.X)):
                if context:
                    await context.warning(
                        "Found NaN or infinite values in data, replacing with 0"
                    )
                adata_subset.X = np.nan_to_num(
                    adata_subset.X, nan=0.0, posinf=0.0, neginf=0.0
                )

        # Report dataset size for LLM awareness
        if adata_subset.n_obs > 3000:
            if context:
                await context.warning(
                    f"Large dataset with {adata_subset.n_obs} spots. "
                    "Processing may take longer. Consider subsampling if needed."
                )

        # Report gene count
        if adata_subset.n_vars > 3000:
            if context:
                await context.warning(
                    f"Processing {adata_subset.n_vars} genes. "
                    "For faster processing, consider selecting highly variable genes."
                )

        # Use pre-selected highly variable genes if available
        if "highly_variable" in adata_subset.var.columns:
            hvg_count = adata_subset.var["highly_variable"].sum()
            if hvg_count > 0:
                adata_subset = adata_subset[
                    :, adata_subset.var["highly_variable"]
                ].copy()
                if context:
                    await context.info(
                        f"Using {hvg_count} pre-selected highly variable genes"
                    )

        # Apply SpaGCN-specific gene filtering (algorithm requirement)
        if context:
            await context.info("Applying SpaGCN-specific gene filtering...")
        try:
            spg.prefilter_genes(adata_subset, min_cells=3)
            spg.prefilter_specialgenes(adata_subset)
        except Exception as e:
            if context:
                await context.warning(
                    f"SpaGCN gene filtering failed: {e}. Continuing without filtering."
                )

        # Identify domains based on method
        if params.method == "spagcn":
            domain_labels, embeddings_key, statistics = await _identify_domains_spagcn(
                adata_subset, params, context
            )
        elif params.method in ["leiden", "louvain"]:
            domain_labels, embeddings_key, statistics = (
                await _identify_domains_clustering(adata_subset, params, context)
            )
        elif params.method == "stagate":
            domain_labels, embeddings_key, statistics = await _identify_domains_stagate(
                adata_subset, params, context
            )
        elif params.method == "graphst":
            domain_labels, embeddings_key, statistics = await _identify_domains_graphst(
                adata_subset, params, context
            )
        else:
            raise ValueError(
                f"Unsupported method: {params.method}. Available methods: spagcn, leiden, louvain, stagate, graphst"
            )

        # Store domain labels in original adata
        domain_key = f"spatial_domains_{params.method}"
        adata.obs[domain_key] = domain_labels
        adata.obs[domain_key] = adata.obs[domain_key].astype("category")

        # Store embeddings if available
        if embeddings_key and embeddings_key in adata_subset.obsm:
            adata.obsm[embeddings_key] = adata_subset.obsm[embeddings_key]

        # Refine domains if requested
        refined_domain_key = None
        if params.refine_domains:
            if context:
                await context.info(
                    "Refining spatial domains using spatial smoothing..."
                )
            try:
                refined_domain_key = f"{domain_key}_refined"
                refined_labels = _refine_spatial_domains(
                    adata,
                    domain_key,
                    refined_domain_key,
                    threshold=params.refinement_threshold,
                )
                adata.obs[refined_domain_key] = refined_labels
                adata.obs[refined_domain_key] = adata.obs[refined_domain_key].astype(
                    "category"
                )
            except Exception as e:
                if context:
                    await context.warning(
                        f"Domain refinement failed: {e}. Proceeding with unrefined domains."
                    )
                refined_domain_key = None  # Reset key if refinement failed

        # Get domain counts
        domain_counts = adata.obs[domain_key].value_counts().to_dict()
        domain_counts = {str(k): int(v) for k, v in domain_counts.items()}

        # ✅ COW FIX: No need to update data_store - changes already reflected via direct reference
        # All modifications to adata.obs/obsm/obsp are in-place and preserved

        # Create result
        result = SpatialDomainResult(
            data_id=data_id,
            method=params.method,
            n_domains=len(domain_counts),
            domain_key=domain_key,
            domain_counts=domain_counts,
            refined_domain_key=refined_domain_key,
            statistics=statistics,
            embeddings_key=embeddings_key,
        )

        if context:
            await context.info(
                f"Successfully identified {len(domain_counts)} spatial domains"
            )

        return result

    except Exception as e:
        error_msg = f"Error in spatial domain identification: {str(e)}"
        if context:
            await context.warning(error_msg)
        raise RuntimeError(error_msg)


async def _identify_domains_spagcn(
    adata: Any, params: SpatialDomainParameters, context: Optional[Context] = None
) -> tuple:
    """
    Identifies spatial domains using the SpaGCN algorithm.

    SpaGCN (Spatial Graph Convolutional Network) constructs a spatial graph where
    each spot is a node. It then uses a graph convolutional network to learn a
    low-dimensional embedding that integrates gene expression, spatial relationships,
    and optionally histology image features. The final domains are obtained by
    clustering these learned embeddings. This method requires the `SpaGCN` package.
    """
    if not SPAGCN_AVAILABLE:
        raise ImportError(
            "SpaGCN is not installed. Please install it with: pip install SpaGCN"
        )

    if context:
        await context.info("Running SpaGCN for spatial domain identification...")

    try:
        # Get spatial coordinates
        if "spatial" in adata.obsm:
            x_array = adata.obsm["spatial"][:, 0].tolist()
            y_array = adata.obsm["spatial"][:, 1].tolist()
        else:
            raise ValueError("Spatial coordinates not found in adata.obsm['spatial']")

        # Validate spatial coordinates
        if len(x_array) == 0 or len(y_array) == 0:
            raise ValueError("Empty spatial coordinates")

        if np.any(np.isnan(x_array)) or np.any(np.isnan(y_array)):
            raise ValueError("NaN values found in spatial coordinates")

        # Check for degenerate spatial coordinates
        if np.std(x_array) == 0 and np.std(y_array) == 0:
            raise ValueError(
                "All spatial coordinates are identical - cannot identify spatial domains"
            )

        # Adaptive parameter adjustment based on data characteristics
        if context:
            await context.info(
                "Adjusting SpaGCN parameters based on data characteristics..."
            )

        # Adjust parameters based on dataset size and spatial spread
        n_spots = len(x_array)
        spatial_spread = np.std(x_array) + np.std(y_array)

        # Report dataset characteristics for LLM awareness
        if n_spots > 2000:
            if context:
                await context.info(
                    f"Large dataset with {n_spots} spots. "
                    f"Using parameters: s={params.spagcn_s}, p={params.spagcn_p}. "
                    "For faster processing, consider s≤0.5, p≤0.3."
                )
        elif n_spots < 100:
            if context:
                await context.info(
                    f"Small dataset with {n_spots} spots. "
                    f"Using parameters: s={params.spagcn_s}, p={params.spagcn_p}. "
                    "For better sensitivity, consider p≥0.4."
                )

        # Report domain-to-spot ratio for LLM awareness
        spots_per_domain = n_spots / params.n_domains
        if spots_per_domain < 10:
            if context:
                await context.warning(
                    f"Requesting {params.n_domains} domains for {n_spots} spots "
                    f"({spots_per_domain:.1f} spots per domain). "
                    "This may result in unstable or noisy domain assignments."
                )
        elif spots_per_domain < 20:
            if context:
                await context.info(
                    f"Note: {params.n_domains} domains for {n_spots} spots "
                    f"({spots_per_domain:.1f} spots per domain)."
                )

        # For SpaGCN, we need pixel coordinates for histology
        # If not available, use array coordinates
        x_pixel = x_array.copy()
        y_pixel = y_array.copy()

        # Create a dummy histology image if not available
        img = None
        scale_factor = 1.0  # Default scale factor

        if params.spagcn_use_histology:
            # Try to get histology image from adata.uns
            if "spatial" in adata.uns and "images" in adata.uns["spatial"]:
                # This is for 10x Visium data
                img_key = list(adata.uns["spatial"].keys())[0]
                if "images" in adata.uns["spatial"][img_key]:
                    img_dict = adata.uns["spatial"][img_key]["images"]
                    # Check for scalefactors
                    if "scalefactors" in adata.uns["spatial"][img_key]:
                        scalefactors = adata.uns["spatial"][img_key]["scalefactors"]
                        if (
                            "hires" in img_dict
                            and "tissue_hires_scalef" in scalefactors
                        ):
                            img = img_dict["hires"]
                            scale_factor = scalefactors["tissue_hires_scalef"]
                        elif (
                            "lowres" in img_dict
                            and "tissue_lowres_scalef" in scalefactors
                        ):
                            img = img_dict["lowres"]
                            scale_factor = scalefactors["tissue_lowres_scalef"]
                    else:
                        # No scalefactors, try to get image anyway
                        if "hires" in img_dict:
                            img = img_dict["hires"]
                        elif "lowres" in img_dict:
                            img = img_dict["lowres"]

        if img is None:
            # Create dummy image or disable histology
            if context:
                await context.info(
                    "No histology image found, running SpaGCN without histology"
                )
            params.spagcn_use_histology = False
            img = np.ones((100, 100, 3), dtype=np.uint8) * 255  # White dummy image
        else:
            # Apply scale factor to get pixel coordinates
            # SpaGCN expects integer pixel coordinates
            x_pixel = [int(x * scale_factor) for x in x_array]
            y_pixel = [int(y * scale_factor) for y in y_array]
            if context:
                await context.info(f"Using image with scale factor {scale_factor}")

        # Run SpaGCN easy mode
        if context:
            await context.info("Running SpaGCN domain detection...")

        # Import and call SpaGCN function directly
        from SpaGCN.ez_mode import detect_spatial_domains_ez_mode

        # Add detailed logging for debugging
        if context:
            await context.info(
                f"SpaGCN parameters: n_clusters={params.n_domains}, s={params.spagcn_s}, b={params.spagcn_b}, p={params.spagcn_p}"
            )
            await context.info(
                f"Data shape: {adata.shape}, spatial coords: {len(x_array)} spots"
            )

        # Call SpaGCN with error handling and timeout protection
        try:
            # Validate input data before calling SpaGCN
            if len(x_array) != adata.shape[0]:
                raise ValueError(
                    f"Spatial coordinates length ({len(x_array)}) doesn't match data ({adata.shape[0]})"
                )

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
                        adata,  # Pass the adata parameter (which is actually adata_subset)
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
                        r_seed=params.spagcn_random_seed,
                    ),
                )

                # Simple, predictable timeout
                timeout_seconds = (
                    params.timeout if params.timeout else 600
                )  # Default 10 minutes

                if context:
                    await context.info(
                        f"Running SpaGCN with {timeout_seconds}s timeout"
                    )

                try:
                    domain_labels = await asyncio.wait_for(
                        future, timeout=timeout_seconds
                    )
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
            error_msg = (
                f"SpaGCN detect_spatial_domains_ez_mode failed: {str(spagcn_error)}"
            )
            if context:
                await context.warning(error_msg)
            raise RuntimeError(error_msg) from spagcn_error

        if context:
            await context.info(
                f"SpaGCN completed, got {len(set(domain_labels))} domains"
            )

        domain_labels = pd.Series(domain_labels, index=adata.obs.index).astype(str)

        statistics = {
            "method": "spagcn",
            "n_clusters": params.n_domains,
            "s_parameter": params.spagcn_s,
            "b_parameter": params.spagcn_b,
            "p_parameter": params.spagcn_p,
            "use_histology": params.spagcn_use_histology,
        }

        return domain_labels, None, statistics

    except Exception as e:
        # Enhanced error reporting
        error_msg = f"SpaGCN execution failed: {str(e)}"
        if context:
            await context.warning(f"Full error details: {traceback.format_exc()}")
        raise RuntimeError(error_msg) from e


async def _identify_domains_clustering(
    adata: Any, params: SpatialDomainParameters, context: Optional[Context] = None
) -> tuple:
    """
    Identifies spatial domains using Leiden or Louvain clustering on a composite graph.

    This function adapts standard graph-based clustering algorithms for spatial data.
    It first constructs a k-nearest neighbor graph based on gene expression (typically
    from PCA embeddings) and another based on spatial coordinates. These two graphs are
    then combined into a single weighted graph. Applying Leiden or Louvain clustering
    to this composite graph partitions the data into domains that are cohesive in both
    expression and physical space.
    """
    if context:
        await context.info(
            f"Running {params.method} clustering for spatial domain identification..."
        )

    try:
        # Get parameters from params, use defaults if not provided
        n_neighbors = params.cluster_n_neighbors or 15
        spatial_weight = params.cluster_spatial_weight or 0.3

        # Validate PCA requirement (Leiden/Louvain clustering official requirement)
        if "X_pca" not in adata.obsm:
            raise ValueError(
                f"{params.method} clustering requires PCA but X_pca not found. "
                "Please run PCA in preprocessing.py: "
                "sc.tl.pca(adata, n_comps=50)"
            )

        # Validate neighborhood graph requirement (Leiden/Louvain clustering official requirement)
        if "neighbors" not in adata.uns:
            raise ValueError(
                f"{params.method} clustering requires neighborhood graph but neighbors not found. "
                "Please compute neighbors in preprocessing.py: "
                f"sc.pp.neighbors(adata, n_neighbors={n_neighbors}, use_rep='X_pca')"
            )

        if context:
            await context.info(
                f"Using pre-computed PCA and neighborhood graph for {params.method} clustering"
            )

        # Add spatial information to the neighborhood graph
        if "spatial" in adata.obsm:
            if context:
                await context.info(
                    "Adding spatial constraints to neighborhood graph..."
                )

            try:
                if not SQUIDPY_AVAILABLE:
                    # Squidpy is a core dependency for spatial analysis
                    raise ImportError(
                        "❌ CRITICAL: squidpy is required for spatial domain analysis but not available.\n\n"
                        "🔬 SCIENTIFIC INTEGRITY NOTICE:\n"
                        "Spatial domain identification requires proper spatial neighbor graphs.\n"
                        "Alternative methods would compromise scientific validity.\n\n"
                        "💡 SOLUTION:\n"
                        "Install squidpy: pip install 'squidpy>=1.2.0'\n\n"
                        "🚫 Cannot proceed with spatial domain analysis without squidpy."
                    )

                # Use squidpy's scientifically validated spatial neighbors
                sq.gr.spatial_neighbors(adata, coord_type="generic")

                # Combine expression and spatial graphs
                expr_weight = 1 - spatial_weight

                if "spatial_connectivities" in adata.obsp:
                    combined_conn = (
                        expr_weight * adata.obsp["connectivities"]
                        + spatial_weight * adata.obsp["spatial_connectivities"]
                    )
                    adata.obsp["connectivities"] = combined_conn

            except Exception as spatial_error:
                # Spatial graph construction failed - fail honestly instead of degrading to expression-only
                error_msg = (
                    f"Spatial graph construction failed: {spatial_error}. "
                    f"Spatial domain identification requires spatial neighbor relationships. "
                    f"This may indicate issues with spatial coordinates, missing dependencies (squidpy), "
                    f"or insufficient spatial data quality. "
                    f"Please check spatial coordinates in adata.obsm['spatial'] or try a different method."
                )
                if context:
                    await context.error(error_msg)
                raise RuntimeError(error_msg)

        # Perform clustering
        if context:
            await context.info(f"Performing {params.method} clustering...")

        # Use a variable to store key_added to ensure consistency
        key_added = (
            f"spatial_{params.method}"  # e.g., "spatial_leiden" or "spatial_louvain"
        )

        if params.method == "leiden":
            sc.tl.leiden(adata, resolution=params.resolution, key_added=key_added)
        else:  # louvain
            try:
                sc.tl.louvain(adata, resolution=params.resolution, key_added=key_added)
            except ImportError as e:
                # Fallback to leiden if louvain is not available
                if context:
                    await context.warning(
                        f"Louvain not available: {e}. Falling back to Leiden clustering."
                    )
                sc.tl.leiden(adata, resolution=params.resolution, key_added=key_added)

        domain_labels = adata.obs[key_added].astype(str)

        statistics = {
            "method": params.method,
            "resolution": params.resolution,
            "n_neighbors": n_neighbors,
            "spatial_weight": spatial_weight if "spatial" in adata.obsm else 0.0,
        }

        return domain_labels, "X_pca", statistics

    except Exception as e:
        raise RuntimeError(f"{params.method} clustering failed: {str(e)}")


def _refine_spatial_domains(
    adata: Any, domain_key: str, refined_key: str, threshold: float = 0.5
) -> pd.Series:
    """
    Refines spatial domain assignments using a spatial smoothing algorithm.

    This post-processing step aims to create more spatially coherent domains by
    reducing noise. It iterates through each spot and re-assigns its domain label
    to the majority label of its k-nearest spatial neighbors, but ONLY if a
    sufficient proportion of neighbors differ from the current label.

    This threshold-based approach follows SpaGCN (Hu et al., Nature Methods 2021),
    which only relabels spots when more than half of their neighbors are assigned
    to a different domain, preventing over-smoothing while still reducing noise.

    Args:
        adata: AnnData object containing spatial data
        domain_key: Column in adata.obs containing domain labels to refine
        refined_key: Name for the refined domain key
        threshold: Minimum proportion of neighbors that must differ to trigger
                  relabeling (default: 0.5, i.e., 50%, following SpaGCN)

    Returns:
        pd.Series: Refined domain labels
    """
    try:
        # Get spatial coordinates - REQUIRED for spatial domain refinement
        if "spatial" not in adata.obsm:
            raise ValueError(
                "CRITICAL ERROR: Cannot perform spatial domain refinement without physical spatial coordinates!\n\n"
                "Spatial domain refinement requires actual physical positions of spots/cells.\n"
                "The data must contain spatial coordinates in adata.obsm['spatial'].\n\n"
                "What you can do:\n"
                "1. Set refine_domains=False to skip refinement\n"
                "2. Ensure your data has spatial coordinates before refinement\n\n"
                "Note: PCA coordinates represent gene expression space, NOT physical space,\n"
                "and cannot be used for spatial smoothing."
            )

        coords = adata.obsm["spatial"]

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
            original_label = labels.iloc[i]
            neighbor_labels = labels.iloc[neighbors]

            # Calculate proportion of neighbors that differ from current label
            different_ratio = (neighbor_labels != original_label).sum() / len(
                neighbor_labels
            )

            # Only relabel if sufficient proportion of neighbors differ (SpaGCN approach)
            if different_ratio >= threshold:
                # Get most common label among neighbors
                most_common = neighbor_labels.mode()
                if len(most_common) > 0:
                    refined_labels.append(most_common.iloc[0])
                else:
                    refined_labels.append(original_label)
            else:
                # Keep original label if not enough neighbors differ
                refined_labels.append(original_label)

        return pd.Series(refined_labels, index=labels.index)

    except Exception as e:
        # Raise error instead of silently failing
        raise RuntimeError(f"Failed to refine spatial domains: {str(e)}") from e


async def _identify_domains_stagate(
    adata: Any, params: SpatialDomainParameters, context: Optional[Context] = None
) -> tuple:
    """
    Identifies spatial domains using the STAGATE algorithm.

    STAGATE (Spatially-aware graph attention network) learns low-dimensional
    embeddings for spots by integrating gene expression with spatial information
    through a graph attention mechanism. This allows the model to weigh the
    importance of neighboring spots adaptively. The resulting embeddings are then
    clustered to define spatial domains. This method requires the `STAGATE_pyG`
    package.
    """
    if not STAGATE_AVAILABLE:
        raise ImportError(
            "STAGATE_pyG is not installed. Please install it from: https://github.com/QIFEIDKN/STAGATE_pyG"
        )

    if context:
        await context.info("Running STAGATE_pyG for spatial domain identification...")

    try:
        import STAGATE_pyG
        import torch

        # STAGATE_pyG works with preprocessed data
        adata_stagate = adata.copy()

        # Calculate spatial graph
        if context:
            await context.info("Calculating STAGATE spatial network...")

        # STAGATE_pyG uses smaller default radius (50 instead of 150)
        rad_cutoff = params.stagate_rad_cutoff or 50
        STAGATE_pyG.Cal_Spatial_Net(adata_stagate, rad_cutoff=rad_cutoff)

        # Optional: Display network statistics
        if context:
            try:
                STAGATE_pyG.Stats_Spatial_Net(adata_stagate)
            except:
                pass  # Stats display is optional

        # Run STAGATE_pyG
        if context:
            await context.info("Training STAGATE_pyG model...")

        # Set device
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        if context:
            await context.info(f"Using device: {device}")

        # Train STAGATE using simple API
        adata_stagate = STAGATE_pyG.train_STAGATE(adata_stagate, device=device)

        # Get embeddings
        embeddings_key = "STAGATE"

        # Perform clustering on STAGATE embeddings (algorithm requirement)
        if context:
            await context.info(
                "Computing clustering on STAGATE embeddings (algorithm requirement)..."
            )

        # STAGATE-specific neighbors computation (algorithm requirement)
        sc.pp.neighbors(
            adata_stagate,
            use_rep="STAGATE",
            n_neighbors=params.cluster_n_neighbors or 15,
        )

        # Use leiden clustering (mclust is optional and has compatibility issues)
        if context:
            await context.info("Performing Leiden clustering on STAGATE embeddings...")
        sc.tl.leiden(adata_stagate, resolution=params.cluster_resolution or 1.0)
        domain_labels = adata_stagate.obs["leiden"].astype(str)

        # Copy embeddings to original adata
        adata.obsm[embeddings_key] = adata_stagate.obsm["STAGATE"]

        statistics = {
            "method": "stagate_pyg",
            "n_clusters": len(domain_labels.unique()),
            "rad_cutoff": rad_cutoff,
            "device": str(device),
            "framework": "PyTorch Geometric",
        }

        return domain_labels, embeddings_key, statistics

    except Exception as e:
        error_msg = f"STAGATE execution failed: {str(e)}"
        if context:
            await context.warning(error_msg)
        raise RuntimeError(error_msg) from e


async def _identify_domains_graphst(
    adata: Any, params: SpatialDomainParameters, context: Optional[Context] = None
) -> tuple:
    """
    Identifies spatial domains using the GraphST algorithm.

    GraphST (Graph Self-supervised Contrastive Learning) learns spatial domain
    representations by combining graph neural networks with self-supervised
    contrastive learning. It constructs a spatial graph based on spot locations
    and learns embeddings that preserve both gene expression patterns and spatial
    relationships. The learned embeddings are then clustered to define spatial
    domains. This method requires the `GraphST` package.
    """
    if not GRAPHST_AVAILABLE:
        raise ImportError(
            "GraphST is not installed. Please install it with: pip install GraphST"
        )

    if context:
        await context.info("Running GraphST for spatial domain identification...")

    try:
        import asyncio
        import concurrent.futures

        import torch

        # GraphST works with preprocessed data
        adata_graphst = adata.copy()

        # Set device (support CUDA, MPS, and CPU)
        if params.graphst_use_gpu:
            if torch.cuda.is_available():
                device = torch.device("cuda:0")
            elif torch.backends.mps.is_available():
                device = torch.device("mps")
            else:
                device = torch.device("cpu")
                if context:
                    await context.warning(
                        "GPU requested but not available. Using CPU instead."
                    )
        else:
            device = torch.device("cpu")
        if context:
            await context.info(f"Using device: {device}")

        # Initialize GraphST model
        if context:
            await context.info("Initializing GraphST model...")

        # Determine number of clusters
        n_clusters = params.graphst_n_clusters or params.n_domains

        # Initialize model (this is fast, no need for async)
        model = GraphST(
            adata_graphst,
            device=device,
            random_seed=params.graphst_random_seed,
        )

        # Train model (this is blocking, run in executor)
        if context:
            await context.info(
                "Training GraphST model (this may take a few minutes)..."
            )

        # Run training in thread pool to avoid blocking
        loop = asyncio.get_running_loop()
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Set timeout
            timeout_seconds = params.timeout or 600

            adata_graphst = await asyncio.wait_for(
                loop.run_in_executor(executor, lambda: model.train()),
                timeout=timeout_seconds,
            )

        if context:
            await context.info("GraphST training completed successfully")

        # Get embeddings key
        embeddings_key = "emb"  # GraphST stores embeddings in adata.obsm['emb']

        # Perform clustering on GraphST embeddings
        if context:
            await context.info(
                f"Performing {params.graphst_clustering_method} clustering on GraphST embeddings..."
            )

        # Run clustering in thread pool
        with concurrent.futures.ThreadPoolExecutor() as executor:

            def run_clustering():
                graphst_clustering(
                    adata_graphst,
                    n_clusters=n_clusters,
                    radius=params.graphst_radius if params.graphst_refinement else None,
                    method=params.graphst_clustering_method,
                    refinement=params.graphst_refinement,
                )

            await loop.run_in_executor(executor, run_clustering)

        # Get domain labels
        domain_labels = adata_graphst.obs["domain"].astype(str)

        # Copy embeddings to original adata
        adata.obsm[embeddings_key] = adata_graphst.obsm["emb"]

        statistics = {
            "method": "graphst",
            "n_clusters": len(domain_labels.unique()),
            "clustering_method": params.graphst_clustering_method,
            "refinement": params.graphst_refinement,
            "device": str(device),
            "framework": "PyTorch",
        }

        if params.graphst_refinement:
            statistics["refinement_radius"] = params.graphst_radius

        return domain_labels, embeddings_key, statistics

    except asyncio.TimeoutError:
        error_msg = f"GraphST training timeout after {params.timeout or 600} seconds"
        if context:
            await context.warning(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        error_msg = f"GraphST execution failed: {str(e)}"
        if context:
            await context.warning(error_msg)
        raise RuntimeError(error_msg) from e
