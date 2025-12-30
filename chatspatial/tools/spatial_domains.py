"""
A module for identifying spatial domains in spatial transcriptomics data.

This module provides an interface to several algorithms designed to partition
spatial data into distinct domains based on gene expression and spatial proximity.
It includes graph-based clustering methods (SpaGCN, STAGATE) and standard clustering
algorithms (Leiden, Louvain) adapted for spatial data. The primary entry point is the `identify_spatial_domains`
function, which handles data preparation and dispatches to the selected method.
"""

from collections import Counter
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import scanpy as sc

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import SpatialDomainResult
from ..models.data import SpatialDomainParameters
from ..utils.adata_utils import (
    ensure_categorical,
    get_spatial_key,
    require_spatial_coords,
)
from ..utils.compute import ensure_neighbors, ensure_pca
from ..utils.dependency_manager import require
from ..utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)


async def identify_spatial_domains(
    data_id: str,
    ctx: "ToolContext",
    params: SpatialDomainParameters = SpatialDomainParameters(),
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
        ctx: The unified ToolContext for data access and logging.
        params: An object containing parameters for the analysis, including the
                method to use and its specific settings.

    Returns:
        A SpatialDomainResult object containing the identified domains and
        associated metadata.
    """
    # COW FIX: Direct reference instead of copy
    # Only add metadata to adata.obs/obsm/obsp, never overwrite entire adata
    adata = await ctx.get_adata(data_id)

    try:
        # Check if spatial coordinates exist
        spatial_key = get_spatial_key(adata)
        if spatial_key is None:
            raise DataNotFoundError("No spatial coordinates found in the dataset")

        # Prepare data for domain identification
        # Use highly variable genes if requested and available
        if params.use_highly_variable and "highly_variable" in adata.var.columns:
            adata_subset = adata[:, adata.var["highly_variable"]].copy()
        else:
            adata_subset = adata.copy()

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

        # Check data preprocessing state
        has_negatives = data_min < 0
        has_large_values = data_max > 100

        # Provide informative warnings without enforcing
        if has_negatives:
            await ctx.warning(
                f"Data contains negative values (min={data_min:.2f}). "
                "This might indicate scaled/z-scored data. "
                "SpaGCN typically works best with normalized, log-transformed data."
            )

            # Use raw data if available for better results
            if adata.raw is not None:
                gene_mask = adata.raw.var_names.isin(adata_subset.var_names)
                adata_subset = adata.raw[:, gene_mask].to_adata()

        elif has_large_values:
            await ctx.warning(
                f"Data contains large values (max={data_max:.2f}). "
                "This might indicate raw count data. "
                "Consider normalizing and log-transforming for better results."
            )

        # Ensure data is float type for SpaGCN compatibility
        if adata_subset.X.dtype != np.float32 and adata_subset.X.dtype != np.float64:
            adata_subset.X = adata_subset.X.astype(np.float32)

        # Check for problematic values that can cause SpaGCN to hang
        # Handle both dense and sparse matrices
        from scipy.sparse import issparse

        if issparse(adata_subset.X):
            # For sparse matrices, check the data attribute
            if np.any(np.isnan(adata_subset.X.data)) or np.any(
                np.isinf(adata_subset.X.data)
            ):
                await ctx.warning(
                    "Found NaN or infinite values in sparse data, replacing with 0"
                )
                adata_subset.X.data = np.nan_to_num(
                    adata_subset.X.data, nan=0.0, posinf=0.0, neginf=0.0
                )
        else:
            # For dense matrices
            if np.any(np.isnan(adata_subset.X)) or np.any(np.isinf(adata_subset.X)):
                await ctx.warning(
                    "Found NaN or infinite values in data, replacing with 0"
                )
                adata_subset.X = np.nan_to_num(
                    adata_subset.X, nan=0.0, posinf=0.0, neginf=0.0
                )

        # Use pre-selected highly variable genes if available
        if "highly_variable" in adata_subset.var.columns:
            hvg_count = adata_subset.var["highly_variable"].sum()
            if hvg_count > 0:
                adata_subset = adata_subset[
                    :, adata_subset.var["highly_variable"]
                ].copy()

        # Identify domains based on method
        if params.method == "spagcn":
            domain_labels, embeddings_key, statistics = await _identify_domains_spagcn(
                adata_subset, params, ctx
            )
        elif params.method in ["leiden", "louvain"]:
            domain_labels, embeddings_key, statistics = (
                await _identify_domains_clustering(adata_subset, params, ctx)
            )
        elif params.method == "stagate":
            domain_labels, embeddings_key, statistics = await _identify_domains_stagate(
                adata_subset, params, ctx
            )
        elif params.method == "graphst":
            domain_labels, embeddings_key, statistics = await _identify_domains_graphst(
                adata_subset, params, ctx
            )
        else:
            raise ParameterError(
                f"Unsupported method: {params.method}. Available methods: spagcn, leiden, louvain, stagate, graphst"
            )

        # Store domain labels in original adata
        domain_key = f"spatial_domains_{params.method}"
        adata.obs[domain_key] = domain_labels
        ensure_categorical(adata, domain_key)

        # Store embeddings if available
        if embeddings_key and embeddings_key in adata_subset.obsm:
            adata.obsm[embeddings_key] = adata_subset.obsm[embeddings_key]

        # Refine domains if requested
        refined_domain_key = None
        if params.refine_domains:
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
                await ctx.warning(
                    f"Domain refinement failed: {e}. Proceeding with unrefined domains."
                )
                refined_domain_key = None  # Reset key if refinement failed

        # Get domain counts
        domain_counts = adata.obs[domain_key].value_counts().to_dict()
        domain_counts = {str(k): int(v) for k, v in domain_counts.items()}

        # COW FIX: No need to update data_store - changes already reflected via direct reference
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

        return result

    except Exception as e:
        raise ProcessingError(
            f"Error in spatial domain identification: {str(e)}"
        ) from e


async def _identify_domains_spagcn(
    adata: Any, params: SpatialDomainParameters, ctx: "ToolContext"
) -> tuple:
    """
    Identifies spatial domains using the SpaGCN algorithm.

    SpaGCN (Spatial Graph Convolutional Network) constructs a spatial graph where
    each spot is a node. It then uses a graph convolutional network to learn a
    low-dimensional embedding that integrates gene expression, spatial relationships,
    and optionally histology image features. The final domains are obtained by
    clustering these learned embeddings. This method requires the `SpaGCN` package.
    """
    spg = require("SpaGCN", ctx, feature="SpaGCN spatial domain identification")

    # Apply SpaGCN-specific gene filtering (algorithm requirement)
    try:
        spg.prefilter_genes(adata, min_cells=3)
        spg.prefilter_specialgenes(adata)
    except Exception as e:
        await ctx.warning(
            f"SpaGCN gene filtering failed: {e}. Continuing without filtering."
        )

    try:
        # Get and validate spatial coordinates (auto-detects key, validates NaN/inf/identical)
        coords = require_spatial_coords(adata)
        n_spots = coords.shape[0]

        # Warn about potentially unstable domain assignments
        spots_per_domain = n_spots / params.n_domains
        if spots_per_domain < 10:
            await ctx.warning(
                f"Requesting {params.n_domains} domains for {n_spots} spots "
                f"({spots_per_domain:.1f} spots per domain). "
                "This may result in unstable or noisy domain assignments."
            )

        # For SpaGCN, we need pixel coordinates for histology
        # If not available, use array coordinates
        x_array = coords[:, 0].tolist()
        y_array = coords[:, 1].tolist()
        x_pixel = x_array.copy()
        y_pixel = y_array.copy()

        # Create a dummy histology image if not available
        img = None
        scale_factor = 1.0  # Default scale factor

        if params.spagcn_use_histology:
            # Try to get histology image from adata.uns
            if "spatial" in adata.uns:
                # This is for 10x Visium data
                # Get the first available library ID
                library_ids = list(adata.uns["spatial"].keys())

                if library_ids:
                    lib_id = library_ids[0]
                    spatial_data = adata.uns["spatial"][lib_id]

                    # Try to get image from spatial data
                    if "images" in spatial_data:
                        img_dict = spatial_data["images"]

                        # Try to get scalefactors
                        scalefactors = spatial_data.get("scalefactors", {})

                        # Prefer high-res image, fall back to low-res
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
                        elif "hires" in img_dict:
                            # Try without scalefactor
                            img = img_dict["hires"]
                        elif "lowres" in img_dict:
                            # Try without scalefactor
                            img = img_dict["lowres"]

        if img is None:
            # Create dummy image or disable histology
            params.spagcn_use_histology = False
            img = np.ones((100, 100, 3), dtype=np.uint8) * 255  # White dummy image
        else:
            # Apply scale factor to get pixel coordinates
            x_pixel = [int(x * scale_factor) for x in x_array]
            y_pixel = [int(y * scale_factor) for y in y_array]

        # Import and call SpaGCN function
        from SpaGCN.ez_mode import detect_spatial_domains_ez_mode

        # Call SpaGCN with error handling and timeout protection
        try:
            # Validate input data before calling SpaGCN
            if len(x_array) != adata.shape[0]:
                raise DataError(
                    f"Spatial coordinates length ({len(x_array)}) doesn't match data ({adata.shape[0]})"
                )

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

                try:
                    domain_labels = await asyncio.wait_for(
                        future, timeout=timeout_seconds
                    )
                except asyncio.TimeoutError as e:
                    error_msg = (
                        f"SpaGCN timed out after {timeout_seconds:.0f} seconds. "
                        f"Dataset: {n_spots} spots, {adata.n_vars} genes. "
                        "Try: 1) Reducing n_domains, 2) Using leiden/louvain instead, "
                        "3) Preprocessing with fewer genes/spots, or 4) Adjusting parameters (s, b, p)."
                    )
                    raise ProcessingError(error_msg) from e
        except Exception as spagcn_error:
            raise ProcessingError(
                f"SpaGCN detect_spatial_domains_ez_mode failed: {str(spagcn_error)}"
            ) from spagcn_error

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
        raise ProcessingError(f"SpaGCN execution failed: {str(e)}") from e


async def _identify_domains_clustering(
    adata: Any, params: SpatialDomainParameters, ctx: "ToolContext"
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
    try:
        # Get parameters from params, use defaults if not provided
        n_neighbors = params.cluster_n_neighbors or 15
        spatial_weight = params.cluster_spatial_weight or 0.3

        # Ensure PCA and neighbors are computed (lazy computation)
        ensure_pca(adata)
        ensure_neighbors(adata, n_neighbors=n_neighbors)

        # Add spatial information to the neighborhood graph
        if "spatial" in adata.obsm:

            try:
                sq = require("squidpy", ctx, feature="spatial neighborhood graph")

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
                raise ProcessingError(
                    f"Spatial graph construction failed: {spatial_error}"
                ) from spatial_error

        # Perform clustering
        # Use a variable to store key_added to ensure consistency
        key_added = (
            f"spatial_{params.method}"  # e.g., "spatial_leiden" or "spatial_louvain"
        )

        if params.method == "leiden":
            sc.tl.leiden(adata, resolution=params.resolution, key_added=key_added)
        else:  # louvain
            # Deprecation notice for louvain
            await ctx.warning(
                "Louvain clustering is deprecated and may not be available on all platforms "
                "(especially macOS due to compilation issues). "
                "Consider using 'leiden' instead, which is an improved algorithm with better performance. "
                "Automatic fallback to Leiden will be used if Louvain is unavailable."
            )
            try:
                sc.tl.louvain(adata, resolution=params.resolution, key_added=key_added)
            except ImportError as e:
                # Fallback to leiden if louvain is not available
                await ctx.warning(
                    f"Louvain not available: {e}. Using Leiden clustering instead."
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
        raise ProcessingError(f"{params.method} clustering failed: {str(e)}") from e


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
        # Get and validate spatial coordinates
        coords = require_spatial_coords(adata)

        # Get domain labels
        labels = adata.obs[domain_key].astype(str)

        if len(labels) == 0:
            raise DataNotFoundError("Dataset is empty, cannot refine domains")

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
            # If nearest neighbors fails, raise error
            raise ProcessingError(
                f"Nearest neighbors computation failed: {nn_error}"
            ) from nn_error

        # Optimized: Pre-extract values and use Counter instead of pandas mode()
        # Counter.most_common() is ~6x faster than pandas Series.mode()
        labels_values = labels.values
        refined_labels = []

        for i, neighbors in enumerate(indices):
            original_label = labels_values[i]
            neighbor_labels = labels_values[neighbors]

            # Calculate proportion of neighbors that differ from current label
            different_count = np.sum(neighbor_labels != original_label)
            different_ratio = different_count / len(neighbor_labels)

            # Only relabel if sufficient proportion of neighbors differ (SpaGCN approach)
            if different_ratio >= threshold:
                # Get most common label using Counter (6x faster than pandas mode)
                counter = Counter(neighbor_labels)
                most_common = counter.most_common(1)[0][0]
                refined_labels.append(most_common)
            else:
                # Keep original label if not enough neighbors differ
                refined_labels.append(original_label)

        return pd.Series(refined_labels, index=labels.index)

    except Exception as e:
        # Raise error instead of silently failing
        raise ProcessingError(f"Failed to refine spatial domains: {str(e)}") from e


async def _identify_domains_stagate(
    adata: Any, params: SpatialDomainParameters, ctx: "ToolContext"
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
    STAGATE_pyG = require(
        "STAGATE_pyG", ctx, feature="STAGATE spatial domain identification"
    )
    import torch

    try:
        # STAGATE_pyG works with preprocessed data
        adata_stagate = adata.copy()

        # Calculate spatial graph
        # STAGATE_pyG uses smaller default radius (50 instead of 150)
        rad_cutoff = params.stagate_rad_cutoff or 50
        STAGATE_pyG.Cal_Spatial_Net(adata_stagate, rad_cutoff=rad_cutoff)

        # Optional: Display network statistics
        try:
            STAGATE_pyG.Stats_Spatial_Net(adata_stagate)
        except Exception:
            pass  # Stats display is optional

        # Set device
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

        # Train STAGATE with timeout protection
        import asyncio
        import concurrent.futures

        loop = asyncio.get_running_loop()
        with concurrent.futures.ThreadPoolExecutor() as executor:
            timeout_seconds = params.timeout or 600

            adata_stagate = await asyncio.wait_for(
                loop.run_in_executor(
                    executor,
                    lambda: STAGATE_pyG.train_STAGATE(adata_stagate, device=device),
                ),
                timeout=timeout_seconds,
            )

        # Get embeddings
        embeddings_key = "STAGATE"

        # Perform clustering on STAGATE embeddings (algorithm requirement)
        # STAGATE-specific neighbors computation (algorithm requirement)
        sc.pp.neighbors(
            adata_stagate,
            use_rep="STAGATE",
            n_neighbors=params.cluster_n_neighbors or 15,
        )

        # Use leiden clustering (mclust is optional and has compatibility issues)
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

    except asyncio.TimeoutError as e:
        raise ProcessingError(
            f"STAGATE training timeout after {params.timeout or 600} seconds"
        ) from e
    except Exception as e:
        raise ProcessingError(f"STAGATE execution failed: {str(e)}") from e


async def _identify_domains_graphst(
    adata: Any, params: SpatialDomainParameters, ctx: "ToolContext"
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
    require("GraphST", ctx, feature="GraphST spatial domain identification")
    import asyncio
    import concurrent.futures

    import torch
    from GraphST import GraphST
    from GraphST.utils import clustering as graphst_clustering

    try:
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
                await ctx.warning("GPU requested but not available. Using CPU instead.")
        else:
            device = torch.device("cpu")

        # Determine number of clusters
        n_clusters = params.graphst_n_clusters or params.n_domains

        # Initialize model
        model = GraphST(
            adata_graphst,
            device=device,
            random_seed=params.graphst_random_seed,
        )

        # Train model (this is blocking, run in executor)
        # Run training in thread pool to avoid blocking
        loop = asyncio.get_running_loop()
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Set timeout
            timeout_seconds = params.timeout or 600

            adata_graphst = await asyncio.wait_for(
                loop.run_in_executor(executor, lambda: model.train()),
                timeout=timeout_seconds,
            )

        # Get embeddings key
        embeddings_key = "emb"  # GraphST stores embeddings in adata.obsm['emb']

        # Perform clustering on GraphST embeddings

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

    except asyncio.TimeoutError as e:
        raise ProcessingError(
            f"GraphST training timeout after {params.timeout or 600} seconds"
        ) from e
    except Exception as e:
        raise ProcessingError(f"GraphST execution failed: {str(e)}") from e
