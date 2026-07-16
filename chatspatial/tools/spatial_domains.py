"""
A module for identifying spatial domains in spatial transcriptomics data.

This module provides an interface to several algorithms designed to partition
spatial data into distinct domains based on gene expression and spatial proximity.
It includes graph-based clustering methods (SpaGCN, STAGATE) and standard clustering
algorithms (Leiden, Louvain) adapted for spatial data. The primary entry point is the `identify_spatial_domains`
function, which handles data preparation and dispatches to the selected method.
"""

import logging
from collections import Counter
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
import scanpy as sc

from ..models.analysis import SpatialDomainResult
from ..models.data import SpatialDomainParameters
from ..utils.adata_utils import (
    ensure_categorical,
    get_spatial_key,
    require_spatial_coords,
    store_analysis_metadata,
)
from ..utils.compute import ensure_neighbors, ensure_pca
from ..utils.dependency_manager import require, require_module
from ..utils.device_utils import resolve_device_async
from ..utils.exceptions import (
    ChatSpatialError,
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)
from ..utils.mcp_utils import suppress_output
from ..utils.results_export import export_analysis_result

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

logger = logging.getLogger(__name__)


def _build_domain_suffix(
    method: str,
    resolution: float,
    n_domains: int,
) -> str:
    """Build a parametric suffix for spatial domain output keys.

    Encodes the parameter that most meaningfully distinguishes runs:
    - leiden/louvain: resolution (e.g. ``leiden_res0_5``)
    - other methods: n_domains (e.g. ``spagcn_n7``)
    """
    if method in ("leiden", "louvain"):
        res_str = f"{resolution:.2f}".replace(".", "_")
        return f"{method}_res{res_str}"
    return f"{method}_n{n_domains}"


async def identify_spatial_domains(
    data_id: str,
    ctx: "ToolContext",
    params: SpatialDomainParameters | None = None,
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
    if params is None:
        params = SpatialDomainParameters()

    # COW FIX: Direct reference instead of copy
    # Only add metadata to adata.obs/obsm/obsp, never overwrite entire adata
    adata = await ctx.get_adata(data_id)

    try:
        # Check if spatial coordinates exist
        spatial_key = get_spatial_key(adata)
        if spatial_key is None:
            raise DataNotFoundError("No spatial coordinates found in the dataset")

        # =================================================================
        # MEMORY OPTIMIZATION: Create working copy exactly once
        # =================================================================
        # Strategy:
        # 1. Determine gene subset (HVG mask) BEFORE copying
        # 2. Check data quality on original (no copy needed for read-only checks)
        # 3. Create single working copy with final gene selection
        # 4. Pass to methods WITHOUT additional copying (they receive independent data)
        # =================================================================

        from scipy.sparse import issparse

        # Step 1: Determine gene subset (no copy yet)
        use_hvg = params.use_highly_variable and "highly_variable" in adata.var.columns
        if params.use_highly_variable and "highly_variable" not in adata.var.columns:
            logger.warning(
                "use_highly_variable=True but 'highly_variable' "
                "column not found in adata.var. "
                "Using all genes instead."
            )
        hvg_mask = adata.var["highly_variable"] if use_hvg else None

        # Step 2: Check data quality on original adata (read-only)
        # Sample a small portion for efficiency
        sample_X = adata.X[:100, :100] if adata.shape[0] > 100 else adata.X
        if issparse(sample_X):
            data_min = sample_X.data.min() if sample_X.data.size > 0 else 0
            data_max = sample_X.data.max() if sample_X.data.size > 0 else 0
        else:
            data_min = float(sample_X.min())
            data_max = float(sample_X.max())

        has_negatives = data_min < 0
        has_large_values = data_max > 100
        use_raw = False

        if has_negatives:
            await ctx.warning(
                f"Data contains negative values (min={data_min:.2f}). "
                "This might indicate scaled/z-scored data. "
                "SpaGCN typically works best with normalized, log-transformed data."
            )
            # Will use raw data if available
            if adata.raw is not None:
                use_raw = True

        elif has_large_values:
            await ctx.warning(
                f"Data contains large values (max={data_max:.2f}). "
                "This might indicate raw count data. "
                "Consider normalizing and log-transforming for better results."
            )

        # Step 3: Create working copy EXACTLY ONCE with final gene selection
        if use_raw:
            # Use raw data (unscaled), subset to HVG if requested
            if hvg_mask is not None:
                # Get genes that are both HVG and in raw
                hvg_genes = adata.var_names[hvg_mask]
                raw_gene_mask = adata.raw.var_names.isin(hvg_genes)
                adata_subset = adata.raw[:, raw_gene_mask].to_adata()
            else:
                adata_subset = adata.raw.to_adata()
            # adata.raw.to_adata() does not carry obsm; restore spatial coords
            # so downstream spatial graph construction works correctly.
            for key in adata.obsm:
                if key not in adata_subset.obsm:
                    adata_subset.obsm[key] = adata.obsm[key]
        elif hvg_mask is not None:
            # Use current X with HVG subset
            adata_subset = adata[:, hvg_mask].copy()
        else:
            # Use full data
            adata_subset = adata.copy()

        # Step 4: In-place data cleaning (no additional copy)
        # Ensure float type for algorithm compatibility
        if adata_subset.X.dtype != np.float32 and adata_subset.X.dtype != np.float64:
            adata_subset.X = adata_subset.X.astype(np.float32)

        # Handle NaN/Inf values in-place
        if issparse(adata_subset.X):
            if adata_subset.X.data.size > 0 and (
                np.any(np.isnan(adata_subset.X.data))
                or np.any(np.isinf(adata_subset.X.data))
            ):
                await ctx.warning(
                    "Found NaN or infinite values in sparse data, replacing with 0"
                )
                adata_subset.X.data = np.nan_to_num(
                    adata_subset.X.data, nan=0.0, posinf=0.0, neginf=0.0
                )
        else:
            if np.any(np.isnan(adata_subset.X)) or np.any(np.isinf(adata_subset.X)):
                await ctx.warning(
                    "Found NaN or infinite values in data, replacing with 0"
                )
                adata_subset.X = np.nan_to_num(
                    adata_subset.X, nan=0.0, posinf=0.0, neginf=0.0
                )

        # NOTE: Removed redundant HVG check (lines 154-159 in original)
        # HVG selection is now handled above in Step 3, avoiding duplicate copy

        # Identify domains based on method
        if params.method == "spagcn":
            domain_labels, embeddings_key, statistics = await _identify_domains_spagcn(
                adata_subset, params, ctx
            )
        elif params.method in ["leiden", "louvain"]:
            (
                domain_labels,
                embeddings_key,
                statistics,
            ) = await _identify_domains_clustering(adata_subset, params, ctx)
        elif params.method == "stagate":
            domain_labels, embeddings_key, statistics = await _identify_domains_stagate(
                adata_subset, params, ctx
            )
        elif params.method == "graphst":
            domain_labels, embeddings_key, statistics = await _identify_domains_graphst(
                adata_subset, params, ctx
            )
        elif params.method == "banksy":
            domain_labels, embeddings_key, statistics = await _identify_domains_banksy(
                adata_subset, params, ctx
            )
        else:
            raise ParameterError(
                f"Unsupported method: {params.method}. Available methods: spagcn, leiden, louvain, stagate, graphst, banksy"
            )

        # Store domain labels in original adata
        # Suffix encodes method + distinguishing param for coexistence
        suffix = _build_domain_suffix(
            params.method, params.resolution, params.n_domains
        )
        domain_key = f"spatial_domains_{suffix}"
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

        # Build results keys for metadata
        results_keys: dict[str, list[str]] = {"obs": [domain_key]}
        if embeddings_key and embeddings_key in adata.obsm:
            results_keys["obsm"] = [embeddings_key]
        if refined_domain_key and refined_domain_key in adata.obs:
            results_keys["obs"].append(refined_domain_key)

        # Store metadata for scientific provenance tracking
        analysis_key = f"spatial_domains_{suffix}"
        store_analysis_metadata(
            adata,
            analysis_name=analysis_key,
            method=params.method,
            parameters={
                "n_domains": params.n_domains,
                "resolution": params.resolution,
                "refine_domains": params.refine_domains,
            },
            results_keys=results_keys,
            statistics=statistics,
        )

        # Export results for reproducibility
        export_analysis_result(adata, data_id, analysis_key)

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

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"Error in spatial domain identification: {e}") from e


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
    # SpaGCN must be imported after restoring SciPy's removed sparse-matrix API.
    from ..utils.compat import ensure_spagcn_compat

    ensure_spagcn_compat()
    spg = require("SpaGCN", ctx, feature="SpaGCN spatial domain identification")
    spagcn_ez = require_module(
        "SpaGCN",
        "SpaGCN.ez_mode",
        ctx,
        feature="SpaGCN spatial domain identification",
    )
    detect_spatial_domains_ez_mode = spagcn_ez.detect_spatial_domains_ez_mode

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

        # Try to get histology image from adata.uns (10x Visium data).
        # Use local flag to avoid mutating the shared params object.
        use_histology = params.spagcn_use_histology

        if use_histology and "spatial" in adata.uns:
            # Select library matching the spatial key, or fall back to
            # the only available library.
            library_ids = list(adata.uns["spatial"].keys())
            lib_id = None
            if len(library_ids) == 1:
                lib_id = library_ids[0]
            elif hasattr(adata, "uns") and "library_id" in adata.obs.columns:
                # Use the most common library in the data
                lib_id = adata.obs["library_id"].mode().iloc[0]
                if lib_id not in library_ids:
                    lib_id = library_ids[0]
            elif library_ids:
                lib_id = library_ids[0]

            # Warn about image-coordinate mismatch in multi-library data
            if len(library_ids) > 1:
                await ctx.warning(
                    f"Multi-library data detected ({len(library_ids)} libraries). "
                    f"SpaGCN will use histology from '{lib_id}' for all spots. "
                    "Image-coordinate mismatch may affect domain boundaries. "
                    "Consider running per-library or set "
                    "spagcn_use_histology=False."
                )

            if lib_id is not None:
                spatial_data = adata.uns["spatial"][lib_id]

                if "images" in spatial_data:
                    img_dict = spatial_data["images"]
                    scalefactors = spatial_data.get("scalefactors", {})

                    if "hires" in img_dict and "tissue_hires_scalef" in scalefactors:
                        img = img_dict["hires"]
                        scale_factor = scalefactors["tissue_hires_scalef"]
                    elif (
                        "lowres" in img_dict and "tissue_lowres_scalef" in scalefactors
                    ):
                        img = img_dict["lowres"]
                        scale_factor = scalefactors["tissue_lowres_scalef"]
                    elif "hires" in img_dict:
                        img = img_dict["hires"]
                    elif "lowres" in img_dict:
                        img = img_dict["lowres"]

        if img is None:
            # No histology available — fall back to expression-only mode
            use_histology = False
            img = np.ones((100, 100, 3), dtype=np.uint8) * 255
        else:
            x_pixel = [int(x * scale_factor) for x in x_array]
            y_pixel = [int(y * scale_factor) for y in y_array]

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

                def _run_spagcn():
                    with suppress_output():
                        return detect_spatial_domains_ez_mode(
                            adata,
                            img,
                            x_array,
                            y_array,
                            x_pixel,
                            y_pixel,
                            n_clusters=params.n_domains,
                            histology=use_histology,
                            s=params.spagcn_s,
                            b=params.spagcn_b,
                            p=params.spagcn_p,
                            r_seed=params.spagcn_random_seed,
                        )

                future = loop.run_in_executor(executor, _run_spagcn)

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
        except ChatSpatialError:
            raise
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
            "use_histology": use_histology,
        }

        return domain_labels, None, statistics

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"SpaGCN execution failed: {e}") from e


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
        n_neighbors = (
            params.cluster_n_neighbors if params.cluster_n_neighbors is not None else 15
        )
        spatial_weight = (
            params.cluster_spatial_weight
            if params.cluster_spatial_weight is not None
            else 0.3
        )

        # Ensure PCA and neighbors are computed (lazy computation)
        ensure_pca(adata)
        ensure_neighbors(adata, n_neighbors=n_neighbors)

        # Add spatial information to the neighborhood graph
        detected_spatial_key = get_spatial_key(adata)
        if detected_spatial_key is not None:
            try:
                sq = require("squidpy", ctx, feature="spatial neighborhood graph")

                # Use squidpy's scientifically validated spatial neighbors
                sq.gr.spatial_neighbors(
                    adata,
                    coord_type="generic",
                    spatial_key=detected_spatial_key,
                )

                # Combine expression and spatial graphs
                expr_weight = 1 - spatial_weight

                if "spatial_connectivities" in adata.obsp:
                    combined_conn = (
                        expr_weight * adata.obsp["connectivities"]
                        + spatial_weight * adata.obsp["spatial_connectivities"]
                    )
                    adata.obsp["connectivities"] = combined_conn

            except ChatSpatialError:
                raise
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
            "spatial_weight": spatial_weight if detected_spatial_key else 0.0,
        }

        return domain_labels, "X_pca", statistics

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"{params.method} clustering failed: {e}") from e


def _refine_spatial_domains(
    adata: Any, domain_key: str, threshold: float = 0.5
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
            # Request k+1 neighbors because kneighbors includes the query point
            # itself as the nearest neighbor (distance=0); we exclude it below.
            nbrs = NearestNeighbors(n_neighbors=k + 1).fit(coords)
            _distances, indices = nbrs.kneighbors(coords)
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
            # Exclude self from neighbor list (first entry is self with dist=0)
            true_neighbors = neighbors[neighbors != i][:k]
            neighbor_labels = labels_values[true_neighbors]

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

    except ChatSpatialError:
        raise
    except Exception as e:
        # Raise error instead of silently failing
        raise ProcessingError(f"Failed to refine spatial domains: {e}") from e


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
    import asyncio

    torch = require("torch", ctx, feature="STAGATE spatial domain identification")

    # Check PyTorch version compatibility with torch_sparse/torch_geometric
    # torch_sparse wheels are only available up to PyTorch 2.8.0
    # See: https://data.pyg.org/whl/
    torch_version = tuple(int(x) for x in torch.__version__.split(".")[:2])
    if torch_version > (2, 8):
        raise ProcessingError(
            f"STAGATE requires PyTorch <= 2.8.0, but found {torch.__version__}. "
            f"torch_sparse/torch_geometric wheels are not available for PyTorch {torch.__version__}. "
            f"Solutions:\n"
            f"  1. Use 'leiden' or 'spagcn' method instead (no PyG dependency)\n"
            f"  2. Downgrade PyTorch: pip install torch==2.8.0\n"
            f"  3. Wait for PyG to support PyTorch {torch.__version__}\n"
            f"See: https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html"
        )

    STAGATE_pyG = require(
        "STAGATE_pyG", ctx, feature="STAGATE spatial domain identification"
    )

    try:
        # MEMORY OPTIMIZATION: adata is already a working copy (adata_subset from caller)
        # No need to copy again - methods receive independent data that can be modified
        adata_stagate = adata

        # Calculate spatial graph
        # STAGATE_pyG uses smaller default radius (50 instead of 150)
        rad_cutoff = (
            params.stagate_rad_cutoff if params.stagate_rad_cutoff is not None else 50
        )
        with suppress_output():
            STAGATE_pyG.Cal_Spatial_Net(adata_stagate, rad_cutoff=rad_cutoff)

        # Optional: Display network statistics
        try:
            with suppress_output():
                STAGATE_pyG.Stats_Spatial_Net(adata_stagate)
        except Exception as exc:
            ctx.debug(f"STAGATE Stats_Spatial_Net skipped: {exc}")

        # Set device (support CUDA, MPS, and CPU)
        device_str = await resolve_device_async(
            prefer_gpu=params.stagate_use_gpu, ctx=ctx
        )
        device = torch.device(device_str)

        # Train STAGATE with timeout protection
        import concurrent.futures

        loop = asyncio.get_running_loop()
        with concurrent.futures.ThreadPoolExecutor() as executor:
            timeout_seconds = params.timeout if params.timeout is not None else 600

            def _train_stagate():
                with suppress_output():
                    return STAGATE_pyG.train_STAGATE(adata_stagate, device=device)

            adata_stagate = await asyncio.wait_for(
                loop.run_in_executor(executor, _train_stagate),
                timeout=timeout_seconds,
            )

        # Get embeddings
        embeddings_key = "STAGATE"
        n_clusters_target = params.n_domains

        # Perform GMM clustering on STAGATE embeddings
        # Using sklearn GaussianMixture with 'tied' covariance (equivalent to mclust EEE)
        # This eliminates R dependency while producing identical results (ARI = 1.0)
        from ..utils.compute import gmm_clustering

        random_seed = (
            params.stagate_random_seed if params.stagate_random_seed is not None else 42
        )
        embedding_data = adata_stagate.obsm[embeddings_key]

        gmm_labels = gmm_clustering(
            data=embedding_data,
            n_clusters=n_clusters_target,
            covariance_type="tied",  # Equivalent to mclust EEE model
            random_state=random_seed,
        )

        # Store in adata - convert to categorical in single operation
        adata_stagate.obs["mclust"] = pd.Categorical(gmm_labels)
        domain_labels = adata_stagate.obs["mclust"].astype(str)
        clustering_method = "gmm_sklearn"  # Updated to reflect actual implementation

        # Copy embeddings to original adata
        adata.obsm[embeddings_key] = adata_stagate.obsm["STAGATE"]

        statistics = {
            "method": "stagate_pyg",
            "n_clusters": len(domain_labels.unique()),
            "target_n_clusters": n_clusters_target,
            "clustering_method": clustering_method,
            "rad_cutoff": rad_cutoff,
            "device": str(device),
            "framework": "PyTorch Geometric",
        }

        return domain_labels, embeddings_key, statistics

    except asyncio.TimeoutError as e:
        raise ProcessingError(
            f"STAGATE training timeout after {params.timeout if params.timeout is not None else 600} seconds"
        ) from e
    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"STAGATE execution failed: {e}") from e


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
    import asyncio
    import concurrent.futures

    torch = require("torch", ctx, feature="GraphST spatial domain identification")

    graphst_module = require_module(
        "GraphST",
        "GraphST.GraphST",
        ctx,
        feature="GraphST spatial domain identification",
    )
    GraphST = graphst_module.GraphST

    try:
        # MEMORY OPTIMIZATION: adata is already a working copy (adata_subset from caller)
        # No need to copy again - methods receive independent data that can be modified
        adata_graphst = adata

        # Set device (support CUDA, MPS, and CPU)
        device_str = await resolve_device_async(
            prefer_gpu=params.graphst_use_gpu, ctx=ctx
        )
        device = torch.device(device_str)

        # Determine number of clusters
        n_clusters = (
            params.graphst_n_clusters
            if params.graphst_n_clusters is not None
            else params.n_domains
        )

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
            timeout_seconds = params.timeout if params.timeout is not None else 600

            def _train_graphst():
                with suppress_output():
                    return model.train()

            adata_graphst = await asyncio.wait_for(
                loop.run_in_executor(executor, _train_graphst),
                timeout=timeout_seconds,
            )

        # Get embeddings key
        embeddings_key = "emb"  # GraphST stores embeddings in adata.obsm['emb']

        refine_label = None
        if params.graphst_refinement:
            graphst_utils = require_module(
                "GraphST",
                "GraphST.utils",
                ctx,
                feature="GraphST spatial domain refinement",
            )
            refine_label = graphst_utils.refine_label

        # Perform clustering on GraphST embeddings
        # OPTIMIZATION: Use binary search instead of GraphST's linear search (290 iterations)
        # GraphST's default search_res uses increment=0.01 from 3.0 to 0.1, which is very slow

        from sklearn.decomposition import PCA

        from ..utils.compute import gmm_clustering

        def run_clustering_optimized():
            with suppress_output():
                # PCA on embeddings (same as GraphST)
                pca = PCA(n_components=20, random_state=42)
                embedding = pca.fit_transform(adata_graphst.obsm["emb"])
                adata_graphst.obsm["emb_pca"] = embedding

                if params.graphst_clustering_method == "mclust":
                    # Use sklearn GMM (equivalent to mclust EEE, eliminates R dependency)
                    gmm_labels = gmm_clustering(
                        data=embedding,
                        n_clusters=n_clusters,
                        covariance_type="tied",  # Equivalent to mclust EEE model
                        random_state=params.graphst_random_seed,
                    )
                    adata_graphst.obs["domain"] = pd.Categorical(gmm_labels)

                    # Apply refinement if requested
                    if params.graphst_refinement:
                        assert refine_label is not None
                        new_type = refine_label(
                            adata_graphst,
                            radius=params.graphst_radius,
                            key="domain",
                        )
                        adata_graphst.obs["domain"] = new_type
                else:
                    # BINARY SEARCH for resolution (replaces GraphST's linear search)
                    # This reduces iterations from 290 to ~10-15
                    sc.pp.neighbors(adata_graphst, n_neighbors=50, use_rep="emb_pca")

                    low, high = 0.1, 3.0
                    best_res, best_diff = 1.0, float("inf")
                    max_iterations = 20  # Binary search converges quickly

                    for _ in range(max_iterations):
                        mid = (low + high) / 2
                        sc.tl.leiden(
                            adata_graphst,
                            resolution=mid,
                            random_state=params.graphst_random_seed,
                        )
                        current_clusters = len(adata_graphst.obs["leiden"].unique())

                        diff = abs(current_clusters - n_clusters)
                        if diff < best_diff:
                            best_diff = diff
                            best_res = mid

                        if current_clusters == n_clusters:
                            break
                        elif current_clusters > n_clusters:
                            high = mid
                        else:
                            low = mid

                        # Early termination if we're close enough
                        if high - low < 0.01:
                            break

                    # Final clustering with best resolution
                    sc.tl.leiden(
                        adata_graphst,
                        resolution=best_res,
                        random_state=params.graphst_random_seed,
                    )
                    adata_graphst.obs["domain"] = adata_graphst.obs["leiden"]

                    # Apply refinement if requested
                    if params.graphst_refinement:
                        assert refine_label is not None
                        new_type = refine_label(
                            adata_graphst,
                            radius=params.graphst_radius,
                            key="domain",
                        )
                        adata_graphst.obs["domain"] = new_type

        # Run clustering in thread pool
        with concurrent.futures.ThreadPoolExecutor() as executor:
            await loop.run_in_executor(executor, run_clustering_optimized)

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
            f"GraphST training timeout after {params.timeout if params.timeout is not None else 600} seconds"
        ) from e
    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"GraphST execution failed: {e}") from e


async def _identify_domains_banksy(
    adata: Any, params: SpatialDomainParameters, ctx: "ToolContext"
) -> tuple:
    """
    Identifies spatial domains using the BANKSY algorithm.

    BANKSY (Building Aggregates with a Neighborhood Kernel and Spatial Yardstick)
    augments gene expression with spatial neighborhood information through:
    1. Neighbor-averaged expression (NBR)
    2. Azimuthal Gabor Filters (AGF) for directional gradients

    Unlike deep learning methods, BANKSY uses explicit mathematical feature
    construction, making it more interpretable and reproducible. This method
    requires the `pybanksy` package.
    """
    import asyncio
    import concurrent.futures

    banksy_embed = require_module(
        "banksy",
        "banksy.embed_banksy",
        ctx,
        feature="BANKSY spatial domain identification",
    )
    banksy_initialize = require_module(
        "banksy",
        "banksy.initialize_banksy",
        ctx,
        feature="BANKSY spatial domain identification",
    )
    generate_banksy_matrix = banksy_embed.generate_banksy_matrix
    initialize_banksy = banksy_initialize.initialize_banksy

    try:
        # MEMORY OPTIMIZATION: adata is already a working copy (adata_subset from caller)
        # No need to copy again - methods receive independent data that can be modified
        adata_banksy = adata

        # Validate and normalize spatial coordinates
        # BANKSY expects coordinates in adata.obsm["spatial"]
        spatial_key = get_spatial_key(adata_banksy)
        if spatial_key is None:
            raise ProcessingError(
                "No spatial coordinates found. Expected in obsm['spatial'], "
                "obsm['X_spatial'], or obsm['coordinates']."
            )

        # Copy coordinates to "spatial" if stored under a different key
        if spatial_key != "spatial":
            adata_banksy.obsm["spatial"] = adata_banksy.obsm[spatial_key]

        # BANKSY coord_keys format: (x_col, y_col, obsm_key)
        # x_col/y_col only used for plotting (disabled), obsm_key is the actual key
        coord_keys = ("x", "y", "spatial")

        # Run BANKSY in thread pool to avoid blocking
        loop = asyncio.get_running_loop()
        timeout_seconds = params.timeout if params.timeout is not None else 600

        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            # Step 1: Initialize BANKSY (compute spatial graphs)
            def init_banksy():
                with suppress_output():
                    return initialize_banksy(
                        adata_banksy,
                        coord_keys=coord_keys,
                        num_neighbours=params.banksy_num_neighbours,
                        max_m=params.banksy_max_m,
                        plt_edge_hist=False,
                        plt_nbr_weights=False,
                        plt_theta=False,
                    )

            banksy_dict = await asyncio.wait_for(
                loop.run_in_executor(executor, init_banksy),
                timeout=timeout_seconds,
            )

            # Step 2: Generate BANKSY matrix (feature augmentation)
            def gen_matrix():
                with suppress_output():
                    return generate_banksy_matrix(
                        adata_banksy,
                        banksy_dict,
                        lambda_list=[params.banksy_lambda],
                        max_m=params.banksy_max_m,
                        verbose=False,
                    )

            _, banksy_matrix = await asyncio.wait_for(
                loop.run_in_executor(executor, gen_matrix),
                timeout=timeout_seconds,
            )

            # Step 3: PCA on BANKSY matrix
            def run_clustering():
                with suppress_output():
                    sc.pp.pca(banksy_matrix, n_comps=params.banksy_pca_dims)
                    sc.pp.neighbors(
                        banksy_matrix,
                        use_rep="X_pca",
                        n_neighbors=params.banksy_num_neighbours,
                    )
                    sc.tl.leiden(
                        banksy_matrix,
                        resolution=params.banksy_cluster_resolution,
                        key_added="banksy_cluster",
                    )
                    return banksy_matrix

            banksy_matrix = await asyncio.wait_for(
                loop.run_in_executor(executor, run_clustering),
                timeout=timeout_seconds,
            )

        # Extract domain labels
        domain_labels = banksy_matrix.obs["banksy_cluster"].astype(str)

        # Store BANKSY embeddings (PCA of augmented matrix)
        embeddings_key = "X_banksy_pca"
        adata.obsm[embeddings_key] = banksy_matrix.obsm["X_pca"]

        # Compute statistics
        n_clusters = len(domain_labels.unique())
        original_features = adata.n_vars
        banksy_features = banksy_matrix.n_vars

        statistics = {
            "method": "banksy",
            "n_clusters": n_clusters,
            "lambda": params.banksy_lambda,
            "num_neighbours": params.banksy_num_neighbours,
            "max_m": params.banksy_max_m,
            "pca_dims": params.banksy_pca_dims,
            "resolution": params.banksy_cluster_resolution,
            "original_features": original_features,
            "banksy_features": banksy_features,
            "feature_expansion": f"{banksy_features / original_features:.1f}x",
        }

        return domain_labels, embeddings_key, statistics

    except asyncio.TimeoutError as e:
        raise ProcessingError(
            f"BANKSY timeout after {params.timeout if params.timeout is not None else 600} seconds"
        ) from e
    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"BANKSY execution failed: {e}") from e
