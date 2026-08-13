"""
A module for identifying spatial domains in spatial transcriptomics data.

This module provides an interface to several algorithms designed to partition
spatial data into distinct domains based on gene expression and spatial proximity.
It includes graph-based clustering methods (SpaGCN, STAGATE) and standard clustering
algorithms (Leiden, Louvain) adapted for spatial data. The primary entry point is the `identify_spatial_domains`
function, which handles data preparation and dispatches to the selected method.
"""

import logging
import sys
import tempfile
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
    to_dense,
)
from ..utils.async_utils import run_sync, run_sync_with_timeout
from ..utils.compute import ensure_neighbors, ensure_pca
from ..utils.dependency_manager import require, require_module
from ..utils.device_utils import resolve_device_async
from ..utils.exceptions import (
    ChatSpatialError,
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)
from ..utils.mcp_utils import suppress_output
from ..utils.results_export import export_analysis_result

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

logger = logging.getLogger(__name__)


def _domain_count_control(params: SpatialDomainParameters) -> str | None:
    """Name the parameter that actually decides how many domains a backend returns.

    Resolution-driven backends cluster by granularity and cannot target an exact
    count, so ``n_domains`` has no effect on them. Naming the real control lets
    callers label outputs honestly and point the user at the knob that works.

    Returns:
        The controlling resolution parameter, or ``None`` when the backend
        honours ``n_domains``.
    """
    if params.method in ("leiden", "louvain"):
        return "resolution"
    if params.method == "banksy":
        return "banksy_cluster_resolution"
    if params.method == "aestetik" and params.aestetik_clustering_method in (
        "leiden",
        "louvain",
    ):
        return "resolution"
    return None


def _build_domain_suffix(params: SpatialDomainParameters) -> str:
    """Build a parametric suffix for spatial domain output keys.

    The suffix names the parameter that distinguishes one run from the next, so
    it must be the parameter the backend actually obeys — encoding ``n_domains``
    for a resolution-driven backend would put a number in the key that the
    result does not honour.
    """
    control = _domain_count_control(params)
    if control is None:
        return f"{params.method}_n{params.n_domains}"
    res_str = f"{getattr(params, control):.2f}".replace(".", "_")
    return f"{params.method}_res{res_str}"


async def _resolve_expression_source(
    adata: Any,
    params: SpatialDomainParameters,
    ctx: "ToolContext",
) -> bool:
    """Decide whether to read ``adata.raw`` and report why.

    Every backend shares this check, so the message names the method actually
    running and states the consequence, rather than guessing which
    normalization produced the values it found.

    Returns:
        True when the backend should read ``adata.raw`` instead of ``adata.X``.
    """
    from scipy.sparse import issparse

    # Sample a small portion for efficiency
    sample_X = adata.X[:100, :100] if adata.shape[0] > 100 else adata.X
    if issparse(sample_X):
        data_min = sample_X.data.min() if sample_X.data.size > 0 else 0
        data_max = sample_X.data.max() if sample_X.data.size > 0 else 0
    else:
        data_min = float(sample_X.min())
        data_max = float(sample_X.max())

    if data_min < 0:
        # Negative values are what variance-stabilizing normalizations produce
        # (Pearson residuals, scaling), so report the fallback, not a guess.
        if adata.raw is not None:
            await ctx.warning(
                f"Expression matrix contains negative values "
                f"(min={data_min:.2f}), which variance-stabilizing "
                f"normalization produces. Using adata.raw for "
                f"{params.method} instead."
            )
            return True
        await ctx.warning(
            f"Expression matrix contains negative values (min={data_min:.2f}) "
            f"and no adata.raw is available, so {params.method} will run on "
            f"them directly."
        )
        return False

    if data_max > 100:
        await ctx.warning(
            f"Expression matrix contains large values (max={data_max:.2f}), "
            f"which suggests unnormalized counts. Normalizing and "
            f"log-transforming before {params.method} usually improves results."
        )

    return False


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

    # Read the managed source without mutating it. Backends receive an isolated
    # workspace below, and the final result is committed only after validation.
    adata = await ctx.get_adata(data_id)

    try:
        # Check if spatial coordinates exist
        spatial_key = get_spatial_key(adata)
        if spatial_key is None:
            raise DataNotFoundError("No spatial coordinates found in the dataset")

        # =================================================================
        # MEMORY OPTIMIZATION: Create the backend working copy exactly once
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
        use_raw = await _resolve_expression_source(adata, params, ctx)

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
        # A resolution-driven backend cannot hit an exact domain count. Say so
        # rather than silently returning a different number than was asked for.
        count_control = _domain_count_control(params)
        if count_control is not None and "n_domains" in params.model_fields_set:
            await ctx.warning(
                f"{params.method} clusters by {count_control}="
                f"{getattr(params, count_control)} and cannot target an exact "
                f"domain count, so n_domains={params.n_domains} is ignored. "
                f"Adjust {count_control} to change how many domains are found."
            )

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
        elif params.method == "aestetik":
            (
                domain_labels,
                embeddings_key,
                statistics,
            ) = await _identify_domains_aestetik(adata_subset, params, ctx)
        else:
            raise ParameterError(
                f"Unsupported method: {params.method}. Available methods: spagcn, leiden, louvain, stagate, graphst, banksy, aestetik"
            )

        # Build the published result on a fresh source copy. Backend workspaces
        # are intentionally isolated above; this second boundary also prevents
        # late metadata, refinement, or export failures from mutating the
        # managed dataset before the complete result is ready.
        adata = adata.copy()

        # Store domain labels in the result candidate
        # Suffix encodes method + distinguishing param for coexistence
        suffix = _build_domain_suffix(params)
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
                if refined_domain_key in adata.obs:
                    del adata.obs[refined_domain_key]
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
        await ctx.set_adata(data_id, adata)
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

            # SpaGCN 1.2.7 reads `adata.X.A`, an attribute scipy removed in
            # 1.14; only its sparse branch does so, and its dense branch is
            # equivalent. Densifying here also spares the two dense copies that
            # branch would otherwise build (one for fit, one for transform).
            adata.X = to_dense(adata.X)

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

            timeout_seconds = params.timeout if params.timeout is not None else 600
            try:
                domain_labels = await run_sync_with_timeout(
                    _run_spagcn,
                    timeout=timeout_seconds,
                    process_name="chatspatial-spagcn",
                )
            except TimeoutError as e:
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

        timeout_seconds = params.timeout if params.timeout is not None else 600

        def _train_stagate():
            with suppress_output():
                return STAGATE_pyG.train_STAGATE(adata_stagate, device=device)

        adata_stagate = await run_sync_with_timeout(
            _train_stagate,
            timeout=timeout_seconds,
            process_name="chatspatial-stagate",
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

    except TimeoutError as e:
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

        timeout_seconds = params.timeout if params.timeout is not None else 600

        def _train_graphst():
            with suppress_output():
                return model.train()

        adata_graphst = await run_sync_with_timeout(
            _train_graphst,
            timeout=timeout_seconds,
            process_name="chatspatial-graphst",
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

        await run_sync(run_clustering_optimized)

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

    except TimeoutError as e:
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

        timeout_seconds = params.timeout if params.timeout is not None else 600

        def _run_banksy():
            with suppress_output():
                banksy_dict = initialize_banksy(
                    adata_banksy,
                    coord_keys=coord_keys,
                    num_neighbours=params.banksy_num_neighbours,
                    max_m=params.banksy_max_m,
                    plt_edge_hist=False,
                    plt_nbr_weights=False,
                    plt_theta=False,
                )
                _, banksy_matrix = generate_banksy_matrix(
                    adata_banksy,
                    banksy_dict,
                    lambda_list=[params.banksy_lambda],
                    max_m=params.banksy_max_m,
                    verbose=False,
                )
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

        banksy_matrix = await run_sync_with_timeout(
            _run_banksy,
            timeout=timeout_seconds,
            process_name="chatspatial-banksy",
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

    except TimeoutError as e:
        raise ProcessingError(
            f"BANKSY timeout after {params.timeout if params.timeout is not None else 600} seconds"
        ) from e
    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"BANKSY execution failed: {e}") from e


# AESTETIK concatenates the two modalities into ``used_obsm_combined``. That
# write happens before its own grid construction reads the transcriptomics
# key, so the combined key must never alias the transcriptomics key.
_AESTETIK_COMBINED_OBSM_KEY = "X_aestetik_combined"
_AESTETIK_TRANSCRIPTOMICS_OBSM_KEY = "X_pca_transcriptomics"
_AESTETIK_MORPHOLOGY_OBSM_KEY = "X_pca_morphology"


def _require_aestetik_representation(
    adata: Any, obsm_key: str, modality: str, guidance: str
) -> np.ndarray:
    """Return a validated 2D representation from ``adata.obsm``."""
    available = sorted(adata.obsm.keys())

    if obsm_key not in adata.obsm:
        raise DataNotFoundError(
            f"AESTETIK requires a precomputed {modality} embedding in "
            f"adata.obsm['{obsm_key}'], which is not present. {guidance} "
            f"Available obsm keys: {available}."
        )

    matrix = np.asarray(adata.obsm[obsm_key])

    if matrix.ndim != 2 or matrix.shape[0] != adata.n_obs or matrix.shape[1] < 1:
        raise DataNotFoundError(
            f"adata.obsm['{obsm_key}'] is not a usable {modality} embedding: "
            f"expected a 2D array with {adata.n_obs} rows and at least one "
            f"column, found shape {matrix.shape}. {guidance}"
        )

    try:
        is_finite = np.all(np.isfinite(matrix))
    except TypeError as exc:
        raise DataNotFoundError(
            f"adata.obsm['{obsm_key}'] is not a numeric {modality} embedding. "
            f"{guidance}"
        ) from exc

    if not is_finite:
        raise DataNotFoundError(
            f"adata.obsm['{obsm_key}'] contains NaN or infinite values and "
            f"cannot be used as the {modality} embedding. {guidance}"
        )

    return matrix


def _managed_aestetik_estimator(
    estimator_cls: type,
    trainer_root: str,
    **kwargs: Any,
) -> Any:
    """Create an AESTETIK estimator whose Lightning trainer writes no artifacts.

    AESTETIK 0.3.x constructs its trainer through a protected factory without
    exposing Lightning configuration. A request-local subclass is safer than a
    process-wide working-directory or monkeypatch guard: concurrent MCP calls
    cannot redirect one another's files or trainer configuration.
    """

    class _ManagedAESTETIK(estimator_cls):
        def _build_trainer(self) -> tuple[Any, list[Any]]:
            trainer, callbacks = super()._build_trainer()
            logger_connector = getattr(trainer, "_logger_connector", None)
            if logger_connector is None or not hasattr(
                logger_connector, "configure_logger"
            ):
                raise DependencyError(
                    "The installed AESTETIK/Lightning combination does not "
                    "support artifact-free server execution. Install the "
                    "supported dependency family with "
                    "pip install 'chatspatial[aestetik]'."
                )
            trainer._default_root_dir = trainer_root
            logger_connector.configure_logger(False)
            return trainer, callbacks

    return _ManagedAESTETIK(**kwargs)


def _install_aestetik_compatibility_inputs(
    adata: Any,
    transcriptomics_key: str,
    morphology_key: str,
) -> dict[str, Any]:
    """Install AESTETIK 0.3.x canonical aliases and return overwritten values.

    AESTETIK exposes configurable input keys, but its training dataset still
    reads cluster labels derived from the two historical key names. The aliases
    live only on ChatSpatial's request-local working copy and are restored after
    fitting, so neither user data nor concurrent requests are affected.
    """
    previous: dict[str, Any] = {}
    aliases = {
        _AESTETIK_TRANSCRIPTOMICS_OBSM_KEY: transcriptomics_key,
        _AESTETIK_MORPHOLOGY_OBSM_KEY: morphology_key,
    }

    for alias, source in aliases.items():
        if alias in adata.obsm and alias != source:
            previous[alias] = adata.obsm[alias]
        elif alias not in adata.obsm:
            previous[alias] = None
        adata.obsm[alias] = adata.obsm[source]

    return previous


def _restore_aestetik_compatibility_inputs(
    adata: Any, previous: dict[str, Any]
) -> None:
    """Restore canonical AESTETIK aliases after request-local fitting."""
    for key, value in previous.items():
        if value is None:
            del adata.obsm[key]
        else:
            adata.obsm[key] = value


def _apply_aestetik_lattice_coordinates(adata: Any) -> str:
    """Populate ``x_array``/``y_array`` on the working copy and report the source.

    AESTETIK matches grid neighbours at an exact integer offset, so pixel
    coordinates are never an acceptable substitute for the array lattice.
    """
    if "x_array" in adata.obs.columns and "y_array" in adata.obs.columns:
        source_columns = ("x_array", "y_array")
    elif "array_row" in adata.obs.columns and "array_col" in adata.obs.columns:
        source_columns = ("array_row", "array_col")
    else:
        raise DataNotFoundError(
            "AESTETIK requires discrete lattice coordinates in adata.obs, as "
            "either 'x_array'/'y_array' or the Visium 'array_row'/'array_col' "
            "columns. Neither pair is present. Pixel coordinates in "
            "obsm['spatial'] cannot be substituted: AESTETIK only accepts grid "
            "neighbours at an exact integer offset, so pixel coordinates would "
            "yield an empty neighborhood for every spot. Reload the dataset "
            "from its Space Ranger output, or add the lattice columns before "
            "running this method. Available obs columns: "
            f"{sorted(adata.obs.columns)}."
        )

    coordinates = []
    for column in source_columns:
        values = pd.to_numeric(adata.obs[column], errors="coerce").to_numpy(dtype=float)
        if not np.all(np.isfinite(values)):
            raise DataNotFoundError(
                f"adata.obs['{column}'] contains non-numeric, NaN, or infinite "
                "values and cannot be used as an AESTETIK lattice coordinate."
            )
        if not np.all(np.equal(np.mod(values, 1), 0)):
            raise DataNotFoundError(
                f"adata.obs['{column}'] holds non-integer values. AESTETIK "
                "matches grid neighbours at an exact integer offset and cannot "
                "use continuous coordinates."
            )
        coordinates.append(values)

    adata.obs["x_array"] = coordinates[0]
    adata.obs["y_array"] = coordinates[1]

    return "/".join(source_columns)


async def _identify_domains_aestetik(
    adata: Any, params: SpatialDomainParameters, ctx: "ToolContext"
) -> tuple:
    """
    Identifies spatial domains using the AESTETIK algorithm.

    AESTETIK is a convolutional autoencoder that learns a representation per
    spot from three inputs: a precomputed transcriptomics embedding, a
    precomputed per-spot morphology embedding, and the spatial neighborhood
    grid around the spot. The learned representations are then clustered to
    define spatial domains. This method requires the `aestetik` package.

    Morphology features are consumed, not produced: ChatSpatial does not
    extract them from tissue images, so ``aestetik_morphology_key`` must
    already be present in the dataset. Cluster refinement is left to the
    shared ``refine_domains`` stage so smoothing is applied exactly once.
    """
    if sys.version_info >= (3, 14):
        raise DependencyError(
            "AESTETIK does not support Python 3.14. Run ChatSpatial with "
            "Python 3.11, 3.12, or 3.13 to use method='aestetik'."
        )

    aestetik = require(
        "aestetik", ctx, feature="AESTETIK spatial domain identification"
    )

    estimator_cls = getattr(aestetik, "AESTETIK", None)
    if estimator_cls is None:
        raise DependencyError(
            "The installed aestetik package does not expose the AESTETIK "
            "estimator. Install a release that provides the fit/predict API: "
            "pip install 'chatspatial[aestetik]'"
        )

    try:
        adata_aestetik = adata

        _require_aestetik_representation(
            adata_aestetik,
            params.aestetik_transcriptomics_key,
            "transcriptomics",
            "Run compute_embeddings first, or point "
            "aestetik_transcriptomics_key at an existing representation.",
        )
        _require_aestetik_representation(
            adata_aestetik,
            params.aestetik_morphology_key,
            "morphology",
            "ChatSpatial does not extract morphology features from tissue "
            "images. Add the embedding to the dataset, or point "
            "aestetik_morphology_key at an existing representation.",
        )

        lattice_source = _apply_aestetik_lattice_coordinates(adata_aestetik)

        # kmeans/bgm take a cluster count; leiden/louvain take a resolution.
        if params.aestetik_clustering_method in ("kmeans", "bgm"):
            n_cluster: float = params.n_domains
        else:
            n_cluster = params.resolution

        timeout_seconds = params.timeout if params.timeout is not None else 600

        def _fit_aestetik():
            previous_inputs = _install_aestetik_compatibility_inputs(
                adata_aestetik,
                params.aestetik_transcriptomics_key,
                params.aestetik_morphology_key,
            )
            try:
                with tempfile.TemporaryDirectory(
                    prefix="chatspatial-aestetik-"
                ) as trainer_root:
                    model = _managed_aestetik_estimator(
                        estimator_cls,
                        trainer_root,
                        n_cluster=n_cluster,
                        morphology_weight=params.aestetik_morphology_weight,
                        window_size=params.aestetik_window_size,
                        clustering_method=params.aestetik_clustering_method,
                        latent_dim=params.aestetik_latent_dim,
                        max_epochs=params.aestetik_max_epochs,
                        random_state=params.aestetik_random_seed,
                        refine_cluster=False,
                        validation_split=0.0,
                        used_obsm_transcriptomics=_AESTETIK_TRANSCRIPTOMICS_OBSM_KEY,
                        used_obsm_morphology=_AESTETIK_MORPHOLOGY_OBSM_KEY,
                        used_obsm_combined=_AESTETIK_COMBINED_OBSM_KEY,
                    )
                    with suppress_output():
                        fitted_model = model.fit(adata_aestetik)
                    return (
                        np.asarray(fitted_model.embedding_),
                        np.asarray(fitted_model.labels_).ravel(),
                    )
            finally:
                _restore_aestetik_compatibility_inputs(adata_aestetik, previous_inputs)

        embedding, labels = await run_sync_with_timeout(
            _fit_aestetik,
            timeout=timeout_seconds,
            process_name="chatspatial-aestetik",
        )

        # Publish the embedding and labels cached by this fit. transform() is
        # stochastic, so recomputing either one here would desynchronize them.
        embedding = np.asarray(embedding)
        labels = np.asarray(labels).ravel()

        if embedding.ndim != 2 or embedding.shape[0] != adata_aestetik.n_obs:
            raise ProcessingError(
                "AESTETIK returned an embedding with unexpected shape "
                f"{embedding.shape} for {adata_aestetik.n_obs} spots."
            )
        try:
            embedding_is_finite = np.all(np.isfinite(embedding))
        except TypeError as exc:
            raise ProcessingError("AESTETIK returned a non-numeric embedding.") from exc
        if not embedding_is_finite:
            raise ProcessingError("AESTETIK returned a non-finite embedding.")
        if labels.shape[0] != adata_aestetik.n_obs:
            raise ProcessingError(
                f"AESTETIK returned {labels.shape[0]} cluster labels for "
                f"{adata_aestetik.n_obs} spots."
            )
        if np.any(pd.isna(labels)):
            raise ProcessingError("AESTETIK returned missing cluster labels.")

        embeddings_key = "X_aestetik"
        adata.obsm[embeddings_key] = embedding
        domain_labels = pd.Series(labels, index=adata_aestetik.obs.index).astype(str)

        statistics = {
            "method": "aestetik",
            "n_clusters": len(domain_labels.unique()),
            "clustering_method": params.aestetik_clustering_method,
            "morphology_weight": params.aestetik_morphology_weight,
            "window_size": params.aestetik_window_size,
            "latent_dim": params.aestetik_latent_dim,
            "max_epochs": params.aestetik_max_epochs,
            "transcriptomics_key": params.aestetik_transcriptomics_key,
            "morphology_key": params.aestetik_morphology_key,
            "lattice_source": lattice_source,
            "framework": "PyTorch Lightning",
        }

        return domain_labels, embeddings_key, statistics

    except TimeoutError as e:
        raise ProcessingError(
            f"AESTETIK timeout after {params.timeout if params.timeout is not None else 600} seconds"
        ) from e
    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"AESTETIK execution failed: {e}") from e
