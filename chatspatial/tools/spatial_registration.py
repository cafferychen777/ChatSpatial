"""
Spatial Registration Tool

Aligns and registers multiple spatial transcriptomics slices using
optimal transport (PASTE) or diffeomorphic mapping (STalign).
"""

import logging
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    import anndata as ad
    from ..spatial_mcp_adapter import ToolContext

from ..models.data import RegistrationParameters
from ..utils.adata_utils import (
    ensure_unique_var_names,
    find_common_genes,
    get_spatial_key,
    require_spatial_coords,
    shallow_copy_adata,
    store_analysis_metadata,
)
from ..utils.dependency_manager import require
from ..utils.device_utils import get_device, get_ot_backend
from ..utils.exceptions import DataError, ParameterError, ProcessingError
from ..utils.results_export import export_analysis_result

logger = logging.getLogger(__name__)

# =============================================================================
# Validation Helpers
# =============================================================================


def _validate_spatial_coords(adata_list: list["ad.AnnData"]) -> str:
    """Validate all slices have spatial coordinates under the same key.

    Delegates per-slice validation (NaN, Inf, dims, degeneracy) to the
    SSOT ``require_spatial_coords``.  Adds cross-slice key consistency.

    Returns the common spatial key.
    """
    spatial_key: str | None = None
    for i, adata in enumerate(adata_list):
        key = get_spatial_key(adata)
        if key is None:
            raise ParameterError(
                f"Slice {i} missing spatial coordinates. "
                f"Expected in adata.obsm['spatial'] or similar."
            )
        if spatial_key is None:
            spatial_key = key
        elif key != spatial_key:
            raise ParameterError(
                f"Spatial key mismatch: slice 0 uses '{spatial_key}' "
                f"but slice {i} uses '{key}'. All slices must store "
                f"coordinates under the same obsm key."
            )
        # Full validation via SSOT (dims, NaN, Inf, degeneracy, dtype)
        require_spatial_coords(adata, spatial_key=key)
    return spatial_key or "spatial"


def _get_common_genes(adata_list: list["ad.AnnData"]) -> list[str]:
    """Get common genes across all slices after making names unique."""
    # Make names unique first
    for adata in adata_list:
        ensure_unique_var_names(adata)

    # Use unified function for intersection
    genes = find_common_genes(*[adata.var_names for adata in adata_list])
    return genes


# =============================================================================
# STalign Image Preparation (module-level, not nested)
# =============================================================================


def _prepare_stalign_image(
    coords: np.ndarray,
    intensity: np.ndarray,
    image_size: tuple,
) -> tuple:
    """
    Convert point cloud to rasterized image for STalign.

    Args:
        coords: Spatial coordinates (N, 2)
        intensity: Intensity values per point (N,)
        image_size: Output image dimensions (height, width)

    Returns:
        Tuple of (xgrid, image_tensor)
    """
    import torch

    coords = np.asarray(coords)
    intensity = np.asarray(intensity, dtype=np.float32)
    padding = 0.1

    def _to_image_indices(values: np.ndarray, size: int) -> np.ndarray:
        values = np.asarray(values, dtype=np.float32)
        vmin, vmax = values.min(), values.max()
        vrange = vmax - vmin
        if vrange > 0:
            target_min = padding * size
            target_max = (1 - padding) * size
            values = ((values - vmin) / vrange) * (target_max - target_min) + target_min
        return np.clip(values.astype(np.int32, copy=False), 0, size - 1)

    # Create coordinate grid
    xgrid = [
        torch.linspace(0, image_size[0], image_size[0], dtype=torch.float32),
        torch.linspace(0, image_size[1], image_size[1], dtype=torch.float32),
    ]

    # Rasterize with Gaussian smoothing (vectorized for 40-500x speedup)
    from scipy.ndimage import gaussian_filter

    # Vectorized coordinate mapping
    x_indices = _to_image_indices(coords[:, 1], image_size[0])
    y_indices = _to_image_indices(coords[:, 0], image_size[1])

    # Accumulate intensities (np.add.at handles duplicate indices correctly)
    image = np.zeros(image_size, dtype=np.float32)
    np.add.at(image, (x_indices, y_indices), intensity)

    # Apply Gaussian filter (sigma=1.0 approximates the original kernel radius 2)
    image = gaussian_filter(image, sigma=1.0)

    # Normalize
    if image.max() > 0:
        image /= image.max()

    return xgrid, torch.tensor(image, dtype=torch.float32)


# =============================================================================
# Core Registration Functions
# =============================================================================


def _prepare_paste_slices(
    adata_list: list["ad.AnnData"],
    common_genes: list[str],
    spatial_key: str,
) -> list["ad.AnnData"]:
    """Subset to common genes, normalize, and ensure ``obsm["spatial"]`` exists.

    PASTE internally accesses ``adata.obsm["spatial"]``.  When the actual
    spatial key is something else (e.g. ``X_spatial``), we copy it to
    ``"spatial"`` so PASTE finds valid coordinates.

    Normalization (``normalize_total`` + ``log1p``) is applied uniformly so
    that optimal-transport costs are computed on the same expression scale
    regardless of how many slices are being aligned.
    """
    import scanpy as sc

    slices: list["ad.AnnData"] = []
    for adata in adata_list:
        s = adata[:, common_genes].copy()
        # Guarantee PASTE can read obsm["spatial"]
        if spatial_key != "spatial":
            s.obsm["spatial"] = s.obsm[spatial_key]
        sc.pp.normalize_total(s, target_sum=1e4)
        sc.pp.log1p(s)
        slices.append(s)
    return slices


def _register_paste(
    adata_list: list["ad.AnnData"],
    params: RegistrationParameters,
    spatial_key: str = "spatial",
) -> list["ad.AnnData"]:
    """Register slices using PASTE optimal transport."""
    import paste as pst

    reference_idx = params.reference_idx or 0
    if reference_idx < 0 or reference_idx >= len(adata_list):
        raise ParameterError(
            f"reference_idx={reference_idx} is out of range for "
            f"{len(adata_list)} slices "
            f"(valid: 0 to {len(adata_list) - 1})"
        )
    # Memory optimization: keep lightweight working copies and avoid duplicating X.
    registered = [shallow_copy_adata(adata) for adata in adata_list]
    common_genes = _get_common_genes(registered)
    if len(common_genes) == 0:
        raise DataError(
            "No common genes found across slices. "
            "Registration requires shared gene space; "
            "check that input datasets use the same gene naming convention."
        )

    # Unified preparation: gene subset + normalize + ensure obsm["spatial"]
    slices = _prepare_paste_slices(registered, common_genes, spatial_key)

    if len(registered) == 2:
        # Pairwise alignment — reference_idx determines which slice
        # is "fixed" (slice1 in PASTE convention).
        ref_idx = reference_idx  # 0 or 1
        other_idx = 1 - ref_idx

        pi = pst.pairwise_align(
            slices[ref_idx],
            slices[other_idx],
            alpha=params.paste_alpha,
            numItermax=params.paste_numItermax,
            verbose=True,
        )

        # Stack and extract aligned coordinates
        aligned = pst.stack_slices_pairwise([slices[ref_idx], slices[other_idx]], [pi])
        registered[ref_idx].obsm["spatial_registered"] = aligned[0].obsm["spatial"]
        registered[other_idx].obsm["spatial_registered"] = aligned[1].obsm["spatial"]

    else:
        # Multi-slice center alignment
        backend = get_ot_backend(params.use_gpu)

        # Initial pairwise alignments to reference
        pis = []
        for i, slice_data in enumerate(slices):
            if i == reference_idx:
                pis.append(np.eye(slices[i].shape[0]))
            else:
                pi = pst.pairwise_align(
                    slices[reference_idx],
                    slice_data,
                    alpha=params.paste_alpha,
                    backend=backend,
                    use_gpu=params.use_gpu,
                    verbose=False,
                    gpu_verbose=False,
                )
                pis.append(pi)

        # Center alignment
        _, pis_new = pst.center_align(
            slices[reference_idx],
            slices,
            pis_init=pis,
            alpha=params.paste_alpha,
            backend=backend,
            use_gpu=params.use_gpu,
            n_components=params.paste_n_components,
            verbose=False,
            gpu_verbose=False,
        )

        # Apply transformations — read coords from prepared slices
        # (which have "spatial" guaranteed by _prepare_paste_slices)
        for i, (adata, pi) in enumerate(zip(registered, pis_new, strict=False)):
            if i == reference_idx:
                adata.obsm["spatial_registered"] = adata.obsm[spatial_key].copy()
            else:
                adata.obsm["spatial_registered"] = _transform_coordinates(
                    pi,
                    slices[reference_idx].obsm["spatial"],
                    query_coords=adata.obsm[spatial_key],
                )

    return registered


def _register_stalign(
    adata_list: list["ad.AnnData"],
    params: RegistrationParameters,
    spatial_key: str = "spatial",
) -> list["ad.AnnData"]:
    """Register slices using STalign diffeomorphic mapping."""
    import STalign.STalign as ST
    import torch

    if len(adata_list) != 2:
        raise ParameterError(
            f"STalign only supports pairwise registration, got {len(adata_list)} slices. "
            f"Use PASTE for multi-slice alignment."
        )

    # Memory optimization: keep lightweight working copies and avoid duplicating X.
    registered = [shallow_copy_adata(adata) for adata in adata_list]
    source, target = registered[0], registered[1]

    # Prepare coordinates
    source_coords = np.asarray(source.obsm[spatial_key], dtype=np.float32)
    target_coords = np.asarray(target.obsm[spatial_key], dtype=np.float32)

    # Prepare intensity
    if params.stalign_use_expression:
        common_genes = _get_common_genes(registered)
        if len(common_genes) == 0:
            raise DataError(
                "No common genes found between the two slices. "
                "Expression-based STalign registration requires "
                "shared gene space; check gene naming conventions."
            )
        if len(common_genes) < 100:
            logger.warning(f"Only {len(common_genes)} common genes found")

        # Compute sum intensity (sparse-aware)
        source_expr = source[:, common_genes].X
        target_expr = target[:, common_genes].X

        def _safe_sum(X):
            return np.asarray(X.sum(axis=1)).flatten().astype(np.float32)

        source_intensity = _safe_sum(source_expr)
        target_intensity = _safe_sum(target_expr)
    else:
        source_intensity = np.ones(len(source_coords), dtype=np.float32)
        target_intensity = np.ones(len(target_coords), dtype=np.float32)

    # Prepare images
    image_size = tuple(params.stalign_image_size)
    source_grid, source_image = _prepare_stalign_image(
        source_coords, source_intensity, image_size
    )
    target_grid, target_image = _prepare_stalign_image(
        target_coords, target_intensity, image_size
    )

    # STalign parameters
    device = get_device(prefer_gpu=params.use_gpu)
    stalign_params = {
        "a": params.stalign_a,
        "p": 2.0,
        "expand": 2.0,
        "nt": 3,
        "niter": params.stalign_niter,
        "diffeo_start": 0,
        "epL": 2e-08,
        "epT": 0.2,
        "epV": 2000.0,
        "sigmaM": 1.0,
        "sigmaB": 2.0,
        "sigmaA": 5.0,
        "sigmaR": 500000.0,
        "sigmaP": 20.0,
        "device": device,
        "dtype": torch.float32,
    }

    try:
        result = ST.LDDMM(
            xI=source_grid,
            I=source_image,
            xJ=target_grid,
            J=target_image,
            **stalign_params,
        )

        A = result.get("A")
        v = result.get("v")
        xv = result.get("xv")

        if A is None or v is None or xv is None:
            raise ProcessingError("STalign did not return valid transformation")

        # Transform coordinates
        source_points = torch.tensor(source_coords, dtype=torch.float32)
        transformed = ST.transform_points_source_to_target(xv, v, A, source_points)

        if isinstance(transformed, torch.Tensor):
            if hasattr(transformed, "detach"):
                transformed = transformed.detach()
            if hasattr(transformed, "cpu"):
                transformed = transformed.cpu()
            if hasattr(transformed, "numpy"):
                transformed = transformed.numpy()

        source.obsm["spatial_registered"] = np.asarray(transformed, dtype=np.float32)
        target.obsm["spatial_registered"] = np.array(target_coords, copy=True)

    except Exception as e:
        raise ProcessingError(
            f"STalign registration failed: {e}. Consider using PASTE method."
        ) from e

    return registered


def _transform_coordinates(
    transport_matrix: np.ndarray,
    reference_coords: np.ndarray,
    query_coords: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Transform coordinates via optimal transport matrix.

    Cells with no transport signal (zero rows) retain their original
    coordinates when *query_coords* is provided, preventing false
    spatial clustering at the origin.
    """
    # Equivalent to (transport_matrix / row_sums) @ reference_coords, but avoids
    # materializing a large normalized matrix in memory.
    row_sums = transport_matrix.sum(axis=1, keepdims=True)
    zero_mask = row_sums.ravel() == 0
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    weighted = transport_matrix @ reference_coords
    result = weighted / row_sums
    if query_coords is not None and zero_mask.any():
        result[zero_mask] = query_coords[zero_mask]
    return result


# =============================================================================
# Public API
# =============================================================================


def register_slices(
    adata_list: list["ad.AnnData"],
    params: Optional[RegistrationParameters] = None,
) -> list["ad.AnnData"]:
    """
    Register multiple spatial transcriptomics slices.

    Args:
        adata_list: List of AnnData objects to register
        params: Registration parameters (uses defaults if None)

    Returns:
        List of registered AnnData objects with 'spatial_registered' in obsm
    """
    if params is None:
        params = RegistrationParameters()

    if len(adata_list) < 2:
        raise ParameterError("Registration requires at least 2 slices")

    # Validate spatial coordinates and get the spatial key
    spatial_key = _validate_spatial_coords(adata_list)

    if params.method == "paste":
        return _register_paste(adata_list, params, spatial_key)
    elif params.method == "stalign":
        return _register_stalign(adata_list, params, spatial_key)
    else:
        raise ParameterError(f"Unknown method: {params.method}")


# =============================================================================
# MCP Tool Wrapper
# =============================================================================


async def register_spatial_slices_mcp(
    source_id: str,
    target_id: str,
    ctx: "ToolContext",
    params: Optional[RegistrationParameters] = None,
) -> dict:
    """
    MCP wrapper for spatial registration.

    Args:
        source_id: Source dataset ID
        target_id: Target dataset ID
        ctx: Tool context for data access
        params: Registration parameters (method, alignment settings, etc.)

    Returns:
        Registration result dictionary
    """
    if params is None:
        params = RegistrationParameters()

    # Check dependencies
    if params.method == "paste":
        require("paste", ctx, feature="PASTE spatial registration")
    elif params.method == "stalign":
        require("STalign", ctx, feature="STalign spatial registration")

    # Get data
    source_adata = await ctx.get_adata(source_id)
    target_adata = await ctx.get_adata(target_id)

    try:
        registered = register_slices([source_adata, target_adata], params)

        # Copy registered coordinates back (in-place modification)
        if "spatial_registered" in registered[0].obsm:
            source_adata.obsm["spatial_registered"] = registered[0].obsm[
                "spatial_registered"
            ]
        if "spatial_registered" in registered[1].obsm:
            target_adata.obsm["spatial_registered"] = registered[1].obsm[
                "spatial_registered"
            ]

        # Store metadata and export results for both datasets
        method = params.method
        results_keys: dict[str, list[str]] = {"obsm": ["spatial_registered"]}

        # Record all parameters that affect registration results
        # so downstream users can reproduce the exact experiment.
        if method == "paste":
            method_params = {
                "paste_alpha": params.paste_alpha,
                "paste_n_components": params.paste_n_components,
                "paste_numItermax": params.paste_numItermax,
            }
        else:  # stalign
            method_params = {
                "stalign_image_size": list(params.stalign_image_size),
                "stalign_niter": params.stalign_niter,
                "stalign_a": params.stalign_a,
                "stalign_use_expression": params.stalign_use_expression,
            }

        parameters = {
            "method": method,
            "target_id": target_id,
            "use_gpu": params.use_gpu,
            **method_params,
        }
        statistics = {
            "n_source_spots": source_adata.n_obs,
            "n_target_spots": target_adata.n_obs,
        }

        # Export source dataset registration results
        store_analysis_metadata(
            source_adata,
            analysis_name=f"registration_{method}",
            method=method,
            parameters=parameters,
            results_keys=results_keys,
            statistics=statistics,
        )
        export_analysis_result(source_adata, source_id, f"registration_{method}")

        # Export target dataset registration results
        store_analysis_metadata(
            target_adata,
            analysis_name=f"registration_{method}",
            method=method,
            parameters={**parameters, "source_id": source_id},
            results_keys=results_keys,
            statistics=statistics,
        )
        export_analysis_result(target_adata, target_id, f"registration_{method}")

        result = {
            "method": method,
            "source_id": source_id,
            "target_id": target_id,
            "n_source_spots": source_adata.n_obs,
            "n_target_spots": target_adata.n_obs,
            "registration_completed": True,
            "spatial_key_registered": "spatial_registered",
        }

        return result

    except Exception as e:
        raise ProcessingError(f"Registration failed: {e}") from e
