"""
Spatial Registration Tool

Aligns and registers multiple spatial transcriptomics slices using
optimal transport (PASTE) or diffeomorphic mapping (STalign).
"""

import logging
from typing import TYPE_CHECKING, List, Optional

import numpy as np

if TYPE_CHECKING:
    import anndata as ad
    from ..spatial_mcp_adapter import ToolContext

from ..models.data import RegistrationParameters
from ..utils.adata_utils import ensure_unique_var_names, get_spatial_key
from ..utils.dependency_manager import require
from ..utils.exceptions import ParameterError, ProcessingError

logger = logging.getLogger(__name__)

# =============================================================================
# Validation Helpers
# =============================================================================


def _validate_spatial_coords(adata_list: List["ad.AnnData"]) -> str:
    """
    Validate all slices have spatial coordinates.

    Returns the spatial key found.
    Raises ParameterError if any slice is missing coordinates.
    """
    spatial_key = None
    for i, adata in enumerate(adata_list):
        key = get_spatial_key(adata)
        if key is None:
            raise ParameterError(
                f"Slice {i} missing spatial coordinates. "
                f"Expected in adata.obsm['spatial'] or similar."
            )
        if spatial_key is None:
            spatial_key = key
    return spatial_key or "spatial"


def _get_common_genes(adata_list: List["ad.AnnData"]) -> List[str]:
    """Get common genes across all slices after making names unique."""
    # Make names unique first
    for adata in adata_list:
        ensure_unique_var_names(adata)

    # Then compute intersection
    common = set(adata_list[0].var_names)
    for adata in adata_list[1:]:
        common &= set(adata.var_names)

    genes = list(common)
    logger.info(f"Found {len(genes)} common genes across {len(adata_list)} slices")
    return genes


def _create_ot_backend(use_gpu: bool):
    """Create optimal transport backend for PASTE."""
    import ot

    if use_gpu:
        try:
            import torch

            if torch.cuda.is_available():
                return ot.backend.TorchBackend()
        except ImportError:
            pass
    return ot.backend.NumpyBackend()


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

    # Normalize coordinates to image space with padding
    coords_norm = coords.copy()
    padding = 0.1

    for dim in range(2):
        cmin, cmax = coords[:, dim].min(), coords[:, dim].max()
        crange = cmax - cmin
        if crange > 0:
            target_min = padding * image_size[dim]
            target_max = (1 - padding) * image_size[dim]
            coords_norm[:, dim] = (coords[:, dim] - cmin) / crange
            coords_norm[:, dim] = (
                coords_norm[:, dim] * (target_max - target_min) + target_min
            )

    # Create coordinate grid
    xgrid = [
        torch.linspace(0, image_size[0], image_size[0], dtype=torch.float32),
        torch.linspace(0, image_size[1], image_size[1], dtype=torch.float32),
    ]

    # Rasterize with Gaussian smoothing
    image = np.zeros(image_size, dtype=np.float32)
    for i in range(len(coords)):
        x_idx = int(np.clip(coords_norm[i, 1], 0, image_size[0] - 1))
        y_idx = int(np.clip(coords_norm[i, 0], 0, image_size[1] - 1))

        # Gaussian kernel (radius 2)
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                xi, yi = x_idx + dx, y_idx + dy
                if 0 <= xi < image_size[0] and 0 <= yi < image_size[1]:
                    weight = np.exp(-(dx * dx + dy * dy) / 2.0)
                    image[xi, yi] += intensity[i] * weight

    # Normalize
    if image.max() > 0:
        image /= image.max()

    return xgrid, torch.tensor(image, dtype=torch.float32)


# =============================================================================
# Core Registration Functions
# =============================================================================


def _register_paste(
    adata_list: List["ad.AnnData"],
    params: RegistrationParameters,
    spatial_key: str = "spatial",
) -> List["ad.AnnData"]:
    """Register slices using PASTE optimal transport."""
    import paste as pst
    import scanpy as sc

    reference_idx = params.reference_idx or 0
    registered = [adata.copy() for adata in adata_list]
    common_genes = _get_common_genes(registered)

    if len(registered) == 2:
        # Pairwise alignment
        logger.info("Performing PASTE pairwise alignment")

        slice1 = registered[0][:, common_genes].copy()
        slice2 = registered[1][:, common_genes].copy()

        # Normalize
        sc.pp.normalize_total(slice1, target_sum=1e4)
        sc.pp.log1p(slice1)
        sc.pp.normalize_total(slice2, target_sum=1e4)
        sc.pp.log1p(slice2)

        # Run PASTE
        pi = pst.pairwise_align(
            slice1,
            slice2,
            alpha=params.paste_alpha,
            numItermax=params.paste_numItermax,
            verbose=True,
        )

        # Stack and extract aligned coordinates
        aligned = pst.stack_slices_pairwise([slice1, slice2], [pi])
        registered[0].obsm["spatial_registered"] = aligned[0].obsm["spatial"]
        registered[1].obsm["spatial_registered"] = aligned[1].obsm["spatial"]

    else:
        # Multi-slice center alignment
        logger.info(f"Performing PASTE center alignment with {len(registered)} slices")

        slices = [adata[:, common_genes] for adata in registered]
        backend = _create_ot_backend(params.use_gpu)

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

        # Apply transformations
        for i, (adata, pi) in enumerate(zip(registered, pis_new, strict=False)):
            if i == reference_idx:
                adata.obsm["spatial_registered"] = adata.obsm[spatial_key].copy()
            else:
                adata.obsm["spatial_registered"] = _transform_coordinates(
                    pi, slices[reference_idx].obsm[spatial_key]
                )

    logger.info("PASTE registration completed")
    return registered


def _register_stalign(
    adata_list: List["ad.AnnData"],
    params: RegistrationParameters,
    spatial_key: str = "spatial",
) -> List["ad.AnnData"]:
    """Register slices using STalign diffeomorphic mapping."""
    import STalign.STalign as ST
    import torch

    if len(adata_list) != 2:
        raise ParameterError(
            f"STalign only supports pairwise registration, got {len(adata_list)} slices. "
            f"Use PASTE for multi-slice alignment."
        )

    logger.info("Performing STalign LDDMM registration")

    registered = [adata.copy() for adata in adata_list]
    source, target = registered[0], registered[1]

    # Prepare coordinates
    source_coords = source.obsm[spatial_key].astype(np.float32)
    target_coords = target.obsm[spatial_key].astype(np.float32)

    # Prepare intensity
    if params.stalign_use_expression:
        common_genes = _get_common_genes(registered)
        if len(common_genes) < 100:
            logger.warning(f"Only {len(common_genes)} common genes found")

        # Compute sum intensity (sparse-aware)
        source_expr = source[:, common_genes].X
        target_expr = target[:, common_genes].X

        def _safe_sum(X):
            if hasattr(X, "toarray"):
                return np.array(X.sum(axis=1)).flatten().astype(np.float32)
            return X.sum(axis=1).astype(np.float32)

        source_intensity = _safe_sum(source_expr)
        target_intensity = _safe_sum(target_expr)
    else:
        source_intensity = np.ones(len(source_coords), dtype=np.float32)
        target_intensity = np.ones(len(target_coords), dtype=np.float32)

    logger.info(f"Registering {len(source_coords)} -> {len(target_coords)} spots")

    # Prepare images
    image_size = params.stalign_image_size
    source_grid, source_image = _prepare_stalign_image(
        source_coords, source_intensity, image_size
    )
    target_grid, target_image = _prepare_stalign_image(
        target_coords, target_intensity, image_size
    )

    # STalign parameters
    device = "cuda:0" if params.use_gpu and torch.cuda.is_available() else "cpu"
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

    logger.info(f"Running STalign LDDMM on {image_size} images")

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
            transformed = transformed.numpy()

        source.obsm["spatial_registered"] = transformed
        target.obsm["spatial_registered"] = target_coords.copy()

        logger.info("STalign registration completed")

    except Exception as e:
        raise ProcessingError(
            f"STalign registration failed: {e}. Consider using PASTE method."
        ) from e

    return registered


def _transform_coordinates(
    transport_matrix: np.ndarray,
    reference_coords: np.ndarray,
) -> np.ndarray:
    """Transform coordinates via optimal transport matrix."""
    # Normalize rows
    row_sums = transport_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    normalized = transport_matrix / row_sums

    # Weighted average of reference coordinates
    return normalized @ reference_coords


# =============================================================================
# Public API
# =============================================================================


def register_slices(
    adata_list: List["ad.AnnData"],
    params: Optional[RegistrationParameters] = None,
) -> List["ad.AnnData"]:
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

    logger.info(f"Registering {len(adata_list)} slices using {params.method}")

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
    method: str = "paste",
) -> dict:
    """
    MCP wrapper for spatial registration.

    Args:
        source_id: Source dataset ID
        target_id: Target dataset ID
        ctx: Tool context for data access
        method: Registration method ('paste' or 'stalign')

    Returns:
        Registration result dictionary
    """
    # Check dependencies
    if method == "paste":
        require("paste", ctx, feature="PASTE spatial registration")
    elif method == "stalign":
        require("STalign", ctx, feature="STalign spatial registration")

    # Get data
    source_adata = await ctx.get_adata(source_id)
    target_adata = await ctx.get_adata(target_id)

    # Create parameters
    params = RegistrationParameters(method=method)

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
