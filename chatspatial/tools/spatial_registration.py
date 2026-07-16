"""
Spatial Registration Tool

Aligns and registers multiple spatial transcriptomics slices using
optimal transport (PASTE) or diffeomorphic mapping (STalign).
"""

import inspect
import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Optional

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
from ..utils.dependency_manager import require, require_module
from ..utils.device_utils import get_device, get_ot_backend
from ..utils.exceptions import (
    ChatSpatialError,
    DataError,
    ParameterError,
    ProcessingError,
)
from ..utils.results_export import export_analysis_result

logger = logging.getLogger(__name__)

# =============================================================================
# Validation Helpers
# =============================================================================


def _patch_paste_line_search_for_pot(
    pst: Any,
    *,
    ot_module: Any | None = None,
) -> None:
    """Patch paste-bio line-search callbacks for POT versions that pass df_G."""
    paste_module = getattr(pst.pairwise_align, "__globals__", {})
    solver = paste_module.get("solve_gromov_linesearch")
    if solver is None:
        return

    solver_params = inspect.signature(solver).parameters
    if "df_G" in solver_params:
        return

    original = paste_module.get("my_fused_gromov_wasserstein")
    if original is None or getattr(original, "_chatspatial_pot_compat", False):
        return
    ot = (
        ot_module
        if ot_module is not None
        else require("ot", feature="PASTE compatibility layer")
    )

    def compatible_my_fused_gromov_wasserstein(
        M,
        C1,
        C2,
        p,
        q,
        G_init=None,
        loss_fun="square_loss",
        alpha=0.5,
        armijo=False,
        log=False,
        numItermax=200,
        tol_rel=1e-9,
        tol_abs=1e-9,
        use_gpu=False,
        **kwargs,
    ):
        p, q = ot.utils.list_to_array(p, q)
        nx = ot.backend.get_backend(p, q, C1, C2, M)
        if G_init is None:
            G0 = p[:, None] * q[None, :]
        else:
            G0 = (1 / nx.sum(G_init)) * G_init
            if use_gpu:
                G0 = G0.cuda()
        constC, hC1, hC2 = ot.gromov.init_matrix(C1, C2, p, q, loss_fun)

        def f(G):
            return ot.gromov.gwloss(constC, hC1, hC2, G)

        def df(G):
            return ot.gromov.gwggrad(constC, hC1, hC2, G)

        armijo_local = True if loss_fun == "kl_loss" else armijo

        if armijo_local:

            def line_search(cost, G, deltaG, Mi, cost_G, df_G=None, **ls_kwargs):
                del df_G
                return ot.optim.line_search_armijo(
                    cost, G, deltaG, Mi, cost_G, nx=nx, **ls_kwargs
                )

        else:

            def line_search(cost, G, deltaG, Mi, cost_G, df_G=None, **ls_kwargs):
                del cost, Mi, df_G
                return solver(
                    G,
                    deltaG,
                    cost_G,
                    C1,
                    C2,
                    M=0.0,
                    reg=1.0,
                    nx=nx,
                    **ls_kwargs,
                )

        if log:
            res, log_dict = ot.optim.cg(
                p,
                q,
                (1 - alpha) * M,
                alpha,
                f,
                df,
                G0,
                line_search,
                log=True,
                numItermax=numItermax,
                stopThr=tol_rel,
                stopThr2=tol_abs,
                **kwargs,
            )
            fgw_dist = log_dict["loss"][-1]
            log_dict["fgw_dist"] = fgw_dist
            return res, log_dict

        return ot.optim.cg(
            p,
            q,
            (1 - alpha) * M,
            alpha,
            f,
            df,
            G0,
            line_search,
            numItermax=numItermax,
            stopThr=tol_rel,
            stopThr2=tol_abs,
            **kwargs,
        )

    compatible_my_fused_gromov_wasserstein.__dict__["_chatspatial_pot_compat"] = True
    paste_module["my_fused_gromov_wasserstein"] = compatible_my_fused_gromov_wasserstein


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


@dataclass(frozen=True)
class _ImageCoordinateTransform:
    """Map spatial x/y coordinates to and from an image row/column grid."""

    x_min: float
    x_max: float
    y_min: float
    y_max: float
    height: int
    width: int
    padding: float = 0.1

    @classmethod
    def from_coords(
        cls,
        coords: np.ndarray,
        image_size: tuple[int, int],
    ) -> "_ImageCoordinateTransform":
        """Build a reversible transform for a non-degenerate 2D point cloud."""
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise DataError(
                f"STalign coordinates must have shape (n_points, 2), got {coords.shape}."
            )
        if coords.shape[0] == 0:
            raise DataError("STalign coordinates must contain at least one point.")
        if not np.isfinite(coords).all():
            raise DataError("STalign coordinates contain non-finite values.")

        height, width = image_size
        if height < 2 or width < 2:
            raise DataError("STalign image dimensions must both be at least 2.")

        x_min = float(coords[:, 0].min())
        x_max = float(coords[:, 0].max())
        y_min = float(coords[:, 1].min())
        y_max = float(coords[:, 1].max())
        if x_max <= x_min or y_max <= y_min:
            raise DataError("STalign coordinates must span both spatial dimensions.")
        return cls(x_min, x_max, y_min, y_max, height, width)

    @staticmethod
    def _scale(
        values: np.ndarray,
        source_min: float,
        source_max: float,
        target_size: int,
        padding: float,
    ) -> np.ndarray:
        target_max = float(target_size - 1)
        lower = padding * target_max
        upper = (1.0 - padding) * target_max
        return lower + (values - source_min) * (upper - lower) / (
            source_max - source_min
        )

    @staticmethod
    def _unscale(
        values: np.ndarray,
        target_min: float,
        target_max: float,
        source_size: int,
        padding: float,
    ) -> np.ndarray:
        source_max_value = float(source_size - 1)
        lower = padding * source_max_value
        upper = (1.0 - padding) * source_max_value
        return target_min + (values - lower) * (target_max - target_min) / (
            upper - lower
        )

    def to_image(self, coords: np.ndarray) -> np.ndarray:
        """Convert spatial x/y points to STalign row/column coordinates."""
        coords = np.asarray(coords, dtype=np.float32)
        rows = self._scale(
            coords[:, 1], self.y_min, self.y_max, self.height, self.padding
        )
        columns = self._scale(
            coords[:, 0], self.x_min, self.x_max, self.width, self.padding
        )
        return np.column_stack((rows, columns)).astype(np.float32, copy=False)

    def from_image(self, row_column: np.ndarray) -> np.ndarray:
        """Convert STalign row/column points back to spatial x/y coordinates."""
        row_column = np.asarray(row_column, dtype=np.float32)
        x = self._unscale(
            row_column[:, 1], self.x_min, self.x_max, self.width, self.padding
        )
        y = self._unscale(
            row_column[:, 0], self.y_min, self.y_max, self.height, self.padding
        )
        return np.column_stack((x, y)).astype(np.float32, copy=False)


def _prepare_stalign_image(
    coords: np.ndarray,
    intensity: np.ndarray,
    image_size: tuple[int, int],
    *,
    torch_module: Any | None = None,
) -> tuple[list[Any], Any, _ImageCoordinateTransform]:
    """
    Convert point cloud to rasterized image for STalign.

    Args:
        coords: Spatial coordinates (N, 2)
        intensity: Intensity values per point (N,)
        image_size: Output image dimensions (height, width)

    Returns:
        Tuple of (row/column grid, channel-first image tensor, coordinate transform)
    """
    coords = np.asarray(coords, dtype=np.float32)
    intensity = np.asarray(intensity, dtype=np.float32)
    transform = _ImageCoordinateTransform.from_coords(coords, image_size)
    if intensity.ndim != 1 or intensity.shape[0] != coords.shape[0]:
        raise DataError(
            "STalign intensity must be a one-dimensional value per spatial point."
        )
    if not np.isfinite(intensity).all():
        raise DataError("STalign intensity contains non-finite values.")
    if np.any(intensity < 0):
        raise DataError("STalign intensity must be non-negative.")
    if not np.any(intensity > 0):
        raise DataError("STalign intensity must contain at least one positive value.")

    if torch_module is None:
        torch_module = require("torch", feature="STalign image preparation")

    # Create coordinate grid
    xgrid = [
        torch_module.arange(image_size[0], dtype=torch_module.float32),
        torch_module.arange(image_size[1], dtype=torch_module.float32),
    ]

    # Rasterize with Gaussian smoothing (vectorized for 40-500x speedup)
    from scipy.ndimage import gaussian_filter

    # Vectorized coordinate mapping
    image_coords = transform.to_image(coords)
    row_indices = np.clip(
        np.rint(image_coords[:, 0]).astype(np.intp), 0, image_size[0] - 1
    )
    column_indices = np.clip(
        np.rint(image_coords[:, 1]).astype(np.intp), 0, image_size[1] - 1
    )

    # Accumulate intensities (np.add.at handles duplicate indices correctly)
    image = np.zeros(image_size, dtype=np.float32)
    np.add.at(image, (row_indices, column_indices), intensity)

    # Apply Gaussian filter (sigma=1.0 approximates the original kernel radius 2)
    image = gaussian_filter(image, sigma=1.0)

    # Normalize
    if image.max() > 0:
        image /= image.max()

    image_tensor = torch_module.tensor(image[None, ...], dtype=torch_module.float32)
    return xgrid, image_tensor, transform


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
    from scipy import sparse

    slices: list["ad.AnnData"] = []
    for adata in adata_list:
        s = adata[:, common_genes].copy()
        if "counts" in s.layers:
            s.X = s.layers["counts"].copy()
        elif s.raw is not None and s.raw.shape == s.shape:
            s.X = s.raw[:, common_genes].X.copy()
        if sparse.issparse(s.X):
            s.X.data = np.clip(s.X.data, 0, None)
        else:
            s.X = np.clip(np.asarray(s.X), 0, None)
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
    *,
    pst: Any | None = None,
) -> list["ad.AnnData"]:
    """Register slices using PASTE optimal transport."""
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

    if pst is None:
        pst = require("paste", feature="PASTE spatial registration")
    _patch_paste_line_search_for_pot(pst)

    # Unified preparation: gene subset + normalize + ensure obsm["spatial"]
    slices = _prepare_paste_slices(registered, common_genes, spatial_key)
    backend = get_ot_backend(params.use_gpu)

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
            backend=backend,
            use_gpu=params.use_gpu,
            verbose=False,
            gpu_verbose=False,
        )

        # Stack and extract aligned coordinates
        aligned = pst.stack_slices_pairwise([slices[ref_idx], slices[other_idx]], [pi])
        registered[ref_idx].obsm["spatial_registered"] = aligned[0].obsm["spatial"]
        registered[other_idx].obsm["spatial_registered"] = aligned[1].obsm["spatial"]

    else:
        # Multi-slice center alignment
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
                    numItermax=params.paste_numItermax,
                    backend=backend,
                    use_gpu=params.use_gpu,
                    verbose=False,
                    gpu_verbose=False,
                )
                pis.append(pi)

        # Center alignment
        center_slice, pis_new = pst.center_align(
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
        if len(pis_new) != len(registered):
            raise ProcessingError(
                "PASTE center alignment returned a different number of transport "
                f"plans than input slices: expected {len(registered)}, "
                f"received {len(pis_new)}."
            )

        # PASTE transport plans map center rows to slice columns. Its own
        # Procrustes projection handles that orientation and unequal spot counts.
        _, aligned_slices = pst.stack_slices_center(
            center_slice,
            slices,
            pis_new,
        )
        if len(aligned_slices) != len(registered):
            raise ProcessingError(
                "PASTE center projection returned a different number of slices "
                f"than requested: expected {len(registered)}, "
                f"received {len(aligned_slices)}."
            )
        for adata, aligned in zip(registered, aligned_slices, strict=True):
            adata.obsm["spatial_registered"] = np.asarray(
                aligned.obsm["spatial"], dtype=np.float32
            )

    return registered


def _register_stalign(
    adata_list: list["ad.AnnData"],
    params: RegistrationParameters,
    spatial_key: str = "spatial",
    *,
    stalign_module: Any | None = None,
    torch_module: Any | None = None,
) -> list["ad.AnnData"]:
    """Register slices using STalign diffeomorphic mapping."""
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

    if stalign_module is None:
        stalign_module = require_module(
            "stalign",
            "STalign.STalign",
            feature="STalign spatial registration",
        )
    if torch_module is None:
        torch_module = require("torch", feature="STalign spatial registration")

    # Prepare images
    image_size = (
        params.stalign_image_size[0],
        params.stalign_image_size[1],
    )
    source_grid, source_image, source_transform = _prepare_stalign_image(
        source_coords,
        source_intensity,
        image_size,
        torch_module=torch_module,
    )
    target_grid, target_image, target_transform = _prepare_stalign_image(
        target_coords,
        target_intensity,
        image_size,
        torch_module=torch_module,
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
        "dtype": torch_module.float32,
    }

    try:
        result = stalign_module.LDDMM(
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
        source_points = torch_module.tensor(
            source_transform.to_image(source_coords),
            dtype=torch_module.float32,
            device=device,
        )
        transformed = stalign_module.transform_points_source_to_target(
            xv, v, A, source_points
        )

        if isinstance(transformed, torch_module.Tensor):
            if hasattr(transformed, "detach"):
                transformed = transformed.detach()
            if hasattr(transformed, "cpu"):
                transformed = transformed.cpu()
            if hasattr(transformed, "numpy"):
                transformed = transformed.numpy()

        transformed = target_transform.from_image(
            np.asarray(transformed, dtype=np.float32)
        )
        source.obsm["spatial_registered"] = transformed
        target.obsm["spatial_registered"] = np.array(target_coords, copy=True)

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(
            f"STalign registration failed: {e}. Consider using PASTE method."
        ) from e

    return registered


# =============================================================================
# Public API
# =============================================================================


def register_slices(
    adata_list: list["ad.AnnData"],
    params: Optional[RegistrationParameters] = None,
    *,
    ctx: Optional["ToolContext"] = None,
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
        pst = require("paste", ctx, feature="PASTE spatial registration")
        return _register_paste(adata_list, params, spatial_key, pst=pst)
    if params.method == "stalign":
        stalign_module = require_module(
            "stalign",
            "STalign.STalign",
            ctx,
            feature="STalign spatial registration",
        )
        torch_module = require("torch", ctx, feature="STalign spatial registration")
        return _register_stalign(
            adata_list,
            params,
            spatial_key,
            stalign_module=stalign_module,
            torch_module=torch_module,
        )
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

    if params.method not in {"paste", "stalign"}:
        raise ParameterError(f"Unknown method: {params.method}")

    # Get data
    source_adata = await ctx.get_adata(source_id)
    target_adata = await ctx.get_adata(target_id)

    try:
        registered = register_slices(
            [source_adata, target_adata],
            params,
            ctx=ctx,
        )

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
        method_params: dict[str, object]
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

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"Registration failed: {e}") from e
