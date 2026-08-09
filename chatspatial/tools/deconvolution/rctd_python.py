"""PyTorch-backed RCTD adapter using the optional rctd-py package."""

import importlib
import os
import shutil
import ssl
import tempfile
import threading
import urllib.request
import warnings
import zipfile
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import certifi
import numpy as np
import pandas as pd

from ...utils.dependency_manager import require
from ...utils.exceptions import (
    ChatSpatialError,
    DataError,
    DependencyError,
    ParameterError,
    ProcessingError,
)
from .base import PreparedDeconvolutionData, create_deconvolution_stats

_VALID_MODES = frozenset({"full", "doublet", "multi"})
_MIN_CELLS_PER_TYPE = 25
_CACHE_DOWNLOAD_LOCK = threading.Lock()


def _is_valid_cache(path: Path) -> bool:
    """Return whether a likelihood cache looks like a complete NPZ archive."""
    return path.is_file() and path.stat().st_size > 0 and zipfile.is_zipfile(path)


def _ensure_likelihood_cache(
    rctd_module: Any,
    cache_path: Path | None = None,
) -> None:
    """Download the rctd-py likelihood cache with a reliable CA bundle."""
    module_name = getattr(rctd_module, "__name__", None)
    if module_name is None:
        return

    target = cache_path or Path.home() / ".cache" / "rctd" / "q_matrices.npz"
    if _is_valid_cache(target):
        return

    with _CACHE_DOWNLOAD_LOCK:
        if _is_valid_cache(target):
            return

        try:
            likelihood = importlib.import_module(f"{module_name}._likelihood")
            download_url = str(likelihood._Q_MATRICES_URL)
        except (ImportError, AttributeError) as exc:
            raise DependencyError(
                "The installed rctd-py release does not expose its likelihood "
                "cache URL. Install a compatible release with: "
                "pip install 'chatspatial[rctd-python]'"
            ) from exc

        target.parent.mkdir(parents=True, exist_ok=True)
        temporary_path: Path | None = None
        try:
            tls_context = ssl.create_default_context(cafile=certifi.where())
            print(f"Downloading rctd-py likelihood cache ({download_url}) ...")
            with (
                urllib.request.urlopen(
                    download_url,
                    timeout=60,
                    context=tls_context,
                ) as response,
                tempfile.NamedTemporaryFile(
                    mode="wb",
                    prefix=".q_matrices.",
                    suffix=".part",
                    dir=target.parent,
                    delete=False,
                ) as temporary_file,
            ):
                temporary_path = Path(temporary_file.name)
                shutil.copyfileobj(response, temporary_file)

            if not _is_valid_cache(temporary_path):
                raise ProcessingError(
                    "Downloaded rctd-py likelihood cache is not a valid NPZ archive."
                )
            os.replace(temporary_path, target)
            temporary_path = None
        except ChatSpatialError:
            raise
        except Exception as exc:
            raise ProcessingError(
                f"Could not download the rctd-py likelihood cache: {exc}"
            ) from exc
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)


def _pixel_mask(result: Any, n_spots: int) -> np.ndarray:
    """Return and validate the upstream mask for analyzed spatial locations."""
    raw_mask = getattr(result, "pixel_mask", None)
    if raw_mask is None:
        return np.ones(n_spots, dtype=bool)

    mask = np.asarray(raw_mask, dtype=bool)
    if mask.shape != (n_spots,):
        raise ProcessingError(
            "rctd-py returned an invalid pixel mask: "
            f"expected {(n_spots,)}, received {mask.shape}."
        )
    return mask


def _expand_rows(
    values: Any,
    mask: np.ndarray,
    *,
    fill_value: Any,
    dtype: Any | None = None,
) -> np.ndarray:
    """Expand values for analyzed spots back to the full spatial observation set."""
    array = np.asarray(values, dtype=dtype)
    expected_rows = int(mask.sum())
    if array.ndim == 0 or array.shape[0] != expected_rows:
        received = array.shape[0] if array.ndim else 0
        raise ProcessingError(
            "rctd-py returned results with an invalid row count: "
            f"expected {expected_rows}, received {received}."
        )
    full_shape = (len(mask), *array.shape[1:])
    expanded = np.full(full_shape, fill_value, dtype=array.dtype)
    expanded[mask] = array
    return expanded


def _type_names(indices: Any, cell_types: list[str]) -> np.ndarray:
    """Map upstream cell-type indices to labels without unsafe negative indexing."""
    raw = np.asarray(indices, dtype=int)
    labels = np.full(raw.shape, "unassigned", dtype=object)
    valid = (raw >= 0) & (raw < len(cell_types))
    if valid.any():
        names = np.asarray(cell_types, dtype=object)
        labels[valid] = names[raw[valid]]
    return labels


def _package_version(rctd_module: Any) -> str:
    """Return the installed rctd-py distribution version when available."""
    try:
        return version("rctd-py")
    except PackageNotFoundError:
        return str(getattr(rctd_module, "__version__", "unknown"))


def _resolved_device(device: str) -> str:
    """Resolve the auto device for provenance without changing backend behavior."""
    if device != "auto":
        return device
    try:
        torch = importlib.import_module("torch")
        return "cuda" if torch.cuda.is_available() else "cpu"
    except (ImportError, AttributeError):
        return "auto"


def _build_backend_outputs(
    result: Any,
    mode: str,
    mask: np.ndarray,
    cell_types: list[str],
    spot_class_names: list[str],
) -> dict[str, dict[str, Any]]:
    """Convert mode-specific rctd-py results into AnnData-ready annotations."""
    outputs: dict[str, dict[str, Any]] = {
        "obs": {
            "rctd_status": np.where(mask, "analyzed", "filtered"),
        },
        "obsm": {},
    }

    if mode == "full":
        outputs["obs"]["rctd_converged"] = _expand_rows(
            result.converged,
            mask,
            fill_value=False,
            dtype=bool,
        )
        return outputs

    if mode == "doublet":
        spot_class = np.asarray(result.spot_class, dtype=int)
        class_labels = np.full(spot_class.shape, "unknown", dtype=object)
        valid_classes = (spot_class >= 0) & (spot_class < len(spot_class_names))
        if valid_classes.any():
            names = np.asarray(spot_class_names, dtype=object)
            class_labels[valid_classes] = names[spot_class[valid_classes]]

        outputs["obs"].update(
            {
                "rctd_spot_class": _expand_rows(
                    class_labels, mask, fill_value="filtered", dtype=object
                ),
                "rctd_first_type": _expand_rows(
                    _type_names(result.first_type, cell_types),
                    mask,
                    fill_value="filtered",
                    dtype=object,
                ),
                "rctd_second_type": _expand_rows(
                    _type_names(result.second_type, cell_types),
                    mask,
                    fill_value="filtered",
                    dtype=object,
                ),
                "rctd_first_class": _expand_rows(
                    result.first_class, mask, fill_value=False, dtype=bool
                ),
                "rctd_second_class": _expand_rows(
                    result.second_class, mask, fill_value=False, dtype=bool
                ),
                "rctd_min_score": _expand_rows(
                    result.min_score, mask, fill_value=np.nan, dtype=float
                ),
                "rctd_singlet_score": _expand_rows(
                    result.singlet_score, mask, fill_value=np.nan, dtype=float
                ),
            }
        )
        outputs["obsm"]["rctd_doublet_weights"] = _expand_rows(
            result.weights_doublet,
            mask,
            fill_value=0.0,
            dtype=float,
        )

        for attribute, key in (
            ("first_class_name", "rctd_first_class_name"),
            ("second_class_name", "rctd_second_class_name"),
        ):
            values = getattr(result, attribute, None)
            if values is not None:
                outputs["obs"][key] = _expand_rows(
                    values, mask, fill_value="filtered", dtype=object
                )
        return outputs

    outputs["obs"].update(
        {
            "rctd_n_types": _expand_rows(result.n_types, mask, fill_value=0, dtype=int),
            "rctd_min_score": _expand_rows(
                result.min_score, mask, fill_value=np.nan, dtype=float
            ),
        }
    )
    outputs["obsm"].update(
        {
            "rctd_sub_weights": _expand_rows(
                result.sub_weights, mask, fill_value=0.0, dtype=float
            ),
            "rctd_cell_type_indices": _expand_rows(
                result.cell_type_indices, mask, fill_value=-1, dtype=int
            ),
            "rctd_confident_assignments": _expand_rows(
                result.conf_list, mask, fill_value=False, dtype=bool
            ),
        }
    )
    return outputs


def deconvolve(
    data: PreparedDeconvolutionData,
    mode: str = "full",
    confidence_threshold: float = 10.0,
    doublet_threshold: float = 25.0,
    max_multi_types: int = 4,
    device: str = "auto",
    batch_size: int | str = "auto",
    dtype: str = "float64",
    sigma_override: int | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Deconvolve spatial counts with the rctd-py PyTorch implementation."""
    if mode not in _VALID_MODES:
        choices = ", ".join(sorted(_VALID_MODES))
        raise ParameterError(f"Invalid RCTD mode '{mode}'. Expected one of: {choices}.")
    if data.spatial_coords is None:
        raise DataError(
            "RCTD requires real spatial coordinates for spatially-informed "
            "deconvolution. No spatial coordinates found. "
            "Use a non-spatial method (e.g., NNLS, DestVI) instead."
        )

    try:
        reference_data = data.reference
        cell_types = reference_data.obs[data.cell_type_key].astype(str)
        type_counts = cell_types.value_counts()
        rare_types = type_counts[type_counts < _MIN_CELLS_PER_TYPE].index.tolist()

        if rare_types:
            warnings.warn(
                f"RCTD requires ≥{_MIN_CELLS_PER_TYPE} cells per cell type. "
                f"Filtering {len(rare_types)} rare types: {rare_types}",
                UserWarning,
                stacklevel=2,
            )
            keep_mask = ~cell_types.isin(rare_types)
            reference_data = reference_data[keep_mask].copy()
            cell_types = cell_types[keep_mask]

        remaining_types = sorted(cell_types.unique().tolist())
        if len(remaining_types) < 2:
            raise DataError(
                f"After filtering rare cell types, only {len(remaining_types)} "
                "cell type(s) remain. RCTD requires at least 2 cell types."
            )
        if mode == "multi" and max_multi_types >= len(remaining_types):
            raise ParameterError(
                f"MAX_MULTI_TYPES ({max_multi_types}) must be less than "
                f"total cell types ({len(remaining_types)})."
            )

        reference_data.obs[data.cell_type_key] = cell_types
        rctd = require("rctd-py", feature="RCTD Python deconvolution")
        required_api = {"Reference", "RCTDConfig", "run_rctd", "SPOT_CLASS_NAMES"}
        missing_api = sorted(name for name in required_api if not hasattr(rctd, name))
        if missing_api:
            raise DependencyError(
                "The installed 'rctd' module is not a compatible rctd-py "
                f"release; missing API: {', '.join(missing_api)}. "
                "Install or upgrade with: pip install 'chatspatial[rctd-python]'"
            )
        _ensure_likelihood_cache(rctd)
        reference = rctd.Reference(
            reference_data,
            cell_type_col=data.cell_type_key,
            cell_min=_MIN_CELLS_PER_TYPE,
            min_UMI=5,
        )
        if mode == "multi" and max_multi_types >= len(reference.cell_type_names):
            raise ParameterError(
                f"MAX_MULTI_TYPES ({max_multi_types}) must be less than "
                f"total cell types ({len(reference.cell_type_names)})."
            )

        config = rctd.RCTDConfig(
            UMI_min_sigma=10,
            MAX_MULTI_TYPES=max_multi_types,
            CONFIDENCE_THRESHOLD=confidence_threshold,
            DOUBLET_THRESHOLD=doublet_threshold,
            dtype=dtype,
            device=device,
        )
        result = rctd.run_rctd(
            data.spatial,
            reference,
            mode=mode,
            config=config,
            batch_size=batch_size,
            sigma_override=sigma_override,
        )

        mask = _pixel_mask(result, data.n_spots)
        cell_type_names = [str(name) for name in result.cell_type_names]
        weights = np.asarray(result.weights, dtype=float)
        expected_shape = (int(mask.sum()), len(cell_type_names))
        if weights.shape != expected_shape:
            raise ProcessingError(
                "rctd-py returned an invalid weights matrix: "
                f"expected {expected_shape}, received {weights.shape}."
            )
        if not np.isfinite(weights).all():
            raise ProcessingError("rctd-py returned non-finite weights.")
        minimum = float(weights.min()) if weights.size else 0.0
        if minimum < -1e-10:
            raise ProcessingError(
                f"rctd-py returned negative weights (minimum={minimum:.3g})."
            )

        full_weights = np.zeros((data.n_spots, len(cell_type_names)), dtype=float)
        full_weights[mask] = np.maximum(weights, 0.0)
        proportions = pd.DataFrame(
            full_weights,
            index=data.spatial.obs_names,
            columns=cell_type_names,
        )
        backend_outputs = _build_backend_outputs(
            result,
            mode,
            mask,
            cell_type_names,
            list(rctd.SPOT_CLASS_NAMES),
        )
        stats = create_deconvolution_stats(
            proportions,
            data.common_genes,
            method=f"RCTD-{mode}",
            device=_resolved_device(device).upper(),
            mode=mode,
            backend="python",
            backend_package="rctd-py",
            backend_version=_package_version(rctd),
            requested_device=device,
            batch_size=batch_size,
            dtype=dtype,
            sigma_override=sigma_override,
            confidence_threshold=confidence_threshold,
            doublet_threshold=doublet_threshold,
            max_multi_types=max_multi_types,
            n_filtered_spots=int((~mask).sum()),
            _backend_outputs=backend_outputs,
        )
        return proportions, stats

    except ChatSpatialError:
        raise
    except Exception as exc:
        raise ProcessingError(f"RCTD Python deconvolution failed: {exc}") from exc
