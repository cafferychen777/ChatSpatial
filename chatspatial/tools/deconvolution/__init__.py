"""
Deconvolution module for spatial transcriptomics data.

This module provides a unified interface for multiple deconvolution methods:
- flashdeconv: Ultra-fast deconvolution with O(N) complexity (recommended)
- cell2location: Bayesian deconvolution with spatial priors
- destvi: Deep learning-based multi-resolution deconvolution
- stereoscope: Two-stage probabilistic deconvolution
- rctd: Robust Cell Type Decomposition (R/spacexr or Python/PyTorch)
- spotlight: NMF-based deconvolution (R-based)
- card: CAR model with spatial correlation (R-based)
- tangram: Deep learning mapping via scvi-tools

Usage:
    from chatspatial.tools.deconvolution import deconvolve_spatial_data
    result = await deconvolve_spatial_data(data_id, ctx, params)
"""

import gc
import importlib
from dataclasses import replace
from typing import TYPE_CHECKING, Any, get_args

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.analysis import DeconvolutionResult
from ...models.data import DeconvolutionParameters
from ...utils.adata_utils import (
    ensure_unique_var_names,
    store_analysis_metadata,
    validate_obs_column,
)
from ...utils.async_utils import run_sync
from ...utils.exceptions import (
    DataError,
    DependencyError,
    ParameterError,
    ProcessingError,
)
from ...utils.results_export import export_analysis_result
from .base import MethodConfig, PreparedDeconvolutionData, prepare_deconvolution


def _build_deconvolution_key(
    method: str,
    reference_data_id: str | None,
) -> str:
    """Build a parametric analysis key for deconvolution results.

    Encodes method + reference so that runs with different reference
    datasets coexist in metadata/export/cache.
    """
    if reference_data_id:
        return f"deconvolution_{method}_{reference_data_id}"
    return f"deconvolution_{method}"


# Export main function and data container
__all__ = ["deconvolve_spatial_data", "PreparedDeconvolutionData", "METHOD_REGISTRY"]


# =============================================================================
# Method Registry - Single Source of Truth
# =============================================================================
#
# Each method is declaratively configured here. To add a new method:
# 1. Add a MethodConfig entry to METHOD_REGISTRY
# 2. Create the method module with a deconvolve() function
# 3. Add parameters to DeconvolutionParameters in models/data.py
#
# The dispatch logic does NOT need to be modified.
#
# param_mapping format: tuple of (DeconvolutionParameters_field, function_arg)

METHOD_REGISTRY: dict[str, MethodConfig] = {
    "flashdeconv": MethodConfig(
        module_name="flashdeconv",
        dependencies=("flashdeconv",),
        param_mapping=(
            ("flashdeconv_sketch_dim", "sketch_dim"),
            ("flashdeconv_lambda_spatial", "lambda_spatial"),
            ("flashdeconv_n_hvg", "n_hvg"),
            ("flashdeconv_n_markers_per_type", "n_markers_per_type"),
        ),
    ),
    "cell2location": MethodConfig(
        module_name="cell2location",
        dependencies=("cell2location", "torch"),
        requires_counts=True,
        supports_gpu=True,
        param_mapping=(
            ("cell2location_ref_model_epochs", "ref_model_epochs"),
            ("cell2location_n_epochs", "n_epochs"),
            ("cell2location_n_cells_per_spot", "n_cells_per_spot"),
            ("cell2location_detection_alpha", "detection_alpha"),
            ("cell2location_batch_key", "batch_key"),
            ("cell2location_categorical_covariate_keys", "categorical_covariate_keys"),
            ("cell2location_ref_model_lr", "ref_model_lr"),
            ("cell2location_lr", "cell2location_lr"),
            ("cell2location_ref_model_train_size", "ref_model_train_size"),
            ("cell2location_train_size", "cell2location_train_size"),
            ("cell2location_early_stopping", "early_stopping"),
            ("cell2location_early_stopping_patience", "early_stopping_patience"),
            ("cell2location_early_stopping_threshold", "early_stopping_threshold"),
            ("cell2location_use_aggressive_training", "use_aggressive_training"),
            ("cell2location_validation_size", "validation_size"),
        ),
    ),
    "destvi": MethodConfig(
        module_name="destvi",
        dependencies=("scvi", "torch"),
        requires_counts=True,
        supports_gpu=True,
        param_mapping=(
            ("destvi_n_epochs", "n_epochs"),
            ("destvi_n_hidden", "n_hidden"),
            ("destvi_n_latent", "n_latent"),
            ("destvi_n_layers", "n_layers"),
            ("destvi_dropout_rate", "dropout_rate"),
            ("destvi_learning_rate", "learning_rate"),
            ("destvi_train_size", "train_size"),
            ("destvi_vamp_prior_p", "vamp_prior_p"),
            ("destvi_l1_reg", "l1_reg"),
        ),
    ),
    "stereoscope": MethodConfig(
        module_name="stereoscope",
        dependencies=("scvi", "torch"),
        requires_counts=True,
        supports_gpu=True,
        param_mapping=(
            ("stereoscope_n_epochs", "n_epochs"),
            ("stereoscope_learning_rate", "learning_rate"),
            ("stereoscope_batch_size", "batch_size"),
        ),
    ),
    "rctd": MethodConfig(
        module_name="rctd",
        dependencies=("rpy2",),
        is_r_based=True,
        requires_counts=True,
        param_mapping=(
            ("rctd_backend", "backend"),
            ("rctd_mode", "mode"),
            ("rctd_max_cores", "max_cores"),
            ("rctd_confidence_threshold", "confidence_threshold"),
            ("rctd_doublet_threshold", "doublet_threshold"),
            ("rctd_max_multi_types", "max_multi_types"),
            ("rctd_device", "device"),
            ("rctd_batch_size", "batch_size"),
            ("rctd_dtype", "dtype"),
            ("rctd_sigma_override", "sigma_override"),
        ),
    ),
    "spotlight": MethodConfig(
        module_name="spotlight",
        dependencies=("rpy2",),
        is_r_based=True,
        requires_counts=True,
        param_mapping=(
            ("spotlight_n_top_genes", "n_top_genes"),
            ("spotlight_nmf_model", "nmf_model"),
            ("spotlight_min_prop", "min_prop"),
            ("spotlight_scale", "scale"),
            ("spotlight_weight_id", "weight_id"),
        ),
    ),
    "card": MethodConfig(
        module_name="card",
        dependencies=("rpy2",),
        is_r_based=True,
        requires_counts=True,
        param_mapping=(
            ("card_sample_key", "sample_key"),
            ("card_minCountGene", "minCountGene"),
            ("card_minCountSpot", "minCountSpot"),
            ("card_imputation", "imputation"),
            ("card_NumGrids", "NumGrids"),
            ("card_ineibor", "ineibor"),
        ),
    ),
    "tangram": MethodConfig(
        module_name="tangram",
        dependencies=("tangram",),
        supports_gpu=True,
        param_mapping=(
            ("tangram_n_epochs", "n_epochs"),
            ("tangram_mode", "mode"),
            ("tangram_learning_rate", "learning_rate"),
            ("tangram_density_prior", "density_prior"),
        ),
    ),
}

# ---------------------------------------------------------------------------
# Structural consistency: registry keys must equal Literal values in the
# parameter model.  Any mismatch crashes at import time, preventing silent
# drift between the two sources.
# ---------------------------------------------------------------------------
_literal_methods = set(
    get_args(DeconvolutionParameters.model_fields["method"].annotation)
)
assert set(METHOD_REGISTRY) == _literal_methods, (
    f"METHOD_REGISTRY keys and DeconvolutionParameters.method Literal are out "
    f"of sync.  Registry-only: {set(METHOD_REGISTRY) - _literal_methods}  "
    f"Literal-only: {_literal_methods - set(METHOD_REGISTRY)}"
)
del _literal_methods  # clean up module namespace


# =============================================================================
# Main Entry Point
# =============================================================================


async def deconvolve_spatial_data(
    data_id: str,
    ctx: "ToolContext",
    params: DeconvolutionParameters,
) -> DeconvolutionResult:
    """Deconvolve spatial transcriptomics data to estimate cell type proportions.

    This is the main entry point for all deconvolution methods. It handles:
    - Data loading and validation
    - Method selection and dependency checking
    - Dispatching to the appropriate method-specific implementation
    - Result storage and formatting

    Args:
        data_id: Dataset ID for spatial data
        ctx: Tool context for data access and logging
        params: Deconvolution parameters (must include method and cell_type_key)

    Returns:
        DeconvolutionResult with cell type proportions and statistics
    """
    # Validate input
    if not data_id:
        raise ParameterError("Dataset ID cannot be empty")

    method = params.method
    if method not in METHOD_REGISTRY:
        raise ParameterError(
            f"Unsupported method: {method}. " f"Supported: {', '.join(METHOD_REGISTRY)}"
        )

    config = _resolve_runtime_config(params, METHOD_REGISTRY[method])

    # Get spatial data
    spatial_adata = await ctx.get_adata(data_id)
    if spatial_adata.n_obs == 0:
        raise DataError(f"Dataset {data_id} contains no observations")

    # Load reference data (all methods require it)
    if not params.reference_data_id:
        raise ParameterError(f"Method '{method}' requires reference_data_id.")

    reference_adata = await ctx.get_adata(params.reference_data_id)
    if reference_adata.n_obs == 0:
        raise DataError(
            f"Reference dataset {params.reference_data_id} contains no observations"
        )

    validate_obs_column(reference_adata, params.cell_type_key, "Cell type")

    # Check method availability
    _check_method_availability(method, config)

    # Prepare data using unified function
    preprocess_hook = (
        _get_preprocess_hook(params) if method == "cell2location" else None
    )

    data = await prepare_deconvolution(
        spatial_adata=spatial_adata,
        reference_adata=reference_adata,
        cell_type_key=params.cell_type_key,
        ctx=ctx,
        require_int_dtype=config.is_r_based,
        require_counts=config.requires_counts,
        preprocess=preprocess_hook,
    )

    # Dispatch to method-specific implementation
    proportions, stats = await run_sync(_dispatch_method, data, params, config)

    # Memory cleanup
    del data
    gc.collect()

    # Publish onto a fresh copy only after the backend has completed. This keeps
    # preparation, backend, metadata, and export failures from contaminating the
    # managed dataset while avoiding another full copy during model execution.
    result_adata = await run_sync(spatial_adata.copy)
    ensure_unique_var_names(result_adata, "spatial data")

    # Store results in AnnData
    result = await _store_results(
        result_adata,
        proportions,
        stats,
        method,
        data_id,
        ctx,
        reference_data_id=params.reference_data_id,
    )

    return result


# =============================================================================
# Internal Functions
# =============================================================================


def _check_method_availability(method: str, config: MethodConfig) -> None:
    """Check if a deconvolution method is available."""
    import importlib.util

    missing = []
    for dep in config.dependencies:
        import_name = _dependency_import_name(dep)
        if importlib.util.find_spec(import_name) is None:
            missing.append(dep)

    if missing:
        # Find available methods for helpful error message
        available = []
        for name, cfg in METHOD_REGISTRY.items():
            check_deps = [_dependency_import_name(x) for x in cfg.dependencies]
            if all(importlib.util.find_spec(x) is not None for x in check_deps):
                available.append(name)

        alt_msg = f"Available: {', '.join(available)}" if available else ""
        if "flashdeconv" in available:
            alt_msg += " (flashdeconv recommended - fastest)"

        raise DependencyError(
            f"Method '{method}' requires: {', '.join(missing)}. {alt_msg}"
        )


def _dependency_import_name(dependency: str) -> str:
    """Return the importable module provided by an optional distribution."""
    return {
        "scvi-tools": "scvi",
        "rctd-py": "rctd",
    }.get(dependency, dependency.replace("-", "_"))


def _resolve_runtime_config(
    params: DeconvolutionParameters,
    config: MethodConfig,
) -> MethodConfig:
    """Select backend-specific dependencies and count preparation settings."""
    if params.method == "rctd" and params.rctd_backend == "python":
        return replace(
            config,
            dependencies=("rctd-py",),
            is_r_based=False,
        )
    return config


def _get_preprocess_hook(params: DeconvolutionParameters):
    """Get cell2location-specific preprocessing hook."""
    if not params.cell2location_apply_gene_filtering:
        return None

    async def cell2location_preprocess(spatial, reference, ctx):
        from .cell2location import apply_gene_filtering

        sp = await apply_gene_filtering(
            spatial,
            ctx,
            cell_count_cutoff=params.cell2location_gene_filter_cell_count_cutoff,
            cell_percentage_cutoff2=params.cell2location_gene_filter_cell_percentage_cutoff2,
            nonz_mean_cutoff=params.cell2location_gene_filter_nonz_mean_cutoff,
        )
        ref = await apply_gene_filtering(
            reference,
            ctx,
            cell_count_cutoff=params.cell2location_gene_filter_cell_count_cutoff,
            cell_percentage_cutoff2=params.cell2location_gene_filter_cell_percentage_cutoff2,
            nonz_mean_cutoff=params.cell2location_gene_filter_nonz_mean_cutoff,
        )
        return sp, ref

    return cell2location_preprocess


def _dispatch_method(
    data: PreparedDeconvolutionData,
    params: DeconvolutionParameters,
    config: MethodConfig,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Dispatch to the appropriate method implementation.

    This function replaces the 124-line if-elif chain with a simple
    registry lookup and dynamic import.
    """
    # Dynamic import of method module
    module = importlib.import_module(f".{config.module_name}", package=__package__)
    deconvolve_func = module.deconvolve

    # Extract method-specific kwargs using declarative mapping
    kwargs = config.extract_kwargs(params)

    return deconvolve_func(data, **kwargs)


def _validate_and_align_proportions(
    proportions: pd.DataFrame,
    obs_names: pd.Index,
    method: str,
) -> pd.DataFrame:
    """Validate backend proportions and align rows to the spatial observation index."""
    if not isinstance(proportions, pd.DataFrame):
        raise ProcessingError(
            f"{method} returned {type(proportions).__name__}, expected a DataFrame."
        )
    if len(obs_names) == 0:
        raise DataError("Cannot store deconvolution results for an empty dataset.")
    if proportions.shape[1] == 0:
        raise ProcessingError(f"{method} returned no cell-type proportion columns.")

    normalized = proportions.copy()
    normalized.index = normalized.index.map(str)
    normalized.columns = normalized.columns.map(str)
    expected_index = pd.Index(obs_names.map(str))

    if not normalized.index.is_unique:
        raise ProcessingError(f"{method} returned duplicate spatial observation IDs.")
    if not normalized.columns.is_unique:
        raise ProcessingError(f"{method} returned duplicate cell-type labels.")

    missing = expected_index.difference(normalized.index)
    unexpected = normalized.index.difference(expected_index)
    if len(missing) or len(unexpected):
        raise ProcessingError(
            f"{method} result observations do not match the spatial dataset: "
            f"{len(missing)} missing and {len(unexpected)} unexpected IDs."
        )

    normalized = normalized.reindex(expected_index)
    raw_values = normalized.to_numpy()
    if np.iscomplexobj(raw_values):
        raise ProcessingError(f"{method} returned complex-valued proportions.")
    try:
        values = np.asarray(raw_values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ProcessingError(f"{method} returned non-numeric proportions.") from exc
    if not np.isfinite(values).all():
        raise ProcessingError(f"{method} returned non-finite proportions.")

    minimum = float(values.min())
    if minimum < -1e-10:
        raise ProcessingError(
            f"{method} returned negative proportions (minimum={minimum:.3g})."
        )
    values = np.maximum(values, 0.0)
    return pd.DataFrame(values, index=expected_index, columns=normalized.columns)


_RCTD_AUX_OBS_KEYS = frozenset(
    {
        "rctd_status",
        "rctd_converged",
        "rctd_spot_class",
        "rctd_first_type",
        "rctd_second_type",
        "rctd_first_class",
        "rctd_second_class",
        "rctd_first_class_name",
        "rctd_second_class_name",
        "rctd_min_score",
        "rctd_singlet_score",
        "rctd_n_types",
    }
)
_RCTD_AUX_OBSM_KEYS = frozenset(
    {
        "rctd_doublet_weights",
        "rctd_sub_weights",
        "rctd_cell_type_indices",
        "rctd_confident_assignments",
    }
)


def _clear_rctd_backend_outputs(spatial_adata: "ad.AnnData") -> None:
    """Remove mode-specific RCTD fields before replacing an existing result."""
    for key in _RCTD_AUX_OBS_KEYS:
        if key in spatial_adata.obs:
            del spatial_adata.obs[key]
    for key in _RCTD_AUX_OBSM_KEYS:
        if key in spatial_adata.obsm:
            del spatial_adata.obsm[key]
    if "rctd_backend" in spatial_adata.uns:
        del spatial_adata.uns["rctd_backend"]


def _store_rctd_backend_outputs(
    spatial_adata: "ad.AnnData",
    outputs: dict[str, dict[str, Any]],
    stats: dict[str, Any],
) -> dict[str, list[str]]:
    """Store compact rctd-py auxiliary results and backend provenance."""
    stored: dict[str, list[str]] = {"obs": [], "obsm": [], "uns": []}
    n_obs = spatial_adata.n_obs

    for key, values in outputs.get("obs", {}).items():
        array = np.asarray(values)
        if array.shape != (n_obs,):
            raise ProcessingError(
                f"RCTD auxiliary obs field '{key}' has shape {array.shape}; "
                f"expected {(n_obs,)}."
            )
        if array.dtype.kind in {"O", "S", "U"}:
            spatial_adata.obs[key] = pd.Categorical(array.astype(str))
        else:
            spatial_adata.obs[key] = array
        stored["obs"].append(key)

    for key, values in outputs.get("obsm", {}).items():
        array = np.asarray(values)
        if array.ndim != 2 or array.shape[0] != n_obs:
            raise ProcessingError(
                f"RCTD auxiliary obsm field '{key}' has shape {array.shape}; "
                f"expected ({n_obs}, n_features)."
            )
        spatial_adata.obsm[key] = array
        stored["obsm"].append(key)

    provenance_keys = (
        "backend",
        "backend_package",
        "backend_version",
        "mode",
        "device",
        "requested_device",
        "batch_size",
        "dtype",
        "sigma_override",
        "confidence_threshold",
        "doublet_threshold",
        "max_multi_types",
        "max_cores",
        "n_filtered_spots",
    )
    provenance = {
        key: stats[key]
        for key in provenance_keys
        if key in stats and stats[key] is not None
    }
    spatial_adata.uns["rctd_backend"] = provenance
    stored["uns"].append("rctd_backend")
    return stored


async def _store_results(
    spatial_adata: "ad.AnnData",
    proportions: pd.DataFrame,
    stats: dict[str, Any],
    method: str,
    data_id: str,
    ctx: "ToolContext",
    *,
    reference_data_id: str | None = None,
) -> DeconvolutionResult:
    """Store deconvolution results in AnnData and return result object."""
    proportions_key = f"deconvolution_{method}"
    aligned = _validate_and_align_proportions(
        proportions, spatial_adata.obs_names, method
    )
    cell_types = aligned.columns.tolist()
    full_proportions = aligned.to_numpy(copy=True)
    result_stats = dict(stats)
    backend_outputs = result_stats.pop("_backend_outputs", None)
    result_stats["n_spots"] = len(full_proportions)
    result_stats["n_cell_types"] = len(cell_types)

    auxiliary_keys: dict[str, list[str]] = {"obs": [], "obsm": [], "uns": []}
    if method == "rctd":
        _clear_rctd_backend_outputs(spatial_adata)

    # Store in obsm
    spatial_adata.obsm[proportions_key] = full_proportions

    # Store cell type names
    spatial_adata.uns[f"{proportions_key}_cell_types"] = cell_types

    # Add individual cell type columns to obs
    for i, ct in enumerate(cell_types):
        spatial_adata.obs[f"{proportions_key}_{ct}"] = full_proportions[:, i]

    # Add dominant cell type annotation (all-zero rows → "unassigned")
    dominant_key = f"dominant_celltype_{method}"
    cell_types_array = np.array(cell_types, dtype=object)
    row_sums = full_proportions.sum(axis=1)
    zero_mask = row_sums == 0
    dominant_types = cell_types_array[np.argmax(full_proportions, axis=1)]
    dominant_types[zero_mask] = "unassigned"
    spatial_adata.obs[dominant_key] = pd.Categorical(dominant_types)

    # Note: "unassigned" is a label in obs[dominant_key] for zero-sum rows,
    # NOT a column in the proportions matrix.  Metadata cell_types must
    # match the proportions matrix columns exactly to avoid downstream
    # DataFrame construction errors.
    has_unassigned = bool(zero_mask.any())

    if method == "rctd":
        auxiliary_keys = _store_rctd_backend_outputs(
            spatial_adata,
            backend_outputs or {},
            result_stats,
        )

    # Store metadata for provenance tracking
    analysis_key = _build_deconvolution_key(method, reference_data_id)
    obs_keys: list[str] = [dominant_key, *auxiliary_keys["obs"]]
    store_analysis_metadata(
        spatial_adata,
        analysis_name=analysis_key,
        method=method,
        parameters={},  # Method-specific params already in stats
        results_keys={
            "obsm": [proportions_key, *auxiliary_keys["obsm"]],
            "obs": obs_keys,
            "uns": [f"{proportions_key}_cell_types", *auxiliary_keys["uns"]],
        },
        statistics={
            "n_cell_types": len(cell_types),
            "n_spots": len(full_proportions),
            "cell_types": cell_types,
            "proportions_key": proportions_key,
            "dominant_type_key": dominant_key,
            "has_unassigned_spots": has_unassigned,
        },
    )

    # Store CARD imputation data if present (bridge analysis → visualization)
    imputation = result_stats.get("imputation")
    if imputation and method == "card":
        n_original = len(full_proportions)
        n_imputed = imputation.get("n_imputed_locations", 0)
        spatial_adata.uns["card_imputation"] = {
            "proportions": imputation["imputed_proportions"],
            "coordinates": imputation["imputed_coordinates"],
            "resolution_increase": n_imputed / max(n_original, 1),
        }

    # Export results to CSV for reproducibility
    export_analysis_result(spatial_adata, data_id, analysis_key)

    result = DeconvolutionResult(
        data_id=data_id,
        method=method,
        dominant_type_key=dominant_key,
        n_cell_types=len(cell_types),
        cell_types=cell_types,
        proportions_key=proportions_key,
        n_spots=len(full_proportions),
        genes_used=result_stats.get("genes_used", result_stats.get("common_genes", 0)),
        statistics=result_stats,
    )
    await ctx.set_adata(data_id, spatial_adata)
    return result
