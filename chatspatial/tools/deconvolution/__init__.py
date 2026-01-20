"""
Deconvolution module for spatial transcriptomics data.

This module provides a unified interface for multiple deconvolution methods:
- flashdeconv: Ultra-fast deconvolution with O(N) complexity (recommended)
- cell2location: Bayesian deconvolution with spatial priors
- destvi: Deep learning-based multi-resolution deconvolution
- stereoscope: Two-stage probabilistic deconvolution
- rctd: Robust Cell Type Decomposition (R-based)
- spotlight: NMF-based deconvolution (R-based)
- card: CAR model with spatial correlation (R-based)
- tangram: Deep learning mapping via scvi-tools

Usage:
    from chatspatial.tools.deconvolution import deconvolve_spatial_data
    result = await deconvolve_spatial_data(data_id, ctx, params)
"""

import gc
import importlib
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.analysis import DeconvolutionResult
from ...models.data import DeconvolutionParameters
from ...utils.adata_utils import (
    ensure_unique_var_names_async,
    store_analysis_metadata,
    validate_obs_column,
)
from ...utils.exceptions import DataError, DependencyError, ParameterError
from ...utils.results_export import export_analysis_result
from .base import MethodConfig, PreparedDeconvolutionData, prepare_deconvolution

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
        param_mapping=(
            ("rctd_mode", "mode"),
            ("max_cores", "max_cores"),
            ("rctd_confidence_threshold", "confidence_threshold"),
            ("rctd_doublet_threshold", "doublet_threshold"),
            ("rctd_max_multi_types", "max_multi_types"),
        ),
    ),
    "spotlight": MethodConfig(
        module_name="spotlight",
        dependencies=("rpy2",),
        is_r_based=True,
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
        dependencies=("scvi", "torch", "tangram", "mudata"),
        supports_gpu=True,
        param_mapping=(
            ("tangram_n_epochs", "n_epochs"),
            ("tangram_mode", "mode"),
            ("tangram_learning_rate", "learning_rate"),
            ("tangram_density_prior", "density_prior"),
        ),
    ),
}


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

    config = METHOD_REGISTRY[method]

    # Get spatial data
    spatial_adata = await ctx.get_adata(data_id)
    if spatial_adata.n_obs == 0:
        raise DataError(f"Dataset {data_id} contains no observations")

    await ensure_unique_var_names_async(spatial_adata, ctx, "spatial data")

    # Load reference data (all methods require it)
    if not params.reference_data_id:
        raise ParameterError(f"Method '{method}' requires reference_data_id.")

    reference_adata = await ctx.get_adata(params.reference_data_id)
    if reference_adata.n_obs == 0:
        raise DataError(
            f"Reference dataset {params.reference_data_id} contains no observations"
        )

    await ensure_unique_var_names_async(reference_adata, ctx, "reference data")
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
        preprocess=preprocess_hook,
    )

    # Dispatch to method-specific implementation
    proportions, stats = _dispatch_method(data, params, config)

    # Memory cleanup
    del data
    gc.collect()

    # Store results in AnnData
    result = await _store_results(
        spatial_adata, proportions, stats, method, data_id, ctx
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
        import_name = "scvi" if dep == "scvi-tools" else dep.replace("-", "_")
        if importlib.util.find_spec(import_name) is None:
            missing.append(dep)

    if missing:
        # Find available methods for helpful error message
        available = []
        for name, cfg in METHOD_REGISTRY.items():
            check_deps = [
                "scvi" if x == "scvi-tools" else x.replace("-", "_")
                for x in cfg.dependencies
            ]
            if all(importlib.util.find_spec(x) is not None for x in check_deps):
                available.append(name)

        alt_msg = f"Available: {', '.join(available)}" if available else ""
        if "flashdeconv" in available:
            alt_msg += " (flashdeconv recommended - fastest)"

        raise DependencyError(
            f"Method '{method}' requires: {', '.join(missing)}. {alt_msg}"
        )


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


async def _store_results(
    spatial_adata: "ad.AnnData",
    proportions: pd.DataFrame,
    stats: dict[str, Any],
    method: str,
    data_id: str,
    ctx: "ToolContext",
) -> DeconvolutionResult:
    """Store deconvolution results in AnnData and return result object."""
    proportions_key = f"deconvolution_{method}"
    cell_types = list(proportions.columns)

    # Align proportions with spatial_adata.obs_names
    full_proportions = proportions.reindex(spatial_adata.obs_names).fillna(0).values

    # Store in obsm
    spatial_adata.obsm[proportions_key] = full_proportions

    # Store cell type names
    spatial_adata.uns[f"{proportions_key}_cell_types"] = cell_types

    # Add individual cell type columns to obs
    for i, ct in enumerate(cell_types):
        spatial_adata.obs[f"{proportions_key}_{ct}"] = full_proportions[:, i]

    # Add dominant cell type annotation
    dominant_key = f"dominant_celltype_{method}"
    cell_types_array = np.array(cell_types)
    dominant_types = cell_types_array[np.argmax(full_proportions, axis=1)]
    spatial_adata.obs[dominant_key] = pd.Categorical(dominant_types)

    # Store metadata for provenance tracking
    store_analysis_metadata(
        spatial_adata,
        analysis_name=f"deconvolution_{method}",
        method=method,
        parameters={},  # Method-specific params already in stats
        results_keys={
            "obsm": [proportions_key],
            "obs": [dominant_key],
            "uns": [f"{proportions_key}_cell_types"],
        },
        statistics={
            "n_cell_types": len(cell_types),
            "n_spots": len(full_proportions),
            "cell_types": cell_types,
            "proportions_key": proportions_key,
            "dominant_type_key": dominant_key,
        },
    )

    # Export results to CSV for reproducibility
    export_analysis_result(spatial_adata, data_id, f"deconvolution_{method}")

    # Save updated data
    await ctx.set_adata(data_id, spatial_adata)

    return DeconvolutionResult(
        data_id=data_id,
        method=method,
        dominant_type_key=dominant_key,
        n_cell_types=len(cell_types),
        cell_types=cell_types,
        proportions_key=proportions_key,
        n_spots=stats.get("n_spots", 0),
        genes_used=stats.get("genes_used", stats.get("common_genes", 0)),
        statistics=stats,
    )
