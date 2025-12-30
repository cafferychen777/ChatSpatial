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
from typing import TYPE_CHECKING, Any, Dict, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext

from ...models.analysis import DeconvolutionResult
from ...models.data import DeconvolutionParameters
from ...utils.adata_utils import ensure_unique_var_names_with_ctx, validate_obs_column
from ...utils.exceptions import DataError, DependencyError, ParameterError
from .base import DeconvolutionContext

# Export main function
__all__ = ["deconvolve_spatial_data", "DeconvolutionContext"]


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

    # Get spatial data
    spatial_adata = await ctx.get_adata(data_id)
    if spatial_adata.n_obs == 0:
        raise DataError(f"Dataset {data_id} contains no observations")

    await ensure_unique_var_names_with_ctx(spatial_adata, ctx, "spatial data")

    # Load reference data for methods that require it
    reference_adata = None
    if params.method in _METHODS_REQUIRING_REFERENCE:
        if not params.reference_data_id:
            raise ParameterError(
                f"Method '{params.method}' requires reference_data_id."
            )

        reference_adata = await ctx.get_adata(params.reference_data_id)
        if reference_adata.n_obs == 0:
            raise DataError(
                f"Reference dataset {params.reference_data_id} contains no observations"
            )

        await ensure_unique_var_names_with_ctx(reference_adata, ctx, "reference data")
        validate_obs_column(reference_adata, params.cell_type_key, "Cell type key")

    # Check method availability
    _check_method_availability(params.method)

    # Create deconvolution context and prepare data
    deconv_ctx = DeconvolutionContext(
        spatial_adata=spatial_adata,
        reference_adata=reference_adata,
        cell_type_key=params.cell_type_key,
        ctx=ctx,
    )

    # Prepare data (R-based methods need int32)
    require_int = params.method in _R_BASED_METHODS
    await deconv_ctx.prepare(require_int_dtype=require_int)

    # Dispatch to method-specific implementation
    proportions, stats = await _dispatch_method(deconv_ctx, params)

    # Memory cleanup: release prepared data from context after dispatch
    # This frees spatial_prepared and reference_prepared copies
    del deconv_ctx
    gc.collect()

    # Store results in AnnData
    result = await _store_results(
        spatial_adata, proportions, stats, params.method, data_id, ctx
    )

    return result


# =============================================================================
# Method Registry
# =============================================================================

_METHODS_REQUIRING_REFERENCE = {
    "flashdeconv",
    "cell2location",
    "rctd",
    "destvi",
    "stereoscope",
    "tangram",
    "spotlight",
    "card",
}

_R_BASED_METHODS = {"rctd", "spotlight", "card"}

_METHOD_DEPENDENCIES = {
    "flashdeconv": ["flashdeconv"],
    "cell2location": ["cell2location", "torch"],
    "destvi": ["scvi", "torch"],
    "stereoscope": ["scvi", "torch"],
    "tangram": ["scvi", "torch", "tangram", "mudata"],
    "rctd": ["rpy2"],
    "spotlight": ["rpy2"],
    "card": ["rpy2"],
}


def _check_method_availability(method: str) -> None:
    """Check if a deconvolution method is available."""
    import importlib.util

    deps = _METHOD_DEPENDENCIES.get(method, [])
    missing = []

    for dep in deps:
        # Handle package name variations
        import_name = "scvi" if dep == "scvi-tools" else dep.replace("-", "_")
        if importlib.util.find_spec(import_name) is None:
            missing.append(dep)

    if missing:
        # Suggest alternatives
        available = []
        for m, d in _METHOD_DEPENDENCIES.items():
            check_deps = [
                "scvi" if x == "scvi-tools" else x.replace("-", "_") for x in d
            ]
            if all(importlib.util.find_spec(x) is not None for x in check_deps):
                available.append(m)

        alt_msg = f"Available: {', '.join(available)}" if available else ""
        if "flashdeconv" in available:
            alt_msg += " (flashdeconv recommended - fastest)"

        raise DependencyError(
            f"Method '{method}' requires: {', '.join(missing)}. {alt_msg}"
        )


async def _dispatch_method(
    deconv_ctx: DeconvolutionContext,
    params: DeconvolutionParameters,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Dispatch to the appropriate method implementation."""
    method = params.method

    if method == "flashdeconv":
        from . import flashdeconv

        return await flashdeconv.deconvolve(
            deconv_ctx,
            sketch_dim=params.flashdeconv_sketch_dim,
            lambda_spatial=params.flashdeconv_lambda_spatial,
            n_hvg=params.flashdeconv_n_hvg,
            n_markers_per_type=params.flashdeconv_n_markers_per_type,
        )

    elif method == "cell2location":
        from . import cell2location

        return await cell2location.deconvolve(
            deconv_ctx,
            ref_model_epochs=params.cell2location_ref_model_epochs,
            n_epochs=params.cell2location_n_epochs,
            n_cells_per_spot=params.cell2location_n_cells_per_spot or 30,
            detection_alpha=params.cell2location_detection_alpha,
            use_gpu=params.use_gpu,
            batch_key=params.cell2location_batch_key,
            categorical_covariate_keys=params.cell2location_categorical_covariate_keys,
            apply_gene_filtering=params.cell2location_apply_gene_filtering,
            gene_filter_cell_count_cutoff=params.cell2location_gene_filter_cell_count_cutoff,
            gene_filter_cell_percentage_cutoff2=params.cell2location_gene_filter_cell_percentage_cutoff2,
            gene_filter_nonz_mean_cutoff=params.cell2location_gene_filter_nonz_mean_cutoff,
            ref_model_lr=params.cell2location_ref_model_lr,
            cell2location_lr=params.cell2location_lr,
            ref_model_train_size=params.cell2location_ref_model_train_size,
            cell2location_train_size=params.cell2location_train_size,
            early_stopping=params.cell2location_early_stopping,
            early_stopping_patience=params.cell2location_early_stopping_patience,
            early_stopping_threshold=params.cell2location_early_stopping_threshold,
            use_aggressive_training=params.cell2location_use_aggressive_training,
            validation_size=params.cell2location_validation_size,
        )

    elif method == "destvi":
        from . import destvi

        return await destvi.deconvolve(
            deconv_ctx,
            n_epochs=params.destvi_n_epochs,
            n_hidden=params.destvi_n_hidden,
            n_latent=params.destvi_n_latent,
            n_layers=params.destvi_n_layers,
            dropout_rate=params.destvi_dropout_rate,
            learning_rate=params.destvi_learning_rate,
            train_size=params.destvi_train_size,
            vamp_prior_p=params.destvi_vamp_prior_p,
            l1_reg=params.destvi_l1_reg,
            use_gpu=params.use_gpu,
        )

    elif method == "stereoscope":
        from . import stereoscope

        return await stereoscope.deconvolve(
            deconv_ctx,
            n_epochs=params.stereoscope_n_epochs,
            learning_rate=params.stereoscope_learning_rate,
            batch_size=params.stereoscope_batch_size,
            use_gpu=params.use_gpu,
        )

    elif method == "rctd":
        from . import rctd

        return await rctd.deconvolve(
            deconv_ctx,
            mode=params.rctd_mode,
            max_cores=params.max_cores,
            confidence_threshold=params.rctd_confidence_threshold,
            doublet_threshold=params.rctd_doublet_threshold,
            max_multi_types=params.rctd_max_multi_types,
        )

    elif method == "spotlight":
        from . import spotlight

        return await spotlight.deconvolve(
            deconv_ctx,
            n_top_genes=params.spotlight_n_top_genes,
            nmf_model=params.spotlight_nmf_model,
            min_prop=params.spotlight_min_prop,
            scale=params.spotlight_scale,
            weight_id=params.spotlight_weight_id,
        )

    elif method == "card":
        from . import card

        return await card.deconvolve(
            deconv_ctx,
            sample_key=params.card_sample_key,
            minCountGene=params.card_minCountGene,
            minCountSpot=params.card_minCountSpot,
            imputation=params.card_imputation,
            NumGrids=params.card_NumGrids,
            ineibor=params.card_ineibor,
        )

    elif method == "tangram":
        from . import tangram

        return await tangram.deconvolve(
            deconv_ctx,
            n_epochs=params.tangram_n_epochs,
            mode=params.tangram_mode,
            learning_rate=params.tangram_learning_rate,
            density_prior=params.tangram_density_prior,
            use_gpu=params.use_gpu,
        )

    else:
        raise ParameterError(
            f"Unsupported method: {params.method}. "
            f"Supported: flashdeconv, cell2location, destvi, stereoscope, "
            f"rctd, spotlight, card, tangram"
        )


async def _store_results(
    spatial_adata: "ad.AnnData",
    proportions: pd.DataFrame,
    stats: Dict[str, Any],
    method: str,
    data_id: str,
    ctx: "ToolContext",
) -> DeconvolutionResult:
    """Store deconvolution results in AnnData and return result object."""
    proportions_key = f"deconvolution_{method}"
    cell_types = list(proportions.columns)

    # Align proportions with spatial_adata.obs_names using vectorized reindex
    # Missing spots (filtered by some methods) will be filled with 0
    # This is O(n) vs the previous O(n²) approach using get_loc() in a loop
    full_proportions = proportions.reindex(spatial_adata.obs_names).fillna(0).values

    # Store in obsm
    spatial_adata.obsm[proportions_key] = full_proportions

    # Store cell type names
    spatial_adata.uns[f"{proportions_key}_cell_types"] = cell_types

    # Add individual cell type columns to obs
    for i, ct in enumerate(cell_types):
        spatial_adata.obs[f"{proportions_key}_{ct}"] = full_proportions[:, i]

    # Add dominant cell type annotation using vectorized numpy operations
    dominant_key = f"dominant_celltype_{method}"
    cell_types_array = np.array(cell_types)
    dominant_types = cell_types_array[np.argmax(full_proportions, axis=1)]
    spatial_adata.obs[dominant_key] = pd.Categorical(dominant_types)

    # Save updated data
    await ctx.set_adata(data_id, spatial_adata)

    return DeconvolutionResult(
        data_id=data_id,
        method=method,
        dominant_type_key=dominant_key,
        n_cell_types=len(cell_types),
        cell_types=cell_types,
        proportions_key=proportions_key,
        statistics=stats,
    )
