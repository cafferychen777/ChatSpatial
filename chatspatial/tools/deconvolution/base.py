"""
Base utilities for deconvolution methods.

Design Philosophy:
- Fixed field bindings around isolated, method-local AnnData workspaces
- Method registry with declarative configuration (MethodConfig)
- Single function API for the common case
- Hook pattern for method-specific preprocessing (e.g., cell2location)
"""

from collections.abc import Awaitable, Callable
from dataclasses import dataclass, field
from typing import (
    TYPE_CHECKING,
    Any,
    Optional,
)

import anndata as ad
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import sparse

if TYPE_CHECKING:
    from ...models.data import DeconvolutionParameters
    from ...spatial_mcp_adapter import ToolContext

from ...utils.adata_utils import (
    ensure_unique_var_names_async,
    find_common_genes,
    get_spatial_key,
    validate_gene_overlap,
    validate_obs_column,
)
from ...utils.exceptions import DataError

# Type aliases
PreprocessHook = Callable[
    [ad.AnnData, ad.AnnData, "ToolContext"],
    Awaitable[tuple[ad.AnnData, ad.AnnData]],
]

DeconvolveFunc = Callable[..., tuple[pd.DataFrame, dict[str, Any]]]


# =============================================================================
# Method Configuration
# =============================================================================


@dataclass(frozen=True)
class MethodConfig:
    """Immutable configuration for a deconvolution method.

    This dataclass follows the same frozen pattern as PreparedDeconvolutionData,
    ensuring configuration cannot be modified after creation.

    The registry pattern replaces the 124-line if-elif chain with declarative
    configuration, adhering to the Open-Closed Principle: new methods can be
    added without modifying the dispatch logic.

    Attributes:
        module_name: Name of the module containing the deconvolve function
        dependencies: Required packages for this method
        is_r_based: Whether method requires R (affects data type conversion)
        requires_counts: Whether method requires integer count data
        supports_gpu: Whether method supports GPU acceleration
        param_mapping: Mapping from DeconvolutionParameters field -> function arg name
                      This is the single source of truth for parameter extraction.

    Example:
        MethodConfig(
            module_name="flashdeconv",
            dependencies=("flashdeconv",),
            param_mapping={
                "flashdeconv_sketch_dim": "sketch_dim",
                "flashdeconv_lambda_spatial": "lambda_spatial",
            },
        )
    """

    module_name: str
    dependencies: tuple[str, ...]
    is_r_based: bool = False
    requires_counts: bool = False
    supports_gpu: bool = False
    param_mapping: tuple[tuple[str, str], ...] = field(default_factory=tuple)

    def extract_kwargs(self, params: "DeconvolutionParameters") -> dict[str, Any]:
        """Extract method-specific kwargs from DeconvolutionParameters.

        Uses param_mapping to transform parameter names from DeconvolutionParameters
        to the function argument names expected by each method's deconvolve function.

        Args:
            params: The full DeconvolutionParameters object

        Returns:
            Dictionary of kwargs to pass to the deconvolve function
        """
        kwargs: dict[str, Any] = {}

        for params_field, func_arg in self.param_mapping:
            value = getattr(params, params_field, None)
            if value is not None:
                kwargs[func_arg] = value

        # Add use_gpu if method supports it
        if self.supports_gpu:
            kwargs["use_gpu"] = params.use_gpu

        return kwargs

    @property
    def requires_reference(self) -> bool:
        """All current deconvolution methods require reference data."""
        return True


# =============================================================================
# Prepared Data Container
# =============================================================================


@dataclass(frozen=True)
class PreparedDeconvolutionData:
    """Prepared deconvolution inputs with immutable field bindings.

    The frozen dataclass prevents methods from replacing fields. The contained
    AnnData objects are independent working copies and intentionally remain
    mutable because third-party backends add annotations and coerce dtypes.

    Attributes:
        spatial: Spatial AnnData subset to common genes (raw counts)
        reference: Reference AnnData subset to common genes (raw counts)
        cell_type_key: Column name for cell types in reference
        cell_types: List of unique cell types
        common_genes: List of genes present in both datasets
        spatial_coords: Spatial coordinates array (n_spots, 2) or None
        ctx: ToolContext for logging/warnings

    Usage:
        data = await prepare_deconvolution(spatial, ref, "cell_type", ctx)
        proportions = run_method(data.spatial, data.reference, data.cell_types)
    """

    spatial: ad.AnnData
    reference: ad.AnnData
    cell_type_key: str
    cell_types: list[str]
    common_genes: list[str]
    spatial_coords: Optional[NDArray[np.floating]]
    ctx: "ToolContext"

    @property
    def n_spots(self) -> int:
        """Number of spatial spots."""
        return self.spatial.n_obs

    @property
    def n_cell_types(self) -> int:
        """Number of cell types."""
        return len(self.cell_types)

    @property
    def n_genes(self) -> int:
        """Number of common genes."""
        return len(self.common_genes)


# =============================================================================
# Single Entry Point
# =============================================================================


async def prepare_deconvolution(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,
    ctx: "ToolContext",
    require_int_dtype: bool = False,
    require_counts: bool = False,
    min_common_genes: int = 100,
    preprocess: Optional[PreprocessHook] = None,
) -> PreparedDeconvolutionData:
    """Prepare data for deconvolution in a single function call.

    This is the primary API for deconvolution data preparation. It handles:
    1. Validation of cell type key
    2. Raw count restoration for both datasets
    3. Optional method-specific preprocessing (via hook)
    4. Common gene identification and validation
    5. Subsetting to common genes

    Args:
        spatial_adata: Spatial transcriptomics AnnData
        reference_adata: Single-cell reference AnnData
        cell_type_key: Column in reference.obs containing cell type labels
        ctx: ToolContext for logging
        require_int_dtype: Convert to int32 (required for R-based methods)
        require_counts: Raise DataError if integer counts unavailable
        min_common_genes: Minimum required gene overlap
        preprocess: Optional async hook for method-specific preprocessing.
                   Signature: async (spatial, reference, ctx) -> (spatial, reference)
                   Called after raw count restoration, before gene finding.

    Returns:
        PreparedDeconvolutionData with all fields populated

    Examples:
        # Standard usage (most methods)
        data = await prepare_deconvolution(spatial, ref, "cell_type", ctx)

        # With custom preprocessing (e.g., cell2location)
        async def custom_filter(sp, ref, ctx):
            sp = await apply_filtering(sp, ctx)
            ref = await apply_filtering(ref, ctx)
            return sp, ref

        data = await prepare_deconvolution(
            spatial, ref, "cell_type", ctx,
            preprocess=custom_filter
        )
    """
    # 1. Extract spatial coordinates from original data (before any processing)
    spatial_coords: Optional[NDArray[np.floating]] = None
    spatial_key = get_spatial_key(spatial_adata)
    if spatial_key:
        spatial_coords = np.asarray(spatial_adata.obsm[spatial_key], dtype=np.float64)

    # 2. Validate cell type key
    validate_obs_column(reference_adata, cell_type_key, "Cell type")

    # 3. Extract cell types
    cell_types = list(reference_adata.obs[cell_type_key].unique())
    if len(cell_types) < 2:
        raise DataError(
            f"Reference data must have at least 2 cell types, found {len(cell_types)}"
        )

    # 4. Restore raw counts
    spatial_prep = await _prepare_counts(
        spatial_adata, "Spatial", ctx, require_int_dtype,
        require_counts=require_counts,
    )
    reference_prep = await _prepare_counts(
        reference_adata, "Reference", ctx, require_int_dtype,
        require_counts=require_counts,
    )

    # Repair names only inside method-local workspaces. The managed spatial
    # and reference datasets are inputs, not scratch buffers for deconvolution.
    await ensure_unique_var_names_async(spatial_prep, ctx, "spatial data")
    await ensure_unique_var_names_async(reference_prep, ctx, "reference data")

    # 5. Optional method-specific preprocessing
    if preprocess is not None:
        spatial_prep, reference_prep = await preprocess(
            spatial_prep, reference_prep, ctx
        )

    # 6. Find common genes
    common_genes = find_common_genes(
        spatial_prep.var_names,
        reference_prep.var_names,
    )

    # 7. Validate gene overlap
    validate_gene_overlap(
        common_genes,
        spatial_prep.n_vars,
        reference_prep.n_vars,
        min_genes=min_common_genes,
        source_name="spatial",
        target_name="reference",
    )

    # 8. Return immutable result with data subset to common genes
    return PreparedDeconvolutionData(
        spatial=spatial_prep[:, common_genes].copy(),
        reference=reference_prep[:, common_genes].copy(),
        cell_type_key=cell_type_key,
        cell_types=cell_types,
        common_genes=common_genes,
        spatial_coords=spatial_coords,
        ctx=ctx,
    )


async def _prepare_counts(
    adata: ad.AnnData,
    label: str,
    ctx: "ToolContext",
    require_int_dtype: bool,
    require_counts: bool = False,
) -> ad.AnnData:
    """Prepare AnnData by restoring raw counts.

    Uses the unified data source priority (layers["counts"] > adata.X)
    to avoid silently using normalized .raw data as counts.

    Args:
        adata: Input AnnData object.
        label: Human-readable label for log messages (e.g. "Spatial").
        ctx: ToolContext for logging.
        require_int_dtype: Convert to int32 (required for R-based methods).
        require_counts: If True, raise DataError when no integer counts
            are available (required for count-based models).
    """
    from ...utils.adata_utils import check_is_integer_counts, get_raw_data_source

    # Use SSOT priority: layers["counts"] > adata.X
    # prefer_complete_genes=False because deconvolution needs alignment
    # with current adata dimensions (obsm, obs, spatial coords).
    result = get_raw_data_source(
        adata,
        prefer_complete_genes=False,
        require_integer_counts=require_counts,
    )

    adata_copy = adata.copy()
    adata_copy.X = result.X

    # Convert to int32 if required (R-based methods like RCTD)
    if require_int_dtype:
        is_int, _, _ = check_is_integer_counts(adata_copy.X)
        if is_int:
            if sparse.issparse(adata_copy.X):
                adata_copy.X = adata_copy.X.astype(np.int32, copy=False)
            else:
                adata_copy.X = np.asarray(adata_copy.X).astype(np.int32, copy=False)

    return adata_copy


# =============================================================================
# Statistics Helper
# =============================================================================


def _dominant_type_counts(proportions: pd.DataFrame) -> dict[str, int]:
    """Count dominant cell types, excluding all-zero rows.

    All-zero rows (no deconvolution signal) are counted under
    "unassigned" instead of being falsely attributed to the first
    cell type column (which is what ``idxmax``/``argmax`` returns
    for tied-at-zero rows).
    """
    row_sums = proportions.sum(axis=1)
    zero_mask = row_sums == 0
    dominant = proportions.loc[~zero_mask].idxmax(axis=1)
    if zero_mask.any():
        unassigned = pd.Series("unassigned", index=proportions.index[zero_mask])
        dominant = pd.concat([dominant, unassigned])
    return dominant.value_counts().to_dict()


def create_deconvolution_stats(
    proportions: pd.DataFrame,
    common_genes: list[str],
    method: str,
    device: str = "CPU",
    **method_specific_params,
) -> dict[str, Any]:
    """Create standardized statistics dictionary for deconvolution results."""
    cell_types = list(proportions.columns)
    stats = {
        "method": method,
        "device": device,
        "n_spots": len(proportions),
        "n_cell_types": len(cell_types),
        "cell_types": cell_types,
        "genes_used": len(common_genes),
        "common_genes": len(common_genes),
        "mean_proportions": proportions.mean().to_dict(),
        "dominant_types": _dominant_type_counts(proportions),
    }
    stats.update(method_specific_params)
    return stats


# =============================================================================
# Convergence Checking
# =============================================================================


def check_model_convergence(
    model,
    model_name: str,
    convergence_threshold: float = 0.001,
    convergence_window: int = 50,
) -> tuple[bool, Optional[str]]:
    """Check if a scvi-tools model has converged based on ELBO history."""
    if not hasattr(model, "history") or model.history is None:
        return True, None

    history = model.history
    elbo_keys = ["elbo_train", "elbo_validation", "train_loss_epoch"]
    elbo_history = None

    for key in elbo_keys:
        if key in history and len(history[key]) > 0:
            elbo_history = history[key]
            break

    if elbo_history is None or len(elbo_history) < convergence_window:
        return True, None

    elbo_arr = np.asarray(elbo_history).ravel()
    recent_elbo = elbo_arr[-convergence_window:]
    elbo_changes = np.abs(np.diff(recent_elbo))

    mean_value = np.abs(np.mean(recent_elbo))
    if mean_value > 0:
        relative_changes = elbo_changes / mean_value
        mean_relative_change = np.mean(relative_changes)

        if mean_relative_change > convergence_threshold:
            return False, (
                f"{model_name} may not have fully converged. "
                f"Mean relative ELBO change: {mean_relative_change:.4f} "
                f"(threshold: {convergence_threshold}). "
                "Consider increasing training epochs."
            )

    return True, None
