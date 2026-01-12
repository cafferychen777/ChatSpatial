"""
Base utilities for deconvolution methods.

Design Philosophy:
- Immutable data container (frozen dataclass) for prepared data
- Single function API for the common case
- Hook pattern for method-specific preprocessing (e.g., cell2location)
"""

from dataclasses import dataclass
from typing import (
    TYPE_CHECKING,
    Any,
    Awaitable,
    Callable,
    Dict,
    List,
    Optional,
    Tuple,
)

import anndata as ad
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ...spatial_mcp_adapter import ToolContext

from ...utils.adata_utils import (
    find_common_genes,
    get_raw_data_source,
    to_dense,
    validate_gene_overlap,
    validate_obs_column,
)
from ...utils.exceptions import DataError

# Type alias for preprocess hook
PreprocessHook = Callable[
    [ad.AnnData, ad.AnnData, "ToolContext"],
    Awaitable[Tuple[ad.AnnData, ad.AnnData]],
]


# =============================================================================
# Immutable Data Container
# =============================================================================


@dataclass(frozen=True)
class PreparedDeconvolutionData:
    """Immutable container for prepared deconvolution data.

    All fields are populated by prepare_deconvolution() and cannot be modified.
    This eliminates state machine complexity and makes data flow explicit.

    Attributes:
        spatial: Spatial AnnData subset to common genes (raw counts)
        reference: Reference AnnData subset to common genes (raw counts)
        cell_type_key: Column name for cell types in reference
        cell_types: List of unique cell types
        common_genes: List of genes present in both datasets
        ctx: ToolContext for logging/warnings

    Usage:
        data = await prepare_deconvolution(spatial, ref, "cell_type", ctx)
        proportions = run_method(data.spatial, data.reference, data.cell_types)
    """

    spatial: ad.AnnData
    reference: ad.AnnData
    cell_type_key: str
    cell_types: List[str]
    common_genes: List[str]
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
    # 1. Validate cell type key
    validate_obs_column(reference_adata, cell_type_key, "Cell type key")

    # 2. Extract cell types
    cell_types = list(reference_adata.obs[cell_type_key].unique())
    if len(cell_types) < 2:
        raise DataError(
            f"Reference data must have at least 2 cell types, found {len(cell_types)}"
        )

    # 3. Restore raw counts
    spatial_prep = await _prepare_counts(
        spatial_adata, "Spatial", ctx, require_int_dtype
    )
    reference_prep = await _prepare_counts(
        reference_adata, "Reference", ctx, require_int_dtype
    )

    # 4. Optional method-specific preprocessing
    if preprocess is not None:
        spatial_prep, reference_prep = await preprocess(
            spatial_prep, reference_prep, ctx
        )

    # 5. Find common genes
    common_genes = find_common_genes(
        spatial_prep.var_names,
        reference_prep.var_names,
    )

    # 6. Validate gene overlap
    validate_gene_overlap(
        common_genes,
        spatial_prep.n_vars,
        reference_prep.n_vars,
        min_genes=min_common_genes,
        source_name="spatial",
        target_name="reference",
    )

    # 7. Return immutable result with data subset to common genes
    return PreparedDeconvolutionData(
        spatial=spatial_prep[:, common_genes].copy(),
        reference=reference_prep[:, common_genes].copy(),
        cell_type_key=cell_type_key,
        cell_types=cell_types,
        common_genes=common_genes,
        ctx=ctx,
    )


async def _prepare_counts(
    adata: ad.AnnData,
    label: str,
    ctx: "ToolContext",
    require_int_dtype: bool,
) -> ad.AnnData:
    """Prepare AnnData by restoring raw counts."""
    result = get_raw_data_source(
        adata,
        prefer_complete_genes=True,
        require_integer_counts=True,
        sample_size=100,
    )

    if not result.is_integer_counts:
        await ctx.warning(
            f"{label}: Using normalized data (no raw counts available). "
            f"This may be acceptable for some reference datasets."
        )

    # Construct copy based on source
    if result.source == "raw":
        adata_copy = adata.raw.to_adata()
    elif result.source == "counts_layer":
        adata_copy = adata.copy()
        adata_copy.X = adata_copy.layers["counts"]
    else:
        adata_copy = adata.copy()

    # Convert to int32 if required (R-based methods)
    if require_int_dtype and result.is_integer_counts:
        dense = to_dense(adata_copy.X)
        adata_copy.X = (
            dense.astype(np.int32, copy=False) if dense.dtype != np.int32 else dense
        )

    return adata_copy


# =============================================================================
# Statistics Helper
# =============================================================================


def create_deconvolution_stats(
    proportions: pd.DataFrame,
    common_genes: List[str],
    method: str,
    device: str = "CPU",
    **method_specific_params,
) -> Dict[str, Any]:
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
        "dominant_types": proportions.idxmax(axis=1).value_counts().to_dict(),
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
) -> Tuple[bool, Optional[str]]:
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
