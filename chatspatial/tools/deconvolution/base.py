"""
Base classes and utilities for deconvolution methods.

This module provides the shared infrastructure for all deconvolution methods,
following the DRY principle by extracting common validation and data preparation logic.
"""

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

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
from ...utils.exceptions import DataError, ProcessingError

# =============================================================================
# Data Preparation Helpers
# =============================================================================


async def prepare_anndata_for_counts(
    adata: ad.AnnData,
    label: str,
    ctx: "ToolContext",
    require_int_dtype: bool = False,
    require_integer_counts: bool = True,
) -> ad.AnnData:
    """Prepare AnnData for deconvolution by ensuring raw counts are available.

    This function:
    1. Gets raw count data from adata.raw, layers['counts'], or X
    2. Optionally converts to integer dtype (required for R-based methods)
    3. Returns a copy to avoid modifying the original

    Args:
        adata: Input AnnData object
        label: Label for logging (e.g., "Spatial", "Reference")
        ctx: ToolContext for warnings
        require_int_dtype: If True, convert to int32 (required for RCTD, SPOTlight, CARD)
        require_integer_counts: If True, raise error on normalized data

    Returns:
        AnnData with raw counts in X
    """
    # Get raw data source - returns RawDataResult object
    result = get_raw_data_source(
        adata,
        prefer_complete_genes=True,
        require_integer_counts=require_integer_counts,
        sample_size=100,
    )

    # Warn if using normalized data (when allowed)
    if not result.is_integer_counts and not require_integer_counts and ctx:
        await ctx.warning(
            f"{label}: Using normalized data (no raw counts available). "
            f"This may be acceptable for some reference datasets."
        )

    # Construct AnnData based on detected source
    if result.source == "raw":
        adata_copy = adata.raw.to_adata()
    elif result.source == "counts_layer":
        # adata.copy() is a deep copy, so layers["counts"] is already independent
        # No need for additional .copy() - this avoids memory overhead
        adata_copy = adata.copy()
        adata_copy.X = adata_copy.layers["counts"]
    else:
        # source == "X": use current X matrix
        adata_copy = adata.copy()

    # Convert to int32 if required (R-based methods need integer counts)
    # This also densifies sparse matrices for R compatibility
    if require_int_dtype and result.is_integer_counts:
        dense = to_dense(adata_copy.X)
        if dense.dtype != np.int32:
            adata_copy.X = dense.astype(np.int32, copy=False)
        else:
            adata_copy.X = dense

    return adata_copy


def create_deconvolution_stats(
    proportions: pd.DataFrame,
    common_genes: List[str],
    method: str,
    device: str = "CPU",
    **method_specific_params,
) -> Dict[str, Any]:
    """Create standardized statistics dictionary for deconvolution results.

    Args:
        proportions: DataFrame with cell type proportions (spots x cell_types)
        common_genes: List of common genes used for deconvolution
        method: Deconvolution method name
        device: Device used for computation (CPU/GPU)
        **method_specific_params: Additional method-specific parameters

    Returns:
        Dictionary with standardized statistics
    """
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

    # Add method-specific parameters
    stats.update(method_specific_params)

    return stats


# =============================================================================
# Deconvolution Context - Encapsulates Common Workflow
# =============================================================================


@dataclass
class DeconvolutionContext:
    """Encapsulates common data preparation for all deconvolution methods.

    This class follows the "prepare once, use many" principle to avoid
    redundant validation and data preparation across deconvolution methods.

    Usage:
        ctx = DeconvolutionContext(spatial, reference, cell_type_key, tool_ctx)
        await ctx.prepare()

        # Now use prepared data in method-specific logic
        ref = ctx.reference_prepared[:, ctx.common_genes]
        sp = ctx.spatial_prepared[:, ctx.common_genes]
    """

    spatial_adata: ad.AnnData
    reference_adata: ad.AnnData
    cell_type_key: str
    ctx: "ToolContext"

    # Prepared data (populated by prepare())
    spatial_prepared: Optional[ad.AnnData] = field(default=None, init=False)
    reference_prepared: Optional[ad.AnnData] = field(default=None, init=False)
    common_genes: Optional[List[str]] = field(default=None, init=False)
    cell_types: Optional[List[str]] = field(default=None, init=False)

    _prepared: bool = field(default=False, init=False)

    async def prepare(
        self,
        require_int_dtype: bool = False,
        min_common_genes: int = 100,
    ) -> "DeconvolutionContext":
        """Execute all common validation and data preparation.

        Args:
            require_int_dtype: If True, convert counts to int32 (for R methods)
            min_common_genes: Minimum required overlapping genes

        Returns:
            self (for method chaining)

        Raises:
            ParameterError: If cell_type_key is invalid
            DataError: If data preparation fails or insufficient gene overlap
        """
        if self._prepared:
            return self

        # 1. Validate cell type key exists in reference
        validate_obs_column(self.reference_adata, self.cell_type_key, "Cell type key")

        # 2. Get cell types
        self.cell_types = list(self.reference_adata.obs[self.cell_type_key].unique())
        if len(self.cell_types) < 2:
            raise DataError(
                f"Reference data must have at least 2 cell types, "
                f"found {len(self.cell_types)}"
            )

        # 3. Prepare spatial data (restore raw counts)
        self.spatial_prepared = await prepare_anndata_for_counts(
            self.spatial_adata, "Spatial", self.ctx, require_int_dtype
        )

        # 4. Prepare reference data (restore raw counts)
        self.reference_prepared = await prepare_anndata_for_counts(
            self.reference_adata, "Reference", self.ctx, require_int_dtype
        )

        # 5. Find and validate common genes
        self.common_genes = find_common_genes(
            self.spatial_prepared.var_names,
            self.reference_prepared.var_names,
        )

        validate_gene_overlap(
            self.common_genes,
            self.spatial_prepared.n_vars,
            self.reference_prepared.n_vars,
            min_genes=min_common_genes,
            source_name="spatial",
            target_name="reference",
        )

        self._prepared = True
        return self

    def get_subset_data(self) -> Tuple[ad.AnnData, ad.AnnData]:
        """Get spatial and reference data subset to common genes.

        Returns:
            Tuple of (spatial_subset, reference_subset)
        """
        if not self._prepared:
            raise ProcessingError("Must call prepare() before get_subset_data()")

        return (
            self.spatial_prepared[:, self.common_genes].copy(),
            self.reference_prepared[:, self.common_genes].copy(),
        )


# =============================================================================
# Convergence Checking
# =============================================================================


def check_model_convergence(
    model,
    model_name: str,
    convergence_threshold: float = 0.001,
    convergence_window: int = 50,
) -> Tuple[bool, Optional[str]]:
    """Check if a scvi-tools model has converged based on ELBO history.

    Args:
        model: Trained model with .history attribute
        model_name: Name for logging
        convergence_threshold: Maximum relative change to consider converged
        convergence_window: Number of epochs to examine

    Returns:
        Tuple of (is_converged, warning_message)
    """
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

    elbo_history = np.atleast_1d(np.array(elbo_history).flatten())
    if len(elbo_history) < convergence_window:
        return True, None

    recent_elbo = elbo_history[-convergence_window:]
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
