"""
Deconvolution tools for spatial transcriptomics data.

This module provides functions for deconvolving spatial transcriptomics data
to estimate cell type proportions in each spatial location.
"""

import traceback
import warnings
from typing import Any, Dict, List, Optional, Tuple

import anndata as ad
import anndata2ri  # For sparse matrix support in R deconvolution methods
import numpy as np
import pandas as pd
import scipy.sparse as sp
from mcp.server.fastmcp import Context

from ..utils.error_handling import suppress_output

# Try importing deep learning dependencies
try:
    import scvi
except ImportError:
    scvi = None

try:
    import torch
except ImportError:
    torch = None

# Import specific scvi modules if available
Stereoscope = None
Tangram = None
DestVI = None

if scvi:
    try:
        from scvi.external import SpatialStereoscope as Stereoscope
        from scvi.external import Tangram
        from scvi.model import DestVI  # noqa: F401
    except ImportError as e:
        # scvi-tools version compatibility issue
        warnings.warn(f"scvi-tools import issue: {e}. Some methods may be unavailable.")

# Import cell2location with graceful fallback
try:
    import cell2location
except ImportError:
    cell2location = None

from ..models.analysis import DeconvolutionResult  # noqa: E402
from ..models.data import DeconvolutionParameters  # noqa: E402

# No longer need local context manager utilities - using centralized version
# Note: cell2location 0.1.4+ includes official one_hot compatibility fix,
# so no monkey patch is needed anymore


# Helper functions to eliminate redundancy


def _validate_deconvolution_inputs(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,
    min_common_genes: int = 100,
) -> List[str]:
    """Validate inputs and return a list of common genes.

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        min_common_genes: Minimum number of common genes required

    Returns:
        List of common genes

    Raises:
        ValueError: If input data is invalid or insufficient common genes
    """
    if spatial_adata is None:
        raise ValueError("Spatial AnnData object cannot be None")
    if reference_adata is None:
        raise ValueError("Reference AnnData object cannot be None")
    if spatial_adata.n_obs == 0:
        raise ValueError("Spatial data contains no observations")
    if reference_adata.n_obs == 0:
        raise ValueError("Reference data contains no observations")
    if cell_type_key not in reference_adata.obs:
        raise ValueError(f"Cell type key '{cell_type_key}' not found in reference data")
    if len(reference_adata.obs[cell_type_key].unique()) < 2:
        raise ValueError(
            f"Reference data must contain at least 2 cell types, found {len(reference_adata.obs[cell_type_key].unique())}"
        )

    # Find common genes
    common_genes = list(set(spatial_adata.var_names) & set(reference_adata.var_names))
    if len(common_genes) < min_common_genes:
        raise ValueError(
            f"Only {len(common_genes)} genes in common between spatial data and reference data. "
            f"Need at least {min_common_genes}. Consider using a different reference dataset or "
            f"reducing the min_common_genes parameter."
        )

    return common_genes


def _get_device(use_gpu: bool, method: str) -> str:
    """Determine the appropriate compute device.

    Args:
        use_gpu: Whether to use GPU for training
        method: Name of the method (for logging purposes)

    Returns:
        Device string ('cuda', 'mps', or 'cpu')
    """
    try:
        import torch
    except ImportError:
        warnings.warn("PyTorch not available. Using CPU.")
        return "cpu"

    if not use_gpu:
        # Using CPU for deconvolution (GPU acceleration disabled)
        return "cpu"

    if torch.cuda.is_available():
        # Using CUDA GPU acceleration
        return "cuda"

    # MPS support - handle differently for different methods
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        # Cell2location has issues with MPS
        warnings.warn(
            f"MPS acceleration is available but disabled for {method} due to numerical instability issues "
            f"with cell2location 0.1.4. This issue persists even with PyTorch 2.10.0 (tested 2025-10-12). "
            f"Using CPU instead. See PyTorch Issue #132605 for status."
        )
        # Using CPU (MPS disabled due to compatibility issues)
        return "cpu"

    warnings.warn(
        f"GPU requested for {method} but neither CUDA nor MPS is available. Using CPU instead."
    )
    # Using CPU for deconvolution
    return "cpu"


def _prepare_anndata_for_counts(
    adata: ad.AnnData, data_name: str, context=None, require_int_dtype: bool = False
) -> ad.AnnData:
    """Ensure AnnData object has raw integer counts in .X with complete gene set

    Checks for raw counts in the following order:
    1. .raw - Returns complete raw data to maximize gene overlap (preferred for deconvolution)
    2. layers["counts"] - Uses counts layer if raw not available
    3. current .X - Only if already integer counts

    IMPORTANT: This function prioritizes returning complete gene sets from .raw
    to ensure maximum gene overlap between spatial and reference datasets during
    deconvolution. This prevents "0 genes in common" errors when datasets have
    different HVG selections after preprocessing.

    Args:
        adata: AnnData object to prepare
        data_name: Name of the data for logging purposes
        context: Optional MCP context for logging
        require_int_dtype: If True, convert float32/64 to int32 for R compatibility
                          If False (default), keep original dtype (for scvi-tools methods)
                          scvi-tools internally uses float32 regardless of input dtype

    Returns:
        AnnData object with integer counts in .X and complete gene set if available

    Raises:
        ValueError: If no raw integer counts can be found
    """
    import logging

    logger = logging.getLogger(__name__)

    # MEMORY OPTIMIZATION: Don't copy adata immediately
    # Only copy if needed (when not using .raw)
    adata_copy = None
    data_source = None

    # Step 1: Check .raw data first (prefer complete gene set for deconvolution)
    if adata.raw is not None:
        logger.info(f"{data_name}: Found .raw data, checking if it contains counts")
        try:
            # Get raw data to check if it contains counts
            raw_adata = adata.raw.to_adata()

            # Validate if raw contains counts (not normalized)
            # Sample first, then convert (memory efficient)
            sample_size = min(100, raw_adata.X.shape[0])
            sample_genes = min(100, raw_adata.X.shape[1])
            raw_X_sample = raw_adata.X[:sample_size, :sample_genes]

            if hasattr(raw_X_sample, "toarray"):
                sample_X = raw_X_sample.toarray()
            else:
                sample_X = raw_X_sample

            has_decimals = not np.allclose(sample_X, np.round(sample_X), atol=1e-6)
            has_negatives = sample_X.min() < 0

            if not has_decimals and not has_negatives:
                # Raw contains counts, use complete gene set
                logger.info(
                    f"{data_name}: .raw contains counts, using complete gene set ({raw_adata.n_vars} genes)"
                )
                if context:
                    context.info(
                        f"Using .raw counts for {data_name} ({raw_adata.n_vars} genes)"
                    )
                adata_copy = raw_adata  # Use raw directly, no copy needed
                data_source = "raw"
            else:
                # Raw is normalized/transformed, skip to layers['counts']
                logger.info(
                    f"{data_name}: .raw is normalized (has_decimals={has_decimals}, "
                    f"has_negatives={has_negatives}), trying counts layer"
                )
                if context:
                    context.info(
                        f"WARNING:.raw for {data_name} is normalized, checking counts layer"
                    )
        except Exception as e:
            logger.warning(
                f"{data_name}: Error accessing .raw data: {e}, trying counts layer"
            )
            data_source = None

    # Step 2: If not using .raw, copy adata now
    if adata_copy is None:
        # Create a copy to avoid modifying original
        adata_copy = adata.copy()

        # Check layers["counts"] if raw not available or not counts
        if "counts" in adata_copy.layers:
            logger.info(f"{data_name}: Found counts layer, using as raw counts")
            if context:
                context.info(f"Using counts layer for {data_name} data")
            adata_copy.X = adata_copy.layers["counts"].copy()
            data_source = "counts_layer"

    # Step 3: Use current X as fallback
    if data_source is None:
        logger.info(f"{data_name}: No counts layer or .raw found, checking current X")
        data_source = "current"

    # Validate data (sparse-aware)
    data_min = adata_copy.X.min()
    data_max = adata_copy.X.max()
    has_negatives = data_min < 0

    # Check for decimals using small sample (avoid converting entire matrix)
    sample_size = min(1000, adata_copy.n_obs)
    sample_genes = min(100, adata_copy.n_vars)
    if hasattr(adata_copy.X, "toarray"):
        X_sample = adata_copy.X[:sample_size, :sample_genes].toarray()
    else:
        X_sample = adata_copy.X[:sample_size, :sample_genes]
    has_decimals = not np.allclose(X_sample, np.round(X_sample), atol=1e-6)

    logger.info(
        f"{data_name} data: source={data_source}, "
        f"range=[{data_min:.2f}, {data_max:.2f}], "
        f"has_negatives={has_negatives}, has_decimals={has_decimals}"
    )

    # MEMORY OPTIMIZATION: Conditional dtype conversion based on downstream method
    # - R methods (RCTD, Spotlight, CARD) require int32 dtype for compatibility
    # - scvi-tools methods (Cell2location, DestVI, Stereoscope) work with float32
    #   (they internally convert to float32 regardless of input dtype)
    if (
        require_int_dtype
        and not has_negatives
        and not has_decimals
        and adata_copy.X.dtype in [np.float32, np.float64]
    ):
        logger.info(
            f"{data_name}: Converting {adata_copy.X.dtype} to int32 for R compatibility"
        )
        if context:
            context.info(
                f"🔄 Converting {data_name} from {adata_copy.X.dtype} to int32 for R method"
            )

        # Direct dtype conversion (works for both sparse and dense matrices)
        # No need for round() since has_decimals=False already verified
        # No need for separate .data handling since .astype() handles both
        adata_copy.X = adata_copy.X.astype(np.int32)

    # Check if data is valid integer counts
    if has_negatives or has_decimals:
        error_msg = (
            f"\n{data_name} data is not raw integer counts:\n"
            f"  • Data source attempted: {data_source}\n"
            f"  • Range: [{data_min:.2f}, {data_max:.2f}]\n"
            f"  • Has negative values: {has_negatives}\n"
            f"  • Has decimal values: {has_decimals}\n\n"
        )

        if has_negatives:
            error_msg += "  WARNING:Data appears to be z-score normalized (contains negative values)\n"
        elif has_decimals and data_max < 20:
            error_msg += "  WARNING:Data appears to be log-transformed\n"
        elif has_decimals:
            error_msg += "  WARNING:Data appears to be normalized (contains decimals)\n"

        error_msg += (
            "\nIMPORTANT: Deconvolution methods (Cell2location, DestVI, RCTD, Stereoscope) "
            "require raw integer counts and CANNOT work with preprocessed data.\n\n"
            "DO NOT use these preprocessing steps before deconvolution:\n"
            "  normalize_total (sc.pp.normalize_total)\n"
            "  log transformation (sc.pp.log1p)\n"
            "  scaling/z-score (sc.pp.scale)\n"
            "  any transformation that creates decimals or negative values\n\n"
            "Solutions:\n"
            "  1. Skip preprocessing before deconvolution:\n"
            "     • Load data → Directly run deconvolution\n"
            "     • Preprocessing can be done AFTER deconvolution if needed\n\n"
            "  2. If you must preprocess first:\n"
            "     • Save counts before preprocessing: adata.layers['counts'] = adata.X.copy()\n"
            "     • Then the deconvolution can use the saved counts\n\n"
            "  3. Use original data files:\n"
            "     • Load fresh data that hasn't been preprocessed\n"
            "     • Ensure the data contains only non-negative integers\n"
        )

        if context:
            context.error(error_msg)
        raise ValueError(error_msg)

    # Add processing info to uns for tracking
    adata_copy.uns[f"{data_name}_data_source"] = {
        "source": data_source,
        "range": [int(data_min), int(data_max)],
    }

    if context:
        context.info(
            f"{data_name} data validated: "
            f"source={data_source}, range=[{int(data_min)}, {int(data_max)}]"
        )

    return adata_copy


def _prepare_reference_for_card(
    adata: ad.AnnData, data_name: str, context=None
) -> ad.AnnData:
    """Prepare reference AnnData for CARD with relaxed validation.

    CARD can accept normalized reference data (e.g., from Smart-seq2 which provides
    TPM/FPKM instead of raw UMI counts). This function tries to get raw counts if
    available, but accepts normalized data if raw counts are not found.

    This matches real-world usage where reference datasets often come from different
    sequencing technologies that don't provide raw UMI counts (e.g., Arora et al. 2023
    used SCTransform-normalized reference data with CARD).

    Args:
        adata: AnnData object to prepare
        data_name: Name of the data for logging purposes
        context: Optional MCP context for logging

    Returns:
        AnnData object prepared for CARD (raw counts if available, normalized otherwise)
    """
    import logging

    logger = logging.getLogger(__name__)

    # Create a copy to avoid modifying original
    adata_copy = adata.copy()

    # Step 1: Try to get raw counts from .raw or layers['counts']
    data_source = None

    if adata_copy.raw is not None:
        try:
            raw_adata = adata_copy.raw.to_adata()

            # Sample first, then convert (memory efficient)
            sample_size = min(100, raw_adata.X.shape[0])
            sample_genes = min(100, raw_adata.X.shape[1])
            raw_X_sample = raw_adata.X[:sample_size, :sample_genes]

            if hasattr(raw_X_sample, "toarray"):
                sample_X = raw_X_sample.toarray()
            else:
                sample_X = raw_X_sample

            has_decimals = not np.allclose(sample_X, np.round(sample_X), atol=1e-6)
            has_negatives = sample_X.min() < 0

            if not has_decimals and not has_negatives:
                logger.info(
                    f"{data_name}: Using raw counts from .raw ({raw_adata.n_vars} genes)"
                )
                if context:
                    context.info(
                        f"Using .raw counts for {data_name} ({raw_adata.n_vars} genes)"
                    )
                adata_copy = raw_adata
                data_source = "raw"
        except Exception as e:
            logger.warning(f"{data_name}: Error accessing .raw data: {e}")

    if data_source is None and "counts" in adata_copy.layers:
        logger.info(f"{data_name}: Using counts from layers['counts']")
        if context:
            context.info(f"Using counts layer for {data_name}")
        adata_copy.X = adata_copy.layers["counts"].copy()
        data_source = "counts_layer"

    # Step 2: If no raw counts found, use current data (may be normalized)
    if data_source is None:
        # Sample first, then convert (more memory efficient)
        sample_size = min(100, adata_copy.n_obs)
        sample_genes = min(100, adata_copy.n_vars)
        X_sample = adata_copy.X[:sample_size, :sample_genes]

        if hasattr(X_sample, "toarray"):
            sample_X = X_sample.toarray()
        else:
            sample_X = X_sample

        has_decimals = not np.allclose(sample_X, np.round(sample_X), atol=1e-6)

        if has_decimals:
            logger.warning(
                f"{data_name}: No raw counts found, using current data (appears normalized). "
                f"Range: [{sample_X.min():.2f}, {sample_X.max():.2f}]"
            )
            if context:
                context.warning(
                    f"WARNING: {data_name}: Using normalized data (no raw counts available). "
                    f"This is acceptable for CARD reference data from technologies like Smart-seq2."
                )
            data_source = "normalized"
        else:
            logger.info(f"{data_name}: Using current data as counts")
            data_source = "current"

    # Log final data info (sparse-aware)
    data_min = adata_copy.X.min()
    data_max = adata_copy.X.max()

    logger.info(
        f"{data_name} prepared: source={data_source}, "
        f"range=[{data_min:.2f}, {data_max:.2f}], "
        f"n_genes={adata_copy.n_vars}"
    )

    if context:
        context.info(
            f"{data_name} prepared: source={data_source}, "
            f"range=[{data_min:.2f}, {data_max:.2f}]"
        )

    return adata_copy


def _create_deconvolution_stats(
    proportions: pd.DataFrame,
    common_genes: List[str],
    method: str,
    device: str,
    **method_specific_params,
) -> Dict[str, Any]:
    """Create standardized statistics dictionary for deconvolution results.

    Args:
        proportions: Cell type proportions DataFrame
        common_genes: List of common genes used
        method: Deconvolution method name
        device: Device used for computation
        **method_specific_params: Additional method-specific parameters

    Returns:
        Statistics dictionary
    """
    stats = {
        "cell_types": list(proportions.columns),
        "n_cell_types": len(proportions.columns),
        "mean_proportions": proportions.mean().to_dict(),
        "genes_used": len(common_genes),
        "common_genes": len(common_genes),
        "method": method,
        "device": device,
    }

    # Add method-specific parameters
    stats.update(method_specific_params)

    return stats


def _validate_and_process_proportions(
    proportions: pd.DataFrame,
    method: str,
    normalize: bool = False,
    context=None,
) -> pd.DataFrame:
    """
    Validate and process deconvolution proportions while preserving data integrity.

    DESIGN NOTE: This function is synchronous (not async) because it needs to be called
    from both sync and async deconvolution methods. The context object supports sync calls,
    so we use sync mode for maximum compatibility.

    SCIENTIFIC INTEGRITY: This function respects the real meaning of special values:
    - NaN means "computation failed" - NOT "zero cells"
    - Negative values indicate algorithm errors - NOT to be silently fixed
    - Non-normalized sums may indicate issues - NOT to be forced to 1

    Args:
        proportions: Raw deconvolution output
        method: Deconvolution method name
        normalize: Whether to normalize to sum to 1 (user must explicitly request)
        context: MCP context for logging

    Returns:
        Processed proportions with transparency

    Raises:
        ValueError: If negative values detected (algorithm error)
    """
    # 1. Check for NaN values - preserve them for transparency
    nan_mask = proportions.isna()
    if nan_mask.any().any():
        nan_count = nan_mask.sum().sum()
        nan_spots = nan_mask.any(axis=1).sum()
        total_values = proportions.shape[0] * proportions.shape[1]
        nan_percentage = (nan_count / total_values) * 100

        if context:
            context.warning(
                f"WARNING:Deconvolution produced {nan_count} NaN values ({nan_percentage:.1f}%) "
                f"in {nan_spots}/{proportions.shape[0]} spots.\n\n"
                "IMPORTANT: NaN indicates computation failure, NOT absence of cell types.\n"
                "These values are preserved for transparency.\n\n"
                "Possible causes:\n"
                "• Algorithm convergence failure\n"
                "• Insufficient gene overlap\n"
                "• Numerical instability\n\n"
                "Consider:\n"
                "1. Checking input data quality\n"
                "2. Adjusting algorithm parameters\n"
                "3. Using a different deconvolution method"
            )

    # 2. Check for negative values - this is a critical error
    if (proportions < 0).any().any():
        neg_mask = proportions < 0
        neg_count = neg_mask.sum().sum()
        neg_spots = neg_mask.any(axis=1).sum()
        min_value = proportions.min().min()

        error_msg = (
            f"CRITICAL: {method} produced {neg_count} negative values "
            f"(minimum: {min_value:.4f}) in {neg_spots} spots.\n\n"
            "This indicates a serious problem:\n"
            "• Algorithm implementation error\n"
            "• Reference-spatial data incompatibility\n"
            "• Invalid input data format\n\n"
            "Negative cell type proportions are biologically impossible.\n"
            "Cannot proceed with invalid results.\n\n"
            "SOLUTIONS:\n"
            "1. Verify input data is raw counts (not normalized)\n"
            "2. Check reference data quality\n"
            "3. Report this as a bug if using standard data"
        )

        if context:
            context.error(error_msg)

        raise ValueError(error_msg)

    # 3. Analyze sum deviation - inform but don't force normalization
    row_sums = proportions.sum(axis=1, skipna=True)

    # Different methods have different expected outputs
    if method.lower() in ["rctd", "spotlight", "tangram"]:
        # These methods should output proportions that sum to ~1
        expected_sum = 1.0
        sum_deviation = abs(row_sums - expected_sum)
        max_deviation = sum_deviation.max()

        if max_deviation > 0.1:  # More than 10% deviation
            spots_affected = (sum_deviation > 0.1).sum()

            if context:
                context.warning(
                    f"WARNING:{method} proportions deviate from expected sum of 1.0:\n"
                    f"• Maximum deviation: {max_deviation:.3f}\n"
                    f"• Spots affected: {spots_affected}/{len(row_sums)}\n"
                    f"• Sum range: [{row_sums.min():.3f}, {row_sums.max():.3f}]\n\n"
                    "This may indicate:\n"
                    "• Incomplete deconvolution\n"
                    "• Missing cell types in reference\n"
                    "• Algorithm convergence issues\n\n"
                    "Original sums are preserved."
                )

    elif method.lower() in ["cell2location", "destvi", "stereoscope"]:
        # These methods may output absolute abundances
        if context:
            context.info(
                f"{method} output statistics:\n"
                f"• Estimated cells per spot: {row_sums.mean():.2f} ± {row_sums.std():.2f}\n"
                f"• Range: [{row_sums.min():.2f}, {row_sums.max():.2f}]\n"
                f"• Zero spots: {(row_sums == 0).sum()}\n\n"
                "Note: These may represent absolute abundances, not proportions."
            )

    # 4. Handle normalization if explicitly requested
    if normalize:
        # Check for zero sums that would cause division errors
        zero_sum_spots = (row_sums == 0).sum()
        if zero_sum_spots > 0:
            if context:
                context.warning(
                    f"WARNING:Cannot normalize {zero_sum_spots} spots with zero total abundance.\n"
                    "These spots will remain as zeros."
                )

        if context:
            context.info(
                "🔄 Normalizing proportions to sum to 1 as requested.\n"
                "Note: This converts absolute abundances to relative proportions."
            )

        # Perform normalization, preserving zeros and NaN
        proportions_normalized = proportions.div(row_sums, axis=0)

        # Don't fill NaN - preserve them
        # Store original sums as metadata
        proportions_normalized.attrs = proportions_normalized.attrs or {}
        proportions_normalized.attrs["original_sums"] = row_sums
        proportions_normalized.attrs["normalization_applied"] = True

        return proportions_normalized

    # Return unmodified data with metadata
    proportions.attrs = proportions.attrs or {}
    proportions.attrs["original_sums"] = row_sums
    proportions.attrs["normalization_applied"] = False

    return proportions


def deconvolve_cell2location(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,  # REQUIRED - LLM will infer from metadata
    ref_model_epochs: int = 250,
    n_epochs: int = 30000,
    n_cells_per_spot: int = 30,
    detection_alpha: float = 200.0,
    use_gpu: bool = False,
    min_common_genes: int = 100,
    context=None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Cell2location

    Cell2location uses a two-stage training process:
    1. Reference model (NB regression) learns cell type gene expression signatures
    2. Cell2location model performs spatial mapping using these signatures

    Args:
        spatial_adata: Spatial transcriptomics AnnData object with raw counts
        reference_adata: Reference single-cell RNA-seq AnnData object with raw counts
        cell_type_key: REQUIRED - Key in reference_adata.obs for cell type information.
                      Common values: 'cell_type', 'celltype', 'annotation'
        ref_model_epochs: Number of training epochs for reference model (NB regression).
                         Default: 250 (recommended by official scvi-tools tutorial)
        n_epochs: Number of training epochs for Cell2location spatial mapping model.
                 Default: 30000 (recommended by official scvi-tools tutorial)
        n_cells_per_spot: Expected number of cells per spatial location.
                         Tissue-dependent hyperparameter. Default: 30
        detection_alpha: Controls platform/technology-specific RNA detection sensitivity.
                        Higher values = less sensitivity variation correction.
                        Default: 200 (recommended by official tutorial)
        use_gpu: Whether to use GPU acceleration for training
        min_common_genes: Minimum number of common genes required between datasets

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Raises:
        ImportError: If cell2location package is not installed
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If cell2location computation fails

    Note:
        Official parameters from scvi-tools tutorial:
        - Reference model: 250 epochs, batch_size=2500
        - Cell2location: 30000 epochs, batch_size=2500
        - N_cells_per_location: 30 (tissue-dependent)
        - detection_alpha: 200
    """
    # Import cell2location
    try:
        import cell2location as cell2location_mod  # noqa: F401
        from cell2location.models import Cell2location, RegressionModel
    except ImportError as e:
        # Provide specific installation guidance
        raise ImportError(
            "cell2location is not installed. "
            "Install with 'pip install chatspatial[full]' or 'pip install cell2location>=0.1.4'. "
            "Note: Requires PyTorch and compatible GPU drivers for optimal performance."
        ) from e

    try:
        # Validate cell type key exists
        if cell_type_key not in reference_adata.obs:
            raise ValueError(
                f"Cell type key '{cell_type_key}' not found in reference data"
            )

        # Unified device selection
        device = _get_device(use_gpu, "Cell2location")

        # Prepare data using helper functions - Cell2location expects raw count data
        # Memory optimization: helper function creates internal copy, no need to copy here
        # Cell2location (scvi-tools) works with float32, no need for int32 conversion
        ref = _prepare_anndata_for_counts(
            reference_adata, "Reference", context, require_int_dtype=False
        )
        sp = _prepare_anndata_for_counts(
            spatial_adata, "Spatial", context, require_int_dtype=False
        )

        # Find common genes AFTER data preparation (uses actual gene set that will be used)
        common_genes = list(set(ref.var_names) & set(sp.var_names))

        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Insufficient common genes after data preparation.\n"
                f"  Found: {len(common_genes)} genes\n"
                f"  Required: {min_common_genes} genes\n"
                f"  Reference data: {ref.n_vars} genes\n"
                f"  Spatial data: {sp.n_vars} genes\n\n"
                f"TIPS:\n"
                f"1. Check gene naming convention (mouse: 'Cd5l', human: 'CD5L')\n"
                f"2. Ensure both datasets are from the same species\n"
                f"3. Try using different reference dataset\n"
                f"4. Reduce min_common_genes parameter (current: {min_common_genes})"
            )

        # Subset to common genes
        ref = ref[:, common_genes]
        sp = sp[:, common_genes]

        # Ensure data is float32 for cell2location compatibility
        if ref.X.dtype != np.float32:
            ref.X = ref.X.astype(np.float32)
        if sp.X.dtype != np.float32:
            sp.X = sp.X.astype(np.float32)

        # Log progress
        # Training cell2location model

        # Check if cell type key has valid values
        if ref.obs[cell_type_key].isna().any():
            warnings.warn(
                f"Reference data contains NaN values in {cell_type_key}. These cells will be excluded."
            )
            ref = ref[~ref.obs[cell_type_key].isna()].copy()

        # Train regression model to get reference cell type signatures
        try:
            RegressionModel.setup_anndata(adata=ref, labels_key=cell_type_key)
        except Exception as e:
            raise RuntimeError(
                f"Failed to setup reference data for RegressionModel: {str(e)}"
            )

        try:
            # Create RegressionModel
            mod = RegressionModel(ref)

            # Use the suppress_output context manager
            with suppress_output():
                # Train Reference Model with official recommended parameters
                # Official tutorial: 250 epochs, batch_size=2500
                if device == "cuda":
                    mod.train(
                        max_epochs=ref_model_epochs, batch_size=2500, accelerator="gpu"
                    )
                else:
                    mod.train(max_epochs=ref_model_epochs, batch_size=2500)

        except Exception as e:
            error_msg = str(e)
            tb = traceback.format_exc()
            raise RuntimeError(f"Failed to train RegressionModel: {error_msg}\n{tb}")

        # Export reference signatures
        try:
            # Export the estimated cell abundance (summary of the posterior distribution)
            ref = mod.export_posterior(
                ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
            )

            # Extract estimated expression in each cluster
            if "means_per_cluster_mu_fg" in ref.varm.keys():
                ref_signatures = ref.varm["means_per_cluster_mu_fg"][
                    [
                        f"means_per_cluster_mu_fg_{i}"
                        for i in ref.uns["mod"]["factor_names"]
                    ]
                ].copy()
            else:
                ref_signatures = ref.var[
                    [
                        f"means_per_cluster_mu_fg_{i}"
                        for i in ref.uns["mod"]["factor_names"]
                    ]
                ].copy()
            ref_signatures.columns = ref.uns["mod"]["factor_names"]
        except Exception as e:
            raise RuntimeError(f"Failed to export reference signatures: {str(e)}")

        # Prepare spatial data for cell2location model
        try:
            Cell2location.setup_anndata(adata=sp, batch_key=None)
        except Exception as e:
            raise RuntimeError(
                f"Failed to setup spatial data for Cell2location: {str(e)}"
            )

        # Run cell2location model
        try:
            # Create Cell2location model (device specification is handled in train method)
            mod = Cell2location(
                sp,
                cell_state_df=ref_signatures,
                N_cells_per_location=n_cells_per_spot,
                detection_alpha=detection_alpha,
            )

            # Use the suppress_output context manager
            with suppress_output():
                # Train Cell2location Model with official recommended parameters
                # Official tutorial: 30000 epochs, batch_size=2500
                if device == "cuda":
                    mod.train(max_epochs=n_epochs, batch_size=2500, accelerator="gpu")
                else:
                    mod.train(max_epochs=n_epochs, batch_size=2500)

        except Exception as e:
            error_msg = str(e)
            tb = traceback.format_exc()
            raise RuntimeError(
                f"Failed to train Cell2location model: {error_msg}\n{tb}"
            )

        # Export results
        try:
            # Export the estimated cell abundance (summary of the posterior distribution)
            sp = mod.export_posterior(
                sp, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
            )
        except Exception as e:
            raise RuntimeError(f"Failed to export Cell2location results: {str(e)}")

        # Get cell abundance - try different possible keys
        cell_abundance = None
        possible_keys = [
            "q05_cell_abundance_w_sf",
            "means_cell_abundance_w_sf",
            "q50_cell_abundance_w_sf",
        ]

        for key in possible_keys:
            if key in sp.obsm:
                cell_abundance = sp.obsm[key]
                # Using cell abundance from stored key
                break

        if cell_abundance is None:
            # List available keys for debugging
            available_keys = list(sp.obsm.keys())
            raise RuntimeError(
                f"Cell2location did not produce expected output. Available keys: {available_keys}"
            )

        # Create DataFrame with results
        proportions = pd.DataFrame(
            cell_abundance, index=sp.obs_names, columns=ref_signatures.columns
        )

        # Process results transparently - Cell2location outputs absolute cell counts
        # so normalization would be converting to proportions (user should decide)
        proportions = _validate_and_process_proportions(
            proportions=proportions,
            method="cell2location",
            normalize=False,  # Cell2location outputs absolute counts, not proportions
            context=context,
        )

        # Create standardized statistics using helper function
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "Cell2location",
            device,
            n_epochs=n_epochs,
            n_cells_per_spot=n_cells_per_spot,
            use_gpu=use_gpu,
        )

        # Add model performance metrics if available
        if hasattr(mod, "history") and mod.history is not None:
            try:
                history = mod.history
                if "elbo_train" in history and len(history["elbo_train"]) > 0:
                    stats["final_elbo"] = float(history["elbo_train"][-1])
                if "elbo_validation" in history and len(history["elbo_validation"]) > 0:
                    stats["final_elbo_validation"] = float(
                        history["elbo_validation"][-1]
                    )
            except Exception as e:
                warnings.warn(f"Failed to extract model history: {str(e)}")

        return proportions, stats

    except Exception as e:
        if not isinstance(e, (ValueError, ImportError, RuntimeError)):
            error_msg = str(e)
            tb = traceback.format_exc()
            # Cell2location deconvolution failed

            # Re-raise the error since we don't have fallback
            raise
        else:
            # Re-raise the error since we don't have fallback
            raise


def is_rctd_available() -> Tuple[bool, str]:
    """Check if RCTD (spacexr) and its dependencies are available

    Returns:
        Tuple of (is_available, error_message)
    """
    try:
        # Try to import rpy2
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter

        # Test R connection with proper converter context
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                ro.r("R.version.string")
        except Exception as e:
            return False, f"R is not accessible: {str(e)}"

        # Try to install/load spacexr package
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                # Check if spacexr is installed
                ro.r("library(spacexr)")
        except Exception as e:
            return (
                False,
                f"spacexr R package is not installed or failed to load: {str(e)}. Install with: install.packages('devtools'); devtools::install_github('dmcable/spacexr', build_vignettes = FALSE)",
            )

        return True, ""

    except ImportError:
        return False, "rpy2 package is not installed. Install with 'pip install rpy2'"


def _detect_gene_format(var_names: pd.Index) -> str:
    """
    Detect gene naming format from var_names.

    Scientific basis: Human/mouse genes have standardized naming conventions:
    - Ensembl IDs: ENSG/ENSMUSG followed by 11 digits
    - Gene symbols: Capitalized names (human) or Title case (mouse)

    Args:
        var_names: Gene names from AnnData.var_names

    Returns:
        'ensembl', 'symbols', or 'mixed'
    """
    import re

    sample_size = min(100, len(var_names))
    sample = var_names[:sample_size]

    ensembl_pattern = re.compile(r"^ENS(G|MUSG)\d{11}$")
    ensembl_count = sum(1 for g in sample if ensembl_pattern.match(str(g)))

    ensembl_ratio = ensembl_count / sample_size

    if ensembl_ratio > 0.9:
        return "ensembl"
    elif ensembl_ratio < 0.1:
        return "symbols"
    else:
        return "mixed"


def _validate_gene_format_compatibility(
    spatial_adata: ad.AnnData, reference_adata: ad.AnnData, context=None
) -> Tuple[bool, str]:
    """
    Validate that spatial and reference data have compatible gene naming.

    Scientific rationale:
    - Gene expression quantification relies on gene ID matching
    - Mismatched formats lead to poor overlap and data loss
    - Statistical power is compromised with <50% gene overlap

    Args:
        spatial_adata: Spatial data
        reference_adata: Reference data
        context: Logging context

    Returns:
        (is_valid, error_message)
    """
    spatial_format = _detect_gene_format(spatial_adata.var_names)
    reference_format = _detect_gene_format(reference_adata.var_names)

    if context:
        context.info("Gene format detection:")
        context.info(
            f"   Spatial: {spatial_format} ({len(spatial_adata.var_names)} genes)"
        )
        context.info(
            f"   Reference: {reference_format} ({len(reference_adata.var_names)} genes)"
        )

    # Check overlap
    common_genes = set(spatial_adata.var_names) & set(reference_adata.var_names)
    overlap_pct = len(common_genes) / len(reference_adata.var_names)

    if overlap_pct < 0.3:  # Statistical threshold: <30% overlap is problematic
        return False, (
            f"Insufficient gene overlap: {overlap_pct:.1%}\n"
            f"   This indicates gene naming format mismatch.\n\n"
            f"   Detected formats:\n"
            f"   - Spatial: {spatial_format}\n"
            f"   - Reference: {reference_format}\n\n"
            f"   SOLUTION:\n"
            f"   If reference uses Ensembl IDs, reload spatial data with:\n"
            f"   sc.read_10x_mtx(path, var_names='gene_ids')\n\n"
            f"   Reference genes (first 5): {list(reference_adata.var_names[:5])}\n"
            f"   Spatial genes (first 5): {list(spatial_adata.var_names[:5])}"
        )

    if reference_format == "mixed" and spatial_format != "mixed":
        if context:
            context.warning(
                "WARNING:Reference has mixed gene naming (symbols + Ensembl IDs).\n"
                "   Using order-preserving gene matching to avoid expression mismatch."
            )

    return True, ""


def _calculate_quality_metrics(adata: ad.AnnData) -> dict:
    """
    Calculate lightweight quality metrics for validation.

    MEMORY OPTIMIZATION: Instead of copying entire AnnData (~4GB),
    compute only required statistics (~400KB for 50k cells).

    Args:
        adata: AnnData object to compute metrics from

    Returns:
        dict with:
        - median_numi: Median UMI counts per cell (float, 8 bytes)
        - n_vars: Number of genes (int, 8 bytes)
        - numi_per_cell: UMI counts per cell (1D array, ~8 bytes * n_obs)

    Memory usage: ~8 * n_obs bytes (e.g., 400KB for 50k cells)
    vs. full copy: ~4GB for typical reference data
    Savings: 99.99% (10000x reduction)
    """
    import numpy as np

    # Calculate total UMI per cell
    if hasattr(adata.X, "toarray"):
        # Sparse matrix
        numi_per_cell = np.array(adata.X.sum(axis=1)).flatten()
    else:
        # Dense matrix
        numi_per_cell = adata.X.sum(axis=1)

    return {
        "median_numi": np.median(numi_per_cell),
        "n_vars": adata.n_vars,
        "numi_per_cell": numi_per_cell,
    }


def _validate_subset_quality(
    before_metrics: dict, subset_adata: ad.AnnData, data_label: str, context=None
) -> None:
    """
    Validate that gene subsetting preserved data quality.

    MEMORY OPTIMIZATION: Uses pre-computed metrics instead of full AnnData copy.
    Saves ~4GB memory per validation (99.99% reduction).

    Statistical rationale:
    - nUMI (UMI counts) is a key quality metric
    - >50% nUMI loss suggests gene name mismatch
    - This is evidence of systematic error, not biological variation

    Args:
        before_metrics: Pre-computed quality metrics from _calculate_quality_metrics()
                       (dict with 'median_numi', 'n_vars', 'numi_per_cell')
        subset_adata: Subsetted AnnData
        data_label: Label for logging (e.g., "Reference")
        context: Logging context

    Raises:
        ValueError: If data quality loss exceeds statistical threshold
    """
    import numpy as np

    # Use pre-computed metrics (no need to recalculate)
    original_median = before_metrics["median_numi"]
    original_n_vars = before_metrics["n_vars"]

    # Calculate subset statistics
    if hasattr(subset_adata.X, "toarray"):
        subset_numi = np.array(subset_adata.X.sum(axis=1)).flatten()
    else:
        subset_numi = subset_adata.X.sum(axis=1)

    # Statistical metrics
    subset_median = np.median(subset_numi)
    loss_ratio = 1 - (subset_median / original_median) if original_median > 0 else 0

    # Count cells with very low UMI (< 5 is RCTD's threshold)
    low_umi_count = (subset_numi < 5).sum()
    low_umi_pct = low_umi_count / len(subset_numi)

    if context:
        context.info(f"📈 {data_label} subset quality check:")
        context.info(f"   Genes: {original_n_vars} → {subset_adata.n_vars}")
        context.info(f"   Median nUMI: {original_median:.0f} → {subset_median:.0f}")
        context.info(f"   nUMI loss: {loss_ratio:.1%}")
        context.info(f"   Cells with nUMI < 5: {low_umi_count} ({low_umi_pct:.1%})")

    # Statistical threshold: >50% median nUMI loss is pathological
    if loss_ratio > 0.5:
        raise ValueError(
            f"{data_label} data quality compromised after gene subsetting!\n\n"
            f"   Statistical evidence of gene name mismatch:\n"
            f"   - Median nUMI loss: {loss_ratio:.1%} (threshold: 50%)\n"
            f"   - Before subsetting: {original_median:.0f} UMI/cell\n"
            f"   - After subsetting: {subset_median:.0f} UMI/cell\n"
            f"   - Cells with nUMI < 5: {low_umi_count}/{len(subset_numi)} ({low_umi_pct:.1%})\n\n"
            f"   This indicates that gene names in var_names don't match\n"
            f"   the actual genes in the expression matrix.\n\n"
            f"   SOLUTION:\n"
            f"   Reload data with matching gene ID format (Ensembl IDs recommended):\n"
            f"   sc.read_10x_mtx(path, var_names='gene_ids')"
        )

    # Warning threshold: >20% loss or >10% low-UMI cells
    if loss_ratio > 0.2 or low_umi_pct > 0.1:
        if context:
            context.warning(
                f"WARNING:Moderate {data_label} data quality loss detected.\n"
                f"   This may impact deconvolution accuracy."
            )


def deconvolve_rctd(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,  # REQUIRED - LLM will infer from metadata
    mode: str = "full",
    max_cores: int = 4,
    confidence_threshold: float = 10.0,
    doublet_threshold: float = 25.0,
    min_common_genes: int = 100,
    context=None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using RCTD (Robust Cell Type Decomposition)

    IMPORTANT: Gene Naming Convention
    ----------------------------------
    Reference and spatial data MUST use the SAME gene identifier format:
    - If reference uses Ensembl IDs (ENSG...), spatial data must too
    - If reference uses gene symbols (Cd5l, CD5L), spatial data must too

    For 10x Visium data, specify gene format when loading:
        # Use Ensembl IDs (recommended - includes all genes)
        adata = sc.read_10x_mtx(path, var_names='gene_ids')

        # Use gene symbols (default - protein-coding genes only)
        adata = sc.read_10x_mtx(path, var_names='gene_symbols')

    Scientific Rationale:
    - Gene expression quantification relies on accurate gene ID matching
    - Mismatched formats cause systematic data loss (>50% nUMI loss)
    - Statistical power is compromised with poor gene overlap (<30%)

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        mode: RCTD mode - 'full', 'doublet', or 'multi'
        max_cores: Maximum number of CPU cores to use (default: 4)
        confidence_threshold: Confidence threshold for cell type assignment (default: 10.0)
        doublet_threshold: Threshold for doublet detection (default: 25.0)
        min_common_genes: Minimum number of common genes required (default: 100)
        context: Logging context for transparent progress tracking

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Raises:
        ImportError: If rpy2 or spacexr package is not available
        ValueError: If gene format mismatch or insufficient data quality
        RuntimeError: If RCTD computation fails

    References:
        Cable et al. (2022) Nat. Biotechnol. "Robust decomposition of cell type
        mixtures in spatial transcriptomics"
    """
    # Validate cell type key exists
    if cell_type_key not in reference_adata.obs:
        raise ValueError(f"Cell type key '{cell_type_key}' not found in reference data")

    # Check if RCTD is available
    is_available, error_message = is_rctd_available()
    if not is_available:
        raise ImportError(f"RCTD is not available: {error_message}")

    # Import rpy2 modules
    import rpy2.robjects as ro
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    try:
        # Load required R packages with proper converter context
        with localconverter(ro.default_converter + pandas2ri.converter):
            rpackages.importr("spacexr")
            rpackages.importr("base")

        # Prepare data first (restore raw counts)
        # Memory optimization: helper function creates internal copy, no need to copy here
        # RCTD (R method) requires int32 dtype for R compatibility
        spatial_data = _prepare_anndata_for_counts(
            spatial_adata, "Spatial", context, require_int_dtype=True
        )
        reference_data = _prepare_anndata_for_counts(
            reference_adata, "Reference", context, require_int_dtype=True
        )

        # === SCIENTIFIC VALIDATION: Gene Format Compatibility ===
        # Validate gene naming format compatibility BEFORE subsetting
        # This prevents systematic data loss from gene name mismatch
        is_valid, error_msg = _validate_gene_format_compatibility(
            spatial_data, reference_data, context
        )
        if not is_valid:
            raise ValueError(error_msg)

        # === SCIENTIFICALLY RIGOROUS GENE MATCHING ===
        # Use order-preserving matching to avoid expression matrix misalignment
        # Rationale: For datasets with mixed gene naming (symbols + Ensembl IDs),
        # set intersection may match symbols while expression is in Ensembl genes
        # This preserves the reference gene order to ensure var_names match expression data
        spatial_genes_set = set(spatial_data.var_names)
        common_genes = [g for g in reference_data.var_names if g in spatial_genes_set]

        if context:
            context.info(
                f"DEBUG:Gene matching: {len(common_genes)} common genes "
                f"({len(common_genes)/len(reference_data.var_names)*100:.1f}% of reference)"
            )

        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Insufficient common genes after data preparation.\n"
                f"  Found: {len(common_genes)} genes\n"
                f"  Required: {min_common_genes} genes\n"
                f"  Reference data: {reference_data.n_vars} genes\n"
                f"  Spatial data: {spatial_data.n_vars} genes\n\n"
                f"TIPS:\n"
                f"1. Check gene naming convention (mouse: 'Cd5l', human: 'CD5L')\n"
                f"2. Ensure both datasets are from the same species\n"
                f"3. Try using different reference dataset\n"
                f"4. Reduce min_common_genes parameter (current: {min_common_genes})"
            )

        # MEMORY OPTIMIZATION: Calculate lightweight metrics instead of full copy
        # Before: spatial_data.copy() + reference_data.copy() = ~4GB
        # After: _calculate_quality_metrics() = ~400KB
        # Savings: 99.99% (10000x reduction)
        spatial_metrics_before = _calculate_quality_metrics(spatial_data)
        reference_metrics_before = _calculate_quality_metrics(reference_data)

        # Subset to common genes
        spatial_data = spatial_data[:, common_genes].copy()
        reference_data = reference_data[:, common_genes].copy()

        # === STATISTICAL VALIDATION: Subset Quality ===
        # Validate that gene subsetting preserved data quality
        # This catches gene name format mismatches that passed initial checks
        _validate_subset_quality(
            reference_metrics_before, reference_data, "Reference", context
        )
        _validate_subset_quality(
            spatial_metrics_before, spatial_data, "Spatial", context
        )

        if context and cell_type_key in reference_data.obs:
            ref_ct_counts = reference_data.obs[cell_type_key].value_counts()
            context.info(
                f"After subsetting: {len(common_genes)} common genes, {len(ref_ct_counts)} cell types"
            )
            min_count = ref_ct_counts.min()
            if min_count < 25:
                context.warning(
                    f"WARNING:Minimum cells per type: {min_count} (< 25 required)"
                )

        # Get spatial coordinates if available
        if "spatial" in spatial_adata.obsm:
            coords = pd.DataFrame(
                spatial_adata.obsm["spatial"],
                index=spatial_adata.obs_names,
                columns=["x", "y"],
            )
        else:
            # Create dummy coordinates
            coords = pd.DataFrame(
                {"x": range(spatial_adata.n_obs), "y": [0] * spatial_adata.n_obs},
                index=spatial_adata.obs_names,
            )

        # Prepare count matrices for R using anndata2ri
        # Note: anndata2ri handles both sparse and dense matrices automatically
        if context:
            matrix_type = "sparse" if sp.issparse(spatial_data.X) else "dense"
            context.info(
                f"Using anndata2ri for matrix transfer to RCTD ({matrix_type} data) "
                f"(spatial: {spatial_data.X.shape}, reference: {reference_data.X.shape})"
            )

        # Prepare cell type information as named factor
        cell_types = reference_data.obs[cell_type_key].copy()
        # Clean cell type names - RCTD doesn't allow special characters
        cell_types = cell_types.str.replace("/", "_", regex=False)
        cell_types = cell_types.str.replace(" ", "_", regex=False)
        cell_types_series = pd.Series(
            cell_types.values, index=reference_data.obs_names, name="cell_type"
        )

        cell_type_counts = cell_types_series.value_counts()
        if context:
            context.info(
                f"Reference: {len(cell_type_counts)} cell types, {len(cell_types_series)} total cells"
            )
            min_count = cell_type_counts.min()
            if min_count < 25:
                context.warning(
                    f"WARNING:WARNING: Minimum cells per type = {min_count} (< 25 required)"
                )
                # List types below threshold
                low_types = [ct for ct, count in cell_type_counts.items() if count < 25]
                context.warning(f"   Types below threshold: {', '.join(low_types)}")

        # Calculate nUMI for both datasets
        spatial_numi = pd.Series(
            np.array(spatial_data.X.sum(axis=1)).flatten(),
            index=spatial_data.obs_names,
            name="nUMI",
        )

        reference_numi = pd.Series(
            np.array(reference_data.X.sum(axis=1)).flatten(),
            index=reference_data.obs_names,
            name="nUMI",
        )

        # anndata2ri handles both sparse and dense matrices automatically
        with localconverter(ro.default_converter + anndata2ri.converter):
            # Transfer matrices directly (genes × spots/cells for R convention)
            ro.globalenv["spatial_counts"] = spatial_data.X.T
            ro.globalenv["reference_counts"] = reference_data.X.T

            # Set row/column names in R
            ro.globalenv["gene_names_spatial"] = ro.StrVector(spatial_data.var_names)
            ro.globalenv["spot_names"] = ro.StrVector(spatial_data.obs_names)
            ro.globalenv["gene_names_ref"] = ro.StrVector(reference_data.var_names)
            ro.globalenv["cell_names"] = ro.StrVector(reference_data.obs_names)

            ro.r(
                """
                rownames(spatial_counts) <- gene_names_spatial
                colnames(spatial_counts) <- spot_names
                rownames(reference_counts) <- gene_names_ref
                colnames(reference_counts) <- cell_names
                """
            )

        # Convert other data using pandas2ri
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_coords = ro.conversion.py2rpy(coords)
            r_spatial_numi = ro.conversion.py2rpy(spatial_numi)
            r_cell_types = ro.conversion.py2rpy(cell_types_series)
            r_reference_numi = ro.conversion.py2rpy(reference_numi)

            # Create SpatialRNA object
            ro.globalenv["coords"] = r_coords
            ro.globalenv["numi_spatial"] = r_spatial_numi

            puck = ro.r(
                """
                SpatialRNA(coords, spatial_counts, numi_spatial)
                """
            )

            # Create Reference object
            ro.globalenv["cell_types_vec"] = r_cell_types
            ro.globalenv["numi_ref"] = r_reference_numi

            reference = ro.r(
                """
                cell_types_factor <- as.factor(cell_types_vec)
                names(cell_types_factor) <- names(cell_types_vec)
                Reference(reference_counts, cell_types_factor, numi_ref, min_UMI = 5)
                """
            )

            # Create RCTD object (must be within localconverter context for async compatibility)
            ro.globalenv["puck"] = puck
            ro.globalenv["reference"] = reference
            ro.globalenv["max_cores_val"] = max_cores

            myRCTD = ro.r(
                f"""
        create.RCTD(puck, reference, max_cores = {max_cores}, UMI_min_sigma = 10)
        """
            )

        # Set RCTD parameters
        ro.globalenv["myRCTD"] = myRCTD
        ro.globalenv["rctd_mode"] = mode
        ro.globalenv["conf_thresh"] = confidence_threshold
        ro.globalenv["doub_thresh"] = doublet_threshold

        ro.r(
            """
        myRCTD@config$CONFIDENCE_THRESHOLD <- conf_thresh
        myRCTD@config$DOUBLET_THRESHOLD <- doub_thresh
        """
        )

        # Run RCTD using the unified run.RCTD function
        myRCTD = ro.r(
            """
        myRCTD <- run.RCTD(myRCTD, doublet_mode = rctd_mode)
        myRCTD
        """
        )

        # Extract results
        if mode == "full":
            # For full mode, get weights matrix and cell type names
            ro.r(
                """
            weights_matrix <- myRCTD@results$weights
            cell_type_names <- myRCTD@cell_type_info$renorm[[2]]
            spot_names <- rownames(weights_matrix)
            """
            )

            weights_r = ro.r("as.matrix(weights_matrix)")
            cell_type_names_r = ro.r("cell_type_names")
            spot_names_r = ro.r("spot_names")

            # Convert back to pandas using proper converter context
            with localconverter(
                ro.default_converter + pandas2ri.converter + numpy2ri.converter
            ):
                weights_array = ro.conversion.rpy2py(weights_r)
                cell_type_names = ro.conversion.rpy2py(cell_type_names_r)
                spot_names = ro.conversion.rpy2py(spot_names_r)

            # Create DataFrame with proper index and column names
            proportions = pd.DataFrame(
                weights_array, index=spot_names, columns=cell_type_names
            )

        else:
            # For doublet/multi mode, use official structure
            try:
                if mode == "doublet":
                    # For doublet mode, use weights_doublet and results_df (official structure)
                    ro.r(
                        """
                        # Extract official doublet mode results
                        if("weights_doublet" %in% names(myRCTD@results) && "results_df" %in% names(myRCTD@results)) {
                            # Get official doublet weights matrix and results dataframe
                            weights_doublet <- myRCTD@results$weights_doublet
                            results_df <- myRCTD@results$results_df
                            
                            # Get all cell type names
                            cell_type_names <- myRCTD@cell_type_info$renorm[[2]]
                            spot_names <- rownames(results_df)
                            n_spots <- length(spot_names)
                            n_cell_types <- length(cell_type_names)
                            
                            # Initialize full weights matrix
                            weights_matrix <- matrix(0, nrow = n_spots, ncol = n_cell_types)
                            rownames(weights_matrix) <- spot_names
                            colnames(weights_matrix) <- cell_type_names
                            
                            # Fill weights based on doublet results
                            for(i in 1:n_spots) {
                                spot_class <- results_df$spot_class[i]
                                
                                if(spot_class %in% c("doublet_certain", "doublet_uncertain")) {
                                    # For doublet spots, use first_type and second_type with doublet weights
                                    first_type <- as.character(results_df$first_type[i])
                                    second_type <- as.character(results_df$second_type[i])
                                    
                                    if(first_type %in% cell_type_names) {
                                        first_idx <- which(cell_type_names == first_type)
                                        weights_matrix[i, first_idx] <- weights_doublet[i, "first_type"]
                                    }
                                    if(second_type %in% cell_type_names && second_type != first_type) {
                                        second_idx <- which(cell_type_names == second_type)
                                        weights_matrix[i, second_idx] <- weights_doublet[i, "second_type"]
                                    }
                                } else if(spot_class == "singlet") {
                                    # For singlet spots, assign 100% to first_type
                                    first_type <- as.character(results_df$first_type[i])
                                    if(first_type %in% cell_type_names) {
                                        first_idx <- which(cell_type_names == first_type)
                                        weights_matrix[i, first_idx] <- 1.0
                                    }
                                }
                                # For "reject" spots, leave as zeros
                            }
                        } else {
                            stop("Official doublet mode structures (weights_doublet, results_df) not found")
                    }
                    """
                    )
                else:
                    # For multi mode, use the old logic (or implement proper multi mode later)
                    ro.r(
                        """
                    # Extract results for multi mode (fallback to old logic)
                    results_list <- myRCTD@results
                    spot_names <- names(results_list)
                    cell_type_names <- myRCTD@cell_type_info$renorm[[2]]
                    n_spots <- length(spot_names)
                    n_cell_types <- length(cell_type_names)

                    # Initialize weights matrix
                    weights_matrix <- matrix(0, nrow = n_spots, ncol = n_cell_types)
                    rownames(weights_matrix) <- spot_names
                    colnames(weights_matrix) <- cell_type_names

                    # Fill in weights for multi mode (implement later if needed)
                    # For now, this will fail gracefully
                    stop("Multi mode not yet properly implemented")
                    """
                    )

                weights_r = ro.r("weights_matrix")
                cell_type_names_r = ro.r("cell_type_names")
                spot_names_r = ro.r("spot_names")

                # Convert back to pandas (already in converter context)
                weights_array = ro.conversion.rpy2py(weights_r)
                cell_type_names = ro.conversion.rpy2py(cell_type_names_r)
                spot_names = ro.conversion.rpy2py(spot_names_r)

                # Create DataFrame with proper index and column names
                proportions = pd.DataFrame(
                    weights_array, index=spot_names, columns=cell_type_names
                )

            except Exception as e:
                # Failed to extract detailed RCTD results
                traceback.print_exc()
                # Re-raise the exception instead of using misleading fallback
                raise RuntimeError(
                    f"RCTD {mode} mode result extraction failed: {str(e)}"
                ) from e

        # SCIENTIFIC INTEGRITY: Process results transparently
        # RCTD should output proportions that sum to ~1, but we won't force it

        # Check for NaN values - preserve them for transparency
        nan_mask = proportions.isna()
        if nan_mask.any().any():
            nan_count = nan_mask.sum().sum()
            nan_spots = nan_mask.any(axis=1).sum()

            if context:
                warnings.warn(
                    f"RCTD produced {nan_count} NaN values in {nan_spots} spots. "
                    "NaN indicates computation failure, NOT absence of cell types. "
                    "These values are preserved for transparency."
                )

        # Check for negative values - this would be a critical error
        if (proportions < 0).any().any():
            neg_count = (proportions < 0).sum().sum()
            min_value = proportions.min().min()

            error_msg = (
                f"CRITICAL: RCTD produced {neg_count} negative values (min: {min_value:.4f}). "
                "This indicates a serious algorithm error. Negative proportions are impossible."
            )
            raise ValueError(error_msg)

        # Analyze sum deviation but don't force normalization
        row_sums = proportions.sum(axis=1, skipna=True)
        sum_deviation = abs(row_sums - 1.0)
        max_deviation = sum_deviation.max()

        if max_deviation > 0.1:  # More than 10% deviation
            spots_affected = (sum_deviation > 0.1).sum()

            if context:
                warnings.warn(
                    f"RCTD proportions deviate from expected sum of 1.0: "
                    f"max deviation {max_deviation:.3f} in {spots_affected} spots. "
                    f"This may indicate incomplete deconvolution or convergence issues."
                )

        # Don't force normalization or fill NaN with 0
        # Users can request normalization if needed
        # Original: proportions = proportions.div(row_sums, axis=0).fillna(0)

        # Store metadata about the results
        proportions.attrs = getattr(proportions, "attrs", {}) or {}
        proportions.attrs["method"] = f"RCTD-{mode}"
        proportions.attrs["original_sums"] = row_sums
        proportions.attrs["has_nan"] = nan_mask.any().any()

        # Create statistics
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            f"RCTD-{mode}",
            "CPU",  # RCTD doesn't use GPU
            mode=mode,
            max_cores=max_cores,
            confidence_threshold=confidence_threshold,
            doublet_threshold=doublet_threshold,
        )

        # RCTD completed successfully

        return proportions, stats

    except Exception as e:
        error_msg = str(e)
        tb = traceback.format_exc()
        raise RuntimeError(f"RCTD deconvolution failed: {error_msg}\n{tb}")
    # Note: No finally block needed with localconverter context managers


async def deconvolve_spatial_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: DeconvolutionParameters,  # No default - must be provided by caller (LLM)
    context: Optional[Context] = None,
) -> DeconvolutionResult:
    """Deconvolve spatial transcriptomics data to estimate cell type proportions

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing datasets
        params: Deconvolution parameters
        context: MCP context

    Returns:
        Deconvolution result with cell type proportions and visualization

    Raises:
        ValueError: If input data is invalid or required parameters are missing
        RuntimeError: If deconvolution computation fails
    """
    try:
        # Validate input parameters
        if not data_id:
            raise ValueError("Dataset ID cannot be empty")
        if not data_store:
            raise ValueError("Data store cannot be empty")

        if context:
            await context.info(
                f"Deconvolving spatial data using {params.method} method"
            )
            await context.info(f"Parameters: {params.model_dump()}")

        # Get spatial data
        if data_id not in data_store:
            raise ValueError(
                f"Dataset {data_id} not found in data store. Available datasets: {list(data_store.keys())}"
            )
        if "adata" not in data_store[data_id]:
            raise ValueError(f"Dataset {data_id} does not contain AnnData object")

        spatial_adata = data_store[data_id]["adata"]
        if spatial_adata.n_obs == 0:
            raise ValueError(f"Dataset {data_id} contains no observations")

        # Ensure spatial data has unique gene names
        if hasattr(spatial_adata, "var_names_make_unique"):
            if context:
                await context.info("Ensuring spatial dataset has unique gene names")
            spatial_adata.var_names_make_unique()

        if context:
            await context.info(f"Spatial dataset shape: {spatial_adata.shape}")

        # Load and validate reference data ONCE for methods that need it
        reference_adata = None
        if params.method in [
            "cell2location",
            "rctd",
            "destvi",
            "stereoscope",
            "tangram",
            "spotlight",
            "card",
        ]:
            if not params.reference_data_id:
                raise ValueError(
                    f"Reference data is required for method '{params.method}'. Please provide reference_data_id."
                )

            if params.reference_data_id not in data_store:
                raise ValueError(
                    f"Reference dataset '{params.reference_data_id}' not found in data store. "
                    f"Available datasets: {list(data_store.keys())}"
                )

            if "adata" not in data_store[params.reference_data_id]:
                raise ValueError(
                    f"Reference dataset {params.reference_data_id} does not contain AnnData object"
                )

            reference_adata = data_store[params.reference_data_id]["adata"]
            if reference_adata.n_obs == 0:
                raise ValueError(
                    f"Reference dataset {params.reference_data_id} contains no observations"
                )

            # Ensure reference data has unique gene names
            if hasattr(reference_adata, "var_names_make_unique"):
                if context:
                    await context.info(
                        "Ensuring reference dataset has unique gene names"
                    )
                reference_adata.var_names_make_unique()

            # Check cell type key
            if params.cell_type_key not in reference_adata.obs:
                raise ValueError(
                    f"Cell type key '{params.cell_type_key}' not found in reference data. "
                    f"Available keys: {list(reference_adata.obs.columns)}"
                )

            if context:
                await context.info(f"Reference dataset shape: {reference_adata.shape}")
                cell_types = reference_adata.obs[params.cell_type_key].unique()
                await context.info(
                    f"Using reference with {len(cell_types)} cell types: {list(cell_types)}"
                )

        # Check method-specific dependencies and provide alternatives
        method_deps = {
            "cell2location": ["cell2location", "torch"],
            "destvi": ["scvi", "torch"],  # Fixed: scvi not scvi-tools
            "stereoscope": ["scvi", "torch"],  # Fixed: scvi not scvi-tools
            "tangram": [
                "scvi",
                "torch",
                "tangram",
            ],  # Fixed: check for both scvi and tangram
            "rctd": ["rpy2"],  # R-based method
            "spotlight": ["rpy2"],  # R-based method
            "card": ["rpy2"],  # R-based method
        }

        # Get available methods
        import importlib.util

        available_methods = []
        for method, deps in method_deps.items():
            # Special handling for package names with different import names
            dep_specs = []
            for d in deps:
                if d == "scvi-tools":
                    # scvi-tools package imports as 'scvi'
                    dep_specs.append(importlib.util.find_spec("scvi"))
                else:
                    dep_specs.append(importlib.util.find_spec(d.replace("-", "_")))

            if all(spec is not None for spec in dep_specs):
                available_methods.append(method)

        if params.method not in available_methods:
            # Suggest alternatives
            if available_methods:
                alt_msg = f"Available alternatives: {', '.join(available_methods)}"
                if "cell2location" in available_methods:
                    alt_msg += " (cell2location is recommended)"
                elif "rctd" in available_methods:
                    alt_msg += " (RCTD provides similar functionality)"
            else:
                alt_msg = "No deconvolution methods available. Install dependencies with 'pip install chatspatial[advanced]'"

            error_msg = f"Method '{params.method}' is not available due to missing dependencies. {alt_msg}"
            if context:
                await context.error(error_msg)
            raise ImportError(error_msg)

        # Run deconvolution with cleaner calls and centralized error handling
        proportions, stats = None, None

        try:
            if params.method == "cell2location":
                if context:
                    await context.info("Running Cell2location deconvolution")
                proportions, stats = deconvolve_cell2location(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    n_cells_per_spot=params.n_cells_per_spot or 10,
                    detection_alpha=params.detection_alpha,
                    use_gpu=params.use_gpu,
                    context=context,
                )

            elif params.method == "rctd":
                if context:
                    await context.info("Running RCTD deconvolution")

                # Check if RCTD is available
                is_available, error_message = is_rctd_available()
                if not is_available:
                    if context:
                        await context.warning(f"RCTD is not available: {error_message}")
                    raise ImportError(f"RCTD is not available: {error_message}")

                # Set RCTD mode - default to 'full' if not specified
                rctd_mode = getattr(params, "rctd_mode", "full")
                max_cores = getattr(params, "max_cores", 4)
                confidence_threshold = getattr(
                    params, "rctd_confidence_threshold", 10.0
                )
                doublet_threshold = getattr(params, "rctd_doublet_threshold", 25.0)

                proportions, stats = deconvolve_rctd(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    mode=rctd_mode,
                    max_cores=max_cores,
                    confidence_threshold=confidence_threshold,
                    doublet_threshold=doublet_threshold,
                    context=context,
                )

            elif params.method == "destvi":
                if context:
                    await context.info("Running DestVI deconvolution")

                if scvi is None:
                    raise ImportError(
                        "scvi-tools package is not installed. Please install it with 'pip install scvi-tools'"
                    )

                proportions, stats = await deconvolve_destvi(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    n_hidden=params.destvi_n_hidden,
                    n_latent=params.destvi_n_latent,
                    n_layers=params.destvi_n_layers,
                    dropout_rate=params.destvi_dropout_rate,
                    learning_rate=params.destvi_learning_rate,
                    use_gpu=params.use_gpu,
                    context=context,
                )

            elif params.method == "stereoscope":
                if context:
                    await context.info("Running Stereoscope deconvolution")

                if Stereoscope is None:
                    raise ImportError(
                        "Stereoscope from scvi-tools package is not installed"
                    )

                proportions, stats = await deconvolve_stereoscope(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.stereoscope_n_epochs,
                    learning_rate=params.stereoscope_learning_rate,
                    batch_size=params.stereoscope_batch_size,
                    use_gpu=params.use_gpu,
                    context=context,
                )

            elif params.method == "spotlight":
                if context:
                    await context.info("Running SPOTlight deconvolution")

                # Check if SPOTlight is available
                is_available, error_message = is_spotlight_available()
                if not is_available:
                    if context:
                        await context.warning(
                            f"SPOTlight is not available: {error_message}"
                        )
                    raise ImportError(f"SPOTlight is not available: {error_message}")

                proportions, stats = deconvolve_spotlight(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_top_genes=params.n_top_genes,
                    context=context,
                )

            elif params.method == "tangram":
                if context:
                    await context.info("Running Tangram deconvolution")

                if scvi is None or Tangram is None:
                    raise ImportError(
                        "scvi-tools and Tangram are required for this method. Install with 'pip install scvi-tools'"
                    )

                proportions, stats = await deconvolve_tangram(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    n_epochs=params.n_epochs,
                    use_gpu=params.use_gpu,
                    context=context,
                )

            elif params.method == "card":
                if context:
                    await context.info("Running CARD deconvolution")

                # Check availability
                is_available, error_message = is_card_available()
                if not is_available:
                    if context:
                        await context.warning(f"CARD is not available: {error_message}")
                    raise ImportError(f"CARD is not available: {error_message}")

                # Get CARD-specific parameters
                minCountGene = getattr(params, "card_minCountGene", 100)
                minCountSpot = getattr(params, "card_minCountSpot", 5)
                sample_key = getattr(params, "card_sample_key", None)
                imputation = getattr(params, "card_imputation", False)
                NumGrids = getattr(params, "card_NumGrids", 2000)
                ineibor = getattr(params, "card_ineibor", 10)

                proportions, stats = deconvolve_card(
                    spatial_adata,
                    reference_adata,
                    cell_type_key=params.cell_type_key,
                    sample_key=sample_key,
                    minCountGene=minCountGene,
                    minCountSpot=minCountSpot,
                    imputation=imputation,
                    NumGrids=NumGrids,
                    ineibor=ineibor,
                    context=context,
                )

            else:
                raise ValueError(
                    f"Unsupported deconvolution method: {params.method}. "
                    f"Supported methods are: cell2location, rctd, destvi, stereoscope, spotlight, tangram, card"
                )

        except Exception as e:
            # Centralized error handling for deconvolution calls
            error_msg = str(e)
            tb = traceback.format_exc()
            if context:
                await context.warning(
                    f"{params.method.capitalize()} failed: {error_msg}"
                )
            raise RuntimeError(
                f"{params.method.capitalize()} deconvolution failed: {error_msg}\n{tb}"
            ) from e

        if context:
            await context.info(
                f"{params.method.capitalize()} completed successfully. Found {stats['n_cell_types']} cell types."
            )

        # Store proportions in AnnData object, handling potential shape mismatch
        proportions_key = f"deconvolution_{params.method}"

        # Create a full-size array filled with zeros for all spots
        full_proportions = np.zeros((spatial_adata.n_obs, proportions.shape[1]))

        # Map the results back to the original spot indices
        spot_mask = spatial_adata.obs_names.isin(proportions.index)
        original_spots_in_results = spatial_adata.obs_names[spot_mask]
        result_indices = [
            proportions.index.get_loc(spot) for spot in original_spots_in_results
        ]
        original_indices = np.where(spot_mask)[0]

        # Fill in the proportions for spots that have results
        full_proportions[original_indices] = proportions.iloc[result_indices].values

        spatial_adata.obsm[proportions_key] = full_proportions

        # Store cell type names in uns
        cell_types_key = f"{proportions_key}_cell_types"
        spatial_adata.uns[cell_types_key] = list(proportions.columns)

        # Also store individual cell type proportions in obs for easier visualization
        for i, cell_type in enumerate(proportions.columns):
            obs_key = f"{proportions_key}_{cell_type}"
            spatial_adata.obs[obs_key] = full_proportions[:, i]

        if context:
            await context.info(
                f"Stored cell type proportions in adata.obsm['{proportions_key}']"
            )
            await context.info(
                f"Stored cell type names in adata.uns['{cell_types_key}']"
            )
            await context.info(
                f"Stored individual cell type proportions in adata.obs with prefix '{proportions_key}_'"
            )

        # Add cell type annotation based on deconvolution results
        if context:
            await context.info(
                "Adding cell type annotation based on deconvolution results"
            )

        # Use method-prefixed output key to avoid overwriting
        # Use 'dominant_celltype' to semantically distinguish from annotation results
        dominant_type_key = f"dominant_celltype_{params.method}"

        # Determine the dominant cell type for each spot using full proportions
        dominant_cell_types = []
        for i in range(full_proportions.shape[0]):
            row = full_proportions[i]
            max_idx = row.argmax()
            dominant_cell_types.append(proportions.columns[max_idx])

        # Add to adata.obs
        spatial_adata.obs[dominant_type_key] = dominant_cell_types

        # Make it categorical
        spatial_adata.obs[dominant_type_key] = spatial_adata.obs[
            dominant_type_key
        ].astype("category")

        if context:
            await context.info(
                f"Added cell type annotation with {len(proportions.columns)} cell types"
            )
            await context.info(
                f"Most common cell type: {spatial_adata.obs[dominant_type_key].value_counts().index[0]}"
            )

        # Store scientific metadata for reproducibility
        from ..utils.metadata_storage import store_analysis_metadata

        # Extract results keys
        results_keys_dict = {
            "obs": [dominant_type_key],
            "obsm": [proportions_key],
            "uns": [cell_types_key],
        }
        # Add individual cell type proportion columns
        for cell_type in proportions.columns:
            results_keys_dict["obs"].append(f"{proportions_key}_{cell_type}")

        # Prepare parameters dict (only scientifically important ones)
        parameters_dict = {}
        if params.method == "cell2location":
            parameters_dict = {
                "n_cells_per_spot": params.n_cells_per_spot,
                "detection_alpha": params.detection_alpha,
            }
        elif params.method == "rctd":
            parameters_dict = {
                "mode": params.rctd_mode,
                "max_cores": params.max_cores,
            }
        elif params.method == "destvi":
            parameters_dict = {
                "n_latent": params.destvi_n_latent,
                "n_hidden": params.destvi_n_hidden,
                "learning_rate": params.destvi_learning_rate,
            }
        elif params.method == "stereoscope":
            parameters_dict = {
                "n_epochs": params.stereoscope_n_epochs,
                "learning_rate": params.stereoscope_learning_rate,
            }
        elif params.method == "spotlight":
            parameters_dict = {
                "hvg": params.hvg,
            }
        elif params.method == "tangram":
            parameters_dict = {
                "use_gpu": params.use_gpu,
                "n_epochs": params.n_epochs,
            }

        # Prepare reference info
        reference_info_dict = None
        if params.reference_data_id:
            reference_info_dict = {
                "reference_data_id": params.reference_data_id,
                "cell_type_key": params.cell_type_key,
            }

        # Store metadata
        store_analysis_metadata(
            spatial_adata,
            analysis_name=f"deconvolution_{params.method}",
            method=params.method,
            parameters=parameters_dict,
            results_keys=results_keys_dict,
            statistics=stats,
            reference_info=reference_info_dict,
        )

        # Return result
        result = DeconvolutionResult(
            data_id=data_id,
            method=params.method,
            dominant_type_key=dominant_type_key,
            cell_types=list(proportions.columns),
            n_cell_types=proportions.shape[1],
            proportions_key=proportions_key,
            statistics=stats,
        )

        if context:
            await context.info(
                f"Deconvolution completed successfully with {result.n_cell_types} cell types"
            )

        return result

    except (
        ValueError,
        ImportError,
        RuntimeError,
        KeyError,
        AttributeError,
        TypeError,
    ) as e:
        # Handle expected business logic errors
        if context:
            await context.warning(f"Deconvolution failed: {str(e)}")
        raise
    except Exception as e:
        # Handle truly unexpected errors while preserving system exceptions
        if isinstance(e, (KeyboardInterrupt, SystemExit, MemoryError)):
            # Don't catch system-level exceptions
            raise

        error_msg = str(e)
        tb = traceback.format_exc()
        if context:
            await context.warning(
                f"Deconvolution failed with unexpected error: {error_msg}"
            )
        raise RuntimeError(
            f"Deconvolution failed with unexpected error: {error_msg}\n{tb}"
        )


async def deconvolve_destvi(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,  # REQUIRED - LLM will infer from metadata
    n_epochs: int = 10000,
    n_hidden: int = 128,
    n_latent: int = 10,
    n_layers: int = 1,
    dropout_rate: float = 0.1,
    learning_rate: float = 1e-3,
    use_gpu: bool = False,
    context: Optional[Context] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using DestVI from scvi-tools (Official API)

    DestVI performs multi-resolution deconvolution by first training a CondSCVI
    model on reference data, then using it to initialize a DestVI model for
    spatial deconvolution.

    This implementation follows the official scvi-tools tutorial and best practices:
    https://docs.scvi-tools.org/en/stable/user_guide/models/destvi.html

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Number of epochs for training (split between CondSCVI and DestVI)
        n_hidden: Number of hidden units in neural networks
        n_latent: Dimensionality of latent space
        n_layers: Number of layers in neural networks
        dropout_rate: Dropout rate
        learning_rate: Learning rate for optimization
        use_gpu: Whether to use GPU for training
        context: MCP context for logging

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    try:
        # Check if scvi-tools is available
        if scvi is None:
            raise ImportError(
                "scvi-tools package is not installed. Please install it with 'pip install scvi-tools'"
            )

        # Validate cell type key exists
        if cell_type_key not in reference_adata.obs:
            raise ValueError(
                f"Cell type key '{cell_type_key}' not found in reference data"
            )

        # Prepare data first (DestVI requires integer counts)
        # Memory optimization: helper function creates internal copy, no need to copy here
        # DestVI (scvi-tools) works with float32, no need for int32 conversion
        ref_data = _prepare_anndata_for_counts(
            reference_adata, "reference", context, require_int_dtype=False
        )
        spatial_data = _prepare_anndata_for_counts(
            spatial_adata, "spatial", context, require_int_dtype=False
        )

        # Find common genes AFTER data preparation
        common_genes = list(set(ref_data.var_names) & set(spatial_data.var_names))
        min_common_genes = 100

        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Insufficient common genes after data preparation.\n"
                f"  Found: {len(common_genes)} genes\n"
                f"  Required: {min_common_genes} genes\n"
                f"  Reference data: {ref_data.n_vars} genes\n"
                f"  Spatial data: {spatial_data.n_vars} genes\n\n"
                f"TIPS:\n"
                f"1. Check gene naming convention (mouse: 'Cd5l', human: 'CD5L')\n"
                f"2. Ensure both datasets are from the same species\n"
                f"3. Try using different reference dataset\n"
                f"4. Reduce min_common_genes parameter (current: {min_common_genes})"
            )

        # Subset to common genes
        ref_data = ref_data[:, common_genes].copy()
        spatial_data = spatial_data[:, common_genes].copy()

        # Validate cell type information
        cell_types = ref_data.obs[cell_type_key].unique()
        if len(cell_types) < 2:
            raise ValueError(
                f"Reference data must contain at least 2 cell types, found {len(cell_types)}"
            )

        if context:
            await context.info(
                f"Training DestVI with {len(common_genes)} genes and {len(cell_types)} cell types: {list(cell_types)}"
            )

        # Calculate optimal epoch distribution (following official tutorials)
        condscvi_epochs = max(400, n_epochs // 5)  # CondSCVI needs sufficient training
        destvi_epochs = max(200, n_epochs // 10)  # DestVI typically needs fewer epochs

        # Step 1: Setup and train CondSCVI model on reference data
        if context:
            await context.info(
                f"Step 1: Training CondSCVI model ({condscvi_epochs} epochs)..."
            )

        # Setup reference data for CondSCVI with proper configuration
        scvi.model.CondSCVI.setup_anndata(
            ref_data,
            labels_key=cell_type_key,
            batch_key=None,  # Explicitly set to None for single-batch data
        )

        # Create CondSCVI model (compatible with scvi-tools 1.3.0)
        condscvi_model = scvi.model.CondSCVI(
            ref_data,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
        )

        # Train CondSCVI model (scvi-tools 1.3.x compatible syntax)
        plan_kwargs = {"lr": learning_rate}
        condscvi_model.train(
            max_epochs=condscvi_epochs,
            accelerator="gpu" if use_gpu else "cpu",  # Correct parameter name for 1.3.x
            train_size=0.9,
            plan_kwargs=plan_kwargs,
        )

        if context:
            await context.info("CondSCVI model training completed")

        # Step 2: Setup spatial data for DestVI
        if context:
            await context.info("Step 2: Setting up DestVI model...")

        # Setup spatial data (no labels needed for spatial data)
        scvi.model.DestVI.setup_anndata(spatial_data)

        # Step 3: Create DestVI model using from_rna_model (official pattern)
        destvi_model = scvi.model.DestVI.from_rna_model(
            spatial_data,
            condscvi_model,
            vamp_prior_p=15,  # Number of VampPrior components (official default)
            l1_reg=10.0,  # L1 regularization for sparsity (official default)
        )

        if context:
            await context.info("DestVI model created successfully")

        # Step 4: Train DestVI model (official training settings)
        if context:
            await context.info(
                f"Step 3: Training DestVI model ({destvi_epochs} epochs)..."
            )

        destvi_model.train(
            max_epochs=destvi_epochs,
            accelerator="gpu" if use_gpu else "cpu",  # Correct parameter name for 1.3.x
            train_size=0.9,
            plan_kwargs=plan_kwargs,
        )

        if context:
            await context.info("DestVI training completed")

        # Step 5: Get results (official API)
        if context:
            await context.info("Extracting cell type proportions...")

        # Get cell type proportions using official method
        proportions_df = destvi_model.get_proportions()
        proportions_df.index = spatial_data.obs_names  # Ensure proper indexing

        # Validate proportions
        if proportions_df.empty or proportions_df.shape[0] != spatial_data.n_obs:
            raise RuntimeError("Failed to extract valid proportions from DestVI model")

        # Process results transparently - DestVI outputs proportions
        proportions_df = _validate_and_process_proportions(
            proportions=proportions_df,
            method="DestVI",
            normalize=False,  # DestVI already outputs proportions
            context=context,
        )

        cell_types_result = list(proportions_df.columns)

        if context:
            await context.info(
                f"Generated proportions for {len(cell_types_result)} cell types: {cell_types_result}"
            )

        # Create enhanced statistics
        stats = _create_deconvolution_stats(
            proportions_df,
            common_genes,
            "DestVI",
            "gpu" if use_gpu else "cpu",
            n_epochs=n_epochs,
            condscvi_epochs=condscvi_epochs,
            destvi_epochs=destvi_epochs,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            vamp_prior_p=15,
            l1_reg=10.0,
            scvi_version="1.3.x_compatible",
        )

        return proportions_df, stats

    except Exception as e:
        error_msg = f"DestVI deconvolution failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg)


async def deconvolve_stereoscope(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,  # REQUIRED - LLM will infer from metadata
    n_epochs: int = 150000,
    learning_rate: float = 0.01,
    batch_size: int = 128,
    use_gpu: bool = False,
    context: Optional[Context] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using official Stereoscope API from scvi-tools

    This implementation follows the official two-stage Stereoscope workflow:
    1. Train RNAStereoscope model on single-cell reference data (75K epochs)
    2. Train SpatialStereoscope model on spatial data using RNA model (75K epochs)

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Total training epochs (default: 150000, matches original Stereoscope)
                  For default 150K: uses fixed 75K RNA + 75K spatial
                  For custom values: splits 50-50 between RNA and spatial models
                  Minimum recommended: 50000 for acceptable results
        learning_rate: Learning rate for optimization (default: 0.01)
                      Matches original Stereoscope default
        batch_size: Minibatch size for training (default: 128)
                   Original Stereoscope used 100
        use_gpu: Whether to use GPU for training
        context: MCP context for logging

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Note:
        Training with fewer than 50000 total epochs may result in uniform
        distribution across cell types. Use 150000 epochs (default) for best results.
    """
    try:
        # Validate cell type key exists
        if cell_type_key not in reference_adata.obs:
            raise ValueError(
                f"Cell type key '{cell_type_key}' not found in reference data"
            )

        # Prepare data first - Stereoscope requires raw counts
        # Memory optimization: helper function creates internal copy, no need to copy here
        # Stereoscope (scvi-tools) works with float32, no need for int32 conversion
        ref_data = _prepare_anndata_for_counts(
            reference_adata, "Reference", context, require_int_dtype=False
        )
        spatial_data = _prepare_anndata_for_counts(
            spatial_adata, "Spatial", context, require_int_dtype=False
        )

        # Find common genes AFTER data preparation
        common_genes = list(set(ref_data.var_names) & set(spatial_data.var_names))
        min_common_genes = 100

        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Insufficient common genes after data preparation.\n"
                f"  Found: {len(common_genes)} genes\n"
                f"  Required: {min_common_genes} genes\n"
                f"  Reference data: {ref_data.n_vars} genes\n"
                f"  Spatial data: {spatial_data.n_vars} genes\n\n"
                f"TIPS:\n"
                f"1. Check gene naming convention (mouse: 'Cd5l', human: 'CD5L')\n"
                f"2. Ensure both datasets are from the same species\n"
                f"3. Try using different reference dataset\n"
                f"4. Reduce min_common_genes parameter (current: {min_common_genes})"
            )

        # Subset to common genes
        ref_data = ref_data[:, common_genes].copy()
        spatial_data = spatial_data[:, common_genes].copy()

        # Ensure cell type key is categorical
        if not ref_data.obs[cell_type_key].dtype.name == "category":
            ref_data.obs[cell_type_key] = ref_data.obs[cell_type_key].astype("category")

        cell_types = list(ref_data.obs[cell_type_key].cat.categories)

        if context:
            await context.info(
                f"Training Stereoscope with {len(common_genes)} genes and {len(cell_types)} cell types"
            )

        # Import official Stereoscope classes
        from scvi.external import RNAStereoscope, SpatialStereoscope

        # Stage 1: Train RNAStereoscope model on reference data
        if context:
            await context.info(
                "Stage 1: Training RNAStereoscope model on reference data..."
            )

        # Setup reference data with cell type labels (no layer specified since X contains counts)
        RNAStereoscope.setup_anndata(ref_data, labels_key=cell_type_key)

        # Create and train RNA model
        rna_model = RNAStereoscope(ref_data)

        # Use fixed epoch values matching original Stereoscope (75K + 75K)
        if n_epochs == 150000:
            rna_epochs = 75000
            spatial_epochs = 75000
        else:
            # Custom n_epochs: split 50-50
            rna_epochs = n_epochs // 2
            spatial_epochs = n_epochs - rna_epochs

        # Prepare training parameters
        plan_kwargs = {"lr": learning_rate}

        if use_gpu:
            rna_model.train(
                max_epochs=rna_epochs,
                batch_size=batch_size,
                accelerator="gpu",
                plan_kwargs=plan_kwargs,
            )
        else:
            rna_model.train(
                max_epochs=rna_epochs, batch_size=batch_size, plan_kwargs=plan_kwargs
            )

        if context:
            await context.info(
                f"RNAStereoscope training completed ({rna_epochs} epochs)"
            )

        # Stage 2: Train SpatialStereoscope model using the RNA model
        if context:
            await context.info(
                "Stage 2: Training SpatialStereoscope model on spatial data..."
            )

        # Setup spatial data
        SpatialStereoscope.setup_anndata(spatial_data)

        # Create spatial model from the trained RNA model
        spatial_model = SpatialStereoscope.from_rna_model(spatial_data, rna_model)

        # Train spatial model (spatial_epochs already calculated above)
        if use_gpu:
            spatial_model.train(
                max_epochs=spatial_epochs,
                batch_size=batch_size,
                accelerator="gpu",
                plan_kwargs=plan_kwargs,
            )
        else:
            spatial_model.train(
                max_epochs=spatial_epochs,
                batch_size=batch_size,
                plan_kwargs=plan_kwargs,
            )

        if context:
            await context.info(
                f"SpatialStereoscope training completed ({spatial_epochs} epochs)"
            )

        # Extract cell type proportions
        if context:
            await context.info("Extracting cell type proportions...")

        proportions_array = spatial_model.get_proportions()

        # Create proportions DataFrame with proper indexing
        proportions = pd.DataFrame(
            proportions_array, index=spatial_data.obs_names, columns=cell_types
        )

        # Process results transparently - Stereoscope outputs proportions
        proportions = _validate_and_process_proportions(
            proportions=proportions,
            method="Stereoscope",
            normalize=False,  # Stereoscope already outputs proportions
            context=context,
        )

        # Create statistics
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "Stereoscope",
            "gpu" if use_gpu else "cpu",
            n_epochs=n_epochs,
            learning_rate=learning_rate,
            batch_size=batch_size,
        )

        # Add Stereoscope-specific statistics
        stats.update(
            {
                "rna_epochs": rna_epochs,
                "spatial_epochs": spatial_epochs,
                "workflow": "RNAStereoscope -> SpatialStereoscope",
            }
        )

        if context:
            await context.info("Stereoscope deconvolution completed successfully")

        return proportions, stats

    except Exception as e:
        error_msg = f"Stereoscope deconvolution failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e


def is_spotlight_available() -> Tuple[bool, str]:
    """Check if SPOTlight (R package) is available through rpy2

    Returns:
        Tuple of (is_available, error_message)
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter

        # Test R connection
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                ro.r("R.version.string")
        except Exception as e:
            return False, f"R is not accessible: {str(e)}"

        # Check if SPOTlight is installed in R
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                ro.r("library(SPOTlight)")
        except Exception as e:
            return (
                False,
                f"SPOTlight R package is not installed: {str(e)}. "
                "Install in R with: BiocManager::install('SPOTlight')",
            )

        return True, ""

    except ImportError:
        return (
            False,
            "rpy2 is not installed. Install with 'pip install rpy2' to use SPOTlight",
        )


def is_card_available() -> Tuple[bool, str]:
    """Check if CARD R package is available through rpy2

    Returns:
        Tuple of (is_available, error_message)
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter

        # Test R connection
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                ro.r("R.version.string")
        except Exception as e:
            return False, f"R is not accessible: {str(e)}"

        # Check if CARD is installed in R
        try:
            with localconverter(ro.default_converter + pandas2ri.converter):
                ro.r("library(CARD)")
        except Exception as e:
            return (
                False,
                f"CARD R package is not installed: {str(e)}. "
                "Install with: devtools::install_github('YingMa0107/CARD')",
            )

        return True, ""

    except ImportError:
        return (
            False,
            "rpy2 is not installed. Install with 'pip install rpy2' to use CARD",
        )


def deconvolve_spotlight(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,  # REQUIRED - LLM will infer from metadata
    n_top_genes: int = 2000,
    min_common_genes: int = 100,
    context=None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using official SPOTlight R package

    This function implements the official SPOTlight workflow using proper
    SingleCellExperiment and SpatialExperiment objects as required by the API.

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_top_genes: Number of top genes to use for deconvolution
        min_common_genes: Minimum number of common genes required

    Returns:
        Tuple of (proportions DataFrame, statistics dict)
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri, pandas2ri
        from rpy2.robjects.conversion import localconverter

        # Validate cell type key exists
        if cell_type_key not in reference_adata.obs:
            raise ValueError(
                f"Cell type key '{cell_type_key}' not found in reference data"
            )

        # Check for spatial coordinates
        if "spatial" not in spatial_adata.obsm:
            raise ValueError(
                "No spatial coordinates found in spatial_adata.obsm['spatial']. "
                "SPOTlight requires spatial coordinates for proper analysis."
            )

        # Prepare full datasets first
        # Memory optimization: helper function creates internal copy, no need to copy here
        # SPOTlight (R method) requires int32 dtype for R compatibility
        spatial_prepared = _prepare_anndata_for_counts(
            spatial_adata, "Spatial", context, require_int_dtype=True
        )
        reference_prepared = _prepare_anndata_for_counts(
            reference_adata, "Reference", context, require_int_dtype=True
        )

        # Find common genes AFTER data preparation
        common_genes = list(
            set(spatial_prepared.var_names) & set(reference_prepared.var_names)
        )

        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Insufficient common genes after data preparation.\n"
                f"  Found: {len(common_genes)} genes\n"
                f"  Required: {min_common_genes} genes\n"
                f"  Reference data: {reference_prepared.n_vars} genes\n"
                f"  Spatial data: {spatial_prepared.n_vars} genes\n\n"
                f"TIPS:\n"
                f"1. Check gene naming convention (mouse: 'Cd5l', human: 'CD5L')\n"
                f"2. Ensure both datasets are from the same species\n"
                f"3. Try using different reference dataset\n"
                f"4. Reduce min_common_genes parameter (current: {min_common_genes})"
            )

        # Subset to common genes
        spatial_subset = spatial_prepared[:, common_genes].copy()
        reference_subset = reference_prepared[:, common_genes].copy()

        # Extract count matrices for R using anndata2ri
        # Note: anndata2ri handles both sparse and dense matrices automatically
        if context:
            matrix_type = "sparse" if sp.issparse(spatial_subset.X) else "dense"
            context.info(
                f"Using anndata2ri for matrix transfer to SPOTlight ({matrix_type} data) "
                f"(spatial: {spatial_subset.X.shape}, reference: {reference_subset.X.shape})"
            )

        # Ensure integer counts (SPOTlight expects integer data)
        spatial_counts = spatial_subset.X.astype(int)
        reference_counts = reference_subset.X.astype(int)

        # Get spatial coordinates
        spatial_coords = spatial_subset.obsm["spatial"]

        # Cell type labels - clean special characters for R compatibility
        cell_types = reference_subset.obs[cell_type_key].astype(str)
        cell_types = cell_types.str.replace("/", "_", regex=False)
        cell_types = cell_types.str.replace(" ", "_", regex=False)

        # Execute SPOTlight using the official API
        # anndata2ri handles both sparse and dense matrices automatically

        # First transfer count matrices using anndata2ri
        with localconverter(ro.default_converter + anndata2ri.converter):
            # Transfer matrices (genes × spots/cells for R convention)
            ro.globalenv["spatial_counts"] = spatial_counts.T
            ro.globalenv["reference_counts"] = reference_counts.T

        # Then transfer other data with numpy2ri + pandas2ri
        with localconverter(
            ro.default_converter + pandas2ri.converter + numpy2ri.converter
        ):
            # Import required R packages
            ro.r("library(SPOTlight)")
            ro.r("library(SingleCellExperiment)")
            ro.r("library(SpatialExperiment)")
            ro.r("library(scran)")
            ro.r("library(scuttle)")

            # Transfer non-matrix data
            ro.globalenv["spatial_coords"] = spatial_coords
            ro.globalenv["gene_names"] = ro.StrVector(common_genes)
            ro.globalenv["spatial_names"] = ro.StrVector(list(spatial_subset.obs_names))
            ro.globalenv["reference_names"] = ro.StrVector(
                list(reference_subset.obs_names)
            )
            ro.globalenv["cell_types"] = ro.StrVector(cell_types.tolist())

        # Create SingleCellExperiment and SpatialExperiment objects
        # SingleCellExperiment and SpatialExperiment fully support sparse matrices
        ro.r(
            """
            # Create SCE object for reference data (sparse matrix supported)
            sce <- SingleCellExperiment(
                assays = list(counts = reference_counts),
                colData = data.frame(
                    cell_type = factor(cell_types),
                    row.names = reference_names
                )
            )
            rownames(sce) <- gene_names

            # Add logcounts
            sce <- logNormCounts(sce)

            # Create SPE object for spatial data (sparse matrix supported)
            spe <- SpatialExperiment(
                assays = list(counts = spatial_counts),
                spatialCoords = spatial_coords,
                colData = data.frame(row.names = spatial_names)
            )
            rownames(spe) <- gene_names
            colnames(spe) <- spatial_names
            
            # Detecting marker genes
            
            # Find marker genes using scran
            markers <- findMarkers(sce, groups = sce$cell_type, test.type = "wilcox")
            
            # Format marker genes for SPOTlight
            cell_type_names <- names(markers)
            mgs_list <- list()
            
            for (ct in cell_type_names) {
                ct_markers <- markers[[ct]]
                # Get top markers for each cell type
                n_markers <- min(50, nrow(ct_markers))
                top_markers <- head(ct_markers[order(ct_markers$p.value), ], n_markers)
                
                mgs_df <- data.frame(
                    gene = rownames(top_markers),
                    cluster = ct,
                    mean.AUC = -log10(top_markers$p.value + 1e-10)  # Use -log10 p-value as weight
                )
                mgs_list[[ct]] <- mgs_df
            }
            
            # Combine all marker genes
            mgs <- do.call(rbind, mgs_list)

            # Run official SPOTlight function
            spotlight_result <- SPOTlight(
                x = sce,                    # SingleCellExperiment object
                y = spe,                    # SpatialExperiment object
                groups = sce$cell_type,     # Cell type labels
                mgs = mgs,                  # Marker genes data frame
                weight_id = "mean.AUC",     # Weight column name
                group_id = "cluster",       # Group column name
                gene_id = "gene",           # Gene column name
                verbose = TRUE              # Verbose output
            )
            
            # SPOTlight deconvolution completed
            """
        )

        # Extract results
        # SPOTlight returns a list with 'mat' field containing deconvolution proportions
        with localconverter(
            ro.default_converter + pandas2ri.converter + numpy2ri.converter
        ):
            proportions_np = np.array(ro.r("spotlight_result$mat"))
            spot_names = list(ro.r("rownames(spotlight_result$mat)"))
            cell_type_names = list(ro.r("colnames(spotlight_result$mat)"))

        # Create proportions DataFrame
        proportions = pd.DataFrame(
            proportions_np, index=spot_names, columns=cell_type_names
        )

        # Process results transparently - SPOTlight should output proportions
        # but check if they already sum to 1
        proportions = _validate_and_process_proportions(
            proportions=proportions,
            method="SPOTlight",
            normalize=False,  # SPOTlight should already output proportions
            context=context,
        )

        # Add results to spatial adata
        for cell_type in proportions.columns:
            col_name = f"SPOTlight_{cell_type}"
            spatial_adata.obs[col_name] = proportions[cell_type].values

        # SPOTlight completed successfully

        # Create statistics
        stats = _create_deconvolution_stats(
            proportions, common_genes, "SPOTlight", "cpu", n_top_genes=n_top_genes
        )

        return proportions, stats

    except Exception as e:
        traceback.format_exc()
        error_msg = f"SPOTlight deconvolution failed: {str(e)}"
        # SPOTlight error occurred
        raise RuntimeError(error_msg) from e


def deconvolve_card(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,  # REQUIRED - LLM will infer from metadata
    sample_key: Optional[str] = None,  # Optional sample/batch info in reference
    minCountGene: int = 100,
    minCountSpot: int = 5,
    min_common_genes: int = 100,
    imputation: bool = False,  # Enable spatial imputation
    NumGrids: int = 2000,  # Number of grids for imputation
    ineibor: int = 10,  # Number of neighbors for imputation
    context=None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using CARD (CAncer Research Deconvolution)

    CARD is a reference-based deconvolution method that accommodates spatial
    correlation in cell type composition across tissue locations. This is a
    unique feature not present in other methods like Cell2location, RCTD, etc.

    Key features of CARD:
    - Models spatial correlation between tissue locations
    - Can create refined high-resolution spatial maps via imputation
    - Designed for cancer tissue analysis but works for all spatial data

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        sample_key: Optional key for sample/subject information in reference data
        minCountGene: Include genes with at least this many counts (default: 100)
        minCountSpot: Include genes expressed in at least this many spots (default: 5)
        min_common_genes: Minimum common genes required (default: 100)
        imputation: Whether to perform spatial imputation for higher resolution (default: False)
        NumGrids: Number of grids for CARD imputation (default: 2000)
        ineibor: Number of neighbors for CARD imputation (default: 10)
        context: MCP context for logging

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Raises:
        ImportError: If rpy2 or CARD package is not available
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If CARD computation fails
    """
    # Validate cell type key exists
    if cell_type_key not in reference_adata.obs:
        raise ValueError(f"Cell type key '{cell_type_key}' not found in reference data")

    # Check availability
    is_available, error_message = is_card_available()
    if not is_available:
        raise ImportError(f"CARD is not available: {error_message}")

    # Import rpy2
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    try:
        # Load CARD package
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.r("library(CARD)")

        # 1. Prepare data
        # IMPORTANT: CARD requires raw counts for spatial data, but can accept
        # normalized reference data (e.g., from Smart-seq2 which doesn't provide raw UMI counts).
        # This matches how Arora et al. 2023 used CARD with SCTransform-normalized reference.

        # Memory optimization: helper functions create internal copies, no need to copy here
        # CARD (R method) requires int32 dtype for R compatibility
        spatial_data = _prepare_anndata_for_counts(
            spatial_adata, "Spatial", context, require_int_dtype=True
        )

        # For reference data: try to get raw counts if available, but accept normalized if not
        reference_data = _prepare_reference_for_card(
            reference_adata, "Reference", context
        )

        # 2. Find common genes AFTER data preparation
        common_genes = list(set(spatial_data.var_names) & set(reference_data.var_names))

        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Insufficient common genes after data preparation.\n"
                f"  Found: {len(common_genes)} genes\n"
                f"  Required: {min_common_genes} genes\n"
                f"  Reference data: {reference_data.n_vars} genes\n"
                f"  Spatial data: {spatial_data.n_vars} genes\n\n"
                f"TIPS:\n"
                f"1. Check gene naming convention (mouse: 'Cd5l', human: 'CD5L')\n"
                f"2. Ensure both datasets are from the same species\n"
                f"3. Try using different reference dataset\n"
                f"4. Reduce min_common_genes parameter (current: {min_common_genes})"
            )

        # 3. Subset to common genes
        spatial_data = spatial_data[:, common_genes].copy()
        reference_data = reference_data[:, common_genes].copy()

        # 4. Get spatial coordinates
        if "spatial" in spatial_adata.obsm:
            spatial_location = pd.DataFrame(
                spatial_adata.obsm["spatial"],
                index=spatial_adata.obs_names,
                columns=["x", "y"],
            )
        else:
            # Create dummy coordinates if not available
            spatial_location = pd.DataFrame(
                {"x": range(spatial_adata.n_obs), "y": [0] * spatial_adata.n_obs},
                index=spatial_adata.obs_names,
            )

        # 5. Prepare count matrices for R using anndata2ri
        # Note: anndata2ri handles both sparse and dense matrices automatically
        if context:
            matrix_type = "sparse" if sp.issparse(spatial_data.X) else "dense"
            context.info(
                f"Using anndata2ri for matrix transfer to CARD ({matrix_type} data) "
                f"(spatial: {spatial_data.X.shape}, reference: {reference_data.X.shape})"
            )

        # 6. Prepare metadata
        sc_meta = reference_data.obs[[cell_type_key]].copy()
        sc_meta.columns = ["cellType"]  # CARD expects 'cellType' column

        # Add sample information
        if sample_key and sample_key in reference_data.obs:
            sc_meta["sampleInfo"] = reference_data.obs[sample_key]
        else:
            sc_meta["sampleInfo"] = "sample1"  # Default single sample

        # 7. Convert to R format and run CARD
        # anndata2ri handles both sparse and dense matrices automatically
        with localconverter(ro.default_converter + anndata2ri.converter):
            # Transfer matrices directly (genes × spots/cells for R convention)
            ro.globalenv["sc_count"] = reference_data.X.T
            ro.globalenv["spatial_count"] = spatial_data.X.T

            # Set row/column names in R
            ro.globalenv["gene_names_ref"] = ro.StrVector(reference_data.var_names)
            ro.globalenv["cell_names"] = ro.StrVector(reference_data.obs_names)
            ro.globalenv["gene_names_spatial"] = ro.StrVector(spatial_data.var_names)
            ro.globalenv["spot_names"] = ro.StrVector(spatial_data.obs_names)

            ro.r(
                """
                rownames(sc_count) <- gene_names_ref
                colnames(sc_count) <- cell_names
                rownames(spatial_count) <- gene_names_spatial
                colnames(spatial_count) <- spot_names
                """
            )

        # Convert metadata and spatial location using pandas2ri
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_sc_meta = ro.conversion.py2rpy(sc_meta)
            r_spatial_location = ro.conversion.py2rpy(spatial_location)

            ro.globalenv["sc_meta"] = r_sc_meta
            ro.globalenv["spatial_location"] = r_spatial_location
            ro.globalenv["minCountGene"] = minCountGene
            ro.globalenv["minCountSpot"] = minCountSpot

        # Create CARD object
        # MCP Protocol: Redirect R stdout to /dev/null to prevent non-JSON output
        # CARD prints progress messages (## QC, ## create) that break MCP JSON-RPC
        ro.r(
            """
        capture.output(
            CARD_obj <- createCARDObject(
                sc_count = sc_count,
                sc_meta = sc_meta,
                spatial_count = spatial_count,
                spatial_location = spatial_location,
                ct.varname = "cellType",
                ct.select = unique(sc_meta$cellType),
                sample.varname = "sampleInfo",
                minCountGene = minCountGene,
                minCountSpot = minCountSpot
            ),
            file = "/dev/null"
        )
        """
        )

        # Run deconvolution
        # MCP Protocol: Suppress stdout to prevent protocol pollution
        ro.r(
            """
        capture.output(
            CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj),
            file = "/dev/null"
        )
        """
        )

        # Extract results
        with localconverter(
            ro.default_converter + pandas2ri.converter + numpy2ri.converter
        ):
            # Get row and column names in R
            row_names = list(ro.r("rownames(CARD_obj@Proportion_CARD)"))
            col_names = list(ro.r("colnames(CARD_obj@Proportion_CARD)"))

            # Get proportions matrix
            proportions_r = ro.r("CARD_obj@Proportion_CARD")
            proportions_array = np.array(proportions_r)

            # Create DataFrame
            proportions = pd.DataFrame(
                proportions_array, index=row_names, columns=col_names
            )

        # Optional: Run CARD imputation for higher resolution
        imputed_proportions = None
        imputed_coordinates = None

        if imputation:
            # MCP Protocol: Suppress stdout
            ro.r(
                f"""
            capture.output(
                CARD_impute <- CARD.imputation(
                    CARD_object = CARD_obj,
                    NumGrids = {NumGrids},
                    ineibor = {ineibor}
                ),
                file = "/dev/null"
            )
            """
            )

            with localconverter(ro.default_converter + pandas2ri.converter):
                # Extract imputed proportions
                imputed_row_names = list(ro.r("rownames(CARD_impute@refined_prop)"))
                imputed_col_names = list(ro.r("colnames(CARD_impute@refined_prop)"))
                imputed_proportions_r = ro.r("CARD_impute@refined_prop")
                imputed_proportions_array = np.array(imputed_proportions_r)

                # Parse coordinates from rownames (format: "x.valuex y.value")
                # Example: "4.1x8.3" -> x=4.1, y=8.3
                coords_list = []
                for name in imputed_row_names:
                    parts = name.split("x")
                    x_val = float(parts[0])
                    y_val = float(parts[1])
                    coords_list.append([x_val, y_val])

                imputed_coords_array = np.array(coords_list)

                # Store imputed results separately (don't replace original proportions)
                imputed_proportions = pd.DataFrame(
                    imputed_proportions_array,
                    index=imputed_row_names,
                    columns=imputed_col_names,
                )

                imputed_coordinates = pd.DataFrame(
                    imputed_coords_array,
                    index=imputed_row_names,
                    columns=["x", "y"],
                )

        # 8. Process results - CARD outputs proportions that sum to 1
        proportions = _validate_and_process_proportions(
            proportions=proportions,
            method="CARD",
            normalize=False,  # CARD already outputs proportions
            context=context,
        )

        # 9. Add results to spatial adata
        # CARD may filter spots during QC, so we need to align results properly
        # Reindex to match spatial_adata.obs, filling missing spots with NaN
        proportions_aligned = proportions.reindex(
            spatial_adata.obs_names, fill_value=np.nan
        )

        for cell_type in proportions.columns:
            col_name = f"CARD_{cell_type}"
            spatial_adata.obs[col_name] = proportions_aligned[cell_type].values

        # Store imputed results if available
        if imputed_proportions is not None and imputed_coordinates is not None:
            spatial_adata.uns["card_imputation"] = {
                "proportions": imputed_proportions,
                "coordinates": imputed_coordinates,
                "n_original_spots": len(row_names),
                "n_imputed_locations": len(imputed_proportions),
                "resolution_increase": len(imputed_proportions) / len(row_names),
                "NumGrids": NumGrids,
                "ineibor": ineibor,
            }
            # CARD imputation stored in uns['card_imputation']

        # CARD completed successfully

        # 10. Create statistics
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "CARD",
            "cpu",
            minCountGene=minCountGene,
            minCountSpot=minCountSpot,
        )

        # Add imputation info to stats if available
        if imputed_proportions is not None:
            stats["imputation"] = {
                "enabled": True,
                "n_imputed_locations": len(imputed_proportions),
                "resolution_increase": f"{len(imputed_proportions) / len(row_names):.1f}x",
            }

        return proportions, stats

    except Exception as e:
        tb = traceback.format_exc()
        error_msg = f"CARD deconvolution failed: {str(e)}\n{tb}"
        # CARD error occurred
        raise RuntimeError(error_msg) from e


async def deconvolve_tangram(
    spatial_adata: ad.AnnData,
    reference_adata: ad.AnnData,
    cell_type_key: str,  # REQUIRED - LLM will infer from metadata
    n_epochs: int = 1000,
    use_gpu: bool = False,
    context: Optional[Context] = None,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using Tangram from scvi-tools

    Tangram maps single-cell RNA-seq data to spatial data, permitting
    deconvolution of cell types in spatial data like Visium.

    Args:
        spatial_adata: Spatial transcriptomics AnnData object
        reference_adata: Reference single-cell RNA-seq AnnData object
        cell_type_key: Key in reference_adata.obs for cell type information
        n_epochs: Number of epochs for training
        use_gpu: Whether to use GPU for training
        context: FastMCP context for logging

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)

    Raises:
        ImportError: If scvi-tools package is not available
        ValueError: If input data is invalid or insufficient common genes
        RuntimeError: If Tangram computation fails
    """
    try:
        if scvi is None or Tangram is None:
            raise ImportError(
                "scvi-tools and Tangram are required for this method. Install with 'pip install scvi-tools'"
            )

        # Import mudata
        try:
            import mudata as md
        except ImportError:
            raise ImportError(
                "mudata package is required for Tangram. Install with 'pip install mudata'"
            )

        # Validate cell type key exists
        if cell_type_key not in reference_adata.obs:
            raise ValueError(
                f"Cell type key '{cell_type_key}' not found in reference data"
            )

        # Find common genes between datasets
        common_genes = list(
            set(spatial_adata.var_names) & set(reference_adata.var_names)
        )
        min_common_genes = 100

        if len(common_genes) < min_common_genes:
            raise ValueError(
                f"Insufficient common genes.\n"
                f"  Found: {len(common_genes)} genes\n"
                f"  Required: {min_common_genes} genes\n"
                f"  Reference data: {reference_adata.n_vars} genes\n"
                f"  Spatial data: {spatial_adata.n_vars} genes\n\n"
                f"TIPS:\n"
                f"1. Check gene naming convention (mouse: 'Cd5l', human: 'CD5L')\n"
                f"2. Ensure both datasets are from the same species\n"
                f"3. Try using different reference dataset\n"
                f"4. Reduce min_common_genes parameter (current: {min_common_genes})"
            )

        # Prepare data - subset to common genes
        ref_data = reference_adata[:, common_genes].copy()
        spatial_data = spatial_adata[:, common_genes].copy()

        if context:
            await context.info(
                f"Training Tangram with {len(common_genes)} genes and {len(ref_data.obs[cell_type_key].unique())} cell types"
            )

        # Check data format - Tangram can work with normalized data but prefers raw counts
        import numpy as np

        # Check spatial data (sparse-aware)
        sp_has_negatives = spatial_data.X.min() < 0

        # Sample for decimal check
        sp_sample_size = min(1000, spatial_data.n_obs)
        sp_sample_genes = min(100, spatial_data.n_vars)
        sp_sample = spatial_data.X[:sp_sample_size, :sp_sample_genes]
        if hasattr(sp_sample, "toarray"):
            sp_sample_dense = sp_sample.toarray()
        else:
            sp_sample_dense = sp_sample
        sp_has_decimals = not np.allclose(
            sp_sample_dense, np.round(sp_sample_dense), atol=1e-6
        )

        # Check reference data (sparse-aware)
        ref_has_negatives = ref_data.X.min() < 0

        # Sample for decimal check
        ref_sample_size = min(1000, ref_data.n_obs)
        ref_sample_genes = min(100, ref_data.n_vars)
        ref_sample = ref_data.X[:ref_sample_size, :ref_sample_genes]
        if hasattr(ref_sample, "toarray"):
            ref_sample_dense = ref_sample.toarray()
        else:
            ref_sample_dense = ref_sample
        ref_has_decimals = not np.allclose(
            ref_sample_dense, np.round(ref_sample_dense), atol=1e-6
        )

        if (
            sp_has_negatives or sp_has_decimals or ref_has_negatives or ref_has_decimals
        ) and context:
            await context.warning(
                "Tangram is using normalized data. While Tangram can handle normalized data, "
                "it performs optimally with raw counts. Consider using raw count data for best results."
            )

        # Create density prior (normalized uniform distribution)
        if context:
            await context.info("Setting up Tangram with MuData...")

        # Create normalized density prior that sums to 1
        density_values = np.ones(spatial_data.n_obs)
        density_values = density_values / density_values.sum()  # Normalize to sum to 1
        spatial_data.obs["density_prior"] = density_values

        # Create MuData object combining spatial and single-cell data
        mdata = md.MuData({"sc_train": ref_data, "sp_train": spatial_data})

        # Setup MuData for Tangram
        Tangram.setup_mudata(
            mdata,
            density_prior_key="density_prior",
            modalities={
                "density_prior_key": "sp_train",
                "sc_layer": "sc_train",
                "sp_layer": "sp_train",
            },
        )

        # Create Tangram model
        target_count = max(1, int(spatial_data.n_obs * 0.1))  # Simple heuristic
        tangram_model = Tangram(mdata, constrained=True, target_count=target_count)

        if context:
            await context.info("Training Tangram model...")

        # Train model
        if use_gpu:
            tangram_model.train(max_epochs=n_epochs, accelerator="gpu")
        else:
            tangram_model.train(max_epochs=n_epochs)

        if context:
            await context.info("Tangram training completed")

        # Get cell type proportions
        if context:
            await context.info("Extracting cell type proportions...")

        # Get mapping matrix (shape: n_cells x n_spots)
        mapping_matrix = tangram_model.get_mapper_matrix()

        # Calculate cell type proportions
        cell_types = ref_data.obs[cell_type_key].unique()
        proportions_list = []

        for cell_type in cell_types:
            # Get cells of this type
            cell_mask = ref_data.obs[cell_type_key] == cell_type
            cell_indices = np.where(cell_mask)[0]

            # Sum mapping weights for this cell type across spots
            # mapping_matrix[cell_indices, :] gives weights for this cell type to all spots
            if len(cell_indices) > 0:
                cell_type_props = mapping_matrix[cell_indices, :].sum(
                    axis=0
                )  # Sum across cells, keep spots
            else:
                cell_type_props = np.zeros(spatial_data.n_obs)

            proportions_list.append(cell_type_props)

        # Create proportions DataFrame (transpose to get spots x cell_types)
        proportions_array = np.column_stack(proportions_list)

        # Create DataFrame without forced normalization
        proportions = pd.DataFrame(
            proportions_array, index=spatial_data.obs_names, columns=cell_types
        )

        # Process results transparently - Tangram mapper matrix returns "equivalent
        # cell counts" (not proportions), so we must normalize to get proportions
        proportions = _validate_and_process_proportions(
            proportions=proportions,
            method="Tangram",
            normalize=True,  # Convert counts to proportions
            context=context,
        )

        # Create statistics dictionary
        stats = _create_deconvolution_stats(
            proportions,
            common_genes,
            "Tangram",
            "GPU" if use_gpu else "CPU",
            n_epochs=n_epochs,
            use_gpu=use_gpu,
        )

        return proportions, stats

    except Exception as e:
        traceback.format_exc()
        error_msg = f"Tangram deconvolution failed: {str(e)}"
        # SPOTlight error occurred
        raise RuntimeError(error_msg) from e
