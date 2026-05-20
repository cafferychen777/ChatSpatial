"""
Cell type annotation tools for spatial transcriptomics data.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from collections import OrderedDict
from pathlib import Path
from typing import TYPE_CHECKING, Any, NamedTuple, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

if TYPE_CHECKING:
    import anndata as ad

    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import AnnotationResult
from ..models.data import AnnotationParameters
from ..utils.adata_utils import (
    ensure_categorical,
    ensure_counts_layer,
    ensure_unique_var_names_async,
    find_common_genes,
    get_cell_type_key,
    get_cluster_key,
    get_raw_data_source,
    get_spatial_key,
    shallow_copy_adata,
    to_dense,
    validate_obs_column,
)
from ..utils.dependency_manager import (
    is_available,
    require,
    validate_r_environment,
    validate_scvi_tools,
)
from ..utils.device_utils import get_device
from ..utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)


class AnnotationMethodOutput(NamedTuple):
    """Unified output from all annotation methods.

    This provides a consistent return type across all annotation methods,
    improving code clarity and preventing positional argument confusion.

    Attributes:
        cell_types: List of unique cell type names identified (deduplicated)
        counts: Mapping of cell type names to number of cells assigned
        confidence: Mapping of cell type names to confidence scores.
                   Empty dict indicates no confidence data available.
        tangram_mapping_score: Tangram-specific mapping quality score (only populated
                              by Tangram method, None for all other methods)
    """

    cell_types: list[str]
    counts: dict[str, int]
    confidence: dict[str, float]
    tangram_mapping_score: Optional[float] = None


# Supported annotation methods
# Confidence behavior by method:
#   - singler/tangram/sctype: Real confidence scores (correlation/probability/scoring)
#   - scanvi/cellassign: Partial confidence (when soft prediction available)
#   - mllmcelltype: No numeric confidence (LLM-based)
SUPPORTED_METHODS = {
    "tangram",
    "scanvi",
    "cellassign",
    "mllmcelltype",
    "sctype",
    "singler",
}


async def _annotate_with_singler(
    adata: "ad.AnnData",
    params: AnnotationParameters,
    ctx: "ToolContext",
    output_key: str,
    confidence_key: str,
    reference_adata: Optional[Any] = None,
) -> AnnotationMethodOutput:
    """Annotate cell types using SingleR reference-based method.

    Memory: SingleR (singler-py) natively supports sparse matrices - no toarray() needed.
    This saves ~1.3 GB for typical 10K cells × 20K genes datasets.

    Confidence scores are transformed to [0, 1] range:
    - Delta scores: 1 - exp(-delta), where delta is best-vs-second-best gap
    - Correlation scores: max(0, r), where negative correlations map to 0
    """
    # Validate and import dependencies
    require("singler", ctx, feature="SingleR annotation")
    require("singlecellexperiment", ctx, feature="SingleR annotation")
    import singler

    # Optional: check for celldex (import to module-level alias to avoid redef)
    celldex_module: Any = None
    if is_available("celldex"):
        import celldex as _celldex_mod

        celldex_module = _celldex_mod

    # Get expression matrix - prefer normalized data
    # IMPORTANT: Ensure test_mat dimensions match adata.var_names (used in test_features)
    if "X_normalized" in adata.layers:
        test_mat = adata.layers["X_normalized"]
    else:
        # Use get_raw_data_source (single source of truth) with prefer_complete_genes=False
        # to ensure dimensions match current adata.var_names
        raw_result = get_raw_data_source(adata, prefer_complete_genes=False)
        test_mat = raw_result.X
        # Always log-normalize raw data for SingleR.
        # get_raw_data_source returns raw counts by definition, so we cannot
        # rely on adata.uns["log1p"] (which only records whether adata.X was
        # log1p'd during preprocessing, not this raw source).
        if raw_result.is_integer_counts:
            await ctx.warning(
                "Data appears to be raw counts. "
                "Applying log1p normalization for SingleR..."
            )
            test_mat = np.log1p(test_mat)
        elif raw_result.has_negatives:
            await ctx.warning(
                "Data contains negative values (likely scaled). "
                "Skipping log1p to avoid NaN for SingleR."
            )
        elif "log1p" not in adata.uns:
            await ctx.warning(
                "Data may not be log-normalized. "
                "Applying log1p for SingleR..."
            )
            test_mat = np.log1p(test_mat)

    # Transpose for SingleR (genes x cells)
    test_mat = test_mat.T

    # Ensure gene names are strings
    test_features = [str(x) for x in adata.var_names]

    # Prepare reference
    reference_name = getattr(params, "singler_reference", None)
    reference_data_id = getattr(params, "reference_data_id", None)
    use_integrated = bool(getattr(params, "singler_integrated", False))
    n_threads = getattr(params, "n_threads", 4)

    ref_data = None
    ref_labels = None
    ref_features_to_use = None  # Only set when using custom reference (not celldex)

    def _extract_reference_labels(ref_obj: Any, ref_id: str) -> Any:
        """Extract labels from celldex reference using prioritized columns."""
        for label_col in ("label.main", "label.fine", "cell_type"):
            try:
                return ref_obj.get_column_data().column(label_col)
            except Exception:
                continue
        raise DataNotFoundError(f"Could not find labels in reference {ref_id}")

    # Priority: reference_name > reference_data_id > default
    if reference_name and celldex_module:
        # Integrated mode supports multiple references via comma-separated names.
        if use_integrated and "," in reference_name:
            ref_names = [x.strip() for x in reference_name.split(",") if x.strip()]
            if len(ref_names) < 2:
                raise ParameterError(
                    "singler_integrated=True requires at least two reference names "
                    "when using comma-separated singler_reference."
                )
            refs: list[Any] = []
            labels: list[Any] = []
            for ref_name in ref_names:
                ref_obj = celldex_module.fetch_reference(
                    ref_name, "2024-02-26", realize_assays=True
                )
                refs.append(ref_obj)
                labels.append(_extract_reference_labels(ref_obj, ref_name))
            ref_data = refs
            ref_labels = labels
        else:
            ref = celldex_module.fetch_reference(
                reference_name, "2024-02-26", realize_assays=True
            )
            ref_labels = _extract_reference_labels(ref, reference_name)
            ref_data = ref

    elif reference_data_id and reference_adata is not None:
        # Use provided reference data (passed from main function via ctx.get_adata())
        # Handle duplicate gene names
        await ensure_unique_var_names_async(reference_adata, ctx, "reference data")
        if await ensure_unique_var_names_async(adata, ctx, "query data") > 0:
            # Update test_features after fixing
            test_features = [str(x) for x in adata.var_names]

        # Get reference expression matrix
        if "X_normalized" in reference_adata.layers:
            ref_mat = reference_adata.layers["X_normalized"]
        else:
            ref_mat = reference_adata.X

        # Ensure log-normalization for reference (skip if negatives present)
        ref_min = ref_mat.min() if hasattr(ref_mat, "min") else np.min(ref_mat)
        if ref_min < 0:
            await ctx.warning(
                "Reference data contains negative values (likely scaled). "
                "Skipping log1p to avoid NaN for SingleR."
            )
        elif "log1p" not in reference_adata.uns:
            await ctx.warning(
                "Reference data may not be log-normalized. Applying log1p..."
            )
            ref_mat = np.log1p(ref_mat)

        # Transpose for SingleR (genes x cells)
        ref_mat = ref_mat.T
        ref_features = [str(x) for x in reference_adata.var_names]

        # Check gene overlap
        common_genes = find_common_genes(test_features, ref_features)

        if len(common_genes) < min(50, len(test_features) * 0.1):
            raise DataError(
                f"Insufficient gene overlap for SingleR: only {len(common_genes)} common genes "
                f"(test: {len(test_features)}, reference: {len(ref_features)})"
            )

        # Get labels from reference - check various common column names
        # cell_type_key is now required (no default value)
        cell_type_key = params.cell_type_key
        if cell_type_key is None:
            raise ParameterError(
                "cell_type_key is required for SingleR annotation with custom reference"
            )

        validate_obs_column(reference_adata, cell_type_key, "Cell type")

        ref_labels = list(reference_adata.obs[cell_type_key])

        # For SingleR, pass the actual expression matrix directly (not SCE)
        # This has been shown to work better in testing
        ref_data = ref_mat
        ref_features_to_use = ref_features  # Keep reference features for gene matching

    elif celldex_module:
        # Use default reference
        ref = celldex_module.fetch_reference(
            "blueprint_encode", "2024-02-26", realize_assays=True
        )
        ref_labels = _extract_reference_labels(ref, "blueprint_encode")
        ref_data = ref
    else:
        raise DataNotFoundError(
            "No reference data. Provide reference_data_id or singler_reference."
        )

    # Run SingleR annotation
    delta_scores = None

    if use_integrated and not isinstance(ref_data, list):
        raise ParameterError(
            "singler_integrated=True requires multiple references. "
            "Provide comma-separated singler_reference names (e.g. 'hpca,blueprint')."
        )

    if use_integrated and isinstance(ref_data, list):
        single_results, integrated = singler.annotate_integrated(
            test_mat,
            ref_data=ref_data,
            ref_labels=ref_labels,
            test_features=test_features,
            n_threads=n_threads,
        )
        best_labels = integrated.column("best_label")
        scores = integrated.column("scores")
    else:
        # Build kwargs for annotate_single
        annotate_kwargs = {
            "test_data": test_mat,
            "test_features": test_features,
            "ref_data": ref_data,
            "ref_labels": ref_labels,
        }

        # Add ref_features if we're using custom reference data (not celldex)
        if ref_features_to_use is not None:
            annotate_kwargs["ref_features"] = ref_features_to_use

        results = singler.annotate_single(**annotate_kwargs)
        best_labels = results.column("best")
        scores = results.column("scores")

        # Try to get delta scores for confidence (higher delta = higher confidence).
        # Wrapped in a single try/except because results.column() may return
        # non-iterable types depending on the SingleR version.
        try:
            delta_scores = results.column("delta")
            n_deltas = len(delta_scores)
            if n_deltas > 0:
                low_delta = sum(1 for d in delta_scores if d is not None and d < 0.05)
                if low_delta > n_deltas * 0.3:
                    await ctx.warning(
                        f"{low_delta}/{n_deltas} cells have low "
                        f"confidence scores (delta < 0.05)"
                    )
        except Exception as exc:
            import logging

            logging.getLogger(__name__).debug(
                "Delta scores unavailable (%s: %s) — confidence from scores only",
                type(exc).__name__,
                exc,
            )

    # Process results
    cell_types = list(best_labels)
    # Preserve deterministic cell-type order based on first appearance.
    unique_types = list(dict.fromkeys(cell_types))
    counts = pd.Series(cell_types).value_counts().to_dict()

    # Calculate per-cell confidence and per-type summary
    n_cells = len(cell_types)
    per_cell_confidence = np.zeros(n_cells, dtype=float)
    confidence_scores: dict[str, float] = {}

    # Prefer delta scores (more meaningful confidence measure)
    # delta→confidence: 1 - exp(-delta), clamped to [0, 1].
    # Negative delta (second-best > best) → confidence 0.
    if delta_scores is not None:
        try:
            for i in range(n_cells):
                if i < len(delta_scores) and delta_scores[i] is not None:
                    raw = 1.0 - np.exp(-delta_scores[i])
                    per_cell_confidence[i] = round(max(0.0, raw), 3)
        except Exception:
            per_cell_confidence[:] = 0.0

    # Fall back to correlation scores if delta not available
    if per_cell_confidence.sum() == 0 and scores is not None:
        try:
            try:
                scores_df = pd.DataFrame(scores.to_dict())
            except AttributeError:
                scores_df = pd.DataFrame(
                    scores.to_numpy() if hasattr(scores, "to_numpy") else scores
                )

            for i, ct in enumerate(cell_types):
                if ct in scores_df.columns:
                    val = float(scores_df.iloc[i][ct])
                    per_cell_confidence[i] = round(max(0.0, val), 3)
        except Exception:
            per_cell_confidence[:] = 0.0

    # Per-type summary (for result metadata, not for obs storage)
    for cell_type in unique_types:
        mask = np.array([ct == cell_type for ct in cell_types])
        if mask.any():
            confidence_scores[cell_type] = round(
                float(per_cell_confidence[mask].mean()), 3
            )

    # Add to AnnData (keys provided by caller for single-point control)
    adata.obs[output_key] = cell_types
    ensure_categorical(adata, output_key)

    has_confidence = per_cell_confidence.sum() > 0
    if has_confidence:
        adata.obs[confidence_key] = per_cell_confidence

    return AnnotationMethodOutput(
        cell_types=unique_types,
        counts=counts,
        confidence=confidence_scores if has_confidence else {},
    )


async def _annotate_with_tangram(
    adata: "ad.AnnData",
    params: AnnotationParameters,
    ctx: "ToolContext",
    output_key: str,
    confidence_key: str,
    reference_adata: Optional[Any] = None,
    tangram_ct_pred_key: str = "tangram_ct_pred",
) -> AnnotationMethodOutput:
    """Annotate cell types using Tangram method"""
    # Validate dependencies with comprehensive error reporting
    require("tangram", ctx, feature="Tangram annotation")
    import tangram as tg

    # Check if reference data is provided
    if reference_adata is None:
        raise ParameterError("Tangram requires reference_data_id parameter.")

    # Use reference single-cell data (passed from main function via ctx.get_adata())
    adata_sc_original = reference_adata

    # ===== CRITICAL FIX: Use raw data for Tangram to preserve gene name case =====
    # Issue: Preprocessed data may have lowercase gene names, while reference has uppercase
    # This causes 0 overlapping genes and Tangram mapping failure (all NaN results)
    # Solution: Use adata.raw which preserves original gene names and full gene set
    if adata.raw is not None:
        # Use raw data which preserves original gene names
        adata_sp = adata.raw.to_adata()
        # Preserve spatial coordinates from preprocessed data
        spatial_key = get_spatial_key(adata)
        if spatial_key:
            adata_sp.obsm[spatial_key] = adata.obsm[spatial_key].copy()
    else:
        adata_sp = adata
        await ctx.warning(
            "Raw data not available - may have gene name mismatches with reference"
        )
    # =============================================================================

    # Handle duplicate gene names
    await ensure_unique_var_names_async(adata_sc_original, ctx, "reference data")
    await ensure_unique_var_names_async(adata_sp, ctx, "spatial data")

    # Determine training genes
    training_genes = params.training_genes

    if training_genes is None:
        # Use marker genes if available
        if params.marker_genes:
            # Flatten marker genes dictionary
            training_genes = []
            for genes in params.marker_genes.values():
                training_genes.extend(genes)
            training_genes = list(dict.fromkeys(training_genes))  # Stable deduplication
        else:
            # Use highly variable genes
            if "highly_variable" not in adata_sc_original.var:
                raise DataNotFoundError(
                    "HVGs not found in reference data. Run preprocessing first."
                )
            training_genes = list(
                adata_sc_original.var_names[adata_sc_original.var.highly_variable]
            )

    # Memory optimization: Use shallow copy (~99% memory savings)
    # Tangram's pp_adatas only adds to uns (training_genes), doesn't modify X
    adata_sc = shallow_copy_adata(adata_sc_original)

    # Preprocess data for Tangram
    tg.pp_adatas(adata_sc, adata_sp, genes=training_genes)

    # Set mapping mode
    mode = params.tangram_mode
    cluster_label = params.cluster_label

    if mode == "clusters" and cluster_label is None:
        await ctx.warning(
            "Cluster label not provided for 'clusters' mode. Using default cell type annotation if available."
        )
        # Try to find a cell type or cluster annotation in the reference data
        cluster_label = get_cell_type_key(adata_sc) or get_cluster_key(adata_sc)

        if cluster_label is None:
            raise ParameterError(
                "No cluster label found. Provide cluster_label parameter."
            )

    # Device selection (supports CUDA, MPS, and CPU)
    device = get_device(prefer_gpu=params.tangram_use_gpu)
    if params.tangram_use_gpu and device == "cpu":
        await ctx.warning("GPU requested but not available - using CPU")

    # Run Tangram mapping with enhanced parameters
    mapping_args: dict[str, str | int | float | None] = {
        "mode": mode,
        "num_epochs": params.n_epochs,
        "device": device,
        "density_prior": params.tangram_density_prior,  # Add density prior
        "learning_rate": params.tangram_learning_rate,  # Add learning rate
    }

    # Add optional regularization parameters
    if params.tangram_lambda_r is not None:
        mapping_args["lambda_r"] = params.tangram_lambda_r

    if params.tangram_lambda_neighborhood is not None:
        mapping_args["lambda_neighborhood"] = params.tangram_lambda_neighborhood

    if mode == "clusters":
        mapping_args["cluster_label"] = cluster_label

    ad_map = tg.map_cells_to_space(adata_sc, adata_sp, **mapping_args)

    # Get mapping score from training history
    tangram_mapping_score = 0.0  # Default score
    try:
        if "training_history" in ad_map.uns:
            history = ad_map.uns["training_history"]

            # Extract score from main_loss (which is actually a similarity score, higher is better)
            if (
                isinstance(history, dict)
                and "main_loss" in history
                and len(history["main_loss"]) > 0
            ):
                import re

                last_value = history["main_loss"][-1]

                # Extract value from tensor string if needed
                if isinstance(last_value, str):
                    # Handle tensor string format: 'tensor(0.9050, grad_fn=...)'
                    match = re.search(r"tensor\(([-\d.]+)", last_value)
                    if match:
                        tangram_mapping_score = float(match.group(1))
                    else:
                        # Try direct conversion
                        try:
                            tangram_mapping_score = float(last_value)
                        except Exception:
                            tangram_mapping_score = 0.0
                else:
                    tangram_mapping_score = float(last_value)

            else:
                error_msg = (
                    f"Tangram history format not recognized: {type(history).__name__}. "
                    f"Upgrade tangram-sc: pip install --upgrade tangram-sc"
                )
                raise ProcessingError(error_msg)
    except Exception as score_error:
        raise ProcessingError(
            f"Tangram mapping completed but score extraction failed: {score_error}"
        ) from score_error

    # Compute validation metrics if requested
    if params.tangram_compute_validation:
        try:
            scores = tg.compare_spatial_geneexp(ad_map, adata_sp)
            adata_sp.uns["tangram_validation_scores"] = scores
        except Exception as val_error:
            await ctx.warning(f"Could not compute validation metrics: {val_error}")

    # Project genes if requested
    if params.tangram_project_genes:
        try:
            ad_ge = tg.project_genes(ad_map, adata_sc)
            adata_sp.obsm["tangram_gene_predictions"] = ad_ge.X
        except Exception as gene_error:
            await ctx.warning(f"Could not project genes: {gene_error}")

    # Project cell annotations to space using proper API function
    try:
        # Determine annotation column
        annotation_col = None
        if mode == "clusters" and cluster_label:
            annotation_col = cluster_label
        else:
            # cell_type_key is now required (no auto-detect)
            if params.cell_type_key not in adata_sc.obs:
                # Improved error message showing available columns
                available_cols = list(adata_sc.obs.columns)
                categorical_cols = [
                    col
                    for col in available_cols
                    if adata_sc.obs[col].dtype.name in ["object", "category"]
                ]

                raise ParameterError(
                    f"Cell type column '{params.cell_type_key}' not found. "
                    f"Available: {categorical_cols[:5]}"
                )

            annotation_col = params.cell_type_key

        # annotation_col is guaranteed to be set (either from cluster_label or cell_type_key)
        tg.project_cell_annotations(ad_map, adata_sp, annotation=annotation_col)
    except Exception as proj_error:
        await ctx.warning(f"Could not project cell annotations: {proj_error}")
        # Continue without projection

    # Get cell type predictions (keys provided by caller for single-point control)
    cell_types: list[str] = []
    counts: dict[str, int] = {}
    confidence_scores: dict[str, float] = {}

    if "tangram_ct_pred" in adata_sp.obsm:
        cell_type_df = adata_sp.obsm["tangram_ct_pred"]

        # Get cell types and counts
        cell_types = list(cell_type_df.columns)

        # ===== Row normalization for proper probability calculation =====
        # tangram_ct_pred contains unnormalized density/abundance values, NOT probabilities
        # Row sums can be != 1.0 and values can exceed 1.0
        # We normalize to convert densities → probability distributions
        row_sums = cell_type_df.sum(axis=1)
        zero_mask = row_sums == 0
        cell_type_prob = cell_type_df.div(row_sums, axis=0)

        # Handle zero-sum rows (spots with no valid mapping)
        if zero_mask.any():
            cell_type_prob.loc[zero_mask] = 0.0

        # Validation: Ensure normalized values are valid probabilities
        if not (cell_type_prob.values >= 0).all():
            await ctx.warning(
                "Some normalized probabilities are negative - data quality issue"
            )
        if not (cell_type_prob.values <= 1.0).all():
            await ctx.warning(
                "Some normalized probabilities exceed 1.0 - normalization failed"
            )
        valid_rows = ~zero_mask
        if valid_rows.any() and not np.allclose(
            cell_type_prob.loc[valid_rows].sum(axis=1), 1.0
        ):
            await ctx.warning(
                "Row sums don't equal 1.0 after normalization - numerical issue"
            )

        # Assign cell type based on highest probability
        adata_sp.obs[output_key] = cell_type_prob.idxmax(axis=1)
        if zero_mask.any():
            adata_sp.obs.loc[zero_mask, output_key] = "unassigned"
        ensure_categorical(adata_sp, output_key)

        # Write per-cell confidence (max probability for assigned type)
        per_cell_conf = cell_type_prob.max(axis=1)
        if zero_mask.any():
            per_cell_conf.loc[zero_mask] = 0.0
        adata_sp.obs[confidence_key] = per_cell_conf.values

        # Include "unassigned" in the cell_types list so that downstream
        # metadata (n_cell_types, cell_types) matches the actual label space.
        if zero_mask.any() and "unassigned" not in cell_types:
            cell_types.append("unassigned")

        # Get counts
        counts = adata_sp.obs[output_key].value_counts().to_dict()

        # Calculate per-cell-type aggregated confidence scores
        confidence_scores = {}
        for cell_type in cell_types:
            if cell_type == "unassigned":
                # Zero-sum rows have confidence 0 by definition
                confidence_scores[cell_type] = 0.0
                continue
            cells_of_type = adata_sp.obs[output_key] == cell_type
            if np.sum(cells_of_type) > 0:
                mean_prob = cell_type_prob.loc[cells_of_type, cell_type].mean()
                confidence_scores[cell_type] = round(float(mean_prob), 3)

    else:
        await ctx.warning("No cell type predictions found in Tangram results")

    # Validate results before returning
    if not cell_types:
        raise ProcessingError(
            "Tangram mapping failed - no cell type predictions generated"
        )

    if tangram_mapping_score <= 0:
        await ctx.warning(
            f"Tangram mapping score is suspiciously low: {tangram_mapping_score}"
        )

    # ===== Copy results from adata_sp back to original adata =====
    # Since adata_sp was created from adata.raw (different object), we need to
    # transfer the Tangram results back to the original adata for downstream use
    if adata_sp is not adata:
        # Copy cell type assignments
        if output_key in adata_sp.obs:
            adata.obs[output_key] = adata_sp.obs[output_key]

        # Copy per-cell confidence scores
        if confidence_key in adata_sp.obs:
            adata.obs[confidence_key] = adata_sp.obs[confidence_key]

        # Copy tangram_ct_pred from obsm (parametrized key)
        if "tangram_ct_pred" in adata_sp.obsm:
            adata.obsm[tangram_ct_pred_key] = adata_sp.obsm["tangram_ct_pred"]

        # Copy tangram_gene_predictions if they exist (parametrized key)
        if "tangram_gene_predictions" in adata_sp.obsm:
            gene_pred_key = tangram_ct_pred_key.replace(
                "tangram_ct_pred", "tangram_gene_predictions"
            )
            adata.obsm[gene_pred_key] = adata_sp.obsm["tangram_gene_predictions"]

    return AnnotationMethodOutput(
        cell_types=cell_types,
        counts=counts,
        confidence=confidence_scores,
        tangram_mapping_score=tangram_mapping_score,
    )


async def _annotate_with_scanvi(
    adata: "ad.AnnData",
    params: AnnotationParameters,
    ctx: "ToolContext",
    output_key: str,
    confidence_key: str,
    reference_adata: Optional[Any] = None,
) -> AnnotationMethodOutput:
    """Annotate cell types using scANVI (semi-supervised variational inference).

    scANVI (single-cell ANnotation using Variational Inference) is a deep learning
    method for transferring cell type labels from reference to query data using
    semi-supervised learning with variational autoencoders.

    Official Implementation: scvi-tools (https://scvi-tools.org)
    Reference: Xu et al. (2021) "Probabilistic harmonization and annotation of
               single-cell transcriptomics data with deep generative models"

    Method Overview:
        1. Trains on reference data with known cell type labels
        2. Learns shared latent representation between reference and query
        3. Transfers labels to query data via probabilistic predictions
        4. Supports batch correction and semi-supervised training

    Requirements:
        - reference_data_id: Must point to preprocessed single-cell reference data
        - cell_type_key: Column in reference data containing cell type labels
        - Both datasets must have 'counts' layer (raw counts, not normalized)
        - Sufficient gene overlap between reference and query data

    Parameters (via AnnotationParameters):
        Core Architecture:
            - scanvi_n_latent (default: 10): Latent space dimensions
            - scanvi_n_hidden (default: 128): Hidden layer units
            - scanvi_n_layers (default: 1): Number of layers
            - scanvi_dropout_rate (default: 0.1): Dropout for regularization

        Training Strategy:
            - scanvi_use_scvi_pretrain (default: True): Use SCVI pretraining
            - scanvi_scvi_epochs (default: 200): SCVI pretraining epochs
            - n_epochs (default: 100): SCANVI training epochs
            - scanvi_query_epochs (default: 100): Query data training epochs

        Advanced:
            - scanvi_unlabeled_category (default: "Unknown"): Label for unlabeled cells
            - scanvi_n_samples_per_label (default: 100): Samples per label
            - batch_key: For batch correction (optional)

    Official Recommendations (scvi-tools):
        For large integration tasks:
        - scanvi_n_layers: 2
        - scanvi_n_latent: 30
        - scanvi_scvi_epochs: 300 (SCVI pretraining)
        - n_epochs: 100 (SCANVI training)
        - scanvi_query_epochs: 100
        - Gene selection: 1000-10000 HVGs recommended

    Empirical Adjustments (not official):
        For small datasets (<1000 genes or <1000 cells):
        - scanvi_n_latent: 3-5 (may prevent NaN/gradient explosion)
        - scanvi_dropout_rate: 0.2-0.3 (may improve regularization)
        - scanvi_use_scvi_pretrain: False (may reduce complexity)
        - n_epochs: 50 (may prevent overfitting)
        - scanvi_query_epochs: 50

    Common Issues:
        - NaN errors during training: Try reducing n_latent or increasing dropout_rate
        - Low confidence scores: Try increasing training epochs or check gene overlap
        - Memory issues: Reduce batch size or use GPU

    Returns:
        Tuple of (cell_types, counts, confidence_scores, None):
        - cell_types: List of predicted cell type categories
        - counts: Dict mapping cell types to number of cells
        - confidence_scores: Dict mapping cell types to mean prediction probability
        - None: (compatibility placeholder)

    Example:
        params = AnnotationParameters(
            method="scanvi",
            reference_data_id="reference_sc",
            cell_type_key="cell_types",
            scanvi_n_latent=5,              # For small dataset
            scanvi_dropout_rate=0.2,        # Better regularization
            scanvi_use_scvi_pretrain=False, # Simpler training
            n_epochs=50,                    # Prevent overfitting
        )
    """

    # Validate dependencies with comprehensive error reporting
    scvi = validate_scvi_tools(ctx, components=["SCANVI"])

    # Check if reference data is provided
    if reference_adata is None:
        raise ParameterError("scANVI requires reference_data_id parameter.")

    # Use reference single-cell data (passed from main function via ctx.get_adata())
    adata_ref_original = reference_adata

    # Handle duplicate gene names
    await ensure_unique_var_names_async(adata_ref_original, ctx, "reference data")
    await ensure_unique_var_names_async(adata, ctx, "query data")

    # Gene alignment
    common_genes = find_common_genes(adata_ref_original.var_names, adata.var_names)

    if len(common_genes) < min(100, adata_ref_original.n_vars * 0.5):
        raise DataError(
            f"Insufficient gene overlap: Only {len(common_genes)} common genes found. "
            f"Reference has {adata_ref_original.n_vars}, query has {adata.n_vars} genes."
        )

    # COW FIX: Operate on temporary copies for gene subsetting
    # This prevents loss of HVG information in the original adata
    if len(common_genes) < adata_ref_original.n_vars:
        await ctx.warning(
            f"Subsetting to {len(common_genes)} common genes for ScanVI training "
            f"(reference: {adata_ref_original.n_vars}, query: {adata.n_vars})"
        )
        # Create subsets for ScanVI (not modifying originals)
        adata_ref = adata_ref_original[:, common_genes].copy()
        adata_subset = adata[:, common_genes].copy()
    else:
        # No subsetting needed - use shallow copy for ~99% memory savings
        # (shares X/layers, copies only obs/uns/var which are modified)
        adata_ref = shallow_copy_adata(adata_ref_original)
        adata_subset = shallow_copy_adata(adata)

    # Data validation
    if "log1p" not in adata_ref.uns:
        await ctx.warning("Reference data may not be log-normalized")
    if "highly_variable" not in adata_ref.var:
        await ctx.warning("No highly variable genes detected in reference")

    # Get parameters
    cell_type_key = getattr(params, "cell_type_key", "cell_type")
    batch_key = getattr(params, "batch_key", None)

    # Optional SCVI Pretraining
    if params.scanvi_use_scvi_pretrain:
        # Ensure counts layer exists (consistent with direct branch)
        ensure_counts_layer(
            adata_ref,
            error_message="scANVI requires raw counts in layers['counts'].",
        )

        # Setup for SCVI with labels (required for SCANVI conversion)
        validate_obs_column(
            adata_ref, cell_type_key, "Cell type column (reference data)"
        )

        # SCVI needs to know about labels for later SCANVI conversion
        scvi.model.SCVI.setup_anndata(
            adata_ref,
            labels_key=cell_type_key,
            batch_key=batch_key,
            layer="counts",
        )

        # Train SCVI
        scvi_model = scvi.model.SCVI(
            adata_ref,
            n_latent=params.scanvi_n_latent,
            n_hidden=params.scanvi_n_hidden,
            n_layers=params.scanvi_n_layers,
            dropout_rate=params.scanvi_dropout_rate,
        )

        scvi_model.train(
            max_epochs=params.scanvi_scvi_epochs,
            early_stopping=True,
            check_val_every_n_epoch=params.scanvi_check_val_every_n_epoch,
        )

        # Convert to SCANVI (no need for setup_anndata, it uses SCVI's setup)
        model = scvi.model.SCANVI.from_scvi_model(
            scvi_model, params.scanvi_unlabeled_category
        )

        # Train SCANVI (fewer epochs needed after pretraining)
        # Use configurable epochs (default: 20, official recommendation after pretraining)
        model.train(
            max_epochs=params.scanvi_scanvi_epochs,
            n_samples_per_label=params.scanvi_n_samples_per_label,
            early_stopping=True,
        )

    else:
        # Direct SCANVI training (existing approach)
        # Ensure counts layer exists (create from adata.raw if needed)
        ensure_counts_layer(
            adata_ref,
            error_message="scANVI requires raw counts in layers['counts'].",
        )

        # Setup AnnData for scANVI
        scvi.model.SCANVI.setup_anndata(
            adata_ref,
            labels_key=cell_type_key,
            unlabeled_category=params.scanvi_unlabeled_category,
            batch_key=batch_key,
            layer="counts",
        )

        # Create scANVI model
        model = scvi.model.SCANVI(
            adata_ref,
            n_hidden=params.scanvi_n_hidden,
            n_latent=params.scanvi_n_latent,
            n_layers=params.scanvi_n_layers,
            dropout_rate=params.scanvi_dropout_rate,
        )

        model.train(
            max_epochs=params.scanvi_scanvi_epochs,
            n_samples_per_label=params.scanvi_n_samples_per_label,
            early_stopping=True,
            check_val_every_n_epoch=params.scanvi_check_val_every_n_epoch,
        )

    # Query data preparation
    adata_subset.obs[cell_type_key] = params.scanvi_unlabeled_category

    # Setup query data (batch handling)
    if batch_key and batch_key not in adata_subset.obs:
        adata_subset.obs[batch_key] = "query_batch"

    # Ensure counts layer exists for query data (create from adata.raw if needed)
    ensure_counts_layer(
        adata_subset,
        error_message="scANVI requires raw counts in layers['counts'].",
    )

    scvi.model.SCANVI.setup_anndata(
        adata_subset,
        labels_key=cell_type_key,
        unlabeled_category=params.scanvi_unlabeled_category,
        batch_key=batch_key,
        layer="counts",
    )

    # Transfer model to spatial data with proper parameters
    spatial_model = scvi.model.SCANVI.load_query_data(adata_subset, model)

    # ===== Improved Query Training (NEW) =====
    spatial_model.train(
        max_epochs=params.scanvi_query_epochs,  # Default: 100 (was 50)
        early_stopping=True,
        plan_kwargs=dict(weight_decay=0.0),  # Critical: preserve reference space
        check_val_every_n_epoch=params.scanvi_check_val_every_n_epoch,
    )

    # COW FIX: Get predictions from adata_subset, then add to original adata
    predictions = spatial_model.predict()
    adata_subset.obs[cell_type_key] = predictions
    ensure_categorical(adata_subset, cell_type_key)

    # Extract results from adata_subset
    cell_types = list(adata_subset.obs[cell_type_key].cat.categories)
    counts = adata_subset.obs[cell_type_key].value_counts().to_dict()

    # Get per-cell confidence from soft predictions
    per_cell_confidence = np.zeros(adata_subset.n_obs, dtype=float)
    confidence_scores: dict[str, float] = {}
    try:
        probs = spatial_model.predict(soft=True)
        predicted_labels = adata_subset.obs[cell_type_key].values

        if isinstance(probs, pd.DataFrame):
            # Each cell's confidence = P(predicted_class) for that cell
            for i, ct in enumerate(predicted_labels):
                if ct in probs.columns:
                    per_cell_confidence[i] = round(float(probs.iloc[i][ct]), 3)
        elif hasattr(probs, "shape") and probs.ndim == 2:
            # Array form: use argmax index from predictions
            pred_idx = probs.argmax(axis=1)
            for i in range(len(pred_idx)):
                per_cell_confidence[i] = round(float(probs[i, pred_idx[i]]), 3)

        # Per-type summary for metadata
        for cell_type in cell_types:
            mask = predicted_labels == cell_type
            if mask.any():
                confidence_scores[cell_type] = round(
                    float(per_cell_confidence[mask].mean()), 2
                )
    except Exception as e:
        await ctx.warning(f"Could not get confidence scores: {e}")
        confidence_scores = {}

    # COW FIX: Add prediction results to original adata.obs using output_key
    adata.obs[output_key] = adata_subset.obs[cell_type_key].values
    ensure_categorical(adata, output_key)

    # Store per-cell confidence if available
    has_confidence = per_cell_confidence.sum() > 0
    if has_confidence:
        adata.obs[confidence_key] = per_cell_confidence

    return AnnotationMethodOutput(
        cell_types=cell_types,
        counts=counts,
        confidence=confidence_scores if has_confidence else {},
    )


async def _annotate_with_mllmcelltype(
    adata: "ad.AnnData",
    params: AnnotationParameters,
    ctx: "ToolContext",
    output_key: str,
    confidence_key: str,
) -> AnnotationMethodOutput:
    """Annotate cell types using mLLMCellType (LLM-based) method.

    Supports both single-model and multi-model consensus annotation.

    Single Model Mode (default):
        - Uses one LLM for annotation
        - Fast and cost-effective
        - Providers: openai, anthropic, gemini, deepseek, qwen, zhipu, stepfun, minimax, grok, openrouter
        - Default models: openai="gpt-5", anthropic="claude-sonnet-4-20250514", gemini="gemini-2.5-pro-preview-03-25"
        - Latest recommended: "gpt-5", "claude-sonnet-4-5-20250929", "claude-opus-4-1-20250805", "gemini-2.5-pro"

    Multi-Model Consensus Mode (set mllm_use_consensus=True):
        - Uses multiple LLMs for collaborative annotation
        - Higher accuracy through consensus
        - Provides uncertainty metrics (consensus proportion, entropy)
        - Structured deliberation for controversial clusters

    Parameters (via AnnotationParameters):
        - cluster_label: Required. Cluster column in adata.obs
        - mllm_species: "human" or "mouse"
        - mllm_tissue: Tissue context (optional but recommended)
        - mllm_provider: LLM provider (single model mode)
        - mllm_model: Model name (None = use default for provider)
        - mllm_use_consensus: Enable multi-model consensus
        - mllm_models: List of models for consensus (e.g., ["gpt-5", "claude-sonnet-4-5-20250929"])
        - mllm_additional_context: Additional context for better annotation
        - mllm_base_urls: Custom API endpoints (useful for proxies)
    """

    # Validate dependencies with comprehensive error reporting
    require("mllmcelltype", ctx, feature="mLLMCellType annotation")
    import mllmcelltype

    # Validate clustering has been performed
    # cluster_label is now required for mLLMCellType (no default value)
    if not params.cluster_label:
        available_cols = list(adata.obs.columns)
        categorical_cols = [
            col
            for col in available_cols
            if adata.obs[col].dtype.name in ["object", "category"]
        ]

        raise ParameterError(
            f"cluster_label parameter is required for mLLMCellType method.\n\n"
            f"Available categorical columns (likely clusters):\n  {', '.join(categorical_cols[:15])}\n"
            f"{f'  ... and {len(categorical_cols)-15} more' if len(categorical_cols) > 15 else ''}\n\n"
            f"Common cluster column names: leiden, louvain, seurat_clusters, phenograph\n\n"
            f"Example: params = {{'cluster_label': 'leiden', ...}}"
        )

    cluster_key = params.cluster_label
    validate_obs_column(adata, cluster_key, "Cluster")

    # Find differentially expressed genes for each cluster

    sc.tl.rank_genes_groups(adata, cluster_key, method="wilcoxon")

    # Extract top marker genes for each cluster
    marker_genes_dict = {}
    n_genes = params.mllm_n_marker_genes

    for cluster in adata.obs[cluster_key].unique():
        # Get top genes for this cluster
        gene_names = adata.uns["rank_genes_groups"]["names"][str(cluster)][:n_genes]
        marker_genes_dict[f"Cluster_{cluster}"] = list(gene_names)

    # Prepare parameters for mllmcelltype
    species = params.mllm_species
    tissue = params.mllm_tissue
    additional_context = params.mllm_additional_context
    use_cache = params.mllm_use_cache
    base_urls = params.mllm_base_urls
    verbose = params.mllm_verbose
    force_rerun = params.mllm_force_rerun
    clusters_to_analyze = params.mllm_clusters_to_analyze

    # Check if using multi-model consensus or single model
    use_consensus = params.mllm_use_consensus

    try:
        if use_consensus:
            # Use interactive_consensus_annotation with multiple models
            models = params.mllm_models
            if not models:
                raise ParameterError(
                    "mllm_models parameter is required when mllm_use_consensus=True. "
                    "Provide a list of model names, e.g., ['gpt-5', 'claude-sonnet-4-5-20250929', 'gemini-2.5-pro']"
                )

            api_keys = params.mllm_api_keys
            consensus_threshold = params.mllm_consensus_threshold
            entropy_threshold = params.mllm_entropy_threshold
            max_discussion_rounds = params.mllm_max_discussion_rounds
            consensus_model = params.mllm_consensus_model

            # Call interactive_consensus_annotation
            consensus_results = mllmcelltype.interactive_consensus_annotation(
                marker_genes=marker_genes_dict,
                species=species,
                models=models,
                api_keys=api_keys,
                tissue=tissue,
                additional_context=additional_context,
                consensus_threshold=consensus_threshold,
                entropy_threshold=entropy_threshold,
                max_discussion_rounds=max_discussion_rounds,
                use_cache=use_cache,
                verbose=verbose,
                consensus_model=consensus_model,
                base_urls=base_urls,
                clusters_to_analyze=clusters_to_analyze,
                force_rerun=force_rerun,
            )

            # Extract consensus annotations
            annotations = consensus_results.get("consensus", {})

        else:
            # Use single model annotation
            provider = params.mllm_provider
            model = params.mllm_model
            api_key = params.mllm_api_key

            # Call annotate_clusters (single model)
            annotations = mllmcelltype.annotate_clusters(
                marker_genes=marker_genes_dict,
                species=species,
                provider=provider,
                model=model,
                api_key=api_key,
                tissue=tissue,
                additional_context=additional_context,
                use_cache=use_cache,
                base_urls=base_urls,
            )
    except ParameterError:
        raise
    except Exception as e:
        raise ProcessingError(f"mLLMCellType annotation failed: {e}") from e

    # Map cluster annotations back to cells
    cluster_to_celltype = {}
    for cluster_name, cell_type in annotations.items():
        # Extract cluster number from "Cluster_X" format
        cluster_id = cluster_name.replace("Cluster_", "")
        cluster_to_celltype[cluster_id] = cell_type

    # Apply cell type annotations to cells (key provided by caller)
    adata.obs[output_key] = adata.obs[cluster_key].astype(str).map(cluster_to_celltype)

    # Handle any unmapped clusters
    unmapped = adata.obs[output_key].isna()
    if unmapped.any():
        await ctx.warning(f"Found {unmapped.sum()} cells in unmapped clusters")
        adata.obs.loc[unmapped, output_key] = "Unknown"

    ensure_categorical(adata, output_key)

    # Get cell types and counts
    cell_types = list(adata.obs[output_key].unique())
    counts = adata.obs[output_key].value_counts().to_dict()

    # LLM-based annotations don't provide numeric confidence scores
    # We intentionally leave this empty rather than assigning misleading values
    return AnnotationMethodOutput(
        cell_types=cell_types,
        counts=counts,
        confidence={},
    )


async def _annotate_with_cellassign(
    adata: "ad.AnnData",
    params: AnnotationParameters,
    ctx: "ToolContext",
    output_key: str,
    confidence_key: str,
) -> AnnotationMethodOutput:
    """Annotate cell types using CellAssign method"""

    # Validate dependencies with comprehensive error reporting
    validate_scvi_tools(ctx, components=["CellAssign"])
    from scvi.external import CellAssign

    # Check if marker genes are provided
    if params.marker_genes is None:
        raise ParameterError(
            "CellAssign requires marker genes to be provided. "
            "Please specify marker_genes parameter with a dictionary of cell types and their marker genes."
        )

    marker_genes = params.marker_genes

    # Use get_raw_data_source (single source of truth) for complete gene coverage
    # Preprocessing filters genes to HVGs, but marker genes may not be in HVGs
    raw_result = get_raw_data_source(
        adata, prefer_complete_genes=True, require_integer_counts=True
    )
    all_genes = set(raw_result.var_names)
    gene_source = raw_result.source
    if raw_result.source != "raw":
        await ctx.warning(
            f"Using {raw_result.source} data for marker gene validation "
            f"({len(all_genes)} genes). Some marker genes may be missing. "
            f"Consider using unpreprocessed data for CellAssign."
        )

    # Validate marker genes exist in dataset
    valid_marker_genes = {}
    total_markers = sum(len(g) for g in marker_genes.values())
    markers_found = 0
    markers_missing = 0

    for cell_type, genes in marker_genes.items():
        existing_genes = [gene for gene in genes if gene in all_genes]
        missing_genes = [gene for gene in genes if gene not in all_genes]

        if existing_genes:
            valid_marker_genes[cell_type] = existing_genes
            markers_found += len(existing_genes)
            if missing_genes and len(missing_genes) > len(existing_genes):
                await ctx.warning(
                    f"Missing most markers for {cell_type}: {len(missing_genes)}/{len(genes)}"
                )
        else:
            markers_missing += len(genes)
            await ctx.warning(
                f"No marker genes found for {cell_type} - all {len(genes)} markers missing!"
            )

    if not valid_marker_genes:
        raise DataError(
            f"No valid marker genes found for any cell type. "
            f"Checked {total_markers} markers against {len(all_genes)} genes in {gene_source}. "
            f"If data was preprocessed, marker genes may have been filtered out. "
            f"Consider using unpreprocessed data or ensure marker genes are highly variable."
        )
    valid_cell_types = list(valid_marker_genes)

    # Create marker gene matrix as DataFrame (required by CellAssign API)
    all_marker_genes = []
    for genes in valid_marker_genes.values():
        all_marker_genes.extend(genes)
    available_marker_genes = list(
        dict.fromkeys(all_marker_genes)
    )  # Stable deduplication
    available_marker_genes_set = set(available_marker_genes)

    # Note: available_marker_genes cannot be empty here because valid_marker_genes
    # is already validated at line 1120 to have at least one cell type with genes

    # Create DataFrame with genes as index, cell types as columns
    marker_gene_matrix = pd.DataFrame(
        np.zeros((len(available_marker_genes), len(valid_cell_types))),
        index=available_marker_genes,
        columns=valid_cell_types,
    )

    # Fill marker matrix
    for cell_type in valid_cell_types:
        for gene in valid_marker_genes[cell_type]:
            if gene in available_marker_genes_set:
                marker_gene_matrix.loc[gene, cell_type] = 1

    # Compute size factors BEFORE subsetting (official CellAssign requirement)
    # Size factors MUST come from the same count source as the modeling matrix
    # (raw_result.X) to satisfy the statistical assumptions of the NB model.
    if "size_factors" not in adata.obs:
        sf_X = raw_result.X
        size_factors = sf_X.sum(axis=1)
        if hasattr(size_factors, "A1"):  # sparse matrix
            size_factors = size_factors.A1
        size_factors = np.asarray(size_factors).flatten()

        # Normalize and ensure positive
        size_factors = np.maximum(size_factors, 1e-6)
        mean_sf = np.mean(size_factors)
        size_factors_normalized = size_factors / mean_sf

        adata.obs["size_factors"] = pd.Series(
            size_factors_normalized, index=adata.obs.index
        )

    # Subset data to marker genes (size factors already computed)
    # Use raw data if available (contains all genes including markers)
    # raw_result already provides X and var_names from get_raw_data_source above
    if gene_source == "raw":
        import anndata as ad_module

        # Use raw_result directly instead of accessing adata.raw again
        gene_indices = [raw_result.var_names.get_loc(g) for g in available_marker_genes]
        adata_subset = ad_module.AnnData(
            X=raw_result.X[:, gene_indices],
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=available_marker_genes),
        )
    else:
        adata_subset = adata[:, available_marker_genes].copy()

    # Check for invalid values in the data
    X_array = to_dense(adata_subset.X)

    # Replace any NaN or Inf values with zeros
    if np.any(np.isnan(X_array)) or np.any(np.isinf(X_array)):
        X_array = np.nan_to_num(X_array, nan=0.0, posinf=0.0, neginf=0.0)
        adata_subset.X = X_array

    # Additional data cleaning for CellAssign compatibility
    # Check for genes with zero variance (which cause numerical issues in CellAssign)
    gene_vars = np.var(X_array, axis=0)
    zero_var_genes = gene_vars == 0
    if np.any(zero_var_genes):
        adata_subset.var_names[zero_var_genes].tolist()
        await ctx.warning(
            f"Found {np.sum(zero_var_genes)} genes with zero variance. "
            f"CellAssign may have numerical issues with these genes."
        )
        # Don't raise error, just warn - CellAssign might handle it

    # Ensure data is non-negative (CellAssign expects count-like data)
    if np.any(X_array < 0):
        X_array = np.maximum(X_array, 0)
        adata_subset.X = X_array

    # Keep size factors aligned to subset observations.
    adata_subset.obs["size_factors"] = (
        adata.obs["size_factors"]
        .reindex(adata_subset.obs.index)
        .fillna(1.0)
        .astype(float)
    )

    # Setup CellAssign on subset data only
    CellAssign.setup_anndata(adata_subset, size_factor_key="size_factors")

    # Train CellAssign model
    try:
        model = CellAssign(
            adata_subset, marker_gene_matrix, n_hidden=params.cellassign_n_hidden
        )
    except TypeError as exc:
        if "n_hidden" not in str(exc):
            raise
        await ctx.warning(
            "Installed scvi-tools CellAssign does not accept n_hidden; "
            "using its default architecture."
        )
        model = CellAssign(adata_subset, marker_gene_matrix)

    model.train(
        max_epochs=params.cellassign_max_iter, lr=params.cellassign_learning_rate
    )

    # Get predictions
    predictions = model.predict()

    # Handle different prediction formats (key provided by caller)
    per_cell_confidence = np.zeros(adata.n_obs, dtype=float)
    confidence_scores: dict[str, float] = {}

    if isinstance(predictions, pd.DataFrame):
        predicted_labels = predictions.idxmax(axis=1).astype(str).tolist()
        unknown_labels = sorted(set(predicted_labels) - set(valid_cell_types))
        if unknown_labels:
            raise DataError(
                "CellAssign returned labels that are not present in the marker matrix: "
                + ", ".join(unknown_labels)
            )
        adata.obs[output_key] = predicted_labels

        for i, label in enumerate(predicted_labels):
            per_cell_confidence[i] = round(float(predictions.iloc[i][label]), 3)

        for cell_type in valid_cell_types:
            mask = np.asarray(predicted_labels) == cell_type
            if mask.any():
                confidence_scores[cell_type] = round(
                    float(per_cell_confidence[mask].mean()), 2
                )
    else:
        # Other models return indices directly — no probability data
        adata.obs[output_key] = [valid_cell_types[i] for i in predictions]

    ensure_categorical(adata, output_key)

    # Store per-cell confidence if available
    has_confidence = per_cell_confidence.sum() > 0
    if has_confidence:
        adata.obs[confidence_key] = per_cell_confidence

    # Get cell types and counts
    counts = adata.obs[output_key].value_counts().to_dict()

    return AnnotationMethodOutput(
        cell_types=valid_cell_types,
        counts=counts,
        confidence=confidence_scores if has_confidence else {},
    )


def _build_annotation_suffix(
    method: str,
    reference_data_id: str | None,
) -> str:
    """Build a parametric suffix for annotation output keys.

    Reference-sensitive methods (tangram, scanvi, singler) encode the
    reference_data_id so that runs with different references coexist.
    Non-reference methods use the method name alone.
    """
    if reference_data_id and method in ("tangram", "scanvi", "singler"):
        return f"{method}_{reference_data_id}"
    return method


async def annotate_cell_types(
    data_id: str,
    ctx: ToolContext,
    params: AnnotationParameters,  # No default - must be provided by caller (LLM)
) -> AnnotationResult:
    """Annotate cell types in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        ctx: Tool context for data access and logging
        params: Annotation parameters

    Returns:
        Annotation result
    """
    # Retrieve the AnnData object via ToolContext
    adata = await ctx.get_adata(data_id)

    # Validate method first - clean and simple
    if params.method not in SUPPORTED_METHODS:
        raise ParameterError(
            f"Unsupported method: {params.method}. Supported: {sorted(SUPPORTED_METHODS)}"
        )

    # Get reference data if needed for methods that require it
    reference_adata = None
    if params.method in ["tangram", "scanvi", "singler"] and params.reference_data_id:
        reference_adata = await ctx.get_adata(params.reference_data_id)

    # Generate output keys in ONE place (single-point control)
    # Suffix encodes reference_data_id for reference-sensitive methods
    # so that runs with different references coexist in obs/uns.
    suffix = _build_annotation_suffix(params.method, params.reference_data_id)
    output_key = f"cell_type_{suffix}"
    confidence_key = f"confidence_{suffix}"

    # Tangram-specific obsm key (parametrized by suffix)
    tangram_ct_pred_key = f"tangram_ct_pred_{suffix}"

    # Route to appropriate annotation method
    try:
        if params.method == "tangram":
            result = await _annotate_with_tangram(
                adata,
                params,
                ctx,
                output_key,
                confidence_key,
                reference_adata,
                tangram_ct_pred_key=tangram_ct_pred_key,
            )
        elif params.method == "scanvi":
            result = await _annotate_with_scanvi(
                adata, params, ctx, output_key, confidence_key, reference_adata
            )
        elif params.method == "cellassign":
            result = await _annotate_with_cellassign(
                adata, params, ctx, output_key, confidence_key
            )
        elif params.method == "mllmcelltype":
            result = await _annotate_with_mllmcelltype(
                adata, params, ctx, output_key, confidence_key
            )
        elif params.method == "singler":
            result = await _annotate_with_singler(
                adata, params, ctx, output_key, confidence_key, reference_adata
            )
        else:  # sctype
            result = await _annotate_with_sctype(
                adata, params, ctx, output_key, confidence_key
            )

    except (ParameterError, DataError, DataNotFoundError, ImportError):
        raise
    except Exception as e:
        raise ProcessingError(f"Annotation failed: {e}") from e

    # Extract values from unified result type
    cell_types = result.cell_types
    counts = result.counts
    confidence_scores = result.confidence
    tangram_mapping_score = result.tangram_mapping_score

    # Determine if confidence_key should be reported (only if we have confidence data)
    confidence_key_for_result = confidence_key if confidence_scores else None

    # Store scientific metadata for reproducibility
    from ..utils.adata_utils import store_analysis_metadata
    from ..utils.results_export import export_analysis_result

    # Store per-cell-type counts summary in uns for export
    counts_key = f"{output_key}_counts"
    if counts:
        adata.uns[counts_key] = dict(counts)

    # Extract results keys
    results_keys_dict: dict[str, list[str]] = {
        "obs": [output_key],
        "obsm": [],
        "uns": [],
    }
    if counts:
        results_keys_dict["uns"].append(counts_key)
    if confidence_key_for_result:
        results_keys_dict["obs"].append(confidence_key)

    # Add method-specific result keys
    # Note: tangram_gene_predictions (n_cells × n_genes) is too large for CSV export
    # Only export tangram_ct_pred (n_cells × n_cell_types) which is reasonably sized
    if params.method == "tangram":
        results_keys_dict["obsm"].append(tangram_ct_pred_key)

    # Prepare parameters dict (only scientifically important ones)
    parameters_dict = {}
    if params.method == "tangram":
        parameters_dict = {
            "use_gpu": params.tangram_use_gpu,
            "n_epochs": params.n_epochs,
            "learning_rate": params.tangram_learning_rate,
        }
    elif params.method == "scanvi":
        parameters_dict = {
            "n_latent": params.scanvi_n_latent,
            "n_hidden": params.scanvi_n_hidden,
            "n_layers": params.scanvi_n_layers,
            "dropout_rate": params.scanvi_dropout_rate,
            "use_scvi_pretrain": params.scanvi_use_scvi_pretrain,
            "scvi_epochs": params.scanvi_scvi_epochs,
            "scanvi_epochs": params.scanvi_scanvi_epochs,
            "query_epochs": params.scanvi_query_epochs,
            "n_samples_per_label": params.scanvi_n_samples_per_label,
        }
    elif params.method == "mllmcelltype":
        parameters_dict = {
            "n_marker_genes": params.mllm_n_marker_genes,
            "species": params.mllm_species,
            "provider": params.mllm_provider,
            "model": params.mllm_model,
            "use_consensus": params.mllm_use_consensus,
        }
    elif params.method == "sctype":
        parameters_dict = {
            "tissue": params.sctype_tissue,
            "scaled": params.sctype_scaled,
            "custom_markers_provided": params.sctype_custom_markers is not None,
            "allow_remote": params.sctype_allow_remote,
            "allow_runtime_r_install": params.sctype_allow_runtime_r_install,
        }
        if params.sctype_db_:
            parameters_dict["db"] = params.sctype_db_
    elif params.method == "singler":
        parameters_dict = {
            "fine_tune": params.singler_fine_tune,
        }

    # Prepare statistics dict (heterogeneous value types)
    statistics_dict: dict[str, int | float] = {"n_cell_types": len(cell_types)}
    if tangram_mapping_score is not None:
        statistics_dict["mapping_score"] = tangram_mapping_score

    # Prepare reference info if applicable
    reference_info_dict = None
    if params.method in ["tangram", "scanvi", "singler"] and params.reference_data_id:
        reference_info_dict = {"reference_data_id": params.reference_data_id}

    # Store metadata (keyed by suffix to allow coexistence of different references)
    analysis_key = f"annotation_{suffix}"
    store_analysis_metadata(
        adata,
        analysis_name=analysis_key,
        method=params.method,
        parameters=parameters_dict,
        results_keys=results_keys_dict,
        statistics=statistics_dict,
        reference_info=reference_info_dict,
    )

    # Export results for reproducibility
    export_analysis_result(adata, data_id, analysis_key)

    # Return result
    return AnnotationResult(
        data_id=data_id,
        method=params.method,
        output_key=output_key,
        confidence_key=confidence_key_for_result,
        cell_types=cell_types,
        counts=counts,
        confidence_scores=confidence_scores,
        tangram_mapping_score=tangram_mapping_score,
    )


# ============================================================================
# SC-TYPE IMPLEMENTATION
# ============================================================================

# Cache for sc-type results (bounded LRU + optional TTL)
_SCTYPE_CACHE: OrderedDict[str, Any] = OrderedDict()
_SCTYPE_CACHE_DIR = Path.home() / ".chatspatial" / "sctype_cache"
_SCTYPE_CACHE_MAX_ITEMS = max(
    1, int(os.getenv("CHATSPATIAL_SCTYPE_CACHE_MAX_ITEMS", "64"))
)
_SCTYPE_CACHE_TTL_SECONDS = max(
    0, int(os.getenv("CHATSPATIAL_SCTYPE_CACHE_TTL_SECONDS", "86400"))
)

# R code constants for sc-type (extracted for clarity)
_R_CHECK_PACKAGES = """
required_packages <- c("dplyr", "openxlsx", "HGNChelper")
missing_packages <- required_packages[!sapply(required_packages, require, character.only = TRUE, quietly = TRUE)]
if (length(missing_packages) > 0) {
    stop(
        paste(
            "Missing required R packages for scType:",
            paste(missing_packages, collapse = ", "),
            ". Install manually or set CHATSPATIAL_ALLOW_RUNTIME_R_INSTALL=1."
        )
    )
}
"""

_R_INSTALL_PACKAGES = """
required_packages <- c("dplyr", "openxlsx", "HGNChelper")
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cran.r-project.org/", quiet = TRUE)
        if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
            stop(paste("Failed to install R package:", pkg))
        }
    }
}
"""

_R_LOAD_SCTYPE_REMOTE = """
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
"""

_SCTYPE_DEFAULT_DB_URL = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

_R_SCTYPE_SCORING = """
# Set row/column names and convert to dense
rownames(scdata) <- gene_names
colnames(scdata) <- cell_names
if (inherits(scdata, 'sparseMatrix')) scdata <- as.matrix(scdata)

# Extract gene sets
gs_positive <- gs_list$gs_positive
gs_negative <- gs_list$gs_negative

if (length(gs_positive) == 0) stop("No valid positive gene sets found")

# Filter gene sets to genes present in data
available_genes <- rownames(scdata)
filtered_gs_positive <- list()
filtered_gs_negative <- list()

for (celltype in names(gs_positive)) {
    pos_genes <- gs_positive[[celltype]]
    neg_genes <- if (celltype %in% names(gs_negative)) gs_negative[[celltype]] else c()
    pos_overlap <- intersect(toupper(pos_genes), toupper(available_genes))
    if (length(pos_overlap) > 0) {
        filtered_gs_positive[[celltype]] <- pos_overlap
        filtered_gs_negative[[celltype]] <- intersect(toupper(neg_genes), toupper(available_genes))
    }
}

if (length(filtered_gs_positive) == 0) {
    stop("No valid cell type gene sets found after filtering.")
}

# Run sc-type scoring
es_max <- sctype_score(
    scRNAseqData = as.matrix(scdata),
    scaled = TRUE,
    gs = filtered_gs_positive,
    gs2 = filtered_gs_negative
)

if (is.null(es_max) || nrow(es_max) == 0) {
    stop("SC-Type scoring failed to produce results.")
}
"""

# Valid tissue types from sc-type database
SCTYPE_VALID_TISSUES = {
    "Adrenal",
    "Brain",
    "Eye",
    "Heart",
    "Hippocampus",
    "Immune system",
    "Intestine",
    "Kidney",
    "Liver",
    "Lung",
    "Muscle",
    "Pancreas",
    "Placenta",
    "Spleen",
    "Stomach",
    "Thymus",
}


def _is_truthy_env(name: str, default: str = "0") -> bool:
    """Interpret standard truthy environment variable values."""
    return os.getenv(name, default).strip().lower() in {"1", "true", "yes", "on"}


def _is_allowed_by_call_or_env(allow_per_call: bool, env_name: str) -> bool:
    return allow_per_call or _is_truthy_env(env_name, "0")


def _is_remote_resource(path: str) -> bool:
    """Whether a resource path points to an HTTP(S) URL."""
    return path.lower().startswith(("http://", "https://"))


def _get_sctype_cache_key(adata, params: AnnotationParameters) -> str:
    """Generate cache key for sc-type results."""
    data_hash = hashlib.md5()

    # Include full shape so any cell/gene count change invalidates
    data_hash.update(f"{adata.n_obs}x{adata.n_vars}".encode())

    # Hash the entire expression matrix for correctness
    X = adata.X
    if sparse.issparse(X):
        csr = X.tocsr()
        data_hash.update(csr.data.tobytes())
        data_hash.update(csr.indices.tobytes())
        data_hash.update(csr.indptr.tobytes())
    else:
        arr = np.asarray(X)
        data_hash.update(np.ascontiguousarray(arr).tobytes())

    # Hash relevant parameters
    params_dict = {
        "tissue": params.sctype_tissue,
        "db": params.sctype_db_,
        "scaled": params.sctype_scaled,
        "custom_markers": params.sctype_custom_markers,
    }
    data_hash.update(json.dumps(params_dict, sort_keys=True, default=str).encode())

    return data_hash.hexdigest()


def _debug_cache(ctx: "ToolContext", msg: str) -> None:
    """Emit cache diagnostics without surfacing noisy user-facing warnings."""
    debug_fn = getattr(ctx, "debug", None)
    if callable(debug_fn):
        debug_fn(msg)


def _pack_cache_entry(
    payload: tuple[Any, Any, Any, Any],
) -> tuple[float, tuple[Any, ...]]:
    """Pack cache payload with timestamp."""
    return (time.time(), payload)


def _unpack_cache_entry(entry: Any) -> Optional[tuple[float, tuple[Any, ...]]]:
    """Unpack modern/legacy cache entry formats.

    Modern format: (timestamp, payload)
    Legacy format: payload only
    """
    if (
        isinstance(entry, tuple)
        and len(entry) == 2
        and isinstance(entry[0], (int, float))
        and isinstance(entry[1], tuple)
        and len(entry[1]) == 4
    ):
        return float(entry[0]), entry[1]

    if isinstance(entry, tuple) and len(entry) == 4:
        # Legacy in-memory entry had no timestamp.
        return time.time(), entry

    return None


def _is_cache_expired(cached_at: float) -> bool:
    """Check TTL-based expiration."""
    if _SCTYPE_CACHE_TTL_SECONDS <= 0:
        return False
    return (time.time() - cached_at) > _SCTYPE_CACHE_TTL_SECONDS


def _evict_sctype_cache_if_needed() -> None:
    """Enforce in-memory cache size bound."""
    if _SCTYPE_CACHE_MAX_ITEMS <= 0:
        _SCTYPE_CACHE.clear()
        return
    while len(_SCTYPE_CACHE) > _SCTYPE_CACHE_MAX_ITEMS:
        popitem_fn = getattr(_SCTYPE_CACHE, "popitem", None)
        if callable(popitem_fn):
            try:
                popitem_fn(last=False)
                continue
            except TypeError:
                pass
        oldest_key = next(iter(_SCTYPE_CACHE))
        _SCTYPE_CACHE.pop(oldest_key, None)


def _touch_sctype_cache_key(cache_key: str) -> None:
    """Mark key as recently used when OrderedDict semantics are available."""
    move_to_end_fn = getattr(_SCTYPE_CACHE, "move_to_end", None)
    if callable(move_to_end_fn):
        move_to_end_fn(cache_key)


def _store_memory_sctype_cache(
    cache_key: str, payload: tuple[Any, Any, Any, Any]
) -> None:
    """Insert entry into bounded LRU cache."""
    _SCTYPE_CACHE[cache_key] = _pack_cache_entry(payload)
    _touch_sctype_cache_key(cache_key)
    _evict_sctype_cache_if_needed()


def _load_sctype_functions(
    ctx: "ToolContext",
    *,
    allow_remote: bool = False,
    allow_runtime_install: bool = False,
) -> None:
    """Load scType R functions with explicit supply-chain controls.

    Security policy (default-safe):
    - Runtime R package auto-install is OFF unless allowed per call or by env
    - Remote script sourcing is OFF unless allowed per call or by env
    - Prefer local script directory via CHATSPATIAL_SCTYPE_R_DIR
    """
    allow_runtime_install = _is_allowed_by_call_or_env(
        allow_runtime_install, "CHATSPATIAL_ALLOW_RUNTIME_R_INSTALL"
    )
    allow_remote_source = _is_allowed_by_call_or_env(
        allow_remote, "CHATSPATIAL_ALLOW_REMOTE_R_SOURCE"
    )
    local_script_dir = os.getenv("CHATSPATIAL_SCTYPE_R_DIR")

    if local_script_dir:
        script_dir = Path(local_script_dir).expanduser()
        gene_sets_prepare = script_dir / "gene_sets_prepare.R"
        sctype_score = script_dir / "sctype_score_.R"
        missing = [str(p) for p in (gene_sets_prepare, sctype_score) if not p.exists()]
        if missing:
            raise ParameterError(
                "Local scType R script directory is configured but files are missing: "
                f"{missing}"
            )
        load_script = (
            f'source("{gene_sets_prepare.as_posix()}")\n'
            f'source("{sctype_score.as_posix()}")'
        )
    else:
        if not allow_remote_source:
            raise ParameterError(
                "scType remote R script sourcing is disabled by default. "
                "For one-off runs, pass sctype_allow_remote=true. "
                "For persistent local use, set CHATSPATIAL_SCTYPE_R_DIR to a local script directory. "
                "For persistent remote use, set CHATSPATIAL_ALLOW_REMOTE_R_SOURCE=1."
            )
        load_script = _R_LOAD_SCTYPE_REMOTE

    robjects, _, _, _, _, default_converter, openrlib, _ = validate_r_environment(ctx)
    from rpy2.robjects import conversion

    with openrlib.rlock:
        with conversion.localconverter(default_converter):
            if allow_runtime_install:
                robjects.r(_R_INSTALL_PACKAGES)
            else:
                robjects.r(_R_CHECK_PACKAGES)
            robjects.r(load_script)


def _prepare_sctype_genesets(
    params: AnnotationParameters,
    ctx: "ToolContext",
    *,
    allow_remote: bool = False,
) -> Any:
    """Prepare gene sets for sc-type."""
    if params.sctype_custom_markers:
        return _convert_custom_markers_to_gs(params.sctype_custom_markers, ctx)

    # Use sc-type database
    tissue = params.sctype_tissue
    if not tissue:
        raise ParameterError("sctype_tissue is required when not using custom markers")

    db_path = params.sctype_db_ or _SCTYPE_DEFAULT_DB_URL
    if _is_remote_resource(db_path):
        allow_remote_db = _is_allowed_by_call_or_env(
            allow_remote, "CHATSPATIAL_ALLOW_REMOTE_SCTYPE_DB"
        )
        if not allow_remote_db:
            raise ParameterError(
                "scType remote database download is disabled by default. "
                "For one-off runs, pass sctype_allow_remote=true. "
                "For persistent local use, provide a local sctype_db_ path. "
                "For persistent remote use, set CHATSPATIAL_ALLOW_REMOTE_SCTYPE_DB=1."
            )
    else:
        local_db = Path(db_path).expanduser()
        if not local_db.exists():
            raise DataNotFoundError(f"scType database file not found: {local_db}")
        db_path = local_db.as_posix()

    robjects, _, _, _, _, default_converter, openrlib, _ = validate_r_environment(ctx)
    from rpy2.robjects import conversion

    with openrlib.rlock:
        with conversion.localconverter(default_converter):
            robjects.r.assign("db_path", db_path)
            robjects.r.assign("tissue_type", tissue)
            robjects.r("gs_list <- gene_sets_prepare(db_path, tissue_type)")
            return robjects.r["gs_list"]


def _convert_custom_markers_to_gs(
    custom_markers: dict[str, dict[str, list[str]]], ctx: "ToolContext"
) -> Any:
    """Convert custom markers to sc-type gene set format"""
    if not custom_markers:
        raise DataError("Custom markers dictionary is empty")

    gs_positive = {}
    gs_negative = {}

    valid_celltypes = 0

    for cell_type, markers in custom_markers.items():
        if not isinstance(markers, dict):
            continue

        positive_genes = []
        negative_genes = []

        if "positive" in markers and isinstance(markers["positive"], list):
            positive_genes = [
                str(g).strip().upper()
                for g in markers["positive"]
                if g and str(g).strip()
            ]

        if "negative" in markers and isinstance(markers["negative"], list):
            negative_genes = [
                str(g).strip().upper()
                for g in markers["negative"]
                if g and str(g).strip()
            ]

        # Only include cell types that have at least some positive markers
        if positive_genes:
            gs_positive[cell_type] = positive_genes
            gs_negative[cell_type] = negative_genes  # Can be empty list
            valid_celltypes += 1

    if valid_celltypes == 0:
        raise DataError(
            "No valid cell types found in custom markers - all cell types need at least one positive marker"
        )

    # Get robjects and converters from validation
    robjects, pandas2ri, _, _, localconverter, default_converter, openrlib, _ = (
        validate_r_environment(ctx)
    )

    # Wrap R calls in conversion context (FIX for contextvars issue)
    with openrlib.rlock:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Convert Python dictionaries to R named lists, handle empty lists properly
            r_gs_positive = robjects.r["list"](
                **{
                    k: robjects.StrVector(v) if v else robjects.StrVector([])
                    for k, v in gs_positive.items()
                }
            )
            r_gs_negative = robjects.r["list"](
                **{
                    k: robjects.StrVector(v) if v else robjects.StrVector([])
                    for k, v in gs_negative.items()
                }
            )

            # Create the final gs_list structure
            gs_list = robjects.r["list"](
                gs_positive=r_gs_positive, gs_negative=r_gs_negative
            )

    return gs_list


def _run_sctype_scoring(
    adata, gs_list, params: AnnotationParameters, ctx: "ToolContext"
) -> pd.DataFrame:
    """Run sc-type scoring algorithm."""
    robjects, pandas2ri, numpy2ri, _, _, default_converter, openrlib, anndata2ri = (
        validate_r_environment(ctx)
    )
    from rpy2.robjects import conversion

    # Prepare expression data
    expr_data = (
        adata.layers["scaled"]
        if params.sctype_scaled and "scaled" in adata.layers
        else adata.X
    )

    with openrlib.rlock:
        with conversion.localconverter(
            default_converter
            + anndata2ri.converter
            + pandas2ri.converter
            + numpy2ri.converter
        ):
            # Transfer data to R (genes × cells for scType)
            robjects.r.assign("scdata", expr_data.T)
            robjects.r.assign("gene_names", list(adata.var_names))
            robjects.r.assign("cell_names", list(adata.obs_names))
            robjects.r.assign("gs_list", gs_list)

            # Run scoring using pre-defined R code
            robjects.r(_R_SCTYPE_SCORING)

            # Get results
            row_names = list(robjects.r("rownames(es_max)"))
            col_names = list(robjects.r("colnames(es_max)"))
            scores_matrix = robjects.r["es_max"]

    # Convert to DataFrame
    if isinstance(scores_matrix, pd.DataFrame):
        scores_df = scores_matrix
        scores_df.index = row_names if row_names else scores_df.index
        scores_df.columns = col_names if col_names else scores_df.columns
    else:
        scores_df = pd.DataFrame(scores_matrix, index=row_names, columns=col_names)

    return scores_df


def _softmax(scores_array: np.ndarray) -> np.ndarray:
    """Compute softmax probabilities from raw scores (numerically stable).

    NaN values are treated as -inf (zero probability) so they don't
    propagate through the entire output.
    """
    clean = np.where(
        np.isnan(scores_array) | np.isinf(scores_array), -np.inf, scores_array
    )
    max_val = clean.max()
    if np.isinf(max_val) and max_val < 0:
        return np.zeros_like(scores_array, dtype=float)
    shifted = clean - max_val
    exp_scores = np.exp(shifted)
    total = np.sum(exp_scores)
    if total == 0:
        # All entries are -inf / NaN → uniform zero
        return np.zeros_like(scores_array, dtype=float)
    return exp_scores / total


def _assign_sctype_celltypes(
    scores_df: pd.DataFrame, ctx: "ToolContext"
) -> tuple[list[str], list[float]]:
    """Assign cell types based on sc-type scores using softmax confidence."""
    if scores_df is None or scores_df.empty:
        raise DataError("Scores DataFrame is empty or None")

    cell_types = []
    confidence_scores = []

    for col_name in scores_df.columns:
        cell_scores = scores_df[col_name]

        # Guard against all-NaN columns (bad marker data)
        if cell_scores.isna().all():
            cell_types.append("Unknown")
            confidence_scores.append(0.0)
            continue

        max_idx = cell_scores.idxmax()
        max_score = cell_scores.loc[max_idx]

        if max_score > 0:
            cell_types.append(str(max_idx))
            # Softmax gives statistically meaningful confidence
            softmax_probs = _softmax(cell_scores.values)
            confidence_scores.append(
                float(softmax_probs[cell_scores.index.get_loc(max_idx)])
            )
        else:
            cell_types.append("Unknown")
            confidence_scores.append(0.0)

    return cell_types, confidence_scores


def _calculate_sctype_stats(cell_types: list[str]) -> dict[str, int]:
    """Calculate cell type counts."""
    from collections import Counter

    return dict(Counter(cell_types))


async def _cache_sctype_results(
    cache_key: str, results: tuple, ctx: "ToolContext"
) -> None:
    """Cache sc-type results to disk as JSON (secure, no pickle)."""
    try:
        _SCTYPE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_file = _SCTYPE_CACHE_DIR / f"{cache_key}.json"

        # Convert tuple to serializable dict
        per_cell_types, counts, confidence_by_celltype, mapping_score = results
        cache_data = {
            "cache_format_version": 2,
            "cached_at": time.time(),
            "per_cell_types": per_cell_types,
            # Backward-compatible alias for older cache readers.
            "cell_types": per_cell_types,
            "counts": counts,
            "confidence_by_celltype": confidence_by_celltype,
            "mapping_score": mapping_score,
        }

        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(cache_data, f)

        _store_memory_sctype_cache(cache_key, results)
    except Exception as e:
        _debug_cache(ctx, f"sc-type cache write skipped: {e}")


def _load_cached_sctype_results(cache_key: str, ctx: "ToolContext") -> Optional[tuple]:
    """Load cached sc-type results from memory or JSON file."""
    # Check memory cache first
    if cache_key in _SCTYPE_CACHE:
        unpacked = _unpack_cache_entry(_SCTYPE_CACHE[cache_key])
        if unpacked is None:
            _SCTYPE_CACHE.pop(cache_key, None)
            _debug_cache(ctx, "sc-type cache entry dropped: invalid in-memory format")
        else:
            cached_at, payload = unpacked
            if _is_cache_expired(cached_at):
                _SCTYPE_CACHE.pop(cache_key, None)
                _debug_cache(ctx, "sc-type cache miss: in-memory entry expired")
            else:
                _touch_sctype_cache_key(cache_key)
                return payload

    # Check disk cache (JSON)
    cache_file = _SCTYPE_CACHE_DIR / f"{cache_key}.json"
    if cache_file.exists():
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                cache_data = json.load(f)

            cached_at = float(cache_data.get("cached_at", time.time()))
            if _is_cache_expired(cached_at):
                cache_file.unlink(missing_ok=True)
                _debug_cache(ctx, "sc-type cache miss: disk entry expired")
                return None

            per_cell_types = cache_data.get(
                "per_cell_types", cache_data.get("cell_types")
            )
            if not isinstance(per_cell_types, list):
                raise ValueError(
                    "Invalid cache payload: missing per-cell sc-type labels"
                )

            results = (
                per_cell_types,
                cache_data["counts"],
                cache_data["confidence_by_celltype"],
                cache_data.get("mapping_score"),
            )
            _store_memory_sctype_cache(cache_key, results)
            return results
        except Exception as e:
            # Cache corrupted or incompatible, will recompute
            _debug_cache(ctx, f"sc-type cache read skipped: {e}")
            try:
                cache_file.unlink(missing_ok=True)
            except Exception as cleanup_error:
                _debug_cache(
                    ctx,
                    f"sc-type cache cleanup skipped for {cache_file.name}: {cleanup_error}",
                )

    return None


async def _annotate_with_sctype(
    adata: sc.AnnData,
    params: AnnotationParameters,
    ctx: "ToolContext",
    output_key: str,
    confidence_key: str,
) -> AnnotationMethodOutput:
    """Annotate cell types using sc-type method."""
    # Validate parameters
    if not params.sctype_tissue and not params.sctype_custom_markers:
        raise ParameterError(
            "Either sctype_tissue or sctype_custom_markers must be specified"
        )

    if params.sctype_tissue and params.sctype_tissue not in SCTYPE_VALID_TISSUES:
        raise ParameterError(
            f"Tissue '{params.sctype_tissue}' not supported. "
            f"Valid: {', '.join(sorted(SCTYPE_VALID_TISSUES))}"
        )

    # Check cache
    cache_key = None
    if params.sctype_use_cache:
        cache_key = _get_sctype_cache_key(adata, params)
        cached = _load_cached_sctype_results(cache_key, ctx)
        if cached:
            per_cell_types, counts, confidence_by_celltype, per_cell_conf = cached
            if len(per_cell_types) == adata.n_obs:
                adata.obs[output_key] = pd.Categorical(per_cell_types)
                # Use per-cell confidence if available (new format),
                # fall back to per-type mean for legacy cache entries.
                # Legacy caches may store a scalar mapping_score instead of
                # a per-cell list; guard with isinstance before len().
                if (
                    isinstance(per_cell_conf, list)
                    and len(per_cell_conf) == adata.n_obs
                ):
                    adata.obs[confidence_key] = per_cell_conf
                else:
                    adata.obs[confidence_key] = [
                        float(confidence_by_celltype.get(ct, 0.0))
                        for ct in per_cell_types
                    ]
                unique_cell_types = list(dict.fromkeys(per_cell_types))
                return AnnotationMethodOutput(
                    cell_types=unique_cell_types,
                    counts=counts,
                    confidence=confidence_by_celltype,
                )

            await ctx.warning(
                "sc-type cache entry is stale (cell count mismatch). Recomputing results."
            )

    # Run sc-type pipeline
    _load_sctype_functions(
        ctx,
        allow_remote=params.sctype_allow_remote,
        allow_runtime_install=params.sctype_allow_runtime_r_install,
    )
    gs_list = _prepare_sctype_genesets(
        params,
        ctx,
        allow_remote=params.sctype_allow_remote,
    )
    scores_df = _run_sctype_scoring(adata, gs_list, params, ctx)
    per_cell_types, per_cell_confidence = _assign_sctype_celltypes(scores_df, ctx)

    # Calculate statistics
    counts = _calculate_sctype_stats(per_cell_types)

    # Average confidence per cell type in O(n) time.
    confidence_sums: dict[str, float] = {}
    confidence_counts: dict[str, int] = {}
    for ct, conf in zip(per_cell_types, per_cell_confidence, strict=False):
        confidence_sums[ct] = confidence_sums.get(ct, 0.0) + float(conf)
        confidence_counts[ct] = confidence_counts.get(ct, 0) + 1
    confidence_by_celltype = {
        ct: (
            (confidence_sums[ct] / confidence_counts[ct])
            if confidence_counts[ct]
            else 0.0
        )
        for ct in dict.fromkeys(per_cell_types)
    }

    # Store in adata.obs (keys provided by caller)
    adata.obs[output_key] = pd.Categorical(per_cell_types)
    adata.obs[confidence_key] = per_cell_confidence

    # Preserve deterministic cell-type order based on first appearance.
    unique_cell_types = list(dict.fromkeys(per_cell_types))

    # Cache results: store per-cell confidence for consistent restore
    if params.sctype_use_cache and cache_key:
        cache_tuple = (
            per_cell_types,
            counts,
            confidence_by_celltype,
            per_cell_confidence,
        )
        await _cache_sctype_results(cache_key, cache_tuple, ctx)

    return AnnotationMethodOutput(
        cell_types=unique_cell_types,
        counts=counts,
        confidence=confidence_by_celltype,
    )
