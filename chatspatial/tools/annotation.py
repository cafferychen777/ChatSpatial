"""
Cell type annotation tools for spatial transcriptomics data.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext

from ..models.analysis import AnnotationResult
from ..models.data import AnnotationParameters
from ..utils.adata_utils import (
    ensure_categorical,
    ensure_counts_layer,
    ensure_unique_var_names_with_ctx,
    find_common_genes,
    get_spatial_key,
    to_dense,
    validate_obs_column,
)
from ..utils.dependency_manager import (
    is_available,
    require,
    validate_r_environment,
    validate_scvi_tools,
)
from ..utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)

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
    adata,
    params: AnnotationParameters,
    ctx: "ToolContext",
    reference_adata: Optional[Any] = None,
):
    """Annotate cell types using SingleR reference-based method"""
    # Validate and import dependencies
    require("singler", ctx, feature="SingleR annotation")
    require("singlecellexperiment", ctx, feature="SingleR annotation")
    import singler

    # Optional: check for celldex
    celldex = None
    if is_available("celldex"):
        import celldex

    # Get expression matrix - prefer normalized data
    # IMPORTANT: Ensure test_mat dimensions match adata.var_names (used in test_features)
    if "X_normalized" in adata.layers:
        test_mat = adata.layers["X_normalized"]
    elif adata.X is not None:
        test_mat = adata.X
    else:
        # Fallback: use raw data, but subset to current var_names to ensure dimension match
        # Note: adata.raw may have full genes while adata has HVG subset
        if adata.raw is not None:
            test_mat = adata.raw[:, adata.var_names].X
        else:
            test_mat = adata.X

    # MEMORY OPTIMIZATION: SingleR (singler-py) natively supports sparse matrices
    # No toarray() needed - both np.log1p() and .T() work with sparse matrices
    # Verified: sparse and dense inputs produce identical results
    # Memory savings: ~1.3 GB for typical 10K cells × 20K genes dataset

    # Ensure log-normalization (SingleR expects log-normalized data)
    if "log1p" not in adata.uns:
        await ctx.warning(
            "Data may not be log-normalized. Applying log1p for SingleR..."
        )
        test_mat = np.log1p(test_mat)

    # Transpose for SingleR (genes x cells)
    test_mat = test_mat.T

    # Ensure gene names are strings
    test_features = [str(x) for x in adata.var_names]

    # Prepare reference
    reference_name = getattr(params, "singler_reference", None)
    reference_data_id = getattr(params, "reference_data_id", None)

    ref_data = None
    ref_labels = None
    ref_features_to_use = None  # Only set when using custom reference (not celldex)

    # Priority: reference_name > reference_data_id > default
    if reference_name and celldex:
        ref = celldex.fetch_reference(reference_name, "2024-02-26", realize_assays=True)
        # Get labels
        for label_col in ["label.main", "label.fine", "cell_type"]:
            try:
                ref_labels = ref.get_column_data().column(label_col)
                break
            except Exception:
                continue  # Try next label column
        if ref_labels is None:
            raise DataNotFoundError(
                f"Could not find labels in reference {reference_name}"
            )
        ref_data = ref

    elif reference_data_id and reference_adata is not None:
        # Use provided reference data (passed from main function via ctx.get_adata())
        # Handle duplicate gene names
        await ensure_unique_var_names_with_ctx(reference_adata, ctx, "reference data")
        if await ensure_unique_var_names_with_ctx(adata, ctx, "query data") > 0:
            # Update test_features after fixing
            test_features = [str(x) for x in adata.var_names]

        # Get reference expression matrix
        if "X_normalized" in reference_adata.layers:
            ref_mat = reference_adata.layers["X_normalized"]
        else:
            ref_mat = reference_adata.X

        # MEMORY OPTIMIZATION: SingleR (singler-py) natively supports sparse matrices
        # No toarray() needed - both np.log1p() and .T() work with sparse matrices
        # Verified: sparse and dense inputs produce identical results
        # Memory savings: ~1.3 GB for typical 10K cells × 20K genes reference dataset

        # Ensure log-normalization for reference
        if "log1p" not in reference_adata.uns:
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

        validate_obs_column(reference_adata, cell_type_key, "Cell type column")

        ref_labels = list(reference_adata.obs[cell_type_key])

        # For SingleR, pass the actual expression matrix directly (not SCE)
        # This has been shown to work better in testing
        ref_data = ref_mat
        ref_features_to_use = ref_features  # Keep reference features for gene matching

    elif celldex:
        # Use default reference
        ref = celldex.fetch_reference(
            "blueprint_encode", "2024-02-26", realize_assays=True
        )
        ref_labels = ref.get_column_data().column("label.main")
        ref_data = ref
    else:
        raise DataNotFoundError(
            "No reference data. Provide reference_data_id or singler_reference."
        )

    # Run SingleR annotation
    use_integrated = getattr(params, "singler_integrated", False)
    num_threads = getattr(params, "num_threads", 4)

    if use_integrated and isinstance(ref_data, list):
        single_results, integrated = singler.annotate_integrated(
            test_mat,
            ref_data=ref_data,
            ref_labels=ref_labels,
            test_features=test_features,
            num_threads=num_threads,
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
            "num_threads": num_threads,
        }

        # Add ref_features if we're using custom reference data (not celldex)
        if ref_features_to_use is not None:
            annotate_kwargs["ref_features"] = ref_features_to_use

        results = singler.annotate_single(**annotate_kwargs)
        best_labels = results.column("best")
        scores = results.column("scores")

        # Try to get delta scores for confidence (higher delta = higher confidence)
        try:
            delta_scores = results.column("delta")
            if delta_scores:
                low_delta = sum(1 for d in delta_scores if d and d < 0.05)
                if low_delta > len(delta_scores) * 0.3:
                    await ctx.warning(
                        f"{low_delta}/{len(delta_scores)} cells have low confidence scores (delta < 0.05)"
                    )
        except Exception:
            delta_scores = None

    # Process results
    cell_types = list(best_labels)
    unique_types = list(set(cell_types))
    counts = pd.Series(cell_types).value_counts().to_dict()

    # Calculate confidence scores - prefer delta scores if available
    confidence_scores = {}

    # First try to use delta scores (more meaningful confidence measure)
    # delta_scores is always defined by the try-except block above (line 429 or 437)
    if delta_scores is not None:
        try:
            for cell_type in unique_types:
                type_indices = [i for i, ct in enumerate(cell_types) if ct == cell_type]
                if type_indices:
                    type_deltas = [
                        delta_scores[i] for i in type_indices if i < len(delta_scores)
                    ]
                    if type_deltas:
                        # Higher delta = higher confidence
                        avg_delta = np.mean([d for d in type_deltas if d is not None])
                        confidence_scores[cell_type] = round(float(avg_delta), 3)
        except Exception:
            # Delta score extraction failed, will fall back to regular scores
            pass

    # Fall back to regular scores if delta not available
    if not confidence_scores and scores is not None:
        try:
            scores_df = pd.DataFrame(scores.to_dict())
        except AttributeError:
            scores_df = pd.DataFrame(
                scores.to_numpy() if hasattr(scores, "to_numpy") else scores
            )

        for cell_type in unique_types:
            mask = [ct == cell_type for ct in cell_types]
            if cell_type in scores_df.columns and any(mask):
                type_scores = scores_df.loc[mask, cell_type]
                avg_score = type_scores.mean()
                confidence = (avg_score + 1) / 2  # Convert correlation to 0-1
                confidence_scores[cell_type] = float(confidence)
            else:
                # No confidence data available - don't assign arbitrary value
                pass  # Cell type won't have confidence score
    else:
        # No confidence information available from this method
        # Don't assign misleading confidence values
        pass

    # Use method-prefixed output keys to avoid overwriting
    output_key = f"cell_type_{params.method}"
    confidence_key = f"confidence_{params.method}"

    # Add to AnnData
    adata.obs[output_key] = cell_types
    ensure_categorical(adata, output_key)

    # Only add confidence column if we have real confidence values
    if confidence_scores:
        # Use 0.0 for cells without confidence (more honest than arbitrary 0.5)
        confidence_array = [confidence_scores.get(ct, 0.0) for ct in cell_types]
        adata.obs[confidence_key] = confidence_array

    return unique_types, counts, confidence_scores, None


async def _annotate_with_tangram(
    adata,
    params: AnnotationParameters,
    ctx: "ToolContext",
    reference_adata: Optional[Any] = None,
):
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
    # Solution: Use adata.raw which preserves original gene names
    if adata.raw is not None:
        adata_sp = adata.raw.to_adata()
        # Preserve spatial coordinates from preprocessed data
        spatial_key = get_spatial_key(adata)
        if spatial_key:
            adata_sp.obsm[spatial_key] = adata.obsm[spatial_key].copy()
    else:
        adata_sp = adata
        await ctx.warning(
            "No raw data available - using preprocessed data (may have gene name mismatches)"
        )
    # =============================================================================

    # Handle duplicate gene names
    await ensure_unique_var_names_with_ctx(adata_sc_original, ctx, "reference data")
    await ensure_unique_var_names_with_ctx(adata_sp, ctx, "spatial data")

    # Determine training genes
    training_genes = params.training_genes

    if training_genes is None:
        # Use marker genes if available
        if params.marker_genes:
            # Flatten marker genes dictionary
            training_genes = []
            for genes in params.marker_genes.values():
                training_genes.extend(genes)
            training_genes = list(set(training_genes))  # Remove duplicates
        else:
            # Use highly variable genes
            if "highly_variable" not in adata_sc_original.var:
                raise DataNotFoundError(
                    "HVGs not found in reference data. Run preprocessing first."
                )
            training_genes = list(
                adata_sc_original.var_names[adata_sc_original.var.highly_variable]
            )

    # COW FIX: Create copy of reference data to avoid modifying original
    # Tangram's pp_adatas adds metadata (uns, obs) but doesn't subset genes
    adata_sc = adata_sc_original.copy()

    # Preprocess data for Tangram
    tg.pp_adatas(adata_sc, adata_sp, genes=training_genes)

    # Set mapping mode
    mode = params.tangram_mode
    cluster_label = params.cluster_label

    if mode == "clusters" and cluster_label is None:
        await ctx.warning(
            "Cluster label not provided for 'clusters' mode. Using default cell type annotation if available."
        )
        # Try to find a cell type annotation in the reference data
        for col in ["cell_type", "celltype", "cell_types", "leiden", "louvain"]:
            if col in adata_sc.obs:
                cluster_label = col
                break

        if cluster_label is None:
            raise ParameterError(
                "No cluster label found. Provide cluster_label parameter."
            )

    # Check GPU availability for device selection
    import torch

    device = params.tangram_device
    if device != "cpu" and not torch.cuda.is_available():
        await ctx.warning("GPU requested but not available - falling back to CPU")
        device = "cpu"

    # Run Tangram mapping with enhanced parameters
    mapping_args = {
        "mode": mode,
        "num_epochs": params.num_epochs,
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

        if annotation_col:
            tg.project_cell_annotations(ad_map, adata_sp, annotation=annotation_col)
        else:
            await ctx.warning("No suitable annotation column found for projection")
    except Exception as proj_error:
        await ctx.warning(f"Could not project cell annotations: {proj_error}")
        # Continue without projection

    # Get cell type predictions
    cell_types = []
    counts = {}
    confidence_scores = {}

    # Use method-prefixed output keys to avoid overwriting
    output_key = f"cell_type_{params.method}"
    confidence_key = f"confidence_{params.method}"

    if "tangram_ct_pred" in adata_sp.obsm:
        cell_type_df = adata_sp.obsm["tangram_ct_pred"]

        # Get cell types and counts
        cell_types = list(cell_type_df.columns)

        # ===== CRITICAL FIX: Row normalization for proper probability calculation =====
        # tangram_ct_pred contains unnormalized density/abundance values, NOT probabilities
        # Row sums can be != 1.0 and values can exceed 1.0
        # We normalize to convert densities → probability distributions
        cell_type_prob = cell_type_df.div(cell_type_df.sum(axis=1), axis=0)

        # Validation: Ensure normalized values are valid probabilities
        if not (cell_type_prob.values >= 0).all():
            await ctx.warning(
                "Some normalized probabilities are negative - data quality issue"
            )
        if not (cell_type_prob.values <= 1.0).all():
            await ctx.warning(
                "Some normalized probabilities exceed 1.0 - normalization failed"
            )
        if not np.allclose(cell_type_prob.sum(axis=1), 1.0):
            await ctx.warning(
                "Row sums don't equal 1.0 after normalization - numerical issue"
            )

        # Assign cell type based on highest probability (argmax is same before/after normalization)
        adata_sp.obs[output_key] = cell_type_prob.idxmax(axis=1)
        ensure_categorical(adata_sp, output_key)

        # Get counts
        counts = adata_sp.obs[output_key].value_counts().to_dict()

        # Calculate confidence scores from NORMALIZED probabilities
        confidence_scores = {}
        for cell_type in cell_types:
            cells_of_type = adata_sp.obs[output_key] == cell_type
            if np.sum(cells_of_type) > 0:
                # Use mean PROBABILITY as confidence (now guaranteed to be in [0, 1])
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

        # Copy confidence scores if they exist
        if confidence_key and confidence_key in adata_sp.obs:
            adata.obs[confidence_key] = adata_sp.obs[confidence_key]

        # Copy tangram_ct_pred from obsm
        if "tangram_ct_pred" in adata_sp.obsm:
            adata.obsm["tangram_ct_pred"] = adata_sp.obsm["tangram_ct_pred"]

        # Copy tangram_gene_predictions if they exist
        if "tangram_gene_predictions" in adata_sp.obsm:
            adata.obsm["tangram_gene_predictions"] = adata_sp.obsm[
                "tangram_gene_predictions"
            ]

    # Final completion log included in main annotate_cell_types function
    return cell_types, counts, confidence_scores, tangram_mapping_score


async def _annotate_with_scanvi(
    adata,
    params: AnnotationParameters,
    ctx: "ToolContext",
    reference_adata: Optional[Any] = None,
):
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
            - num_epochs (default: 100): SCANVI training epochs
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
        - num_epochs: 100 (SCANVI training)
        - scanvi_query_epochs: 100
        - Gene selection: 1000-10000 HVGs recommended

    Empirical Adjustments (not official):
        For small datasets (<1000 genes or <1000 cells):
        - scanvi_n_latent: 3-5 (may prevent NaN/gradient explosion)
        - scanvi_dropout_rate: 0.2-0.3 (may improve regularization)
        - scanvi_use_scvi_pretrain: False (may reduce complexity)
        - num_epochs: 50 (may prevent overfitting)
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
            num_epochs=50,                  # Prevent overfitting
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
    await ensure_unique_var_names_with_ctx(adata_ref_original, ctx, "reference data")
    await ensure_unique_var_names_with_ctx(adata, ctx, "query data")

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
        # No subsetting needed
        adata_ref = adata_ref_original.copy()
        adata_subset = adata.copy()

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
        # Setup for SCVI with labels (required for SCANVI conversion)
        # First ensure the reference has the cell type labels
        validate_obs_column(
            adata_ref, cell_type_key, "Cell type column (reference data)"
        )

        # SCVI needs to know about labels for later SCANVI conversion
        scvi.model.SCVI.setup_anndata(
            adata_ref,
            labels_key=cell_type_key,  # Important: include labels_key
            batch_key=batch_key,
            layer=params.layer,
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
            max_epochs=params.num_epochs,
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

    # Get prediction probabilities as confidence scores
    try:
        probs = spatial_model.predict(soft=True)
        confidence_scores = {}
        for i, cell_type in enumerate(cell_types):
            cells_of_type = adata_subset.obs[cell_type_key] == cell_type
            if np.sum(cells_of_type) > 0 and isinstance(probs, pd.DataFrame):
                if cell_type in probs.columns:
                    mean_prob = probs.loc[cells_of_type, cell_type].mean()
                    confidence_scores[cell_type] = round(float(mean_prob), 2)
                else:
                    # No probability column for this cell type - skip confidence
                    pass
            elif (
                np.sum(cells_of_type) > 0
                and hasattr(probs, "shape")
                and probs.shape[1] > i
            ):
                mean_prob = probs[cells_of_type, i].mean()
                confidence_scores[cell_type] = round(float(mean_prob), 2)
            # else: No cells of this type or no probability data - skip confidence
    except Exception as e:
        await ctx.warning(f"Could not get confidence scores: {e}")
        # Could not extract probabilities - return empty confidence dict
        confidence_scores = (
            {}
        )  # Empty dict clearly indicates no confidence data available

    # COW FIX: Add prediction results to original adata.obs
    # This will be picked up by the main function to store in adata.obs[output_key]
    adata.obs[cell_type_key] = adata_subset.obs[cell_type_key].values

    return cell_types, counts, confidence_scores, None


async def _annotate_with_mllmcelltype(
    adata, params: AnnotationParameters, ctx: "ToolContext"
):
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
    validate_obs_column(adata, cluster_key, "cluster_key")

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
    except Exception as e:
        raise ProcessingError(f"mLLMCellType annotation failed: {str(e)}") from e

    # Map cluster annotations back to cells
    cluster_to_celltype = {}
    for cluster_name, cell_type in annotations.items():
        # Extract cluster number from "Cluster_X" format
        cluster_id = cluster_name.replace("Cluster_", "")
        cluster_to_celltype[cluster_id] = cell_type

    # Use method-prefixed output keys to avoid overwriting
    output_key = f"cell_type_{params.method}"

    # Apply cell type annotations to cells
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
    confidence_scores = {}

    return cell_types, counts, confidence_scores, None


async def _annotate_with_cellassign(
    adata, params: AnnotationParameters, ctx: "ToolContext"
):
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

    # CRITICAL FIX: Use adata.raw for marker gene validation if available
    # Preprocessing filters genes to HVGs, but marker genes may not be in HVGs
    # adata.raw contains all original genes and should be checked first
    if adata.raw is not None:
        all_genes = set(adata.raw.var_names)
        gene_source = "adata.raw"
    else:
        all_genes = set(adata.var_names)
        gene_source = "adata.var_names"
        await ctx.warning(
            f"Using filtered gene set for marker gene validation "
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
    valid_cell_types = list(valid_marker_genes.keys())

    # Create marker gene matrix as DataFrame (required by CellAssign API)
    all_marker_genes = []
    for genes in valid_marker_genes.values():
        all_marker_genes.extend(genes)
    available_marker_genes = list(set(all_marker_genes))  # Remove duplicates

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
            if gene in available_marker_genes:
                marker_gene_matrix.loc[gene, cell_type] = 1

    # Compute size factors BEFORE subsetting (official CellAssign requirement)
    if "size_factors" not in adata.obs:
        # Calculate size factors from FULL dataset
        if hasattr(adata.X, "sum"):
            size_factors = adata.X.sum(axis=1)
            if hasattr(size_factors, "A1"):  # sparse matrix
                size_factors = size_factors.A1
        else:
            size_factors = np.sum(adata.X, axis=1)

        # Normalize and ensure positive
        size_factors = np.maximum(size_factors, 1e-6)
        mean_sf = np.mean(size_factors)
        size_factors_normalized = size_factors / mean_sf

        adata.obs["size_factors"] = pd.Series(
            size_factors_normalized, index=adata.obs.index
        )

    # Subset data to marker genes (size factors already computed)
    # Use adata.raw if available (contains all genes including markers)
    if adata.raw is not None:
        import anndata as ad_module

        adata_subset = ad_module.AnnData(
            X=adata.raw[:, available_marker_genes].X,
            obs=adata.obs.copy(),
            var=adata.raw.var.loc[available_marker_genes].copy(),
        )
    else:
        adata_subset = adata[:, available_marker_genes].copy()

    # Check for invalid values in the data
    X_array = to_dense(adata_subset.X)

    # Replace any NaN or Inf values with zeros
    if np.any(np.isnan(X_array)) or np.any(np.isinf(X_array)):
        await ctx.warning("Found NaN or Inf values in data, replacing with zeros")
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
        await ctx.warning("Found negative values in data, clipping to zero")
        X_array = np.maximum(X_array, 0)
        adata_subset.X = X_array

    # Verify size factors were transferred to subset
    if "size_factors" not in adata_subset.obs:
        raise ProcessingError(
            "Size factors not found in adata.obs. This should not happen - "
            "they should have been computed before subsetting. Please report this bug."
        )

    # Setup CellAssign on subset data only
    CellAssign.setup_anndata(adata_subset, size_factor_key="size_factors")

    # Train CellAssign model
    model = CellAssign(adata_subset, marker_gene_matrix)

    model.train(
        max_epochs=params.cellassign_max_iter, lr=params.cellassign_learning_rate
    )

    # Get predictions
    predictions = model.predict()

    # Use method-prefixed output keys to avoid overwriting
    output_key = f"cell_type_{params.method}"

    # Handle different prediction formats
    if isinstance(predictions, pd.DataFrame):
        # CellAssign returns DataFrame with probabilities
        predicted_indices = predictions.values.argmax(axis=1)
        adata.obs[output_key] = [valid_cell_types[i] for i in predicted_indices]

        # Get confidence scores from probabilities DataFrame
        confidence_scores = {}
        for i, cell_type in enumerate(valid_cell_types):
            cells_of_type = adata.obs[output_key] == cell_type
            if np.sum(cells_of_type) > 0:
                # Use iloc with boolean indexing properly
                cell_indices = np.where(cells_of_type)[0]
                mean_prob = predictions.iloc[cell_indices, i].mean()
                confidence_scores[cell_type] = round(float(mean_prob), 2)
            # else: No cells of this type - skip confidence
    else:
        # Other models return indices directly
        adata.obs[output_key] = [valid_cell_types[i] for i in predictions]
        # CellAssign returned indices, not probabilities - no confidence available
        confidence_scores = {}  # Empty dict indicates no confidence data

    ensure_categorical(adata, output_key)

    # Get cell types and counts
    cell_types = valid_cell_types
    counts = adata.obs[output_key].value_counts().to_dict()

    return cell_types, counts, confidence_scores, None


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

    # Route to appropriate annotation method
    try:
        if params.method == "tangram":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_tangram(adata, params, ctx, reference_adata)
            )
        elif params.method == "scanvi":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_scanvi(adata, params, ctx, reference_adata)
            )
        elif params.method == "cellassign":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_cellassign(adata, params, ctx)
            )
        elif params.method == "mllmcelltype":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_mllmcelltype(adata, params, ctx)
            )
        elif params.method == "singler":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_singler(adata, params, ctx, reference_adata)
            )
        else:  # sctype
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_sctype(adata, params, ctx)
            )

    except Exception as e:
        raise ProcessingError(f"Annotation failed: {str(e)}") from e

    # COW FIX: Each annotation method handles adata.obs assignment internally
    # This maintains consistency and prevents duplicate assignments
    # Construct output keys for result reporting
    output_key = f"cell_type_{params.method}"
    confidence_key = f"confidence_{params.method}" if confidence_scores else None

    # Store scientific metadata for reproducibility
    from ..utils.adata_utils import store_analysis_metadata

    # Extract results keys
    results_keys_dict = {"obs": [output_key], "obsm": [], "uns": []}
    if confidence_key:
        results_keys_dict["obs"].append(confidence_key)

    # Add method-specific result keys
    if params.method == "tangram":
        results_keys_dict["obsm"].extend(
            ["tangram_ct_pred", "tangram_gene_predictions"]
        )

    # Prepare parameters dict (only scientifically important ones)
    parameters_dict = {}
    if params.method == "tangram":
        parameters_dict = {
            "device": params.tangram_device,
            "n_epochs": params.num_epochs,  # Fixed: use num_epochs instead of tangram_num_epochs
            "learning_rate": params.tangram_learning_rate,
        }
    elif params.method == "scanvi":
        parameters_dict = {
            "n_latent": params.scanvi_n_latent,
            "n_hidden": params.scanvi_n_hidden,
            "dropout_rate": params.scanvi_dropout_rate,
            "use_scvi_pretrain": params.scanvi_use_scvi_pretrain,
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
        }
    elif params.method == "singler":
        parameters_dict = {
            "fine_tune": params.singler_fine_tune,
        }

    # Prepare statistics dict
    statistics_dict = {"n_cell_types": len(cell_types)}
    if tangram_mapping_score is not None:
        statistics_dict["mapping_score"] = tangram_mapping_score

    # Prepare reference info if applicable
    reference_info_dict = None
    if params.method in ["tangram", "scanvi", "singler"] and params.reference_data_id:
        reference_info_dict = {"reference_data_id": params.reference_data_id}

    # Store metadata
    store_analysis_metadata(
        adata,
        analysis_name=f"annotation_{params.method}",
        method=params.method,
        parameters=parameters_dict,
        results_keys=results_keys_dict,
        statistics=statistics_dict,
        reference_info=reference_info_dict,
    )

    # Return result
    return AnnotationResult(
        data_id=data_id,
        method=params.method,
        output_key=output_key,
        confidence_key=confidence_key,
        cell_types=cell_types,
        counts=counts,
        confidence_scores=confidence_scores,
        tangram_mapping_score=tangram_mapping_score,
    )


# ============================================================================
# SC-TYPE IMPLEMENTATION
# ============================================================================

# Cache for sc-type results (memory only, no pickle)
_SCTYPE_CACHE: Dict[str, Any] = {}
_SCTYPE_CACHE_DIR = Path.home() / ".chatspatial" / "sctype_cache"

# R code constants for sc-type (extracted for clarity)
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

_R_LOAD_SCTYPE = """
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
"""

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


def _get_sctype_cache_key(adata, params: AnnotationParameters) -> str:
    """Generate cache key for sc-type results"""
    # Create a hash based on data and parameters
    data_hash = hashlib.md5()

    # Hash expression data (sample first 1000 cells and 500 genes for efficiency)
    sample_slice = adata.X[: min(1000, adata.n_obs), : min(500, adata.n_vars)]
    sample_data = to_dense(sample_slice)
    data_hash.update(sample_data.tobytes())

    # Hash relevant parameters
    params_dict = {
        "tissue": params.sctype_tissue,
        "db": params.sctype_db_,
        "scaled": params.sctype_scaled,
        "custom_markers": params.sctype_custom_markers,
    }
    data_hash.update(str(params_dict).encode())

    return data_hash.hexdigest()


def _load_sctype_functions(ctx: "ToolContext") -> None:
    """Load sc-type R functions and auto-install R packages if needed."""
    robjects, _, _, _, _, default_converter, openrlib, _ = validate_r_environment(ctx)
    from rpy2.robjects import conversion

    with openrlib.rlock:
        with conversion.localconverter(default_converter):
            robjects.r(_R_INSTALL_PACKAGES)
            robjects.r(_R_LOAD_SCTYPE)


def _prepare_sctype_genesets(params: AnnotationParameters, ctx: "ToolContext"):
    """Prepare gene sets for sc-type."""
    if params.sctype_custom_markers:
        return _convert_custom_markers_to_gs(params.sctype_custom_markers, ctx)

    # Use sc-type database
    tissue = params.sctype_tissue
    if not tissue:
        raise ParameterError("sctype_tissue is required when not using custom markers")

    robjects, _, _, _, _, default_converter, openrlib, _ = validate_r_environment(ctx)
    from rpy2.robjects import conversion

    db_path = (
        params.sctype_db_
        or "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    )

    with openrlib.rlock:
        with conversion.localconverter(default_converter):
            robjects.r.assign("db_path", db_path)
            robjects.r.assign("tissue_type", tissue)
            robjects.r("gs_list <- gene_sets_prepare(db_path, tissue_type)")
            return robjects.r["gs_list"]


def _convert_custom_markers_to_gs(
    custom_markers: Dict[str, Dict[str, List[str]]], ctx: "ToolContext"
):
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
                gene.upper().strip()
                for gene in markers["positive"]
                if gene and str(gene).strip()
            ]

        if "negative" in markers and isinstance(markers["negative"], list):
            negative_genes = [
                gene.upper().strip()
                for gene in markers["negative"]
                if gene and str(gene).strip()
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
    """Compute softmax probabilities from raw scores (numerically stable)."""
    shifted = scores_array - np.max(scores_array)
    exp_scores = np.exp(shifted)
    return exp_scores / np.sum(exp_scores)


def _assign_sctype_celltypes(
    scores_df: pd.DataFrame, ctx: "ToolContext"
) -> tuple[List[str], List[float]]:
    """Assign cell types based on sc-type scores using softmax confidence."""
    if scores_df is None or scores_df.empty:
        raise DataError("Scores DataFrame is empty or None")

    cell_types = []
    confidence_scores = []

    for col_idx in range(len(scores_df.columns)):
        cell_scores = scores_df.iloc[:, col_idx]
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


def _calculate_sctype_stats(cell_types: List[str]) -> Dict[str, int]:
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
        cell_types, counts, confidence_by_celltype, mapping_score = results
        cache_data = {
            "cell_types": cell_types,
            "counts": counts,
            "confidence_by_celltype": confidence_by_celltype,
            "mapping_score": mapping_score,
        }

        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(cache_data, f)

        _SCTYPE_CACHE[cache_key] = results
    except Exception as e:
        await ctx.warning(f"Failed to cache results: {e}")


def _load_cached_sctype_results(cache_key: str, ctx: "ToolContext") -> Optional[tuple]:
    """Load cached sc-type results from memory or JSON file."""
    # Check memory cache first
    if cache_key in _SCTYPE_CACHE:
        return _SCTYPE_CACHE[cache_key]

    # Check disk cache (JSON)
    cache_file = _SCTYPE_CACHE_DIR / f"{cache_key}.json"
    if cache_file.exists():
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                cache_data = json.load(f)

            results = (
                cache_data["cell_types"],
                cache_data["counts"],
                cache_data["confidence_by_celltype"],
                cache_data.get("mapping_score"),
            )
            _SCTYPE_CACHE[cache_key] = results
            return results
        except Exception:
            # Cache corrupted or incompatible, will recompute
            pass

    return None


async def _annotate_with_sctype(
    adata: sc.AnnData, params: AnnotationParameters, ctx: "ToolContext"
) -> tuple[List[str], Dict[str, int], Dict[str, float], Optional[float]]:
    """Annotate cell types using sc-type method."""
    # Validate R environment
    validate_r_environment(ctx)

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
            return cached

    # Run sc-type pipeline
    _load_sctype_functions(ctx)
    gs_list = _prepare_sctype_genesets(params, ctx)
    scores_df = _run_sctype_scoring(adata, gs_list, params, ctx)
    cell_types, confidence_scores = _assign_sctype_celltypes(scores_df, ctx)

    # Calculate statistics
    counts = _calculate_sctype_stats(cell_types)

    # Average confidence per cell type
    confidence_by_celltype = {}
    for ct in set(cell_types):
        ct_confs = [c for i, c in enumerate(confidence_scores) if cell_types[i] == ct]
        confidence_by_celltype[ct] = sum(ct_confs) / len(ct_confs) if ct_confs else 0.0

    # Store in adata.obs
    adata.obs[f"cell_type_{params.method}"] = pd.Categorical(cell_types)
    adata.obs[f"confidence_{params.method}"] = confidence_scores

    results = (list(set(cell_types)), counts, confidence_by_celltype, None)

    # Cache results
    if params.sctype_use_cache and cache_key:
        await _cache_sctype_results(cache_key, results, ctx)

    return results
