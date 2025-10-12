"""
Cell type annotation tools for spatial transcriptomics data.
"""

import hashlib
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from mcp.server.fastmcp import Context

# Import Tangram for cell type annotation
try:
    import tangram as tg
except ImportError:
    tg = None

# Optional imports - actual validation happens at runtime
# This allows the module to load even if optional dependencies are missing

# R interface validation is now handled by _validate_rpy2_and_r() function

from ..models.analysis import AnnotationResult
from ..models.data import AnnotationParameters

# ============================================================================
# DEPENDENCY VALIDATION SYSTEM
# ============================================================================


def _validate_scvi_tools(context: Optional[Context] = None):
    """Validate scvi-tools availability and return the module"""
    try:
        import scvi
    except ImportError:
        scvi = None

    if scvi is None:
        error_msg = "scvi-tools is required for scANVI and CellAssign methods. Install with: pip install scvi-tools"
        if context:
            context.error(f"scvi-tools not available: {error_msg}")
        raise ImportError(error_msg)

    try:
        from scvi.external import CellAssign

        # Log version information if context is provided
        if context:
            try:
                version = getattr(scvi, "__version__", "unknown")
                context.info(f"Using scvi-tools version {version}")
            except Exception:
                pass

        return scvi, CellAssign
    except ImportError as e:
        error_msg = (
            f"scvi-tools is installed but CellAssign is not available. "
            f"This may be due to version incompatibility. Error: {e}"
        )
        if context:
            context.error(error_msg)
        raise ImportError(error_msg) from e


def _validate_tangram(context: Optional[Context] = None):
    """Validate tangram availability and return the module"""
    try:
        import tangram as tg

        if context and hasattr(tg, "__version__"):
            context.info(f"Using tangram version {tg.__version__}")

        return tg
    except ImportError as e:
        raise ImportError(
            "tangram-sc is required for Tangram method. Install with: pip install tangram-sc"
        ) from e


def _validate_mllmcelltype(context: Optional[Context] = None):
    """Validate mllmcelltype availability and return the module"""
    try:
        import mllmcelltype

        if context and hasattr(mllmcelltype, "__version__"):
            context.info(f"Using mllmcelltype version {mllmcelltype.__version__}")

        return mllmcelltype
    except ImportError as e:
        raise ImportError(
            "mllmcelltype is required for mLLMCellType method. Install with: pip install mllmcelltype"
        ) from e


def _validate_rpy2_and_r(context: Optional[Context] = None):
    """Validate R and rpy2 availability with detailed error reporting"""
    try:
        # First check rpy2
        import rpy2.robjects as robjects
        from rpy2.rinterface_lib import openrlib  # For thread safety
        from rpy2.robjects import conversion, default_converter, numpy2ri, pandas2ri
        from rpy2.robjects.conversion import localconverter
        from rpy2.robjects.packages import importr

        # Test R availability with proper conversion context (FIX for contextvars issue)
        # Same pattern as SPARK-X to prevent "Conversion rules missing" errors
        with openrlib.rlock:
            with conversion.localconverter(default_converter):
                robjects.r("R.version")

                if context:
                    r_version = robjects.r("R.version.string")[0]
                    context.info(f"Using R: {r_version}")

        return (
            robjects,
            pandas2ri,
            numpy2ri,
            importr,
            localconverter,
            default_converter,
            openrlib,
        )
    except ImportError as e:
        raise ImportError(
            "rpy2 is required for sc-type method. Install with: pip install rpy2 (requires R installation)"
        ) from e
    except Exception as e:
        error_msg = f"""
R environment setup failed: {str(e)}

Common solutions:
  • Install R from https://www.r-project.org/
  • Set R_HOME environment variable
  • Install required R packages: install.packages(c('dplyr', 'openxlsx', 'HGNChelper'))
  • macOS: brew install r
  • Ubuntu: sudo apt install r-base
  • Windows: Download from CRAN
"""
        raise ImportError(
            f"R environment setup failed for sc-type method: {str(e)}"
        ) from e


def _validate_singler(context: Optional[Context] = None):
    """Validate SingleR and its dependencies"""
    try:
        import singlecellexperiment as sce
        import singler

        # Optional: check for celldex
        celldex = None
        try:
            import celldex

            if context:
                context.info("celldex available for pre-built references")
        except ImportError:
            if context:
                context.info("celldex not available - will use custom references")

        if context and hasattr(singler, "__version__"):
            context.info(f"Using SingleR version {singler.__version__}")

        return singler, sce, celldex
    except ImportError as e:
        raise ImportError(
            "singler is required for SingleR method. Install with: pip install singler singlecellexperiment"
        ) from e


# Constants for annotation
DEFAULT_HVG_COUNT = 2000
# DEFAULT_SCANVI_EPOCHS is now replaced by scanvi_query_epochs parameter (default: 100)
CONFIDENCE_MAX = 0.99

# Documentation of which methods provide real confidence scores
# IMPORTANT: Methods should only return confidence scores when they represent
# real statistical measures (e.g., probabilities, correlations, overlap ratios).
# Empty confidence_scores dict indicates no confidence data available.
METHODS_WITH_REAL_CONFIDENCE = {
    "singler",  # Correlation scores from reference matching
    "tangram",  # Mapping probabilities from spatial alignment
    "sctype",  # Scoring based on marker gene expression levels
}
METHODS_WITHOUT_CONFIDENCE = {
    "mllmcelltype"  # LLM-based annotations don't provide numeric confidence
}
METHODS_WITH_PARTIAL_CONFIDENCE = {
    "scanvi",  # Provides probabilities when soft=True prediction works
    "cellassign",  # Provides probabilities when prediction returns DataFrame
}

# Supported annotation methods
SUPPORTED_METHODS = {
    "tangram",
    "scanvi",
    "cellassign",
    "mllmcelltype",
    "sctype",
    "singler",
}


async def _handle_annotation_error(
    error: Exception, method: str, context: Optional[Context] = None
) -> None:
    """Handle annotation errors consistently"""
    error_msg = f"{method} annotation failed: {str(error)}"
    if context:
        await context.error(f"Error in {method} annotation: {str(error)}")
    raise ValueError(error_msg)


async def _annotate_with_singler(
    adata,
    params: AnnotationParameters,
    data_store: Dict[str, Any] = None,
    context: Optional[Context] = None,
):
    """Annotate cell types using SingleR reference-based method"""
    try:
        singler, sce, celldex = _validate_singler(context)

        if context:
            await context.info("🔬 Starting SingleR annotation...")
            await context.info(f"   Cells: {adata.n_obs}, Genes: {adata.n_vars}")

        # Convert AnnData to SingleCellExperiment
        test_sce = sce.SingleCellExperiment.from_anndata(adata)

        # Get expression matrix - prefer normalized data
        if "X_normalized" in adata.layers:
            test_mat = adata.layers["X_normalized"]
        elif adata.X is not None:
            test_mat = adata.X
        else:
            test_mat = adata.raw.X if adata.raw else adata.X

        # Convert sparse to dense if needed
        if hasattr(test_mat, "toarray"):
            test_mat = test_mat.toarray()

        # Ensure log-normalization (SingleR expects log-normalized data)
        if "log1p" not in adata.uns:
            if context:
                await context.warning(
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

        # Priority: reference_name > reference_data_id > default
        if reference_name and celldex:
            if context:
                await context.info(f"Loading reference: {reference_name}")
            ref = celldex.fetch_reference(
                reference_name, "2024-02-26", realize_assays=True
            )
            # Get labels
            for label_col in ["label.main", "label.fine", "cell_type"]:
                try:
                    ref_labels = ref.get_column_data().column(label_col)
                    break
                except:
                    continue
            if ref_labels is None:
                raise ValueError(f"Could not find labels in reference {reference_name}")
            ref_data = ref

        elif reference_data_id and data_store and reference_data_id in data_store:
            # Use provided reference data
            if context:
                await context.info("Using provided reference data")
            reference_adata = data_store[reference_data_id]["adata"]

            # ===== Handle Duplicate Gene Names (CRITICAL FIX) =====
            if not reference_adata.var_names.is_unique:
                n_duplicates = (
                    len(reference_adata.var_names)
                    - len(set(reference_adata.var_names))
                )
                if context:
                    await context.warning(
                        f"⚠️  Found {n_duplicates} duplicate gene names in reference data, fixing..."
                    )
                reference_adata.var_names_make_unique()

            if not adata.var_names.is_unique:
                n_duplicates = len(adata.var_names) - len(set(adata.var_names))
                if context:
                    await context.warning(
                        f"⚠️  Found {n_duplicates} duplicate gene names in query data, fixing..."
                    )
                adata.var_names_make_unique()
                # Update test_features after fixing
                test_features = [str(x) for x in adata.var_names]

            # Get reference expression matrix
            if "X_normalized" in reference_adata.layers:
                ref_mat = reference_adata.layers["X_normalized"]
            else:
                ref_mat = reference_adata.X

            # Convert sparse to dense if needed
            if hasattr(ref_mat, "toarray"):
                ref_mat = ref_mat.toarray()

            # Ensure log-normalization for reference
            if "log1p" not in reference_adata.uns:
                if context:
                    await context.warning(
                        "Reference data may not be log-normalized. Applying log1p..."
                    )
                ref_mat = np.log1p(ref_mat)

            # Transpose for SingleR (genes x cells)
            ref_mat = ref_mat.T
            ref_features = [str(x) for x in reference_adata.var_names]

            # Check gene overlap
            common_genes = set(test_features) & set(ref_features)
            if context:
                await context.info(
                    f"Gene overlap: {len(common_genes)}/{len(test_features)} test genes match reference"
                )

            if len(common_genes) < min(50, len(test_features) * 0.1):
                # Too few common genes
                if context:
                    await context.error(
                        f"Insufficient gene overlap: Only {len(common_genes)} common genes. "
                        f"Test has {len(test_features)}, Reference has {len(ref_features)}"
                    )
                    # Show sample of gene names for debugging
                    test_sample = list(test_features)[:5]
                    ref_sample = list(ref_features)[:5]
                    await context.info(f"Test gene sample: {test_sample}")
                    await context.info(f"Reference gene sample: {ref_sample}")
                raise ValueError(
                    f"Insufficient gene overlap for SingleR: only {len(common_genes)} common genes"
                )

            # Get labels from reference - check various common column names
            # cell_type_key is now required (no default value)
            cell_type_key = params.cell_type_key

            if cell_type_key not in reference_adata.obs.columns:
                # Improved error message showing available columns
                available_cols = list(reference_adata.obs.columns)

                # Categorize columns by type for better guidance
                categorical_cols = [
                    col
                    for col in available_cols
                    if reference_adata.obs[col].dtype.name in ["object", "category"]
                ]

                raise ValueError(
                    f"Cell type column '{cell_type_key}' not found in reference data.\n\n"
                    f"Available categorical columns (likely cell types):\n  {', '.join(categorical_cols[:15])}\n"
                    f"{f'  ... and {len(categorical_cols)-15} more' if len(categorical_cols) > 15 else ''}\n\n"
                    f"All columns ({len(available_cols)} total):\n  {', '.join(available_cols[:20])}\n"
                    f"{f'  ... and {len(available_cols)-20} more' if len(available_cols) > 20 else ''}\n\n"
                    f"Please specify the correct column using: cell_type_key='your_column_name'"
                )

            ref_labels = list(reference_adata.obs[cell_type_key])
            if context:
                await context.info(
                    f"Using '{cell_type_key}' column for reference labels"
                )

            # For SingleR, pass the actual expression matrix directly (not SCE)
            # This has been shown to work better in testing
            ref_data = ref_mat
            ref_features_to_use = (
                ref_features  # Keep reference features for gene matching
            )

        elif celldex:
            # Use default reference
            if context:
                await context.info("Using default BlueprintEncode reference")
            ref = celldex.fetch_reference(
                "blueprint_encode", "2024-02-26", realize_assays=True
            )
            ref_labels = ref.get_column_data().column("label.main")
            ref_data = ref
        else:
            raise ValueError(
                "No reference data available. Please either:\n"
                "1. Provide reference_data_id\n"
                "2. Provide singler_reference name\n"
                "3. Install celldex for pre-built references"
            )

        # Run SingleR annotation
        if context:
            await context.info("Running SingleR classification...")
            await context.info(f"Test features: {len(test_features)} genes")
            await context.info(f"Reference labels: {len(set(ref_labels))} unique types")

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
            if "ref_features_to_use" in locals() and ref_features_to_use is not None:
                annotate_kwargs["ref_features"] = ref_features_to_use
                if context:
                    await context.info(
                        f"Using reference features: {len(ref_features_to_use)} genes"
                    )

            results = singler.annotate_single(**annotate_kwargs)
            best_labels = results.column("best")
            scores = results.column("scores")

            # Try to get delta scores for confidence (higher delta = higher confidence)
            try:
                delta_scores = results.column("delta")
                if delta_scores and context:
                    low_delta = sum(1 for d in delta_scores if d and d < 0.05)
                    if low_delta > len(delta_scores) * 0.3:
                        await context.warning(
                            f"⚠️ {low_delta}/{len(delta_scores)} cells have low confidence scores (delta < 0.05)"
                        )
            except:
                delta_scores = None

        # Process results
        cell_types = list(best_labels)
        unique_types = list(set(cell_types))
        counts = pd.Series(cell_types).value_counts().to_dict()

        # Calculate confidence scores - prefer delta scores if available
        confidence_scores = {}

        # First try to use delta scores (more meaningful confidence measure)
        if "delta_scores" in locals() and delta_scores is not None:
            try:
                for cell_type in unique_types:
                    type_indices = [
                        i for i, ct in enumerate(cell_types) if ct == cell_type
                    ]
                    if type_indices:
                        type_deltas = [
                            delta_scores[i]
                            for i in type_indices
                            if i < len(delta_scores)
                        ]
                        if type_deltas:
                            # Higher delta = higher confidence
                            avg_delta = np.mean(
                                [d for d in type_deltas if d is not None]
                            )
                            confidence_scores[cell_type] = round(float(avg_delta), 3)
                if context and confidence_scores:
                    await context.info(
                        f"Using delta scores for confidence (avg: {np.mean(list(confidence_scores.values())):.3f})"
                    )
            except:
                pass  # Fall back to regular scores

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
        adata.obs[output_key] = adata.obs[output_key].astype("category")

        # Only add confidence column if we have real confidence values
        if confidence_scores:
            # Use 0.0 for cells without confidence (more honest than arbitrary 0.5)
            confidence_array = [confidence_scores.get(ct, 0.0) for ct in cell_types]
            adata.obs[confidence_key] = confidence_array

        if context:
            await context.info("✅ SingleR annotation completed!")
            await context.info(f"   Found {len(unique_types)} cell types")
            top_types = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:3]
            await context.info(
                f"   Top types: {', '.join([f'{t}({c})' for t, c in top_types])}"
            )
            await context.info(
                f"Cell type annotation complete. Cell types stored in '{output_key}' column"
            )
            await context.info(
                f"Use visualize_data tool with feature='{output_key}' to visualize results"
            )

        return unique_types, counts, confidence_scores, None

    except Exception as e:
        await _handle_annotation_error(e, "singler", context)


async def _annotate_with_tangram(
    adata,
    params: AnnotationParameters,
    data_store: Dict[str, Any],
    context: Optional[Context] = None,
):
    """Annotate cell types using Tangram method"""
    try:
        if context:
            await context.info("Using Tangram method for annotation")

        # Validate dependencies with comprehensive error reporting
        tg = _validate_tangram(context)

        # Check if reference data is provided
        if params.reference_data_id is None:
            raise ValueError("Reference data ID is required for Tangram method")

        if params.reference_data_id not in data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")

        # Get reference single-cell data
        adata_sc_original = data_store[params.reference_data_id]["adata"]
        adata_sp = adata  # Spatial data (will be modified in-place for results)

        # ===== Handle Duplicate Gene Names (CRITICAL FIX) =====
        if not adata_sc_original.var_names.is_unique:
            n_duplicates = (
                len(adata_sc_original.var_names)
                - len(set(adata_sc_original.var_names))
            )
            if context:
                await context.warning(
                    f"⚠️  Found {n_duplicates} duplicate gene names in reference data, fixing..."
                )
            adata_sc_original.var_names_make_unique()

        if not adata_sp.var_names.is_unique:
            n_duplicates = len(adata_sp.var_names) - len(set(adata_sp.var_names))
            if context:
                await context.warning(
                    f"⚠️  Found {n_duplicates} duplicate gene names in spatial data, fixing..."
                )
            adata_sp.var_names_make_unique()

        if context:
            await context.info(
                f"Using reference dataset {params.reference_data_id} with {adata_sc_original.n_obs} cells"
            )

        # Determine training genes
        training_genes = params.training_genes

        if training_genes is None:
            # Use marker genes if available
            if params.marker_genes:
                if context:
                    await context.info(
                        "Using provided marker genes for Tangram mapping"
                    )
                # Flatten marker genes dictionary
                training_genes = []
                for genes in params.marker_genes.values():
                    training_genes.extend(genes)
                training_genes = list(set(training_genes))  # Remove duplicates
            else:
                # Use highly variable genes
                if "highly_variable" not in adata_sc_original.var:
                    raise ValueError(
                        "Highly variable genes not found in reference data. "
                        "Tangram mapping requires HVG selection. "
                        "Please run HVG selection in preprocessing.py or use: sc.pp.highly_variable_genes(adata)"
                    )
                training_genes = list(
                    adata_sc_original.var_names[adata_sc_original.var.highly_variable]
                )

        if context:
            await context.info(f"Using {len(training_genes)} genes for Tangram mapping")

        # ✅ COW FIX: Create copy of reference data to avoid modifying original
        # Tangram's pp_adatas adds metadata (uns, obs) but doesn't subset genes
        adata_sc = adata_sc_original.copy()
        if context:
            await context.info(
                "Created working copy of reference data to preserve original"
            )

        # Preprocess data for Tangram
        tg.pp_adatas(adata_sc, adata_sp, genes=training_genes)

        if context:
            await context.info(
                f"Preprocessed data for Tangram mapping. {len(adata_sc.uns['training_genes'])} training genes selected."
            )

        # Set mapping mode
        mode = params.mode
        cluster_label = params.cluster_label

        if mode == "clusters" and cluster_label is None:
            if context:
                await context.warning(
                    "Cluster label not provided for 'clusters' mode. Using default cell type annotation if available."
                )
            # Try to find a cell type annotation in the reference data
            for col in ["cell_type", "celltype", "cell_types", "leiden", "louvain"]:
                if col in adata_sc.obs:
                    cluster_label = col
                    break

            if cluster_label is None:
                raise ValueError(
                    "No cluster label found in reference data. Please provide a cluster_label parameter."
                )

            if context:
                await context.info(
                    f"Using '{cluster_label}' as cluster label for Tangram mapping"
                )

        # Check GPU availability for device selection
        import torch

        device = params.tangram_device
        if device != "cpu" and torch.cuda.is_available():
            if context:
                await context.info(f"GPU detected - using device: {device}")
        elif device != "cpu":
            if context:
                await context.warning(
                    "GPU requested but not available - falling back to CPU"
                )
            device = "cpu"

        # Run Tangram mapping with enhanced parameters
        if context:
            await context.info(
                f"Running Tangram mapping in '{mode}' mode for {params.num_epochs} epochs on {device}"
            )

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
                            except:
                                tangram_mapping_score = 0.0
                    else:
                        tangram_mapping_score = float(last_value)

                    if context:
                        await context.info(
                            f"Extracted Tangram mapping score: {tangram_mapping_score:.4f}"
                        )

                else:
                    # NO FALLBACK: Require modern Tangram format for scientific integrity
                    error_msg = (
                        "❌ Tangram training history format not recognized.\n\n"
                        "Expected dictionary with 'main_loss' key, but got: "
                        f"{type(history).__name__ if history else 'None'}\n\n"
                        "🔧 SOLUTIONS:\n"
                        "1. Ensure using modern Tangram-sc version:\n"
                        "   pip install --upgrade tangram-sc\n\n"
                        "2. Verify Tangram completed training successfully\n\n"
                        "3. Check that mapping_result has valid training_history\n\n"
                        "📋 SCIENTIFIC INTEGRITY: We require consistent Tangram output format "
                        "to ensure reproducible cell type mapping scores."
                    )
                    if context:
                        await context.error(error_msg)
                    raise ValueError(error_msg)
        except Exception as score_error:
            if context:
                await context.error(
                    f"Failed to extract Tangram mapping score: {score_error}"
                )
            # Don't return a misleading score - fail gracefully
            raise RuntimeError(
                f"Tangram mapping completed but score extraction failed: {score_error}"
            )

        # Compute validation metrics if requested
        if params.tangram_compute_validation:
            try:
                if context:
                    await context.info(
                        "Computing validation metrics for spatial gene expression"
                    )
                # Compare predicted vs actual spatial gene expression
                scores = tg.compare_spatial_geneexp(ad_map, adata_sp)
                # Store validation scores in adata
                adata_sp.uns["tangram_validation_scores"] = scores
                if context:
                    await context.info(
                        "Validation completed - scores stored in adata.uns['tangram_validation_scores']"
                    )
            except Exception as val_error:
                if context:
                    await context.warning(
                        f"Could not compute validation metrics: {val_error}"
                    )

        # Project genes if requested
        if params.tangram_project_genes:
            try:
                if context:
                    await context.info("Projecting gene expression to spatial data")
                ad_ge = tg.project_genes(ad_map, adata_sc)
                adata_sp.obsm["tangram_gene_predictions"] = ad_ge.X
                if context:
                    await context.info(
                        "Gene expression projections stored in adata.obsm['tangram_gene_predictions']"
                    )
            except Exception as gene_error:
                if context:
                    await context.warning(f"Could not project genes: {gene_error}")

        if context:
            await context.info(
                f"Tangram mapping completed with score: {tangram_mapping_score}"
            )

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

                    raise ValueError(
                        f"Cell type column '{params.cell_type_key}' not found in reference data.\n\n"
                        f"Available categorical columns:\n  {', '.join(categorical_cols[:15])}\n"
                        f"{f'  ... and {len(categorical_cols)-15} more' if len(categorical_cols) > 15 else ''}\n\n"
                        f"All columns ({len(available_cols)} total):\n  {', '.join(available_cols[:20])}\n"
                        f"{f'  ... and {len(available_cols)-20} more' if len(available_cols) > 20 else ''}"
                    )

                annotation_col = params.cell_type_key
                if context:
                    await context.info(
                        f"Using cell type column: '{params.cell_type_key}'"
                    )

            if annotation_col:
                if context:
                    await context.info(
                        f"Projecting cell annotations using '{annotation_col}' column"
                    )
                # Use project_cell_annotations instead of plot_cell_annotation
                tg.project_cell_annotations(ad_map, adata_sp, annotation=annotation_col)
            else:
                if context:
                    await context.warning(
                        "No suitable annotation column found for projection"
                    )
        except Exception as proj_error:
            if context:
                await context.warning(
                    f"Could not project cell annotations: {proj_error}"
                )
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

            # Assign cell type based on highest probability
            adata_sp.obs[output_key] = cell_type_df.idxmax(axis=1)
            adata_sp.obs[output_key] = adata_sp.obs[output_key].astype("category")

            # Get counts
            counts = adata_sp.obs[output_key].value_counts().to_dict()

            # Calculate confidence scores (use max probability as confidence)
            confidence_scores = {}
            for cell_type in cell_types:
                cells_of_type = adata_sp.obs[output_key] == cell_type
                if np.sum(cells_of_type) > 0:
                    # Use mean probability as confidence
                    mean_prob = cell_type_df.loc[cells_of_type, cell_type].mean()
                    confidence_scores[cell_type] = round(float(mean_prob), 2)
                # Note: If no cells assigned to this type, we don't report a confidence
                # This is more scientifically honest than assigning an arbitrary value

            # Note: Visualizations should be created using the separate visualize_data tool
            # This maintains clean separation between analysis and visualization
            if context:
                await context.info(
                    f"Cell type mapping complete. Cell types stored in '{output_key}' column"
                )
                await context.info(
                    f"Use visualize_data tool with feature='{output_key}' to visualize results"
                )

        else:
            if context:
                await context.warning(
                    "No cell type predictions found in Tangram results"
                )

        # Validate results before returning
        if not cell_types:
            raise RuntimeError(
                "Tangram mapping failed - no cell type predictions generated"
            )

        if tangram_mapping_score <= 0:
            if context:
                await context.warning(
                    f"Tangram mapping score is suspiciously low: {tangram_mapping_score}"
                )

        if context:
            await context.info(
                f"Cell type annotation complete. Cell types stored in '{output_key}' column"
            )
            await context.info(
                f"Use visualize_data tool with feature='{output_key}' to visualize results"
            )

        return cell_types, counts, confidence_scores, tangram_mapping_score

    except Exception as e:
        await _handle_annotation_error(e, "tangram", context)


async def _annotate_with_scanvi(
    adata,
    params: AnnotationParameters,
    data_store: Dict[str, Any],
    context: Optional[Context] = None,
):
    """Annotate cell types using scANVI method with official best practices"""
    try:
        if context:
            await context.info("Using scANVI method for annotation")

        # Validate dependencies with comprehensive error reporting
        scvi, CellAssign = _validate_scvi_tools(context)

        # Check if reference data is provided
        if params.reference_data_id is None:
            raise ValueError("Reference data ID is required for scANVI method")

        if params.reference_data_id not in data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")

        # Get reference single-cell data
        adata_ref_original = data_store[params.reference_data_id]["adata"]

        # ===== Handle Duplicate Gene Names (CRITICAL FIX) =====
        if context:
            await context.info("Checking for duplicate gene names...")

        # Fix duplicate gene names in reference data
        if not adata_ref_original.var_names.is_unique:
            n_duplicates_ref = (
                len(adata_ref_original.var_names)
                - len(set(adata_ref_original.var_names))
            )
            if context:
                await context.warning(
                    f"⚠️  Found {n_duplicates_ref} duplicate gene names in reference data"
                )
                await context.info("Fixing duplicate gene names with unique suffixes...")
            adata_ref_original.var_names_make_unique()
            if context:
                await context.info("✓ Reference gene names fixed")

        # Fix duplicate gene names in query data
        if not adata.var_names.is_unique:
            n_duplicates_query = len(adata.var_names) - len(set(adata.var_names))
            if context:
                await context.warning(
                    f"⚠️  Found {n_duplicates_query} duplicate gene names in query data"
                )
                await context.info("Fixing duplicate gene names with unique suffixes...")
            adata.var_names_make_unique()
            if context:
                await context.info("✓ Query gene names fixed")

        # ===== Gene Alignment (NEW) =====
        if context:
            await context.info("Aligning genes between reference and query data...")

        common_genes = adata_ref_original.var_names.intersection(adata.var_names)

        if len(common_genes) < min(100, adata_ref_original.n_vars * 0.5):
            raise ValueError(
                f"Insufficient gene overlap: Only {len(common_genes)} common genes found. "
                f"Reference has {adata_ref_original.n_vars}, query has {adata.n_vars} genes."
            )

        # ✅ COW FIX: Operate on temporary copies for gene subsetting
        # This prevents loss of HVG information in the original adata
        if len(common_genes) < adata_ref_original.n_vars:
            if context:
                await context.warning(
                    f"Subsetting to {len(common_genes)} common genes for ScanVI training "
                    f"(reference: {adata_ref_original.n_vars}, query: {adata.n_vars})"
                )
                await context.info(
                    "Note: Original gene set and HVG information will be preserved"
                )
            # Create subsets for ScanVI (not modifying originals)
            adata_ref = adata_ref_original[:, common_genes].copy()
            adata_subset = adata[:, common_genes].copy()
        else:
            # No subsetting needed
            adata_ref = adata_ref_original.copy()
            adata_subset = adata.copy()

        # ===== Data Validation (NEW) =====
        if context:
            await context.info("Validating data preprocessing...")

        # ✅ DEBUG: Check for counts layer existence
        if context:
            await context.info(f"Reference data has 'counts' layer: {'counts' in adata_ref.layers}")
            await context.info(f"Query data has 'counts' layer: {'counts' in adata_subset.layers}")
            if 'counts' in adata_ref.layers:
                await context.info(f"  Reference counts max: {adata_ref.layers['counts'].max():.2f}")
            if 'counts' in adata_subset.layers:
                await context.info(f"  Query counts max: {adata_subset.layers['counts'].max():.2f}")

        # Check if reference data is normalized
        if "log1p" not in adata_ref.uns:
            if context:
                await context.warning("Reference data may not be log-normalized")

        # Check for HVGs
        if "highly_variable" not in adata_ref.var:
            if context:
                await context.warning("No highly variable genes detected in reference")

        # Get parameters
        cell_type_key = getattr(params, "cell_type_key", "cell_type")
        batch_key = getattr(params, "batch_key", None)

        # ===== Optional SCVI Pretraining (NEW) =====
        if params.scanvi_use_scvi_pretrain:
            if context:
                await context.info(
                    f"Step 1/3: Pretraining SCVI model for {params.scanvi_scvi_epochs} epochs..."
                )

            # Setup for SCVI with labels (required for SCANVI conversion)
            # First ensure the reference has the cell type labels
            if cell_type_key not in adata_ref.obs.columns:
                raise ValueError(f"Reference data missing '{cell_type_key}' column")

            # SCVI needs to know about labels for later SCANVI conversion
            scvi.model.SCVI.setup_anndata(
                adata_ref,
                labels_key=cell_type_key,  # Important: include labels_key
                batch_key=batch_key,
                layer=params.layer if hasattr(params, "layer") else None,
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

            if context:
                await context.info("Step 2/3: Converting SCVI to SCANVI model...")

            # Convert to SCANVI (no need for setup_anndata, it uses SCVI's setup)
            model = scvi.model.SCANVI.from_scvi_model(
                scvi_model, params.scanvi_unlabeled_category
            )

            # Train SCANVI (fewer epochs needed after pretraining)
            model.train(
                max_epochs=20,
                n_samples_per_label=params.scanvi_n_samples_per_label,
                early_stopping=True,
            )

        else:
            # Direct SCANVI training (existing approach)
            if context:
                await context.info("Training SCANVI model directly...")

            # Setup AnnData for scANVI
            # ✅ FIX: Use raw counts from layers['counts'] instead of normalized adata.X
            scvi.model.SCANVI.setup_anndata(
                adata_ref,
                labels_key=cell_type_key,
                unlabeled_category=params.scanvi_unlabeled_category,
                batch_key=batch_key,
                layer="counts",  # CRITICAL: Use raw counts, not normalized data
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

        # ===== Query Data Preparation (IMPROVED) =====
        if context:
            await context.info("Step 3/3: Preparing and training on query data...")

        # ✅ COW FIX: Work on adata_subset for query data preparation
        # Add unlabeled category for all cells (on subset, not original)
        adata_subset.obs[cell_type_key] = params.scanvi_unlabeled_category

        # Setup query data (batch handling) - TRANSPARENT TEMPORARY METADATA
        if batch_key and batch_key not in adata_subset.obs:
            # ScANVI requires batch information for technical reasons
            adata_subset.obs[batch_key] = "query_batch"
            if context:
                await context.info(
                    f"🔧 Added temporary batch label '{batch_key}' = 'query_batch' for ScANVI compatibility.\n"
                    f"   This is TECHNICAL METADATA, not real batch information from your experiment.\n"
                    f"   ScANVI algorithm requires batch labels to function properly."
                )

        # ✅ FIX: Use raw counts for query data as well
        scvi.model.SCANVI.setup_anndata(
            adata_subset,
            labels_key=cell_type_key,
            unlabeled_category=params.scanvi_unlabeled_category,
            batch_key=batch_key,
            layer="counts",  # CRITICAL: Use raw counts, not normalized data
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

        # ✅ COW FIX: Get predictions from adata_subset, then add to original adata
        predictions = spatial_model.predict()
        adata_subset.obs[cell_type_key] = predictions
        adata_subset.obs[cell_type_key] = adata_subset.obs[cell_type_key].astype(
            "category"
        )

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
            if context:
                await context.warning(f"Could not get confidence scores: {e}")
            # Could not extract probabilities - return empty confidence dict
            confidence_scores = (
                {}
            )  # Empty dict clearly indicates no confidence data available

        # ✅ COW FIX: Add prediction results to original adata.obs
        # This will be picked up by the main function to store in adata.obs[output_key]
        adata.obs[cell_type_key] = adata_subset.obs[cell_type_key].values

        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization

        if context:
            await context.info(
                "Cell type annotation complete. Cell types stored in 'cell_type' column"
            )
            await context.info(
                "Use visualize_data tool with feature='cell_type' to visualize results"
            )

        return cell_types, counts, confidence_scores, None

    except Exception as e:
        await _handle_annotation_error(e, "scanvi", context)


async def _annotate_with_mllmcelltype(
    adata, params: AnnotationParameters, context: Optional[Context] = None
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
    try:
        if context:
            await context.info("Using mLLMCellType (LLM-based) method for annotation")

        # Validate dependencies with comprehensive error reporting
        mllmcelltype = _validate_mllmcelltype(context)

        # Validate clustering has been performed
        # cluster_label is now required for mLLMCellType (no default value)
        if not params.cluster_label:
            available_cols = list(adata.obs.columns)
            categorical_cols = [
                col
                for col in available_cols
                if adata.obs[col].dtype.name in ["object", "category"]
            ]

            raise ValueError(
                f"cluster_label parameter is required for mLLMCellType method.\n\n"
                f"Available categorical columns (likely clusters):\n  {', '.join(categorical_cols[:15])}\n"
                f"{f'  ... and {len(categorical_cols)-15} more' if len(categorical_cols) > 15 else ''}\n\n"
                f"Common cluster column names: leiden, louvain, seurat_clusters, phenograph\n\n"
                f"Example: params = {{'cluster_label': 'leiden', ...}}"
            )

        cluster_key = params.cluster_label

        if cluster_key not in adata.obs:
            available_cols = list(adata.obs.columns)
            categorical_cols = [
                col
                for col in available_cols
                if adata.obs[col].dtype.name in ["object", "category"]
            ]

            raise ValueError(
                f"Clustering key '{cluster_key}' not found in adata.obs.\n\n"
                f"Available categorical columns:\n  {', '.join(categorical_cols[:15])}\n"
                f"{f'  ... and {len(categorical_cols)-15} more' if len(categorical_cols) > 15 else ''}\n\n"
                f"mLLMCellType annotation requires clustering information.\n"
                f"Please run clustering in preprocessing.py first."
            )

        # Find differentially expressed genes for each cluster
        if context:
            await context.info("Finding marker genes for each cluster")

        sc.tl.rank_genes_groups(adata, cluster_key, method="wilcoxon")

        # Extract top marker genes for each cluster
        marker_genes_dict = {}
        n_genes = (
            params.mllm_n_marker_genes if hasattr(params, "mllm_n_marker_genes") else 20
        )

        for cluster in adata.obs[cluster_key].unique():
            # Get top genes for this cluster
            gene_names = adata.uns["rank_genes_groups"]["names"][str(cluster)][:n_genes]
            marker_genes_dict[f"Cluster_{cluster}"] = list(gene_names)

        if context:
            await context.info(
                f"Found marker genes for {len(marker_genes_dict)} clusters"
            )

        # Prepare parameters for mllmcelltype
        species = params.mllm_species if hasattr(params, "mllm_species") else "human"
        tissue = params.mllm_tissue if hasattr(params, "mllm_tissue") else None
        additional_context = params.mllm_additional_context if hasattr(params, "mllm_additional_context") else None
        use_cache = params.mllm_use_cache if hasattr(params, "mllm_use_cache") else True
        base_urls = params.mllm_base_urls if hasattr(params, "mllm_base_urls") else None
        verbose = params.mllm_verbose if hasattr(params, "mllm_verbose") else False
        force_rerun = params.mllm_force_rerun if hasattr(params, "mllm_force_rerun") else False
        clusters_to_analyze = params.mllm_clusters_to_analyze if hasattr(params, "mllm_clusters_to_analyze") else None

        # Check if using multi-model consensus or single model
        use_consensus = params.mllm_use_consensus if hasattr(params, "mllm_use_consensus") else False

        try:
            if use_consensus:
                # Use interactive_consensus_annotation with multiple models
                models = params.mllm_models if hasattr(params, "mllm_models") else None
                if not models:
                    raise ValueError(
                        "mllm_models parameter is required when mllm_use_consensus=True. "
                        "Provide a list of model names, e.g., ['gpt-5', 'claude-sonnet-4-5-20250929', 'gemini-2.5-pro']"
                    )

                api_keys = params.mllm_api_keys if hasattr(params, "mllm_api_keys") else None
                consensus_threshold = params.mllm_consensus_threshold if hasattr(params, "mllm_consensus_threshold") else 0.7
                entropy_threshold = params.mllm_entropy_threshold if hasattr(params, "mllm_entropy_threshold") else 1.0
                max_discussion_rounds = params.mllm_max_discussion_rounds if hasattr(params, "mllm_max_discussion_rounds") else 3
                consensus_model = params.mllm_consensus_model if hasattr(params, "mllm_consensus_model") else None

                if context:
                    await context.info(
                        f"Calling interactive consensus annotation with {len(models)} models"
                    )

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

                # Store additional metadata for reporting
                consensus_proportions = consensus_results.get("consensus_proportions", {})
                entropy = consensus_results.get("entropy", {})

                if context:
                    await context.info(
                        f"Consensus reached for {len(annotations)} clusters. "
                        f"Mean consensus proportion: {sum(consensus_proportions.values()) / len(consensus_proportions):.2f}"
                    )
            else:
                # Use single model annotation
                provider = params.mllm_provider if hasattr(params, "mllm_provider") else "openai"
                model = params.mllm_model if hasattr(params, "mllm_model") else None
                api_key = params.mllm_api_key if hasattr(params, "mllm_api_key") else None

                if context:
                    await context.info(
                        f"Calling LLM ({provider}/{model or 'default'}) for cell type annotation"
                    )

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
            if context:
                await context.error(f"mLLMCellType annotation failed: {str(e)}")
            raise

        if context:
            await context.info(f"Received annotations for {len(annotations)} clusters")

        # Map cluster annotations back to cells
        cluster_to_celltype = {}
        for cluster_name, cell_type in annotations.items():
            # Extract cluster number from "Cluster_X" format
            cluster_id = cluster_name.replace("Cluster_", "")
            cluster_to_celltype[cluster_id] = cell_type

        # Use method-prefixed output keys to avoid overwriting
        output_key = f"cell_type_{params.method}"

        # Apply cell type annotations to cells
        adata.obs[output_key] = (
            adata.obs[cluster_key].astype(str).map(cluster_to_celltype)
        )

        # Handle any unmapped clusters
        unmapped = adata.obs[output_key].isna()
        if unmapped.any():
            if context:
                await context.warning(
                    f"Found {unmapped.sum()} cells in unmapped clusters"
                )
            adata.obs.loc[unmapped, output_key] = "Unknown"

        adata.obs[output_key] = adata.obs[output_key].astype("category")

        # Get cell types and counts
        cell_types = list(adata.obs[output_key].unique())
        counts = adata.obs[output_key].value_counts().to_dict()

        # Calculate confidence scores based on cluster homogeneity
        confidence_scores = {}
        for cell_type in cell_types:
            if cell_type != "Unknown":
                # LLM-based annotations don't provide numeric confidence scores
                # Don't assign arbitrary values that could be misleading
                pass  # No numeric confidence available
            else:
                # Unknown cells have no confidence
                pass

        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization
        if context:
            await context.info(
                f"Cell type annotation complete. Cell types stored in '{output_key}' column"
            )
            await context.info(
                f"Use visualize_data tool with feature='{output_key}' to visualize results"
            )

        return cell_types, counts, confidence_scores, None, None

    except Exception as e:
        await _handle_annotation_error(e, "mllmcelltype", context)


async def _annotate_with_cellassign(
    adata, params: AnnotationParameters, context: Optional[Context] = None
):
    """Annotate cell types using CellAssign method"""
    try:
        if context:
            await context.info("Using CellAssign method for annotation")

        # Validate dependencies with comprehensive error reporting
        scvi, CellAssign = _validate_scvi_tools(context)

        # Check if marker genes are provided
        if params.marker_genes is None:
            raise ValueError(
                "CellAssign requires marker genes to be provided. "
                "Please specify marker_genes parameter with a dictionary of cell types and their marker genes."
            )

        marker_genes = params.marker_genes

        # Validate marker genes exist in dataset
        all_genes = set(adata.var_names)
        valid_marker_genes = {}

        for cell_type, genes in marker_genes.items():
            existing_genes = [gene for gene in genes if gene in all_genes]
            if existing_genes:
                valid_marker_genes[cell_type] = existing_genes
                if context:
                    await context.info(
                        f"Found {len(existing_genes)} marker genes for {cell_type}"
                    )
            elif context:
                await context.warning(f"No marker genes found for {cell_type}")

        if not valid_marker_genes:
            raise ValueError("No valid marker genes found for any cell type")
        valid_cell_types = list(valid_marker_genes.keys())

        # Create marker gene matrix as DataFrame (required by CellAssign API)
        all_marker_genes = []
        for genes in valid_marker_genes.values():
            all_marker_genes.extend(genes)
        available_marker_genes = list(set(all_marker_genes))  # Remove duplicates

        if not available_marker_genes:
            raise ValueError("No marker genes found in the dataset")

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

        # Subset data to only marker genes first
        adata_subset = adata[:, available_marker_genes].copy()

        # Check for invalid values in the data
        if hasattr(adata_subset.X, "toarray"):
            X_array = adata_subset.X.toarray()
        else:
            X_array = adata_subset.X

        # Replace any NaN or Inf values with zeros
        if np.any(np.isnan(X_array)) or np.any(np.isinf(X_array)):
            if context:
                await context.warning(
                    "Found NaN or Inf values in data, replacing with zeros"
                )
            X_array = np.nan_to_num(X_array, nan=0.0, posinf=0.0, neginf=0.0)
            adata_subset.X = X_array

        # Additional data cleaning for CellAssign compatibility
        # Check for genes with zero variance (which cause numerical issues in CellAssign)
        gene_vars = np.var(X_array, axis=0)
        zero_var_genes = gene_vars == 0
        if np.any(zero_var_genes):
            zero_var_gene_names = adata_subset.var_names[zero_var_genes].tolist()
            raise ValueError(
                f"Found {np.sum(zero_var_genes)} genes with zero variance: {zero_var_gene_names[:10]}... "
                f"CellAssign requires all genes to have non-zero variance. "
                f"Please filter these genes in preprocessing: adata = adata[:, adata.var['std'] > 0] "
                f"or use sc.pp.highly_variable_genes() to select informative genes."
            )

        # Ensure data is non-negative (CellAssign expects count-like data)
        if np.any(X_array < 0):
            if context:
                await context.warning("Found negative values in data, clipping to zero")
            X_array = np.maximum(X_array, 0)
            adata_subset.X = X_array

        # Add size factors if not present
        if "size_factors" not in adata_subset.obs:
            # Calculate size factors from subset data
            if hasattr(adata_subset.X, "sum"):
                size_factors = adata_subset.X.sum(axis=1)
                if hasattr(size_factors, "A1"):  # sparse matrix
                    size_factors = size_factors.A1
            else:
                size_factors = np.sum(adata_subset.X, axis=1)

            # Ensure size factors are positive (avoid division by zero)
            size_factors = np.maximum(size_factors, 1e-6)
            adata_subset.obs["size_factors"] = pd.Series(
                size_factors, index=adata_subset.obs.index
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
        confidence_key = f"confidence_{params.method}"

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

        adata.obs[output_key] = adata.obs[output_key].astype("category")

        # Get cell types and counts
        cell_types = valid_cell_types
        counts = adata.obs[output_key].value_counts().to_dict()

        # Note: Visualizations should be created using the separate visualize_data tool
        # This maintains clean separation between analysis and visualization

        if context:
            await context.info(
                f"Cell type annotation complete. Cell types stored in '{output_key}' column"
            )
            await context.info(
                f"Use visualize_data tool with feature='{output_key}' to visualize results"
            )

        return cell_types, counts, confidence_scores, None

    except Exception as e:
        await _handle_annotation_error(e, "cellassign", context)


async def annotate_cell_types(
    data_id: str,
    data_store: Dict[str, Any],
    params: AnnotationParameters,  # No default - must be provided by caller (LLM)
    context: Optional[Context] = None,
) -> AnnotationResult:
    """Annotate cell types in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Annotation parameters
        context: MCP context

    Returns:
        Annotation result
    """
    if context:
        await context.info(f"Annotating cell types using {params.method} method")

    # Retrieve the AnnData object from data store
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    adata = data_store[data_id]["adata"]

    # Validate method first - clean and simple
    if params.method not in SUPPORTED_METHODS:
        raise ValueError(
            f"Unsupported method: {params.method}. Supported: {sorted(SUPPORTED_METHODS)}"
        )

    # Route to appropriate annotation method
    try:
        if params.method == "tangram":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_tangram(adata, params, data_store, context)
            )
        elif params.method == "scanvi":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_scanvi(adata, params, data_store, context)
            )
        elif params.method == "cellassign":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_cellassign(adata, params, context)
            )
        elif params.method == "mllmcelltype":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_mllmcelltype(adata, params, context)
            )
        elif params.method == "singler":
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_singler(adata, params, data_store, context)
            )
        else:  # sctype
            cell_types, counts, confidence_scores, tangram_mapping_score = (
                await _annotate_with_sctype(adata, params, context)
            )

    except Exception as e:
        if context:
            await context.error(f"Annotation failed: {str(e)}")
        raise

    # ✅ COW FIX: Each annotation method handles adata.obs assignment internally
    # This maintains consistency and prevents duplicate assignments
    # Construct output keys for result reporting
    output_key = f"cell_type_{params.method}"
    confidence_key = f"confidence_{params.method}" if confidence_scores else None

    # ❌ REMOVED: data_store[data_id]["adata"] = adata
    # No longer overwrite adata to preserve HVG information and original gene set

    # Inform user about visualization options
    if context:
        await context.info(
            f"Cell type annotation complete. Cell types stored in '{output_key}' column"
        )
        await context.info(
            f"Use visualize_data tool with feature='{output_key}' to visualize results"
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

# Cache for sc-type results to avoid repeated R calls
_SCTYPE_CACHE = {}
_SCTYPE_CACHE_DIR = Path.home() / ".chatspatial" / "sctype_cache"


def _get_sctype_cache_key(adata, params: AnnotationParameters) -> str:
    """Generate cache key for sc-type results"""
    # Create a hash based on data and parameters
    data_hash = hashlib.md5()

    # Hash expression data (sample first 1000 cells and 500 genes for efficiency)
    if hasattr(adata.X, "toarray"):
        sample_data = adata.X[
            : min(1000, adata.n_obs), : min(500, adata.n_vars)
        ].toarray()
    else:
        sample_data = adata.X[: min(1000, adata.n_obs), : min(500, adata.n_vars)]
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


async def _load_sctype_functions(context: Optional[Context] = None) -> None:
    """Load sc-type R functions and auto-install R packages if needed"""
    if context:
        await context.info("📚 Loading sc-type R functions...")

    try:
        # Get robjects from validation
        robjects, _, _, _, _, default_converter, openrlib = _validate_rpy2_and_r(
            context
        )

        # Import conversion module
        from rpy2.robjects import conversion

        # Wrap R calls in conversion context (FIX for contextvars issue)
        with openrlib.rlock:
            with conversion.localconverter(default_converter):
                # Auto-install required R packages
                robjects.r(
                    """
                    # Install packages automatically if not present
                    required_packages <- c("dplyr", "openxlsx", "HGNChelper")
                    for (pkg in required_packages) {
                        if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
                            cat("Installing R package:", pkg, "\\n")
                            install.packages(pkg, repos = "https://cran.r-project.org/", quiet = TRUE)
                            if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
                                stop(paste("Failed to install required R package:", pkg))
                            }
                        }
                    }
                """
                )

                if context:
                    await context.info("✅ R packages loaded/installed successfully")

                # Load sc-type functions from GitHub
                robjects.r(
                    """
                    # Load gene sets preparation function
                    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")

                    # Load scoring function
                    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
                """
                )

        if context:
            await context.info("✅ sc-type R functions loaded successfully")

    except Exception as e:
        error_msg = f"Failed to load sc-type R functions: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


async def _prepare_sctype_genesets(
    params: AnnotationParameters, context: Optional[Context] = None
):
    """Prepare gene sets for sc-type"""
    if context:
        await context.info("🧬 Preparing sc-type gene sets...")

    try:
        # Get robjects from validation
        robjects, _, _, _, _, default_converter, openrlib = _validate_rpy2_and_r(
            context
        )

        if params.sctype_custom_markers:
            # Use custom markers
            if context:
                await context.info("Using custom marker gene sets")
            return _convert_custom_markers_to_gs(params.sctype_custom_markers, context)
        else:
            # Use sc-type database
            tissue = params.sctype_tissue
            if not tissue:
                raise ValueError(
                    "sctype_tissue parameter is required when not using custom markers"
                )

            if context:
                await context.info(f"Using sc-type database for tissue: {tissue}")

            # Download and use ScTypeDB
            db_path = (
                params.sctype_db_
                or "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
            )

            # Import conversion module
            from rpy2.robjects import conversion

            # Wrap R calls in conversion context (FIX for contextvars issue)
            with openrlib.rlock:
                with conversion.localconverter(default_converter):
                    robjects.r.assign("db_path", db_path)
                    robjects.r.assign("tissue_type", tissue)

                    robjects.r(
                        """
                        # Load gene sets
                        gs_list <- gene_sets_prepare(db_path, tissue_type)
                    """
                    )

                    return robjects.r["gs_list"]

    except Exception as e:
        error_msg = f"Failed to prepare gene sets: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


def _convert_custom_markers_to_gs(
    custom_markers: Dict[str, Dict[str, List[str]]], context: Optional[Context] = None
):
    """Convert custom markers to sc-type gene set format"""
    if not custom_markers:
        raise ValueError("Custom markers dictionary is empty")

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
        raise ValueError(
            "No valid cell types found in custom markers - all cell types need at least one positive marker"
        )

    # Get robjects and converters from validation
    robjects, pandas2ri, _, _, localconverter, default_converter, openrlib = (
        _validate_rpy2_and_r(context)
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


async def _run_sctype_scoring(
    adata, gs_list, params: AnnotationParameters, context: Optional[Context] = None
):
    """Run sc-type scoring algorithm"""
    if context:
        await context.info("🔬 Running sc-type scoring...")

    try:
        # Get robjects and converters from validation
        robjects, pandas2ri, _, _, localconverter, default_converter, openrlib = (
            _validate_rpy2_and_r(context)
        )

        # Import conversion module
        from rpy2.robjects import conversion

        # Prepare expression data
        if params.sctype_scaled and "scaled" in adata.layers:
            expr_data = adata.layers["scaled"]
            if context:
                await context.info(
                    "Using scaled expression data from adata.layers['scaled']"
                )
        else:
            expr_data = adata.X
            if context:
                await context.info("Using raw expression data from adata.X")

        # Convert to dense array if sparse
        if hasattr(expr_data, "toarray"):
            expr_data = expr_data.toarray()

        # Convert to DataFrame with proper gene and cell names
        expr_df = pd.DataFrame(
            expr_data.T,  # sc-type expects genes as rows, cells as columns
            index=adata.var_names,
            columns=adata.obs_names,
        )

        # Wrap ALL R calls in conversion context (FIX for contextvars issue)
        with openrlib.rlock:
            with conversion.localconverter(default_converter + pandas2ri.converter):
                # Transfer to R
                robjects.r.assign("scdata", expr_df)
                robjects.r.assign("gs_list", gs_list)

                # Extract gs_positive and gs_negative from gs_list in R
                robjects.r(
                    """
                    gs_positive <- gs_list$gs_positive
                    gs_negative <- gs_list$gs_negative
                """
                )

                # Run sc-type scoring with comprehensive error handling
                robjects.r(
                    """
                    # Check if gene sets are valid
                    if (length(gs_positive) == 0) {
                        stop("No valid positive gene sets found")
                    }

                    # Get available genes in the dataset
                    available_genes <- rownames(scdata)

                    # Check each cell type for overlapping genes and filter empty ones
                    filtered_gs_positive <- list()
                    filtered_gs_negative <- list()

                    for (celltype in names(gs_positive)) {
                        pos_genes <- gs_positive[[celltype]]
                        neg_genes <- if (celltype %in% names(gs_negative)) gs_negative[[celltype]] else c()

                        # Find overlapping genes
                        pos_overlap <- intersect(toupper(pos_genes), toupper(available_genes))
                        neg_overlap <- intersect(toupper(neg_genes), toupper(available_genes))

                        # Only keep cell types with at least one positive marker gene
                        if (length(pos_overlap) > 0) {
                            filtered_gs_positive[[celltype]] <- pos_overlap
                            filtered_gs_negative[[celltype]] <- neg_overlap
                        }
                    }

                    # Check if we have any valid cell types after filtering
                    if (length(filtered_gs_positive) == 0) {
                        # Fail explicitly when no valid gene sets are available
                        stop("No valid cell type gene sets found after filtering. Available tissues: ",
                             paste(unique(tissue_df$tissueType), collapse=", "),
                             ". Please check your tissue parameter or provide custom markers.")
                    }

                    # Run sc-type scoring with filtered gene sets
                    tryCatch({
                        es_max <- sctype_score(
                            scRNAseqData = as.matrix(scdata),
                            scaled = TRUE,
                            gs = filtered_gs_positive,
                            gs2 = filtered_gs_negative
                        )

                        # Check for valid results
                        if (is.null(es_max) || nrow(es_max) == 0 || ncol(es_max) == 0) {
                            stop("SC-Type analysis failed to generate valid results. This may be due to: ",
                                 "1) Insufficient gene overlap between data and markers, ",
                                 "2) Poor data quality, or 3) Inappropriate tissue selection.")
                        }
                    }, error = function(e) {
                        # Propagate the actual error instead of masking it
                        stop("SC-Type scoring failed: ", e$message)
                    })
                    """
                )

                # Get results back to Python with row and column names preserved
                # Extract row names (cell type names) and column names from R
                row_names = list(robjects.r("rownames(es_max)"))
                col_names = list(robjects.r("colnames(es_max)"))
                scores_matrix = robjects.r["es_max"]

        # Convert to DataFrame with proper index and columns
        if isinstance(scores_matrix, pd.DataFrame):
            scores_df = scores_matrix
            if row_names:
                scores_df.index = row_names
            if col_names:
                scores_df.columns = col_names
        else:
            scores_df = pd.DataFrame(
                scores_matrix,
                index=row_names if row_names else None,
                columns=col_names if col_names else None,
            )

        if context:
            await context.info(
                f"✅ Scoring completed - DataFrame shape: {scores_df.shape}"
            )
            await context.info(f"   Rows (should be cell types): {scores_df.shape[0]}")
            await context.info(f"   Cols (should be cells): {scores_df.shape[1]}")
            await context.info(
                f"   Expected cells: {len(col_names) if col_names else 'unknown'}"
            )
            if scores_df.index is not None and len(scores_df.index) > 0:
                cell_type_names = list(scores_df.index)[:5]
                await context.info(
                    f"📋 Cell types detected: {', '.join(str(ct) for ct in cell_type_names)}"
                    + ("..." if len(scores_df.index) > 5 else "")
                )

        return scores_df

    except Exception as e:
        error_msg = f"Failed to run sc-type scoring: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


async def _assign_sctype_celltypes(scores_df, context: Optional[Context] = None):
    """Assign cell types based on sc-type scores"""
    if context:
        await context.info("🏷️  Assigning cell types based on scores...")

    try:
        # Handle empty or invalid scores
        if scores_df is None:
            raise ValueError("Scores data is None")

        # Convert to DataFrame if it's a numpy array
        if isinstance(scores_df, np.ndarray):
            # Check for empty array
            if scores_df.size == 0:
                raise ValueError("Scores array is empty")

            # Convert to DataFrame - need to get column and row names from R
            import pandas as pd

            scores_df = pd.DataFrame(scores_df)
            if context:
                await context.info(
                    f"Converted numpy array to DataFrame: {scores_df.shape}"
                )

        # Check DataFrame has data
        if hasattr(scores_df, "empty") and scores_df.empty:
            raise ValueError("Scores DataFrame is empty")

        # Get the cell type with highest score for each cell
        cell_types = []
        confidence_scores = []

        # Handle both DataFrame and numpy array cases
        if hasattr(scores_df, "columns"):
            # DataFrame case
            if len(scores_df.columns) == 0:
                raise ValueError("DataFrame has no columns")

            n_cells = len(scores_df.columns)
            n_celltypes = len(scores_df.index)

            # Handle single cell type case
            if n_celltypes == 1:
                # Only one cell type available - assign to all cells based on scores
                single_celltype = str(scores_df.index[0])

                for col_idx in range(n_cells):
                    cell_score = scores_df.iloc[0, col_idx]
                    if cell_score > 0:
                        cell_types.append(single_celltype)
                        confidence_scores.append(min(cell_score / 10.0, 1.0))
                    else:
                        cell_types.append("Unknown")
                        confidence_scores.append(0.0)
            else:
                # Multiple cell types available

                for col_idx in range(n_cells):
                    cell_scores = scores_df.iloc[
                        :, col_idx
                    ]  # All cell types for this cell

                    # Find cell type with maximum score
                    max_score_idx = cell_scores.idxmax()
                    max_score = cell_scores.loc[max_score_idx]

                    # If max score is positive, assign cell type, otherwise "Unknown"
                    if max_score > 0:
                        cell_types.append(str(max_score_idx))
                        # Normalize confidence score to 0-1 range
                        confidence_scores.append(min(max_score / 10.0, 1.0))
                    else:
                        cell_types.append("Unknown")
                        confidence_scores.append(0.0)
        else:
            # Numpy array case
            if len(scores_df.shape) == 0 or scores_df.size == 0:
                raise ValueError("Empty scores array")

            n_cells = scores_df.shape[1] if len(scores_df.shape) > 1 else 1
            n_celltypes = (
                scores_df.shape[0] if len(scores_df.shape) > 1 else scores_df.shape[0]
            )

            if n_celltypes == 0:
                raise ValueError("No cell types in scores array")

            for cell_idx in range(n_cells):
                if len(scores_df.shape) > 1:
                    cell_scores = scores_df[:, cell_idx]
                else:
                    cell_scores = (
                        scores_df
                        if n_cells == 1
                        else scores_df[cell_idx : cell_idx + 1]
                    )

                # Handle case where cell_scores is empty
                if len(cell_scores) == 0:
                    cell_types.append("Unknown")
                    confidence_scores.append(0.0)
                    continue

                # Find cell type with maximum score
                max_score_idx = np.argmax(cell_scores)
                max_score = cell_scores[max_score_idx]

                # If max score is positive, assign cell type, otherwise "Unknown"
                if max_score > 0:
                    cell_types.append(f"CellType_{max_score_idx}")
                    # Normalize confidence score to 0-1 range
                    confidence_scores.append(min(max_score / 10.0, 1.0))
                else:
                    cell_types.append("Unknown")
                    confidence_scores.append(0.0)

        if context:
            unique_types = list(set(cell_types))
            await context.info(
                f"✅ Assigned {len(unique_types)} unique cell types: {', '.join(unique_types[:10])}"
                + ("..." if len(unique_types) > 10 else "")
            )

        return cell_types, confidence_scores

    except Exception as e:
        error_msg = f"Failed to assign cell types: {str(e)}"
        if context:
            await context.error(error_msg)
        raise ValueError(error_msg)


def _calculate_sctype_stats(cell_types):
    """Calculate statistics for sc-type results"""
    from collections import Counter

    counts = Counter(cell_types)
    return dict(counts)


async def _cache_sctype_results(
    cache_key: str, results, context: Optional[Context] = None
) -> None:
    """Cache sc-type results to disk"""
    if context:
        await context.info("💾 Caching sc-type results...")

    try:
        # Create cache directory
        _SCTYPE_CACHE_DIR.mkdir(parents=True, exist_ok=True)

        # Save to cache file
        cache_file = _SCTYPE_CACHE_DIR / f"{cache_key}.pkl"
        with open(cache_file, "wb") as f:
            pickle.dump(results, f)

        # Also store in memory cache
        _SCTYPE_CACHE[cache_key] = results

        if context:
            await context.info("✅ Results cached successfully")

    except Exception as e:
        if context:
            await context.warning(f"Failed to cache results: {str(e)}")


async def _load_cached_sctype_results(
    cache_key: str, context: Optional[Context] = None
):
    """Load cached sc-type results"""
    # Check memory cache first
    if cache_key in _SCTYPE_CACHE:
        if context:
            await context.info("📂 Using cached results from memory")
        return _SCTYPE_CACHE[cache_key]

    # Check disk cache
    cache_file = _SCTYPE_CACHE_DIR / f"{cache_key}.pkl"
    if cache_file.exists():
        try:
            with open(cache_file, "rb") as f:
                results = pickle.load(f)

            # Store in memory cache for next time
            _SCTYPE_CACHE[cache_key] = results

            if context:
                await context.info("📂 Using cached results from disk")
            return results

        except Exception as e:
            if context:
                await context.warning(f"Failed to load cached results: {str(e)}")

    return None


async def _annotate_with_sctype(
    adata: sc.AnnData, params: AnnotationParameters, context: Optional[Context] = None
) -> tuple[List[str], Dict[str, int], List[float], Optional[float]]:
    """
    Annotate cell types using sc-type method

    Args:
        adata: AnnData object
        params: Annotation parameters
        context: MCP context

    Returns:
        Tuple of (cell_types, counts, confidence_scores, mapping_score)
    """

    # Validate dependencies with comprehensive error reporting
    (
        robjects,
        pandas2ri,
        numpy2ri,
        importr,
        localconverter,
        default_converter,
        openrlib,
    ) = _validate_rpy2_and_r(context)

    # Define supported tissue types from sc-type database
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

    # Validate required parameters
    if not params.sctype_tissue and not params.sctype_custom_markers:
        raise ValueError(
            "Either sctype_tissue or sctype_custom_markers must be specified"
        )

    # Validate tissue type if provided
    if params.sctype_tissue and params.sctype_tissue not in SCTYPE_VALID_TISSUES:
        valid_tissues = sorted(SCTYPE_VALID_TISSUES)
        raise ValueError(
            f"Tissue type '{params.sctype_tissue}' is not supported by sc-type database.\n"
            f"Supported tissues: {', '.join(valid_tissues)}\n"
            f"Alternatively, use sctype_custom_markers to define custom cell type markers."
        )

    try:
        if context:
            await context.info("🧬 Starting sc-type cell type annotation...")
            await context.info(
                f"📊 Analyzing {adata.n_obs} cells with {adata.n_vars} genes"
            )
            if params.sctype_tissue:
                await context.info(f"🔬 Using tissue type: {params.sctype_tissue}")
            else:
                await context.info(
                    f"🔬 Using custom markers for {len(params.sctype_custom_markers)} cell types"
                )

        # Check cache if enabled
        cache_key = None
        if params.sctype_use_cache:
            cache_key = _get_sctype_cache_key(adata, params)
            cached_results = await _load_cached_sctype_results(cache_key, context)
            if cached_results:
                return cached_results

        # Step 1: Load sc-type functions
        await _load_sctype_functions(context)

        # Step 2: Prepare gene sets
        gs_list = await _prepare_sctype_genesets(params, context)

        # Step 3: Run sc-type scoring
        scores_df = await _run_sctype_scoring(adata, gs_list, params, context)

        # Step 4: Assign cell types
        cell_types, confidence_scores = await _assign_sctype_celltypes(
            scores_df, context
        )

        # Step 5: Calculate statistics
        counts = _calculate_sctype_stats(cell_types)

        # Step 6: Calculate average confidence scores per cell type (to match AnnotationResult interface)
        confidence_by_celltype = {}
        for cell_type in set(cell_types):
            # Get confidence scores for this cell type
            celltype_confidences = [
                conf
                for i, conf in enumerate(confidence_scores)
                if cell_types[i] == cell_type
            ]
            if celltype_confidences:
                confidence_by_celltype[cell_type] = sum(celltype_confidences) / len(
                    celltype_confidences
                )
            else:
                confidence_by_celltype[cell_type] = 0.0

        # ✅ COW FIX: Assign results to adata.obs (like tangram and other methods)
        # This maintains consistency across all annotation methods
        output_key = f"cell_type_{params.method}"
        confidence_key = f"confidence_{params.method}"

        adata.obs[output_key] = cell_types
        adata.obs[output_key] = adata.obs[output_key].astype("category")
        adata.obs[confidence_key] = confidence_scores

        # Return unique cell types and per-cell-type confidence dict (for AnnotationResult)
        unique_cell_types = list(set(cell_types))
        results = (unique_cell_types, counts, confidence_by_celltype, None)

        # Cache results if enabled
        if params.sctype_use_cache and cache_key:
            await _cache_sctype_results(cache_key, results, context)

        if context:
            await context.info("🎉 sc-type annotation completed successfully!")
            await context.info(
                f"📈 Found {len(counts)} cell types: {', '.join(list(counts.keys())[:5])}"
                + ("..." if len(counts) > 5 else "")
            )

        return results

    except Exception as e:
        await _handle_annotation_error(e, "sctype", context)
