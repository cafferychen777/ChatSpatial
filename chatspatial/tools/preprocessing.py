"""
Preprocessing tools for spatial transcriptomics data.
"""

import logging
import traceback
from typing import Any, Dict, Optional

import numpy as np
import scanpy as sc
import scipy.sparse
import squidpy as sq
from anndata import AnnData
from mcp.server.fastmcp import Context

from ..models.analysis import PreprocessingResult
from ..models.data import (PreprocessingParameters,
                           ResolVIPreprocessingParameters)
from ..utils.data_adapter import standardize_adata
from ..utils.tool_error_handling import mcp_tool_error_handler

# Import scvi-tools for advanced preprocessing
try:
    import scvi
    import torch
    from scvi.external import RESOLVI
except ImportError:
    torch = None
    scvi = None
    RESOLVI = None

# Import scvelo for RNA velocity preprocessing
try:
    import scvelo as scv
except ImportError:
    scv = None

# Import SpaGCN for spatial domain-specific preprocessing
try:
    import SpaGCN as spg
except ImportError:
    spg = None

# Check for R sctransform availability via rpy2
def _is_r_sctransform_available() -> tuple[bool, str]:
    """Check if R sctransform package is available via rpy2.

    Returns:
        Tuple of (is_available, error_message)
    """
    try:
        import rpy2.robjects as ro
        import rpy2.robjects.packages as rpackages
        from rpy2.robjects.conversion import localconverter

        with localconverter(ro.default_converter):
            rpackages.importr("sctransform")
            rpackages.importr("Matrix")
        return True, ""
    except ImportError:
        return False, (
            "rpy2 is not installed. Install with: pip install 'rpy2>=3.5.0'"
        )
    except Exception as e:
        if "sctransform" in str(e).lower():
            return False, (
                "R package 'sctransform' is not installed.\n"
                "Install in R: install.packages('sctransform')"
            )
        elif "matrix" in str(e).lower():
            return False, (
                "R package 'Matrix' is not installed.\n"
                "Install in R: install.packages('Matrix')"
            )
        return False, f"R sctransform check failed: {e}"


# Lazy check - only validate when actually used
R_SCTRANSFORM_AVAILABLE: Optional[bool] = None

# Setup logger
logger = logging.getLogger(__name__)
MIN_KMEANS_CLUSTERS = 2
MAX_TSNE_PCA_COMPONENTS = 50


def _should_use_all_genes_for_hvg(adata) -> bool:
    """
    Check if we should use all genes for HVG selection.

    Only applies to very small gene sets (e.g., MERFISH with <100 genes)
    where statistical HVG selection is not meaningful.
    """
    return adata.n_vars < 100


@mcp_tool_error_handler()
async def preprocess_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: PreprocessingParameters = PreprocessingParameters(),
    context: Optional[Context] = None,
) -> PreprocessingResult:
    """Preprocess spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Preprocessing parameters
        context: MCP context

    Returns:
        Preprocessing result summary
    """
    try:
        if context:
            await context.info(
                f"Preprocessing dataset {data_id} with {params.normalization} normalization"
            )

        # Retrieve the AnnData object from data store
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found in data store")

        # Memory optimization: Preprocessing modifies dataset in-place (by design)
        # The preprocessed result replaces the original at line 982, so no copy needed
        # Original counts are preserved in adata.raw (line 311-316)
        adata = data_store[data_id]["adata"]

        # LINUS FIX: Standardize data format at the entry point
        # This eliminates all downstream special cases for data format handling
        if context:
            await context.info("Standardizing data structure to ChatSpatial format...")
        try:
            adata = standardize_adata(
                adata, copy=False, strict=False, preserve_original=True
            )
            if context:
                await context.info("Data structure standardized successfully")
        except Exception as e:
            if context:
                await context.warning(
                    f"Data standardization failed: {e}. Proceeding with original data."
                )
            # Continue with original data if standardization fails

        # Validate input data
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError(
                f"Dataset {data_id} is empty: {adata.n_obs} cells, {adata.n_vars} genes"
            )

        # ===== Handle Duplicate Gene Names (CRITICAL FIX) =====
        # Must be done BEFORE any gene-based operations (QC, HVG selection, etc.)
        if not adata.var_names.is_unique:
            n_duplicates = len(adata.var_names) - len(set(adata.var_names))
            if context:
                await context.warning(
                    f"Found {n_duplicates} duplicate gene names in data"
                )
                await context.info(
                    "Fixing duplicate gene names with unique suffixes..."
                )
            adata.var_names_make_unique()
            if context:
                await context.info(f"Fixed {n_duplicates} duplicate gene names")

        # 1. Calculate QC metrics
        if context:
            await context.info("Calculating QC metrics...")
        try:
            sc.pp.calculate_qc_metrics(adata, inplace=True)
        except Exception as e:
            # Don't fallback with fake data - fail fast with actionable error message
            # calculate_qc_metrics failures are typically dependency issues, not data issues
            error_msg = (
                f"QC metrics calculation failed: {str(e)}\n\n"
                "COMMON CAUSES:\n"
                "1. Numba/llvmlite version incompatibility\n"
                "   → Fix: pip install --upgrade numba llvmlite\n"
                "2. Corrupted or empty expression matrix\n"
                "   → Check: adata.X should contain valid numeric data\n"
                "3. Memory issues with large sparse matrices\n"
                "   → Try: Reduce dataset size or increase available memory\n\n"
                f"DATASET INFO: {adata.n_obs} cells × {adata.n_vars} genes\n"
                f"MATRIX TYPE: {type(adata.X).__name__}\n\n"
                "SCIENTIFIC INTEGRITY: We refuse to create fake QC metrics.\n"
                "Please fix the underlying issue before proceeding."
            )
            if context:
                await context.error(error_msg)
            raise RuntimeError(error_msg) from e

        # Store original QC metrics before filtering
        qc_metrics = {
            "n_cells_before_filtering": int(adata.n_obs),
            "n_genes_before_filtering": int(adata.n_vars),
            "median_genes_per_cell": float(np.median(adata.obs.n_genes_by_counts)),
            "median_umi_per_cell": float(np.median(adata.obs.total_counts)),
        }

        # 2. Apply user-controlled data filtering and subsampling
        if context:
            await context.info("Applying data filtering and subsampling...")

        # Apply gene filtering using LLM-controlled parameters
        min_cells = params.filter_genes_min_cells
        if min_cells is not None and min_cells > 0:
            if context:
                await context.info(f"Filtering genes: min_cells={min_cells}")
            sc.pp.filter_genes(adata, min_cells=min_cells)

        # Apply cell filtering using LLM-controlled parameters
        min_genes = params.filter_cells_min_genes
        if min_genes is not None and min_genes > 0:
            if context:
                await context.info(f"Filtering cells: min_genes={min_genes}")
            sc.pp.filter_cells(adata, min_genes=min_genes)

        # Apply spot subsampling if requested
        if params.subsample_spots is not None and params.subsample_spots < adata.n_obs:
            if context:
                await context.info(
                    f"Subsampling from {adata.n_obs} to {params.subsample_spots} spots"
                )
            sc.pp.subsample(
                adata,
                n_obs=params.subsample_spots,
                random_state=params.subsample_random_seed,
            )

        # Apply gene subsampling if requested (after HVG selection)
        gene_subsample_requested = params.subsample_genes is not None

        # Save raw data before normalization (required for some analysis methods)
        if context:
            await context.info("Saving raw data for downstream analysis...")

        # IMPORTANT: Create a proper frozen copy for .raw to preserve counts
        # Using `adata.raw = adata` creates a view that gets modified during normalization
        # We need to create an independent AnnData object to truly preserve counts
        import anndata as ad_module

        # Memory optimization: AnnData.raw internally copies var, so no need for .copy()
        # obs MUST be copied to prevent contamination from later preprocessing steps
        # uns can be empty dict as raw doesn't need metadata
        adata.raw = ad_module.AnnData(
            X=adata.X.copy(),  # Must copy - will be modified during normalization
            var=adata.var,  # No copy needed - AnnData internally creates independent copy
            obs=adata.obs.copy(),  # Must copy - will be modified by clustering/annotation
            uns={},  # Empty dict - raw doesn't need uns metadata
        )

        # Update QC metrics after filtering
        qc_metrics.update(
            {
                "n_cells_after_filtering": int(adata.n_obs),
                "n_genes_after_filtering": int(adata.n_vars),
            }
        )

        # 3. Advanced preprocessing with ResolVI (if enabled)
        resolvi_used = False
        if params.use_resolvi_preprocessing and params.resolvi_params:
            if context:
                await context.info(
                    "Using ResolVI for advanced molecular reassignment correction..."
                )

            try:
                # Run ResolVI preprocessing
                await preprocess_with_resolvi(adata, params.resolvi_params, context)

                # Verify that ResolVI actually created the expected outputs
                if "X_resolvi" not in adata.obsm:
                    raise RuntimeError(
                        "ResolVI preprocessing completed but did not create X_resolvi. "
                        "This indicates a critical error in ResolVI execution."
                    )

                # ResolVI creates a latent representation in adata.obsm['X_resolvi']
                # and corrected counts in adata.layers['resolvi_corrected']
                if context:
                    await context.info("ResolVI preprocessing completed successfully")
                    await context.info(
                        f"Created latent representation with shape {adata.obsm['X_resolvi'].shape}"
                    )

                # Skip standard normalization since ResolVI handles it
                resolvi_used = True
            except Exception as e:
                if context:
                    await context.warning(
                        f"ResolVI preprocessing failed: {e}. Falling back to standard normalization."
                    )
                resolvi_used = False

        # 3. Normalize data (skip if ResolVI was used)
        if not resolvi_used:
            # Log normalization configuration
            logger.info("=" * 50)
            logger.info("Normalization Configuration:")
            logger.info(f"  Method: {params.normalization}")
            if params.normalize_target_sum is not None:
                logger.info(f"  Target sum: {params.normalize_target_sum:.0f}")
            else:
                logger.info("  Target sum: ADAPTIVE (using median counts)")
            if params.scale:
                if params.scale_max_value is not None:
                    logger.info(f"  Scale clipping: ±{params.scale_max_value} SD")
                else:
                    logger.info("  Scale clipping: NONE (preserving all outliers)")
            logger.info("=" * 50)

            if context:
                await context.info(
                    f"Normalizing data using {params.normalization} method..."
                )

            if params.normalization == "log":
                # Standard log normalization
                # Check if data appears to be already normalized
                if scipy.sparse.issparse(adata.X):
                    X_sample = (
                        adata.X.data[: min(1000, len(adata.X.data))]
                        if hasattr(adata.X, "data")
                        else adata.X[:1000].toarray().flatten()
                    )
                else:
                    X_sample = adata.X.flatten()[: min(1000, adata.X.size)]

                # Check for negative values (indicates already log-normalized data)
                if np.any(X_sample < 0):
                    error_msg = (
                        "Log normalization requires non-negative data (raw or normalized counts). "
                        "Data contains negative values, suggesting it has already been log-normalized. "
                        "Options:\n"
                        "• Use normalization='none' if data is already pre-processed\n"
                        "• Load raw count data instead of processed data\n"
                        "• Remove the log transformation from your data before re-processing"
                    )
                    if context:
                        await context.error(error_msg)
                    raise ValueError(error_msg)

                if params.normalize_target_sum is not None:
                    sc.pp.normalize_total(adata, target_sum=params.normalize_target_sum)
                    logger.info(
                        f"✓ Normalized to target_sum={params.normalize_target_sum:.0f}"
                    )
                else:
                    # Calculate median and inform user transparently
                    calculated_median = np.median(
                        np.array(adata.X.sum(axis=1)).flatten()
                    )
                    if context:
                        await context.info(
                            f"normalize_target_sum not specified. Using adaptive normalization:\n"
                            f"   • Calculated median counts: {calculated_median:.0f}\n"
                            f"   • This value was automatically determined from your data\n"
                            f"   • For reproducible results, use: normalize_target_sum={calculated_median:.0f}"
                        )
                    sc.pp.normalize_total(adata, target_sum=calculated_median)
                    logger.info(
                        f"✓ Normalized to adaptive target_sum={calculated_median:.0f}"
                    )
                sc.pp.log1p(adata)
                logger.info("Applied log1p transformation")
            elif params.normalization == "sct":
                # SCTransform v2 variance-stabilizing normalization via R's sctransform
                # Check R sctransform availability
                is_available, error_msg = _is_r_sctransform_available()
                if not is_available:
                    full_error = (
                        f"SCTransform requires R and the sctransform package.\n\n"
                        f"ERROR: {error_msg}\n\n"
                        "INSTALLATION:\n"
                        "  1. Install R (https://cran.r-project.org/)\n"
                        "  2. In R: install.packages('sctransform')\n"
                        "  3. pip install 'rpy2>=3.5.0'\n\n"
                        "ALTERNATIVES:\n"
                        "• Use normalization='pearson_residuals' (built-in, similar results)\n"
                        "• Use normalization='log' (standard method)"
                    )
                    if context:
                        await context.error(full_error)
                    raise ImportError(full_error)

                # Check if data appears to be raw counts (required for SCTransform)
                if scipy.sparse.issparse(adata.X):
                    X_sample = (
                        adata.X.data[: min(1000, len(adata.X.data))]
                        if hasattr(adata.X, "data")
                        else adata.X[:1000].toarray().flatten()
                    )
                else:
                    X_sample = adata.X.flatten()[: min(1000, adata.X.size)]

                # Check for non-integer values (indicates normalized data)
                if np.any((X_sample % 1) != 0):
                    error_msg = (
                        "SCTransform requires raw count data (integers).\n"
                        "Your data contains non-integer values.\n\n"
                        "SOLUTIONS:\n"
                        "• Load raw count data instead of normalized data\n"
                        "• Use normalization='none' if data is pre-normalized\n"
                        "• Use normalization='log' for already-processed data"
                    )
                    if context:
                        await context.error(error_msg)
                    raise ValueError(error_msg)

                # Map method parameter to vst.flavor
                vst_flavor = "v2" if params.sct_method == "fix-slope" else "v1"
                if context:
                    await context.info(
                        f"Applying SCTransform (vst.flavor={vst_flavor}, "
                        f"var_features={params.sct_var_features_n})"
                    )

                try:
                    # Import rpy2 modules
                    import rpy2.robjects as ro
                    from rpy2.robjects import numpy2ri
                    from rpy2.robjects.conversion import localconverter

                    # Store raw counts before replacing X
                    adata.layers["counts"] = adata.X.copy()

                    # Convert to sparse CSC matrix (genes × cells) for R's dgCMatrix
                    if scipy.sparse.issparse(adata.X):
                        counts_sparse = scipy.sparse.csc_matrix(adata.X.T)
                    else:
                        counts_sparse = scipy.sparse.csc_matrix(adata.X.T)

                    # Transfer sparse matrix components to R
                    with localconverter(ro.default_converter + numpy2ri.converter):
                        ro.globalenv["sp_data"] = counts_sparse.data.astype(np.float64)
                        ro.globalenv["sp_indices"] = counts_sparse.indices.astype(
                            np.int32
                        )
                        ro.globalenv["sp_indptr"] = counts_sparse.indptr.astype(
                            np.int32
                        )
                        ro.globalenv["n_genes"] = counts_sparse.shape[0]
                        ro.globalenv["n_cells"] = counts_sparse.shape[1]
                        ro.globalenv["gene_names"] = ro.StrVector(
                            adata.var_names.tolist()
                        )
                        ro.globalenv["cell_names"] = ro.StrVector(
                            adata.obs_names.tolist()
                        )
                        ro.globalenv["vst_flavor"] = vst_flavor
                        ro.globalenv["n_cells_param"] = (
                            params.sct_n_cells if params.sct_n_cells else ro.NULL
                        )

                    # Reconstruct sparse matrix and run SCTransform in R
                    ro.r(
                        """
                        library(Matrix)
                        library(sctransform)

                        # Create dgCMatrix from components
                        umi_matrix <- new(
                            "dgCMatrix",
                            x = as.numeric(sp_data),
                            i = as.integer(sp_indices),
                            p = as.integer(sp_indptr),
                            Dim = as.integer(c(n_genes, n_cells)),
                            Dimnames = list(gene_names, cell_names)
                        )

                        # Run SCTransform
                        suppressWarnings({
                            vst_result <- sctransform::vst(
                                umi = umi_matrix,
                                vst.flavor = vst_flavor,
                                return_gene_attr = TRUE,
                                return_cell_attr = TRUE,
                                n_cells = n_cells_param,
                                verbosity = 0
                            )
                        })

                        # Convert output to dense matrix for transfer
                        pearson_residuals <- as.matrix(vst_result$y)
                        residual_variance <- vst_result$gene_attr$residual_variance
                    """
                    )

                    # Extract results from R
                    with localconverter(ro.default_converter + numpy2ri.converter):
                        pearson_residuals = np.array(ro.r("pearson_residuals"))
                        residual_variance = np.array(ro.r("residual_variance"))

                    # Transpose back to cells × genes for AnnData format
                    adata.X = pearson_residuals.T

                    # Store SCTransform metadata
                    adata.uns["sctransform"] = {
                        "method": params.sct_method,
                        "vst_flavor": vst_flavor,
                        "var_features_n": params.sct_var_features_n,
                        "exclude_poisson": params.sct_exclude_poisson,
                        "n_cells": params.sct_n_cells,
                    }

                    # Mark highly variable genes based on residual variance
                    adata.var["sct_residual_variance"] = residual_variance

                    # Select top N genes by residual variance
                    n_hvg = min(params.sct_var_features_n, len(residual_variance))
                    top_hvg_indices = np.argsort(residual_variance)[-n_hvg:]
                    adata.var["highly_variable"] = False
                    adata.var.iloc[top_hvg_indices, adata.var.columns.get_loc(
                        "highly_variable"
                    )] = True

                    if context:
                        await context.info(
                            f"✓ SCTransform: {n_hvg} highly variable genes identified"
                        )

                    logger.info("Applied SCTransform normalization via R")

                except MemoryError:
                    error_msg = (
                        f"Insufficient memory for SCTransform on "
                        f"{adata.n_obs}×{adata.n_vars} matrix.\n"
                        "SOLUTIONS:\n"
                        "• Reduce dataset size with subsample_spots parameter\n"
                        "• Use normalization='pearson_residuals' (more memory efficient)\n"
                        "• Use normalization='log' (minimal memory usage)"
                    )
                    if context:
                        await context.error(error_msg)
                    raise MemoryError(error_msg)
                except Exception as e:
                    error_msg = f"SCTransform failed: {str(e)}"
                    if context:
                        await context.error(error_msg)
                        await context.info(f"Error details: {traceback.format_exc()}")
                    raise RuntimeError(error_msg) from e
            elif params.normalization == "pearson_residuals":
                # Modern Pearson residuals normalization (recommended for UMI data)
                if context:
                    await context.info("Using Pearson residuals normalization...")

                # Check if method is available
                if not hasattr(sc.experimental.pp, "normalize_pearson_residuals"):
                    error_msg = (
                        "Pearson residuals normalization not available (requires scanpy>=1.9.0).\n"
                        "Options:\n"
                        "• Install newer scanpy: pip install 'scanpy>=1.9.0'\n"
                        "• Use log normalization instead: params.normalization='log'\n"
                        "• Skip normalization if data is pre-processed: params.normalization='none'"
                    )
                    if context:
                        await context.error(error_msg)
                    raise ValueError(error_msg)

                # Check if data appears to be raw counts
                if scipy.sparse.issparse(adata.X):
                    # Sample first 1000 values for efficiency
                    X_sample = (
                        adata.X.data[: min(1000, len(adata.X.data))]
                        if hasattr(adata.X, "data")
                        else adata.X[:1000].toarray().flatten()
                    )
                else:
                    X_sample = adata.X.flatten()[: min(1000, adata.X.size)]

                # Check for non-integer values (indicates normalized data)
                if np.any((X_sample % 1) != 0):
                    error_msg = (
                        "Pearson residuals requires raw count data (integers). "
                        "Data contains non-integer values. "
                        "Use params.normalization='none' if data is already normalized, "
                        "or params.normalization='log' for standard normalization."
                    )
                    if context:
                        await context.error(error_msg)
                    raise ValueError(error_msg)

                # Execute normalization
                try:
                    # Apply Pearson residuals normalization (to all genes)
                    # Note: High variable gene selection happens later in the pipeline
                    sc.experimental.pp.normalize_pearson_residuals(adata)
                    logger.info("Applied Pearson residuals normalization")
                except MemoryError:
                    raise MemoryError(
                        f"Insufficient memory for Pearson residuals on {adata.n_obs}×{adata.n_vars} matrix. "
                        "Try reducing n_hvgs or use 'log' normalization."
                    )
                except Exception as e:
                    raise RuntimeError(
                        f"Pearson residuals normalization failed: {e}. "
                        "Consider using 'log' normalization instead."
                    )
            elif params.normalization == "none":
                # Explicitly skip normalization
                if context:
                    await context.info(
                        "Skipping normalization (data assumed to be pre-normalized)"
                    )
                logger.info("Normalization skipped (using pre-normalized data)")

                # CRITICAL: Check if data appears to be raw counts
                # HVG selection requires normalized data for statistical validity
                if scipy.sparse.issparse(adata.X):
                    X_sample = (
                        adata.X.data[: min(1000, len(adata.X.data))]
                        if hasattr(adata.X, "data")
                        else adata.X[:1000].toarray().flatten()
                    )
                else:
                    X_sample = adata.X.flatten()[: min(1000, adata.X.size)]

                # Check if data looks raw (all integers and high values)
                if np.all((X_sample % 1) == 0) and np.max(X_sample) > 100:
                    error_msg = (
                        "STATISTICAL ERROR: Cannot perform HVG selection on raw counts with normalization='none'\n\n"
                        "Your data appears to be raw counts (integer values with max > 100), but you specified "
                        "normalization='none'. Highly variable gene (HVG) selection requires normalized data "
                        "for statistical validity because:\n"
                        "• Raw count variance scales non-linearly with expression level\n"
                        "• This prevents accurate comparison of variability across genes\n"
                        "• Scanpy's HVG algorithm will fail with 'infinity' errors\n\n"
                        "REQUIRED ACTIONS:\n"
                        "Option 1 (Recommended): Use normalization='log' for standard log-normalization\n"
                        "Option 2: Use normalization='pearson_residuals' for variance-stabilizing normalization\n"
                        "Option 3: Pre-normalize your data externally, then reload with normalized values\n\n"
                        "WARNING: If your data is already normalized but appears raw, verify data integrity."
                    )
                    if context:
                        await context.error(error_msg)
                    raise ValueError(error_msg)
            elif params.normalization == "scvi":
                # scVI deep learning-based normalization
                # Uses variational autoencoder to learn latent representation
                if scvi is None:
                    error_msg = (
                        "scVI normalization requires scvi-tools package.\n\n"
                        "INSTALLATION:\n"
                        "  pip install scvi-tools\n\n"
                        "ALTERNATIVES:\n"
                        "• Use normalization='log' (standard method)\n"
                        "• Use normalization='pearson_residuals' (variance-stabilizing)"
                    )
                    if context:
                        await context.error(error_msg)
                    raise ImportError(error_msg)

                # Check if data appears to be raw counts (required for scVI)
                if scipy.sparse.issparse(adata.X):
                    X_sample = (
                        adata.X.data[: min(1000, len(adata.X.data))]
                        if hasattr(adata.X, "data")
                        else adata.X[:1000].toarray().flatten()
                    )
                else:
                    X_sample = adata.X.flatten()[: min(1000, adata.X.size)]

                # Check for negative values (indicates already normalized data)
                if np.any(X_sample < 0):
                    error_msg = (
                        "scVI requires non-negative count data.\n"
                        "Data contains negative values, suggesting log-normalization.\n\n"
                        "SOLUTIONS:\n"
                        "• Load raw count data instead of normalized data\n"
                        "• Use normalization='none' if data is pre-normalized"
                    )
                    if context:
                        await context.error(error_msg)
                    raise ValueError(error_msg)

                if context:
                    await context.info(
                        f"Applying scVI normalization "
                        f"(n_latent={params.scvi_n_latent}, "
                        f"n_hidden={params.scvi_n_hidden}, "
                        f"gene_likelihood={params.scvi_gene_likelihood})"
                    )

                try:
                    # Store raw counts before any modification
                    adata.layers["counts"] = adata.X.copy()

                    # Setup AnnData for scVI
                    # Use counts layer since we just saved it
                    scvi.model.SCVI.setup_anndata(
                        adata,
                        layer="counts",
                        batch_key=(
                            params.batch_key
                            if params.batch_key in adata.obs.columns
                            else None
                        ),
                    )

                    # Create scVI model with user-specified parameters
                    scvi_model = scvi.model.SCVI(
                        adata,
                        n_hidden=params.scvi_n_hidden,
                        n_latent=params.scvi_n_latent,
                        n_layers=params.scvi_n_layers,
                        dropout_rate=params.scvi_dropout_rate,
                        gene_likelihood=params.scvi_gene_likelihood,
                    )

                    if context:
                        await context.info("Training scVI model...")

                    # Train the model
                    # Use reasonable defaults for training
                    scvi_model.train(
                        max_epochs=400,
                        early_stopping=True,
                        early_stopping_patience=20,
                        early_stopping_monitor="elbo_validation",
                        train_size=0.9,
                    )

                    if context:
                        await context.info("scVI training completed")

                    # Get latent representation (replaces PCA)
                    adata.obsm["X_scvi"] = scvi_model.get_latent_representation()

                    # Get normalized expression for downstream analysis
                    # This is the denoised, batch-corrected expression
                    normalized_expr = scvi_model.get_normalized_expression(
                        library_size=1e4  # Normalize to 10k counts
                    )
                    # Store as dense array (normalized expression is typically dense)
                    if hasattr(normalized_expr, "values"):
                        adata.X = normalized_expr.values
                    else:
                        adata.X = np.array(normalized_expr)

                    # Apply log1p for downstream compatibility
                    adata.X = np.log1p(adata.X)

                    # Store scVI metadata
                    adata.uns["scvi"] = {
                        "n_hidden": params.scvi_n_hidden,
                        "n_latent": params.scvi_n_latent,
                        "n_layers": params.scvi_n_layers,
                        "dropout_rate": params.scvi_dropout_rate,
                        "gene_likelihood": params.scvi_gene_likelihood,
                        "training_completed": True,
                    }

                    if context:
                        await context.info(
                            f"✓ scVI: Latent representation stored in X_scvi "
                            f"(shape: {adata.obsm['X_scvi'].shape})"
                        )
                        await context.info(
                            "✓ scVI: Normalized expression stored in adata.X"
                        )

                    logger.info("Applied scVI normalization")

                except Exception as e:
                    error_msg = f"scVI normalization failed: {str(e)}"
                    if context:
                        await context.error(error_msg)
                        await context.info(f"Error details: {traceback.format_exc()}")
                    raise RuntimeError(error_msg) from e
            else:
                # Catch unknown normalization methods
                valid_methods = ["log", "sct", "pearson_residuals", "none", "scvi"]
                raise ValueError(
                    f"Unknown normalization method: '{params.normalization}'. "
                    f"Valid options are: {', '.join(valid_methods)}"
                )

        # 4. Find highly variable genes and apply gene subsampling
        if context:
            await context.info("Finding highly variable genes...")

        # Determine number of HVGs to select
        if gene_subsample_requested:
            # User wants to subsample genes
            n_hvgs = min(params.subsample_genes, adata.n_vars - 1, params.n_hvgs)
            if context:
                await context.info(
                    f"User requested {params.subsample_genes} genes, selecting {n_hvgs} highly variable genes"
                )
        else:
            # Use standard HVG selection
            n_hvgs = min(params.n_hvgs, adata.n_vars - 1)
            if context:
                await context.info(f"Using {n_hvgs} highly variable genes...")

        # Statistical warning: Very low HVG count may lead to unstable clustering
        # Based on literature consensus: 500-5000 genes recommended, 1000-2000 typical
        # References:
        # - Bioconductor OSCA: "any value from 500 to 5000 is reasonable"
        # - Single-cell best practices: typical range 1000-2000
        if n_hvgs < 500 and context:
            await context.warning(
                f"Using only {n_hvgs} HVGs is below the recommended minimum of 500 genes.\n"
                f"   • Literature consensus: 500-5000 genes (typical: 1000-2000)\n"
                f"   • Low gene counts may lead to unstable clustering results\n"
                f"   • Recommended: Use n_hvgs=1000-2000 for most analyses\n"
                f"   • Current dataset: {adata.n_obs} cells × {adata.n_vars} total genes"
            )

        # Check if we should use all genes (for very small gene sets)
        if _should_use_all_genes_for_hvg(adata):
            if context:
                await context.info(
                    f"Small gene set detected ({adata.n_vars} genes), using all genes for analysis"
                )
            adata.var["highly_variable"] = True
        else:
            # Attempt HVG selection - no fallback for failures
            try:
                sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs)
                logger.info(f"Selected {n_hvgs} highly variable genes")
            except Exception as e:
                # Provide clear error message without overly complex diagnostics
                # scanpy's highly_variable_genes is robust - if it fails, the data has issues
                error_msg = (
                    f"Highly variable gene selection failed: {e}\n\n"
                    "COMMON CAUSES:\n"
                    "1. Data not properly normalized\n"
                    "   → Ensure normalization completed successfully before HVG selection\n"
                    "2. Too few genes remaining after filtering\n"
                    f"   → Current: {adata.n_vars} genes, requested: {n_hvgs} HVGs\n"
                    "3. All genes have zero or constant expression\n"
                    "   → Check data quality and filtering parameters\n\n"
                    f"DATASET INFO: {adata.n_obs} cells × {adata.n_vars} genes\n\n"
                    "SOLUTIONS:\n"
                    "• If n_hvgs > n_genes: reduce n_hvgs parameter\n"
                    "• Check that normalization step completed without errors\n"
                    "• Verify input data contains meaningful expression variation"
                )
                if context:
                    await context.error(error_msg)
                raise RuntimeError(error_msg) from e

        # Apply gene subsampling if requested
        if gene_subsample_requested and params.subsample_genes < adata.n_vars:
            # Ensure HVG selection was successful
            if "highly_variable" not in adata.var:
                error_msg = (
                    "Gene subsampling requested but no highly variable genes were identified. "
                    "This indicates a failure in the HVG selection step."
                )
                if context:
                    await context.error(error_msg)
                raise RuntimeError(error_msg)

            if not adata.var["highly_variable"].any():
                error_msg = (
                    "Gene subsampling requested but no genes were marked as highly variable. "
                    "Check HVG selection parameters or data quality."
                )
                if context:
                    await context.error(error_msg)
                raise RuntimeError(error_msg)

            # Use properly identified HVGs
            adata = adata[:, adata.var["highly_variable"]].copy()
            if context:
                await context.info(
                    f"Subsampled to {adata.n_vars} highly variable genes"
                )
            logger.info(f"Gene subsampling: kept {adata.n_vars} HVGs")

        # 5. Batch effect correction (if applicable)
        if (
            params.batch_key in adata.obs
            and len(adata.obs[params.batch_key].unique()) > 1
        ):
            if context:
                await context.info(
                    "Detected batch information. Applying batch effect correction with Harmony..."
                )
            try:
                # Use Harmony for batch correction (modern standard, works on PCA space)
                # Harmony is more robust than ComBat for single-cell/spatial data
                import scanpy.external as sce

                # Harmony requires PCA to be computed first
                if "X_pca" not in adata.obsm:
                    sc.tl.pca(adata, n_comps=min(50, adata.n_vars - 1))

                sce.pp.harmony_integrate(adata, key=params.batch_key)
                if context:
                    await context.info(
                        "Batch effect correction completed using Harmony"
                    )
            except ImportError:
                if context:
                    await context.warning(
                        "Harmony not available (pip install harmonypy). "
                        "Skipping batch correction."
                    )
            except Exception as e:
                if context:
                    await context.warning(
                        f"Batch effect correction failed: {e}. "
                        "Continuing without correction."
                    )

        # 6. Scale data (if requested)
        if params.scale:
            if context:
                await context.info("Scaling data...")
            try:
                # Trust scanpy's internal zero-variance handling and sparse matrix optimization
                sc.pp.scale(adata, max_value=params.scale_max_value)

                # Log scaling results
                if params.scale_max_value is not None:
                    logger.info(
                        f"✓ Scaled to unit variance with clipping at ±{params.scale_max_value} SD"
                    )
                else:
                    logger.info(
                        "✓ Scaled to unit variance without clipping (preserving all outliers)"
                    )

                # Clean up any NaN/Inf values that might remain (sparse-matrix safe)
                # Only apply if we have a max_value for clipping
                if params.scale_max_value is not None:
                    if hasattr(adata.X, "data"):
                        # Sparse matrix - only modify the data array
                        adata.X.data = np.nan_to_num(
                            adata.X.data,
                            nan=0.0,
                            posinf=params.scale_max_value,
                            neginf=-params.scale_max_value,
                        )
                    else:
                        # Dense matrix
                        adata.X = np.nan_to_num(
                            adata.X,
                            nan=0.0,
                            posinf=params.scale_max_value,
                            neginf=-params.scale_max_value,
                        )

            except Exception as e:
                if context:
                    await context.warning(
                        f"Scaling failed: {e}. Continuing without scaling."
                    )

        # 7. Run PCA (skip if ResolVI or scVI was used)
        scvi_used = params.normalization == "scvi" and "X_scvi" in adata.obsm
        if resolvi_used:
            # ResolVI already created a latent representation
            if context:
                await context.info("Skipping PCA (using ResolVI latent representation)")
            # Use ResolVI latent dimensions for downstream analysis
            n_pcs = (
                params.resolvi_params.n_latent
                if params.resolvi_params
                else params.n_pcs
            )
        elif scvi_used:
            # scVI already created a latent representation
            if context:
                await context.info("Skipping PCA (using scVI latent representation)")
            # Use scVI latent dimensions for downstream analysis
            n_pcs = params.scvi_n_latent
        else:
            if context:
                await context.info("Running PCA...")

            # Adjust n_pcs based on dataset size
            n_pcs = min(
                params.n_pcs, adata.n_vars - 1, adata.n_obs - 1
            )  # Ensure we don't use more PCs than possible

            if context:
                await context.info(f"Using {n_pcs} principal components...")

            try:
                sc.tl.pca(adata, n_comps=n_pcs)
            except Exception as e:
                # PCA failed - provide detailed error information without arbitrary fallback
                error_msg = (
                    f"PCA computation failed with {n_pcs} components: {str(e)}. "
                    f"This indicates a potential data quality issue. "
                    f"Current dataset: {adata.n_obs} cells × {adata.n_vars} genes. "
                    f"Possible solutions: "
                    f"1) Reduce n_pcs parameter manually (current: {n_pcs}) "
                    f"2) Check data normalization and filtering "
                    f"3) Ensure sufficient gene expression variation "
                    f"4) Check for memory limitations with large datasets"
                )
                if context:
                    await context.error(error_msg)
                raise RuntimeError(error_msg)

        # 8. Compute neighbors graph
        if context:
            await context.info("Computing neighbors graph...")

        # Use scientifically-informed n_neighbors with validation
        n_neighbors = params.n_neighbors

        # Validate n_neighbors against dataset constraints
        if n_neighbors >= adata.n_obs:
            raise ValueError(
                f"n_neighbors ({n_neighbors}) must be less than dataset size ({adata.n_obs})"
            )

        # Statistical warning: n_neighbors > sqrt(N) may reduce clustering resolution
        # Based on:
        # - KNN rule of thumb: k ≈ sqrt(N) or sqrt(N)/2
        # - UMAP developers: "range 5 to 50, with a choice of 10 to 15 being sensible default"
        # - GitHub discussion (scverse/scanpy#223): small datasets (~1000 cells) use n_neighbors=5
        # Examples: 100 cells → k≈10, 50 cells → k≈7, 20 cells → k≈4-5
        recommended_max = int(np.sqrt(adata.n_obs))
        if n_neighbors > recommended_max and context:
            await context.warning(
                f"n_neighbors={n_neighbors} exceeds the sqrt(N) rule of thumb (√{adata.n_obs} ≈ {recommended_max}).\n"
                f"   • High neighbor counts may reduce clustering resolution for small datasets\n"
                f"   • Rule of thumb: k ≈ √N (based on KNN literature)\n"
                f"   • Recommended: n_neighbors={recommended_max} or lower for this dataset\n"
                f"   • General guidelines: 5-50 (UMAP), 10-15 default (Scanpy)"
            )

        if context:
            await context.info(
                f"Using n_neighbors: {n_neighbors} (Scanpy industry standard)"
            )

        if context:
            if resolvi_used:
                await context.info(
                    f"Using {n_neighbors} neighbors with ResolVI latent representation for graph construction..."
                )
            elif scvi_used:
                await context.info(
                    f"Using {n_neighbors} neighbors with scVI latent representation for graph construction..."
                )
            else:
                await context.info(
                    f"Using {n_neighbors} neighbors and {n_pcs} PCs for graph construction..."
                )

        try:
            if resolvi_used:
                # Use ResolVI latent representation for neighbors
                sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_resolvi")
            elif scvi_used:
                # Use scVI latent representation for neighbors
                sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_scvi")
            else:
                sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

            # 9. Run UMAP for visualization
            if context:
                await context.info("Running UMAP...")
            sc.tl.umap(adata)

            # 10. Run clustering
            if context:
                await context.info("Running Leiden clustering...")

            # Use explicit clustering resolution
            resolution = params.clustering_resolution
            if context:
                await context.info(f"Using clustering resolution: {resolution}")

            if context:
                await context.info(
                    f"Using Leiden clustering with resolution {resolution}..."
                )

            sc.tl.leiden(adata, resolution=resolution, key_added=params.cluster_key)

            # Count clusters
            n_clusters = len(adata.obs[params.cluster_key].unique())

            # Compute diffusion map for trajectory analysis
            if context:
                await context.info("Computing diffusion map for trajectory analysis...")
            sc.tl.diffmap(adata)
        except Exception as e:
            if context:
                await context.error(f"Clustering failed: {str(e)}")
            # Don't silently create fallback clustering - let the LLM handle the failure
            raise RuntimeError(
                f"Clustering failed: {str(e)}. "
                "Please check data quality or try different preprocessing parameters."
            )

        # 11. Add RNA velocity preprocessing if requested
        if getattr(params, "enable_rna_velocity", False):
            if context:
                await context.info("Adding RNA velocity preprocessing...")
            try:
                # Check for velocity data (spliced/unspliced layers)
                has_velocity_data = (
                    "spliced" in adata.layers and "unspliced" in adata.layers
                )

                if has_velocity_data:
                    # Use unified velocity computation from trajectory module
                    from ..models.data import RNAVelocityParameters
                    from .trajectory import compute_rna_velocity

                    # Use velocity_params or create with defaults - TRANSPARENT
                    velocity_params = params.velocity_params or RNAVelocityParameters()

                    if params.velocity_params is None and context:
                        default_params = RNAVelocityParameters()
                        await context.info(
                            f"velocity_params not specified. Using default parameters:\n"
                            f"   • Method: {default_params.method}\n"
                            f"   • Mode: {default_params.scvelo_mode}\n"
                            f"   • N top genes: {default_params.n_top_genes}\n"
                            f"   • N neighbors: {default_params.n_neighbors}\n"
                            f"   • Min shared counts: {default_params.min_shared_counts}\n"
                            f"   For custom settings, provide velocity_params parameter"
                        )

                    # Compute velocity with unified function
                    adata = compute_rna_velocity(adata, params=velocity_params)

                    if context:
                        await context.info(
                            f"RNA velocity preprocessing completed using {velocity_params.scvelo_mode} mode"
                        )
                else:
                    if context:
                        await context.warning(
                            "RNA velocity requested but no spliced/unspliced layers found"
                        )
            except Exception as e:
                if context:
                    await context.warning(f"RNA velocity preprocessing failed: {e}")

        # 12. Add trajectory analysis preprocessing if requested
        if getattr(params, "enable_trajectory_analysis", False):
            if context:
                await context.info("Adding trajectory analysis preprocessing...")
            try:
                # Compute diffusion map for pseudotime analysis
                sc.tl.diffmap(adata)

                # Optionally compute DPT if root cell is specified
                if (
                    hasattr(params, "dpt_root_cell")
                    and params.dpt_root_cell is not None
                ):
                    if params.dpt_root_cell in adata.obs_names:
                        adata.uns["iroot"] = np.where(
                            adata.obs_names == params.dpt_root_cell
                        )[0][0]
                        sc.tl.dpt(adata)
                        if context:
                            await context.info(
                                "Diffusion pseudotime computed with specified root cell"
                            )
                    else:
                        if context:
                            await context.warning(
                                f"Root cell {params.dpt_root_cell} not found"
                            )

                if context:
                    await context.info("Trajectory analysis preprocessing completed")
            except Exception as e:
                if context:
                    await context.warning(
                        f"Trajectory analysis preprocessing failed: {e}"
                    )

        # 13. Add spatial domain-specific preprocessing if requested
        if getattr(params, "enable_spatial_domains", False) and spg is not None:
            if context:
                await context.info("Adding spatial domain-specific preprocessing...")
            try:
                # Apply SpaGCN-specific gene filtering
                spg.prefilter_genes(adata, min_cells=3)
                spg.prefilter_specialgenes(adata)
                if context:
                    await context.info("SpaGCN gene filtering completed")
            except Exception as e:
                if context:
                    await context.warning(f"Spatial domain preprocessing failed: {e}")

        # 14. For spatial data, compute spatial neighbors if not already done
        if params.spatial_key and (
            params.spatial_key in adata.uns
            or any(params.spatial_key in key for key in adata.obsm.keys())
        ):
            if context:
                await context.info("Computing spatial neighbors...")
            try:
                # Check if spatial neighbors already exist
                if "spatial_connectivities" not in adata.obsp:
                    # For MERFISH data, we need to ensure the spatial coordinates are correctly formatted
                    if params.spatial_key in adata.obsm:
                        # Check if this is MERFISH data by looking at the shape and content
                        if adata.obsm[params.spatial_key].shape[1] == 2:
                            # This is likely MERFISH or other single-cell resolution spatial data
                            # Use delaunay method which works better for single-cell data
                            sq.gr.spatial_neighbors(
                                adata, coord_type="generic", delaunay=True
                            )
                        else:
                            # Use default parameters for other spatial data types
                            sq.gr.spatial_neighbors(adata)
                    else:
                        # Use default parameters if spatial key is in uns but not in obsm
                        sq.gr.spatial_neighbors(adata)
            except Exception as e:
                if context:
                    await context.warning(
                        f"Could not compute spatial neighbors: {str(e)}"
                    )
                    await context.info("Continuing without spatial neighbors...")

        # Store the processed AnnData object back in the data store
        data_store[data_id]["adata"] = adata

        # Return preprocessing result
        return PreprocessingResult(
            data_id=data_id,
            n_cells=adata.n_obs,
            n_genes=adata.n_vars,
            n_hvgs=(
                int(sum(adata.var.highly_variable))
                if "highly_variable" in adata.var
                else n_hvgs
            ),
            clusters=n_clusters,
            qc_metrics=qc_metrics,
        )

    except Exception as e:
        error_msg = f"Error in preprocessing: {str(e)}"
        tb = traceback.format_exc()
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {tb}")
        raise RuntimeError(f"{error_msg}\n{tb}")


def detect_spatial_key(adata: AnnData, user_key: Optional[str] = None) -> Optional[str]:
    """
    Intelligently detect spatial coordinate key in the AnnData object

    Priority:
    1. User-specified key (if exists)
    2. Common spatial key names in obsm
    3. Create from obs columns if x,y coordinates exist

    Args:
        adata: AnnData object
        user_key: User-specified spatial key name

    Returns:
        Spatial key name if found, None otherwise
    """
    # If user specified a key, validate it exists
    if user_key:
        if user_key in adata.obsm:
            coords = adata.obsm[user_key]
            if coords.shape[1] >= 2:
                logger.info(f"Using user-specified spatial key: '{user_key}'")
                return user_key
        else:
            logger.warning(f"User-specified spatial key '{user_key}' not found in obsm")

    # Auto-detect common spatial key names
    common_keys = [
        "X_spatial",  # ResolVI standard
        "spatial",  # Common default
        "coordinates",  # Some formats
        "spatial_coords",  # Variant
        "X_coordinates",  # Another variant
        "xy_coords",  # Additional variant
    ]

    for key in common_keys:
        if key in adata.obsm:
            coords = adata.obsm[key]
            # Validate it's proper spatial coordinates (at least 2D)
            if coords.shape[1] >= 2:
                logger.info(f"Auto-detected spatial coordinates in '{key}'")
                return key

    # Check if coordinates exist in obs columns and create obsm entry
    if "x" in adata.obs.columns and "y" in adata.obs.columns:
        adata.obsm["X_spatial"] = adata.obs[["x", "y"]].values.astype(np.float32)
        logger.info("Created X_spatial from x,y columns in obs")
        return "X_spatial"

    if "x_centroid" in adata.obs.columns and "y_centroid" in adata.obs.columns:
        adata.obsm["X_spatial"] = adata.obs[["x_centroid", "y_centroid"]].values.astype(
            np.float32
        )
        logger.info("Created X_spatial from x_centroid,y_centroid columns in obs")
        return "X_spatial"

    # Check for other coordinate column naming conventions
    if "X" in adata.obs.columns and "Y" in adata.obs.columns:
        adata.obsm["X_spatial"] = adata.obs[["X", "Y"]].values.astype(np.float32)
        logger.info("Created X_spatial from X,Y columns in obs")
        return "X_spatial"

    return None


async def preprocess_with_resolvi(
    adata,
    params: Optional["ResolVIPreprocessingParameters"] = None,
    context: Optional[Context] = None,
) -> Dict[str, Any]:
    """Preprocess spatial transcriptomics data using ResolVI

    ResolVI addresses noise and bias in single-cell resolved spatial
    transcriptomics data through molecular reassignment correction.

    Scientific Requirements:
    - High-resolution cellular-resolved spatial data (Xenium, MERFISH, etc.)
    - Real spatial coordinates (cannot work without spatial information)
    - Preferably raw count data (not log-transformed)

    Args:
        adata: Spatial transcriptomics AnnData object
        params: ResolVI parameters (uses defaults if None)
        context: MCP context for logging

    Returns:
        Dictionary containing ResolVI preprocessing results

    Raises:
        ImportError: If scvi-tools package is not available
        ValueError: If data doesn't meet ResolVI requirements
        RuntimeError: If RESOLVI computation fails
    """
    try:
        # Import the parameter model if needed
        from ..models.data import ResolVIPreprocessingParameters

        if scvi is None or RESOLVI is None:
            raise ImportError(
                "scvi-tools package is required for ResolVI preprocessing. "
                "Install with 'pip install scvi-tools'"
            )

        # Use default parameters if none provided
        if params is None:
            params = ResolVIPreprocessingParameters()

        if context:
            await context.info("Starting ResolVI spatial data preprocessing...")

        # 1. Validate data size
        if adata.n_obs < params.min_cells:
            raise ValueError(
                f"ResolVI requires at least {params.min_cells} cells for meaningful analysis. "
                f"Found only {adata.n_obs} cells. "
                "For small datasets, consider using standard denoising methods instead."
            )

        # 2. Detect and validate spatial coordinates (CRITICAL)
        spatial_key = detect_spatial_key(adata, params.spatial_key)

        if spatial_key is None:
            if params.require_spatial:
                # No spatial coordinates found - this is a critical error for ResolVI
                raise ValueError(
                    "ResolVI REQUIRES real spatial coordinates for molecular reassignment.\n"
                    "No spatial coordinates found in adata.obsm.\n"
                    "ResolVI is specifically designed for high-resolution spatial transcriptomics:\n"
                    "  - Xenium (10X Genomics)\n"
                    "  - MERFISH (Vizgen)\n"
                    "  - CosMx (NanoString)\n"
                    "  - Visium HD (10X Genomics)\n"
                    "\nFor non-spatial data, consider alternatives:\n"
                    "  - scVI for general denoising\n"
                    "  - MAGIC for imputation\n"
                    "  - DCA for count data denoising"
                )
            else:
                # require_spatial is False but no spatial coordinates found
                raise ValueError(
                    "No spatial coordinates found in dataset. "
                    "ResolVI requires spatial information. "
                    "Please check that your data has spatial coordinates in adata.obsm."
                )

        if context:
            await context.info(f"Using spatial coordinates from '{spatial_key}'")

        # MEMORY OPTIMIZATION: ResolVI natively supports sparse matrices
        # Verified with scvi-tools v1.3.0 - both CSR and dense formats work correctly
        #
        # Memory savings for typical spatial transcriptomics datasets:
        #   - 10K cells × 5K genes: ~181 MB saved (95% reduction)
        #   - 50K cells × 10K genes: ~1.8 GB saved (95% reduction)
        #   - Xenium/MERFISH (100K+ cells): ~5-10 GB saved
        #
        # ResolVI will handle both sparse and dense formats correctly during training.
        # We only need to ensure integer counts (conversion below preserves sparsity).
        import scipy.sparse as sp

        # RESOLVI expects count data (integers)
        # Check if data appears to be normalized (has decimals)
        # This check works for both sparse and dense matrices
        if sp.issparse(adata.X):
            # For sparse matrices, check only non-zero values
            has_decimals = np.any(adata.X.data != np.round(adata.X.data))
        else:
            # For dense matrices, check all values
            has_decimals = np.any(adata.X != adata.X.astype(int))

        if has_decimals:
            if context:
                await context.warning(
                    "RESOLVI expects integer count data, but found decimals. Converting to integer counts."
                )
            # Round to nearest integer and ensure non-negative
            # Sparse-aware conversion preserves memory efficiency
            if sp.issparse(adata.X):
                # For sparse: only modify non-zero data values
                adata.X.data = np.maximum(0, np.round(adata.X.data)).astype(np.int32)
                adata.X.eliminate_zeros()  # Remove any zeros created by rounding
            else:
                # For dense: standard numpy operations
                adata.X = np.maximum(0, np.round(adata.X)).astype(np.int32)
        else:
            # Ensure integer type (preserves sparsity if already sparse)
            adata.X = adata.X.astype(np.int32)

        if context:
            await context.info(
                f"Preprocessing {adata.n_obs} spots and {adata.n_vars} genes with RESOLVI"
            )

        # Setup RESOLVI - at this point spatial_key is guaranteed to be not None
        # and should exist in adata.obsm due to detect_spatial_key validation
        if spatial_key not in adata.obsm:
            raise RuntimeError(
                f"Internal error: Spatial key '{spatial_key}' not found in adata.obsm. "
                f"Available keys: {list(adata.obsm.keys())}. "
                "This should not happen - please report this issue."
            )

        # Ensure spatial coordinates are in the correct format for RESOLVI
        spatial_coords = adata.obsm[spatial_key]
        if spatial_coords.shape[1] != 2:
            raise ValueError(
                f"ResolVI requires 2D spatial coordinates, but got shape {spatial_coords.shape}. "
                "Expected (n_cells, 2) for x,y coordinates."
            )

        # Add spatial coordinates as X_spatial for RESOLVI (required by the model)
        adata.obsm["X_spatial"] = spatial_coords

        # RESOLVI setup without spatial_key parameter (not supported in current version)
        RESOLVI.setup_anndata(adata)

        # Create RESOLVI model with proper parameter types
        # Disable downsample_counts to avoid the LogNormal parameter issue
        resolvi_model = RESOLVI(
            adata,
            n_hidden=int(params.n_hidden),
            n_latent=int(params.n_latent),
            downsample_counts=False,  # Disable to avoid torch tensor type issues
        )

        if context:
            await context.info("Training RESOLVI model...")

        # Train model
        if params.use_gpu and torch.cuda.is_available():
            resolvi_model.train(max_epochs=params.n_epochs, accelerator="gpu")
        else:
            resolvi_model.train(max_epochs=params.n_epochs)

        if context:
            await context.info("RESOLVI training completed")

        # Get results
        if context:
            await context.info("Extracting RESOLVI denoised data...")

        # Get denoised expression
        denoised_expression = resolvi_model.get_normalized_expression()

        # Get latent representation (without distribution parameters)
        latent = resolvi_model.get_latent_representation(return_dist=False)

        # Store results in adata
        adata.layers["resolvi_denoised"] = denoised_expression
        adata.obsm["X_resolvi"] = latent  # Store as X_resolvi for downstream use

        # Calculate preprocessing statistics
        original_mean = np.mean(adata.X)
        denoised_mean = np.mean(denoised_expression)

        # Calculate noise reduction metrics
        # get_normalized_expression() returns pandas DataFrame with dense values
        # Convert original data to dense only for this statistical computation
        if sp.issparse(adata.X):
            orig_data = adata.X.toarray()
        else:
            orig_data = adata.X

        if hasattr(denoised_expression, "toarray"):
            denoised_data = denoised_expression.toarray()
        else:
            # DataFrame case - extract underlying numpy array
            denoised_data = (
                denoised_expression.values
                if hasattr(denoised_expression, "values")
                else denoised_expression
            )

        noise_reduction = np.mean(np.abs(orig_data - denoised_data))

        # Calculate summary statistics
        results = {
            "method": "RESOLVI",
            "n_latent_dims": params.n_latent,
            "n_epochs": params.n_epochs,
            "denoising_completed": True,
            "original_mean_expression": float(original_mean),
            "denoised_mean_expression": float(denoised_mean),
            "noise_reduction_metric": float(noise_reduction),
            "latent_shape": latent.shape,
            "denoised_shape": denoised_expression.shape,
            "training_completed": True,
            "device": (
                "GPU" if (params.use_gpu and torch.cuda.is_available()) else "CPU"
            ),
        }

        if context:
            await context.info("RESOLVI preprocessing completed successfully")
            await context.info(
                "Stored denoised data in adata.layers['resolvi_denoised']"
            )
            await context.info(
                "Stored latent representation in adata.obsm['X_resolvi']"
            )
            await context.info(f"Noise reduction metric: {noise_reduction:.4f}")

        return results

    except Exception as e:
        error_msg = f"RESOLVI preprocessing failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e
