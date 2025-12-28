"""
Preprocessing tools for spatial transcriptomics data.
"""

import traceback

import numpy as np
import scanpy as sc
import scipy.sparse

from ..models.analysis import PreprocessingResult
from ..models.data import PreprocessingParameters
from ..spatial_mcp_adapter import ToolContext
from ..utils.adata_utils import (ensure_unique_var_names_with_ctx,
                                 sample_expression_values, standardize_adata)
from ..utils.compute import ensure_pca
from ..utils.dependency_manager import require, validate_r_package
from ..utils.mcp_utils import mcp_tool_error_handler


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
    ctx: ToolContext,
    params: PreprocessingParameters = PreprocessingParameters(),
) -> PreprocessingResult:
    """Preprocess spatial transcriptomics data

    Args:
        data_id: Dataset ID
        ctx: Tool context for data access and logging
        params: Preprocessing parameters

    Returns:
        Preprocessing result summary
    """
    try:
        await ctx.info(
            f"Preprocessing dataset {data_id} with {params.normalization} normalization"
        )

        # Get AnnData directly via ToolContext (no redundant dict wrapping)
        # Memory optimization: Preprocessing modifies dataset in-place (by design)
        # Original counts are preserved in adata.raw
        adata = await ctx.get_adata(data_id)

        # Standardize data format at the entry point
        # This eliminates all downstream special cases for data format handling
        await ctx.info("Standardizing data structure to ChatSpatial format...")
        try:
            adata = standardize_adata(
                adata, copy=False, strict=False, preserve_original=True
            )
            await ctx.info("Data structure standardized successfully")
        except Exception as e:
            await ctx.warning(
                f"Data standardization failed: {e}. Proceeding with original data."
            )
            # Continue with original data if standardization fails

        # Validate input data
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError(
                f"Dataset {data_id} is empty: {adata.n_obs} cells, {adata.n_vars} genes"
            )

        # Handle duplicate gene names (must be done before gene-based operations)
        await ensure_unique_var_names_with_ctx(adata, ctx, "data")

        # 1. Calculate QC metrics (including mitochondrial percentage)
        await ctx.info("Calculating QC metrics...")
        try:
            # Identify mitochondrial genes (MT-* for human, mt-* for mouse)
            # This is needed for both QC metrics and optional filtering
            adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
            n_mt_genes = adata.var["mt"].sum()
            if n_mt_genes > 0:
                await ctx.info(
                    f"Identified {n_mt_genes} mitochondrial genes (MT-*/mt-*)"
                )

            # Identify ribosomal genes (RPS*, RPL* for human, Rps*, Rpl* for mouse)
            adata.var["ribo"] = adata.var_names.str.startswith(
                ("RPS", "RPL", "Rps", "Rpl")
            )
            n_ribo_genes = adata.var["ribo"].sum()
            if n_ribo_genes > 0:
                await ctx.info(
                    f"Identified {n_ribo_genes} ribosomal genes (RPS*/RPL*/Rps*/Rpl*)"
                )

            # FIX: Adjust percent_top for small datasets
            #
            # Problem: sc.pp.calculate_qc_metrics() uses default percent_top=[50, 100, 200, 500]
            # to calculate "percentage of counts in top N genes". When n_genes < 500,
            # scanpy raises IndexError: "Positions outside range of features"
            # (see scanpy/preprocessing/_qc.py line 392: check_ns decorator)
            #
            # Solution: Dynamically adjust percent_top to only include values < n_genes
            n_genes = adata.n_vars
            default_percent_top = [50, 100, 200, 500]

            # Filter to only include values that are valid for this dataset
            safe_percent_top = [p for p in default_percent_top if p < n_genes]

            # For very small datasets (n_genes < 50), create proportional values
            if not safe_percent_top:
                safe_percent_top = []
                for fraction in [0.1, 0.25, 0.5]:
                    val = max(1, int(n_genes * fraction))
                    if val < n_genes and val not in safe_percent_top:
                        safe_percent_top.append(val)

            # Add the largest possible value (n_genes - 1) if reasonable
            if n_genes > 1 and (n_genes - 1) not in safe_percent_top:
                safe_percent_top.append(n_genes - 1)

            safe_percent_top = (
                sorted(set(safe_percent_top)) if safe_percent_top else None
            )

            if safe_percent_top != default_percent_top[: len(safe_percent_top)]:
                await ctx.info(
                    f"Small dataset ({n_genes} genes): using percent_top={safe_percent_top}"
                )

            # Calculate QC metrics including mitochondrial and ribosomal percentages
            sc.pp.calculate_qc_metrics(
                adata,
                qc_vars=["mt", "ribo"],
                percent_top=safe_percent_top,
                inplace=True,
            )
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
            await ctx.error(error_msg)
            raise RuntimeError(error_msg) from e

        # Store original QC metrics before filtering (including mito stats)
        mito_pct_col = "pct_counts_mt" if "pct_counts_mt" in adata.obs else None
        qc_metrics = {
            "n_cells_before_filtering": int(adata.n_obs),
            "n_genes_before_filtering": int(adata.n_vars),
            "median_genes_per_cell": float(np.median(adata.obs.n_genes_by_counts)),
            "median_umi_per_cell": float(np.median(adata.obs.total_counts)),
        }
        # Add mitochondrial stats if available
        if mito_pct_col:
            qc_metrics["median_mito_pct"] = float(np.median(adata.obs[mito_pct_col]))
            qc_metrics["max_mito_pct"] = float(np.max(adata.obs[mito_pct_col]))
            qc_metrics["n_mt_genes"] = int(adata.var["mt"].sum())

        # 2. Apply user-controlled data filtering and subsampling
        await ctx.info("Applying data filtering and subsampling...")

        # Apply gene filtering using LLM-controlled parameters
        min_cells = params.filter_genes_min_cells
        if min_cells is not None and min_cells > 0:
            await ctx.info(f"Filtering genes: min_cells={min_cells}")
            sc.pp.filter_genes(adata, min_cells=min_cells)

        # Apply cell filtering using LLM-controlled parameters
        min_genes = params.filter_cells_min_genes
        if min_genes is not None and min_genes > 0:
            await ctx.info(f"Filtering cells: min_genes={min_genes}")
            sc.pp.filter_cells(adata, min_genes=min_genes)

        # Apply mitochondrial percentage filtering (BEST PRACTICE for spatial data)
        # High mito% indicates damaged cells that have lost cytoplasmic mRNA
        if params.filter_mito_pct is not None and mito_pct_col:
            n_before = adata.n_obs
            high_mito_mask = adata.obs[mito_pct_col] > params.filter_mito_pct
            n_high_mito = high_mito_mask.sum()

            if n_high_mito > 0:
                adata = adata[~high_mito_mask].copy()
                await ctx.info(
                    f"Filtered {n_high_mito} spots with mito% > {params.filter_mito_pct}% "
                    f"({n_before} → {adata.n_obs} spots)"
                )
                # Update qc_metrics with mito filtering info
                qc_metrics["n_spots_filtered_mito"] = int(n_high_mito)
            else:
                await ctx.info(
                    f"No spots exceeded mito% threshold of {params.filter_mito_pct}%"
                )
        elif params.filter_mito_pct is not None and not mito_pct_col:
            await ctx.warning(
                "Mitochondrial filtering requested but no mito genes detected. "
                "This may indicate non-standard gene naming or imaging-based data."
            )

        # Apply spot subsampling if requested
        if params.subsample_spots is not None and params.subsample_spots < adata.n_obs:
            await ctx.info(
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
        await ctx.info("Saving raw data for downstream analysis...")

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

        # Store counts layer for scVI-tools compatibility (Cell2location, scANVI, DestVI)
        # Note: This layer follows adata through HVG subsetting, complementing adata.raw
        # - adata.raw: Full gene set (for cell communication needing complete L-R coverage)
        # - adata.layers["counts"]: HVG subset after filtering (for scVI-tools alignment)
        adata.layers["counts"] = adata.X.copy()

        # Store preprocessing metadata following scanpy/anndata conventions
        # This metadata enables downstream tools to reuse gene annotations
        adata.uns["preprocessing"] = {
            "normalization": params.normalization,
            "raw_preserved": True,
            "counts_layer": True,
            "n_genes_before_norm": adata.n_vars,
            # Gene type annotations - downstream tools should reuse these
            "gene_annotations": {
                "mt_column": "mt" if "mt" in adata.var.columns else None,
                "ribo_column": "ribo" if "ribo" in adata.var.columns else None,
                "n_mt_genes": (
                    int(adata.var["mt"].sum()) if "mt" in adata.var.columns else 0
                ),
                "n_ribo_genes": (
                    int(adata.var["ribo"].sum()) if "ribo" in adata.var.columns else 0
                ),
            },
        }

        # Update QC metrics after filtering
        qc_metrics.update(
            {
                "n_cells_after_filtering": int(adata.n_obs),
                "n_genes_after_filtering": int(adata.n_vars),
            }
        )

        # 3. Normalize data
        # Log normalization configuration (developer log)
        norm_config = {
            "Method": params.normalization,
            "Target sum": (
                f"{params.normalize_target_sum:.0f}"
                if params.normalize_target_sum is not None
                else "ADAPTIVE (using median counts)"
            ),
        }
        if params.scale:
            norm_config["Scale clipping"] = (
                f"±{params.scale_max_value} SD"
                if params.scale_max_value is not None
                else "NONE (preserving all outliers)"
            )
        ctx.log_config("Normalization Configuration", norm_config)

        await ctx.info(f"Normalizing data using {params.normalization} method...")

        if params.normalization == "log":
            # Standard log normalization
            # Check if data appears to be already normalized
            X_sample = sample_expression_values(adata)

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
                await ctx.error(error_msg)
                raise ValueError(error_msg)

            if params.normalize_target_sum is not None:
                sc.pp.normalize_total(adata, target_sum=params.normalize_target_sum)
            else:
                # Calculate median and inform user transparently
                calculated_median = np.median(np.array(adata.X.sum(axis=1)).flatten())
                await ctx.info(
                    f"normalize_target_sum not specified. Using adaptive normalization:\n"
                    f"   • Calculated median counts: {calculated_median:.0f}\n"
                    f"   • This value was automatically determined from your data\n"
                    f"   • For reproducible results, use: normalize_target_sum={calculated_median:.0f}"
                )
                sc.pp.normalize_total(adata, target_sum=calculated_median)
            sc.pp.log1p(adata)
        elif params.normalization == "sct":
            # SCTransform v2 variance-stabilizing normalization via R's sctransform
            # Check R sctransform availability using centralized dependency manager
            try:
                validate_r_package("sctransform", ctx)
                validate_r_package("Matrix", ctx)
            except ImportError as e:
                full_error = (
                    f"SCTransform requires R and the sctransform package.\n\n"
                    f"ERROR: {str(e)}\n\n"
                    "INSTALLATION:\n"
                    "  1. Install R (https://cran.r-project.org/)\n"
                    "  2. In R: install.packages('sctransform')\n"
                    "  3. pip install 'rpy2>=3.5.0'\n\n"
                    "ALTERNATIVES:\n"
                    "• Use normalization='pearson_residuals' (built-in, similar results)\n"
                    "• Use normalization='log' (standard method)"
                )
                await ctx.error(full_error)
                raise ImportError(full_error) from e

            # Check if data appears to be raw counts (required for SCTransform)
            X_sample = sample_expression_values(adata)

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
                await ctx.error(error_msg)
                raise ValueError(error_msg)

            # Map method parameter to vst.flavor
            vst_flavor = "v2" if params.sct_method == "fix-slope" else "v1"
            await ctx.info(
                f"Applying SCTransform (vst.flavor={vst_flavor}, "
                f"var_features={params.sct_var_features_n})"
            )

            try:
                # Import rpy2 modules
                import rpy2.robjects as ro
                from rpy2.robjects import numpy2ri
                from rpy2.robjects.conversion import localconverter

                # Note: counts layer already saved in unified preprocessing step (line 338)
                # It will be properly subsetted if SCT filters genes
                # Convert to sparse CSC matrix (genes × cells) for R's dgCMatrix
                if scipy.sparse.issparse(adata.X):
                    counts_sparse = scipy.sparse.csc_matrix(adata.X.T)
                else:
                    counts_sparse = scipy.sparse.csc_matrix(adata.X.T)

                # Transfer sparse matrix components to R
                with localconverter(ro.default_converter + numpy2ri.converter):
                    ro.globalenv["sp_data"] = counts_sparse.data.astype(np.float64)
                    ro.globalenv["sp_indices"] = counts_sparse.indices.astype(np.int32)
                    ro.globalenv["sp_indptr"] = counts_sparse.indptr.astype(np.int32)
                    ro.globalenv["n_genes"] = counts_sparse.shape[0]
                    ro.globalenv["n_cells"] = counts_sparse.shape[1]
                    ro.globalenv["gene_names"] = ro.StrVector(adata.var_names.tolist())
                    ro.globalenv["cell_names"] = ro.StrVector(adata.obs_names.tolist())
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
                    # Extract gene names that survived SCTransform filtering
                    kept_genes <- rownames(vst_result$y)
                """
                )

                # Extract results from R
                with localconverter(ro.default_converter + numpy2ri.converter):
                    pearson_residuals = np.array(ro.r("pearson_residuals"))
                    residual_variance = np.array(ro.r("residual_variance"))
                    kept_genes = list(ro.r("kept_genes"))

                # CRITICAL FIX: Subset adata to match genes returned by SCTransform
                # R's sctransform internally filters genes, so we need to subset
                n_genes_before_sct = adata.n_vars
                if len(kept_genes) != adata.n_vars:
                    n_filtered = adata.n_vars - len(kept_genes)
                    await ctx.info(
                        f"SCTransform filtered {n_filtered} additional genes internally. "
                        f"Keeping {len(kept_genes)} genes for analysis."
                    )
                    # Subset adata to keep only genes returned by SCTransform
                    adata = adata[:, kept_genes].copy()
                else:
                    n_filtered = 0

                # Transpose back to cells × genes for AnnData format
                adata.X = pearson_residuals.T

                # Store SCTransform metadata
                adata.uns["sctransform"] = {
                    "method": params.sct_method,
                    "vst_flavor": vst_flavor,
                    "var_features_n": params.sct_var_features_n,
                    "exclude_poisson": params.sct_exclude_poisson,
                    "n_cells": params.sct_n_cells,
                    "n_genes_before": n_genes_before_sct,
                    "n_genes_after": len(kept_genes),
                    "n_genes_filtered_by_sct": n_filtered,
                }

                # Mark highly variable genes based on residual variance
                # Now adata has been subset, so residual_variance should match adata.n_vars
                if len(residual_variance) != adata.n_vars:
                    error_msg = (
                        f"Dimension mismatch after SCTransform: "
                        f"residual_variance has {len(residual_variance)} values "
                        f"but adata has {adata.n_vars} genes"
                    )
                    raise ValueError(error_msg)

                adata.var["sct_residual_variance"] = residual_variance

                # Select top N genes by residual variance
                n_hvg = min(params.sct_var_features_n, len(residual_variance))
                top_hvg_indices = np.argsort(residual_variance)[-n_hvg:]
                adata.var["highly_variable"] = False
                adata.var.iloc[
                    top_hvg_indices, adata.var.columns.get_loc("highly_variable")
                ] = True

                await ctx.info(
                    f"[OK]SCTransform: {n_hvg} highly variable genes identified"
                )

            except MemoryError:
                error_msg = (
                    f"Insufficient memory for SCTransform on "
                    f"{adata.n_obs}×{adata.n_vars} matrix.\n"
                    "SOLUTIONS:\n"
                    "• Reduce dataset size with subsample_spots parameter\n"
                    "• Use normalization='pearson_residuals' (more memory efficient)\n"
                    "• Use normalization='log' (minimal memory usage)"
                )
                await ctx.error(error_msg)
                raise MemoryError(error_msg)
            except Exception as e:
                error_msg = f"SCTransform failed: {str(e)}"
                await ctx.error(error_msg)
                await ctx.info(f"Error details: {traceback.format_exc()}")
                raise RuntimeError(error_msg) from e
        elif params.normalization == "pearson_residuals":
            # Modern Pearson residuals normalization (recommended for UMI data)
            await ctx.info("Using Pearson residuals normalization...")

            # Check if method is available
            if not hasattr(sc.experimental.pp, "normalize_pearson_residuals"):
                error_msg = (
                    "Pearson residuals normalization not available (requires scanpy>=1.9.0).\n"
                    "Options:\n"
                    "• Install newer scanpy: pip install 'scanpy>=1.9.0'\n"
                    "• Use log normalization instead: params.normalization='log'\n"
                    "• Skip normalization if data is pre-processed: params.normalization='none'"
                )
                await ctx.error(error_msg)
                raise ValueError(error_msg)

            # Check if data appears to be raw counts
            X_sample = sample_expression_values(adata)

            # Check for non-integer values (indicates normalized data)
            if np.any((X_sample % 1) != 0):
                error_msg = (
                    "Pearson residuals requires raw count data (integers). "
                    "Data contains non-integer values. "
                    "Use params.normalization='none' if data is already normalized, "
                    "or params.normalization='log' for standard normalization."
                )
                await ctx.error(error_msg)
                raise ValueError(error_msg)

            # Execute normalization
            try:
                # Apply Pearson residuals normalization (to all genes)
                # Note: High variable gene selection happens later in the pipeline
                sc.experimental.pp.normalize_pearson_residuals(adata)
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
            await ctx.info("Skipping normalization (data assumed to be pre-normalized)")

            # CRITICAL: Check if data appears to be raw counts
            # HVG selection requires normalized data for statistical validity
            X_sample = sample_expression_values(adata)

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
                await ctx.error(error_msg)
                raise ValueError(error_msg)
        elif params.normalization == "scvi":
            # scVI deep learning-based normalization
            # Uses variational autoencoder to learn latent representation
            require("scvi", feature="scVI normalization")
            import scvi

            # Check if data appears to be raw counts (required for scVI)
            X_sample = sample_expression_values(adata)

            # Check for negative values (indicates already normalized data)
            if np.any(X_sample < 0):
                error_msg = (
                    "scVI requires non-negative count data.\n"
                    "Data contains negative values, suggesting log-normalization.\n\n"
                    "SOLUTIONS:\n"
                    "• Load raw count data instead of normalized data\n"
                    "• Use normalization='none' if data is pre-normalized"
                )
                await ctx.error(error_msg)
                raise ValueError(error_msg)

            await ctx.info(
                f"Applying scVI normalization "
                f"(n_latent={params.scvi_n_latent}, "
                f"n_hidden={params.scvi_n_hidden}, "
                f"gene_likelihood={params.scvi_gene_likelihood})"
            )

            try:
                # Note: counts layer already saved in unified preprocessing step (line 338)
                # scVI requires this layer for proper count-based modeling

                # Setup AnnData for scVI using the pre-saved counts layer
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

                await ctx.info(
                    f"Training scVI model (max_epochs={params.scvi_max_epochs}, "
                    f"early_stopping={params.scvi_early_stopping}, "
                    f"train_size={params.scvi_train_size})..."
                )

                # Train the model with user-configurable parameters
                scvi_model.train(
                    max_epochs=params.scvi_max_epochs,
                    early_stopping=params.scvi_early_stopping,
                    early_stopping_patience=params.scvi_early_stopping_patience,
                    early_stopping_monitor="elbo_validation",
                    train_size=params.scvi_train_size,
                )

                await ctx.info("scVI training completed")

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

                await ctx.info(
                    f"[OK]scVI: Latent representation stored in X_scvi "
                    f"(shape: {adata.obsm['X_scvi'].shape})"
                )
                await ctx.info("[OK]scVI: Normalized expression stored in adata.X")

            except Exception as e:
                error_msg = f"scVI normalization failed: {str(e)}"
                await ctx.error(error_msg)
                await ctx.info(f"Error details: {traceback.format_exc()}")
                raise RuntimeError(error_msg) from e
        else:
            # Catch unknown normalization methods
            valid_methods = ["log", "sct", "pearson_residuals", "none", "scvi"]
            raise ValueError(
                f"Unknown normalization method: '{params.normalization}'. "
                f"Valid options are: {', '.join(valid_methods)}"
            )

        # 4. Find highly variable genes and apply gene subsampling
        await ctx.info("Finding highly variable genes...")

        # Determine number of HVGs to select
        if gene_subsample_requested:
            # User wants to subsample genes
            n_hvgs = min(params.subsample_genes, adata.n_vars - 1, params.n_hvgs)
            await ctx.info(
                f"User requested {params.subsample_genes} genes, selecting {n_hvgs} highly variable genes"
            )
        else:
            # Use standard HVG selection
            n_hvgs = min(params.n_hvgs, adata.n_vars - 1)
            await ctx.info(f"Using {n_hvgs} highly variable genes...")

        # Statistical warning: Very low HVG count may lead to unstable clustering
        # Based on literature consensus: 500-5000 genes recommended, 1000-2000 typical
        # References:
        # - Bioconductor OSCA: "any value from 500 to 5000 is reasonable"
        # - Single-cell best practices: typical range 1000-2000
        if n_hvgs < 500:
            await ctx.warning(
                f"Using only {n_hvgs} HVGs is below the recommended minimum of 500 genes.\n"
                f"   • Literature consensus: 500-5000 genes (typical: 1000-2000)\n"
                f"   • Low gene counts may lead to unstable clustering results\n"
                f"   • Recommended: Use n_hvgs=1000-2000 for most analyses\n"
                f"   • Current dataset: {adata.n_obs} cells × {adata.n_vars} total genes"
            )

        # Check if we should use all genes (for very small gene sets)
        if _should_use_all_genes_for_hvg(adata):
            await ctx.info(
                f"Small gene set detected ({adata.n_vars} genes), using all genes for analysis"
            )
            adata.var["highly_variable"] = True
        else:
            # Attempt HVG selection - no fallback for failures
            try:
                sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs)
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
                await ctx.error(error_msg)
                raise RuntimeError(error_msg) from e

        # Exclude mitochondrial genes from HVG selection (BEST PRACTICE)
        # Mito genes can dominate HVG due to high expression and technical variation
        if params.remove_mito_genes and "mt" in adata.var.columns:
            n_mito_hvg = (adata.var["highly_variable"] & adata.var["mt"]).sum()
            if n_mito_hvg > 0:
                adata.var.loc[adata.var["mt"], "highly_variable"] = False
                await ctx.info(
                    f"Excluded {n_mito_hvg} mitochondrial genes from HVG selection "
                    f"(genes retained in adata.raw for downstream analysis)"
                )

        # Exclude ribosomal genes from HVG selection (optional)
        if params.remove_ribo_genes and "ribo" in adata.var.columns:
            n_ribo_hvg = (adata.var["highly_variable"] & adata.var["ribo"]).sum()
            if n_ribo_hvg > 0:
                adata.var.loc[adata.var["ribo"], "highly_variable"] = False
                await ctx.info(
                    f"Excluded {n_ribo_hvg} ribosomal genes from HVG selection "
                    f"(genes retained in adata.raw for downstream analysis)"
                )

        # Log final HVG count after exclusions
        final_hvg_count = adata.var["highly_variable"].sum()
        await ctx.info(f"Final HVG count after exclusions: {final_hvg_count}")

        # Apply gene subsampling if requested
        if gene_subsample_requested and params.subsample_genes < adata.n_vars:
            # Ensure HVG selection was successful
            if "highly_variable" not in adata.var:
                error_msg = (
                    "Gene subsampling requested but no highly variable genes were identified. "
                    "This indicates a failure in the HVG selection step."
                )
                await ctx.error(error_msg)
                raise RuntimeError(error_msg)

            if not adata.var["highly_variable"].any():
                error_msg = (
                    "Gene subsampling requested but no genes were marked as highly variable. "
                    "Check HVG selection parameters or data quality."
                )
                await ctx.error(error_msg)
                raise RuntimeError(error_msg)

            # Use properly identified HVGs
            adata = adata[:, adata.var["highly_variable"]].copy()
            await ctx.info(f"Subsampled to {adata.n_vars} highly variable genes")

        # 5. Batch effect correction (if applicable)
        if (
            params.batch_key in adata.obs
            and len(adata.obs[params.batch_key].unique()) > 1
        ):
            await ctx.info(
                "Detected batch information. Applying batch effect correction with Harmony..."
            )
            try:
                # Use Harmony for batch correction (modern standard, works on PCA space)
                # Harmony is more robust than ComBat for single-cell/spatial data
                # Use centralized dependency manager for consistent error handling
                require(
                    "harmonypy"
                )  # Raises ImportError with install instructions if missing
                import scanpy.external as sce

                # Harmony requires PCA - use lazy computation
                ensure_pca(adata, n_comps=min(50, adata.n_vars - 1))

                sce.pp.harmony_integrate(adata, key=params.batch_key)
                await ctx.info("Batch effect correction completed using Harmony")
            except Exception as e:
                # Harmony failed - raise error, don't silently continue
                raise RuntimeError(
                    f"Batch effect correction with Harmony failed: {e}\n\n"
                    "Batch effects can severely impact downstream analyses. "
                    "Continuing without correction would produce unreliable results.\n\n"
                    "POSSIBLE CAUSES:\n"
                    "1. Insufficient cells per batch (need at least 30-50 per batch)\n"
                    "2. Too many batches relative to cells\n"
                    "3. PCA dimensionality issues\n\n"
                    "SOLUTIONS:\n"
                    "1. Check batch sizes: adata.obs['batch'].value_counts()\n"
                    "2. Merge small batches before preprocessing\n"
                    "3. Try alternative integration methods (scVI, BBKNN)"
                ) from e

        # 6. Scale data (if requested)
        if params.scale:
            await ctx.info("Scaling data...")
            try:
                # Trust scanpy's internal zero-variance handling and sparse matrix optimization
                sc.pp.scale(adata, max_value=params.scale_max_value)

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
                await ctx.warning(f"Scaling failed: {e}. Continuing without scaling.")

        # Store preprocessing metadata for downstream tools
        # PCA, UMAP, clustering, and spatial neighbors are computed lazily
        # by analysis tools using ensure_* functions from utils.compute
        adata.uns["preprocessing"]["completed"] = True
        adata.uns["preprocessing"]["n_pcs"] = params.n_pcs
        adata.uns["preprocessing"]["n_neighbors"] = params.n_neighbors
        adata.uns["preprocessing"][
            "clustering_resolution"
        ] = params.clustering_resolution

        await ctx.info(
            "Core preprocessing complete. "
            "PCA, UMAP, and clustering will be computed on-demand by analysis tools."
        )

        # Store the processed AnnData object back via ToolContext
        await ctx.set_adata(data_id, adata)

        # Return preprocessing result
        # Note: clusters=0 indicates clustering not yet performed
        # Analysis tools will compute clustering lazily when needed
        return PreprocessingResult(
            data_id=data_id,
            n_cells=adata.n_obs,
            n_genes=adata.n_vars,
            n_hvgs=(
                int(sum(adata.var.highly_variable))
                if "highly_variable" in adata.var
                else 0
            ),
            clusters=0,  # Clustering computed lazily by analysis tools
            qc_metrics=qc_metrics,
        )

    except Exception as e:
        error_msg = f"Error in preprocessing: {str(e)}"
        tb = traceback.format_exc()
        await ctx.warning(error_msg)
        await ctx.info(f"Error details: {tb}")
        raise RuntimeError(f"{error_msg}\n{tb}")
