"""
Preprocessing tools for spatial transcriptomics data.
"""

from typing import Any

import numpy as np
import scanpy as sc
import scipy.sparse

from ..models.analysis import PreprocessingResult
from ..models.data import PreprocessingParameters
from ..spatial_mcp_adapter import ToolContext
from ..utils.adata_utils import (
    check_is_integer_counts,
    ensure_unique_var_names_async,
    sample_expression_values,
    standardize_adata,
)
from ..utils.dependency_manager import require, validate_r_environment
from ..utils.exceptions import (
    ChatSpatialError,
    DataError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


def _select_hvgs_by_variance(adata, n_hvgs: int) -> None:
    """Select HVGs by finite per-gene variance when Scanpy binning is invalid."""
    if scipy.sparse.issparse(adata.X):
        means = np.asarray(adata.X.mean(axis=0)).ravel()
        sq_means = np.asarray(adata.X.power(2).mean(axis=0)).ravel()
        var = sq_means - np.square(means)
    else:
        var = np.var(np.asarray(adata.X), axis=0)
    var = np.nan_to_num(var, nan=-np.inf, posinf=-np.inf, neginf=-np.inf)
    n_select = min(max(int(n_hvgs), 1), adata.n_vars)
    top_idx = np.argpartition(var, -n_select)[-n_select:]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[top_idx] = True
    adata.var["highly_variable"] = mask


def _compute_safe_percent_top(n_genes: int) -> list[int] | None:
    """Compute valid percent_top values for scanpy QC metrics.

    scanpy's calculate_qc_metrics requires all percent_top values < n_genes,
    otherwise raises IndexError. This function adapts the standard defaults
    [50, 100, 200, 500] to work with any dataset size.
    """
    if n_genes <= 1:
        return None

    # Standard scanpy defaults, filtered to valid range
    result = [p for p in [50, 100, 200, 500] if p < n_genes]

    # For small datasets (< 50 genes), use proportional values instead
    if not result:
        result = [
            max(1, int(n_genes * f))
            for f in [0.1, 0.25, 0.5]
            if int(n_genes * f) < n_genes
        ]

    # Include n_genes - 1 as maximum coverage point
    result.append(n_genes - 1)

    return sorted(set(result)) or None


def _calculate_qc_metrics(adata) -> tuple[str | None, dict[str, int | float]]:
    """Annotate gene classes and return pre-filtering QC statistics."""
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))
    try:
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt", "ribo"],
            percent_top=_compute_safe_percent_top(adata.n_vars),
            inplace=True,
        )
    except Exception as exc:
        raise ProcessingError(
            f"QC metrics failed: {exc}. "
            f"Data: {adata.n_obs}×{adata.n_vars}, type: {type(adata.X).__name__}"
        ) from exc

    mito_pct_col = "pct_counts_mt" if "pct_counts_mt" in adata.obs else None
    metrics: dict[str, int | float] = {
        "n_cells_before_filtering": int(adata.n_obs),
        "n_genes_before_filtering": int(adata.n_vars),
        "median_genes_per_cell": float(np.median(adata.obs.n_genes_by_counts)),
        "median_umi_per_cell": float(np.median(adata.obs.total_counts)),
    }
    if mito_pct_col:
        metrics.update(
            {
                "median_mito_pct": float(np.median(adata.obs[mito_pct_col])),
                "max_mito_pct": float(np.max(adata.obs[mito_pct_col])),
                "n_mt_genes": int(adata.var["mt"].sum()),
            }
        )
    return mito_pct_col, metrics


async def _filter_and_subsample(
    adata,
    params: PreprocessingParameters,
    ctx: ToolContext,
    qc_metrics: dict[str, int | float | bool | str],
    mito_pct_col: str | None,
):
    """Apply filtering and spot subsampling without mutating the source object."""
    if params.filter_genes_min_cells is not None and params.filter_genes_min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=params.filter_genes_min_cells)
    if params.filter_cells_min_genes is not None and params.filter_cells_min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=params.filter_cells_min_genes)

    if params.filter_mito_pct is not None and mito_pct_col:
        high_mito_mask = adata.obs[mito_pct_col] > params.filter_mito_pct
        n_high_mito = int(high_mito_mask.sum())
        if n_high_mito:
            adata = adata[~high_mito_mask].copy()
            qc_metrics["n_spots_filtered_mito"] = n_high_mito
    elif params.filter_mito_pct is not None:
        await ctx.warning(
            "Mitochondrial filtering requested but no mito genes detected. "
            "This may indicate non-standard gene naming or imaging-based data."
        )

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise DataError(
            f"No data remaining after filtering: {adata.n_obs} cells, "
            f"{adata.n_vars} genes. Relax filtering parameters "
            "(filter_genes_min_cells, filter_cells_min_genes, filter_mito_pct)."
        )
    if params.subsample_spots is not None and params.subsample_spots < adata.n_obs:
        sc.pp.subsample(
            adata,
            n_obs=params.subsample_spots,
            random_state=params.subsample_random_seed,
        )
    return adata


async def _run_scrublet(
    adata,
    params: PreprocessingParameters,
    ctx: ToolContext,
    qc_metrics: dict[str, int | float | bool | str],
):
    """Run optional doublet detection and return the accepted workspace."""
    qc_metrics["scrublet_requested"] = params.use_scrublet
    qc_metrics["use_scrublet"] = False
    if not params.use_scrublet:
        return adata

    try:
        min_cells_for_scrublet = 100
        if adata.n_obs < min_cells_for_scrublet:
            await ctx.warning(
                f"Scrublet requires at least {min_cells_for_scrublet} cells, "
                f"but only {adata.n_obs} present. Skipping doublet detection."
            )
            qc_metrics["scrublet_skip_reason"] = "insufficient_cells"
            return adata

        candidate = adata.copy()
        sc.pp.scrublet(
            candidate,
            expected_doublet_rate=params.scrublet_expected_doublet_rate,
            threshold=params.scrublet_threshold,
            sim_doublet_ratio=params.scrublet_sim_doublet_ratio,
            n_prin_comps=min(params.scrublet_n_prin_comps, adata.n_vars - 1),
            batch_key=(
                params.batch_key if params.batch_key in candidate.obs.columns else None
            ),
        )
        n_doublets = int(candidate.obs["predicted_doublet"].sum())
        doublet_rate = n_doublets / candidate.n_obs
        qc_metrics.update(
            {
                "use_scrublet": True,
                "n_doublets_detected": n_doublets,
                "doublet_rate": float(doublet_rate),
                "scrublet_threshold": float(
                    params.scrublet_threshold
                    if params.scrublet_threshold is not None
                    else candidate.uns.get("scrublet", {}).get("threshold", 0.0)
                ),
                "median_doublet_score": float(
                    np.median(candidate.obs["doublet_score"])
                ),
            }
        )
        action = "kept in dataset"
        if params.scrublet_filter_doublets and n_doublets > 0:
            candidate = candidate[~candidate.obs["predicted_doublet"]].copy()
            qc_metrics["n_cells_after_doublet_filter"] = int(candidate.n_obs)
            action = "removed from dataset"
        await ctx.info(
            f"Scrublet: Detected {n_doublets} doublets "
            f"({doublet_rate:.1%}), {action}."
        )
        return candidate
    except Exception as exc:
        await ctx.warning(
            f"Scrublet doublet detection failed: {exc}. "
            "Continuing without doublet filtering."
        )
        qc_metrics["scrublet_error"] = str(exc)
        return adata


def _preserve_raw_counts(adata, params: PreprocessingParameters) -> None:
    """Freeze raw counts and record the preprocessing input contract."""
    import anndata as ad_module

    if adata.raw is None:
        adata.raw = ad_module.AnnData(
            X=adata.X.copy(),
            var=adata.var,
            obs=adata.obs.copy(),
            uns={},
        )
    if adata.raw is not None and adata.raw.shape == adata.shape:
        adata.layers["counts"] = adata.raw.X
    else:
        adata.layers["counts"] = adata.X.copy()

    adata.uns["preprocessing"] = {
        "normalization": params.normalization,
        "raw_preserved": True,
        "counts_layer": True,
        "n_genes_before_norm": adata.n_vars,
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


async def _select_and_subsample_genes(
    adata,
    params: PreprocessingParameters,
    ctx: ToolContext,
):
    """Select HVGs, apply exclusions, and honor explicit gene subsampling."""
    requested_gene_count = params.subsample_genes
    use_sct_hvgs = (
        params.normalization == "sct" and "highly_variable" in adata.var.columns
    )
    if use_sct_hvgs:
        current_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
        current_count = int(current_mask.sum())
        n_hvgs = (
            min(requested_gene_count, current_count)
            if requested_gene_count is not None
            else current_count
        )
        if 0 < n_hvgs < current_count:
            selected_idx = np.flatnonzero(current_mask)
            if "sct_residual_variance" in adata.var.columns:
                residual_var = adata.var["sct_residual_variance"].to_numpy()
                selected_scores = residual_var[selected_idx]
                top_local = np.argpartition(selected_scores, -n_hvgs)[-n_hvgs:]
                selected_idx = selected_idx[top_local]
            else:
                selected_idx = selected_idx[:n_hvgs]
            new_mask = np.zeros(adata.n_vars, dtype=bool)
            new_mask[selected_idx] = True
            adata.var["highly_variable"] = new_mask
    else:
        n_hvgs = (
            min(requested_gene_count, adata.n_vars - 1, params.n_hvgs)
            if requested_gene_count is not None
            else min(params.n_hvgs, adata.n_vars - 1)
        )

    if n_hvgs <= 0:
        n_hvgs = 1
        await ctx.warning(
            "Computed n_hvgs=0; forcing to 1 to avoid selecting all genes."
        )

    if not use_sct_hvgs:
        if adata.n_vars < 100:
            if requested_gene_count is not None and n_hvgs < adata.n_vars:
                _select_hvgs_by_variance(adata, n_hvgs)
            else:
                adata.var["highly_variable"] = True
        else:
            try:
                sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs)
            except KeyError as exc:
                if "nan" not in str(exc).lower():
                    raise ProcessingError(
                        f"HVG selection failed: {exc}. "
                        f"Data: {adata.n_obs}×{adata.n_vars}, requested: {n_hvgs} HVGs."
                    ) from exc
                _select_hvgs_by_variance(adata, n_hvgs)
                await ctx.warning(
                    "Scanpy HVG binning produced NaN bins; "
                    "selected HVGs by finite variance instead."
                )
            except Exception as exc:
                raise ProcessingError(
                    f"HVG selection failed: {exc}. "
                    f"Data: {adata.n_obs}×{adata.n_vars}, requested: {n_hvgs} HVGs."
                ) from exc

    if params.remove_mito_genes and "mt" in adata.var.columns:
        adata.var.loc[adata.var["mt"], "highly_variable"] = False
    if params.remove_ribo_genes and "ribo" in adata.var.columns:
        adata.var.loc[adata.var["ribo"], "highly_variable"] = False
    if "highly_variable" not in adata.var:
        raise ProcessingError(
            "HVG selection failed: no highly_variable annotations were produced."
        )

    selected_count = int(adata.var["highly_variable"].sum())
    if selected_count == 0:
        raise DataError(
            "HVG selection produced no usable genes. Check the normalization, "
            "HVG, and mitochondrial/ribosomal exclusion parameters."
        )
    if selected_count < adata.n_vars and adata.n_vars < 500:
        await ctx.warning(
            f"Using {selected_count} of {adata.n_vars} available genes in a small panel.\n"
            "   • Small panels are already curated and often benefit from retaining all genes\n"
            "   • Consider increasing n_hvgs or removing subsample_genes\n"
            f"   • Current dataset: {adata.n_obs} cells × {adata.n_vars} total genes"
        )
    elif adata.n_vars >= 500 and selected_count < 500:
        await ctx.warning(
            f"Using only {selected_count} HVGs is below the recommended minimum of 500 genes.\n"
            "   • Literature consensus: 500-5000 genes (typical: 1000-2000)\n"
            "   • Low gene counts may lead to unstable clustering results\n"
            "   • Recommended: Use n_hvgs=1000-2000 for most analyses\n"
            f"   • Current dataset: {adata.n_obs} cells × {adata.n_vars} total genes"
        )
    if requested_gene_count is not None and requested_gene_count < adata.n_vars:
        adata = adata[:, adata.var["highly_variable"]].copy()
    return adata


async def _scale_preprocessed_data(
    adata,
    params: PreprocessingParameters,
    ctx: ToolContext,
) -> tuple[Any, bool, str | None]:
    """Apply scaling transactionally so failure leaves expression unchanged."""
    if not params.scale:
        return adata, False, None
    try:
        candidate = adata.copy()
        sc.pp.scale(candidate, max_value=params.scale_max_value)
        values = (
            candidate.X.data
            if scipy.sparse.issparse(candidate.X)
            else np.asarray(candidate.X)
        )
        if not np.isfinite(values).all():
            raise ProcessingError("Scaling produced non-finite expression values")
        return candidate, True, None
    except Exception as exc:
        await ctx.warning(f"Scaling failed: {exc}. Continuing without scaling.")
        return adata, False, str(exc)


async def preprocess_data(
    data_id: str,
    ctx: ToolContext,
    params: PreprocessingParameters | None = None,
) -> PreprocessingResult:
    """Preprocess spatial transcriptomics data

    Args:
        data_id: Dataset ID
        ctx: Tool context for data access and logging
        params: Preprocessing parameters

    Returns:
        Preprocessing result summary
    """
    if params is None:
        params = PreprocessingParameters()

    # Work on a copy and commit via ctx.set_adata only after all steps succeed.
    source_adata = await ctx.get_adata(data_id)

    # Standardize data format at the entry point
    try:
        adata = standardize_adata(source_adata, copy=True)
    except Exception as e:
        await ctx.warning(
            f"Data standardization failed: {e}. Proceeding with original data."
        )
        adata = source_adata.copy()

    # Validate input data
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise DataError(
            f"Dataset {data_id} is empty: {adata.n_obs} cells, {adata.n_vars} genes"
        )

    # Handle duplicate gene names (must be done before gene-based operations)
    await ensure_unique_var_names_async(adata, ctx, "data")

    mito_pct_col, base_qc_metrics = _calculate_qc_metrics(adata)
    qc_metrics: dict[str, int | float | bool | str] = dict(base_qc_metrics)
    adata = await _filter_and_subsample(
        adata,
        params,
        ctx,
        qc_metrics,
        mito_pct_col,
    )
    adata = await _run_scrublet(adata, params, ctx, qc_metrics)
    _preserve_raw_counts(adata, params)

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
            raise DataError(error_msg)

        if params.normalize_target_sum is not None:
            sc.pp.normalize_total(adata, target_sum=params.normalize_target_sum)
        else:
            # Calculate median for adaptive normalization
            calculated_median = np.median(np.array(adata.X.sum(axis=1)).flatten())
            sc.pp.normalize_total(adata, target_sum=calculated_median)
        sc.pp.log1p(adata)
    elif params.normalization == "sct":
        # SCTransform v2 variance-stabilizing normalization via R's sctransform
        # Validate the count contract before loading the optional R runtime.
        X_sample = sample_expression_values(adata)
        if np.any((X_sample % 1) != 0):
            raise DataError(
                "SCTransform requires raw count data (integers). "
                "Use normalization='log' for normalized data."
            )

        try:
            r_env = validate_r_environment(
                ctx, required_packages=["sctransform", "Matrix"]
            )
        except DependencyError as e:
            full_error = (
                f"SCTransform requires R and the sctransform package.\n\n"
                f"ERROR: {e}\n\n"
                "INSTALLATION:\n"
                "  1. Install R (https://cran.r-project.org/)\n"
                "  2. In R: install.packages('sctransform')\n"
                "  3. pip install 'rpy2>=3.5.0'\n\n"
                "ALTERNATIVES:\n"
                "• Use normalization='pearson_residuals' (built-in, similar results)\n"
                "• Use normalization='log' (standard method)"
            )
            raise DependencyError(full_error) from e

        # Map method parameter to vst.flavor
        vst_flavor = "v2" if params.sct_method == "fix-slope" else "v1"

        try:
            # Note: counts layer is already created earlier in this preprocessing workflow.
            # It will be properly subsetted if SCT filters genes
            # Convert to sparse CSC matrix (genes × cells) for R's dgCMatrix
            counts_sparse = scipy.sparse.csc_matrix(adata.X.T)

            # Keep the complete R session under one lock and conversion context.
            with r_env.local_context(numpy=True) as r_context:
                ro = r_env.robjects
                r_context["sp_data"] = counts_sparse.data.astype(np.float64)
                r_context["sp_indices"] = counts_sparse.indices.astype(np.int32)
                r_context["sp_indptr"] = counts_sparse.indptr.astype(np.int32)
                r_context["n_genes"] = counts_sparse.shape[0]
                r_context["n_cells"] = counts_sparse.shape[1]
                r_context["gene_names"] = ro.StrVector(adata.var_names.tolist())
                r_context["cell_names"] = ro.StrVector(adata.obs_names.tolist())
                r_context["vst_flavor"] = vst_flavor
                r_context["n_cells_param"] = (
                    params.sct_n_cells if params.sct_n_cells else ro.NULL
                )

                # Reconstruct sparse matrix and run SCTransform in R.
                ro.r(
                    """
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

                pearson_residuals = np.array(ro.r("pearson_residuals"))
                residual_variance = np.array(ro.r("residual_variance"))
                kept_genes = list(ro.r("kept_genes"))

            # CRITICAL FIX: Subset adata to match genes returned by SCTransform
            # R's sctransform internally filters genes, so we need to subset
            n_genes_before_sct = adata.n_vars
            if len(kept_genes) != adata.n_vars:
                n_filtered = adata.n_vars - len(kept_genes)
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
                raise ProcessingError(error_msg)

            adata.var["sct_residual_variance"] = residual_variance

            # Select top N genes by residual variance
            n_hvg = min(params.sct_var_features_n, len(residual_variance))
            adata.var["highly_variable"] = False
            if n_hvg >= len(residual_variance):
                adata.var["highly_variable"] = True
            elif n_hvg > 0:
                top_hvg_indices = np.argpartition(residual_variance, -n_hvg)[-n_hvg:]
                adata.var.iloc[
                    top_hvg_indices, adata.var.columns.get_loc("highly_variable")
                ] = True

        except MemoryError as e:
            raise MemoryError(
                f"Memory error for SCTransform on {adata.n_obs}×{adata.n_vars} matrix. "
                f"Use normalization='log' or subsample data."
            ) from e
        except ChatSpatialError:
            raise
        except Exception as e:
            raise ProcessingError(f"SCTransform failed: {e}") from e
    elif params.normalization == "pearson_residuals":
        # Modern Pearson residuals normalization (recommended for UMI data)

        # Check if method is available
        if not hasattr(sc.experimental.pp, "normalize_pearson_residuals"):
            error_msg = (
                "Pearson residuals normalization not available (requires scanpy>=1.9.0).\n"
                "Options:\n"
                "• Install newer scanpy: pip install 'scanpy>=1.9.0'\n"
                "• Use log normalization instead: params.normalization='log'\n"
                "• Skip normalization if data is pre-processed: params.normalization='none'"
            )
            raise DependencyError(error_msg)

        # Check if data appears to be raw counts
        X_sample = sample_expression_values(adata)

        # Check for non-integer values (indicates normalized data)
        if np.any((X_sample % 1) != 0):
            raise DataError(
                "Pearson residuals requires raw count data (integers). "
                "Data contains non-integer values. "
                "Use params.normalization='none' if data is already normalized, "
                "or params.normalization='log' for standard normalization."
            )

        # Execute normalization
        try:
            # Apply Pearson residuals normalization (to all genes)
            # Note: High variable gene selection happens later in the pipeline
            sc.experimental.pp.normalize_pearson_residuals(adata)
        except MemoryError as e:
            raise MemoryError(
                f"Insufficient memory for Pearson residuals on {adata.n_obs}×{adata.n_vars} matrix. "
                "Try reducing n_hvgs or use 'log' normalization."
            ) from e
        except Exception as e:
            raise ProcessingError(
                f"Pearson residuals normalization failed: {e}. "
                "Consider using 'log' normalization instead."
            ) from e
    elif params.normalization == "none":
        # Explicitly skip normalization

        # CRITICAL: Raw integer counts MUST NOT go through standard HVG selection.
        # Use check_is_integer_counts (same validator as scVI) instead of
        # ad-hoc heuristic — low-depth platforms (MERFISH, Xenium) have
        # max values well below 100 but are still raw counts.
        is_int, _, _ = check_is_integer_counts(adata.X)
        if is_int:
            raise DataError(
                "STATISTICAL ERROR: Cannot perform HVG selection on raw counts "
                "with normalization='none'.\n\n"
                "Your data contains integer counts, indicating raw (unnormalized) "
                "data. HVG selection on raw counts is statistically invalid "
                "because count variance scales non-linearly with expression "
                "level.\n\n"
                "REQUIRED ACTIONS:\n"
                "Option 1 (Recommended): Use normalization='log' for standard "
                "log-normalization\n"
                "Option 2: Use normalization='pearson_residuals' for "
                "variance-stabilizing normalization\n"
                "Option 3: Pre-normalize your data externally, then reload"
            )
    elif params.normalization == "scvi":
        # scVI deep learning-based normalization
        # Uses variational autoencoder to learn latent representation
        # Validate counts layer (created earlier) for scVI's count-based model
        is_int, has_neg, _ = check_is_integer_counts(adata.layers["counts"])
        if has_neg:
            raise DataError(
                "scVI requires non-negative count data. "
                "The counts layer contains negative values, "
                "indicating already-normalized data."
            )
        if not is_int:
            import logging

            logging.getLogger(__name__).warning(
                "scVI expects integer count data but the counts layer "
                "contains non-integer values. This may indicate "
                "pre-normalized data; results may be suboptimal."
            )

        scvi = require("scvi", ctx, feature="scVI normalization")

        try:
            # Note: counts layer is already created earlier in this preprocessing workflow.
            # scVI requires this layer for proper count-based modeling

            # Setup AnnData for scVI using the pre-saved counts layer
            scvi.model.SCVI.setup_anndata(
                adata,
                layer="counts",
                batch_key=(
                    params.batch_key if params.batch_key in adata.obs.columns else None
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

            # Train the model with user-configurable parameters
            scvi_model.train(
                max_epochs=params.scvi_max_epochs,
                early_stopping=params.scvi_early_stopping,
                early_stopping_patience=params.scvi_early_stopping_patience,
                early_stopping_monitor="elbo_validation",
                train_size=params.scvi_train_size,
            )

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

        except Exception as e:
            raise ProcessingError(f"scVI normalization failed: {e}") from e
    else:
        # Catch unknown normalization methods
        valid_methods = ["log", "sct", "pearson_residuals", "none", "scvi"]
        raise ParameterError(
            f"Unknown normalization method: '{params.normalization}'. "
            f"Valid options are: {', '.join(valid_methods)}"
        )

    adata = await _select_and_subsample_genes(adata, params, ctx)
    adata, scale_applied, scale_error = await _scale_preprocessed_data(
        adata,
        params,
        ctx,
    )

    qc_metrics["scale_requested"] = params.scale
    qc_metrics["scale_applied"] = scale_applied
    if scale_error is not None:
        qc_metrics["scale_error"] = scale_error

    # Store only metadata describing work completed by preprocessing.
    # Embedding and clustering parameters belong to compute_embeddings.
    adata.uns["preprocessing"]["completed"] = True
    adata.uns["preprocessing"]["scale_requested"] = params.scale
    adata.uns["preprocessing"]["scale_applied"] = scale_applied
    if scale_error is not None:
        adata.uns["preprocessing"]["scale_error"] = scale_error

    # clusters=0 indicates that compute_embeddings has not run yet.
    result = PreprocessingResult(
        data_id=data_id,
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        n_hvgs=(
            int(sum(adata.var.highly_variable)) if "highly_variable" in adata.var else 0
        ),
        clusters=0,
        qc_metrics=qc_metrics,
    )
    await ctx.set_adata(data_id, adata)
    return result
