"""
RCTD (Robust Cell Type Decomposition) deconvolution method.

RCTD is an R-based deconvolution method that performs robust
decomposition of cell type mixtures via the spacexr package.
"""

import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ...utils.dependency_manager import validate_r_environment
from ...utils.exceptions import (
    ChatSpatialError,
    DataError,
    ParameterError,
    ProcessingError,
)
from .base import PreparedDeconvolutionData, create_deconvolution_stats

_VALID_MODES = frozenset({"full", "doublet", "multi"})
_RCTD_SUBPROCESS_CODE = r"""
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
    stop("Expected input path, output path, and RCTD mode")
}
if (length(args) > 3L) {
    .libPaths(unique(c(args[4:length(args)], .libPaths())))
}
suppressPackageStartupMessages(library(spacexr))

input_bundle <- readRDS(args[[1L]])
if (!is.null(input_bundle$random_seed)) {
    assign(".Random.seed", input_bundle$random_seed, envir = .GlobalEnv)
}

myRCTD <- run.RCTD(input_bundle$rctd, doublet_mode = args[[3L]])
random_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
saveRDS(list(rctd = myRCTD, random_seed = random_seed), args[[2L]])
"""


def _run_rctd_subprocess(
    robjects: Any,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    mode: str,
) -> None:
    """Run noisy RCTD workers without redirecting process-global file descriptors."""
    rscript = Path(
        str(
            robjects.r(
                'file.path(R.home("bin"), '
                'if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")'
            )[0]
        )
    )
    if not rscript.is_file():
        raise ProcessingError(f"Rscript executable not found: {rscript}")

    library_paths = [str(path) for path in robjects.r(".libPaths()")]
    command = [
        str(rscript),
        "--vanilla",
        "-e",
        _RCTD_SUBPROCESS_CODE,
        str(input_path),
        str(output_path),
        mode,
        *library_paths,
    ]

    try:
        with log_path.open("w", encoding="utf-8") as log_file:
            completed = subprocess.run(
                command,
                stdin=subprocess.DEVNULL,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                check=False,
            )
    except OSError as exc:
        raise ProcessingError(f"Could not start the RCTD R subprocess: {exc}") from exc

    if completed.returncode != 0:
        try:
            details = log_path.read_text(encoding="utf-8", errors="replace")[-4000:]
        except OSError:
            details = ""
        message = f"RCTD R subprocess exited with status {completed.returncode}"
        if details := details.strip():
            message = f"{message}:\n{details}"
        raise ProcessingError(message)

    if not output_path.is_file():
        raise ProcessingError("RCTD R subprocess completed without producing a result")


def deconvolve(
    data: PreparedDeconvolutionData,
    mode: str = "full",
    max_cores: int = 4,
    confidence_threshold: float = 10.0,
    doublet_threshold: float = 25.0,
    max_multi_types: int = 4,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Deconvolve spatial data using RCTD from spacexr R package.

    Args:
        data: Prepared deconvolution data (immutable, includes spatial coordinates)
        mode: RCTD mode - 'full', 'doublet', or 'multi'
        max_cores: Maximum CPU cores
        confidence_threshold: Confidence threshold
        doublet_threshold: Doublet detection threshold
        max_multi_types: Max cell types per spot in multi mode

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    ctx = data.ctx

    if mode not in _VALID_MODES:
        choices = ", ".join(sorted(_VALID_MODES))
        raise ParameterError(f"Invalid RCTD mode '{mode}'. Expected one of: {choices}.")

    # Validate mode-specific parameters
    if mode == "multi" and max_multi_types >= data.n_cell_types:
        raise ParameterError(
            f"MAX_MULTI_TYPES ({max_multi_types}) must be less than "
            f"total cell types ({data.n_cell_types})."
        )

    try:
        # Data already copied in prepare_deconvolution
        spatial_data = data.spatial
        reference_data = data.reference

        # Get spatial coordinates from prepared data
        if data.spatial_coords is None:
            raise DataError(
                "RCTD requires real spatial coordinates for spatially-informed "
                "deconvolution. No spatial coordinates found. "
                "Use a non-spatial method (e.g., NNLS, DestVI) instead."
            )
        coords = pd.DataFrame(
            data.spatial_coords[:, :2],
            index=spatial_data.obs_names,
            columns=["x", "y"],
        )

        # Prepare cell type information
        cell_types = reference_data.obs[data.cell_type_key].astype(str)

        # RCTD requires minimum 25 cells per cell type
        MIN_CELLS_PER_TYPE = 25
        cell_type_counts = cell_types.value_counts()
        rare_types = cell_type_counts[
            cell_type_counts < MIN_CELLS_PER_TYPE
        ].index.tolist()

        if rare_types:
            warnings.warn(
                f"RCTD requires ≥{MIN_CELLS_PER_TYPE} cells per cell type. "
                f"Filtering {len(rare_types)} rare types: {rare_types}",
                UserWarning,
                stacklevel=2,
            )
            keep_mask = ~cell_types.isin(rare_types)
            reference_data = reference_data[keep_mask].copy()
            cell_types = cell_types[keep_mask]

            remaining_types = cell_types.unique()
            if len(remaining_types) < 2:
                raise DataError(
                    f"After filtering rare cell types, only {len(remaining_types)} "
                    f"cell type(s) remain. RCTD requires at least 2 cell types."
                )

        cell_types_series = pd.Series(
            cell_types.values, index=reference_data.obs_names, name="cell_type"
        )

        # Calculate nUMI
        spatial_numi = pd.Series(
            np.asarray(spatial_data.X.sum(axis=1)).ravel(),
            index=spatial_data.obs_names,
            name="nUMI",
        )
        reference_numi = pd.Series(
            np.asarray(reference_data.X.sum(axis=1)).ravel(),
            index=reference_data.obs_names,
            name="nUMI",
        )

        r_env = validate_r_environment(
            ctx,
            required_packages=["spacexr"],
            require_anndata2ri=True,
            package_install_commands={
                "spacexr": (
                    "devtools::install_github('dmcable/spacexr', "
                    "build_vignettes = FALSE)"
                )
            },
        )
        ro = r_env.robjects

        # spacexr starts PSOCK workers with outfile="", so their output bypasses
        # rpy2 callbacks. A request-private R subprocess contains those file
        # descriptors without mutating stdout or stderr for concurrent MCP work.
        with tempfile.TemporaryDirectory(prefix="chatspatial-rctd-") as temp_dir:
            temp_root = Path(temp_dir)
            input_path = temp_root / "input.rds"
            output_path = temp_root / "output.rds"
            log_path = temp_root / "rctd.log"

            with r_env.local_context(
                anndata=True, pandas=True, numpy=True
            ) as r_context:
                r_context["spatial_counts"] = spatial_data.X.T
                r_context["reference_counts"] = reference_data.X.T

                r_context["gene_names_spatial"] = ro.StrVector(spatial_data.var_names)
                r_context["spot_names"] = ro.StrVector(spatial_data.obs_names)
                r_context["gene_names_ref"] = ro.StrVector(reference_data.var_names)
                r_context["cell_names"] = ro.StrVector(reference_data.obs_names)
                r_context["coords"] = ro.conversion.py2rpy(coords)
                r_context["numi_spatial"] = ro.conversion.py2rpy(spatial_numi)
                r_context["cell_types_vec"] = ro.conversion.py2rpy(cell_types_series)
                r_context["numi_ref"] = ro.conversion.py2rpy(reference_numi)
                r_context["max_cores_val"] = max_cores
                r_context["conf_thresh"] = confidence_threshold
                r_context["doub_thresh"] = doublet_threshold
                r_context["max_multi_types_val"] = max_multi_types
                r_context["rctd_input_path"] = str(input_path)

                ro.r(
                    """
                    invisible(suppressMessages(suppressWarnings(capture.output({
                        rownames(spatial_counts) <- gene_names_spatial
                        colnames(spatial_counts) <- spot_names
                        rownames(reference_counts) <- gene_names_ref
                        colnames(reference_counts) <- cell_names
                        puck <- SpatialRNA(coords, spatial_counts, numi_spatial)
                        cell_types_factor <- as.factor(cell_types_vec)
                        names(cell_types_factor) <- names(cell_types_vec)
                        reference <- Reference(
                            reference_counts,
                            cell_types_factor,
                            numi_ref,
                            min_UMI = 5
                        )
                        myRCTD <- create.RCTD(
                            puck,
                            reference,
                            max_cores = max_cores_val,
                            MAX_MULTI_TYPES = max_multi_types_val,
                            UMI_min_sigma = 10
                        )
                        myRCTD@config$CONFIDENCE_THRESHOLD <- conf_thresh
                        myRCTD@config$DOUBLET_THRESHOLD <- doub_thresh
                        # run.RCTD samples sigma candidates. Transfer the active
                        # RNG state so process isolation preserves R semantics.
                        random_seed <- get0(
                            ".Random.seed",
                            envir = .GlobalEnv,
                            inherits = FALSE
                        )
                        saveRDS(
                            list(rctd = myRCTD, random_seed = random_seed),
                            rctd_input_path
                        )
                    }))))
                """
                )

                _run_rctd_subprocess(
                    ro,
                    input_path,
                    output_path,
                    log_path,
                    mode,
                )

                r_context["rctd_output_path"] = str(output_path)
                ro.r(
                    """
                    rctd_bundle <- readRDS(rctd_output_path)
                    myRCTD <- rctd_bundle$rctd
                    if (!is.null(rctd_bundle$random_seed)) {
                        assign(
                            ".Random.seed",
                            rctd_bundle$random_seed,
                            envir = .GlobalEnv
                        )
                    }
                """
                )
                proportions = _extract_rctd_results(mode, ro)

        # Validate results
        if proportions.isna().any().any():
            nan_count = proportions.isna().sum().sum()
            warnings.warn(
                f"RCTD produced {nan_count} NaN values", UserWarning, stacklevel=2
            )

        if (proportions < 0).any().any():
            neg_count = (proportions < 0).sum().sum()
            raise ProcessingError(f"RCTD error: {neg_count} negative values")

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            data.common_genes,
            method=f"RCTD-{mode}",
            device="CPU",
            mode=mode,
            max_cores=max_cores,
            confidence_threshold=confidence_threshold,
            doublet_threshold=doublet_threshold,
        )

        return proportions, stats

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"RCTD deconvolution failed: {e}") from e


def _extract_rctd_results(mode: str, robjects: Any) -> pd.DataFrame:
    """Extract RCTD results; caller must hold the R lock/converter context."""
    ro = robjects

    if mode == "full":
        ro.r(
            """
            weights_matrix <- myRCTD@results$weights
            cell_type_names <- myRCTD@cell_type_info$renorm[[2]]
            spot_names <- rownames(weights_matrix)
        """
        )
    elif mode == "doublet":
        ro.r(
            """
            if("weights_doublet" %in% names(myRCTD@results) && "results_df" %in% names(myRCTD@results)) {
                weights_doublet <- myRCTD@results$weights_doublet
                results_df <- myRCTD@results$results_df
                cell_type_names <- myRCTD@cell_type_info$renorm[[2]]
                spot_names <- rownames(results_df)
                n_spots <- length(spot_names)
                n_cell_types <- length(cell_type_names)
                weights_matrix <- matrix(0, nrow = n_spots, ncol = n_cell_types)
                rownames(weights_matrix) <- spot_names
                colnames(weights_matrix) <- cell_type_names
                for(i in 1:n_spots) {
                    spot_class <- results_df$spot_class[i]
                    if(spot_class %in% c("doublet_certain", "doublet_uncertain")) {
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
                        first_type <- as.character(results_df$first_type[i])
                        if(first_type %in% cell_type_names) {
                            first_idx <- which(cell_type_names == first_type)
                            weights_matrix[i, first_idx] <- 1.0
                        }
                    }
                }
            } else {
                stop("Official doublet mode structures not found")
            }
        """
        )
    else:  # multi mode
        ro.r(
            """
            results_list <- myRCTD@results
            spot_names <- colnames(myRCTD@spatialRNA@counts)
            cell_type_names <- myRCTD@cell_type_info$renorm[[2]]
            n_spots <- length(spot_names)
            n_cell_types <- length(cell_type_names)
            weights_matrix <- matrix(0, nrow = n_spots, ncol = n_cell_types)
            rownames(weights_matrix) <- spot_names
            colnames(weights_matrix) <- cell_type_names
            for(i in 1:n_spots) {
                spot_result <- results_list[[i]]
                predicted_types <- spot_result$cell_type_list
                proportions <- spot_result$sub_weights
                for(j in seq_along(predicted_types)) {
                    cell_type <- predicted_types[j]
                    if(cell_type %in% cell_type_names) {
                        col_idx <- which(cell_type_names == cell_type)
                        weights_matrix[i, col_idx] <- proportions[j]
                    }
                }
            }
        """
        )

    weights_r = ro.r("as.matrix(weights_matrix)")
    cell_type_names_r = ro.r("cell_type_names")
    spot_names_r = ro.r("spot_names")

    weights_array = ro.conversion.rpy2py(weights_r)
    cell_type_names = ro.conversion.rpy2py(cell_type_names_r)
    spot_names = ro.conversion.rpy2py(spot_names_r)

    return pd.DataFrame(weights_array, index=spot_names, columns=cell_type_names)
