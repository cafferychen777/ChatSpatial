"""
RCTD (Robust Cell Type Decomposition) deconvolution method.

RCTD is an R-based deconvolution method that performs robust
decomposition of cell type mixtures via the spacexr package.
"""

from typing import TYPE_CHECKING, Any, Dict, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    pass

from ...utils.adata_utils import get_spatial_key
from ...utils.dependency_manager import validate_r_package
from ...utils.exceptions import DataError, ParameterError, ProcessingError
from .base import DeconvolutionContext, create_deconvolution_stats


async def deconvolve(
    deconv_ctx: DeconvolutionContext,
    mode: str = "full",
    max_cores: int = 4,
    confidence_threshold: float = 10.0,
    doublet_threshold: float = 25.0,
    max_multi_types: int = 4,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using RCTD from spacexr R package.

    Args:
        deconv_ctx: Prepared DeconvolutionContext
        mode: RCTD mode - 'full', 'doublet', or 'multi'
        max_cores: Maximum CPU cores
        confidence_threshold: Confidence threshold
        doublet_threshold: Doublet detection threshold
        max_multi_types: Max cell types per spot in multi mode

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    import anndata2ri
    import rpy2.robjects as ro
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    ctx = deconv_ctx.ctx
    cell_type_key = deconv_ctx.cell_type_key

    # Validate mode-specific parameters
    if mode == "multi" and max_multi_types >= len(deconv_ctx.cell_types):
        raise ParameterError(
            f"MAX_MULTI_TYPES ({max_multi_types}) must be less than "
            f"total cell types ({len(deconv_ctx.cell_types)})."
        )

    # Validate R package
    validate_r_package(
        "spacexr",
        ctx,
        install_cmd="devtools::install_github('dmcable/spacexr', build_vignettes = FALSE)",
    )

    try:
        # Load R packages
        with localconverter(ro.default_converter + pandas2ri.converter):
            rpackages.importr("spacexr")
            rpackages.importr("base")

        # Get subset data
        spatial_data, reference_data = deconv_ctx.get_subset_data()
        common_genes = deconv_ctx.common_genes

        # Get spatial coordinates
        spatial_key = get_spatial_key(deconv_ctx.spatial_adata)
        if spatial_key:
            coords = pd.DataFrame(
                deconv_ctx.spatial_adata.obsm[spatial_key],
                index=deconv_ctx.spatial_adata.obs_names,
                columns=["x", "y"],
            )
        else:
            coords = pd.DataFrame(
                {"x": range(spatial_data.n_obs), "y": [0] * spatial_data.n_obs},
                index=spatial_data.obs_names,
            )

        # Prepare cell type information
        cell_types = reference_data.obs[cell_type_key].copy()
        cell_types = cell_types.str.replace("/", "_", regex=False)
        cell_types = cell_types.str.replace(" ", "_", regex=False)

        # RCTD requires minimum 25 cells per cell type
        # Filter out rare cell types to avoid R-side errors
        MIN_CELLS_PER_TYPE = 25
        cell_type_counts = cell_types.value_counts()
        rare_types = cell_type_counts[
            cell_type_counts < MIN_CELLS_PER_TYPE
        ].index.tolist()

        if rare_types:
            await ctx.warning(
                f"RCTD requires ≥{MIN_CELLS_PER_TYPE} cells per cell type. "
                f"Filtering {len(rare_types)} rare types: {rare_types}"
            )
            # Filter reference data to exclude rare cell types
            keep_mask = ~cell_types.isin(rare_types)
            reference_data = reference_data[keep_mask].copy()
            cell_types = cell_types[keep_mask]

            # Check if enough cell types remain
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
        # Memory optimization: ravel() returns view when possible, flatten() always copies
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

        # Transfer matrices to R
        with localconverter(ro.default_converter + anndata2ri.converter):
            ro.globalenv["spatial_counts"] = spatial_data.X.T
            ro.globalenv["reference_counts"] = reference_data.X.T

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

        # Transfer other data
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.globalenv["coords"] = ro.conversion.py2rpy(coords)
            ro.globalenv["numi_spatial"] = ro.conversion.py2rpy(spatial_numi)
            ro.globalenv["cell_types_vec"] = ro.conversion.py2rpy(cell_types_series)
            ro.globalenv["numi_ref"] = ro.conversion.py2rpy(reference_numi)
            ro.globalenv["max_cores_val"] = max_cores
            ro.globalenv["rctd_mode"] = mode
            ro.globalenv["conf_thresh"] = confidence_threshold
            ro.globalenv["doub_thresh"] = doublet_threshold
            ro.globalenv["max_multi_types_val"] = max_multi_types

        # Run RCTD in R
        ro.r(
            """
            puck <- SpatialRNA(coords, spatial_counts, numi_spatial)
            cell_types_factor <- as.factor(cell_types_vec)
            names(cell_types_factor) <- names(cell_types_vec)
            reference <- Reference(reference_counts, cell_types_factor, numi_ref, min_UMI = 5)
            myRCTD <- create.RCTD(puck, reference, max_cores = max_cores_val,
                                  MAX_MULTI_TYPES = max_multi_types_val, UMI_min_sigma = 10)
            myRCTD@config$CONFIDENCE_THRESHOLD <- conf_thresh
            myRCTD@config$DOUBLET_THRESHOLD <- doub_thresh
            myRCTD <- run.RCTD(myRCTD, doublet_mode = rctd_mode)
        """
        )

        # Extract results
        proportions = _extract_rctd_results(mode)

        # Validate results
        if proportions.isna().any().any():
            nan_count = proportions.isna().sum().sum()
            await ctx.warning(f"RCTD produced {nan_count} NaN values")

        if (proportions < 0).any().any():
            neg_count = (proportions < 0).sum().sum()
            raise ProcessingError(f"RCTD error: {neg_count} negative values")

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            common_genes,
            f"RCTD-{mode}",
            "CPU",
            mode=mode,
            max_cores=max_cores,
            confidence_threshold=confidence_threshold,
            doublet_threshold=doublet_threshold,
        )

        # Clean up R global environment to free memory
        ro.r(
            """
            rm(list = c("spatial_counts", "reference_counts", "gene_names_spatial",
                        "spot_names", "gene_names_ref", "cell_names", "coords",
                        "numi_spatial", "cell_types_vec", "numi_ref", "max_cores_val",
                        "rctd_mode", "conf_thresh", "doub_thresh", "max_multi_types_val",
                        "puck", "cell_types_factor", "reference", "myRCTD",
                        "weights_matrix", "cell_type_names"),
                   envir = .GlobalEnv)
            gc()
        """
        )

        return proportions, stats

    except Exception as e:
        if isinstance(e, (ParameterError, ProcessingError)):
            raise
        raise ProcessingError(f"RCTD deconvolution failed: {str(e)}") from e


def _extract_rctd_results(mode: str) -> pd.DataFrame:
    """Extract RCTD results from R environment."""
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    with localconverter(
        ro.default_converter + pandas2ri.converter + numpy2ri.converter
    ):
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
