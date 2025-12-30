"""
SPOTlight deconvolution method.

SPOTlight is an R-based deconvolution method that uses NMF
(Non-negative Matrix Factorization) for cell type decomposition.
"""

from typing import TYPE_CHECKING, Any, Dict, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    pass

from ...utils.adata_utils import require_spatial_coords, to_dense
from ...utils.dependency_manager import validate_r_package
from ...utils.exceptions import ProcessingError
from .base import DeconvolutionContext, create_deconvolution_stats


async def deconvolve(
    deconv_ctx: DeconvolutionContext,
    n_top_genes: int = 2000,
    nmf_model: str = "ns",
    min_prop: float = 0.01,
    scale: bool = True,
    weight_id: str = "mean.AUC",
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using SPOTlight R package.

    Args:
        deconv_ctx: Prepared DeconvolutionContext
        n_top_genes: Number of top HVGs to use
        nmf_model: NMF model type - 'ns' (non-smooth) or 'std' (standard)
        min_prop: Minimum proportion threshold
        scale: Whether to scale data
        weight_id: Column name for marker gene weights

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    import anndata2ri
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter

    ctx = deconv_ctx.ctx
    cell_type_key = deconv_ctx.cell_type_key

    # Validate R package
    validate_r_package(
        "SPOTlight",
        ctx,
        install_cmd="BiocManager::install('SPOTlight')",
    )

    try:
        # Validate spatial coordinates
        spatial_coords = require_spatial_coords(deconv_ctx.spatial_adata)

        # Get subset data
        spatial_data, reference_data = deconv_ctx.get_subset_data()
        common_genes = deconv_ctx.common_genes

        # Ensure integer counts for R interface
        # Note: DeconvolutionContext.prepare() already converted to int32 for R-based methods
        # but we convert here explicitly to handle any edge cases
        dense = to_dense(spatial_data.X)
        spatial_counts = (
            dense.astype(np.int32, copy=False) if dense.dtype != np.int32 else dense
        )

        dense = to_dense(reference_data.X)
        reference_counts = (
            dense.astype(np.int32, copy=False) if dense.dtype != np.int32 else dense
        )

        # Clean cell type labels
        cell_types = reference_data.obs[cell_type_key].astype(str)
        cell_types = cell_types.str.replace("/", "_", regex=False)
        cell_types = cell_types.str.replace(" ", "_", regex=False)

        # Transfer matrices to R using anndata2ri
        with localconverter(ro.default_converter + anndata2ri.converter):
            ro.globalenv["spatial_counts"] = spatial_counts.T
            ro.globalenv["reference_counts"] = reference_counts.T

        # Transfer other data
        with localconverter(
            ro.default_converter + pandas2ri.converter + numpy2ri.converter
        ):
            ro.r("library(SPOTlight)")
            ro.r("library(SingleCellExperiment)")
            ro.r("library(SpatialExperiment)")
            ro.r("library(scran)")
            ro.r("library(scuttle)")

            ro.globalenv["spatial_coords"] = spatial_coords
            ro.globalenv["gene_names"] = ro.StrVector(common_genes)
            ro.globalenv["spatial_names"] = ro.StrVector(list(spatial_data.obs_names))
            ro.globalenv["reference_names"] = ro.StrVector(
                list(reference_data.obs_names)
            )
            ro.globalenv["cell_types"] = ro.StrVector(cell_types.tolist())
            ro.globalenv["nmf_model"] = nmf_model
            ro.globalenv["min_prop"] = min_prop
            ro.globalenv["scale_data"] = scale
            ro.globalenv["weight_id"] = weight_id

        # Create SCE and SPE objects, run SPOTlight
        ro.r(
            """
            # Create SingleCellExperiment for reference
            sce <- SingleCellExperiment(
                assays = list(counts = reference_counts),
                colData = data.frame(
                    cell_type = factor(cell_types),
                    row.names = reference_names
                )
            )
            rownames(sce) <- gene_names
            sce <- logNormCounts(sce)

            # Create SpatialExperiment for spatial data
            spe <- SpatialExperiment(
                assays = list(counts = spatial_counts),
                spatialCoords = spatial_coords,
                colData = data.frame(row.names = spatial_names)
            )
            rownames(spe) <- gene_names
            colnames(spe) <- spatial_names

            # Find marker genes using scran
            markers <- findMarkers(sce, groups = sce$cell_type, test.type = "wilcox")

            # Format marker genes for SPOTlight
            cell_type_names <- names(markers)
            mgs_list <- list()

            for (ct in cell_type_names) {
                ct_markers <- markers[[ct]]
                n_markers <- min(50, nrow(ct_markers))
                top_markers <- head(ct_markers[order(ct_markers$p.value), ], n_markers)

                mgs_df <- data.frame(
                    gene = rownames(top_markers),
                    cluster = ct,
                    mean.AUC = -log10(top_markers$p.value + 1e-10)
                )
                mgs_list[[ct]] <- mgs_df
            }

            mgs <- do.call(rbind, mgs_list)

            # Run SPOTlight
            spotlight_result <- SPOTlight(
                x = sce,
                y = spe,
                groups = sce$cell_type,
                mgs = mgs,
                weight_id = weight_id,
                group_id = "cluster",
                gene_id = "gene",
                model = nmf_model,
                min_prop = min_prop,
                scale = scale_data,
                verbose = TRUE
            )
        """
        )

        # Extract results
        with localconverter(
            ro.default_converter + pandas2ri.converter + numpy2ri.converter
        ):
            proportions_np = np.array(ro.r("spotlight_result$mat"))
            spot_names = list(ro.r("rownames(spotlight_result$mat)"))
            cell_type_names = list(ro.r("colnames(spotlight_result$mat)"))

        proportions = pd.DataFrame(
            proportions_np, index=spot_names, columns=cell_type_names
        )

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            common_genes,
            "SPOTlight",
            "CPU",
            n_top_genes=n_top_genes,
            nmf_model=nmf_model,
            min_prop=min_prop,
        )

        # Clean up R global environment to free memory
        ro.r(
            """
            rm(list = c("spatial_counts", "reference_counts", "spatial_coords",
                        "gene_names", "spatial_names", "reference_names", "cell_types",
                        "nmf_model", "min_prop", "scale_data", "weight_id",
                        "sce", "spe", "markers", "mgs", "spotlight_result"),
                   envir = .GlobalEnv)
            gc()
        """
        )

        return proportions, stats

    except Exception as e:
        if isinstance(e, ProcessingError):
            raise
        raise ProcessingError(f"SPOTlight deconvolution failed: {str(e)}") from e
