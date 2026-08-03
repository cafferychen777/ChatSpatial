"""
SPOTlight deconvolution method.

SPOTlight is an R-based deconvolution method that uses NMF
(Non-negative Matrix Factorization) for cell type decomposition.
"""

from typing import Any

import numpy as np
import pandas as pd

from ...utils.adata_utils import to_dense
from ...utils.dependency_manager import validate_r_environment
from ...utils.exceptions import ChatSpatialError, DataError, ProcessingError
from .base import PreparedDeconvolutionData, create_deconvolution_stats


def deconvolve(
    data: PreparedDeconvolutionData,
    n_top_genes: int = 50,
    nmf_model: str = "ns",
    min_prop: float = 0.01,
    scale: bool = True,
    weight_id: str = "mean.AUC",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Deconvolve spatial data using SPOTlight R package.

    Args:
        data: Prepared deconvolution data (immutable, includes spatial coordinates)
        n_top_genes: Number of marker genes per cell type
        nmf_model: NMF model type - 'ns' (non-smooth) or 'std' (standard)
        min_prop: Minimum proportion threshold
        scale: Whether to scale data
        weight_id: Column name for marker gene weights

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    ctx = data.ctx

    try:
        # Validate spatial coordinates from prepared data
        if data.spatial_coords is None:
            raise DataError(
                "SPOTlight requires spatial coordinates. "
                "Ensure spatial data has 'spatial' key in obsm."
            )
        spatial_coords = data.spatial_coords

        # Data already copied in prepare_deconvolution
        spatial_data = data.spatial
        reference_data = data.reference

        # Ensure integer counts for R interface
        dense = to_dense(spatial_data.X)
        spatial_counts = (
            dense.astype(np.int32, copy=False) if dense.dtype != np.int32 else dense
        )

        dense = to_dense(reference_data.X)
        reference_counts = (
            dense.astype(np.int32, copy=False) if dense.dtype != np.int32 else dense
        )

        # R factors support arbitrary string values; preserve biological labels.
        cell_types = reference_data.obs[data.cell_type_key].astype(str)

        required_packages = [
            "SPOTlight",
            "SingleCellExperiment",
            "SpatialExperiment",
            "scran",
            "scuttle",
        ]
        r_env = validate_r_environment(
            ctx,
            required_packages=required_packages,
            package_install_commands={
                package: f"BiocManager::install('{package}')"
                for package in required_packages
            },
        )
        ro = r_env.robjects

        with r_env.local_context(pandas=True, numpy=True) as r_context:
            r_context["spatial_counts"] = spatial_counts.T
            r_context["reference_counts"] = reference_counts.T
            r_context["spatial_coords"] = spatial_coords
            r_context["gene_names"] = ro.StrVector(data.common_genes)
            r_context["spatial_names"] = ro.StrVector(list(spatial_data.obs_names))
            r_context["reference_names"] = ro.StrVector(list(reference_data.obs_names))
            r_context["cell_types"] = ro.StrVector(cell_types.tolist())
            r_context["nmf_model"] = nmf_model
            r_context["min_prop"] = min_prop
            r_context["scale_data"] = scale
            r_context["weight_id"] = weight_id
            r_context["n_top_genes"] = n_top_genes

            ro.r(
                """
                    sce <- SingleCellExperiment(
                        assays = list(counts = reference_counts),
                        colData = data.frame(
                            cell_type = factor(cell_types),
                            row.names = reference_names
                        )
                    )
                    rownames(sce) <- gene_names
                    sce <- logNormCounts(sce)

                    spe <- SpatialExperiment(
                        assays = list(counts = spatial_counts),
                        spatialCoords = spatial_coords,
                        colData = data.frame(row.names = spatial_names)
                    )
                    rownames(spe) <- gene_names
                    colnames(spe) <- spatial_names

                    markers <- findMarkers(sce, groups = sce$cell_type, test.type = "wilcox")
                    cell_type_names <- names(markers)
                    mgs_list <- list()

                    for (ct in cell_type_names) {
                        ct_markers <- markers[[ct]]
                        n_markers <- min(n_top_genes, nrow(ct_markers))
                        top_markers <- head(ct_markers[order(ct_markers$p.value), ], n_markers)

                        mgs_df <- data.frame(
                            gene = rownames(top_markers),
                            cluster = ct,
                            mean.AUC = -log10(top_markers$p.value + 1e-10)
                        )
                        mgs_list[[ct]] <- mgs_df
                    }

                    mgs <- do.call(rbind, mgs_list)

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
                        verbose = FALSE
                    )
                """
            )

            proportions_np = np.array(ro.r("spotlight_result$mat"))
            spot_names = list(ro.r("rownames(spotlight_result$mat)"))
            cell_type_names = list(ro.r("colnames(spotlight_result$mat)"))

            proportions = pd.DataFrame(
                proportions_np, index=spot_names, columns=cell_type_names
            )

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            data.common_genes,
            method="SPOTlight",
            device="CPU",
            n_top_genes=n_top_genes,
            nmf_model=nmf_model,
            min_prop=min_prop,
        )

        return proportions, stats

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"SPOTlight deconvolution failed: {e}") from e
