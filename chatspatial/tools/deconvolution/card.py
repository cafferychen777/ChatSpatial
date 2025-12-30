"""
CARD (Conditional AutoRegressive-based Deconvolution) method.

CARD models spatial correlation in cell type composition using a
CAR (Conditional AutoRegressive) model. Unique features:
- Spatial correlation modeling
- Optional high-resolution imputation
"""

from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    pass

from ...utils.adata_utils import get_spatial_key
from ...utils.dependency_manager import validate_r_package
from ...utils.exceptions import ProcessingError
from .base import DeconvolutionContext, create_deconvolution_stats


async def deconvolve(
    deconv_ctx: DeconvolutionContext,
    sample_key: Optional[str] = None,
    minCountGene: int = 100,
    minCountSpot: int = 5,
    imputation: bool = False,
    NumGrids: int = 2000,
    ineibor: int = 10,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Deconvolve spatial data using CARD R package.

    Args:
        deconv_ctx: Prepared DeconvolutionContext
        sample_key: Optional sample/batch key in reference data
        minCountGene: Include genes with at least this many counts
        minCountSpot: Include genes expressed in at least this many spots
        imputation: Whether to perform spatial imputation
        NumGrids: Number of grids for imputation
        ineibor: Number of neighbors for imputation

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
        "CARD",
        ctx,
        install_cmd="devtools::install_github('YingMa0107/CARD')",
    )

    try:
        # Load CARD
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.r("library(CARD)")

        # Get subset data
        spatial_data, reference_data = deconv_ctx.get_subset_data()
        common_genes = deconv_ctx.common_genes

        # Get spatial coordinates
        spatial_key = get_spatial_key(deconv_ctx.spatial_adata)
        if spatial_key:
            spatial_location = pd.DataFrame(
                deconv_ctx.spatial_adata.obsm[spatial_key],
                index=deconv_ctx.spatial_adata.obs_names,
                columns=["x", "y"],
            )
        else:
            spatial_location = pd.DataFrame(
                {"x": range(spatial_data.n_obs), "y": [0] * spatial_data.n_obs},
                index=spatial_data.obs_names,
            )

        # Prepare metadata
        sc_meta = reference_data.obs[[cell_type_key]].copy()
        sc_meta.columns = ["cellType"]

        if sample_key and sample_key in reference_data.obs:
            sc_meta["sampleInfo"] = reference_data.obs[sample_key]
        else:
            sc_meta["sampleInfo"] = "sample1"

        # Transfer matrices to R
        with localconverter(ro.default_converter + anndata2ri.converter):
            ro.globalenv["sc_count"] = reference_data.X.T
            ro.globalenv["spatial_count"] = spatial_data.X.T

            ro.globalenv["gene_names_ref"] = ro.StrVector(reference_data.var_names)
            ro.globalenv["cell_names"] = ro.StrVector(reference_data.obs_names)
            ro.globalenv["gene_names_spatial"] = ro.StrVector(spatial_data.var_names)
            ro.globalenv["spot_names"] = ro.StrVector(spatial_data.obs_names)

            ro.r(
                """
                rownames(sc_count) <- gene_names_ref
                colnames(sc_count) <- cell_names
                rownames(spatial_count) <- gene_names_spatial
                colnames(spatial_count) <- spot_names
            """
            )

        # Transfer metadata
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.globalenv["sc_meta"] = ro.conversion.py2rpy(sc_meta)
            ro.globalenv["spatial_location"] = ro.conversion.py2rpy(spatial_location)
            ro.globalenv["minCountGene"] = minCountGene
            ro.globalenv["minCountSpot"] = minCountSpot

        # Create CARD object and run deconvolution (suppress R output)
        ro.r(
            """
            capture.output(
                CARD_obj <- createCARDObject(
                    sc_count = sc_count,
                    sc_meta = sc_meta,
                    spatial_count = spatial_count,
                    spatial_location = spatial_location,
                    ct.varname = "cellType",
                    ct.select = unique(sc_meta$cellType),
                    sample.varname = "sampleInfo",
                    minCountGene = minCountGene,
                    minCountSpot = minCountSpot
                ),
                file = "/dev/null"
            )
            capture.output(
                CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj),
                file = "/dev/null"
            )
        """
        )

        # Extract results
        with localconverter(
            ro.default_converter + pandas2ri.converter + numpy2ri.converter
        ):
            row_names = list(ro.r("rownames(CARD_obj@Proportion_CARD)"))
            col_names = list(ro.r("colnames(CARD_obj@Proportion_CARD)"))
            proportions_r = ro.r("CARD_obj@Proportion_CARD")
            proportions_array = np.array(proportions_r)

            proportions = pd.DataFrame(
                proportions_array, index=row_names, columns=col_names
            )

        # Optional imputation
        imputed_proportions = None
        imputed_coordinates = None

        if imputation:
            ro.r(
                f"""
                capture.output(
                    CARD_impute <- CARD.imputation(
                        CARD_object = CARD_obj,
                        NumGrids = {NumGrids},
                        ineibor = {ineibor}
                    ),
                    file = "/dev/null"
                )
            """
            )

            with localconverter(ro.default_converter + pandas2ri.converter):
                imputed_row_names = list(ro.r("rownames(CARD_impute@refined_prop)"))
                imputed_col_names = list(ro.r("colnames(CARD_impute@refined_prop)"))
                imputed_proportions_r = ro.r("CARD_impute@refined_prop")
                imputed_proportions_array = np.array(imputed_proportions_r)

                # Parse coordinates from rownames
                coords_list = []
                for name in imputed_row_names:
                    parts = name.split("x")
                    coords_list.append([float(parts[0]), float(parts[1])])

                imputed_proportions = pd.DataFrame(
                    imputed_proportions_array,
                    index=imputed_row_names,
                    columns=imputed_col_names,
                )
                imputed_coordinates = pd.DataFrame(
                    coords_list, index=imputed_row_names, columns=["x", "y"]
                )

        # Create statistics
        stats = create_deconvolution_stats(
            proportions,
            common_genes,
            "CARD",
            "CPU",
            minCountGene=minCountGene,
            minCountSpot=minCountSpot,
        )

        if imputation and imputed_proportions is not None:
            stats["imputation"] = {
                "enabled": True,
                "n_imputed_locations": len(imputed_proportions),
                "resolution_increase": f"{len(imputed_proportions) / len(row_names):.1f}x",
                "imputed_proportions": imputed_proportions,
                "imputed_coordinates": imputed_coordinates,
            }

        # Clean up R global environment to free memory
        cleanup_vars = [
            "sc_count",
            "spatial_count",
            "gene_names_ref",
            "cell_names",
            "gene_names_spatial",
            "spot_names",
            "sc_meta",
            "spatial_location",
            "minCountGene",
            "minCountSpot",
            "CARD_obj",
        ]
        if imputation:
            cleanup_vars.append("CARD_impute")

        ro.r(
            f"""
            rm(list = c({', '.join(f'"{v}"' for v in cleanup_vars)}),
               envir = .GlobalEnv)
            gc()
        """
        )

        return proportions, stats

    except Exception as e:
        if isinstance(e, ProcessingError):
            raise
        raise ProcessingError(f"CARD deconvolution failed: {str(e)}") from e
