"""
CARD (Conditional AutoRegressive-based Deconvolution) method.

CARD models spatial correlation in cell type composition using a
CAR (Conditional AutoRegressive) model. Unique features:
- Spatial correlation modeling
- Optional high-resolution imputation
"""

from typing import Any, Optional

import numpy as np
import pandas as pd

from ...utils.dependency_manager import validate_r_environment
from ...utils.exceptions import ChatSpatialError, DataError, ProcessingError
from .base import PreparedDeconvolutionData, create_deconvolution_stats


def deconvolve(
    data: PreparedDeconvolutionData,
    sample_key: Optional[str] = None,
    minCountGene: int = 100,
    minCountSpot: int = 5,
    imputation: bool = False,
    NumGrids: int = 2000,
    ineibor: int = 10,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Deconvolve spatial data using CARD R package.

    Args:
        data: Prepared deconvolution data (immutable, includes spatial coordinates)
        sample_key: Optional sample/batch key in reference data
        minCountGene: Include genes with at least this many counts
        minCountSpot: Include genes expressed in at least this many spots
        imputation: Whether to perform spatial imputation
        NumGrids: Number of grids for imputation
        ineibor: Number of neighbors for imputation

    Returns:
        Tuple of (proportions DataFrame, statistics dictionary)
    """
    ctx = data.ctx

    try:
        # Data already copied in prepare_deconvolution
        spatial_data = data.spatial
        reference_data = data.reference

        # Get spatial coordinates from prepared data
        if data.spatial_coords is None:
            raise DataError(
                "CARD requires real spatial coordinates for spatially-informed "
                "deconvolution. No spatial coordinates found. "
                "Use a non-spatial method (e.g., NNLS, DestVI) instead."
            )
        spatial_location = pd.DataFrame(
            data.spatial_coords[:, :2],
            index=spatial_data.obs_names,
            columns=["x", "y"],
        )

        # Prepare metadata
        sc_meta = reference_data.obs[[data.cell_type_key]].copy()
        sc_meta.columns = ["cellType"]

        if sample_key and sample_key in reference_data.obs:
            sc_meta["sampleInfo"] = reference_data.obs[sample_key]
        else:
            sc_meta["sampleInfo"] = "sample1"

        r_env = validate_r_environment(
            ctx,
            required_packages=["CARD"],
            require_anndata2ri=True,
            package_install_commands={
                "CARD": "devtools::install_github('YingMa0107/CARD')"
            },
        )
        ro = r_env.robjects

        imputed_proportions = None
        imputed_coordinates = None

        with r_env.local_context(
            anndata=True, pandas=True, numpy=True
        ) as r_context:
            r_context["sc_count"] = reference_data.X.T
            r_context["spatial_count"] = spatial_data.X.T

            r_context["gene_names_ref"] = ro.StrVector(reference_data.var_names)
            r_context["cell_names"] = ro.StrVector(reference_data.obs_names)
            r_context["gene_names_spatial"] = ro.StrVector(spatial_data.var_names)
            r_context["spot_names"] = ro.StrVector(spatial_data.obs_names)

            ro.r("""
                    rownames(sc_count) <- gene_names_ref
                    colnames(sc_count) <- cell_names
                    rownames(spatial_count) <- gene_names_spatial
                    colnames(spatial_count) <- spot_names
                """)

            r_context["sc_meta"] = ro.conversion.py2rpy(sc_meta)
            r_context["spatial_location"] = ro.conversion.py2rpy(spatial_location)
            r_context["minCountGene"] = minCountGene
            r_context["minCountSpot"] = minCountSpot

            ro.r("""
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
                """)

            row_names = list(ro.r("rownames(CARD_obj@Proportion_CARD)"))
            col_names = list(ro.r("colnames(CARD_obj@Proportion_CARD)"))
            proportions_r = ro.r("CARD_obj@Proportion_CARD")
            proportions_array = np.array(proportions_r)

            proportions = pd.DataFrame(
                proportions_array, index=row_names, columns=col_names
            )

            if imputation:
                ro.r(f"""
                        capture.output(
                            CARD_impute <- CARD.imputation(
                                CARD_object = CARD_obj,
                                NumGrids = {NumGrids},
                                ineibor = {ineibor}
                            ),
                            file = "/dev/null"
                        )
                """)

                imputed_row_names = list(ro.r("rownames(CARD_impute@refined_prop)"))
                imputed_col_names = list(ro.r("colnames(CARD_impute@refined_prop)"))
                imputed_proportions_r = ro.r("CARD_impute@refined_prop")
                imputed_proportions_array = np.array(imputed_proportions_r)

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
            data.common_genes,
            method="CARD",
            device="CPU",
            minCountGene=minCountGene,
            minCountSpot=minCountSpot,
        )

        if imputation and imputed_proportions is not None:
            stats["imputation"] = {
                "enabled": True,
                "n_imputed_locations": len(imputed_proportions),
                "resolution_increase": (
                    f"{len(imputed_proportions) / len(row_names):.1f}x"
                ),
                "imputed_proportions": imputed_proportions,
                "imputed_coordinates": imputed_coordinates,
            }

        return proportions, stats

    except ChatSpatialError:
        raise
    except Exception as e:
        raise ProcessingError(f"CARD deconvolution failed: {e}") from e
