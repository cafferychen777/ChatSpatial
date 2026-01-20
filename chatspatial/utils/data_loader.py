"""
Data loading utilities for spatial transcriptomics data.

Handles loading various spatial data formats:
- H5AD files (AnnData format)
- H5 files (10x Genomics format)
- MTX directories (10x Visium structure)
- Visium directories with spatial information
- Xenium directories with cell-level spatial data

For data persistence, see persistence.py.
"""

import logging
import os
from typing import TYPE_CHECKING, Any, Optional, cast

from ..models.data import SpatialPlatform

if TYPE_CHECKING:
    from zarr.core.array import Array
    from zarr.core.group import Group
from .adata_utils import ensure_unique_var_names, get_adata_profile
from .dependency_manager import is_available
from .exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)

logger = logging.getLogger(__name__)


def _load_xenium_zarr(data_path: str) -> Any:
    """Load Xenium data from zarr format.

    Args:
        data_path: Path to Xenium output directory containing zarr files

    Returns:
        AnnData object with expression data and spatial coordinates
    """
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import zarr
    from zarr.storage import ZipStore

    matrix_zarr = os.path.join(data_path, "cell_feature_matrix.zarr.zip")
    cells_zarr = os.path.join(data_path, "cells.zarr.zip")

    # Load cell_feature_matrix.zarr
    # Note: zarr.open() returns AnyArray | Group, we cast to Group for type safety
    matrix_store = ZipStore(matrix_zarr, mode="r")
    matrix_root = cast("Group", zarr.open(matrix_store, mode="r"))
    cf = cast("Group", matrix_root["cell_features"])

    # Get sparse matrix components (CSC format)
    # cf["data"], cf["indices"], cf["indptr"] are zarr Arrays
    data = np.asarray(cast("Array", cf["data"])[:])
    indices = np.asarray(cast("Array", cf["indices"])[:])
    indptr = np.asarray(cast("Array", cf["indptr"])[:])

    # attrs values are JSON type; Xenium format guarantees these are integers
    n_cells = cast(int, cf.attrs["number_cells"])
    n_features = cast(int, cf.attrs["number_features"])

    # Create CSC matrix then convert to CSR (cells x genes)
    X_csc = sp.csc_matrix((data, indices, indptr), shape=(n_cells, n_features))
    X = X_csc.tocsr()

    # Get feature names from zarr attrs (JSON type -> list[str])
    # Xenium format guarantees these are string lists
    feature_keys = cast(list[str], cf.attrs.get("feature_keys", []))
    feature_ids_raw = cast(list[str], cf.attrs.get("feature_ids", []))
    feature_names: list[str] = list(feature_keys)
    feature_ids: list[str] = list(feature_ids_raw)

    # Load cells.zarr for spatial coordinates
    cells_store = ZipStore(cells_zarr, mode="r")
    cells_root = cast("Group", zarr.open(cells_store, mode="r"))

    # cells_root["cell_summary"] and cells_root["cell_id"] are zarr Arrays
    cell_summary = np.asarray(cast("Array", cells_root["cell_summary"])[:])
    cell_id = np.asarray(cast("Array", cells_root["cell_id"])[:])
    cell_summary_arr = cast("Array", cells_root["cell_summary"])
    column_names_raw = cast(list[str], cell_summary_arr.attrs.get("column_names", []))
    column_names: list[str] = list(column_names_raw)

    # Create AnnData
    obs_names = [str(cid[0]) for cid in cell_id]

    var = pd.DataFrame({"gene_ids": feature_ids}, index=feature_names)

    obs = pd.DataFrame(index=obs_names)
    for i, col in enumerate(column_names):
        obs[col] = cell_summary[:, i]

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Set spatial coordinates (cell_centroid_x, cell_centroid_y)
    adata.obsm["spatial"] = cell_summary[:, :2]

    return adata


async def load_spatial_data(
    data_path: str,
    data_type: SpatialPlatform,
    name: Optional[str] = None,
) -> dict[str, Any]:
    """Load spatial transcriptomics data.

    Args:
        data_path: Path to the data file or directory
        data_type: Type of spatial data (visium, xenium, slide_seq, merfish, seqfish, generic)
        name: Optional name for the dataset

    Returns:
        Dictionary with dataset information and AnnData object

    Note:
        Async interface for consistency with data management layer, enabling
        future async I/O (aiofiles, remote storage) without API changes.
    """
    # Validate path
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data path not found: {data_path}")

    platform = data_type

    # Import dependencies
    import scanpy as sc

    # Load data based on platform type
    if platform == "visium":
        # For 10x Visium, we need to provide the path to the directory containing the data
        try:
            # Check if it's a directory or an h5ad file
            if os.path.isdir(data_path):
                # Check if the directory has the expected structure
                if os.path.exists(
                    os.path.join(data_path, "filtered_feature_bc_matrix.h5")
                ):
                    # H5 file based 10x Visium directory structure
                    adata = sc.read_visium(data_path)
                elif os.path.exists(
                    os.path.join(data_path, "filtered_feature_bc_matrix")
                ):
                    # Check if it contains MTX files (compressed or uncompressed)
                    mtx_dir = os.path.join(data_path, "filtered_feature_bc_matrix")
                    if os.path.exists(
                        os.path.join(mtx_dir, "matrix.mtx.gz")
                    ) or os.path.exists(os.path.join(mtx_dir, "matrix.mtx")):
                        # Matrix files based 10x Visium directory structure
                        # Use scanpy's read_10x_mtx function
                        adata = sc.read_10x_mtx(
                            os.path.join(data_path, "filtered_feature_bc_matrix"),
                            var_names="gene_symbols",
                            cache=False,
                        )
                        # Try to load spatial coordinates if available
                        spatial_dir = os.path.join(data_path, "spatial")
                        if os.path.exists(spatial_dir):
                            try:
                                # Add spatial information manually
                                import json

                                import pandas as pd

                                # Load tissue positions
                                positions_path = os.path.join(
                                    spatial_dir, "tissue_positions_list.csv"
                                )
                                if os.path.exists(positions_path):
                                    # Try to detect if file has header
                                    with open(positions_path, "r") as f:
                                        first_line = f.readline().strip()

                                    if first_line.startswith("barcode"):
                                        # File has header
                                        positions = pd.read_csv(positions_path)
                                    else:
                                        # File has no header
                                        positions = pd.read_csv(
                                            positions_path, header=None
                                        )
                                        positions.columns = [
                                            "barcode",
                                            "in_tissue",
                                            "array_row",
                                            "array_col",
                                            "pxl_row_in_fullres",
                                            "pxl_col_in_fullres",
                                        ]

                                    positions.set_index("barcode", inplace=True)

                                    # Filter for spots in tissue
                                    positions = positions[positions["in_tissue"] == 1]

                                    # Add spatial coordinates to adata
                                    adata.obsm["spatial"] = positions.loc[
                                        adata.obs_names,
                                        ["pxl_col_in_fullres", "pxl_row_in_fullres"],
                                    ].values

                                    # Load scalefactors
                                    scalefactors_path = os.path.join(
                                        spatial_dir, "scalefactors_json.json"
                                    )
                                    if os.path.exists(scalefactors_path):
                                        with open(scalefactors_path, "r") as f:
                                            scalefactors = json.load(f)

                                        # Add scalefactors to adata
                                        adata.uns["spatial"] = {
                                            "scalefactors": scalefactors
                                        }
                            except Exception as e:
                                logger.warning(
                                    f"Could not load spatial information: {e}"
                                )
                else:
                    raise DataCompatibilityError(
                        f"Directory {data_path} does not have the expected 10x Visium structure"
                    )
            elif os.path.isfile(data_path) and data_path.endswith(".h5"):
                # Single H5 file - new support for 10x H5 format
                adata = sc.read_10x_h5(data_path)

                # Try to find and add spatial information
                spatial_path = _find_spatial_folder(data_path)
                if spatial_path:
                    try:
                        adata = _add_spatial_info_to_adata(adata, spatial_path)
                    except Exception as e:
                        logger.warning(f"Could not add spatial information: {e}")
            elif os.path.isfile(data_path) and data_path.endswith(".h5ad"):
                # If it's an h5ad file but marked as visium, read it as h5ad
                adata = sc.read_h5ad(data_path)
                # Check if it has the necessary spatial information
                if "spatial" not in adata.uns and not any(
                    "spatial" in key for key in adata.obsm.keys()
                ):
                    logger.warning(
                        "The h5ad file does not contain spatial information typically required for Visium data"
                    )
            else:
                raise ParameterError(
                    f"Unsupported file format for visium: {data_path}. Supported formats: directory with Visium structure, .h5 file, or .h5ad file"
                )

        except FileNotFoundError as e:
            raise DataNotFoundError(f"File not found: {e}") from e
        except Exception as e:
            # Provide more detailed error information
            error_msg = f"Error loading Visium data from {data_path}: {e}"

            # Add helpful suggestions based on error type
            if "No matching barcodes" in str(e):
                error_msg += "\n\nPossible solutions:"
                error_msg += "\n1. Check if the H5 file and spatial coordinates are from the same sample"
                error_msg += "\n2. Verify barcode format (with or without -1 suffix)"
                error_msg += "\n3. Ensure the spatial folder contains the correct tissue_positions_list.csv file"
            elif ".h5" in data_path and "read_10x_h5" in str(e):
                error_msg += "\n\nThis might not be a valid 10x H5 file. Try:"
                error_msg += (
                    "\n1. Set data_type='generic' if this is an AnnData H5AD file"
                )
                error_msg += (
                    "\n2. Verify the file is from 10x Genomics Cell Ranger output"
                )
            elif "spatial" in str(e).lower():
                error_msg += "\n\nSpatial data issue detected. Try:"
                error_msg += (
                    "\n1. Loading without spatial data by using data_type='generic'"
                )
                error_msg += "\n2. Ensuring spatial folder contains: tissue_positions_list.csv and scalefactors_json.json"

            raise ProcessingError(error_msg) from e
    elif platform == "xenium":
        # For 10x Xenium data - supports two formats:
        # 1. Zarr format: cell_feature_matrix.zarr.zip + cells.zarr.zip
        # 2. Standard format: cell_feature_matrix.h5 + cells.parquet/csv.gz
        try:
            import pandas as pd

            # Check for zarr format
            matrix_zarr = os.path.join(data_path, "cell_feature_matrix.zarr.zip")
            cells_zarr = os.path.join(data_path, "cells.zarr.zip")
            has_zarr = os.path.exists(matrix_zarr) and os.path.exists(cells_zarr)

            # Check for standard format
            cell_matrix_h5 = os.path.join(data_path, "cell_feature_matrix.h5")
            cell_matrix_dir = os.path.join(data_path, "cell_feature_matrix")
            cells_parquet = os.path.join(data_path, "cells.parquet")
            cells_csv = os.path.join(data_path, "cells.csv.gz")
            has_standard = (
                os.path.exists(cell_matrix_h5) or os.path.exists(cell_matrix_dir)
            ) and (os.path.exists(cells_parquet) or os.path.exists(cells_csv))

            if has_zarr:
                # Load zarr format
                logger.info("Loading Xenium data from zarr format")
                adata = _load_xenium_zarr(data_path)

            elif has_standard:
                # Load standard format
                logger.info("Loading Xenium data from standard format")
                if os.path.exists(cell_matrix_h5):
                    adata = sc.read_10x_h5(cell_matrix_h5)
                else:
                    adata = sc.read_10x_mtx(cell_matrix_dir, var_names="gene_symbols")

                # Load cell metadata
                if os.path.exists(cells_parquet):
                    cells = pd.read_parquet(cells_parquet)
                else:
                    cells = pd.read_csv(cells_csv, compression="gzip")

                # Set cell_id as index for alignment
                if "cell_id" in cells.columns:
                    cells = cells.set_index("cell_id")

                # Align cells metadata with adata
                common_cells = adata.obs_names.intersection(cells.index)
                if len(common_cells) == 0:
                    cells.index = cells.index.astype(str)
                    common_cells = adata.obs_names.intersection(cells.index)

                if len(common_cells) == 0:
                    raise DataCompatibilityError(
                        "No matching cell IDs between count matrix and cells metadata. "
                        f"Count matrix has {adata.n_obs} cells, "
                        f"cells metadata has {len(cells)} cells."
                    )

                # Filter to common cells
                if len(common_cells) < adata.n_obs:
                    logger.info(
                        f"Filtering to {len(common_cells)} cells with spatial coordinates "
                        f"(from {adata.n_obs} total)"
                    )
                    adata = adata[common_cells, :].copy()
                    cells = cells.loc[common_cells]

                # Add spatial coordinates
                if "x_centroid" in cells.columns and "y_centroid" in cells.columns:
                    adata.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].values
                else:
                    raise DataCompatibilityError(
                        "Xenium cells metadata missing x_centroid/y_centroid columns. "
                        f"Available columns: {list(cells.columns)}"
                    )

                # Add other useful cell metadata
                metadata_cols = [
                    "transcript_counts",
                    "control_probe_counts",
                    "control_codeword_counts",
                    "cell_area",
                    "nucleus_area",
                ]
                for col in metadata_cols:
                    if col in cells.columns:
                        adata.obs[col] = cells[col].values

            else:
                raise DataNotFoundError(
                    f"No valid Xenium data found in {data_path}. "
                    "Expected either zarr format (cell_feature_matrix.zarr.zip + cells.zarr.zip) "
                    "or standard format (cell_feature_matrix.h5 + cells.parquet/csv.gz)"
                )

        except (DataNotFoundError, DataCompatibilityError):
            raise
        except Exception as e:
            raise ProcessingError(
                f"Error loading Xenium data from {data_path}: {e}"
            ) from e

    elif platform in ("slide_seq", "merfish", "seqfish", "generic"):
        # For h5ad files or other spatial platforms
        try:
            adata = sc.read_h5ad(data_path)
        except Exception as e:
            raise ProcessingError(f"Error loading {platform} data: {e}") from e
    else:
        raise ParameterError(f"Unsupported platform type: {platform}")

    # Set dataset name
    dataset_name = name or os.path.basename(data_path).split(".")[0]

    # Calculate basic statistics
    n_cells = adata.n_obs
    n_genes = adata.n_vars

    # Check if spatial coordinates are available
    # Priority: obsm["spatial"] is the actual coordinate storage location
    # uns["spatial"] only contains metadata (scalefactors, images) not coordinates
    spatial_coordinates_available = (
        hasattr(adata, "obsm")
        and "spatial" in adata.obsm
        and adata.obsm["spatial"] is not None
        and len(adata.obsm["spatial"]) > 0
    )

    # Check if tissue image is available (for Visium data)
    # Structure: adata.uns["spatial"][library_id]["images"]["hires"/"lowres"]
    # Must check for actual hires or lowres images, not just non-empty dict
    tissue_image_available = False
    if "spatial" in adata.uns and isinstance(adata.uns["spatial"], dict):
        for _sample_key, sample_data in adata.uns["spatial"].items():
            # Each sample_data should be a dict with "images" key
            if isinstance(sample_data, dict) and "images" in sample_data:
                images_dict = sample_data["images"]
                # Check if images dict has actual hires or lowres images
                if isinstance(images_dict, dict) and (
                    "hires" in images_dict or "lowres" in images_dict
                ):
                    tissue_image_available = True
                    break

    # Make variable names unique to avoid reindexing issues
    ensure_unique_var_names(adata)

    # Preserve raw data for downstream analysis (C2 strategy)
    # Only save if .raw doesn't already exist - respect user's existing .raw
    import anndata as ad

    if adata.raw is None:
        # Save current data state to .raw
        # This ensures downstream tools always have access to original loaded data
        # Note: Raw only stores X, var, varm - obs is NOT stored in raw
        adata.raw = ad.AnnData(X=adata.X.copy(), var=adata.var)

    # Also ensure layers["counts"] exists for scVI-tools compatibility
    # Reference raw.X instead of creating another copy - safe because:
    # 1. adata.copy() creates deep copies of layers (won't affect original raw.X)
    # 2. Preprocessing modifies adata.X, not layers["counts"]
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.raw.X

    # Get metadata profile for LLM understanding
    profile = get_adata_profile(adata)

    # Return dataset info and AnnData object with comprehensive metadata
    return {
        "name": dataset_name,
        "type": platform,  # Always a valid SpatialPlatform value
        "path": data_path,
        "adata": adata,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "spatial_coordinates_available": spatial_coordinates_available,
        "tissue_image_available": tissue_image_available,
        # Metadata profile from adata_utils
        **profile,
    }


def _find_spatial_folder(h5_path: str) -> Optional[str]:
    """
    Intelligently find spatial information folder for a given H5 file.

    Search strategy:
    1. Same directory 'spatial' folder
    2. Parent directory 'spatial' folder
    3. Same name prefix spatial folder
    4. Common variations

    Args:
        h5_path: Path to the H5 file

    Returns:
        Path to spatial folder if found, None otherwise
    """
    base_dir = os.path.dirname(h5_path)
    base_name = os.path.splitext(os.path.basename(h5_path))[0]

    # Candidate paths to check
    candidates = [
        os.path.join(base_dir, "spatial"),
        os.path.join(base_dir, "..", "spatial"),
        os.path.join(base_dir, f"{base_name}_spatial"),
        os.path.join(base_dir, "spatial_data"),
        # Check for sample-specific spatial folders
        os.path.join(
            base_dir, base_name.replace("_filtered_feature_bc_matrix", "_spatial")
        ),
        os.path.join(base_dir, base_name.replace("_matrix", "_spatial")),
    ]

    for candidate in candidates:
        candidate = os.path.normpath(candidate)
        if os.path.exists(candidate) and os.path.isdir(candidate):
            # Verify it contains required spatial files
            required_files = ["tissue_positions_list.csv", "scalefactors_json.json"]
            if all(os.path.exists(os.path.join(candidate, f)) for f in required_files):
                return candidate

    logger.warning(f"No spatial folder found for {h5_path}")
    return None


def _add_spatial_info_to_adata(adata: Any, spatial_path: str) -> Any:
    """
    Add spatial information to an AnnData object.

    Args:
        adata: AnnData object with expression data
        spatial_path: Path to spatial information folder

    Returns:
        AnnData object with spatial information added
    """
    import json

    import numpy as np
    import pandas as pd

    try:
        # Load tissue positions
        positions_file = os.path.join(spatial_path, "tissue_positions_list.csv")

        # Try to detect if file has header
        with open(positions_file, "r") as f:
            first_line = f.readline().strip()

        if first_line.startswith("barcode"):
            # File has header
            positions = pd.read_csv(positions_file)
        else:
            # File has no header
            positions = pd.read_csv(positions_file, header=None)

            # Handle different formats of tissue positions file
            if len(positions.columns) == 6:
                positions.columns = [
                    "barcode",
                    "in_tissue",
                    "array_row",
                    "array_col",
                    "pxl_row_in_fullres",
                    "pxl_col_in_fullres",
                ]
            elif len(positions.columns) == 5:
                # Some datasets don't have the 'in_tissue' column
                positions.columns = [
                    "barcode",
                    "array_row",
                    "array_col",
                    "pxl_row_in_fullres",
                    "pxl_col_in_fullres",
                ]
                positions["in_tissue"] = 1  # Assume all spots are in tissue
            else:
                raise DataCompatibilityError(
                    f"Unexpected tissue positions format with {len(positions.columns)} columns"
                )

        positions.set_index("barcode", inplace=True)

        # Find common barcodes between expression data and spatial coordinates
        common_barcodes = adata.obs_names.intersection(positions.index)

        if len(common_barcodes) == 0:
            # Try with modified barcode format (sometimes -1 suffix is added/removed)
            if all("-1" in bc for bc in adata.obs_names[:10]):
                # Expression data has -1 suffix, spatial doesn't
                positions.index = positions.index + "-1"
            elif all("-1" not in bc for bc in adata.obs_names[:10]) and all(
                "-1" in bc for bc in positions.index[:10]
            ):
                # Spatial has -1 suffix, expression doesn't
                positions.index = positions.index.str.replace("-1", "")

            # Try again
            common_barcodes = adata.obs_names.intersection(positions.index)

        if len(common_barcodes) == 0:
            raise DataCompatibilityError(
                "No matching barcodes between expression data and spatial coordinates"
            )

        # Filter to common barcodes
        adata = adata[common_barcodes, :].copy()
        positions = positions.loc[common_barcodes]

        # Add spatial coordinates
        adata.obsm["spatial"] = positions[
            ["pxl_col_in_fullres", "pxl_row_in_fullres"]
        ].values.astype(float)

        # Add tissue information
        if "in_tissue" in positions.columns:
            adata.obs["in_tissue"] = positions["in_tissue"].values

        # Load scalefactors
        scalefactors_file = os.path.join(spatial_path, "scalefactors_json.json")
        with open(scalefactors_file, "r") as f:
            scalefactors = json.load(f)

        # Generate meaningful library_id from spatial_path
        # Priority: parent directory name (usually sample name) > "sample_1" default
        # Avoid using "spatial" as library_id to prevent confusing adata.uns["spatial"]["spatial"] nesting
        parent_dir = os.path.dirname(spatial_path.rstrip(os.sep))
        if parent_dir and os.path.basename(parent_dir) != "":
            library_id = os.path.basename(parent_dir)
        else:
            library_id = "sample_1"  # Fallback to clear default name

        # Create spatial uns structure (scanpy expects nested structure)
        adata.uns["spatial"] = {
            library_id: {"scalefactors": scalefactors, "images": {}}
        }

        # Try to load images if available (using centralized dependency manager)
        if is_available("Pillow"):
            from PIL import Image

            for img_name in ["tissue_hires_image.png", "tissue_lowres_image.png"]:
                img_path = os.path.join(spatial_path, img_name)
                if os.path.exists(img_path):
                    try:
                        img = np.array(Image.open(img_path))

                        img_key = "hires" if "hires" in img_name else "lowres"
                        adata.uns["spatial"][library_id]["images"][img_key] = img
                    except Exception as e:
                        logger.warning(f"Could not load image {img_name}: {e}")
        else:
            logger.warning("Pillow not available, skipping tissue image loading")

        return adata

    except Exception as e:
        logger.error(f"Failed to add spatial information: {e}")
        raise
