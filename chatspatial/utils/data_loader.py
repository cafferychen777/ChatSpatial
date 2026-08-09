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
from .adata_utils import (
    check_is_integer_counts,
    ensure_unique_var_names,
    get_adata_profile,
    get_spatial_key,
    has_tissue_image,
)
from .async_utils import run_sync
from .dependency_manager import is_available
from .exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)

logger = logging.getLogger(__name__)

_TISSUE_POSITIONS_FILENAMES = (
    "tissue_positions.csv",
    "tissue_positions_list.csv",
)

# Space Ranger metadata that belongs in ``adata.obs``. Pixel coordinates are
# represented canonically by ``adata.obsm["spatial"]`` below, while the array
# coordinates preserve the discrete Visium lattice required by grid-based
# spatial methods.
_TISSUE_POSITION_OBS_COLUMNS = (
    "in_tissue",
    "array_row",
    "array_col",
)


def _find_tissue_positions_file(spatial_path: str) -> Optional[str]:
    """Return the supported Space Ranger tissue positions file, if present."""
    for filename in _TISSUE_POSITIONS_FILENAMES:
        candidate = os.path.join(spatial_path, filename)
        if os.path.isfile(candidate):
            return candidate
    return None


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
    # Use context manager to ensure ZipStore file handles are released
    with ZipStore(matrix_zarr, mode="r") as matrix_store:
        matrix_root = cast("Group", zarr.open(matrix_store, mode="r"))
        cf = cast("Group", matrix_root["cell_features"])

        # Get sparse matrix components (CSC format)
        # All data is fully materialized via [:] so store can close after
        data = np.asarray(cast("Array", cf["data"])[:])
        indices = np.asarray(cast("Array", cf["indices"])[:])
        indptr = np.asarray(cast("Array", cf["indptr"])[:])

        # attrs values are JSON type; Xenium format guarantees these are integers
        n_cells = cast(int, cf.attrs["number_cells"])
        n_features = cast(int, cf.attrs["number_features"])

        # Get feature names from zarr attrs (JSON type -> list[str])
        # Xenium format guarantees these are string lists
        feature_keys = cast(list[str], cf.attrs.get("feature_keys", []))
        feature_ids_raw = cast(list[str], cf.attrs.get("feature_ids", []))
        feature_names: list[str] = list(feature_keys)
        feature_ids: list[str] = list(feature_ids_raw)

    # Create CSC matrix then convert to CSR (cells x genes)
    X_csc = sp.csc_matrix((data, indices, indptr), shape=(n_cells, n_features))
    X = X_csc.tocsr()

    # Load cells.zarr for spatial coordinates
    with ZipStore(cells_zarr, mode="r") as cells_store:
        cells_root = cast("Group", zarr.open(cells_store, mode="r"))

        # All data is fully materialized via [:] so store can close after
        cell_summary = np.asarray(cast("Array", cells_root["cell_summary"])[:])
        cell_id = np.asarray(cast("Array", cells_root["cell_id"])[:])
        cell_summary_arr = cast("Array", cells_root["cell_summary"])
        column_names_raw = cast(
            list[str], cell_summary_arr.attrs.get("column_names", [])
        )
        column_names: list[str] = list(column_names_raw)

    # Create AnnData
    obs_names = [str(cid[0]) for cid in cell_id]

    var = pd.DataFrame({"gene_ids": feature_ids}, index=feature_names)

    obs = pd.DataFrame(index=obs_names)
    for i, col in enumerate(column_names):
        obs[col] = cell_summary[:, i]

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Set spatial coordinates using column names instead of hardcoded indices
    x_col = None
    y_col = None
    for i, name in enumerate(column_names):
        lower = name.lower()
        if x_col is None and ("centroid_x" in lower or lower == "x_centroid"):
            x_col = i
        elif y_col is None and ("centroid_y" in lower or lower == "y_centroid"):
            y_col = i
    if x_col is not None and y_col is not None:
        adata.obsm["spatial"] = np.column_stack(
            [cell_summary[:, x_col], cell_summary[:, y_col]]
        )
    else:
        # Fallback: assume first two columns are x, y
        adata.obsm["spatial"] = cell_summary[:, :2]

    return adata


async def load_spatial_data(
    data_path: str,
    data_type: SpatialPlatform,
    name: Optional[str] = None,
) -> dict[str, Any]:
    """Load spatial data without blocking concurrent MCP requests."""
    return await run_sync(_load_spatial_data_sync, data_path, data_type, name)


def _load_visium_directory(data_path: str, scanpy: Any) -> Any:
    """Load either supported Space Ranger directory layout."""
    h5_matrix = os.path.join(data_path, "filtered_feature_bc_matrix.h5")
    matrix_dir = os.path.join(data_path, "filtered_feature_bc_matrix")
    if os.path.exists(h5_matrix):
        return scanpy.read_visium(data_path)
    if not os.path.isdir(matrix_dir):
        raise DataCompatibilityError(
            f"Directory {data_path} does not have the expected 10x Visium structure"
        )
    if not any(
        os.path.exists(os.path.join(matrix_dir, filename))
        for filename in ("matrix.mtx.gz", "matrix.mtx")
    ):
        raise DataCompatibilityError(
            f"Directory {matrix_dir} is missing matrix.mtx or matrix.mtx.gz"
        )

    adata = scanpy.read_10x_mtx(
        matrix_dir,
        var_names="gene_symbols",
        cache=False,
    )
    spatial_dir = os.path.join(data_path, "spatial")
    if not os.path.isdir(spatial_dir):
        raise DataNotFoundError(
            f"Visium spatial directory not found: {spatial_dir}. "
            "Use data_type='generic' to load expression-only data."
        )
    try:
        return _add_spatial_info_to_adata(adata, spatial_dir)
    except (DataNotFoundError, DataCompatibilityError):
        raise
    except Exception as exc:
        raise DataCompatibilityError(
            f"Visium spatial information failed to load: {exc}. "
            "Spatial coordinates are required for Visium data. "
            "If you only need expression data, use data_type='generic'."
        ) from exc


def _load_visium_h5(data_path: str, scanpy: Any) -> Any:
    """Load a standalone 10x H5 matrix and its adjacent spatial metadata."""
    adata = scanpy.read_10x_h5(data_path)
    spatial_path = _find_spatial_folder(data_path)
    if spatial_path is None:
        raise DataNotFoundError(
            f"Visium spatial metadata not found for {data_path}. "
            "Expected a spatial directory containing a tissue positions "
            "file and scalefactors_json.json. Use data_type='generic' "
            "to load expression-only data."
        )
    try:
        return _add_spatial_info_to_adata(adata, spatial_path)
    except (DataNotFoundError, DataCompatibilityError):
        raise
    except Exception as exc:
        raise DataCompatibilityError(
            f"Visium spatial information failed to load: {exc}. "
            "Spatial coordinates are required for Visium data. "
            "If you only need expression data, use data_type='generic'."
        ) from exc


def _visium_error_message(data_path: str, error: Exception) -> str:
    """Build actionable guidance for unexpected Visium loader failures."""
    message = f"Error loading Visium data from {data_path}: {error}"
    error_text = str(error)
    if "No matching barcodes" in error_text:
        return message + (
            "\n\nPossible solutions:"
            "\n1. Check if the H5 file and spatial coordinates are from the same sample"
            "\n2. Verify barcode format (with or without -1 suffix)"
            "\n3. Ensure the spatial folder contains tissue_positions.csv "
            "or tissue_positions_list.csv"
        )
    if ".h5" in data_path and "read_10x_h5" in error_text:
        return message + (
            "\n\nThis might not be a valid 10x H5 file. Try:"
            "\n1. Set data_type='generic' if this is an AnnData H5AD file"
            "\n2. Verify the file is from 10x Genomics Cell Ranger output"
        )
    if "spatial" in error_text.lower():
        return message + (
            "\n\nSpatial data issue detected. Try:"
            "\n1. Loading without spatial data by using data_type='generic'"
            "\n2. Ensuring the spatial folder contains: tissue_positions.csv "
            "(or tissue_positions_list.csv) and scalefactors_json.json"
        )
    return message


def _load_visium_data(data_path: str, scanpy: Any) -> Any:
    """Load a Visium dataset while preserving domain-specific errors."""
    try:
        if os.path.isdir(data_path):
            return _load_visium_directory(data_path, scanpy)
        if os.path.isfile(data_path) and data_path.endswith(".h5"):
            return _load_visium_h5(data_path, scanpy)
        if os.path.isfile(data_path) and data_path.endswith(".h5ad"):
            adata = scanpy.read_h5ad(data_path)
            if "spatial" not in adata.uns and get_spatial_key(adata) is None:
                logger.warning(
                    "The h5ad file does not contain spatial information typically "
                    "required for Visium data"
                )
            return adata
        raise ParameterError(
            f"Unsupported file format for visium: {data_path}. Supported formats: "
            "directory with Visium structure, .h5 file, or .h5ad file"
        )
    except (DataNotFoundError, DataCompatibilityError, ParameterError):
        raise
    except FileNotFoundError as exc:
        raise DataNotFoundError(f"File not found: {exc}") from exc
    except Exception as exc:
        raise ProcessingError(_visium_error_message(data_path, exc)) from exc


def _align_xenium_metadata(adata: Any, cells: Any) -> Any:
    """Align Xenium cell metadata to count-matrix observation identity."""
    if "cell_id" in cells.columns:
        cells = cells.set_index("cell_id")
    cells.index = cells.index.astype(str)
    if not adata.obs_names.is_unique:
        raise DataCompatibilityError(
            "Xenium count matrix contains duplicate cell IDs; "
            "cell-level metadata cannot be aligned unambiguously."
        )
    if not cells.index.is_unique:
        examples = cells.index[cells.index.duplicated()].unique()[:5]
        raise DataCompatibilityError(
            "Xenium cells metadata contains duplicate cell IDs: " f"{examples.tolist()}"
        )

    common_cells = adata.obs_names[adata.obs_names.isin(cells.index)]
    if len(common_cells) == 0:
        raise DataCompatibilityError(
            "No matching cell IDs between count matrix and cells metadata. "
            f"Count matrix has {adata.n_obs} cells, cells metadata has {len(cells)} cells."
        )
    if len(common_cells) < adata.n_obs:
        logger.info(
            "Filtering to %d cells with spatial coordinates (from %d total)",
            len(common_cells),
            adata.n_obs,
        )
        adata = adata[common_cells, :].copy()

    cells = cells.reindex(adata.obs_names)
    required_coordinates = ["x_centroid", "y_centroid"]
    if not set(required_coordinates).issubset(cells.columns):
        raise DataCompatibilityError(
            "Xenium cells metadata missing x_centroid/y_centroid columns. "
            f"Available columns: {list(cells.columns)}"
        )
    adata.obsm["spatial"] = cells[required_coordinates].to_numpy()
    for column in (
        "transcript_counts",
        "control_probe_counts",
        "control_codeword_counts",
        "cell_area",
        "nucleus_area",
    ):
        if column in cells:
            adata.obs[column] = cells[column]
    return adata


def _load_xenium_data(data_path: str, scanpy: Any) -> Any:
    """Load one supported Xenium directory layout."""
    import pandas as pd

    matrix_zarr = os.path.join(data_path, "cell_feature_matrix.zarr.zip")
    cells_zarr = os.path.join(data_path, "cells.zarr.zip")
    matrix_h5 = os.path.join(data_path, "cell_feature_matrix.h5")
    matrix_dir = os.path.join(data_path, "cell_feature_matrix")
    cells_parquet = os.path.join(data_path, "cells.parquet")
    cells_csv = os.path.join(data_path, "cells.csv.gz")

    try:
        if os.path.exists(matrix_zarr) and os.path.exists(cells_zarr):
            logger.info("Loading Xenium data from zarr format")
            return _load_xenium_zarr(data_path)

        has_matrix = os.path.exists(matrix_h5) or os.path.exists(matrix_dir)
        has_cells = os.path.exists(cells_parquet) or os.path.exists(cells_csv)
        if not has_matrix or not has_cells:
            raise DataNotFoundError(
                f"No valid Xenium data found in {data_path}. Expected either zarr "
                "format (cell_feature_matrix.zarr.zip + cells.zarr.zip) or standard "
                "format (cell_feature_matrix.h5 + cells.parquet/csv.gz)"
            )

        logger.info("Loading Xenium data from standard format")
        adata = (
            scanpy.read_10x_h5(matrix_h5)
            if os.path.exists(matrix_h5)
            else scanpy.read_10x_mtx(matrix_dir, var_names="gene_symbols")
        )
        cells = (
            pd.read_parquet(cells_parquet)
            if os.path.exists(cells_parquet)
            else pd.read_csv(cells_csv, compression="gzip")
        )
        return _align_xenium_metadata(adata, cells)
    except (DataNotFoundError, DataCompatibilityError):
        raise
    except Exception as exc:
        raise ProcessingError(
            f"Error loading Xenium data from {data_path}: {exc}"
        ) from exc


def _ensure_raw_count_sources(adata: Any) -> None:
    """Populate raw/count containers only when integer counts are proven."""
    import anndata as ad

    if adata.raw is None:
        is_integer, has_negative, _ = check_is_integer_counts(adata.X)
        if is_integer and not has_negative:
            adata.raw = ad.AnnData(X=adata.X.copy(), var=adata.var)
    if "counts" in adata.layers:
        return

    is_integer, has_negative, _ = check_is_integer_counts(adata.X)
    if is_integer and not has_negative:
        adata.layers["counts"] = adata.X.copy()
        return
    if adata.raw is None:
        logger.info(
            "Loaded data does not contain proven raw integer counts. "
            "Run preprocess_data() to generate a counts layer."
        )
        return

    raw_adata = adata.raw.to_adata()
    same_genes = len(raw_adata.var_names) == adata.n_vars and set(
        raw_adata.var_names
    ) == set(adata.var_names)
    if not same_genes:
        logger.info(
            "adata.raw has different genes than current adata (%d vs %d). "
            "Skipping layers['counts'] creation.",
            raw_adata.n_vars,
            adata.n_vars,
        )
        return

    is_raw_integer, raw_has_negative, _ = check_is_integer_counts(raw_adata.X)
    if not is_raw_integer or raw_has_negative:
        logger.info(
            "Neither adata.X nor adata.raw contain raw integer counts. "
            "Skipping layers['counts'] creation."
        )
        return
    aligned_raw = (
        raw_adata
        if raw_adata.var_names.equals(adata.var_names)
        else raw_adata[:, adata.var_names]
    )
    adata.layers["counts"] = aligned_raw.X.copy()


def _finalize_loaded_data(
    adata: Any,
    *,
    data_path: str,
    platform: SpatialPlatform,
    name: str | None,
) -> dict[str, Any]:
    """Apply shared invariants and build the loader response."""
    dataset_name = name or os.path.basename(data_path).split(".")[0]
    n_cells, n_genes = adata.n_obs, adata.n_vars
    spatial_coordinates_available = get_spatial_key(adata) is not None
    tissue_image_available = has_tissue_image(adata)
    ensure_unique_var_names(adata)
    _ensure_raw_count_sources(adata)
    return {
        "name": dataset_name,
        "type": platform,
        "path": data_path,
        "adata": adata,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "spatial_coordinates_available": spatial_coordinates_available,
        "tissue_image_available": tissue_image_available,
        **get_adata_profile(adata),
    }


def _load_spatial_data_sync(
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

    This private implementation is synchronous and is called through the public
    async boundary above.
    """
    # Validate path
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data path not found: {data_path}")

    platform = data_type

    # Import dependencies
    import scanpy as sc

    if platform == "visium":
        adata = _load_visium_data(data_path, sc)
    elif platform == "xenium":
        adata = _load_xenium_data(data_path, sc)

    elif platform in ("slide_seq", "merfish", "seqfish", "generic"):
        # For h5ad files or other spatial platforms
        try:
            adata = sc.read_h5ad(data_path)
        except Exception as e:
            raise ProcessingError(f"Error loading {platform} data: {e}") from e
    else:
        raise ParameterError(f"Unsupported platform type: {platform}")

    return _finalize_loaded_data(
        adata,
        data_path=data_path,
        platform=platform,
        name=name,
    )


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
        normalized_candidate = os.path.normpath(candidate)
        if os.path.isdir(normalized_candidate):
            # Verify it contains a supported positions file and scale factors.
            has_positions = (
                _find_tissue_positions_file(normalized_candidate) is not None
            )
            has_scalefactors = os.path.isfile(
                os.path.join(normalized_candidate, "scalefactors_json.json")
            )
            if has_positions and has_scalefactors:
                return normalized_candidate

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
        positions_file = _find_tissue_positions_file(spatial_path)
        if positions_file is None:
            supported = ", ".join(_TISSUE_POSITIONS_FILENAMES)
            raise DataNotFoundError(
                f"No tissue positions file found in {spatial_path}. "
                f"Expected one of: {supported}."
            )

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
                positions.index = positions.index.str.replace(r"-1$", "", regex=True)

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

        # Preserve canonical Space Ranger spot metadata in observation order.
        # Assign positionally after aligning ``positions`` to ``common_barcodes``
        # so pandas cannot silently reindex values against barcode labels.
        for column in _TISSUE_POSITION_OBS_COLUMNS:
            if column in positions.columns:
                adata.obs[column] = positions[column].to_numpy(copy=True)

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
