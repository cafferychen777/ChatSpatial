"""
Data loading utilities for spatial transcriptomics data.

Handles loading various spatial data formats:
- H5AD files (AnnData format)
- H5 files (10x Genomics format)
- MTX directories (10x Visium structure)
- Visium directories with spatial information

For data persistence, see persistence.py.
"""

import logging
import os
from typing import Any, Dict, Literal, Optional

from .adata_utils import ensure_unique_var_names, get_adata_profile
from .dependency_manager import is_available
from .exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)

logger = logging.getLogger(__name__)


async def load_spatial_data(
    data_path: str,
    data_type: Literal[
        "10x_visium", "slide_seq", "merfish", "seqfish", "other", "auto", "h5ad"
    ] = "auto",
    name: Optional[str] = None,
) -> Dict[str, Any]:
    """Load spatial transcriptomics data

    Args:
        data_path: Path to the data file or directory
        data_type: Type of spatial data. If 'auto', will try to determine the type from the file extension or directory structure.
        name: Optional name for the dataset

    Returns:
        Dictionary with dataset information and AnnData object
    """
    # Validate path
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data path not found: {data_path}")

    # Auto-detect data type if set to 'auto'
    if data_type == "auto":
        if os.path.isfile(data_path):
            if data_path.endswith(".h5ad"):
                # It's an h5ad file
                data_type = "h5ad"
            elif data_path.endswith(".h5"):
                # It's likely a 10x H5 file
                data_type = "10x_visium"
                logger.info("Auto-detected as 10x H5 file, using 10x_visium loader")
            else:
                # Default to other for unknown file types
                data_type = "other"
        elif os.path.isdir(data_path):
            # Check if it has the structure of a 10x Visium dataset
            if os.path.exists(
                os.path.join(data_path, "filtered_feature_bc_matrix")
            ) or os.path.exists(
                os.path.join(data_path, "filtered_feature_bc_matrix.h5")
            ):
                data_type = "10x_visium"
            else:
                # Default to other if we can't determine
                data_type = "other"
        else:
            # Default to other for unknown file types
            data_type = "other"

    # Convert h5ad to other for backward compatibility
    if data_type == "h5ad":
        data_type = "other"

    # Import dependencies
    import scanpy as sc
    import squidpy as sq

    # Load data based on data_type
    if data_type == "10x_visium":
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
                                    f"Could not load spatial information: {str(e)}"
                                )
                else:
                    raise DataCompatibilityError(
                        f"Directory {data_path} does not have the expected 10x Visium structure"
                    )
            elif os.path.isfile(data_path) and data_path.endswith(".h5"):
                # Single H5 file - new support for 10x H5 format
                logger.info(f"Loading 10x H5 file: {data_path}")
                adata = sc.read_10x_h5(data_path)

                # Try to find and add spatial information
                spatial_path = _find_spatial_folder(data_path)
                if spatial_path:
                    try:
                        adata = _add_spatial_info_to_adata(adata, spatial_path)
                        logger.info(
                            f"Successfully added spatial information from {spatial_path}"
                        )
                    except Exception as e:
                        logger.warning(f"Could not add spatial information: {str(e)}")
                        logger.info("Proceeding with expression data only")
                else:
                    logger.info(
                        "No spatial folder found. Loading expression data only."
                    )
                    logger.info(
                        "Tip: Place spatial files in a 'spatial' folder in the same directory as the H5 file"
                    )
            elif os.path.isfile(data_path) and data_path.endswith(".h5ad"):
                # If it's an h5ad file but marked as 10x_visium, read it as h5ad
                adata = sc.read_h5ad(data_path)
                # Check if it has the necessary spatial information
                if "spatial" not in adata.uns and not any(
                    "spatial" in key for key in adata.obsm.keys()
                ):
                    logger.warning(
                        "The h5ad file does not contain spatial information typically required for 10x Visium data"
                    )
            else:
                raise ParameterError(
                    f"Unsupported file format for 10x_visium: {data_path}. Supported formats: directory with Visium structure, .h5 file, or .h5ad file"
                )

            # Add spatial neighborhood graph if not already present
            if "spatial_connectivities" not in adata.obsp and "spatial" in adata.obsm:
                try:
                    sq.gr.spatial_neighbors(adata)
                except Exception as e:
                    logger.warning(f"Could not compute spatial neighbors: {str(e)}")
        except FileNotFoundError as e:
            raise DataNotFoundError(f"File not found: {str(e)}") from e
        except Exception as e:
            # Provide more detailed error information
            error_msg = f"Error loading 10x Visium data from {data_path}: {str(e)}"

            # Add helpful suggestions based on error type
            if "No matching barcodes" in str(e):
                error_msg += "\n\nPossible solutions:"
                error_msg += "\n1. Check if the H5 file and spatial coordinates are from the same sample"
                error_msg += "\n2. Verify barcode format (with or without -1 suffix)"
                error_msg += "\n3. Ensure the spatial folder contains the correct tissue_positions_list.csv file"
            elif ".h5" in data_path and "read_10x_h5" in str(e):
                error_msg += "\n\nThis might not be a valid 10x H5 file. Try:"
                error_msg += "\n1. Set data_type='h5ad' if this is an AnnData H5AD file"
                error_msg += (
                    "\n2. Verify the file is from 10x Genomics Cell Ranger output"
                )
            elif "spatial" in str(e).lower():
                error_msg += "\n\nSpatial data issue detected. Try:"
                error_msg += (
                    "\n1. Loading without spatial data by using data_type='other'"
                )
                error_msg += "\n2. Ensuring spatial folder contains: tissue_positions_list.csv and scalefactors_json.json"

            raise ProcessingError(error_msg) from e
    elif data_type == "h5ad" or data_type in [
        "slide_seq",
        "merfish",
        "seqfish",
        "other",
    ]:
        # For h5ad files or other data types
        try:
            adata = sc.read_h5ad(data_path)
        except Exception as e:
            raise ProcessingError(f"Error loading {data_type} data: {str(e)}") from e
    else:
        raise ParameterError(f"Unsupported data type: {data_type}")

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
        adata.raw = ad.AnnData(
            X=adata.X.copy(),
            var=adata.var,
            obs=adata.obs.copy(),
            uns={},
        )
        logger.info("Preserved loaded data state in adata.raw")

    # Also ensure layers["counts"] exists for scVI-tools compatibility
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        logger.info("Created counts layer from current data")

    # Get metadata profile for LLM understanding
    profile = get_adata_profile(adata)

    # Return dataset info and AnnData object with comprehensive metadata
    return {
        "name": dataset_name,
        "type": data_type,
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
                logger.info(f"Found spatial folder at: {candidate}")
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

        logger.info(
            f"Found {len(common_barcodes)} matching barcodes out of {len(adata)} cells"
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

        logger.info(f"Using library_id: {library_id}")

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
                        logger.info(f"Loaded {img_key} tissue image")
                    except Exception as e:
                        logger.warning(f"Could not load image {img_name}: {str(e)}")
        else:
            logger.warning("Pillow not available, skipping tissue image loading")

        return adata

    except Exception as e:
        logger.error(f"Failed to add spatial information: {str(e)}")
        raise
