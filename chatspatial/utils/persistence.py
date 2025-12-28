"""
Data persistence utilities for spatial transcriptomics data.

Handles saving AnnData objects to disk with proper path management.
"""

import os
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def get_save_path(data_id: str, original_path: str) -> Path:
    """
    Get save path for adata, supports environment variable configuration.

    Priority:
    1. CHATSPATIAL_DATA_DIR environment variable
    2. .chatspatial_saved/ directory next to original data (default)

    Args:
        data_id: Dataset identifier
        original_path: Original data file path

    Returns:
        Directory path for saving
    """
    env_dir = os.getenv("CHATSPATIAL_DATA_DIR")
    if env_dir:
        save_dir = Path(env_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        return save_dir

    # Default: use directory next to original data
    path_obj = Path(original_path)

    # Determine parent directory based on whether path looks like a file
    # Check if path has a file extension or ends with a known data format
    if path_obj.suffix in [".h5ad", ".h5", ".csv", ".txt", ".mtx", ".gz"]:
        # It's a file path, use parent directory
        parent_dir = path_obj.parent
    elif path_obj.is_dir():
        # It's an existing directory
        parent_dir = path_obj
    else:
        # Assume it's a file path (even if doesn't exist yet)
        parent_dir = path_obj.parent

    save_dir = parent_dir / ".chatspatial_saved"
    save_dir.mkdir(parents=True, exist_ok=True)
    return save_dir


def save_adata(data_id: str, adata: "AnnData", original_path: str) -> Path:
    """
    Save AnnData object to disk.

    Args:
        data_id: Dataset identifier
        adata: AnnData object to save
        original_path: Original data file path

    Returns:
        Path where data was saved

    Raises:
        IOError: If save fails
    """
    save_dir = get_save_path(data_id, original_path)
    save_path = save_dir / f"{data_id}.h5ad"

    try:
        adata.write_h5ad(save_path, compression="gzip", compression_opts=4)
        return save_path
    except Exception as e:
        raise IOError(f"Failed to save data to {save_path}: {str(e)}") from e
