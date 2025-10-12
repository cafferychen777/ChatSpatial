"""
Manual data persistence utilities for AnnData objects.

Provides simple manual save functionality without automatic triggers.
"""

import os
from pathlib import Path

import anndata as ad
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
    original_path = Path(original_path)

    # Determine parent directory based on whether path looks like a file
    # Check if path has a file extension or ends with a known data format
    if original_path.suffix in [".h5ad", ".h5", ".csv", ".txt", ".mtx", ".gz"]:
        # It's a file path, use parent directory
        parent_dir = original_path.parent
    elif original_path.is_dir():
        # It's an existing directory
        parent_dir = original_path
    else:
        # Assume it's a file path (even if doesn't exist yet)
        parent_dir = original_path.parent

    save_dir = parent_dir / ".chatspatial_saved"
    save_dir.mkdir(parents=True, exist_ok=True)
    return save_dir


def save_adata(data_id: str, adata: AnnData, original_path: str) -> Path:
    """
    Manually save AnnData object to disk.

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


def load_saved_adata(data_id: str, original_path: str) -> AnnData:
    """
    Load previously saved AnnData object.

    Args:
        data_id: Dataset identifier
        original_path: Original data file path

    Returns:
        Loaded AnnData object

    Raises:
        FileNotFoundError: If saved file doesn't exist
    """
    save_dir = get_save_path(data_id, original_path)
    save_path = save_dir / f"{data_id}.h5ad"

    if not save_path.exists():
        raise FileNotFoundError(f"No saved data found at {save_path}")

    return ad.read_h5ad(save_path)


def check_saved_exists(data_id: str, original_path: str) -> bool:
    """
    Check if saved version exists for a dataset.

    Args:
        data_id: Dataset identifier
        original_path: Original data file path

    Returns:
        True if saved version exists
    """
    save_dir = get_save_path(data_id, original_path)
    save_path = save_dir / f"{data_id}.h5ad"
    return save_path.exists()
