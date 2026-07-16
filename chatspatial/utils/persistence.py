"""
Data persistence utilities for spatial transcriptomics data.

Provides a fixed active directory for MCP-script data sharing:
~/.chatspatial/active/{data_id}.h5ad

Design Principles:
- Convention over configuration: Fixed path, no environment variables
- Predictable: Users always know where data is
- Symmetric: export_data and reload_data use same paths
"""

from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
from pandas.api.types import is_object_dtype, is_string_dtype


def _coerce_object_series_for_h5ad(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    non_null = series.notna()
    if non_null.any() and numeric[non_null].notna().all():
        return numeric
    return series.where(non_null, "").map(str)


if TYPE_CHECKING:
    from anndata import AnnData


def get_active_dir() -> Path:
    """
    Get the fixed active directory for MCP-script data sharing.

    Returns:
        Path to ~/.chatspatial/active/
    """
    active_dir = Path.home() / ".chatspatial" / "active"
    active_dir.mkdir(parents=True, exist_ok=True)
    return active_dir


def get_active_path(data_id: str) -> Path:
    """
    Get the active file path for a dataset.

    Args:
        data_id: Dataset identifier

    Returns:
        Path to ~/.chatspatial/active/{data_id}.h5ad
    """
    return get_active_dir() / f"{data_id}.h5ad"


def _sanitize_dataframe_for_h5ad(df: pd.DataFrame) -> pd.DataFrame:
    sanitized = df.copy()
    sanitized.index = sanitized.index.map(str)
    sanitized.columns = sanitized.columns.map(str)
    for column in sanitized.columns:
        series = sanitized[column]
        if is_object_dtype(series.dtype) or is_string_dtype(series.dtype):
            sanitized[column] = _coerce_object_series_for_h5ad(series)
    return sanitized


def _sanitize_value_for_h5ad(value: Any) -> Any:
    if isinstance(value, pd.DataFrame):
        return _sanitize_dataframe_for_h5ad(value)
    if isinstance(value, dict):
        return {str(key): _sanitize_value_for_h5ad(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_sanitize_value_for_h5ad(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_sanitize_value_for_h5ad(item) for item in value)
    return value


def _prepare_export_adata(adata: "AnnData") -> "AnnData":
    export_adata_obj = adata.copy()
    for key in list(export_adata_obj.uns.keys()):
        export_adata_obj.uns[key] = _sanitize_value_for_h5ad(export_adata_obj.uns[key])
    return export_adata_obj


def export_adata(data_id: str, adata: "AnnData", path: Path | None = None) -> Path:
    """
    Export AnnData object to disk.

    Args:
        data_id: Dataset identifier
        adata: AnnData object to export
        path: Optional custom path. If None, uses active directory.

    Returns:
        Path where data was exported

    Raises:
        IOError: If export fails
    """
    if path is None:
        export_path = get_active_path(data_id)
    else:
        export_path = Path(path)
        export_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        export_adata_obj = _prepare_export_adata(adata)
        export_adata_obj.write_h5ad(export_path, compression="gzip", compression_opts=4)
        return export_path
    except Exception as e:
        raise IOError(f"Failed to export data to {export_path}: {e}") from e


def load_adata_from_active(data_id: str, path: Path | None = None) -> "AnnData":
    """
    Load AnnData object from active directory or custom path.

    Args:
        data_id: Dataset identifier
        path: Optional custom path. If None, uses active directory.

    Returns:
        Loaded AnnData object

    Raises:
        FileNotFoundError: If file doesn't exist
        IOError: If load fails
    """
    import anndata

    load_path = get_active_path(data_id) if path is None else Path(path)

    if not load_path.exists():
        raise FileNotFoundError(f"Data file not found: {load_path}")

    try:
        return anndata.read_h5ad(load_path)
    except Exception as e:
        raise IOError(f"Failed to load data from {load_path}: {e}") from e
