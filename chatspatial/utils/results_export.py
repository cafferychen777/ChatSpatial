"""
Analysis results export utilities for reproducibility.

Provides automatic CSV export of analysis results:
~/.chatspatial/results/{data_id}/{analysis_type}/{method}_{key}.csv

Design Principles:
- Convention over configuration: Fixed path structure
- Metadata-driven: Uses results_keys from store_analysis_metadata
- Reproducibility: Every analysis generates traceable output files
"""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

if TYPE_CHECKING:
    from anndata import AnnData

logger = logging.getLogger(__name__)


# =============================================================================
# Directory Management
# =============================================================================


def get_results_dir(data_id: str) -> Path:
    """
    Get the results directory for a dataset.

    Args:
        data_id: Dataset identifier

    Returns:
        Path to ~/.chatspatial/results/{data_id}/
    """
    results_dir = Path.home() / ".chatspatial" / "results" / data_id
    results_dir.mkdir(parents=True, exist_ok=True)
    return results_dir


def get_analysis_dir(data_id: str, analysis_name: str) -> Path:
    """
    Get the directory for a specific analysis type.

    Args:
        data_id: Dataset identifier
        analysis_name: Analysis name (e.g., "differential_expression")

    Returns:
        Path to ~/.chatspatial/results/{data_id}/{analysis_name}/
    """
    analysis_dir = get_results_dir(data_id) / analysis_name
    analysis_dir.mkdir(parents=True, exist_ok=True)
    return analysis_dir


# =============================================================================
# Core Export Function
# =============================================================================


def export_analysis_result(
    adata: "AnnData",
    data_id: str,
    analysis_name: str,
) -> list[Path]:
    """
    Export analysis results to CSV based on stored metadata.

    Automatically discovers result locations from {analysis_name}_metadata
    and exports to appropriate CSV files.

    Args:
        adata: AnnData object containing analysis results
        data_id: Dataset identifier
        analysis_name: Analysis name (must match store_analysis_metadata)

    Returns:
        List of exported file paths

    Example:
        >>> export_analysis_result(adata, "my_data", "differential_expression")
        [Path('~/.chatspatial/results/my_data/differential_expression/wilcoxon_rank_genes_groups.csv')]
    """
    metadata_key = f"{analysis_name}_metadata"
    if metadata_key not in adata.uns:
        logger.debug(f"No metadata found for {analysis_name}, skipping export")
        return []

    metadata = adata.uns[metadata_key]
    results_keys = metadata.get("results_keys", {})
    method = metadata.get("method", "unknown")

    if not results_keys:
        logger.debug(f"No results_keys in metadata for {analysis_name}")
        return []

    analysis_dir = get_analysis_dir(data_id, analysis_name)
    exported_files: list[Path] = []

    # Export based on storage location
    for location, keys in results_keys.items():
        for key in keys:
            try:
                df = _extract_as_dataframe(adata, location, key, analysis_name)
                if df is not None and len(df) > 0:
                    # Sanitize key for filename
                    safe_key = key.replace("/", "_").replace("\\", "_")
                    filename = f"{method}_{safe_key}.csv"
                    filepath = analysis_dir / filename
                    df.to_csv(filepath, index=True)
                    exported_files.append(filepath)
                    logger.info(f"Exported {analysis_name} result to {filepath}")
            except Exception as e:
                logger.warning(f"Failed to export {key} from {location}: {e}")

    # Update index
    if exported_files:
        _update_index(data_id, analysis_name, method, metadata, exported_files)

    return exported_files


# =============================================================================
# Data Extraction
# =============================================================================


def _extract_as_dataframe(
    adata: "AnnData",
    location: str,
    key: str,
    analysis_name: str,
) -> pd.DataFrame | None:
    """
    Extract data from adata and convert to DataFrame.

    Handles different storage locations and data formats:
    - uns: dict, DataFrame, structured array
    - obs: categorical/numeric columns
    - var: gene-level statistics
    - obsm: cell-level matrices
    """
    # Special case: scanpy rank_genes_groups
    if key == "rank_genes_groups" and location == "uns":
        return _extract_rank_genes_groups(adata)

    # uns location
    if location == "uns":
        return _extract_from_uns(adata, key)

    # obs location
    if location == "obs":
        return _extract_from_obs(adata, key)

    # var location
    if location == "var":
        return _extract_from_var(adata, key)

    # obsm location
    if location == "obsm":
        return _extract_from_obsm(adata, key)

    logger.warning(f"Unknown location: {location}")
    return None


def _extract_rank_genes_groups(adata: "AnnData") -> pd.DataFrame | None:
    """Extract differential expression results using scanpy's built-in function."""
    try:
        import scanpy as sc

        if "rank_genes_groups" not in adata.uns:
            return None

        # Get results for all groups
        df = sc.get.rank_genes_groups_df(adata, group=None)
        return df
    except Exception as e:
        logger.warning(f"Failed to extract rank_genes_groups: {e}")
        return None


def _extract_from_uns(adata: "AnnData", key: str) -> pd.DataFrame | None:
    """Extract data from adata.uns."""
    if key not in adata.uns:
        return None

    data = adata.uns[key]

    # Already a DataFrame
    if isinstance(data, pd.DataFrame):
        return data.copy()

    # Special case: squidpy nhood_enrichment or co_occurrence
    # Format: {"zscore": np.ndarray, "count": np.ndarray}
    if key.endswith("_nhood_enrichment") or key.endswith("_co_occurrence"):
        return _extract_squidpy_spatial_result(adata, key, data)

    # Special case: squidpy ripley
    # Format: {"L_stat": DataFrame, "sims_stat": DataFrame, "bins": array, "pvalues": array}
    if "_ripley_" in key and isinstance(data, dict) and "L_stat" in data:
        # Export the main L_stat DataFrame
        return (
            data["L_stat"].copy() if isinstance(data["L_stat"], pd.DataFrame) else None
        )

    # Dictionary - convert to DataFrame
    if isinstance(data, dict):
        return _dict_to_dataframe(data)

    # Structured array (numpy recarray)
    if hasattr(data, "dtype") and data.dtype.names is not None:
        return pd.DataFrame({name: data[name] for name in data.dtype.names})

    logger.debug(f"Unsupported data type in uns[{key}]: {type(data)}")
    return None


def _extract_squidpy_spatial_result(
    adata: "AnnData", key: str, data: dict[str, Any]
) -> pd.DataFrame | None:
    """
    Extract squidpy spatial analysis results (nhood_enrichment, co_occurrence).

    Formats:
    - nhood_enrichment: {"zscore": 2D array, "count": 2D array}
    - co_occurrence: {"occ": 3D array (clusters, clusters, intervals), "interval": 1D array}

    Returns a DataFrame with cluster labels as index.
    """
    import numpy as np

    if not isinstance(data, dict):
        return None

    # Extract cluster_key from key name
    if key.endswith("_nhood_enrichment"):
        cluster_key = key[: -len("_nhood_enrichment")]
    elif key.endswith("_co_occurrence"):
        cluster_key = key[: -len("_co_occurrence")]
    else:
        return None

    # Get cluster labels from adata.obs
    if cluster_key not in adata.obs.columns:
        logger.warning(f"Cluster key '{cluster_key}' not found in adata.obs")
        return None

    cluster_col = adata.obs[cluster_key]
    if hasattr(cluster_col, "cat"):
        labels = list(cluster_col.cat.categories)
    else:
        labels = list(cluster_col.unique())

    # Handle co_occurrence special case (3D array with intervals)
    if "occ" in data and isinstance(data["occ"], np.ndarray) and data["occ"].ndim == 3:
        return _extract_co_occurrence(data, labels, cluster_key)

    # Handle nhood_enrichment and similar 2D matrix results
    result_dfs = []

    for metric_name, metric_data in data.items():
        if isinstance(metric_data, np.ndarray) and metric_data.ndim == 2:
            # Create DataFrame with cluster labels
            n_rows, n_cols = metric_data.shape

            # Ensure labels match dimensions
            row_labels = labels[:n_rows] if len(labels) >= n_rows else labels
            col_labels = labels[:n_cols] if len(labels) >= n_cols else labels

            df = pd.DataFrame(
                metric_data,
                index=row_labels,
                columns=[f"{metric_name}_{c}" for c in col_labels],
            )
            df.index.name = cluster_key
            result_dfs.append(df)

    if not result_dfs:
        return None

    # Concatenate all metrics horizontally
    return pd.concat(result_dfs, axis=1)


def _extract_co_occurrence(
    data: dict[str, Any], labels: list, cluster_key: str
) -> pd.DataFrame | None:
    """
    Extract co_occurrence results with 3D array format.

    Format: {"occ": (n_clusters, n_clusters, n_intervals), "interval": (n_intervals+1,)}

    Exports the mean co-occurrence across all distance intervals as a 2D matrix.
    """
    import numpy as np

    occ = data["occ"]  # Shape: (n_clusters, n_clusters, n_intervals)

    n_clusters_row, n_clusters_col, n_intervals = occ.shape

    # Ensure labels match dimensions
    row_labels = labels[:n_clusters_row] if len(labels) >= n_clusters_row else labels
    col_labels = labels[:n_clusters_col] if len(labels) >= n_clusters_col else labels

    # Export mean co-occurrence across intervals
    mean_occ = np.mean(occ, axis=2)
    df_mean = pd.DataFrame(
        mean_occ,
        index=row_labels,
        columns=[f"occ_mean_{c}" for c in col_labels],
    )

    # Also export first interval (nearest neighbor) for comparison
    first_occ = occ[:, :, 0]
    df_first = pd.DataFrame(
        first_occ,
        index=row_labels,
        columns=[f"occ_nearest_{c}" for c in col_labels],
    )

    # Combine both metrics
    result = pd.concat([df_mean, df_first], axis=1)
    result.index.name = cluster_key

    return result


def _extract_from_obs(adata: "AnnData", key: str) -> pd.DataFrame | None:
    """Extract columns from adata.obs matching the key pattern."""
    # Find columns containing the key
    matching_cols = [c for c in adata.obs.columns if key in c]

    if not matching_cols:
        # Try exact match
        if key in adata.obs.columns:
            matching_cols = [key]
        else:
            return None

    df = adata.obs[matching_cols].copy()
    df.index.name = "cell_id"
    return df


def _extract_from_var(adata: "AnnData", key: str) -> pd.DataFrame | None:
    """Extract columns from adata.var matching the key pattern."""
    # Find columns containing the key
    matching_cols = [c for c in adata.var.columns if key in c]

    if not matching_cols:
        # Try exact match
        if key in adata.var.columns:
            matching_cols = [key]
        else:
            return None

    df = adata.var[matching_cols].copy()
    df.index.name = "gene"
    return df


def _extract_from_obsm(adata: "AnnData", key: str) -> pd.DataFrame | None:
    """Extract matrix from adata.obsm and convert to DataFrame."""
    if key not in adata.obsm:
        return None

    data = adata.obsm[key]

    # Handle different data types
    if isinstance(data, pd.DataFrame):
        df = data.copy()
        df.index = adata.obs_names
        return df

    # Check for CellRank Lineage object
    if hasattr(data, "names") and hasattr(data, "X"):
        # Lineage object
        return pd.DataFrame(
            data.X if hasattr(data, "X") else data,
            index=adata.obs_names,
            columns=data.names,
        )

    # Regular numpy array
    import numpy as np

    if isinstance(data, np.ndarray):
        n_cols = data.shape[1] if len(data.shape) > 1 else 1

        # Try to get column names from related metadata
        col_names = _infer_obsm_columns(adata, key, n_cols)

        return pd.DataFrame(data, index=adata.obs_names, columns=col_names)

    return None


def _infer_obsm_columns(adata: "AnnData", key: str, n_cols: int) -> list[str]:
    """Infer column names for obsm matrices."""
    # Check for cell type proportions
    if "cell_type_proportions" in key:
        # Try to find cell types from deconvolution result
        method = key.replace("cell_type_proportions_", "")
        result_key = f"deconvolution_result_{method}"
        if result_key in adata.uns and "cell_types" in adata.uns[result_key]:
            return list(adata.uns[result_key]["cell_types"])

    # Check for spatial scores with interactions
    if "spatial_scores" in key or "spatial_pvals" in key:
        # LIANA spatial results
        if "liana_spatial_interactions" in adata.uns:
            interactions = adata.uns["liana_spatial_interactions"]
            if len(interactions) == n_cols:
                return list(interactions)

    # Default: numeric indices
    return [f"{key}_{i}" for i in range(n_cols)]


# =============================================================================
# Dictionary Conversion
# =============================================================================


def _dict_to_dataframe(data: dict[str, Any]) -> pd.DataFrame:
    """Convert dictionary to DataFrame with appropriate structure."""
    if not data:
        return pd.DataFrame()

    # Check first value type
    first_val = next(iter(data.values()))

    # Nested dict: {pathway: {score: x, pval: y}}
    if isinstance(first_val, dict):
        return pd.DataFrame.from_dict(data, orient="index")

    # List of values: {gene: [val1, val2, ...]}
    if isinstance(first_val, (list, tuple)):
        return pd.DataFrame.from_dict(data, orient="index")

    # Flat dict: {pathway: score}
    return pd.DataFrame(
        {"key": list(data.keys()), "value": list(data.values())}
    ).set_index("key")


# =============================================================================
# Index Management
# =============================================================================


def _update_index(
    data_id: str,
    analysis_name: str,
    method: str,
    metadata: dict[str, Any],
    exported_files: list[Path],
) -> None:
    """Update the _index.json with export information."""
    index_path = get_results_dir(data_id) / "_index.json"

    # Load existing index or create new
    if index_path.exists():
        try:
            with open(index_path) as f:
                index = json.load(f)
        except json.JSONDecodeError:
            index = {"data_id": data_id, "analyses": {}, "created_at": _now()}
    else:
        index = {"data_id": data_id, "analyses": {}, "created_at": _now()}

    # Update with new analysis
    index["analyses"][analysis_name] = {
        "method": method,
        "parameters": _sanitize_for_json(metadata.get("parameters", {})),
        "statistics": _sanitize_for_json(metadata.get("statistics", {})),
        "exported_at": _now(),
        "files": [f.name for f in exported_files],
    }
    index["updated_at"] = _now()

    # Write index
    with open(index_path, "w") as f:
        json.dump(index, f, indent=2, default=str)


def _sanitize_for_json(obj: Any) -> Any:
    """Sanitize object for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _sanitize_for_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_for_json(v) for v in obj]
    if isinstance(obj, (int, float, str, bool, type(None))):
        return obj
    # Convert other types to string
    return str(obj)


def _now() -> str:
    """Get current timestamp as ISO string."""
    return datetime.now().isoformat()


# =============================================================================
# Utility Functions
# =============================================================================


def list_exported_results(data_id: str) -> dict[str, list[str]]:
    """
    List all exported results for a dataset.

    Args:
        data_id: Dataset identifier

    Returns:
        Dictionary mapping analysis names to file lists
    """
    index_path = get_results_dir(data_id) / "_index.json"
    if not index_path.exists():
        return {}

    with open(index_path) as f:
        index = json.load(f)

    return {
        name: info.get("files", []) for name, info in index.get("analyses", {}).items()
    }


def get_result_path(data_id: str, analysis_name: str, filename: str) -> Path:
    """
    Get the full path to an exported result file.

    Args:
        data_id: Dataset identifier
        analysis_name: Analysis name
        filename: Result filename

    Returns:
        Full path to the file
    """
    return get_analysis_dir(data_id, analysis_name) / filename
