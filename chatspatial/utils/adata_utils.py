"""
AnnData utilities for ChatSpatial.

This module provides:
1. Standard field name constants
2. Field discovery functions (get_*_key)
3. Data access functions (get_*)
4. Validation functions (validate_*)

One file for all AnnData-related utilities. No duplication.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

from .exceptions import DataError

# =============================================================================
# Constants: Standard Field Names
# =============================================================================
SPATIAL_KEY = "spatial"
CELL_TYPE_KEY = "cell_type"
CLUSTER_KEY = "leiden"
BATCH_KEY = "batch"

# Alternative names for compatibility
ALTERNATIVE_SPATIAL_KEYS: Set[str] = {
    "spatial",
    "X_spatial",
    "coordinates",
    "coords",
    "spatial_coords",
    "positions",
}
ALTERNATIVE_CELL_TYPE_KEYS: Set[str] = {
    "cell_type",
    "celltype",
    "cell_types",
    "annotation",
    "cell_annotation",
    "predicted_celltype",
}
ALTERNATIVE_CLUSTER_KEYS: Set[str] = {
    "leiden",
    "louvain",
    "clusters",
    "cluster",
    "clustering",
    "cluster_labels",
    "spatial_domains",
}
ALTERNATIVE_BATCH_KEYS: Set[str] = {
    "batch",
    "sample",
    "dataset",
    "experiment",
    "replicate",
    "batch_id",
    "sample_id",
}


# =============================================================================
# Field Discovery: Find keys in AnnData
# =============================================================================
def get_spatial_key(adata: "ad.AnnData") -> Optional[str]:
    """Find spatial coordinate key in adata.obsm."""
    for key in ALTERNATIVE_SPATIAL_KEYS:
        if key in adata.obsm:
            return key
    return None


def get_cell_type_key(adata: "ad.AnnData") -> Optional[str]:
    """Find cell type column in adata.obs."""
    for key in ALTERNATIVE_CELL_TYPE_KEYS:
        if key in adata.obs:
            return key
    return None


def get_cluster_key(adata: "ad.AnnData") -> Optional[str]:
    """Find cluster column in adata.obs."""
    for key in ALTERNATIVE_CLUSTER_KEYS:
        if key in adata.obs:
            return key
    return None


def get_batch_key(adata: "ad.AnnData") -> Optional[str]:
    """Find batch/sample column in adata.obs."""
    for key in ALTERNATIVE_BATCH_KEYS:
        if key in adata.obs:
            return key
    return None


# =============================================================================
# Data Access: Get data from AnnData
# =============================================================================
def get_spatial_coordinates(adata: "ad.AnnData") -> Tuple[np.ndarray, np.ndarray]:
    """
    Get spatial coordinates (x, y) from AnnData.

    Checks obsm['spatial'], alternative keys, and obs['x'/'y'].

    Returns:
        Tuple of (x_coords, y_coords)

    Raises:
        DataError: If no spatial coordinates found
    """
    # Check obsm
    for key in ALTERNATIVE_SPATIAL_KEYS:
        if key in adata.obsm:
            coords = adata.obsm[key]
            return coords[:, 0], coords[:, 1]

    # Check obs x/y
    if "x" in adata.obs and "y" in adata.obs:
        x = pd.to_numeric(adata.obs["x"], errors="coerce").values
        y = pd.to_numeric(adata.obs["y"], errors="coerce").values
        if not (np.any(np.isnan(x)) or np.any(np.isnan(y))):
            return x, y

    raise DataError(
        "No spatial coordinates found. Expected in adata.obsm['spatial'] "
        "or adata.obs['x'/'y']"
    )


# =============================================================================
# Validation: Check and validate AnnData
# =============================================================================
def validate_obs_column(
    adata: "ad.AnnData",
    column: str,
    friendly_name: Optional[str] = None,
) -> None:
    """
    Validate that a column exists in adata.obs.

    Raises:
        DataError: If column not found
    """
    if column not in adata.obs.columns:
        name = friendly_name or f"Column '{column}'"
        available = ", ".join(list(adata.obs.columns)[:10])
        suffix = "..." if len(adata.obs.columns) > 10 else ""
        raise DataError(
            f"{name} not found in adata.obs. Available: {available}{suffix}"
        )


def validate_var_column(
    adata: "ad.AnnData",
    column: str,
    friendly_name: Optional[str] = None,
) -> None:
    """
    Validate that a column exists in adata.var.

    Raises:
        DataError: If column not found
    """
    if column not in adata.var.columns:
        name = friendly_name or f"Column '{column}'"
        available = ", ".join(list(adata.var.columns)[:10])
        suffix = "..." if len(adata.var.columns) > 10 else ""
        raise DataError(
            f"{name} not found in adata.var. Available: {available}{suffix}"
        )


def validate_obs_columns(adata: "ad.AnnData", columns: List[str]) -> None:
    """Validate that multiple columns exist in adata.obs."""
    missing = [col for col in columns if col not in adata.obs.columns]
    if missing:
        available = ", ".join(list(adata.obs.columns)[:10])
        suffix = "..." if len(adata.obs.columns) > 10 else ""
        raise DataError(
            f"Columns not found in adata.obs: {', '.join(missing)}. Available: {available}{suffix}"
        )


def validate_var_columns(adata: "ad.AnnData", columns: List[str]) -> None:
    """Validate that multiple columns exist in adata.var."""
    missing = [col for col in columns if col not in adata.var.columns]
    if missing:
        available = ", ".join(list(adata.var.columns)[:10])
        suffix = "..." if len(adata.var.columns) > 10 else ""
        raise DataError(
            f"Columns not found in adata.var: {', '.join(missing)}. Available: {available}{suffix}"
        )


def validate_adata_basics(
    adata: "ad.AnnData",
    min_obs: int = 1,
    min_vars: int = 1,
) -> None:
    """Validate basic AnnData structure."""
    if adata is None:
        raise DataError("AnnData object cannot be None")
    if adata.n_obs < min_obs:
        raise DataError(f"Dataset has {adata.n_obs} observations, need {min_obs}")
    if adata.n_vars < min_vars:
        raise DataError(f"Dataset has {adata.n_vars} variables, need {min_vars}")


def ensure_categorical(adata: "ad.AnnData", column: str) -> None:
    """Ensure a column is categorical dtype, converting if needed."""
    if column not in adata.obs.columns:
        raise DataError(f"Column '{column}' not found in adata.obs")
    if not pd.api.types.is_categorical_dtype(adata.obs[column]):
        adata.obs[column] = adata.obs[column].astype("category")


# =============================================================================
# Ensure: Make sure something exists
# =============================================================================
def ensure_spatial_key(adata: "ad.AnnData") -> str:
    """
    Ensure spatial coordinates exist and return the key.

    Creates from obs['x'/'y'] if needed.

    Returns:
        The spatial key name

    Raises:
        DataError: If no spatial coordinates found
    """
    key = get_spatial_key(adata)
    if key:
        return key

    # Try to create from obs x/y
    if "x" in adata.obs and "y" in adata.obs:
        x = pd.to_numeric(adata.obs["x"], errors="coerce").values
        y = pd.to_numeric(adata.obs["y"], errors="coerce").values
        if not (np.any(np.isnan(x)) or np.any(np.isnan(y))):
            adata.obsm[SPATIAL_KEY] = np.column_stack([x, y]).astype("float64")
            return SPATIAL_KEY

    raise DataError("No spatial coordinates found")


# =============================================================================
# Standardization
# =============================================================================
def standardize_adata(adata: "ad.AnnData", copy: bool = True) -> "ad.AnnData":
    """
    Standardize AnnData to ChatSpatial conventions.

    Does:
    1. Move spatial coordinates to obsm['spatial']
    2. Make gene names unique
    3. Convert known categorical columns to category dtype

    Does NOT:
    - Compute HVGs (use preprocessing)
    - Compute spatial neighbors (computed by analysis tools)
    """
    if copy:
        adata = adata.copy()

    # Standardize spatial coordinates
    _move_spatial_to_standard(adata)

    # Make gene names unique
    ensure_unique_var_names(adata)

    # Ensure categorical columns
    all_categorical_keys = (
        ALTERNATIVE_CELL_TYPE_KEYS | ALTERNATIVE_CLUSTER_KEYS | ALTERNATIVE_BATCH_KEYS
    )
    for key in adata.obs.columns:
        if key in all_categorical_keys:
            if not pd.api.types.is_categorical_dtype(adata.obs[key]):
                adata.obs[key] = adata.obs[key].astype("category")

    return adata


def _move_spatial_to_standard(adata: "ad.AnnData") -> None:
    """Move spatial coordinates to standard obsm['spatial'] location."""
    if SPATIAL_KEY in adata.obsm:
        return

    # Check alternative obsm keys
    for key in ALTERNATIVE_SPATIAL_KEYS:
        if key in adata.obsm and key != SPATIAL_KEY:
            adata.obsm[SPATIAL_KEY] = adata.obsm[key]
            return

    # Check obs x/y
    if "x" in adata.obs and "y" in adata.obs:
        try:
            x = pd.to_numeric(adata.obs["x"], errors="coerce").values
            y = pd.to_numeric(adata.obs["y"], errors="coerce").values
            if not (np.any(np.isnan(x)) or np.any(np.isnan(y))):
                adata.obsm[SPATIAL_KEY] = np.column_stack([x, y]).astype("float64")
        except Exception:
            pass


# =============================================================================
# Advanced Validation: validate_adata with optional checks
# =============================================================================
def validate_adata(
    adata: "ad.AnnData",
    required_keys: dict,
    check_spatial: bool = False,
    check_velocity: bool = False,
    spatial_key: str = "spatial",
) -> None:
    """
    Validate AnnData object has required keys and optional data integrity checks.

    Args:
        adata: AnnData object to validate
        required_keys: Dict of required keys by category (obs, var, obsm, etc.)
        check_spatial: Whether to validate spatial coordinates
        check_velocity: Whether to validate velocity data layers
        spatial_key: Key for spatial coordinates in adata.obsm

    Raises:
        DataError: If required keys are missing or validation fails
    """
    missing = []

    for category, keys in required_keys.items():
        if isinstance(keys, str):
            keys = [keys]

        attr = getattr(adata, category, None)
        if attr is None:
            missing.extend([f"{category}.{k}" for k in keys])
            continue

        for key in keys:
            if hasattr(attr, "columns"):  # DataFrame
                if key not in attr.columns:
                    missing.append(f"{category}.{key}")
            elif hasattr(attr, "keys"):  # Dict-like
                if key not in attr.keys():
                    missing.append(f"{category}.{key}")
            else:
                missing.append(f"{category}.{key}")

    if missing:
        raise DataError(f"Missing required keys: {', '.join(missing)}")

    # Enhanced validation checks
    if check_spatial:
        _validate_spatial_data(adata, spatial_key, missing)

    if check_velocity:
        _validate_velocity_data(adata, missing)

    if missing:
        raise DataError(f"Validation failed: {', '.join(missing)}")


def _validate_spatial_data(
    adata: "ad.AnnData", spatial_key: str, issues: List[str]
) -> None:
    """Internal helper for spatial data validation."""
    if spatial_key not in adata.obsm:
        issues.append(f"Missing '{spatial_key}' coordinates in adata.obsm")
        return

    spatial_coords = adata.obsm[spatial_key]

    if spatial_coords.shape[1] < 2:
        issues.append(
            f"Spatial coordinates should have at least 2 dimensions, "
            f"found {spatial_coords.shape[1]}"
        )

    if np.any(np.isnan(spatial_coords)):
        issues.append("Spatial coordinates contain NaN values")

    if np.std(spatial_coords[:, 0]) == 0 and np.std(spatial_coords[:, 1]) == 0:
        issues.append("All spatial coordinates are identical")


def _validate_velocity_data(adata: "ad.AnnData", issues: List[str]) -> None:
    """Internal helper for velocity data validation."""
    if "spliced" not in adata.layers:
        issues.append("Missing 'spliced' layer required for RNA velocity")
    if "unspliced" not in adata.layers:
        issues.append("Missing 'unspliced' layer required for RNA velocity")

    if "spliced" in adata.layers and "unspliced" in adata.layers:
        for layer_name in ["spliced", "unspliced"]:
            layer_data = adata.layers[layer_name]

            if hasattr(layer_data, "nnz"):  # Sparse matrix
                if layer_data.nnz == 0:
                    issues.append(f"'{layer_name}' layer is empty (all zeros)")
            else:  # Dense matrix
                if np.all(layer_data == 0):
                    issues.append(f"'{layer_name}' layer is empty (all zeros)")

            if hasattr(layer_data, "data"):  # Sparse matrix
                if np.any(np.isnan(layer_data.data)):
                    issues.append(f"'{layer_name}' layer contains NaN values")
            else:  # Dense matrix
                if np.any(np.isnan(layer_data)):
                    issues.append(f"'{layer_name}' layer contains NaN values")


# =============================================================================
# Metadata Storage: Scientific Provenance Tracking
# =============================================================================
def store_analysis_metadata(
    adata: "ad.AnnData",
    analysis_name: str,
    method: str,
    parameters: Dict[str, Any],
    results_keys: Dict[str, List[str]],
    statistics: Optional[Dict[str, Any]] = None,
    species: Optional[str] = None,
    database: Optional[str] = None,
    reference_info: Optional[Dict[str, Any]] = None,
) -> None:
    """Store analysis metadata in adata.uns for scientific provenance tracking.

    This function stores ONLY scientifically important metadata:
    - Method name (required for reproducibility)
    - Parameters (required for reproducibility)
    - Results locations (required for data access)
    - Statistics (required for quality assessment)
    - Species/Database (required for biological interpretation)
    - Reference info (required for reference-based methods)

    Args:
        adata: AnnData object to store metadata in
        analysis_name: Name of the analysis (e.g., "annotation_tangram")
        method: Method name (e.g., "tangram", "liana", "cellrank")
        parameters: Dictionary of analysis parameters
        results_keys: Dictionary mapping storage location to list of keys
            Example: {"obs": ["cell_type_tangram"], "obsm": ["tangram_ct_pred"]}
        statistics: Optional dictionary of quality/summary statistics
        species: Optional species identifier (critical for communication/enrichment)
        database: Optional database/resource name (critical for communication/enrichment)
        reference_info: Optional reference dataset information
    """
    # Build metadata dictionary - only scientifically important information
    metadata = {
        "method": method,
        "parameters": parameters,
        "results_keys": results_keys,
    }

    # Add optional scientific metadata
    if statistics is not None:
        metadata["statistics"] = statistics

    if species is not None:
        metadata["species"] = species

    if database is not None:
        metadata["database"] = database

    if reference_info is not None:
        metadata["reference_info"] = reference_info

    # Store in adata.uns with unique key
    metadata_key = f"{analysis_name}_metadata"
    adata.uns[metadata_key] = metadata


# =============================================================================
# Gene Selection Utilities
# =============================================================================
def get_highly_variable_genes(
    adata: "ad.AnnData",
    max_genes: int = 500,
    fallback_to_variance: bool = True,
) -> List[str]:
    """
    Get highly variable genes from AnnData.

    Priority order:
    1. Use precomputed HVG from adata.var['highly_variable']
    2. If fallback enabled, compute variance and return top variable genes

    Args:
        adata: AnnData object
        max_genes: Maximum number of genes to return
        fallback_to_variance: If True, compute variance when HVG not available

    Returns:
        List of gene names (may be shorter than max_genes if fewer available)
    """
    # Try precomputed HVG first
    if "highly_variable" in adata.var.columns:
        hvg_genes = adata.var_names[adata.var["highly_variable"]].tolist()
        return hvg_genes[:max_genes]

    # Fallback to variance calculation
    if fallback_to_variance:
        from scipy import sparse

        if sparse.issparse(adata.X):
            var_scores = np.array(adata.X.toarray().var(axis=0)).flatten()
        else:
            var_scores = np.array(adata.X.var(axis=0)).flatten()

        top_indices = np.argsort(var_scores)[-max_genes:]
        return adata.var_names[top_indices].tolist()

    return []


# =============================================================================
# Gene Name Utilities
# =============================================================================
def ensure_unique_var_names(
    adata: "ad.AnnData",
    label: str = "data",
) -> int:
    """
    Ensure gene names are unique, fixing duplicates if needed.

    Args:
        adata: AnnData object (modified in-place)
        label: Label for logging (not used in sync version, for API consistency)

    Returns:
        Number of duplicate gene names that were fixed (0 if already unique)
    """
    if adata.var_names.is_unique:
        return 0

    n_duplicates = len(adata.var_names) - len(set(adata.var_names))
    adata.var_names_make_unique()
    return n_duplicates


async def ensure_unique_var_names_with_ctx(
    adata: "ad.AnnData",
    ctx: Any,  # ToolContext, use Any to avoid circular import
    label: str = "data",
) -> int:
    """
    Ensure gene names are unique with user feedback via ctx.

    Args:
        adata: AnnData object (modified in-place)
        ctx: ToolContext for logging warnings to user
        label: Descriptive label for the data (e.g., "reference data", "query data")

    Returns:
        Number of duplicate gene names that were fixed (0 if already unique)
    """
    n_fixed = ensure_unique_var_names(adata, label)
    if n_fixed > 0:
        await ctx.warning(f"Found {n_fixed} duplicate gene names in {label}, fixed")
    return n_fixed
