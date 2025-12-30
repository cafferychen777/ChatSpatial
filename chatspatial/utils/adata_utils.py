"""
AnnData utilities for ChatSpatial.

This module provides:
1. Standard field name constants
2. Field discovery functions (get_*_key)
3. Data access functions (get_*)
4. Validation functions (validate_*)

One file for all AnnData-related utilities. No duplication.
"""

from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional, Set, Tuple

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

from scipy import sparse

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


def sample_expression_values(
    adata: "ad.AnnData",
    n_samples: int = 1000,
    layer: Optional[str] = None,
) -> np.ndarray:
    """
    Sample expression values from data matrix for validation checks.

    Efficiently samples values from sparse or dense matrices without
    materializing the full matrix. Used for data type detection
    (integer vs float, negative values, etc.).

    Args:
        adata: AnnData object
        n_samples: Maximum number of values to sample (default: 1000)
        layer: Optional layer name. If None, uses adata.X

    Returns:
        1D numpy array of sampled expression values

    Examples:
        # Check for negative values (indicates log-normalized data)
        sample = sample_expression_values(adata)
        if np.any(sample < 0):
            raise ValueError("Log normalization requires non-negative data")

        # Check for non-integer values (indicates normalized data)
        sample = sample_expression_values(adata)
        if np.any((sample % 1) != 0):
            raise ValueError("Method requires raw count data (integers)")
    """
    # Get the data matrix
    X = adata.layers[layer] if layer is not None else adata.X

    # Handle sparse matrices efficiently
    if sparse.issparse(X):
        # For sparse matrices, sample from .data array (non-zero values only)
        # This is efficient as it doesn't require converting to dense
        if hasattr(X, "data") and len(X.data) > 0:
            return X.data[: min(n_samples, len(X.data))]
        else:
            # Fallback for empty or unusual sparse matrix
            return X[:n_samples].toarray().flatten()
    else:
        # For dense matrices, flatten and sample
        return X.flatten()[: min(n_samples, X.size)]


def require_spatial_coords(
    adata: "ad.AnnData",
    spatial_key: Optional[str] = None,
    validate: bool = True,
) -> np.ndarray:
    """
    Get validated spatial coordinates array from AnnData.

    Unlike get_spatial_coordinates() which returns (x, y) tuple, this
    returns the full coordinates array and provides optional validation.

    Args:
        adata: AnnData object
        spatial_key: Optional key in obsm. If None, auto-detects using
            ALTERNATIVE_SPATIAL_KEYS
        validate: If True (default), validates coordinates for:
            - At least 2 dimensions
            - No NaN values
            - Not all identical

    Returns:
        Spatial coordinates as 2D numpy array (n_cells, n_dims)

    Raises:
        DataError: If spatial coordinates not found or validation fails

    Examples:
        # Auto-detect spatial key
        coords = require_spatial_coords(adata)

        # Use specific key without validation
        coords = require_spatial_coords(adata, spatial_key="X_spatial", validate=False)
    """
    # Find spatial key if not specified
    if spatial_key is None:
        spatial_key = get_spatial_key(adata)
        if spatial_key is None:
            # Also check obs for x/y columns
            if "x" in adata.obs and "y" in adata.obs:
                x = pd.to_numeric(adata.obs["x"], errors="coerce").values
                y = pd.to_numeric(adata.obs["y"], errors="coerce").values
                coords = np.column_stack([x, y])
                if validate and np.any(np.isnan(coords)):
                    raise DataError("Spatial coordinates in obs['x'/'y'] contain NaN")
                return coords

            raise DataError(
                "No spatial coordinates found. Expected in adata.obsm['spatial'] "
                "or similar key. Available obsm keys: "
                f"{list(adata.obsm.keys()) if adata.obsm else 'none'}"
            )

    # Check if key exists
    if spatial_key not in adata.obsm:
        raise DataError(
            f"Spatial coordinates '{spatial_key}' not found in adata.obsm. "
            f"Available keys: {list(adata.obsm.keys())}"
        )

    coords = adata.obsm[spatial_key]

    # Validate if requested
    if validate:
        if coords.shape[1] < 2:
            raise DataError(
                f"Spatial coordinates should have at least 2 dimensions, "
                f"found {coords.shape[1]}"
            )
        if np.any(np.isnan(coords)):
            raise DataError("Spatial coordinates contain NaN values")
        if np.any(np.isinf(coords)):
            raise DataError("Spatial coordinates contain infinite values")
        if np.std(coords[:, 0]) == 0 and np.std(coords[:, 1]) == 0:
            raise DataError("All spatial coordinates are identical")

    return coords


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

    # Ensure categorical columns for known key types
    all_categorical_keys = (
        ALTERNATIVE_CELL_TYPE_KEYS | ALTERNATIVE_CLUSTER_KEYS | ALTERNATIVE_BATCH_KEYS
    )
    for key in adata.obs.columns:
        if key in all_categorical_keys:
            ensure_categorical(adata, key)

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
            # Compute variance on sparse matrix without converting to dense
            # Var(X) = E[X^2] - E[X]^2 (memory efficient, ~5x faster)
            mean = np.array(adata.X.mean(axis=0)).flatten()
            mean_sq = np.array(adata.X.power(2).mean(axis=0)).flatten()
            var_scores = mean_sq - mean**2
        else:
            var_scores = np.array(adata.X.var(axis=0)).flatten()

        top_indices = np.argsort(var_scores)[-max_genes:]
        return adata.var_names[top_indices].tolist()

    return []


def select_genes_for_analysis(
    adata: "ad.AnnData",
    genes: Optional[List[str]] = None,
    n_genes: int = 20,
    require_hvg: bool = True,
    analysis_name: str = "analysis",
) -> List[str]:
    """
    Select genes for spatial/statistical analysis.

    Unified gene selection logic for all analysis tools. Replaces duplicated
    code across spatial_statistics.py and other tools.

    Priority:
        1. User-specified genes (filtered to existing genes)
        2. Highly variable genes (HVG) from preprocessing

    Args:
        adata: AnnData object
        genes: User-specified gene list. If provided, filters to genes in adata.
        n_genes: Maximum number of genes to return when using HVG.
        require_hvg: If True (default), raise error when HVG not found.
                    If False, return empty list when HVG not found.
        analysis_name: Name of analysis for error messages (e.g., "Moran's I").

    Returns:
        List of gene names to analyze.

    Raises:
        DataError: If genes specified but none found, or HVG required but missing.

    Examples:
        # Use user-specified genes
        genes = select_genes_for_analysis(adata, genes=["CD4", "CD8A"])

        # Use top 50 HVGs
        genes = select_genes_for_analysis(adata, n_genes=50)

        # For analysis that can work without HVG
        genes = select_genes_for_analysis(adata, require_hvg=False)
    """
    # Case 1: User specified genes
    if genes is not None:
        valid_genes = [g for g in genes if g in adata.var_names]
        if not valid_genes:
            # Find closest matches for better error message
            from difflib import get_close_matches

            suggestions = []
            for g in genes[:3]:  # Check first 3 genes
                matches = get_close_matches(
                    g, adata.var_names.tolist(), n=1, cutoff=0.6
                )
                if matches:
                    suggestions.append(f"'{g}' → '{matches[0]}'?")

            suggestion_str = (
                f" Did you mean: {', '.join(suggestions)}" if suggestions else ""
            )
            raise DataError(
                f"None of the specified genes found in data: {genes[:5]}..."
                f"{suggestion_str}"
            )
        return valid_genes

    # Case 2: Use HVG
    if "highly_variable" in adata.var.columns and adata.var["highly_variable"].any():
        hvg_genes = adata.var_names[adata.var["highly_variable"]].tolist()
        return hvg_genes[:n_genes]

    # Case 3: HVG not available
    if require_hvg:
        raise DataError(
            f"Highly variable genes (HVG) required for {analysis_name}.\n\n"
            "Solutions:\n"
            "1. Run preprocess_data() first to compute HVGs\n"
            "2. Specify genes explicitly via 'genes' parameter"
        )

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


# =============================================================================
# Raw Counts Data Access: Unified interface for accessing raw data
# =============================================================================


def check_is_integer_counts(X: Any, sample_size: int = 100) -> Tuple[bool, bool, bool]:
    """Check if a matrix contains integer counts.

    This is a lightweight utility for checking data format without
    going through the full data source detection logic.

    Args:
        X: Data matrix (sparse or dense)
        sample_size: Number of rows/cols to sample for efficiency

    Returns:
        Tuple of (is_integer, has_negatives, has_decimals)
    """
    n_rows = min(sample_size, X.shape[0])
    n_cols = min(sample_size, X.shape[1])
    sample = X[:n_rows, :n_cols]

    if sparse.issparse(sample):
        sample = sample.toarray()

    has_negatives = float(sample.min()) < 0
    has_decimals = not np.allclose(sample, np.round(sample), atol=1e-6)
    is_integer = not has_negatives and not has_decimals

    return is_integer, has_negatives, has_decimals


def ensure_counts_layer(
    adata: "ad.AnnData",
    layer_name: str = "counts",
    error_message: Optional[str] = None,
) -> bool:
    """Ensure a counts layer exists in AnnData, creating from raw if needed.

    This is the single source of truth for counts layer preparation.
    Used by scVI-tools methods (scANVI, Cell2location, etc.) that require
    raw counts in a specific layer.

    Args:
        adata: AnnData object (modified in-place)
        layer_name: Name of the layer to ensure (default: "counts")
        error_message: Custom error message if counts cannot be created

    Returns:
        True if layer was created, False if already existed

    Raises:
        DataNotFoundError: If counts layer cannot be created

    Note:
        When adata has been subsetted to HVGs, this function correctly
        subsets adata.raw to match the current var_names.

    Examples:
        # Ensure counts layer exists before scANVI setup
        ensure_counts_layer(adata_ref)
        scvi.model.SCANVI.setup_anndata(adata_ref, layer="counts", ...)

        # With custom error message
        ensure_counts_layer(adata, error_message="scANVI requires raw counts")
    """
    from .exceptions import DataNotFoundError

    if layer_name in adata.layers:
        return False

    if adata.raw is not None:
        # Get raw counts, subsetting to current var_names
        # Note: adata.raw may have full genes while adata has HVG subset
        adata.layers[layer_name] = adata.raw[:, adata.var_names].X
        return True

    # Cannot create counts layer
    default_error = (
        f"Cannot create '{layer_name}' layer: adata.raw is None. "
        "Load unpreprocessed data or ensure adata.raw is preserved during preprocessing."
    )
    raise DataNotFoundError(error_message or default_error)


class RawDataResult:
    """Result of raw data extraction."""

    def __init__(
        self,
        X: Any,  # sparse or dense matrix
        var_names: pd.Index,
        source: str,
        is_integer_counts: bool,
        has_negatives: bool = False,
        has_decimals: bool = False,
    ):
        self.X = X
        self.var_names = var_names
        self.source = source
        self.is_integer_counts = is_integer_counts
        self.has_negatives = has_negatives
        self.has_decimals = has_decimals


def get_raw_data_source(
    adata: "ad.AnnData",
    prefer_complete_genes: bool = True,
    require_integer_counts: bool = False,
    sample_size: int = 100,
) -> RawDataResult:
    """
    Get raw count data from AnnData using a unified priority order.

    This is THE single source of truth for accessing raw counts data.
    All tools should use this function instead of implementing their own logic.

    Priority order (when prefer_complete_genes=True):
        1. adata.raw - Complete gene set, preserved before HVG filtering
        2. adata.layers["counts"] - Raw counts layer
        3. adata.X - Current expression matrix

    Priority order (when prefer_complete_genes=False):
        1. adata.layers["counts"] - Raw counts layer
        2. adata.X - Current expression matrix
        (adata.raw is skipped as it may have different dimensions)

    Args:
        adata: AnnData object
        prefer_complete_genes: If True, prefer adata.raw for complete gene coverage.
            Set to False when you need data aligned with current adata dimensions.
        require_integer_counts: If True, validate that data contains integer counts.
            Raises DataError if only normalized data is found.
        sample_size: Number of cells/genes to sample for validation.

    Returns:
        RawDataResult with data matrix, var_names, source name, and validation info.

    Raises:
        DataError: If require_integer_counts=True and no integer counts found.

    Example:
        result = get_raw_data_source(adata, prefer_complete_genes=True)
        print(f"Using {result.source}: {len(result.var_names)} genes")
        if result.is_integer_counts:
            # Safe to use for deconvolution/velocity
            pass
    """
    sources_tried = []

    # Source 1: adata.raw (complete gene set)
    if prefer_complete_genes and adata.raw is not None:
        try:
            raw_adata = adata.raw.to_adata()
            is_int, has_neg, has_dec = check_is_integer_counts(raw_adata.X, sample_size)

            if is_int or not require_integer_counts:
                return RawDataResult(
                    X=raw_adata.X,
                    var_names=raw_adata.var_names,
                    source="raw",
                    is_integer_counts=is_int,
                    has_negatives=has_neg,
                    has_decimals=has_dec,
                )
            sources_tried.append("raw (normalized, skipped)")
        except Exception:
            sources_tried.append("raw (error, skipped)")

    # Source 2: layers["counts"]
    if "counts" in adata.layers:
        X_counts = adata.layers["counts"]
        is_int, has_neg, has_dec = check_is_integer_counts(X_counts, sample_size)

        if is_int or not require_integer_counts:
            return RawDataResult(
                X=X_counts,
                var_names=adata.var_names,
                source="counts_layer",
                is_integer_counts=is_int,
                has_negatives=has_neg,
                has_decimals=has_dec,
            )
        sources_tried.append("counts_layer (normalized, skipped)")

    # Source 3: current X
    is_int, has_neg, has_dec = check_is_integer_counts(adata.X, sample_size)

    if is_int or not require_integer_counts:
        return RawDataResult(
            X=adata.X,
            var_names=adata.var_names,
            source="current",
            is_integer_counts=is_int,
            has_negatives=has_neg,
            has_decimals=has_dec,
        )

    # If we get here with require_integer_counts=True, no valid source found
    if require_integer_counts:
        raise DataError(
            f"No raw integer counts found. Sources tried: {sources_tried + ['current (normalized)']}. "
            f"Data appears to be normalized (has_negatives={has_neg}, has_decimals={has_dec}). "
            "Deconvolution and velocity methods require raw integer counts. "
            "Solutions: (1) Load unpreprocessed data, (2) Ensure adata.layers['counts'] "
            "contains raw counts, or (3) Re-run preprocessing with adata.raw preservation."
        )

    # Fallback: return normalized data with warning flags
    return RawDataResult(
        X=adata.X,
        var_names=adata.var_names,
        source="current",
        is_integer_counts=False,
        has_negatives=has_neg,
        has_decimals=has_dec,
    )


# =============================================================================
# Expression Data Extraction: Unified sparse/dense handling
# =============================================================================
def to_dense(X: Any, copy: bool = False) -> np.ndarray:
    """
    Convert sparse matrix to dense numpy array.

    Handles both scipy sparse matrices and already-dense arrays uniformly.
    This is THE single function for sparse-to-dense conversion across ChatSpatial.

    Args:
        X: Expression matrix (sparse or dense)
        copy: If True, always return a copy (safe for modification).
              If False (default), may return view for dense input (read-only use).

    Returns:
        Dense numpy array

    Note:
        - Sparse input: Always returns a new array (toarray() creates copy)
        - Dense input with copy=False: May return view (no memory overhead)
        - Dense input with copy=True: Always returns independent copy

    Examples:
        # Read-only use (default, memory efficient)
        dense_X = to_dense(adata.X)

        # When you need to modify the result
        dense_X = to_dense(adata.X, copy=True)
        dense_X[0, 0] = 999  # Safe, won't affect original
    """
    if sparse.issparse(X):
        return X.toarray()
    # For dense: np.array with copy=False may still copy if needed (e.g., non-contiguous)
    # np.array with copy=True always copies
    return np.array(X, copy=copy)


def get_gene_expression(
    adata: "ad.AnnData",
    gene: str,
    layer: Optional[str] = None,
) -> np.ndarray:
    """
    Get expression values of a single gene as 1D array.

    This is THE single function for extracting single-gene expression.
    Replaces 12+ duplicated code patterns across the codebase.

    Args:
        adata: AnnData object
        gene: Gene name (must exist in adata.var_names)
        layer: Optional layer name. If None, uses adata.X

    Returns:
        1D numpy array of expression values (length = n_obs)

    Raises:
        DataError: If gene not found in adata

    Examples:
        # Basic usage
        cd4_expr = get_gene_expression(adata, "CD4")

        # From specific layer
        counts = get_gene_expression(adata, "CD4", layer="counts")

        # Use in visualization
        adata.obs["_temp_expr"] = get_gene_expression(adata, gene)
    """
    if gene not in adata.var_names:
        raise DataError(
            f"Gene '{gene}' not found in data. "
            f"Available genes (first 5): {adata.var_names[:5].tolist()}"
        )

    if layer is not None:
        if layer not in adata.layers:
            raise DataError(
                f"Layer '{layer}' not found. Available: {list(adata.layers.keys())}"
            )
        gene_idx = adata.var_names.get_loc(gene)
        X = adata.layers[layer][:, gene_idx]
    else:
        X = adata[:, gene].X

    return to_dense(X).flatten()


def get_genes_expression(
    adata: "ad.AnnData",
    genes: List[str],
    layer: Optional[str] = None,
) -> np.ndarray:
    """
    Get expression values of multiple genes as 2D array.

    Args:
        adata: AnnData object
        genes: List of gene names (must exist in adata.var_names)
        layer: Optional layer name. If None, uses adata.X

    Returns:
        2D numpy array of shape (n_obs, n_genes)

    Raises:
        DataError: If any gene not found in adata

    Examples:
        # Get expression matrix for heatmap
        expr_matrix = get_genes_expression(adata, ["CD4", "CD8A", "CD3D"])

        # From counts layer
        counts = get_genes_expression(adata, marker_genes, layer="counts")
    """
    # Validate genes
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise DataError(
            f"Genes not found: {missing[:5]}{'...' if len(missing) > 5 else ''}. "
            f"Available genes (first 5): {adata.var_names[:5].tolist()}"
        )

    if layer is not None:
        if layer not in adata.layers:
            raise DataError(
                f"Layer '{layer}' not found. Available: {list(adata.layers.keys())}"
            )
        gene_indices = [adata.var_names.get_loc(g) for g in genes]
        X = adata.layers[layer][:, gene_indices]
    else:
        X = adata[:, genes].X

    return to_dense(X)


# =============================================================================
# Metadata Profiling: Extract structure information for LLM understanding
# =============================================================================
def get_column_profile(
    adata: "ad.AnnData", layer: Literal["obs", "var"] = "obs"
) -> List[Dict[str, Any]]:
    """
    Get metadata column profile for obs or var.

    Returns detailed information about each column to help LLM understand the data.

    Args:
        adata: AnnData object
        layer: Which layer to profile ("obs" or "var")

    Returns:
        List of column information dictionaries with keys:
        - name: Column name
        - dtype: "numerical" or "categorical"
        - n_unique: Number of unique values
        - range: (min, max) for numerical columns, None for categorical
        - sample_values: Sample values for categorical columns, None for numerical
    """
    df = adata.obs if layer == "obs" else adata.var
    profiles = []

    for col in df.columns:
        col_data = df[col]

        # Determine if numeric
        is_numeric = pd.api.types.is_numeric_dtype(col_data)

        if is_numeric:
            # Numerical column
            profiles.append(
                {
                    "name": col,
                    "dtype": "numerical",
                    "n_unique": int(col_data.nunique()),
                    "range": (float(col_data.min()), float(col_data.max())),
                    "sample_values": None,
                }
            )
        else:
            # Categorical column
            unique_vals = col_data.unique()
            n_unique = len(unique_vals)

            # Take first 5 sample values (or 3 if too many unique values)
            if n_unique <= 100:
                sample_vals = unique_vals[:5].tolist()
            else:
                sample_vals = unique_vals[:3].tolist()

            profiles.append(
                {
                    "name": col,
                    "dtype": "categorical",
                    "n_unique": n_unique,
                    "sample_values": [str(v) for v in sample_vals],
                    "range": None,
                }
            )

    return profiles


def get_gene_profile(
    adata: "ad.AnnData",
) -> Tuple[Optional[List[str]], List[str]]:
    """
    Get gene expression profile including HVGs and top expressed genes.

    Args:
        adata: AnnData object

    Returns:
        Tuple of (top_highly_variable_genes, top_expressed_genes)
        - top_highly_variable_genes: List of HVG names or None if not computed
        - top_expressed_genes: List of top 10 expressed gene names
    """
    # Highly variable genes (no fallback - only return if precomputed)
    hvg_list = get_highly_variable_genes(
        adata, max_genes=10, fallback_to_variance=False
    )
    top_hvg = hvg_list if hvg_list else None

    # Top expressed genes
    try:
        mean_expr = np.array(adata.X.mean(axis=0)).flatten()
        top_idx = np.argsort(mean_expr)[-10:][::-1]  # Descending order
        top_expr = adata.var_names[top_idx].tolist()
    except Exception:
        top_expr = adata.var_names[:10].tolist()  # Fallback

    return top_hvg, top_expr


def get_adata_profile(adata: "ad.AnnData") -> Dict[str, Any]:
    """
    Get comprehensive metadata profile for LLM understanding.

    This is the main function for extracting dataset information that helps
    LLM make informed analysis decisions.

    Args:
        adata: AnnData object

    Returns:
        Dictionary containing:
        - obs_columns: Profile of observation metadata columns
        - var_columns: Profile of variable metadata columns
        - obsm_keys: List of keys in obsm (embeddings, coordinates)
        - uns_keys: List of keys in uns (unstructured annotations)
        - top_highly_variable_genes: Top HVGs if computed
        - top_expressed_genes: Top expressed genes
    """
    # Get column profiles
    obs_profile = get_column_profile(adata, layer="obs")
    var_profile = get_column_profile(adata, layer="var")

    # Get gene profiles
    top_hvg, top_expr = get_gene_profile(adata)

    # Get multi-dimensional data keys
    obsm_keys = list(adata.obsm.keys()) if hasattr(adata, "obsm") else []
    uns_keys = list(adata.uns.keys()) if hasattr(adata, "uns") else []

    return {
        "obs_columns": obs_profile,
        "var_columns": var_profile,
        "obsm_keys": obsm_keys,
        "uns_keys": uns_keys,
        "top_highly_variable_genes": top_hvg,
        "top_expressed_genes": top_expr,
    }


# =============================================================================
# Gene Overlap: Find and validate common genes between datasets
# =============================================================================
def find_common_genes(
    genes_a: List[str],
    genes_b: List[str],
) -> List[str]:
    """
    Find common genes between two gene lists.

    This is THE single function for computing gene intersections across ChatSpatial.

    Args:
        genes_a: First gene list (e.g., adata1.var_names)
        genes_b: Second gene list (e.g., adata2.var_names)

    Returns:
        List of common gene names (order not guaranteed)

    Examples:
        # Between two AnnData objects
        common = find_common_genes(adata1.var_names.tolist(), adata2.var_names.tolist())

        # With Index objects (will be converted)
        common = find_common_genes(list(adata.var_names), list(ref.var_names))
    """
    return list(set(genes_a) & set(genes_b))


def validate_gene_overlap(
    common_genes: List[str],
    source_n_genes: int,
    target_n_genes: int,
    min_genes: int = 100,
    source_name: str = "source",
    target_name: str = "target",
) -> None:
    """
    Validate that gene overlap meets minimum requirements.

    This is THE single validation function for gene overlap across ChatSpatial.
    Moved from deconvolution.py._validate_common_genes for reuse.

    Args:
        common_genes: List of common gene names
        source_n_genes: Number of genes in source data
        target_n_genes: Number of genes in target data
        min_genes: Minimum required common genes (default: 100)
        source_name: Name of source data for error messages
        target_name: Name of target data for error messages

    Raises:
        DataError: If insufficient common genes

    Examples:
        # Basic validation
        common = find_common_genes(spatial.var_names, reference.var_names)
        validate_gene_overlap(common, spatial.n_vars, reference.n_vars)

        # With custom threshold and names
        validate_gene_overlap(
            common, spatial.n_vars, reference.n_vars,
            min_genes=50, source_name="spatial", target_name="reference"
        )
    """
    if len(common_genes) < min_genes:
        raise DataError(
            f"Insufficient gene overlap: {len(common_genes)} < {min_genes} required. "
            f"{source_name}: {source_n_genes} genes, {target_name}: {target_n_genes} genes. "
            f"Check species/gene naming convention match."
        )
