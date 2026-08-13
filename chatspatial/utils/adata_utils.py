"""
AnnData utilities for ChatSpatial.

This module provides:
1. Standard field name constants
2. Field discovery functions (get_*_key)
3. Data access functions (get_*)
4. Validation functions (validate_*)
5. Ensure functions (ensure_*)

One file for all AnnData-related utilities. No duplication.

Naming Conventions (MUST follow across codebase):
-------------------------------------------------
- validate_*(adata, ...) -> None
    Check-only. Raises exception if validation fails.
    Does NOT modify data. Use for precondition checks.
    Example: validate_obs_column(adata, "leiden")

- ensure_*(adata, ...) -> bool
    Check-and-fix. Returns True if action was taken, False if already OK.
    MAY modify data in-place. Idempotent (safe to call multiple times).
    Example: ensure_categorical(adata, "leiden")

- require(name, ctx, feature) -> module
    Dependency check. Raises DependencyError with install or repair guidance.
    Used by optional analysis backends.

Async variants: Add '_async' suffix (e.g., ensure_unique_var_names_async).
"""

from typing import TYPE_CHECKING, Any, Literal, Optional

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad

from scipy import sparse

from .exceptions import DataError, DataNotFoundError, ParameterError

# =============================================================================
# Constants: Standard Field Names
# =============================================================================
SPATIAL_KEY = "spatial"
CELL_TYPE_KEY = "cell_type"
CLUSTER_KEY = "leiden"
BATCH_KEY = "batch"

# Alternative names for compatibility
# Ordered priority lists for deterministic key resolution.
# Do NOT convert these to sets: set iteration order is intentionally unstable.
ALTERNATIVE_SPATIAL_KEYS: tuple[str, ...] = (
    "spatial",
    "X_spatial",
    "coordinates",
    "coords",
    "spatial_coords",
    "positions",
)
ALTERNATIVE_CELL_TYPE_KEYS: tuple[str, ...] = (
    "cell_type",
    "celltype",
    "cell_types",
    "annotation",
    "cell_annotation",
    "predicted_celltype",
)
ALTERNATIVE_CLUSTER_KEYS: tuple[str, ...] = (
    "leiden",
    "louvain",
    "clusters",
    "cluster",
    "clustering",
    "cluster_labels",
    "spatial_domains",
)
ALTERNATIVE_BATCH_KEYS: tuple[str, ...] = (
    "batch",
    "sample",
    "dataset",
    "experiment",
    "replicate",
    "batch_id",
    "sample_id",
)


# =============================================================================
# Field Discovery: Find keys in AnnData
# =============================================================================
def get_spatial_key(adata: "ad.AnnData") -> Optional[str]:
    """Find spatial coordinate key in adata.obsm."""
    for key in ALTERNATIVE_SPATIAL_KEYS:
        if key in adata.obsm:
            return key
    return None


def has_tissue_image(adata: Any) -> bool:
    """Return whether AnnData contains a usable Visium-like tissue image."""
    uns = getattr(adata, "uns", None)
    if uns is None:
        return False
    spatial = uns.get("spatial") if hasattr(uns, "get") else None
    if not isinstance(spatial, dict):
        return False

    for sample_data in spatial.values():
        if not isinstance(sample_data, dict):
            continue
        images = sample_data.get("images")
        if not isinstance(images, dict):
            continue
        if any(images.get(key) is not None for key in ("hires", "lowres")):
            return True
    return False


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

    # Deterministic seed from matrix metadata for reproducibility
    nnz = X.nnz if sparse.issparse(X) else int(np.count_nonzero(X))
    seed = (
        X.shape[0] * 2654435761 ^ X.shape[1] * 40503 ^ nnz * 65537 ^ n_samples
    ) & 0xFFFFFFFF
    rng = np.random.default_rng(seed)

    # Handle sparse matrices efficiently
    if sparse.issparse(X):
        # Sample from .data array (non-zero values only)
        data = X.data
        if len(data) > 0:
            k = min(n_samples, len(data))
            idx = rng.integers(0, len(data), size=k)
            return data[idx]
        else:
            # Empty sparse matrix
            total_values = int(X.shape[0]) * int(X.shape[1])
            return np.zeros(
                min(n_samples, total_values),
                dtype=getattr(X, "dtype", float),
            )
    else:
        # For dense matrices, sample random positions
        arr = np.asarray(X).ravel()
        k = min(n_samples, arr.size)
        idx = rng.integers(0, arr.size, size=k)
        return arr[idx]


def require_spatial_coords(
    adata: "ad.AnnData",
    spatial_key: Optional[str] = None,
    validate: bool = True,
) -> np.ndarray:
    """
    Get validated spatial coordinates array from AnnData.

    Pure read-only accessor and single source of truth for spatial coordinate
    validation. Format unification (obs x/y → obsm['spatial']) is handled
    upstream by _move_spatial_to_standard() at load time. This function
    never mutates adata; it reads from obsm or falls back to obs x/y.

    Args:
        adata: AnnData object
        spatial_key: Optional key in obsm. If None, auto-detects using
            ALTERNATIVE_SPATIAL_KEYS
        validate: If True (default), validates coordinates for:
            - At least 2 dimensions
            - No NaN or Inf values
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
    # Find spatial key if not specified.
    # Normally obs x/y is already in obsm['spatial'] via _move_spatial_to_standard()
    # at load time. If not, check obs x/y as a read-only fallback (no mutation).
    if spatial_key is None:
        spatial_key = get_spatial_key(adata)
        if spatial_key is None:
            # Read-only fallback: parse obs x/y without writing to obsm
            if "x" in adata.obs and "y" in adata.obs:
                try:
                    x = pd.to_numeric(adata.obs["x"], errors="coerce").values
                    y = pd.to_numeric(adata.obs["y"], errors="coerce").values
                    coords = np.column_stack([x, y]).astype("float64")
                    if validate:
                        if np.any(np.isnan(coords)):
                            raise DataError(
                                "Spatial coordinates in obs['x'/'y'] contain NaN"
                            )
                        if np.any(np.isinf(coords)):
                            raise DataError(
                                "Spatial coordinates in obs['x'/'y'] contain "
                                "infinite values"
                            )
                        if np.std(coords[:, 0]) == 0 and np.std(coords[:, 1]) == 0:
                            raise DataError("All spatial coordinates are identical")
                    return coords
                except (TypeError, ValueError) as exc:
                    raise DataError(
                        f"Failed to parse obs['x'/'y'] as coordinates: {exc}"
                    ) from exc
            raise DataNotFoundError(
                "No spatial coordinates found. Expected in adata.obsm['spatial'] "
                "or similar key. Available obsm keys: "
                f"{list(adata.obsm.keys()) if adata.obsm else 'none'}"
            )

    # Check if key exists
    if spatial_key not in adata.obsm:
        raise DataNotFoundError(
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
    check_empty_ratio: bool = False,
    max_empty_obs_ratio: float = 0.1,
    max_empty_vars_ratio: float = 0.5,
) -> None:
    """Validate basic AnnData structure.

    Args:
        adata: AnnData object to validate
        min_obs: Minimum number of observations (cells/spots) required
        min_vars: Minimum number of variables (genes) required
        check_empty_ratio: If True, also check for empty cells/genes
        max_empty_obs_ratio: Max fraction of cells with zero expression (default 10%)
        max_empty_vars_ratio: Max fraction of genes with zero expression (default 50%)

    Raises:
        DataError: If validation fails
    """
    if adata is None:
        raise DataError("AnnData object cannot be None")
    if adata.n_obs < min_obs:
        raise DataError(f"Dataset has {adata.n_obs} observations, need {min_obs}")
    if adata.n_vars < min_vars:
        raise DataError(f"Dataset has {adata.n_vars} variables, need {min_vars}")

    if check_empty_ratio:
        # Count non-zero entries per cell/gene (sparse-aware)
        if sparse.issparse(adata.X):
            cell_nnz = np.array(adata.X.getnnz(axis=1)).flatten()
            gene_nnz = np.array(adata.X.getnnz(axis=0)).flatten()
        else:
            cell_nnz = np.sum(adata.X > 0, axis=1)
            gene_nnz = np.sum(adata.X > 0, axis=0)

        empty_cells = np.sum(cell_nnz == 0)
        empty_genes = np.sum(gene_nnz == 0)

        if empty_cells > adata.n_obs * max_empty_obs_ratio:
            pct = empty_cells / adata.n_obs * 100
            raise DataError(
                f"{empty_cells} cells ({pct:.1f}%) have zero expression. "
                f"Check data quality and consider filtering."
            )

        if empty_genes > adata.n_vars * max_empty_vars_ratio:
            pct = empty_genes / adata.n_vars * 100
            raise DataError(
                f"{empty_genes} genes ({pct:.1f}%) have zero expression. "
                f"Consider gene filtering."
            )


def ensure_categorical(adata: "ad.AnnData", column: str) -> None:
    """Ensure a column is categorical dtype, converting if needed."""
    if column not in adata.obs.columns:
        raise DataError(f"Column '{column}' not found in adata.obs")
    if not isinstance(adata.obs[column].dtype, pd.CategoricalDtype):
        adata.obs[column] = adata.obs[column].astype("category")


def is_categorical_or_string(values: Any) -> bool:
    """Return whether values carry labels rather than continuous magnitudes.

    Pandas 3 uses ``StringDtype`` for inferred text columns instead of object
    dtype. Keep semantic type detection independent of that storage change and
    avoid the deprecated ``is_categorical_dtype`` predicate.
    """
    dtype = getattr(values, "dtype", values)
    return (
        isinstance(dtype, pd.CategoricalDtype)
        or pd.api.types.is_object_dtype(dtype)
        or pd.api.types.is_string_dtype(dtype)
    )


# =============================================================================
# Memory Optimization: Efficient AnnData operations
# =============================================================================
def shallow_copy_adata(adata: "ad.AnnData") -> "ad.AnnData":
    """Create a memory-efficient shallow copy of AnnData.

    This function shares X and existing layers (read-only) but copies obs/uns/var.
    Provides ~99% memory savings compared to full .copy().

    Safe to use when:
    - X will not be modified (e.g., scvi-tools setup_anndata only modifies uns)
    - Existing layers will not be modified (only new layers added)
    - Only obs, uns, var need to be independent

    Memory comparison (4000 cells × 20000 genes):
    - Full .copy():     ~916 MB
    - shallow_copy():   ~0.7 MB
    - Savings:          ~99.9%

    Args:
        adata: Source AnnData object

    Returns:
        New AnnData with shared X/layers, independent obs/uns/var

    Example:
        # Safe for scvi-tools workflows (annotation.py scANVI)
        adata_work = shallow_copy_adata(adata_original)
        adata_work.obs["label"] = "Unknown"  # Doesn't affect original
        scvi.model.SCVI.setup_anndata(adata_work, ...)  # Adds to uns only
    """
    import anndata as ad

    # Create new AnnData with shared X
    adata_new = ad.AnnData(
        X=adata.X,  # Share X (not copied)
        obs=adata.obs.copy(),  # Copy obs (may be modified)
        var=adata.var.copy(),  # Copy var (for safety)
        uns=adata.uns.copy(),  # Copy uns (scvi adds to it)
    )

    # Share obsm arrays (direct assignment shares memory)
    for key in adata.obsm:
        adata_new.obsm[key] = adata.obsm[key]

    # Share varm arrays
    for key in adata.varm:
        adata_new.varm[key] = adata.varm[key]

    # Share layers arrays (direct assignment shares memory)
    # NOTE: adata.layers.copy() does deep copy! Must assign individually.
    for key in adata.layers:
        adata_new.layers[key] = adata.layers[key]

    # Share raw if present
    if adata.raw is not None:
        # AnnData.raw setter only accepts AnnData, but we intentionally keep
        # shallow-copy semantics here to avoid duplicating raw matrices.
        adata_new._raw = adata.raw

    return adata_new


_VELOVI_REQUIRED_KEYS: tuple[str, ...] = (
    "velovi_obs_names",
    "velovi_gene_names",
    "velovi_velocity",
    "velovi_Ms",
    "velovi_connectivities",
)


def _velovi_axis(values: Any, axis_name: str) -> pd.Index:
    """Normalize and validate a persisted VELOVI identity axis."""
    axis = pd.Index([str(value) for value in values])
    if len(axis) == 0:
        raise DataError(f"VELOVI {axis_name} axis is empty.")
    if not axis.is_unique:
        raise DataError(f"VELOVI {axis_name} axis contains duplicate IDs.")
    return axis


def _validate_velovi_matrix(
    matrix: Any, expected_shape: tuple[int, int], key: str
) -> Any:
    """Validate a persisted VELOVI matrix without changing its representation."""
    if getattr(matrix, "shape", None) != expected_shape:
        raise DataError(
            f"{key} has shape {getattr(matrix, 'shape', None)}; "
            f"expected {expected_shape}."
        )
    raw_values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if np.iscomplexobj(raw_values):
        raise DataError(f"{key} contains complex values.")
    try:
        finite = np.isfinite(raw_values)
    except TypeError as exc:
        raise DataError(f"{key} must be numeric.") from exc
    if not finite.all():
        raise DataError(f"{key} contains non-finite values.")
    return matrix


def store_velovi_essential_data(
    adata: "ad.AnnData", adata_velovi: "ad.AnnData"
) -> None:
    """Store the labeled VELOVI payload required to reconstruct CellRank input.

    Stores in adata.uns:
        - velovi_obs_names: cell identities for safe reconstruction
        - velovi_gene_names: filtered gene names
        - velovi_velocity: velocity matrix
        - velovi_Ms: smoothed spliced (required by VelocityKernel)
        - velovi_Mu: smoothed unspliced
        - velovi_connectivities: neighbors graph (sparse)
        - velovi_distances: neighbors distances (sparse)

    Args:
        adata: Original AnnData to store data in
        adata_velovi: Preprocessed velovi AnnData with velocity results
    """
    source_cells = _velovi_axis(adata_velovi.obs_names, "observation")
    target_cells = _velovi_axis(adata.obs_names, "observation")
    if not source_cells.equals(target_cells):
        raise DataError(
            "VELOVI workspace observations do not match the target dataset order."
        )
    gene_names = _velovi_axis(adata_velovi.var_names, "gene")
    matrix_shape = (len(source_cells), len(gene_names))
    graph_shape = (len(source_cells), len(source_cells))

    if "velocity_velovi" in adata_velovi.layers:
        velocity = adata_velovi.layers["velocity_velovi"]
    elif "velocity" in adata_velovi.layers:
        velocity = adata_velovi.layers["velocity"]
    else:
        raise DataNotFoundError(
            "VELOVI workspace is missing the velocity result layer."
        )
    if "Ms" not in adata_velovi.layers:
        raise DataNotFoundError("VELOVI workspace is missing the required Ms layer.")
    if "connectivities" not in adata_velovi.obsp:
        raise DataNotFoundError(
            "VELOVI workspace is missing the required connectivities graph."
        )

    payload: dict[str, Any] = {
        "velovi_obs_names": source_cells.tolist(),
        "velovi_gene_names": gene_names.tolist(),
        "velovi_velocity": _validate_velovi_matrix(
            velocity, matrix_shape, "velovi_velocity"
        ),
        "velovi_Ms": _validate_velovi_matrix(
            adata_velovi.layers["Ms"], matrix_shape, "velovi_Ms"
        ),
        "velovi_connectivities": _validate_velovi_matrix(
            adata_velovi.obsp["connectivities"],
            graph_shape,
            "velovi_connectivities",
        ),
    }
    if "Mu" in adata_velovi.layers:
        payload["velovi_Mu"] = _validate_velovi_matrix(
            adata_velovi.layers["Mu"], matrix_shape, "velovi_Mu"
        )
    if "distances" in adata_velovi.obsp:
        payload["velovi_distances"] = _validate_velovi_matrix(
            adata_velovi.obsp["distances"], graph_shape, "velovi_distances"
        )

    for optional_key in ("velovi_Mu", "velovi_distances"):
        if optional_key not in payload:
            adata.uns.pop(optional_key, None)
    adata.uns.update(payload)


def reconstruct_velovi_adata(adata: "ad.AnnData") -> "ad.AnnData":
    """Reconstruct velovi AnnData from stored essential data.

    This is the inverse of store_velovi_essential_data(). It creates a minimal
    AnnData with all data required by CellRank's VelocityKernel.

    Args:
        adata: AnnData with stored velovi essential data in uns

    Returns:
        Reconstructed velovi AnnData suitable for CellRank

    Raises:
        DataError: If essential velovi data is missing
    """
    import anndata as ad

    missing_keys = [key for key in _VELOVI_REQUIRED_KEYS if key not in adata.uns]
    if missing_keys:
        raise DataNotFoundError(
            "Missing VELOVI reconstruction data: "
            f"{', '.join(missing_keys)}. Run velocity analysis again."
        )

    stored_cells = _velovi_axis(adata.uns["velovi_obs_names"], "observation")
    current_cells = _velovi_axis(adata.obs_names, "observation")
    missing_cells = current_cells.difference(stored_cells)
    if len(missing_cells):
        raise DataError(
            "Current dataset contains observations that were not present in the "
            f"VELOVI analysis: {missing_cells[:5].tolist()}."
        )
    cell_positions = stored_cells.get_indexer(current_cells)

    gene_names = _velovi_axis(adata.uns["velovi_gene_names"], "gene")
    stored_matrix_shape = (len(stored_cells), len(gene_names))
    stored_graph_shape = (len(stored_cells), len(stored_cells))
    velocity = _validate_velovi_matrix(
        adata.uns["velovi_velocity"], stored_matrix_shape, "velovi_velocity"
    )
    ms_data = _validate_velovi_matrix(
        adata.uns["velovi_Ms"], stored_matrix_shape, "velovi_Ms"
    )
    connectivities = _validate_velovi_matrix(
        adata.uns["velovi_connectivities"],
        stored_graph_shape,
        "velovi_connectivities",
    )
    mu_data = (
        _validate_velovi_matrix(
            adata.uns["velovi_Mu"], stored_matrix_shape, "velovi_Mu"
        )
        if "velovi_Mu" in adata.uns
        else None
    )
    distances = (
        _validate_velovi_matrix(
            adata.uns["velovi_distances"],
            stored_graph_shape,
            "velovi_distances",
        )
        if "velovi_distances" in adata.uns
        else None
    )

    n_cells = len(current_cells)
    n_genes = len(gene_names)

    # CellRank's velocity kernel does not consume the expression matrix itself.
    X_placeholder = np.zeros((n_cells, n_genes), dtype=np.float32)

    # Create reconstructed AnnData
    adata_velovi = ad.AnnData(
        X=X_placeholder,
        obs=adata.obs.copy(),
    )
    adata_velovi.var_names = gene_names
    adata_velovi.obs_names = adata.obs_names.copy()

    # Add velocity layer (essential)
    adata_velovi.layers["velocity"] = velocity[cell_positions, :]

    # Add Ms/Mu layers (required by VelocityKernel)
    adata_velovi.layers["Ms"] = ms_data[cell_positions, :]
    if mu_data is not None:
        adata_velovi.layers["Mu"] = mu_data[cell_positions, :]

    # Add neighbors graph (essential for kernels)
    adata_velovi.obsp["connectivities"] = connectivities[cell_positions, :][
        :, cell_positions
    ]
    if distances is not None:
        adata_velovi.obsp["distances"] = distances[cell_positions, :][:, cell_positions]

    # Add neighbors metadata (required by CellRank)
    neighbors: dict[str, Any] = {
        "connectivities_key": "connectivities",
        "params": {"method": "umap"},
    }
    if distances is not None:
        neighbors["distances_key"] = "distances"
    adata_velovi.uns["neighbors"] = neighbors

    # Add spatial coordinates if available
    spatial_key = get_spatial_key(adata)
    if spatial_key and spatial_key in adata.obsm:
        adata_velovi.obsm["spatial"] = adata.obsm[spatial_key]

    return adata_velovi


def has_velovi_essential_data(adata: "ad.AnnData") -> bool:
    """Check if AnnData has essential velovi data for reconstruction."""
    return all(key in adata.uns for key in _VELOVI_REQUIRED_KEYS)


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
    all_categorical_keys = set(
        ALTERNATIVE_CELL_TYPE_KEYS + ALTERNATIVE_CLUSTER_KEYS + ALTERNATIVE_BATCH_KEYS
    )
    for key in adata.obs.columns:
        if key in all_categorical_keys:
            ensure_categorical(adata, key)

    return adata


def _move_spatial_to_standard(adata: "ad.AnnData") -> None:
    """Move spatial coordinates to standard obsm['spatial'] location.

    This function only handles format unification (moving coordinates to
    the canonical obsm['spatial'] key). It does NOT validate coordinate
    values — validation is the sole responsibility of require_spatial_coords().
    """
    if SPATIAL_KEY in adata.obsm:
        return

    # Check alternative obsm keys
    for key in ALTERNATIVE_SPATIAL_KEYS:
        if key in adata.obsm and key != SPATIAL_KEY:
            adata.obsm[SPATIAL_KEY] = adata.obsm[key]
            return

    # Check obs x/y — move unconditionally so downstream validation
    # (require_spatial_coords) can report precise errors instead of
    # "no spatial coordinates found".
    if "x" in adata.obs and "y" in adata.obs:
        try:
            x = pd.to_numeric(adata.obs["x"], errors="coerce").values
            y = pd.to_numeric(adata.obs["y"], errors="coerce").values
            adata.obsm[SPATIAL_KEY] = np.column_stack([x, y]).astype("float64")
        except (TypeError, ValueError):
            return


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

    for category, configured_keys in required_keys.items():
        keys = (
            [configured_keys] if isinstance(configured_keys, str) else configured_keys
        )

        attr = getattr(adata, category, None)
        if attr is None:
            missing.extend([f"{category}.{k}" for k in keys])
            continue

        for key in keys:
            if hasattr(attr, "columns"):  # DataFrame
                if key not in attr.columns:
                    missing.append(f"{category}.{key}")
            elif hasattr(attr, "keys"):  # Dict-like
                if key not in attr:
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
    adata: "ad.AnnData", spatial_key: str, issues: list[str]
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


def _validate_velocity_data(adata: "ad.AnnData", issues: list[str]) -> None:
    """Internal helper for velocity data validation."""
    if "spliced" not in adata.layers:
        issues.append("Missing 'spliced' layer required for RNA velocity")
    if "unspliced" not in adata.layers:
        issues.append("Missing 'unspliced' layer required for RNA velocity")

    if "spliced" in adata.layers and "unspliced" in adata.layers:
        for layer_name in ["spliced", "unspliced"]:
            layer_data = adata.layers[layer_name]
            is_sparse_layer = sparse.issparse(layer_data)

            if is_sparse_layer:
                if layer_data.nnz == 0:
                    issues.append(f"'{layer_name}' layer is empty (all zeros)")
            else:
                if np.all(layer_data == 0):
                    issues.append(f"'{layer_name}' layer is empty (all zeros)")

            if is_sparse_layer:
                if np.any(np.isnan(layer_data.data)):
                    issues.append(f"'{layer_name}' layer contains NaN values")
            else:
                if np.any(np.isnan(layer_data)):
                    issues.append(f"'{layer_name}' layer contains NaN values")


# =============================================================================
# Metadata Storage: Scientific Provenance Tracking
# =============================================================================
def store_analysis_metadata(
    adata: "ad.AnnData",
    analysis_name: str,
    method: str,
    parameters: dict[str, Any],
    results_keys: dict[str, list[str]],
    statistics: Optional[dict[str, Any]] = None,
    species: Optional[str] = None,
    database: Optional[str] = None,
    reference_info: Optional[dict[str, Any]] = None,
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


def get_analysis_metadata_field(
    adata: "ad.AnnData",
    analysis_name: str,
    field_name: str,
    default: Any = None,
) -> Any:
    """Get any top-level field from stored analysis metadata.

    Reads directly from the metadata dict stored by store_analysis_metadata().
    Use this for top-level fields like ``results_keys``, ``statistics``,
    ``method``, ``species``, etc.

    For fields nested under ``parameters``, prefer ``get_analysis_parameter()``.

    Args:
        adata: AnnData object
        analysis_name: Name of the analysis (e.g., "deconvolution_flashdeconv")
        field_name: Top-level key (e.g., "results_keys", "statistics")
        default: Default value if not found

    Returns:
        Field value or default
    """
    metadata_key = f"{analysis_name}_metadata"
    if metadata_key not in adata.uns:
        return default
    return adata.uns[metadata_key].get(field_name, default)


def get_analysis_parameter(
    adata: "ad.AnnData",
    analysis_name: str,
    parameter_name: str,
    default: Any = None,
) -> Any:
    """Get a parameter from stored analysis metadata.

    Retrieves parameters stored by store_analysis_metadata(). Use this to
    access analysis parameters (like cluster_key) without re-inferring them.

    Args:
        adata: AnnData object
        analysis_name: Name of the analysis (e.g., "spatial_stats_neighborhood")
        parameter_name: Name of the parameter (e.g., "cluster_key")
        default: Default value if parameter not found

    Returns:
        Parameter value or default

    Example:
        # Get cluster_key used in neighborhood analysis
        cluster_key = get_analysis_parameter(
            adata, "spatial_stats_neighborhood", "cluster_key"
        )
    """
    metadata_key = f"{analysis_name}_metadata"
    if metadata_key not in adata.uns:
        return default

    metadata = adata.uns[metadata_key]
    if "parameters" not in metadata:
        return default

    return metadata["parameters"].get(parameter_name, default)


# =============================================================================
# Gene Selection Utilities
# =============================================================================
def get_highly_variable_genes(
    adata: "ad.AnnData",
    max_genes: int = 500,
    fallback_to_variance: bool = True,
) -> list[str]:
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

        top_k = min(max_genes, len(var_scores))
        if top_k <= 0:
            return []

        top_indices = np.argpartition(var_scores, -top_k)[-top_k:]
        top_indices = top_indices[np.argsort(var_scores[top_indices])[::-1]]
        return adata.var_names[top_indices].tolist()

    return []


def select_genes_for_analysis(
    adata: "ad.AnnData",
    genes: Optional[list[str]] = None,
    n_genes: int = 20,
    require_hvg: bool = True,
    analysis_name: str = "analysis",
) -> list[str]:
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
def make_unique_names(names: list[str]) -> list[str]:
    """
    Make names unique by adding suffixes (-1, -2, etc.) to duplicates.

    Uses the same suffix convention as scanpy/pandas for consistency across
    the entire codebase. This is THE single source of truth for name deduplication.

    Args:
        names: List of names (may contain duplicates)

    Returns:
        List of unique names with suffixes added to duplicates.
        First occurrence is unchanged, subsequent get -1, -2, etc.

    Examples:
        >>> make_unique_names(["A", "B", "A", "C", "A"])
        ["A", "B", "A-1", "C", "A-2"]

        >>> make_unique_names(["gene1", "gene2"])  # Already unique
        ["gene1", "gene2"]
    """
    if len(names) == len(set(names)):
        return names

    from collections import Counter

    name_counts = Counter(names)
    unique_names = []
    seen_counts: dict[str, int] = {}

    for name in names:
        if name_counts[name] > 1:
            if name not in seen_counts:
                seen_counts[name] = 0
                unique_names.append(name)
            else:
                seen_counts[name] += 1
                # Use "-" suffix to match scanpy convention
                unique_names.append(f"{name}-{seen_counts[name]}")
        else:
            unique_names.append(name)

    return unique_names


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


async def ensure_unique_var_names_async(
    adata: "ad.AnnData",
    ctx: Any,  # ToolContext, use Any to avoid circular import
    label: str = "data",
) -> int:
    """
    Ensure gene names are unique with user feedback via ctx.

    Async variant of ensure_unique_var_names with context logging.

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


def check_is_integer_counts(X: Any, sample_size: int = 1000) -> tuple[bool, bool, bool]:
    """Check if a matrix contains integer counts.

    For sparse matrices, checks ALL non-zero values (typically a small
    fraction of the matrix) to avoid missing rare decimals.
    For dense matrices, uses random sampling with a larger default.

    Args:
        X: Data matrix (sparse or dense)
        sample_size: Number of elements to sample for dense matrices

    Returns:
        Tuple of (is_integer, has_negatives, has_decimals)
    """
    if sparse.issparse(X):
        data = X.data
        if data.size == 0:
            return True, False, False
        # For sparse matrices, non-zero values are typically a small
        # fraction — check them all for reliable detection.
        sample = data
    else:
        arr = np.asarray(X)
        total = arr.size
        if total == 0:
            return True, False, False
        if sample_size >= total:
            sample = arr.ravel()
        else:
            # Use larger sample for big matrices to reduce false negatives
            effective_size = min(max(sample_size, total // 100), total)
            nnz = getattr(X, "nnz", None)
            seed = (
                X.shape[0] * 2654435761
                ^ X.shape[1] * 40503
                ^ int(nnz or 0) * 65537
                ^ effective_size
            ) & 0xFFFFFFFF
            rng = np.random.default_rng(seed)
            flat = arr.ravel()
            idx = rng.integers(0, total, size=effective_size)
            sample = flat[idx]

    has_negatives = bool(np.any(sample < 0))
    has_decimals = not np.allclose(sample, np.round(sample), atol=1e-6, rtol=0)
    is_integer = not has_negatives and not has_decimals

    return is_integer, has_negatives, bool(has_decimals)


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
        # Use intersection to avoid KeyError when gene names don't fully overlap
        shared_genes = adata.var_names.intersection(adata.raw.var_names)
        if len(shared_genes) == len(adata.var_names):
            raw_X = adata.raw[:, adata.var_names].X
        elif len(shared_genes) > 0:
            # Partial overlap: subset to shared genes, then reindex to full var_names
            raw_subset = adata.raw[:, shared_genes].X
            if sparse.issparse(raw_subset):
                raw_X = sparse.lil_matrix(
                    (adata.n_obs, adata.n_vars), dtype=raw_subset.dtype
                )
                col_idx = [adata.var_names.get_loc(g) for g in shared_genes]
                raw_X[:, col_idx] = raw_subset
                raw_X = raw_X.tocsr()
            else:
                raw_X = np.zeros((adata.n_obs, adata.n_vars), dtype=raw_subset.dtype)
                col_idx = [adata.var_names.get_loc(g) for g in shared_genes]
                raw_X[:, col_idx] = raw_subset
        else:
            raw_X = None  # No overlap — skip .raw

        if raw_X is not None:
            is_int, has_neg, _ = check_is_integer_counts(raw_X)
            if is_int and not has_neg:
                adata.layers[layer_name] = raw_X
                return True
        # .raw exists but doesn't contain valid integer counts — fall through

    # Fallback: check if adata.X itself is valid integer counts
    is_int_x, has_neg_x, _ = check_is_integer_counts(adata.X)
    if is_int_x and not has_neg_x:
        adata.layers[layer_name] = adata.X.copy()
        return True

    # Cannot create counts layer — neither .raw, X, nor layers had valid counts
    if adata.raw is not None:
        default_error = (
            f"Cannot create '{layer_name}' layer: adata.raw and adata.X both "
            "contain normalized (non-integer) data. Load unpreprocessed data "
            "or ensure adata.layers['counts'] is set before preprocessing."
        )
    else:
        default_error = (
            f"Cannot create '{layer_name}' layer: adata.raw is None and "
            "adata.X contains normalized (non-integer) data. Load "
            "unpreprocessed data or ensure adata.layers['counts'] is set "
            "before preprocessing."
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
    sample_size: int = 1000,
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
        except Exception as exc:
            sources_tried.append(f"raw (error: {exc.__class__.__name__}, skipped)")

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

    # If we reach here, require_integer_counts=True but no valid source found
    # (line 1012 would have returned if require_integer_counts=False)
    raise DataError(
        f"No raw integer counts found. Sources tried: {sources_tried + ['current (normalized)']}. "
        f"Data appears to be normalized (has_negatives={has_neg}, has_decimals={has_dec}). "
        "Deconvolution and velocity methods require raw integer counts. "
        "Solutions: (1) Load unpreprocessed data, (2) Ensure adata.layers['counts'] "
        "contains raw counts, or (3) Re-run preprocessing with adata.raw preservation."
    )


# =============================================================================
# Expression Data Extraction: Unified sparse/dense handling
# =============================================================================
# Key under which the loader records which platform produced a dataset.
SPATIAL_PLATFORM_KEY = "spatial_platform"

# Platforms whose observations pool several cells. A label assigned to one of
# these describes a mixture, not a cell identity.
_SPOT_BASED_PLATFORMS = frozenset({"visium", "slide_seq"})

# Platforms that resolve individual cells.
_SINGLE_CELL_PLATFORMS = frozenset({"xenium", "merfish", "seqfish"})


def get_spatial_platform(adata: Any) -> Optional[str]:
    """Return the platform recorded when the dataset was loaded.

    Returns None for datasets loaded before this was recorded, or loaded as
    "generic", where the platform genuinely is not known.
    """
    platform = adata.uns.get(SPATIAL_PLATFORM_KEY)
    return str(platform) if platform is not None else None


def is_spot_based_platform(adata: Any) -> Optional[bool]:
    """Whether each observation pools several cells.

    Returns None when the platform is unknown, so callers can distinguish "this
    is a spot platform" from "we cannot tell" — absence of evidence is not
    evidence of single-cell resolution.
    """
    platform = get_spatial_platform(adata)
    if platform is None:
        return None
    if platform in _SPOT_BASED_PLATFORMS:
        return True
    if platform in _SINGLE_CELL_PLATFORMS:
        return False
    return None


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
    genes: list[str],
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
) -> list[dict[str, Any]]:
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
        # Complex values have no meaningful one-dimensional min/max range and
        # JSON cannot represent them as numbers. Expose them as categorical
        # strings instead of silently discarding their imaginary component.
        is_numeric = pd.api.types.is_numeric_dtype(
            col_data
        ) and not pd.api.types.is_complex_dtype(col_data)

        if is_numeric:
            # JSON has no representation for NaN or infinity. Profile the
            # finite observations only and omit the range when a column has no
            # finite values, keeping MCP structured output schema-valid.
            numeric_values = pd.to_numeric(col_data, errors="coerce").to_numpy(
                dtype=float,
                na_value=np.nan,
            )
            finite_values = numeric_values[np.isfinite(numeric_values)]
            value_range = (
                (float(finite_values.min()), float(finite_values.max()))
                if finite_values.size
                else None
            )
            profiles.append(
                {
                    "name": col,
                    "dtype": "numerical",
                    "n_unique": int(np.unique(finite_values).size),
                    "range": value_range,
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
) -> tuple[Optional[list[str]], list[str]]:
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

    if adata.n_vars == 0:
        return top_hvg, []

    try:
        mean_expr = np.asarray(adata.X.mean(axis=0), dtype=float).ravel()
    except (TypeError, ValueError) as exc:
        raise DataError("Expression matrix means must be numeric") from exc

    if mean_expr.size != adata.n_vars:
        raise DataError(
            "Expression matrix mean length does not match adata.n_vars: "
            f"{mean_expr.size} != {adata.n_vars}"
        )

    # Missing/non-finite means are not evidence of high expression. Rank only
    # finite values, deterministically preserving gene order for exact ties.
    finite_indices = np.flatnonzero(np.isfinite(mean_expr))
    top_k = min(10, finite_indices.size)
    if top_k == 0:
        return top_hvg, []

    finite_means = mean_expr[finite_indices]
    ranked = np.argsort(-finite_means, kind="stable")[:top_k]
    top_expr = adata.var_names[finite_indices[ranked]].tolist()

    return top_hvg, top_expr


def get_adata_profile(adata: "ad.AnnData") -> dict[str, Any]:
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
def find_common_genes(*gene_collections: Any) -> list[str]:
    """
    Find common genes across multiple gene collections.

    This is THE single function for computing gene intersections across ChatSpatial.
    Supports any number of gene collections (2 or more).

    Args:
        *gene_collections: Two or more gene collections. Each can be:
            - List[str]: Gene name list
            - pd.Index: AnnData var_names
            - Any Iterable[str]: Will be converted to set

    Returns:
        List of common gene names in first-collection order (stable/deterministic)

    Raises:
        ValueError: If fewer than 2 collections provided

    Examples:
        # Between two AnnData objects
        common = find_common_genes(adata1.var_names, adata2.var_names)

        # Multiple datasets (e.g., spatial registration)
        common = find_common_genes(
            adata1.var_names, adata2.var_names, adata3.var_names
        )

        # With explicit lists
        common = find_common_genes(["GeneA", "GeneB"], ["GeneB", "GeneC"])
    """
    if len(gene_collections) < 2:
        raise ParameterError("find_common_genes requires at least 2 gene collections")

    first_collection = list(gene_collections[0])
    result = set(first_collection)

    # Intersect with remaining collections
    for genes in gene_collections[1:]:
        result &= set(genes)

    if not result:
        return []

    # Preserve first collection order while deduplicating.
    seen: set[str] = set()
    ordered_common: list[str] = []
    for gene in first_collection:
        if gene in result and gene not in seen:
            ordered_common.append(gene)
            seen.add(gene)
    return ordered_common


def validate_gene_overlap(
    common_genes: list[str],
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
