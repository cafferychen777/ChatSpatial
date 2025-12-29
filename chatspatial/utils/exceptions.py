"""
Exception classes for ChatSpatial.

Exception Hierarchy Design Principles:
======================================
1. SEMANTIC CLARITY: Exception type immediately indicates error category
2. MINIMAL HIERARCHY: Only add subclasses when semantically distinct
3. CONSISTENCY: Same semantic error uses same exception type everywhere
4. DEBUGGABILITY: Always preserve exception chain with `raise X from e`

Exception Hierarchy:
====================
ChatSpatialError (base)
├── DataError (data access/availability/format issues)
│   ├── DataNotFoundError (required data missing)
│   └── DataCompatibilityError (format/species mismatch)
├── ParameterError (invalid user input/parameters)
├── ProcessingError (algorithm/computation failures)
└── DependencyError (missing packages/R environment)

Usage Guidelines:
=================

INSTEAD OF ValueError, USE:
---------------------------
- ParameterError: Invalid parameter values, invalid combinations
  Example: "n_clusters must be > 0", "method must be 'leiden' or 'louvain'"

- DataError: Data validation failures (missing columns, wrong format)
  Example: "Column 'cell_type' not found", "Data contains NaN values"

INSTEAD OF RuntimeError, USE:
-----------------------------
- ProcessingError: Algorithm failures, computation errors
  Example: "Clustering failed to converge", "Model training failed"

- DependencyError: Package/environment issues
  Example: "scvi-tools not installed", "R environment not configured"

ALWAYS use exception chaining:
------------------------------
    try:
        result = compute()
    except Exception as e:
        raise ProcessingError(f"Computation failed: {e}") from e

NEVER silently swallow exceptions:
----------------------------------
    # BAD
    except Exception:
        pass

    # GOOD
    except ExpectedError as e:
        logger.debug(f"Expected error: {e}")
        # Handle gracefully
"""


class ChatSpatialError(Exception):
    """Base exception for all ChatSpatial errors.

    All custom exceptions inherit from this class, allowing:
    - Catch all ChatSpatial errors: `except ChatSpatialError`
    - Distinguish from Python built-in exceptions
    """

    pass


# =============================================================================
# Data Errors - Issues with data access, availability, or format
# =============================================================================


class DataError(ChatSpatialError):
    """Data-related errors (missing, invalid format, validation failures).

    Use when:
    - Required data column/key is missing
    - Data format is invalid (wrong dtype, contains NaN/Inf)
    - Data validation fails (too few cells, no spatial coordinates)

    Examples:
        raise DataError("Spatial coordinates contain NaN values")
        raise DataError("Expression matrix is empty")
        raise DataError(f"Column '{col}' not found in adata.obs")
    """

    pass


class DataNotFoundError(DataError):
    """Required data not found in the expected location.

    Use when:
    - Dataset ID not found in registry
    - Required analysis results missing (e.g., deconvolution not run)
    - Expected file/resource doesn't exist

    Examples:
        raise DataNotFoundError(f"Dataset '{data_id}' not found")
        raise DataNotFoundError("Deconvolution results not found")
        raise DataNotFoundError("Velocity results not found in adata.obsm")
    """

    pass


class DataCompatibilityError(DataError):
    """Data format or compatibility issues between datasets.

    Use when:
    - Gene naming convention mismatch (Ensembl vs Symbol)
    - Species mismatch between spatial and reference data
    - Incompatible data formats for integration

    Examples:
        raise DataCompatibilityError("Gene naming mismatch: symbols vs Ensembl")
        raise DataCompatibilityError("Species mismatch: mouse vs human")
    """

    pass


# =============================================================================
# Parameter Errors - Invalid user input or parameter combinations
# =============================================================================


class ParameterError(ChatSpatialError):
    """Invalid parameter values or combinations.

    Use when:
    - Parameter value out of valid range
    - Invalid parameter combination
    - Required parameter missing
    - Parameter type is wrong

    REPLACES: ValueError for parameter validation

    Examples:
        raise ParameterError("n_clusters must be > 0")
        raise ParameterError("method must be one of: 'leiden', 'louvain'")
        raise ParameterError("Cannot use both 'n_clusters' and 'resolution'")
        raise ParameterError("species parameter is required")
    """

    pass


# =============================================================================
# Processing Errors - Algorithm or computation failures
# =============================================================================


class ProcessingError(ChatSpatialError):
    """Errors during analysis processing or computation.

    Use when:
    - Algorithm fails to converge
    - Numerical computation errors
    - Model training failures
    - Unexpected processing failures

    REPLACES: RuntimeError for processing failures

    Examples:
        raise ProcessingError("Leiden clustering failed to converge")
        raise ProcessingError(f"Model training failed: {e}") from e
        raise ProcessingError("PCA computation failed: matrix is singular")
    """

    pass


# =============================================================================
# Dependency Errors - Missing packages or environment issues
# =============================================================================


class DependencyError(ChatSpatialError):
    """Missing or incompatible dependency.

    Use when:
    - Required Python package not installed
    - R environment not configured
    - R package missing
    - Version incompatibility

    REPLACES: ImportError for optional dependencies

    Examples:
        raise DependencyError("scvi-tools required: pip install scvi-tools")
        raise DependencyError("R environment not configured")
        raise DependencyError("CellChat R package not installed")
    """

    pass
