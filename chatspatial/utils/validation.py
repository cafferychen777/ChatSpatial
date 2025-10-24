"""
Validation utilities for spatial transcriptomics data.

Centralized validation functions to reduce code duplication across tools.
"""

from typing import List, Optional

import anndata as ad


def validate_obs_column(
    adata: ad.AnnData,
    column: str,
    friendly_name: Optional[str] = None,
) -> None:
    """
    Validate that a column exists in adata.obs.

    Args:
        adata: AnnData object
        column: Column name to check
        friendly_name: Optional friendly name for error message

    Raises:
        ValueError: If column not found

    Examples:
        >>> validate_obs_column(adata, "leiden", "Cluster labels")
        >>> validate_obs_column(adata, "cell_type")
    """
    if column not in adata.obs.columns:
        name = friendly_name or f"Column '{column}'"
        raise ValueError(
            f"{name} not found in observation data (adata.obs).\n"
            f"Available columns: {', '.join(list(adata.obs.columns)[:10])}"
            f"{'...' if len(adata.obs.columns) > 10 else ''}"
        )


def validate_var_column(
    adata: ad.AnnData,
    column: str,
    friendly_name: Optional[str] = None,
) -> None:
    """
    Validate that a column exists in adata.var.

    Args:
        adata: AnnData object
        column: Column name to check
        friendly_name: Optional friendly name for error message

    Raises:
        ValueError: If column not found

    Examples:
        >>> validate_var_column(adata, "highly_variable")
        >>> validate_var_column(adata, "gene_symbols", "Gene symbols")
    """
    if column not in adata.var.columns:
        name = friendly_name or f"Column '{column}'"
        raise ValueError(
            f"{name} not found in variable data (adata.var).\n"
            f"Available columns: {', '.join(list(adata.var.columns)[:10])}"
            f"{'...' if len(adata.var.columns) > 10 else ''}"
        )


def validate_obs_columns(
    adata: ad.AnnData,
    columns: List[str],
) -> None:
    """
    Validate that multiple columns exist in adata.obs.

    Args:
        adata: AnnData object
        columns: List of column names to check

    Raises:
        ValueError: If any column not found

    Examples:
        >>> validate_obs_columns(adata, ["leiden", "cell_type", "batch"])
    """
    missing = [col for col in columns if col not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Columns not found in observation data: {', '.join(missing)}\n"
            f"Available columns: {', '.join(list(adata.obs.columns)[:10])}"
            f"{'...' if len(adata.obs.columns) > 10 else ''}"
        )


def validate_var_columns(
    adata: ad.AnnData,
    columns: List[str],
) -> None:
    """
    Validate that multiple columns exist in adata.var.

    Args:
        adata: AnnData object
        columns: List of column names to check

    Raises:
        ValueError: If any column not found

    Examples:
        >>> validate_var_columns(adata, ["highly_variable", "gene_symbols"])
    """
    missing = [col for col in columns if col not in adata.var.columns]
    if missing:
        raise ValueError(
            f"Columns not found in variable data: {', '.join(missing)}\n"
            f"Available columns: {', '.join(list(adata.var.columns)[:10])}"
            f"{'...' if len(adata.var.columns) > 10 else ''}"
        )


def validate_adata_basics(adata: ad.AnnData, min_obs: int = 1, min_vars: int = 1) -> None:
    """
    Validate basic AnnData structure and size requirements.

    Args:
        adata: AnnData object to validate
        min_obs: Minimum number of observations required
        min_vars: Minimum number of variables required

    Raises:
        ValueError: If basic requirements not met

    Examples:
        >>> validate_adata_basics(adata)
        >>> validate_adata_basics(adata, min_obs=10, min_vars=50)
    """
    if adata is None:
        raise ValueError("AnnData object cannot be None")

    if adata.n_obs < min_obs:
        raise ValueError(
            f"Dataset has {adata.n_obs} observations, but {min_obs} required"
        )

    if adata.n_vars < min_vars:
        raise ValueError(
            f"Dataset has {adata.n_vars} variables, but {min_vars} required"
        )


def validate_categorical(
    adata: ad.AnnData,
    column: str,
    auto_convert: bool = False,
) -> None:
    """
    Validate that a column is categorical dtype.

    Args:
        adata: AnnData object
        column: Column name in adata.obs
        auto_convert: If True, automatically convert to categorical

    Raises:
        ValueError: If column is not categorical and auto_convert=False

    Examples:
        >>> validate_categorical(adata, "leiden")
        >>> validate_categorical(adata, "cell_type", auto_convert=True)
    """
    import pandas as pd

    if column not in adata.obs.columns:
        raise ValueError(f"Column '{column}' not found in adata.obs")

    if not pd.api.types.is_categorical_dtype(adata.obs[column]):
        if auto_convert:
            adata.obs[column] = adata.obs[column].astype("category")
        else:
            raise ValueError(
                f"Column '{column}' must be categorical dtype. "
                f"Current dtype: {adata.obs[column].dtype}"
            )


def ensure_categorical(adata: ad.AnnData, column: str) -> None:
    """
    Ensure a column is categorical, converting if needed.

    Convenience wrapper for validate_categorical with auto_convert=True.

    Args:
        adata: AnnData object
        column: Column name in adata.obs

    Raises:
        ValueError: If column not found

    Examples:
        >>> ensure_categorical(adata, "leiden")
        >>> ensure_categorical(adata, "cell_type")
    """
    validate_categorical(adata, column, auto_convert=True)
