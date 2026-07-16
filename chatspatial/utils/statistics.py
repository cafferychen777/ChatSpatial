"""Shared statistical post-processing utilities."""

from typing import TYPE_CHECKING, Any, Optional

import numpy as np

from .dependency_manager import require_module

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


def adjust_pvalues(
    pvalues: Any,
    *,
    alpha: float = 0.05,
    method: str = "fdr_bh",
    ctx: Optional["ToolContext"] = None,
    feature: str = "multiple-testing correction",
) -> tuple[np.ndarray, np.ndarray]:
    """Correct finite p-values while preserving invalid positions as non-results.

    The returned arrays always match the one-dimensional input shape. NaN and
    infinite inputs are excluded from the correction family, remain rejected
    as ``False``, and receive ``NaN`` adjusted p-values.
    """
    values = np.asarray(pvalues, dtype=float)
    if values.ndim != 1:
        raise ValueError(f"p-values must be one-dimensional, got shape {values.shape}.")

    reject = np.zeros(values.shape, dtype=bool)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return reject, adjusted

    multitest = require_module(
        "statsmodels",
        "statsmodels.stats.multitest",
        ctx,
        feature=feature,
    )
    valid_reject, valid_adjusted, _, _ = multitest.multipletests(
        values[valid],
        alpha=alpha,
        method=method,
        is_sorted=False,
        returnsorted=False,
    )
    reject[valid] = valid_reject
    adjusted[valid] = valid_adjusted
    return reject, adjusted
