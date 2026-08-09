"""Regression tests for squidpy matrix export label matching."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from chatspatial.utils.results_export import _match_labels


@pytest.mark.unit
def test_match_labels_prefers_all_categories():
    """Matrix sized by all categories should use the full list."""
    all_labs = ["A", "B", "C"]
    obs_labs = ["A", "B"]
    assert _match_labels(all_labs, obs_labs, 3) == all_labs


@pytest.mark.unit
def test_match_labels_falls_back_to_observed():
    """Matrix sized by observed categories should use observed list."""
    all_labs = ["A", "B", "C"]
    obs_labs = ["A", "B"]
    assert _match_labels(all_labs, obs_labs, 2) == obs_labs


@pytest.mark.unit
def test_match_labels_integer_fallback():
    """When neither list matches, integer indices should be returned."""
    all_labs = ["A", "B", "C"]
    obs_labs = ["A", "B"]
    result = _match_labels(all_labs, obs_labs, 5)
    assert result == [0, 1, 2, 3, 4]


@pytest.mark.unit
def test_extract_squidpy_nhood_matrix_with_unused_categories(
    minimal_spatial_adata,
):
    """2x2 matrix with 1 observed category + 1 unused must not ValueError."""
    from chatspatial.utils.results_export import _extract_squidpy_spatial_result

    adata = minimal_spatial_adata.copy()
    # Categorical with 2 levels, only 1 observed
    adata.obs["leiden"] = pd.Categorical(["A"] * adata.n_obs, categories=["A", "B"])
    # Squidpy builds 2x2 matrix (all categories)
    adata.uns["leiden_nhood_enrichment"] = {
        "zscore": np.array([[1.0, 0.5], [0.5, 2.0]]),
    }

    result = _extract_squidpy_spatial_result(
        adata, "leiden_nhood_enrichment", adata.uns["leiden_nhood_enrichment"]
    )
    assert result is not None
    assert result.shape == (2, 2)
    assert list(result.index) == ["A", "B"]


@pytest.mark.unit
def test_extract_squidpy_nhood_matrix_observed_only(
    minimal_spatial_adata,
):
    """1x1 matrix matching only observed categories should also work."""
    from chatspatial.utils.results_export import _extract_squidpy_spatial_result

    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["A"] * adata.n_obs, categories=["A", "B"])
    # Hypothetical 1x1 matrix (observed-only)
    adata.uns["leiden_nhood_enrichment"] = {
        "zscore": np.array([[3.0]]),
    }

    result = _extract_squidpy_spatial_result(
        adata, "leiden_nhood_enrichment", adata.uns["leiden_nhood_enrichment"]
    )
    assert result is not None
    assert result.shape == (1, 1)
    assert list(result.index) == ["A"]
