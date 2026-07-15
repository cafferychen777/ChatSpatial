"""Contract tests for core AnnData utility helpers."""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from chatspatial.utils.adata_utils import (
    get_raw_data_source,
    get_spatial_key,
    require_spatial_coords,
    to_dense,
)
from chatspatial.utils.exceptions import DataError, DataNotFoundError


def test_get_spatial_key_returns_known_key(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    assert get_spatial_key(adata) in adata.obsm


def test_require_spatial_coords_falls_back_to_obs_xy(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    adata.obs["x"] = np.arange(adata.n_obs, dtype=float)
    adata.obs["y"] = np.arange(adata.n_obs, dtype=float) + 1.0

    coords = require_spatial_coords(adata)
    assert coords.shape == (adata.n_obs, 2)
    assert np.allclose(coords[:, 1], coords[:, 0] + 1.0)


def test_require_spatial_coords_rejects_identical_points(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obsm["spatial"] = np.ones((adata.n_obs, 2), dtype=float)

    with pytest.raises(DataError, match="identical"):
        require_spatial_coords(adata)


def test_require_spatial_coords_uses_not_found_error_for_missing_data(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]

    with pytest.raises(DataNotFoundError, match="No spatial coordinates found"):
        require_spatial_coords(adata)


def test_get_raw_data_source_prefers_raw_for_complete_genes(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    adata = adata[:, :5].copy()

    out = get_raw_data_source(adata, prefer_complete_genes=True)

    assert out.source == "raw"
    assert len(out.var_names) == 24


def test_get_raw_data_source_uses_counts_when_not_preferring_raw(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    adata.layers["counts"] = np.asarray(adata.X).astype(np.int64)

    out = get_raw_data_source(adata, prefer_complete_genes=False)

    assert out.source == "counts_layer"
    assert len(out.var_names) == adata.n_vars


def test_get_raw_data_source_requires_integer_counts(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.X = (np.asarray(adata.X, dtype=float) + 0.25).astype(np.float64)
    if "counts" in adata.layers:
        del adata.layers["counts"]

    with pytest.raises(DataError, match="No raw integer counts found"):
        get_raw_data_source(
            adata,
            prefer_complete_genes=False,
            require_integer_counts=True,
        )


def test_to_dense_handles_sparse_and_dense_copy_semantics():
    dense = np.arange(6, dtype=float).reshape(2, 3)
    sparse_x = sparse.csr_matrix(dense)

    out_sparse = to_dense(sparse_x)
    assert isinstance(out_sparse, np.ndarray)
    assert np.array_equal(out_sparse, dense)

    out_copy = to_dense(dense, copy=True)
    out_copy[0, 0] = -1
    assert dense[0, 0] == 0
