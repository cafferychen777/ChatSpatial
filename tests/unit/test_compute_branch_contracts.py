"""Branch-focused contracts for compute ensure/has utility functions."""

from __future__ import annotations

import numpy as np
import pytest

from chatspatial.utils import compute


def test_ensure_neighbors_skips_when_graph_already_available(
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    adata.uns["neighbors"] = {"params": {}}
    adata.obsp["connectivities"] = np.eye(adata.n_obs)
    assert compute.ensure_neighbors(adata) is False


def test_ensure_umap_and_leiden_and_louvain_skip_if_present(
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2))
    adata.obs["leiden"] = ["0"] * adata.n_obs
    adata.obs["louvain"] = ["1"] * adata.n_obs

    assert compute.ensure_umap(adata) is False
    assert compute.ensure_leiden(adata, key_added="leiden") is False
    assert compute.ensure_louvain(adata, key_added="louvain") is False


def test_ensure_louvain_computes_and_casts_categorical(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
) -> None:
    adata = minimal_spatial_adata.copy()
    called = {"neighbors": False, "categorical": False}

    def _fake_neighbors(_adata):
        called["neighbors"] = True
        return True

    def _fake_louvain(_adata, **kwargs):
        _adata.obs[kwargs["key_added"]] = ["1"] * _adata.n_obs

    def _fake_categorical(_adata, key):
        called["categorical"] = True
        assert key == "louvain"
        return True

    monkeypatch.setattr(compute, "ensure_neighbors", _fake_neighbors)
    monkeypatch.setattr(compute.sc.tl, "louvain", _fake_louvain)
    monkeypatch.setattr(compute, "ensure_categorical", _fake_categorical)

    assert compute.ensure_louvain(adata, key_added="louvain") is True
    assert called["neighbors"] and called["categorical"]


def test_ensure_diffmap_compute_and_skip_paths(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
) -> None:
    adata = minimal_spatial_adata.copy()
    called = {"neighbors": False, "diffmap": False}

    def _fake_neighbors(_adata):
        called["neighbors"] = True
        _adata.uns["neighbors"] = {}
        _adata.obsp["connectivities"] = np.eye(_adata.n_obs)
        return True

    def _fake_diffmap(_adata, n_comps):
        called["diffmap"] = True
        _adata.obsm["X_diffmap"] = np.zeros((_adata.n_obs, n_comps))

    monkeypatch.setattr(compute, "ensure_neighbors", _fake_neighbors)
    monkeypatch.setattr(compute.sc.tl, "diffmap", _fake_diffmap)

    assert compute.ensure_diffmap(adata, n_comps=5) is True
    assert called["neighbors"] and called["diffmap"]
    assert compute.ensure_diffmap(adata, n_comps=5) is False


def test_ensure_spatial_neighbors_skips_when_precomputed(minimal_spatial_adata) -> None:
    adata = minimal_spatial_adata.copy()
    adata.obsp["spatial_connectivities"] = np.eye(adata.n_obs)
    assert compute.ensure_spatial_neighbors(adata) is False


def test_has_helpers_reflect_current_data_state(minimal_spatial_adata) -> None:
    adata = minimal_spatial_adata.copy()
    assert not compute.has_pca(adata)
    assert not compute.has_neighbors(adata)
    assert not compute.has_umap(adata)
    assert not compute.has_clustering(adata, key="leiden")
    assert not compute.has_spatial_neighbors(adata)
    assert not compute.has_hvg(adata)

    adata.obsm["X_pca"] = np.zeros((adata.n_obs, 3))
    adata.uns["neighbors"] = {"params": {}}
    adata.obsp["connectivities"] = np.eye(adata.n_obs)
    adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2))
    adata.obs["leiden"] = ["0"] * adata.n_obs
    adata.obsp["spatial_connectivities"] = np.eye(adata.n_obs)
    adata.var["highly_variable"] = [True] + [False] * (adata.n_vars - 1)

    assert compute.has_pca(adata)
    assert compute.has_neighbors(adata)
    assert compute.has_umap(adata)
    assert compute.has_clustering(adata, key="leiden")
    assert compute.has_spatial_neighbors(adata)
    assert compute.has_hvg(adata)
