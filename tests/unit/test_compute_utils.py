"""Unit contracts for compute.ensure_* utilities and GMM clustering."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from chatspatial.utils import compute
from chatspatial.utils.exceptions import DataNotFoundError


def test_ensure_pca_skips_when_exists(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_pca"] = np.zeros((adata.n_obs, 2), dtype=float)

    monkeypatch.setattr(
        compute.sc.tl,
        "pca",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("should not run")),
    )
    assert compute.ensure_pca(adata) is False


def test_ensure_pca_computes_with_safe_n_comps(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()[:, :3].copy()
    captured: dict[str, object] = {}

    def _fake_pca(_adata, **kwargs):
        captured.update(kwargs)
        _adata.obsm["X_pca"] = np.zeros((_adata.n_obs, kwargs["n_comps"]))

    monkeypatch.setattr(compute.sc.tl, "pca", _fake_pca)
    out = compute.ensure_pca(adata, n_comps=50)

    assert out is True
    assert captured["n_comps"] == min(50, min(adata.n_obs, adata.n_vars) - 1)
    assert "X_pca" in adata.obsm


def test_ensure_pca_uses_current_scanpy_mask_var_api(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.var["highly_variable"] = [True] * adata.n_vars
    captured: dict[str, object] = {}

    def _fake_pca(
        _adata,
        *,
        n_comps,
        mask_var,
        random_state,
    ):
        captured.update(
            n_comps=n_comps,
            mask_var=mask_var,
            random_state=random_state,
        )
        _adata.obsm["X_pca"] = np.zeros((_adata.n_obs, n_comps))

    monkeypatch.setattr(compute.sc.tl, "pca", _fake_pca)

    assert compute.ensure_pca(adata, n_comps=6, random_state=13) is True
    assert captured == {
        "n_comps": 6,
        "mask_var": "highly_variable",
        "random_state": 13,
    }


def test_ensure_neighbors_calls_prerequisites(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    called = {"ensure_pca": False, "neighbors": False}

    def _fake_ensure_pca(_adata, **_kwargs):
        called["ensure_pca"] = True
        _adata.obsm["X_pca"] = np.zeros((_adata.n_obs, 4))
        return True

    def _fake_neighbors(_adata, **kwargs):
        called["neighbors"] = True
        _adata.uns["neighbors"] = {}
        _adata.obsp["connectivities"] = np.eye(_adata.n_obs)
        assert kwargs["use_rep"] == "X_pca"

    monkeypatch.setattr(compute, "ensure_pca", _fake_ensure_pca)
    monkeypatch.setattr(compute.sc.pp, "neighbors", _fake_neighbors)

    assert compute.ensure_neighbors(adata, use_rep="X_pca") is True
    assert called["ensure_pca"] and called["neighbors"]


def test_ensure_neighbors_clamps_n_pcs_to_available_representation(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    captured: dict[str, object] = {}

    def _fake_ensure_pca(_adata, **kwargs):
        captured["pca_kwargs"] = kwargs
        _adata.obsm["X_pca"] = np.zeros((_adata.n_obs, 4))
        return True

    def _fake_neighbors(_adata, **kwargs):
        captured["neighbors_kwargs"] = kwargs
        _adata.uns["neighbors"] = {}
        _adata.obsp["connectivities"] = np.eye(_adata.n_obs)

    monkeypatch.setattr(compute, "ensure_pca", _fake_ensure_pca)
    monkeypatch.setattr(compute.sc.pp, "neighbors", _fake_neighbors)

    assert compute.ensure_neighbors(adata, n_pcs=30, random_state=17) is True
    assert captured["pca_kwargs"] == {"n_comps": 30, "random_state": 17}
    assert captured["neighbors_kwargs"]["n_pcs"] == 4


def test_ensure_umap_calls_neighbors_then_umap(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    called = {"neighbors": False, "umap": False}

    def _fake_neighbors(_adata):
        called["neighbors"] = True
        _adata.uns["neighbors"] = {}
        _adata.obsp["connectivities"] = np.eye(_adata.n_obs)
        return True

    def _fake_umap(_adata, **_kwargs):
        called["umap"] = True
        _adata.obsm["X_umap"] = np.zeros((_adata.n_obs, 2))

    monkeypatch.setattr(compute, "ensure_neighbors", _fake_neighbors)
    monkeypatch.setattr(compute.sc.tl, "umap", _fake_umap)

    assert compute.ensure_umap(adata) is True
    assert called["neighbors"] and called["umap"]


def test_ensure_leiden_sets_categorical(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    called = {"neighbors": False, "categorical": False}

    def _fake_neighbors(_adata):
        called["neighbors"] = True
        return True

    def _fake_leiden(_adata, **kwargs):
        _adata.obs[kwargs["key_added"]] = ["0"] * _adata.n_obs

    def _fake_categorical(_adata, key):
        called["categorical"] = True
        assert key == "leiden"
        return True

    monkeypatch.setattr(compute, "ensure_neighbors", _fake_neighbors)
    monkeypatch.setattr(compute.sc.tl, "leiden", _fake_leiden)
    monkeypatch.setattr(compute, "ensure_categorical", _fake_categorical)

    assert compute.ensure_leiden(adata, key_added="leiden") is True
    assert called["neighbors"] and called["categorical"]


def test_ensure_spatial_neighbors_requires_spatial_key(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    with pytest.raises(DataNotFoundError, match="Spatial coordinates"):
        compute.ensure_spatial_neighbors(adata)


def test_ensure_spatial_neighbors_grid_and_generic_dispatch(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    calls: list[dict[str, object]] = []

    fake_sq = ModuleType("squidpy")
    fake_sq.gr = SimpleNamespace(
        spatial_neighbors=lambda _adata, **kwargs: calls.append(kwargs)
        or _adata.obsp.__setitem__("spatial_connectivities", np.eye(_adata.n_obs))
    )
    monkeypatch.setitem(sys.modules, "squidpy", fake_sq)

    assert compute.ensure_spatial_neighbors(adata, coord_type="grid", n_rings=2)
    assert calls[-1] == {"coord_type": "grid", "n_rings": 2, "spatial_key": "spatial"}

    del adata.obsp["spatial_connectivities"]
    assert compute.ensure_spatial_neighbors(adata, coord_type="generic", n_neighs=8)
    assert calls[-1] == {
        "coord_type": "generic",
        "n_neighs": 8,
        "spatial_key": "spatial",
    }


def test_top_n_desc_indices_returns_descending_top_k():
    values = np.array([1.0, 9.0, 3.0, 7.0, 5.0])
    out = compute.top_n_desc_indices(values, 3)
    assert out.tolist() == [1, 3, 4]


def test_top_n_desc_indices_sanitizes_non_finite_values():
    values = np.array([1.0, np.nan, np.inf, 2.0])
    out = compute.top_n_desc_indices(values, 2, sanitize_nonfinite=True)
    assert out.tolist() == [3, 0]


def test_gmm_clustering_validates_input_shape_and_cluster_count():
    with pytest.raises(ValueError, match="2D array"):
        compute.gmm_clustering(np.array([1.0, 2.0]), n_clusters=2)
    with pytest.raises(ValueError, match=">= 1"):
        compute.gmm_clustering(np.zeros((5, 2)), n_clusters=0)
    with pytest.raises(ValueError, match="cannot exceed"):
        compute.gmm_clustering(np.zeros((5, 2)), n_clusters=6)


def test_gmm_clustering_returns_one_indexed_labels():
    data = np.vstack(
        [
            np.random.default_rng(0).normal(loc=0, scale=0.1, size=(10, 2)),
            np.random.default_rng(1).normal(loc=3, scale=0.1, size=(10, 2)),
        ]
    )
    labels = compute.gmm_clustering(data, n_clusters=2, random_state=0)
    assert labels.min() >= 1
    assert labels.max() <= 2
    assert len(labels) == 20


class TestEnsureHighlyVariableGenes:
    """Scanpy's mean-binned dispersion is undefined on degenerate inputs.

    Genes that are all zero — routine after merging samples with different gene
    sets — make the bins unusable, and pandas signals that as a KeyError or a
    ValueError depending on where it gives up. Both mean the same thing.
    """

    @pytest.mark.parametrize(
        "binning_error",
        [
            KeyError("[nan] not in index"),
            ValueError(
                "cannot specify integer `bins` when input data contains infinity"
            ),
        ],
        ids=["nan-bin-edges", "infinite-bins"],
    )
    def test_unusable_bins_fall_back_to_variance(
        self, minimal_spatial_adata, monkeypatch, binning_error: Exception
    ):
        adata = minimal_spatial_adata.copy()

        def _raise(*_args, **_kwargs):
            raise binning_error

        monkeypatch.setattr(compute.sc.pp, "highly_variable_genes", _raise)

        used_fallback = compute.ensure_highly_variable_genes(adata, n_top_genes=5)

        assert used_fallback is True
        assert int(adata.var["highly_variable"].sum()) == 5

    def test_reports_no_fallback_when_scanpy_succeeds(
        self, minimal_spatial_adata, monkeypatch
    ):
        adata = minimal_spatial_adata.copy()

        def _mark(a, **_kwargs):
            a.var["highly_variable"] = True

        monkeypatch.setattr(compute.sc.pp, "highly_variable_genes", _mark)

        assert compute.ensure_highly_variable_genes(adata, n_top_genes=5) is False

    def test_unrelated_failures_still_propagate(
        self, minimal_spatial_adata, monkeypatch
    ):
        adata = minimal_spatial_adata.copy()

        def _raise(*_args, **_kwargs):
            raise RuntimeError("scanpy is broken")

        monkeypatch.setattr(compute.sc.pp, "highly_variable_genes", _raise)

        with pytest.raises(RuntimeError, match="scanpy is broken"):
            compute.ensure_highly_variable_genes(adata, n_top_genes=5)

    def test_scanpy_kwargs_are_forwarded(self, minimal_spatial_adata, monkeypatch):
        adata = minimal_spatial_adata.copy()
        seen: dict[str, object] = {}

        def _capture(a, **kwargs):
            seen.update(kwargs)
            a.var["highly_variable"] = True

        monkeypatch.setattr(compute.sc.pp, "highly_variable_genes", _capture)

        compute.ensure_highly_variable_genes(adata, n_top_genes=7, batch_key="batch")

        assert seen["n_top_genes"] == 7
        assert seen["batch_key"] == "batch"

    def test_pearson_residual_flavor_uses_the_experimental_implementation(
        self, minimal_spatial_adata, monkeypatch
    ):
        """Only the experimental entry point accepts this flavor.

        The stable one rejects it with a ValueError that the degenerate-input
        fallback would swallow, answering with plain variance ranking instead.
        """
        adata = minimal_spatial_adata.copy()
        adata.layers["counts"] = adata.X.copy()
        seen: dict[str, object] = {}

        def _stable(*_args, **_kwargs):
            raise ValueError('`flavor` needs to be "seurat" or "cell_ranger"')

        def _experimental(a, **kwargs):
            seen.update(kwargs)
            a.var["highly_variable"] = True

        monkeypatch.setattr(compute.sc.pp, "highly_variable_genes", _stable)
        monkeypatch.setattr(
            compute.sc.experimental.pp, "highly_variable_genes", _experimental
        )

        used_fallback = compute.ensure_highly_variable_genes(
            adata, n_top_genes=5, flavor="pearson_residuals", layer="counts"
        )

        assert used_fallback is False
        assert seen["flavor"] == "pearson_residuals"
        assert seen["layer"] == "counts"


def test_select_hvgs_by_variance_picks_the_most_variable(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.X = np.zeros((adata.n_obs, adata.n_vars), dtype=np.float32)
    # Give two genes real variance; everything else is constant.
    adata.X[:, 3] = np.arange(adata.n_obs, dtype=np.float32)
    adata.X[:, 7] = np.arange(adata.n_obs, dtype=np.float32) * 2

    compute.select_hvgs_by_variance(adata, n_hvgs=2)

    assert list(np.flatnonzero(adata.var["highly_variable"].to_numpy())) == [3, 7]
