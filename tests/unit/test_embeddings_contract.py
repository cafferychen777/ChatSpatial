"""Unit contracts for embedding computation tool."""

from __future__ import annotations

import numpy as np
import pytest

from chatspatial.tools import embeddings as emb
from chatspatial.utils.exceptions import DataNotFoundError


class DummyCtx:
    def __init__(self, adata):
        self._adata = adata
        self.warnings: list[str] = []
        self.set_calls = 0

    async def get_adata(self, _data_id: str):
        return self._adata

    async def set_adata(self, _data_id: str, adata) -> None:
        self._adata = adata
        self.set_calls += 1

    async def warning(self, msg: str):
        self.warnings.append(msg)


@pytest.mark.asyncio
async def test_compute_embeddings_happy_path_leiden_with_metadata_export(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)
    calls: dict[str, object] = {}

    monkeypatch.setattr(emb, "ensure_pca", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(emb, "ensure_neighbors", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(emb, "ensure_diffmap", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(emb, "ensure_spatial_neighbors", lambda *_args, **_kwargs: True)

    def _ensure_leiden(adata_obj, *, key_added, **_kwargs):
        adata_obj.obs[key_added] = ["0"] * (adata_obj.n_obs // 2) + ["1"] * (
            adata_obj.n_obs - adata_obj.n_obs // 2
        )
        return True

    monkeypatch.setattr(emb, "ensure_leiden", _ensure_leiden)
    monkeypatch.setattr(emb, "ensure_louvain", lambda *_args, **_kwargs: False)

    def _store(*_args, **kwargs):
        calls["metadata"] = kwargs

    def _export(_adata, data_id, name):
        calls["export"] = (data_id, name)

    monkeypatch.setattr(emb, "store_analysis_metadata", _store)
    monkeypatch.setattr(emb, "export_analysis_result", _export)

    adata.uns["pca"] = {"variance_ratio": np.array([0.4, 0.2])}
    params = emb.EmbeddingParameters(
        compute_diffmap=True,
        clustering_method="leiden",
        clustering_key="cluster_x",
    )

    out = await emb.compute_embeddings("d1", ctx, params)

    assert out.data_id == "d1"
    assert "PCA" in out.computed
    assert "neighbors" in out.computed
    assert "UMAP" in out.computed
    assert "diffusion map" in out.computed
    assert "spatial neighbors" in out.computed
    assert out.n_clusters == 2
    assert out.pca_variance_ratio == pytest.approx(0.6)
    assert calls["export"] == ("d1", "embeddings_leiden")
    assert calls["metadata"]["analysis_name"] == "embeddings_leiden"
    assert calls["metadata"]["statistics"]["n_clusters"] == 2


@pytest.mark.asyncio
async def test_compute_embeddings_louvain_and_skip_paths(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["louvain_key"] = ["0"] * adata.n_obs
    ctx = DummyCtx(adata)

    monkeypatch.setattr(emb, "ensure_pca", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_neighbors", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_diffmap", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_spatial_neighbors", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_leiden", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_louvain", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "store_analysis_metadata", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(emb, "export_analysis_result", lambda *_args, **_kwargs: None)

    out = await emb.compute_embeddings(
        "d2",
        ctx,
        emb.EmbeddingParameters(
            clustering_method="louvain",
            clustering_key="louvain_key",
            compute_diffmap=True,
        ),
    )

    assert out.computed == []
    assert any("PCA (already exists)" in s for s in out.skipped)
    assert any("neighbors (already exists)" in s for s in out.skipped)
    assert any("UMAP (already exists)" in s for s in out.skipped)
    assert any("diffusion map (already exists)" in s for s in out.skipped)
    assert any("spatial neighbors (already exists)" in s for s in out.skipped)
    assert any("louvain_key (already exists)" in s for s in out.skipped)
    assert out.n_clusters == 1


@pytest.mark.asyncio
async def test_compute_embeddings_force_removes_existing_artifacts(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_pca"] = np.ones((adata.n_obs, 2))
    adata.uns["pca"] = {"variance_ratio": np.array([1.0])}
    adata.uns["neighbors"] = {"params": {}}
    adata.obsp["connectivities"] = np.eye(adata.n_obs)
    adata.obsp["distances"] = np.eye(adata.n_obs)
    adata.obsm["X_umap"] = np.ones((adata.n_obs, 2))
    adata.obs["leiden"] = ["0"] * adata.n_obs
    adata.obsm["X_diffmap"] = np.ones((adata.n_obs, 2))
    adata.obsp["spatial_connectivities"] = np.eye(adata.n_obs)
    adata.obsp["spatial_distances"] = np.eye(adata.n_obs)

    ctx = DummyCtx(adata)
    monkeypatch.setattr(emb, "ensure_pca", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_neighbors", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_diffmap", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_spatial_neighbors", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_leiden", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_louvain", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "store_analysis_metadata", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(emb, "export_analysis_result", lambda *_args, **_kwargs: None)

    await emb.compute_embeddings(
        "d3",
        ctx,
        emb.EmbeddingParameters(
            force=True,
            clustering_key="leiden",
            compute_diffmap=True,
            compute_spatial_neighbors=True,
        ),
    )

    committed = ctx._adata
    assert committed is not adata
    assert "X_pca" not in committed.obsm
    assert "pca" not in committed.uns
    assert "neighbors" not in committed.uns
    assert "connectivities" not in committed.obsp
    assert "distances" not in committed.obsp
    assert "X_umap" not in committed.obsm
    assert "leiden" not in committed.obs
    assert "X_diffmap" not in committed.obsm
    assert "spatial_connectivities" not in committed.obsp
    assert "spatial_distances" not in committed.obsp
    assert ctx.set_calls == 1


@pytest.mark.asyncio
async def test_compute_embeddings_force_failure_preserves_source_artifacts(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_pca"] = np.ones((adata.n_obs, 2))
    adata.uns["pca"] = {"variance_ratio": np.array([1.0])}
    ctx = DummyCtx(adata)

    def _fail_after_force_cleanup(*_args, **_kwargs):
        raise RuntimeError("PCA failed")

    monkeypatch.setattr(emb, "ensure_pca", _fail_after_force_cleanup)

    with pytest.raises(RuntimeError, match="PCA failed"):
        await emb.compute_embeddings(
            "d3",
            ctx,
            emb.EmbeddingParameters(
                force=True,
                compute_neighbors=False,
                compute_umap=False,
                compute_clustering=False,
                compute_spatial_neighbors=False,
            ),
        )

    assert ctx._adata is adata
    assert "X_pca" in adata.obsm
    assert "pca" in adata.uns
    assert ctx.set_calls == 0


@pytest.mark.asyncio
async def test_compute_embeddings_force_removes_dependency_artifacts(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_pca"] = np.ones((adata.n_obs, 2))
    adata.uns["pca"] = {"variance_ratio": np.array([1.0])}
    adata.uns["neighbors"] = {"params": {"n_neighbors": 99}}
    adata.obsp["connectivities"] = np.eye(adata.n_obs)
    adata.obsp["distances"] = np.eye(adata.n_obs)
    ctx = DummyCtx(adata)

    def _ensure_neighbors(adata_obj, **_kwargs):
        assert "X_pca" not in adata_obj.obsm
        assert "pca" not in adata_obj.uns
        assert "neighbors" not in adata_obj.uns
        assert "connectivities" not in adata_obj.obsp
        assert "distances" not in adata_obj.obsp
        return True

    monkeypatch.setattr(emb, "ensure_neighbors", _ensure_neighbors)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_args, **_kwargs: True)

    out = await emb.compute_embeddings(
        "d3",
        ctx,
        emb.EmbeddingParameters(
            force=True,
            compute_pca=False,
            compute_neighbors=False,
            compute_umap=True,
            compute_clustering=False,
            compute_diffmap=False,
            compute_spatial_neighbors=False,
        ),
    )

    assert "neighbors (dependency)" in out.computed
    assert "UMAP" in out.computed


@pytest.mark.asyncio
async def test_compute_embeddings_spatial_neighbor_error_is_non_fatal(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    previous_connectivities = np.eye(adata.n_obs)
    previous_distances = np.eye(adata.n_obs) * 2
    previous_metadata = {"params": {"n_neighs": 4}}
    adata.obsp["spatial_connectivities"] = previous_connectivities
    adata.obsp["spatial_distances"] = previous_distances
    adata.uns["spatial_neighbors"] = previous_metadata
    ctx = DummyCtx(adata)

    monkeypatch.setattr(emb, "ensure_pca", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_neighbors", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_diffmap", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_leiden", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_louvain", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "store_analysis_metadata", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(emb, "export_analysis_result", lambda *_args, **_kwargs: None)

    def _spatial_fail(candidate, *_args, **_kwargs):
        candidate.obsp["spatial_connectivities"] = np.zeros((adata.n_obs, adata.n_obs))
        candidate.uns["spatial_neighbors"] = {"partial": True}
        raise ValueError("no spatial coordinates")

    monkeypatch.setattr(emb, "ensure_spatial_neighbors", _spatial_fail)

    out = await emb.compute_embeddings(
        "d4", ctx, emb.EmbeddingParameters(force=True)
    )
    assert any("spatial neighbors (error: no spatial coordinates)" in s for s in out.skipped)
    assert any("Could not compute spatial neighbors" in w for w in ctx.warnings)
    np.testing.assert_array_equal(
        ctx._adata.obsp["spatial_connectivities"], previous_connectivities
    )
    np.testing.assert_array_equal(
        ctx._adata.obsp["spatial_distances"], previous_distances
    )
    assert ctx._adata.uns["spatial_neighbors"] == previous_metadata


@pytest.mark.asyncio
async def test_compute_embeddings_missing_spatial_coordinates_is_non_fatal(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    ctx = DummyCtx(adata)

    monkeypatch.setattr(emb, "ensure_pca", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_neighbors", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_diffmap", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_leiden", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "ensure_louvain", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(emb, "store_analysis_metadata", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(emb, "export_analysis_result", lambda *_args, **_kwargs: None)

    def _spatial_fail(*_args, **_kwargs):
        raise DataNotFoundError("no spatial coordinates")

    monkeypatch.setattr(emb, "ensure_spatial_neighbors", _spatial_fail)

    out = await emb.compute_embeddings("d4", ctx, emb.EmbeddingParameters())
    assert any(
        "spatial neighbors (error: no spatial coordinates)" in item
        for item in out.skipped
    )
    assert any("Could not compute spatial neighbors" in item for item in ctx.warnings)


@pytest.mark.asyncio
async def test_compute_embeddings_uses_requested_neighbor_params_for_dependencies(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm.clear()
    adata.uns.clear()
    ctx = DummyCtx(adata)
    captured: dict[str, object] = {}

    def _ensure_neighbors(adata_obj, **kwargs):
        captured.update(kwargs)
        adata_obj.uns["neighbors"] = {}
        adata_obj.obsp["connectivities"] = np.eye(adata_obj.n_obs)
        return True

    monkeypatch.setattr(emb, "ensure_neighbors", _ensure_neighbors)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_args, **_kwargs: True)

    out = await emb.compute_embeddings(
        "d4",
        ctx,
        emb.EmbeddingParameters(
            compute_pca=False,
            compute_neighbors=False,
            compute_umap=True,
            compute_clustering=False,
            compute_spatial_neighbors=False,
            n_pcs=12,
            n_neighbors=7,
            random_state=19,
        ),
    )

    assert captured == {"n_neighbors": 7, "n_pcs": 12, "random_state": 19}
    assert "neighbors (dependency)" in out.computed
    assert "UMAP" in out.computed


@pytest.mark.asyncio
async def test_compute_embeddings_louvain_computed_branch_sets_cluster_count(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    monkeypatch.setattr(emb, "ensure_pca", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_neighbors", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_diffmap", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_spatial_neighbors", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_leiden", lambda *_a, **_k: False)

    def _ensure_louvain(adata_obj, *, key_added, **_kwargs):
        adata_obj.obs[key_added] = ["0"] * (adata_obj.n_obs // 2) + ["1"] * (
            adata_obj.n_obs - adata_obj.n_obs // 2
        )
        return True

    monkeypatch.setattr(emb, "ensure_louvain", _ensure_louvain)
    monkeypatch.setattr(emb, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(emb, "export_analysis_result", lambda *_a, **_k: None)

    out = await emb.compute_embeddings(
        "d5",
        ctx,
        emb.EmbeddingParameters(
            clustering_method="louvain",
            clustering_key="louvain_auto",
            compute_diffmap=False,
        ),
    )

    assert out.n_clusters == 2
    assert any("Louvain clustering (2 clusters)" in s for s in out.computed)


@pytest.mark.asyncio
async def test_compute_embeddings_louvain_missing_key_reports_skip_message(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    monkeypatch.setattr(emb, "ensure_pca", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_neighbors", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_umap", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_diffmap", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_spatial_neighbors", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_leiden", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "ensure_louvain", lambda *_a, **_k: False)
    monkeypatch.setattr(emb, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(emb, "export_analysis_result", lambda *_a, **_k: None)

    out = await emb.compute_embeddings(
        "d6",
        ctx,
        emb.EmbeddingParameters(
            clustering_method="louvain",
            clustering_key="not_written",
            compute_diffmap=False,
        ),
    )

    assert any("not_written (missing; clustering not computed)" in s for s in out.skipped)
