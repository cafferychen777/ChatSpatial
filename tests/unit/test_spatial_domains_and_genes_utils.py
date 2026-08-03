"""Unit tests for lightweight helpers in spatial_domains and spatial_genes."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from chatspatial.models.data import SpatialDomainParameters
from chatspatial.tools import spatial_domains as sd
from chatspatial.tools import spatial_genes as sg
from chatspatial.tools.spatial_domains import _refine_spatial_domains
from chatspatial.tools.spatial_genes import _calculate_sparse_gene_stats
from chatspatial.utils.exceptions import (
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


@pytest.fixture(autouse=True)
def _resolve_fake_spatial_domain_submodules(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        sd,
        "require_module",
        lambda _name, module_name, *_args, **_kwargs: __import__("sys").modules[
            module_name
        ],
    )


def _required_dependency(name: str, *_args, **_kwargs):
    return __import__("sys").modules.get(name, object())


def _patch_stagate_dependencies(
    monkeypatch: pytest.MonkeyPatch,
    torch_module: object,
    stagate_module: object,
) -> None:
    dependencies = {"torch": torch_module, "STAGATE_pyG": stagate_module}

    def _require(name: str, *_args, **_kwargs):
        try:
            return dependencies[name]
        except KeyError as exc:
            raise AssertionError(f"Unexpected dependency request: {name}") from exc

    monkeypatch.setattr(sd, "require", _require)


class DummyCtx:
    def __init__(self, adata):
        self.adata = adata
        self.warnings: list[str] = []
        self.debug_logs: list[str] = []
        self.set_calls = 0

    async def get_adata(self, _data_id: str):
        return self.adata

    async def set_adata(self, _data_id: str, adata) -> None:
        self.adata = adata
        self.set_calls += 1

    async def warning(self, msg: str):
        self.warnings.append(msg)

    def debug(self, msg: str):
        self.debug_logs.append(msg)


def test_calculate_sparse_gene_stats_handles_dense_and_sparse_equivalently(
    minimal_spatial_adata,
):
    X_dense = np.asarray(minimal_spatial_adata.X)
    X_sparse = sp.csr_matrix(X_dense)

    dense_totals, dense_expr = _calculate_sparse_gene_stats(X_dense)
    sparse_totals, sparse_expr = _calculate_sparse_gene_stats(X_sparse)

    np.testing.assert_allclose(dense_totals, sparse_totals)
    np.testing.assert_array_equal(dense_expr, sparse_expr)


def test_top_n_indices_returns_descending_and_handles_bounds():
    values = np.array([1.0, 7.0, 4.0, 7.5], dtype=float)

    top2 = sg._top_n_indices(values, 2)
    assert list(top2) == [3, 1]

    top_all = sg._top_n_indices(values, 10)
    assert list(top_all) == [3, 1, 2, 0]

    top_zero = sg._top_n_indices(values, 0)
    assert top_zero.size == 0


def test_refine_spatial_domains_returns_same_length_series(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["domain"] = ["0"] * 20 + ["1"] * 20 + ["2"] * 20

    refined = _refine_spatial_domains(adata, "domain", threshold=0.5)
    assert len(refined) == adata.n_obs
    assert set(refined.unique()).issubset({"0", "1", "2"})


@pytest.mark.asyncio
async def test_identify_spatial_domains_missing_spatial_coordinates_raises_processing_error(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    params = SpatialDomainParameters(method="leiden", refine_domains=False)

    with pytest.raises(DataNotFoundError, match="No spatial coordinates found"):
        await sd.identify_spatial_domains("d1", DummyCtx(adata), params)


@pytest.mark.asyncio
async def test_identify_spatial_domains_leiden_happy_path_stores_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)
    captured: dict[str, object] = {}

    async def _fake_clustering(adata_subset, params, _ctx):
        adata_subset.obsm["X_fake"] = np.zeros((adata_subset.n_obs, 2))
        labels = ["0"] * (adata_subset.n_obs // 2) + ["1"] * (
            adata_subset.n_obs - adata_subset.n_obs // 2
        )
        return labels, "X_fake", {"silhouette": 0.5}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(
        sd,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_args, **_kwargs: [])

    out = await sd.identify_spatial_domains(
        "d1",
        ctx,
        SpatialDomainParameters(method="leiden", refine_domains=False),
    )
    assert out.method == "leiden"
    assert out.domain_key == "spatial_domains_leiden_res0_50"
    assert out.refined_domain_key is None
    assert captured["analysis_name"] == "spatial_domains_leiden_res0_50"
    assert captured["results_keys"] == {
        "obs": ["spatial_domains_leiden_res0_50"],
        "obsm": ["X_fake"],
    }
    assert ctx.adata is not adata
    assert out.domain_key in ctx.adata.obs
    assert ctx.set_calls == 1


@pytest.mark.asyncio
async def test_identify_spatial_domains_refinement_failure_does_not_abort(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def _fake_clustering(adata_subset, params, _ctx):
        labels = ["0"] * (adata_subset.n_obs // 2) + ["1"] * (
            adata_subset.n_obs - adata_subset.n_obs // 2
        )
        return labels, None, {"silhouette": 0.4}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)

    def _partially_refine_then_fail(candidate, domain_key, **_kwargs):
        candidate.obs[f"{domain_key}_refined"] = "partial"
        raise RuntimeError("refine boom")

    monkeypatch.setattr(sd, "_refine_spatial_domains", _partially_refine_then_fail)
    monkeypatch.setattr(sd, "store_analysis_metadata", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_args, **_kwargs: [])

    out = await sd.identify_spatial_domains(
        "d1",
        ctx,
        SpatialDomainParameters(method="leiden", refine_domains=True),
    )
    assert out.refined_domain_key is None
    assert any("Domain refinement failed" in w for w in ctx.warnings)
    assert "spatial_domains_leiden_res0_50_refined" not in ctx.adata.obs


@pytest.mark.asyncio
async def test_identify_spatial_domains_late_failure_preserves_managed_source(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def _fake_clustering(adata_subset, params, _ctx):
        del params, _ctx
        return ["0"] * adata_subset.n_obs, None, {"score": 1.0}

    def _partially_store_then_fail(candidate, **_kwargs):
        candidate.uns["partial_domain_metadata"] = True
        raise RuntimeError("metadata failed")

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(sd, "store_analysis_metadata", _partially_store_then_fail)

    with pytest.raises(ProcessingError, match="metadata failed"):
        await sd.identify_spatial_domains(
            "d1",
            ctx,
            SpatialDomainParameters(method="leiden", refine_domains=False),
        )

    assert ctx.adata is adata
    assert "spatial_domains_leiden_res0_50" not in adata.obs
    assert "partial_domain_metadata" not in adata.uns
    assert ctx.set_calls == 0


@pytest.mark.asyncio
async def test_identify_spatial_genes_routes_to_spatialde(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def _fake_spatialde(data_id, adata_obj, params, _ctx):
        return sg.SpatialVariableGenesResult(
            data_id=data_id,
            method="spatialde",
            n_genes_analyzed=10,
            n_significant_genes=2,
            spatial_genes=["gene_0", "gene_1"],
            results_key="spatialde_results_d1",
        )

    monkeypatch.setattr(sg, "_identify_spatial_genes_spatialde", _fake_spatialde)

    out = await sg.identify_spatial_genes(
        "d1",
        ctx,
        sg.SpatialVariableGenesParameters(method="spatialde", spatial_key="spatial"),
    )
    assert out.method == "spatialde"
    assert out.n_significant_genes == 2


@pytest.mark.asyncio
async def test_identify_spatial_genes_routes_to_sparkx(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def _fake_sparkx(data_id, adata_obj, params, _ctx):
        return sg.SpatialVariableGenesResult(
            data_id=data_id,
            method="sparkx",
            n_genes_analyzed=8,
            n_significant_genes=1,
            spatial_genes=["gene_2"],
            results_key="sparkx_results_d1",
        )

    monkeypatch.setattr(sg, "_identify_spatial_genes_sparkx", _fake_sparkx)

    out = await sg.identify_spatial_genes(
        "d1",
        ctx,
        sg.SpatialVariableGenesParameters(method="sparkx", spatial_key="spatial"),
    )
    assert out.method == "sparkx"
    assert out.spatial_genes == ["gene_2"]


@pytest.mark.asyncio
async def test_identify_spatial_genes_routes_to_flashs(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def _fake_flashs(data_id, adata_obj, params, _ctx):
        return sg.SpatialVariableGenesResult(
            data_id=data_id,
            method="flashs",
            n_genes_analyzed=12,
            n_significant_genes=3,
            spatial_genes=["gene_1", "gene_4", "gene_8"],
            results_key="flashs_results_d1",
        )

    monkeypatch.setattr(sg, "_identify_spatial_genes_flashs", _fake_flashs)

    out = await sg.identify_spatial_genes(
        "d1",
        ctx,
        sg.SpatialVariableGenesParameters(method="flashs", spatial_key="spatial"),
    )
    assert out.method == "flashs"
    assert out.n_significant_genes == 3


@pytest.mark.asyncio
async def test_identify_spatial_genes_late_failure_preserves_managed_source(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def _partially_publish_then_fail(data_id, candidate, params, _ctx):
        del data_id, params, _ctx
        candidate.var["partial_spatial_score"] = 1.0
        candidate.uns["partial_spatial_metadata"] = True
        raise RuntimeError("spatial export failed")

    monkeypatch.setattr(
        sg,
        "_identify_spatial_genes_flashs",
        _partially_publish_then_fail,
    )

    with pytest.raises(RuntimeError, match="spatial export failed"):
        await sg.identify_spatial_genes(
            "d1",
            ctx,
            sg.SpatialVariableGenesParameters(method="flashs"),
        )

    assert ctx.adata is adata
    assert "partial_spatial_score" not in adata.var
    assert "partial_spatial_metadata" not in adata.uns
    assert ctx.set_calls == 0


@pytest.mark.asyncio
async def test_identify_spatial_genes_rejects_unknown_method_via_runtime_guard(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)
    params = sg.SpatialVariableGenesParameters(method="sparkx").model_copy(
        update={"method": "unknown_method"}
    )

    with pytest.raises(ParameterError, match="Unsupported method"):
        await sg.identify_spatial_genes("d1", ctx, params)


@pytest.mark.asyncio
async def test_identify_spatial_genes_requires_spatial_coordinates(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    ctx = DummyCtx(adata)

    with pytest.raises(Exception, match="Spatial coordinates"):
        await sg.identify_spatial_genes(
            "d1", ctx, sg.SpatialVariableGenesParameters(method="sparkx")
        )


@pytest.mark.asyncio
async def test_identify_spatial_domains_unknown_method_is_wrapped(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    params = SpatialDomainParameters(method="leiden").model_copy(
        update={"method": "unknown"}
    )

    with pytest.raises(ParameterError, match="Unsupported method"):
        await sd.identify_spatial_domains("d1", DummyCtx(adata), params)


@pytest.mark.asyncio
async def test_identify_domains_clustering_louvain_falls_back_to_leiden(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)
    # _identify_domains_clustering expects expression graph to exist.
    adata.obsp["connectivities"] = np.eye(adata.n_obs, dtype=float)

    monkeypatch.setattr(sd, "ensure_pca", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "ensure_neighbors", lambda *_a, **_k: None)

    class _FakeSq:
        class gr:
            @staticmethod
            def spatial_neighbors(a, **_kw):
                a.obsp["spatial_connectivities"] = np.eye(a.n_obs, dtype=float)

    monkeypatch.setattr(sd, "require", lambda *_a, **_k: _FakeSq)
    monkeypatch.setattr(
        sd.sc.tl,
        "louvain",
        lambda *_a, **_k: (_ for _ in ()).throw(ImportError("missing louvain")),
    )

    def _fake_leiden(a, resolution=1.0, key_added="spatial_leiden"):
        del resolution
        a.obs[key_added] = ["0"] * a.n_obs

    monkeypatch.setattr(sd.sc.tl, "leiden", _fake_leiden)

    labels, emb_key, stats = await sd._identify_domains_clustering(
        adata,
        SpatialDomainParameters(method="louvain", resolution=0.5),
        ctx,
    )

    assert len(labels) == adata.n_obs
    assert emb_key == "X_pca"
    assert stats["method"] == "louvain"
    assert any("deprecated" in w.lower() for w in ctx.warnings)
    assert any("using leiden clustering instead" in w.lower() for w in ctx.warnings)


@pytest.mark.asyncio
async def test_identify_domains_clustering_spatial_graph_failure_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    monkeypatch.setattr(sd, "ensure_pca", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "ensure_neighbors", lambda *_a, **_k: None)

    class _FakeSq:
        class gr:
            @staticmethod
            def spatial_neighbors(_a, **_kw):
                raise RuntimeError("sq boom")

    monkeypatch.setattr(sd, "require", lambda *_a, **_k: _FakeSq)

    with pytest.raises(ProcessingError, match="Spatial graph construction failed"):
        await sd._identify_domains_clustering(
            adata,
            SpatialDomainParameters(method="leiden"),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_clustering_preserves_dependency_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(sd, "ensure_pca", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(sd, "ensure_neighbors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        sd,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("squidpy is unavailable")
        ),
    )

    with pytest.raises(DependencyError, match="squidpy is unavailable"):
        await sd._identify_domains_clustering(
            adata,
            SpatialDomainParameters(method="leiden"),
            DummyCtx(adata),
        )


def test_refine_spatial_domains_single_spot_preserves_coordinate_error(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata[:1].copy()
    adata.obs["domain"] = ["A"]

    with pytest.raises(DataError, match="All spatial coordinates are identical"):
        _refine_spatial_domains(adata, "domain", threshold=0.5)


@pytest.mark.asyncio
async def test_identify_spatial_domains_cleans_dense_nan_inf_before_clustering(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    dense = np.asarray(adata.X, dtype=float)
    dense[0, 0] = np.nan
    dense[1, 1] = np.inf
    adata.X = dense

    async def _fake_clustering(adata_subset, params, _ctx):
        del params
        assert np.isfinite(np.asarray(adata_subset.X)).all()
        return ["0"] * adata_subset.n_obs, None, {"ok": 1}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(sd, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_a, **_k: [])

    ctx = DummyCtx(adata)
    out = await sd.identify_spatial_domains(
        "d1", ctx, SpatialDomainParameters(method="leiden", refine_domains=False)
    )

    assert out.n_domains == 1
    assert any("NaN or infinite values" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_identify_spatial_domains_uses_raw_when_current_has_negatives(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    adata.X = np.asarray(adata.X, dtype=float)
    adata.X[0, 0] = -1.0

    async def _fake_clustering(adata_subset, params, _ctx):
        del params
        assert np.min(np.asarray(adata_subset.X)) >= 0.0
        return ["0"] * adata_subset.n_obs, None, {"ok": 1}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(sd, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_a, **_k: [])

    ctx = DummyCtx(adata)
    await sd.identify_spatial_domains(
        "d1", ctx, SpatialDomainParameters(method="leiden", refine_domains=False)
    )

    assert any("negative values" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_spagcn_applies_compatibility_patch_before_import(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    events: list[str] = []

    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spagcn_compat",
        lambda: events.append("compat"),
    )

    def _stop_after_import(name: str, *_args, **_kwargs):
        events.append(name)
        raise DependencyError("stop after import order check")

    monkeypatch.setattr(sd, "require", _stop_after_import)

    with pytest.raises(DependencyError, match="import order check"):
        await sd._identify_domains_spagcn(
            minimal_spatial_adata.copy(),
            SpatialDomainParameters(method="spagcn"),
            DummyCtx(minimal_spatial_adata),
        )

    assert events == ["compat", "SpaGCN"]


@pytest.mark.asyncio
async def test_identify_domains_spagcn_success_with_dummy_histology_fallback(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    class _FakeSpg:
        @staticmethod
        def prefilter_genes(_adata, min_cells=3):
            del _adata, min_cells

        @staticmethod
        def prefilter_specialgenes(_adata):
            del _adata

    monkeypatch.setattr(sd, "require", lambda *_a, **_k: _FakeSpg)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spagcn_compat", lambda *_a, **_k: None
    )

    import sys
    import types

    ez_mod = types.ModuleType("SpaGCN.ez_mode")
    ez_mod.detect_spatial_domains_ez_mode = lambda ad, *args, **kwargs: ["0"] * ad.n_obs
    pkg = types.ModuleType("SpaGCN")
    pkg.ez_mode = ez_mod
    sys.modules["SpaGCN"] = pkg
    sys.modules["SpaGCN.ez_mode"] = ez_mod

    params = SpatialDomainParameters(
        method="spagcn",
        n_domains=20,
        spagcn_use_histology=True,
    )

    labels, emb_key, stats = await sd._identify_domains_spagcn(adata, params, ctx)

    assert len(labels) == adata.n_obs
    assert emb_key is None
    assert stats["method"] == "spagcn"
    assert stats["use_histology"] is False
    assert any("unstable or noisy" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_identify_domains_spagcn_timeout_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class _FakeSpg:
        @staticmethod
        def prefilter_genes(_adata, min_cells=3):
            del _adata, min_cells

        @staticmethod
        def prefilter_specialgenes(_adata):
            del _adata

    monkeypatch.setattr(sd, "require", lambda *_a, **_k: _FakeSpg)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spagcn_compat", lambda *_a, **_k: None
    )

    import asyncio as _asyncio
    import sys
    import types

    ez_mod = types.ModuleType("SpaGCN.ez_mode")
    ez_mod.detect_spatial_domains_ez_mode = lambda ad, *args, **kwargs: ["0"] * ad.n_obs
    pkg = types.ModuleType("SpaGCN")
    pkg.ez_mode = ez_mod
    sys.modules["SpaGCN"] = pkg
    sys.modules["SpaGCN.ez_mode"] = ez_mod

    async def _raise_timeout(_future, timeout=None):
        del _future, timeout
        raise _asyncio.TimeoutError()

    monkeypatch.setattr("asyncio.wait_for", _raise_timeout)

    with pytest.raises(ProcessingError, match="SpaGCN timed out"):
        await sd._identify_domains_spagcn(
            adata,
            SpatialDomainParameters(method="spagcn", timeout=1),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_spagcn_rejects_mismatched_coordinate_length(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class _FakeSpg:
        @staticmethod
        def prefilter_genes(_adata, min_cells=3):
            del _adata, min_cells

        @staticmethod
        def prefilter_specialgenes(_adata):
            del _adata

    monkeypatch.setattr(sd, "require", lambda *_a, **_k: _FakeSpg)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spagcn_compat", lambda *_a, **_k: None
    )
    monkeypatch.setattr(
        sd,
        "require_spatial_coords",
        lambda _adata: np.asarray(_adata.obsm["spatial"])[:-1],
    )

    import sys
    import types

    ez_mod = types.ModuleType("SpaGCN.ez_mode")
    ez_mod.detect_spatial_domains_ez_mode = lambda ad, *args, **kwargs: ["0"] * ad.n_obs
    pkg = types.ModuleType("SpaGCN")
    pkg.ez_mode = ez_mod
    sys.modules["SpaGCN"] = pkg
    sys.modules["SpaGCN.ez_mode"] = ez_mod

    with pytest.raises(DataError, match="doesn't match data"):
        await sd._identify_domains_spagcn(
            adata,
            SpatialDomainParameters(method="spagcn"),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_stagate_rejects_unsupported_torch_version(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import types

    adata = minimal_spatial_adata.copy()
    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.9.0"
    torch_mod.device = lambda x: x
    _patch_stagate_dependencies(monkeypatch, torch_mod, object())

    with pytest.raises(ProcessingError, match="requires PyTorch <= 2.8.0"):
        await sd._identify_domains_stagate(
            adata,
            SpatialDomainParameters(method="stagate"),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_stagate_success_returns_embeddings_and_stats(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: f"device:{x}"
    sys.modules["torch"] = torch_mod

    class _FakeSTAGATE:
        @staticmethod
        def Cal_Spatial_Net(a, rad_cutoff=50):
            a.uns["rad_cutoff"] = rad_cutoff

        @staticmethod
        def Stats_Spatial_Net(_a):
            return None

        @staticmethod
        def train_STAGATE(a, device=None):
            del device
            a.obsm["STAGATE"] = np.ones((a.n_obs, 4), dtype=float)
            return a

    _patch_stagate_dependencies(monkeypatch, torch_mod, _FakeSTAGATE)
    monkeypatch.setattr(
        "chatspatial.utils.compute.gmm_clustering",
        lambda data, n_clusters, **_kwargs: np.arange(data.shape[0]) % n_clusters,
    )

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    labels, emb_key, stats = await sd._identify_domains_stagate(
        adata,
        SpatialDomainParameters(
            method="stagate",
            n_domains=3,
            stagate_use_gpu=False,
            stagate_rad_cutoff=40,
        ),
        DummyCtx(adata),
    )

    assert len(labels) == adata.n_obs
    assert emb_key == "STAGATE"
    assert "STAGATE" in adata.obsm
    assert stats["method"] == "stagate_pyg"
    assert stats["target_n_clusters"] == 3
    assert stats["rad_cutoff"] == 40


@pytest.mark.asyncio
async def test_identify_domains_stagate_stats_failure_logs_debug_and_continues(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: f"device:{x}"
    sys.modules["torch"] = torch_mod

    class _FakeSTAGATE:
        @staticmethod
        def Cal_Spatial_Net(a, rad_cutoff=50):
            a.uns["rad_cutoff"] = rad_cutoff

        @staticmethod
        def Stats_Spatial_Net(_a):
            raise RuntimeError("stats boom")

        @staticmethod
        def train_STAGATE(a, device=None):
            del device
            a.obsm["STAGATE"] = np.ones((a.n_obs, 2), dtype=float)
            return a

    _patch_stagate_dependencies(monkeypatch, torch_mod, _FakeSTAGATE)
    monkeypatch.setattr(
        "chatspatial.utils.compute.gmm_clustering",
        lambda data, n_clusters, **_kwargs: np.arange(data.shape[0]) % n_clusters,
    )

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    labels, emb_key, stats = await sd._identify_domains_stagate(
        adata,
        SpatialDomainParameters(method="stagate", n_domains=2, stagate_use_gpu=False),
        ctx,
    )

    assert len(labels) == adata.n_obs
    assert emb_key == "STAGATE"
    assert stats["method"] == "stagate_pyg"
    assert any("stats boom" in msg for msg in ctx.debug_logs)


@pytest.mark.asyncio
async def test_identify_domains_stagate_timeout_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import asyncio as _asyncio
    import sys
    import types

    adata = minimal_spatial_adata.copy()

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: x
    sys.modules["torch"] = torch_mod

    class _FakeSTAGATE:
        @staticmethod
        def Cal_Spatial_Net(_a, rad_cutoff=50):
            del rad_cutoff

        @staticmethod
        def train_STAGATE(a, device=None):
            del a, device
            return None

    _patch_stagate_dependencies(monkeypatch, torch_mod, _FakeSTAGATE)

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    async def _raise_timeout(_future, timeout=None):
        del _future, timeout
        raise _asyncio.TimeoutError()

    monkeypatch.setattr("asyncio.wait_for", _raise_timeout)

    with pytest.raises(ProcessingError, match="STAGATE training timeout"):
        await sd._identify_domains_stagate(
            adata,
            SpatialDomainParameters(method="stagate", timeout=1),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_graphst_mclust_path_success(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types
    import warnings

    adata = minimal_spatial_adata.copy()

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: f"device:{x}"
    sys.modules["torch"] = torch_mod

    class _FakeGraphST:
        def __init__(self, adata_graphst, device=None, random_seed=0):
            del device, random_seed
            self._adata = adata_graphst

        def train(self):
            self._adata.obsm["emb"] = np.ones((self._adata.n_obs, 24), dtype=float)
            return self._adata

    graphst_sub = types.ModuleType("GraphST.GraphST")
    graphst_sub.GraphST = _FakeGraphST
    graphst_pkg = types.ModuleType("GraphST")
    graphst_pkg.GraphST = graphst_sub
    sys.modules["GraphST"] = graphst_pkg
    sys.modules["GraphST.GraphST"] = graphst_sub

    monkeypatch.setattr(sd, "require", _required_dependency)
    monkeypatch.setattr(
        "chatspatial.utils.compute.gmm_clustering",
        lambda data, n_clusters, **_kwargs: np.arange(data.shape[0]) % n_clusters,
    )

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in divide",
            category=RuntimeWarning,
        )
        labels, emb_key, stats = await sd._identify_domains_graphst(
            adata,
            SpatialDomainParameters(
                method="graphst",
                graphst_clustering_method="mclust",
                graphst_refinement=False,
                graphst_n_clusters=4,
                n_domains=4,
            ),
            DummyCtx(adata),
        )

    assert len(labels) == adata.n_obs
    assert emb_key == "emb"
    assert "emb" in adata.obsm
    assert stats["method"] == "graphst"
    assert stats["clustering_method"] == "mclust"
    assert stats["n_clusters"] == 4


@pytest.mark.asyncio
async def test_identify_domains_graphst_timeout_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import asyncio as _asyncio
    import sys
    import types

    adata = minimal_spatial_adata.copy()

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: x
    sys.modules["torch"] = torch_mod

    class _FakeGraphST:
        def __init__(self, adata_graphst, device=None, random_seed=0):
            del adata_graphst, device, random_seed

        def train(self):
            return adata

    graphst_sub = types.ModuleType("GraphST.GraphST")
    graphst_sub.GraphST = _FakeGraphST
    graphst_pkg = types.ModuleType("GraphST")
    graphst_pkg.GraphST = graphst_sub
    sys.modules["GraphST"] = graphst_pkg
    sys.modules["GraphST.GraphST"] = graphst_sub

    monkeypatch.setattr(sd, "require", _required_dependency)

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    async def _raise_timeout(_future, timeout=None):
        del _future, timeout
        raise _asyncio.TimeoutError()

    monkeypatch.setattr("asyncio.wait_for", _raise_timeout)

    with pytest.raises(ProcessingError, match="GraphST training timeout"):
        await sd._identify_domains_graphst(
            adata,
            SpatialDomainParameters(method="graphst", timeout=1),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_banksy_success_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()
    banksy_matrix = adata.copy()

    init_mod = types.ModuleType("banksy.initialize_banksy")
    init_mod.initialize_banksy = lambda *_a, **_k: {"ok": True}
    embed_mod = types.ModuleType("banksy.embed_banksy")
    embed_mod.generate_banksy_matrix = lambda *_a, **_k: (None, banksy_matrix)

    sys.modules["banksy"] = types.ModuleType("banksy")
    sys.modules["banksy.initialize_banksy"] = init_mod
    sys.modules["banksy.embed_banksy"] = embed_mod

    monkeypatch.setattr(sd, "require", _required_dependency)
    monkeypatch.setattr(
        sd.sc.pp,
        "pca",
        lambda a, n_comps=20: a.obsm.__setitem__("X_pca", np.ones((a.n_obs, n_comps))),
    )
    monkeypatch.setattr(sd.sc.pp, "neighbors", lambda *_a, **_k: None)
    monkeypatch.setattr(
        sd.sc.tl,
        "leiden",
        lambda a, resolution=0.5, key_added="banksy_cluster": a.obs.__setitem__(
            key_added, ["0"] * (a.n_obs // 2) + ["1"] * (a.n_obs - a.n_obs // 2)
        ),
    )

    labels, emb_key, stats = await sd._identify_domains_banksy(
        adata,
        SpatialDomainParameters(method="banksy", banksy_pca_dims=5, n_domains=2),
        DummyCtx(adata),
    )

    assert len(labels) == adata.n_obs
    assert emb_key == "X_banksy_pca"
    assert "X_banksy_pca" in adata.obsm
    assert stats["method"] == "banksy"
    assert stats["n_clusters"] == 2
    assert stats["feature_expansion"].endswith("x")


@pytest.mark.asyncio
async def test_identify_domains_banksy_timeout_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import asyncio as _asyncio
    import sys
    import types

    adata = minimal_spatial_adata.copy()

    init_mod = types.ModuleType("banksy.initialize_banksy")
    init_mod.initialize_banksy = lambda *_a, **_k: {"ok": True}
    embed_mod = types.ModuleType("banksy.embed_banksy")
    embed_mod.generate_banksy_matrix = lambda *_a, **_k: (None, adata.copy())

    sys.modules["banksy"] = types.ModuleType("banksy")
    sys.modules["banksy.initialize_banksy"] = init_mod
    sys.modules["banksy.embed_banksy"] = embed_mod

    monkeypatch.setattr(sd, "require", _required_dependency)

    async def _raise_timeout(_future, timeout=None):
        del _future, timeout
        raise _asyncio.TimeoutError()

    monkeypatch.setattr("asyncio.wait_for", _raise_timeout)

    with pytest.raises(ProcessingError, match="BANKSY timeout"):
        await sd._identify_domains_banksy(
            adata,
            SpatialDomainParameters(method="banksy", timeout=1),
            DummyCtx(adata),
        )


def _install_fake_spagcn_modules(monkeypatch: pytest.MonkeyPatch, detect_fn) -> None:
    import sys
    import types

    ez_mod = types.ModuleType("SpaGCN.ez_mode")
    ez_mod.detect_spatial_domains_ez_mode = detect_fn
    pkg = types.ModuleType("SpaGCN")
    pkg.ez_mode = ez_mod

    monkeypatch.setitem(sys.modules, "SpaGCN", pkg)
    monkeypatch.setitem(sys.modules, "SpaGCN.ez_mode", ez_mod)


def _install_fake_banksy_modules(
    monkeypatch: pytest.MonkeyPatch,
    initialize_fn,
    generate_fn,
) -> None:
    import sys
    import types

    init_mod = types.ModuleType("banksy.initialize_banksy")
    init_mod.initialize_banksy = initialize_fn
    embed_mod = types.ModuleType("banksy.embed_banksy")
    embed_mod.generate_banksy_matrix = generate_fn

    monkeypatch.setitem(sys.modules, "banksy", types.ModuleType("banksy"))
    monkeypatch.setitem(sys.modules, "banksy.initialize_banksy", init_mod)
    monkeypatch.setitem(sys.modules, "banksy.embed_banksy", embed_mod)


@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("method", "helper_name"),
    [
        ("spagcn", "_identify_domains_spagcn"),
        ("stagate", "_identify_domains_stagate"),
        ("graphst", "_identify_domains_graphst"),
        ("banksy", "_identify_domains_banksy"),
    ],
)
async def test_identify_spatial_domains_dispatches_non_clustering_methods(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
    method: str,
    helper_name: str,
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)
    labels = ["0"] * adata.n_obs
    called = {"value": False}

    async def _fake_method(*_args, **_kwargs):
        called["value"] = True
        return labels, None, {"method": method}

    monkeypatch.setattr(sd, helper_name, _fake_method)
    monkeypatch.setattr(sd, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_a, **_k: [])

    out = await sd.identify_spatial_domains(
        "d1",
        ctx,
        SpatialDomainParameters(method=method, refine_domains=False),
    )

    assert called["value"] is True
    assert out.method == method
    # Non-clustering methods use n_domains suffix (default n_domains=7)
    assert out.domain_key == f"spatial_domains_{method}_n7"


@pytest.mark.asyncio
async def test_identify_spatial_domains_refinement_success_updates_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)
    captured: dict[str, object] = {}

    async def _fake_clustering(adata_subset, params, _ctx):
        del params, _ctx
        return ["0"] * adata_subset.n_obs, None, {"score": 1.0}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(
        sd,
        "_refine_spatial_domains",
        lambda a, domain_key, threshold=0.5: np.array(["1"] * a.n_obs, dtype=object),
    )
    monkeypatch.setattr(
        sd,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_a, **_k: [])

    out = await sd.identify_spatial_domains(
        "d1",
        ctx,
        SpatialDomainParameters(method="leiden", refine_domains=True),
    )

    assert out.refined_domain_key == "spatial_domains_leiden_res0_50_refined"
    assert captured["results_keys"] == {
        "obs": [
            "spatial_domains_leiden_res0_50",
            "spatial_domains_leiden_res0_50_refined",
        ]
    }


@pytest.mark.asyncio
async def test_identify_spatial_domains_sparse_zero_matrix_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.X = sp.csr_matrix((adata.n_obs, adata.n_vars), dtype=np.float32)

    async def _fake_clustering(adata_subset, params, _ctx):
        del params, _ctx
        assert sp.issparse(adata_subset.X)
        assert adata_subset.X.data.size == 0
        return ["0"] * adata_subset.n_obs, None, {"ok": 1}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(sd, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_a, **_k: [])

    out = await sd.identify_spatial_domains(
        "d1",
        DummyCtx(adata),
        SpatialDomainParameters(method="leiden", refine_domains=False),
    )

    assert out.n_domains == 1


@pytest.mark.asyncio
async def test_identify_spatial_domains_hvg_subset_casts_to_float_and_warns_large_values(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    mask = np.arange(adata.n_vars) % 2 == 0
    adata.var["highly_variable"] = mask
    adata.X = np.asarray(adata.X, dtype=np.int16)
    adata.X[0, 0] = 120

    async def _fake_clustering(adata_subset, params, _ctx):
        del params, _ctx
        assert adata_subset.n_vars == int(mask.sum())
        assert adata_subset.X.dtype in (np.float32, np.float64)
        return ["0"] * adata_subset.n_obs, None, {"ok": 1}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(sd, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_a, **_k: [])

    ctx = DummyCtx(adata)
    await sd.identify_spatial_domains(
        "d1",
        ctx,
        SpatialDomainParameters(
            method="leiden",
            refine_domains=False,
            use_highly_variable=True,
        ),
    )

    assert any("large values" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_identify_spatial_domains_uses_raw_hvg_subset_with_negative_current_matrix(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.raw = adata.copy()
    mask = np.arange(adata.n_vars) % 2 == 1
    adata.var["highly_variable"] = mask
    adata.X = np.asarray(adata.X, dtype=float)
    adata.X[0, 0] = -2.0

    expected_hvg = list(adata.var_names[mask])

    async def _fake_clustering(adata_subset, params, _ctx):
        del params, _ctx
        assert list(adata_subset.var_names) == expected_hvg
        assert float(np.min(np.asarray(adata_subset.X))) >= 0.0
        return ["0"] * adata_subset.n_obs, None, {"ok": 1}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(sd, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_a, **_k: [])

    ctx = DummyCtx(adata)
    await sd.identify_spatial_domains(
        "d1",
        ctx,
        SpatialDomainParameters(
            method="leiden",
            refine_domains=False,
            use_highly_variable=True,
        ),
    )

    assert any("negative values" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_identify_spatial_domains_cleans_sparse_nan_inf(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    dense = np.asarray(adata.X, dtype=float)
    dense[0, 0] = np.nan
    dense[1, 1] = np.inf
    adata.X = sp.csr_matrix(dense)

    async def _fake_clustering(adata_subset, params, _ctx):
        del params, _ctx
        assert sp.issparse(adata_subset.X)
        assert np.isfinite(adata_subset.X.data).all()
        return ["0"] * adata_subset.n_obs, None, {"ok": 1}

    monkeypatch.setattr(sd, "_identify_domains_clustering", _fake_clustering)
    monkeypatch.setattr(sd, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "export_analysis_result", lambda *_a, **_k: [])

    ctx = DummyCtx(adata)
    await sd.identify_spatial_domains(
        "d1",
        ctx,
        SpatialDomainParameters(method="leiden", refine_domains=False),
    )

    assert any("sparse data" in w for w in ctx.warnings)


@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("images", "scalefactors", "expected_scale"),
    [
        (
            {"hires": np.ones((5, 5, 3), dtype=np.uint8)},
            {"tissue_hires_scalef": 2.0},
            2.0,
        ),
        (
            {"lowres": np.ones((5, 5, 3), dtype=np.uint8)},
            {"tissue_lowres_scalef": 3.0},
            3.0,
        ),
        ({"hires": np.ones((5, 5, 3), dtype=np.uint8)}, {}, 1.0),
        ({"lowres": np.ones((5, 5, 3), dtype=np.uint8)}, {}, 1.0),
    ],
)
async def test_identify_domains_spagcn_histology_image_selection_and_scaling(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
    images: dict,
    scalefactors: dict,
    expected_scale: float,
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)
    adata.uns["spatial"] = {"lib1": {"images": images, "scalefactors": scalefactors}}
    captured: dict[str, object] = {}

    class _FakeSpg:
        @staticmethod
        def prefilter_genes(_adata, min_cells=3):
            del _adata, min_cells

        @staticmethod
        def prefilter_specialgenes(_adata):
            del _adata

    def _detect_fn(_adata, _img, x_array, _y_array, x_pixel, _y_pixel, **_kwargs):
        captured["x_array"] = x_array
        captured["x_pixel"] = x_pixel
        return ["0"] * _adata.n_obs

    monkeypatch.setattr(sd, "require", lambda *_a, **_k: _FakeSpg)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spagcn_compat", lambda *_a, **_k: None
    )
    _install_fake_spagcn_modules(monkeypatch, _detect_fn)

    labels, emb_key, stats = await sd._identify_domains_spagcn(
        adata,
        SpatialDomainParameters(method="spagcn", spagcn_use_histology=True),
        ctx,
    )

    expected_pixels = [int(x * expected_scale) for x in captured["x_array"]]
    assert list(captured["x_pixel"]) == expected_pixels
    assert len(labels) == adata.n_obs
    assert emb_key is None
    assert stats["use_histology"] is True


@pytest.mark.asyncio
async def test_identify_domains_spagcn_prefilter_failure_warns_and_continues(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    class _FakeSpg:
        @staticmethod
        def prefilter_genes(_adata, min_cells=3):
            del _adata, min_cells
            raise RuntimeError("prefilter failed")

        @staticmethod
        def prefilter_specialgenes(_adata):
            del _adata

    monkeypatch.setattr(sd, "require", lambda *_a, **_k: _FakeSpg)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spagcn_compat", lambda *_a, **_k: None
    )
    _install_fake_spagcn_modules(monkeypatch, lambda ad, *_a, **_k: ["0"] * ad.n_obs)

    labels, _, _ = await sd._identify_domains_spagcn(
        adata,
        SpatialDomainParameters(method="spagcn"),
        ctx,
    )

    assert len(labels) == adata.n_obs
    assert any("gene filtering failed" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_identify_domains_clustering_leiden_without_spatial_sets_zero_spatial_weight(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]

    monkeypatch.setattr(sd, "ensure_pca", lambda *_a, **_k: None)
    monkeypatch.setattr(sd, "ensure_neighbors", lambda *_a, **_k: None)

    def _fake_leiden(a, resolution=1.0, key_added="spatial_leiden"):
        del resolution
        a.obs[key_added] = ["0"] * a.n_obs

    monkeypatch.setattr(sd.sc.tl, "leiden", _fake_leiden)

    labels, _, stats = await sd._identify_domains_clustering(
        adata,
        SpatialDomainParameters(method="leiden"),
        DummyCtx(adata),
    )

    assert len(labels) == adata.n_obs
    assert stats["spatial_weight"] == 0.0


def test_refine_spatial_domains_empty_dataset_preserves_data_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata[:0].copy()
    adata.obs["domain"] = []
    monkeypatch.setattr(sd, "require_spatial_coords", lambda _adata: np.empty((0, 2)))

    with pytest.raises(DataNotFoundError, match="Dataset is empty"):
        _refine_spatial_domains(adata, "domain", threshold=0.5)


def test_refine_spatial_domains_single_spot_returns_original_when_coords_are_mocked(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata[:1].copy()
    adata.obs["domain"] = ["A"]
    monkeypatch.setattr(
        sd, "require_spatial_coords", lambda _adata: np.array([[1.0, 2.0]])
    )

    refined = _refine_spatial_domains(adata, "domain", threshold=0.5)
    assert refined.tolist() == ["A"]


def test_refine_spatial_domains_nn_failure_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["domain"] = ["A"] * adata.n_obs

    class _BrokenNN:
        def __init__(self, n_neighbors):
            del n_neighbors

        def fit(self, coords):
            del coords
            raise RuntimeError("nn boom")

    monkeypatch.setattr("sklearn.neighbors.NearestNeighbors", _BrokenNN)

    with pytest.raises(ProcessingError, match="Nearest neighbors computation failed"):
        _refine_spatial_domains(adata, "domain", threshold=0.5)


@pytest.mark.asyncio
async def test_identify_domains_stagate_generic_error_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()
    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: x
    monkeypatch.setitem(sys.modules, "torch", torch_mod)

    class _FakeSTAGATE:
        @staticmethod
        def Cal_Spatial_Net(_adata, rad_cutoff=50):
            del _adata, rad_cutoff
            raise RuntimeError("cal graph failed")

    _patch_stagate_dependencies(monkeypatch, torch_mod, _FakeSTAGATE)

    with pytest.raises(ProcessingError, match="STAGATE execution failed"):
        await sd._identify_domains_stagate(
            adata,
            SpatialDomainParameters(method="stagate"),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_graphst_leiden_refinement_branch_and_radius_stats(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: x
    monkeypatch.setitem(sys.modules, "torch", torch_mod)

    class _FakeGraphST:
        def __init__(self, adata_graphst, device=None, random_seed=0):
            del device, random_seed
            self._adata = adata_graphst

        def train(self):
            rng = np.random.default_rng(7)
            self._adata.obsm["emb"] = rng.normal(size=(self._adata.n_obs, 24))
            return self._adata

    graphst_sub = types.ModuleType("GraphST.GraphST")
    graphst_sub.GraphST = _FakeGraphST
    graphst_pkg = types.ModuleType("GraphST")
    graphst_pkg.GraphST = graphst_sub
    graphst_utils = types.ModuleType("GraphST.utils")
    graphst_utils.refine_label = lambda a, radius=50, key="domain": [
        f"r_{v}" for v in a.obs[key].astype(str).tolist()
    ]

    monkeypatch.setitem(sys.modules, "GraphST", graphst_pkg)
    monkeypatch.setitem(sys.modules, "GraphST.GraphST", graphst_sub)
    monkeypatch.setitem(sys.modules, "GraphST.utils", graphst_utils)

    monkeypatch.setattr(sd, "require", _required_dependency)
    monkeypatch.setattr(sd.sc.pp, "neighbors", lambda *_a, **_k: None)

    calls = {"n": 0}

    def _fake_leiden(a, resolution=1.0, random_state=0):
        del resolution, random_state
        calls["n"] += 1
        if calls["n"] == 1:
            n = 2  # < target
        elif calls["n"] == 2:
            n = 4  # > target
        else:
            n = 3  # == target
        a.obs["leiden"] = [str(i % n) for i in range(a.n_obs)]

    monkeypatch.setattr(sd.sc.tl, "leiden", _fake_leiden)

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    labels, emb_key, stats = await sd._identify_domains_graphst(
        adata,
        SpatialDomainParameters(
            method="graphst",
            n_domains=3,
            graphst_n_clusters=3,
            graphst_clustering_method="leiden",
            graphst_refinement=True,
            graphst_radius=77,
        ),
        DummyCtx(adata),
    )

    assert len(labels) == adata.n_obs
    assert all(label.startswith("r_") for label in labels[:5])
    assert emb_key == "emb"
    assert stats["refinement"] is True
    assert stats["refinement_radius"] == 77


@pytest.mark.asyncio
async def test_identify_domains_banksy_missing_spatial_coordinates_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    adata.uns.pop("spatial", None)

    _install_fake_banksy_modules(
        monkeypatch,
        initialize_fn=lambda *_a, **_k: {"ok": True},
        generate_fn=lambda *_a, **_k: (None, adata.copy()),
    )
    monkeypatch.setattr(sd, "require", _required_dependency)

    with pytest.raises(ProcessingError, match="No spatial coordinates found"):
        await sd._identify_domains_banksy(
            adata,
            SpatialDomainParameters(method="banksy"),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_banksy_uses_alternative_spatial_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_spatial"] = adata.obsm.pop("spatial")
    banksy_matrix = adata.copy()

    _install_fake_banksy_modules(
        monkeypatch,
        initialize_fn=lambda *_a, **_k: {"ok": True},
        generate_fn=lambda *_a, **_k: (None, banksy_matrix),
    )
    monkeypatch.setattr(sd, "require", _required_dependency)
    monkeypatch.setattr(
        sd.sc.pp,
        "pca",
        lambda a, n_comps=20: a.obsm.__setitem__("X_pca", np.ones((a.n_obs, n_comps))),
    )
    monkeypatch.setattr(sd.sc.pp, "neighbors", lambda *_a, **_k: None)
    monkeypatch.setattr(
        sd.sc.tl,
        "leiden",
        lambda a, resolution=0.5, key_added="banksy_cluster": a.obs.__setitem__(
            key_added, ["0"] * a.n_obs
        ),
    )

    labels, emb_key, _stats = await sd._identify_domains_banksy(
        adata,
        SpatialDomainParameters(method="banksy", banksy_pca_dims=5),
        DummyCtx(adata),
    )

    assert len(labels) == adata.n_obs
    assert emb_key == "X_banksy_pca"
    assert "spatial" in adata.obsm


@pytest.mark.asyncio
async def test_identify_domains_banksy_generic_error_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    _install_fake_banksy_modules(
        monkeypatch,
        initialize_fn=lambda *_a, **_k: (_ for _ in ()).throw(
            RuntimeError("init boom")
        ),
        generate_fn=lambda *_a, **_k: (None, adata.copy()),
    )
    monkeypatch.setattr(sd, "require", _required_dependency)

    with pytest.raises(ProcessingError, match="BANKSY execution failed"):
        await sd._identify_domains_banksy(
            adata,
            SpatialDomainParameters(method="banksy"),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_identify_domains_graphst_mclust_refinement_branch(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: x
    monkeypatch.setitem(sys.modules, "torch", torch_mod)

    class _FakeGraphST:
        def __init__(self, adata_graphst, device=None, random_seed=0):
            del device, random_seed
            self._adata = adata_graphst

        def train(self):
            rng = np.random.default_rng(11)
            self._adata.obsm["emb"] = rng.normal(size=(self._adata.n_obs, 24))
            return self._adata

    graphst_sub = types.ModuleType("GraphST.GraphST")
    graphst_sub.GraphST = _FakeGraphST
    graphst_pkg = types.ModuleType("GraphST")
    graphst_pkg.GraphST = graphst_sub
    graphst_utils = types.ModuleType("GraphST.utils")
    graphst_utils.refine_label = lambda a, radius=50, key="domain": [
        f"m_{v}" for v in a.obs[key].astype(str).tolist()
    ]

    monkeypatch.setitem(sys.modules, "GraphST", graphst_pkg)
    monkeypatch.setitem(sys.modules, "GraphST.GraphST", graphst_sub)
    monkeypatch.setitem(sys.modules, "GraphST.utils", graphst_utils)

    monkeypatch.setattr(sd, "require", _required_dependency)
    monkeypatch.setattr(
        "chatspatial.utils.compute.gmm_clustering",
        lambda data, n_clusters, **_kwargs: np.arange(data.shape[0]) % n_clusters,
    )

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    labels, _, stats = await sd._identify_domains_graphst(
        adata,
        SpatialDomainParameters(
            method="graphst",
            graphst_clustering_method="mclust",
            graphst_refinement=True,
            graphst_radius=66,
            n_domains=3,
            graphst_n_clusters=3,
        ),
        DummyCtx(adata),
    )

    assert len(labels) == adata.n_obs
    assert all(label.startswith("m_") for label in labels[:5])
    assert stats["refinement_radius"] == 66


@pytest.mark.asyncio
async def test_identify_domains_graphst_leiden_early_stopping_branch(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: x
    monkeypatch.setitem(sys.modules, "torch", torch_mod)

    class _FakeGraphST:
        def __init__(self, adata_graphst, device=None, random_seed=0):
            del device, random_seed
            self._adata = adata_graphst

        def train(self):
            rng = np.random.default_rng(13)
            self._adata.obsm["emb"] = rng.normal(size=(self._adata.n_obs, 24))
            return self._adata

    graphst_sub = types.ModuleType("GraphST.GraphST")
    graphst_sub.GraphST = _FakeGraphST
    graphst_pkg = types.ModuleType("GraphST")
    graphst_pkg.GraphST = graphst_sub
    monkeypatch.setitem(sys.modules, "GraphST", graphst_pkg)
    monkeypatch.setitem(sys.modules, "GraphST.GraphST", graphst_sub)

    monkeypatch.setattr(sd, "require", _required_dependency)
    monkeypatch.setattr(sd.sc.pp, "neighbors", lambda *_a, **_k: None)

    def _fake_leiden(a, resolution=1.0, random_state=0):
        del resolution, random_state
        a.obs["leiden"] = [str(i % 2) for i in range(a.n_obs)]

    monkeypatch.setattr(sd.sc.tl, "leiden", _fake_leiden)

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    labels, _, stats = await sd._identify_domains_graphst(
        adata,
        SpatialDomainParameters(
            method="graphst",
            graphst_clustering_method="leiden",
            graphst_refinement=False,
            n_domains=9,  # Unreachable by fake leiden => triggers early-stop branch
            graphst_n_clusters=9,
        ),
        DummyCtx(adata),
    )

    assert len(labels) == adata.n_obs
    assert stats["clustering_method"] == "leiden"


@pytest.mark.asyncio
async def test_identify_domains_graphst_generic_error_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import sys
    import types

    adata = minimal_spatial_adata.copy()

    torch_mod = types.ModuleType("torch")
    torch_mod.__version__ = "2.8.0"
    torch_mod.device = lambda x: x
    monkeypatch.setitem(sys.modules, "torch", torch_mod)

    class _FakeGraphST:
        def __init__(self, adata_graphst, device=None, random_seed=0):
            del adata_graphst, device, random_seed

        def train(self):
            raise RuntimeError("train failed")

    graphst_sub = types.ModuleType("GraphST.GraphST")
    graphst_sub.GraphST = _FakeGraphST
    graphst_pkg = types.ModuleType("GraphST")
    graphst_pkg.GraphST = graphst_sub
    monkeypatch.setitem(sys.modules, "GraphST", graphst_pkg)
    monkeypatch.setitem(sys.modules, "GraphST.GraphST", graphst_sub)

    monkeypatch.setattr(sd, "require", _required_dependency)

    async def _fake_resolve_device(prefer_gpu: bool, ctx):
        del prefer_gpu, ctx
        return "cpu"

    monkeypatch.setattr(sd, "resolve_device_async", _fake_resolve_device)

    with pytest.raises(ProcessingError, match="GraphST execution failed"):
        await sd._identify_domains_graphst(
            adata,
            SpatialDomainParameters(method="graphst"),
            DummyCtx(adata),
        )
