"""Unit tests for lightweight spatial_statistics helper contracts."""

from __future__ import annotations

import importlib
from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import SpatialStatisticsParameters
from chatspatial.tools import spatial_statistics as ss
from chatspatial.tools.spatial_statistics import (
    _build_results_keys,
    _extract_result_summary,
)
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


class DummyCtx:
    def __init__(self, adata=None, error: Exception | None = None):
        self.adata = adata
        self.error = error

    async def get_adata(self, _data_id: str):
        if self.error is not None:
            raise self.error
        return self.adata

    async def info(self, _msg: str):
        return None


def _required_dependency(name: str, *_args, **_kwargs):
    return __import__("sys").modules.get(name) or importlib.import_module(name)


def _patch_required_modules(
    monkeypatch: pytest.MonkeyPatch,
    modules: dict[str, object],
) -> None:
    def _require_module(
        _dependency_name: str,
        module_name: str,
        *_args,
        **_kwargs,
    ):
        try:
            return modules[module_name]
        except KeyError as exc:
            raise AssertionError(f"Unexpected module request: {module_name}") from exc

    monkeypatch.setattr(ss, "require_module", _require_module)


def test_build_results_keys_neighborhood_uses_cluster_dynamic_uns_key():
    keys = _build_results_keys("neighborhood", genes=None, cluster_key="leiden")
    assert keys["uns"] == ["leiden_nhood_enrichment"]
    assert keys["obs"] == []


def test_build_results_keys_local_moran_expands_gene_specific_obs_keys():
    keys = _build_results_keys("local_moran", genes=["gene_a", "gene_b"])
    assert "gene_a_local_morans" in keys["obs"]
    assert "gene_a_lisa_cluster" in keys["obs"]
    assert "gene_b_lisa_pvalue" in keys["obs"]


def test_extract_result_summary_moran_compacts_expected_fields():
    summary = _extract_result_summary(
        {
            "n_genes_analyzed": 12,
            "n_significant": 3,
            "top_highest_autocorrelation": ["g1", "g2", "g3"],
            "mean_morans_i": 0.31,
            "analysis_key": "moranI",
        },
        "moran",
    )
    assert summary["n_features_analyzed"] == 12
    assert summary["n_significant"] == 3
    assert summary["top_features"] == ["g1", "g2", "g3"]
    assert summary["results_key"] == "moranI"


def test_extract_result_summary_local_join_count_aggregates_per_category_stats():
    summary = _extract_result_summary(
        {
            "n_categories": 2,
            "categories": ["A", "B"],
            "per_category_stats": {
                "A": {"n_significant": 4, "n_hotspots": 2},
                "B": {"n_significant": 1, "n_hotspots": 1},
            },
        },
        "local_join_count",
    )
    assert summary["n_features_analyzed"] == 2
    assert summary["n_significant"] == 5
    assert summary["top_features"] == ["A", "B"]
    assert summary["summary_metrics"]["total_significant_clusters"] == 5


def test_build_results_keys_unknown_analysis_returns_empty_key_structure():
    keys = _build_results_keys("unknown_analysis", genes=["g1"], cluster_key="x")
    assert keys == {"obs": [], "var": [], "obsm": [], "uns": []}


def test_extract_result_summary_getis_ord_uses_corrected_hot_and_cold_counts():
    summary = _extract_result_summary(
        {
            "genes_analyzed": ["GeneA", "GeneB"],
            "parameters": {"correction": "fdr_bh"},
            "results": {
                "GeneA": {
                    "n_hot_spots": 50,
                    "n_cold_spots": 60,
                    "n_hot_spots_corrected": 3,
                    "n_cold_spots_corrected": 4,
                    "n_significant_corrected": 7,
                },
                "GeneB": {
                    "n_hot_spots": 70,
                    "n_cold_spots": 80,
                    "n_hot_spots_corrected": 5,
                    "n_cold_spots_corrected": 6,
                    "n_significant_corrected": 11,
                },
            },
        },
        "getis_ord",
    )

    assert summary["n_features_analyzed"] == 2
    assert summary["n_significant"] == 18
    assert summary["summary_metrics"] == {
        "total_hotspots": 8,
        "total_coldspots": 10,
    }


def test_extract_result_summary_join_count_and_network_properties_branches():
    join_summary = _extract_result_summary(
        {"bb": 10, "ww": 20, "bw": 5, "J": 35, "p_value": 0.03},
        "join_count",
    )
    assert join_summary["n_features_analyzed"] == 2
    assert join_summary["n_significant"] == 1
    assert join_summary["summary_metrics"]["total_joins"] == 35
    assert join_summary["summary_metrics"]["p_value"] == 0.03

    no_pvalue_summary = _extract_result_summary(
        {"bb": 10, "ww": 20, "bw": 5, "J": 35, "p_value": None},
        "join_count",
    )
    assert no_pvalue_summary["n_significant"] == 0
    assert "p_value" not in no_pvalue_summary["summary_metrics"]

    net_summary = _extract_result_summary(
        {"analysis_key": "network_properties", "density": 0.12, "n_cells": 60},
        "network_properties",
    )
    assert net_summary["results_key"] == "network_properties"
    assert net_summary["summary_metrics"] == {"density": 0.12}


def test_extract_result_summary_spatial_centrality_reports_cluster_and_global_metrics():
    summary = _extract_result_summary(
        {
            "analysis_key": "spatial_centrality",
            "centrality_computed": True,
            "cluster_centrality": {
                "A": {
                    "mean_degree": 0.2,
                    "mean_closeness": 0.3,
                    "mean_betweenness": 0.01,
                },
                "B": {
                    "mean_degree": 0.8,
                    "mean_closeness": 0.4,
                    "mean_betweenness": 0.02,
                },
            },
            "global_stats": {
                "mean_degree": 0.5,
                "mean_closeness": 0.35,
                "mean_betweenness": 0.015,
            },
        },
        "spatial_centrality",
    )

    assert summary["results_key"] == "spatial_centrality"
    assert summary["n_features_analyzed"] == 2
    assert summary["top_features"] == ["B", "A"]
    assert summary["summary_metrics"]["n_centrality_metrics"] == 3.0
    assert summary["summary_metrics"]["top_centrality_score"] == 0.8
    assert summary["summary_metrics"]["top_degree"] == 0.8
    assert summary["summary_metrics"]["average_degree"] == 0.5
    assert summary["summary_metrics"]["global_mean_degree"] == 0.5
    assert "centrality_computed" not in summary["summary_metrics"]


def test_get_optimal_n_jobs_respects_requested_and_auto_rules(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr("os.cpu_count", lambda: 8)
    assert ss._get_optimal_n_jobs(500, requested_n_jobs=-1) == 8
    assert ss._get_optimal_n_jobs(500, requested_n_jobs=3) == 3
    assert ss._get_optimal_n_jobs(500, requested_n_jobs=None) == 1
    assert ss._get_optimal_n_jobs(3000, requested_n_jobs=None) == 2
    assert ss._get_optimal_n_jobs(9000, requested_n_jobs=None) == 4


def test_registry_handler_direct_call_with_unified_signature(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """Registry handlers all use unified (adata, params, ctx) signature."""
    params = SpatialStatisticsParameters()
    ctx = DummyCtx(minimal_spatial_adata)

    # Create a test handler matching the unified signature
    def test_handler(adata, p, c):
        return {"n_obs": adata.n_obs, "n_neighbors": p.n_neighbors}

    config = ss._AnalysisConfig(
        handler=test_handler, needs_cluster=False, metadata_keys={}
    )
    monkeypatch.setitem(ss._ANALYSIS_REGISTRY, "_test", config)

    # Dispatch via registry (same as main function does)
    result = ss._ANALYSIS_REGISTRY["_test"].handler(minimal_spatial_adata, params, ctx)
    assert result["n_obs"] == 60
    assert result["n_neighbors"] == params.n_neighbors


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_rejects_invalid_analysis_type():
    params = SpatialStatisticsParameters().model_copy(
        update={"analysis_type": "invalid"}
    )
    with pytest.raises(ParameterError, match="Unsupported analysis type"):
        await ss.analyze_spatial_statistics("d1", DummyCtx(None), params)


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_rejects_non_positive_neighbors():
    params = SpatialStatisticsParameters().model_copy(update={"n_neighbors": 0})
    with pytest.raises(ParameterError, match="n_neighbors must be positive"):
        await ss.analyze_spatial_statistics("d1", DummyCtx(None), params)


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_requires_cluster_key_for_cluster_analysis(
    minimal_spatial_adata,
):
    params = SpatialStatisticsParameters(analysis_type="neighborhood", cluster_key=None)
    with pytest.raises(ParameterError, match="cluster_key is required"):
        await ss.analyze_spatial_statistics(
            "d1", DummyCtx(minimal_spatial_adata.copy()), params
        )


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_passes_through_data_not_found():
    params = SpatialStatisticsParameters(analysis_type="moran")
    with pytest.raises(DataNotFoundError, match="missing dataset"):
        await ss.analyze_spatial_statistics(
            "missing",
            DummyCtx(None, error=DataNotFoundError("missing dataset")),
            params,
        )


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_moran_success_path_updates_metadata_and_summary(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    params = SpatialStatisticsParameters(
        analysis_type="moran", genes=["gene_0"], n_neighbors=7
    )
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        ss,
        "ensure_spatial_neighbors",
        lambda _a, n_neighs=6, **_kw: captured.setdefault("n_neighs", n_neighs),
    )

    fake_result = {
        "n_genes_analyzed": 1,
        "n_significant": 1,
        "top_highest_autocorrelation": ["gene_0"],
        "mean_morans_i": 0.42,
        "analysis_key": "moranI",
    }
    monkeypatch.setitem(
        ss._ANALYSIS_REGISTRY,
        "moran",
        ss._AnalysisConfig(
            handler=lambda *_a, **_kw: fake_result,
            needs_cluster=False,
            metadata_keys={"uns": ["moranI"]},
        ),
    )

    def _store(
        _adata,
        *,
        analysis_name,
        method,
        parameters,
        results_keys,
        statistics,
        **_kwargs,
    ):
        captured["analysis_name"] = analysis_name
        captured["method"] = method
        captured["parameters"] = parameters
        captured["results_keys"] = results_keys
        captured["statistics"] = statistics

    monkeypatch.setattr(ss, "store_analysis_metadata", _store)
    monkeypatch.setattr(ss, "export_analysis_result", lambda *_args, **_kwargs: [])

    out = await ss.analyze_spatial_statistics("d1", DummyCtx(adata), params)
    assert out.analysis_type == "moran"
    assert out.n_features_analyzed == 1
    assert out.n_significant == 1
    assert out.top_features == ["gene_0"]
    assert out.summary_metrics["mean_morans_i"] == 0.42
    assert out.results_key == "moranI"
    assert captured["n_neighs"] == 7
    assert captured["analysis_name"] == "spatial_stats_moran"
    assert captured["method"] == "moran"
    assert captured["parameters"] == {
        "n_neighbors": 7,
        "genes": ["gene_0"],
        "n_perms": 10,
    }
    assert captured["results_keys"] == {
        "obs": [],
        "var": [],
        "obsm": [],
        "uns": ["moranI"],
    }
    assert captured["statistics"] == {"n_cells": 60, "n_significant": 1}


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_accepts_result_with_dict_method(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["0"] * 30 + ["1"] * 30)
    params = SpatialStatisticsParameters(
        analysis_type="neighborhood", cluster_key="leiden"
    )
    monkeypatch.setattr(ss, "ensure_spatial_neighbors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(ss, "store_analysis_metadata", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(ss, "export_analysis_result", lambda *_args, **_kwargs: [])
    monkeypatch.setitem(
        ss._ANALYSIS_REGISTRY,
        "neighborhood",
        ss._AnalysisConfig(
            handler=lambda *_a, **_kw: SimpleNamespace(
                dict=lambda: {
                    "n_clusters": 2,
                    "max_enrichment": 1.2,
                    "min_enrichment": -0.7,
                    "n_enriched_pairs": 3,
                    "n_depleted_pairs": 1,
                    "top_enriched_pairs": ["0-1", "1-1", "0-0"],
                    "z_threshold": 1.96,
                    "analysis_key": "leiden_nhood_enrichment",
                }
            ),
            needs_cluster=True,
            metadata_keys={"uns": ["neighborhood"]},
        ),
    )

    out = await ss.analyze_spatial_statistics("d1", DummyCtx(adata), params)
    assert out.analysis_type == "neighborhood"
    assert out.n_features_analyzed == 2
    assert out.n_significant == 3
    assert out.top_features == ["0-1", "1-1", "0-0"]
    assert out.summary_metrics == {
        "max_enrichment": 1.2,
        "min_enrichment": -0.7,
        "z_threshold": 1.96,
        "n_depleted_pairs": 1,
    }
    assert out.results_key == "leiden_nhood_enrichment"


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_passes_alternative_spatial_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """ensure_spatial_neighbors receives detected key when obsm uses X_spatial."""
    adata = minimal_spatial_adata.copy()
    coords = adata.obsm.pop("spatial")
    adata.obsm["X_spatial"] = coords

    captured: dict[str, object] = {}

    def _capture_key(_a, n_neighs=6, spatial_key="spatial"):
        captured["spatial_key"] = spatial_key

    params = SpatialStatisticsParameters(analysis_type="moran")
    monkeypatch.setattr(ss, "ensure_spatial_neighbors", _capture_key)
    monkeypatch.setattr(ss, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(ss, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setitem(
        ss._ANALYSIS_REGISTRY,
        "moran",
        ss._AnalysisConfig(
            handler=lambda *_a, **_kw: {
                "n_genes_analyzed": 0,
                "n_significant": 0,
                "top_highest_autocorrelation": [],
                "mean_morans_i": 0.0,
                "analysis_key": "moranI",
            },
            needs_cluster=False,
            metadata_keys={"uns": ["moranI"]},
        ),
    )

    await ss.analyze_spatial_statistics("d1", DummyCtx(adata), params)
    assert captured["spatial_key"] == "X_spatial"


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_wraps_invalid_result_format(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    params = SpatialStatisticsParameters(analysis_type="moran")
    monkeypatch.setattr(ss, "ensure_spatial_neighbors", lambda *_args, **_kwargs: None)
    monkeypatch.setitem(
        ss._ANALYSIS_REGISTRY,
        "moran",
        ss._AnalysisConfig(
            handler=lambda *_a, **_kw: object(),
            needs_cluster=False,
            metadata_keys={"uns": ["moranI"]},
        ),
    )

    with pytest.raises(ProcessingError, match="Invalid result format"):
        await ss.analyze_spatial_statistics("d1", DummyCtx(adata), params)


def test_analyze_join_count_preserves_parameter_error_for_non_binary_clusters(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cluster3"] = pd.Categorical(["a"] * 20 + ["b"] * 20 + ["c"] * 20)

    monkeypatch.setattr(
        ss,
        "require_module",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("Dependencies should not load before category validation")
        ),
    )

    with pytest.raises(ParameterError, match="requires binary data"):
        ss._analyze_join_count(
            adata,
            SpatialStatisticsParameters(
                analysis_type="join_count",
                cluster_key="cluster3",
                n_neighbors=5,
            ),
            DummyCtx(adata),
        )


def test_analyze_neighborhood_enrichment_uses_nan_safe_extrema(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["0"] * 30 + ["1"] * 30)

    captured: dict[str, object] = {}

    def _fake_nhood(a, cluster_key, **kwargs):
        captured["cluster_key"] = cluster_key
        captured.update(kwargs)
        a.uns["leiden_nhood_enrichment"] = {
            "zscore": np.array([[2.5, 3.1], [3.1, -2.2]])
        }

    monkeypatch.setattr(ss.sq.gr, "nhood_enrichment", _fake_nhood)

    out = ss._analyze_neighborhood_enrichment(
        adata,
        SpatialStatisticsParameters(
            analysis_type="neighborhood",
            cluster_key="leiden",
            neighborhood_n_perms=25,
            n_jobs=1,
            backend="threading",
        ),
        DummyCtx(adata),
    )
    assert captured == {
        "cluster_key": "leiden",
        "n_perms": 25,
        "n_jobs": 1,
        "backend": "threading",
        "seed": 0,
        "show_progress_bar": False,
    }
    assert out["n_clusters"] == 2
    assert out["n_enriched_pairs"] == 2
    assert out["n_depleted_pairs"] == 1
    assert out["top_enriched_pairs"] == ["0-1", "0-0"]
    assert out["z_threshold"] == 1.96
    assert out["max_enrichment"] == 3.1
    assert out["min_enrichment"] == -2.2
    assert out["analysis_key"] == "leiden_nhood_enrichment"


def test_analyze_co_occurrence_uses_interval_from_params_and_returns_distance_range(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["0"] * 30 + ["1"] * 30)

    captured: dict[str, object] = {}

    def _fake_co_occurrence(a, cluster_key, interval):
        captured["interval"] = interval
        captured["cluster_key"] = cluster_key
        occ = np.zeros((2, 2, 4))
        occ[0, 1, 2] = 3.5
        occ[1, 0, 1] = 2.0
        a.uns["leiden_co_occurrence"] = {
            "occ": occ,
            # bin edges: 5 edges → 4 intervals (matches occ shape)
            "interval": np.array([0.0, 5.0, 10.0, 15.0, 20.0]),
        }

    monkeypatch.setattr(ss.sq.gr, "co_occurrence", _fake_co_occurrence)

    params = SpatialStatisticsParameters(
        analysis_type="co_occurrence",
        cluster_key="leiden",
        co_occurrence_interval=42,
    )
    out = ss._analyze_co_occurrence(adata, params, DummyCtx(adata))

    assert captured == {"interval": 42, "cluster_key": "leiden"}
    assert out["n_clusters"] == 2
    assert out["n_intervals"] == 4  # 5 bin edges → 4 intervals
    assert out["distance_range"] == (0.0, 20.0)
    assert out["top_features"] == ["0-1", "1-0"]
    assert out["summary_metrics"]["peak_co_occurrence"] == 3.5
    assert out["summary_metrics"]["peak_interval_index"] == 2.0
    assert out["summary_metrics"]["distance_min"] == 0.0
    assert out["summary_metrics"]["distance_max"] == 20.0
    assert out["summary_metrics"]["peak_distance_midpoint"] == 12.5


def test_analyze_co_occurrence_uses_default_interval_when_not_set(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["0"] * 30 + ["1"] * 30)

    captured: dict[str, object] = {}

    def _fake_co_occurrence(a, cluster_key, interval):
        del cluster_key
        captured["interval"] = interval
        occ = np.zeros((2, 2, 3))
        occ[1, 0, 2] = 1.7
        a.uns["leiden_co_occurrence"] = {"occ": occ}

    monkeypatch.setattr(ss.sq.gr, "co_occurrence", _fake_co_occurrence)

    params = SpatialStatisticsParameters(
        analysis_type="co_occurrence",
        cluster_key="leiden",
        co_occurrence_interval=None,
    )
    out = ss._analyze_co_occurrence(adata, params, DummyCtx(adata))

    assert captured["interval"] == 50
    assert out["n_intervals"] == 50
    assert out["top_features"] == ["1-0", "0-1"]
    assert out["summary_metrics"]["n_intervals"] == 3.0
    assert out["summary_metrics"]["peak_co_occurrence"] == 1.7
    assert "distance_min" not in out["summary_metrics"]


def test_analyze_gearys_c_raises_when_results_missing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    monkeypatch.setattr(
        ss, "select_genes_for_analysis", lambda *_a, **_k: ["gene_0", "gene_1"]
    )
    monkeypatch.setattr(ss.sq.gr, "spatial_autocorr", lambda *_a, **_k: None)

    with pytest.raises(
        ProcessingError, match="Geary's C computation did not produce results"
    ):
        ss._analyze_gearys_c(
            adata,
            SpatialStatisticsParameters(analysis_type="geary", genes=["gene_0"]),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_wraps_unexpected_get_adata_errors():
    params = SpatialStatisticsParameters(analysis_type="moran")
    with pytest.raises(ProcessingError, match="Error in moran analysis"):
        await ss.analyze_spatial_statistics(
            "d-bad",
            DummyCtx(None, error=RuntimeError("boom")),
            params,
        )


def test_analyze_getis_ord_fdr_branch_adds_corrected_pvalues(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    monkeypatch.setattr(
        ss,
        "select_genes_for_analysis",
        lambda *_args, **_kwargs: ["gene_0", "gene_1"],
    )

    class _FakeGLocal:
        def __init__(self, y, w, transform="R", star=True):
            del w, transform, star
            y = np.asarray(y)
            self.Zs = np.where(y >= np.median(y), 2.5, -2.5)
            self.p_sim = np.full(y.shape[0], 0.01, dtype=float)

    fake_esda_getisord = ModuleType("esda.getisord")
    fake_esda_getisord.G_Local = _FakeGLocal

    class _KNNObj:
        def __init__(self):
            self.transform = None

    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.KNN = SimpleNamespace(
        from_array=lambda *_args, **_kwargs: _KNNObj()
    )
    _patch_required_modules(
        monkeypatch,
        {
            "esda.getisord": fake_esda_getisord,
            "libpysal.weights": fake_libpysal_weights,
        },
    )

    params = SpatialStatisticsParameters(
        analysis_type="getis_ord",
        genes=["gene_0", "gene_1"],
        n_neighbors=5,
        getis_ord_correction="fdr_bh",
        getis_ord_alpha=0.05,
    )
    out = ss._analyze_getis_ord(adata, params, DummyCtx(adata))

    assert out["n_genes_analyzed"] == 2
    assert "gene_0" in out["results"]
    assert "n_significant_corrected" in out["results"]["gene_0"]
    assert "gene_0_getis_ord_p_corrected" in adata.obs.columns
    assert "gene_1_getis_ord_p_corrected" in adata.obs.columns


def test_analyze_bivariate_moran_computes_and_persists_results(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    from scipy.sparse import csr_matrix

    adata = minimal_spatial_adata.copy()
    n = adata.n_obs
    dense_w = np.ones((n, n), dtype=float)
    np.fill_diagonal(dense_w, 0.0)

    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.KNN = SimpleNamespace(
        from_array=lambda *_args, **_kwargs: SimpleNamespace(
            sparse=csr_matrix(dense_w), transform=None
        )
    )
    _patch_required_modules(monkeypatch, {"libpysal.weights": fake_libpysal_weights})

    params = SpatialStatisticsParameters(
        analysis_type="bivariate_moran",
        gene_pairs=[("gene_0", "gene_1"), ("gene_0", "not_exist")],
        n_neighbors=4,
    )
    out = ss._analyze_bivariate_moran(adata, params, DummyCtx(adata))

    assert out["n_pairs_analyzed"] == 1
    assert "gene_0_vs_gene_1" in out["bivariate_morans_i"]
    assert "bivariate_moran" in adata.uns


def test_analyze_network_properties_disconnected_component_branch(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    from scipy.sparse import csr_matrix

    adata = minimal_spatial_adata.copy()
    n = adata.n_obs

    # Two disconnected components: first 30, last 30
    conn = np.zeros((n, n), dtype=float)
    for i in range(29):
        conn[i, i + 1] = 1.0
        conn[i + 1, i] = 1.0
    for i in range(30, 59):
        conn[i, i + 1] = 1.0
        conn[i + 1, i] = 1.0
    adata.obsp["spatial_connectivities"] = csr_matrix(conn)

    monkeypatch.setattr(ss, "require", _required_dependency)

    out = ss._analyze_network_properties(
        adata,
        SpatialStatisticsParameters(
            analysis_type="network_properties",
            cluster_key="group",
            n_neighbors=4,
        ),
        DummyCtx(adata),
    )

    assert out["is_connected"] is False
    assert out["n_components"] == 2
    assert out["largest_component_size"] == 30
    assert "network_properties" in adata.uns


def test_analyze_spatial_centrality_handles_missing_node_keys(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    from scipy.sparse import eye

    adata = minimal_spatial_adata.copy()
    adata.obs["group"] = pd.Categorical(["a"] * 30 + ["b"] * 30)
    adata.obsp["spatial_connectivities"] = eye(adata.n_obs, format="csr")

    fake_nx = ModuleType("networkx")
    fake_nx.from_scipy_sparse_array = lambda _arr: object()
    fake_nx.degree_centrality = lambda _g: {0: 1.0}  # missing all other nodes
    fake_nx.closeness_centrality = lambda _g: {0: 0.5}
    fake_nx.betweenness_centrality = lambda _g: {0: 0.2}

    monkeypatch.setitem(__import__("sys").modules, "networkx", fake_nx)
    monkeypatch.setattr(ss, "require", _required_dependency)

    with pytest.warns(UserWarning, match="Centrality computation incomplete"):
        out = ss._analyze_spatial_centrality(
            adata,
            SpatialStatisticsParameters(
                analysis_type="spatial_centrality",
                cluster_key="group",
                n_neighbors=4,
            ),
            DummyCtx(adata),
        )

    assert out["analysis_key"] == "spatial_centrality"
    assert out["centrality_computed"] is True
    assert adata.uns["spatial_centrality"] is out
    assert "degree_centrality" in adata.obs.columns
    assert adata.obs["degree_centrality"].shape[0] == adata.n_obs
    # Only node 0 from fake centrality has non-zero; others should be safely backfilled with 0
    assert float(adata.obs["degree_centrality"].iloc[0]) == 1.0
    assert float(adata.obs["degree_centrality"].iloc[1]) == 0.0


def test_build_results_keys_dynamic_cluster_variants_and_getis_ord_obs_keys():
    co_keys = _build_results_keys("co_occurrence", genes=None, cluster_key="leiden")
    ripley_keys = _build_results_keys("ripley", genes=None, cluster_key="celltype")
    cent_keys = _build_results_keys("centrality", genes=None, cluster_key="group")
    getis_keys = _build_results_keys("getis_ord", genes=["g1", "g2"])
    spatial_centrality_keys = _build_results_keys(
        "spatial_centrality", genes=None, cluster_key="group"
    )

    assert co_keys["uns"] == ["leiden_co_occurrence"]
    assert ripley_keys["uns"] == ["celltype_ripley_L"]
    assert cent_keys["uns"] == ["group_centrality_scores"]
    assert set(getis_keys["obs"]) == {
        "g1_getis_ord_z",
        "g1_getis_ord_p",
        "g1_getis_ord_p_corrected",
        "g2_getis_ord_z",
        "g2_getis_ord_p",
        "g2_getis_ord_p_corrected",
    }
    assert spatial_centrality_keys["uns"] == ["spatial_centrality"]


def test_extract_result_summary_covers_remaining_analysis_branches():
    geary = _extract_result_summary(
        {
            "n_genes_analyzed": 5,
            "n_significant": 2,
            "top_positive_autocorrelation": ["g1", "g2", "g3"],
            "mean_gearys_c": 0.88,
            "min_gearys_c": 0.2,
            "max_gearys_c": 1.4,
            "analysis_key": "gearyC",
        },
        "geary",
    )
    assert geary["n_features_analyzed"] == 5
    assert geary["n_significant"] == 2
    assert geary["top_features"] == ["g1", "g2", "g3"]
    assert geary["summary_metrics"] == {
        "mean_gearys_c": 0.88,
        "min_gearys_c": 0.2,
        "max_gearys_c": 1.4,
    }
    assert geary["results_key"] == "gearyC"

    local_moran = _extract_result_summary(
        {
            "genes_analyzed": ["g1", "g2"],
            "results": {
                "g1": {"n_significant": 2, "n_hotspots": 1, "n_coldspots": 3},
                "g2": {"n_significant": 4, "n_hotspots": 5, "n_coldspots": 1},
            },
        },
        "local_moran",
    )
    assert local_moran["n_features_analyzed"] == 2
    assert local_moran["n_significant"] == 6
    assert local_moran["summary_metrics"]["mean_hotspots_per_gene"] == 3.0
    assert local_moran["summary_metrics"]["mean_coldspots_per_gene"] == 2.0

    getis = _extract_result_summary(
        {
            "genes_analyzed": ["g1"],
            "results": {"g1": {"n_hot_spots": 7, "n_cold_spots": 2}},
        },
        "getis_ord",
    )
    assert getis["n_features_analyzed"] == 1
    assert getis["summary_metrics"] == {"total_hotspots": 7, "total_coldspots": 2}

    co = _extract_result_summary(
        {
            "n_clusters": 3,
            "analysis_key": "x_co",
            "top_features": ["0-1"],
            "summary_metrics": {"peak_co_occurrence": 2.5},
        },
        "co_occurrence",
    )
    rip = _extract_result_summary({"n_clusters": 4, "analysis_key": "x_rip"}, "ripley")
    cen = _extract_result_summary(
        {
            "n_clusters": 2,
            "analysis_key": "x_cent",
            "top_features": ["1"],
            "summary_metrics": {"top_centrality_score": 0.9},
        },
        "centrality",
    )
    assert co["n_features_analyzed"] == 3 and co["results_key"] == "x_co"
    assert co["top_features"] == ["0-1"]
    assert co["summary_metrics"] == {"peak_co_occurrence": 2.5}
    assert rip["n_features_analyzed"] == 4 and rip["results_key"] == "x_rip"
    assert cen["n_features_analyzed"] == 2 and cen["results_key"] == "x_cent"
    assert cen["top_features"] == ["1"]
    assert cen["summary_metrics"] == {"top_centrality_score": 0.9}

    biv = _extract_result_summary(
        {
            "n_pairs_analyzed": 3,
            "bivariate_morans_i": {"a_vs_b": 0.1, "c_vs_d": 0.4, "e_vs_f": -0.35},
            "mean_bivariate_i": 0.05,
        },
        "bivariate_moran",
    )
    assert biv["n_features_analyzed"] == 3
    assert biv["n_significant"] == 2
    assert biv["top_features"] == ["c_vs_d", "e_vs_f"]
    assert biv["summary_metrics"]["mean_bivariate_i"] == 0.05
    assert biv["summary_metrics"]["max_abs_bivariate_i"] == 0.4
    assert biv["summary_metrics"]["notable_abs_i_threshold"] == 0.3

    ljc_empty = _extract_result_summary(
        {"n_categories": 0, "categories": [], "per_category_stats": {}},
        "local_join_count",
    )
    assert ljc_empty["summary_metrics"] == {
        "mean_hotspots_per_category": 0.0,
        "total_significant_clusters": 0,
    }


def test_extract_result_summary_bivariate_moran_reports_strongest_non_notable_pair():
    summary = _extract_result_summary(
        {
            "n_pairs_analyzed": 2,
            "bivariate_morans_i": {"a_vs_b": 0.1, "c_vs_d": -0.25},
            "mean_bivariate_i": -0.075,
        },
        "bivariate_moran",
    )

    assert summary["n_significant"] == 0
    assert summary["top_features"] == []
    assert summary["summary_metrics"]["max_abs_bivariate_i"] == 0.25
    assert summary["summary_metrics"]["max_bivariate_i"] == 0.1
    assert summary["summary_metrics"]["min_bivariate_i"] == -0.25


def test_extract_result_summary_local_moran_empty_results_uses_zero_defaults():
    summary = _extract_result_summary(
        {"genes_analyzed": [], "results": {}},
        "local_moran",
    )
    assert summary["n_features_analyzed"] == 0
    assert summary["n_significant"] == 0
    assert summary["summary_metrics"] == {
        "mean_hotspots_per_gene": 0.0,
        "mean_coldspots_per_gene": 0.0,
    }


def test_analyze_morans_i_handles_small_gene_set_without_duplicate_rank_lists(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(
        ss,
        "select_genes_for_analysis",
        lambda *_args, **_kwargs: ["gene_0", "gene_1", "gene_2", "gene_3"],
    )

    def _fake_autocorr(a, **_kwargs):
        a.uns["moranI"] = pd.DataFrame(
            # p-values chosen so 2 genes survive FDR correction (BH):
            # sorted [0.001, 0.02, 0.2, 0.9] → FDR [0.004, 0.04, 0.267, 0.9]
            {"I": [0.7, 0.2, -0.1, 0.4], "pval_norm": [0.001, 0.2, 0.02, 0.9]},
            index=["gene_0", "gene_1", "gene_2", "gene_3"],
        )

    monkeypatch.setattr(ss.sq.gr, "spatial_autocorr", _fake_autocorr)
    out = ss._analyze_morans_i(
        adata,
        SpatialStatisticsParameters(analysis_type="moran", n_top_genes=4),
        DummyCtx(adata),
    )
    assert out["n_genes_analyzed"] == 4
    assert out["n_significant"] == 2  # gene_0 and gene_2 survive FDR
    assert out["top_highest_autocorrelation"] == []
    assert out["top_lowest_autocorrelation"] == []
    assert out["analysis_key"] == "moranI"


def test_analyze_morans_i_raises_when_results_missing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(
        ss, "select_genes_for_analysis", lambda *_args, **_kwargs: ["gene_0"]
    )
    monkeypatch.setattr(ss.sq.gr, "spatial_autocorr", lambda *_args, **_kwargs: None)
    with pytest.raises(
        ProcessingError, match="Moran's I computation did not produce results"
    ):
        ss._analyze_morans_i(
            adata,
            SpatialStatisticsParameters(analysis_type="moran", genes=["gene_0"]),
            DummyCtx(adata),
        )


def test_analyze_gearys_c_success_dataframe_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(
        ss, "select_genes_for_analysis", lambda *_args, **_kwargs: ["gene_0", "gene_1"]
    )

    def _fake_geary(a, **_kwargs):
        a.uns["gearyC"] = pd.DataFrame(
            {"C": [0.2, 1.0], "pval_norm": [0.001, 0.2]},
            index=["gene_0", "gene_1"],
        )

    monkeypatch.setattr(ss.sq.gr, "spatial_autocorr", _fake_geary)
    out = ss._analyze_gearys_c(
        adata,
        SpatialStatisticsParameters(analysis_type="geary", genes=["gene_0", "gene_1"]),
        DummyCtx(adata),
    )
    assert out["n_genes_analyzed"] == 2
    assert out["n_significant"] == 1
    assert out["top_positive_autocorrelation"] == ["gene_0", "gene_1"]
    assert out["top_weak_autocorrelation"] == ["gene_1", "gene_0"]
    assert out["mean_gearys_c"] == pytest.approx(0.6)
    assert out["min_gearys_c"] == pytest.approx(0.2)
    assert out["max_gearys_c"] == pytest.approx(1.0)
    assert out["analysis_key"] == "gearyC"
    assert adata.uns["gearyC"].loc["gene_0", "pval_norm_fdr"] == pytest.approx(0.002)


def test_analyze_neighborhood_and_co_occurrence_raise_when_results_missing(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["0"] * 30 + ["1"] * 30)

    monkeypatch.setattr(ss.sq.gr, "nhood_enrichment", lambda *_args, **_kwargs: None)
    nhood_params = SpatialStatisticsParameters(
        analysis_type="neighborhood", cluster_key="leiden"
    )
    with pytest.raises(
        ProcessingError, match="Neighborhood enrichment did not produce results"
    ):
        ss._analyze_neighborhood_enrichment(adata, nhood_params, DummyCtx(adata))

    monkeypatch.setattr(ss.sq.gr, "co_occurrence", lambda *_args, **_kwargs: None)
    with pytest.raises(
        ProcessingError, match="Co-occurrence analysis did not produce results"
    ):
        ss._analyze_co_occurrence(
            adata,
            SpatialStatisticsParameters(
                analysis_type="co_occurrence", cluster_key="leiden"
            ),
            DummyCtx(adata),
        )


def test_analyze_ripleys_k_success_and_error_paths(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["0"] * 20 + ["1"] * 20 + ["2"] * 20)

    monkeypatch.setattr(ss.sq.gr, "ripley", lambda *_args, **_kwargs: None)
    ripley_params = SpatialStatisticsParameters(
        analysis_type="ripley", cluster_key="leiden"
    )
    out = ss._analyze_ripleys_k(adata, ripley_params, DummyCtx(adata))
    assert out["analysis_completed"] is True
    assert out["analysis_key"] == "leiden_ripley_L"
    assert out["n_clusters"] == 3

    monkeypatch.setattr(
        ss.sq.gr,
        "ripley",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("bad ripley")),
    )
    with pytest.raises(ProcessingError, match="Ripley's K analysis failed"):
        ss._analyze_ripleys_k(adata, ripley_params, DummyCtx(adata))


def test_analyze_getis_ord_bonferroni_branch_adds_corrected_columns(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    monkeypatch.setattr(
        ss,
        "select_genes_for_analysis",
        lambda *_args, **_kwargs: ["gene_0", "gene_1"],
    )

    class _FakeGLocal:
        def __init__(self, y, w, transform="R", star=True):
            del w, transform, star
            n = len(np.asarray(y))
            self.Zs = np.tile(np.array([2.2, -2.2, 0.0, 1.2]), int(np.ceil(n / 4)))[:n]
            self.p_sim = np.full(n, 0.02, dtype=float)

    fake_esda_getisord = ModuleType("esda.getisord")
    fake_esda_getisord.G_Local = _FakeGLocal

    class _KNNObj:
        def __init__(self):
            self.transform = None

    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.KNN = SimpleNamespace(
        from_array=lambda *_args, **_kwargs: _KNNObj()
    )
    _patch_required_modules(
        monkeypatch,
        {
            "esda.getisord": fake_esda_getisord,
            "libpysal.weights": fake_libpysal_weights,
        },
    )

    params = SpatialStatisticsParameters(
        analysis_type="getis_ord",
        genes=["gene_0", "gene_1"],
        n_neighbors=6,
        getis_ord_correction="bonferroni",
        getis_ord_alpha=0.05,
    )
    out = ss._analyze_getis_ord(adata, params, DummyCtx(adata))

    assert out["n_genes_analyzed"] == 2
    assert "n_hot_spots_corrected" in out["results"]["gene_0"]
    assert "n_cold_spots_corrected" in out["results"]["gene_1"]
    assert "gene_0_getis_ord_p_corrected" in adata.obs.columns
    assert "gene_1_getis_ord_p_corrected" in adata.obs.columns


def test_analyze_getis_ord_wraps_runtime_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    monkeypatch.setattr(
        ss, "select_genes_for_analysis", lambda *_args, **_kwargs: ["gene_0"]
    )

    class _BrokenGLocal:
        def __init__(self, *_args, **_kwargs):
            raise RuntimeError("cannot compute")

    fake_esda_getisord = ModuleType("esda.getisord")
    fake_esda_getisord.G_Local = _BrokenGLocal
    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.KNN = SimpleNamespace(
        from_array=lambda *_args, **_kwargs: SimpleNamespace(transform=None)
    )
    _patch_required_modules(
        monkeypatch,
        {
            "esda.getisord": fake_esda_getisord,
            "libpysal.weights": fake_libpysal_weights,
        },
    )

    with pytest.raises(ProcessingError, match="Getis-Ord analysis failed"):
        ss._analyze_getis_ord(
            adata,
            SpatialStatisticsParameters(analysis_type="getis_ord", genes=["gene_0"]),
            DummyCtx(adata),
        )


def test_analyze_centrality_success_and_missing_key_failure(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["0"] * 30 + ["1"] * 30)

    def _fake_centrality_scores(a, cluster_key):
        a.uns[f"{cluster_key}_centrality_scores"] = pd.DataFrame(
            {
                "degree_centrality": [0.2, 0.9],
                "closeness_centrality": [0.3, 0.6],
                "label": ["low", "high"],
            },
            index=["0", "1"],
        )

    monkeypatch.setattr(ss.sq.gr, "centrality_scores", _fake_centrality_scores)
    cent_params = SpatialStatisticsParameters(
        analysis_type="centrality", cluster_key="leiden"
    )
    out = ss._analyze_centrality(adata, cent_params, DummyCtx(adata))
    assert out["analysis_completed"] is True
    assert out["analysis_key"] == "leiden_centrality_scores"
    assert out["n_clusters"] == 2
    assert out["top_features"] == ["1", "0"]
    assert out["summary_metrics"]["n_centrality_metrics"] == 2.0
    assert out["summary_metrics"]["top_centrality_score"] == 0.9
    assert out["summary_metrics"]["top_degree_centrality"] == 0.9
    assert out["summary_metrics"]["average_degree_centrality"] == 0.55

    adata.uns["legacy_centrality_scores"] = {
        "0": {"degree_centrality": 0.1},
        "1": {"degree_centrality": 0.8},
    }
    legacy_summary = ss._summarize_centrality_scores(
        adata.uns["legacy_centrality_scores"]
    )
    assert legacy_summary["top_features"] == ["1", "0"]

    monkeypatch.setattr(ss.sq.gr, "centrality_scores", lambda *_args, **_kwargs: None)
    adata.uns.pop("leiden_centrality_scores", None)
    with pytest.raises(
        ProcessingError, match="Centrality analysis did not produce results"
    ):
        ss._analyze_centrality(adata, cent_params, DummyCtx(adata))


def test_analyze_bivariate_moran_requires_gene_pairs_and_handles_zero_variance(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    from scipy.sparse import csr_matrix

    adata = minimal_spatial_adata.copy()
    adata.X[:, 0] = 1.0
    adata.X[:, 1] = 1.0
    with pytest.raises(ParameterError, match="requires gene_pairs"):
        ss._analyze_bivariate_moran(
            adata,
            SpatialStatisticsParameters(
                analysis_type="bivariate_moran", gene_pairs=None
            ),
            DummyCtx(adata),
        )

    n = adata.n_obs
    dense_w = np.ones((n, n), dtype=float)
    np.fill_diagonal(dense_w, 0.0)
    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.KNN = SimpleNamespace(
        from_array=lambda *_args, **_kwargs: SimpleNamespace(
            sparse=csr_matrix(dense_w), transform=None
        )
    )
    _patch_required_modules(monkeypatch, {"libpysal.weights": fake_libpysal_weights})

    out = ss._analyze_bivariate_moran(
        adata,
        SpatialStatisticsParameters(
            analysis_type="bivariate_moran",
            gene_pairs=[("gene_0", "gene_1")],
            n_neighbors=5,
        ),
        DummyCtx(adata),
    )
    assert out["bivariate_morans_i"]["gene_0_vs_gene_1"] == 0.0


def test_analyze_bivariate_moran_wraps_internal_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    _patch_required_modules(
        monkeypatch,
        {"libpysal.weights": ModuleType("libpysal.weights")},
    )
    monkeypatch.setattr(
        ss,
        "require_spatial_coords",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("coords missing")),
    )

    with pytest.raises(
        ProcessingError, match="Bivariate Moran's I failed: coords missing"
    ):
        ss._analyze_bivariate_moran(
            adata,
            SpatialStatisticsParameters(
                analysis_type="bivariate_moran",
                gene_pairs=[("gene_0", "gene_1")],
                n_neighbors=4,
            ),
            DummyCtx(adata),
        )


@pytest.mark.parametrize(
    "analysis_type",
    [
        "getis_ord",
        "bivariate_moran",
        "join_count",
        "local_join_count",
        "network_properties",
        "spatial_centrality",
    ],
)
def test_spatial_helpers_preserve_coordinate_domain_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, analysis_type
):
    adata = minimal_spatial_adata.copy()
    adata.obs["binary"] = pd.Categorical(["A"] * 30 + ["B"] * 30)
    adata.obs["celltype"] = pd.Categorical(["A"] * 20 + ["B"] * 40)
    adata.obs["group"] = pd.Categorical(["A"] * 30 + ["B"] * 30)
    adata.obsp.pop("spatial_connectivities", None)

    monkeypatch.setattr(
        ss, "require_module", lambda *_args, **_kwargs: SimpleNamespace()
    )
    monkeypatch.setattr(ss, "require", lambda *_args, **_kwargs: SimpleNamespace())
    monkeypatch.setattr(
        ss,
        "select_genes_for_analysis",
        lambda *_args, **_kwargs: ["gene_0"],
    )
    monkeypatch.setattr(
        ss,
        "require_spatial_coords",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DataCompatibilityError("coordinates are invalid")
        ),
    )

    cases = {
        "getis_ord": (
            ss._analyze_getis_ord,
            SpatialStatisticsParameters(analysis_type="getis_ord", genes=["gene_0"]),
        ),
        "bivariate_moran": (
            ss._analyze_bivariate_moran,
            SpatialStatisticsParameters(
                analysis_type="bivariate_moran",
                gene_pairs=[("gene_0", "gene_1")],
            ),
        ),
        "join_count": (
            ss._analyze_join_count,
            SpatialStatisticsParameters(
                analysis_type="join_count", cluster_key="binary"
            ),
        ),
        "local_join_count": (
            ss._analyze_local_join_count,
            SpatialStatisticsParameters(
                analysis_type="local_join_count", cluster_key="celltype"
            ),
        ),
        "network_properties": (
            ss._analyze_network_properties,
            SpatialStatisticsParameters(
                analysis_type="network_properties", cluster_key="group"
            ),
        ),
        "spatial_centrality": (
            ss._analyze_spatial_centrality,
            SpatialStatisticsParameters(
                analysis_type="spatial_centrality", cluster_key="group"
            ),
        ),
    }
    analysis, params = cases[analysis_type]

    with pytest.raises(DataCompatibilityError, match="coordinates are invalid"):
        analysis(adata, params, DummyCtx(adata))


def test_analyze_join_count_success_and_wraps_runtime_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["binary"] = pd.Categorical(["x"] * 30 + ["y"] * 30)

    monkeypatch.setattr(ss, "require_spatial_coords", lambda a: a.obsm["spatial"])

    class _FakeJoinCounts:
        def __init__(self, y, w):
            del y, w
            self.bb = 11.0
            self.ww = 7.0
            self.bw = 5.0
            self.J = 23.0
            self.p_sim = 0.04

    fake_esda_join_counts = ModuleType("esda.join_counts")
    fake_esda_join_counts.Join_Counts = _FakeJoinCounts
    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.KNN = SimpleNamespace(
        from_array=lambda *_args, **_kwargs: object()
    )
    _patch_required_modules(
        monkeypatch,
        {
            "esda.join_counts": fake_esda_join_counts,
            "libpysal.weights": fake_libpysal_weights,
        },
    )

    jc_params = SpatialStatisticsParameters(
        analysis_type="join_count", cluster_key="binary", n_neighbors=4
    )
    out = ss._analyze_join_count(adata, jc_params, DummyCtx(adata))
    assert out == {"bb": 11.0, "ww": 7.0, "bw": 5.0, "J": 23.0, "p_value": 0.04}

    class _BrokenJoinCounts:
        def __init__(self, *_args, **_kwargs):
            raise RuntimeError("join count failed")

    fake_esda_join_counts.Join_Counts = _BrokenJoinCounts
    with pytest.raises(ProcessingError, match="Join Count analysis failed"):
        ss._analyze_join_count(adata, jc_params, DummyCtx(adata))


def test_analyze_local_join_count_success_and_failure_paths(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["celltype"] = pd.Categorical(["A"] * 20 + ["B"] * 20 + ["C"] * 20)

    monkeypatch.setattr(ss, "require_spatial_coords", lambda a: a.obsm["spatial"])

    class _FakeLJC:
        def __init__(self, connectivity):
            self.connectivity = connectivity

        def fit(self, y):
            del self.connectivity
            y = np.asarray(y)
            self.LJC = y.astype(float)
            self.p_sim = np.where(y > 0, 0.01, 0.5)
            return self

    fake_esda_ljc = ModuleType("esda.join_counts_local")
    fake_esda_ljc.Join_Counts_Local = _FakeLJC
    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.KNN = SimpleNamespace(
        from_array=lambda *_args, **_kwargs: object()
    )
    _patch_required_modules(
        monkeypatch,
        {
            "esda.join_counts_local": fake_esda_ljc,
            "libpysal.weights": fake_libpysal_weights,
        },
    )

    ljc_params = SpatialStatisticsParameters(
        analysis_type="local_join_count", cluster_key="celltype", n_neighbors=4
    )
    out = ss._analyze_local_join_count(adata, ljc_params, DummyCtx(adata))
    assert out["n_categories"] == 3
    assert "A" in out["per_category_stats"]
    assert "local_join_count" in adata.uns
    assert "ljc_A" in adata.obs.columns
    assert "ljc_B_pvalue" in adata.obs.columns

    class _BrokenLJC:
        def __init__(self, connectivity):
            self.connectivity = connectivity

        def fit(self, *_args, **_kwargs):
            raise RuntimeError("local join count failed")

    fake_esda_ljc.Join_Counts_Local = _BrokenLJC
    with pytest.raises(ProcessingError, match="Local Join Count analysis failed"):
        ss._analyze_local_join_count(adata, ljc_params, DummyCtx(adata))


def test_analyze_network_properties_connected_fallback_and_avg_clustering_guard(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import networkx as nx

    adata = minimal_spatial_adata.copy()
    adata.obsp.pop("spatial_connectivities", None)

    monkeypatch.setattr(ss, "require", lambda *_args, **_kwargs: nx)
    monkeypatch.setattr(
        nx,
        "average_clustering",
        lambda _g: (_ for _ in ()).throw(RuntimeError("clustering failed")),
    )

    out = ss._analyze_network_properties(
        adata,
        SpatialStatisticsParameters(
            analysis_type="network_properties",
            cluster_key="group",
            n_neighbors=59,
        ),
        DummyCtx(adata),
    )
    assert out["is_connected"] is True
    assert "diameter" in out and "radius" in out
    assert out["avg_clustering"] is None
    assert "network_properties" in adata.uns


def test_analyze_spatial_centrality_fallback_kneighbors_graph_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    from scipy.sparse import eye

    adata = minimal_spatial_adata.copy()
    adata.obs["group"] = pd.Categorical(["a"] * 30 + ["b"] * 30)
    adata.obsp.pop("spatial_connectivities", None)
    n_nodes = adata.n_obs

    fake_nx = ModuleType("networkx")
    fake_nx.from_scipy_sparse_array = lambda _arr: object()
    fake_nx.degree_centrality = lambda _g: {i: 0.1 for i in range(n_nodes)}
    fake_nx.closeness_centrality = lambda _g: {i: 0.2 for i in range(n_nodes)}
    fake_nx.betweenness_centrality = lambda _g: {i: 0.3 for i in range(n_nodes)}
    monkeypatch.setattr(ss, "require", lambda *_args, **_kwargs: fake_nx)

    fake_sklearn_neighbors = ModuleType("sklearn.neighbors")
    fake_sklearn_neighbors.kneighbors_graph = lambda *_args, **_kwargs: eye(
        n_nodes, format="csr"
    )
    monkeypatch.setitem(
        __import__("sys").modules, "sklearn.neighbors", fake_sklearn_neighbors
    )

    out = ss._analyze_spatial_centrality(
        adata,
        SpatialStatisticsParameters(
            analysis_type="spatial_centrality",
            cluster_key="group",
            n_neighbors=5,
        ),
        DummyCtx(adata),
    )
    assert out["analysis_key"] == "spatial_centrality"
    assert out["centrality_computed"] is True
    assert adata.uns["spatial_centrality"] is out
    assert out["global_stats"]["mean_degree"] == pytest.approx(0.1)
    assert out["global_stats"]["mean_closeness"] == pytest.approx(0.2)
    assert out["global_stats"]["mean_betweenness"] == pytest.approx(0.3)


def test_analyze_local_moran_success_with_fdr_and_persists_obs_uns(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    from scipy.sparse import csr_matrix

    adata = minimal_spatial_adata.copy()
    adata.obsp["spatial_connectivities"] = csr_matrix(np.eye(adata.n_obs))

    monkeypatch.setattr(
        ss, "select_genes_for_analysis", lambda *_args, **_kwargs: ["gene_0"]
    )

    class _FakePySALWeights:
        def __init__(self, neighbors, weights):
            self.neighbors = neighbors
            self.weights = weights

    class _FakeMoranLocal:
        def __init__(self, expr, w, permutations):
            del w, permutations
            n = len(expr)
            self.Is = np.linspace(0.0, 1.0, n)
            self.p_sim = np.tile(np.array([0.01, 0.2]), int(np.ceil(n / 2)))[:n]
            self.q = np.tile(np.array([1, 3, 4, 2]), int(np.ceil(n / 4)))[:n]

    fake_esda_moran = ModuleType("esda.moran")
    fake_esda_moran.Moran_Local = _FakeMoranLocal
    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.W = _FakePySALWeights

    _patch_required_modules(
        monkeypatch,
        {
            "esda.moran": fake_esda_moran,
            "libpysal.weights": fake_libpysal_weights,
        },
    )
    monkeypatch.setattr(
        ss,
        "adjust_pvalues",
        lambda p_values, **_kwargs: (
            np.asarray(p_values) < 0.05,
            np.asarray(p_values) * 0.5,
        ),
    )

    out = ss._analyze_local_moran(
        adata,
        SpatialStatisticsParameters(
            analysis_type="local_moran",
            genes=["gene_0"],
            local_moran_permutations=99,
            local_moran_alpha=0.05,
            local_moran_fdr_correction=True,
            n_neighbors=5,
        ),
        DummyCtx(adata),
    )

    assert out["analysis_type"] == "local_moran"
    assert out["genes_analyzed"] == ["gene_0"]
    assert "gene_0_local_morans" in adata.obs.columns
    assert "gene_0_lisa_cluster" in adata.obs.columns
    assert "gene_0_lisa_pvalue" in adata.obs.columns
    assert "local_moran" in adata.uns
    assert out["results"]["gene_0"]["fdr_corrected"] is True


def test_analyze_local_moran_dense_connectivity_and_no_fdr_branch(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsp["spatial_connectivities"] = np.eye(adata.n_obs, dtype=float)

    monkeypatch.setattr(
        ss, "select_genes_for_analysis", lambda *_args, **_kwargs: ["gene_0"]
    )

    class _FakePySALWeights:
        def __init__(self, neighbors, weights):
            self.neighbors = neighbors
            self.weights = weights

    class _FakeMoranLocal:
        def __init__(self, expr, w, permutations):
            del expr, w, permutations
            n = adata.n_obs
            self.Is = np.linspace(0.1, 1.0, n)
            self.p_sim = np.tile(np.array([0.01, 0.5]), int(np.ceil(n / 2)))[:n]
            self.q = np.tile(np.array([1, 3, 4, 2]), int(np.ceil(n / 4)))[:n]

    fake_esda_moran = ModuleType("esda.moran")
    fake_esda_moran.Moran_Local = _FakeMoranLocal
    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.W = _FakePySALWeights
    _patch_required_modules(
        monkeypatch,
        {
            "esda.moran": fake_esda_moran,
            "libpysal.weights": fake_libpysal_weights,
        },
    )

    out = ss._analyze_local_moran(
        adata,
        SpatialStatisticsParameters(
            analysis_type="local_moran",
            genes=["gene_0"],
            local_moran_permutations=25,
            local_moran_alpha=0.05,
            local_moran_fdr_correction=False,
            n_neighbors=5,
        ),
        DummyCtx(adata),
    )

    assert out["results"]["gene_0"]["fdr_corrected"] is False
    assert out["results"]["gene_0"]["n_significant"] > 0
    assert "gene_0_lisa_cluster" in adata.obs.columns


def test_analyze_local_moran_preserves_fdr_dependency_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    from scipy.sparse import csr_matrix

    adata = minimal_spatial_adata.copy()
    adata.obsp["spatial_connectivities"] = csr_matrix(np.eye(adata.n_obs))

    class _FakePySALWeights:
        def __init__(self, _neighbors, _weights):
            pass

    class _FakeMoranLocal:
        def __init__(self, expr, _weights, permutations):
            del permutations
            n = len(expr)
            self.Is = np.zeros(n)
            self.p_sim = np.full(n, 0.1)
            self.q = np.ones(n, dtype=int)

    fake_moran = ModuleType("esda.moran")
    fake_moran.Moran_Local = _FakeMoranLocal
    fake_weights = ModuleType("libpysal.weights")
    fake_weights.W = _FakePySALWeights

    def _require_module(_dependency, module_name, *_args, **_kwargs):
        if module_name == "esda.moran":
            return fake_moran
        if module_name == "libpysal.weights":
            return fake_weights
        raise AssertionError(f"Unexpected module request: {module_name}")

    def _raise_dependency_error(*_args, **_kwargs):
        raise DependencyError("statsmodels submodule is broken")

    monkeypatch.setattr(ss, "require_module", _require_module)
    monkeypatch.setattr(ss, "adjust_pvalues", _raise_dependency_error)

    with pytest.raises(DependencyError, match="statsmodels submodule is broken"):
        ss._analyze_local_moran(
            adata,
            SpatialStatisticsParameters(
                analysis_type="local_moran",
                genes=["gene_0"],
                local_moran_permutations=99,
                local_moran_fdr_correction=True,
            ),
            DummyCtx(adata),
        )


def test_analyze_local_moran_wraps_runtime_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    from scipy.sparse import csr_matrix

    adata = minimal_spatial_adata.copy()
    adata.obsp["spatial_connectivities"] = csr_matrix(np.eye(adata.n_obs))

    monkeypatch.setattr(
        ss, "select_genes_for_analysis", lambda *_args, **_kwargs: ["gene_0"]
    )

    class _FakePySALWeights:
        def __init__(self, neighbors, weights):
            self.neighbors = neighbors
            self.weights = weights

    class _BrokenMoranLocal:
        def __init__(self, *_args, **_kwargs):
            raise RuntimeError("local moran boom")

    fake_esda_moran = ModuleType("esda.moran")
    fake_esda_moran.Moran_Local = _BrokenMoranLocal
    fake_libpysal_weights = ModuleType("libpysal.weights")
    fake_libpysal_weights.W = _FakePySALWeights
    _patch_required_modules(
        monkeypatch,
        {
            "esda.moran": fake_esda_moran,
            "libpysal.weights": fake_libpysal_weights,
        },
    )

    with pytest.raises(
        ProcessingError, match="Local Moran's I analysis failed: local moran boom"
    ):
        ss._analyze_local_moran(
            adata,
            SpatialStatisticsParameters(
                analysis_type="local_moran",
                genes=["gene_0"],
                local_moran_fdr_correction=False,
            ),
            DummyCtx(adata),
        )


def test_analyze_network_properties_wraps_internal_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    fake_nx = ModuleType("networkx")
    fake_nx.from_scipy_sparse_array = lambda *_args, **_kwargs: (_ for _ in ()).throw(
        RuntimeError("graph build failed")
    )
    monkeypatch.setattr(ss, "require", lambda *_args, **_kwargs: fake_nx)

    with pytest.raises(
        ProcessingError, match="Network properties analysis failed: graph build failed"
    ):
        ss._analyze_network_properties(
            adata,
            SpatialStatisticsParameters(
                analysis_type="network_properties",
                cluster_key="group",
                n_neighbors=4,
            ),
            DummyCtx(adata),
        )


def test_analyze_spatial_centrality_wraps_internal_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["group"] = pd.Categorical(["a"] * 30 + ["b"] * 30)
    fake_nx = ModuleType("networkx")
    fake_nx.from_scipy_sparse_array = lambda *_args, **_kwargs: object()
    fake_nx.degree_centrality = lambda *_args, **_kwargs: (_ for _ in ()).throw(
        RuntimeError("centrality failed")
    )
    fake_nx.closeness_centrality = lambda *_args, **_kwargs: {}
    fake_nx.betweenness_centrality = lambda *_args, **_kwargs: {}
    monkeypatch.setattr(ss, "require", lambda *_args, **_kwargs: fake_nx)

    with pytest.raises(
        ProcessingError, match="Spatial centrality analysis failed: centrality failed"
    ):
        ss._analyze_spatial_centrality(
            adata,
            SpatialStatisticsParameters(
                analysis_type="spatial_centrality",
                cluster_key="group",
                n_neighbors=4,
            ),
            DummyCtx(adata),
        )


@pytest.mark.asyncio
async def test_analyze_spatial_statistics_metadata_includes_mean_score_when_present(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = pd.Categorical(["0"] * 30 + ["1"] * 30)
    params = SpatialStatisticsParameters(
        analysis_type="neighborhood", cluster_key="leiden"
    )
    captured: dict[str, object] = {}

    monkeypatch.setattr(ss, "ensure_spatial_neighbors", lambda *_args, **_kwargs: None)
    fake_result = {
        "n_clusters": 2,
        "max_enrichment": 1.1,
        "min_enrichment": -0.2,
        "analysis_key": "leiden_nhood_enrichment",
        "mean_score": 0.45,
    }
    monkeypatch.setitem(
        ss._ANALYSIS_REGISTRY,
        "neighborhood",
        ss._AnalysisConfig(
            handler=lambda *_a, **_kw: fake_result,
            needs_cluster=True,
            metadata_keys={"uns": ["neighborhood"]},
        ),
    )
    monkeypatch.setattr(ss, "export_analysis_result", lambda *_args, **_kwargs: [])

    def _store(_adata, *, statistics, **_kwargs):
        captured["statistics"] = statistics

    monkeypatch.setattr(ss, "store_analysis_metadata", _store)

    out = await ss.analyze_spatial_statistics("d1", DummyCtx(adata), params)
    assert out.analysis_type == "neighborhood"
    assert captured["statistics"] == {
        "n_cells": adata.n_obs,
        "mean_score": 0.45,
    }
