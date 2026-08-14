"""Integration contract tests for core analysis endpoints in server layer."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from chatspatial.models.analysis import (
    EnrichmentResult,
    IntegrationResult,
    SpatialStatisticsResult,
)
from chatspatial.models.data import (
    EnrichmentParameters,
    IntegrationParameters,
    SpatialStatisticsParameters,
)
from chatspatial.server import (
    analyze_enrichment,
    analyze_spatial_statistics,
    integrate_samples,
)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_spatial_statistics_passes_inputs_and_returns_public_result(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    calls: dict[str, object] = {}

    async def fake_analyze(data_id, ctx, params):
        calls["data_id"] = data_id
        calls["ctx"] = ctx
        calls["params"] = params
        return SpatialStatisticsResult(
            data_id=data_id,
            analysis_type="moran",
            n_features_analyzed=3,
            n_significant=1,
            top_features=["gene_0"],
            summary_metrics={"moran_i_mean": 0.42},
            results_key="spatial_stats_moran",
        )

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.spatial_statistics",
        SimpleNamespace(analyze_spatial_statistics=fake_analyze),
    )

    params = SpatialStatisticsParameters(analysis_type="moran", genes=["gene_0"])
    result = await analyze_spatial_statistics("d1", params=params)

    assert isinstance(result, SpatialStatisticsResult)
    assert calls["data_id"] == "d1"
    assert calls["params"] is params
    assert result.analysis_type == "moran"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_integrate_samples_passes_inputs_and_returns_integrated_id(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    calls: dict[str, object] = {}

    async def fake_integrate(data_ids, ctx, params):
        calls["data_ids"] = data_ids
        calls["ctx"] = ctx
        calls["params"] = params
        return IntegrationResult(
            data_id="integrated_7",
            n_samples=len(data_ids),
            integration_method=params.method,
            n_cells=120,
            batch_key=params.batch_key,
            embedding_key="X_pca_harmony",
        )

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.integration",
        SimpleNamespace(integrate_samples=fake_integrate),
    )

    params = IntegrationParameters(method="harmony")
    result = await integrate_samples(["d1", "d2"], params=params)

    assert isinstance(result, IntegrationResult)
    assert calls["data_ids"] == ["d1", "d2"]
    assert calls["params"] is params
    assert result.data_id == "integrated_7"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_enrichment_passes_inputs_and_returns_public_result(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    calls: dict[str, object] = {}

    async def fake_enrichment(data_id, ctx, params):
        calls["data_id"] = data_id
        calls["ctx"] = ctx
        calls["params"] = params
        return EnrichmentResult(
            method="pathway_ora",
            n_gene_sets=4,
            n_significant=2,
            top_gene_sets=["GO_A", "GO_B"],
            top_depleted_sets=["GO_C"],
        )

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.enrichment",
        SimpleNamespace(analyze_enrichment=fake_enrichment),
    )

    params = EnrichmentParameters(species="human", method="pathway_ora")
    result = await analyze_enrichment("d5", params=params)

    assert isinstance(result, EnrichmentResult)
    assert calls["data_id"] == "d5"
    assert calls["params"] is params


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_spatial_statistics_materializes_default_params(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_analyze(data_id, ctx, params):
        assert isinstance(params, SpatialStatisticsParameters)
        return SpatialStatisticsResult(
            data_id=data_id,
            analysis_type=params.analysis_type,
            n_features_analyzed=0,
            n_significant=0,
            top_features=[],
            summary_metrics={},
            results_key="spatial_stats_default",
        )

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.spatial_statistics",
        SimpleNamespace(analyze_spatial_statistics=fake_analyze),
    )

    out = await analyze_spatial_statistics("d_defaults")
    assert isinstance(out, SpatialStatisticsResult)
    assert out.analysis_type == "neighborhood"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_integrate_samples_materializes_default_params(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_integrate(data_ids, ctx, params):
        assert isinstance(params, IntegrationParameters)
        return IntegrationResult(
            data_id="integrated_defaults",
            n_samples=len(data_ids),
            integration_method=params.method,
            n_cells=64,
            batch_key=params.batch_key,
            embedding_key="X_pca_harmony",
        )

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.integration",
        SimpleNamespace(integrate_samples=fake_integrate),
    )

    out = await integrate_samples(["d1", "d2"])
    assert isinstance(out, IntegrationResult)
    assert out.integration_method == "harmony"
