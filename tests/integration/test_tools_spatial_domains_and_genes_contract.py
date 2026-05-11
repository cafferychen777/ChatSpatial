"""Integration contracts for spatial_domains and spatial_genes tool entrypoints."""

from __future__ import annotations

import pandas as pd
import pytest

from chatspatial.models.analysis import SpatialDomainResult, SpatialVariableGenesResult
from chatspatial.models.data import SpatialDomainParameters, SpatialVariableGenesParameters
from chatspatial.tools import spatial_domains as domains_module
from chatspatial.tools import spatial_genes as genes_module
from chatspatial.tools.spatial_domains import identify_spatial_domains
from chatspatial.tools.spatial_genes import identify_spatial_genes
from chatspatial.utils.exceptions import DataError, DataNotFoundError


class DummyCtx:
    def __init__(self, adata):
        self._adata = adata
        self.warnings: list[str] = []

    async def get_adata(self, data_id: str):
        return self._adata

    async def warning(self, msg: str):
        self.warnings.append(msg)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_identify_spatial_domains_leiden_contract_with_mocked_backend(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def fake_clustering(adata_subset, params, ctx):
        labels = pd.Series(["0"] * 30 + ["1"] * 30, index=adata_subset.obs.index)
        return labels, "X_pca", {"method": "leiden", "resolution": params.resolution}

    monkeypatch.setattr(
        domains_module, "_identify_domains_clustering", fake_clustering
    )
    monkeypatch.setattr(domains_module, "store_analysis_metadata", lambda *a, **k: None)
    monkeypatch.setattr(domains_module, "export_analysis_result", lambda *a, **k: None)

    params = SpatialDomainParameters(method="leiden", refine_domains=False)
    result = await identify_spatial_domains("d1", ctx, params)

    assert isinstance(result, SpatialDomainResult)
    assert result.data_id == "d1"
    assert result.method == "leiden"
    assert result.domain_key == "spatial_domains_leiden_res0_50"
    assert sum(result.domain_counts.values()) == adata.n_obs
    assert "spatial_domains_leiden_res0_50" in adata.obs


@pytest.mark.integration
@pytest.mark.asyncio
async def test_identify_spatial_domains_requires_spatial_coords(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    ctx = DummyCtx(adata)

    params = SpatialDomainParameters(method="leiden", refine_domains=False)
    with pytest.raises(DataNotFoundError, match="No spatial coordinates found"):
        await identify_spatial_domains("d2", ctx, params)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_identify_spatial_genes_sparkx_dispatch_contract(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def fake_sparkx(data_id, adata_obj, params, ctx_obj):
        return SpatialVariableGenesResult(
            data_id=data_id,
            method="sparkx",
            n_genes_analyzed=24,
            n_significant_genes=2,
            spatial_genes=["gene_1", "gene_2"],
            results_key=f"sparkx_results_{data_id}",
        )

    monkeypatch.setattr(genes_module, "_identify_spatial_genes_sparkx", fake_sparkx)

    params = SpatialVariableGenesParameters(method="sparkx")
    result = await identify_spatial_genes("d3", ctx, params)

    assert isinstance(result, SpatialVariableGenesResult)
    assert result.method == "sparkx"
    assert result.results_key == "sparkx_results_d3"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_identify_spatial_genes_flashs_dispatch_contract(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    async def fake_flashs(data_id, adata_obj, params, ctx_obj):
        return SpatialVariableGenesResult(
            data_id=data_id,
            method="flashs",
            n_genes_analyzed=18,
            n_significant_genes=3,
            spatial_genes=["gene_0", "gene_5", "gene_9"],
            results_key=f"flashs_results_{data_id}",
        )

    monkeypatch.setattr(genes_module, "_identify_spatial_genes_flashs", fake_flashs)

    params = SpatialVariableGenesParameters(method="flashs")
    result = await identify_spatial_genes("d3f", ctx, params)

    assert isinstance(result, SpatialVariableGenesResult)
    assert result.method == "flashs"
    assert result.results_key == "flashs_results_d3f"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_identify_spatial_genes_requires_spatial_coordinates(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    ctx = DummyCtx(adata)

    with pytest.raises(DataError):
        await identify_spatial_genes("d4", ctx, SpatialVariableGenesParameters(method="sparkx"))
