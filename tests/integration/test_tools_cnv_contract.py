"""Integration contracts for tools.cnv_analysis.infer_cnv."""

from __future__ import annotations

import pytest

from chatspatial.models.analysis import CNVResult
from chatspatial.models.data import CNVParameters
from chatspatial.tools import cnv_analysis as cnv_module
from chatspatial.tools.cnv_analysis import infer_cnv
from chatspatial.utils.exceptions import ParameterError


class DummyCtx:
    def __init__(self, adata):
        self._adata = adata

    async def get_adata(self, data_id: str):
        return self._adata

    async def set_adata(self, data_id: str, adata):
        self._adata = adata


@pytest.mark.integration
@pytest.mark.asyncio
async def test_infer_cnv_rejects_missing_reference_category(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    ctx = DummyCtx(adata)

    params = CNVParameters(
        method="infercnvpy",
        reference_key="cell_type",
        reference_categories=["C"],  # absent
    )

    with pytest.raises(ParameterError, match="not found"):
        await infer_cnv("d1", ctx, params)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_infer_cnv_dispatches_to_infercnvpy_handler(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    ctx = DummyCtx(adata)

    async def fake_infercnvpy(data_id, adata, params, ctx):
        return CNVResult(
            data_id=data_id,
            method="infercnvpy",
            reference_key=params.reference_key,
            reference_categories=list(params.reference_categories),
            n_chromosomes=1,
            n_genes_analyzed=adata.n_vars,
            cnv_score_key="X_cnv",
            visualization_available=True,
        )

    monkeypatch.setattr(cnv_module, "_infer_cnv_infercnvpy", fake_infercnvpy)

    params = CNVParameters(
        method="infercnvpy",
        reference_key="cell_type",
        reference_categories=["A"],
    )
    result = await infer_cnv("d2", ctx, params)

    assert isinstance(result, CNVResult)
    assert result.method == "infercnvpy"
    assert result.data_id == "d2"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_infer_cnv_dispatches_to_numbat_handler(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    ctx = DummyCtx(adata)

    def fake_numbat(data_id, adata, params, ctx):
        return CNVResult(
            data_id=data_id,
            method="numbat",
            reference_key=params.reference_key,
            reference_categories=list(params.reference_categories),
            n_chromosomes=0,
            n_genes_analyzed=adata.n_vars,
            cnv_score_key="X_cnv_numbat",
            visualization_available=True,
        )

    monkeypatch.setattr(cnv_module, "_infer_cnv_numbat", fake_numbat)

    params = CNVParameters(
        method="numbat",
        reference_key="cell_type",
        reference_categories=["A"],
    )
    result = await infer_cnv("d3", ctx, params)

    assert isinstance(result, CNVResult)
    assert result.method == "numbat"
    assert result.data_id == "d3"
