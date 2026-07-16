"""Integration contracts for tools.cell_communication entrypoints."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import chatspatial.tools.cell_communication as ccc_module
from chatspatial.models.analysis import CellCommunicationResult
from chatspatial.models.data import CellCommunicationParameters
from chatspatial.tools.cell_communication import (
    CCCStorage,
    _run_ccc_analysis,
    _validate_ccc_params,
    analyze_cell_communication,
)
from chatspatial.utils.exceptions import DataNotFoundError, ParameterError, ProcessingError


class DummyCtx:
    def __init__(self, adata):
        self._adata = adata
        self.warnings: list[str] = []

    async def get_adata(self, data_id: str):
        return self._adata

    async def set_adata(self, data_id: str, adata):
        self._adata = adata

    async def warning(self, msg: str):
        self.warnings.append(msg)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_validate_ccc_params_requires_spatial_connectivity_for_liana_spatial(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["T"] * 30 + ["B"] * 30
    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
        perform_spatial_analysis=True,
    )

    with pytest.raises(DataNotFoundError, match="Spatial connectivity required"):
        await _validate_ccc_params(adata, params, DummyCtx(adata))


@pytest.mark.integration
@pytest.mark.asyncio
async def test_validate_ccc_params_warns_for_mouse_consensus_choice(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["T"] * 30 + ["B"] * 30
    adata.obsp["spatial_connectivities"] = np.eye(adata.n_obs)
    ctx = DummyCtx(adata)

    params = CellCommunicationParameters(
        method="liana",
        species="mouse",
        cell_type_key="cell_type",
        liana_resource="consensus",
        perform_spatial_analysis=True,
    )
    await _validate_ccc_params(adata, params, ctx)

    assert any("mouseconsensus" in msg for msg in ctx.warnings)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_run_ccc_analysis_rejects_unsupported_method(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx(adata)

    class BadParams:
        method = "not_supported"

    with pytest.raises(ParameterError, match="Unsupported method"):
        await _run_ccc_analysis(adata, BadParams(), ctx)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_cell_communication_stores_results_and_returns_contract(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["T"] * 30 + ["B"] * 30
    ctx = DummyCtx(adata)

    async def fake_validate(_adata, _params, _ctx):
        return None

    async def fake_run(_adata, _params, _ctx):
        return CCCStorage(
            method="liana",
            analysis_type="spatial",
            species="human",
            database="consensus",
            lr_pairs=["L1^R1", "L2^R2"],
            top_lr_pairs=["L1^R1"],
            n_pairs=2,
            n_significant=1,
            results=pd.DataFrame({"score": [0.5]}, index=["L1_R1"]),
            method_data={
                "spatial_scores": np.ones((adata.n_obs, 2), dtype=float),
                "spatial_pvals": np.zeros((adata.n_obs, 2), dtype=float),
            },
            statistics={"demo": True},
        )

    monkeypatch.setattr(ccc_module, "_validate_ccc_params", fake_validate)
    monkeypatch.setattr(ccc_module, "_run_ccc_analysis", fake_run)
    monkeypatch.setattr("chatspatial.utils.adata_utils.store_analysis_metadata", lambda *a, **k: None)
    monkeypatch.setattr("chatspatial.utils.results_export.export_analysis_result", lambda *a, **k: None)

    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
        perform_spatial_analysis=True,
    )
    result = await analyze_cell_communication("d1", ctx, params)

    assert isinstance(result, CellCommunicationResult)
    assert result.method == "liana"
    assert result.analysis_type == "spatial"
    assert result.n_lr_pairs == 2
    assert ctx._adata is not adata
    assert "ccc" not in adata.uns
    assert "ccc" in ctx._adata.uns
    assert "ccc_spatial_scores" in ctx._adata.obsm
    assert "ccc_spatial_pvals" in ctx._adata.obsm


@pytest.mark.integration
@pytest.mark.asyncio
async def test_analyze_cell_communication_wraps_errors_as_processing_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["cell_type"] = ["T"] * 30 + ["B"] * 30
    ctx = DummyCtx(adata)

    async def fake_validate(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(ccc_module, "_validate_ccc_params", fake_validate)

    params = CellCommunicationParameters(
        method="liana",
        species="human",
        cell_type_key="cell_type",
    )
    with pytest.raises(ProcessingError, match="Error in cell communication analysis"):
        await analyze_cell_communication("d2", ctx, params)
