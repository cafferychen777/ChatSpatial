"""Integration contracts for deconvolution module entrypoints."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import chatspatial.tools.deconvolution as deconv_module
from chatspatial.models.analysis import DeconvolutionResult
from chatspatial.models.data import DeconvolutionParameters
from chatspatial.tools.deconvolution import (
    _store_results,
    _validate_and_align_proportions,
    deconvolve_spatial_data,
)
from chatspatial.tools.deconvolution.base import PreparedDeconvolutionData
from chatspatial.utils.exceptions import DataError, ParameterError, ProcessingError


class DummyCtx:
    def __init__(self, datasets: dict[str, object]):
        self.datasets = datasets
        self.updated: dict[str, object] = {}
        self.warnings: list[str] = []

    async def get_adata(self, data_id: str):
        return self.datasets[data_id]

    async def set_adata(self, data_id: str, adata):
        self.updated[data_id] = adata

    async def warning(self, message: str):
        self.warnings.append(message)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_deconvolve_spatial_data_requires_nonempty_data_id():
    params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id="ref",
        cell_type_key="cell_type",
    )
    with pytest.raises(ParameterError, match="Dataset ID cannot be empty"):
        await deconvolve_spatial_data("", DummyCtx({}), params)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_deconvolve_spatial_data_requires_reference_data_id(
    minimal_spatial_adata,
):
    spatial = minimal_spatial_adata.copy()
    spatial.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    ctx = DummyCtx({"d1": spatial})

    params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id=None,
        cell_type_key="cell_type",
    )
    with pytest.raises(ParameterError, match="requires reference_data_id"):
        await deconvolve_spatial_data("d1", ctx, params)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_deconvolve_spatial_data_rejects_empty_spatial_dataset(
    minimal_spatial_adata,
):
    spatial = minimal_spatial_adata[:0, :].copy()
    reference = minimal_spatial_adata.copy()
    reference.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    ctx = DummyCtx({"d1": spatial, "ref": reference})

    params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id="ref",
        cell_type_key="cell_type",
    )
    with pytest.raises(DataError, match="contains no observations"):
        await deconvolve_spatial_data("d1", ctx, params)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_deconvolve_spatial_data_dispatch_contract_with_mocks(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()
    reference.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    ctx = DummyCtx({"d1": spatial, "ref": reference})

    monkeypatch.setattr(
        deconv_module, "_check_method_availability", lambda *a, **k: None
    )

    prepared = PreparedDeconvolutionData(
        spatial=spatial,
        reference=reference,
        cell_type_key="cell_type",
        cell_types=["A", "B"],
        common_genes=list(spatial.var_names),
        spatial_coords=spatial.obsm["spatial"],
        ctx=ctx,
    )
    calls: dict[str, object] = {}

    async def fake_prepare_deconvolution(**kwargs):
        calls["prepared"] = True
        return prepared

    def fake_dispatch(data, params, config):
        calls["method"] = params.method
        props = pd.DataFrame(
            {"A": [0.8] * data.spatial.n_obs, "B": [0.2] * data.spatial.n_obs},
            index=data.spatial.obs_names,
        )
        return props, {"n_spots": data.spatial.n_obs, "genes_used": data.n_genes}

    async def fake_store(
        spatial_adata, proportions, stats, method, data_id, ctx, **_kw
    ):
        calls["stored"] = (method, data_id, len(proportions))
        calls["result_adata"] = spatial_adata
        return DeconvolutionResult(
            data_id=data_id,
            method=method,
            dominant_type_key=f"dominant_celltype_{method}",
            n_cell_types=2,
            cell_types=["A", "B"],
            proportions_key=f"deconvolution_{method}",
            n_spots=stats["n_spots"],
            genes_used=stats["genes_used"],
        )

    monkeypatch.setattr(
        deconv_module, "prepare_deconvolution", fake_prepare_deconvolution
    )
    monkeypatch.setattr(deconv_module, "_dispatch_method", fake_dispatch)
    monkeypatch.setattr(deconv_module, "_store_results", fake_store)

    params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id="ref",
        cell_type_key="cell_type",
    )
    result = await deconvolve_spatial_data("d1", ctx, params)

    assert isinstance(result, DeconvolutionResult)
    assert result.method == "flashdeconv"
    assert calls["prepared"] is True
    assert calls["method"] == "flashdeconv"
    assert calls["stored"] == ("flashdeconv", "d1", spatial.n_obs)
    assert calls["result_adata"] is not spatial


@pytest.mark.integration
@pytest.mark.asyncio
async def test_deconvolution_late_failure_does_not_publish_partial_results(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()
    reference.obs["cell_type"] = ["A"] * 30 + ["B"] * 30
    ctx = DummyCtx({"d1": spatial, "ref": reference})

    monkeypatch.setattr(
        deconv_module, "_check_method_availability", lambda *_args, **_kwargs: None
    )

    async def _prepare(**kwargs):
        return PreparedDeconvolutionData(
            spatial=kwargs["spatial_adata"].copy(),
            reference=kwargs["reference_adata"].copy(),
            cell_type_key="cell_type",
            cell_types=["A", "B"],
            common_genes=list(spatial.var_names),
            spatial_coords=spatial.obsm["spatial"],
            ctx=ctx,
        )

    def _dispatch(data, _params, _config):
        proportions = pd.DataFrame(
            {"A": np.full(data.n_spots, 0.7), "B": np.full(data.n_spots, 0.3)},
            index=data.spatial.obs_names,
        )
        return proportions, {"genes_used": data.n_genes}

    monkeypatch.setattr(deconv_module, "prepare_deconvolution", _prepare)
    monkeypatch.setattr(deconv_module, "_dispatch_method", _dispatch)
    monkeypatch.setattr(
        deconv_module, "store_analysis_metadata", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        deconv_module,
        "export_analysis_result",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("export failed")),
    )

    params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id="ref",
        cell_type_key="cell_type",
    )
    with pytest.raises(RuntimeError, match="export failed"):
        await deconvolve_spatial_data("d1", ctx, params)

    assert "deconvolution_flashdeconv" not in spatial.obsm
    assert "dominant_celltype_flashdeconv" not in spatial.obs
    assert ctx.updated == {}


@pytest.mark.integration
@pytest.mark.asyncio
async def test_store_results_persists_expected_keys_and_calls_set_adata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})

    monkeypatch.setattr(deconv_module, "store_analysis_metadata", lambda *a, **k: None)
    monkeypatch.setattr(deconv_module, "export_analysis_result", lambda *a, **k: None)

    proportions = pd.DataFrame(
        {"T": np.full(adata.n_obs, 0.6), "B": np.full(adata.n_obs, 0.4)},
        index=adata.obs_names,
    )
    stats = {"n_spots": adata.n_obs, "genes_used": adata.n_vars}

    result = await _store_results(
        spatial_adata=adata,
        proportions=proportions,
        stats=stats,
        method="flashdeconv",
        data_id="d1",
        ctx=ctx,
    )

    assert isinstance(result, DeconvolutionResult)
    assert result.proportions_key == "deconvolution_flashdeconv"
    assert "deconvolution_flashdeconv" in adata.obsm
    assert "dominant_celltype_flashdeconv" in adata.obs
    assert "d1" in ctx.updated


def test_validate_and_align_proportions_reorders_without_inventing_rows(
    minimal_spatial_adata,
):
    obs_names = minimal_spatial_adata.obs_names
    proportions = pd.DataFrame(
        {"T": [0.25, 0.75], "B": [0.75, 0.25]},
        index=[obs_names[1], obs_names[0]],
    )

    aligned = _validate_and_align_proportions(
        proportions, pd.Index(obs_names[:2]), "test_method"
    )

    assert aligned.index.tolist() == obs_names[:2].tolist()
    np.testing.assert_array_equal(aligned["T"], [0.75, 0.25])


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (lambda df: df.iloc[:-1], "1 missing"),
        (lambda df: df.assign(T=np.nan), "non-finite"),
        (lambda df: df.assign(T=-0.1), "negative proportions"),
        (lambda df: df.rename(columns={"T": "B"}), "duplicate cell-type"),
    ],
)
def test_validate_and_align_proportions_rejects_invalid_backend_output(
    minimal_spatial_adata, mutator, message
):
    obs_names = minimal_spatial_adata.obs_names[:3]
    valid = pd.DataFrame(
        {"T": [0.6, 0.6, 0.6], "B": [0.4, 0.4, 0.4]},
        index=obs_names,
    )

    with pytest.raises(ProcessingError, match=message):
        _validate_and_align_proportions(mutator(valid), pd.Index(obs_names), "backend")


@pytest.mark.integration
@pytest.mark.asyncio
async def test_store_results_handles_zero_sum_rows(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """When some spots have all-zero proportions, 'unassigned' appears in
    obs[dominant_key] but NOT in cell_types/n_cell_types (which must match
    the proportions matrix columns exactly)."""
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})

    captured: dict[str, object] = {}
    monkeypatch.setattr(
        deconv_module,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(deconv_module, "export_analysis_result", lambda *a, **k: None)

    # First row all-zero → "unassigned" in dominant type annotation
    proportions = pd.DataFrame(
        {
            "T": np.concatenate([[0.0], np.full(adata.n_obs - 1, 0.6)]),
            "B": np.concatenate([[0.0], np.full(adata.n_obs - 1, 0.4)]),
        },
        index=adata.obs_names,
    )
    stats = {"n_spots": adata.n_obs, "genes_used": adata.n_vars}

    result = await _store_results(
        spatial_adata=adata,
        proportions=proportions,
        stats=stats,
        method="flashdeconv",
        data_id="d1",
        ctx=ctx,
    )

    # cell_types/n_cell_types must match proportions matrix columns (no unassigned)
    assert result.cell_types == ["T", "B"]
    assert result.n_cell_types == 2

    # Metadata cell_types must also match matrix columns
    assert captured["statistics"]["cell_types"] == ["T", "B"]
    assert captured["statistics"]["n_cell_types"] == 2
    assert captured["statistics"]["has_unassigned_spots"] is True

    # But obs[dominant_key] labels zero-sum rows as "unassigned"
    dominant = adata.obs["dominant_celltype_flashdeconv"]
    assert "unassigned" in dominant.values

    # A spot without proportions carries no result, so the caller has to hear
    # about it instead of reading an empty map as a finding.
    assert any("unassigned" in warning for warning in ctx.warnings)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_store_results_warns_loudly_when_every_spot_is_unassigned(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """An all-zero proportion matrix means the backend produced nothing."""
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})

    monkeypatch.setattr(
        deconv_module, "store_analysis_metadata", lambda _adata, **kwargs: None
    )
    monkeypatch.setattr(deconv_module, "export_analysis_result", lambda *a, **k: None)

    proportions = pd.DataFrame(
        {
            "T": np.zeros(adata.n_obs),
            "B": np.zeros(adata.n_obs),
        },
        index=adata.obs_names,
    )

    await _store_results(
        spatial_adata=adata,
        proportions=proportions,
        stats={"n_spots": adata.n_obs, "genes_used": adata.n_vars},
        method="flashdeconv",
        data_id="d1",
        ctx=ctx,
    )

    assert len(ctx.warnings) == 1
    assert "100.0%" in ctx.warnings[0]
    assert "No spot received any signal" in ctx.warnings[0]


@pytest.mark.integration
@pytest.mark.asyncio
async def test_store_results_no_unassigned_when_all_rows_have_signal(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """When no zero-sum rows exist, 'unassigned' should NOT be in cell_types."""
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx({"d1": adata})

    monkeypatch.setattr(deconv_module, "store_analysis_metadata", lambda *a, **k: None)
    monkeypatch.setattr(deconv_module, "export_analysis_result", lambda *a, **k: None)

    proportions = pd.DataFrame(
        {"T": np.full(adata.n_obs, 0.6), "B": np.full(adata.n_obs, 0.4)},
        index=adata.obs_names,
    )

    result = await _store_results(
        spatial_adata=adata,
        proportions=proportions,
        stats={"n_spots": adata.n_obs, "genes_used": adata.n_vars},
        method="flashdeconv",
        data_id="d1",
        ctx=ctx,
    )

    assert "unassigned" not in result.cell_types
    assert result.n_cell_types == 2
