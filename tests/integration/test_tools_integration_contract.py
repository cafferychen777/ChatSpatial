"""Integration tests for chatspatial.tools.integration contracts."""

from __future__ import annotations

import numpy as np
import pytest

from chatspatial.models.analysis import IntegrationResult
from chatspatial.models.data import IntegrationParameters
from chatspatial.tools import integration as integration_module
from chatspatial.tools.integration import integrate_multiple_samples, integrate_samples
from chatspatial.utils.exceptions import DataError, ParameterError


class DummyCtx:
    def __init__(self, datasets: dict[str, object]):
        self.datasets = datasets
        self.added: dict[str, object] = {}

    async def get_adata(self, data_id: str):
        return self.datasets[data_id]

    async def add_dataset(self, adata, prefix: str = "data"):
        data_id = f"{prefix}_1"
        self.added[data_id] = adata
        return data_id


@pytest.mark.integration
def test_integrate_multiple_samples_requires_at_least_two_datasets(minimal_spatial_adata):
    with pytest.raises(ParameterError, match="at least 2 datasets"):
        integrate_multiple_samples([minimal_spatial_adata], method="harmony")


@pytest.mark.integration
def test_integrate_multiple_samples_rejects_raw_count_like_input(minimal_spatial_adata):
    adata1 = minimal_spatial_adata.copy()
    adata2 = minimal_spatial_adata.copy()
    adata2.obs_names = [f"raw2_{i}" for i in range(adata2.n_obs)]

    # Construct obvious raw-count-like matrices (non-negative and very large max values).
    adata1.X = np.full(adata1.X.shape, 200.0, dtype=np.float32)
    adata2.X = np.full(adata2.X.shape, 220.0, dtype=np.float32)

    with pytest.raises(DataError, match="raw counts"):
        integrate_multiple_samples([adata1, adata2], method="harmony")


@pytest.mark.integration
@pytest.mark.asyncio
async def test_integrate_samples_preserves_scvi_data_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata1 = minimal_spatial_adata.copy()
    adata2 = minimal_spatial_adata.copy()
    adata2.obs_names = [f"batch2_{i}" for i in range(adata2.n_obs)]
    normalized = np.log1p(np.random.default_rng(0).poisson(4, adata1.X.shape)).astype(
        np.float32
    )
    for adata in (adata1, adata2):
        adata.X = normalized.copy()
        adata.layers["counts"] = np.full(adata.X.shape, 0.25, dtype=np.float32)
        adata.var["highly_variable"] = True

    ctx = DummyCtx({"d1": adata1, "d2": adata2})

    def fake_integrate_with_scvi(*_args, **_kwargs):
        raise DataError("scVI requires raw integer counts")

    monkeypatch.setattr(
        integration_module, "integrate_with_scvi", fake_integrate_with_scvi
    )

    params = IntegrationParameters(method="scvi", batch_key="batch", align_spatial=False)
    with pytest.raises(DataError, match="raw integer counts"):
        await integrate_samples(["d1", "d2"], ctx, params)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_integrate_samples_adds_integrated_dataset_and_exports(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata1 = minimal_spatial_adata.copy()
    adata2 = minimal_spatial_adata.copy()

    ctx = DummyCtx({"d1": adata1, "d2": adata2})

    integrated = minimal_spatial_adata.copy()
    integrated.obsm["spatial"] = integrated.obsm["spatial"].copy()
    integrated.obsm["spatial_aligned"] = integrated.obsm["spatial"].copy()

    called: dict[str, object] = {}

    def fake_integrate_multiple_samples(adatas, batch_key, method, n_pcs, params):
        called["integrate"] = {
            "n_adatas": len(adatas),
            "batch_key": batch_key,
            "method": method,
            "n_pcs": n_pcs,
        }
        return integrated

    def fake_rescale_spatial_coordinates(combined_adata, batch_key, reference_batch):
        called["align"] = {
            "batch_key": batch_key,
            "reference_batch": reference_batch,
            "has_spatial": "spatial" in combined_adata.obsm,
        }
        return combined_adata

    def fake_export_analysis_result(adata, data_id, analysis_name):
        called.setdefault("exports", []).append((data_id, analysis_name))

    monkeypatch.setattr(
        integration_module, "integrate_multiple_samples", fake_integrate_multiple_samples
    )
    monkeypatch.setattr(
        integration_module, "rescale_spatial_coordinates", fake_rescale_spatial_coordinates
    )
    monkeypatch.setattr(
        integration_module, "export_analysis_result", fake_export_analysis_result
    )

    params = IntegrationParameters(
        method="harmony",
        batch_key="batch",
        n_pcs=20,
        align_spatial=True,
    )
    result = await integrate_samples(["d1", "d2"], ctx, params)

    assert isinstance(result, IntegrationResult)
    assert result.data_id == "integrated_1"
    assert result.n_samples == 2
    assert result.integration_method == "harmony"
    assert called["integrate"]["method"] == "harmony"
    assert called["align"]["has_spatial"] is True
    assert ("integrated_1", "integration_harmony") in called["exports"]
    assert ("integrated_1", "spatial_alignment") in called["exports"]
    assert "integrated_1" in ctx.added

