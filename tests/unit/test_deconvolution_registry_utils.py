"""Unit tests for deconvolution registry/dispatch helper logic."""

from __future__ import annotations

import importlib.util
from types import SimpleNamespace

import pandas as pd
import pytest

import chatspatial.tools.deconvolution as deconv_module
from chatspatial.models.data import DeconvolutionParameters
from chatspatial.tools.deconvolution.base import MethodConfig
from chatspatial.utils.exceptions import DataError, DependencyError, ParameterError


class _Ctx:
    def __init__(self, datasets: dict[str, object]) -> None:
        self.datasets = datasets

    async def get_adata(self, data_id: str):
        return self.datasets[data_id]

    async def set_adata(self, data_id: str, adata):
        self.datasets[data_id] = adata


@pytest.mark.asyncio
async def test_deconvolve_spatial_data_rejects_unsupported_method() -> None:
    params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id="ref",
        cell_type_key="cell_type",
    )
    params.method = "not_real"

    with pytest.raises(ParameterError, match="Unsupported method: not_real"):
        await deconv_module.deconvolve_spatial_data("d1", _Ctx({}), params)


@pytest.mark.asyncio
async def test_deconvolve_spatial_data_rejects_empty_reference_dataset(
    minimal_spatial_adata,
) -> None:
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata[:0, :].copy()
    ctx = _Ctx({"d1": spatial, "ref": reference})

    params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id="ref",
        cell_type_key="cell_type",
    )

    with pytest.raises(
        DataError, match="Reference dataset ref contains no observations"
    ):
        await deconv_module.deconvolve_spatial_data("d1", ctx, params)


def test_check_method_availability_maps_scvi_tools_to_scvi(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[str] = []

    def _find_spec(name: str):
        calls.append(name)
        return object() if name == "scvi" else None

    monkeypatch.setattr(importlib.util, "find_spec", _find_spec)

    deconv_module._check_method_availability(
        "dummy",
        MethodConfig(module_name="dummy", dependencies=("scvi-tools",)),
    )

    assert calls == ["scvi"]


def test_tangram_availability_depends_only_on_its_runtime_package(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[str] = []

    def _find_spec(name: str):
        calls.append(name)
        return object() if name == "tangram" else None

    monkeypatch.setattr(importlib.util, "find_spec", _find_spec)

    deconv_module._check_method_availability(
        "tangram", deconv_module.METHOD_REGISTRY["tangram"]
    )

    assert calls == ["tangram"]


def test_check_method_availability_reports_available_alternatives(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _find_spec(name: str):
        # Keep only flashdeconv available to validate fallback recommendation text.
        return object() if name == "flashdeconv" else None

    monkeypatch.setattr(importlib.util, "find_spec", _find_spec)

    with pytest.raises(DependencyError) as exc_info:
        deconv_module._check_method_availability(
            "rctd",
            deconv_module.METHOD_REGISTRY["rctd"],
        )

    message = str(exc_info.value)
    assert "requires: rpy2" in message
    assert "Available: flashdeconv" in message
    assert "flashdeconv recommended - fastest" in message


def test_check_method_availability_handles_no_available_methods(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: None)

    with pytest.raises(DependencyError) as exc_info:
        deconv_module._check_method_availability(
            "rctd",
            deconv_module.METHOD_REGISTRY["rctd"],
        )

    message = str(exc_info.value)
    assert "requires: rpy2" in message
    assert "Available:" not in message


def test_dispatch_method_dynamic_import_and_kwarg_extraction(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    params = DeconvolutionParameters(
        method="flashdeconv",
        reference_data_id="ref",
        cell_type_key="cell_type",
        use_gpu=True,
        cell2location_n_epochs=123,
    )
    config = MethodConfig(
        module_name="fake_backend",
        dependencies=(),
        supports_gpu=True,
        param_mapping=(("cell2location_n_epochs", "n_epochs"),),
    )

    captured: dict[str, object] = {}

    def _fake_deconvolve(data, **kwargs):
        captured["data"] = data
        captured["kwargs"] = kwargs
        return pd.DataFrame({"A": [1.0]}, index=["spot_1"]), {"ok": True}

    def _fake_import(name: str, package: str | None = None):
        captured["import"] = (name, package)
        return SimpleNamespace(deconvolve=_fake_deconvolve)

    monkeypatch.setattr(deconv_module.importlib, "import_module", _fake_import)

    dummy_data = object()
    proportions, stats = deconv_module._dispatch_method(dummy_data, params, config)

    assert proportions.index.tolist() == ["spot_1"]
    assert stats == {"ok": True}
    assert captured["import"] == (".fake_backend", deconv_module.__package__)
    assert captured["data"] is dummy_data
    assert captured["kwargs"] == {"n_epochs": 123, "use_gpu": True}


def test_get_preprocess_hook_returns_none_when_filtering_disabled() -> None:
    params = DeconvolutionParameters(
        method="cell2location",
        reference_data_id="ref",
        cell_type_key="cell_type",
        cell2location_apply_gene_filtering=False,
    )

    assert deconv_module._get_preprocess_hook(params) is None


@pytest.mark.asyncio
async def test_get_preprocess_hook_applies_filtering_to_spatial_and_reference(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from chatspatial.tools.deconvolution import cell2location as c2l_module

    params = DeconvolutionParameters(
        method="cell2location",
        reference_data_id="ref",
        cell_type_key="cell_type",
        cell2location_apply_gene_filtering=True,
        cell2location_gene_filter_cell_count_cutoff=9,
        cell2location_gene_filter_cell_percentage_cutoff2=0.25,
        cell2location_gene_filter_nonz_mean_cutoff=1.75,
    )

    calls: list[dict[str, float | int]] = []

    async def _fake_apply_gene_filtering(adata, _ctx, **kwargs):
        calls.append(kwargs)
        out = adata.copy()
        out.uns["filtered"] = True
        return out

    monkeypatch.setattr(c2l_module, "apply_gene_filtering", _fake_apply_gene_filtering)

    hook = deconv_module._get_preprocess_hook(params)
    assert hook is not None

    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()

    spatial_out, reference_out = await hook(spatial, reference, object())

    assert spatial_out is not spatial
    assert reference_out is not reference
    assert spatial_out.uns["filtered"] is True
    assert reference_out.uns["filtered"] is True
    assert len(calls) == 2
    assert calls[0] == {
        "cell_count_cutoff": 9,
        "cell_percentage_cutoff2": 0.25,
        "nonz_mean_cutoff": 1.75,
    }
    assert calls[1] == calls[0]
