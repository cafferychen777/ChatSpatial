"""Lightweight contracts for R-based deconvolution modules via mocked rpy2."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import replace
from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from chatspatial.tools.deconvolution import card as card_module
from chatspatial.tools.deconvolution import rctd as rctd_module
from chatspatial.tools.deconvolution import spotlight as spotlight_module
from chatspatial.tools.deconvolution.base import PreparedDeconvolutionData
from chatspatial.utils.dependency_manager import REnvironment
from chatspatial.utils.exceptions import (
    DataError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


class DummyCtx:
    async def warning(self, _msg: str):
        return None


def _prepared_data(minimal_spatial_adata) -> PreparedDeconvolutionData:
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()
    reference.obs["cell_type"] = ["A"] * (reference.n_obs // 2) + ["B"] * (
        reference.n_obs - reference.n_obs // 2
    )
    return PreparedDeconvolutionData(
        spatial=spatial,
        reference=reference,
        cell_type_key="cell_type",
        cell_types=["A", "B"],
        common_genes=list(spatial.var_names),
        spatial_coords=spatial.obsm["spatial"],
        ctx=DummyCtx(),
    )


def _install_fake_r_modules(monkeypatch: pytest.MonkeyPatch, ro_r):
    modules = __import__("sys").modules

    class _Converter:
        def __add__(self, _other):
            return self

    @contextmanager
    def _localconverter(_conv):
        yield

    ro_mod = ModuleType("rpy2.robjects")
    ro_mod.globalenv = {}
    ro_mod.default_converter = _Converter()
    ro_mod.StrVector = lambda x: list(x)
    ro_mod.r = ro_r
    ro_mod.conversion = SimpleNamespace(
        py2rpy=lambda x: x,
        rpy2py=lambda x: x,
        localconverter=_localconverter,
    )

    pandas2ri_mod = ModuleType("rpy2.robjects.pandas2ri")
    pandas2ri_mod.converter = _Converter()
    numpy2ri_mod = ModuleType("rpy2.robjects.numpy2ri")
    numpy2ri_mod.converter = _Converter()

    conversion_mod = ModuleType("rpy2.robjects.conversion")
    conversion_mod.localconverter = _localconverter

    class _RLock:
        def __enter__(self):
            return None

        def __exit__(self, *_args):
            return False

    rinterface_lib_mod = ModuleType("rpy2.rinterface_lib")
    rinterface_lib_mod.openrlib = SimpleNamespace(rlock=_RLock())

    rpy2_mod = ModuleType("rpy2")
    anndata2ri_mod = ModuleType("anndata2ri")
    anndata2ri_mod.converter = _Converter()

    monkeypatch.setitem(modules, "rpy2", rpy2_mod)
    monkeypatch.setitem(modules, "rpy2.robjects", ro_mod)
    monkeypatch.setitem(modules, "rpy2.robjects.pandas2ri", pandas2ri_mod)
    monkeypatch.setitem(modules, "rpy2.robjects.numpy2ri", numpy2ri_mod)
    monkeypatch.setitem(modules, "rpy2.robjects.conversion", conversion_mod)
    monkeypatch.setitem(modules, "rpy2.rinterface_lib", rinterface_lib_mod)
    monkeypatch.setitem(modules, "anndata2ri", anndata2ri_mod)

    environment = REnvironment(
        robjects=ro_mod,
        pandas2ri=pandas2ri_mod,
        numpy2ri=numpy2ri_mod,
        packages=SimpleNamespace(importr=lambda _package: object()),
        conversion=conversion_mod,
        openrlib=rinterface_lib_mod.openrlib,
        anndata2ri=anndata2ri_mod,
    )
    for target_module in (rctd_module, card_module, spotlight_module):
        monkeypatch.setattr(
            target_module,
            "validate_r_environment",
            lambda *_args, **_kwargs: environment,
        )

    class _Logger:
        level = 0

        def setLevel(self, level):
            self.level = level

    callbacks = SimpleNamespace(
        consolewrite_print=lambda _text: None,
        consolewrite_warnerror=lambda _text: None,
        logger=_Logger(),
    )
    monkeypatch.setattr(
        rctd_module,
        "require_module",
        lambda *_args, **_kwargs: callbacks,
    )
    return environment


@pytest.mark.parametrize(
    "target_module",
    [rctd_module, card_module, spotlight_module],
)
def test_r_deconvolution_methods_preserve_dependency_errors(
    target_module,
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    data = _prepared_data(minimal_spatial_adata)
    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)
    monkeypatch.setattr(
        target_module,
        "validate_r_environment",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("R runtime unavailable")
        ),
    )

    with pytest.raises(DependencyError, match="R runtime unavailable"):
        target_module.deconvolve(data)


def test_rctd_mode_multi_parameter_guard_before_heavy_execution(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)

    with pytest.raises(ParameterError, match="MAX_MULTI_TYPES"):
        rctd_module.deconvolve(data, mode="multi", max_multi_types=2)


def test_rctd_extract_results_full_mode_with_fake_r(monkeypatch: pytest.MonkeyPatch):
    def _ro_r(code: str):
        if code.strip() == "as.matrix(weights_matrix)":
            return np.array([[0.7, 0.3], [0.2, 0.8]])
        if code.strip() == "cell_type_names":
            return ["A", "B"]
        if code.strip() == "spot_names":
            return ["s1", "s2"]
        return None

    environment = _install_fake_r_modules(monkeypatch, ro_r=_ro_r)
    out = rctd_module._extract_rctd_results("full", environment.robjects)
    assert list(out.index) == ["s1", "s2"]
    assert list(out.columns) == ["A", "B"]
    assert out.shape == (2, 2)


def test_rctd_extract_results_doublet_mode_with_fake_r(monkeypatch: pytest.MonkeyPatch):
    def _ro_r(code: str):
        if code.strip() == "as.matrix(weights_matrix)":
            return np.array([[1.0, 0.0], [0.3, 0.7]])
        if code.strip() == "cell_type_names":
            return ["A", "B"]
        if code.strip() == "spot_names":
            return ["s1", "s2"]
        return None

    environment = _install_fake_r_modules(monkeypatch, ro_r=_ro_r)
    out = rctd_module._extract_rctd_results("doublet", environment.robjects)
    assert out.shape == (2, 2)
    assert np.isclose(out.loc["s1", "A"], 1.0)


def test_rctd_extract_results_multi_mode_with_fake_r(monkeypatch: pytest.MonkeyPatch):
    def _ro_r(code: str):
        if code.strip() == "as.matrix(weights_matrix)":
            return np.array([[0.4, 0.6]])
        if code.strip() == "cell_type_names":
            return ["A", "B"]
        if code.strip() == "spot_names":
            return ["s1"]
        return None

    environment = _install_fake_r_modules(monkeypatch, ro_r=_ro_r)
    out = rctd_module._extract_rctd_results("multi", environment.robjects)
    assert out.shape == (1, 2)
    assert np.isclose(out.loc["s1", "B"], 0.6)


def test_spotlight_wraps_runtime_errors_as_processing_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    def _ro_r(code: str):
        if "library(SPOTlight)" in code:
            raise RuntimeError("R init failed")
        return None

    _install_fake_r_modules(monkeypatch, ro_r=_ro_r)

    with pytest.raises(ProcessingError, match="SPOTlight deconvolution failed"):
        spotlight_module.deconvolve(data)


def test_spotlight_missing_spatial_coords_preserves_data_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = replace(_prepared_data(minimal_spatial_adata), spatial_coords=None)

    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)

    with pytest.raises(DataError, match="requires spatial coordinates"):
        spotlight_module.deconvolve(data)


def test_spotlight_success_casts_counts_and_returns_stats(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.spatial.X = data.spatial.X.astype(np.float32)
    data.reference.X = data.reference.X.astype(np.float64)
    data.reference.obs["cell_type"] = ["A/B"] * (data.reference.n_obs // 2) + [
        "B C"
    ] * (data.reference.n_obs - data.reference.n_obs // 2)

    def _ro_r(code: str):
        text = code.strip()
        if "SPOTlight(" in text:
            ro_mod = __import__("sys").modules["rpy2.robjects"]
            ro_mod.globalenv["last_spotlight_code"] = code
        if text == "spotlight_result$mat":
            return np.array([[0.9, 0.1], [0.4, 0.6]], dtype=float)
        if text == "rownames(spotlight_result$mat)":
            return ["s1", "s2"]
        if text == "colnames(spotlight_result$mat)":
            return ["A/B", "B C"]
        return None

    _install_fake_r_modules(monkeypatch, ro_r=_ro_r)

    proportions, stats = spotlight_module.deconvolve(
        data,
        n_top_genes=1234,
        nmf_model="std",
        min_prop=0.05,
        scale=False,
        weight_id="weight_col",
    )

    import rpy2.robjects as ro

    assert ro.globalenv["spatial_counts"].dtype == np.int32
    assert ro.globalenv["reference_counts"].dtype == np.int32
    assert ro.globalenv["cell_types"][0] == "A/B"
    assert ro.globalenv["cell_types"][-1] == "B C"
    assert ro.globalenv["n_top_genes"] == 1234
    assert "verbose = FALSE" in ro.globalenv["last_spotlight_code"]
    assert proportions.shape == (2, 2)
    assert list(proportions.columns) == ["A/B", "B C"]
    assert stats["method"] == "SPOTlight"
    assert stats["n_top_genes"] == 1234
    assert stats["nmf_model"] == "std"
    assert stats["min_prop"] == pytest.approx(0.05)


def test_spotlight_passthrough_processing_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)
    monkeypatch.setattr(
        spotlight_module,
        "to_dense",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            ProcessingError("dense failed")
        ),
    )

    with pytest.raises(ProcessingError, match="dense failed"):
        spotlight_module.deconvolve(data)


def test_rctd_deconvolve_filters_rare_types_and_raises_when_insufficient_types(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs["cell_type"] = ["A"] * (data.reference.n_obs - 1) + ["B"]

    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)

    with pytest.warns(UserWarning, match="Filtering 1 rare types"):
        with pytest.raises(DataError, match="RCTD requires at least 2 cell types"):
            rctd_module.deconvolve(data, mode="full")


def test_rctd_deconvolve_raises_for_negative_proportions(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)

    bad = pd.DataFrame(
        {"A": [0.5, -0.1], "B": [0.5, 1.1]},
        index=["s1", "s2"],
    )
    monkeypatch.setattr(
        rctd_module, "_extract_rctd_results", lambda _mode, _robjects: bad
    )

    with pytest.raises(ProcessingError, match="negative values"):
        rctd_module.deconvolve(data, mode="full")


def test_rctd_preserves_distinct_cell_type_labels_supported_by_r(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs["cell_type"] = ["A/B"] * 30 + ["A B"] * 30
    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)
    monkeypatch.setattr(
        rctd_module,
        "_extract_rctd_results",
        lambda _mode, _robjects: pd.DataFrame(
            np.tile([0.6, 0.4], (data.spatial.n_obs, 1)),
            index=data.spatial.obs_names,
            columns=["A/B", "A B"],
        ),
    )

    proportions, _stats = rctd_module.deconvolve(data, mode="full")

    import rpy2.robjects as ro

    assert set(ro.globalenv["cell_types_vec"]) == {"A/B", "A B"}
    assert proportions.columns.tolist() == ["A/B", "A B"]


def test_rctd_deconvolve_warns_on_nan_and_returns_stats(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)

    out_df = pd.DataFrame(
        {"A": [0.7, np.nan], "B": [0.3, 0.6]},
        index=["s1", "s2"],
    )
    monkeypatch.setattr(
        rctd_module, "_extract_rctd_results", lambda _mode, _robjects: out_df
    )

    with pytest.warns(UserWarning, match="NaN values"):
        proportions, stats = rctd_module.deconvolve(data, mode="full")

    assert proportions.shape == (2, 2)
    assert stats["method"] == "RCTD-full"


def test_rctd_deconvolve_without_spatial_coords_raises_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """RCTD must fail-fast when spatial coordinates are missing."""
    data = replace(_prepared_data(minimal_spatial_adata), spatial_coords=None)

    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)

    with pytest.raises(DataError, match="RCTD requires real spatial coordinates"):
        rctd_module.deconvolve(data, mode="full")


def test_card_wraps_runtime_errors_as_processing_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    def _ro_r(code: str):
        if "library(CARD)" in code:
            raise RuntimeError("CARD load failed")
        return None

    _install_fake_r_modules(monkeypatch, ro_r=_ro_r)

    with pytest.raises(ProcessingError, match="CARD deconvolution failed"):
        card_module.deconvolve(data)


def test_card_success_with_fake_r_outputs(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    def _ro_r(code: str):
        text = code.strip()
        if text == "rownames(CARD_obj@Proportion_CARD)":
            return ["s1", "s2"]
        if text == "colnames(CARD_obj@Proportion_CARD)":
            return ["A", "B"]
        if text == "CARD_obj@Proportion_CARD":
            return np.array([[0.7, 0.3], [0.2, 0.8]])
        return None

    _install_fake_r_modules(monkeypatch, ro_r=_ro_r)

    proportions, stats = card_module.deconvolve(data)

    assert proportions.shape == (2, 2)
    assert list(proportions.columns) == ["A", "B"]
    assert stats["method"] == "CARD"
    assert stats["device"] == "CPU"
    assert stats["common_genes"] == len(data.common_genes)


def test_card_raises_error_without_spatial_coords(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """CARD must fail-fast when spatial coordinates are missing."""
    data = replace(_prepared_data(minimal_spatial_adata), spatial_coords=None)

    _install_fake_r_modules(monkeypatch, ro_r=lambda _code: None)

    with pytest.raises(DataError, match="CARD requires real spatial coordinates"):
        card_module.deconvolve(data)


def test_card_uses_reference_sample_info(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs["sample_id"] = [f"s{i % 3}" for i in range(data.reference.n_obs)]

    def _ro_r(code: str):
        text = code.strip()
        if text == "rownames(CARD_obj@Proportion_CARD)":
            return ["s1", "s2"]
        if text == "colnames(CARD_obj@Proportion_CARD)":
            return ["A", "B"]
        if text == "CARD_obj@Proportion_CARD":
            return np.array([[0.7, 0.3], [0.2, 0.8]])
        return None

    _install_fake_r_modules(monkeypatch, ro_r=_ro_r)

    proportions, _stats = card_module.deconvolve(data, sample_key="sample_id")

    import rpy2.robjects as ro

    sc_meta = ro.globalenv["sc_meta"]
    assert sc_meta["sampleInfo"].tolist() == data.reference.obs["sample_id"].tolist()
    assert proportions.shape == (2, 2)


def test_card_re_raises_processing_error_without_wrapping(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    def _ro_r(code: str):
        text = code.strip()
        if text == "rownames(CARD_obj@Proportion_CARD)":
            return ["s1", "s2"]
        if text == "colnames(CARD_obj@Proportion_CARD)":
            return ["A", "B"]
        if text == "CARD_obj@Proportion_CARD":
            return np.array([[0.7, 0.3], [0.2, 0.8]])
        return None

    _install_fake_r_modules(monkeypatch, ro_r=_ro_r)
    monkeypatch.setattr(
        card_module,
        "create_deconvolution_stats",
        lambda *_a, **_k: (_ for _ in ()).throw(ProcessingError("stats failed")),
    )

    with pytest.raises(ProcessingError, match="stats failed"):
        card_module.deconvolve(data)


def test_card_success_with_imputation_adds_imputation_statistics(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    def _ro_r(code: str):
        text = code.strip()
        if text == "rownames(CARD_obj@Proportion_CARD)":
            return ["s1", "s2"]
        if text == "colnames(CARD_obj@Proportion_CARD)":
            return ["A", "B"]
        if text == "CARD_obj@Proportion_CARD":
            return np.array([[0.6, 0.4], [0.1, 0.9]])
        if text == "rownames(CARD_impute@refined_prop)":
            return ["1x2", "3x4"]
        if text == "colnames(CARD_impute@refined_prop)":
            return ["A", "B"]
        if text == "CARD_impute@refined_prop":
            return np.array([[0.5, 0.5], [0.8, 0.2]])
        return None

    _install_fake_r_modules(monkeypatch, ro_r=_ro_r)

    proportions, stats = card_module.deconvolve(
        data,
        imputation=True,
        NumGrids=500,
        ineibor=5,
    )

    assert proportions.shape == (2, 2)
    assert "imputation" in stats
    assert stats["imputation"]["enabled"] is True
    assert stats["imputation"]["n_imputed_locations"] == 2
    assert stats["imputation"]["resolution_increase"] == "1.0x"
