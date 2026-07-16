"""Unit contracts for deconvolution method submodules (lightweight, mocked)."""

from __future__ import annotations

from contextlib import nullcontext
from dataclasses import replace
from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import torch

from chatspatial.tools.deconvolution import cell2location as c2l_module
from chatspatial.tools.deconvolution import destvi as destvi_module
from chatspatial.tools.deconvolution import flashdeconv as flash_module
from chatspatial.tools.deconvolution import stereoscope as stereo_module
from chatspatial.tools.deconvolution import tangram as tangram_module
from chatspatial.tools.deconvolution.base import PreparedDeconvolutionData
from chatspatial.utils.exceptions import DataError, DependencyError, ProcessingError


class DummyCtx:
    def __init__(self):
        self.warnings: list[str] = []

    async def warning(self, msg: str):
        self.warnings.append(msg)


def _prepared_data(minimal_spatial_adata, *, n_types: int = 2) -> PreparedDeconvolutionData:
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()
    labels = ["A", "B"] if n_types == 2 else ["A"]
    reference.obs["cell_type"] = [labels[i % len(labels)] for i in range(reference.n_obs)]
    return PreparedDeconvolutionData(
        spatial=spatial,
        reference=reference,
        cell_type_key="cell_type",
        cell_types=labels,
        common_genes=list(spatial.var_names),
        spatial_coords=spatial.obsm["spatial"],
        ctx=DummyCtx(),
    )


def _patch_required_dependency(
    monkeypatch: pytest.MonkeyPatch,
    target_module: object,
    dependency_name: str,
    dependency: object,
) -> None:
    def _require(requested_name: str, *_args, **_kwargs):
        if requested_name != dependency_name:
            raise AssertionError(f"Unexpected dependency request: {requested_name}")
        return dependency

    monkeypatch.setattr(target_module, "require", _require)


def _patch_destvi_dependencies(
    monkeypatch: pytest.MonkeyPatch,
    scvi_module: object,
) -> None:
    dependencies = {"scvi": scvi_module, "torch": torch}

    def _require(requested_name: str, *_args, **_kwargs):
        try:
            return dependencies[requested_name]
        except KeyError as exc:
            raise AssertionError(
                f"Unexpected dependency request: {requested_name}"
            ) from exc

    monkeypatch.setattr(destvi_module, "require", _require)


def _patch_required_module(
    monkeypatch: pytest.MonkeyPatch,
    target_module: object,
    module_name: str,
    module: object,
) -> None:
    def _require_module(
        _dependency_name: str,
        requested_module: str,
        *_args,
        **_kwargs,
    ):
        if requested_module != module_name:
            raise AssertionError(f"Unexpected module request: {requested_module}")
        return module

    monkeypatch.setattr(target_module, "require_module", _require_module)


def test_flashdeconv_dependency_error(minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch):
    data = _prepared_data(minimal_spatial_adata)
    monkeypatch.setattr(
        flash_module,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("FlashDeconv is not available")
        ),
    )
    with pytest.raises(DependencyError, match="FlashDeconv is not available"):
        flash_module.deconvolve(data)


def test_flashdeconv_success_with_fake_backend(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    def _fake_run(adata_st, _reference, **_kwargs):
        adata_st.obsm["flashdeconv"] = np.tile(np.array([0.7, 0.3]), (adata_st.n_obs, 1))

    fake_mod = ModuleType("flashdeconv")
    fake_mod.tl = SimpleNamespace(deconvolve=_fake_run)
    _patch_required_dependency(monkeypatch, flash_module, "flashdeconv", fake_mod)

    proportions, stats = flash_module.deconvolve(data)
    assert proportions.shape == (data.n_spots, 2)
    assert list(proportions.columns) == data.cell_types
    assert stats["method"] == "FlashDeconv"


def test_flashdeconv_missing_output_raises_processing_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    fake_mod = ModuleType("flashdeconv")
    fake_mod.tl = SimpleNamespace(deconvolve=lambda *_args, **_kwargs: None)
    _patch_required_dependency(monkeypatch, flash_module, "flashdeconv", fake_mod)

    with pytest.raises(ProcessingError, match="did not produce output"):
        flash_module.deconvolve(data)


def test_flashdeconv_dataframe_output_uses_spatial_observation_labels(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    def _fake_run(adata_st, _reference, **_kwargs):
        adata_st.obsm["flashdeconv"] = pd.DataFrame(
            np.tile([0.6, 0.4], (adata_st.n_obs, 1)),
            index=adata_st.obs_names,
            columns=data.cell_types,
        )

    fake_mod = ModuleType("flashdeconv")
    fake_mod.tl = SimpleNamespace(deconvolve=_fake_run)
    _patch_required_dependency(monkeypatch, flash_module, "flashdeconv", fake_mod)

    proportions, _stats = flash_module.deconvolve(data)
    assert list(proportions.index) == list(data.spatial.obs_names)


def test_flashdeconv_wraps_unexpected_backend_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    def _fake_run(*_args, **_kwargs):
        raise RuntimeError("backend crashed")

    fake_mod = ModuleType("flashdeconv")
    fake_mod.tl = SimpleNamespace(deconvolve=_fake_run)
    _patch_required_dependency(monkeypatch, flash_module, "flashdeconv", fake_mod)

    with pytest.raises(ProcessingError, match="FlashDeconv deconvolution failed"):
        flash_module.deconvolve(data)


def test_tangram_dependency_error_is_wrapped(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    monkeypatch.setattr(
        tangram_module,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("tangram-sc is required")
        ),
    )
    with pytest.raises(DependencyError, match="tangram-sc is required"):
        tangram_module.deconvolve(data)


def test_tangram_success_with_fake_module(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    class _Map:
        def __init__(self, X, labels):
            self.X = X
            self.obs = pd.DataFrame({data.cell_type_key: labels})

    def _pp_adatas(_ref, _spatial, genes):
        assert len(genes) > 0

    def _map_cells_to_space(_ref_data, spatial_data, **kwargs):
        assert kwargs["mode"] == "clusters"
        return _Map(
            np.tile(np.array([[0.8], [0.2]]), (1, spatial_data.n_obs)).astype(float),
            ["B", "A"],
        )

    fake_mod = ModuleType("tangram")
    fake_mod.pp_adatas = _pp_adatas
    fake_mod.map_cells_to_space = _map_cells_to_space
    _patch_required_dependency(monkeypatch, tangram_module, "tangram", fake_mod)

    proportions, stats = tangram_module.deconvolve(data, mode="clusters", n_epochs=5)
    assert proportions.shape == (data.n_spots, 2)
    assert list(proportions.columns) == ["B", "A"]
    assert proportions.iloc[0].to_dict() == pytest.approx({"B": 0.8, "A": 0.2})
    assert np.allclose(proportions.sum(axis=1).values, 1.0)
    assert stats["method"] == "Tangram"


def test_tangram_cells_mode_uses_fallback_cell_type_column(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs["ct"] = data.reference.obs["cell_type"].copy()
    data.reference.obs = data.reference.obs.drop(columns=["cell_type"])
    data = replace(data, cell_type_key="ct")

    class _Map:
        def __init__(self, X):
            self.X = X

    def _pp_adatas(_ref, _spatial, genes):
        assert len(genes) > 0

    def _map_cells_to_space(ref_data, spatial_data, **kwargs):
        assert kwargs["mode"] == "cells"
        assert "cluster_label" not in kwargs
        assert "cell_type" in ref_data.obs.columns
        return _Map(np.ones((ref_data.n_obs, spatial_data.n_obs), dtype=float))

    fake_mod = ModuleType("tangram")
    fake_mod.pp_adatas = _pp_adatas
    fake_mod.map_cells_to_space = _map_cells_to_space
    _patch_required_dependency(monkeypatch, tangram_module, "tangram", fake_mod)
    monkeypatch.setattr(tangram_module, "get_device", lambda **_kwargs: "cpu")

    proportions, stats = tangram_module.deconvolve(data, mode="cells", n_epochs=3)
    assert proportions.shape == (data.n_spots, 2)
    assert set(proportions.columns) == {"A", "B"}
    assert stats["device"] == "CPU"


def test_tangram_cluster_mode_rejects_cell_level_mapping_shape(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    class _Map:
        def __init__(self, X):
            self.X = X

    fake_mod = ModuleType("tangram")
    fake_mod.pp_adatas = lambda *_a, **_k: None
    fake_mod.map_cells_to_space = lambda ref_data, spatial_data, **_k: _Map(
        np.ones((ref_data.n_obs, spatial_data.n_obs), dtype=float)
    )
    _patch_required_dependency(monkeypatch, tangram_module, "tangram", fake_mod)

    with pytest.raises(ProcessingError, match="Unexpected Tangram cluster mapping shape"):
        tangram_module.deconvolve(data, mode="clusters", n_epochs=5)


def test_tangram_wraps_unexpected_errors(minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch):
    data = _prepared_data(minimal_spatial_adata)

    def _pp_adatas(_ref, _spatial, genes):
        assert len(genes) > 0

    def _map_cells_to_space(*_args, **_kwargs):
        raise RuntimeError("bad mapping")

    fake_mod = ModuleType("tangram")
    fake_mod.pp_adatas = _pp_adatas
    fake_mod.map_cells_to_space = _map_cells_to_space
    _patch_required_dependency(monkeypatch, tangram_module, "tangram", fake_mod)

    with pytest.raises(ProcessingError, match="Tangram deconvolution failed"):
        tangram_module.deconvolve(data, mode="cells")


def test_tangram_preserves_processing_error_from_internal_checks(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    class _Map:
        def __init__(self, X):
            self.X = X
            self.obs = pd.DataFrame({data.cell_type_key: ["A", "B"]})

    fake_mod = ModuleType("tangram")
    fake_mod.pp_adatas = lambda *_a, **_k: None
    fake_mod.map_cells_to_space = lambda _ref_data, spatial_data, **_k: _Map(
        np.ones((2, spatial_data.n_obs), dtype=float)
    )
    _patch_required_dependency(monkeypatch, tangram_module, "tangram", fake_mod)
    monkeypatch.setattr(
        tangram_module,
        "create_deconvolution_stats",
        lambda *_a, **_k: (_ for _ in ()).throw(ProcessingError("stats failed tangram")),
    )

    with pytest.raises(ProcessingError, match="stats failed tangram"):
        tangram_module.deconvolve(data, mode="clusters")


def test_destvi_dependency_error(minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch):
    data = _prepared_data(minimal_spatial_adata)
    monkeypatch.setattr(
        destvi_module,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("scvi is required for DestVI deconvolution")
        ),
    )
    with pytest.raises(DependencyError, match="scvi is required"):
        destvi_module.deconvolve(data)


def test_destvi_validates_minimum_cell_types(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata, n_types=1)
    monkeypatch.setattr(
        destvi_module,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("Dependency should not load before data validation")
        ),
    )

    with pytest.raises(DataError, match="at least 2 cell types"):
        destvi_module.deconvolve(data)


def test_destvi_success_with_fake_scvi(minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch):
    data = _prepared_data(minimal_spatial_adata)
    data.spatial.X = data.spatial.X.astype(np.float64)
    data.reference.X = data.reference.X.astype(np.int64)

    train_calls: list[dict[str, object]] = []

    class FakeCondSCVI:
        @staticmethod
        def setup_anndata(ref, **_kwargs):
            assert ref.X.dtype == np.float32

        def __init__(self, ref, **kwargs):
            assert ref.X.dtype == np.float32
            assert kwargs["prior"] == "normal"
            self.history = None

        def train(self, **kwargs):
            train_calls.append({"stage": "condscvi", **kwargs})

        def get_vamp_prior(self, *_args, **_kwargs):
            return {
                "mean_vprior": np.zeros((2, 3, 5), dtype=np.float32),
                "var_vprior": np.ones((2, 3, 5), dtype=np.float32),
                "weights_vprior": np.full((2, 3), 1 / 3, dtype=np.float32),
            }

    class FakeDestVI:
        @staticmethod
        def setup_anndata(spatial, **_kwargs):
            assert spatial.X.dtype == np.float32

        @classmethod
        def from_rna_model(cls, spatial_data, cond_model, **kwargs):
            assert spatial_data.X.dtype == np.float32
            assert "l1_reg" not in kwargs
            assert kwargs["vamp_prior_p"] == 15
            prior = cond_model.get_vamp_prior()
            assert all(isinstance(value, torch.Tensor) for value in prior.values())
            assert all(value.dtype == torch.float32 for value in prior.values())
            inst = cls()
            inst._n = spatial_data.n_obs
            return inst

        def train(self, **kwargs):
            train_calls.append({"stage": "destvi", **kwargs})

        def get_proportions(self):
            return pd.DataFrame(
                np.tile([0.55, 0.45], (self._n, 1)),
                columns=["A", "B"],
            )

    fake_scvi = ModuleType("scvi")
    fake_scvi.model = SimpleNamespace(CondSCVI=FakeCondSCVI, DestVI=FakeDestVI)
    _patch_destvi_dependencies(monkeypatch, fake_scvi)

    proportions, stats = destvi_module.deconvolve(data, n_epochs=20, use_gpu=False)
    assert proportions.shape == (data.n_spots, 2)
    assert stats["method"] == "DestVI"
    assert stats["n_cell_types"] == 2

    condscvi_call = next(c for c in train_calls if c["stage"] == "condscvi")
    destvi_call = next(c for c in train_calls if c["stage"] == "destvi")
    assert condscvi_call["plan_kwargs"] == {"lr": 1e-3}
    assert destvi_call["plan_kwargs"] == {"lr": 1e-3, "ct_sparsity_weight": 10.0}


def test_destvi_raises_processing_error_for_invalid_proportion_shape(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    class FakeCondSCVI:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, _ref, **_kwargs):
            return None

        def train(self, **_kwargs):
            return None

    class FakeDestVI:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        @classmethod
        def from_rna_model(cls, _spatial_data, _cond_model, **_kwargs):
            return cls()

        def train(self, **_kwargs):
            return None

        def get_proportions(self):
            return pd.DataFrame()

    fake_scvi = ModuleType("scvi")
    fake_scvi.model = SimpleNamespace(CondSCVI=FakeCondSCVI, DestVI=FakeDestVI)
    _patch_destvi_dependencies(monkeypatch, fake_scvi)

    with pytest.raises(ProcessingError, match="Failed to extract valid proportions from DestVI"):
        destvi_module.deconvolve(data)


def test_destvi_wraps_unexpected_runtime_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    class FakeCondSCVI:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, _ref, **_kwargs):
            return None

        def train(self, **_kwargs):
            raise RuntimeError("condscvi failed")

    class FakeDestVI:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        @classmethod
        def from_rna_model(cls, _spatial_data, _cond_model, **_kwargs):
            return cls()

        def train(self, **_kwargs):
            return None

        def get_proportions(self):
            return pd.DataFrame()

    fake_scvi = ModuleType("scvi")
    fake_scvi.model = SimpleNamespace(CondSCVI=FakeCondSCVI, DestVI=FakeDestVI)
    _patch_destvi_dependencies(monkeypatch, fake_scvi)

    with pytest.raises(ProcessingError, match="DestVI deconvolution failed: condscvi failed"):
        destvi_module.deconvolve(data)


def test_stereoscope_success_with_fake_scvi_external(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs["cell_type"] = data.reference.obs["cell_type"].astype(str)

    class FakeRNAStereoscope:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, _ref):
            return None

        def train(self, **_kwargs):
            return None

    class FakeSpatialStereoscope:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        @classmethod
        def from_rna_model(cls, spatial_data, _rna_model):
            inst = cls()
            inst._n = spatial_data.n_obs
            return inst

        def train(self, **_kwargs):
            return None

        def get_proportions(self):
            return np.tile(np.array([0.8, 0.2]), (self._n, 1))

    external_mod = ModuleType("scvi.external")
    external_mod.RNAStereoscope = FakeRNAStereoscope
    external_mod.SpatialStereoscope = FakeSpatialStereoscope
    _patch_required_module(
        monkeypatch, stereo_module, "scvi.external", external_mod
    )

    proportions, stats = stereo_module.deconvolve(data, n_epochs=100, use_gpu=False)
    assert proportions.shape == (data.n_spots, 2)
    assert stats["method"] == "Stereoscope"
    assert stats["rna_epochs"] == 50
    assert stats["spatial_epochs"] == 50


def test_stereoscope_default_epochs_and_gpu_kwargs(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs["cell_type"] = data.reference.obs["cell_type"].astype(str)
    monkeypatch.setattr(stereo_module, "get_accelerator", lambda **_kwargs: "gpu")

    train_calls: list[dict[str, object]] = []

    class FakeRNAStereoscope:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, _ref):
            return None

        def train(self, **kwargs):
            train_calls.append({"stage": "rna", **kwargs})

    class FakeSpatialStereoscope:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        @classmethod
        def from_rna_model(cls, spatial_data, _rna_model):
            inst = cls()
            inst._n = spatial_data.n_obs
            return inst

        def train(self, **kwargs):
            train_calls.append({"stage": "spatial", **kwargs})

        def get_proportions(self):
            return np.tile(np.array([0.7, 0.3]), (self._n, 1))

    external_mod = ModuleType("scvi.external")
    external_mod.RNAStereoscope = FakeRNAStereoscope
    external_mod.SpatialStereoscope = FakeSpatialStereoscope
    _patch_required_module(
        monkeypatch, stereo_module, "scvi.external", external_mod
    )

    proportions, stats = stereo_module.deconvolve(data, n_epochs=150000, use_gpu=True)
    assert proportions.shape == (data.n_spots, 2)
    assert stats["rna_epochs"] == 75000
    assert stats["spatial_epochs"] == 75000
    assert stats["device"] == "gpu"

    rna_call = next(c for c in train_calls if c["stage"] == "rna")
    spatial_call = next(c for c in train_calls if c["stage"] == "spatial")
    assert rna_call["max_epochs"] == 75000
    assert spatial_call["max_epochs"] == 75000
    assert rna_call["accelerator"] == "gpu"
    assert spatial_call["accelerator"] == "gpu"


def test_stereoscope_wraps_unexpected_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs["cell_type"] = data.reference.obs["cell_type"].astype(str)

    class FakeRNAStereoscope:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, _ref):
            return None

        def train(self, **_kwargs):
            return None

    class FakeSpatialStereoscope:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        @classmethod
        def from_rna_model(cls, _spatial_data, _rna_model):
            raise RuntimeError("spatial model failed")

    external_mod = ModuleType("scvi.external")
    external_mod.RNAStereoscope = FakeRNAStereoscope
    external_mod.SpatialStereoscope = FakeSpatialStereoscope
    _patch_required_module(
        monkeypatch, stereo_module, "scvi.external", external_mod
    )

    with pytest.raises(ProcessingError, match="Stereoscope deconvolution failed"):
        stereo_module.deconvolve(data, n_epochs=10, use_gpu=False)


def test_stereoscope_preserves_processing_error_from_internal_checks(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs["cell_type"] = data.reference.obs["cell_type"].astype(str)

    class FakeRNAStereoscope:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, _ref):
            return None

        def train(self, **_kwargs):
            return None

    class FakeSpatialStereoscope:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        @classmethod
        def from_rna_model(cls, spatial_data, _rna_model):
            inst = cls()
            inst._n = spatial_data.n_obs
            return inst

        def train(self, **_kwargs):
            return None

        def get_proportions(self):
            return np.tile(np.array([0.8, 0.2]), (self._n, 1))

    external_mod = ModuleType("scvi.external")
    external_mod.RNAStereoscope = FakeRNAStereoscope
    external_mod.SpatialStereoscope = FakeSpatialStereoscope
    _patch_required_module(
        monkeypatch, stereo_module, "scvi.external", external_mod
    )
    monkeypatch.setattr(
        stereo_module,
        "create_deconvolution_stats",
        lambda *_a, **_k: (_ for _ in ()).throw(ProcessingError("stats failed")),
    )

    with pytest.raises(ProcessingError, match="stats failed"):
        stereo_module.deconvolve(data, n_epochs=20, use_gpu=False)


@pytest.mark.asyncio
async def test_cell2location_apply_gene_filtering_unavailable_returns_copy(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx()
    monkeypatch.setattr(c2l_module, "is_available", lambda *_: False)

    out = await c2l_module.apply_gene_filtering(adata, ctx)
    assert out is not adata
    assert out.n_obs == adata.n_obs
    assert out.n_vars == adata.n_vars
    assert any("Skipping gene filtering" in m for m in ctx.warnings)


def test_cell2location_build_train_kwargs_variants():
    aggressive = c2l_module._build_train_kwargs(
        epochs=10,
        lr=0.01,
        train_size=0.8,
        accelerator="gpu",
        early_stopping=True,
        early_stopping_patience=5,
        early_stopping_threshold=0.01,
        validation_size=0.2,
        use_aggressive=True,
    )
    assert aggressive["accelerator"] == "gpu"
    assert aggressive["check_val_every_n_epoch"] == 1
    assert aggressive["train_size"] == pytest.approx(0.8)
    assert aggressive["early_stopping_min_delta"] == 0.01

    standard = c2l_module._build_train_kwargs(
        epochs=10,
        lr=0.01,
        train_size=0.7,
        accelerator="cpu",
        early_stopping=False,
        early_stopping_patience=5,
        early_stopping_threshold=0.0,
        validation_size=0.1,
        use_aggressive=False,
    )
    assert standard["max_epochs"] == 10
    assert standard["train_size"] == pytest.approx(0.7)
    assert "accelerator" not in standard


def test_cell2location_build_train_kwargs_standard_gpu_and_aggressive_no_es():
    standard_gpu = c2l_module._build_train_kwargs(
        epochs=8,
        lr=0.02,
        train_size=0.6,
        accelerator="gpu",
        early_stopping=False,
        early_stopping_patience=5,
        early_stopping_threshold=0.0,
        validation_size=0.1,
        use_aggressive=False,
    )
    assert standard_gpu["accelerator"] == "gpu"
    assert standard_gpu["train_size"] == pytest.approx(0.6)

    aggressive_no_es = c2l_module._build_train_kwargs(
        epochs=12,
        lr=0.01,
        train_size=0.75,
        accelerator="cpu",
        early_stopping=False,
        early_stopping_patience=5,
        early_stopping_threshold=0.0,
        validation_size=0.2,
        use_aggressive=True,
    )
    assert aggressive_no_es["train_size"] == pytest.approx(0.75)
    assert "check_val_every_n_epoch" not in aggressive_no_es


def test_cell2location_extract_reference_signatures_and_abundance(minimal_spatial_adata):
    ref = minimal_spatial_adata.copy()
    ref.uns["mod"] = {"factor_names": ["A", "B"]}
    ref.var["means_per_cluster_mu_fg_A"] = np.ones(ref.n_vars)
    ref.var["means_per_cluster_mu_fg_B"] = np.full(ref.n_vars, 2.0)

    sig = c2l_module._extract_reference_signatures(ref)
    assert list(sig.columns) == ["A", "B"]
    assert sig.shape[0] == ref.n_vars

    sp = minimal_spatial_adata.copy()
    sp.obsm["q05_cell_abundance_w_sf"] = pd.DataFrame(
        np.tile([0.4, 0.6], (sp.n_obs, 1)),
        index=sp.obs_names,
        columns=[
            "q05cell_abundance_w_sf_B",
            "q05cell_abundance_w_sf_A",
        ],
    )
    abundance = c2l_module._extract_cell_abundance(sp, ["A", "B"])
    assert abundance.shape == (sp.n_obs, 2)
    assert abundance.columns.tolist() == ["A", "B"]
    np.testing.assert_allclose(abundance.iloc[0], [0.6, 0.4])


def test_cell2location_extract_cell_abundance_missing_key_raises(minimal_spatial_adata):
    with pytest.raises(ProcessingError, match="did not produce expected output"):
        c2l_module._extract_cell_abundance(
            minimal_spatial_adata.copy(), ["A", "B"]
        )


@pytest.mark.asyncio
async def test_cell2location_apply_gene_filtering_uses_filter_and_subsets(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx()

    monkeypatch.setattr(c2l_module, "is_available", lambda *_: True)
    monkeypatch.setattr(c2l_module, "non_interactive_backend", nullcontext)

    filtering_mod = ModuleType("cell2location.utils.filtering")

    def _filter_genes(_adata, **_kwargs):
        mask = np.zeros(_adata.n_vars, dtype=bool)
        mask[:5] = True
        return mask

    filtering_mod.filter_genes = _filter_genes
    _patch_required_module(
        monkeypatch,
        c2l_module,
        "cell2location.utils.filtering",
        filtering_mod,
    )

    out = await c2l_module.apply_gene_filtering(adata, ctx)

    assert out.n_obs == adata.n_obs
    assert out.n_vars == 5


def test_cell2location_extract_reference_signatures_prefers_varm(minimal_spatial_adata):
    ref = minimal_spatial_adata.copy()
    ref.uns["mod"] = {"factor_names": ["A", "B"]}
    ref.varm["means_per_cluster_mu_fg"] = pd.DataFrame(
        {
            "means_per_cluster_mu_fg_A": np.ones(ref.n_vars),
            "means_per_cluster_mu_fg_B": np.full(ref.n_vars, 2.0),
        },
        index=ref.var_names,
    )

    sig = c2l_module._extract_reference_signatures(ref)
    assert sig.shape == (ref.n_vars, 2)
    assert list(sig.columns) == ["A", "B"]


def test_cell2location_extract_cell_abundance_fallback_key(minimal_spatial_adata):
    sp = minimal_spatial_adata.copy()
    sp.obsm["means_cell_abundance_w_sf"] = np.tile([0.1, 0.9], (sp.n_obs, 1))
    abundance = c2l_module._extract_cell_abundance(sp, ["A", "B"])
    assert abundance.shape == (sp.n_obs, 2)


def test_cell2location_extract_cell_abundance_rejects_mislabeled_columns(
    minimal_spatial_adata,
):
    sp = minimal_spatial_adata.copy()
    sp.obsm["q05_cell_abundance_w_sf"] = pd.DataFrame(
        np.ones((sp.n_obs, 2)),
        index=sp.obs_names,
        columns=["unknown_A", "unknown_B"],
    )

    with pytest.raises(ProcessingError, match="columns do not match"):
        c2l_module._extract_cell_abundance(sp, ["A", "B"])


@pytest.mark.parametrize(
    ("values", "message"),
    [
        (lambda n: np.ones((n, 1)), "expected"),
        (lambda n: np.full((n, 2), np.nan), "non-finite"),
        (lambda n: np.full((n, 2), -0.1), "negative"),
    ],
)
def test_cell2location_extract_cell_abundance_rejects_invalid_values(
    minimal_spatial_adata, values, message
):
    sp = minimal_spatial_adata.copy()
    sp.obsm["means_cell_abundance_w_sf"] = values(sp.n_obs)

    with pytest.raises(ProcessingError, match=message):
        c2l_module._extract_cell_abundance(sp, ["A", "B"])


def test_cell2location_extract_cell_abundance_validates_exported_factor_order(
    minimal_spatial_adata,
):
    sp = minimal_spatial_adata.copy()
    sp.uns["mod"] = {"factor_names": ["B", "A"]}
    sp.obsm["means_cell_abundance_w_sf"] = np.ones((sp.n_obs, 2))

    with pytest.raises(ProcessingError, match="factor names do not match"):
        c2l_module._extract_cell_abundance(sp, ["A", "B"])


def test_cell2location_deconvolve_success_with_fake_models(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.obs.loc[data.reference.obs_names[0], "cell_type"] = np.nan

    monkeypatch.setattr(c2l_module, "get_device", lambda **_kwargs: "cpu")
    monkeypatch.setattr(c2l_module, "get_accelerator", lambda **_kwargs: "cpu")
    monkeypatch.setattr(c2l_module, "suppress_output", nullcontext)
    monkeypatch.setattr(
        c2l_module,
        "check_model_convergence",
        lambda *_args, **_kwargs: (False, "model not converged"),
    )

    class _FakeRegressionModel:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, ref):
            self.ref = ref
            self.history = {"elbo_train": pd.Series([1.0, 0.9])}

        def train(self, **_kwargs):
            return None

        def export_posterior(self, ref, **_kwargs):
            ref.uns["mod"] = {"factor_names": ["A", "B"]}
            ref.var["means_per_cluster_mu_fg_A"] = np.ones(ref.n_vars)
            ref.var["means_per_cluster_mu_fg_B"] = np.full(ref.n_vars, 2.0)
            return ref

    class _FakeCell2location:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, sp, **_kwargs):
            self.sp = sp
            self.history = {"elbo_train": pd.Series([9.0, 8.5])}

        def train(self, **_kwargs):
            return None

        def export_posterior(self, sp, **_kwargs):
            sp.obsm["q05_cell_abundance_w_sf"] = pd.DataFrame(
                np.tile([0.4, 0.6], (sp.n_obs, 1)),
                index=sp.obs_names,
                columns=[
                    "q05cell_abundance_w_sf_B",
                    "q05cell_abundance_w_sf_A",
                ],
            )
            return sp

    models_mod = ModuleType("cell2location.models")
    models_mod.Cell2location = _FakeCell2location
    models_mod.RegressionModel = _FakeRegressionModel
    _patch_required_module(
        monkeypatch, c2l_module, "cell2location.models", models_mod
    )

    with pytest.warns(UserWarning) as record:
        proportions, stats = c2l_module.deconvolve(data, n_epochs=5, ref_model_epochs=5)
    assert any("not converged" in str(w.message) for w in record)

    assert proportions.shape == (data.n_spots, 2)
    assert set(proportions.columns) == {"A", "B"}
    np.testing.assert_allclose(proportions.iloc[0], [0.6, 0.4])
    assert stats["method"] == "Cell2location"
    assert stats["device"] == "cpu"
    assert "final_elbo" in stats


def test_cell2location_deconvolve_casts_ref_and_spatial_to_float32(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)
    data.reference.X = data.reference.X.astype(np.int64)
    data.spatial.X = data.spatial.X.astype(np.float64)

    monkeypatch.setattr(c2l_module, "get_device", lambda **_kwargs: "cpu")
    monkeypatch.setattr(c2l_module, "get_accelerator", lambda **_kwargs: "cpu")
    monkeypatch.setattr(c2l_module, "suppress_output", nullcontext)
    monkeypatch.setattr(
        c2l_module,
        "check_model_convergence",
        lambda *_args, **_kwargs: (True, None),
    )

    class _FakeRegressionModel:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, ref):
            assert ref.X.dtype == np.float32
            self.history = {"elbo_train": pd.Series([1.0])}

        def train(self, **_kwargs):
            return None

        def export_posterior(self, ref, **_kwargs):
            ref.uns["mod"] = {"factor_names": ["A", "B"]}
            ref.var["means_per_cluster_mu_fg_A"] = np.ones(ref.n_vars)
            ref.var["means_per_cluster_mu_fg_B"] = np.full(ref.n_vars, 2.0)
            return ref

    class _FakeCell2location:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

        def __init__(self, sp, **_kwargs):
            assert sp.X.dtype == np.float32
            self._n = sp.n_obs
            self.history = {"elbo_train": pd.Series([2.0])}

        def train(self, **_kwargs):
            return None

        def export_posterior(self, sp, **_kwargs):
            sp.obsm["q05_cell_abundance_w_sf"] = np.tile([0.5, 0.5], (sp.n_obs, 1))
            return sp

    models_mod = ModuleType("cell2location.models")
    models_mod.Cell2location = _FakeCell2location
    models_mod.RegressionModel = _FakeRegressionModel
    _patch_required_module(
        monkeypatch, c2l_module, "cell2location.models", models_mod
    )

    proportions, stats = c2l_module.deconvolve(data, n_epochs=2, ref_model_epochs=2)
    assert proportions.shape == (data.n_spots, 2)
    assert stats["method"] == "Cell2location"


def test_cell2location_deconvolve_passthrough_data_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    monkeypatch.setattr(c2l_module, "get_device", lambda **_kwargs: "cpu")
    monkeypatch.setattr(c2l_module, "get_accelerator", lambda **_kwargs: "cpu")

    class _FakeRegressionModel:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            raise DataError("invalid reference")

    class _FakeCell2location:
        @staticmethod
        def setup_anndata(*_args, **_kwargs):
            return None

    models_mod = ModuleType("cell2location.models")
    models_mod.Cell2location = _FakeCell2location
    models_mod.RegressionModel = _FakeRegressionModel
    _patch_required_module(
        monkeypatch, c2l_module, "cell2location.models", models_mod
    )

    with pytest.raises(DataError, match="invalid reference"):
        c2l_module.deconvolve(data)


def test_cell2location_deconvolve_wraps_unexpected_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    data = _prepared_data(minimal_spatial_adata)

    monkeypatch.setattr(
        c2l_module,
        "get_device",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("broken env")),
    )

    models_mod = ModuleType("cell2location.models")
    models_mod.Cell2location = object
    models_mod.RegressionModel = object
    _patch_required_module(
        monkeypatch, c2l_module, "cell2location.models", models_mod
    )

    with pytest.raises(ProcessingError, match="Cell2location deconvolution failed"):
        c2l_module.deconvolve(data)
