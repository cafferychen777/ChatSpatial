"""Unit tests for integration tool utility contracts and regressions."""

from __future__ import annotations

import warnings
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData, ImplicitModificationWarning

from chatspatial.tools.integration import (
    _orient_harmony_embedding,
    _reassemble_scanorama_embeddings,
    integrate_multiple_samples,
    integrate_with_scvi,
    rescale_spatial_coordinates,
)
from chatspatial.utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
    ProcessingError,
)


def _install_fake_scvi(monkeypatch: pytest.MonkeyPatch, calls: dict[str, object]) -> None:
    class FakeSCVI:
        @staticmethod
        def setup_anndata(combined, batch_key, layer):
            calls["setup"] = {
                "batch_key": batch_key,
                "layer": layer,
                "n_obs": combined.n_obs,
            }

        def __init__(self, combined, **kwargs):
            calls["init"] = {"kwargs": kwargs, "n_obs": combined.n_obs}
            self._combined = combined

        def train(self, **kwargs):
            calls["train"] = kwargs

        def get_latent_representation(self):
            return np.zeros((self._combined.n_obs, 3), dtype=np.float32)

    fake_scvi = ModuleType("scvi")
    fake_scvi.model = SimpleNamespace(SCVI=FakeSCVI)

    monkeypatch.setattr("chatspatial.tools.integration.require", lambda *_args, **_kwargs: None)
    monkeypatch.setitem(__import__("sys").modules, "scvi", fake_scvi)


def test_rescale_spatial_coordinates_preserves_row_mapping_for_interleaved_batches(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) % 2 == 0, "A", "B")

    ref_mask = adata.obs["batch"] == "A"
    b_mask = adata.obs["batch"] == "B"
    ref_coords = adata.obsm["spatial"][ref_mask.to_numpy()]
    b_coords = adata.obsm["spatial"][b_mask.to_numpy()]

    mean = ref_coords.mean(axis=0)
    std = ref_coords.std(axis=0)
    expected = np.zeros_like(adata.obsm["spatial"], dtype=float)
    expected[ref_mask.to_numpy()] = (ref_coords - mean) / std
    expected[b_mask.to_numpy()] = (b_coords - mean) / std

    out = rescale_spatial_coordinates(adata, batch_key="batch", reference_batch="A")

    np.testing.assert_allclose(out.obsm["spatial_aligned"], expected, atol=1e-8)


def test_orient_harmony_embedding_accepts_supported_axis_orders():
    row_major = np.arange(12, dtype=np.float32).reshape(4, 3)
    np.testing.assert_array_equal(_orient_harmony_embedding(row_major, 4), row_major)
    np.testing.assert_array_equal(
        _orient_harmony_embedding(row_major.T, 4), row_major
    )


@pytest.mark.parametrize(
    ("embedding", "message"),
    [
        (np.ones(4), "two-dimensional"),
        (np.ones((3, 2)), "does not match the dataset"),
        (np.ones((4, 0)), "zero embedding components"),
        (np.full((4, 2), np.nan), "non-finite"),
        ([["not-numeric"]], "non-numeric"),
    ],
)
def test_orient_harmony_embedding_rejects_invalid_output(embedding, message):
    with pytest.raises(ProcessingError, match=message):
        _orient_harmony_embedding(embedding, 4)


@pytest.mark.integration
def test_integrate_multiple_samples_scvi_metadata_uses_x_scvi_key(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True

    captured: dict[str, object] = {}

    monkeypatch.setattr(
        "chatspatial.tools.integration.check_is_integer_counts",
        lambda *_a, **_k: (False, False, False),
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.integrate_with_scvi",
        lambda combined, **_kwargs: combined,
    )
    monkeypatch.setattr("scanpy.tl.umap", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "chatspatial.tools.integration.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )

    out = integrate_multiple_samples(adata, method="scvi")

    assert out is adata
    assert captured["analysis_name"] == "integration_scvi"
    assert captured["results_keys"] == {"obsm": ["X_scvi"]}


def test_integrate_with_scvi_rejects_unpreprocessed_values(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.X = np.full(adata.X.shape, 80.0, dtype=np.float32)

    _install_fake_scvi(monkeypatch, calls={})

    with pytest.raises(DataError, match="preprocessed"):
        integrate_with_scvi(adata, batch_key="batch")


def test_integrate_with_scvi_requires_batch_key(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()

    _install_fake_scvi(monkeypatch, calls={})

    with pytest.raises(ParameterError, match="Batch key 'batch'"):
        integrate_with_scvi(adata, batch_key="batch")


def test_integrate_with_scvi_requires_multiple_batches(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = "single"

    _install_fake_scvi(monkeypatch, calls={})

    with pytest.raises(DataError, match="at least 2 batches"):
        integrate_with_scvi(adata, batch_key="batch")


def test_integrate_with_scvi_happy_path_sets_latent_and_neighbors(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.layers["counts"] = np.abs(
        np.round(np.asarray(adata.X))
    ).astype(float)
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    calls: dict[str, object] = {}
    _install_fake_scvi(monkeypatch, calls)

    monkeypatch.setattr(
        "chatspatial.tools.integration.get_accelerator",
        lambda prefer_gpu=False: "cpu",
    )
    monkeypatch.setattr(
        "scanpy.pp.neighbors",
        lambda _adata, use_rep: calls.update({"neighbors": {"use_rep": use_rep}}),
    )

    out = integrate_with_scvi(adata, batch_key="batch", n_epochs=7, use_gpu=False)

    assert out is adata
    assert out.obsm["X_scvi"].shape == (adata.n_obs, 3)
    assert calls["setup"]["batch_key"] == "batch"
    assert calls["train"]["max_epochs"] == 7
    assert calls["train"]["accelerator"] == "cpu"
    assert calls["neighbors"]["use_rep"] == "X_scvi"


def test_integrate_with_scvi_uses_counts_layer_when_available(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    """scVI should pass layer='counts' to setup_anndata when counts layer exists."""
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(
        np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2"
    )
    # Normalized X (small values, passes the max_val < 50 check)
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)
    # Add a counts layer with raw integer counts
    adata.layers["counts"] = np.abs(
        np.round(np.asarray(minimal_spatial_adata.X))
    ).astype(float)

    calls: dict[str, object] = {}
    _install_fake_scvi(monkeypatch, calls)

    monkeypatch.setattr(
        "chatspatial.tools.integration.get_accelerator",
        lambda prefer_gpu=False: "cpu",
    )
    monkeypatch.setattr(
        "scanpy.pp.neighbors",
        lambda _adata, use_rep: None,
    )

    out = integrate_with_scvi(adata, batch_key="batch", n_epochs=5, use_gpu=False)

    assert out is adata
    assert calls["setup"]["layer"] == "counts"


def test_integrate_with_scvi_salvages_counts_from_raw(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    """When X is not integer counts and no layers['counts'] exists, but .raw
    has integer counts, ensure_counts_layer should salvage them."""
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(
        np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2"
    )
    # Store raw integer counts in .raw before normalizing.
    # AnnData.raw is created from the current state of adata, so we first
    # ensure X contains integer counts, snapshot .raw, then normalize X.
    adata.X = np.abs(np.round(np.asarray(adata.X))).astype(np.float32)
    adata.raw = adata.copy()
    # Normalized X — non-integer, no counts layer
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)
    assert "counts" not in adata.layers

    calls: dict[str, object] = {}
    _install_fake_scvi(monkeypatch, calls)

    monkeypatch.setattr(
        "chatspatial.tools.integration.get_accelerator",
        lambda prefer_gpu=False: "cpu",
    )
    monkeypatch.setattr(
        "scanpy.pp.neighbors",
        lambda _adata, use_rep: None,
    )

    integrate_with_scvi(adata, batch_key="batch", n_epochs=5, use_gpu=False)

    # ensure_counts_layer should have created layers["counts"] from .raw
    assert "counts" in adata.layers
    assert calls["setup"]["layer"] == "counts"


def test_integrate_with_scvi_errors_when_no_counts_available(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    """scVI should raise DataError when only normalized data is available
    (no counts, no .raw)."""
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(
        np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2"
    )
    # Normalized X — non-integer, no counts layer, no .raw
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    calls: dict[str, object] = {}
    _install_fake_scvi(monkeypatch, calls)

    monkeypatch.setattr(
        "chatspatial.tools.integration.get_accelerator",
        lambda prefer_gpu=False: "cpu",
    )
    monkeypatch.setattr(
        "scanpy.pp.neighbors",
        lambda _adata, use_rep: None,
    )

    with pytest.raises(DataError, match="scVI requires raw count data"):
        integrate_with_scvi(adata, batch_key="batch", n_epochs=5, use_gpu=False)


def test_integrate_with_scvi_uses_x_directly_when_x_is_counts(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    """scVI should use layer=None when X itself contains integer counts."""
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(
        np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2"
    )
    # X is integer counts and small enough to pass max_val check
    adata.X = np.clip(
        np.abs(np.round(np.asarray(adata.X))).astype(np.float32), 0, 40
    )

    calls: dict[str, object] = {}
    _install_fake_scvi(monkeypatch, calls)

    monkeypatch.setattr(
        "chatspatial.tools.integration.get_accelerator",
        lambda prefer_gpu=False: "cpu",
    )
    monkeypatch.setattr(
        "scanpy.pp.neighbors",
        lambda _adata, use_rep: None,
    )

    out = integrate_with_scvi(adata, batch_key="batch", n_epochs=5, use_gpu=False)

    assert out is adata
    assert calls["setup"]["layer"] is None


def test_integrate_multiple_samples_autofills_missing_batch_labels_with_concat(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata1 = minimal_spatial_adata.copy()
    adata2 = minimal_spatial_adata.copy()
    adata2.obs_names = [f"b_{i}" for i in range(adata2.n_obs)]
    adata1.var_names = [f"a_{g}" for g in adata1.var_names]
    adata2.var_names = [f"b_{g}" for g in adata2.var_names]

    captured: dict[str, object] = {}

    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.integrate_with_scvi",
        lambda combined, **_kwargs: combined,
    )
    monkeypatch.setattr("scanpy.tl.umap", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "chatspatial.tools.integration.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )

    out = integrate_multiple_samples([adata1, adata2], method="scvi", batch_key="batch")

    assert out.n_obs == adata1.n_obs + adata2.n_obs
    assert set(out.obs["batch"].astype(str).unique()) == {"sample_0", "sample_1"}
    assert captured["analysis_name"] == "integration_scvi"


def test_integrate_multiple_samples_requires_at_least_two_datasets(minimal_spatial_adata):
    with pytest.raises(ParameterError, match="at least 2 datasets"):
        integrate_multiple_samples([minimal_spatial_adata.copy()], method="harmony")


def test_integrate_multiple_samples_requires_batch_key_for_merged_adata(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    if "batch" in adata.obs:
        del adata.obs["batch"]

    with pytest.raises(ParameterError, match="missing batch information key 'batch'"):
        integrate_multiple_samples(adata, method="harmony", batch_key="batch")


def test_integrate_with_scvi_auto_epochs_for_small_dataset(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.layers["counts"] = np.abs(
        np.round(np.asarray(adata.X))
    ).astype(float)
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    calls: dict[str, object] = {}
    _install_fake_scvi(monkeypatch, calls)

    monkeypatch.setattr(
        "chatspatial.tools.integration.get_accelerator",
        lambda prefer_gpu=False: "cpu",
    )
    monkeypatch.setattr("scanpy.pp.neighbors", lambda *_args, **_kwargs: None)

    out = integrate_with_scvi(adata, batch_key="batch", n_epochs=None, use_gpu=False)

    assert out is adata
    assert calls["train"]["max_epochs"] == 400


@pytest.mark.asyncio
async def test_integrate_samples_scvi_exports_and_adds_dataset(monkeypatch: pytest.MonkeyPatch):
    from chatspatial.models.data import IntegrationParameters
    from chatspatial.tools.integration import integrate_samples

    adata = AnnData(np.ones((12, 8), dtype=np.float32))
    adata.obs["batch"] = ["a"] * 6 + ["b"] * 6
    adata.obsm["spatial"] = np.ones((12, 2), dtype=float)

    class _Ctx:
        def __init__(self):
            self.added: list[tuple[str, AnnData]] = []

        async def get_adata(self, data_id: str):
            del data_id
            return adata

        async def add_dataset(self, ds: AnnData, prefix: str = "data"):
            data_id = f"{prefix}_1"
            self.added.append((data_id, ds))
            return data_id

    ctx = _Ctx()
    exported: list[tuple[str, str]] = []

    monkeypatch.setattr(
        "chatspatial.tools.integration.integrate_multiple_samples",
        lambda adatas, **kwargs: adatas[0],
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.rescale_spatial_coordinates",
        lambda combined_adata, **kwargs: (combined_adata.obsm.__setitem__("spatial_aligned", combined_adata.obsm["spatial"]), combined_adata)[1],
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.export_analysis_result",
        lambda adata_obj, dataset_id, analysis_name: exported.append((dataset_id, analysis_name)),
    )

    params = IntegrationParameters(method="scvi", align_spatial=True, batch_key="batch")
    result = await integrate_samples(["d1", "d2"], ctx, params)

    assert result.data_id == "integrated_1"
    assert result.integration_method == "scvi"
    assert len(ctx.added) == 1
    assert exported == [
        ("integrated_1", "integration_scvi"),
        ("integrated_1", "spatial_alignment"),
    ]


def _install_classical_integration_mocks(
    monkeypatch: pytest.MonkeyPatch,
    captured: dict[str, object],
) -> None:
    monkeypatch.setattr(
        "chatspatial.tools.integration.check_is_integer_counts",
        lambda *_a, **_k: (False, False, False),
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    def _set_hvg_true(adata, **_kwargs):
        adata.var.loc[:, "highly_variable"] = True

    monkeypatch.setattr("scanpy.pp.highly_variable_genes", _set_hvg_true)
    monkeypatch.setattr("scanpy.pp.scale", lambda *_args, **_kwargs: None)

    def _fake_pca(adata, n_comps, svd_solver, zero_center=False):
        del svd_solver, zero_center
        _set_obsm(adata, "X_pca", np.zeros((adata.n_obs, n_comps), dtype=np.float32))

    monkeypatch.setattr("scanpy.tl.pca", _fake_pca)
    monkeypatch.setattr("scanpy.tl.umap", lambda *_args, **_kwargs: None)

    def _fake_neighbors(adata, **kwargs):
        captured.setdefault("neighbors", []).append(kwargs)
        adata.uns["neighbors_called"] = True

    monkeypatch.setattr("scanpy.pp.neighbors", _fake_neighbors)
    monkeypatch.setattr(
        "chatspatial.tools.integration.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.require",
        lambda *_args, **_kwargs: None,
    )


def _set_hvg_value(adata, value: bool) -> None:
    adata.var.loc[:, "highly_variable"] = bool(value)


def _set_obsm(adata, key: str, value) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ImplicitModificationWarning)
        adata.obsm[key] = value


def _set_x(adata, value) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.simplefilter("ignore", category=ImplicitModificationWarning)
        adata.X = value



def test_integrate_multiple_samples_harmony_uses_transposed_output_when_needed(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)

    fake_harmonypy = ModuleType("harmonypy")
    # Old harmonypy shape: (n_pcs, n_cells) -> should be transposed by code
    fake_harmonypy.run_harmony = lambda **_kwargs: SimpleNamespace(
        Z_corr=np.full((3, adata.n_obs), 1.0, dtype=np.float32)
    )
    monkeypatch.setitem(__import__("sys").modules, "harmonypy", fake_harmonypy)

    out = integrate_multiple_samples(adata, method="harmony", batch_key="batch", n_pcs=3)

    assert out.n_obs == adata.n_obs
    assert out.obsm["X_pca_harmony"].shape == (adata.n_obs, 3)
    assert captured["analysis_name"] == "integration_harmony"
    assert captured["results_keys"] == {"obsm": ["X_pca_harmony"]}
    assert captured["neighbors"][0]["use_rep"] == "X_pca_harmony"



def test_integrate_multiple_samples_bbknn_dispatches_module_call(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)

    fake_bbknn = ModuleType("bbknn")
    fake_bbknn.bbknn = lambda adata, batch_key, neighbors_within_batch: captured.update(
        {
            "bbknn": {
                "batch_key": batch_key,
                "neighbors_within_batch": neighbors_within_batch,
                "n_obs": adata.n_obs,
            }
        }
    )
    monkeypatch.setitem(__import__("sys").modules, "bbknn", fake_bbknn)

    out = integrate_multiple_samples(adata, method="bbknn", batch_key="batch", n_pcs=3)

    assert out.n_obs == adata.n_obs
    assert captured["bbknn"] == {
        "batch_key": "batch",
        "neighbors_within_batch": 3,
        "n_obs": adata.n_obs,
    }
    assert captured["analysis_name"] == "integration_bbknn"
    assert captured["results_keys"] == {}



def test_integrate_multiple_samples_scanorama_wrapper_path_uses_x_scanorama(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)

    def _scanorama_integrate(combined, key, basis, adjusted_basis):
        del key, basis
        _set_obsm(combined, adjusted_basis, np.zeros((combined.n_obs, 2), dtype=np.float32))

    monkeypatch.setattr("scanpy.external.pp.scanorama_integrate", _scanorama_integrate)

    out = integrate_multiple_samples(adata, method="scanorama", batch_key="batch", n_pcs=3)

    assert out.n_obs == adata.n_obs
    assert "X_scanorama" in out.obsm
    assert captured["analysis_name"] == "integration_scanorama"
    assert captured["results_keys"] == {"obsm": ["X_scanorama"]}
    assert captured["neighbors"][0]["use_rep"] == "X_scanorama"



def test_integrate_multiple_samples_unknown_method_falls_back_to_default_neighbors(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)

    out = integrate_multiple_samples(adata, method="not_real", batch_key="batch", n_pcs=3)

    assert out.n_obs == adata.n_obs
    assert captured["analysis_name"] == "integration_not_real"
    assert captured["results_keys"] == {"obsm": ["X_pca"]}
    # default neighbors call should not force use_rep
    assert captured["neighbors"][0] == {}



def test_integrate_multiple_samples_wraps_harmony_errors(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)

    fake_harmonypy = ModuleType("harmonypy")
    fake_harmonypy.run_harmony = lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("hboom"))
    monkeypatch.setitem(__import__("sys").modules, "harmonypy", fake_harmonypy)

    with pytest.raises(ProcessingError, match="Harmony integration failed"):
        integrate_multiple_samples(adata, method="harmony", batch_key="batch", n_pcs=3)


def test_integrate_multiple_samples_cleans_var_na_and_diffmap_artifacts(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ad1 = minimal_spatial_adata[:, :4].copy()
    ad2 = minimal_spatial_adata[:, 2:6].copy()
    ad1.obs["batch"] = "b1"
    ad2.obs["batch"] = "b2"

    ad1.var["flag"] = [True, False, True, False]
    ad2.var["flag"] = [True, True, False, False]
    ad1.var["symbol"] = ["A", "B", "C", "D"]
    ad2.var["symbol"] = ["C", "D", "E", "F"]
    ad1.obs_names = [f"{name}_b1" for name in ad1.obs_names]
    ad2.obs_names = [f"{name}_b2" for name in ad2.obs_names]

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)
    monkeypatch.setattr("scanpy.tl.umap", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "scanpy.pp.highly_variable_genes",
        lambda adata, **_kwargs: _set_hvg_value(adata, True),
    )

    def _fake_pca(adata, n_comps, svd_solver, zero_center=False):
        del svd_solver, zero_center
        _set_obsm(adata, "X_pca", np.zeros((adata.n_obs, min(n_comps, 3)), dtype=np.float32))

    monkeypatch.setattr("scanpy.tl.pca", _fake_pca)
    monkeypatch.setattr("scanpy.pp.neighbors", lambda *_args, **_kwargs: None)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Downcasting object dtype arrays on .fillna",
            category=FutureWarning,
        )
        out = integrate_multiple_samples(
            [ad1, ad2], method="not_real", batch_key="batch", n_pcs=3
        )

    assert out.var["flag"].dtype == bool
    assert out.var["symbol"].dtype == object
    assert not out.var["symbol"].isna().any()
    assert "X_diffmap" not in out.obsm
    assert "diffmap_evals" not in out.uns


def test_integrate_multiple_samples_raises_when_hvg_recalc_still_empty(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = False
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "scanpy.pp.highly_variable_genes",
        lambda adata, **_kwargs: _set_hvg_value(adata, False),
    )

    with pytest.raises(DataError, match="No highly variable genes found"):
        integrate_multiple_samples(adata, method="not_real", batch_key="batch", n_pcs=3)


def test_integrate_multiple_samples_sparse_zero_variance_and_scale_fallback(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    X = np.asarray(adata.X, dtype=np.float32)
    X[:, 0] = 0.0
    adata.X = sp.csr_matrix(X / 10.0)
    adata.var["highly_variable"] = True

    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "scanpy.pp.highly_variable_genes",
        lambda adata, **_kwargs: _set_hvg_value(adata, True),
    )

    scale_calls: list[bool] = []

    def _scale(_adata, zero_center=True, max_value=10):
        del max_value
        scale_calls.append(bool(zero_center))
        if zero_center:
            raise RuntimeError("zero-center fail")

    monkeypatch.setattr("scanpy.pp.scale", _scale)
    monkeypatch.setattr(
        "scanpy.tl.pca",
        lambda adata, n_comps, svd_solver, zero_center=False: _set_obsm(
            adata, "X_pca", np.zeros((adata.n_obs, min(n_comps, 3)), dtype=np.float32)
        ),
    )
    monkeypatch.setattr("scanpy.pp.neighbors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("scanpy.tl.umap", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "chatspatial.tools.integration.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )

    out = integrate_multiple_samples(adata, method="not_real", batch_key="batch", n_pcs=3)
    assert out.n_obs == adata.n_obs
    assert scale_calls == [True, False]


def test_integrate_multiple_samples_raises_when_both_scale_strategies_fail(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "scanpy.pp.scale",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("scale fail")),
    )

    with pytest.raises(ProcessingError, match="Data scaling failed completely"):
        integrate_multiple_samples(adata, method="not_real", batch_key="batch", n_pcs=3)


def test_integrate_multiple_samples_pca_nan_inf_and_solver_failures(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "scanpy.pp.highly_variable_genes",
        lambda adata, **_kwargs: _set_hvg_value(adata, True),
    )

    def _scale_nan(adata_obj, **_kwargs):
        arr = np.asarray(adata_obj.X, dtype=np.float32).copy()
        arr[0, 0] = np.nan
        _set_x(adata_obj, arr)

    monkeypatch.setattr("scanpy.pp.scale", _scale_nan)
    adata_nan = adata.copy()
    with pytest.raises(DataError, match="NaN values after scaling"):
        integrate_multiple_samples(adata_nan, method="not_real", batch_key="batch", n_pcs=3)

    def _scale_inf(adata_obj, **_kwargs):
        arr = np.asarray(adata_obj.X, dtype=np.float32).copy()
        arr[0, 0] = np.inf
        _set_x(adata_obj, arr)

    monkeypatch.setattr("scanpy.pp.scale", _scale_inf)
    adata_inf = adata.copy()
    with pytest.raises(DataError, match="infinite values after scaling"):
        integrate_multiple_samples(adata_inf, method="not_real", batch_key="batch", n_pcs=3)

    monkeypatch.setattr("scanpy.pp.scale", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "scanpy.tl.pca",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("pca fail")),
    )
    with pytest.raises(ProcessingError, match="PCA failed"):
        integrate_multiple_samples(adata, method="not_real", batch_key="batch", n_pcs=3)


def test_integrate_multiple_samples_harmony_keeps_new_shape_without_transpose(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)

    fake_harmonypy = ModuleType("harmonypy")
    fake_harmonypy.run_harmony = lambda **_kwargs: SimpleNamespace(
        Z_corr=np.full((adata.n_obs, 4), 2.0, dtype=np.float32)
    )
    monkeypatch.setitem(__import__("sys").modules, "harmonypy", fake_harmonypy)

    out = integrate_multiple_samples(adata, method="harmony", batch_key="batch", n_pcs=4)
    assert out.obsm["X_pca_harmony"].shape == (adata.n_obs, 4)
    assert captured["neighbors"][0]["use_rep"] == "X_pca_harmony"


def test_integrate_multiple_samples_scanorama_raw_fallback_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)

    import scanpy.external as sce

    monkeypatch.setattr(sce, "pp", SimpleNamespace(), raising=False)

    fake_scanorama = ModuleType("scanorama")
    fake_scanorama.integrate = lambda datasets, genes_list, dimred=100: (
        [np.zeros((ds.shape[0], 3), dtype=np.float32) for ds in datasets],
        genes_list[0],
    )
    monkeypatch.setitem(__import__("sys").modules, "scanorama", fake_scanorama)

    out = integrate_multiple_samples(adata, method="scanorama", batch_key="batch", n_pcs=3)
    assert "X_scanorama" in out.obsm
    assert captured["neighbors"][0]["use_rep"] == "X_scanorama"


def test_reassemble_scanorama_embeddings_restores_interleaved_rows():
    out = _reassemble_scanorama_embeddings(
        [np.array([0, 2]), np.array([1, 3])],
        [
            np.array([[10.0, 11.0], [20.0, 21.0]]),
            np.array([[30.0, 31.0], [40.0, 41.0]]),
        ],
        n_obs=4,
    )

    np.testing.assert_array_equal(
        out,
        np.array(
            [[10.0, 11.0], [30.0, 31.0], [20.0, 21.0], [40.0, 41.0]],
            dtype=np.float32,
        ),
    )


@pytest.mark.parametrize(
    ("batch_indices", "embeddings", "message"),
    [
        ([np.array([0]), np.array([1])], [np.ones((1, 2))], "number of embeddings"),
        ([np.array([0, 1])], [np.ones((1, 2))], "rows for 2 observations"),
        (
            [np.array([0]), np.array([1])],
            [np.ones((1, 2)), np.ones((1, 3))],
            "inconsistent embedding dimensions",
        ),
        (
            [np.array([0]), np.array([0])],
            [np.ones((1, 2)), np.ones((1, 2))],
            "duplicated or outside",
        ),
        ([np.array([0])], [np.ones((1, 2))], "missing embeddings for 1"),
    ],
)
def test_reassemble_scanorama_embeddings_rejects_inconsistent_backend_output(
    batch_indices, embeddings, message
):
    with pytest.raises(ProcessingError, match=message):
        _reassemble_scanorama_embeddings(batch_indices, embeddings, n_obs=2)


def test_rescale_spatial_coordinates_error_branches(minimal_spatial_adata):
    adata_missing = minimal_spatial_adata.copy()
    adata_missing.obs["batch"] = "only"
    del adata_missing.obsm["spatial"]
    with pytest.raises(DataNotFoundError, match="spatial coordinates"):
        rescale_spatial_coordinates(adata_missing, batch_key="batch")

    adata_empty = minimal_spatial_adata.copy()[:0].copy()
    adata_empty.obsm["spatial"] = np.zeros((0, 2), dtype=float)
    adata_empty.obs["batch"] = []
    with pytest.raises(DataError, match="Dataset is empty"):
        rescale_spatial_coordinates(adata_empty, batch_key="batch")

    adata_missing_batch = minimal_spatial_adata.copy()
    with pytest.raises(ParameterError, match="Batch key 'batch' not found"):
        rescale_spatial_coordinates(adata_missing_batch, batch_key="batch")

    adata_missing_label = minimal_spatial_adata.copy()
    adata_missing_label.obs["batch"] = np.array(
        ["A"] * (adata_missing_label.n_obs - 1) + [None], dtype=object
    )
    with pytest.raises(DataError, match="contains missing values"):
        rescale_spatial_coordinates(adata_missing_label, batch_key="batch")

    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = "only"
    with pytest.raises(ParameterError, match="Reference batch 'missing' not found"):
        rescale_spatial_coordinates(adata, batch_key="batch", reference_batch="missing")


def test_rescale_spatial_coordinates_preserves_coordinate_dimensions(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) % 2 == 0, 1, 2)
    third_dimension = np.linspace(0.0, 1.0, adata.n_obs)[:, None]
    adata.obsm["spatial"] = np.column_stack(
        [adata.obsm["spatial"], third_dimension]
    )

    out = rescale_spatial_coordinates(
        adata,
        batch_key="batch",
        reference_batch="1",
    )

    assert out.obsm["spatial_aligned"].shape == (adata.n_obs, 3)
    assert out.uns["spatial_alignment_metadata"]["parameters"]["reference_batch"] == "1"


def test_integrate_with_scvi_auto_epochs_for_medium_and_large_datasets(
    monkeypatch: pytest.MonkeyPatch,
):
    calls_medium: dict[str, object] = {}
    calls_large: dict[str, object] = {}

    _install_fake_scvi(monkeypatch, calls_medium)
    monkeypatch.setattr("scanpy.pp.neighbors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("chatspatial.tools.integration.get_accelerator", lambda **_kwargs: "cpu")

    adata_medium = AnnData(np.full((1001, 2), 1.0, dtype=np.float32))
    adata_medium.obs["batch"] = np.where(np.arange(1001) < 500, "a", "b")
    integrate_with_scvi(adata_medium, batch_key="batch", n_epochs=None, use_gpu=False)
    assert calls_medium["train"]["max_epochs"] == 200

    _install_fake_scvi(monkeypatch, calls_large)
    adata_large = AnnData(np.full((10001, 2), 1.0, dtype=np.float32))
    adata_large.obs["batch"] = np.where(np.arange(10001) < 5000, "a", "b")
    integrate_with_scvi(adata_large, batch_key="batch", n_epochs=None, use_gpu=False)
    assert calls_large["train"]["max_epochs"] == 100


def test_integrate_multiple_samples_scvi_wraps_integration_errors(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.integrate_with_scvi",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
    )

    with pytest.raises(ProcessingError, match="scVI integration failed: boom"):
        integrate_multiple_samples(adata, method="scvi", batch_key="batch")


def test_integrate_multiple_samples_raises_when_too_few_pca_components_possible(
    monkeypatch: pytest.MonkeyPatch,
):
    X = np.vstack(
        [
            np.zeros((1, 60), dtype=np.float32),
            np.ones((1, 60), dtype=np.float32),
        ]
    )
    adata = AnnData(X)
    adata.obs["batch"] = ["b1", "b2"]
    adata.var["highly_variable"] = True

    monkeypatch.setattr(
        "chatspatial.tools.integration.check_is_integer_counts",
        lambda *_a, **_k: (False, False, False),
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr("scanpy.pp.scale", lambda *_args, **_kwargs: None)

    with pytest.raises(DataError, match="Cannot perform PCA: only 1 components possible"):
        integrate_multiple_samples(adata, method="not_real", batch_key="batch", n_pcs=30)


def test_integrate_multiple_samples_sparse_nan_values_raise_data_error(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = sp.csr_matrix(np.asarray(adata.X, dtype=np.float32) / 10.0)

    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "scanpy.pp.highly_variable_genes",
        lambda adata, **_kwargs: _set_hvg_value(adata, True),
    )

    def _scale_with_nan(adata_obj, **_kwargs):
        arr = np.asarray(adata_obj.X.toarray(), dtype=np.float32)
        arr[0, 0] = np.nan
        _set_x(adata_obj, sp.csr_matrix(arr))

    monkeypatch.setattr("scanpy.pp.scale", _scale_with_nan)

    with pytest.raises(DataError, match="NaN values after scaling"):
        integrate_multiple_samples(adata, method="not_real", batch_key="batch", n_pcs=3)


def test_integrate_multiple_samples_sparse_inf_values_raise_data_error(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = sp.csr_matrix(np.asarray(adata.X, dtype=np.float32) / 10.0)

    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "scanpy.pp.highly_variable_genes",
        lambda adata, **_kwargs: _set_hvg_value(adata, True),
    )

    def _scale_with_inf(adata_obj, **_kwargs):
        arr = np.asarray(adata_obj.X.toarray(), dtype=np.float32)
        arr[0, 0] = np.inf
        _set_x(adata_obj, sp.csr_matrix(arr))

    monkeypatch.setattr("scanpy.pp.scale", _scale_with_inf)

    with pytest.raises(DataError, match="infinite values after scaling"):
        integrate_multiple_samples(adata, method="not_real", batch_key="batch", n_pcs=3)


def test_integrate_multiple_samples_scanorama_wraps_exceptions(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)
    monkeypatch.setattr(
        "scanpy.external.pp.scanorama_integrate",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("sboom")),
    )

    with pytest.raises(ProcessingError, match="Scanorama integration failed: sboom"):
        integrate_multiple_samples(adata, method="scanorama", batch_key="batch", n_pcs=3)


def test_integrate_multiple_samples_harmony_metadata_falls_back_to_x_harmony_key(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "b1", "b2")
    adata.var["highly_variable"] = True
    adata.X = np.clip(adata.X.astype(np.float32) / 10.0, 0, 10)

    captured: dict[str, object] = {}
    _install_classical_integration_mocks(monkeypatch, captured)

    fake_harmonypy = ModuleType("harmonypy")
    fake_harmonypy.run_harmony = lambda **_kwargs: SimpleNamespace(
        Z_corr=np.full((adata.n_obs, 3), 1.0, dtype=np.float32)
    )
    monkeypatch.setitem(__import__("sys").modules, "harmonypy", fake_harmonypy)

    def _neighbors(adata_obj, **kwargs):
        captured.setdefault("neighbors", []).append(kwargs)
        adata_obj.obsm["X_harmony"] = adata_obj.obsm["X_pca_harmony"]
        del adata_obj.obsm["X_pca_harmony"]

    monkeypatch.setattr("scanpy.pp.neighbors", _neighbors)

    out = integrate_multiple_samples(adata, method="harmony", batch_key="batch", n_pcs=3)

    assert out.n_obs == adata.n_obs
    assert captured["results_keys"] == {"obsm": ["X_harmony"]}


def test_integrate_multiple_samples_deletes_incomplete_diffmap_artifacts_with_concat(
    monkeypatch: pytest.MonkeyPatch,
):
    adata_1 = AnnData(np.vstack([np.zeros((5, 60), dtype=np.float32), np.ones((5, 60), dtype=np.float32)]))
    adata_2 = AnnData(np.vstack([np.ones((5, 60), dtype=np.float32), np.zeros((5, 60), dtype=np.float32)]))
    adata_1.obs["batch"] = "b1"
    adata_2.obs["batch"] = "b2"

    combined = AnnData(np.vstack([np.zeros((10, 60), dtype=np.float32), np.ones((10, 60), dtype=np.float32)]))
    combined.obs["batch"] = ["b1"] * 10 + ["b2"] * 10
    combined.var["highly_variable"] = True
    _set_obsm(combined, "X_diffmap", np.zeros((20, 2), dtype=np.float32))
    combined.uns["diffmap_evals"] = np.array([1.0, 0.5], dtype=np.float32)

    monkeypatch.setattr(
        "chatspatial.tools.integration.check_is_integer_counts",
        lambda *_a, **_k: (False, False, False),
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.ad.concat",
        lambda *_args, **_kwargs: combined,
    )
    monkeypatch.setattr(
        "chatspatial.tools.integration.validate_adata_basics",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr("scanpy.pp.scale", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "scanpy.tl.pca",
        lambda adata, n_comps, svd_solver, zero_center=False: _set_obsm(
            adata, "X_pca", np.zeros((adata.n_obs, min(n_comps, 3)), dtype=np.float32)
        ),
    )
    monkeypatch.setattr("scanpy.pp.neighbors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr("scanpy.tl.umap", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        "chatspatial.tools.integration.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )

    out = integrate_multiple_samples([adata_1, adata_2], method="not_real", batch_key="batch")

    assert "X_diffmap" not in out.obsm
    assert "diffmap_evals" not in out.uns


def test_rescale_spatial_coordinates_uses_first_batch_when_reference_is_none(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = np.where(np.arange(adata.n_obs) < adata.n_obs // 2, "first", "second")

    captured: dict[str, object] = {}
    monkeypatch.setattr(
        "chatspatial.tools.integration.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )

    out = rescale_spatial_coordinates(adata, batch_key="batch", reference_batch=None)

    assert out.obsm["spatial_aligned"].shape == (adata.n_obs, 2)
    assert captured["parameters"]["reference_batch"] == "first"
