"""Unit tests for deconvolution base utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from chatspatial.tools.deconvolution.base import (
    MethodConfig,
    _prepare_counts,
    check_model_convergence,
    create_deconvolution_stats,
    prepare_deconvolution,
)
from chatspatial.utils.exceptions import DataError


class DummyCtx:
    async def warning(self, msg: str):
        return None


def test_method_config_extract_kwargs_maps_and_adds_gpu_flag():
    class Params:
        flashdeconv_sketch_dim = 256
        flashdeconv_lambda_spatial = 1234.0
        use_gpu = True

    cfg = MethodConfig(
        module_name="flashdeconv",
        dependencies=("flashdeconv",),
        supports_gpu=True,
        param_mapping=(
            ("flashdeconv_sketch_dim", "sketch_dim"),
            ("flashdeconv_lambda_spatial", "lambda_spatial"),
        ),
    )

    kwargs = cfg.extract_kwargs(Params())
    assert kwargs["sketch_dim"] == 256
    assert kwargs["lambda_spatial"] == 1234.0
    assert kwargs["use_gpu"] is True


def test_create_deconvolution_stats_has_consistent_summary_fields():
    proportions = pd.DataFrame(
        {
            "T": [0.9, 0.2, 0.1],
            "B": [0.1, 0.8, 0.9],
        },
        index=["spot_1", "spot_2", "spot_3"],
    )
    stats = create_deconvolution_stats(
        proportions=proportions,
        common_genes=["g1", "g2", "g3"],
        method="flashdeconv",
        device="CPU",
    )

    assert stats["method"] == "flashdeconv"
    assert stats["n_spots"] == 3
    assert stats["n_cell_types"] == 2
    assert stats["genes_used"] == 3
    assert sum(stats["dominant_types"].values()) == 3


@pytest.mark.asyncio
async def test_prepare_counts_uses_counts_layer_over_raw(minimal_spatial_adata):
    """_prepare_counts should prefer layers['counts'] over .raw."""
    adata = minimal_spatial_adata.copy()
    counts_val = np.full((adata.n_obs, adata.n_vars), 42.0)
    adata.layers["counts"] = counts_val
    raw = adata.copy()
    raw.X = np.full_like(counts_val, 99.0)
    adata.raw = raw
    adata.obsm["extra"] = np.ones((adata.n_obs, 2))

    out = await _prepare_counts(adata, "Spatial", DummyCtx(), require_int_dtype=False)
    assert np.allclose(out.X, 42.0), "Should use layers['counts'], not .raw"
    assert "extra" in out.obsm


@pytest.mark.asyncio
async def test_prepare_counts_falls_back_to_x_when_no_counts_layer(
    minimal_spatial_adata,
):
    """_prepare_counts should fall back to adata.X when no counts layer."""
    adata = minimal_spatial_adata.copy()
    # Set .raw to something different to verify it is NOT used
    raw = adata.copy()
    raw.X = np.full((adata.n_obs, adata.n_vars), 99.0)
    adata.raw = raw
    if "counts" in adata.layers:
        del adata.layers["counts"]

    original_x = np.asarray(adata.X).copy()
    out = await _prepare_counts(adata, "Spatial", DummyCtx(), require_int_dtype=False)
    assert np.allclose(np.asarray(out.X), original_x), "Should use adata.X, not .raw"


@pytest.mark.asyncio
async def test_prepare_counts_converts_integer_like_data_to_int32(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.X = np.rint(np.asarray(adata.X)).astype(np.float64)

    out = await _prepare_counts(adata, "Spatial", DummyCtx(), require_int_dtype=True)
    assert out.X.dtype == np.int32


@pytest.mark.asyncio
async def test_prepare_counts_prefers_counts_layer_when_raw_absent(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    if adata.raw is not None:
        adata.raw = None
    counts = np.full((adata.n_obs, adata.n_vars), 7.0, dtype=np.float64)
    adata.layers["counts"] = counts
    adata.X = np.zeros_like(counts)

    out = await _prepare_counts(adata, "Spatial", DummyCtx(), require_int_dtype=False)
    assert np.allclose(out.X, counts)


@pytest.mark.asyncio
async def test_prepare_counts_does_not_force_int_for_decimal_values(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.X = np.asarray(adata.X, dtype=np.float64) + 0.123

    out = await _prepare_counts(adata, "Spatial", DummyCtx(), require_int_dtype=True)
    assert out.X.dtype == np.float64


@pytest.mark.asyncio
async def test_prepare_counts_handles_empty_sparse_matrix_sampling(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.X = sparse.csr_matrix((adata.n_obs, adata.n_vars), dtype=np.float64)

    out = await _prepare_counts(adata, "Spatial", DummyCtx(), require_int_dtype=True)
    assert out.X.dtype == np.int32
    assert out.X.shape == (adata.n_obs, adata.n_vars)


@pytest.mark.asyncio
async def test_prepare_counts_integer_sparse_input_stays_sparse(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    dense = np.rint(np.asarray(adata.X)).astype(np.float64)
    adata.X = sparse.csr_matrix(dense)

    out = await _prepare_counts(adata, "Spatial", DummyCtx(), require_int_dtype=True)
    assert sparse.issparse(out.X)
    assert out.X.dtype == np.int32


@pytest.mark.asyncio
async def test_prepare_deconvolution_with_preprocess_hook_and_no_spatial_key(
    minimal_spatial_adata,
):
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()
    reference.obs["cell_type"] = ["A"] * (reference.n_obs // 2) + ["B"] * (
        reference.n_obs - reference.n_obs // 2
    )
    del spatial.obsm["spatial"]

    calls = {"preprocess": 0}

    async def _preprocess(sp, ref, _ctx):
        calls["preprocess"] += 1
        sp = sp.copy()
        ref = ref.copy()
        sp.uns["prep"] = "ok"
        ref.uns["prep"] = "ok"
        return sp, ref

    data = await prepare_deconvolution(
        spatial_adata=spatial,
        reference_adata=reference,
        cell_type_key="cell_type",
        ctx=DummyCtx(),
        min_common_genes=5,
        preprocess=_preprocess,
    )

    assert calls["preprocess"] == 1
    assert data.spatial_coords is None
    assert data.n_spots == spatial.n_obs
    assert data.n_cell_types == 2
    assert data.n_genes == spatial.n_vars
    assert data.spatial.uns["prep"] == "ok"


@pytest.mark.asyncio
async def test_prepare_deconvolution_repairs_names_only_in_workspaces(
    minimal_spatial_adata,
):
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()
    reference.obs["cell_type"] = ["A"] * (reference.n_obs // 2) + ["B"] * (
        reference.n_obs - reference.n_obs // 2
    )
    duplicate_names = spatial.var_names.tolist()
    duplicate_names[1] = duplicate_names[0]
    spatial.var_names = duplicate_names
    reference.var_names = duplicate_names

    with pytest.warns(UserWarning, match="Variable names are not unique"):
        data = await prepare_deconvolution(
            spatial_adata=spatial,
            reference_adata=reference,
            cell_type_key="cell_type",
            ctx=DummyCtx(),
            min_common_genes=5,
        )

    assert not spatial.var_names.is_unique
    assert not reference.var_names.is_unique
    assert data.spatial.var_names.is_unique
    assert data.reference.var_names.is_unique


@pytest.mark.asyncio
async def test_prepare_deconvolution_rejects_single_cell_type(minimal_spatial_adata):
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()
    reference.obs["cell_type"] = ["A"] * reference.n_obs

    with pytest.raises(DataError, match="at least 2 cell types"):
        await prepare_deconvolution(
            spatial_adata=spatial,
            reference_adata=reference,
            cell_type_key="cell_type",
            ctx=DummyCtx(),
            min_common_genes=5,
        )


def test_method_config_requires_reference_is_true():
    cfg = MethodConfig(module_name="x", dependencies=("dep",))
    assert cfg.requires_reference is True


def test_check_model_convergence_returns_true_without_history():
    model = object()
    converged, message = check_model_convergence(model, "Dummy")
    assert converged is True
    assert message is None


def test_check_model_convergence_returns_true_for_short_history():
    model = type("Model", (), {"history": {"elbo_train": [1.0, 0.9, 0.8]}})()
    converged, message = check_model_convergence(model, "Dummy", convergence_window=10)
    assert converged is True
    assert message is None


def test_check_model_convergence_detects_unstable_training():
    history = [1.0, 3.0, 1.2, 3.1, 1.1, 3.2, 1.0, 3.3, 1.2, 3.0]
    model = type("Model", (), {"history": {"elbo_validation": history}})()
    converged, message = check_model_convergence(
        model,
        "DemoModel",
        convergence_threshold=0.05,
        convergence_window=10,
    )
    assert converged is False
    assert "DemoModel may not have fully converged" in (message or "")


def test_check_model_convergence_handles_zero_mean_without_warning():
    history = [0.0] * 10
    model = type("Model", (), {"history": {"train_loss_epoch": history}})()
    converged, message = check_model_convergence(
        model,
        "ZeroModel",
        convergence_window=10,
    )
    assert converged is True
    assert message is None


@pytest.mark.asyncio
async def test_prepare_counts_recovers_expression_from_raw_when_x_is_scaled(
    minimal_spatial_adata,
):
    """Scaled data cannot be mixed; every backend would return all-zero proportions.

    A scanpy-tutorial reference keeps z-scored values in X and the expression
    it was scaled from in .raw, so .raw is what deconvolution must use.
    """
    adata = minimal_spatial_adata.copy()
    raw = adata.copy()
    raw.X = np.full((adata.n_obs, adata.n_vars), 7.0)
    adata.raw = raw
    adata.X = np.full((adata.n_obs, adata.n_vars), -1.5)
    adata.layers.pop("counts", None)

    warnings: list[str] = []

    class RecordingCtx:
        async def warning(self, msg: str):
            warnings.append(msg)

    out = await _prepare_counts(
        adata, "Reference", RecordingCtx(), require_int_dtype=False
    )

    assert np.allclose(out.X, 7.0)
    assert out.n_obs == adata.n_obs
    assert any("scaled" in message for message in warnings)


@pytest.mark.asyncio
async def test_prepare_counts_rejects_scaled_data_without_recoverable_source(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.raw = None
    adata.layers.pop("counts", None)
    adata.X = np.full((adata.n_obs, adata.n_vars), -0.5)

    with pytest.raises(DataError, match="negative values"):
        await _prepare_counts(adata, "Reference", DummyCtx(), require_int_dtype=False)


@pytest.mark.unit
def test_missing_counts_name_the_dataset_that_lacks_them():
    """Deconvolution holds two datasets, so an unlabelled error is a guess.

    RCTD on a scanpy-processed reference failed with "No raw integer counts
    found" and no indication of whether the spatial data or the reference was
    the one to fix.
    """
    import anndata as ad

    from chatspatial.tools.deconvolution import base
    from chatspatial.utils.exceptions import DataError

    normalized = ad.AnnData(np.array([[0.5, -1.25], [2.5, 0.75]], dtype=np.float32))

    with pytest.raises(DataError) as raised:
        base._counts_for(
            normalized, "Reference", prefer_complete_genes=False, require_counts=True
        )

    assert str(raised.value).startswith("Reference:")
    # The original diagnosis is kept, not replaced.
    assert "integer counts" in str(raised.value)


@pytest.mark.unit
def test_available_counts_are_returned_unlabelled():
    import anndata as ad

    from chatspatial.tools.deconvolution import base

    counts = ad.AnnData(np.array([[3.0, 1.0], [2.0, 4.0]], dtype=np.float32))

    result = base._counts_for(
        counts, "Spatial", prefer_complete_genes=False, require_counts=True
    )

    assert result.is_integer_counts
    assert list(result.var_names) == list(counts.var_names)
