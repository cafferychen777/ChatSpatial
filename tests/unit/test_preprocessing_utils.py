"""Unit tests for preprocessing helper and preprocess_data contracts."""

from __future__ import annotations

from contextlib import contextmanager
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData

from chatspatial.models.data import PreprocessingParameters
from chatspatial.tools import preprocessing as preprocessing_mod
from chatspatial.tools.preprocessing import _compute_safe_percent_top, preprocess_data
from chatspatial.utils.exceptions import DataError, DependencyError, ParameterError, ProcessingError


class DummyCtx:
    def __init__(self, adata: AnnData):
        self._adata = adata
        self.saved_adata: AnnData | None = None
        self.warnings: list[str] = []
        self.infos: list[str] = []
        self.config_logs: list[tuple[str, dict]] = []

    async def get_adata(self, _data_id: str) -> AnnData:
        return self._adata

    async def set_adata(self, _data_id: str, adata: AnnData) -> None:
        self.saved_adata = adata

    async def warning(self, msg: str) -> None:
        self.warnings.append(msg)

    async def info(self, msg: str) -> None:
        self.infos.append(msg)

    def log_config(self, title: str, config: dict) -> None:
        self.config_logs.append((title, config))


def _make_adata(n_obs: int = 30, n_vars: int = 120) -> AnnData:
    rng = np.random.default_rng(7)
    X = rng.poisson(5, size=(n_obs, n_vars)).astype(np.float32)
    adata = AnnData(X=X)
    # Include mito/ribo-like names to exercise annotation columns.
    var_names = [f"gene_{i}" for i in range(n_vars)]
    if n_vars > 2:
        var_names[0] = "MT-ND1"
        var_names[1] = "RPS3"
    adata.var_names = var_names
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    return adata


def _install_lightweight_preprocess_mocks(monkeypatch: pytest.MonkeyPatch) -> None:
    async def _noop_ensure_unique(_adata, _ctx, _label):
        return None

    def _identity_standardize(adata, copy=False):
        return adata.copy() if copy else adata

    def _sample_values(adata):
        arr = adata.X
        flat = np.asarray(arr).reshape(-1)
        return flat[: min(len(flat), 128)]

    def _calc_qc_metrics(adata, qc_vars=None, percent_top=None, inplace=True):
        del qc_vars, percent_top, inplace
        counts = np.asarray(adata.X)
        adata.obs["n_genes_by_counts"] = (counts > 0).sum(axis=1)
        adata.obs["total_counts"] = counts.sum(axis=1)
        if "mt" in adata.var.columns:
            mt_mask = adata.var["mt"].to_numpy()
            mt_counts = counts[:, mt_mask].sum(axis=1)
            total = np.clip(counts.sum(axis=1), 1e-9, None)
            adata.obs["pct_counts_mt"] = mt_counts / total * 100.0

    def _hvg(adata, n_top_genes=2000):
        flags = np.zeros(adata.n_vars, dtype=bool)
        flags[: min(max(int(n_top_genes), 0), adata.n_vars)] = True
        adata.var["highly_variable"] = flags

    def _no_op(*args, **kwargs):
        del args, kwargs
        return None

    monkeypatch.setattr(preprocessing_mod, "ensure_unique_var_names_async", _noop_ensure_unique)
    monkeypatch.setattr(preprocessing_mod, "standardize_adata", _identity_standardize)
    monkeypatch.setattr(preprocessing_mod, "sample_expression_values", _sample_values)
    monkeypatch.setattr(preprocessing_mod.sc.pp, "calculate_qc_metrics", _calc_qc_metrics)
    monkeypatch.setattr(preprocessing_mod.sc.pp, "filter_genes", _no_op)
    monkeypatch.setattr(preprocessing_mod.sc.pp, "filter_cells", _no_op)
    monkeypatch.setattr(preprocessing_mod.sc.pp, "subsample", _no_op)
    monkeypatch.setattr(preprocessing_mod.sc.pp, "normalize_total", _no_op)
    monkeypatch.setattr(preprocessing_mod.sc.pp, "log1p", _no_op)
    monkeypatch.setattr(preprocessing_mod.sc.pp, "highly_variable_genes", _hvg)
    monkeypatch.setattr(preprocessing_mod.sc.pp, "scale", _no_op)


def _install_fake_rpy2_for_sct(
    monkeypatch: pytest.MonkeyPatch,
    *,
    pearson_residuals: np.ndarray,
    residual_variance: np.ndarray,
    kept_genes: list[str],
) -> None:
    class _Converter:
        def __add__(self, _other):
            return self

    @contextmanager
    def _localconverter(_converter):
        yield

    def _ro_r(code: str):
        text = code.strip()
        if text == "pearson_residuals":
            return pearson_residuals
        if text == "residual_variance":
            return residual_variance
        if text == "kept_genes":
            return kept_genes
        return None

    ro_mod = ModuleType("rpy2.robjects")
    ro_mod.default_converter = _Converter()
    ro_mod.globalenv = {}
    ro_mod.NULL = object()
    ro_mod.StrVector = lambda vals: list(vals)
    ro_mod.r = _ro_r

    conversion_mod = ModuleType("rpy2.robjects.conversion")
    conversion_mod.localconverter = _localconverter

    numpy2ri_mod = ModuleType("rpy2.robjects.numpy2ri")
    numpy2ri_mod.converter = _Converter()

    rpy2_mod = ModuleType("rpy2")
    monkeypatch.setitem(sys.modules, "rpy2", rpy2_mod)
    monkeypatch.setitem(sys.modules, "rpy2.robjects", ro_mod)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.numpy2ri", numpy2ri_mod)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.conversion", conversion_mod)


@pytest.mark.asyncio
async def test_preprocess_data_success_persists_core_artifacts(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    adata = _make_adata(n_obs=24, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(
        normalization="log",
        n_hvgs=20,
        subsample_genes=12,
        filter_genes_min_cells=1,
        filter_cells_min_genes=1,
        remove_mito_genes=False,
        remove_ribo_genes=False,
        filter_mito_pct=None,
        scale=False,
    )

    result = await preprocess_data("d1", ctx, params)

    assert result.data_id == "d1"
    assert result.n_cells == 24
    assert result.n_genes == 12
    assert result.n_hvgs == 12

    assert ctx.saved_adata is not None
    assert "counts" in ctx.saved_adata.layers
    assert ctx.saved_adata.layers["counts"].shape == ctx.saved_adata.X.shape
    assert ctx.saved_adata.raw is not None
    assert ctx.saved_adata.uns["preprocessing"]["completed"] is True


@pytest.mark.asyncio
async def test_preprocess_data_reuses_raw_matrix_for_counts_layer_when_aligned(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    adata = _make_adata(n_obs=24, n_vars=120)
    adata.raw = adata.copy()
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(
        normalization="log",
        n_hvgs=20,
        filter_genes_min_cells=1,
        filter_cells_min_genes=1,
        remove_mito_genes=False,
        remove_ribo_genes=False,
        filter_mito_pct=None,
        scale=False,
    )

    await preprocess_data("d1_raw_reuse", ctx, params)

    assert ctx.saved_adata is not None
    counts = ctx.saved_adata.layers["counts"]
    raw_x = ctx.saved_adata.raw.X
    assert counts.shape == raw_x.shape
    if sp.issparse(counts):
        assert counts is raw_x
    else:
        assert np.shares_memory(np.asarray(counts), np.asarray(raw_x))


@pytest.mark.asyncio
async def test_preprocess_data_warns_when_hvgs_too_low(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    adata = _make_adata(n_obs=20, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(
        normalization="log",
        n_hvgs=30,
        filter_genes_min_cells=1,
        filter_cells_min_genes=1,
        remove_mito_genes=False,
        filter_mito_pct=None,
    )

    await preprocess_data("d2", ctx, params)

    assert any("recommended minimum of 500" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_preprocess_data_rejects_none_normalization_for_raw_counts(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    # _make_adata produces integer counts (Poisson draws) — check_is_integer_counts
    # detects these regardless of magnitude, catching low-depth platforms too.
    adata = _make_adata(n_obs=10, n_vars=120)
    original_shape = adata.shape
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(normalization="none", filter_mito_pct=None)

    with pytest.raises(DataError, match="Cannot perform HVG selection on raw counts"):
        await preprocess_data("d3", ctx, params)

    assert adata.shape == original_shape
    assert ctx.saved_adata is None


@pytest.mark.asyncio
async def test_preprocess_data_rejects_none_normalization_low_depth_counts(monkeypatch: pytest.MonkeyPatch):
    """Low-depth integer data (max < 100) is still detected as raw counts."""
    _install_lightweight_preprocess_mocks(monkeypatch)

    rng = np.random.default_rng(42)
    # Simulate targeted panel (MERFISH-like): integer counts, max ~ 10
    X = rng.poisson(2, size=(20, 50)).astype(np.float32)
    adata = AnnData(X=X)
    adata.var_names = [f"gene_{i}" for i in range(50)]
    adata.obs_names = [f"cell_{i}" for i in range(20)]

    ctx = DummyCtx(adata)
    params = PreprocessingParameters(normalization="none", filter_mito_pct=None)

    with pytest.raises(DataError, match="Cannot perform HVG selection on raw counts"):
        await preprocess_data("d3_low", ctx, params)


@pytest.mark.asyncio
async def test_preprocess_data_failed_validation_does_not_mutate_source(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    adata = _make_adata(n_obs=20, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(
        normalization="none",
        subsample_spots=5,
        filter_mito_pct=None,
    )

    with pytest.raises(DataError, match="Cannot perform HVG selection on raw counts"):
        await preprocess_data("d3_subsample", ctx, params)

    assert adata.n_obs == 20
    assert ctx.saved_adata is None


@pytest.mark.asyncio
async def test_preprocess_data_unknown_normalization_raises(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    adata = _make_adata(n_obs=10, n_vars=120)
    ctx = DummyCtx(adata)
    # bypass Literal validation to test runtime defensive branch
    params = PreprocessingParameters.model_construct(normalization="invalid", filter_mito_pct=None)

    with pytest.raises(ParameterError, match="Unknown normalization method"):
        await preprocess_data("d4", ctx, params)


@pytest.mark.asyncio
async def test_preprocess_data_gene_subsample_requires_nonempty_hvgs(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _all_false_hvg(adata, n_top_genes=2000):
        del n_top_genes
        adata.var["highly_variable"] = False

    monkeypatch.setattr(preprocessing_mod.sc.pp, "highly_variable_genes", _all_false_hvg)

    adata = _make_adata(n_obs=16, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(
        normalization="log",
        subsample_genes=20,
        n_hvgs=20,
        filter_genes_min_cells=1,
        filter_cells_min_genes=1,
        filter_mito_pct=None,
    )

    with pytest.raises(DataError, match="no genes were marked as highly variable"):
        await preprocess_data("d5", ctx, params)


@pytest.mark.asyncio
async def test_preprocess_data_pearson_residuals_requires_scanpy_support(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(
        preprocessing_mod.sc,
        "experimental",
        SimpleNamespace(pp=SimpleNamespace()),
        raising=False,
    )

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(normalization="pearson_residuals", filter_mito_pct=None)

    with pytest.raises(DependencyError, match="Pearson residuals normalization not available"):
        await preprocess_data("d6", ctx, params)


@pytest.mark.asyncio
async def test_preprocess_data_pearson_residuals_rejects_non_integer_input(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    monkeypatch.setattr(
        preprocessing_mod.sc,
        "experimental",
        SimpleNamespace(pp=SimpleNamespace(normalize_pearson_residuals=lambda _adata: None)),
        raising=False,
    )
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.1, 1.5, 2.0]),
    )

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(normalization="pearson_residuals", filter_mito_pct=None)

    with pytest.raises(DataError, match="requires raw count data"):
        await preprocess_data("d7", ctx, params)


@pytest.mark.asyncio
async def test_preprocess_data_scvi_success_writes_latent_and_metadata(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "require", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 3.0]),
    )

    class _FakeSCVI:
        @staticmethod
        def setup_anndata(adata, layer=None, batch_key=None):
            del adata, layer, batch_key
            return None

        def __init__(self, adata, **kwargs):
            self._adata = adata
            self._kwargs = kwargs

        def train(self, **kwargs):
            del kwargs
            return None

        def get_latent_representation(self):
            return np.ones((self._adata.n_obs, 2), dtype=float)

        def get_normalized_expression(self, library_size=1e4):
            del library_size
            return np.full((self._adata.n_obs, self._adata.n_vars), 2.0, dtype=float)

    fake_scvi = SimpleNamespace(model=SimpleNamespace(SCVI=_FakeSCVI))
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    adata = _make_adata(n_obs=14, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(
        normalization="scvi",
        n_hvgs=30,
        filter_genes_min_cells=1,
        filter_cells_min_genes=1,
        filter_mito_pct=None,
        remove_mito_genes=False,
    )

    result = await preprocess_data("d8", ctx, params)

    assert result.n_cells == 14
    assert ctx.saved_adata is not None
    assert "X_scvi" in ctx.saved_adata.obsm
    assert ctx.saved_adata.obsm["X_scvi"].shape == (14, 2)
    assert ctx.saved_adata.uns["scvi"]["training_completed"] is True


@pytest.mark.asyncio
async def test_preprocess_data_scvi_failure_is_wrapped_as_processing_error(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "require", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 3.0]),
    )

    class _FakeSCVI:
        @staticmethod
        def setup_anndata(adata, layer=None, batch_key=None):
            del adata, layer, batch_key
            return None

        def __init__(self, adata, **kwargs):
            del adata, kwargs

        def train(self, **kwargs):
            del kwargs
            raise RuntimeError("boom")

    fake_scvi = SimpleNamespace(model=SimpleNamespace(SCVI=_FakeSCVI))
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)

    adata = _make_adata(n_obs=14, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(normalization="scvi", filter_mito_pct=None)

    with pytest.raises(ProcessingError, match="scVI normalization failed: boom"):
        await preprocess_data("d9", ctx, params)


def test_compute_safe_percent_top_small_gene_set():
    values = _compute_safe_percent_top(10)
    assert values is not None
    assert all(v < 10 for v in values)
    assert values[-1] == 9


def test_compute_safe_percent_top_standard_gene_set():
    values = _compute_safe_percent_top(1000)
    assert values == [50, 100, 200, 500, 999]


def test_compute_safe_percent_top_degenerate_case():
    assert _compute_safe_percent_top(1) is None


@pytest.mark.asyncio
async def test_preprocess_data_sct_missing_dependency_raises_dependency_error(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _raise_import_error(_pkg, _ctx):
        raise ImportError("r package missing")

    monkeypatch.setattr(preprocessing_mod, "validate_r_package", _raise_import_error)

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(normalization="sct", filter_mito_pct=None)

    with pytest.raises(DependencyError, match="SCTransform requires R and the sctransform package"):
        await preprocess_data("d10", ctx, params)


@pytest.mark.asyncio
async def test_preprocess_data_sct_rejects_non_integer_input(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "validate_r_package", lambda *_a, **_k: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.1, 1.2, 3.0]),
    )

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(normalization="sct", filter_mito_pct=None)

    with pytest.raises(DataError, match="SCTransform requires raw count data"):
        await preprocess_data("d11", ctx, params)


@pytest.mark.asyncio
async def test_preprocess_data_sct_import_failure_wrapped_as_processing_error(
    monkeypatch: pytest.MonkeyPatch,
):
    import builtins

    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "validate_r_package", lambda *_a, **_k: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0]),
    )

    real_import = builtins.__import__

    def _wrapped_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name.startswith("rpy2"):
            raise ModuleNotFoundError("forced missing rpy2")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", _wrapped_import)

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    params = PreprocessingParameters(normalization="sct", filter_mito_pct=None)

    with pytest.raises(ProcessingError, match="SCTransform failed"):
        await preprocess_data("d12", ctx, params)


@pytest.mark.asyncio
async def test_preprocess_data_standardize_failure_warns_and_continues(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _raise_standardize(_adata, copy=False):
        del copy
        raise RuntimeError("std boom")

    monkeypatch.setattr(preprocessing_mod, "standardize_adata", _raise_standardize)

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)

    result = await preprocess_data(
        "d13",
        ctx,
        PreprocessingParameters(normalization="log", filter_mito_pct=None),
    )

    assert result.n_cells == 12
    assert any("Data standardization failed" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_preprocess_data_log_normalization_rejects_negative_values(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([1.0, -0.2, 3.0]),
    )

    adata = _make_adata(n_obs=10, n_vars=120)
    ctx = DummyCtx(adata)

    with pytest.raises(DataError, match="requires non-negative data"):
        await preprocess_data(
            "d14",
            ctx,
            PreprocessingParameters(normalization="log", filter_mito_pct=None),
        )


@pytest.mark.asyncio
async def test_preprocess_data_scrublet_warns_when_too_few_cells(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    adata = _make_adata(n_obs=40, n_vars=120)
    ctx = DummyCtx(adata)

    await preprocess_data(
        "d15",
        ctx,
        PreprocessingParameters(
            normalization="log",
            filter_mito_pct=None,
            use_scrublet=True,
        ),
    )

    assert any("Scrublet requires at least 100 cells" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_preprocess_data_rejects_empty_dataset(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)
    adata = AnnData(np.zeros((0, 10), dtype=np.float32))
    ctx = DummyCtx(adata)

    with pytest.raises(DataError, match="Dataset d16 is empty"):
        await preprocess_data("d16", ctx, PreprocessingParameters(normalization="log"))


@pytest.mark.asyncio
async def test_preprocess_data_qc_failure_is_wrapped(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _boom(*_args, **_kwargs):
        raise RuntimeError("qc boom")

    monkeypatch.setattr(preprocessing_mod.sc.pp, "calculate_qc_metrics", _boom)
    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)

    with pytest.raises(ProcessingError, match="QC metrics failed: qc boom"):
        await preprocess_data("d17", ctx, PreprocessingParameters(normalization="log"))


@pytest.mark.asyncio
async def test_preprocess_data_filters_high_mito_cells(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)
    adata = _make_adata(n_obs=20, n_vars=120)
    adata.X[:4, :] = 0.0
    adata.X[:4, 0] = 2000.0
    adata.X[4:, 0] = 1.0
    adata.X[4:, 2:] = 40.0
    ctx = DummyCtx(adata)

    result = await preprocess_data(
        "d18",
        ctx,
        PreprocessingParameters(
            normalization="log",
            filter_mito_pct=20.0,
            remove_mito_genes=False,
            filter_genes_min_cells=1,
            filter_cells_min_genes=1,
        ),
    )

    assert result.n_cells == 16
    assert result.qc_metrics["n_spots_filtered_mito"] == 4


@pytest.mark.asyncio
async def test_preprocess_data_warns_when_mito_filter_has_no_pct_column(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _calc_qc_without_mito_pct(adata, qc_vars=None, percent_top=None, inplace=True):
        del qc_vars, percent_top, inplace
        counts = np.asarray(adata.X)
        adata.obs["n_genes_by_counts"] = (counts > 0).sum(axis=1)
        adata.obs["total_counts"] = counts.sum(axis=1)

    monkeypatch.setattr(
        preprocessing_mod.sc.pp, "calculate_qc_metrics", _calc_qc_without_mito_pct
    )

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)

    await preprocess_data(
        "d19",
        ctx,
        PreprocessingParameters(normalization="log", filter_mito_pct=15.0),
    )

    assert any("Mitochondrial filtering requested" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_preprocess_data_subsample_spots_invokes_scanpy(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)
    called: dict[str, object] = {}

    def _capture_subsample(adata, n_obs=None, random_state=None):
        called["n_obs"] = n_obs
        called["random_state"] = random_state
        adata._inplace_subset_obs(np.arange(int(n_obs)))

    monkeypatch.setattr(preprocessing_mod.sc.pp, "subsample", _capture_subsample)

    adata = _make_adata(n_obs=18, n_vars=120)
    ctx = DummyCtx(adata)
    result = await preprocess_data(
        "d20",
        ctx,
        PreprocessingParameters(
            normalization="log",
            subsample_spots=7,
            subsample_random_seed=123,
            filter_mito_pct=None,
        ),
    )

    assert called["n_obs"] == 7
    assert called["random_state"] == 123
    assert result.n_cells == 7


@pytest.mark.asyncio
async def test_preprocess_data_scrublet_detects_and_filters_doublets(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _fake_scrublet(adata, **_kwargs):
        predicted = np.zeros(adata.n_obs, dtype=bool)
        predicted[:12] = True
        adata.obs["predicted_doublet"] = predicted
        adata.obs["doublet_score"] = np.linspace(0.0, 1.0, adata.n_obs)
        adata.uns["scrublet"] = {"threshold": 0.33}

    monkeypatch.setattr(preprocessing_mod.sc.pp, "scrublet", _fake_scrublet)
    adata = _make_adata(n_obs=120, n_vars=120)
    ctx = DummyCtx(adata)

    result = await preprocess_data(
        "d21",
        ctx,
        PreprocessingParameters(
            normalization="log",
            filter_mito_pct=None,
            use_scrublet=True,
            scrublet_filter_doublets=True,
        ),
    )

    assert result.n_cells == 108
    assert result.qc_metrics["use_scrublet"] is True
    assert result.qc_metrics["n_doublets_detected"] == 12
    assert result.qc_metrics["n_cells_after_doublet_filter"] == 108
    assert any("removed from dataset" in msg for msg in ctx.infos)


@pytest.mark.asyncio
async def test_preprocess_data_scrublet_detects_and_keeps_doublets(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _fake_scrublet(adata, **_kwargs):
        predicted = np.zeros(adata.n_obs, dtype=bool)
        predicted[:5] = True
        adata.obs["predicted_doublet"] = predicted
        adata.obs["doublet_score"] = np.linspace(0.1, 0.9, adata.n_obs)

    monkeypatch.setattr(preprocessing_mod.sc.pp, "scrublet", _fake_scrublet)
    adata = _make_adata(n_obs=120, n_vars=120)
    ctx = DummyCtx(adata)

    result = await preprocess_data(
        "d22",
        ctx,
        PreprocessingParameters(
            normalization="log",
            filter_mito_pct=None,
            use_scrublet=True,
            scrublet_filter_doublets=False,
            scrublet_threshold=0.42,
        ),
    )

    assert result.n_cells == 120
    assert result.qc_metrics["n_doublets_detected"] == 5
    assert result.qc_metrics["scrublet_threshold"] == pytest.approx(0.42)
    assert any("kept in dataset" in msg for msg in ctx.infos)


@pytest.mark.asyncio
async def test_preprocess_data_scrublet_failure_records_warning(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(
        preprocessing_mod.sc.pp,
        "scrublet",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("scrub boom")),
    )

    adata = _make_adata(n_obs=120, n_vars=120)
    ctx = DummyCtx(adata)
    result = await preprocess_data(
        "d23",
        ctx,
        PreprocessingParameters(normalization="log", filter_mito_pct=None, use_scrublet=True),
    )

    assert result.qc_metrics["use_scrublet"] is False
    assert "scrub boom" in result.qc_metrics["scrublet_error"]
    assert any("Scrublet doublet detection failed" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_preprocess_data_log_normalization_uses_explicit_target_sum(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    captured: dict[str, float] = {}

    def _capture_normalize_total(_adata, target_sum=None):
        captured["target_sum"] = float(target_sum)

    monkeypatch.setattr(preprocessing_mod.sc.pp, "normalize_total", _capture_normalize_total)
    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)

    await preprocess_data(
        "d24",
        ctx,
        PreprocessingParameters(
            normalization="log",
            normalize_target_sum=1234.0,
            filter_mito_pct=None,
        ),
    )

    assert captured["target_sum"] == pytest.approx(1234.0)


@pytest.mark.asyncio
async def test_preprocess_data_sct_success_subsets_genes_and_stores_metadata(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "validate_r_package", lambda *_a, **_k: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0, 3.0]),
    )

    adata = _make_adata(n_obs=16, n_vars=120)
    kept_genes = list(adata.var_names[:110])
    pearson_residuals = np.full((110, adata.n_obs), 0.25, dtype=np.float32)
    residual_variance = np.linspace(0.0, 2.0, 110, dtype=np.float32)
    _install_fake_rpy2_for_sct(
        monkeypatch,
        pearson_residuals=pearson_residuals,
        residual_variance=residual_variance,
        kept_genes=kept_genes,
    )

    ctx = DummyCtx(adata)
    result = await preprocess_data(
        "d25",
        ctx,
        PreprocessingParameters(
            normalization="sct",
            filter_mito_pct=None,
            remove_mito_genes=False,
        ),
    )

    assert result.n_genes == 110
    assert ctx.saved_adata is not None
    assert ctx.saved_adata.uns["sctransform"]["n_genes_before"] == 120
    assert ctx.saved_adata.uns["sctransform"]["n_genes_after"] == 110
    assert "sct_residual_variance" in ctx.saved_adata.var.columns


@pytest.mark.asyncio
async def test_preprocess_data_sct_dimension_mismatch_raises_processing_error(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "validate_r_package", lambda *_a, **_k: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0]),
    )

    adata = _make_adata(n_obs=10, n_vars=120)
    kept_genes = list(adata.var_names[:110])
    _install_fake_rpy2_for_sct(
        monkeypatch,
        pearson_residuals=np.ones((110, adata.n_obs), dtype=np.float32),
        residual_variance=np.ones(100, dtype=np.float32),
        kept_genes=kept_genes,
    )
    ctx = DummyCtx(adata)

    with pytest.raises(ProcessingError, match="Dimension mismatch after SCTransform"):
        await preprocess_data(
            "d26",
            ctx,
            PreprocessingParameters(normalization="sct", filter_mito_pct=None),
        )


@pytest.mark.asyncio
async def test_preprocess_data_sct_memory_error_has_actionable_message(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "validate_r_package", lambda *_a, **_k: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0]),
    )
    monkeypatch.setattr(
        preprocessing_mod.scipy.sparse,
        "csc_matrix",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(MemoryError("oom")),
    )

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    with pytest.raises(MemoryError, match="Memory error for SCTransform"):
        await preprocess_data(
            "d27",
            ctx,
            PreprocessingParameters(normalization="sct", filter_mito_pct=None),
        )


@pytest.mark.asyncio
async def test_preprocess_data_pearson_residuals_memory_error_has_guidance(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(
        preprocessing_mod.sc,
        "experimental",
        SimpleNamespace(
            pp=SimpleNamespace(
                normalize_pearson_residuals=lambda _adata: (_ for _ in ()).throw(
                    MemoryError("oom")
                )
            )
        ),
        raising=False,
    )
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0]),
    )

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    with pytest.raises(MemoryError, match="Insufficient memory for Pearson residuals"):
        await preprocess_data(
            "d28",
            ctx,
            PreprocessingParameters(normalization="pearson_residuals", filter_mito_pct=None),
        )


@pytest.mark.asyncio
async def test_preprocess_data_pearson_residuals_runtime_error_is_wrapped(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(
        preprocessing_mod.sc,
        "experimental",
        SimpleNamespace(
            pp=SimpleNamespace(
                normalize_pearson_residuals=lambda _adata: (_ for _ in ()).throw(
                    RuntimeError("pearson boom")
                )
            )
        ),
        raising=False,
    )
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0]),
    )

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    with pytest.raises(
        ProcessingError,
        match="Pearson residuals normalization failed: pearson boom",
    ):
        await preprocess_data(
            "d29",
            ctx,
            PreprocessingParameters(normalization="pearson_residuals", filter_mito_pct=None),
        )


@pytest.mark.asyncio
async def test_preprocess_data_scvi_rejects_negative_values(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "require", lambda *_args, **_kwargs: None)
    monkeypatch.setitem(
        sys.modules,
        "scvi",
        SimpleNamespace(model=SimpleNamespace(SCVI=object)),
    )

    # Inject negative values so counts layer (created from X) triggers rejection
    adata = _make_adata(n_obs=12, n_vars=120)
    adata.X[0, 0] = -0.5
    ctx = DummyCtx(adata)
    with pytest.raises(DataError, match="requires non-negative count data"):
        await preprocess_data(
            "d30",
            ctx,
            PreprocessingParameters(normalization="scvi", filter_mito_pct=None),
        )


@pytest.mark.asyncio
async def test_preprocess_data_scale_dense_cleans_nan_and_inf(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _inject_dense_scale(adata, max_value=None):
        del max_value
        adata.X = np.array(
            [[np.nan, np.inf, -np.inf], [1.0, -2.0, 3.0]], dtype=np.float32
        )

    monkeypatch.setattr(preprocessing_mod.sc.pp, "scale", _inject_dense_scale)
    adata = _make_adata(n_obs=2, n_vars=3)
    ctx = DummyCtx(adata)
    result = await preprocess_data(
        "d31",
        ctx,
        PreprocessingParameters(
            normalization="log",
            filter_mito_pct=None,
            scale=True,
            scale_max_value=5.0,
            n_hvgs=2,
        ),
    )

    assert result.n_cells == 2
    assert ctx.saved_adata is not None
    assert np.isfinite(ctx.saved_adata.X).all()
    assert np.max(ctx.saved_adata.X) <= 5.0
    assert np.min(ctx.saved_adata.X) >= -5.0


@pytest.mark.asyncio
async def test_preprocess_data_scale_sparse_cleans_nan_and_inf(monkeypatch: pytest.MonkeyPatch):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _inject_sparse_scale(adata, max_value=None):
        del max_value
        mat = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 2.0]], dtype=np.float32))
        mat.data = np.array([np.nan, np.inf], dtype=np.float32)
        adata.X = mat

    monkeypatch.setattr(preprocessing_mod.sc.pp, "scale", _inject_sparse_scale)
    adata = _make_adata(n_obs=2, n_vars=2)
    ctx = DummyCtx(adata)
    await preprocess_data(
        "d32",
        ctx,
        PreprocessingParameters(
            normalization="log",
            filter_mito_pct=None,
            scale=True,
            scale_max_value=4.0,
            n_hvgs=1,
        ),
    )

    assert ctx.saved_adata is not None
    assert sp.issparse(ctx.saved_adata.X)
    assert np.isfinite(ctx.saved_adata.X.data).all()
    assert np.max(ctx.saved_adata.X.data) <= 4.0


@pytest.mark.asyncio
async def test_preprocess_data_scale_failure_warns_and_continues(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(
        preprocessing_mod.sc.pp,
        "scale",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("scale boom")),
    )

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    result = await preprocess_data(
        "d33",
        ctx,
        PreprocessingParameters(normalization="log", filter_mito_pct=None, scale=True),
    )

    assert result.n_cells == 12
    assert any("Scaling failed: scale boom" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_preprocess_data_sct_sparse_input_without_gene_filtering(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    
    def _sparse_safe_qc(adata, qc_vars=None, percent_top=None, inplace=True):
        del qc_vars, percent_top, inplace
        counts = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        adata.obs["n_genes_by_counts"] = (counts > 0).sum(axis=1)
        adata.obs["total_counts"] = counts.sum(axis=1)
        if "mt" in adata.var.columns:
            mt_mask = adata.var["mt"].to_numpy()
            mt_counts = counts[:, mt_mask].sum(axis=1)
            total = np.clip(counts.sum(axis=1), 1e-9, None)
            adata.obs["pct_counts_mt"] = mt_counts / total * 100.0

    monkeypatch.setattr(preprocessing_mod.sc.pp, "calculate_qc_metrics", _sparse_safe_qc)
    monkeypatch.setattr(preprocessing_mod, "validate_r_package", lambda *_a, **_k: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0]),
    )

    adata = _make_adata(n_obs=14, n_vars=120)
    adata.X = sp.csr_matrix(adata.X)
    kept_genes = list(adata.var_names)
    pearson_residuals = np.full((adata.n_vars, adata.n_obs), 0.2, dtype=np.float32)
    residual_variance = np.linspace(0.0, 1.0, adata.n_vars, dtype=np.float32)
    _install_fake_rpy2_for_sct(
        monkeypatch,
        pearson_residuals=pearson_residuals,
        residual_variance=residual_variance,
        kept_genes=kept_genes,
    )

    ctx = DummyCtx(adata)
    result = await preprocess_data(
        "d34",
        ctx,
        PreprocessingParameters(
            normalization="sct",
            filter_mito_pct=None,
            remove_mito_genes=False,
        ),
    )

    assert result.n_genes == 120
    assert ctx.saved_adata is not None
    assert ctx.saved_adata.uns["sctransform"]["n_genes_filtered_by_sct"] == 0


@pytest.mark.asyncio
async def test_preprocess_data_sct_selects_expected_top_hvgs(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "validate_r_package", lambda *_a, **_k: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0]),
    )

    adata = _make_adata(n_obs=14, n_vars=120)
    kept_genes = list(adata.var_names[:110])
    pearson_residuals = np.full((110, adata.n_obs), 0.3, dtype=np.float32)
    residual_variance = np.linspace(0.0, 1.0, 110, dtype=np.float32)
    _install_fake_rpy2_for_sct(
        monkeypatch,
        pearson_residuals=pearson_residuals,
        residual_variance=residual_variance,
        kept_genes=kept_genes,
    )

    ctx = DummyCtx(adata)
    await preprocess_data(
        "d34_top_hvg",
        ctx,
        PreprocessingParameters(
            normalization="sct",
            filter_mito_pct=None,
            remove_mito_genes=False,
            sct_var_features_n=100,
        ),
    )

    assert ctx.saved_adata is not None
    hvg_mask = ctx.saved_adata.var["highly_variable"].to_numpy()
    assert int(hvg_mask.sum()) == 100
    assert not np.any(hvg_mask[:10])
    assert np.all(hvg_mask[-100:])


@pytest.mark.asyncio
async def test_preprocess_data_scvi_normalized_expression_uses_values_attribute(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(preprocessing_mod, "require", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        preprocessing_mod,
        "sample_expression_values",
        lambda _adata: np.array([0.0, 1.0, 2.0]),
    )

    class _FakeSCVI:
        @staticmethod
        def setup_anndata(adata, layer=None, batch_key=None):
            del adata, layer, batch_key
            return None

        def __init__(self, adata, **kwargs):
            self._adata = adata
            del kwargs

        def train(self, **kwargs):
            del kwargs
            return None

        def get_latent_representation(self):
            return np.ones((self._adata.n_obs, 2), dtype=float)

        def get_normalized_expression(self, library_size=1e4):
            del library_size
            return SimpleNamespace(
                values=np.full((self._adata.n_obs, self._adata.n_vars), 3.0, dtype=float)
            )

    monkeypatch.setitem(sys.modules, "scvi", SimpleNamespace(model=SimpleNamespace(SCVI=_FakeSCVI)))

    adata = _make_adata(n_obs=12, n_vars=120)
    ctx = DummyCtx(adata)
    await preprocess_data(
        "d35",
        ctx,
        PreprocessingParameters(normalization="scvi", filter_mito_pct=None),
    )

    assert ctx.saved_adata is not None
    assert np.isfinite(ctx.saved_adata.X).all()


@pytest.mark.asyncio
async def test_preprocess_data_hvg_selection_failure_is_wrapped(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    monkeypatch.setattr(
        preprocessing_mod.sc.pp,
        "highly_variable_genes",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("hvg boom")),
    )

    adata = _make_adata(n_obs=14, n_vars=120)
    ctx = DummyCtx(adata)

    with pytest.raises(ProcessingError, match="HVG selection failed: hvg boom"):
        await preprocess_data(
            "d36",
            ctx,
            PreprocessingParameters(normalization="log", filter_mito_pct=None),
        )


@pytest.mark.asyncio
async def test_preprocess_data_hvg_nan_bin_keyerror_falls_back_to_variance(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _raise_nan_bin_keyerror(*_args, **_kwargs):
        raise KeyError("[nan] not in index")

    monkeypatch.setattr(
        preprocessing_mod.sc.pp,
        "highly_variable_genes",
        _raise_nan_bin_keyerror,
    )

    adata = _make_adata(n_obs=14, n_vars=120)
    ctx = DummyCtx(adata)

    await preprocess_data(
        "d36b",
        ctx,
        PreprocessingParameters(
            normalization="log",
            filter_mito_pct=None,
            n_hvgs=20,
            remove_mito_genes=False,
            remove_ribo_genes=False,
        ),
    )

    assert ctx.saved_adata is not None
    assert int(ctx.saved_adata.var["highly_variable"].sum()) == 20
    assert any("NaN bins" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_preprocess_data_remove_ribo_genes_drops_ribo_from_hvgs(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)
    adata = _make_adata(n_obs=14, n_vars=120)
    ctx = DummyCtx(adata)

    await preprocess_data(
        "d37",
        ctx,
        PreprocessingParameters(
            normalization="log",
            filter_mito_pct=None,
            remove_mito_genes=False,
            remove_ribo_genes=True,
            n_hvgs=20,
        ),
    )

    assert ctx.saved_adata is not None
    assert bool(ctx.saved_adata.var.loc["RPS3", "ribo"]) is True
    assert bool(ctx.saved_adata.var.loc["RPS3", "highly_variable"]) is False


@pytest.mark.asyncio
async def test_preprocess_data_gene_subsample_requires_hvg_column_presence(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_lightweight_preprocess_mocks(monkeypatch)

    def _drop_hvg_column(adata, n_top_genes=2000):
        del n_top_genes
        if "highly_variable" in adata.var.columns:
            del adata.var["highly_variable"]

    monkeypatch.setattr(preprocessing_mod.sc.pp, "highly_variable_genes", _drop_hvg_column)

    adata = _make_adata(n_obs=16, n_vars=120)
    ctx = DummyCtx(adata)

    with pytest.raises(ProcessingError, match="Gene subsampling failed: no HVGs identified"):
        await preprocess_data(
            "d38",
            ctx,
            PreprocessingParameters(
                normalization="log",
                filter_mito_pct=None,
                subsample_genes=20,
                n_hvgs=20,
                remove_mito_genes=False,
                remove_ribo_genes=False,
            ),
        )
