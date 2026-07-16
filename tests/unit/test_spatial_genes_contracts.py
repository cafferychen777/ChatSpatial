"""Unit tests for spatial_genes core contracts with lightweight dependency stubs."""

from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import SpatialVariableGenesParameters
from chatspatial.tools import spatial_genes as sg
from chatspatial.utils.dependency_manager import REnvironment
from chatspatial.utils.exceptions import (
    DataError,
    DataNotFoundError,
    DependencyError,
    ProcessingError,
)


def _required_module(name: str, *_args, **_kwargs):
    """Return the fake backend installed by each test."""
    module_name = {
        "spatialde": "SpatialDE",
        "naivede": "NaiveDE",
    }.get(name, name)
    return sys.modules.get(module_name, object())


@pytest.fixture(autouse=True)
def _resolve_fake_dependency_submodules(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        sg,
        "require_module",
        lambda _name, module_name, *_args, **_kwargs: sys.modules[module_name],
    )


class DummyCtx:
    def __init__(self):
        self.warnings: list[str] = []

    async def warning(self, msg: str):
        self.warnings.append(msg)


def _install_fake_rpy2(monkeypatch: pytest.MonkeyPatch) -> None:
    class _Lock:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    fake_ro = ModuleType("rpy2.robjects")
    fake_ro.r = {}
    fake_ro.default_converter = object()

    fake_conversion = ModuleType("rpy2.robjects.conversion")

    class _LocalConverter:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return None

        def __exit__(self, exc_type, exc, tb):
            return False

    fake_conversion.localconverter = _LocalConverter
    fake_conversion.default_converter = object()

    fake_packages = ModuleType("rpy2.robjects.packages")
    fake_packages.importr = lambda *_args, **_kwargs: None

    fake_rinterface_lib = ModuleType("rpy2.rinterface_lib")
    fake_rinterface_lib.openrlib = SimpleNamespace(rlock=_Lock())

    monkeypatch.setitem(__import__("sys").modules, "rpy2", ModuleType("rpy2"))
    monkeypatch.setitem(__import__("sys").modules, "rpy2.robjects", fake_ro)
    monkeypatch.setitem(
        __import__("sys").modules, "rpy2.robjects.conversion", fake_conversion
    )
    monkeypatch.setitem(
        __import__("sys").modules, "rpy2.robjects.packages", fake_packages
    )
    monkeypatch.setitem(
        __import__("sys").modules, "rpy2.rinterface_lib", fake_rinterface_lib
    )

    environment = REnvironment(
        robjects=fake_ro,
        pandas2ri=SimpleNamespace(converter=object()),
        numpy2ri=SimpleNamespace(converter=object()),
        packages=fake_packages,
        conversion=fake_conversion,
        openrlib=fake_rinterface_lib.openrlib,
        anndata2ri=None,
    )
    monkeypatch.setattr(
        sg,
        "validate_r_environment",
        lambda *_args, **_kwargs: environment,
    )


def _install_fake_rpy2_runtime(monkeypatch: pytest.MonkeyPatch, spark_factory) -> None:
    class _Lock:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class _FakeRMatrix:
        def __init__(self, vec, nrow: int, ncol: int, byrow: bool = True):
            self.values = list(vec)
            self.nrow = nrow
            self.ncol = ncol
            self.byrow = byrow
            self.rownames = []
            self.colnames = []

    class _FakeR:
        def __init__(self):
            self._funcs = {
                "data.frame": lambda **kwargs: kwargs,
                "is.data.frame": lambda _obj: [True],
                "$": lambda obj, key: obj.get(key) if isinstance(obj, dict) else None,
            }
            self.matrix = lambda vec, nrow, ncol, byrow=True: _FakeRMatrix(
                vec, nrow, ncol, byrow
            )

        def __getitem__(self, key: str):
            return self._funcs[key]

    fake_r = _FakeR()
    fake_ro = ModuleType("rpy2.robjects")
    fake_ro.r = fake_r
    fake_ro.NULL = None
    fake_ro.default_converter = object()
    fake_ro.IntVector = lambda x: list(x)
    fake_ro.FloatVector = lambda x: list(x)
    fake_ro.StrVector = lambda x: list(x)

    fake_conversion = ModuleType("rpy2.robjects.conversion")

    class _LocalConverter:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return None

        def __exit__(self, exc_type, exc, tb):
            return False

    fake_conversion.localconverter = _LocalConverter
    fake_conversion.default_converter = object()
    fake_ro.conversion = fake_conversion

    fake_packages = ModuleType("rpy2.robjects.packages")
    fake_packages.importr = lambda _name: spark_factory()

    fake_rinterface_lib = ModuleType("rpy2.rinterface_lib")
    fake_rinterface_lib.openrlib = SimpleNamespace(rlock=_Lock())

    monkeypatch.setitem(__import__("sys").modules, "rpy2", ModuleType("rpy2"))
    monkeypatch.setitem(__import__("sys").modules, "rpy2.robjects", fake_ro)
    monkeypatch.setitem(
        __import__("sys").modules, "rpy2.robjects.conversion", fake_conversion
    )
    monkeypatch.setitem(
        __import__("sys").modules, "rpy2.robjects.packages", fake_packages
    )
    monkeypatch.setitem(
        __import__("sys").modules, "rpy2.rinterface_lib", fake_rinterface_lib
    )

    environment = REnvironment(
        robjects=fake_ro,
        pandas2ri=SimpleNamespace(converter=object()),
        numpy2ri=SimpleNamespace(converter=object()),
        packages=fake_packages,
        conversion=fake_conversion,
        openrlib=fake_rinterface_lib.openrlib,
        anndata2ri=None,
        loaded_packages={"SPARK": spark_factory()},
    )
    monkeypatch.setattr(
        sg,
        "validate_r_environment",
        lambda *_args, **_kwargs: environment,
    )


@pytest.mark.asyncio
async def test_spatialde_applies_compatibility_patch_before_import(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    events: list[str] = []

    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spatialde_compat",
        lambda: events.append("compat"),
    )

    def _stop_after_first_import(name: str, *_args, **_kwargs):
        events.append(name)
        raise DependencyError("stop after import order check")

    monkeypatch.setattr(sg, "require", _stop_after_first_import)

    with pytest.raises(DependencyError, match="import order check"):
        await sg._identify_spatial_genes_spatialde(
            "d1",
            minimal_spatial_adata.copy(),
            SpatialVariableGenesParameters(method="spatialde"),
            DummyCtx(),
        )

    assert events == ["compat", "spatialde"]


@pytest.mark.asyncio
async def test_spatialde_success_stores_var_outputs_and_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    captured: dict[str, object] = {}

    fake_naivede = ModuleType("NaiveDE")
    fake_naivede.stabilize = lambda x: x
    fake_naivede.regress_out = lambda _tc, expr_t, _formula: expr_t

    fake_spatialde = ModuleType("SpatialDE")
    fake_spatialde.run = lambda _coords, _expr: pd.DataFrame(
        {"g": ["gene_0", "gene_1"], "pval": [0.001, 0.02], "l": [1.2, 0.8]}
    )
    fake_spatialde_util = ModuleType("SpatialDE.util")
    fake_spatialde_util.qvalue = lambda pvals, pi0=None: np.array([0.01, 0.04])

    monkeypatch.setitem(__import__("sys").modules, "NaiveDE", fake_naivede)
    monkeypatch.setitem(__import__("sys").modules, "SpatialDE", fake_spatialde)
    monkeypatch.setitem(
        __import__("sys").modules, "SpatialDE.util", fake_spatialde_util
    )

    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spatialde_compat",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_spatialde(
        "d1",
        adata,
        SpatialVariableGenesParameters(method="spatialde", spatial_key="spatial"),
        DummyCtx(),
    )
    assert out.method == "spatialde"
    assert out.n_significant_genes == 2
    assert "spatialde_pval" in adata.var.columns
    assert "spatialde_qval" in adata.var.columns
    assert captured["analysis_name"] == "spatial_genes_spatialde"
    assert captured["results_keys"]["var"] == [
        "spatialde_pval",
        "spatialde_qval",
        "spatialde_l",
    ]
    # results_key must point to a real metadata key, not a phantom uns key
    assert out.results_key == "spatial_genes_spatialde_metadata"


@pytest.mark.asyncio
async def test_flashs_success_stores_var_outputs_and_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    captured: dict[str, object] = {}

    class _FakeResult:
        def __init__(self, genes: list[str]):
            n = len(genes)
            self.gene_names = genes
            self.pvalues = np.linspace(0.001, 0.8, n)
            self.qvalues = np.linspace(0.01, 0.9, n)
            self.statistics = np.linspace(1.0, 2.0, n)
            self.effect_size = np.linspace(0.5, 1.5, n)
            self.pvalues_binary = np.linspace(0.002, 0.7, n)
            self.pvalues_rank = np.linspace(0.003, 0.6, n)
            self.n_expressed = np.arange(10, 10 + n)
            self.tested_mask = np.array([True] * (n - 1) + [False], dtype=bool)
            self.n_tested = n - 1
            self.n_significant = 1

    class _FakeFlashS:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def fit_test(self, coords, X, gene_names):
            del coords, X
            return _FakeResult(gene_names)

    fake_flashs = ModuleType("flashs")
    fake_flashs.FlashS = _FakeFlashS
    monkeypatch.setitem(__import__("sys").modules, "flashs", fake_flashs)

    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="current",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_flashs(
        "d_flashs",
        adata,
        SpatialVariableGenesParameters(method="flashs", spatial_key="spatial"),
        DummyCtx(),
    )

    assert out.method == "flashs"
    assert out.n_genes_analyzed == adata.n_vars - 1
    assert "flashs_pval" in adata.var.columns
    assert "flashs_qval" in adata.var.columns
    assert "flashs_statistic" in adata.var.columns
    assert "flashs_effect_size" in adata.var.columns
    assert "flashs_pval_binary" in adata.var.columns
    assert "flashs_pval_rank" in adata.var.columns
    assert "flashs_n_expressed" in adata.var.columns
    assert "flashs_tested" in adata.var.columns
    assert captured["analysis_name"] == "spatial_genes_flashs"
    assert captured["method"] == "flashs"
    assert "flashs_qval" in captured["results_keys"]["var"]
    assert out.results_key == "spatial_genes_flashs_metadata"


@pytest.mark.asyncio
async def test_flashs_missing_dependency_raises_dependency_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    def _raise_dependency(*_args, **_kwargs):
        raise DependencyError("flashs dependency missing")

    monkeypatch.setattr(sg, "require", _raise_dependency)

    with pytest.raises(DependencyError, match="flashs dependency missing"):
        await sg._identify_spatial_genes_flashs(
            "d_missing",
            adata,
            SpatialVariableGenesParameters(method="flashs", spatial_key="spatial"),
            DummyCtx(),
        )


@pytest.mark.asyncio
async def test_sparkx_requires_hvg_column_when_test_only_hvg_enabled(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    _install_fake_rpy2(monkeypatch)
    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )

    with pytest.raises(DataError, match="Highly variable genes marker"):
        await sg._identify_spatial_genes_sparkx(
            "d1",
            adata,
            SpatialVariableGenesParameters(method="sparkx", test_only_hvg=True),
            DummyCtx(),
        )


@pytest.mark.asyncio
async def test_sparkx_raises_when_no_hvgs_found(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.var["highly_variable"] = False
    _install_fake_rpy2(monkeypatch)
    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )

    with pytest.raises(DataNotFoundError, match="No HVGs found"):
        await sg._identify_spatial_genes_sparkx(
            "d1",
            adata,
            SpatialVariableGenesParameters(method="sparkx", test_only_hvg=True),
            DummyCtx(),
        )


@pytest.mark.asyncio
async def test_sparkx_missing_r_package_raises_dependency_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    _install_fake_rpy2(monkeypatch)
    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )

    monkeypatch.setattr(
        sg,
        "validate_r_environment",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("R package 'SPARK' not installed")
        ),
    )

    with pytest.raises(DependencyError, match="R package 'SPARK' not installed"):
        await sg._identify_spatial_genes_sparkx(
            "d2",
            adata,
            SpatialVariableGenesParameters(
                method="sparkx", spatial_key="spatial", test_only_hvg=False
            ),
            DummyCtx(),
        )


@pytest.mark.asyncio
async def test_spatialde_warns_for_large_gene_set_runtime(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    import anndata as ad

    base = minimal_spatial_adata.copy()
    ctx = DummyCtx()

    X_big = np.tile(np.asarray(base.X), (1, 260)).astype(np.float32)
    adata = ad.AnnData(X_big)
    adata.obs_names = base.obs_names.copy()
    adata.var_names = [f"gene_{i}" for i in range(adata.n_vars)]
    adata.obsm["spatial"] = np.asarray(base.obsm["spatial"]).copy()

    fake_naivede = ModuleType("NaiveDE")
    fake_naivede.stabilize = lambda x: x
    fake_naivede.regress_out = lambda _tc, expr_t, _formula: expr_t

    fake_spatialde = ModuleType("SpatialDE")
    fake_spatialde.run = lambda _coords, _expr: pd.DataFrame(
        {"g": ["gene_0", "gene_1"], "pval": [0.001, 0.02], "l": [1.0, 0.5]}
    )
    fake_spatialde_util = ModuleType("SpatialDE.util")
    fake_spatialde_util.qvalue = lambda pvals, pi0=None: np.array([0.01, 0.04])

    monkeypatch.setitem(__import__("sys").modules, "NaiveDE", fake_naivede)
    monkeypatch.setitem(__import__("sys").modules, "SpatialDE", fake_spatialde)
    monkeypatch.setitem(
        __import__("sys").modules, "SpatialDE.util", fake_spatialde_util
    )

    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spatialde_compat",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_spatialde(
        "d3",
        adata,
        SpatialVariableGenesParameters(method="spatialde", spatial_key="spatial"),
        ctx,
    )

    assert out.method == "spatialde"
    assert any("may take" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_spatialde_prefers_hvgs_and_passes_pi0_to_qvalue(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata[:, :6].copy()
    adata.var["highly_variable"] = [True, False, True, True, False, True]
    captured: dict[str, object] = {}

    fake_naivede = ModuleType("NaiveDE")
    fake_naivede.stabilize = lambda x: x
    fake_naivede.regress_out = lambda _tc, expr_t, _formula: expr_t

    def _fake_run(_coords, expr):
        captured["genes_in_run"] = list(expr.columns)
        genes = list(expr.columns)
        return pd.DataFrame(
            {
                "g": genes,
                "pval": np.linspace(0.001, 0.04, len(genes)),
                "l": np.linspace(0.1, 1.0, len(genes)),
            }
        )

    fake_spatialde = ModuleType("SpatialDE")
    fake_spatialde.run = _fake_run
    fake_spatialde_util = ModuleType("SpatialDE.util")

    def _fake_qvalue(pvals, pi0=None):
        captured["pi0"] = pi0
        return np.linspace(0.01, 0.04, len(pvals))

    fake_spatialde_util.qvalue = _fake_qvalue

    monkeypatch.setitem(__import__("sys").modules, "NaiveDE", fake_naivede)
    monkeypatch.setitem(__import__("sys").modules, "SpatialDE", fake_spatialde)
    monkeypatch.setitem(
        __import__("sys").modules, "SpatialDE.util", fake_spatialde_util
    )

    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spatialde_compat",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_spatialde(
        "spatialde_hvg",
        adata,
        SpatialVariableGenesParameters(
            method="spatialde",
            spatial_key="spatial",
            n_top_genes=3,
            spatialde_pi0=0.7,
        ),
        DummyCtx(),
    )

    assert out.n_genes_analyzed == 3
    assert captured["pi0"] == 0.7
    assert captured["genes_in_run"] == ["gene_0", "gene_2", "gene_3"]


@pytest.mark.asyncio
async def test_spatialde_falls_back_to_expression_when_hvgs_insufficient(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata[:, :5].copy()
    adata.var["highly_variable"] = [True, False, False, False, False]
    adata.X = np.tile(
        np.asarray([[10, 9, 8, 1, 0]], dtype=np.float32), (adata.n_obs, 1)
    )
    captured: dict[str, object] = {}

    fake_naivede = ModuleType("NaiveDE")
    fake_naivede.stabilize = lambda x: x
    fake_naivede.regress_out = lambda _tc, expr_t, _formula: expr_t

    def _fake_run(_coords, expr):
        captured["genes_in_run"] = list(expr.columns)
        genes = list(expr.columns)
        return pd.DataFrame(
            {
                "g": genes,
                "pval": [0.001, 0.002, 0.1],
                "l": [0.2, 0.3, 0.4],
            }
        )

    fake_spatialde = ModuleType("SpatialDE")
    fake_spatialde.run = _fake_run
    fake_spatialde_util = ModuleType("SpatialDE.util")
    fake_spatialde_util.qvalue = lambda pvals, pi0=None: np.asarray([0.01, 0.02, 0.2])

    monkeypatch.setitem(__import__("sys").modules, "NaiveDE", fake_naivede)
    monkeypatch.setitem(__import__("sys").modules, "SpatialDE", fake_spatialde)
    monkeypatch.setitem(
        __import__("sys").modules, "SpatialDE.util", fake_spatialde_util
    )
    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spatialde_compat",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_spatialde(
        "spatialde_fallback",
        adata,
        SpatialVariableGenesParameters(
            method="spatialde",
            spatial_key="spatial",
            n_top_genes=3,
        ),
        DummyCtx(),
    )

    assert out.n_genes_analyzed == 3
    assert captured["genes_in_run"] == ["gene_0", "gene_1", "gene_2"]


@pytest.mark.asyncio
async def test_spatialde_selects_by_expression_without_hvg_column(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata[:, :4].copy()
    if "highly_variable" in adata.var:
        adata.var = adata.var.drop(columns=["highly_variable"])
    adata.X = np.tile(np.asarray([[7, 6, 1, 0]], dtype=np.float32), (adata.n_obs, 1))
    captured: dict[str, object] = {}

    fake_naivede = ModuleType("NaiveDE")
    fake_naivede.stabilize = lambda x: x
    fake_naivede.regress_out = lambda _tc, expr_t, _formula: expr_t

    def _fake_run(_coords, expr):
        captured["genes_in_run"] = list(expr.columns)
        genes = list(expr.columns)
        return pd.DataFrame(
            {
                "g": genes,
                "pval": [0.001, 0.03],
                "l": [0.2, 0.8],
            }
        )

    fake_spatialde = ModuleType("SpatialDE")
    fake_spatialde.run = _fake_run
    fake_spatialde_util = ModuleType("SpatialDE.util")
    fake_spatialde_util.qvalue = lambda pvals, pi0=None: np.asarray([0.01, 0.04])

    monkeypatch.setitem(__import__("sys").modules, "NaiveDE", fake_naivede)
    monkeypatch.setitem(__import__("sys").modules, "SpatialDE", fake_spatialde)
    monkeypatch.setitem(
        __import__("sys").modules, "SpatialDE.util", fake_spatialde_util
    )
    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        "chatspatial.utils.compat.ensure_spatialde_compat",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X,
            var_names=_adata.var_names,
            source="raw",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_spatialde(
        "spatialde_no_hvg",
        adata,
        SpatialVariableGenesParameters(
            method="spatialde",
            spatial_key="spatial",
            n_top_genes=2,
        ),
        DummyCtx(),
    )

    assert out.n_genes_analyzed == 2
    assert captured["genes_in_run"] == ["gene_0", "gene_1"]


@pytest.mark.asyncio
async def test_flashs_reindexs_to_adata_var_and_defaults_tested_mask(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class _FakeResult:
        def __init__(self):
            self.pvalues = np.asarray([0.01, 0.2, 0.03], dtype=float)
            self.qvalues = np.asarray([0.02, 0.3, 0.04], dtype=float)
            self.statistics = np.asarray([1.0, 0.5, 0.8], dtype=float)
            self.effect_size = np.asarray([0.9, 0.1, 0.7], dtype=float)
            self.pvalues_binary = np.asarray([0.01, 0.2, 0.03], dtype=float)
            self.pvalues_rank = np.asarray([0.01, 0.2, 0.03], dtype=float)
            self.n_expressed = np.asarray([10, 20, 30], dtype=int)
            self.tested_mask = None
            self.n_tested = 3
            self.n_significant = 2

    class _FakeFlashS:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def fit_test(self, coords, X, gene_names):
            del coords, X, gene_names
            return _FakeResult()

    fake_flashs = ModuleType("flashs")
    fake_flashs.FlashS = _FakeFlashS
    monkeypatch.setitem(__import__("sys").modules, "flashs", fake_flashs)

    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=_adata.X[:, :3],
            var_names=pd.Index(["gene_1", "gene_0", "not_in_var"]),
            source="current",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_flashs(
        "d_flashs_reindex",
        adata,
        SpatialVariableGenesParameters(method="flashs", spatial_key="spatial"),
        DummyCtx(),
    )

    assert out.n_genes_analyzed == 3
    assert float(adata.var.loc["gene_2", "flashs_pval"]) == 1.0
    assert bool(adata.var.loc["gene_1", "flashs_tested"]) is True
    assert bool(adata.var.loc["gene_2", "flashs_tested"]) is False
    assert int(adata.var.loc["gene_0", "flashs_n_expressed"]) == 20


@pytest.mark.asyncio
async def test_sparkx_success_covers_filtering_and_housekeeping_warning(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata[:, :6].copy()
    adata.var_names = ["MT-A", "RPS1", "ACTB", "GAPDH", "GENE1", "GENE2"]
    ctx = DummyCtx()
    captured: dict[str, object] = {}

    class _SparkPkg:
        @staticmethod
        def sparkx(count_in, locus_in, X_in, numCores, option, verbose):
            del locus_in, X_in, numCores, option, verbose
            n = count_in.nrow
            return SimpleNamespace(
                rx2=lambda key: {
                    "combinedPval": np.linspace(0.001, 0.02, n),
                    "adjustedPval": np.linspace(0.002, 0.04, n),
                }
                if key == "res_mtest"
                else None
            )

    _install_fake_rpy2_runtime(monkeypatch, spark_factory=lambda: _SparkPkg())

    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=np.asarray(_adata.X, dtype=float),
            var_names=_adata.var_names,
            source="raw",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_sparkx(
        "spark_ok",
        adata,
        SpatialVariableGenesParameters(
            method="sparkx",
            spatial_key="spatial",
            test_only_hvg=False,
            filter_mt_genes=True,
            filter_ribo_genes=True,
            warn_housekeeping=True,
            sparkx_percentage=0.01,
            sparkx_min_total_counts=1,
        ),
        ctx,
    )

    assert out.method == "sparkx"
    assert out.n_genes_analyzed == 4
    assert "sparkx_pval" in adata.var.columns
    assert "sparkx_qval" in adata.var.columns
    assert captured["statistics"]["n_genes_analyzed"] == 4
    assert any("Housekeeping gene dominance detected" in w for w in ctx.warnings)
    assert out.results_key == "spatial_genes_sparkx_metadata"


@pytest.mark.asyncio
async def test_sparkx_success_with_hvg_only_branch_and_low_result_warning(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata[:, :4].copy()
    adata.var_names = ["G1", "G2", "G3", "G4"]
    adata.var["highly_variable"] = [True, True, True, False]
    ctx = DummyCtx()

    class _SparkPkg:
        @staticmethod
        def sparkx(count_in, locus_in, X_in, numCores, option, verbose):
            del locus_in, X_in, numCores, option, verbose
            # Return fewer rows than input genes to trigger quality warning branch.
            return SimpleNamespace(
                rx2=lambda key: {
                    "combinedPval": [0.01],
                    "adjustedPval": [0.02],
                }
                if key == "res_mtest"
                else None
            )

    _install_fake_rpy2_runtime(monkeypatch, spark_factory=lambda: _SparkPkg())
    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=np.asarray(_adata.X, dtype=float),
            var_names=_adata.var_names,
            source="raw",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    out = await sg._identify_spatial_genes_sparkx(
        "spark_hvg",
        adata,
        SpatialVariableGenesParameters(
            method="sparkx",
            spatial_key="spatial",
            test_only_hvg=True,
            filter_mt_genes=False,
            filter_ribo_genes=False,
            sparkx_percentage=0.01,
            sparkx_min_total_counts=1,
        ),
        ctx,
    )

    assert out.n_genes_analyzed == 1
    assert any("returned results for only" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_sparkx_hvg_no_overlap_raises_data_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata[:, :4].copy()
    adata.var["highly_variable"] = [True, False, True, False]
    _install_fake_rpy2(monkeypatch)
    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=np.asarray(_adata.X[:, :4], dtype=float),
            var_names=pd.Index(["X1", "X2", "X3", "X4"]),
            source="raw",
        ),
    )

    with pytest.raises(DataError, match="no overlap found"):
        await sg._identify_spatial_genes_sparkx(
            "spark_no_overlap",
            adata,
            SpatialVariableGenesParameters(
                method="sparkx",
                spatial_key="spatial",
                test_only_hvg=True,
                sparkx_percentage=0.01,
                sparkx_min_total_counts=1,
            ),
            DummyCtx(),
        )


@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("res_mtest", "is_dataframe", "error_match"),
    [
        (None, True, "returned None for res_mtest"),
        (
            {"combinedPval": [0.01], "adjustedPval": [0.02]},
            False,
            "output format error",
        ),
        ({"adjustedPval": [0.02]}, True, "missing 'combinedPval'"),
        ({"combinedPval": [0.01]}, True, "missing 'adjustedPval'"),
    ],
)
async def test_sparkx_invalid_output_formats_raise_processing_error(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
    res_mtest,
    is_dataframe: bool,
    error_match: str,
):
    import sys

    adata = minimal_spatial_adata[:, :4].copy()

    class _SparkPkg:
        @staticmethod
        def sparkx(count_in, locus_in, X_in, numCores, option, verbose):
            del count_in, locus_in, X_in, numCores, option, verbose
            return SimpleNamespace(
                rx2=lambda key: res_mtest if key == "res_mtest" else None
            )

    _install_fake_rpy2_runtime(monkeypatch, spark_factory=lambda: _SparkPkg())
    fake_ro = sys.modules["rpy2.robjects"]
    fake_ro.r._funcs["is.data.frame"] = lambda _obj: [is_dataframe]

    monkeypatch.setattr(sg, "require", _required_module)
    monkeypatch.setattr(
        sg,
        "get_raw_data_source",
        lambda _adata, **_kw: SimpleNamespace(
            X=np.asarray(_adata.X, dtype=float),
            var_names=_adata.var_names,
            source="raw",
        ),
    )
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    with pytest.raises(ProcessingError, match=error_match):
        await sg._identify_spatial_genes_sparkx(
            "spark_bad_format",
            adata,
            SpatialVariableGenesParameters(
                method="sparkx",
                spatial_key="spatial",
                test_only_hvg=False,
                sparkx_percentage=0.01,
                sparkx_min_total_counts=1,
            ),
            DummyCtx(),
        )
