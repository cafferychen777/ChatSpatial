"""Unit tests for lightweight sc-type helper utilities in annotation module."""

from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import AnnotationParameters
from chatspatial.tools import annotation as ann
from chatspatial.utils.dependency_manager import REnvironment
from chatspatial.utils.exceptions import DataError


class DummyCtx:
    async def warning(self, _msg: str):
        return None


class _FakeRLock:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


class _FakeLocalConverter:
    def __enter__(self):
        return None

    def __exit__(self, exc_type, exc, tb):
        return False


class _FakeConverter:
    def __add__(self, _other):
        return self


def _mock_validate_r_environment(
    monkeypatch: pytest.MonkeyPatch,
    r_obj: object,
    *,
    converter: object | None = None,
) -> None:
    base_converter = converter or _FakeConverter()
    robjects = SimpleNamespace(
        r=r_obj,
        default_converter=base_converter,
        StrVector=lambda values: list(values),
    )
    conversion = SimpleNamespace(
        localconverter=lambda _converter: _FakeLocalConverter()
    )
    environment = REnvironment(
        robjects=robjects,
        pandas2ri=SimpleNamespace(converter=base_converter),
        numpy2ri=SimpleNamespace(converter=base_converter),
        packages=SimpleNamespace(importr=lambda _package: object()),
        conversion=conversion,
        openrlib=SimpleNamespace(rlock=_FakeRLock()),
        anndata2ri=SimpleNamespace(converter=base_converter),
    )
    monkeypatch.setattr(
        ann,
        "validate_r_environment",
        lambda _ctx, **_kwargs: environment,
    )


def test_softmax_is_stable_and_normalized():
    arr = np.array([1000.0, 1001.0, 999.0], dtype=float)
    out = ann._softmax(arr)
    assert np.isfinite(out).all()
    assert np.isclose(out.sum(), 1.0)
    assert out.argmax() == 1


def test_assign_sctype_celltypes_and_unknown_handling():
    scores = pd.DataFrame(
        {
            "cell_1": [2.0, 0.1],
            "cell_2": [-1.0, -0.2],
        },
        index=["T", "B"],
    )

    types, conf = ann._assign_sctype_celltypes(scores, DummyCtx())
    assert types[0] == "T"
    assert conf[0] > 0
    assert types[1] == "Unknown"
    assert conf[1] == 0.0


def test_assign_sctype_celltypes_rejects_empty_scores():
    with pytest.raises(DataError, match="Scores DataFrame is empty or None"):
        ann._assign_sctype_celltypes(pd.DataFrame(), DummyCtx())


def test_calculate_sctype_stats_counts_labels():
    out = ann._calculate_sctype_stats(["T", "B", "T", "Unknown"])
    assert out == {"T": 2, "B": 1, "Unknown": 1}


def test_get_sctype_cache_key_changes_with_params(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    p1 = AnnotationParameters(method="sctype", sctype_tissue="Liver")
    p2 = AnnotationParameters(method="sctype", sctype_tissue="Brain")
    k1 = ann._get_sctype_cache_key(adata, p1)
    k2 = ann._get_sctype_cache_key(adata, p2)
    assert k1 != k2


def test_load_cached_sctype_results_reads_json(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_DIR", tmp_path)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE", {})
    cache_key = "abc"
    payload = {
        "cell_types": ["T", "B"],
        "counts": {"T": 1, "B": 1},
        "confidence_by_celltype": {"T": 0.7, "B": 0.6},
        "mapping_score": None,
    }
    (tmp_path / f"{cache_key}.json").write_text(json.dumps(payload), encoding="utf-8")

    out = ann._load_cached_sctype_results(cache_key, DummyCtx())
    assert out is not None
    assert out[0] == ["T", "B"]
    assert out[1]["T"] == 1


@pytest.mark.asyncio
async def test_cache_sctype_results_writes_json(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_DIR", tmp_path)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE", {})
    cache_key = "k1"
    results = (["T"], {"T": 1}, {"T": 0.9}, 0.5)

    await ann._cache_sctype_results(cache_key, results, DummyCtx())

    f = tmp_path / f"{cache_key}.json"
    assert f.exists()
    data = json.loads(f.read_text(encoding="utf-8"))
    assert data["cell_types"] == ["T"]
    assert data["mapping_score"] == 0.5


@pytest.mark.asyncio
async def test_annotate_with_sctype_cache_hit_short_circuits_pipeline(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    params = AnnotationParameters(
        method="sctype", sctype_tissue="Brain", sctype_use_cache=True
    )

    cached_cell_types = ["T"] * (adata.n_obs // 2) + ["B"] * (
        adata.n_obs - adata.n_obs // 2
    )
    cached_counts = {"T": adata.n_obs // 2, "B": adata.n_obs - adata.n_obs // 2}

    def _validate_should_not_run(*_args, **_kwargs):
        raise AssertionError("R validation should not run on sc-type cache hit")

    monkeypatch.setattr(ann, "validate_r_environment", _validate_should_not_run)
    monkeypatch.setattr(ann, "_get_sctype_cache_key", lambda *_args, **_kwargs: "k1")
    monkeypatch.setattr(
        ann,
        "_load_cached_sctype_results",
        lambda *_args, **_kwargs: (
            cached_cell_types,
            cached_counts,
            {"T": 0.8, "B": 0.7},
            None,
        ),
    )

    def _should_not_run(*_args, **_kwargs):
        raise AssertionError("pipeline should not run on cache hit")

    monkeypatch.setattr(ann, "_load_sctype_functions", _should_not_run)
    monkeypatch.setattr(ann, "_prepare_sctype_genesets", _should_not_run)
    monkeypatch.setattr(ann, "_run_sctype_scoring", _should_not_run)

    out = await ann._annotate_with_sctype(
        adata,
        params,
        DummyCtx(),
        output_key="cell_type_sctype",
        confidence_key="confidence_sctype",
    )

    assert out.cell_types == ["T", "B"]
    assert out.counts == cached_counts
    assert adata.obs["cell_type_sctype"].dtype.name == "category"
    assert "confidence_sctype" in adata.obs


@pytest.mark.asyncio
async def test_annotate_with_sctype_cache_miss_preserves_cell_type_order_and_caches(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    params = AnnotationParameters(
        method="sctype", sctype_tissue="Brain", sctype_use_cache=True
    )

    captured: dict[str, object] = {}

    per_cell_types = ["B", "T", "B", "Unknown"] * (adata.n_obs // 4)
    per_cell_conf = [0.9, 0.8, 0.7, 0.0] * (adata.n_obs // 4)

    monkeypatch.setattr(
        ann, "validate_r_environment", lambda *_args, **_kwargs: object()
    )
    monkeypatch.setattr(ann, "_get_sctype_cache_key", lambda *_args, **_kwargs: "k2")
    monkeypatch.setattr(
        ann, "_load_cached_sctype_results", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(ann, "_load_sctype_functions", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(ann, "_prepare_sctype_genesets", lambda *_args, **_kwargs: "GS")
    monkeypatch.setattr(
        ann,
        "_run_sctype_scoring",
        lambda *_args, **_kwargs: pd.DataFrame({"c1": [1.0], "c2": [1.0]}, index=["T"]),
    )
    monkeypatch.setattr(
        ann,
        "_assign_sctype_celltypes",
        lambda *_args, **_kwargs: (per_cell_types, per_cell_conf),
    )

    async def _fake_cache(cache_key, results, _ctx):
        captured["cache_key"] = cache_key
        captured["results"] = results

    monkeypatch.setattr(ann, "_cache_sctype_results", _fake_cache)

    out = await ann._annotate_with_sctype(
        adata,
        params,
        DummyCtx(),
        output_key="cell_type_sctype",
        confidence_key="confidence_sctype",
    )

    assert out.cell_types == ["B", "T", "Unknown"]
    assert out.counts == {"B": 30, "T": 15, "Unknown": 15}
    assert captured["cache_key"] == "k2"
    assert captured["results"][0] == per_cell_types
    # Cache stores per-cell confidence as 4th element (not None)
    assert captured["results"][3] == per_cell_conf


@pytest.mark.asyncio
async def test_annotate_with_sctype_forwards_remote_options(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    adata = minimal_spatial_adata.copy()
    params = AnnotationParameters(
        method="sctype",
        sctype_tissue="Brain",
        sctype_use_cache=False,
        sctype_allow_remote=True,
        sctype_allow_runtime_r_install=True,
    )
    captured: dict[str, object] = {}

    def _fake_load(_ctx, *, allow_remote=False, allow_runtime_install=False):
        captured["load_allow_remote"] = allow_remote
        captured["load_allow_runtime_install"] = allow_runtime_install

    def _fake_prepare(_params, _ctx, *, allow_remote=False):
        captured["prepare_allow_remote"] = allow_remote
        return "GS"

    monkeypatch.setattr(ann, "_load_sctype_functions", _fake_load)
    monkeypatch.setattr(ann, "_prepare_sctype_genesets", _fake_prepare)
    monkeypatch.setattr(
        ann,
        "_run_sctype_scoring",
        lambda *_args, **_kwargs: pd.DataFrame({"c1": [1.0]}, index=["T"]),
    )
    monkeypatch.setattr(
        ann,
        "_assign_sctype_celltypes",
        lambda *_args, **_kwargs: (["T"] * adata.n_obs, [0.9] * adata.n_obs),
    )

    out = await ann._annotate_with_sctype(
        adata,
        params,
        DummyCtx(),
        output_key="cell_type_sctype",
        confidence_key="confidence_sctype",
    )

    assert out.cell_types == ["T"]
    assert captured == {
        "load_allow_remote": True,
        "load_allow_runtime_install": True,
        "prepare_allow_remote": True,
    }


def test_prepare_sctype_genesets_requires_tissue_without_custom_markers():
    with pytest.raises(
        ann.ParameterError,
        match="sctype_tissue is required when not using custom markers",
    ):
        ann._prepare_sctype_genesets(
            AnnotationParameters(
                method="sctype", sctype_tissue=None, sctype_custom_markers=None
            ),
            DummyCtx(),
        )


def test_convert_custom_markers_rejects_empty_dict():
    with pytest.raises(DataError, match="Custom markers dictionary is empty"):
        ann._convert_custom_markers_to_gs({}, DummyCtx())


def test_convert_custom_markers_rejects_without_positive_markers():
    bad = {"T": {"negative": ["MALAT1"]}}
    with pytest.raises(DataError, match="No valid cell types found"):
        ann._convert_custom_markers_to_gs(bad, DummyCtx())


def test_convert_custom_markers_normalizes_and_filters(monkeypatch: pytest.MonkeyPatch):
    class _Conv:
        def __add__(self, _other):
            return self

    class _Lock:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class _LCtx:
        def __enter__(self):
            return None

        def __exit__(self, exc_type, exc, tb):
            return False

    class _R:
        def __getitem__(self, name: str):
            if name == "list":
                return lambda **kwargs: kwargs
            raise KeyError(name)

    converter = _Conv()
    environment = REnvironment(
        robjects=SimpleNamespace(
            default_converter=converter,
            StrVector=lambda values: list(values),
            r=_R(),
        ),
        pandas2ri=SimpleNamespace(converter=converter),
        numpy2ri=SimpleNamespace(converter=converter),
        packages=SimpleNamespace(importr=lambda _package: object()),
        conversion=SimpleNamespace(localconverter=lambda _converter: _LCtx()),
        openrlib=SimpleNamespace(rlock=_Lock()),
        anndata2ri=None,
    )
    monkeypatch.setattr(
        ann, "validate_r_environment", lambda _ctx, **_kwargs: environment
    )

    markers = {
        "T": {"positive": [" cd3d ", "", None, "CD3E"], "negative": [" malat1 "]},
        "B": {"positive": ["MS4A1"], "negative": []},
        "ignored": ["not-a-dict"],
    }

    out = ann._convert_custom_markers_to_gs(markers, DummyCtx())

    assert set(out.keys()) == {"gs_positive", "gs_negative"}
    assert out["gs_positive"]["T"] == ["CD3D", "CD3E"]
    assert out["gs_negative"]["T"] == ["MALAT1"]
    assert out["gs_positive"]["B"] == ["MS4A1"]


def test_load_sctype_functions_runs_install_and_load_scripts(
    monkeypatch: pytest.MonkeyPatch,
):
    calls: list[str] = []
    monkeypatch.setenv("CHATSPATIAL_ALLOW_RUNTIME_R_INSTALL", "1")
    monkeypatch.setenv("CHATSPATIAL_ALLOW_REMOTE_R_SOURCE", "1")

    class _R:
        def __call__(self, script: str):
            calls.append(script)
            return None

    _mock_validate_r_environment(monkeypatch, _R())

    ann._load_sctype_functions(DummyCtx())

    assert any("required_packages" in script for script in calls)
    assert any("gene_sets_prepare.R" in script for script in calls)


def test_load_sctype_functions_rejects_remote_by_default(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.delenv("CHATSPATIAL_ALLOW_REMOTE_R_SOURCE", raising=False)
    monkeypatch.delenv("CHATSPATIAL_SCTYPE_R_DIR", raising=False)

    with pytest.raises(ann.ParameterError, match="sctype_allow_remote=true"):
        ann._load_sctype_functions(DummyCtx())


def test_load_sctype_functions_allows_remote_per_call(
    monkeypatch: pytest.MonkeyPatch,
):
    calls: list[str] = []
    monkeypatch.delenv("CHATSPATIAL_ALLOW_REMOTE_R_SOURCE", raising=False)
    monkeypatch.delenv("CHATSPATIAL_ALLOW_RUNTIME_R_INSTALL", raising=False)
    monkeypatch.delenv("CHATSPATIAL_SCTYPE_R_DIR", raising=False)

    class _R:
        def __call__(self, script: str):
            calls.append(script)
            return None

    _mock_validate_r_environment(monkeypatch, _R())

    ann._load_sctype_functions(DummyCtx(), allow_remote=True)

    assert any("gene_sets_prepare.R" in script for script in calls)
    assert not any("install.packages" in script for script in calls)


def test_prepare_sctype_genesets_uses_custom_markers_short_circuit(
    monkeypatch: pytest.MonkeyPatch,
):
    sentinel = {"gs_positive": {"T": ["CD3D"]}, "gs_negative": {"T": []}}
    monkeypatch.setattr(
        ann, "_convert_custom_markers_to_gs", lambda *_a, **_k: sentinel
    )

    out = ann._prepare_sctype_genesets(
        AnnotationParameters(
            method="sctype",
            sctype_custom_markers={"T": {"positive": ["CD3D"], "negative": []}},
        ),
        DummyCtx(),
    )
    assert out is sentinel


def test_prepare_sctype_genesets_loads_database_and_returns_gs_list(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    assigns: dict[str, object] = {}
    executed: list[str] = []

    class _R:
        def assign(self, key: str, value):
            assigns[key] = value

        def __call__(self, code: str):
            executed.append(code)
            return None

        def __getitem__(self, name: str):
            if name == "gs_list":
                return {"ok": True}
            raise KeyError(name)

    _mock_validate_r_environment(monkeypatch, _R(), converter=_FakeConverter())

    db_path = tmp_path / "db.xlsx"
    db_path.write_bytes(b"dummy")

    out = ann._prepare_sctype_genesets(
        AnnotationParameters(
            method="sctype",
            sctype_tissue="Brain",
            sctype_db_=str(db_path),
        ),
        DummyCtx(),
    )

    assert out == {"ok": True}
    assert assigns["db_path"] == db_path.as_posix()
    assert assigns["tissue_type"] == "Brain"
    assert any("gene_sets_prepare" in code for code in executed)


def test_prepare_sctype_genesets_rejects_remote_db_by_default():
    with pytest.raises(ann.ParameterError, match="sctype_allow_remote=true"):
        ann._prepare_sctype_genesets(
            AnnotationParameters(method="sctype", sctype_tissue="Brain"),
            DummyCtx(),
        )


def test_prepare_sctype_genesets_allows_remote_db_per_call(
    monkeypatch: pytest.MonkeyPatch,
):
    assigns: dict[str, object] = {}

    class _R:
        def assign(self, key: str, value):
            assigns[key] = value

        def __call__(self, _code: str):
            return None

        def __getitem__(self, name: str):
            if name == "gs_list":
                return {"ok": True}
            raise KeyError(name)

    _mock_validate_r_environment(monkeypatch, _R(), converter=_FakeConverter())
    monkeypatch.delenv("CHATSPATIAL_ALLOW_REMOTE_SCTYPE_DB", raising=False)

    out = ann._prepare_sctype_genesets(
        AnnotationParameters(method="sctype", sctype_tissue="Brain"),
        DummyCtx(),
        allow_remote=True,
    )

    assert out == {"ok": True}
    assert assigns["db_path"] == ann._SCTYPE_DEFAULT_DB_URL


def test_run_sctype_scoring_converts_r_matrix_to_dataframe(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    params = AnnotationParameters(
        method="sctype", sctype_tissue="Brain", sctype_scaled=False
    )
    assigned: dict[str, object] = {}

    class _Conv:
        def __add__(self, _other):
            return self

    class _Lock:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class _LCtx:
        def __enter__(self):
            return None

        def __exit__(self, exc_type, exc, tb):
            return False

    class _R:
        def assign(self, key: str, value):
            assigned[key] = value

        def __call__(self, code: str):
            if code == "rownames(es_max)":
                return ["T", "B"]
            if code == "colnames(es_max)":
                return ["cell_1", "cell_2"]
            return None

        def __getitem__(self, name: str):
            if name == "es_max":
                return np.array([[1.0, 0.2], [0.1, 0.9]])
            raise KeyError(name)

    converter = _Conv()
    environment = REnvironment(
        robjects=SimpleNamespace(r=_R(), default_converter=converter),
        pandas2ri=SimpleNamespace(converter=converter),
        numpy2ri=SimpleNamespace(converter=converter),
        packages=SimpleNamespace(importr=lambda _package: object()),
        conversion=SimpleNamespace(localconverter=lambda _converter: _LCtx()),
        openrlib=SimpleNamespace(rlock=_Lock()),
        anndata2ri=SimpleNamespace(converter=converter),
    )
    monkeypatch.setattr(
        ann,
        "validate_r_environment",
        lambda _ctx, **_kwargs: environment,
    )

    out = ann._run_sctype_scoring(
        adata, gs_list={"ok": True}, params=params, ctx=DummyCtx()
    )

    assert list(out.index) == ["T", "B"]
    assert list(out.columns) == ["cell_1", "cell_2"]
    assert assigned["gs_list"] == {"ok": True}


def test_run_sctype_scoring_preserves_dataframe_and_relabels_axes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    params = AnnotationParameters(
        method="sctype", sctype_tissue="Brain", sctype_scaled=False
    )

    class _Conv:
        def __add__(self, _other):
            return self

    class _Lock:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class _LCtx:
        def __enter__(self):
            return None

        def __exit__(self, exc_type, exc, tb):
            return False

    class _R:
        def assign(self, key: str, value):
            del key, value

        def __call__(self, code: str):
            if code == "rownames(es_max)":
                return ["TypeA", "TypeB"]
            if code == "colnames(es_max)":
                return ["cell_a", "cell_b"]
            return None

        def __getitem__(self, name: str):
            if name == "es_max":
                return pd.DataFrame(
                    [[3.0, 1.0], [0.5, 2.0]],
                    index=["x", "y"],
                    columns=["u", "v"],
                )
            raise KeyError(name)

    converter = _Conv()
    environment = REnvironment(
        robjects=SimpleNamespace(r=_R(), default_converter=converter),
        pandas2ri=SimpleNamespace(converter=converter),
        numpy2ri=SimpleNamespace(converter=converter),
        packages=SimpleNamespace(importr=lambda _package: object()),
        conversion=SimpleNamespace(localconverter=lambda _converter: _LCtx()),
        openrlib=SimpleNamespace(rlock=_Lock()),
        anndata2ri=SimpleNamespace(converter=converter),
    )
    monkeypatch.setattr(
        ann,
        "validate_r_environment",
        lambda _ctx, **_kwargs: environment,
    )

    out = ann._run_sctype_scoring(
        adata, gs_list={"ok": True}, params=params, ctx=DummyCtx()
    )

    assert list(out.index) == ["TypeA", "TypeB"]
    assert list(out.columns) == ["cell_a", "cell_b"]


def test_load_cached_sctype_results_returns_memory_cache_hit(
    monkeypatch: pytest.MonkeyPatch,
):
    expected = (["T"], {"T": 1}, {"T": 0.9}, None)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE", {"hit": expected})
    out = ann._load_cached_sctype_results("hit", DummyCtx())
    assert out == expected


def test_load_cached_sctype_results_returns_none_on_corrupted_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_DIR", tmp_path)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE", {})
    (tmp_path / "bad.json").write_text("{not valid json", encoding="utf-8")

    out = ann._load_cached_sctype_results("bad", DummyCtx())
    assert out is None


def test_load_cached_sctype_results_expires_stale_memory_entry(
    monkeypatch: pytest.MonkeyPatch,
):
    expected = (["T"], {"T": 1}, {"T": 0.9}, None)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE", {"hit": (0.0, expected)})
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_TTL_SECONDS", 1)
    monkeypatch.setattr(ann.time, "time", lambda: 100.0)

    out = ann._load_cached_sctype_results("hit", DummyCtx())
    assert out is None
    assert "hit" not in ann._SCTYPE_CACHE


def test_load_cached_sctype_results_expires_stale_disk_entry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_DIR", tmp_path)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE", {})
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_TTL_SECONDS", 1)
    monkeypatch.setattr(ann.time, "time", lambda: 100.0)

    cache_key = "stale"
    payload = {
        "cached_at": 0.0,
        "per_cell_types": ["T"],
        "counts": {"T": 1},
        "confidence_by_celltype": {"T": 0.9},
        "mapping_score": None,
    }
    cache_file = tmp_path / f"{cache_key}.json"
    cache_file.write_text(json.dumps(payload), encoding="utf-8")

    out = ann._load_cached_sctype_results(cache_key, DummyCtx())
    assert out is None
    assert not cache_file.exists()


@pytest.mark.asyncio
async def test_cache_sctype_results_failure_is_non_fatal(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    blocked = tmp_path / "not_a_dir"
    blocked.write_text("file", encoding="utf-8")
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_DIR", blocked)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE", {})

    await ann._cache_sctype_results(
        "k", (["T"], {"T": 1}, {"T": 0.9}, None), DummyCtx()
    )
    assert ann._SCTYPE_CACHE == {}


@pytest.mark.asyncio
async def test_cache_sctype_results_respects_lru_capacity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_DIR", tmp_path)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE", OrderedDict())
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_MAX_ITEMS", 1)
    monkeypatch.setattr(ann, "_SCTYPE_CACHE_TTL_SECONDS", 0)

    await ann._cache_sctype_results(
        "k1", (["T"], {"T": 1}, {"T": 0.9}, None), DummyCtx()
    )
    await ann._cache_sctype_results(
        "k2", (["B"], {"B": 1}, {"B": 0.8}, None), DummyCtx()
    )

    assert "k1" not in ann._SCTYPE_CACHE
    assert "k2" in ann._SCTYPE_CACHE


@pytest.mark.asyncio
async def test_annotate_with_sctype_requires_tissue_or_custom_markers(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(ann, "validate_r_environment", lambda *_a, **_k: object())

    with pytest.raises(
        ann.ParameterError, match="Either sctype_tissue or sctype_custom_markers"
    ):
        await ann._annotate_with_sctype(
            minimal_spatial_adata.copy(),
            AnnotationParameters(
                method="sctype", sctype_tissue=None, sctype_custom_markers=None
            ),
            DummyCtx(),
            "cell_type_sctype",
            "confidence_sctype",
        )


@pytest.mark.asyncio
async def test_annotate_with_sctype_rejects_invalid_tissue(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(ann, "validate_r_environment", lambda *_a, **_k: object())

    with pytest.raises(ann.ParameterError, match="not supported"):
        await ann._annotate_with_sctype(
            minimal_spatial_adata.copy(),
            AnnotationParameters(method="sctype", sctype_tissue="InvalidTissue"),
            DummyCtx(),
            "cell_type_sctype",
            "confidence_sctype",
        )


# =============================================================================
# Issue 1 regression: legacy cache with scalar mapping_score
# =============================================================================


@pytest.mark.asyncio
async def test_annotate_with_sctype_cache_hit_scalar_mapping_score_no_crash(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    """Legacy cache entries store mapping_score as a scalar float, not a
    per-cell list. The code must not call len() on a float."""
    adata = minimal_spatial_adata.copy()
    params = AnnotationParameters(
        method="sctype", sctype_tissue="Brain", sctype_use_cache=True
    )

    cached_cell_types = ["T"] * (adata.n_obs // 2) + ["B"] * (
        adata.n_obs - adata.n_obs // 2
    )
    cached_counts = {
        "T": adata.n_obs // 2,
        "B": adata.n_obs - adata.n_obs // 2,
    }

    # Legacy cache: 4th element is a scalar float (old mapping_score), not a list
    monkeypatch.setattr(ann, "_get_sctype_cache_key", lambda *_a, **_k: "legacy")
    monkeypatch.setattr(
        ann,
        "_load_cached_sctype_results",
        lambda *_a, **_k: (
            cached_cell_types,
            cached_counts,
            {"T": 0.8, "B": 0.7},
            0.95,  # scalar float — legacy format
        ),
    )

    # Pipeline should NOT run
    monkeypatch.setattr(
        ann,
        "_load_sctype_functions",
        lambda *_a, **_k: (_ for _ in ()).throw(
            AssertionError("pipeline should not run")
        ),
    )

    out = await ann._annotate_with_sctype(
        adata, params, DummyCtx(), "cell_type_sctype", "confidence_sctype"
    )

    assert out.cell_types == ["T", "B"]
    assert "confidence_sctype" in adata.obs
    # Confidence should fall back to per-type mean, not crash
    assert adata.obs["confidence_sctype"].iloc[0] == 0.8  # T-cell confidence
