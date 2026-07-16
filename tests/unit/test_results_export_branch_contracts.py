"""Branch-focused tests for results export extraction and index behavior."""

from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from chatspatial.utils import results_export as re


def _patch_home(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(re.Path, "home", classmethod(lambda _cls: tmp_path))


def test_export_analysis_result_returns_empty_when_results_keys_missing(
    minimal_spatial_adata, monkeypatch, tmp_path: Path
) -> None:
    _patch_home(monkeypatch, tmp_path)
    adata = minimal_spatial_adata.copy()
    adata.uns["demo_metadata"] = {"method": "demo", "results_keys": {}}
    assert re.export_analysis_result(adata, "d_meta", "demo") == []


def test_get_export_limit_handles_invalid_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CHATSPATIAL_EXPORT_MAX_ROWS", "not-an-int")
    assert re._get_export_limit("CHATSPATIAL_EXPORT_MAX_ROWS", 123) == 123


def test_export_analysis_result_can_be_disabled_by_env(
    minimal_spatial_adata, monkeypatch, tmp_path: Path
) -> None:
    _patch_home(monkeypatch, tmp_path)
    monkeypatch.setenv("CHATSPATIAL_EXPORT_RESULTS", "0")

    adata = minimal_spatial_adata.copy()
    adata.uns["demo_metadata"] = {
        "method": "demo",
        "results_keys": {"obs": ["group"]},
    }

    assert re.export_analysis_result(adata, "d_disabled", "demo") == []
    assert not (tmp_path / ".chatspatial" / "results" / "d_disabled").exists()


def test_warn_large_export_does_not_truncate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reproducibility exports must be lossless — warn but never truncate."""
    monkeypatch.setenv("CHATSPATIAL_EXPORT_WARN_ROWS", "3")
    monkeypatch.setenv("CHATSPATIAL_EXPORT_WARN_COLS", "2")
    df = pd.DataFrame(
        np.arange(20, dtype=float).reshape(5, 4),
        columns=["c0", "c1", "c2", "c3"],
    )

    # Should not raise or modify df
    re._warn_large_export(df, location="obs", key="demo")
    assert df.shape == (5, 4)  # unchanged


def test_export_analysis_result_unknown_location_is_skipped(
    minimal_spatial_adata, monkeypatch, tmp_path: Path
) -> None:
    _patch_home(monkeypatch, tmp_path)
    adata = minimal_spatial_adata.copy()
    adata.uns["demo_metadata"] = {
        "method": "demo",
        "results_keys": {"unknown_loc": ["abc"]},
    }
    assert re.export_analysis_result(adata, "d_loc", "demo") == []


def test_extract_rank_genes_groups_missing_or_failing_returns_none(
    minimal_spatial_adata, monkeypatch
) -> None:
    adata = minimal_spatial_adata.copy()
    fake_scanpy = SimpleNamespace(
        get=SimpleNamespace(rank_genes_groups_df=lambda *_a, **_k: pd.DataFrame())
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    assert re._extract_rank_genes_groups(adata) is None

    adata.uns["rank_genes_groups"] = {"names": np.array([])}
    fake_scanpy_bad = SimpleNamespace(
        get=SimpleNamespace(
            rank_genes_groups_df=lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("boom"))
        )
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy_bad)
    assert re._extract_rank_genes_groups(adata) is None


def test_extract_from_uns_covers_dataframe_ripley_dict_structured_and_unsupported(
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    assert re._extract_from_uns(adata, "missing_key") is None

    adata.uns["as_df"] = pd.DataFrame({"v": [1.0, 2.0]})
    out_df = re._extract_from_uns(adata, "as_df")
    assert isinstance(out_df, pd.DataFrame)
    assert list(out_df.columns) == ["v"]

    adata.uns["my_ripley_domain"] = {"L_stat": pd.DataFrame({"L": [0.1, 0.2]})}
    out_ripley = re._extract_from_uns(adata, "my_ripley_domain")
    assert isinstance(out_ripley, pd.DataFrame)
    assert list(out_ripley.columns) == ["L"]

    adata.uns["flat_dict"] = {"a": 1.0, "b": 2.0}
    out_flat = re._extract_from_uns(adata, "flat_dict")
    assert "value" in out_flat.columns

    rec = np.array([(1, 2.0)], dtype=[("x", "i4"), ("y", "f8")])
    adata.uns["recarray"] = rec
    out_rec = re._extract_from_uns(adata, "recarray")
    assert list(out_rec.columns) == ["x", "y"]

    adata.uns["unsupported"] = 123
    assert re._extract_from_uns(adata, "unsupported") is None


def test_extract_from_uns_routes_squidpy_keys_to_spatial_extractor(
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    adata.obs["domain"] = pd.Categorical([f"d{i % 3}" for i in range(adata.n_obs)])
    adata.uns["domain_nhood_enrichment"] = {"zscore": np.ones((3, 3))}

    out = re._extract_from_uns(adata, "domain_nhood_enrichment")
    assert isinstance(out, pd.DataFrame)
    assert any(c.startswith("zscore_") for c in out.columns)


def test_extract_squidpy_result_handles_non_dict_missing_cluster_and_empty_metrics(
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    assert re._extract_squidpy_spatial_result(adata, "domain_co_occurrence", ["bad"]) is None
    assert re._extract_squidpy_spatial_result(adata, "nope_nhood_enrichment", {"zscore": np.eye(2)}) is None

    adata.obs["domain"] = [f"d{i % 3}" for i in range(adata.n_obs)]  # non-categorical branch
    out_empty = re._extract_squidpy_spatial_result(
        adata,
        "domain_nhood_enrichment",
        {"not_matrix": np.array([1, 2, 3])},
    )
    assert out_empty is None

    out_2d = re._extract_squidpy_spatial_result(
        adata,
        "domain_nhood_enrichment",
        {"zscore": np.ones((3, 3)), "count": np.ones((3, 3))},
    )
    assert isinstance(out_2d, pd.DataFrame)
    assert any(c.startswith("zscore_") for c in out_2d.columns)
    assert any(c.startswith("count_") for c in out_2d.columns)


def test_extract_squidpy_result_rejects_unrecognized_key_suffix(
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    adata.obs["domain"] = pd.Categorical([f"d{i % 2}" for i in range(adata.n_obs)])
    assert re._extract_squidpy_spatial_result(adata, "domain", {"zscore": np.eye(2)}) is None


def test_extract_from_obs_and_var_return_none_for_non_matching_key(minimal_spatial_adata) -> None:
    adata = minimal_spatial_adata.copy()
    assert re._extract_from_obs(adata, "not_present") is None
    assert re._extract_from_var(adata, "not_present") is None


def test_extract_from_obs_and_var_support_prefix_wildcard(minimal_spatial_adata) -> None:
    adata = minimal_spatial_adata.copy()
    adata.obs["score_alpha"] = np.arange(adata.n_obs)
    adata.obs["score_beta"] = np.arange(adata.n_obs)
    adata.var["gene_metric_a"] = np.arange(adata.n_vars)
    adata.var["gene_metric_b"] = np.arange(adata.n_vars)

    out_obs = re._extract_from_obs(adata, "score_*")
    out_var = re._extract_from_var(adata, "gene_metric_*")

    assert out_obs is not None and set(out_obs.columns) == {"score_alpha", "score_beta"}
    assert out_var is not None and set(out_var.columns) == {
        "gene_metric_a",
        "gene_metric_b",
    }

    # No implicit substring matching without wildcard.
    assert re._extract_from_obs(adata, "score_") is None
    assert re._extract_from_var(adata, "gene_metric_") is None


def test_extract_from_obs_and_var_exact_match_fallback_for_non_substring_keys() -> None:
    class _Key:
        def __contains__(self, _item) -> bool:
            return False

    obs_key = _Key()
    var_key = _Key()
    adata = SimpleNamespace(
        obs=pd.DataFrame({obs_key: [1, 2]}, index=["c0", "c1"]),
        var=pd.DataFrame({var_key: [0.1, 0.2]}, index=["g0", "g1"]),
    )

    out_obs = re._extract_from_obs(adata, obs_key)
    out_var = re._extract_from_var(adata, var_key)

    assert out_obs is not None and list(out_obs.columns) == [obs_key]
    assert out_var is not None and list(out_var.columns) == [var_key]


def test_extract_from_obsm_covers_missing_dataframe_array_and_unknown_object(
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    assert re._extract_from_obsm(adata, "not_present") is None

    adata.obsm["df_like"] = pd.DataFrame(
        {"a": np.arange(adata.n_obs)},
        index=adata.obs_names,
    )
    out_df = re._extract_from_obsm(adata, "df_like")
    assert isinstance(out_df, pd.DataFrame)
    assert list(out_df.index) == list(adata.obs_names)

    adata.obsm["array_like"] = np.ones((adata.n_obs, 2))
    out_arr = re._extract_from_obsm(adata, "array_like")
    assert list(out_arr.columns) == ["array_like_0", "array_like_1"]

    class _FakeAdata:
        def __init__(self):
            self.obsm = {"unknown": object()}
            self.obs_names = ["c0", "c1"]
            self.uns = {}

    assert re._extract_from_obsm(_FakeAdata(), "unknown") is None


def test_extract_as_dataframe_routes_obsm_location(minimal_spatial_adata) -> None:
    adata = minimal_spatial_adata.copy()
    adata.obsm["latent"] = np.ones((adata.n_obs, 2), dtype=float)
    out = re._extract_as_dataframe(adata, "obsm", "latent", "demo")
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == ["latent_0", "latent_1"]


def test_infer_obsm_columns_prefers_metadata_and_falls_back(minimal_spatial_adata) -> None:
    adata = minimal_spatial_adata.copy()

    # Deconvolution: uns["deconvolution_{method}_cell_types"]
    adata.uns["deconvolution_spotlight_cell_types"] = ["T", "B"]
    cols_deconv = re._infer_obsm_columns(adata, "deconvolution_spotlight", 2)
    assert cols_deconv == ["T", "B"]

    # CCC spatial (shared key): uns["ccc"]["lr_pairs"]
    adata.uns["ccc"] = {"lr_pairs": ["A^B", "B^C"]}
    cols_scores = re._infer_obsm_columns(adata, "ccc_spatial_scores", 2)
    assert cols_scores == ["A^B", "B^C"]
    cols_pvals = re._infer_obsm_columns(adata, "ccc_spatial_pvals", 2)
    assert cols_pvals == ["A^B", "B^C"]

    # CCC spatial (per-method key): uns["ccc_liana"]["lr_pairs"] takes priority
    adata.uns["ccc_liana"] = {"lr_pairs": ["X^Y", "Y^Z"]}
    cols_method = re._infer_obsm_columns(
        adata, "ccc_spatial_scores_liana", 2
    )
    assert cols_method == ["X^Y", "Y^Z"]
    cols_method_pvals = re._infer_obsm_columns(
        adata, "ccc_spatial_pvals_liana", 2
    )
    assert cols_method_pvals == ["X^Y", "Y^Z"]

    # CCC spatial (per-method key missing → falls back to shared)
    cols_fallback = re._infer_obsm_columns(
        adata, "ccc_spatial_scores_cellphonedb", 2
    )
    assert cols_fallback == ["A^B", "B^C"]

    # Fallback: numeric indices
    cols_default = re._infer_obsm_columns(adata, "x_embed", 3)
    assert cols_default == ["x_embed_0", "x_embed_1", "x_embed_2"]


def test_dict_to_dataframe_covers_empty_nested_list_and_flat() -> None:
    out_empty = re._dict_to_dataframe({})
    assert out_empty.empty

    out_nested = re._dict_to_dataframe({"p1": {"score": 1.0}, "p2": {"score": 2.0}})
    assert "score" in out_nested.columns

    out_list = re._dict_to_dataframe({"g1": [1, 2], "g2": [3, 4]})
    assert out_list.shape == (2, 2)

    out_flat = re._dict_to_dataframe({"k1": 1.0, "k2": 2.0})
    assert list(out_flat.columns) == ["value"]


def test_update_index_recovers_from_invalid_json_and_list_exported_results_missing_index(
    monkeypatch, tmp_path: Path
) -> None:
    _patch_home(monkeypatch, tmp_path)
    assert re.list_exported_results("no_idx") == {}

    results_dir = re.get_results_dir("d_idx")
    index_path = results_dir / "_index.json"
    index_path.write_text("{invalid json", encoding="utf-8")

    export_file = re.get_analysis_dir("d_idx", "demo") / "demo_key.csv"
    export_file.parent.mkdir(parents=True, exist_ok=True)
    export_file.write_text("a,b\n1,2\n", encoding="utf-8")

    re._update_index(
        data_id="d_idx",
        analysis_name="demo",
        method="demo_method",
        metadata={"parameters": {"path": Path("x/y")}, "statistics": {"n": 1}},
        exported_files=[export_file],
    )

    loaded = json.loads(index_path.read_text(encoding="utf-8"))
    assert loaded["analyses"]["demo"]["method"] == "demo_method"
    assert loaded["analyses"]["demo"]["files"] == ["demo_key.csv"]
    assert loaded["analyses"]["demo"]["parameters"]["path"] == "x/y"


def test_update_index_failure_preserves_previous_complete_index(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_home(monkeypatch, tmp_path)
    index_path = re.get_results_dir("dataset") / "_index.json"
    previous = {
        "data_id": "dataset",
        "analyses": {"existing": {"method": "kept", "files": []}},
    }
    index_path.write_text(json.dumps(previous), encoding="utf-8")

    def _partial_dump(_value, file, **_kwargs) -> None:
        file.write('{"partial"')
        raise OSError("disk full")

    monkeypatch.setattr(re.json, "dump", _partial_dump)

    with pytest.raises(OSError, match="disk full"):
        re._update_index("dataset", "new", "method", {}, [])

    assert json.loads(index_path.read_text(encoding="utf-8")) == previous
    assert list(index_path.parent.iterdir()) == [index_path]


def test_concurrent_index_updates_preserve_every_analysis(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    _patch_home(monkeypatch, tmp_path)

    def _update(name: str) -> None:
        re._update_index("dataset", name, name, {}, [])

    names = [f"analysis_{index}" for index in range(8)]
    with ThreadPoolExecutor(max_workers=8) as executor:
        list(executor.map(_update, names))

    index_path = re.get_results_dir("dataset") / "_index.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    assert set(index["analyses"]) == set(names)


def test_result_timestamps_are_explicit_utc() -> None:
    timestamp = datetime.fromisoformat(re._now())
    assert timestamp.utcoffset() == timedelta(0)


def test_sanitize_for_json_preserves_numeric_types_and_rejects_nonfinite_values() -> None:
    sanitized = re._sanitize_for_json(
        {
            1: np.int64(3),
            "ratio": np.float32(0.5),
            "values": np.array([1.0, np.nan, np.inf]),
            "path": Path("data/result.csv"),
        }
    )

    assert sanitized == {
        "1": 3,
        "ratio": 0.5,
        "values": [1.0, None, None],
        "path": "data/result.csv",
    }


def test_extract_squidpy_matrix_filters_unused_categories(
    minimal_spatial_adata,
) -> None:
    """Squidpy matrix export uses only observed categories, not all cat levels."""
    adata = minimal_spatial_adata.copy()

    # Create a categorical with an unused level "C"
    adata.obs["cluster"] = pd.Categorical(
        ["A", "B"] * (adata.n_obs // 2) + ["A"] * (adata.n_obs % 2),
        categories=["A", "B", "C"],  # "C" is never observed
    )

    # Key format: {cluster_key}_nhood_enrichment (squidpy convention)
    # Create a 2x2 matrix matching only observed categories (A, B)
    matrix_2x2 = np.array([[1.0, 0.5], [0.5, 1.0]])
    adata.uns["cluster_nhood_enrichment"] = {"zscore": matrix_2x2}

    df = re._extract_from_uns(adata, "cluster_nhood_enrichment")
    assert df is not None
    # Should have 2 labels (A, B), not 3 (A, B, C)
    # Columns are prefixed with metric name: zscore_A, zscore_B
    assert list(df.columns) == ["zscore_A", "zscore_B"]
    assert list(df.index) == ["A", "B"]


def test_extract_rank_genes_groups_per_run_key(
    minimal_spatial_adata,
) -> None:
    """Per-run rank_genes_groups_* key should export correctly, not crash."""
    import scanpy as sc

    adata = minimal_spatial_adata.copy()
    sc.tl.rank_genes_groups(adata, groupby="group", method="wilcoxon")

    # Create per-run copy (simulates what find_markers does)
    per_run_key = "rank_genes_groups_wilcoxon_all_groups"
    adata.uns[per_run_key] = adata.uns["rank_genes_groups"]

    df = re._extract_as_dataframe(adata, "uns", per_run_key, "test")
    assert df is not None
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_extract_rank_genes_groups_per_run_key_restores_original(
    minimal_spatial_adata,
) -> None:
    """Per-run key export must restore the original rank_genes_groups entry."""
    import scanpy as sc

    adata = minimal_spatial_adata.copy()
    sc.tl.rank_genes_groups(adata, groupby="group", method="wilcoxon")
    original_rgg = adata.uns["rank_genes_groups"]

    per_run_key = "rank_genes_groups_wilcoxon_all_groups"
    adata.uns[per_run_key] = adata.uns["rank_genes_groups"]

    re._extract_as_dataframe(adata, "uns", per_run_key, "test")
    # Original key must still be intact
    assert adata.uns["rank_genes_groups"] is original_rgg


def test_extract_rank_genes_groups_passes_custom_key_without_mutating_uns(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
) -> None:
    import scanpy as sc

    adata = minimal_spatial_adata.copy()
    custom_result = {"sentinel": True}
    adata.uns["rank_genes_groups_custom"] = custom_result
    captured: dict[str, object] = {}

    def _extract(_adata, group, *, key):
        captured.update(group=group, key=key)
        return pd.DataFrame({"names": ["G1"]})

    monkeypatch.setattr(sc.get, "rank_genes_groups_df", _extract)

    result = re._extract_rank_genes_groups(adata, "rank_genes_groups_custom")

    assert result is not None
    assert captured == {"group": None, "key": "rank_genes_groups_custom"}
    assert "rank_genes_groups" not in adata.uns
    assert adata.uns["rank_genes_groups_custom"] is custom_result


def test_extract_ccc_per_method_key(
    minimal_spatial_adata,
) -> None:
    """Per-method ccc_<method> key should extract results DataFrame."""
    adata = minimal_spatial_adata.copy()

    results_df = pd.DataFrame(
        {"source": ["A", "B"], "target": ["C", "D"], "score": [0.5, 0.8]}
    )
    ccc_storage = {
        "results": results_df,
        "lr_pairs": ["A_C", "B_D"],
        "statistics": {"n_pairs": 2},
    }
    adata.uns["ccc_liana"] = ccc_storage

    df = re._extract_as_dataframe(
        adata, "uns", "ccc_liana", "test"
    )
    assert df is not None
    assert isinstance(df, pd.DataFrame)
    assert "source" in df.columns
    assert len(df) == 2


def test_export_analysis_result_updates_index_when_all_exports_fail(
    minimal_spatial_adata, monkeypatch, tmp_path: Path
) -> None:
    """Audit trail must be recorded even when every export raises an exception."""
    _patch_home(monkeypatch, tmp_path)
    monkeypatch.setenv("CHATSPATIAL_EXPORT_RESULTS", "1")
    adata = minimal_spatial_adata.copy()
    adata.uns["fail_metadata"] = {
        "method": "fail_method",
        "results_keys": {"obs": ["nonexistent_key_a", "nonexistent_key_b"]},
    }

    result = re.export_analysis_result(adata, "d_fail", "fail")
    assert result == []

    # Index should still have been written for the audit trail
    index_path = re.get_results_dir("d_fail") / "_index.json"
    assert index_path.exists()
    loaded = json.loads(index_path.read_text(encoding="utf-8"))
    assert "fail" in loaded["analyses"]
    assert loaded["analyses"]["fail"]["method"] == "fail_method"
