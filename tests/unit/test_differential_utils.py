"""Unit tests for differential expression tool contracts."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from chatspatial.models.analysis import DifferentialExpressionResult
from chatspatial.models.data import DifferentialExpressionParameters
from chatspatial.tools import differential as differential_mod
from chatspatial.tools.differential import _run_pydeseq2, differential_expression
from chatspatial.utils.exceptions import (
    DataError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


class DummyCtx:
    def __init__(self, adata: AnnData):
        self._adata = adata
        self.warnings: list[str] = []
        self.infos: list[str] = []

    async def get_adata(self, _data_id: str) -> AnnData:
        return self._adata

    async def set_adata(self, _data_id: str, adata: AnnData) -> None:
        self._adata = adata

    async def warning(self, msg: str) -> None:
        self.warnings.append(msg)

    async def info(self, msg: str) -> None:
        self.infos.append(msg)


def _make_de_adata() -> AnnData:
    X = np.array(
        [
            [20, 1, 0],
            [18, 2, 0],
            [21, 1, 0],
            [19, 3, 0],
            [17, 2, 0],
            [1, 15, 0],
            [2, 14, 0],
            [1, 16, 0],
            [2, 13, 0],
            [1, 12, 0],
            [10, 10, 0],
            [11, 9, 0],
        ],
        dtype=np.float32,
    )
    adata = AnnData(X=X)
    adata.var_names = ["gene_0", "gene_1", "gene_2"]
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.obs["cluster"] = ["A"] * 5 + ["B"] * 5 + ["tiny"] * 2
    # Preserve raw for log2FC path
    adata.raw = adata.copy()
    return adata


def _fake_rank_genes_groups_factory(names_by_group: dict[str, list[str]]):
    def _fake_rank_genes_groups(
        adata, groupby, method, n_genes, reference, groups=None
    ):
        del groupby, method, reference, groups
        fields = [(g, "U64") for g in names_by_group]
        out = np.zeros((n_genes,), dtype=fields)
        pvals = np.zeros((n_genes,), dtype=[(g, "f8") for g in names_by_group])
        for g, genes in names_by_group.items():
            padded = list(genes) + [genes[-1]] * max(0, n_genes - len(genes))
            out[g] = np.array(padded[:n_genes], dtype="U64")
            pvals[g] = np.linspace(0.001, 0.05, n_genes)
        adata.uns["rank_genes_groups"] = {"names": out, "pvals_adj": pvals}

    return _fake_rank_genes_groups


@pytest.mark.asyncio
async def test_differential_expression_dispatches_pydeseq2(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_de_adata()
    ctx = DummyCtx(adata)

    expected = DifferentialExpressionResult(
        data_id="d1",
        comparison="A vs B",
        n_genes=1,
        top_genes=["gene_0"],
        statistics={"method": "pydeseq2"},
    )

    async def _fake_run(_data_id, _ctx, _params):
        return expected

    monkeypatch.setattr(differential_mod, "_run_pydeseq2", _fake_run)

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
    )
    result = await differential_expression("d1", ctx, params)

    assert result == expected


@pytest.mark.asyncio
async def test_run_pydeseq2_requires_sample_key():
    adata = _make_de_adata()
    ctx = DummyCtx(adata)
    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key=None,
    )

    with pytest.raises(ParameterError, match="sample_key is required"):
        await _run_pydeseq2("d2", ctx, params)


@pytest.mark.asyncio
async def test_differential_all_groups_skips_tiny_groups(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_de_adata()
    ctx = DummyCtx(adata)

    monkeypatch.setattr(
        differential_mod.sc.tl,
        "rank_genes_groups",
        _fake_rank_genes_groups_factory(
            {
                "A": ["gene_0", "gene_1", "gene_2"],
                "B": ["gene_1", "gene_0", "gene_2"],
            }
        ),
    )
    monkeypatch.setattr(
        differential_mod, "store_analysis_metadata", lambda *a, **k: None
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *a, **k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1=None,
        method="wilcoxon",
        n_top_genes=3,
        min_cells=3,
    )
    result = await differential_expression("d3", ctx, params)

    assert result.comparison == "All groups in cluster"
    assert result.n_genes == 3
    assert result.top_genes == ["gene_0", "gene_1", "gene_2"]
    assert any("Skipped 1 group(s)" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_differential_late_failure_does_not_publish_partial_results(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_de_adata()
    ctx = DummyCtx(adata)
    monkeypatch.setattr(
        differential_mod.sc.tl,
        "rank_genes_groups",
        _fake_rank_genes_groups_factory({"A": ["gene_0"], "B": ["gene_1"]}),
    )
    monkeypatch.setattr(
        differential_mod,
        "export_analysis_result",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("export failed")),
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1=None,
        method="wilcoxon",
        n_top_genes=1,
        min_cells=3,
    )
    with pytest.raises(RuntimeError, match="export failed"):
        await differential_expression("d3", ctx, params)

    assert ctx._adata is adata
    assert "rank_genes_groups" not in adata.uns


@pytest.mark.asyncio
async def test_differential_specific_group_accepts_any_raw_source(
    monkeypatch: pytest.MonkeyPatch,
):
    """Fold change calculation accepts any source from get_raw_data_source (raw, counts, X)."""
    adata = _make_de_adata()
    ctx = DummyCtx(adata)

    monkeypatch.setattr(
        differential_mod.sc.tl,
        "rank_genes_groups",
        _fake_rank_genes_groups_factory({"A": ["gene_0", "gene_1"]}),
    )
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(
            source="X", X=adata.X, var_names=adata.var_names, is_integer_counts=False
        ),
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1="A",
        group2="B",
        method="wilcoxon",
        n_top_genes=2,
    )

    result = await differential_expression("d4", ctx, params)
    assert result is not None
    assert result.n_genes > 0


@pytest.mark.asyncio
async def test_differential_specific_group_warns_on_missing_genes(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_de_adata()
    ctx = DummyCtx(adata)

    monkeypatch.setattr(
        differential_mod.sc.tl,
        "rank_genes_groups",
        _fake_rank_genes_groups_factory({"A": ["gene_0", "missing_gene"]}),
    )
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(
            source="raw",
            X=adata.raw.X,
            var_names=adata.raw.var_names,
            is_integer_counts=True,
        ),
    )
    monkeypatch.setattr(
        differential_mod, "store_analysis_metadata", lambda *a, **k: None
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *a, **k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1="A",
        group2="B",
        method="wilcoxon",
        n_top_genes=2,
    )
    result = await differential_expression("d5", ctx, params)

    assert result.comparison == "A vs B"
    assert result.top_genes == ["gene_0", "missing_gene"]
    assert result.statistics["mean_log2fc"] is not None
    assert any("genes not found in raw data" in msg for msg in ctx.warnings)


def _patch_pydeseq2_modules(
    monkeypatch: pytest.MonkeyPatch,
    dds_module=None,
    ds_module=None,
):
    dds_module = dds_module or SimpleNamespace(DeseqDataSet=object)
    ds_module = ds_module or SimpleNamespace(DeseqStats=object)
    modules = {
        "pydeseq2.dds": dds_module,
        "pydeseq2.ds": ds_module,
    }
    monkeypatch.setattr(
        differential_mod,
        "require_module",
        lambda _name, module_name, *_args, **_kwargs: modules[module_name],
    )


def _install_fake_pydeseq2(
    monkeypatch: pytest.MonkeyPatch, results_df, captured: dict[str, object]
):
    from types import ModuleType

    class _DDS:
        def __init__(self, counts, metadata, design_factors):
            captured["counts_shape"] = tuple(counts.shape)
            captured["metadata_shape"] = tuple(metadata.shape)
            captured["design_factors"] = design_factors

        def deseq2(self):
            captured["deseq2_called"] = True

    class _Stats:
        def __init__(self, dds, contrast):
            del dds
            captured["contrast"] = contrast
            self.results_df = results_df.copy()

        def summary(self):
            captured["summary_called"] = True

    m_pkg = ModuleType("pydeseq2")
    m_dds = ModuleType("pydeseq2.dds")
    m_ds = ModuleType("pydeseq2.ds")
    m_dds.DeseqDataSet = _DDS
    m_ds.DeseqStats = _Stats

    monkeypatch.setitem(__import__("sys").modules, "pydeseq2", m_pkg)
    monkeypatch.setitem(__import__("sys").modules, "pydeseq2.dds", m_dds)
    monkeypatch.setitem(__import__("sys").modules, "pydeseq2.ds", m_ds)
    _patch_pydeseq2_modules(monkeypatch, m_dds, m_ds)


def _make_pydeseq2_adata() -> AnnData:
    # 8 cells, 3 genes, two conditions across two biological samples
    X = np.array(
        [
            [20, 2, 1],
            [18, 1, 1],
            [2, 21, 1],
            [1, 19, 1],
            [22, 3, 1],
            [17, 2, 1],
            [3, 20, 1],
            [2, 18, 1],
        ],
        dtype=np.int64,
    )
    adata = AnnData(X=X)
    adata.var_names = ["gene_0", "gene_1", "gene_2"]
    adata.obs_names = [f"c{i}" for i in range(adata.n_obs)]
    adata.obs["cluster"] = ["A", "A", "B", "B", "A", "A", "B", "B"]
    adata.obs["sample"] = ["s1", "s1", "s1", "s1", "s2", "s2", "s2", "s2"]
    return adata


@pytest.mark.asyncio
async def test_run_pydeseq2_preserves_managed_submodule_error(
    monkeypatch: pytest.MonkeyPatch,
):
    def _missing_submodule(*_args, **_kwargs):
        raise DependencyError("pydeseq2.dds is unavailable")

    monkeypatch.setattr(differential_mod, "require_module", _missing_submodule)

    with pytest.raises(DependencyError, match="pydeseq2.dds is unavailable"):
        await _run_pydeseq2(
            "d_missing",
            DummyCtx(_make_pydeseq2_adata()),
            DifferentialExpressionParameters(
                group_key="cluster",
                method="pydeseq2",
                sample_key="sample",
            ),
        )


@pytest.mark.asyncio
async def test_run_pydeseq2_success_auto_group_selection_and_persistence(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_pydeseq2_adata()
    ctx = DummyCtx(adata)
    captured: dict[str, object] = {}

    results_df = pd.DataFrame(
        {
            "padj": [0.01, 0.2],
            "log2FoldChange": [1.3, -0.4],
        },
        index=["gene_0", "gene_1"],
    )

    _install_fake_pydeseq2(monkeypatch, results_df, captured)

    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(X=adata.X, var_names=adata.var_names),
    )
    monkeypatch.setattr(
        differential_mod, "check_is_integer_counts", lambda _x: (True, None, None)
    )
    monkeypatch.setattr(
        differential_mod,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *_a, **_k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
        group1=None,
        group2=None,
        n_top_genes=5,
    )

    out = await _run_pydeseq2("d6", ctx, params)

    assert out.comparison == "A vs B"
    assert out.n_genes == 2
    assert out.top_genes == ["gene_0", "gene_1"]
    assert out.statistics["n_pseudobulk_samples"] == 4
    assert any("No group specified" in msg for msg in ctx.infos)
    assert captured["analysis_name"] == "de_pydeseq2_A_vs_B"
    assert captured["method"] == "pydeseq2"
    assert captured["statistics"]["n_pseudobulk_samples"] == 4
    assert captured["contrast"] == ["condition", "A", "B"]
    assert "pydeseq2_results_A_vs_B" not in adata.uns
    assert "pydeseq2_results_A_vs_B" in ctx._adata.uns
    assert "_de_condition" not in adata.obs.columns
    assert "_pseudobulk_id" not in adata.obs.columns


@pytest.mark.asyncio
async def test_run_pydeseq2_warns_for_non_integer_counts(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_pydeseq2_adata()
    ctx = DummyCtx(adata)
    captured: dict[str, object] = {}

    results_df = pd.DataFrame(
        {"padj": [0.04], "log2FoldChange": [0.8]}, index=["gene_0"]
    )
    _install_fake_pydeseq2(monkeypatch, results_df, captured)

    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(
            X=adata.X.astype(float), var_names=adata.var_names
        ),
    )
    monkeypatch.setattr(
        differential_mod, "check_is_integer_counts", lambda _x: (False, None, None)
    )
    monkeypatch.setattr(
        differential_mod,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *_a, **_k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
        group1="A",
        group2="B",
    )

    out = await _run_pydeseq2("d7", ctx, params)

    assert out.comparison == "A vs B"
    assert any("requires raw integer counts" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_run_pydeseq2_rejects_insufficient_total_pseudobulk_samples(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_pydeseq2_adata()
    adata.obs["sample"] = "s1"  # only two pseudobulk samples: s1_A, s1_B
    ctx = DummyCtx(adata)

    _patch_pydeseq2_modules(monkeypatch)
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(X=adata.X, var_names=adata.var_names),
    )
    monkeypatch.setattr(
        differential_mod, "check_is_integer_counts", lambda _x: (True, None, None)
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
        group1="A",
        group2="B",
    )

    with pytest.raises(DataError, match="Found only 2 total pseudobulk samples"):
        await _run_pydeseq2("d8", ctx, params)


@pytest.mark.asyncio
async def test_run_pydeseq2_rejects_condition_with_single_sample(
    monkeypatch: pytest.MonkeyPatch,
):
    # 4 pseudobulk groups total, but B has only one biological sample
    X = np.array([[10, 1, 0], [1, 10, 0], [11, 1, 0], [12, 1, 0]], dtype=np.int64)
    adata = AnnData(X=X)
    adata.var_names = ["gene_0", "gene_1", "gene_2"]
    adata.obs_names = ["c1", "c2", "c3", "c4"]
    adata.obs["cluster"] = ["A", "B", "A", "A"]
    adata.obs["sample"] = ["s1", "s1", "s2", "s3"]
    ctx = DummyCtx(adata)

    _patch_pydeseq2_modules(monkeypatch)
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(X=adata.X, var_names=adata.var_names),
    )
    monkeypatch.setattr(
        differential_mod, "check_is_integer_counts", lambda _x: (True, None, None)
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
        group1="A",
        group2="B",
    )

    with pytest.raises(DataError, match="at least 2 samples per condition"):
        await _run_pydeseq2("d9", ctx, params)


@pytest.mark.asyncio
async def test_differential_all_groups_raises_when_all_groups_below_min_cells():
    adata = _make_de_adata()
    ctx = DummyCtx(adata)

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1=None,
        method="wilcoxon",
        min_cells=10,
    )

    with pytest.raises(DataError, match="All groups have <10 cells"):
        await differential_expression("d10", ctx, params)


@pytest.mark.asyncio
async def test_differential_all_groups_converts_float16_after_filtering(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_de_adata()
    adata.X = adata.X.astype(np.float16)
    adata.raw = adata.copy()
    ctx = DummyCtx(adata)
    captured: dict[str, np.dtype] = {}

    def _fake_rank(adata, groupby, method, n_genes, reference, groups=None):
        del groupby, method, reference, groups
        captured["dtype"] = adata.X.dtype
        names = np.zeros((n_genes,), dtype=[("A", "U64"), ("B", "U64")])
        pvals = np.zeros((n_genes,), dtype=[("A", "f8"), ("B", "f8")])
        names["A"] = np.array(["gene_0"] * n_genes, dtype="U64")
        names["B"] = np.array(["gene_1"] * n_genes, dtype="U64")
        pvals["A"] = np.linspace(0.001, 0.01, n_genes)
        pvals["B"] = np.linspace(0.001, 0.01, n_genes)
        adata.uns["rank_genes_groups"] = {"names": names, "pvals_adj": pvals}

    monkeypatch.setattr(differential_mod.sc.tl, "rank_genes_groups", _fake_rank)
    monkeypatch.setattr(
        differential_mod, "store_analysis_metadata", lambda *a, **k: None
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *a, **k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1=None,
        method="wilcoxon",
        n_top_genes=2,
        min_cells=3,
    )
    await differential_expression("d11", ctx, params)

    assert captured["dtype"] == np.float32


@pytest.mark.asyncio
async def test_differential_specific_group_rejects_unknown_group2():
    adata = _make_de_adata()
    ctx = DummyCtx(adata)

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1="A",
        group2="missing",
        method="wilcoxon",
    )

    with pytest.raises(ParameterError, match="Group 'missing' not found"):
        await differential_expression("d12", ctx, params)


@pytest.mark.asyncio
async def test_differential_specific_group_uses_fallback_name_column_and_float16_cast(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_de_adata()
    adata.X = adata.X.astype(np.float16)
    adata.raw = adata.copy()
    ctx = DummyCtx(adata)
    captured: dict[str, np.dtype] = {}

    def _fake_rank(adata, groupby, groups, reference, method, n_genes):
        del groupby, groups, reference, method
        captured["dtype"] = adata.X.dtype
        names = np.zeros((n_genes,), dtype=[("fallback", "U64")])
        pvals = np.zeros((n_genes,), dtype=[("fallback", "f8")])
        names["fallback"] = np.array(["gene_0"] * n_genes, dtype="U64")
        pvals["fallback"] = np.linspace(0.01, 0.02, n_genes)
        adata.uns["rank_genes_groups"] = {"names": names, "pvals_adj": pvals}

    monkeypatch.setattr(differential_mod.sc.tl, "rank_genes_groups", _fake_rank)
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(
            source="raw",
            X=adata.raw.X,
            var_names=adata.raw.var_names,
            is_integer_counts=True,
        ),
    )
    monkeypatch.setattr(
        differential_mod, "store_analysis_metadata", lambda *a, **k: None
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *a, **k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1="A",
        group2="B",
        method="wilcoxon",
        n_top_genes=1,
    )

    out = await differential_expression("d13", ctx, params)

    assert captured["dtype"] == np.float32
    assert out.top_genes == ["gene_0"]


@pytest.mark.asyncio
async def test_differential_specific_group_raises_when_rank_results_have_no_names(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_de_adata()
    ctx = DummyCtx(adata)

    def _fake_rank(adata, **_kwargs):
        adata.uns["rank_genes_groups"] = {
            "pvals_adj": np.zeros((1,), dtype=[("A", "f8")])
        }

    monkeypatch.setattr(differential_mod.sc.tl, "rank_genes_groups", _fake_rank)
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(
            source="raw",
            X=adata.raw.X,
            var_names=adata.raw.var_names,
            is_integer_counts=True,
        ),
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1="A",
        group2="B",
        method="wilcoxon",
        n_top_genes=1,
    )

    with pytest.raises(ProcessingError, match="No DE genes found"):
        await differential_expression("d14", ctx, params)


@pytest.mark.asyncio
async def test_differential_specific_group_sets_none_mean_log2fc_when_all_genes_missing_in_raw(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_de_adata()
    ctx = DummyCtx(adata)

    monkeypatch.setattr(
        differential_mod.sc.tl,
        "rank_genes_groups",
        _fake_rank_genes_groups_factory({"A": ["missing_1", "missing_2"]}),
    )
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(
            source="raw",
            X=adata.raw.X,
            var_names=adata.raw.var_names,
            is_integer_counts=True,
        ),
    )
    monkeypatch.setattr(
        differential_mod, "store_analysis_metadata", lambda *a, **k: None
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *a, **k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1="A",
        group2="B",
        method="wilcoxon",
        n_top_genes=2,
    )
    out = await differential_expression("d15", ctx, params)

    assert out.statistics["mean_log2fc"] is None
    assert any("genes not found in raw data" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_differential_specific_group_warns_on_extreme_fold_change(
    monkeypatch: pytest.MonkeyPatch,
):
    X = np.array(
        [
            [10000, 1],
            [9000, 1],
            [8000, 1],
            [7000, 1],
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
        ],
        dtype=np.float32,
    )
    adata = AnnData(X=X)
    adata.var_names = ["gene_0", "gene_1"]
    adata.obs["cluster"] = ["A"] * 4 + ["B"] * 4
    adata.raw = adata.copy()
    ctx = DummyCtx(adata)

    monkeypatch.setattr(
        differential_mod.sc.tl,
        "rank_genes_groups",
        _fake_rank_genes_groups_factory({"A": ["gene_0"]}),
    )
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(
            source="raw",
            X=adata.raw.X,
            var_names=adata.raw.var_names,
            is_integer_counts=True,
        ),
    )
    monkeypatch.setattr(
        differential_mod, "store_analysis_metadata", lambda *a, **k: None
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *a, **k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        group1="A",
        group2="B",
        method="wilcoxon",
        n_top_genes=1,
        pseudocount=1e-6,
    )
    await differential_expression("d16", ctx, params)

    assert any("Extreme fold change" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_run_pydeseq2_rejects_single_unique_group_when_group1_not_provided(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_pydeseq2_adata()
    adata.obs["cluster"] = ["A"] * adata.n_obs
    ctx = DummyCtx(adata)

    _patch_pydeseq2_modules(monkeypatch)
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(X=adata.X, var_names=adata.var_names),
    )
    monkeypatch.setattr(
        differential_mod, "check_is_integer_counts", lambda _x: (True, None, None)
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
        group1=None,
    )

    with pytest.raises(DataError, match="Need at least 2 groups"):
        await _run_pydeseq2("d17", ctx, params)


@pytest.mark.asyncio
async def test_run_pydeseq2_group1_vs_rest_path(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_pydeseq2_adata()
    ctx = DummyCtx(adata)
    captured: dict[str, object] = {}
    results_df = pd.DataFrame(
        {"padj": [0.01, 0.2], "log2FoldChange": [1.1, -0.2]},
        index=["gene_0", "gene_1"],
    )
    _install_fake_pydeseq2(monkeypatch, results_df, captured)

    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(X=adata.X, var_names=adata.var_names),
    )
    monkeypatch.setattr(
        differential_mod, "check_is_integer_counts", lambda _x: (True, None, None)
    )
    monkeypatch.setattr(
        differential_mod, "store_analysis_metadata", lambda *a, **k: None
    )
    monkeypatch.setattr(
        differential_mod, "export_analysis_result", lambda *a, **k: None
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
        group1="A",
        group2=None,
    )
    out = await _run_pydeseq2("d18", ctx, params)

    assert out.comparison == "A vs rest"
    assert captured["contrast"] == ["condition", "A", "rest"]
    assert "_de_condition" not in adata.obs.columns
    assert "_pseudobulk_id" not in adata.obs.columns


@pytest.mark.asyncio
async def test_run_pydeseq2_wraps_runtime_failure_as_processing_error(
    monkeypatch: pytest.MonkeyPatch,
):
    from types import ModuleType

    adata = _make_pydeseq2_adata()
    ctx = DummyCtx(adata)

    class _DDS:
        def __init__(self, counts, metadata, design_factors):
            del counts, metadata, design_factors

        def deseq2(self):
            raise RuntimeError("boom")

    class _Stats:
        def __init__(self, dds, contrast):
            del dds, contrast
            self.results_df = pd.DataFrame()

        def summary(self):
            return None

    m_pkg = ModuleType("pydeseq2")
    m_dds = ModuleType("pydeseq2.dds")
    m_ds = ModuleType("pydeseq2.ds")
    m_dds.DeseqDataSet = _DDS
    m_ds.DeseqStats = _Stats
    monkeypatch.setitem(__import__("sys").modules, "pydeseq2", m_pkg)
    monkeypatch.setitem(__import__("sys").modules, "pydeseq2.dds", m_dds)
    monkeypatch.setitem(__import__("sys").modules, "pydeseq2.ds", m_ds)

    _patch_pydeseq2_modules(monkeypatch, m_dds, m_ds)
    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(X=adata.X, var_names=adata.var_names),
    )
    monkeypatch.setattr(
        differential_mod, "check_is_integer_counts", lambda _x: (True, None, None)
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
        group1="A",
        group2="B",
    )

    with pytest.raises(ProcessingError, match="PyDESeq2 analysis failed: boom"):
        await _run_pydeseq2("d19", ctx, params)


@pytest.mark.asyncio
async def test_run_pydeseq2_raises_when_no_top_genes_after_dropna(
    monkeypatch: pytest.MonkeyPatch,
):
    adata = _make_pydeseq2_adata()
    ctx = DummyCtx(adata)
    captured: dict[str, object] = {}
    results_df = pd.DataFrame(
        {"padj": [np.nan, np.nan], "log2FoldChange": [0.3, -0.1]},
        index=["gene_0", "gene_1"],
    )
    _install_fake_pydeseq2(monkeypatch, results_df, captured)

    monkeypatch.setattr(
        differential_mod,
        "get_raw_data_source",
        lambda *_a, **_k: SimpleNamespace(X=adata.X, var_names=adata.var_names),
    )
    monkeypatch.setattr(
        differential_mod, "check_is_integer_counts", lambda _x: (True, None, None)
    )

    params = DifferentialExpressionParameters(
        group_key="cluster",
        method="pydeseq2",
        sample_key="sample",
        group1="A",
        group2="B",
    )

    with pytest.raises(ProcessingError, match="No DE genes found"):
        await _run_pydeseq2("d20", ctx, params)
