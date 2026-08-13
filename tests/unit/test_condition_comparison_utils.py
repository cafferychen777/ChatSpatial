"""Unit tests for condition_comparison low-level utilities."""

from __future__ import annotations

from types import ModuleType

import numpy as np
import pandas as pd
import pytest

from chatspatial.models.analysis import DEGene
from chatspatial.models.data import ConditionComparisonParameters
from chatspatial.tools import condition_comparison as cc_module
from chatspatial.tools.condition_comparison import (
    _create_pseudobulk,
    _run_deseq2,
    _run_global_comparison,
    _run_stratified_comparison,
    _validate_sample_condition_mapping,
    compare_conditions,
)
from chatspatial.utils.exceptions import (
    DataError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


class DummyCtx:
    def __init__(self, adata=None):
        self._adata = adata
        self.info_logs: list[str] = []
        self.warn_logs: list[str] = []

    async def get_adata(self, _data_id: str):
        return self._adata

    async def set_adata(self, _data_id: str, adata):
        self._adata = adata

    async def info(self, msg: str):
        self.info_logs.append(msg)

    async def warning(self, msg: str):
        self.warn_logs.append(msg)


def _patch_pydeseq2_modules(
    monkeypatch: pytest.MonkeyPatch,
    dds_module: object,
    ds_module: object,
) -> None:
    modules = {
        "pydeseq2.dds": dds_module,
        "pydeseq2.ds": ds_module,
    }
    monkeypatch.setattr(
        cc_module,
        "require_module",
        lambda _dependency, module_name, *_args, **_kwargs: modules[module_name],
    )


def test_create_pseudobulk_aggregates_per_sample(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["condition"] = ["treated"] * 30 + ["control"] * 30
    adata.obs["sample"] = ["s1"] * 15 + ["s2"] * 15 + ["s3"] * 15 + ["s4"] * 15

    counts_df, metadata_df, cell_counts = _create_pseudobulk(
        adata=adata,
        raw_X=adata.X,
        var_names=adata.var_names,
        sample_key="sample",
        condition_key="condition",
        min_cells_per_sample=10,
    )

    assert counts_df.shape[0] == 4
    assert list(metadata_df.columns) == ["condition"]
    assert metadata_df.loc["s1", "condition"] == "treated"
    assert metadata_df.loc["s4", "condition"] == "control"
    assert cell_counts == {"s1": 15, "s2": 15, "s3": 15, "s4": 15}


def test_create_pseudobulk_raises_when_no_sample_meets_min_cells(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["condition"] = ["treated"] * 30 + ["control"] * 30
    adata.obs["sample"] = [f"s{i}" for i in range(60)]  # one cell per sample

    with pytest.raises(DataError, match="No samples have >="):
        _create_pseudobulk(
            adata=adata,
            raw_X=adata.X,
            var_names=adata.var_names,
            sample_key="sample",
            condition_key="condition",
            min_cells_per_sample=2,
        )


def test_create_pseudobulk_cell_type_filter_only_keeps_selected_type(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["condition"] = ["treated"] * 30 + ["control"] * 30
    adata.obs["sample"] = ["s1"] * 15 + ["s2"] * 15 + ["s3"] * 15 + ["s4"] * 15
    adata.obs["cell_type"] = ["T"] * 20 + ["B"] * 40

    counts_df, metadata_df, cell_counts = _create_pseudobulk(
        adata=adata,
        raw_X=adata.X,
        var_names=adata.var_names,
        sample_key="sample",
        condition_key="condition",
        cell_type="T",
        cell_type_key="cell_type",
        min_cells_per_sample=5,
    )

    # Only samples containing enough T cells survive.
    assert counts_df.shape[0] >= 1
    assert set(metadata_df["condition"]).issubset({"treated", "control"})
    assert all(v >= 5 for v in cell_counts.values())


def test_run_deseq2_filters_thresholds_and_nan(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
):
    counts_df = pd.DataFrame(
        [[10, 5, 2, 1], [12, 3, 1, 2], [6, 10, 2, 2], [8, 11, 1, 3]],
        index=["s1", "s2", "s3", "s4"],
        columns=["gene_up", "gene_down", "gene_nonsig", "gene_nan"],
    )
    metadata_df = pd.DataFrame(
        {"condition": ["treated", "treated", "control", "control"]},
        index=counts_df.index,
    )

    captured: dict[str, object] = {}

    class FakeDeseqDataSet:
        def __init__(self, counts, metadata, design_factors):
            print("dataset init output")
            captured["counts"] = counts
            captured["metadata"] = metadata
            captured["design_factors"] = design_factors

        def deseq2(self):
            print("deseq2 output")
            captured["deseq2_called"] = True

    class FakeDeseqStats:
        def __init__(self, dds, contrast):
            print("stats init output")
            captured["dds"] = dds
            captured["contrast"] = contrast
            self.results_df = pd.DataFrame(
                {
                    "log2FoldChange": [1.2, -1.5, 0.2, 0.9],
                    "pvalue": [0.001, 0.002, 0.3, 0.01],
                    "padj": [0.01, 0.03, 0.2, np.nan],
                },
                index=["gene_up", "gene_down", "gene_nonsig", "gene_nan"],
            )

        def summary(self):
            print("summary output")
            captured["summary_called"] = True

    dds_mod = ModuleType("pydeseq2.dds")
    dds_mod.DeseqDataSet = FakeDeseqDataSet
    ds_mod = ModuleType("pydeseq2.ds")
    ds_mod.DeseqStats = FakeDeseqStats
    _patch_pydeseq2_modules(monkeypatch, dds_mod, ds_mod)

    top_up, top_down, n_significant, results_df, n_up, n_down = _run_deseq2(
        counts_df,
        metadata_df,
        condition1="treated",
        condition2="control",
        n_top_genes=1,
        padj_threshold=0.05,
        log2fc_threshold=0.5,
    )

    assert captured["design_factors"] == "condition"
    assert captured["contrast"] == ["condition", "treated", "control"]
    assert captured["deseq2_called"] is True
    assert captured["summary_called"] is True
    assert capsys.readouterr().out == ""
    assert n_significant == 2
    assert len(results_df) == 3
    assert top_up[0].gene == "gene_up"
    assert top_down[0].gene == "gene_down"


def test_run_deseq2_preserves_managed_submodule_error(
    monkeypatch: pytest.MonkeyPatch,
):
    def _missing_submodule(*_args, **_kwargs):
        raise DependencyError("pydeseq2.dds is unavailable")

    monkeypatch.setattr(cc_module, "require_module", _missing_submodule)

    with pytest.raises(DependencyError, match="pydeseq2.dds is unavailable"):
        _run_deseq2(
            pd.DataFrame([[1]], index=["s1"], columns=["g1"]),
            pd.DataFrame({"condition": ["treated"]}, index=["s1"]),
            condition1="treated",
            condition2="control",
        )


@pytest.mark.asyncio
async def test_run_global_comparison_requires_two_samples_per_condition(monkeypatch):
    adata = type("A", (), {"obs": pd.DataFrame()})()
    ctx = DummyCtx()
    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        min_cells_per_sample=1,
    )

    counts_df = pd.DataFrame([[1], [2]], index=["s1", "s2"], columns=["g1"])
    metadata_df = pd.DataFrame(
        {"condition": ["treated", "control"]}, index=["s1", "s2"]
    )
    monkeypatch.setattr(
        cc_module,
        "_create_pseudobulk",
        lambda *_args, **_kwargs: (counts_df, metadata_df, {"s1": 5, "s2": 5}),
    )

    with pytest.raises(DataError, match="at least 2 samples per condition"):
        await _run_global_comparison(
            adata,
            np.array([[1], [2]]),
            pd.Index(["g1"]),
            ctx,
            params,
        )


@pytest.mark.asyncio
async def test_run_global_comparison_wraps_deseq2_errors(monkeypatch):
    adata = type("A", (), {"obs": pd.DataFrame()})()
    ctx = DummyCtx()
    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        min_cells_per_sample=1,
    )

    counts_df = pd.DataFrame(
        [[1], [2], [3], [4]],
        index=["s1", "s2", "s3", "s4"],
        columns=["g1"],
    )
    metadata_df = pd.DataFrame(
        {
            "condition": ["treated", "treated", "control", "control"],
        },
        index=counts_df.index,
    )
    monkeypatch.setattr(
        cc_module,
        "_create_pseudobulk",
        lambda *_args, **_kwargs: (counts_df, metadata_df, {}),
    )
    monkeypatch.setattr(
        cc_module,
        "_run_deseq2",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("deseq failed")),
    )

    with pytest.raises(ProcessingError, match="DESeq2 analysis failed"):
        await _run_global_comparison(
            adata,
            np.array([[1], [2], [3], [4]]),
            pd.Index(["g1"]),
            ctx,
            params,
        )


@pytest.mark.asyncio
async def test_run_global_comparison_preserves_dependency_errors(monkeypatch):
    adata = type("A", (), {"obs": pd.DataFrame()})()
    ctx = DummyCtx()
    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        min_cells_per_sample=1,
    )
    counts_df = pd.DataFrame(
        [[1], [2], [3], [4]],
        index=["s1", "s2", "s3", "s4"],
        columns=["g1"],
    )
    metadata_df = pd.DataFrame(
        {"condition": ["treated", "treated", "control", "control"]},
        index=counts_df.index,
    )
    monkeypatch.setattr(
        cc_module,
        "_create_pseudobulk",
        lambda *_args, **_kwargs: (counts_df, metadata_df, {}),
    )
    monkeypatch.setattr(
        cc_module,
        "_run_deseq2",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("pydeseq2.dds is unavailable")
        ),
    )

    with pytest.raises(DependencyError, match="pydeseq2.dds is unavailable"):
        await _run_global_comparison(
            adata,
            np.array([[1], [2], [3], [4]]),
            pd.Index(["g1"]),
            ctx,
            params,
        )


@pytest.mark.asyncio
async def test_run_global_comparison_success_builds_contract(monkeypatch):
    adata = type("A", (), {"obs": pd.DataFrame()})()
    ctx = DummyCtx()
    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        min_cells_per_sample=1,
        n_top_genes=3,
        padj_threshold=0.05,
    )

    counts_df = pd.DataFrame(
        [[1], [2], [3], [4]],
        index=["s1", "s2", "s3", "s4"],
        columns=["g1"],
    )
    metadata_df = pd.DataFrame(
        {"condition": ["treated", "treated", "control", "control"]},
        index=counts_df.index,
    )
    monkeypatch.setattr(
        cc_module,
        "_create_pseudobulk",
        lambda *_args, **_kwargs: (counts_df, metadata_df, {}),
    )

    top_up = [DEGene(gene="g_up", log2fc=1.0, pvalue=0.001, padj=0.01)]
    top_down = [DEGene(gene="g_down", log2fc=-1.0, pvalue=0.002, padj=0.02)]
    monkeypatch.setattr(
        cc_module,
        "_run_deseq2",
        lambda *_args, **_kwargs: (top_up, top_down, 2, pd.DataFrame(), 1, 1),
    )

    result, results_df = await _run_global_comparison(
        adata,
        np.array([[1], [2], [3], [4]]),
        pd.Index(["g1"]),
        ctx,
        params,
    )

    assert result.method == "pseudobulk"
    assert result.comparison == "treated vs control"
    assert result.global_n_significant == 2
    assert result.statistics["analysis_type"] == "global"
    assert result.statistics["n_pseudobulk_samples"] == 4
    assert result.statistics["n_upregulated"] == 1
    assert result.statistics["n_downregulated"] == 1
    assert any("Created 4 pseudobulk samples" in m for m in ctx.info_logs)
    assert any("Found 2 significant DE genes" in m for m in ctx.info_logs)


@pytest.mark.asyncio
async def test_run_stratified_comparison_no_valid_cell_types_raises(monkeypatch):
    n_cells = 12
    adata = type(
        "A",
        (),
        {
            "obs": pd.DataFrame(
                {
                    "condition": ["treated"] * 6 + ["control"] * 6,
                    "sample": ["s1", "s1", "s2", "s2", "s3", "s3"] * 2,
                    "cell_type": ["small"] * 3 + ["large"] * 9,
                }
            )
        },
    )()
    ctx = DummyCtx()
    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        cell_type_key="cell_type",
        min_cells_per_sample=2,
    )
    raw_X = np.ones((n_cells, 2), dtype=np.int64)
    var_names = pd.Index(["g1", "g2"])

    def _fake_create(*_args, **kwargs):
        if kwargs.get("cell_type") == "large":
            counts_df = pd.DataFrame([[1], [2]], index=["s1", "s3"], columns=["g1"])
            metadata_df = pd.DataFrame(
                {"condition": ["treated", "control"]}, index=counts_df.index
            )
            return counts_df, metadata_df, {}
        raise AssertionError("small should be skipped before pseudobulk")

    monkeypatch.setattr(cc_module, "_create_pseudobulk", _fake_create)

    with pytest.raises(ProcessingError, match="No cell types had sufficient samples"):
        await _run_stratified_comparison(adata, raw_X, var_names, ctx, params)

    assert any("Skipping small: only 3 cells" in m for m in ctx.warn_logs)
    assert any("Skipping large: insufficient samples" in m for m in ctx.warn_logs)


@pytest.mark.asyncio
async def test_run_stratified_comparison_continues_after_cell_type_failure(monkeypatch):
    n_cells = 24
    obs = pd.DataFrame(
        {
            "condition": ["treated"] * 12 + ["control"] * 12,
            "sample": (["s1"] * 6 + ["s2"] * 6 + ["s3"] * 6 + ["s4"] * 6),
            "cell_type": ["bad"] * 4 + ["good"] * 20,
        }
    )
    adata = type(
        "A",
        (),
        {
            "obs": obs,
            "__getitem__": lambda self, idx: type("B", (), {"obs": self.obs[idx]})(),
        },
    )()
    ctx = DummyCtx()
    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        cell_type_key="cell_type",
        min_cells_per_sample=2,
    )
    raw_X = np.ones((n_cells, 2), dtype=np.int64)
    var_names = pd.Index(["g1", "g2"])

    counts_df = pd.DataFrame(
        [[1, 2], [2, 1], [3, 4], [4, 3]],
        index=["s1", "s2", "s3", "s4"],
        columns=var_names,
    )
    metadata_df = pd.DataFrame(
        {"condition": ["treated", "treated", "control", "control"]},
        index=counts_df.index,
    )

    def _fake_create(*_args, **kwargs):
        if kwargs.get("cell_type") == "bad":
            raise RuntimeError("bad type failed")
        return counts_df, metadata_df, {}

    monkeypatch.setattr(cc_module, "_create_pseudobulk", _fake_create)
    monkeypatch.setattr(
        cc_module,
        "_run_deseq2",
        lambda *_args, **_kwargs: (
            [DEGene(gene="g_up", log2fc=1.2, pvalue=0.001, padj=0.01)],
            [DEGene(gene="g_down", log2fc=-1.1, pvalue=0.002, padj=0.02)],
            2,
            pd.DataFrame(),
            1,
            1,
        ),
    )

    result, combined_df = await _run_stratified_comparison(
        adata, raw_X, var_names, ctx, params
    )

    assert result.statistics["analysis_type"] == "cell_type_stratified"
    assert result.statistics["n_cell_types_analyzed"] == 1
    assert result.statistics["total_significant_hits"] == 2
    assert result.cell_type_results is not None
    assert result.cell_type_results[0].cell_type == "good"
    assert any("Analysis failed for bad: bad type failed" in m for m in ctx.warn_logs)
    assert any("good: 2 significant genes" in m for m in ctx.info_logs)


@pytest.mark.asyncio
async def test_run_stratified_comparison_preserves_dependency_errors(monkeypatch):
    obs = pd.DataFrame(
        {
            "condition": ["treated"] * 4 + ["control"] * 4,
            "sample": ["s1"] * 2 + ["s2"] * 2 + ["s3"] * 2 + ["s4"] * 2,
            "cell_type": ["T"] * 8,
        }
    )
    adata = type("A", (), {"obs": obs})()
    ctx = DummyCtx()
    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        cell_type_key="cell_type",
        min_cells_per_sample=1,
    )
    counts_df = pd.DataFrame(
        [[1], [2], [3], [4]],
        index=["s1", "s2", "s3", "s4"],
        columns=["g1"],
    )
    metadata_df = pd.DataFrame(
        {"condition": ["treated", "treated", "control", "control"]},
        index=counts_df.index,
    )
    monkeypatch.setattr(
        cc_module,
        "_create_pseudobulk",
        lambda *_args, **_kwargs: (counts_df, metadata_df, {}),
    )
    monkeypatch.setattr(
        cc_module,
        "_run_deseq2",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("pydeseq2.ds is unavailable")
        ),
    )

    with pytest.raises(DependencyError, match="pydeseq2.ds is unavailable"):
        await _run_stratified_comparison(
            adata,
            np.ones((8, 1), dtype=np.int64),
            pd.Index(["g1"]),
            ctx,
            params,
        )


@pytest.mark.asyncio
async def test_compare_conditions_missing_condition2_raises_parameter_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["condition"] = ["treated"] * 30 + ["control"] * 30
    adata.obs["sample"] = ["s1"] * 15 + ["s2"] * 15 + ["s3"] * 15 + ["s4"] * 15
    ctx = DummyCtx(adata)
    monkeypatch.setattr(
        cc_module,
        "require_module",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("dependency should load only after input validation")
        ),
    )

    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="missing",
        sample_key="sample",
    )
    with pytest.raises(ParameterError, match="Condition 'missing'"):
        await compare_conditions("d1", ctx, params)


@pytest.mark.asyncio
async def test_compare_conditions_min_samples_guard_for_condition2(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["condition"] = ["treated"] * 40 + ["control"] * 20
    adata.obs["sample"] = ["s1"] * 20 + ["s2"] * 20 + ["s3"] * 20
    ctx = DummyCtx(adata)

    class _RawStub:
        def __init__(self, X, var_names):
            self.X = X
            self.var_names = var_names

    def _raw_source(_adata, *, prefer_complete_genes, require_integer_counts):
        assert prefer_complete_genes is True
        assert require_integer_counts is True
        return _RawStub(adata.X, adata.var_names)

    monkeypatch.setattr(cc_module, "get_raw_data_source", _raw_source)

    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        min_samples_per_condition=2,
    )
    with pytest.raises(DataError, match="Insufficient samples for control: 1"):
        await compare_conditions("d1", ctx, params)


@pytest.mark.asyncio
async def test_compare_conditions_min_samples_guard_for_condition1(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["condition"] = ["treated"] * 20 + ["control"] * 40
    adata.obs["sample"] = ["s1"] * 20 + ["s2"] * 20 + ["s3"] * 20
    ctx = DummyCtx(adata)

    class _RawStub:
        def __init__(self, X, var_names):
            self.X = X
            self.var_names = var_names

    def _raw_source(_adata, *, prefer_complete_genes, require_integer_counts):
        assert prefer_complete_genes is True
        assert require_integer_counts is True
        return _RawStub(adata.X, adata.var_names)

    monkeypatch.setattr(cc_module, "get_raw_data_source", _raw_source)

    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
        min_samples_per_condition=2,
    )
    with pytest.raises(DataError, match="Insufficient samples for treated: 1"):
        await compare_conditions("d1", ctx, params)


def test_validate_sample_condition_mapping_clean():
    obs = pd.DataFrame(
        {"sample": ["S1", "S1", "S2", "S2"], "cond": ["A", "A", "B", "B"]}
    )
    result = _validate_sample_condition_mapping(obs, "sample", "cond")
    assert result["S1"] == "A"
    assert result["S2"] == "B"


def test_validate_sample_condition_mapping_rejects_mixed():
    obs = pd.DataFrame(
        {"sample": ["S1", "S1", "S2", "S2"], "cond": ["A", "B", "B", "B"]}
    )
    with pytest.raises(DataError, match="multiple conditions"):
        _validate_sample_condition_mapping(obs, "sample", "cond")


def test_create_pseudobulk_rejects_mixed_condition_sample(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["condition"] = ["treated"] * 30 + ["control"] * 30
    adata.obs["sample"] = ["s1"] * 15 + ["s2"] * 15 + ["s3"] * 15 + ["s4"] * 15
    # Corrupt s1: mix in a "control" cell
    adata.obs.iloc[0, adata.obs.columns.get_loc("condition")] = "control"

    with pytest.raises(DataError, match="multiple conditions"):
        _create_pseudobulk(
            adata=adata,
            raw_X=adata.X,
            var_names=adata.var_names,
            sample_key="sample",
            condition_key="condition",
            min_cells_per_sample=1,
        )


def test_run_deseq2_reports_expression_behind_each_fold_change(
    monkeypatch: pytest.MonkeyPatch,
):
    """A fold change alone cannot be judged.

    Doubling a gene at 5 counts and one at 5000 are different findings, so each
    result carries the per-condition means the ratio was taken over.
    """
    counts_df = pd.DataFrame(
        [[300, 5], [340, 7], [100, 6], [120, 4]],
        index=["s1", "s2", "s3", "s4"],
        columns=["gene_high", "gene_low"],
    )
    metadata_df = pd.DataFrame(
        {"condition": ["treated", "treated", "control", "control"]},
        index=counts_df.index,
    )

    class FakeDeseqDataSet:
        def __init__(self, counts, metadata, design_factors):
            pass

        def deseq2(self):
            pass

    class FakeDeseqStats:
        def __init__(self, dds, contrast):
            self.results_df = pd.DataFrame(
                {
                    "log2FoldChange": [1.6, 1.6],
                    "pvalue": [0.001, 0.001],
                    "padj": [0.01, 0.01],
                },
                index=["gene_high", "gene_low"],
            )

        def summary(self):
            pass

    dds_mod = ModuleType("pydeseq2.dds")
    dds_mod.DeseqDataSet = FakeDeseqDataSet
    ds_mod = ModuleType("pydeseq2.ds")
    ds_mod.DeseqStats = FakeDeseqStats
    _patch_pydeseq2_modules(monkeypatch, dds_mod, ds_mod)

    top_up, _down, _n_sig, _df, _up, _down_n = _run_deseq2(
        counts_df,
        metadata_df,
        condition1="treated",
        condition2="control",
        n_top_genes=5,
        padj_threshold=0.05,
        log2fc_threshold=0.5,
    )

    by_gene = {g.gene: g for g in top_up}
    # Identical fold changes, three orders of magnitude apart in expression.
    assert by_gene["gene_high"].mean_expr_condition1 == pytest.approx(320.0)
    assert by_gene["gene_high"].mean_expr_condition2 == pytest.approx(110.0)
    assert by_gene["gene_low"].mean_expr_condition1 == pytest.approx(6.0)
    assert by_gene["gene_low"].mean_expr_condition2 == pytest.approx(5.0)
