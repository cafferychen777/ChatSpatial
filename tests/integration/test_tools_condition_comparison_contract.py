"""Integration contracts for tools.condition_comparison.compare_conditions."""

from __future__ import annotations

import pytest

import pandas as pd

from chatspatial.models.analysis import ConditionComparisonResult
from chatspatial.models.data import ConditionComparisonParameters
from chatspatial.tools import condition_comparison as cc_module
from chatspatial.tools.condition_comparison import compare_conditions
from chatspatial.utils.exceptions import ParameterError


class DummyCtx:
    def __init__(self, adata):
        self._adata = adata
        self.info_logs: list[str] = []
        self.warn_logs: list[str] = []

    async def get_adata(self, data_id: str):
        return self._adata

    async def info(self, msg: str):
        self.info_logs.append(msg)

    async def warning(self, msg: str):
        self.warn_logs.append(msg)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_compare_conditions_rejects_missing_condition_value(
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
        condition1="missing_condition",
        condition2="control",
        sample_key="sample",
    )
    with pytest.raises(ParameterError, match="missing_condition"):
        await compare_conditions("d1", ctx, params)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_compare_conditions_uses_global_branch_and_returns_contract(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["condition"] = ["treated"] * 30 + ["control"] * 30
    adata.obs["sample"] = ["s1"] * 15 + ["s2"] * 15 + ["s3"] * 15 + ["s4"] * 15
    ctx = DummyCtx(adata)

    class _RawStub:
        def __init__(self, X, var_names):
            self.X = X
            self.var_names = var_names

    monkeypatch.setattr(
        cc_module,
        "get_raw_data_source",
        lambda *args, **kwargs: _RawStub(adata.X, adata.var_names),
    )

    fake_de_results = pd.DataFrame(
        {"log2FoldChange": [1.5, -0.8], "padj": [0.01, 0.05]},
        index=["gene1", "gene2"],
    )

    async def fake_run_global(*args, **kwargs):
        data_id = kwargs.get("data_id", "")
        results_key = kwargs.get("results_key", "")
        return (
            ConditionComparisonResult(
                data_id=data_id,
                method="pseudobulk",
                comparison="treated vs control",
                condition_key="condition",
                condition1="treated",
                condition2="control",
                sample_key="sample",
                cell_type_key=None,
                n_samples_condition1=2,
                n_samples_condition2=2,
                global_n_significant=3,
                global_top_upregulated=[],
                global_top_downregulated=[],
                cell_type_results=None,
                results_key=results_key,
                statistics={"analysis_type": "global", "n_significant_genes": 3},
            ),
            fake_de_results,
        )

    monkeypatch.setattr(cc_module, "_run_global_comparison", fake_run_global)
    captured_meta: dict[str, object] = {}
    monkeypatch.setattr(
        cc_module,
        "store_analysis_metadata",
        lambda _a, **kw: captured_meta.update(kw),
    )
    monkeypatch.setattr(
        cc_module, "export_analysis_result", lambda *args, **kwargs: None
    )

    params = ConditionComparisonParameters(
        condition_key="condition",
        condition1="treated",
        condition2="control",
        sample_key="sample",
    )
    result = await compare_conditions("d1", ctx, params)

    assert isinstance(result, ConditionComparisonResult)
    assert result.data_id == "d1"
    assert result.results_key == "condition_comparison_treated_vs_control"
    assert result.n_samples_condition1 == 2
    assert result.n_samples_condition2 == 2
    assert "condition_comparison_treated_vs_control" in adata.uns
    # Regression: gene-level DE results must be stored for export
    de_key = "condition_comparison_treated_vs_control_de_results"
    assert de_key in adata.uns
    assert isinstance(adata.uns[de_key], pd.DataFrame)
    assert list(adata.uns[de_key].columns) == ["log2FoldChange", "padj"]
    # Regression: analysis_name must be comparison-specific (not generic)
    # so multiple comparisons don't overwrite each other's provenance
    assert captured_meta["analysis_name"] == "condition_comparison_treated_vs_control"
