"""Unit tests for the rule-based annotation methods marker_score and threshold."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from pydantic import ValidationError

from chatspatial.models.analysis import AnnotationResult
from chatspatial.models.data import AnnotationParameters
from chatspatial.tools import annotation as ann
from chatspatial.utils.exceptions import (
    DataError,
    DataNotFoundError,
    ParameterError,
)

pytestmark = pytest.mark.unit

TC = ["KRT14", "KRT5", "MKI67"]
LE = ["VIM", "FN1", "COL1A1"]
OTHER = ["G1", "G2", "G3", "G4"]
DOMAINS = ["0", "1", "2", "3", "4", "5"]


class DummyCtx:
    def __init__(self, adata):
        self.adata = adata
        self.warnings: list[str] = []

    async def get_adata(self, _data_id: str):
        return self.adata

    async def set_adata(self, _data_id: str, adata):
        self.adata = adata

    async def warning(self, msg: str):
        self.warnings.append(msg)


@pytest.fixture
def no_export(monkeypatch: pytest.MonkeyPatch):
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )
    return captured


def _domain_adata() -> ad.AnnData:
    """Six domains of 10 spots with integer counts.

    Tumor-core markers are highest in domains 0 and 1. Leading-edge markers are
    highest in domain 1 (a tumor-core domain), then 2, then 3.
    """
    rng = np.random.default_rng(0)
    genes = TC + LE + OTHER
    tc_level = {"0": 30, "1": 25, "2": 2, "3": 2, "4": 1, "5": 1}
    le_level = {"0": 1, "1": 40, "2": 30, "3": 20, "4": 2, "5": 1}
    rows = []
    for domain in DOMAINS:
        for _ in range(10):
            row = [tc_level[domain]] * len(TC) + [le_level[domain]] * len(LE)
            row += list(rng.integers(200, 220, size=len(OTHER)))
            rows.append(row)
    adata = ad.AnnData(np.asarray(rows, dtype=np.float32))
    adata.var_names = genes
    adata.obs_names = [f"spot_{i}" for i in range(adata.n_obs)]
    adata.obs["domain"] = pd.Categorical(np.repeat(DOMAINS, 10), categories=DOMAINS)
    return adata


def _params(**overrides) -> AnnotationParameters:
    values = {
        "method": "marker_score",
        "cluster_label": "domain",
        "marker_genes": {"TC": TC, "LE": LE},
    }
    values.update(overrides)
    return AnnotationParameters(**values)


@pytest.mark.asyncio
async def test_marker_score_rank_then_drop_reproduces_the_tc_le_rule(no_export):
    ctx = DummyCtx(_domain_adata())

    out = await ann.annotate_cell_types("d1", ctx, _params(output_key="tc_le_region"))

    # Top two LE domains are 1 and 2; domain 1 is already TC, so LE keeps only 2.
    assert out.group_labels == {
        "0": "TC",
        "1": "TC",
        "2": "LE",
        "3": "Other",
        "4": "Other",
        "5": "Other",
    }
    assert out.output_key == "tc_le_region"
    assert out.cell_types == ["TC", "LE", "Other"]
    assert out.counts == {"TC": 20, "LE": 10, "Other": 30}
    assert out.markers_used == {"TC": TC, "LE": LE}
    labels = ctx.adata.obs["tc_le_region"]
    assert list(labels.cat.categories) == ["TC", "LE", "Other"]
    assert labels.iloc[0] == "TC" and labels.iloc[25] == "LE"
    table = ctx.adata.uns["tc_le_region_group_scores"]
    assert table.loc["2", "assigned_label"] == "LE"
    assert int(table.loc["2", "n_obs"]) == 10
    assert no_export["analysis_name"] == "annotation_marker_score_tc_le_region"
    assert "tc_le_region_group_scores" in no_export["results_keys"]["uns"]


@pytest.mark.asyncio
async def test_marker_score_drop_then_rank_fills_each_label(no_export):
    ctx = DummyCtx(_domain_adata())

    out = await ann.annotate_cell_types(
        "d1", ctx, _params(marker_score_exclusion="drop_then_rank")
    )

    assert out.output_key == "cell_type_marker_score"
    assert {g for g, label in out.group_labels.items() if label == "LE"} == {"2", "3"}
    assert {g for g, label in out.group_labels.items() if label == "TC"} == {"0", "1"}


@pytest.mark.asyncio
async def test_marker_score_is_mean_log_normalized_expression(no_export):
    adata = _domain_adata()
    ctx = DummyCtx(adata.copy())

    out = await ann.annotate_cell_types("d1", ctx, _params(marker_score_n_top=1))

    expected = adata.copy()
    sc.pp.normalize_total(expected, target_sum=1e4)
    sc.pp.log1p(expected)
    for domain in DOMAINS:
        mask = (expected.obs["domain"] == domain).to_numpy()
        for label, genes in {"TC": TC, "LE": LE}.items():
            value = float(np.asarray(expected[mask, genes].X).mean())
            assert out.group_scores[domain][label] == pytest.approx(value, abs=1e-5)
    assert ctx.adata.uns["cell_type_marker_score_rule"]["expression_source"] == (
        "log1p-normalized current"
    )


@pytest.mark.asyncio
async def test_marker_score_reads_counts_from_raw_when_x_is_normalized(no_export):
    adata = _domain_adata()
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    ctx = DummyCtx(adata)

    out = await ann.annotate_cell_types("d1", ctx, _params())

    tc_mask = (adata.obs["domain"] == "0").to_numpy()
    assert out.group_scores["0"]["TC"] == pytest.approx(
        float(np.asarray(adata[tc_mask, TC].X).mean()), abs=1e-5
    )


@pytest.mark.asyncio
async def test_marker_score_skips_missing_markers_with_warning(no_export):
    ctx = DummyCtx(_domain_adata())

    out = await ann.annotate_cell_types(
        "d1", ctx, _params(marker_genes={"TC": [*TC, "PCNA"], "LE": LE})
    )

    assert out.markers_used["TC"] == TC
    assert any("PCNA" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_marker_score_rejects_label_without_any_present_gene(no_export):
    ctx = DummyCtx(_domain_adata())

    with pytest.raises(DataError, match="None of the marker genes for 'LE'"):
        await ann.annotate_cell_types(
            "d1", ctx, _params(marker_genes={"TC": TC, "LE": ["ABSENT1"]})
        )


@pytest.mark.asyncio
async def test_marker_score_leaves_unlabeled_groups_missing(no_export):
    adata = _domain_adata()
    domain = adata.obs["domain"].astype(object).to_numpy()
    domain[:3] = None
    adata.obs["domain"] = pd.Categorical(domain, categories=DOMAINS)
    ctx = DummyCtx(adata)

    out = await ann.annotate_cell_types("d1", ctx, _params())

    assert ctx.adata.obs["cell_type_marker_score"].iloc[:3].isna().all()
    assert sum(out.counts.values()) == adata.n_obs - 3


@pytest.mark.asyncio
async def test_marker_score_rejects_output_key_equal_to_cluster_label(no_export):
    ctx = DummyCtx(_domain_adata())
    with pytest.raises(ParameterError, match="output_key"):
        await ann.annotate_cell_types("d1", ctx, _params(output_key="domain"))


def test_marker_score_parameter_requirements():
    with pytest.raises(ValidationError, match="cluster_label"):
        AnnotationParameters(method="marker_score", marker_genes={"TC": TC})
    with pytest.raises(ValidationError, match="marker_genes"):
        AnnotationParameters(method="marker_score", cluster_label="domain")
    with pytest.raises(ValidationError, match="without genes"):
        _params(marker_genes={"TC": TC, "LE": []})
    with pytest.raises(ValidationError, match="unassigned"):
        _params(marker_score_unassigned_label="TC")
    with pytest.raises(ValidationError):
        _params(marker_score_n_top=0)


def _proportion_adata() -> ad.AnnData:
    adata = ad.AnnData(np.ones((6, 3), dtype=np.float32))
    adata.obs_names = [f"spot_{i}" for i in range(6)]
    tumor = np.array([0.0, 0.1, 0.15, 0.2, 0.9, 0.149])
    adata.obsm["deconvolution_rctd"] = np.column_stack([tumor, 1 - tumor])
    adata.uns["deconvolution_rctd_cell_types"] = ["Tumor_cells", "Fibroblasts"]
    adata.obs["deconvolution_rctd_Tumor_cells"] = tumor
    adata.obs["deconvolution_rctd_Fibroblasts"] = 1 - tumor
    return adata


def _threshold_params(**overrides) -> AnnotationParameters:
    values = {
        "method": "threshold",
        "threshold_key": "deconvolution_rctd_Tumor_cells",
        "threshold_value": 0.15,
        "threshold_pass_label": "tumor_enriched",
        "threshold_fail_label": "low_tumor",
        "output_key": "rctd_tumor_class",
    }
    values.update(overrides)
    return AnnotationParameters(**values)


@pytest.mark.asyncio
async def test_threshold_labels_deconvolved_proportions(no_export):
    ctx = DummyCtx(_proportion_adata())

    out = await ann.annotate_cell_types("d1", ctx, _threshold_params())

    labels = ctx.adata.obs["rctd_tumor_class"]
    assert labels.tolist() == [
        "low_tumor",
        "low_tumor",
        "tumor_enriched",
        "tumor_enriched",
        "tumor_enriched",
        "low_tumor",
    ]
    assert list(labels.cat.categories) == ["tumor_enriched", "low_tumor"]
    assert out.counts == {"tumor_enriched": 3, "low_tumor": 3}
    assert out.cell_types == ["tumor_enriched", "low_tumor"]
    assert no_export["parameters"]["comparison"] == ">="
    assert no_export["analysis_name"] == "annotation_threshold_rctd_tumor_class"


@pytest.mark.parametrize(
    ("comparison", "n_pass"),
    [(">=", 3), (">", 2), ("<=", 4), ("<", 3)],
)
@pytest.mark.asyncio
async def test_threshold_comparisons(no_export, comparison, n_pass):
    ctx = DummyCtx(_proportion_adata())

    out = await ann.annotate_cell_types(
        "d1",
        ctx,
        AnnotationParameters(
            method="threshold",
            threshold_key="deconvolution_rctd_Tumor_cells",
            threshold_value=0.15,
            threshold_comparison=comparison,
        ),
    )

    assert out.output_key == "cell_type_threshold"
    assert out.counts["pass"] == n_pass
    assert out.counts["fail"] == 6 - n_pass


@pytest.mark.asyncio
async def test_threshold_leaves_missing_values_unlabeled(no_export):
    adata = _proportion_adata()
    adata.obs.loc["spot_0", "deconvolution_rctd_Tumor_cells"] = np.nan
    ctx = DummyCtx(adata)

    out = await ann.annotate_cell_types("d1", ctx, _threshold_params())

    assert pd.isna(ctx.adata.obs["rctd_tumor_class"].iloc[0])
    assert out.counts == {"tumor_enriched": 3, "low_tumor": 2}
    assert any("no value" in w for w in ctx.warnings)


@pytest.mark.asyncio
async def test_threshold_missing_key_lists_deconvolution_columns(no_export):
    ctx = DummyCtx(_proportion_adata())

    with pytest.raises(DataNotFoundError, match="deconvolution_rctd_Tumor_cells"):
        await ann.annotate_cell_types(
            "d1", ctx, _threshold_params(threshold_key="tumor_fraction")
        )


@pytest.mark.asyncio
async def test_threshold_rejects_non_numeric_column(no_export):
    adata = _proportion_adata()
    adata.obs["label"] = pd.Categorical(["a", "b"] * 3)
    ctx = DummyCtx(adata)

    with pytest.raises(ParameterError, match="numeric"):
        await ann.annotate_cell_types(
            "d1", ctx, _threshold_params(threshold_key="label")
        )


def test_threshold_parameter_requirements():
    with pytest.raises(ValidationError, match="threshold_key"):
        AnnotationParameters(method="threshold", threshold_value=0.1)
    with pytest.raises(ValidationError, match="threshold_value"):
        AnnotationParameters(method="threshold", threshold_key="x")
    with pytest.raises(ValidationError, match="must differ"):
        _threshold_params(threshold_pass_label="a", threshold_fail_label="a")
    with pytest.raises(ValidationError):
        _threshold_params(threshold_comparison="==")
    with pytest.raises(ValidationError):
        _threshold_params(threshold_value=float("nan"))


def test_marker_score_result_serializes_group_scores():
    result = AnnotationResult(
        data_id="d1",
        method="marker_score",
        output_key="region",
        cell_types=["TC", "Other"],
        counts={"TC": 1, "Other": 1},
        group_labels={"0": "TC", "1": "Other"},
        group_scores={"0": {"TC": 1.5}, "1": {"TC": 0.5}},
        markers_used={"TC": ["KRT5"]},
    )
    payload = result.model_dump(mode="json")
    assert payload["group_scores"]["0"]["TC"] == 1.5
    assert payload["group_labels"]["1"] == "Other"
