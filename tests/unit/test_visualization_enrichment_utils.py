"""Unit tests for enrichment visualization utility contracts."""

from __future__ import annotations

import sys
import warnings
from types import ModuleType

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import enrichment as viz_enrich
from chatspatial.utils.exceptions import (
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


class DummyCtx:
    def __init__(self):
        self.infos: list[str] = []

    async def info(self, msg: str):
        self.infos.append(msg)


class _FakeAxes:
    def __init__(self):
        self._fig, _ = plt.subplots()

    def get_figure(self):
        return self._fig


def test_get_score_columns_fallback_suffix_search(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["Wnt_score"] = 0.1
    adata.obs["NotScore"] = 0.2
    assert viz_enrich._get_score_columns(adata) == ["Wnt_score"]


def test_resolve_score_column_validates_and_defaults(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = 0.2
    cols = ["A_score"]
    assert viz_enrich._resolve_score_column(adata, None, cols) == "A_score"
    assert viz_enrich._resolve_score_column(adata, "A_score", cols) == "A_score"
    assert viz_enrich._resolve_score_column(adata, "A", cols) == "A_score"

    with pytest.raises(DataNotFoundError, match="Score column 'missing' not found"):
        viz_enrich._resolve_score_column(adata, "missing", cols)

    with pytest.raises(DataNotFoundError, match="No enrichment scores found"):
        viz_enrich._resolve_score_column(adata, None, [])


@pytest.mark.asyncio
async def test_create_enrichment_visualization_routes_violin(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = 0.3
    sentinel = object()

    def _fake_violin(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_enrich, "_create_enrichment_violin", _fake_violin)
    out = await viz_enrich._create_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="violin"),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_enrichment_visualization_requires_scores(minimal_spatial_adata):
    with pytest.raises(DataNotFoundError, match="No enrichment scores found"):
        await viz_enrich._create_enrichment_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="enrichment", subtype="violin"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
@pytest.mark.parametrize("subtype", ["spatial_score", "spatial", "violin"])
async def test_create_pathway_enrichment_visualization_routes_score_subtypes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, subtype: str
):
    sentinel = object()

    async def _fake_router(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_enrich, "_create_enrichment_visualization", _fake_router)
    out = await viz_enrich.create_pathway_enrichment_visualization(
        minimal_spatial_adata,
        VisualizationParameters(plot_type="enrichment", subtype=subtype),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_pathway_enrichment_visualization_requires_results_key(
    minimal_spatial_adata,
):
    with pytest.raises(DataNotFoundError, match="Enrichment results not found"):
        await viz_enrich.create_pathway_enrichment_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="enrichment", subtype="barplot"),
            context=DummyCtx(),
        )


def test_gsea_results_to_dataframe_and_term_resolution():
    df = viz_enrich._gsea_results_to_dataframe(
        {"PathA": {"NES": 1.2, "pval": 0.01}, "PathB": {"NES": -1.1, "pval": 0.2}}
    )
    assert set(df["Term"]) == {"PathA", "PathB"}

    df2 = pd.DataFrame({"pathway": ["A"], "pval": [0.1]})
    viz_enrich._ensure_term_column(df2)
    assert "Term" in df2.columns


def test_find_pvalue_column_priority():
    df = pd.DataFrame({"Adjusted P-value": [0.1], "pval": [0.2]})
    assert viz_enrich._find_pvalue_column(df) == "Adjusted P-value"
    df2 = pd.DataFrame({"fdr": [0.1]})
    assert viz_enrich._find_pvalue_column(df2) == "fdr"


def test_gsea_barplot_requires_pvalue_column():
    with pytest.raises(DataNotFoundError, match="requires tabular GSEA/ORA results"):
        viz_enrich._create_gsea_barplot(
            pd.DataFrame({"Term": ["PathA"], "method": ["spatial_enrichmap"]}),
            VisualizationParameters(plot_type="enrichment", subtype="barplot"),
        )


def test_ensure_term_column_raises_when_missing():
    df = pd.DataFrame({"x": [1], "y": [2]})
    with pytest.raises(DataNotFoundError, match="No pathway/term column found"):
        viz_enrich._ensure_term_column(df)


def test_create_enrichmap_single_score_requires_feature(minimal_spatial_adata):
    with pytest.raises(DataNotFoundError, match="Feature parameter required"):
        viz_enrich._create_enrichmap_single_score(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="enrichment", subtype="spatial_score"),
            "sample_1",
            em=object(),
            context=None,
        )


def test_create_enrichmap_spatial_requires_dependency(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    def _require_missing(*_args, **_kwargs):
        raise DependencyError("EnrichMap is required")

    monkeypatch.setattr(viz_enrich, "require", _require_missing)

    with pytest.raises(DependencyError, match="EnrichMap is required"):
        viz_enrich._create_enrichmap_spatial(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="enrichment", subtype="spatial_score"),
            score_cols=["A_score"],
            context=None,
        )


def test_gsea_results_to_dataframe_restores_numeric_object_columns():
    df = viz_enrich._gsea_results_to_dataframe(
        pd.DataFrame(
            {
                "Term": ["PathA", "PathB"],
                "FDR q-val": pd.Series(["0.05", "0.20"], dtype=object),
                "Lead_genes": ["GeneA;GeneB", "GeneC"],
            }
        )
    )

    assert df["FDR q-val"].tolist() == [0.05, 0.2]
    assert np.issubdtype(df["FDR q-val"].dtype, np.number)
    assert df["Lead_genes"].tolist() == ["GeneA;GeneB", "GeneC"]


def test_create_gsea_barplot_wraps_gseapy_errors(monkeypatch: pytest.MonkeyPatch):
    fake_gp = ModuleType("gseapy")

    def _boom(**_kwargs):
        raise RuntimeError("bad data")

    fake_gp.barplot = _boom
    monkeypatch.setitem(__import__("sys").modules, "gseapy", fake_gp)

    with pytest.raises(ProcessingError, match="gseapy.barplot failed"):
        viz_enrich._create_gsea_barplot(
            {"PathA": {"Adjusted P-value": 0.1}},
            VisualizationParameters(plot_type="enrichment", subtype="barplot"),
        )

    plt.close("all")


def test_ensure_enrichmap_compatibility_adds_minimum_metadata(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    if "spatial" in adata.uns:
        del adata.uns["spatial"]
    if "library_id" in adata.obs.columns:
        del adata.obs["library_id"]

    result = viz_enrich._ensure_enrichmap_compatibility(adata)

    assert "library_id" not in adata.obs.columns
    assert "spatial" not in adata.uns
    assert "library_id" in result.obs.columns
    assert "spatial" in result.uns
    assert "sample_1" in result.uns["spatial"]


def test_ensure_enrichmap_compatibility_uses_existing_spatial_library(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    if "library_id" in adata.obs.columns:
        del adata.obs["library_id"]
    adata.uns["spatial"] = {
        "V1_Human_Lymph_Node": {
            "images": {},
            "scalefactors": {"spot_diameter_fullres": 1.0},
        }
    }

    result = viz_enrich._ensure_enrichmap_compatibility(adata)

    assert "library_id" not in adata.obs.columns
    assert result.obs["library_id"].unique().tolist() == ["V1_Human_Lymph_Node"]
    assert list(result.uns["spatial"].keys()) == ["V1_Human_Lymph_Node"]
    assert "sample_1" not in result.uns["spatial"]


def test_get_score_columns_prefers_metadata(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = 0.1
    adata.obs["ssgsea_B"] = 0.2
    adata.uns["enrichment_spatial_metadata"] = {
        "parameters": {},
        "results_keys": {"obs": ["A_score", "missing_score"]},
    }
    adata.uns["enrichment_ssgsea_metadata"] = {
        "parameters": {},
        "results_keys": {"obs": ["ssgsea_B"]},
    }

    out = viz_enrich._get_score_columns(adata)
    assert out == ["A_score", "ssgsea_B"]


def test_get_score_columns_finds_parametrized_metadata_keys(minimal_spatial_adata):
    """Parametrized enrichment keys like enrichment_ssgsea_KEGG_Pathways must be found."""
    adata = minimal_spatial_adata.copy()
    adata.obs["ssgsea_PathA"] = 0.3
    adata.obs["Wnt_score"] = 0.5

    adata.uns["enrichment_ssgsea_KEGG_Pathways_metadata"] = {
        "parameters": {},
        "results_keys": {"obs": ["ssgsea_PathA"]},
    }
    adata.uns["enrichment_spatial_MSigDB_Hallmark_metadata"] = {
        "parameters": {},
        "results_keys": {"obs": ["Wnt_score"]},
    }

    out = viz_enrich._get_score_columns(adata)
    assert "ssgsea_PathA" in out
    assert "Wnt_score" in out
    assert len(out) == 2


def test_get_score_columns_prioritizes_latest_score_key(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["ssgsea_Old"] = 0.1
    adata.obs["ssgsea_New"] = 0.2
    adata.uns["enrichment_ssgsea_OldDB_metadata"] = {
        "parameters": {},
        "results_keys": {"obs": ["ssgsea_Old"]},
    }
    adata.uns["enrichment_ssgsea_NewDB_metadata"] = {
        "parameters": {},
        "results_keys": {"obs": ["ssgsea_New"]},
    }
    adata.uns["enrichment_latest_score_key"] = "enrichment_ssgsea_NewDB"

    out = viz_enrich._get_score_columns(adata)

    assert out == ["ssgsea_New", "ssgsea_Old"]


@pytest.mark.asyncio
async def test_create_enrichment_visualization_routes_spatial_prefix(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = 0.1
    sentinel = object()

    def _fake_enrichmap(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_enrich, "_create_enrichmap_spatial", _fake_enrichmap)

    out = await viz_enrich._create_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="spatial_score"),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_enrichment_visualization_default_routes_spatial_scatter(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = 0.1
    sentinel = object()

    async def _fake_spatial(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_enrich, "_create_enrichment_spatial", _fake_spatial)
    out = await viz_enrich._create_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="spatial"),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_pathway_enrichment_visualization_uses_alternate_result_keys(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["pathway_enrichment"] = {"PathA": {"Adjusted P-value": 0.03}}
    sentinel = object()

    def _fake_barplot(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_enrich, "_create_gsea_barplot", _fake_barplot)
    out = await viz_enrich.create_pathway_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="barplot"),
        context=DummyCtx(),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_pathway_enrichment_visualization_prefers_latest_results_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["gsea_results"] = pd.DataFrame(
        {"Term": ["Old_GSEA"], "FDR q-val": [0.001]}
    )
    adata.uns["enrichr_results"] = pd.DataFrame(
        {"Term": ["Latest_Enrichr"], "Adjusted P-value": [0.9]}
    )
    adata.uns[viz_enrich._LATEST_ENRICHMENT_RESULTS_KEY] = "enrichr_results"
    captured: dict[str, object] = {}

    def _fake_barplot(gsea_results, *_args, **_kwargs):
        captured["terms"] = list(gsea_results["Term"])
        return object()

    monkeypatch.setattr(viz_enrich, "_create_gsea_barplot", _fake_barplot)
    await viz_enrich.create_pathway_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="barplot"),
        context=DummyCtx(),
    )

    assert captured["terms"] == ["Latest_Enrichr"]


@pytest.mark.asyncio
async def test_create_pathway_enrichment_visualization_defaults_to_spatial_scores(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = 0.1
    adata.uns["enrichment_spatial_MSigDB_Hallmark_metadata"] = {
        "results_keys": {"obs": ["A_score"], "uns": []}
    }
    sentinel = object()

    async def _fake_spatial(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_enrich, "_create_enrichment_visualization", _fake_spatial)

    out = await viz_enrich.create_pathway_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="barplot"),
        context=DummyCtx(),
    )

    assert out is sentinel


@pytest.mark.asyncio
async def test_create_pathway_enrichment_visualization_uses_rank_genes_when_available(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = 0.1
    adata.uns["enrichment_spatial_MSigDB_Hallmark_metadata"] = {
        "results_keys": {"obs": ["A_score"], "uns": []}
    }
    adata.uns["rank_genes_groups"] = {"names": {"0": ["gene_0"]}}
    sentinel = object()

    def _fake_barplot(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(viz_enrich, "_create_gsea_barplot", _fake_barplot)

    out = await viz_enrich.create_pathway_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="barplot"),
        context=DummyCtx(),
    )

    assert out is sentinel


@pytest.mark.asyncio
async def test_create_pathway_enrichment_visualization_unknown_subtype_error(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.uns["gsea_results"] = {"PathA": {"Adjusted P-value": 0.03}}
    with pytest.raises(ParameterError, match="Unknown enrichment visualization type"):
        await viz_enrich.create_pathway_enrichment_visualization(
            adata,
            VisualizationParameters(plot_type="enrichment", subtype="weird"),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_create_pathway_enrichment_visualization_routes_enrichment_plot_and_dotplot(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.uns["gsea_results"] = {"PathA": {"Adjusted P-value": 0.03, "RES": [0.1, 0.2]}}
    sentinel_enrichment = object()
    sentinel_dot = object()

    monkeypatch.setattr(viz_enrich, "_create_gsea_enrichment_plot", lambda *_a, **_k: sentinel_enrichment)
    monkeypatch.setattr(viz_enrich, "_create_gsea_dotplot", lambda *_a, **_k: sentinel_dot)

    out_enrich = await viz_enrich.create_pathway_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="enrichment_plot"),
        context=DummyCtx(),
    )
    out_dot = await viz_enrich.create_pathway_enrichment_visualization(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="dotplot"),
        context=DummyCtx(),
    )
    assert out_enrich is sentinel_enrichment
    assert out_dot is sentinel_dot


def test_create_enrichment_violin_requires_cluster_key(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = 0.3
    with pytest.raises(ParameterError, match="requires 'cluster_key'"):
        viz_enrich._create_enrichment_violin(
            adata,
            VisualizationParameters(plot_type="enrichment", subtype="violin"),
            score_cols=["A_score"],
        )


def test_create_enrichment_violin_multifeature_layout(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    adata.obs["B_score"] = np.linspace(1.0, 0.0, adata.n_obs)
    fig = viz_enrich._create_enrichment_violin(
        adata,
        VisualizationParameters(
            plot_type="enrichment",
            subtype="violin",
            cluster_key="group",
            feature=["A_score", "B_score"],
        ),
        score_cols=["A_score", "B_score"],
    )
    assert len(fig.axes) == 2
    fig.clf()


def test_create_enrichment_violin_defaults_to_score_cols_and_single_axis(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    fig = viz_enrich._create_enrichment_violin(
        adata,
        VisualizationParameters(
            plot_type="enrichment",
            subtype="violin",
            cluster_key="group",
        ),
        score_cols=["A_score"],
    )
    assert len(fig.axes) == 1
    fig.clf()


@pytest.mark.asyncio
async def test_create_enrichment_spatial_multifeature_and_single_feature_paths(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    adata.obs["B_score"] = np.linspace(1.0, 0.0, adata.n_obs)

    seen: list[str] = []

    def _fake_plot_spatial_feature(_adata, feature, ax, params):
        seen.append(feature)
        ax.scatter([0], [0], c=[1.0])

    monkeypatch.setattr(viz_enrich, "plot_spatial_feature", _fake_plot_spatial_feature)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=".*not compatible with tight_layout.*",
            category=UserWarning,
        )
        fig_multi = await viz_enrich._create_enrichment_spatial(
            adata,
            VisualizationParameters(
                plot_type="enrichment",
                subtype="spatial",
                feature=["A_score", "B_score"],
            ),
            score_cols=["A_score", "B_score"],
            context=None,
        )
    assert seen[:2] == ["A_score", "B_score"]
    fig_multi.clf()

    ctx = DummyCtx()
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=".*not compatible with tight_layout.*",
            category=UserWarning,
        )
        fig_single = await viz_enrich._create_enrichment_spatial(
            adata,
            VisualizationParameters(
                plot_type="enrichment",
                subtype="spatial",
                feature="A",
                show_colorbar=True,
            ),
            score_cols=["A_score", "B_score"],
            context=ctx,
        )
    assert any("Using score column: A_score" in msg for msg in ctx.infos)
    fig_single.clf()


@pytest.mark.asyncio
async def test_create_enrichment_spatial_multifeature_suffix_resolution_and_missing_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    adata.obs["B_score"] = np.linspace(1.0, 0.0, adata.n_obs)

    seen: list[str] = []
    monkeypatch.setattr(
        viz_enrich,
        "plot_spatial_feature",
        lambda _adata, feature, ax, params: seen.append(feature) or ax.scatter([0], [0], c=[1.0]),
    )

    fig = await viz_enrich._create_enrichment_spatial(
        adata,
        VisualizationParameters(
            plot_type="enrichment",
            subtype="spatial",
            feature=["A", "B"],
        ),
        score_cols=["A_score", "B_score"],
        context=None,
    )
    assert seen[:2] == ["A_score", "B_score"]
    fig.clf()

    with pytest.raises(DataNotFoundError, match="None of the specified scores found"):
        await viz_enrich._create_enrichment_spatial(
            adata,
            VisualizationParameters(
                plot_type="enrichment",
                subtype="spatial",
                feature=["X", "Y"],
            ),
            score_cols=["A_score", "B_score"],
            context=None,
        )


def test_create_enrichmap_spatial_uses_visium_spatial_library_id(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    if "library_id" in adata.obs.columns:
        del adata.obs["library_id"]
    adata.uns["spatial"] = {
        "V1_Human_Lymph_Node": {
            "images": {},
            "scalefactors": {"spot_diameter_fullres": 1.0},
        }
    }
    adata.obs["A_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    captured: dict[str, str] = {}

    class _PL:
        @staticmethod
        def spatial_enrichmap(*_args, **kwargs):
            captured["library_id"] = kwargs["library_id"]
            plt.figure()

    fake_em = ModuleType("enrichmap")
    fake_em.pl = _PL()
    monkeypatch.setitem(sys.modules, "enrichmap", fake_em)

    fig = viz_enrich._create_enrichmap_spatial(
        adata,
        VisualizationParameters(
            plot_type="enrichment",
            subtype="spatial_score",
            feature="A",
        ),
        score_cols=["A_score"],
    )

    assert captured["library_id"] == "V1_Human_Lymph_Node"
    fig.clf()


def test_create_enrichmap_spatial_routes_cross_and_wraps_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    sentinel = object()

    fake_em = ModuleType("enrichmap")
    monkeypatch.setitem(sys.modules, "enrichmap", fake_em)

    monkeypatch.setattr(
        viz_enrich,
        "_create_enrichmap_cross_correlation",
        lambda *_a, **_k: sentinel,
    )
    out = viz_enrich._create_enrichmap_spatial(
        adata,
        VisualizationParameters(plot_type="enrichment", subtype="spatial_cross_correlation"),
        score_cols=["A_score"],
    )
    assert out is sentinel

    def _boom(*_args, **_kwargs):
        raise RuntimeError("bad")

    monkeypatch.setattr(viz_enrich, "_create_enrichmap_single_score", _boom)
    with pytest.raises(ProcessingError, match="EnrichMap spatial_score visualization failed"):
        viz_enrich._create_enrichmap_spatial(
            adata,
            VisualizationParameters(plot_type="enrichment", subtype="spatial_score", feature="A"),
            score_cols=["A_score"],
        )


def test_create_enrichmap_spatial_reraises_data_not_found(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    fake_em = ModuleType("enrichmap")
    monkeypatch.setitem(sys.modules, "enrichmap", fake_em)
    monkeypatch.setattr(
        viz_enrich,
        "_create_enrichmap_single_score",
        lambda *_a, **_k: (_ for _ in ()).throw(DataNotFoundError("missing score")),
    )
    with pytest.raises(DataNotFoundError, match="missing score"):
        viz_enrich._create_enrichmap_spatial(
            adata,
            VisualizationParameters(
                plot_type="enrichment",
                subtype="spatial_score",
                feature="A",
            ),
            score_cols=["A_score"],
        )


def test_create_enrichmap_cross_correlation_prefers_per_run_gene_sets(
    minimal_spatial_adata,
):
    """Per-run spatial gene-set key should be found as last-resort fallback;
    non-spatial keys must be ignored.  Priority order:
    1. Per-run key matching analysis_method hint
    2. Shared key (enrichment_spatial_gene_sets)
    3. Any per-run key (fallback)
    """
    adata = minimal_spatial_adata.copy()
    adata.obs["Wnt_score"] = 0.1
    adata.obs["Notch_score"] = 0.2

    # Only per-run spatial key present (no shared key)
    adata.uns["enrichment_spatial_gene_sets_spatial_KEGG_Pathways"] = {
        "Wnt": ["G1", "G2"],
        "Notch": ["G3", "G4"],
    }

    class _PL:
        @staticmethod
        def cross_moran_scatter(*_args, **_kwargs):
            plt.figure()

    class _EM:
        pl = _PL()

    params = VisualizationParameters(
        plot_type="enrichment",
        subtype="spatial_cross_correlation",
    )
    # Should NOT raise DataNotFoundError — per-run spatial key is found
    fig = viz_enrich._create_enrichmap_cross_correlation(
        adata, params, "sample_1", em=_EM()
    )
    fig.clf()

    # When shared key exists, it takes priority over per-run (no hint)
    adata.uns["enrichment_spatial_gene_sets"] = {
        "SharedA": ["G10"],
        "SharedB": ["G11"],
    }
    captured_pathways: list[str] = []

    class _PL2:
        @staticmethod
        def cross_moran_scatter(*_args, **_kwargs):
            captured_pathways.extend([_kwargs["score_x"], _kwargs["score_y"]])
            plt.figure()

    class _EM2:
        pl = _PL2()

    fig2 = viz_enrich._create_enrichmap_cross_correlation(
        adata, params, "sample_1", em=_EM2()
    )
    fig2.clf()
    # Shared key pathways (SharedA/SharedB) should be used, not KEGG ones
    assert captured_pathways[-2:] == ["SharedA_score", "SharedB_score"]

    # Non-spatial gene-set keys must NOT be picked up
    adata2 = minimal_spatial_adata.copy()
    adata2.uns["enrichment_gsea_gene_sets"] = {
        "PathA": ["G1"],
        "PathB": ["G2"],
    }
    adata2.uns["enrichment_ora_gene_sets"] = {
        "PathC": ["G3"],
        "PathD": ["G4"],
    }
    adata2.uns["enrichment_ssgsea_gene_sets"] = {
        "PathE": ["G5"],
        "PathF": ["G6"],
    }
    with pytest.raises(DataNotFoundError, match="Spatial enrichment gene sets not found"):
        viz_enrich._create_enrichmap_cross_correlation(
            adata2, params, "sample_1", em=_EM()
        )


def test_create_enrichmap_cross_correlation_analysis_method_selects_kegg(
    minimal_spatial_adata,
):
    """When analysis_method='KEGG', the KEGG per-run key should be selected
    over the Reactome per-run key."""
    adata = minimal_spatial_adata.copy()
    adata.obs["Wnt_score"] = 0.1
    adata.obs["Notch_score"] = 0.2
    adata.obs["PathA_score"] = 0.3
    adata.obs["PathB_score"] = 0.4

    # Two per-run keys: KEGG and Reactome
    adata.uns["enrichment_spatial_gene_sets_spatial_KEGG_Pathways"] = {
        "Wnt": ["G1", "G2"],
        "Notch": ["G3", "G4"],
    }
    adata.uns["enrichment_spatial_gene_sets_spatial_Reactome"] = {
        "PathA": ["G5", "G6"],
        "PathB": ["G7", "G8"],
    }

    captured_kwargs: list[dict] = []

    class _PL:
        @staticmethod
        def cross_moran_scatter(*_args, **_kwargs):
            captured_kwargs.append(_kwargs)
            plt.figure()

    class _EM:
        pl = _PL()

    # With analysis_method="KEGG", should pick KEGG key
    params_kegg = VisualizationParameters(
        plot_type="enrichment",
        subtype="spatial_cross_correlation",
        analysis_method="KEGG",
    )
    fig = viz_enrich._create_enrichmap_cross_correlation(
        adata, params_kegg, "sample_1", em=_EM()
    )
    fig.clf()
    # The pathways should come from the KEGG key (Wnt, Notch)
    assert captured_kwargs[-1]["score_x"] == "Wnt_score"
    assert captured_kwargs[-1]["score_y"] == "Notch_score"

    # With analysis_method="Reactome", should pick Reactome key
    params_reactome = VisualizationParameters(
        plot_type="enrichment",
        subtype="spatial_cross_correlation",
        analysis_method="Reactome",
    )
    fig2 = viz_enrich._create_enrichmap_cross_correlation(
        adata, params_reactome, "sample_1", em=_EM()
    )
    fig2.clf()
    assert captured_kwargs[-1]["score_x"] == "PathA_score"
    assert captured_kwargs[-1]["score_y"] == "PathB_score"


def test_create_enrichmap_cross_correlation_validation_and_success(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    params = VisualizationParameters(
        plot_type="enrichment",
        subtype="spatial_cross_correlation",
        figure_size=(8, 5),
        dpi=140,
    )

    with pytest.raises(DataNotFoundError, match="Spatial enrichment gene sets not found"):
        viz_enrich._create_enrichmap_cross_correlation(adata, params, "sample_1", em=object())

    adata.uns["enrichment_spatial_gene_sets"] = {"PathA": {"G1"}}
    with pytest.raises(DataNotFoundError, match="Need at least 2 pathways"):
        viz_enrich._create_enrichmap_cross_correlation(adata, params, "sample_1", em=object())

    adata.uns["enrichment_spatial_gene_sets"] = {"PathA": {"G1"}, "PathB": {"G2"}}

    class _PL:
        @staticmethod
        def cross_moran_scatter(*_args, **_kwargs):
            plt.figure()

    class _EM:
        pl = _PL()

    fig = viz_enrich._create_enrichmap_cross_correlation(adata, params, "sample_1", em=_EM())
    assert tuple(fig.get_size_inches()) == pytest.approx((8.0, 5.0))
    assert fig.get_dpi() == 140
    fig.clf()


def test_create_enrichmap_single_score_routes_all_supported_subtypes(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obs["A_score"] = np.linspace(0.0, 1.0, adata.n_obs)
    calls: list[str] = []

    class _PL:
        @staticmethod
        def morans_correlogram(*_args, **_kwargs):
            calls.append("correlogram")
            plt.figure()

        @staticmethod
        def variogram(*_args, **_kwargs):
            calls.append("variogram")
            plt.figure()

        @staticmethod
        def spatial_enrichmap(*_args, **_kwargs):
            calls.append("spatial_score")
            plt.figure()

    class _EM:
        pl = _PL()

    for subtype in ["spatial_correlogram", "spatial_variogram", "spatial_score"]:
        fig = viz_enrich._create_enrichmap_single_score(
            adata,
            VisualizationParameters(plot_type="enrichment", subtype=subtype, feature="A"),
            library_id="sample_1",
            em=_EM(),
            context=None,
        )
        fig.clf()

    assert calls == ["correlogram", "variogram", "spatial_score"]

    fig2 = viz_enrich._create_enrichmap_single_score(
        adata,
        VisualizationParameters(
            plot_type="enrichment",
            subtype="spatial_score",
            feature="A",
            figure_size=(7, 4),
            dpi=150,
        ),
        library_id="sample_1",
        em=_EM(),
        context=None,
    )
    assert tuple(fig2.get_size_inches()) == pytest.approx((7.0, 4.0))
    assert fig2.get_dpi() == 150
    fig2.clf()


def test_create_gsea_enrichment_plot_validation_and_success(monkeypatch: pytest.MonkeyPatch):
    with pytest.raises(DataNotFoundError, match="requires running enrichment scores"):
        viz_enrich._create_gsea_enrichment_plot(
            pd.DataFrame({"Term": ["PathA"]}),
            VisualizationParameters(plot_type="enrichment", subtype="enrichment_plot"),
        )

    with pytest.raises(DataNotFoundError, match="requires 'RES'"):
        viz_enrich._create_gsea_enrichment_plot(
            {"PathA": {"NES": 1.0}},
            VisualizationParameters(plot_type="enrichment", subtype="enrichment_plot"),
        )

    fake_gp = ModuleType("gseapy")
    fake_gp.gseaplot = lambda **_kwargs: plt.figure(figsize=(6, 4))
    monkeypatch.setitem(sys.modules, "gseapy", fake_gp)

    fig = viz_enrich._create_gsea_enrichment_plot(
        {"PathA": {"RES": [0.1, 0.2], "NES": 1.1, "pval": 0.02}},
        VisualizationParameters(plot_type="enrichment", subtype="enrichment_plot"),
    )
    assert fig is not None
    fig.clf()

    with pytest.raises(ParameterError, match="Unsupported GSEA results format"):
        viz_enrich._create_gsea_enrichment_plot(
            ["bad"],
            VisualizationParameters(plot_type="enrichment", subtype="enrichment_plot"),
        )

    fig_path = viz_enrich._create_gsea_enrichment_plot(
        {
            "PathA": {"RES": [0.1, 0.2], "NES": 1.1, "pval": 0.02},
            "PathB": {"RES": [0.2, 0.3], "NES": 1.3, "pval": 0.01},
        },
        VisualizationParameters(
            plot_type="enrichment",
            subtype="enrichment_plot",
            feature="PathB",
        ),
    )
    assert fig_path is not None
    fig_path.clf()


def test_create_gsea_dotplot_nested_and_error_wrap(monkeypatch: pytest.MonkeyPatch):
    fake_gp = ModuleType("gseapy")
    fake_gp.dotplot = lambda **_kwargs: _FakeAxes()
    monkeypatch.setitem(sys.modules, "gseapy", fake_gp)

    fig = viz_enrich._create_gsea_dotplot(
        {
            "ConditionA": {"PathA": {"Adjusted P-value": 0.01}},
            "ConditionB": {"PathB": {"Adjusted P-value": 0.02}},
        },
        VisualizationParameters(plot_type="enrichment", subtype="dotplot"),
    )
    assert fig is not None
    fig.clf()

    def _boom(**_kwargs):
        raise RuntimeError("dotplot failed")

    fake_gp.dotplot = _boom
    with pytest.raises(ProcessingError, match="gseapy.dotplot failed"):
        viz_enrich._create_gsea_dotplot(
            pd.DataFrame({"Term": ["PathA"], "Adjusted P-value": [0.1]}),
            VisualizationParameters(plot_type="enrichment", subtype="dotplot"),
        )

    with pytest.raises(DataNotFoundError, match="No enrichment results found"):
        viz_enrich._create_gsea_dotplot(
            {},
            VisualizationParameters(plot_type="enrichment", subtype="dotplot"),
        )


def test_create_gsea_barplot_success_empty_and_figure_size(monkeypatch: pytest.MonkeyPatch):
    captured: dict[str, object] = {}
    fake_gp = ModuleType("gseapy")

    def _fake_barplot(**kwargs):
        captured["figsize"] = kwargs.get("figsize")
        return _FakeAxes()

    fake_gp.barplot = _fake_barplot
    monkeypatch.setitem(sys.modules, "gseapy", fake_gp)

    fig = viz_enrich._create_gsea_barplot(
        {"PathA": {"Adjusted P-value": 0.01}},
        VisualizationParameters(
            plot_type="enrichment",
            subtype="barplot",
            figure_size=(9, 5),
        ),
    )
    assert captured["figsize"] == (9.0, 5.0)
    fig.clf()

    with pytest.raises(DataNotFoundError, match="No enrichment results found"):
        viz_enrich._create_gsea_barplot(
            {},
            VisualizationParameters(plot_type="enrichment", subtype="barplot"),
        )


def test_utility_branches_for_dataframe_conversion_and_feature_resolution():
    assert viz_enrich._resolve_feature_list(None, pd.Index(["A"]), ["A_score"]) == []
    assert viz_enrich._resolve_feature_list("A", pd.Index(["A"]), ["A_score"]) == ["A"]
    assert viz_enrich._resolve_feature_list(["A", "B"], pd.Index(["A"]), ["A_score"]) == [
        "A",
        "B",
    ]

    with pytest.raises(ParameterError, match="Unsupported GSEA results format"):
        viz_enrich._gsea_results_to_dataframe(["bad"])

    df, x_col = viz_enrich._nested_dict_to_dataframe(
        {
            "A": {"Path1": {"pval": 0.1}},
            "B": {"Path2": {"pval": 0.2}},
        }
    )
    assert x_col == "Group"
    assert set(df["Group"]) == {"A", "B"}

    assert viz_enrich._find_pvalue_column(pd.DataFrame({"x": [1]})) == "Adjusted P-value"

    df3 = pd.DataFrame({"score": [1, 2]}, index=pd.Index(["P1", "P2"], name="path"))
    viz_enrich._ensure_term_column(df3)
    assert list(df3["Term"]) == ["P1", "P2"]
