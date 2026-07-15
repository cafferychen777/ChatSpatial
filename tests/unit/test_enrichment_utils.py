"""Unit tests for enrichment helper and dispatch contracts."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

# Keep tests lightweight: if gseapy is unavailable in local env, inject a tiny stub.
if "gseapy" not in sys.modules:
    sys.modules["gseapy"] = SimpleNamespace()

from chatspatial.models.analysis import EnrichmentResult
from chatspatial.models.data import EnrichmentParameters
from chatspatial.tools import enrichment as enrichment_module
from chatspatial.tools.enrichment import (
    _build_enrichment_key,
    _compute_std_sparse_compatible,
    _convert_gene_format_for_matching,
    _filter_gene_sets_by_size,
    _filter_significant_statistics,
    _get_organism_name,
    analyze_enrichment,
    load_gene_sets,
    map_gene_set_database_to_enrichr_library,
)
from chatspatial.utils.exceptions import DependencyError, ParameterError, ProcessingError


class DummyCtx:
    def __init__(self, adata):
        self._adata = adata
        self.errors: list[str] = []
        self.set_adata_calls: list[tuple[str, object]] = []

    async def get_adata(self, data_id: str):
        return self._adata

    async def set_adata(self, data_id: str, adata):
        self.set_adata_calls.append((data_id, adata))
        self._adata = adata

    async def error(self, msg: str):
        self.errors.append(msg)

    async def warning(self, _msg: str):
        return None


def test_filter_significant_statistics_uses_gsea_threshold_025():
    stats, scores, pvals, adj = _filter_significant_statistics(
        gene_set_statistics={"A": {"x": 1}, "B": {"x": 2}},
        enrichment_scores={"A": 0.3, "B": 0.2},
        pvalues={"A": 0.01, "B": 0.2},
        adjusted_pvalues={"A": 0.24, "B": 0.26},
        method="gsea",
    )
    assert set(stats) == {"A"}
    assert set(scores) == {"A"}
    assert set(pvals) == {"A"}
    assert set(adj) == {"A"}


def test_map_gene_set_database_to_enrichr_library_maps_species_specific_kegg():
    human = map_gene_set_database_to_enrichr_library("KEGG_Pathways", "human")
    mouse = map_gene_set_database_to_enrichr_library("KEGG_Pathways", "mouse")
    assert human == "KEGG_2021_Human"
    assert mouse == "KEGG_2019_Mouse"


def test_load_gene_sets_unknown_database_raises_parameter_error():
    with pytest.raises(ParameterError, match="Unknown database"):
        load_gene_sets("NOT_A_DB")


def test_filter_gene_sets_by_size_keeps_only_in_range():
    gene_sets = {"small": ["A"], "ok": ["A", "B", "C"], "large": list("ABCDEFG")}
    out = _filter_gene_sets_by_size(gene_sets, min_size=2, max_size=4)
    assert out == {"ok": ["A", "B", "C"]}


def test_compute_std_sparse_compatible_matches_dense():
    X = np.array([[1.0, 0.0, 3.0], [2.0, 1.0, 1.0], [3.0, 1.0, 2.0]])
    X_sp = sparse.csr_matrix(X)

    dense_std = _compute_std_sparse_compatible(X, axis=0, ddof=1)
    sparse_std = _compute_std_sparse_compatible(X_sp, axis=0, ddof=1)
    assert np.allclose(dense_std, sparse_std)


def test_convert_gene_format_for_matching_mouse_and_human_rules():
    # Mouse: uppercase pathway genes should map to title/lowercase dataset genes
    mouse_genes, mouse_map = _convert_gene_format_for_matching(
        pathway_genes=["CD5L", "GM42418"],
        dataset_genes={"Cd5l", "Gm42418"},
        species="mouse",
    )
    assert set(mouse_genes) == {"Cd5l", "Gm42418"}
    assert mouse_map["Cd5l"] == "CD5L"

    # Human: lowercase pathway genes should map to uppercase dataset genes
    human_genes, human_map = _convert_gene_format_for_matching(
        pathway_genes=["hes1"],
        dataset_genes={"HES1"},
        species="human",
    )
    assert human_genes == ["HES1"]
    assert human_map["HES1"] == "hes1"


def test_get_organism_name_maps_species():
    assert _get_organism_name("human") == "Homo sapiens"
    assert _get_organism_name("mouse") == "Mus musculus"


@pytest.mark.asyncio
async def test_analyze_enrichment_pathway_ora_dispatch_with_loaded_gene_sets(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    async def _fake_get_adata(_):
        return minimal_spatial_adata

    ctx = DummyCtx(minimal_spatial_adata)
    ctx.get_adata = _fake_get_adata  # type: ignore[method-assign]

    monkeypatch.setattr(
        enrichment_module,
        "load_gene_sets",
        lambda **kwargs: {"set_a": ["gene_0", "gene_1"]},
    )

    def fake_perform_ora(**kwargs):
        assert kwargs["gene_sets"] == {"set_a": ["gene_0", "gene_1"]}
        return EnrichmentResult(
            method="pathway_ora",
            n_gene_sets=1,
            n_significant=1,
            top_gene_sets=["set_a"],
            top_depleted_sets=[],
        )

    monkeypatch.setattr(enrichment_module, "perform_ora", fake_perform_ora)

    params = EnrichmentParameters(species="human", method="pathway_ora")
    result = await analyze_enrichment("d1", ctx, params)

    assert result.method == "pathway_ora"
    assert result.n_significant == 1
    assert ctx.set_adata_calls == [("d1", minimal_spatial_adata)]


@pytest.mark.asyncio
async def test_analyze_enrichment_wraps_gene_set_loading_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx(minimal_spatial_adata)

    def boom_loader(**kwargs):
        raise RuntimeError("network down")

    monkeypatch.setattr(enrichment_module, "load_gene_sets", boom_loader)

    params = EnrichmentParameters(species="human", method="pathway_ora")
    with pytest.raises(ProcessingError, match="Failed to load gene sets"):
        await analyze_enrichment("d1", ctx, params)

    assert any("Gene set database loading failed" in msg for msg in ctx.errors)


@pytest.mark.asyncio
async def test_analyze_enrichment_requires_params(minimal_spatial_adata):
    ctx = DummyCtx(minimal_spatial_adata)
    with pytest.raises(ParameterError, match="params parameter is required"):
        await analyze_enrichment("d1", ctx, None)


@pytest.mark.asyncio
async def test_analyze_enrichment_rejects_empty_gene_sets_when_no_database(
    minimal_spatial_adata,
):
    ctx = DummyCtx(minimal_spatial_adata)
    params = EnrichmentParameters(
        species="human",
        method="pathway_ora",
        gene_sets={},
        gene_set_database=None,
    )
    with pytest.raises(ProcessingError, match="No valid gene sets available"):
        await analyze_enrichment("d1", ctx, params)


@pytest.mark.asyncio
async def test_analyze_enrichment_normalizes_list_gene_sets_for_ora(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx(minimal_spatial_adata)
    captured: dict[str, object] = {}

    def _fake_ora(**kwargs):
        captured["gene_sets"] = kwargs["gene_sets"]
        return EnrichmentResult(
            method="pathway_ora",
            n_gene_sets=1,
            n_significant=0,
            top_gene_sets=[],
            top_depleted_sets=[],
        )

    monkeypatch.setattr(enrichment_module, "perform_ora", _fake_ora)
    params = EnrichmentParameters(
        species="human",
        method="pathway_ora",
        gene_sets=["gene_0", "gene_1"],
        gene_set_database=None,
    )
    out = await analyze_enrichment("d1", ctx, params)
    assert out.method == "pathway_ora"
    assert captured["gene_sets"] == {"user_genes": ["gene_0", "gene_1"]}


@pytest.mark.asyncio
async def test_analyze_enrichment_limits_database_spatial_enrichmap_to_best_overlaps(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx(minimal_spatial_adata)
    captured: dict[str, object] = {}
    many_gene_sets = {
        f"set_{i:02d}": [f"gene_{i % 24}", f"gene_{(i + 1) % 24}"]
        for i in range(60)
    }
    many_gene_sets["best_match"] = ["gene_0", "gene_1", "gene_2", "gene_3"]

    monkeypatch.setattr(enrichment_module, "load_gene_sets", lambda **_kwargs: many_gene_sets)

    async def _fake_spatial(**kwargs):
        captured.update(kwargs)
        return EnrichmentResult(
            method="spatial_enrichmap",
            n_gene_sets=len(kwargs["gene_sets"]),
            n_significant=0,
            top_gene_sets=[],
            top_depleted_sets=[],
        )

    monkeypatch.setattr(enrichment_module, "perform_spatial_enrichment", _fake_spatial)

    params = EnrichmentParameters(
        species="human",
        method="spatial_enrichmap",
        gene_set_database="GO_Biological_Process",
    )
    await analyze_enrichment("d1", ctx, params)

    assert len(captured["gene_sets"]) == 50
    assert "best_match" in captured["gene_sets"]


@pytest.mark.asyncio
async def test_analyze_enrichment_uses_first_score_key_for_gsea(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx(minimal_spatial_adata)
    captured: dict[str, object] = {}

    def _fake_gsea(**kwargs):
        captured["ranking_key"] = kwargs["ranking_key"]
        return EnrichmentResult(
            method="pathway_gsea",
            n_gene_sets=1,
            n_significant=1,
            top_gene_sets=["set_a"],
            top_depleted_sets=[],
        )

    monkeypatch.setattr(enrichment_module, "perform_gsea", _fake_gsea)
    minimal_spatial_adata.var["score_a"] = np.linspace(1.0, -1.0, minimal_spatial_adata.n_vars)
    params = EnrichmentParameters(
        species="human",
        method="pathway_gsea",
        gene_sets={"set_a": ["gene_0", "gene_1"]},
        gene_set_database=None,
        score_keys=["score_a", "score_b"],
    )
    out = await analyze_enrichment("d1", ctx, params)
    assert out.method == "pathway_gsea"
    assert captured["ranking_key"] == "score_a"


@pytest.mark.asyncio
async def test_analyze_enrichment_rejects_missing_gsea_score_key(
    minimal_spatial_adata,
):
    ctx = DummyCtx(minimal_spatial_adata)
    params = EnrichmentParameters(
        species="human",
        method="pathway_gsea",
        gene_sets={"set_a": ["gene_0", "gene_1"]},
        gene_set_database=None,
        score_keys="spatial_domains_spagcn_n7_refined",
    )

    with pytest.raises(ParameterError, match="numeric gene-level ranking"):
        await analyze_enrichment("d1", ctx, params)


@pytest.mark.asyncio
async def test_analyze_enrichment_rejects_nonnumeric_gsea_score_key(
    minimal_spatial_adata,
):
    ctx = DummyCtx(minimal_spatial_adata)
    minimal_spatial_adata.var["rank_label"] = ["high"] * minimal_spatial_adata.n_vars
    params = EnrichmentParameters(
        species="human",
        method="pathway_gsea",
        gene_sets={"set_a": ["gene_0", "gene_1"]},
        gene_set_database=None,
        score_keys="rank_label",
    )

    with pytest.raises(ParameterError, match="numeric gene-level ranking"):
        await analyze_enrichment("d1", ctx, params)


@pytest.mark.asyncio
async def test_analyze_enrichment_enrichr_uses_highly_variable_genes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx(minimal_spatial_adata)
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.get_highly_variable_genes",
        lambda _adata, max_genes=500: ["gene_2", "gene_3"],
    )

    def _fake_enrichr(**kwargs):
        captured.update(kwargs)
        return EnrichmentResult(
            method="pathway_enrichr",
            n_gene_sets=2,
            n_significant=1,
            top_gene_sets=["KEGG_A"],
            top_depleted_sets=[],
        )

    monkeypatch.setattr(enrichment_module, "perform_enrichr", _fake_enrichr)
    params = EnrichmentParameters(
        species="human",
        method="pathway_enrichr",
        gene_set_database="KEGG_Pathways",
        gene_sets={"dummy_set": ["gene_0"]},
    )
    out = await analyze_enrichment("d1", ctx, params)
    assert out.method == "pathway_enrichr"
    assert captured["gene_list"] == ["gene_2", "gene_3"]
    assert captured["gene_sets"] == "KEGG_Pathways"
    assert captured["organism"] == "human"
    assert captured["adata"] is minimal_spatial_adata
    assert captured["species"] == "human"
    assert captured["database"] == "KEGG_Pathways"
    assert captured["data_id"] == "d1"


@pytest.mark.asyncio
async def test_analyze_enrichment_unknown_method_raises_parameter_error(
    minimal_spatial_adata,
):
    ctx = DummyCtx(minimal_spatial_adata)
    params = EnrichmentParameters(
        species="human",
        method="pathway_ora",
        gene_sets={"set_a": ["gene_0"]},
        gene_set_database=None,
    ).model_copy(update={"method": "not_real"})
    with pytest.raises(ParameterError, match="Unknown enrichment method"):
        await analyze_enrichment("d1", ctx, params)


@pytest.mark.asyncio
async def test_analyze_enrichment_preserves_database_dependency_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx(minimal_spatial_adata)
    params = EnrichmentParameters(
        species="human",
        method="pathway_ora",
        gene_sets=None,
        gene_set_database="GO_Biological_Process",
    )

    def _missing_backend(*_args, **_kwargs):
        raise DependencyError("gseapy is required")

    monkeypatch.setattr(enrichment_module, "load_gene_sets", _missing_backend)

    with pytest.raises(DependencyError, match="gseapy is required"):
        await analyze_enrichment("d1", ctx, params)


@pytest.mark.asyncio
async def test_perform_spatial_enrichment_partial_failure_still_returns_success(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class CtxWithLogs(DummyCtx):
        def __init__(self, ad):
            super().__init__(ad)
            self.warnings: list[str] = []
            self.infos: list[str] = []

        async def warning(self, msg: str):
            self.warnings.append(msg)

        async def info(self, msg: str):
            self.infos.append(msg)

    ctx = CtxWithLogs(adata)

    fake_enrichmap = SimpleNamespace()

    def _score(*, adata, gene_set, score_key, **_kwargs):
        if score_key == "sig_bad":
            raise RuntimeError("simulated failure")
        adata.obs[f"{score_key}_score"] = np.linspace(0.0, 1.0, adata.n_obs)

    fake_enrichmap.tl = SimpleNamespace(score=_score)

    captured: dict[str, object] = {}
    monkeypatch.setattr(
        enrichment_module, "require", lambda *_args, **_kwargs: fake_enrichmap
    )
    monkeypatch.setitem(__import__("sys").modules, "enrichmap", fake_enrichmap)
    monkeypatch.setattr(
        "chatspatial.tools.enrichment.store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        "chatspatial.tools.enrichment.export_analysis_result",
        lambda *_args, **_kwargs: [],
    )

    result = await enrichment_module.perform_spatial_enrichment(
        data_id="d1",
        ctx=ctx,
        gene_sets={"sig_ok": ["gene_0", "gene_1"], "sig_bad": ["gene_2", "gene_3"]},
        species="human",
        database="KEGG_Pathways",
    )

    assert result.method == "spatial_enrichmap"
    assert result.n_gene_sets == 1
    # Spatial enrichment has no significance testing — n_significant is always 0
    assert result.n_significant == 0
    assert result.n_successful_signatures == 1
    assert "sig_ok" in result.enrichment_scores
    assert "sig_ok_score" in adata.obs.columns
    # Gene sets keyed by method to avoid cross-analysis overwrite
    assert "enrichment_spatial_gene_sets" in adata.uns
    assert list(adata.uns["enrichment_spatial_gene_sets"].keys()) == ["sig_ok"]
    # Parametrized: includes database when provided
    assert captured["analysis_name"] == "enrichment_spatial_KEGG_Pathways"
    assert captured["results_keys"]["obs"] == ["sig_ok_score"]
    assert any("Failed to process 1 gene sets" in msg for msg in ctx.warnings)


@pytest.mark.asyncio
async def test_perform_spatial_enrichment_raises_when_all_signatures_fail(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()

    class CtxWithLogs(DummyCtx):
        async def warning(self, _msg: str):
            return None

        async def info(self, _msg: str):
            return None

    ctx = CtxWithLogs(adata)

    fake_enrichmap = SimpleNamespace(
        tl=SimpleNamespace(score=lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("boom")))
    )
    monkeypatch.setattr(
        enrichment_module, "require", lambda *_args, **_kwargs: fake_enrichmap
    )
    monkeypatch.setitem(__import__("sys").modules, "enrichmap", fake_enrichmap)

    with pytest.raises(ProcessingError, match="All EnrichMap scoring failed"):
        await enrichment_module.perform_spatial_enrichment(
            data_id="d1",
            ctx=ctx,
            gene_sets={"sig1": ["gene_0", "gene_1"], "sig2": ["gene_2", "gene_3"]},
        )


@pytest.mark.asyncio
async def test_perform_spatial_enrichment_requires_spatial_coordinates(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]

    class CtxWithLogs(DummyCtx):
        async def warning(self, _msg: str):
            return None

        async def info(self, _msg: str):
            return None

    ctx = CtxWithLogs(adata)

    fake_enrichmap = SimpleNamespace(tl=SimpleNamespace(score=lambda **_kwargs: None))
    monkeypatch.setattr(
        enrichment_module, "require", lambda *_args, **_kwargs: fake_enrichmap
    )
    monkeypatch.setitem(__import__("sys").modules, "enrichmap", fake_enrichmap)

    with pytest.raises(ProcessingError, match="Spatial coordinates 'spatial' not found"):
        await enrichment_module.perform_spatial_enrichment(
            data_id="d1",
            ctx=ctx,
            gene_sets={"sig1": ["gene_0", "gene_1"]},
        )


def test_perform_enrichr_maps_library_and_filters_significant(monkeypatch: pytest.MonkeyPatch):
    captured: dict[str, object] = {}

    class _EnrResult:
        def __init__(self):
            self.results = pd.DataFrame(
                {
                    "Term": ["Path_B", "Path_A"],
                    "Combined Score": [8.0, 5.0],
                    "P-value": [0.2, 0.001],
                    "Adjusted P-value": [0.2, 0.01],
                    "Z-score": [1.0, 2.5],
                    "Overlap": ["1/10", "2/10"],
                    "Genes": ["GENE3", "GENE1;GENE2"],
                    "Odds Ratio": [1.2, 2.0],
                }
            )

    def _fake_enrichr(**kwargs):
        captured.update(kwargs)
        return _EnrResult()

    monkeypatch.setattr(enrichment_module.gp, "enrichr", _fake_enrichr, raising=False)

    out = enrichment_module.perform_enrichr(
        gene_list=["GENE1", "GENE2"],
        gene_sets="KEGG_Pathways",
        organism="mouse",
    )

    assert captured["gene_sets"] == ["KEGG_2019_Mouse"]
    assert captured["organism"] == "Mouse"
    assert out.method == "enrichr"
    assert out.n_gene_sets == 2
    assert out.n_significant == 1
    assert set(out.gene_set_statistics.keys()) == {"Path_A"}
    assert out.top_gene_sets == ["Path_B", "Path_A"]


def test_perform_enrichr_stores_results_for_visualization(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    class _EnrResult:
        def __init__(self):
            self.results = pd.DataFrame(
                {
                    "Term": ["Path_B", "Path_A"],
                    "Combined Score": [8.0, 5.0],
                    "P-value": [0.2, 0.001],
                    "Adjusted P-value": [0.2, 0.01],
                    "Z-score": [1.0, 2.5],
                    "Overlap": ["1/10", "2/10"],
                    "Genes": ["GENE3", "GENE1;GENE2"],
                    "Odds Ratio": [1.2, 2.0],
                }
            )

    monkeypatch.setattr(
        enrichment_module.gp,
        "enrichr",
        lambda **_kwargs: _EnrResult(),
        raising=False,
    )

    out = enrichment_module.perform_enrichr(
        gene_list=["GENE1", "GENE2"],
        gene_sets="MSigDB_Hallmark",
        organism="human",
        adata=minimal_spatial_adata,
        species="human",
        database="MSigDB_Hallmark",
        data_id=None,
    )

    assert out.method == "enrichr"
    assert "enrichr_results" in minimal_spatial_adata.uns
    assert "enrichr_results_enrichr_MSigDB_Hallmark" in minimal_spatial_adata.uns
    assert "enrichment_enrichr_MSigDB_Hallmark_metadata" in minimal_spatial_adata.uns
    assert (
        minimal_spatial_adata.uns[enrichment_module._LATEST_ENRICHMENT_RESULTS_KEY]
        == "enrichr_results"
    )


def test_perform_gsea_stores_latest_result_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    minimal_spatial_adata.var["rank_score"] = np.linspace(
        1.0, -1.0, minimal_spatial_adata.n_vars
    )

    class _PreRankResult:
        res2d = pd.DataFrame(
            {
                "Term": ["Path_A"],
                "ES": [0.5],
                "NES": [1.2],
                "NOM p-val": [0.01],
                "FDR q-val": [0.05],
                "Matched_size": [2],
                "Lead_genes": ["gene_0;gene_1"],
            }
        )

    monkeypatch.setattr(
        enrichment_module.gp,
        "prerank",
        lambda **_kwargs: _PreRankResult(),
        raising=False,
    )

    enrichment_module.perform_gsea(
        adata=minimal_spatial_adata,
        gene_sets={"Path_A": ["gene_0", "gene_1"]},
        ranking_key="rank_score",
        database="MSigDB_Hallmark",
        data_id=None,
    )

    assert (
        minimal_spatial_adata.uns[enrichment_module._LATEST_ENRICHMENT_RESULTS_KEY]
        == "gsea_results"
    )


def test_perform_ora_stores_latest_result_key(minimal_spatial_adata):
    minimal_spatial_adata.var["highly_variable"] = True
    enrichment_module.perform_ora(
        adata=minimal_spatial_adata,
        gene_sets={"Path_A": ["gene_0", "gene_1"], "Path_B": ["missing"]},
        database="MSigDB_Hallmark",
        data_id=None,
    )

    assert (
        minimal_spatial_adata.uns[enrichment_module._LATEST_ENRICHMENT_RESULTS_KEY]
        == "ora_results"
    )


def test_load_gene_sets_dispatches_to_expected_loader(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(enrichment_module, "load_go_gene_sets", lambda *args, **kwargs: {"go": ["G1", "G2"]})
    monkeypatch.setattr(enrichment_module, "load_kegg_gene_sets", lambda *args, **kwargs: {"kegg": ["G1", "G2"]})

    out_go = load_gene_sets("GO_Biological_Process", species="human")
    out_kegg = load_gene_sets("KEGG_Pathways", species="mouse")

    assert out_go == {"go": ["G1", "G2"]}
    assert out_kegg == {"kegg": ["G1", "G2"]}


def test_load_go_gene_sets_invalid_aspect_raises_parameter_error():
    with pytest.raises(ParameterError, match="Invalid GO aspect"):
        enrichment_module.load_go_gene_sets("human", aspect="BAD")


@pytest.mark.asyncio
async def test_analyze_enrichment_spatial_dispatches_to_spatial_function(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx(minimal_spatial_adata)
    captured: dict[str, object] = {}

    async def _fake_spatial(**kwargs):
        captured.update(kwargs)
        return EnrichmentResult(
            method="spatial_enrichmap",
            n_gene_sets=1,
            n_significant=1,
            top_gene_sets=["sig_a"],
            top_depleted_sets=[],
        )

    monkeypatch.setattr(enrichment_module, "perform_spatial_enrichment", _fake_spatial)

    out = await analyze_enrichment(
        "d1",
        ctx,
        EnrichmentParameters(
            species="human",
            method="spatial_enrichmap",
            gene_sets={"sig_a": ["gene_0", "gene_1"]},
            score_keys=["sig_a"],
            spatial_key="spatial",
            n_neighbors=8,
            smoothing=True,
            correct_spatial_covariates=True,
        ),
    )

    assert out.method == "spatial_enrichmap"
    assert captured["data_id"] == "d1"
    assert captured["gene_sets"] == {"sig_a": ["gene_0", "gene_1"]}
    assert captured["score_keys"] == ["sig_a"]


def test_perform_enrichr_uses_default_libraries_and_handles_missing_optional_columns(
    monkeypatch: pytest.MonkeyPatch,
):
    captured: dict[str, object] = {}

    class _EnrResult:
        def __init__(self):
            # No Z-score or Odds Ratio columns to exercise defaults.
            self.results = pd.DataFrame(
                {
                    "Term": ["Term_A"],
                    "Combined Score": [2.5],
                    "P-value": [0.001],
                    "Adjusted P-value": [0.01],
                    "Overlap": ["1/10"],
                    "Genes": ["G1;G2"],
                }
            )

    def _fake_enrichr(**kwargs):
        captured.update(kwargs)
        return _EnrResult()

    monkeypatch.setattr(enrichment_module.gp, "enrichr", _fake_enrichr, raising=False)

    out = enrichment_module.perform_enrichr(
        gene_list=["G1", "G2"],
        gene_sets=None,
        organism="human",
    )

    assert isinstance(captured["gene_sets"], list)
    assert "GO_Biological_Process_2025" in captured["gene_sets"]
    assert out.n_gene_sets == 1
    assert out.gene_set_statistics["Term_A"]["z_score"] != float("inf")
    assert out.gene_set_statistics["Term_A"]["odds_ratio"] == 1.0


@pytest.mark.asyncio
async def test_analyze_enrichment_dispatches_to_ssgsea(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ctx = DummyCtx(minimal_spatial_adata)
    captured: dict[str, object] = {}

    def _fake_ssgsea(**kwargs):
        captured.update(kwargs)
        return EnrichmentResult(
            method="pathway_ssgsea",
            n_gene_sets=1,
            n_significant=0,
            top_gene_sets=["set_a"],
            top_depleted_sets=[],
        )

    monkeypatch.setattr(enrichment_module, "perform_ssgsea", _fake_ssgsea)

    out = await analyze_enrichment(
        "d1",
        ctx,
        EnrichmentParameters(
            species="human",
            method="pathway_ssgsea",
            gene_sets={"set_a": ["gene_0", "gene_1"]},
            gene_set_database=None,
            min_genes=2,
            max_genes=100,
        ),
    )

    assert out.method == "pathway_ssgsea"
    assert captured["gene_sets"] == {"set_a": ["gene_0", "gene_1"]}
    assert captured["min_size"] == 2
    assert captured["max_size"] == 100


def test_load_msigdb_hallmark_filters_by_size(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        enrichment_module.gp,
        "get_library",
        lambda name, organism: {
            "small": ["A"],
            "ok": ["A", "B", "C"],
            "large": [str(i) for i in range(10)],
        },
        raising=False,
    )

    out = enrichment_module.load_msigdb_gene_sets(
        species="human",
        collection="H",
        min_size=2,
        max_size=5,
    )
    assert out == {"ok": ["A", "B", "C"]}


def test_load_library_falls_back_when_gseapy_leaves_lines_as_bytes(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(
        enrichment_module.gp,
        "get_library",
        lambda *_a, **_k: (_ for _ in ()).throw(
            TypeError("a bytes-like object is required, not 'str'")
        ),
        raising=False,
    )
    monkeypatch.setattr(
        enrichment_module,
        "_download_enrichr_library",
        lambda name, organism: {b"PathA": [b"Gene1", "Gene2"]},
    )

    out = enrichment_module._load_library_first_available(
        "MSigDB_Hallmark", "Mus musculus"
    )

    assert out == {"PathA": ["Gene1", "Gene2"]}


def test_load_go_gene_sets_calls_gseapy_with_expected_library(monkeypatch: pytest.MonkeyPatch):
    captured: dict[str, object] = {}

    def _fake_get_library(name, organism):
        captured["name"] = name
        captured["organism"] = organism
        return {"go_ok": ["A", "B", "C"]}

    monkeypatch.setattr(enrichment_module.gp, "get_library", _fake_get_library, raising=False)

    out = enrichment_module.load_go_gene_sets("mouse", aspect="BP", min_size=2, max_size=10)

    assert out == {"go_ok": ["A", "B", "C"]}
    assert captured["name"] == "GO_Biological_Process_2025"
    assert captured["organism"] == "Mus musculus"


def test_load_kegg_gene_sets_uses_species_specific_library(monkeypatch: pytest.MonkeyPatch):
    calls: list[tuple[str, str]] = []

    def _fake_get_library(name, organism):
        calls.append((name, organism))
        return {"k": ["A", "B", "C"]}

    monkeypatch.setattr(enrichment_module.gp, "get_library", _fake_get_library, raising=False)

    out_h = enrichment_module.load_kegg_gene_sets("human", min_size=2, max_size=10)
    out_m = enrichment_module.load_kegg_gene_sets("mouse", min_size=2, max_size=10)

    assert out_h == {"k": ["A", "B", "C"]}
    assert out_m == {"k": ["A", "B", "C"]}
    assert calls[0][0] == "KEGG_2021_Human"
    assert calls[1][0] == "KEGG_2019_Mouse"


def test_load_reactome_and_cell_marker_gene_sets(monkeypatch: pytest.MonkeyPatch):
    calls: list[tuple[str, str]] = []

    def _fake_get_library(name, organism):
        calls.append((name, organism))
        return {
            "ok": ["A", "B", "C"],
            "small": ["A"],
        }

    monkeypatch.setattr(enrichment_module.gp, "get_library", _fake_get_library, raising=False)

    reactome = enrichment_module.load_reactome_gene_sets("human", min_size=2, max_size=10)
    cellm = enrichment_module.load_cell_marker_gene_sets("mouse", min_size=2, max_size=10)

    assert reactome == {"ok": ["A", "B", "C"]}
    assert cellm == {"ok": ["A", "B", "C"]}
    assert calls[0] == ("Reactome_Pathways_2024", "Homo sapiens")
    assert calls[1] == ("CellMarker_Augmented_2021", "Mus musculus")


def test_perform_gsea_with_ranking_key_persists_results(monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.var["rank_metric"] = np.linspace(1.0, 0.1, adata.n_vars)

    class _Res:
        def __init__(self):
            self.res2d = pd.DataFrame(
                {
                    "Term": ["GS_A", "GS_B"],
                    "ES": [0.5, -0.2],
                    "NES": [1.7, -1.1],
                    "NOM p-val": [0.01, 0.2],
                    "FDR q-val": [0.2, 0.3],
                    "Matched_size": [12, 10],
                    "Lead_genes": ["gene_0;gene_1", "gene_2"],
                }
            )

    captured: dict[str, object] = {}

    monkeypatch.setattr(enrichment_module.gp, "prerank", lambda **kwargs: _Res(), raising=False)
    monkeypatch.setattr(
        enrichment_module,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(enrichment_module, "export_analysis_result", lambda *_a, **_k: None)

    out = enrichment_module.perform_gsea(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1"], "GS_B": ["gene_2", "gene_3"]},
        ranking_key="rank_metric",
        data_id="d1",
        species="human",
        database="KEGG_Pathways",
    )

    assert out.method == "gsea"
    assert out.n_gene_sets == 2
    # GSEA uses FDR < 0.25 threshold (Subramanian et al. 2005); GS_A has FDR=0.2
    assert out.n_significant == 1
    assert set(out.gene_set_statistics) == {"GS_A"}
    assert "gsea_results" in adata.uns
    assert "enrichment_gsea_gene_sets" in adata.uns
    # Parametrized: includes database when provided
    assert captured["analysis_name"] == "enrichment_gsea_KEGG_Pathways"
    # Per-run results key coexists with shared key
    assert "gsea_results_gsea_KEGG_Pathways" in adata.uns


def test_perform_ssgsea_success_populates_obs_and_uns(monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()

    class _Res:
        def __init__(self, obs_names):
            rows = []
            for i, sample in enumerate(obs_names):
                rows.extend(
                    [
                        {"Name": sample, "Term": "GS_A", "ES": 0.4 + (i % 3) * 0.1},
                        {"Name": sample, "Term": "GS_B", "ES": 0.2},
                    ]
                )
            self.res2d = pd.DataFrame(rows)
            self.results = {
                sample: {
                    "GS_A": {"es": 0.4 + (i % 3) * 0.1},
                    "GS_B": {"es": 0.2},
                }
                for i, sample in enumerate(obs_names)
            }

    captured: dict[str, object] = {}

    monkeypatch.setattr(
        enrichment_module.gp,
        "ssgsea",
        lambda **kwargs: _Res(adata.obs_names),
        raising=False,
    )
    monkeypatch.setattr(
        enrichment_module,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(enrichment_module, "export_analysis_result", lambda *_a, **_k: None)

    out = enrichment_module.perform_ssgsea(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1"], "GS_B": ["gene_2", "gene_3"]},
        data_id="d1",
        species="mouse",
    )

    assert out.method == "ssgsea"
    assert out.n_gene_sets == 2
    assert out.n_significant == 0
    assert set(out.top_gene_sets) == {"GS_A", "GS_B"}
    assert "ssgsea_GS_A" in adata.obs.columns
    assert "ssgsea_GS_B" in adata.obs.columns
    assert "enrichment_ssgsea_gene_sets" in adata.uns
    assert adata.uns["enrichment_latest_score_key"] == "enrichment_ssgsea"
    # Parametrized: no database → falls back to method-only key
    assert captured["analysis_name"] == "enrichment_ssgsea"

    # With database, the key includes database name
    # (tested separately in test_enrichment_branch_contracts.py)


def test_perform_ssgsea_invalid_result_format_raises_processing_error(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
):
    adata = minimal_spatial_adata.copy()

    class _BadRes:
        def __init__(self):
            self.results = "bad-format"

    monkeypatch.setattr(
        enrichment_module.gp,
        "ssgsea",
        lambda **kwargs: _BadRes(),
        raising=False,
    )

    with pytest.raises(ProcessingError, match="ssGSEA results format not recognized"):
        enrichment_module.perform_ssgsea(
            adata=adata,
            gene_sets={"GS_A": ["gene_0", "gene_1"]},
        )


def test_build_enrichment_key_with_database():
    key = _build_enrichment_key("gsea", "KEGG_Pathways")
    assert key == "enrichment_gsea_KEGG_Pathways"


def test_build_enrichment_key_without_database():
    key = _build_enrichment_key("ora", None)
    assert key == "enrichment_ora"


def test_build_enrichment_key_sanitizes_spaces_and_slashes():
    key = _build_enrichment_key("gsea", "GO Biological/Process")
    assert key == "enrichment_gsea_GO_Biological_Process"


def test_perform_gsea_without_database_uses_method_only_key(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
):
    """When database is None, analysis_name falls back to enrichment_gsea."""
    adata = minimal_spatial_adata.copy()
    adata.var["rank_metric"] = np.linspace(1.0, 0.1, adata.n_vars)

    class _Res:
        def __init__(self):
            self.res2d = pd.DataFrame(
                {
                    "Term": ["GS_A"],
                    "ES": [0.5],
                    "NES": [1.7],
                    "NOM p-val": [0.01],
                    "FDR q-val": [0.2],
                    "Matched_size": [12],
                    "Lead_genes": ["gene_0;gene_1"],
                }
            )

    captured: dict[str, object] = {}
    monkeypatch.setattr(
        enrichment_module.gp,
        "prerank",
        lambda **kwargs: _Res(),
        raising=False,
    )
    monkeypatch.setattr(
        enrichment_module,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        enrichment_module,
        "export_analysis_result",
        lambda *_a, **_k: None,
    )

    enrichment_module.perform_gsea(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1"]},
        ranking_key="rank_metric",
        data_id="d1",
        database=None,
    )

    assert captured["analysis_name"] == "enrichment_gsea"
    # Per-run key without database suffix
    assert "gsea_results_gsea" in adata.uns


def test_perform_ora_parametrized_key_with_database(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
):
    """ORA with database produces parametrized analysis key."""
    adata = minimal_spatial_adata.copy()
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        enrichment_module,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.update(kwargs),
    )
    monkeypatch.setattr(
        enrichment_module,
        "export_analysis_result",
        lambda *_a, **_k: None,
    )

    enrichment_module.perform_ora(
        adata=adata,
        gene_sets={"GS_A": adata.var_names[:3].tolist()},
        gene_list=adata.var_names[:3].tolist(),
        min_size=1,
        max_size=100,
        database="GO_Biological_Process",
        data_id="d1",
    )

    assert captured["analysis_name"] == "enrichment_ora_GO_Biological_Process"
    assert "ora_results_ora_GO_Biological_Process" in adata.uns


# ============================================================================
# Bug-fix regression tests
# ============================================================================


def test_check_is_integer_counts_default_sample_size_catches_low_frequency_decimals():
    """With sample_size=1000, a matrix with 1% decimals is reliably detected."""
    from chatspatial.utils.adata_utils import check_is_integer_counts

    rng = np.random.default_rng(42)
    # 10_000 integer values with 1% replaced by decimals
    data = rng.integers(0, 100, size=10_000).astype(float)
    n_decimal = int(len(data) * 0.01)
    decimal_indices = rng.choice(len(data), size=n_decimal, replace=False)
    data[decimal_indices] += 0.5

    mat = sparse.csr_matrix(data.reshape(100, 100))
    is_int, _has_neg, has_dec = check_is_integer_counts(mat)
    # Default sample_size=1000 should detect the decimals
    assert not is_int, "Expected non-integer detection with 1% decimals"
    assert has_dec, "Expected has_decimals=True"


def test_map_gene_set_database_to_enrichr_library_raises_for_zebrafish_kegg():
    """Zebrafish KEGG should raise ParameterError, not silently fall back to mouse."""
    with pytest.raises(ParameterError, match="KEGG_Pathways is only available"):
        map_gene_set_database_to_enrichr_library("KEGG_Pathways", "zebrafish")


def test_load_kegg_gene_sets_raises_for_unsupported_species():
    """load_kegg_gene_sets should raise for non-human/mouse species."""
    with pytest.raises(ParameterError, match="KEGG_Pathways is only available"):
        enrichment_module.load_kegg_gene_sets("zebrafish")


def test_perform_enrichr_passes_pvalue_cutoff_through(
    monkeypatch: pytest.MonkeyPatch,
):
    """perform_enrichr should forward pvalue_cutoff to gp.enrichr cutoff."""
    captured: dict[str, object] = {}

    class _EnrResult:
        def __init__(self):
            self.results = pd.DataFrame(
                {
                    "Term": ["Path_A", "Path_B"],
                    "Combined Score": [5.0, 8.0],
                    "P-value": [0.001, 0.001],
                    "Adjusted P-value": [0.005, 0.02],
                    "Z-score": [2.0, 1.5],
                    "Overlap": ["2/10", "1/10"],
                    "Genes": ["GENE1;GENE2", "GENE3"],
                    "Odds Ratio": [2.0, 1.5],
                }
            )

    def _fake_enrichr(**kwargs):
        captured.update(kwargs)
        return _EnrResult()

    monkeypatch.setattr(
        enrichment_module.gp, "enrichr", _fake_enrichr, raising=False
    )

    # Use a stricter cutoff of 0.01
    out = enrichment_module.perform_enrichr(
        gene_list=["GENE1", "GENE2"],
        gene_sets="GO_Biological_Process",
        organism="human",
        pvalue_cutoff=0.01,
    )

    assert captured["cutoff"] == 0.01, (
        "gp.enrichr should receive the user-specified cutoff"
    )
    # Path_A has adj-pval 0.005 < 0.01, Path_B has 0.02 >= 0.01
    assert out.n_significant == 1
