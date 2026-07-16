"""Focused branch/contract tests for enrichment analysis utilities."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from chatspatial.tools import enrichment as enrichment_module
from chatspatial.tools.enrichment import (
    _compute_variance_ranking,
    _convert_gene_format_for_matching,
    _filter_significant_statistics,
    _top_n_desc_indices,
    load_cell_marker_gene_sets,
    load_go_gene_sets,
    load_kegg_gene_sets,
    load_msigdb_gene_sets,
    load_reactome_gene_sets,
    map_gene_set_database_to_enrichr_library,
    perform_enrichr,
    perform_gsea,
    perform_ora,
    perform_spatial_enrichment,
    perform_ssgsea,
)
from chatspatial.utils.exceptions import DataNotFoundError, ParameterError, ProcessingError


def _patch_gseapy(monkeypatch: pytest.MonkeyPatch, **methods):
    backend = SimpleNamespace(**methods)
    monkeypatch.setattr(
        enrichment_module,
        "_get_gseapy",
        lambda *_args, **_kwargs: backend,
    )
    return backend


class _LogCtx:
    def __init__(self, adata):
        self._adata = adata
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.infos: list[str] = []

    async def get_adata(self, _data_id: str):
        return self._adata

    async def set_adata(self, _data_id: str, adata):
        self._adata = adata

    async def error(self, msg: str):
        self.errors.append(msg)

    async def warning(self, msg: str):
        self.warnings.append(msg)

    async def info(self, msg: str):
        self.infos.append(msg)


def _gsea_res_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Term": ["GS_A", "GS_B"],
            "ES": [0.5, -0.3],
            "NES": [1.8, -1.2],
            "NOM p-val": [0.01, 0.2],
            "FDR q-val": [0.2, 0.3],
        }
    )


def _patch_metadata_noop(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(enrichment_module, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(enrichment_module, "export_analysis_result", lambda *_a, **_k: None)


def test_filter_significant_statistics_no_pvalues_and_default_method_threshold() -> None:
    stats = {"A": {"x": 1}, "B": {"x": 2}}
    scores = {"A": 0.4, "B": 0.2}
    pvals = {"A": 0.01, "B": 0.2}

    out_no_p = _filter_significant_statistics(
        gene_set_statistics=stats,
        enrichment_scores=scores,
        pvalues=pvals,
        adjusted_pvalues={},
        method="ssgsea",
    )
    assert out_no_p == (stats, scores, pvals, {})

    out_default = _filter_significant_statistics(
        gene_set_statistics=stats,
        enrichment_scores=scores,
        pvalues=pvals,
        adjusted_pvalues={"A": 0.049, "B": 0.051},
        method="unknown_method",
    )
    assert out_default[0] == {"A": {"x": 1}}
    assert out_default[1] == {"A": 0.4}


def test_convert_gene_format_direct_match_keeps_original_mapping() -> None:
    genes, conv = _convert_gene_format_for_matching(
        pathway_genes=["Cd5l"],
        dataset_genes={"Cd5l"},
        species="mouse",
    )
    assert genes == ["Cd5l"]
    assert conv == {"Cd5l": "Cd5l"}


def test_map_gene_set_database_to_enrichr_library_rejects_unknown_option() -> None:
    with pytest.raises(ParameterError, match="Unknown gene set database"):
        map_gene_set_database_to_enrichr_library("NotARealDb", "human")


@pytest.mark.parametrize(
    "mode, ranking_method, expected_top",
    [
        ("binary", "signal_to_noise", None),
        ("multigroup", "signal_to_noise", None),
        ("hvg", "highly_variable_rank", "gene_0"),
        ("dispersion", "dispersions_norm", "gene_0"),
        ("cv", "coefficient_of_variation", None),
    ],
)
def test_perform_gsea_covers_ranking_strategy_branches(
    mode: str,
    ranking_method: str,
    expected_top: str | None,
    monkeypatch: pytest.MonkeyPatch,
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)

    # Use deterministic values so branch-specific ranking is stable.
    adata.X = np.asarray(adata.X, dtype=np.float64)
    adata.X[:, 0] = np.linspace(1.0, 100.0, adata.n_obs)
    adata.X[:, 1] = np.linspace(100.0, 1.0, adata.n_obs)

    if mode == "multigroup":
        adata.obs["group"] = ["A"] * 20 + ["B"] * 20 + ["C"] * 20
    elif mode in ("hvg", "dispersion", "cv"):
        if "group" in adata.obs:
            del adata.obs["group"]
        if mode == "hvg":
            adata.var["highly_variable_rank"] = np.linspace(10.0, 1.0, adata.n_vars)
        elif mode == "dispersion":
            adata.var["dispersions_norm"] = np.linspace(5.0, 0.1, adata.n_vars)

    monkeypatch.setattr(
        enrichment_module,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: SimpleNamespace(
            X=_adata.X, var_names=_adata.var_names
        ),
    )

    captured: dict[str, object] = {}

    class _Res:
        def __init__(self):
            self.res2d = _gsea_res_df()

    def _fake_prerank(**kwargs):
        captured["rnk"] = kwargs["rnk"]
        return _Res()

    _patch_gseapy(monkeypatch, prerank=_fake_prerank)

    out = perform_gsea(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1"], "GS_B": ["gene_2", "gene_3"]},
        ranking_key=None,
        method=ranking_method,
        permutation_num=10,
        min_size=1,
        max_size=1000,
        data_id="d1",
    )

    assert out.method == "gsea"
    assert out.n_gene_sets == 2
    rnk_df = captured["rnk"]
    assert isinstance(rnk_df, pd.DataFrame)
    assert rnk_df.shape[0] == adata.n_vars
    assert np.isfinite(rnk_df["score"]).all()

    if expected_top is not None:
        assert rnk_df.index[0] == expected_top


def test_perform_gsea_rejects_unknown_ranking_method(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)

    with pytest.raises(ParameterError, match="Unsupported GSEA ranking method"):
        perform_gsea(
            adata=adata,
            gene_sets={"GS_A": ["gene_0", "gene_1"]},
            method="not_a_real_method",
        )


def test_compute_variance_ranking_covers_sparse_and_dense_paths(
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    dense_ranking = _compute_variance_ranking(np.asarray(adata.X), adata.var_names)
    assert len(dense_ranking) == adata.n_vars

    from scipy import sparse as sp

    sparse_ranking = _compute_variance_ranking(
        sp.csr_matrix(np.asarray(adata.X)), adata.var_names
    )
    assert len(sparse_ranking) == adata.n_vars


def test_perform_gsea_uses_variance_ranking_mode(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)

    monkeypatch.setattr(
        enrichment_module,
        "get_raw_data_source",
        lambda _adata, prefer_complete_genes=True: SimpleNamespace(
            X=_adata.X, var_names=_adata.var_names
        ),
    )

    class _Res:
        def __init__(self):
            self.res2d = _gsea_res_df()

    _patch_gseapy(monkeypatch, prerank=lambda **_kwargs: _Res())

    out = perform_gsea(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1"]},
        method="variance",
        permutation_num=10,
    )
    assert out.method == "gsea"


@pytest.mark.parametrize("method_name", ["highly_variable_rank", "dispersions_norm"])
def test_perform_gsea_raises_when_requested_ranking_column_missing(
    method_name: str, minimal_spatial_adata
) -> None:
    with pytest.raises(DataNotFoundError, match="requested but adata.var"):
        perform_gsea(
            adata=minimal_spatial_adata.copy(),
            gene_sets={"GS_A": ["gene_0", "gene_1"]},
            method=method_name,
        )


def test_perform_gsea_propagates_prerank_failures(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
) -> None:
    adata = minimal_spatial_adata.copy()
    adata.var["rank_metric"] = np.linspace(1.0, 0.0, adata.n_vars)
    _patch_metadata_noop(monkeypatch)

    _patch_gseapy(
        monkeypatch,
        prerank=lambda **_kwargs: (_ for _ in ()).throw(
            RuntimeError("prerank failed")
        ),
    )
    with pytest.raises(RuntimeError, match="prerank failed"):
        perform_gsea(
            adata=adata,
            gene_sets={"GS": ["gene_0", "gene_1"]},
            ranking_key="rank_metric",
        )


def test_perform_ora_builds_gene_list_from_rank_genes_groups_without_pvals(
    monkeypatch: pytest.MonkeyPatch,
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)

    n = 120  # forces i >= 100 branch for no-pvals path
    names = np.zeros(n, dtype=[("A", "U16"), ("B", "U16")])
    for i in range(n):
        names["A"][i] = f"gene_{i % adata.n_vars}"
        names["B"][i] = f"gene_{(i + 1) % adata.n_vars}"
    adata.uns["rank_genes_groups"] = {"names": names}

    out = perform_ora(
        adata=adata,
        gene_sets={"GS_A": adata.var_names[:12].tolist()},
        min_size=1,
        max_size=200,
        data_id="d1",
    )

    assert out.method == "ora"
    assert out.n_gene_sets == 1


def test_perform_ora_uses_pvals_adj_filtering_when_available(
    monkeypatch: pytest.MonkeyPatch,
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)

    names = np.zeros(3, dtype=[("A", "U16")])
    names["A"] = ["gene_0", "gene_1", "gene_2"]
    pvals_adj = np.zeros(3, dtype=[("A", "f8")])
    pvals_adj["A"] = [0.01, 0.2, 0.03]
    adata.uns["rank_genes_groups"] = {"names": names, "pvals_adj": pvals_adj}

    out = perform_ora(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_2"], "GS_B": ["gene_1"]},
        gene_list=None,
        pvalue_threshold=0.05,
        min_size=1,
        max_size=100,
    )
    assert out.gene_set_statistics["GS_A"]["query_size"] == 2
    assert "GS_B" not in out.gene_set_statistics


def test_perform_ora_uses_pvals_filtering_when_adjusted_not_available(
    monkeypatch: pytest.MonkeyPatch,
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)

    names = np.zeros(3, dtype=[("A", "U16")])
    names["A"] = ["gene_0", "gene_1", "gene_2"]
    pvals = np.zeros(3, dtype=[("A", "f8")])
    pvals["A"] = [0.01, 0.04, 0.2]
    adata.uns["rank_genes_groups"] = {"names": names, "pvals": pvals}

    out = perform_ora(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1"], "GS_B": ["gene_2"]},
        gene_list=None,
        pvalue_threshold=0.05,
        min_size=1,
        max_size=100,
    )
    assert out.gene_set_statistics["GS_A"]["query_size"] == 2
    assert "GS_B" not in out.gene_set_statistics


def test_perform_ora_uses_highly_variable_genes_when_rank_genes_missing(
    monkeypatch: pytest.MonkeyPatch,
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)
    adata.var["highly_variable"] = False
    adata.var.loc[adata.var_names[:3], "highly_variable"] = True

    out = perform_ora(
        adata=adata,
        gene_sets={"GS_A": adata.var_names[:3].tolist(), "GS_B": adata.var_names[5:8].tolist()},
        gene_list=None,
        min_size=1,
        max_size=100,
    )
    assert out.method == "ora"
    assert out.gene_set_statistics["GS_A"]["query_size"] == 3


def test_perform_ora_fallback_cv_and_case_insensitive_gene_matching(
    monkeypatch: pytest.MonkeyPatch,
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)

    # CV fallback: no rank_genes_groups and no highly_variable in var.
    out_cv = perform_ora(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1", "gene_2"]},
        gene_list=None,
        min_size=1,
        max_size=100,
    )
    assert out_cv.method == "ora"

    # Case-insensitive query matching (GENE_0 -> gene_0).
    out_case = perform_ora(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1"]},
        gene_list=["GENE_0", "GENE_1"],
        min_size=1,
        max_size=100,
    )
    assert out_case.gene_set_statistics["GS_A"]["query_size"] == 2


def test_top_n_desc_indices_handles_non_finite_and_bounds() -> None:
    values = np.array([0.2, np.nan, 1.5, 0.7, -1.0], dtype=float)

    top3 = _top_n_desc_indices(values, 3)
    assert list(top3) == [2, 3, 0]

    top_all = _top_n_desc_indices(values, 10)
    assert list(top_all) == [2, 3, 0, 4, 1]

    top_zero = _top_n_desc_indices(values, 0)
    assert top_zero.size == 0


def test_perform_ora_handles_empty_pvalues_after_size_filter(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
) -> None:
    adata = minimal_spatial_adata.copy()
    _patch_metadata_noop(monkeypatch)
    out = perform_ora(
        adata=adata,
        gene_sets={"too_small": ["gene_0"]},
        gene_list=["gene_0"],
        min_size=10,
        max_size=500,
    )
    assert out.enrichment_scores == {}
    assert out.pvalues == {}
    assert out.adjusted_pvalues == {}


def test_perform_ssgsea_large_batch_processing_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    ad = pytest.importorskip("anndata")
    n_obs, n_vars = 620, 8  # >500 triggers batch path
    X = np.arange(n_obs * n_vars, dtype=np.float32).reshape(n_obs, n_vars)
    adata = ad.AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"gene_{i}" for i in range(n_vars)]
    _patch_metadata_noop(monkeypatch)

    class _BatchRes:
        def __init__(self, sample_names: list[str]):
            self.results = {
                name: pd.DataFrame({"Term": ["GS_A", "GS_B"], "ES": [0.1, 0.2]})
                for name in sample_names
            }

    backend = SimpleNamespace(
        ssgsea=lambda **kwargs: _BatchRes(list(kwargs["data"].columns))
    )
    resolved_contexts: list[object] = []
    monkeypatch.setattr(
        enrichment_module,
        "_get_gseapy",
        lambda ctx=None: resolved_contexts.append(ctx) or backend,
    )

    out = perform_ssgsea(
        adata=adata,
        gene_sets={"GS_A": ["gene_0", "gene_1"], "GS_B": ["gene_2", "gene_3"]},
        min_size=1,
        max_size=100,
    )
    assert out.method == "ssgsea"
    assert "ssgsea_GS_A" in adata.obs.columns
    assert "ssgsea_GS_B" in adata.obs.columns
    assert resolved_contexts == [None]


def test_perform_enrichr_propagates_service_failures(monkeypatch: pytest.MonkeyPatch) -> None:
    _patch_gseapy(
        monkeypatch,
        enrichr=lambda **_kwargs: (_ for _ in ()).throw(
            RuntimeError("service unavailable")
        ),
    )
    with pytest.raises(RuntimeError, match="service unavailable"):
        perform_enrichr(gene_list=["G1", "G2"], gene_sets="KEGG_Pathways", organism="human")


@pytest.mark.asyncio
@pytest.mark.parametrize(
    "score_keys, expected_name",
    [
        (None, "enrichmap_signature"),
        ("my_sig", "my_sig"),
        (["first_sig", "second_sig"], "first_sig"),
    ],
)
async def test_perform_spatial_enrichment_normalizes_list_input_and_score_keys(
    score_keys,
    expected_name: str,
    monkeypatch: pytest.MonkeyPatch,
    minimal_spatial_adata,
) -> None:
    adata = minimal_spatial_adata.copy()
    # Force conversion path by making dataset genes title-case while request is uppercase.
    adata.var_names = [f"Gene_{i}" for i in range(adata.n_vars)]
    ctx = _LogCtx(adata)

    def _fake_score(*, adata, score_key, **_kwargs):
        adata.obs[f"{score_key}_score"] = np.linspace(0.0, 1.0, adata.n_obs)

    fake_enrichmap = SimpleNamespace(tl=SimpleNamespace(score=_fake_score))
    monkeypatch.setattr(enrichment_module, "require", lambda *_a, **_k: fake_enrichmap)
    monkeypatch.setitem(sys.modules, "enrichmap", fake_enrichmap)
    monkeypatch.setattr(enrichment_module, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(enrichment_module, "export_analysis_result", lambda *_a, **_k: None)

    out = await perform_spatial_enrichment(
        data_id="d1",
        ctx=ctx,
        gene_sets=["GENE_0", "GENE_1"],
        score_keys=score_keys,
        species="mouse",
    )

    assert out.method == "spatial_enrichmap"
    assert f"{expected_name}_score" not in adata.obs.columns
    assert f"{expected_name}_score" in ctx._adata.obs.columns
    assert expected_name in ctx._adata.uns["enrichment_spatial_gene_sets"]


@pytest.mark.asyncio
async def test_perform_spatial_enrichment_warns_and_raises_when_all_signatures_invalid(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
) -> None:
    adata = minimal_spatial_adata.copy()
    ctx = _LogCtx(adata)

    fake_enrichmap = SimpleNamespace(tl=SimpleNamespace(score=lambda **_kwargs: None))
    monkeypatch.setattr(enrichment_module, "require", lambda *_a, **_k: fake_enrichmap)
    monkeypatch.setitem(
        sys.modules,
        "enrichmap",
        fake_enrichmap,
    )

    with pytest.raises(ProcessingError, match="No valid gene signatures found"):
        await perform_spatial_enrichment(
            data_id="d1",
            ctx=ctx,
            gene_sets={"too_small": ["missing_gene"]},
            species="human",
        )
    assert any("Skipping." in msg for msg in ctx.warnings)


def test_load_msigdb_additional_collections_and_species_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[str] = []

    def _fake_get_library(name: str, organism: str):
        calls.append(name)
        return {"ok": ["A", "B", "C"]}

    _patch_gseapy(monkeypatch, get_library=_fake_get_library)

    out_kegg_mouse = load_msigdb_gene_sets(
        species="mouse",
        collection="C2",
        subcollection="CP:KEGG",
        min_size=1,
        max_size=10,
    )
    out_kegg_human = load_msigdb_gene_sets(
        species="human",
        collection="C2",
        subcollection="CP:KEGG",
        min_size=1,
        max_size=10,
    )
    out_reactome = load_msigdb_gene_sets(
        species="human",
        collection="C2",
        subcollection="CP:REACTOME",
        min_size=1,
        max_size=10,
    )
    out_go_all = load_msigdb_gene_sets(
        species="human",
        collection="C5",
        subcollection=None,
        min_size=1,
        max_size=10,
    )
    out_c8 = load_msigdb_gene_sets(
        species="human",
        collection="C8",
        subcollection=None,
        min_size=1,
        max_size=10,
    )

    assert out_kegg_mouse == {"ok": ["A", "B", "C"]}
    assert out_kegg_human == {"ok": ["A", "B", "C"]}
    assert out_reactome == {"ok": ["A", "B", "C"]}
    assert out_go_all == {"ok": ["A", "B", "C"]}
    assert out_c8 == {"ok": ["A", "B", "C"]}
    assert "KEGG_2019_Mouse" in calls
    assert "KEGG_2021_Human" in calls
    assert "Reactome_Pathways_2024" in calls
    assert "GO_Biological_Process_2025" in calls
    assert "GO_Molecular_Function_2025" in calls
    assert "GO_Cellular_Component_2025" in calls
    assert "CellMarker_Augmented_2021" in calls


def test_loader_functions_wrap_external_failures_in_processing_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_gseapy(
        monkeypatch,
        get_library=lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
    )
    with pytest.raises(ProcessingError, match="Failed to load MSigDB gene sets"):
        load_msigdb_gene_sets("human")
    with pytest.raises(ProcessingError, match="Failed to load GO gene sets"):
        load_go_gene_sets("human", aspect="BP")
    with pytest.raises(ProcessingError, match="Failed to load KEGG pathways"):
        load_kegg_gene_sets("human")
    with pytest.raises(ProcessingError, match="Failed to load Reactome pathways"):
        load_reactome_gene_sets("human")
    with pytest.raises(ProcessingError, match="Failed to load cell markers"):
        load_cell_marker_gene_sets("human")
