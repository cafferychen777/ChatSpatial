"""Unit tests for analyze_cnv gene positions, gene exclusion and query clustering."""

from __future__ import annotations

from types import ModuleType, SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError

from chatspatial.models.data import CNVParameters
from chatspatial.tools import cnv_analysis as cnv
from chatspatial.utils.exceptions import DataCompatibilityError, ParameterError
from chatspatial.utils.gene_classes import gene_class_mask

pytestmark = pytest.mark.unit

HUMAN_GENES = ["TP53", "EGFR", "MYC", "IGHG1", "RPL13", "MT-CO1", "HLA-A", "NOTAGENE"]


class DummyCtx:
    def __init__(self, adata):
        self.adata = adata
        self.warnings: list[str] = []
        self.committed = None

    async def get_adata(self, _data_id: str):
        return self.adata

    async def set_adata(self, _data_id: str, adata):
        self.adata = adata
        self.committed = adata

    async def warning(self, msg: str):
        self.warnings.append(msg)


def _adata(genes: list[str], n_obs: int = 20) -> ad.AnnData:
    adata = ad.AnnData(np.ones((n_obs, len(genes)), dtype=np.float32))
    adata.var_names = genes
    adata.obs_names = [f"spot_{i}" for i in range(n_obs)]
    half = n_obs // 2
    adata.obs["tumor_class"] = pd.Categorical(
        ["low_tumor"] * half + ["tumor_enriched"] * (n_obs - half)
    )
    return adata


@pytest.fixture
def fake_infercnvpy(monkeypatch: pytest.MonkeyPatch):
    """Replace infercnvpy with a stub that records the genes it receives."""
    seen: dict[str, object] = {}
    module = ModuleType("infercnvpy")

    def _infercnv(adata_obj, **kwargs):
        seen["genes"] = list(adata_obj.var_names)
        seen["var"] = adata_obj.var.copy()
        seen["kwargs"] = kwargs
        adata_obj.obsm["X_cnv"] = np.arange(adata_obj.n_obs * 2, dtype=float).reshape(
            adata_obj.n_obs, 2
        )
        adata_obj.uns["cnv"] = {"chr_pos": {"chr1": 0}}

    module.tl = SimpleNamespace(infercnv=_infercnv)
    monkeypatch.setitem(__import__("sys").modules, "infercnvpy", module)
    monkeypatch.setattr(cnv, "require", lambda name, *_a, **_k: module)
    monkeypatch.setattr(cnv, "export_analysis_result", lambda *_a, **_k: [])
    monkeypatch.setattr(cnv, "store_analysis_metadata", lambda *_a, **_k: None)
    return seen


def _params(**overrides) -> CNVParameters:
    values = {
        "method": "infercnvpy",
        "reference_key": "tumor_class",
        "reference_categories": ["low_tumor"],
    }
    values.update(overrides)
    return CNVParameters(**values)


def test_bundled_tables_are_compact_and_positioned():
    for species in ("human", "mouse"):
        table = cnv._load_gene_position_table(species)
        assert list(table.columns) == [
            "gene_id",
            "gene_name",
            "chromosome",
            "start",
            "end",
        ]
        assert len(table) > 50_000
        assert table["chromosome"].str.startswith("chr").all()
        assert (table["end"] >= table["start"]).all()
    resource_dir = cnv.files("chatspatial").joinpath("resources", "gene_positions")
    for filename, _ in cnv._GENE_POSITION_TABLES.values():
        with resource_dir.joinpath(filename).open("rb") as handle:
            assert len(handle.read()) < 2 * 1024 * 1024


def test_annotate_gene_positions_matches_symbols():
    adata = _adata(HUMAN_GENES)

    positions, summary = cnv._annotate_gene_positions(adata, "human")

    assert positions.loc["TP53", "chromosome"] == "chr17"
    assert positions.loc["TP53", "start"] == 7661779
    assert positions.loc["EGFR", "chromosome"] == "chr7"
    assert pd.isna(positions.loc["NOTAGENE", "chromosome"])
    # The bundled table covers chromosomes 1-22, X and Y only.
    assert pd.isna(positions.loc["MT-CO1", "chromosome"])
    assert summary["id_source"] is None
    assert summary["n_matched_by_id"] == 0
    assert summary["n_matched_by_symbol"] == 6
    assert summary["n_unmatched"] == 2


def test_annotate_gene_positions_prefers_ensembl_ids_over_symbols():
    adata = _adata(["renamed_tp53", "EGFR", "MYC"])
    adata.var["gene_ids"] = [
        "ENSG00000141510.18",  # TP53, with a version suffix
        "ENSG00000146648",
        "ENSG99999999999",  # unknown ID: falls back to the symbol
    ]

    positions, summary = cnv._annotate_gene_positions(adata, "human")

    assert positions.loc["renamed_tp53", "chromosome"] == "chr17"
    assert positions.loc["MYC", "chromosome"] == "chr8"
    assert summary["id_source"] == "var['gene_ids']"
    assert summary["n_matched_by_id"] == 2
    assert summary["n_matched_by_symbol"] == 1


def test_annotate_gene_positions_mouse():
    adata = _adata(["Trp53", "Actb", "TP53"])

    positions, summary = cnv._annotate_gene_positions(adata, "mouse")

    assert positions.loc["Trp53", "chromosome"] == "chr11"
    assert positions.loc["Actb", "chromosome"] == "chr5"
    assert summary["n_unmatched"] == 1


@pytest.mark.asyncio
async def test_infer_cnv_auto_annotates_missing_positions(fake_infercnvpy):
    ctx = DummyCtx(_adata(HUMAN_GENES))

    out = await cnv.infer_cnv("d1", ctx, _params())

    assert out.gene_positions["source"] == "Ensembl 110 GRCh38"
    assert out.gene_positions["n_matched_by_symbol"] == 6
    assert "NOTAGENE" not in fake_infercnvpy["genes"]
    assert "TP53" in fake_infercnvpy["genes"]
    assert ctx.committed.var.loc["TP53", "chromosome"] == "chr17"
    assert pd.isna(ctx.committed.var.loc["NOTAGENE", "start"])
    assert ctx.committed.uns["cnv_analysis_infercnvpy_low_tumor"]["gene_positions"] == (
        "Ensembl 110 GRCh38"
    )


@pytest.mark.asyncio
async def test_infer_cnv_auto_keeps_provided_positions(fake_infercnvpy):
    adata = _adata(["g1", "g2", "g3"])
    adata.var["chromosome"] = ["chr1", "chr1", "chr2"]
    adata.var["start"] = [1, 100, 1]
    adata.var["end"] = [50, 150, 50]
    ctx = DummyCtx(adata)

    out = await cnv.infer_cnv("d1", ctx, _params())

    assert out.gene_positions == {"source": "adata.var"}
    assert fake_infercnvpy["genes"] == ["g1", "g2", "g3"]


@pytest.mark.asyncio
async def test_infer_cnv_annotate_rejects_unmatched_species(fake_infercnvpy):
    ctx = DummyCtx(_adata(["Trp53", "Actb"]))

    with pytest.raises(DataCompatibilityError, match="species"):
        await cnv.infer_cnv("d1", ctx, _params(gene_positions="annotate"))


def test_gene_class_masks():
    names = pd.Index(
        [
            "IGHG1",
            "IGKC",
            "IGLC2",
            "IGHV3-23",
            "JCHAIN",
            "IGLL5",
            "IGHMBP2",
            "Ighm",
            "RPL13",
            "RPS6KA1",
            "RPS19BP1",
            "Rpl13",
            "MT-CO1",
            "mt-Co1",
            "HLA-A",
            "H2-K1",
            "TP53",
        ]
    )
    selected = {
        cls: set(names[gene_class_mask(names, cls)])
        for cls in ("immunoglobulin", "ribosomal", "mitochondrial", "hla")
    }
    assert selected["immunoglobulin"] == {
        "IGHG1",
        "IGKC",
        "IGLC2",
        "IGHV3-23",
        "JCHAIN",
        "IGLL5",
        "Ighm",
    }
    assert selected["ribosomal"] == {"RPL13", "Rpl13"}
    assert selected["mitochondrial"] == {"MT-CO1", "mt-Co1"}
    assert selected["hla"] == {"HLA-A", "H2-K1"}


@pytest.mark.asyncio
async def test_infer_cnv_excludes_gene_classes_and_patterns(fake_infercnvpy):
    ctx = DummyCtx(_adata(HUMAN_GENES))

    out = await cnv.infer_cnv(
        "d1",
        ctx,
        _params(
            exclude_gene_classes=["immunoglobulin", "ribosomal", "hla"],
            exclude_gene_patterns=["MY*"],
        ),
    )

    assert sorted(fake_infercnvpy["genes"]) == ["EGFR", "TP53"]
    assert out.n_genes_excluded == {
        "immunoglobulin": 1,
        "ribosomal": 1,
        "hla": 1,
        "pattern:MY*": 1,
    }
    assert out.statistics["n_genes_excluded"] == 4
    assert out.n_genes_analyzed == 2
    # Excluded genes leave the CNV input only; the dataset keeps them.
    assert "IGHG1" in ctx.committed.var_names


@pytest.mark.asyncio
async def test_infer_cnv_clusters_only_query_categories(
    fake_infercnvpy, monkeypatch: pytest.MonkeyPatch
):
    clustered: dict[str, object] = {}

    def _neighbors(adata_obj, **_kwargs):
        clustered["obs_names"] = list(adata_obj.obs_names)

    def _leiden(adata_obj, key_added):
        adata_obj.obs[key_added] = pd.Categorical(["0", "1"] * (adata_obj.n_obs // 2))

    monkeypatch.setattr(cnv.sc.pp, "neighbors", _neighbors)
    monkeypatch.setattr(cnv.sc.tl, "leiden", _leiden)
    ctx = DummyCtx(_adata(HUMAN_GENES))

    out = await cnv.infer_cnv(
        "d1",
        ctx,
        _params(cluster_cells=True, query_categories=["tumor_enriched"]),
    )

    assert clustered["obs_names"] == [f"spot_{i}" for i in range(10, 20)]
    labels = ctx.committed.obs["cnv_clusters"]
    assert labels.iloc[:10].isna().all()
    assert labels.iloc[10:].notna().all()
    assert list(labels.cat.categories) == ["0", "1"]
    assert out.n_query_cells == 10
    # Inference itself still sees every observation.
    assert fake_infercnvpy["kwargs"]["reference_cat"] == ["low_tumor"]
    assert ctx.committed.obsm["X_cnv"].shape[0] == 20


@pytest.mark.asyncio
async def test_infer_cnv_rejects_unknown_query_category(fake_infercnvpy):
    ctx = DummyCtx(_adata(HUMAN_GENES))

    with pytest.raises(ParameterError, match="Query categories"):
        await cnv.infer_cnv(
            "d1", ctx, _params(cluster_cells=True, query_categories=["tumour"])
        )


def test_cnv_parameter_validation():
    with pytest.raises(ValidationError, match="cluster_cells"):
        _params(query_categories=["tumor_enriched"])
    with pytest.raises(ValidationError, match="overlap"):
        _params(cluster_cells=True, query_categories=["low_tumor"])
    with pytest.raises(ValidationError, match="infercnvpy only"):
        _params(method="numbat", exclude_gene_classes=["ribosomal"])
    with pytest.raises(ValidationError):
        _params(exclude_gene_classes=["ribosome"])
    with pytest.raises(ValidationError):
        _params(species="zebrafish")
    defaults = _params()
    assert defaults.gene_positions == "auto"
    assert defaults.species == "human"
    assert defaults.exclude_gene_classes == []
    assert defaults.exclude_gene_patterns == []
    assert defaults.query_categories is None
