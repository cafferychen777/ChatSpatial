"""Unit tests for gene-class exclusion in find_markers."""

from __future__ import annotations

from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError

from chatspatial.models.data import DifferentialExpressionParameters
from chatspatial.tools import differential as de
from chatspatial.utils.exceptions import DataError

pytestmark = pytest.mark.unit

# Group A over-expresses one gene of every class plus two ordinary markers.
CLASS_GENES = ["IGKC", "RPL13", "MT-CO1", "HLA-A", "MALAT1"]
MARKERS = ["KRT5", "EPCAM"]
BACKGROUND = [f"G{i}" for i in range(13)]
GENES = CLASS_GENES + MARKERS + BACKGROUND


class DummyCtx:
    def __init__(self, adata):
        self.adata = adata
        self.warnings: list[str] = []

    async def get_adata(self, _data_id):
        return self.adata

    async def set_adata(self, _data_id, adata):
        self.adata = adata

    async def warning(self, msg):
        self.warnings.append(msg)


@pytest.fixture(autouse=True)
def no_export(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(de, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(de, "export_analysis_result", lambda *_a, **_k: [])


def _adata() -> ad.AnnData:
    rng = np.random.default_rng(0)
    n = 40
    counts = rng.poisson(5, size=(n, len(GENES))).astype(np.float32)
    group_a = np.arange(n) < n // 2
    # Class genes are the strongest group-A signal, then the ordinary markers.
    counts[np.ix_(group_a, range(len(CLASS_GENES)))] += 200
    counts[np.ix_(group_a, range(len(CLASS_GENES), len(CLASS_GENES) + 2))] += 60
    adata = ad.AnnData(counts)
    adata.var_names = GENES
    adata.obs_names = [f"spot_{i}" for i in range(n)]
    adata.obs["group"] = pd.Categorical(np.where(group_a, "A", "B"))
    adata.obs["sample"] = [f"s{i % 4}" for i in range(n)]
    return adata


def _params(**overrides) -> DifferentialExpressionParameters:
    values = {
        "group_key": "group",
        "n_top_genes": 5,
        "exclude_gene_classes": [
            "immunoglobulin",
            "ribosomal",
            "mitochondrial",
            "hla",
        ],
        "exclude_gene_patterns": ["MALAT*"],
    }
    values.update(overrides)
    return DifferentialExpressionParameters(**values)


@pytest.mark.asyncio
async def test_without_exclusion_class_genes_lead_the_markers():
    out = await de.differential_expression(
        "d",
        DummyCtx(_adata()),
        DifferentialExpressionParameters(group_key="group", n_top_genes=5),
    )
    assert set(out.top_genes_by_group["A"]) == set(CLASS_GENES)
    assert out.n_genes_excluded is None


@pytest.mark.asyncio
async def test_one_vs_rest_never_reports_excluded_genes():
    out = await de.differential_expression("d", DummyCtx(_adata()), _params())

    reported = {g for genes in out.top_genes_by_group.values() for g in genes}
    assert not reported & set(CLASS_GENES)
    assert set(out.top_genes_by_group["A"][:2]) == set(MARKERS)
    assert out.n_genes_excluded == {
        "immunoglobulin": 1,
        "ribosomal": 1,
        "mitochondrial": 1,
        "hla": 1,
        "pattern:MALAT*": 1,
    }


@pytest.mark.asyncio
async def test_pairwise_comparison_excludes_genes_before_ranking():
    ctx = DummyCtx(_adata())
    out = await de.differential_expression(
        "d", ctx, _params(group1="A", group2="B", n_top_genes=20)
    )

    assert not set(out.top_genes_by_group["A"]) & set(CLASS_GENES)
    names = ctx.adata.uns["rank_genes_groups"]["names"]["A"]
    assert not set(names) & set(CLASS_GENES)
    assert len(names) == len(GENES) - len(CLASS_GENES)


@pytest.mark.asyncio
async def test_pydeseq2_receives_counts_without_excluded_genes(
    monkeypatch: pytest.MonkeyPatch,
):
    seen: dict[str, pd.DataFrame] = {}

    class FakeDataSet:
        def __init__(self, counts, metadata, design_factors):
            seen["counts"] = counts

        def deseq2(self):
            return None

    class FakeStats:
        def __init__(self, _dds, contrast):
            genes = seen["counts"].columns
            self.results_df = pd.DataFrame(
                {
                    "log2FoldChange": np.ones(len(genes)),
                    "padj": np.linspace(0.001, 0.5, len(genes)),
                },
                index=genes,
            )

        def summary(self):
            return None

    modules = {
        "pydeseq2.dds": SimpleNamespace(DeseqDataSet=FakeDataSet),
        "pydeseq2.ds": SimpleNamespace(DeseqStats=FakeStats),
    }
    monkeypatch.setattr(
        de, "require_module", lambda _package, module, *_a, **_k: modules[module]
    )
    adata = _adata()
    adata.obs["sample"] = [f"s{i % 4}_{g}" for i, g in enumerate(adata.obs["group"])]

    out = await de.differential_expression(
        "d",
        DummyCtx(adata),
        _params(method="pydeseq2", sample_key="sample", group1="A", group2="B"),
    )

    assert not set(seen["counts"].columns) & set(CLASS_GENES)
    assert not set(out.top_genes_by_group["A"]) & set(CLASS_GENES)
    assert out.n_genes_excluded["immunoglobulin"] == 1


@pytest.mark.asyncio
async def test_excluding_every_gene_is_an_error():
    with pytest.raises(DataError, match="Every gene"):
        await de.differential_expression(
            "d",
            DummyCtx(_adata()),
            _params(exclude_gene_classes=[], exclude_gene_patterns=["*"]),
        )


def test_exclusion_parameter_defaults_and_validation():
    defaults = DifferentialExpressionParameters(group_key="group")
    assert defaults.exclude_gene_classes == []
    assert defaults.exclude_gene_patterns == []
    with pytest.raises(ValidationError):
        DifferentialExpressionParameters(
            group_key="group", exclude_gene_classes=["ribosome"]
        )
