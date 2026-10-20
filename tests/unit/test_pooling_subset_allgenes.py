"""Unit tests for concatenate-only integration, CCC observation subsets and all-gene Moran's I."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from pydantic import ValidationError

from chatspatial.models.data import (
    CellCommunicationParameters,
    IntegrationParameters,
    SpatialStatisticsParameters,
)
from chatspatial.tools import cell_communication as ccc
from chatspatial.tools import spatial_statistics as spstats
from chatspatial.tools.cell_communication import CCCStorage
from chatspatial.tools.integration import integrate_samples
from chatspatial.utils.exceptions import DataError, ParameterError

pytestmark = pytest.mark.unit


class DummyCtx:
    def __init__(self, datasets: dict[str, ad.AnnData]):
        self.datasets = datasets
        self.added: dict[str, ad.AnnData] = {}
        self.warnings: list[str] = []

    async def get_adata(self, data_id: str):
        return self.datasets[data_id]

    async def set_adata(self, data_id: str, adata):
        self.datasets[data_id] = adata

    def reserve_dataset_id(self, prefix: str = "data"):
        return f"{prefix}_1"

    async def add_dataset(self, adata, prefix: str = "data", data_id=None):
        self.added[data_id or f"{prefix}_1"] = adata
        return data_id

    async def warning(self, msg: str):
        self.warnings.append(msg)


@pytest.fixture
def no_export(monkeypatch: pytest.MonkeyPatch):
    for target in (
        "chatspatial.tools.integration.export_analysis_result",
        "chatspatial.tools.spatial_statistics.export_analysis_result",
        "chatspatial.utils.results_export.export_analysis_result",
    ):
        monkeypatch.setattr(target, lambda *_a, **_k: [])


def _section(name: str, genes: list[str], seed: int, n_obs: int = 30) -> ad.AnnData:
    """A preprocessed section: counts in raw and layers, log-normalized X."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(3, size=(n_obs, len(genes))).astype(np.float32)
    adata = ad.AnnData(counts)
    adata.var_names = genes
    adata.obs_names = [f"AAAC-{i}" for i in range(n_obs)]  # same barcodes per section
    adata.obs["region"] = pd.Categorical(["TC", "LE", "Other"] * (n_obs // 3))
    adata.obs["dominant_celltype_flashdeconv"] = pd.Categorical(
        ["Malignant", "Fibroblast"] * (n_obs // 2)
    )
    adata.obsm["spatial"] = rng.uniform(0, 100, size=(n_obs, 2))
    adata.uns["spatial"] = {name: {"scalefactors": {"spot_diameter_fullres": 1.0}}}
    adata.raw = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


@pytest.mark.asyncio
async def test_concatenate_only_pools_counts_of_shared_genes(no_export):
    s1 = _section("lib1", [f"G{i}" for i in range(10)], seed=1)
    s2 = _section("lib2", [f"G{i}" for i in range(2, 12)], seed=2)
    ctx = DummyCtx({"s1": s1, "s2": s2})

    result = await integrate_samples(
        ["s1", "s2"],
        ctx,
        IntegrationParameters(batch_correction=False, batch_key="sample"),
    )

    pooled = ctx.added["integrated_1"]
    assert result.integration_method == "concatenate"
    assert result.embedding_key is None
    assert result.n_cells == 60
    assert result.n_genes == 8  # G2..G9 are measured in both sections
    assert list(pooled.var_names) == [f"G{i}" for i in range(2, 10)]
    assert pooled.obs["sample"].astype(str).tolist() == ["s1"] * 30 + ["s2"] * 30
    assert pooled.obs_names.is_unique
    expected = np.vstack(
        [
            s1.raw.to_adata()[:, pooled.var_names].X,
            s2.raw.to_adata()[:, pooled.var_names].X,
        ]
    )
    np.testing.assert_allclose(pooled.layers["counts"].toarray(), expected)
    totals = expected.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(
        pooled.X.toarray(), np.log1p(expected / totals * 1e4), rtol=1e-5
    )
    assert pooled.obs["dominant_celltype_flashdeconv"].notna().all()
    assert pooled.obsm["spatial"].shape == (60, 2)
    assert set(pooled.uns["spatial"]) == {"lib1", "lib2"}
    assert pooled.uns["integration_has_counts"] is True
    assert "integration_concatenate_metadata" in pooled.uns
    assert "X_pca" not in pooled.obsm


@pytest.mark.asyncio
async def test_concatenate_only_reads_raw_when_no_counts_layer(no_export):
    sections = {}
    for name, seed in (("s1", 1), ("s2", 2)):
        adata = _section(name, [f"G{i}" for i in range(8)], seed=seed)
        del adata.layers["counts"]
        sections[name] = adata
    ctx = DummyCtx(sections)

    await integrate_samples(
        ["s1", "s2"], ctx, IntegrationParameters(batch_correction=False)
    )

    pooled = ctx.added["integrated_1"]
    np.testing.assert_allclose(
        pooled.layers["counts"].toarray(),
        np.vstack([sections["s1"].raw.X, sections["s2"].raw.X]),
    )


@pytest.mark.asyncio
async def test_concatenate_only_refuses_mixed_counts_and_normalized(no_export):
    with_counts = _section("lib1", [f"G{i}" for i in range(8)], seed=1)
    normalized = _section("lib2", [f"G{i}" for i in range(8)], seed=2)
    del normalized.layers["counts"]
    normalized.raw = None
    ctx = DummyCtx({"s1": with_counts, "s2": normalized})

    with pytest.raises(DataError, match="s2 hold no integer counts"):
        await integrate_samples(
            ["s1", "s2"], ctx, IntegrationParameters(batch_correction=False)
        )
    assert not ctx.added


@pytest.mark.asyncio
async def test_concatenate_only_pools_normalized_inputs_as_they_are(no_export):
    sections = {}
    for name, seed in (("s1", 1), ("s2", 2)):
        adata = _section(name, [f"G{i}" for i in range(8)], seed=seed)
        del adata.layers["counts"]
        adata.raw = None
        sections[name] = adata
    ctx = DummyCtx(sections)

    await integrate_samples(
        ["s1", "s2"], ctx, IntegrationParameters(batch_correction=False)
    )

    pooled = ctx.added["integrated_1"]
    assert "counts" not in pooled.layers
    assert pooled.uns["integration_has_counts"] is False
    np.testing.assert_allclose(
        pooled.X.toarray(), np.vstack([sections["s1"].X, sections["s2"].X]), rtol=1e-6
    )


@pytest.mark.asyncio
async def test_concatenate_only_accepts_unpreprocessed_counts(no_export):
    rng = np.random.default_rng(0)
    sections = {}
    for name in ("a", "b"):
        adata = ad.AnnData(rng.poisson(5, size=(20, 6)).astype(np.float32))
        adata.var_names = [f"G{i}" for i in range(6)]
        adata.obsm["spatial"] = rng.uniform(size=(20, 2))
        sections[name] = adata
    ctx = DummyCtx(sections)

    result = await integrate_samples(
        ["a", "b"], ctx, IntegrationParameters(batch_correction=False)
    )

    pooled = ctx.added["integrated_1"]
    assert result.n_cells == 40
    np.testing.assert_allclose(
        pooled.layers["counts"].toarray(),
        np.vstack([sections["a"].X, sections["b"].X]),
    )


def _ccc_adata() -> ad.AnnData:
    rng = np.random.default_rng(3)
    adata = ad.AnnData(rng.poisson(3, size=(30, 5)).astype(np.float32))
    adata.obs_names = [f"spot_{i}" for i in range(30)]
    adata.obs["tc_le_region"] = pd.Categorical(["TC", "LE", "Other"] * 10)
    labels = np.array(["Malignant", "Fibroblast", "Macrophage"] * 10, dtype=object)
    labels[2::3] = "Endothelial"  # only in "Other" spots
    adata.obs["dominant_celltype_flashdeconv"] = pd.Categorical(labels)
    adata.obsm["spatial"] = rng.uniform(size=(30, 2))
    return adata


def _ccc_params(**overrides) -> CellCommunicationParameters:
    values = {
        "method": "cellchat_r",
        "species": "human",
        "cell_type_key": "dominant_celltype_flashdeconv",
        "subset_key": "tc_le_region",
        "subset_categories": ["TC", "LE"],
    }
    values.update(overrides)
    return CellCommunicationParameters(**values)


def _storage(n_obs: int, spatial: bool = False) -> CCCStorage:
    return CCCStorage(
        method="cellchat_r",
        analysis_type="spatial" if spatial else "cluster",
        species="human",
        database="CellChatDB.human",
        lr_pairs=["A^B"],
        top_lr_pairs=["A^B"],
        n_pairs=1,
        n_significant=1,
        results=pd.DataFrame({"score": [1.0]}, index=["A^B"]),
        method_data=(
            {
                "spatial_scores": np.ones((n_obs, 1)),
                "spatial_pvals": np.zeros((n_obs, 1)),
            }
            if spatial
            else {}
        ),
    )


@pytest.fixture
def no_ccc_metadata(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        "chatspatial.utils.adata_utils.store_analysis_metadata", lambda *a, **k: None
    )
    monkeypatch.setattr(
        "chatspatial.utils.results_export.export_analysis_result",
        lambda *a, **k: None,
    )


@pytest.mark.parametrize("method", ["cellchat_r", "liana"])
@pytest.mark.asyncio
async def test_ccc_backends_receive_only_the_subset(
    monkeypatch: pytest.MonkeyPatch, no_ccc_metadata, method
):
    seen: dict[str, ad.AnnData] = {}

    def _fake_cellchat(adata, _params, _ctx):
        seen["adata"] = adata
        return _storage(adata.n_obs)

    async def _fake_liana(adata, _params, _ctx):
        seen["adata"] = adata
        return _storage(adata.n_obs)

    monkeypatch.setattr(ccc, "_analyze_communication_cellchat_r", _fake_cellchat)
    monkeypatch.setattr(ccc, "_analyze_communication_liana", _fake_liana)
    source = _ccc_adata()
    ctx = DummyCtx({"pooled": source})

    result = await ccc.analyze_cell_communication(
        "pooled",
        ctx,
        _ccc_params(method=method, perform_spatial_analysis=False),
    )

    analyzed = seen["adata"]
    assert analyzed.n_obs == 20
    assert set(analyzed.obs["tc_le_region"].astype(str)) == {"TC", "LE"}
    assert list(analyzed.obs["dominant_celltype_flashdeconv"].cat.categories) == [
        "Fibroblast",
        "Malignant",
    ]
    assert result.n_obs_analyzed == 20
    published = ctx.datasets["pooled"]
    assert published.n_obs == 30
    assert published.uns["ccc"]["statistics"]["n_obs_analyzed"] == 20


@pytest.mark.asyncio
async def test_ccc_subset_expands_per_spot_scores(
    monkeypatch: pytest.MonkeyPatch, no_ccc_metadata
):
    async def _fake_run(adata, _params, _ctx):
        return _storage(adata.n_obs, spatial=True)

    monkeypatch.setattr(ccc, "_run_ccc_analysis", _fake_run)
    monkeypatch.setattr(ccc, "_validate_ccc_params", lambda *a, **k: _noop())
    ctx = DummyCtx({"pooled": _ccc_adata()})

    await ccc.analyze_cell_communication(
        "pooled", ctx, _ccc_params(method="liana", perform_spatial_analysis=True)
    )

    scores = ctx.datasets["pooled"].obsm["ccc_spatial_scores_liana"]
    assert scores.shape == (30, 1)
    assert np.isnan(scores[2::3]).all()
    assert (scores[0::3] == 1).all()
    pvals = ctx.datasets["pooled"].obsm["ccc_spatial_pvals_liana"]
    assert np.isnan(pvals[2::3]).all()


async def _noop():
    return None


@pytest.mark.asyncio
async def test_ccc_subset_rejects_unknown_category(no_ccc_metadata):
    ctx = DummyCtx({"pooled": _ccc_adata()})
    with pytest.raises(ParameterError, match="subset_categories"):
        await ccc.analyze_cell_communication(
            "pooled", ctx, _ccc_params(subset_categories=["TC", "Edge"])
        )


def test_ccc_subset_parameters_come_in_pairs():
    with pytest.raises(ValidationError, match="both"):
        _ccc_params(subset_categories=None)
    with pytest.raises(ValidationError, match="both"):
        _ccc_params(subset_key=None)
    with pytest.raises(ValidationError, match="at least one"):
        _ccc_params(subset_categories=[])
    plain = _ccc_params(subset_key=None, subset_categories=None)
    assert plain.subset_key is None


def _moran_adata() -> ad.AnnData:
    rng = np.random.default_rng(5)
    n = 100
    coords = np.column_stack([np.repeat(np.arange(10), 10), np.tile(np.arange(10), 10)])
    gradient = coords[:, 0].astype(float)
    X = np.column_stack(
        [gradient + rng.normal(0, 0.1, n)] + [rng.normal(size=n) for _ in range(29)]
    ).astype(np.float32)
    adata = ad.AnnData(X)
    adata.var_names = ["GRAD"] + [f"N{i}" for i in range(29)]
    adata.obsm["spatial"] = coords.astype(float)
    adata.var["highly_variable"] = [True] * 5 + [False] * 25
    return adata


@pytest.mark.asyncio
async def test_moran_gene_selection_all_tables_every_gene(no_export):
    ctx = DummyCtx({"d": _moran_adata()})

    result = await spstats.analyze_spatial_statistics(
        "d",
        ctx,
        SpatialStatisticsParameters(
            analysis_type="moran", gene_selection="all", moran_n_perms=10
        ),
    )

    table = ctx.datasets["d"].uns["moranI"]
    assert result.n_features_analyzed == 30
    assert len(table) == 30
    assert result.results_key == "moranI"
    assert len(result.top_features) <= 10
    assert result.top_features[0] == "GRAD"
    assert table.loc["GRAD", "rank"] == 1
    assert sorted(table["rank"]) == list(range(1, 31))
    assert {"I", "pval_norm", "pval_norm_fdr"} <= set(table.columns)
    payload = result.model_dump(mode="json")
    assert "statistics" not in payload


@pytest.mark.asyncio
async def test_moran_default_still_uses_hvgs(no_export):
    ctx = DummyCtx({"d": _moran_adata()})

    result = await spstats.analyze_spatial_statistics(
        "d", ctx, SpatialStatisticsParameters(analysis_type="moran")
    )

    assert result.n_features_analyzed == 5


@pytest.mark.asyncio
async def test_geary_gene_selection_all_ranks_by_smallest_c(no_export):
    ctx = DummyCtx({"d": _moran_adata()})

    await spstats.analyze_spatial_statistics(
        "d",
        ctx,
        SpatialStatisticsParameters(analysis_type="geary", gene_selection="all"),
    )

    table = ctx.datasets["d"].uns["gearyC"]
    assert len(table) == 30
    assert table.loc["GRAD", "rank"] == 1


def test_gene_selection_all_validation():
    with pytest.raises(ValidationError, match="moran and geary"):
        SpatialStatisticsParameters(analysis_type="getis_ord", gene_selection="all")
    with pytest.raises(ValidationError, match="not both"):
        SpatialStatisticsParameters(
            analysis_type="moran", gene_selection="all", genes=["A"]
        )
    assert SpatialStatisticsParameters().gene_selection == "hvg"


def test_integration_batch_correction_default():
    assert IntegrationParameters().batch_correction is True
