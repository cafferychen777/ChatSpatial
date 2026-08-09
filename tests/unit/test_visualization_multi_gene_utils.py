"""Unit tests for multi-gene and LR visualization contracts."""

from __future__ import annotations

import sys
from types import ModuleType

import matplotlib.pyplot as plt
import numpy as np
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import multi_gene as viz_mg
from chatspatial.utils.exceptions import DataNotFoundError, ProcessingError


class DummyCtx:
    def __init__(self):
        self.infos: list[str] = []
        self.warnings: list[str] = []

    async def info(self, msg: str):
        self.infos.append(msg)

    async def warning(self, msg: str):
        self.warnings.append(msg)


@pytest.mark.asyncio
async def test_create_multi_gene_visualization_requires_requested_basis(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    with pytest.raises(DataNotFoundError, match="UMAP embedding not found"):
        await viz_mg.create_multi_gene_visualization(
            adata,
            VisualizationParameters(
                plot_type="expression",
                subtype="heatmap",
                feature=["gene_0"],
                basis="umap",
            ),
        )
    adata.obsm.pop("spatial", None)
    with pytest.raises(DataNotFoundError, match="Spatial coordinates not found"):
        await viz_mg.create_multi_gene_visualization(
            adata,
            VisualizationParameters(
                plot_type="expression",
                subtype="heatmap",
                feature=["gene_0"],
                basis="spatial",
            ),
        )


@pytest.mark.asyncio
async def test_create_multi_gene_visualization_spatial_success_and_cleanup(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_spatial"] = adata.obsm.pop("spatial")

    async def _validated(*_a, **_k):
        return ["gene_0", "gene_1"]

    monkeypatch.setattr(viz_mg, "get_validated_features", _validated)
    monkeypatch.setattr(
        viz_mg,
        "get_gene_expression",
        lambda _adata, gene: (
            np.linspace(0.0, 1.0, _adata.n_obs)
            if gene == "gene_0"
            else np.linspace(1.0, 0.0, _adata.n_obs)
        ),
    )

    def _plot_spatial(_adata, ax, values=None, params=None, show_colorbar=False, **_kw):
        vals = values if values is not None else np.zeros(_adata.n_obs)
        spatial_key = viz_mg.get_spatial_key(_adata)
        assert spatial_key is not None
        ax.scatter(
            _adata.obsm[spatial_key][:, 0],
            _adata.obsm[spatial_key][:, 1],
            c=vals,
        )

    monkeypatch.setattr(viz_mg, "plot_spatial_feature", _plot_spatial)

    params = VisualizationParameters(
        plot_type="expression",
        subtype="heatmap",
        basis="spatial",
        color_scale="sqrt",
        show_colorbar=True,
        add_gene_labels=True,
        feature=["gene_0", "gene_1"],
    )
    ctx = DummyCtx()
    fig = await viz_mg.create_multi_gene_visualization(adata, params, context=ctx)

    # No temporary obs columns should be created (values= is used directly)
    assert all(not col.startswith("_") for col in adata.obs.columns)
    assert any("Visualizing 2 genes on spatial" in msg for msg in ctx.infos)
    assert "gene_0" in fig.axes[0].get_title()
    fig.clf()


@pytest.mark.asyncio
async def test_create_multi_gene_visualization_umap_branch(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.random.default_rng(0).normal(size=(adata.n_obs, 2))

    async def _validated(*_a, **_k):
        return ["gene_2"]

    monkeypatch.setattr(viz_mg, "get_validated_features", _validated)
    monkeypatch.setattr(
        viz_mg, "get_gene_expression", lambda _adata, _gene: np.arange(_adata.n_obs)
    )

    captured: dict[str, object] = {}

    def _umap(*_args, **kwargs):
        captured["color"] = kwargs.get("color")
        captured["vmin"] = kwargs.get("vmin")
        captured["vmax"] = kwargs.get("vmax")
        kwargs["ax"].scatter([0], [0], c=[1.0])

    monkeypatch.setattr(viz_mg.sc.pl, "umap", _umap)

    fig = await viz_mg.create_multi_gene_visualization(
        adata,
        VisualizationParameters(
            plot_type="expression",
            subtype="heatmap",
            basis="umap",
            feature=["gene_2"],
            color_scale="log",
        ),
        context=DummyCtx(),
    )
    assert captured["color"] == "_multi_gene_umap_temp"
    assert captured["vmax"] >= captured["vmin"]
    # try/finally ensures temp column is cleaned up even on error
    assert "_multi_gene_umap_temp" not in adata.obs.columns
    fig.clf()


def test_parse_lr_pairs_supports_multiple_sources(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    params = VisualizationParameters(
        plot_type="interaction",
        feature=["L1^R1", "L2^R2"],
    )
    assert viz_mg._parse_lr_pairs(adata, params) == [("L1", "R1"), ("L2", "R2")]

    # Unified CCC storage contract: adata.uns["ccc"]["top_lr_pairs"]
    adata.uns["ccc"] = {"top_lr_pairs": ["A^B", "C^D"]}
    out = viz_mg._parse_lr_pairs(
        adata,
        VisualizationParameters(plot_type="interaction"),
    )
    assert out == [("A", "B"), ("C", "D")]

    explicit = viz_mg._parse_lr_pairs(
        minimal_spatial_adata.copy(),
        VisualizationParameters(
            plot_type="interaction",
            lr_pairs=[("LIG", "REC")],
        ),
    )
    assert explicit == [("LIG", "REC")]


def test_parse_lr_pairs_raises_when_missing(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.uns["ccc"] = "foreign metadata"
    with pytest.raises(
        DataNotFoundError, match="No ligand-receptor pairs to visualize"
    ):
        viz_mg._parse_lr_pairs(
            adata,
            VisualizationParameters(plot_type="interaction"),
        )


@pytest.mark.asyncio
async def test_create_lr_pairs_visualization_requires_available_pairs(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(viz_mg, "_parse_lr_pairs", lambda *_a, **_k: [("x", "y")])
    with pytest.raises(DataNotFoundError, match="None of the specified LR pairs found"):
        await viz_mg.create_lr_pairs_visualization(
            adata, VisualizationParameters(plot_type="interaction"), context=DummyCtx()
        )


@pytest.mark.asyncio
async def test_create_lr_pairs_visualization_limits_pairs_and_cleans_temp(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    pairs = [
        ("gene_0", "gene_1"),
        ("gene_2", "gene_3"),
        ("gene_4", "gene_5"),
        ("gene_6", "gene_7"),
        ("gene_8", "gene_9"),
    ]
    monkeypatch.setattr(viz_mg, "_parse_lr_pairs", lambda *_a, **_k: pairs)
    monkeypatch.setattr(viz_mg, "ensure_unique_var_names", lambda _a: None)
    monkeypatch.setattr(
        viz_mg,
        "get_gene_expression",
        lambda _adata, gene: np.linspace(0.0, 1.0, _adata.n_obs)
        + (int(gene.split("_")[1]) * 0.01),
    )
    monkeypatch.setattr(
        viz_mg,
        "plot_spatial_feature",
        lambda _adata, ax, values=None, params=None, show_colorbar=False, **_kw: ax.scatter(
            _adata.obsm["spatial"][:, 0],
            _adata.obsm["spatial"][:, 1],
            c=values if values is not None else np.zeros(_adata.n_obs),
        ),
    )

    ctx = DummyCtx()
    fig = await viz_mg.create_lr_pairs_visualization(
        adata,
        VisualizationParameters(
            plot_type="interaction",
            show_colorbar=False,
            show_correlation_stats=True,
            correlation_method="spearman",
        ),
        context=ctx,
    )
    assert any("Limiting to first 4" in msg for msg in ctx.warnings)
    # No temporary obs columns created (values= passed directly)
    assert all(not col.startswith("_") for col in adata.obs.columns)
    fig.clf()


@pytest.mark.asyncio
@pytest.mark.parametrize("color_scale", ["log", "sqrt"])
async def test_create_lr_pairs_visualization_covers_scaling_colorbar_and_pearson(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch, color_scale: str
):
    adata = minimal_spatial_adata.copy()
    monkeypatch.setattr(
        viz_mg, "_parse_lr_pairs", lambda *_a, **_k: [("gene_0", "gene_1")]
    )
    monkeypatch.setattr(viz_mg, "ensure_unique_var_names", lambda _a: None)

    def _expr(_adata, gene):
        base = np.linspace(0.0, 2.0, _adata.n_obs)
        return base if gene == "gene_0" else base * 0.6 + 0.2

    monkeypatch.setattr(viz_mg, "get_gene_expression", _expr)
    import scipy.stats as _stats

    monkeypatch.setattr(_stats, "pearsonr", lambda *_a, **_k: (0.5, 0.01))
    monkeypatch.setattr(
        viz_mg,
        "plot_spatial_feature",
        lambda _adata, ax, values=None, params=None, show_colorbar=False, **_kw: ax.scatter(
            _adata.obsm["spatial"][:, 0],
            _adata.obsm["spatial"][:, 1],
            c=values if values is not None else np.zeros(_adata.n_obs),
        ),
    )

    fig = await viz_mg.create_lr_pairs_visualization(
        adata,
        VisualizationParameters(
            plot_type="interaction",
            show_colorbar=True,
            show_correlation_stats=True,
            correlation_method="pearson",
            color_scale=color_scale,
        ),
        context=DummyCtx(),
    )
    assert "Correlation:" in fig.axes[2].get_title()
    assert "lr_expr_temp_viz_99_unique" not in adata.obs.columns
    fig.clf()


@pytest.mark.asyncio
async def test_create_lr_pairs_visualization_without_spatial_plots_correlation_only(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    monkeypatch.setattr(
        viz_mg, "_parse_lr_pairs", lambda *_a, **_k: [("gene_0", "gene_1")]
    )
    monkeypatch.setattr(viz_mg, "ensure_unique_var_names", lambda _a: None)
    monkeypatch.setattr(
        viz_mg, "get_gene_expression", lambda _adata, _gene: np.arange(_adata.n_obs)
    )

    fig = await viz_mg.create_lr_pairs_visualization(
        adata,
        VisualizationParameters(
            plot_type="interaction",
            show_correlation_stats=False,
            correlation_method="kendall",
        ),
    )
    assert len(fig.axes) >= 1
    assert "vs" in fig.axes[0].get_title()
    fig.clf()


@pytest.mark.asyncio
async def test_create_gene_correlation_visualization_success(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata
):
    adata = minimal_spatial_adata.copy()

    async def _validated(*_a, **_k):
        return ["gene_0", "gene_1", "gene_2"]

    monkeypatch.setattr(viz_mg, "get_validated_features", _validated)
    monkeypatch.setattr(
        viz_mg,
        "get_genes_expression",
        lambda _adata, genes: np.vstack(
            [
                np.linspace(0.0, 1.0, _adata.n_obs),
                np.linspace(1.0, 0.0, _adata.n_obs),
                np.ones(_adata.n_obs),
            ]
        ).T,
    )

    class _FakeCluster:
        def __init__(self):
            self.fig = plt.figure()

    fake_sns = ModuleType("seaborn")
    fake_sns.clustermap = lambda *args, **kwargs: _FakeCluster()
    monkeypatch.setitem(sys.modules, "seaborn", fake_sns)

    fig = await viz_mg.create_gene_correlation_visualization(
        adata,
        VisualizationParameters(plot_type="expression", correlation_method="pearson"),
        context=DummyCtx(),
    )
    assert "Gene Correlation" in fig._suptitle.get_text()
    fig.clf()


@pytest.mark.asyncio
@pytest.mark.parametrize("color_scale", ["log", "sqrt"])
async def test_create_gene_correlation_visualization_color_scale_branches(
    monkeypatch: pytest.MonkeyPatch, minimal_spatial_adata, color_scale: str
):
    adata = minimal_spatial_adata.copy()

    async def _validated(*_a, **_k):
        return ["gene_0", "gene_1"]

    monkeypatch.setattr(viz_mg, "get_validated_features", _validated)
    monkeypatch.setattr(
        viz_mg,
        "get_genes_expression",
        lambda _adata, genes: np.vstack(
            [np.linspace(0.0, 2.0, _adata.n_obs), np.linspace(1.0, 3.0, _adata.n_obs)]
        ).T,
    )

    class _FakeCluster:
        def __init__(self):
            self.fig = plt.figure()

    fake_sns = ModuleType("seaborn")
    fake_sns.clustermap = lambda *args, **kwargs: _FakeCluster()
    monkeypatch.setitem(sys.modules, "seaborn", fake_sns)

    fig = await viz_mg.create_gene_correlation_visualization(
        adata,
        VisualizationParameters(
            plot_type="expression",
            correlation_method="pearson",
            color_scale=color_scale,
        ),
        context=DummyCtx(),
    )
    assert fig is not None
    fig.clf()


@pytest.mark.asyncio
async def test_spatial_interaction_visualization_validates_pairs_and_wraps_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    with pytest.raises(ProcessingError, match="No ligand-receptor pairs provided"):
        await viz_mg.create_spatial_interaction_visualization(
            adata,
            VisualizationParameters(plot_type="interaction", lr_pairs=None),
            context=DummyCtx(),
        )

    monkeypatch.setattr(
        viz_mg,
        "require_spatial_coords",
        lambda _a: (_ for _ in ()).throw(ValueError("bad spatial")),
    )
    with pytest.raises(ValueError, match="bad spatial"):
        await viz_mg.create_spatial_interaction_visualization(
            adata,
            VisualizationParameters(
                plot_type="interaction", lr_pairs=[("gene_0", "gene_1")]
            ),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_spatial_interaction_visualization_success_and_missing_gene_warning(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    adata = minimal_spatial_adata.copy()
    rng = np.random.default_rng(42)
    monkeypatch.setattr(viz_mg, "require_spatial_coords", lambda _a: _a.obsm["spatial"])

    def _expr(_adata, gene):
        if gene == "gene_0":
            out = np.zeros(_adata.n_obs)
            out[:15] = rng.uniform(2.0, 3.0, 15)
            return out
        if gene == "gene_1":
            out = np.zeros(_adata.n_obs)
            out[10:35] = rng.uniform(2.0, 3.0, 25)
            return out
        return np.zeros(_adata.n_obs)

    monkeypatch.setattr(viz_mg, "get_gene_expression", _expr)
    ctx = DummyCtx()
    fig = await viz_mg.create_spatial_interaction_visualization(
        adata,
        VisualizationParameters(
            plot_type="interaction",
            lr_pairs=[("gene_0", "gene_1"), ("missing", "gene_1")],
        ),
        context=ctx,
    )
    assert "Ligand-Receptor" in fig.axes[0].get_title()
    assert any("not found in expression data" in msg for msg in ctx.warnings)
    fig.clf()
