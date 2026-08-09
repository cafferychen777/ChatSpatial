"""Unit tests for expression visualization routing and contracts."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import expression as expr
from chatspatial.utils.exceptions import ParameterError


class DummyCtx:
    def __init__(self):
        self.infos: list[str] = []

    async def info(self, msg: str):
        self.infos.append(msg)


@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("subtype", "target_attr"),
    [
        ("violin", "_create_violin"),
        ("dotplot", "_create_dotplot"),
        ("correlation", "_create_correlation"),
    ],
)
async def test_create_expression_visualization_routes_other_subtypes(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
    subtype: str,
    target_attr: str,
):
    sentinel = object()
    ctx = DummyCtx()

    async def _handler(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(expr, target_attr, _handler)
    out = await expr.create_expression_visualization(
        minimal_spatial_adata,
        VisualizationParameters(plot_type="expression", subtype=subtype),
        context=ctx,
    )
    assert out is sentinel
    assert any(subtype in m for m in ctx.infos)


@pytest.mark.asyncio
async def test_create_expression_visualization_routes_heatmap_by_subtype(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    sentinel = object()

    async def _fake_heatmap(*_args, **_kwargs):
        return sentinel

    monkeypatch.setattr(expr, "_create_heatmap", _fake_heatmap)
    out = await expr.create_expression_visualization(
        minimal_spatial_adata,
        VisualizationParameters(plot_type="expression", subtype="heatmap"),
    )
    assert out is sentinel


@pytest.mark.asyncio
async def test_create_expression_visualization_rejects_invalid_subtype(
    minimal_spatial_adata,
):
    with pytest.raises(ParameterError, match="Invalid expression subtype"):
        await expr.create_expression_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="expression", subtype="bad"),
        )


@pytest.mark.asyncio
async def test_heatmap_requires_cluster_key(minimal_spatial_adata):
    with pytest.raises(ParameterError, match="Heatmap requires cluster_key"):
        await expr._create_heatmap(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="expression", subtype="heatmap"),
            context=None,
        )


@pytest.mark.asyncio
async def test_heatmap_requires_valid_features(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    monkeypatch.setattr(expr, "validate_obs_column", lambda *_args, **_kwargs: None)

    async def _empty_features(*_args, **_kwargs):
        return []

    monkeypatch.setattr(expr, "get_validated_features", _empty_features)

    with pytest.raises(ParameterError, match="No valid gene features"):
        await expr._create_heatmap(
            adata,
            VisualizationParameters(
                plot_type="expression", subtype="heatmap", cluster_key="leiden"
            ),
            context=None,
        )


@pytest.mark.asyncio
async def test_heatmap_success_calls_scanpy_and_logs(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    ctx = DummyCtx()
    captured: dict[str, object] = {}

    monkeypatch.setattr(expr, "validate_obs_column", lambda *_args, **_kwargs: None)

    async def _features(*_args, **_kwargs):
        return ["gene_0", "gene_1"]

    monkeypatch.setattr(expr, "get_validated_features", _features)

    def _heatmap(*_args, **kwargs):
        captured.update(kwargs)
        plt.figure()
        return None

    monkeypatch.setattr(expr.sc.pl, "heatmap", _heatmap)

    fig = await expr._create_heatmap(
        adata,
        VisualizationParameters(
            plot_type="expression",
            subtype="heatmap",
            cluster_key="leiden",
            dotplot_dendrogram=True,
            dotplot_swap_axes=True,
            dotplot_standard_scale="var",
        ),
        context=ctx,
    )

    assert fig is not None
    assert captured["groupby"] == "leiden"
    assert captured["dendrogram"] is True
    assert captured["swap_axes"] is True
    assert captured["standard_scale"] == "var"
    assert any("Creating heatmap for 2 genes grouped by leiden" in m for m in ctx.infos)


@pytest.mark.asyncio
async def test_violin_requires_cluster_key(minimal_spatial_adata):
    with pytest.raises(ParameterError, match="Violin plot requires cluster_key"):
        await expr._create_violin(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="expression", subtype="violin"),
            context=None,
        )


@pytest.mark.asyncio
async def test_violin_requires_valid_features(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    monkeypatch.setattr(expr, "validate_obs_column", lambda *_args, **_kwargs: None)

    async def _empty_features(*_args, **_kwargs):
        return []

    monkeypatch.setattr(expr, "get_validated_features", _empty_features)

    with pytest.raises(
        ParameterError, match="No valid gene features provided for violin plot"
    ):
        await expr._create_violin(
            adata,
            VisualizationParameters(
                plot_type="expression", subtype="violin", cluster_key="leiden"
            ),
            context=None,
        )


@pytest.mark.asyncio
async def test_violin_success_calls_scanpy_and_logs(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    ctx = DummyCtx()
    captured: dict[str, object] = {}

    monkeypatch.setattr(expr, "validate_obs_column", lambda *_args, **_kwargs: None)

    async def _features(*_args, **_kwargs):
        return ["gene_0"]

    monkeypatch.setattr(expr, "get_validated_features", _features)

    def _violin(*_args, **kwargs):
        captured.update(kwargs)
        plt.figure()
        return None

    monkeypatch.setattr(expr.sc.pl, "violin", _violin)

    fig = await expr._create_violin(
        adata,
        VisualizationParameters(
            plot_type="expression", subtype="violin", cluster_key="leiden"
        ),
        context=ctx,
    )

    assert fig is not None
    assert captured["groupby"] == "leiden"
    assert any(
        "Creating violin plot for 1 genes grouped by leiden" in m for m in ctx.infos
    )


@pytest.mark.asyncio
async def test_dotplot_requires_cluster_key(minimal_spatial_adata):
    with pytest.raises(ParameterError, match="Dot plot requires cluster_key"):
        await expr._create_dotplot(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="expression", subtype="dotplot"),
            context=None,
        )


@pytest.mark.asyncio
async def test_dotplot_requires_valid_features(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    monkeypatch.setattr(expr, "validate_obs_column", lambda *_args, **_kwargs: None)

    async def _empty_features(*_args, **_kwargs):
        return []

    monkeypatch.setattr(expr, "get_validated_features", _empty_features)

    with pytest.raises(
        ParameterError, match="No valid gene features provided for dot plot"
    ):
        await expr._create_dotplot(
            adata,
            VisualizationParameters(
                plot_type="expression", subtype="dotplot", cluster_key="leiden"
            ),
            context=None,
        )


@pytest.mark.asyncio
async def test_dotplot_passes_optional_kwargs(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    called: dict[str, object] = {}

    monkeypatch.setattr(expr, "validate_obs_column", lambda *_args, **_kwargs: None)

    async def _features(*_args, **_kwargs):
        return ["gene_0", "gene_1"]

    monkeypatch.setattr(expr, "get_validated_features", _features)

    def _dotplot(**kwargs):
        called.update(kwargs)
        plt.figure()
        return None

    monkeypatch.setattr(expr.sc.pl, "dotplot", _dotplot)

    params = VisualizationParameters(
        plot_type="expression",
        subtype="dotplot",
        cluster_key="leiden",
        dotplot_dendrogram=True,
        dotplot_swap_axes=True,
        dotplot_standard_scale="var",
        dotplot_dot_min=0.1,
        dotplot_dot_max=0.9,
        dotplot_smallest_dot=2.0,
    )

    fig = await expr._create_dotplot(adata, params, context=None)
    assert fig is not None
    assert called["groupby"] == "leiden"
    assert called["dendrogram"] is True
    assert called["swap_axes"] is True
    assert called["standard_scale"] == "var"
    assert called["dot_min"] == pytest.approx(0.1)
    assert called["dot_max"] == pytest.approx(0.9)
    assert called["smallest_dot"] == pytest.approx(2.0)


@pytest.mark.asyncio
async def test_dotplot_maps_var_groups_and_logs(minimal_spatial_adata, monkeypatch):
    adata = minimal_spatial_adata.copy()
    adata.obs["leiden"] = ["0"] * adata.n_obs
    called: dict[str, object] = {}
    ctx = DummyCtx()

    monkeypatch.setattr(expr, "validate_obs_column", lambda *_args, **_kwargs: None)

    async def _features(*_args, **_kwargs):
        return ["gene_0", "gene_1"]

    monkeypatch.setattr(expr, "get_validated_features", _features)

    def _dotplot(**kwargs):
        called.update(kwargs)
        plt.figure()
        return None

    monkeypatch.setattr(expr.sc.pl, "dotplot", _dotplot)

    fig = await expr._create_dotplot(
        adata,
        VisualizationParameters(
            plot_type="expression",
            subtype="dotplot",
            cluster_key="leiden",
            dotplot_var_groups={"T cells": ["gene_0"], "B cells": ["gene_1"]},
        ),
        context=ctx,
    )

    assert fig is not None
    # var_groups dict is passed directly as var_names (scanpy native grouping)
    assert called["var_names"] == {"T cells": ["gene_0"], "B cells": ["gene_1"]}
    assert "var_group_positions" not in called
    assert "var_group_labels" not in called
    assert any(
        "Creating dot plot for 2 genes grouped by leiden" in m for m in ctx.infos
    )


@pytest.mark.asyncio
async def test_correlation_rejects_single_gene(minimal_spatial_adata, monkeypatch):
    """Regression: correlation requires >= 2 genes; 1 gene must fail early."""
    adata = minimal_spatial_adata.copy()

    async def _single_gene(*_args, **_kwargs):
        return ["gene_0"]

    monkeypatch.setattr(expr, "get_validated_features", _single_gene)

    with pytest.raises(ParameterError, match="at least 2 valid genes"):
        await expr._create_correlation(
            adata,
            VisualizationParameters(plot_type="expression", feature=["gene_0"]),
            context=DummyCtx(),
        )


@pytest.mark.asyncio
async def test_correlation_uses_requested_method_and_logs(
    minimal_spatial_adata, monkeypatch
):
    adata = minimal_spatial_adata.copy()
    ctx = DummyCtx()
    captured: dict[str, object] = {}

    async def _features(*_args, **_kwargs):
        return ["gene_0", "gene_1", "gene_2"]

    monkeypatch.setattr(expr, "get_validated_features", _features)
    monkeypatch.setattr(
        expr,
        "get_genes_expression",
        lambda *_args, **_kwargs: np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]]),
    )

    class _Grid:
        def __init__(self):
            self.fig = plt.figure()

    def _clustermap(corr_df, **kwargs):
        captured["corr"] = corr_df
        captured["kwargs"] = kwargs
        return _Grid()

    monkeypatch.setattr(expr.sns, "clustermap", _clustermap)

    fig = await expr._create_correlation(
        adata,
        VisualizationParameters(
            plot_type="expression",
            subtype="correlation",
            correlation_method="spearman",
        ),
        context=ctx,
    )

    assert fig is not None
    assert captured["kwargs"]["fmt"] == ".2f"
    assert any("gene correlation visualization" in m.lower() for m in ctx.infos)


@pytest.mark.asyncio
@pytest.mark.parametrize("color_scale", ["log", "sqrt"])
async def test_correlation_applies_color_scale_transform(
    minimal_spatial_adata, monkeypatch, color_scale: str
):
    adata = minimal_spatial_adata.copy()
    called = {"log": 0, "sqrt": 0}

    async def _features(*_args, **_kwargs):
        return ["gene_0", "gene_1"]

    monkeypatch.setattr(expr, "get_validated_features", _features)
    monkeypatch.setattr(
        expr,
        "get_genes_expression",
        lambda *_args, **_kwargs: np.array([[1.0, 4.0], [2.0, 8.0]]),
    )

    orig_log1p = expr.np.log1p
    orig_sqrt = expr.np.sqrt

    def _log1p(arr):
        called["log"] += 1
        return orig_log1p(arr)

    def _sqrt(arr):
        called["sqrt"] += 1
        return orig_sqrt(arr)

    monkeypatch.setattr(expr.np, "log1p", _log1p)
    monkeypatch.setattr(expr.np, "sqrt", _sqrt)

    class _Grid:
        def __init__(self):
            self.fig = plt.figure()

    monkeypatch.setattr(expr.sns, "clustermap", lambda *_args, **_kwargs: _Grid())

    fig = await expr._create_correlation(
        adata,
        VisualizationParameters(
            plot_type="expression", subtype="correlation", color_scale=color_scale
        ),
        context=None,
    )

    assert fig is not None
    if color_scale == "log":
        assert called["log"] == 1
        assert called["sqrt"] == 0
    else:
        assert called["sqrt"] == 1
        assert called["log"] == 0
