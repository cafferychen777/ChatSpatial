"""Unit tests for integration visualization branch contracts."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import integration as viz_integ
from chatspatial.utils.exceptions import DataError, DataNotFoundError, ParameterError


class DummyCtx:
    def __init__(self):
        self.infos: list[str] = []

    async def info(self, msg: str):
        self.infos.append(msg)


def _with_umap_and_batch(adata, *, n_batches: int = 2, categorical: bool = True):
    out = adata.copy()
    out.obsm["X_umap"] = np.column_stack([np.arange(out.n_obs), np.arange(out.n_obs)])
    labels = [f"b{i % n_batches}" for i in range(out.n_obs)]
    out.obs["batch"] = pd.Categorical(labels) if categorical else labels
    return out


@pytest.mark.asyncio
async def test_create_batch_integration_visualization_routes_subtypes(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    sentinel_batch = object()
    sentinel_cluster = object()
    sentinel_highlight = object()

    async def _batch(*_args, **_kwargs):
        return sentinel_batch

    async def _cluster(*_args, **_kwargs):
        return sentinel_cluster

    async def _highlight(*_args, **_kwargs):
        return sentinel_highlight

    monkeypatch.setattr(viz_integ, "_create_umap_by_batch", _batch)
    monkeypatch.setattr(viz_integ, "_create_umap_by_cluster", _cluster)
    monkeypatch.setattr(viz_integ, "_create_batch_highlight", _highlight)

    assert (
        await viz_integ.create_batch_integration_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="integration", subtype="batch"),
        )
        is sentinel_batch
    )
    assert (
        await viz_integ.create_batch_integration_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="integration", subtype="cluster"),
        )
        is sentinel_cluster
    )
    assert (
        await viz_integ.create_batch_integration_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="integration", subtype="highlight"),
        )
        is sentinel_highlight
    )


@pytest.mark.asyncio
async def test_create_batch_integration_visualization_rejects_unknown_subtype(
    minimal_spatial_adata,
):
    with pytest.raises(ParameterError, match="Unknown integration subtype"):
        await viz_integ.create_batch_integration_visualization(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="integration", subtype="bad"),
        )


@pytest.mark.asyncio
async def test_umap_by_batch_requires_batch_column(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2))

    with pytest.raises(DataError, match="Batch not found"):
        await viz_integ._create_umap_by_batch(
            adata, VisualizationParameters(plot_type="integration"), context=None
        )


@pytest.mark.asyncio
async def test_umap_by_batch_requires_umap_coordinates(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = ["b1"] * adata.n_obs

    with pytest.raises(DataNotFoundError, match="UMAP coordinates not found"):
        await viz_integ._create_umap_by_batch(
            adata,
            VisualizationParameters(plot_type="integration", batch_key="batch"),
            context=None,
        )


@pytest.mark.asyncio
@pytest.mark.parametrize("categorical", [True, False])
async def test_umap_by_batch_success_logs_and_renders_legend(
    minimal_spatial_adata, categorical: bool
):
    adata = _with_umap_and_batch(minimal_spatial_adata, categorical=categorical)
    ctx = DummyCtx()
    fig = await viz_integ._create_umap_by_batch(
        adata,
        VisualizationParameters(
            plot_type="integration",
            batch_key="batch",
            show_legend=True,
            title="Batch panel",
        ),
        context=ctx,
    )
    assert fig.axes[0].get_title() == "Batch panel"
    assert fig.axes[0].get_legend() is not None
    expected_colors = viz_integ.get_category_colors(2)
    actual_colors = [
        collection.get_facecolors()[0] for collection in fig.axes[0].collections
    ]
    np.testing.assert_allclose(
        [color[:3] for color in actual_colors],
        [color[:3] for color in expected_colors],
    )
    assert any("colored by batch" in m for m in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_umap_by_cluster_requires_umap_coordinates(minimal_spatial_adata):
    with pytest.raises(DataNotFoundError, match="UMAP coordinates not found"):
        await viz_integ._create_umap_by_cluster(
            minimal_spatial_adata,
            VisualizationParameters(plot_type="integration", subtype="cluster"),
            context=None,
        )


@pytest.mark.asyncio
async def test_umap_by_cluster_requires_cluster_column(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2))

    with pytest.raises(DataNotFoundError, match="Cluster column not found"):
        await viz_integ._create_umap_by_cluster(
            adata,
            VisualizationParameters(plot_type="integration", subtype="cluster"),
            context=None,
        )


@pytest.mark.asyncio
async def test_umap_by_cluster_auto_detects_cluster_and_uses_two_column_legend(
    minimal_spatial_adata,
):
    adata = minimal_spatial_adata.copy()
    adata.obsm["X_umap"] = np.column_stack(
        [np.arange(adata.n_obs), np.arange(adata.n_obs)]
    )
    # 12 clusters -> legend ncol path (2)
    adata.obs["cluster"] = pd.Categorical([f"c{i % 12}" for i in range(adata.n_obs)])
    ctx = DummyCtx()

    fig = await viz_integ._create_umap_by_cluster(
        adata,
        VisualizationParameters(
            plot_type="integration", subtype="cluster", show_legend=True
        ),
        context=ctx,
    )

    assert "UMAP by cluster" in fig.axes[0].get_title()
    legend = fig.axes[0].get_legend()
    assert legend is not None
    assert getattr(legend, "_ncols", 1) == 2
    assert any("colored by cluster" in m for m in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_batch_highlight_requires_batch_and_umap(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    with pytest.raises(DataError, match="Batch not found"):
        await viz_integ._create_batch_highlight(
            adata, VisualizationParameters(plot_type="integration"), context=None
        )

    adata.obs["batch"] = ["b1"] * adata.n_obs
    with pytest.raises(DataNotFoundError, match="UMAP coordinates not found"):
        await viz_integ._create_batch_highlight(
            adata,
            VisualizationParameters(plot_type="integration", batch_key="batch"),
            context=None,
        )


@pytest.mark.asyncio
async def test_batch_highlight_creates_grid_and_hides_unused_axes(
    minimal_spatial_adata,
):
    adata = _with_umap_and_batch(minimal_spatial_adata, n_batches=5, categorical=True)
    ctx = DummyCtx()
    fig = await viz_integ._create_batch_highlight(
        adata,
        VisualizationParameters(
            plot_type="integration",
            subtype="highlight",
            batch_key="batch",
            title="Highlight panel",
        ),
        context=ctx,
    )

    # 5 batches -> 2x4 grid = 8 axes, 3 hidden
    assert len(fig.axes) == 8
    hidden_axes = [ax for ax in fig.axes if not ax.axison]
    assert len(hidden_axes) == 3
    assert fig._suptitle is not None
    assert fig._suptitle.get_text() == "Highlight panel"
    assert any("per-batch highlight visualization" in m for m in ctx.infos)
    fig.clf()


@pytest.mark.asyncio
async def test_batch_highlight_respects_explicit_figure_size(minimal_spatial_adata):
    adata = _with_umap_and_batch(minimal_spatial_adata, n_batches=2, categorical=True)
    fig = await viz_integ._create_batch_highlight(
        adata,
        VisualizationParameters(
            plot_type="integration",
            subtype="highlight",
            batch_key="batch",
            figure_size=[9, 4],
        ),
        context=None,
    )
    width, height = fig.get_size_inches()
    assert width == pytest.approx(9.0)
    assert height == pytest.approx(4.0)
    fig.clf()


# ---------------------------------------------------------------------------
# Regression: all-NaN batch column must raise DataError, not ZeroDivisionError
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_batch_highlight_raises_on_all_nan_batch(minimal_spatial_adata):
    """batch_key with only NaN values must not ZeroDivisionError."""
    adata = minimal_spatial_adata.copy()
    adata.obs["batch"] = pd.array([pd.NA] * adata.n_obs, dtype=object)
    adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2))

    with pytest.raises(DataError, match="only missing values"):
        await viz_integ._create_batch_highlight(
            adata,
            VisualizationParameters(
                plot_type="integration",
                subtype="highlight",
                batch_key="batch",
            ),
            context=DummyCtx(),
        )
