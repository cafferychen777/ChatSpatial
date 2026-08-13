"""Contract tests for visualization main entrypoint behavior."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib.pyplot as plt
import pytest

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import main as viz_main
from chatspatial.utils.exceptions import (
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)


class _Ctx:
    def __init__(self, adata):
        self._adata = adata

    async def get_adata(self, _data_id: str):
        return self._adata


def _params(**overrides):
    base = {
        "plot_type": "feature",
        "subtype": None,
        "dpi": 120,
        "output_path": None,
        "output_format": "png",
    }
    base.update(overrides)
    return SimpleNamespace(**base)


async def test_visualize_data_rejects_unknown_plot_type():
    with pytest.raises(ParameterError, match="Invalid plot_type"):
        await viz_main.visualize_data("d1", _Ctx(None), _params(plot_type="unknown"))


async def test_visualize_data_rejects_dataset_with_no_observations():
    adata = SimpleNamespace(n_obs=0, n_vars=8)
    with pytest.raises(DataNotFoundError, match="no observations"):
        await viz_main.visualize_data("d1", _Ctx(adata), _params())


async def test_visualize_data_rejects_dataset_with_too_few_genes():
    adata = SimpleNamespace(n_obs=8, n_vars=4)
    with pytest.raises(DataNotFoundError, match="too few genes"):
        await viz_main.visualize_data("d1", _Ctx(adata), _params())


async def test_visualize_data_wraps_unexpected_handler_errors(monkeypatch):
    adata = SimpleNamespace(n_obs=8, n_vars=8)
    existing = plt.figure()
    created_number: int | None = None

    async def _boom(_adata, _params, _ctx):
        nonlocal created_number
        created_number = plt.figure().number
        raise RuntimeError("handler failed")

    monkeypatch.setattr(viz_main, "PLOT_HANDLERS", {"feature": _boom})

    try:
        with pytest.raises(
            ProcessingError, match="Failed to create feature visualization"
        ):
            await viz_main.visualize_data("d1", _Ctx(adata), _params())

        assert existing.number in plt.get_fignums()
        assert created_number not in plt.get_fignums()
    finally:
        plt.close(existing)


async def test_visualize_data_preserves_dependency_errors(monkeypatch):
    adata = SimpleNamespace(n_obs=8, n_vars=8)

    async def _missing_dependency(_adata, _params, _ctx):
        raise DependencyError("plot backend is unavailable")

    monkeypatch.setattr(
        viz_main,
        "PLOT_HANDLERS",
        {"feature": _missing_dependency},
    )

    with pytest.raises(DependencyError, match="plot backend is unavailable"):
        await viz_main.visualize_data("d1", _Ctx(adata), _params())


async def test_visualize_data_builds_plot_type_key_with_subtype(monkeypatch):
    adata = SimpleNamespace(n_obs=8, n_vars=8)
    captured = {"plot_type": None}

    async def _handler(_adata, _params, _ctx):
        return object()

    async def _optimize(fig, params, ctx, data_id=None, plot_type=None):
        assert fig is not None
        assert params.plot_type == "feature"
        assert data_id == "d1"
        captured["plot_type"] = plot_type
        return "saved"

    monkeypatch.setattr(viz_main, "PLOT_HANDLERS", {"feature": _handler})
    monkeypatch.setattr(viz_main, "optimize_fig_to_image_with_cache", _optimize)

    result = await viz_main.visualize_data("d1", _Ctx(adata), _params(subtype="umap"))

    assert result == "saved"
    assert captured["plot_type"] == "feature_umap"


# ---------------------------------------------------------------------------
# Regression: non-gene-dependent plot types must accept low n_vars
# ---------------------------------------------------------------------------


async def test_visualize_data_allows_few_genes_for_non_gene_plots(monkeypatch):
    """integration/deconvolution/cnv etc. should not be rejected for low n_vars."""
    adata = SimpleNamespace(n_obs=60, n_vars=3)
    sentinel = object()

    async def _handler(_adata, _params, _ctx):
        return sentinel

    async def _optimize(fig, params, ctx, data_id=None, plot_type=None):
        return "ok"

    monkeypatch.setattr(viz_main, "PLOT_HANDLERS", {"integration": _handler})
    monkeypatch.setattr(viz_main, "optimize_fig_to_image_with_cache", _optimize)

    # Should NOT raise DataNotFoundError for n_vars=3
    result = await viz_main.visualize_data(
        "d1",
        _Ctx(adata),
        _params(plot_type="integration"),
    )
    assert result == "ok"


async def test_visualize_data_still_rejects_few_genes_for_feature_type():
    """feature/expression types must still enforce the gene count check."""
    adata = SimpleNamespace(n_obs=60, n_vars=3)
    with pytest.raises(DataNotFoundError, match="too few genes"):
        await viz_main.visualize_data("d1", _Ctx(adata), _params(plot_type="feature"))


class _WarningCtx(_Ctx):
    def __init__(self, adata):
        super().__init__(adata)
        self.warnings: list[str] = []

    async def warning(self, msg: str):
        self.warnings.append(msg)


async def test_batch_key_is_not_ignored_in_silence(monkeypatch):
    """batch_key steers only the integration views.

    Elsewhere it reads like a request to split panels per sample; without a
    word, a multi-slice dataset comes back overplotted on a single panel.
    """
    adata = SimpleNamespace(n_obs=8, n_vars=8)
    ctx = _WarningCtx(adata)

    async def _handler(*_args, **_kwargs):
        return plt.figure()

    monkeypatch.setitem(viz_main.PLOT_HANDLERS, "feature", _handler)

    async def _export(*_a, **_k):
        return "out.png"

    monkeypatch.setattr(viz_main, "optimize_fig_to_image_with_cache", _export)

    params = VisualizationParameters(plot_type="feature", batch_key="batch")
    await viz_main.visualize_data("d1", ctx, params)

    assert [w for w in ctx.warnings if "batch_key is only used by" in w]


async def test_no_batch_warning_when_left_at_its_default(monkeypatch):
    adata = SimpleNamespace(n_obs=8, n_vars=8)
    ctx = _WarningCtx(adata)

    async def _handler(*_args, **_kwargs):
        return plt.figure()

    monkeypatch.setitem(viz_main.PLOT_HANDLERS, "feature", _handler)

    async def _export(*_a, **_k):
        return "out.png"

    monkeypatch.setattr(viz_main, "optimize_fig_to_image_with_cache", _export)

    await viz_main.visualize_data(
        "d1", ctx, VisualizationParameters(plot_type="feature")
    )

    assert not [w for w in ctx.warnings if "batch_key" in w]
