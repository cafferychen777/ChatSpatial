"""Contract tests for visualization main entrypoint behavior."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib.pyplot as plt
import pytest

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
