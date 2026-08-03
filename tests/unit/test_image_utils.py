"""Unit tests for image utility branch behavior."""

from __future__ import annotations

import asyncio
from pathlib import Path
from types import SimpleNamespace

import matplotlib
import matplotlib.pyplot as plt
import pytest

from chatspatial.utils import image_utils


class _FakeFigure:
    """Minimal figure stub that captures save parameters and writes output."""

    def __init__(self):
        self.saved_path: Path | None = None
        self.saved_kwargs: dict | None = None
        self._chatspatial_extra_artists = ("artist",)

    def savefig(self, path: str, **kwargs):
        self.saved_path = Path(path)
        self.saved_kwargs = kwargs
        self.saved_path.parent.mkdir(parents=True, exist_ok=True)
        self.saved_path.write_bytes(b"img")


def test_non_interactive_plotting_preserves_backend_and_restores_state(monkeypatch):
    original_backend = matplotlib.get_backend()
    was_interactive = plt.isinteractive()
    show_calls = {"count": 0}

    def _show(*_args, **_kwargs) -> None:
        show_calls["count"] += 1

    monkeypatch.setattr(plt, "show", _show)
    original_show = plt.show

    monkeypatch.setattr(
        matplotlib,
        "use",
        lambda *_args, **_kwargs: pytest.fail("backend switching is destructive"),
    )

    try:
        plt.ion()
        with image_utils.non_interactive_plotting():
            assert matplotlib.get_backend() == original_backend
            assert plt.isinteractive() is False
            assert plt.show is not original_show
            plt.show()

        assert matplotlib.get_backend() == original_backend
        assert plt.isinteractive() is True
        assert plt.show is original_show
        assert show_calls["count"] == 0
    finally:
        plt.interactive(was_interactive)


def test_non_interactive_plotting_restores_state_after_failure():
    original_show = plt.show
    was_interactive = plt.isinteractive()

    try:
        plt.ioff()
        with pytest.raises(RuntimeError, match="plot failed"):
            with image_utils.non_interactive_plotting():
                raise RuntimeError("plot failed")

        assert plt.isinteractive() is False
        assert plt.show is original_show
    finally:
        plt.interactive(was_interactive)


def test_isolated_figure_scope_closes_only_figures_created_inside():
    existing = plt.figure()
    created_number: int | None = None

    try:
        with pytest.raises(RuntimeError, match="plot failed"):
            with image_utils.isolated_figure_scope():
                created_number = plt.figure().number
                raise RuntimeError("plot failed")

        assert existing.number in plt.get_fignums()
        assert created_number not in plt.get_fignums()
    finally:
        plt.close(existing)


@pytest.mark.asyncio
async def test_serialized_figure_scope_excludes_concurrent_coroutines():
    active = 0
    max_active = 0

    async def render() -> None:
        nonlocal active, max_active
        async with image_utils.serialized_figure_scope():
            active += 1
            max_active = max(max_active, active)
            await asyncio.sleep(0.01)
            active -= 1

    await asyncio.gather(render(), render())

    assert max_active == 1


async def test_optimize_fig_to_image_with_cache_jpg_sets_quality_and_suffix(
    monkeypatch, tmp_path: Path
):
    fig = _FakeFigure()
    closed = {"fig": None}
    params = SimpleNamespace(
        dpi=144,
        output_format="jpg",
        output_path=str(tmp_path / "exports" / "figure_without_suffix"),
    )

    monkeypatch.setattr(plt, "close", lambda figure: closed.__setitem__("fig", figure))

    result = await image_utils.optimize_fig_to_image_with_cache(
        fig, params, plot_type="umap"
    )

    assert fig.saved_path is not None
    assert fig.saved_path.suffix == ".jpg"
    assert fig.saved_kwargs is not None
    assert fig.saved_kwargs["format"] == "jpg"
    assert fig.saved_kwargs["dpi"] == 144
    assert fig.saved_kwargs["pil_kwargs"] == {"quality": 95}
    assert "metadata" not in fig.saved_kwargs
    assert "Format: JPG" in result
    assert "Resolution: 144 DPI" in result
    assert closed["fig"] is fig


@pytest.mark.parametrize("output_format", ["jpeg", "svg", "eps", "tif", "tiff"])
async def test_optimize_fig_to_image_with_cache_removes_unsupported_metadata(
    monkeypatch, tmp_path: Path, output_format: str
):
    fig = _FakeFigure()
    params = SimpleNamespace(dpi=None, output_format=output_format, output_path=None)
    monkeypatch.setattr(image_utils, "get_default_output_dir", lambda: tmp_path)
    monkeypatch.setattr(
        image_utils.uuid,
        "uuid4",
        lambda: SimpleNamespace(hex="deadbeefcafebabe"),
    )
    monkeypatch.setattr(plt, "close", lambda _figure: None)

    result = await image_utils.optimize_fig_to_image_with_cache(
        fig, params, plot_type="scatter"
    )

    assert (
        fig.saved_path
        == tmp_path / "visualizations" / f"scatter_deadbeef.{output_format}"
    )
    assert fig.saved_kwargs is not None
    assert fig.saved_kwargs["format"] == output_format
    assert fig.saved_kwargs["dpi"] == 100
    assert fig.saved_kwargs["bbox_extra_artists"] == ("artist",)
    assert "metadata" not in fig.saved_kwargs
    assert f"Format: {output_format.upper()}" in result
