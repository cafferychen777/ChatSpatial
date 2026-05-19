"""Unit tests for image utility branch behavior."""

from __future__ import annotations

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


def test_ensure_non_interactive_backend_switches_from_interactive_backend(
    monkeypatch,
):
    use_calls: list[str] = []
    ioff_calls = {"count": 0}

    monkeypatch.setattr(matplotlib, "get_backend", lambda: "MacOSX")
    monkeypatch.setattr(matplotlib, "use", lambda backend: use_calls.append(backend))
    monkeypatch.setattr(plt, "ioff", lambda: ioff_calls.__setitem__("count", 1))

    image_utils._ensure_non_interactive_backend()

    assert use_calls == ["Agg"]
    assert ioff_calls["count"] == 1


def test_ensure_non_interactive_backend_is_noop_when_already_agg(monkeypatch):
    use_calls: list[str] = []
    ioff_calls = {"count": 0}

    monkeypatch.setattr(matplotlib, "get_backend", lambda: "Agg")
    monkeypatch.setattr(matplotlib, "use", lambda backend: use_calls.append(backend))
    monkeypatch.setattr(plt, "ioff", lambda: ioff_calls.__setitem__("count", 1))

    image_utils._ensure_non_interactive_backend()

    assert use_calls == []
    assert ioff_calls["count"] == 0


def test_non_interactive_backend_restores_original_backend(monkeypatch):
    use_calls: list[str] = []
    ioff_calls = {"count": 0}

    monkeypatch.setattr(matplotlib, "get_backend", lambda: "QtAgg")
    monkeypatch.setattr(matplotlib, "use", lambda backend: use_calls.append(backend))
    monkeypatch.setattr(plt, "ioff", lambda: ioff_calls.__setitem__("count", 1))

    with image_utils.non_interactive_backend():
        pass

    assert use_calls == ["Agg", "QtAgg"]
    assert ioff_calls["count"] == 1


def test_non_interactive_backend_is_noop_when_backend_is_agg(monkeypatch):
    use_calls: list[str] = []
    ioff_calls = {"count": 0}

    monkeypatch.setattr(matplotlib, "get_backend", lambda: "Agg")
    monkeypatch.setattr(matplotlib, "use", lambda backend: use_calls.append(backend))
    monkeypatch.setattr(plt, "ioff", lambda: ioff_calls.__setitem__("count", 1))

    with image_utils.non_interactive_backend():
        pass

    assert use_calls == []
    assert ioff_calls["count"] == 0


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

    monkeypatch.setattr(image_utils, "_ensure_non_interactive_backend", lambda: None)
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
    monkeypatch.setattr(image_utils, "_ensure_non_interactive_backend", lambda: None)
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

    assert fig.saved_path == tmp_path / "visualizations" / f"scatter_deadbeef.{output_format}"
    assert fig.saved_kwargs is not None
    assert fig.saved_kwargs["format"] == output_format
    assert fig.saved_kwargs["dpi"] == 100
    assert fig.saved_kwargs["bbox_extra_artists"] == ("artist",)
    assert "metadata" not in fig.saved_kwargs
    assert f"Format: {output_format.upper()}" in result
