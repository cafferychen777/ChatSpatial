"""
Image utilities for spatial transcriptomics MCP.

This module provides:
- Matplotlib backend management (prevent GUI popups)
- Figure export to file with user-specified format/DPI
"""

import threading
import uuid
from collections.abc import AsyncIterator, Generator
from contextlib import asynccontextmanager, contextmanager
from typing import TYPE_CHECKING, Any, Optional

import anyio

from ..config import get_default_output_dir

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


# =============================================================================
# Matplotlib State Management
# =============================================================================

_PLOTTING_STATE_LOCK = threading.RLock()
_ASYNC_PLOTTING_LOCK = anyio.Lock()


@contextmanager
def non_interactive_plotting() -> Generator[None, None, None]:
    """Temporarily suppress interactive display without switching backend.

    Use this when calling external plotting functions (e.g., cellrank, scvelo)
    that may call ``pyplot.show()``. Backend switching is intentionally avoided:
    Matplotlib 3.5 closes every open figure when switching backends. The
    caller's interactive mode and ``show`` function are restored on exit.

    Usage:
        with non_interactive_plotting():
            cr.pl.circular_projection(adata, ...)
            fig = plt.gcf()

    Yields:
        None
    """
    import matplotlib.pyplot as plt

    with _PLOTTING_STATE_LOCK:
        was_interactive = plt.isinteractive()
        original_show = plt.show
        try:
            plt.ioff()
            plt.show = lambda *_args, **_kwargs: None
            yield
        finally:
            plt.show = original_show
            if was_interactive:
                plt.ion()
            else:
                plt.ioff()


@contextmanager
def isolated_figure_scope() -> Generator[None, None, None]:
    """Close only Matplotlib figures created inside this context.

    This prevents a failing visualization or third-party plotting call from
    leaking figures without closing figures owned by another task.
    """
    import matplotlib.pyplot as plt

    existing_figures = set(plt.get_fignums())
    try:
        yield
    finally:
        created_figures = set(plt.get_fignums()) - existing_figures
        for figure_number in created_figures:
            plt.close(figure_number)


@asynccontextmanager
async def serialized_figure_scope() -> AsyncIterator[None]:
    """Serialize complete Matplotlib lifecycles across async tool calls.

    Matplotlib state is process-global. A thread lock cannot distinguish
    coroutines running on the same event-loop thread, so async callers must
    hold this lock from figure creation through export and cleanup.
    """
    async with _ASYNC_PLOTTING_LOCK:
        with isolated_figure_scope():
            yield


if TYPE_CHECKING:
    import matplotlib.pyplot as plt


# Standard savefig parameters for consistent figure output
SAVEFIG_PARAMS: dict[str, Any] = {
    "bbox_inches": "tight",
    "transparent": False,
    "facecolor": "white",
    "edgecolor": "none",
    "pad_inches": 0.1,
    "metadata": {"Software": "spatial-transcriptomics-mcp"},
}


# ============ Figure Export (Unified in visualize_data) ============


async def optimize_fig_to_image_with_cache(
    fig: "plt.Figure",
    params: Any,
    ctx: Optional["ToolContext"] = None,
    data_id: Optional[str] = None,
    plot_type: Optional[str] = None,
) -> str:
    """Save figure to file and return path.

    Supports user-specified output_path and output_format from VisualizationParameters.
    No separate save_visualization tool needed - all export options unified here.

    Args:
        fig: Matplotlib figure
        params: VisualizationParameters with dpi, output_path, output_format
        ctx: ToolContext (unused, kept for API compatibility)
        data_id: Dataset ID (unused, kept for API compatibility)
        plot_type: Plot type (used for filename generation)

    Returns:
        str with file path
    """
    import matplotlib.pyplot as plt

    # Extract parameters
    target_dpi = getattr(params, "dpi", None) or 100
    output_format = getattr(params, "output_format", None) or "png"
    output_path = getattr(params, "output_path", None)
    extra_artists = getattr(fig, "_chatspatial_extra_artists", None)

    # Determine output file path
    if output_path:
        # User specified path
        from pathlib import Path

        filepath = Path(output_path)
        # Ensure parent directory exists
        filepath.parent.mkdir(parents=True, exist_ok=True)
        # Ensure extension matches output_format
        expected_suffix = f".{output_format}"
        if not filepath.suffix:
            filepath = filepath.with_suffix(expected_suffix)
        elif filepath.suffix.lower() != expected_suffix.lower():
            # Override format to match user-specified extension
            output_format = filepath.suffix.lstrip(".").lower()
    else:
        # Default path
        output_dir = get_default_output_dir() / "visualizations"
        output_dir.mkdir(parents=True, exist_ok=True)
        filename = (
            f"{plot_type}_{uuid.uuid4().hex[:8]}.{output_format}"
            if plot_type
            else f"viz_{uuid.uuid4().hex[:8]}.{output_format}"
        )
        filepath = output_dir / filename

    # Prepare save parameters
    save_params = {
        **SAVEFIG_PARAMS,
        "dpi": target_dpi,
        "format": output_format,
        "bbox_extra_artists": extra_artists,
    }

    # Format-specific settings
    if output_format == "pdf":
        import matplotlib

        save_params["metadata"] = {
            "Title": f"{plot_type or 'visualization'}",
            "Creator": "ChatSpatial",
            "Producer": f"matplotlib {matplotlib.__version__}",
        }
    elif output_format in ["jpg", "jpeg"]:
        save_params["pil_kwargs"] = {"quality": 95}
        save_params.pop("metadata", None)
    elif output_format in ["svg", "eps", "tif", "tiff"]:
        save_params.pop("metadata", None)

    # Save figure
    fig.savefig(str(filepath), **save_params)
    actual_size = filepath.stat().st_size
    plt.close(fig)

    return (
        f"Visualization saved: {filepath}\n"
        f"Type: {plot_type or 'visualization'}\n"
        f"Format: {output_format.upper()}\n"
        f"Size: {actual_size // 1024}KB\n"
        f"Resolution: {target_dpi} DPI"
    )
