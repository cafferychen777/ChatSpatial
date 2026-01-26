"""
Image utilities for spatial transcriptomics MCP.

This module provides:
- Matplotlib backend management (prevent GUI popups)
- Figure export to file with user-specified format/DPI
"""

import uuid
from collections.abc import Generator
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any, Optional

from .path_utils import get_safe_output_path

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


# =============================================================================
# Matplotlib Backend Management
# =============================================================================


def _ensure_non_interactive_backend() -> None:
    """Ensure matplotlib uses non-interactive backend to prevent GUI popups on macOS."""
    import matplotlib

    current_backend = matplotlib.get_backend()
    if current_backend != "Agg":
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.ioff()


@contextmanager
def non_interactive_backend() -> Generator[None, None, None]:
    """Context manager for temporary non-interactive matplotlib backend.

    Use this when calling external plotting functions (e.g., cellrank, scvelo)
    that may trigger interactive backend behavior. The original backend is
    restored after the context exits.

    For internal plotting code that doesn't need backend restoration,
    use _ensure_non_interactive_backend() instead.

    Usage:
        with non_interactive_backend():
            cr.pl.circular_projection(adata, ...)
            fig = plt.gcf()

    Yields:
        None
    """
    import matplotlib

    original_backend = matplotlib.get_backend()
    needs_restore = original_backend != "Agg"

    try:
        if needs_restore:
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            plt.ioff()
        yield
    finally:
        if needs_restore:
            matplotlib.use(original_backend)


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
    _ensure_non_interactive_backend()
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
        # Ensure correct extension
        if not filepath.suffix:
            filepath = filepath.with_suffix(f".{output_format}")
    else:
        # Default path
        output_dir = get_safe_output_path("./visualizations")
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
    elif output_format in ["svg", "eps"]:
        # SVG and EPS writers don't support "Software" metadata key, remove it
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
