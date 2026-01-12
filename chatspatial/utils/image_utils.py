"""
Image utilities for spatial transcriptomics MCP.

This module provides standardized functions for handling images in the MCP.
All functions return Image objects that can be directly used in MCP tools.
"""

import base64
import io
import os
import uuid
from typing import TYPE_CHECKING, Any, Dict, Optional

from mcp.types import ImageContent

from .exceptions import ProcessingError

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


# Function to ensure matplotlib uses non-interactive backend
def _ensure_non_interactive_backend():
    """Ensure matplotlib uses non-interactive backend to prevent GUI popups on macOS."""
    import matplotlib

    current_backend = matplotlib.get_backend()
    if current_backend != "Agg":
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.ioff()  # Turn off interactive mode


if TYPE_CHECKING:
    import matplotlib.pyplot as plt


# Standard savefig parameters for consistent figure output
SAVEFIG_PARAMS: Dict[str, Any] = {
    "bbox_inches": "tight",
    "transparent": False,
    "facecolor": "white",
    "edgecolor": "none",
    "pad_inches": 0.1,
    "metadata": {"Software": "spatial-transcriptomics-mcp"},
}


def bytes_to_image_content(data: bytes, format: str = "png") -> ImageContent:
    """Convert raw image bytes to MCP ImageContent.

    This unified utility function handles the conversion from raw image bytes
    to the MCP-compatible ImageContent type with proper MIME type mapping.

    Args:
        data: Raw image bytes
        format: Image format (png, jpg, jpeg, gif, webp)

    Returns:
        ImageContent object ready for MCP tool return

    Examples:
        >>> img_bytes = fig.savefig(buf, format='png')
        >>> content = bytes_to_image_content(img_bytes, format='png')
    """
    # MIME type mapping for common image formats
    format_to_mime = {
        "png": "image/png",
        "jpg": "image/jpeg",
        "jpeg": "image/jpeg",
        "gif": "image/gif",
        "webp": "image/webp",
    }

    # Get MIME type, default to PNG if format is unknown
    mime_type = format_to_mime.get(format.lower(), "image/png")

    # Encode to base64 string as required by ImageContent
    encoded_data = base64.b64encode(data).decode("utf-8")

    return ImageContent(type="image", data=encoded_data, mimeType=mime_type)


def fig_to_image(
    fig: "plt.Figure",
    dpi: int = 100,
    format: str = "png",
    close_fig: bool = True,
) -> ImageContent:
    """Convert matplotlib figure to ImageContent

    This function respects user's DPI and format settings without any
    automatic compression or quality reduction. Large images are handled
    by optimize_fig_to_image_with_cache which saves them to disk.

    Args:
        fig: Matplotlib figure
        dpi: Resolution in dots per inch (user's setting is always respected)
        format: Image format (png or jpg)
        close_fig: Whether to close the figure after conversion

    Returns:
        ImageContent object ready for MCP tool return
    """
    _ensure_non_interactive_backend()  # Prevent GUI popups on macOS
    import matplotlib.pyplot as plt

    buf = io.BytesIO()

    # Save figure with user's exact settings - no compromise
    # Check for extra artists (e.g., legends positioned outside plot area)
    extra_artists = getattr(fig, "_chatspatial_extra_artists", None)

    try:
        if format == "jpg":
            try:
                # Try with quality parameter first (newer matplotlib)
                fig.savefig(
                    buf,
                    format=format,
                    dpi=dpi,
                    bbox_extra_artists=extra_artists,
                    quality=85,
                    **SAVEFIG_PARAMS,
                )
            except TypeError:
                # Fallback for older matplotlib without quality parameter
                fig.savefig(
                    buf,
                    format=format,
                    dpi=dpi,
                    bbox_extra_artists=extra_artists,
                    **SAVEFIG_PARAMS,
                )
        else:  # PNG
            fig.savefig(
                buf,
                format=format,
                dpi=dpi,
                bbox_extra_artists=extra_artists,
                **SAVEFIG_PARAMS,
            )

        buf.seek(0)
        img_data = buf.read()

        if close_fig:
            plt.close(fig)

        # Convert to ImageContent using unified utility
        return bytes_to_image_content(img_data, format=format)

    except Exception as e:
        if close_fig:
            plt.close(fig)
        raise ProcessingError(f"Failed to convert figure to image: {str(e)}") from e


# ============ Token Optimization and Publication Export Support ============


async def optimize_fig_to_image_with_cache(
    fig: "plt.Figure",
    params: Any,
    ctx: Optional["ToolContext"] = None,
    data_id: Optional[str] = None,
    plot_type: Optional[str] = None,
    mode: str = "auto",  # Kept for API compatibility, ignored
) -> str:
    """Save figure to file and return path for MCP protocol efficiency.

    MCP protocol has ~3.3x overhead for embedded ImageContent, so we always
    save to file and return path as text (following MCP best practice:
    "Prefer using URIs over embedded content for large files").

    Args:
        fig: Matplotlib figure
        params: Visualization parameters (used for DPI extraction)
        ctx: ToolContext (unused, kept for API compatibility)
        data_id: Dataset ID (unused, kept for API compatibility)
        plot_type: Plot type (used for filename generation)
        mode: Ignored (always saves to file)

    Returns:
        str with file path (FastMCP auto-converts to TextContent)
    """
    _ensure_non_interactive_backend()
    import matplotlib.pyplot as plt

    target_dpi = params.dpi if hasattr(params, "dpi") and params.dpi else 100
    extra_artists = getattr(fig, "_chatspatial_extra_artists", None)

    # Generate image
    buf = io.BytesIO()
    fig.savefig(
        buf,
        format="png",
        dpi=target_dpi,
        bbox_extra_artists=extra_artists,
        **SAVEFIG_PARAMS,
    )
    actual_size = buf.tell()

    # Save to file
    os.makedirs("/tmp/chatspatial/visualizations", exist_ok=True)
    filename = (
        f"{plot_type}_{uuid.uuid4().hex[:8]}.png"
        if plot_type
        else f"viz_{uuid.uuid4().hex[:8]}.png"
    )
    filepath = f"/tmp/chatspatial/visualizations/{filename}"

    buf.seek(0)
    with open(filepath, "wb") as f:
        f.write(buf.read())
    buf.close()
    plt.close(fig)

    return (
        f"Visualization saved: {filepath}\n"
        f"Type: {plot_type or 'visualization'}\n"
        f"Size: {actual_size // 1024}KB\n"
        f"Resolution: {target_dpi} DPI"
    )
