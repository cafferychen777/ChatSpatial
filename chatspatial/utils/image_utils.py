"""
Image utilities for spatial transcriptomics MCP.

This module provides standardized functions for handling images in the MCP.
All functions return Image objects that can be directly used in MCP tools.
"""

import base64
import io
import json
import os
import pickle
import sys
import uuid
import weakref
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple, Union

from mcp.types import EmbeddedResource, ImageContent, TextResourceContents


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
    try:
        if format == "jpg":
            try:
                # Try with quality parameter first (newer matplotlib)
                fig.savefig(
                    buf,
                    format=format,
                    dpi=dpi,
                    bbox_inches="tight",
                    transparent=False,
                    facecolor="white",
                    edgecolor="none",
                    pad_inches=0.1,
                    quality=85,
                    metadata={"Software": "spatial-transcriptomics-mcp"},
                )
            except TypeError:
                # Fallback for older matplotlib without quality parameter
                fig.savefig(
                    buf,
                    format=format,
                    dpi=dpi,
                    bbox_inches="tight",
                    transparent=False,
                    facecolor="white",
                    edgecolor="none",
                    pad_inches=0.1,
                    metadata={"Software": "spatial-transcriptomics-mcp"},
                )
        else:  # PNG
            fig.savefig(
                buf,
                format=format,
                dpi=dpi,
                bbox_inches="tight",
                transparent=False,
                facecolor="white",
                edgecolor="none",
                pad_inches=0.1,
                metadata={"Software": "spatial-transcriptomics-mcp"},
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
        raise RuntimeError(f"Failed to convert figure to image: {str(e)}") from e


def create_placeholder_image(
    message: str = "No visualization available",
    figsize: Tuple[int, int] = (6, 6),
    format: str = "png",
) -> ImageContent:
    """Create a placeholder image with a message

    Args:
        message: Message to display in the placeholder image
        figsize: Figure size in inches
        format: Image format (png, jpg)

    Returns:
        ImageContent object ready for MCP tool return
    """
    _ensure_non_interactive_backend()  # Prevent GUI popups on macOS
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize)
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Convert to Image object
    return fig_to_image(fig, format=format)


# ============ Token Optimization and Publication Export Support ============

# Global Figure cache (using weak references to avoid memory leaks)
_figure_cache: Dict[str, weakref.ReferenceType] = {}


def cache_figure(key: str, fig: "plt.Figure"):
    """Cache matplotlib figure object for high-quality export

    Args:
        key: Cache key (usually data_id_plot_type)
        fig: Matplotlib figure to cache
    """
    _figure_cache[key] = weakref.ref(fig)


def get_cached_figure(key: str) -> Optional["plt.Figure"]:
    """Get cached figure object

    Args:
        key: Cache key

    Returns:
        Cached figure or None if not found/expired
    """
    if key in _figure_cache:
        fig_ref = _figure_cache[key]
        fig = fig_ref()
        if fig is not None:
            return fig
    return None


def save_figure_pickle(fig: "plt.Figure", path: str):
    """Save figure object to pickle file (preserves all information)

    Args:
        fig: Matplotlib figure
        path: Path to save pickle file
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(fig, f)


def load_figure_pickle(path: str) -> "plt.Figure":
    """Load figure object from pickle file

    Args:
        path: Path to pickle file

    Returns:
        Loaded figure object
    """
    with open(path, "rb") as f:
        return pickle.load(f)


async def optimize_fig_to_image_with_cache(
    fig: "plt.Figure",
    params: Any,
    context: Optional[Any] = None,
    data_id: str = None,
    plot_type: str = None,
    mode: str = "auto",
) -> Union[ImageContent, str]:
    """Optimized image conversion with Figure caching for high-quality export

    This function implements MCP 2025 best practice token optimization:
    - Small images (<70KB): Direct embedding as ImageContent
    - Large images (≥70KB): Save to file, return path as text (URI over embedded content)
    - Caches Figure object for later high-quality export

    Following MCP specification recommendation:
    "Prefer using URIs over embedded content for large files"

    Args:
        fig: Matplotlib figure
        params: Visualization parameters
        context: MCP context for logging
        data_id: Dataset ID (for cache key)
        plot_type: Plot type (for cache key)
        mode: Optimization mode - "auto" or "direct"

    Returns:
        Small images: ImageContent object (embedded)
        Large images: str with file path (FastMCP auto-converts to TextContent)
    """
    _ensure_non_interactive_backend()  # Prevent GUI popups on macOS
    import matplotlib.pyplot as plt

    # Initialize variables
    cache_key = None
    pickle_path = None

    # Cache Figure object for high-quality export
    if data_id and plot_type:
        cache_key = f"{data_id}_{plot_type}"
        cache_figure(cache_key, fig)

        # Also save pickle file for persistence
        os.makedirs("/tmp/chatspatial/figures", exist_ok=True)
        pickle_path = f"/tmp/chatspatial/figures/{cache_key}.pkl"
        save_figure_pickle(fig, pickle_path)

        if context:
            await context.info(
                f"Cached figure object for high-quality export: {cache_key}"
            )

    # Estimate original image size
    test_buf = io.BytesIO()
    target_dpi = params.dpi if hasattr(params, "dpi") and params.dpi else 100
    fig.savefig(test_buf, format="png", dpi=target_dpi, bbox_inches="tight")
    estimated_size = test_buf.tell()
    test_buf.close()

    # MCP 2025 best practice: prefer URIs over embedded content for large files
    # Threshold: 70KB (safe for MCP 25K token limit)
    # Calculation: 70KB × 1.33 (base64) ÷ 4 (chars/token) ≈ 23K tokens
    DIRECT_EMBED_THRESHOLD = 70 * 1024  # 70KB

    # Small images: Direct embedding
    if mode == "direct" or (mode == "auto" and estimated_size < DIRECT_EMBED_THRESHOLD):
        if context:
            await context.info(
                f"Small image ({estimated_size//1024}KB), embedding directly"
            )
        return fig_to_image(fig, dpi=target_dpi, format="png")

    # Large images: Save to file, return path as text
    # This follows MCP best practice and avoids token limits
    if context:
        await context.info(
            f"Large image ({estimated_size//1024}KB), saving to file "
            f"(following MCP best practice: URI over embedded content)"
        )

    # Save high-quality version
    os.makedirs("/tmp/chatspatial/visualizations", exist_ok=True)
    hq_filename = (
        f"{plot_type}_{uuid.uuid4().hex[:8]}.png"
        if plot_type
        else f"viz_{uuid.uuid4().hex[:8]}.png"
    )
    hq_path = f"/tmp/chatspatial/visualizations/{hq_filename}"

    fig.savefig(
        hq_path,
        dpi=target_dpi if target_dpi else 300,
        format="png",
        bbox_inches="tight",
        facecolor="white",
    )

    # Close figure
    plt.close(fig)

    # NOTE: Metadata is now managed by server.py, not here
    # Pickle file serves as source of truth for figure reconstruction

    # Return text message with file path (MCP best practice for large files)
    # FastMCP will auto-convert str to TextContent
    message = (
        f"✓ Visualization created successfully!\n\n"
        f"📊 **{plot_type if plot_type else 'Visualization'}** "
        f"({estimated_size//1024}KB estimated)\n\n"
        f"📁 **High-quality image saved to:**\n`{hq_path}`\n\n"
        f"🎨 Resolution: {target_dpi if target_dpi else 300} DPI\n"
        f"📦 Format: PNG\n\n"
        f"💡 Following MCP 2025 best practice: large images returned as file paths "
        f"to avoid token limits and reduce conversation costs.\n\n"
        f"You can view this image using your system's image viewer, "
        f"or use the export/save tools to convert to other formats."
    )

    if context:
        await context.info(f"✓ Saved high-quality image to: {hq_path}")

    return message
