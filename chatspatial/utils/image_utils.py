"""
Image utilities for spatial transcriptomics MCP.

This module provides standardized functions for handling images in the MCP.
All functions return Image objects that can be directly used in MCP tools.
"""

import base64
import io
import os
import pickle
import uuid
import weakref
from typing import TYPE_CHECKING, Any, Dict, Optional, Union

from mcp.types import ImageContent


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


def save_figure_pickle(fig: "plt.Figure", path: str) -> bool:
    """Save figure object to pickle file (preserves all information)

    Args:
        fig: Matplotlib figure
        path: Path to save pickle file

    Returns:
        True if saved successfully, False if pickle failed (e.g., networkx figures)
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    try:
        with open(path, "wb") as f:
            pickle.dump(fig, f)
        return True
    except (AttributeError, TypeError) as e:
        # Some figures (e.g., networkx circle_plot) cannot be pickled
        # This is not a fatal error - just skip caching
        import logging
        logging.debug(f"Figure pickle failed (non-fatal): {e}")
        return False


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
    data_id: Optional[str] = None,
    plot_type: Optional[str] = None,
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

    # Generate image once with ACTUAL parameters (not estimation)
    target_dpi = params.dpi if hasattr(params, "dpi") and params.dpi else 100

    # Use fig_to_image to get actual image with full parameters
    actual_buf = io.BytesIO()
    fig.savefig(
        actual_buf,
        format="png",
        dpi=target_dpi,
        bbox_inches="tight",
        transparent=False,
        facecolor="white",
        edgecolor="none",
        pad_inches=0.1,
        metadata={"Software": "spatial-transcriptomics-mcp"},
    )
    actual_size = actual_buf.tell()

    # MCP 2025 best practice: prefer URIs over embedded content for large files
    # Root cause identified: MCP protocol has ~3.3x overhead for ImageContent
    # Even 16K token images become 53K tokens in MCP responses (beyond our control)
    # Solution: Always use file URIs to avoid MCP ImageContent overhead
    DIRECT_EMBED_THRESHOLD = 0  # 0 = Always save to file (avoids MCP overhead)

    # Small images: Direct embedding (use already-generated image)
    if mode == "direct" or (mode == "auto" and actual_size < DIRECT_EMBED_THRESHOLD):
        if context:
            await context.info(
                f"Small image ({actual_size//1024}KB), embedding directly"
            )
        # Convert buffer to ImageContent
        actual_buf.seek(0)
        img_data = actual_buf.read()
        actual_buf.close()
        import matplotlib.pyplot as plt

        plt.close(fig)
        return bytes_to_image_content(img_data, format="png")

    # Large images: Save to file, return path as text
    # This follows MCP best practice and avoids token limits
    if context:
        await context.info(
            f"Large image ({actual_size//1024}KB), saving to file "
            f"(following MCP best practice: URI over embedded content)"
        )

    # Save the already-generated image to file (avoid regenerating)
    os.makedirs("/tmp/chatspatial/visualizations", exist_ok=True)
    hq_filename = (
        f"{plot_type}_{uuid.uuid4().hex[:8]}.png"
        if plot_type
        else f"viz_{uuid.uuid4().hex[:8]}.png"
    )
    hq_path = f"/tmp/chatspatial/visualizations/{hq_filename}"

    # Write the image data we already generated
    actual_buf.seek(0)
    with open(hq_path, "wb") as f:
        f.write(actual_buf.read())
    actual_buf.close()

    # Close figure
    import matplotlib.pyplot as plt

    plt.close(fig)

    # NOTE: Metadata is now managed by server.py, not here
    # Pickle file serves as source of truth for figure reconstruction

    # Return text message with file path (MCP best practice for large files)
    # FastMCP will auto-convert str to TextContent
    message = (
        f"Visualization saved: {hq_path}\n"
        f"Type: {plot_type if plot_type else 'visualization'}\n"
        f"Size: {actual_size//1024}KB\n"
        f"Resolution: {target_dpi if target_dpi else 300} DPI"
    )

    if context:
        await context.info(f"Saved visualization to: {hq_path}")

    return message
