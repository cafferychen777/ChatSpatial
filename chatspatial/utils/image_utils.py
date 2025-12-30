"""
Image utilities for spatial transcriptomics MCP.

This module provides standardized functions for handling images in the MCP.
All functions return Image objects that can be directly used in MCP tools.
"""

import base64
import io
import json
import os
import uuid
import weakref
from typing import TYPE_CHECKING, Any, Dict, Optional, Union

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


def save_visualization_metadata(
    path: str,
    data_id: str,
    plot_type: str,
    params: Any,
) -> bool:
    """Save visualization metadata as JSON for later regeneration.

    This replaces unsafe pickle serialization with a secure JSON-based approach.
    Instead of storing the matplotlib figure object, we store the parameters
    needed to regenerate it. This is fundamentally more secure because:
    1. JSON cannot contain executable code (unlike pickle)
    2. Parameters are human-readable and auditable
    3. Regeneration uses the original trusted codebase

    Args:
        path: Path to save JSON metadata file
        data_id: Dataset identifier
        plot_type: Type of plot (with optional subtype suffix)
        params: VisualizationParameters object or dict

    Returns:
        True if saved successfully
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)

    # Convert params to dict if it's a Pydantic model
    if hasattr(params, "model_dump"):
        params_dict = params.model_dump()
    elif hasattr(params, "dict"):
        params_dict = params.dict()  # Pydantic v1 fallback
    elif isinstance(params, dict):
        params_dict = params
    else:
        params_dict = {}

    metadata = {
        "data_id": data_id,
        "plot_type": plot_type,
        "params": params_dict,
        "version": "1.0",  # For future compatibility
    }

    with open(path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, default=str)

    return True


def load_visualization_metadata(path: str) -> Dict[str, Any]:
    """Load visualization metadata from JSON file.

    Args:
        path: Path to JSON metadata file

    Returns:
        Dictionary containing data_id, plot_type, and params

    Raises:
        FileNotFoundError: If metadata file doesn't exist
        json.JSONDecodeError: If file is not valid JSON
    """
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


async def optimize_fig_to_image_with_cache(
    fig: "plt.Figure",
    params: Any,
    ctx: Optional["ToolContext"] = None,
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
        ctx: ToolContext for logging and data access
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

    # Cache Figure object for in-session high-quality export (WeakRef only)
    # For cross-session persistence, we save JSON metadata (not pickle)
    if data_id and plot_type:
        cache_key = f"{data_id}_{plot_type}"
        cache_figure(cache_key, fig)

        # Save visualization metadata as JSON (secure, no pickle)
        # This allows regenerating the figure in future sessions
        os.makedirs("/tmp/chatspatial/figures", exist_ok=True)
        metadata_path = f"/tmp/chatspatial/figures/{cache_key}.json"
        save_visualization_metadata(metadata_path, data_id, plot_type, params)

    # Generate image once with ACTUAL parameters (not estimation)
    target_dpi = params.dpi if hasattr(params, "dpi") and params.dpi else 100

    # Check for extra artists (e.g., legends positioned outside plot area)
    extra_artists = getattr(fig, "_chatspatial_extra_artists", None)

    # Use fig_to_image to get actual image with full parameters
    actual_buf = io.BytesIO()
    fig.savefig(
        actual_buf,
        format="png",
        dpi=target_dpi,
        bbox_extra_artists=extra_artists,
        **SAVEFIG_PARAMS,
    )
    actual_size = actual_buf.tell()

    # MCP 2025 best practice: prefer URIs over embedded content for large files
    # Root cause identified: MCP protocol has ~3.3x overhead for ImageContent
    # Even 16K token images become 53K tokens in MCP responses (beyond our control)
    # Solution: Always use file URIs to avoid MCP ImageContent overhead
    DIRECT_EMBED_THRESHOLD = 0  # 0 = Always save to file (avoids MCP overhead)

    # Small images: Direct embedding (use already-generated image)
    if mode == "direct" or (mode == "auto" and actual_size < DIRECT_EMBED_THRESHOLD):
        # Convert buffer to ImageContent
        actual_buf.seek(0)
        img_data = actual_buf.read()
        actual_buf.close()
        import matplotlib.pyplot as plt

        plt.close(fig)
        return bytes_to_image_content(img_data, format="png")

    # Large images: Save to file, return path as text
    # This follows MCP best practice and avoids token limits
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

    # NOTE: JSON metadata saved above enables figure regeneration
    # Data + params are the source of truth, not serialized figure objects

    # Return text message with file path (MCP best practice for large files)
    # FastMCP will auto-convert str to TextContent
    message = (
        f"Visualization saved: {hq_path}\n"
        f"Type: {plot_type if plot_type else 'visualization'}\n"
        f"Size: {actual_size//1024}KB\n"
        f"Resolution: {target_dpi if target_dpi else 300} DPI"
    )

    return message
