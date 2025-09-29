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
from pathlib import Path
from typing import TYPE_CHECKING, Tuple, Optional, Union, Dict, Any

from mcp.server.fastmcp.utilities.types import Image
from mcp.types import EmbeddedResource, TextResourceContents

# Function to ensure matplotlib uses non-interactive backend
def _ensure_non_interactive_backend():
    """Ensure matplotlib uses non-interactive backend to prevent GUI popups on macOS."""
    import matplotlib
    current_backend = matplotlib.get_backend()
    if current_backend != 'Agg':
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.ioff()  # Turn off interactive mode

if TYPE_CHECKING:
    import matplotlib.pyplot as plt


def fig_to_image(
    fig: "plt.Figure",
    dpi: int = 100,
    format: str = "png",
    max_size_kb: int = 900,
    close_fig: bool = True,
) -> Image:
    """Convert matplotlib figure to Image object with size control

    Args:
        fig: Matplotlib figure
        dpi: Resolution in dots per inch (lower = smaller file)
        format: Image format (png, jpg)
        max_size_kb: Maximum size in KB for the image (MCP limit is 1MB)
        close_fig: Whether to close the figure after conversion

    Returns:
        Image object
    """
    _ensure_non_interactive_backend()  # Prevent GUI popups on macOS
    import matplotlib.pyplot as plt

    # Try different compression settings until we get a small enough image
    current_dpi = dpi
    min_dpi = 40  # Lower minimum DPI to allow smaller files

    # For heatmaps and other complex plots, start with a lower DPI
    if fig.get_size_inches()[0] > 8 or fig.get_size_inches()[1] > 8:
        current_dpi = min(current_dpi, 80)  # Start with lower DPI for large figures

    # Try JPEG format for potentially smaller file size if PNG is too large
    current_format = format

    for attempt in range(4):  # Try up to 4 times with different settings
        buf = io.BytesIO()

        try:
            # Save figure with appropriate settings
            # Remove quality parameter which is not supported in some matplotlib versions
            if current_format == "jpg":
                try:
                    fig.savefig(
                        buf,
                        format=current_format,
                        dpi=current_dpi,
                        bbox_inches="tight",
                        transparent=False,
                        facecolor="white",
                        edgecolor="none",
                        pad_inches=0.1,
                        quality=85,  # Try with quality parameter
                        metadata={"Software": "spatial-transcriptomics-mcp"},
                    )
                except TypeError:
                    # If quality parameter is not supported, try without it
                    fig.savefig(
                        buf,
                        format=current_format,
                        dpi=current_dpi,
                        bbox_inches="tight",
                        transparent=False,
                        facecolor="white",
                        edgecolor="none",
                        pad_inches=0.1,
                        metadata={"Software": "spatial-transcriptomics-mcp"},
                    )
            else:
                fig.savefig(
                    buf,
                    format=current_format,
                    dpi=current_dpi,
                    bbox_inches="tight",
                    transparent=False,
                    facecolor="white",
                    edgecolor="none",
                    pad_inches=0.1,
                    metadata={"Software": "spatial-transcriptomics-mcp"},
                )

            buf.seek(0)
            img_data = buf.read()
            size_kb = len(img_data) / 1024

            # If the image is small enough, return it
            if size_kb <= max_size_kb or (
                current_dpi <= min_dpi and current_format == "jpg"
            ):
                if close_fig:
                    plt.close(fig)
                # Return Image object with bytes data
                return Image(data=img_data, format=current_format)

            # Try different strategies to reduce size
            if current_format == "png" and attempt == 1:
                # Switch to jpg format which is usually smaller
                current_format = "jpg"
            else:
                # Reduce DPI more aggressively for subsequent attempts
                current_dpi = max(min_dpi, int(current_dpi * 0.6))  # Reduce DPI by 40%

        except Exception as e:
            print(f"Error saving figure: {str(e)}", file=sys.stderr)
            # Try with a lower DPI
            current_dpi = max(min_dpi, int(current_dpi * 0.6))

    # If we still can't get a small enough image, return a simplified version
    if close_fig:
        plt.close(fig)

    # Create a simplified figure with minimal content
    simple_fig = plt.figure(figsize=(4, 3), dpi=72)
    ax = simple_fig.add_subplot(111)
    ax.text(
        0.5,
        0.5,
        "Image too large - simplified version\nPlease try with fewer features or groups",
        ha="center",
        va="center",
        fontsize=10,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    buf = io.BytesIO()
    # Use PNG format instead of JPG for Claude compatibility
    simple_fig.savefig(buf, format="png", dpi=72, bbox_inches="tight")

    buf.seek(0)
    img_data = buf.read()
    plt.close(simple_fig)

    # Return Image object with bytes data
    return Image(data=img_data, format="png")


def create_placeholder_image(
    message: str = "No visualization available",
    figsize: Tuple[int, int] = (6, 6),
    format: str = "png",
) -> Image:
    """Create a placeholder image with a message

    Args:
        message: Message to display in the placeholder image
        figsize: Figure size in inches
        format: Image format (png, jpg)

    Returns:
        Image object
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
    with open(path, 'wb') as f:
        pickle.dump(fig, f)


def load_figure_pickle(path: str) -> "plt.Figure":
    """Load figure object from pickle file
    
    Args:
        path: Path to pickle file
        
    Returns:
        Loaded figure object
    """
    with open(path, 'rb') as f:
        return pickle.load(f)




async def optimize_fig_to_image_with_cache(
    fig: "plt.Figure",
    params: Any,
    context: Optional[Any] = None,
    data_id: str = None,
    plot_type: str = None,
    mode: str = "auto"
) -> Union[Image, Tuple[Image, EmbeddedResource]]:
    """Optimized image conversion with Figure caching for high-quality export
    
    This function implements the token optimization strategy:
    - Small images (<100KB): Direct embedding
    - Large images (>=100KB): Preview + resource reference
    - Caches Figure object for later high-quality export
    
    Args:
        fig: Matplotlib figure
        params: Visualization parameters
        context: MCP context for logging
        data_id: Dataset ID (for cache key)
        plot_type: Plot type (for cache key) 
        mode: Optimization mode - "auto", "preview", or "direct"
        
    Returns:
        Small images: Image object
        Large images: Tuple[Image preview, EmbeddedResource with high quality]
    """
    _ensure_non_interactive_backend()  # Prevent GUI popups on macOS
    import matplotlib.pyplot as plt
    
    # Cache Figure object for high-quality export
    if data_id and plot_type:
        cache_key = f"{data_id}_{plot_type}"
        cache_figure(cache_key, fig)
        
        # Also save pickle file for persistence
        os.makedirs("/tmp/chatspatial/figures", exist_ok=True)
        pickle_path = f"/tmp/chatspatial/figures/{cache_key}.pkl"
        save_figure_pickle(fig, pickle_path)
        
        if context:
            await context.info(f"Cached figure object for high-quality export: {cache_key}")
    
    # Estimate original image size
    test_buf = io.BytesIO()
    target_dpi = params.dpi if hasattr(params, 'dpi') and params.dpi else 100
    fig.savefig(test_buf, format='png', dpi=target_dpi, bbox_inches='tight')
    estimated_size = test_buf.tell()
    test_buf.close()
    
    # Decision thresholds
    DIRECT_EMBED_THRESHOLD = 100 * 1024  # 100KB
    PREVIEW_THRESHOLD = 500 * 1024       # 500KB
    
    # Mode: direct - force direct embedding
    if mode == "direct" or (mode == "auto" and estimated_size < DIRECT_EMBED_THRESHOLD):
        if context:
            await context.info(f"Small image ({estimated_size//1024}KB), embedding directly")
        return fig_to_image(fig, dpi=target_dpi, format="png", max_size_kb=900)
    
    # Mode: preview - use preview+resource strategy
    if mode == "preview" or (mode == "auto" and estimated_size > PREVIEW_THRESHOLD):
        if context:
            await context.info(f"Large image ({estimated_size//1024}KB), using preview+resource strategy")
        
        # 1. Create low-quality preview (target: 50KB)
        preview_buf = io.BytesIO()
        try:
            # Try with quality parameter
            fig.savefig(
                preview_buf,
                format='jpeg',
                dpi=60,
                quality=40,
                bbox_inches='tight',
                facecolor='white'
            )
        except TypeError:
            # If quality parameter is not supported, try without it
            fig.savefig(
                preview_buf,
                format='jpeg',
                dpi=60,
                bbox_inches='tight',
                facecolor='white'
            )
        preview_buf.seek(0)
        preview_image = Image(data=preview_buf.read(), format="jpeg")
        
        # 2. Save high-quality version
        os.makedirs("/tmp/chatspatial/visualizations", exist_ok=True)
        hq_filename = f"{plot_type}_{uuid.uuid4().hex[:8]}.png" if plot_type else f"viz_{uuid.uuid4().hex[:8]}.png"
        hq_path = f"/tmp/chatspatial/visualizations/{hq_filename}"
        
        fig.savefig(
            hq_path,
            dpi=target_dpi if target_dpi else 300,
            format='png',
            bbox_inches='tight',
            facecolor='white'
        )
        
        # 3. Create resource reference with metadata
        metadata = {
            "description": f"High-quality {plot_type if plot_type else 'visualization'} at {target_dpi} DPI",
            "figure_pickle": pickle_path if data_id and plot_type else None,
            "can_export_pdf": True,
            "can_export_svg": True,
            "cache_key": cache_key if data_id and plot_type else None
        }
        
        resource = EmbeddedResource(
            type="resource",
            resource=TextResourceContents(
                uri=f"file://{hq_path}",
                mimeType="image/png",
                text=json.dumps(metadata)
            )
        )
        
        if context:
            await context.info(
                f"✅ Preview: {len(preview_image.data)//1024}KB | "
                f"High-quality: {hq_path} | "
                f"Figure cached for PDF/SVG export"
            )
        
        return preview_image, resource
    
    # Mode: auto with medium size - compress but don't use preview
    if context:
        await context.info(f"Medium image ({estimated_size//1024}KB), using compression")
    return fig_to_image(fig, dpi=80, format="png", max_size_kb=200)


