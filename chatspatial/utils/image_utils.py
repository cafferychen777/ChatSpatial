"""
Image utilities for spatial transcriptomics MCP.

This module provides standardized functions for handling images in the MCP.
All functions return Image objects that can be directly used in MCP tools.
"""

import io
import base64
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional, Union, Tuple
from mcp.server.fastmcp.utilities.types import Image


def fig_to_image(
    fig: plt.Figure,
    dpi: int = 100,
    format: str = 'png',
    max_size_kb: int = 900,
    close_fig: bool = True
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
            if current_format == 'jpg':
                try:
                    fig.savefig(
                        buf,
                        format=current_format,
                        dpi=current_dpi,
                        bbox_inches='tight',
                        transparent=False,
                        facecolor='white',
                        edgecolor='none',
                        pad_inches=0.1,
                        quality=85,  # Try with quality parameter
                        metadata={'Software': 'spatial-transcriptomics-mcp'}
                    )
                except TypeError:
                    # If quality parameter is not supported, try without it
                    fig.savefig(
                        buf,
                        format=current_format,
                        dpi=current_dpi,
                        bbox_inches='tight',
                        transparent=False,
                        facecolor='white',
                        edgecolor='none',
                        pad_inches=0.1,
                        metadata={'Software': 'spatial-transcriptomics-mcp'}
                    )
            else:
                fig.savefig(
                    buf,
                    format=current_format,
                    dpi=current_dpi,
                    bbox_inches='tight',
                    transparent=False,
                    facecolor='white',
                    edgecolor='none',
                    pad_inches=0.1,
                    metadata={'Software': 'spatial-transcriptomics-mcp'}
                )

            buf.seek(0)
            img_data = buf.read()
            size_kb = len(img_data) / 1024

            # If the image is small enough, return it
            if size_kb <= max_size_kb or (current_dpi <= min_dpi and current_format == 'jpg'):
                if close_fig:
                    plt.close(fig)
                # Return Image object with bytes data
                return Image(data=img_data, format=current_format)

            # Try different strategies to reduce size
            if current_format == 'png' and attempt == 1:
                # Switch to jpg format which is usually smaller
                current_format = 'jpg'
            else:
                # Reduce DPI more aggressively for subsequent attempts
                current_dpi = max(min_dpi, int(current_dpi * 0.6))  # Reduce DPI by 40%

        except Exception as e:
            print(f"Error saving figure: {str(e)}")
            # Try with a lower DPI
            current_dpi = max(min_dpi, int(current_dpi * 0.6))

    # If we still can't get a small enough image, return a simplified version
    if close_fig:
        plt.close(fig)

    # Create a simplified figure with minimal content
    simple_fig = plt.figure(figsize=(4, 3), dpi=72)
    ax = simple_fig.add_subplot(111)
    ax.text(0.5, 0.5, "Image too large - simplified version\nPlease try with fewer features or groups",
            ha='center', va='center', fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    buf = io.BytesIO()
    # Use PNG format instead of JPG for Claude compatibility
    simple_fig.savefig(buf, format='png', dpi=72, bbox_inches='tight')

    buf.seek(0)
    img_data = buf.read()
    plt.close(simple_fig)

    # Return Image object with bytes data
    return Image(data=img_data, format='png')


def create_placeholder_image(
    message: str = "No visualization available",
    figsize: Tuple[int, int] = (6, 6),
    format: str = 'png'
) -> Image:
    """Create a placeholder image with a message

    Args:
        message: Message to display in the placeholder image
        figsize: Figure size in inches
        format: Image format (png, jpg)

    Returns:
        Image object
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # Convert to Image object
    return fig_to_image(fig, format=format)


def fig_to_base64(
    fig: plt.Figure,
    dpi: int = 100,
    format: str = 'png',
    max_size_mb: float = 5,
    close_fig: bool = True
) -> str:
    """Convert matplotlib figure to base64 encoded string

    This function is provided for backward compatibility.
    New code should use fig_to_image instead.

    Args:
        fig: Matplotlib figure
        dpi: Resolution in dots per inch
        format: Image format (png, jpg)
        max_size_mb: Maximum size in MB for the image
        close_fig: Whether to close the figure after conversion

    Returns:
        Base64 encoded string of the image
    """
    # Convert to Image object first
    image = fig_to_image(fig, dpi=dpi, format=format, max_size_kb=max_size_mb*1024, close_fig=close_fig)

    # Extract base64 string from Image object
    if isinstance(image.data, bytes):
        return base64.b64encode(image.data).decode('utf-8')
    else:
        # This should not happen, but just in case
        buf = io.BytesIO()
        fig.savefig(buf, format=format, dpi=dpi)
        buf.seek(0)
        if close_fig:
            plt.close(fig)
        return base64.b64encode(buf.read()).decode('utf-8')


