"""
Visualization persistence functions.

This module contains:
- save_visualization: Save single visualization to disk
- export_all_visualizations: Export all cached visualizations
- clear_visualization_cache: Clear visualization cache
"""

from datetime import datetime
from typing import TYPE_CHECKING, Any, Dict, List, Optional

import matplotlib
import matplotlib.pyplot as plt

from ...models.data import VisualizationParameters
from ...utils.exceptions import DataNotFoundError, ParameterError, ProcessingError
from ...utils.path_utils import get_output_dir_from_config, get_safe_output_path

if TYPE_CHECKING:
    import anndata as ad

    from ...spatial_mcp_adapter import ToolContext


# =============================================================================
# Internal Helper Functions
# =============================================================================


async def _regenerate_figure_for_export(
    adata: "ad.AnnData",
    params: VisualizationParameters,
    context: Optional["ToolContext"] = None,
) -> plt.Figure:
    """Regenerate a matplotlib figure from saved parameters for high-quality export.

    This is an internal helper function used by save_visualization to recreate
    figures from JSON metadata. It directly returns the matplotlib Figure object
    (instead of ImageContent) so it can be exported at arbitrary DPI/format.

    Args:
        adata: AnnData object containing the data
        params: VisualizationParameters reconstructed from saved metadata
        context: MCP context for logging

    Returns:
        Matplotlib Figure object ready for export

    Raises:
        ValueError: If plot_type is unknown
    """
    # Import here to avoid circular imports
    from . import PLOT_HANDLERS

    plot_type = params.plot_type

    if plot_type not in PLOT_HANDLERS:
        raise ParameterError(f"Unknown plot type: {plot_type}")

    # Get the appropriate visualization function
    viz_func = PLOT_HANDLERS[plot_type]

    # Call the visualization function to get the figure
    fig = await viz_func(adata, params, context)

    return fig


# =============================================================================
# Public API Functions
# =============================================================================


async def save_visualization(
    data_id: str,
    ctx: "ToolContext",
    plot_type: str,
    subtype: Optional[str] = None,
    output_dir: str = "./outputs",
    filename: Optional[str] = None,
    format: str = "png",
    dpi: Optional[int] = None,
) -> str:
    """Save a visualization to disk at publication quality by regenerating from metadata.

    This function regenerates visualizations from stored metadata (JSON) and the original
    data, then exports at the requested quality. This approach is more secure than
    loading serialized figure objects (pickle) because:
    1. JSON metadata cannot contain executable code
    2. Regeneration uses the trusted visualization codebase
    3. All parameters are human-readable and auditable

    Supports multiple formats including vector (PDF, SVG, EPS) and raster (PNG, JPEG, TIFF)
    with publication-ready metadata.

    Args:
        data_id: Dataset ID
        ctx: ToolContext for unified data access and logging
        plot_type: Type of plot to save (e.g., 'spatial', 'deconvolution', 'spatial_statistics')
        subtype: Optional subtype for plot types with variants
                 - For deconvolution: 'spatial_multi', 'dominant_type', 'diversity', etc.
                 - For spatial_statistics: 'neighborhood', 'co_occurrence', 'ripley', etc.
        output_dir: Directory to save the file (default: ./outputs)
        filename: Custom filename (optional, auto-generated if not provided)
        format: Image format (png, jpg, jpeg, pdf, svg, eps, ps, tiff)
        dpi: DPI for raster formats (default: 300 for publication quality)
              Vector formats (PDF, SVG, EPS, PS) ignore DPI

    Returns:
        Path to the saved file

    Raises:
        DataNotFoundError: If visualization metadata not found
        ProcessingError: If regeneration or saving fails
    """
    try:
        # Use environment variable for output_dir if default value was passed
        if output_dir == "./outputs":
            output_dir = get_output_dir_from_config(default="./outputs")

        # Validate format
        valid_formats = ["png", "jpg", "jpeg", "pdf", "svg", "eps", "ps", "tiff"]
        if format.lower() not in valid_formats:
            raise ParameterError(
                f"Invalid format: {format}. Must be one of {valid_formats}"
            )

        # Generate cache key with subtype if provided
        cache_key = (
            f"{data_id}_{plot_type}_{subtype}" if subtype else f"{data_id}_{plot_type}"
        )

        # Check if visualization exists in registry
        viz_entry = ctx.get_visualization(cache_key)

        # Set default DPI based on format
        if dpi is None:
            dpi = 300  # High quality for all formats (publication-ready)

        # Create output directory using safe path handling
        try:
            output_path = get_safe_output_path(
                output_dir, fallback_to_tmp=True, create_if_missing=True
            )
        except PermissionError as e:
            raise ProcessingError(
                f"Cannot save to {output_dir}: {e}. Check permissions."
            ) from e

        # Generate filename if not provided
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            plot_name = f"{plot_type}_{subtype}" if subtype else plot_type

            if dpi != 100:
                filename = f"{data_id}_{plot_name}_{dpi}dpi_{timestamp}.{format}"
            else:
                filename = f"{data_id}_{plot_name}_{timestamp}.{format}"
        else:
            # Ensure filename has correct extension
            if not filename.endswith(f".{format}"):
                filename = f"{filename}.{format}"

        # Full path for the file
        file_path = output_path / filename

        # Get visualization params from registry and regenerate figure
        if viz_entry is None:
            raise DataNotFoundError(
                f"Visualization '{plot_type}' not found. Use visualize_data() first."
            )

        # Regenerate figure from stored params (first principles: params + data = output)
        try:
            # Get params from registry entry and override DPI
            viz_params_dict = viz_entry.params.model_dump()
            viz_params_dict["dpi"] = dpi
            viz_params = VisualizationParameters(**viz_params_dict)

            # Regenerate the figure
            adata = await ctx.get_adata(data_id)
            cached_fig = await _regenerate_figure_for_export(
                adata, viz_params, ctx._mcp_context
            )

        except Exception as e:
            raise ProcessingError(
                f"Failed to regenerate '{cache_key}': {str(e)}"
            ) from e

        try:
            # Prepare save parameters
            save_params: Dict[str, Any] = {
                "bbox_inches": "tight",
                "facecolor": "white",
                "edgecolor": "none",
                "transparent": False,
                "pad_inches": 0.1,
            }

            # Format-specific settings
            if format.lower() == "pdf":
                save_params["dpi"] = dpi
                save_params["format"] = "pdf"
                save_params["metadata"] = {
                    "Title": f"{plot_type} visualization of {data_id}",
                    "Author": "ChatSpatial MCP",
                    "Subject": "Spatial Transcriptomics Analysis",
                    "Keywords": f"{plot_type}, {data_id}, spatial transcriptomics",
                    "Creator": "ChatSpatial with matplotlib",
                    "Producer": f"matplotlib {matplotlib.__version__}",
                }
            elif format.lower() == "svg":
                save_params["format"] = "svg"
            elif format.lower() in ["eps", "ps"]:
                save_params["format"] = format.lower()
            elif format.lower() in ["png", "jpg", "jpeg", "tiff"]:
                save_params["dpi"] = dpi
                save_params["format"] = format.lower()
                if format.lower() in ["jpg", "jpeg"]:
                    save_params["pil_kwargs"] = {"quality": 95}

            # Save the figure
            cached_fig.savefig(str(file_path), **save_params)

            return str(file_path)

        except Exception as e:
            raise ProcessingError(f"Failed to export visualization: {str(e)}") from e

    except (DataNotFoundError, ParameterError):
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to save visualization: {str(e)}") from e


async def export_all_visualizations(
    data_id: str,
    ctx: "ToolContext",
    output_dir: str = "./exports",
    format: str = "png",
    dpi: Optional[int] = None,
) -> List[str]:
    """Export all cached visualizations for a dataset to disk.

    Args:
        data_id: Dataset ID to export visualizations for
        ctx: ToolContext for unified data access and logging
        output_dir: Directory to save files
        format: Image format (png, jpg, pdf, svg)
        dpi: DPI for saved images (default: 300 for publication quality)

    Returns:
        List of paths to saved files
    """
    try:
        # Get visualization keys from registry (single source of truth)
        relevant_keys = ctx.list_visualizations(data_id)

        if not relevant_keys:
            await ctx.warning(f"No visualizations found for dataset '{data_id}'")
            return []

        saved_files = []

        for cache_key in relevant_keys:
            # Extract plot_type and subtype from cache key
            remainder = cache_key.replace(f"{data_id}_", "")

            # Known plot types that support subtypes
            known_plot_types_with_subtype = ["deconvolution", "spatial_statistics"]

            plot_type = None
            subtype = None

            # Try to match known plot types with subtypes
            for known_type in known_plot_types_with_subtype:
                if remainder.startswith(f"{known_type}_"):
                    plot_type = known_type
                    subtype = remainder[len(known_type) + 1 :]
                    break

            # If no match, treat the entire remainder as plot_type
            if plot_type is None:
                plot_type = remainder
                subtype = None

            try:
                saved_path = await save_visualization(
                    data_id=data_id,
                    ctx=ctx,
                    plot_type=plot_type,
                    subtype=subtype,
                    output_dir=output_dir,
                    format=format,
                    dpi=dpi,
                )
                saved_files.append(saved_path)
            except Exception as e:
                await ctx.warning(f"Failed to export {cache_key}: {str(e)}")

        return saved_files

    except ProcessingError:
        raise
    except Exception as e:
        raise ProcessingError(f"Failed to export visualizations: {str(e)}") from e


async def clear_visualization_cache(
    ctx: "ToolContext",
    data_id: Optional[str] = None,
) -> int:
    """Clear visualization cache to free memory.

    Args:
        ctx: ToolContext for unified data access and logging
        data_id: Optional dataset ID to clear specific visualizations

    Returns:
        Number of visualizations cleared
    """
    try:
        if data_id:
            # Clear specific dataset visualizations using prefix
            cleared_count = ctx.clear_visualizations(prefix=f"{data_id}_")
        else:
            # Clear all visualizations
            cleared_count = ctx.clear_visualizations()

        return cleared_count

    except Exception as e:
        raise ProcessingError(f"Failed to clear cache: {str(e)}") from e
