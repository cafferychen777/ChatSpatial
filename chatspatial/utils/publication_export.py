"""
High-quality scientific publication export utilities.

This module provides functions for exporting visualizations at publication quality,
including 300 DPI PDFs and vector formats (SVG, EPS) for scientific journals.
"""

import os
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any, List
# Set non-interactive backend for matplotlib to prevent GUI popups on macOS
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend - prevents Dock popup
import matplotlib.pyplot as plt
plt.ioff()  # Turn off interactive mode

from .image_utils import get_cached_figure, load_figure_pickle


async def export_for_publication(
    data_id: str,
    plot_type: str,
    output_path: str,
    format: str = "pdf",
    dpi: int = 300,
    params: Optional[Dict[str, Any]] = None,
    context: Optional[Any] = None
) -> str:
    """Export high-quality image for scientific publication
    
    This function exports cached visualizations at publication quality,
    suitable for submission to scientific journals like Nature, Science, Cell.
    
    Args:
        data_id: Dataset identifier
        plot_type: Type of plot (e.g., 'spatial', 'heatmap', 'umap')
        output_path: Full path for output file
        format: Output format ('pdf', 'svg', 'png', 'eps', 'ps')
        dpi: DPI for raster formats (default: 300 for publication)
        params: Additional export parameters
        context: MCP context for logging
        
    Returns:
        Path to the saved file
        
    Raises:
        FileNotFoundError: If figure not found in cache
        ValueError: If invalid format specified
        
    Examples:
        # For Nature/Science submission (vector PDF)
        await export_for_publication("mouse_brain", "spatial", "fig1.pdf", format="pdf")
        
        # For web publication (scalable SVG)
        await export_for_publication("mouse_brain", "heatmap", "fig2.svg", format="svg")
        
        # For high-res raster (PowerPoint, posters)
        await export_for_publication("mouse_brain", "umap", "fig3.png", dpi=600)
    """
    
    # Try to get Figure from memory cache first
    cache_key = f"{data_id}_{plot_type}"
    fig = get_cached_figure(cache_key)
    
    # If not in memory, try loading from pickle file
    if fig is None:
        pickle_path = f"/tmp/chatspatial/figures/{cache_key}.pkl"
        if os.path.exists(pickle_path):
            if context:
                await context.info(f"Loading figure from pickle: {pickle_path}")
            try:
                fig = load_figure_pickle(pickle_path)
            except Exception as e:
                if context:
                    await context.error(f"Failed to load pickle: {e}")
                raise FileNotFoundError(
                    f"Could not load figure from pickle. Error: {e}"
                )
        else:
            raise FileNotFoundError(
                f"Figure not found in cache for {cache_key}. "
                f"Please regenerate the visualization first."
            )
    
    # Validate format
    valid_formats = ['pdf', 'svg', 'png', 'eps', 'ps', 'jpg', 'jpeg', 'tiff']
    if format.lower() not in valid_formats:
        raise ValueError(
            f"Invalid format: {format}. Must be one of {valid_formats}"
        )
    
    # Prepare save parameters
    save_params = {
        'bbox_inches': 'tight',
        'facecolor': 'white',
        'edgecolor': 'none',
        'transparent': False,
        'pad_inches': 0.1
    }
    
    # Format-specific settings
    if format.lower() == 'pdf':
        # PDF specific settings for publication
        save_params['dpi'] = dpi
        save_params['format'] = 'pdf'  # CRITICAL: Must set format, not backend!
        save_params['metadata'] = {
            'Title': f'{plot_type} visualization of {data_id}',
            'Author': 'ChatSpatial MCP',
            'Subject': 'Spatial Transcriptomics Analysis',
            'Keywords': f'{plot_type}, {data_id}, spatial transcriptomics',
            'Creator': 'ChatSpatial with matplotlib',
            'Producer': f'matplotlib {matplotlib.__version__}'
        }
    elif format.lower() == 'svg':
        # SVG is vector, doesn't need DPI
        save_params.pop('dpi', None)
        save_params['format'] = 'svg'
    elif format.lower() in ['eps', 'ps']:
        # PostScript formats
        save_params['format'] = format.lower()
        save_params.pop('dpi', None)  # Vector format
    elif format.lower() in ['png', 'jpg', 'jpeg', 'tiff']:
        # Raster formats need DPI
        save_params['dpi'] = dpi
        save_params['format'] = format.lower()
        if format.lower() in ['jpg', 'jpeg']:
            save_params['quality'] = 95  # High quality for publication
    
    # Apply user-provided parameters
    if params:
        save_params.update(params)
    
    # Ensure output directory exists
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Log export details
    if context:
        format_info = f"{format.upper()}"
        if format.lower() in ['pdf', 'svg', 'eps', 'ps']:
            format_info += " (vector)"
        else:
            format_info += f" ({dpi} DPI)"
        
        await context.info(
            f"Exporting {plot_type} as {format_info} for publication..."
        )
    
    try:
        # Save the figure
        fig.savefig(output_path, **save_params)
        
        # Get file size
        file_size = Path(output_path).stat().st_size
        
        # Log success
        if context:
            size_str = f"{file_size/1024:.1f} KB" if file_size < 1024*1024 else f"{file_size/(1024*1024):.1f} MB"
            
            await context.info(
                f"✅ Exported publication-quality {format.upper()}: {output_path} ({size_str})"
            )
            
            # Provide format-specific advice
            if format.lower() == 'pdf':
                await context.info(
                    "📊 PDF ready for journal submission (Nature/Science/Cell compatible)"
                )
            elif format.lower() == 'svg':
                await context.info(
                    "📊 SVG ready for web publication or further editing in Illustrator/Inkscape"
                )
            elif format.lower() in ['eps', 'ps']:
                await context.info(
                    "📊 PostScript format ready for LaTeX inclusion or professional printing"
                )
            elif dpi >= 300:
                await context.info(
                    f"📊 High-resolution ({dpi} DPI) raster image suitable for publication"
                )
        
        return str(output_path)
        
    except Exception as e:
        error_msg = f"Failed to export {format.upper()}: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e


async def batch_export_for_publication(
    data_id: str,
    output_dir: str = "./publication_figures",
    formats: List[str] = None,
    dpi: int = 300,
    context: Optional[Any] = None
) -> Dict[str, List[str]]:
    """Batch export all cached visualizations for publication
    
    Export all visualizations for a dataset in multiple formats suitable
    for scientific publication. This is useful when preparing figures for
    a manuscript submission.
    
    Args:
        data_id: Dataset identifier
        output_dir: Directory for output files
        formats: List of formats (default: ['pdf', 'png'])
        dpi: DPI for raster formats (default: 300)
        context: MCP context for logging
        
    Returns:
        Dictionary mapping format to list of exported file paths
        
    Examples:
        # Export all figures for a paper
        results = await batch_export_for_publication(
            "mouse_brain",
            output_dir="./manuscript/figures",
            formats=["pdf", "svg", "png"],
            dpi=300
        )
    """
    import glob
    
    if formats is None:
        formats = ["pdf", "png"]
    
    # Find all pickle files for this dataset
    pickle_pattern = f"/tmp/chatspatial/figures/{data_id}_*.pkl"
    pickle_files = glob.glob(pickle_pattern)
    
    if not pickle_files:
        raise FileNotFoundError(
            f"No cached figures found for dataset: {data_id}. "
            f"Please generate visualizations first."
        )
    
    # Log batch export start
    if context:
        await context.info(
            f"Starting batch export for {len(pickle_files)} figures "
            f"in {len(formats)} format(s)"
        )
    
    # Initialize results dictionary
    results = {fmt: [] for fmt in formats}
    export_errors = []
    
    # Process each cached figure
    for pickle_path in pickle_files:
        # Extract plot_type from filename
        filename = Path(pickle_path).stem
        plot_type = filename.replace(f"{data_id}_", "")
        
        # Export in each requested format
        for fmt in formats:
            # Generate output filename with timestamp
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_filename = f"{data_id}_{plot_type}_{timestamp}.{fmt}"
            output_path = os.path.join(output_dir, output_filename)
            
            try:
                # Export the figure
                exported_path = await export_for_publication(
                    data_id=data_id,
                    plot_type=plot_type,
                    output_path=output_path,
                    format=fmt,
                    dpi=dpi,
                    context=context
                )
                results[fmt].append(exported_path)
                
            except Exception as e:
                error_msg = f"Failed to export {plot_type} as {fmt}: {e}"
                export_errors.append(error_msg)
                if context:
                    await context.warning(error_msg)
    
    # Log summary
    if context:
        total_exported = sum(len(files) for files in results.values())
        await context.info(
            f"✅ Batch export complete: {total_exported} files exported"
        )
        
        if export_errors:
            await context.warning(
                f"⚠️ {len(export_errors)} export(s) failed. Check warnings above."
            )
    
    return results


