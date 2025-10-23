"""
Main server implementation for ChatSpatial using the Spatial MCP Adapter.
"""

import logging
import os
import sys
import warnings
from typing import Any, Dict, List, Optional, Union

# Suppress warnings to speed up startup
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# CRITICAL: Disable progress bars to prevent stdout pollution
# This protects against accidental stdout usage if server is imported directly
os.environ["TQDM_DISABLE"] = "1"

# Suppress scanpy/squidpy verbosity
try:
    import scanpy as sc

    sc.settings.verbosity = 0
except ImportError:
    pass

from mcp.server.fastmcp import Context  # noqa: E402
from mcp.types import ImageContent  # noqa: E402

from .models.analysis import (
    AnnotationResult,
    CellCommunicationResult,  # noqa: E402
    CNVResult,
    DeconvolutionResult,
    DifferentialExpressionResult,
    EnrichmentResult,
    IntegrationResult,
    PreprocessingResult,
    RNAVelocityResult,
    SpatialStatisticsResult,
    SpatialDomainResult,
    SpatialVariableGenesResult,
    TrajectoryResult,
)
from .models.data import (
    CellCommunicationParameters,
    CNVParameters,  # noqa: E402
    ColumnInfo,
    DeconvolutionParameters,
    EnrichmentParameters,
    IntegrationParameters,
    RNAVelocityParameters,
    SpatialDataset,
    SpatialDomainParameters,
    SpatialVariableGenesParameters,
    TrajectoryParameters,
)
from .spatial_mcp_adapter import (
    MCPToolMetadata,
    create_spatial_mcp_server,
)  # noqa: E402
from .utils.error_handling import ProcessingError  # noqa: E402
from .utils.mcp_parameter_handler import (
    manual_parameter_validation,  # noqa: E402
    validate_analysis_params,
    validate_annotation_params,
    validate_cell_communication_params,
    validate_spatial_analysis_params,
    validate_visualization_params,
)
from .utils.tool_error_handling import mcp_tool_error_handler  # noqa: E402

# Lazy imports for heavy dependencies - imported when first used
# This significantly speeds up server startup time
# from .tools.annotation import annotate_cell_types  # Lazy loaded
# from .tools.cnv_analysis import infer_cnv  # Lazy loaded
# from .tools.deconvolution import deconvolve_spatial_data  # Lazy loaded
# from .tools.differential import differential_expression  # Lazy loaded
# from .tools.spatial_genes import identify_spatial_genes  # Lazy loaded
# from .tools.spatial_statistics import analyze_spatial_statistics  # Lazy loaded (squidpy is slow)
# from .tools.trajectory import analyze_rna_velocity  # Lazy loaded


logger = logging.getLogger(__name__)

# Create MCP server and adapter
mcp, adapter = create_spatial_mcp_server("ChatSpatial")

# Get data manager from adapter
data_manager = adapter.data_manager


def validate_dataset(data_id: str) -> None:
    """Validate that a dataset exists in the data store

    Args:
        data_id: Dataset ID

    Raises:
        ValueError: If the dataset is not found
    """
    if data_id not in data_manager.data_store:
        raise ValueError(f"Dataset {data_id} not found")


@mcp.tool()
@mcp_tool_error_handler()
async def load_data(
    data_path: str,
    data_type: str = "auto",
    name: Optional[str] = None,
    context: Context = None,
) -> SpatialDataset:
    """Load spatial transcriptomics data with comprehensive metadata profile

    Returns detailed information about the dataset structure to help with analysis:
    - Cell and gene counts
    - Available metadata columns with types and sample values
    - Multi-dimensional data (spatial coordinates, dimensionality reduction, etc.)
    - Gene expression profiles

    Args:
        data_path: Path to the data file or directory
        data_type: Type of spatial data (auto, 10x_visium, slide_seq, merfish, seqfish, other, h5ad).
                  If 'auto', will try to determine the type from the file extension or directory structure.
        name: Optional name for the dataset

    Returns:
        Comprehensive dataset information including metadata profiles
    """
    if context:
        await context.info(f"Loading data from {data_path} (type: {data_type})")

    # Load data using data manager
    data_id = await data_manager.load_dataset(data_path, data_type, name)
    dataset_info = await data_manager.get_dataset(data_id)

    if context:
        await context.info(
            f"Successfully loaded {dataset_info['type']} data with {dataset_info['n_cells']} cells and {dataset_info['n_genes']} genes"
        )

    # Create resource for the dataset
    await adapter.resource_manager.create_dataset_resource(data_id, dataset_info)

    # Convert column info from dict to ColumnInfo objects
    obs_columns = (
        [ColumnInfo(**col) for col in dataset_info.get("obs_columns", [])]
        if dataset_info.get("obs_columns")
        else None
    )
    var_columns = (
        [ColumnInfo(**col) for col in dataset_info.get("var_columns", [])]
        if dataset_info.get("var_columns")
        else None
    )

    # Return comprehensive dataset information
    return SpatialDataset(
        id=data_id,
        name=dataset_info["name"],
        data_type=dataset_info["type"],  # Use normalized type from dataset_info
        description=f"Spatial data: {dataset_info['n_cells']} cells × {dataset_info['n_genes']} genes",
        n_cells=dataset_info["n_cells"],
        n_genes=dataset_info["n_genes"],
        spatial_coordinates_available=dataset_info["spatial_coordinates_available"],
        tissue_image_available=dataset_info["tissue_image_available"],
        obs_columns=obs_columns,
        var_columns=var_columns,
        obsm_keys=dataset_info.get("obsm_keys"),
        uns_keys=dataset_info.get("uns_keys"),
        top_highly_variable_genes=dataset_info.get("top_highly_variable_genes"),
        top_expressed_genes=dataset_info.get("top_expressed_genes"),
    )


@mcp.tool()
@mcp_tool_error_handler()
@manual_parameter_validation(("params", validate_analysis_params))
async def preprocess_data(
    data_id: str, params: Any = None, context: Context = None
) -> PreprocessingResult:
    """Preprocess spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Preprocessing parameters

    Returns:
        Preprocessing result

    Notes:
        Available normalization methods:
        - log: Standard log normalization (default)
        - sct: SCTransform normalization
        - pearson_residuals: Modern Pearson residuals normalization (recommended for UMI data)
        - none: No normalization
        - scvi: Use scVI for normalization and dimensionality reduction

        When use_scvi_preprocessing=True, scVI will be used for advanced preprocessing
        including denoising and batch effect correction.

        Advanced configuration options:
        - n_neighbors: Number of neighbors for graph construction (default: 15)
        - clustering_resolution: Leiden clustering resolution (default: 1.0)
        - clustering_key: Key name for storing clustering results (default: "leiden")
        - spatial_key: Key name for spatial coordinates in obsm (default: None, auto-detected)
        - batch_key: Key name for batch information in obs (default: "batch")

        IMPORTANT: This preprocessing creates a filtered gene set for analysis efficiency.
        Raw data is automatically preserved in adata.raw for downstream analyses requiring
        comprehensive gene coverage (e.g., cell communication analysis with LIANA+).

        For cell communication analysis, you can later use data_source="raw" parameter
        to access the full unfiltered gene set.
    """
    # Import to avoid name conflict
    from .tools.preprocessing import preprocess_data as preprocess_func

    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Call preprocessing function
    result = await preprocess_func(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save preprocessing result
    await data_manager.save_result(data_id, "preprocessing", result)

    return result


@mcp.tool()
@mcp_tool_error_handler()  # CRITICAL: This decorator has special Image handling - see /docs/CRITICAL_IMAGE_DISPLAY_BUG.md
@manual_parameter_validation(("params", validate_visualization_params))
async def visualize_data(
    data_id: str, params: Any = None, context: Context = None
) -> Union[
    ImageContent, str
]:  # Simplified: ImageContent or str (MCP 2025 best practice)
    """Visualize spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Visualization parameters including:
            - plot_type: Type of visualization. Available types:
                        * Basic plots: spatial, heatmap, violin, umap
                        * Analysis results: spatial_domains, cell_communication, deconvolution,
                          trajectory, rna_velocity, spatial_statistics
                        * Multi-gene/correlation: multi_gene, lr_pairs, gene_correlation
                        * Enrichment: pathway_enrichment, spatial_enrichment
                        * Integration/QC: spatial_interaction, batch_integration
                        * CNV analysis: cnv_heatmap, spatial_cnv
                        * High-resolution: card_imputation
            - feature: Gene or feature to visualize (single gene as string or multiple genes as list)
                      For cell types, use the actual column name created by annotation methods:
                      - After Tangram: 'cell_type_tangram'
                      - After scANVI: 'cell_type_scanvi'
                      - After CellAssign: 'cell_type_cellassign'
                      - Or use clustering results: 'leiden', 'louvain', etc.
                      For lr_pairs plot_type: Can pass L-R pairs as ["Ligand^Receptor"] format
            - lr_pairs: (Optional) For lr_pairs plot_type, explicit list of (ligand, receptor) tuples
                       Example: [("Fn1", "Cd79a"), ("Vegfa", "Nrp2")]
            - colormap: Color scheme for visualization
            - figure_size: Size of the output figure

    Returns:
        Visualization image

    Examples:
        # Visualize cell types after annotation
        params = {
            "plot_type": "spatial",
            "feature": "cell_type_tangram",  # Use method-specific column name (e.g., cell_type_tangram, cell_type_scanvi)
            "colormap": "tab20"
        }

        # Visualize gene expression
        params = {
            "plot_type": "spatial",
            "feature": "Cd7",  # Gene name
            "colormap": "viridis"
        }

        # Visualize L-R pairs (Method 1: Using lr_pairs parameter)
        params = {
            "plot_type": "lr_pairs",
            "lr_pairs": [("Fn1", "Cd79a"), ("Vegfa", "Nrp2")]
        }

        # Visualize L-R pairs (Method 2: Using feature with ^ format)
        params = {
            "plot_type": "lr_pairs",
            "feature": ["Fn1^Cd79a", "Vegfa^Nrp2"]  # Will be parsed automatically
        }

        # Visualize spatial statistics results (Ripley, neighborhood enrichment, etc.)
        params = {
            "plot_type": "spatial_statistics",
            "subtype": "ripley",  # neighborhood, co_occurrence, ripley, moran, centrality, getis_ord
            "cluster_key": "leiden"
        }

        # Visualize deconvolution results
        params = {
            "plot_type": "deconvolution",
            "subtype": "dominant_type"  # spatial_multi, dominant_type, diversity, stacked_bar, scatterpie, umap
        }

        # NOTE: This tool does NOT provide demo/default data for scientific integrity.
        # All visualizations must be based on actual experimental results.
    """
    # Import to avoid name conflict
    from .tools.visualization import visualize_data as visualize_func

    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Parameter preprocessing and validation is now handled by the @manual_parameter_validation decorator
    # params is already a validated VisualizationParameters instance

    # Call visualization function
    image = await visualize_func(data_id, data_store, params, context)

    # Create visualization resource and return the image
    if image is not None:
        import time

        # Generate cache key with subtype if applicable
        # This handles plot types with subtypes (e.g., deconvolution, spatial_statistics)
        subtype = (
            params.subtype if hasattr(params, "subtype") and params.subtype else None
        )

        if subtype:
            cache_key = f"{data_id}_{params.plot_type}_{subtype}"
            viz_id = f"{data_id}_{params.plot_type}_{subtype}_{int(time.time())}"
        else:
            cache_key = f"{data_id}_{params.plot_type}"
            viz_id = f"{data_id}_{params.plot_type}_{int(time.time())}"

        metadata = {
            "data_id": data_id,
            "plot_type": params.plot_type,
            "subtype": subtype if subtype else "N/A",
            "feature": getattr(params, "feature", "N/A"),
            "timestamp": int(time.time()),
            "name": f"{params.plot_type} - {getattr(params, 'feature', 'N/A')}",
            "description": f"Visualization of {data_id}",
        }

        # Handle two return types: str (large images) or ImageContent (small images)
        if isinstance(image, str):
            # Large image: file path returned as text (MCP 2025 best practice)
            # Store a marker in cache indicating file path return
            adapter.resource_manager._visualization_cache[cache_key] = {
                "type": "file_path",
                "message": image,
                "timestamp": int(time.time()),
            }

            if context:
                await context.info(
                    "Large visualization saved to disk (following MCP best practice: URI over embedded content)"
                )
                await context.info(
                    f"Visualization type: {params.plot_type}, feature: {getattr(params, 'feature', 'N/A')}"
                )

            return image  # Return str directly (FastMCP converts to TextContent)

        # ImageContent (small images <70KB)
        image_data = image.data

        # Decode base64 string to bytes before caching
        import base64

        image_bytes = base64.b64decode(image_data)

        # Store in cache with simple key for easy retrieval
        adapter.resource_manager._visualization_cache[cache_key] = image_bytes

        # Also create resource with timestamped ID for uniqueness
        await adapter.resource_manager.create_visualization_resource(
            viz_id, image_bytes, metadata
        )

        file_size_kb = len(image_data) / 1024

        if context:
            await context.info(
                f"Image saved as resource: spatial://visualizations/{viz_id} ({file_size_kb:.1f}KB)"
            )
            await context.info(
                f"Visualization type: {params.plot_type}, feature: {getattr(params, 'feature', 'N/A')}"
            )

        # !!!!!!!!!! CRITICAL WARNING - DO NOT MODIFY !!!!!!!!!!
        # MUST return the raw Image object here!
        # DO NOT wrap in dictionary, DO NOT call to_image_content()!
        # The error handler will pass it through to FastMCP unchanged.
        # See /docs/CRITICAL_IMAGE_DISPLAY_BUG.md for why this is critical.
        # This bug took 2 WEEKS to find - DO NOT CHANGE!
        # !!!!!!!!!! CRITICAL WARNING - DO NOT MODIFY !!!!!!!!!!
        return image

    else:
        # Return error message if no image was generated
        return "Visualization generation failed, please check the data and parameter settings."


@mcp.tool()
@mcp_tool_error_handler()
async def save_visualization(
    data_id: str,
    plot_type: str,
    subtype: Optional[str] = None,
    output_dir: str = "./outputs",
    filename: Optional[str] = None,
    format: str = "png",
    dpi: Optional[int] = None,
    context: Context = None,
) -> str:
    """Save a visualization from cache to disk

    Args:
        data_id: Dataset ID
        plot_type: Type of plot to save (e.g., 'spatial', 'umap', 'deconvolution', 'spatial_statistics')
        subtype: Optional subtype for plot types with variants (e.g., 'neighborhood', 'scatterpie')
                 - For pathway_enrichment: 'enrichment_plot', 'barplot', 'dotplot', 'spatial'
                 - For deconvolution: 'spatial_multi', 'dominant_type', 'diversity', 'stacked_bar', 'scatterpie', 'umap'
                 - For spatial_statistics: 'neighborhood', 'co_occurrence', 'ripley', 'moran', 'centrality', 'getis_ord'
        output_dir: Directory to save the file (default: ./outputs)
        filename: Custom filename (optional, auto-generated if not provided)
        format: Image format (png, jpg, pdf, svg)
        dpi: DPI for saved image (default: 300 for publication quality)
             For publication quality, use 300+ DPI

    Returns:
        Path to the saved file

    Examples:
        Save a spatial plot: save_visualization("data1", "spatial")
        Save with subtype: save_visualization("data1", "spatial_statistics", subtype="neighborhood")
        Save deconvolution: save_visualization("data1", "deconvolution", subtype="scatterpie", format="pdf")
        Save for publication: save_visualization("data1", "spatial", dpi=300, format="png")
    """
    from .tools.visualization import save_visualization as save_func

    # Get the visualization cache from the adapter's resource manager
    visualization_cache = adapter.resource_manager._visualization_cache

    result = await save_func(
        data_id=data_id,
        plot_type=plot_type,
        subtype=subtype,
        output_dir=output_dir,
        filename=filename,
        format=format,
        dpi=dpi,
        visualization_cache=visualization_cache,
        context=context,
    )

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def list_saved_visualizations(
    output_dir: str = "./outputs",
    pattern: Optional[str] = None,
    context: Context = None,
) -> List[Dict[str, Any]]:
    """List all saved visualizations in the output directory

    Args:
        output_dir: Directory to search for saved files (default: ./outputs)
        pattern: Optional glob pattern to filter files (e.g., "*spatial*")

    Returns:
        List of saved visualization files with metadata

    Examples:
        List all saved files: list_saved_visualizations()
        List spatial plots: list_saved_visualizations(pattern="*spatial*")
        List from custom dir: list_saved_visualizations(output_dir="./my_outputs")
    """
    from .tools.visualization import list_saved_visualizations as list_func

    result = await list_func(
        output_dir=output_dir,
        pattern=pattern,
        context=context,
    )

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def export_all_visualizations(
    data_id: str,
    output_dir: str = "./exports",
    format: str = "png",
    dpi: Optional[int] = None,
    context: Context = None,
) -> List[str]:
    """Export all cached visualizations for a dataset to disk

    Args:
        data_id: Dataset ID to export visualizations for
        output_dir: Directory to save files (default: ./exports)
        format: Image format (png, jpg, jpeg, pdf, svg, eps, ps, tiff) (default: png)
        dpi: DPI for raster formats (default: 300 for publication quality)

    Returns:
        List of paths to saved files

    Examples:
        # Export all visualizations as PNG
        export_all_visualizations("data1")

        # Export all as PDF for publication
        export_all_visualizations("data1", format="pdf", dpi=300)

        # Export to custom directory as SVG
        export_all_visualizations("data1", "./my_exports", format="svg")
    """
    from .tools.visualization import export_all_visualizations as export_func

    # Get the visualization cache from the adapter's resource manager
    visualization_cache = adapter.resource_manager._visualization_cache

    result = await export_func(
        data_id=data_id,
        output_dir=output_dir,
        format=format,
        dpi=dpi,
        visualization_cache=visualization_cache,
        context=context,
    )

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def clear_visualization_cache(
    data_id: Optional[str] = None,
    context: Context = None,
) -> int:
    """Clear visualization cache to free memory

    Args:
        data_id: Optional dataset ID to clear specific visualizations (if None, clears all)

    Returns:
        Number of visualizations cleared

    Examples:
        Clear all visualizations: clear_visualization_cache()
        Clear for specific dataset: clear_visualization_cache("data1")
    """
    from .tools.visualization import clear_visualization_cache as clear_func

    # Get the visualization cache from the adapter's resource manager
    visualization_cache = adapter.resource_manager._visualization_cache

    result = await clear_func(
        data_id=data_id,
        visualization_cache=visualization_cache,
        context=context,
    )

    return result


@mcp.tool()
@mcp_tool_error_handler()
@manual_parameter_validation(("params", validate_annotation_params))
async def annotate_cell_types(
    data_id: str, params: Any = None, context: Context = None
) -> AnnotationResult:
    """Annotate cell types in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Annotation parameters

    Returns:
        Annotation result with cell type information and optional visualization

    Notes:
        Annotation methods (status):
        - tangram: Implemented (requires reference_data_id and PREPROCESSED reference data with HVGs)
        - scanvi: Implemented (deep learning label transfer via scvi-tools, requires reference_data_id)
        - cellassign: Implemented (via scvi-tools, requires marker_genes parameter)
        - mllmcelltype: Implemented (multimodal LLM classifier)
        - sctype: Implemented (requires R and rpy2)
        - singler: Implemented (Python-based via singler/celldex packages, requires singler_reference parameter)

        For methods requiring reference data (tangram, scanvi, singler):
        - tangram/scanvi: reference_data_id must point to a loaded AND PREPROCESSED single-cell dataset
        - IMPORTANT: Reference data MUST be preprocessed with preprocess_data() before use!
        - cell_type_key: Leave as None for auto-detection. Only set if you know the exact column name in reference data
        - Common cell type column names: 'cell_type', 'cell_types', 'celltype'
        - singler: Can use either reference_data_id OR singler_reference (celldex built-in references)

        Tangram-specific notes:
        - Method: Deep learning-based spatial mapping of single-cell to spatial transcriptomics
        - Requires: reference_data_id with PREPROCESSED single-cell data
        - Mapping modes (mode parameter):
          * mode="cells" (default): Maps individual cells to spatial locations
            - Preserves single-cell heterogeneity and fine-grained resolution
            - More computationally intensive (GPU recommended for large datasets)
            - Best for: Same specimen data, when cell-level detail is critical
          * mode="clusters" (recommended for cross-specimen): Aggregates cells by type before mapping
            - Dramatically improves performance, runs on standard laptop
            - Official recommendation: "Our choice when scRNAseq and spatial data come from different specimens"
            - Requires: cluster_label parameter (e.g., "cell_type")
            - Best for: Different specimens, limited resources, cell type distributions
            - Trades single-cell resolution for stability and speed
        - Confidence scores: Automatically normalized to [0, 1] probability range
        - GPU acceleration: Set tangram_device='cuda:0' if GPU available
        - Other parameters: tangram_density_prior, tangram_learning_rate, tangram_lambda_r

        scANVI-specific notes:
        - Method: Semi-supervised variational inference for label transfer
        - Requires: Both datasets must have 'counts' layer (raw counts)
        - Architecture: Configurable via scanvi_n_latent, scanvi_n_hidden, scanvi_dropout_rate
        - Small datasets (<1000 genes/cells): Use scanvi_n_latent=3-5, scanvi_dropout_rate=0.2,
          scanvi_use_scvi_pretrain=False, num_epochs=50 to prevent NaN errors
        - Returns probabilistic cell type predictions with confidence scores
        - GPU acceleration available (set tangram_device='cuda:0' if available)

        SingleR-specific notes:
        - Method: Reference-based correlation matching for cell type annotation
        - Reference options:
          * Built-in celldex references (via singler_reference parameter):
            - Human: 'hpca' (recommended), 'blueprint_encode', 'dice', 'monaco_immune', 'novershtern_hematopoietic'
            - Mouse: 'immgen' (recommended), 'mouse_rnaseq'
          * Custom reference (via reference_data_id parameter)
        - Common mistakes:
          * ✗ 'HumanPrimaryCellAtlasData' → ✓ use 'hpca'
          * ✗ 'ImmGenData' → ✓ use 'immgen'
        - Returns correlation-based confidence scores for cell type assignments
        - No GPU required (Python-based implementation via singler/celldex packages)
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Validate reference data for methods that require it
    if params.method in ["tangram", "scanvi", "singler"] and params.reference_data_id:
        if params.reference_data_id not in data_manager.data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")
        ref_info = await data_manager.get_dataset(params.reference_data_id)
        data_store[params.reference_data_id] = ref_info

    # Lazy import annotation tool (avoids slow startup)
    from .tools.annotation import annotate_cell_types

    # Call annotation function
    result = await annotate_cell_types(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save annotation result
    await data_manager.save_result(data_id, "annotation", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "annotation", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
@manual_parameter_validation(("params", validate_spatial_analysis_params))
async def analyze_spatial_statistics(
    data_id: str, params: Any = None, context: Context = None
) -> SpatialStatisticsResult:
    """Analyze spatial statistics and autocorrelation patterns

    Args:
        data_id: Dataset ID
        params: Analysis parameters

    Returns:
        Spatial statistics analysis result with statistics and optional visualization

    Notes:
        Available analysis types (implemented):
        - moran: Global Moran's I spatial autocorrelation (squidpy)
        - local_moran: Local Moran's I (LISA) for spatial clustering detection
        - geary: Geary's C spatial autocorrelation (squidpy)
        - getis_ord: Getis-Ord Gi* hot/cold spot detection (esda/PySAL)
          * Detects statistically significant spatial clusters of high/low values
          * Parameters: getis_ord_alpha (significance level), getis_ord_correction (FDR/Bonferroni)
          * Returns raw and corrected hotspot/coldspot counts
        - neighborhood: Neighborhood enrichment (squidpy)
        - co_occurrence: Co-occurrence analysis (squidpy)
        - centrality: Graph centrality scores (squidpy)
        - ripley: Ripley's K/L spatial point patterns
        - bivariate_moran: Bivariate Moran's I for gene pair correlation

        **Categorical Data Analysis (Choose based on number of categories):**
        - join_count: Traditional Join Count for BINARY data (exactly 2 categories)
          * Use for: Binary presence/absence, case/control, treated/untreated
          * Returns: Global statistics (BB/WW/BW joins, p-value)
          * Reference: Cliff & Ord (1981)

        - local_join_count: Local Join Count for MULTI-CATEGORY data (>2 categories)
          * Use for: Cell types, tissue domains, multi-class categorical variables
          * Returns: Per-category local clustering statistics with p-values
          * Identifies WHERE each category spatially clusters
          * Reference: Anselin & Li (2019)

        - network_properties: Spatial network analysis
        - spatial_centrality: Spatial-specific centrality measures
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Lazy import spatial_statistics (squidpy is slow to import)
    from .tools.spatial_statistics import (
        analyze_spatial_statistics as _analyze_spatial_statistics,
    )

    # Call spatial statistics analysis function
    result = await _analyze_spatial_statistics(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save spatial statistics result
    await data_manager.save_result(data_id, "spatial_statistics", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(
        data_id, "spatial_statistics", result
    )

    # Note: Visualization should be created separately using create_visualization tool
    # This maintains clean separation between analysis and visualization

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def find_markers(
    data_id: str,
    group_key: str,
    group1: Optional[str] = None,
    group2: Optional[str] = None,
    method: str = "wilcoxon",
    n_top_genes: int = 25,  # Number of top differentially expressed genes to return
    context: Context = None,
) -> DifferentialExpressionResult:
    """Find differentially expressed genes between groups

    Args:
        data_id: Dataset ID
        group_key: Column name defining groups
        group1: First group (if None, compare against all others)
        group2: Second group (if None, compare group1 against all others)
        method: Statistical test method
        n_top_genes: Number of top differentially expressed genes to return

    Returns:
        Differential expression result with top marker genes
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Lazy import differential expression tool
    from .tools.differential import differential_expression

    # Call differential expression function
    result = await differential_expression(
        data_id=data_id,
        data_store=data_store,
        group_key=group_key,
        group1=group1,
        group2=group2,
        method=method,
        n_top_genes=n_top_genes,  # Direct parameter pass-through - no conversion needed
        context=context,
    )

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save differential expression result
    await data_manager.save_result(data_id, "differential_expression", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "markers", result)

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def analyze_cnv(
    data_id: str,
    reference_key: str,
    reference_categories: List[str],
    method: str = "infercnvpy",
    window_size: int = 100,
    step: int = 10,
    exclude_chromosomes: Optional[List[str]] = None,
    dynamic_threshold: Optional[float] = 1.5,
    cluster_cells: bool = False,
    dendrogram: bool = False,
    numbat_genome: str = "hg38",
    numbat_allele_data_key: str = "allele_counts",
    numbat_t: float = 0.15,
    numbat_max_entropy: float = 0.8,
    numbat_min_cells: int = 10,
    numbat_ncores: int = 1,
    numbat_skip_nj: bool = False,
    context: Context = None,
) -> CNVResult:
    """Analyze copy number variations (CNVs) in spatial transcriptomics data

    Supports two CNV analysis methods:
    - infercnvpy: Expression-based CNV inference (default, fast)
    - Numbat: Haplotype-aware CNV analysis (requires allele data, more accurate)

    Args:
        data_id: Dataset identifier
        reference_key: Column name in adata.obs for cell type labels
        reference_categories: List of cell types to use as reference (normal cells)
        method: CNV analysis method ("infercnvpy" or "numbat", default: "infercnvpy")
        window_size: Number of genes for CNV averaging window (default: 100)
        step: Step size for sliding window (default: 10)
        exclude_chromosomes: Chromosomes to exclude (e.g., ['chrX', 'chrY'])
        dynamic_threshold: Threshold for dynamic CNV calling (default: 1.5)
        cluster_cells: Whether to cluster cells by CNV pattern
        dendrogram: Whether to compute hierarchical clustering dendrogram
        context: MCP context

    Returns:
        CNV analysis result with statistics and visualization availability

    Notes:
        CNV analysis methods:
        - infercnvpy: Expression-based (implemented, no allele data required)
        - numbat: Haplotype-aware (implemented when rpy2 installed, requires allele data)

        Numbat-specific notes:
        - Method: Haplotype-aware CNV analysis with phylogeny reconstruction
        - Requires: Allele-specific counts in adata.layers or adata.obsm
        - Allele data preparation: Use cellSNP-lite, pileup_and_phase, or similar tools
        - Genome options: hg38, hg19, mm10, mm39
        - Returns: CNV matrix, clone assignments, phylogeny tree
        - GPU acceleration: Not applicable (R-based method)

    Examples:
        # Basic infercnvpy analysis
        analyze_cnv("data1", "cell_type", ["T cells", "B cells"])

        # Numbat analysis (requires allele data)
        analyze_cnv("data1", "cell_type", ["T cells", "B cells"],
                   method="numbat", numbat_genome="hg38")

        # With clustering
        analyze_cnv("data1", "leiden", ["0", "1"], cluster_cells=True)
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Create CNVParameters object
    # Type: ignore needed for Literal parameters validated at runtime by Pydantic
    params = CNVParameters(
        method=method,  # type: ignore[arg-type]
        reference_key=reference_key,
        reference_categories=reference_categories,
        window_size=window_size,
        step=step,
        exclude_chromosomes=exclude_chromosomes,
        dynamic_threshold=dynamic_threshold,
        cluster_cells=cluster_cells,
        dendrogram=dendrogram,
        numbat_genome=numbat_genome,  # type: ignore[arg-type]
        numbat_allele_data_key=numbat_allele_data_key,
        numbat_t=numbat_t,
        numbat_max_entropy=numbat_max_entropy,
        numbat_min_cells=numbat_min_cells,
        numbat_ncores=numbat_ncores,
        numbat_skip_nj=numbat_skip_nj,
    )

    # Lazy import CNV analysis tool
    from .tools.cnv_analysis import infer_cnv

    # Call CNV inference function
    result = await infer_cnv(
        data_id=data_id, data_store=data_store, params=params, context=context
    )

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save CNV result
    await data_manager.save_result(data_id, "cnv_analysis", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "cnv", result)

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def analyze_velocity_data(
    data_id: str,
    params: RNAVelocityParameters = RNAVelocityParameters(),
    context: Context = None,
) -> RNAVelocityResult:
    """Analyze RNA velocity to understand cellular dynamics

    Args:
        data_id: Dataset ID
        params: RNA velocity parameters

    Returns:
        RNA velocity analysis result

    Notes:
        Velocity methods (status):
        - scvelo: scVelo with three modes (implemented, tested)
          - deterministic: Deterministic rate model
          - stochastic: Stochastic rate model (default)
          - dynamical: Dynamical model with ODE fitting
        - velovi: VeloVI deep learning method (implemented, requires scvi-tools, tested)
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Lazy import trajectory analysis tool
    from .tools.trajectory import analyze_rna_velocity

    # Call RNA velocity function
    result = await analyze_rna_velocity(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save velocity result
    await data_manager.save_result(data_id, "rna_velocity", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "velocity", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def analyze_trajectory_data(
    data_id: str,
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Context = None,
) -> TrajectoryResult:
    """Infer cellular trajectories and pseudotime

    Args:
        data_id: Dataset ID
        params: Trajectory analysis parameters

    Returns:
        Trajectory analysis result

    Notes:
        Trajectory methods (status):
        - dpt: Diffusion pseudotime (implemented)
        - palantir: Probabilistic trajectory inference (implemented when palantir installed)
        - cellrank: RNA velocity-based trajectory inference (implemented when cellrank installed)
        - velovi: scvi-tools VeloVI (implemented when scvi-tools available)
        - paga: Not implemented in tools/trajectory.py (planned)
    """
    # Import trajectory function
    from .tools.trajectory import analyze_trajectory

    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Call trajectory function
    result = await analyze_trajectory(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save trajectory result
    await data_manager.save_result(data_id, "trajectory", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "trajectory", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def integrate_samples(
    data_ids: List[str],
    params: IntegrationParameters = IntegrationParameters(),
    context: Context = None,
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples

    Args:
        data_ids: List of dataset IDs to integrate
        params: Integration parameters

    Returns:
        Integration result with integrated dataset ID

    Notes:
        Integration methods (status):
        - harmony, bbknn, scanorama: Classical methods (implemented)
        - scvi: Deep learning method (implemented, requires scvi-tools)

        Removed methods:
        - multivi: Requires MuData format (not compatible with current workflow)
        - contrastivevi: Not integrated (designed for Perturb-seq use cases)
    """
    # Import integration function
    from .tools.integration import integrate_samples

    # Validate all datasets
    for data_id in data_ids:
        validate_dataset(data_id)

    # Get all datasets from data manager
    data_store = {}
    for data_id in data_ids:
        dataset_info = await data_manager.get_dataset(data_id)
        data_store[data_id] = dataset_info

    # Call integration function
    result = await integrate_samples(data_ids, data_store, params, context)

    # Save integrated dataset
    integrated_id = result.data_id
    if integrated_id and integrated_id in data_store:
        data_manager.data_store[integrated_id] = data_store[integrated_id]

        # Create resource for integrated dataset
        await adapter.resource_manager.create_dataset_resource(
            integrated_id, data_store[integrated_id]
        )

    # Save integration result
    await data_manager.save_result(integrated_id, "integration", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(
        integrated_id, "integration", result
    )

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def deconvolve_data(
    data_id: str,
    params: DeconvolutionParameters,  # No default - LLM must provide parameters
    context: Context = None,
) -> DeconvolutionResult:
    """Deconvolve spatial spots to estimate cell type proportions

    Args:
        data_id: Dataset ID
        params: Deconvolution parameters including:
                - method: Deconvolution method to use
                - cell_type_key: Key in reference data for cell types (REQUIRED)
                - reference_data_id: Reference single-cell dataset ID (required for most methods)

                Cell2location-specific parameters (official scvi-tools recommendations):
                - ref_model_epochs: Reference model training epochs (default: 250)
                - n_epochs: Cell2location model training epochs (default: 30000)
                - n_cells_per_spot: Expected cells per location (default: 30, tissue-dependent)
                - detection_alpha: RNA detection sensitivity (default: 200)

    Returns:
        Deconvolution result with cell type proportions

    Notes:
        Deconvolution methods (status):
        - cell2location, destvi, stereoscope, tangram: Implemented when scvi-tools available
        - rctd, spotlight: Implemented via rpy2/R when R packages are installed
        - card: Implemented (CARD deconvolution method)

        Cell2location uses two-stage training:
        1. Reference model (NB regression): Learns cell type signatures (250 epochs)
        2. Cell2location model: Maps cell types to spatial locations (30000 epochs)
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Validate reference data if provided
    if params.reference_data_id:
        if params.reference_data_id not in data_manager.data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")
        ref_info = await data_manager.get_dataset(params.reference_data_id)
        data_store[params.reference_data_id] = ref_info

    # Lazy import deconvolution tool
    from .tools.deconvolution import deconvolve_spatial_data

    # Call deconvolution function
    result = await deconvolve_spatial_data(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save deconvolution result
    await data_manager.save_result(data_id, "deconvolution", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(
        data_id, "deconvolution", result
    )

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def identify_spatial_domains(
    data_id: str,
    params: SpatialDomainParameters = SpatialDomainParameters(),
    context: Context = None,
) -> SpatialDomainResult:
    """Identify spatial domains and tissue architecture

    Args:
        data_id: Dataset ID
        params: Spatial domain parameters

    Returns:
        Spatial domain result with identified domains

    Notes:
        Spatial domain methods (status):
        - spagcn: SpaGCN graph convolutional network (implemented; optional dependency SpaGCN)
        - leiden / louvain: clustering-based (implemented; no extra deps)
        - stagate: STAGATE (implemented; optional dependency STAGATE)
        - graphst: GraphST graph self-supervised contrastive learning (implemented; optional dependency GraphST)
        - stlearn / sedr / bayesspace: not implemented in this server; planned/experimental
    """
    # Import spatial domains function
    from .tools.spatial_domains import identify_spatial_domains as identify_domains_func

    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Call spatial domains function
    result = await identify_domains_func(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save spatial domains result
    await data_manager.save_result(data_id, "spatial_domains", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "domains", result)

    return result


@mcp.tool()
@mcp_tool_error_handler()
@manual_parameter_validation(("params", validate_cell_communication_params))
async def analyze_cell_communication(
    data_id: str,
    params: CellCommunicationParameters,  # No default - LLM must provide parameters
    context: Context = None,
) -> CellCommunicationResult:
    """Analyze cell-cell communication patterns

    Args:
        data_id: Dataset ID
        params: Cell communication parameters

    Returns:
        Cell communication analysis result

    Notes:
        Cell communication methods (status):
        - liana: Implemented (global/cluster and spatial bivariate modes; requires liana)
        - cellphonedb: Implemented (statistical analysis with spatial microenvironments; requires cellphonedb)
        - cellchat_liana: Implemented (CellChat algorithm via LIANA framework; requires liana)
        - nichenet / connectome / cytotalk / squidpy: Not implemented in this server

        IMPORTANT: For comprehensive cell communication analysis:

        **Species-specific configuration:**
        - species="mouse" + liana_resource="mouseconsensus" for mouse data
        - species="human" + liana_resource="consensus" for human data
        - species="zebrafish" for zebrafish data

        **Available LIANA resources (liana_resource parameter):**
        - "consensus" (default, recommended): Consensus of multiple databases
        - "mouseconsensus": Mouse-specific consensus database
        - "cellphonedb": CellPhoneDB database (curated, stringent)
        - "celltalkdb": CellTalkDB database (large, comprehensive)
        - "icellnet": iCellNet database (immune cell focus)
        - "cellchatdb": CellChat database
        - "connectomedb2020": Connectome database 2020
        - "baccin2019", "cellcall", "cellinker", "embrace", "guide2pharma",
          "hpmr", "italk", "kirouac2010", "lrdb", "ramilowski2015": Additional resources

        **Data source selection:**
        - data_source="raw" - Use raw unfiltered data (recommended for comprehensive gene coverage)
        - data_source="current" - Use current processed data (may have limited genes)

        **Common failure scenarios and solutions:**
        1. "Too few features from resource found in data":
           - Use data_source="raw" to access full gene set
           - Ensure species matches data (mouse vs human)
           - Use species-appropriate resource (mouseconsensus for mouse)

        2. Missing spatial connectivity:
           - Run spatial neighbor computation in preprocessing step (see below)

        3. Missing cell type annotations:
           - Ensure cell_type_key column exists or run annotation first

        **Spatial connectivity computation (preprocessing step):**

        The spatial neighborhood definition profoundly impacts cell communication analysis results.
        Choose parameters based on your spatial transcriptomics platform and biological question:

        **Platform-specific recommendations:**

        10x Visium (hexagonal grid, 55µm spots, 100µm center-to-center spacing):
          • coord_type: "grid" (for hexagonal layout) or "generic" (for custom)
          • n_neighs: 6 (direct neighbors in hexagonal grid)
          • n_rings: 1-2 (for grid mode: 1=first ring only, 2=first+second ring)
          • radius: 150-200 pixels (for distance-based, ~captures first neighbor ring)
          ├─ Local interactions (paracrine signaling): n_neighs=6 or n_rings=1
          ├─ Microenvironment analysis: n_neighs=12-18 or n_rings=2
          └─ Broader spatial context: radius=300-500 pixels

        Slide-seq/Slide-seqV2 (10µm beads, high density):
          • coord_type: "generic"
          • n_neighs: 10-30 (higher density requires more neighbors)
          • radius: 50-100 µm (typical cell-cell signaling range)
          ├─ Dense regions: n_neighs=20-30
          ├─ Sparse regions: n_neighs=10-15
          └─ Distance-based: radius=50-100 µm (matches biological signaling range)

        MERFISH/seqFISH+ (single-cell resolution, <1µm precision):
          • coord_type: "generic"
          • n_neighs: 3-10 (nearest cell neighbors)
          • radius: 20-50 µm (direct cell-cell contact to short-range paracrine)
          ├─ Direct contact: n_neighs=3-5 or radius=10-20 µm
          ├─ Paracrine signaling: n_neighs=5-10 or radius=30-50 µm
          └─ Microenvironment: radius=50-100 µm

        **Biological considerations:**

        Cell communication distance ranges (from literature):
          • Juxtacrine signaling: 0-10 µm (direct contact)
          • Paracrine signaling: 10-100 µm (e.g., Wnt/Wg: ~50-100 µm)
          • Broader microenvironment: 100-500 µm

        Analysis goal-based selection:
          • Identify direct cell-cell interactions → Use smaller neighborhoods (n_neighs=6-10, radius=50-100 µm)
          • Study tissue microenvironments → Use larger neighborhoods (n_neighs=15-30, radius=200-500 µm)
          • Rare cell type interactions → Use adaptive/larger k to avoid missing signals
          • Abundant cell types → Use smaller k to avoid spurious connections

        **Parameter tradeoffs:**
          • Larger neighborhoods: Capture long-range signals but lose spatial specificity
          • Smaller neighborhoods: High spatial precision but may miss important interactions
          • Fixed k (n_neighs): Same number for all spots, may overcluster dense regions
          • Distance-based (radius): More biologically meaningful but varying neighbor counts

        **Examples:**

        Visium - local paracrine signaling:
          # Step 1: Compute spatial neighbors (preprocessing)
          import squidpy as sq
          sq.gr.spatial_neighbors(adata, coord_type='grid', n_rings=1)

          # Step 2: Analyze communication
          params = {
              "species": "human",
              "liana_resource": "consensus",
              "data_source": "raw"
          }

        Visium - microenvironment analysis:
          # Step 1: Compute spatial neighbors (preprocessing)
          import squidpy as sq
          sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=18)

          # Step 2: Analyze communication
          params = {
              "species": "human",
              "data_source": "raw"
          }

        MERFISH - direct cell-cell contact:
          # Step 1: Compute spatial neighbors (preprocessing)
          import squidpy as sq
          sq.gr.spatial_neighbors(adata, coord_type='generic', radius=20)

          # Step 2: Analyze communication
          params = {
              "species": "mouse",
              "liana_resource": "mouseconsensus",
              "data_source": "raw"
          }

        **References:**
          • Squidpy framework: Palla et al., Nat Methods 2022
          • LIANA+: Dimitrov et al., Nat Cell Biol 2024
          • Visium resolution: 10x Genomics Technical Note
          • Signaling ranges: Literature-based (Wnt/Wg: ~50-100 µm)
    """
    # Import cell communication function
    from .tools.cell_communication import (
        analyze_cell_communication as analyze_comm_func,
    )

    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Call cell communication function
    result = await analyze_comm_func(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save communication result
    await data_manager.save_result(data_id, "cell_communication", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(
        data_id, "communication", result
    )

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def analyze_enrichment(
    data_id: str,
    params: Optional[EnrichmentParameters] = None,
    context: Context = None,
) -> EnrichmentResult:
    """Perform gene set enrichment analysis

    Args:
        data_id: Dataset ID
        params: Enrichment analysis parameters (REQUIRED: species must be specified)

    Returns:
        Enrichment analysis result

    IMPORTANT - Species and Database Selection:
    You MUST specify 'species' parameter explicitly. No default species is assumed.

    Recommended database combinations by species:

    FOR MOUSE DATA (species="mouse"):
    - "KEGG_Pathways" (recommended, uses KEGG_2019_Mouse internally)
    - "Reactome_Pathways" (comprehensive pathway database)
    - "MSigDB_Hallmark" (curated hallmark gene sets)
    - "GO_Biological_Process" (works but may have fewer matches)

    FOR HUMAN DATA (species="human"):
    - "KEGG_Pathways" (recommended, uses KEGG_2021_Human internally)
    - "Reactome_Pathways" (comprehensive pathway database)
    - "MSigDB_Hallmark" (curated hallmark gene sets)
    - "GO_Biological_Process" (standard GO terms)

    Available gene_set_database options:
    - "GO_Biological_Process" (default, auto-adapts to species)
    - "GO_Molecular_Function" (GO molecular function terms)
    - "GO_Cellular_Component" (GO cellular component terms)
    - "KEGG_Pathways" (species-specific: KEGG_2021_Human or KEGG_2019_Mouse)
    - "Reactome_Pathways" (Reactome_2022 pathway database)
    - "MSigDB_Hallmark" (MSigDB_Hallmark_2020 curated gene sets)
    - "Cell_Type_Markers" (cell type marker genes)
    - Custom gene sets via gene_sets parameter

    Methods available:
    - "pathway_ora": Over-representation analysis (recommended)
    - "pathway_enrichr": Enrichr web service
    - "pathway_gsea": Gene Set Enrichment Analysis
    - "pathway_ssgsea": Single-sample GSEA
    - "spatial_enrichmap": Spatial enrichment mapping

    Result Optimization:
    This tool automatically limits returned results to the 50 most significant pathways
    for efficient LLM communication. Complete results are preserved in adata.uns for
    downstream visualization and analysis. This optimization is transparent to users
    and maintains scientific validity while reducing token usage.

    Example usage:
    For mouse data:  params={"species": "mouse", "gene_set_database": "KEGG_Pathways"}
    For human data:  params={"species": "human", "gene_set_database": "KEGG_Pathways"}
    """
    # Import enrichment analysis function
    import time

    from .tools.enrichment import (
        perform_spatial_enrichment as perform_enrichment_analysis,
    )

    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Start timing
    start_time = time.time()

    # Check if params is None (parameter is required now)
    if params is None:
        raise ValueError(
            "params parameter is required for enrichment analysis.\n"
            "You must provide EnrichmentParameters with at least 'species' specified.\n"
            "Example: params={'species': 'mouse', 'method': 'pathway_ora'}"
        )

    # Get adata for gene set handling
    adata = data_store[data_id]["adata"]

    # Handle gene sets - either user-provided or from database
    gene_sets = params.gene_sets

    # If no gene sets provided, load from database
    if gene_sets is None and params.gene_set_database:
        if context:
            await context.info(f"Loading gene sets from {params.gene_set_database}")

        # Load gene sets based on database name
        try:
            from .tools.enrichment import load_gene_sets
        except ImportError:
            raise ProcessingError(
                "gseapy package is required for gene set loading. Install with: pip install gseapy"
            )
        try:
            # Use species from explicit parameter (no defaults, no auto-detection)
            if not hasattr(params, "species") or params.species is None:
                raise ValueError(
                    "species parameter is required for enrichment analysis.\n"
                    "Please specify:\n"
                    "  - 'human' for human data (genes like CD5L, PTPRC)\n"
                    "  - 'mouse' for mouse data (genes like Cd5l, Ptprc)\n"
                    "  - 'zebrafish' for zebrafish data\n"
                    "Example: params={'species': 'mouse', 'method': 'pathway_ora'}"
                )

            species = params.species  # Use explicit value from LLM

            gene_sets = await load_gene_sets(
                database=params.gene_set_database,
                species=species,
                min_genes=params.min_genes,
                max_genes=params.max_genes if hasattr(params, "max_genes") else 500,
                context=context,
            )

            if context:
                await context.info(
                    f"Loaded {len(gene_sets)} gene sets from {params.gene_set_database}"
                )

            # === CRITICAL TOKEN OPTIMIZATION ===
            # Limit gene sets to prevent token bloat for insights-based analysis
            if len(gene_sets) > 50:  # Limit to top 50 pathways for token efficiency
                if context:
                    await context.info(
                        f"Limiting to top 50 pathways from {len(gene_sets)} total for token optimization"
                    )

                # Sort pathways by gene count (larger pathways often more important)
                sorted_pathways = sorted(
                    gene_sets.items(), key=lambda x: len(x[1]), reverse=True
                )
                gene_sets = dict(sorted_pathways[:50])

                if context:
                    await context.info(
                        f"Using {len(gene_sets)} top pathways for analysis"
                    )

        except Exception as e:
            # NO FALLBACK: Enrichment analysis requires specific gene sets for scientific validity
            error_msg = (
                f"Failed to load gene sets from {params.gene_set_database}: {e}\n\n"
                f"ENRICHMENT ANALYSIS REQUIRES SPECIFIC GENE SETS\n\n"
                f"Gene set enrichment analysis cannot proceed with arbitrary gene substitutions.\n"
                f"This preserves scientific integrity and prevents misleading results.\n\n"
                f"SOLUTIONS:\n"
                f"1. Check your internet connection (required for database access)\n"
                f"2. Verify species parameter: '{params.species}' (use 'human' or 'mouse')\n"
                f"3. Try a different database:\n"
                f"   - 'KEGG_Pathways' (recommended for pathway analysis)\n"
                f"   - 'GO_Biological_Process' (for biological processes)\n"
                f"   - 'Reactome_Pathways' (for molecular pathways)\n"
                f"   - 'MSigDB_Hallmark' (for hallmark gene sets)\n"
                f"4. Provide custom gene sets via 'gene_sets' parameter\n"
                f"5. Use spatial analysis tools for data-driven insights without predefined pathways\n\n"
                f"WHY NO FALLBACK:\n"
                f"Using different gene sets (like highly variable genes) would produce\n"
                f"scientifically different results while appearing to be pathway analysis."
            )

            if context:
                await context.error(f"Gene set database loading failed: {e}")
                await context.error("No fallback - preserving scientific integrity")

            raise ProcessingError(error_msg)

    # Verify we have valid gene sets (should not be None after proper error handling above)
    if gene_sets is None or len(gene_sets) == 0:
        # This should not happen with proper error handling above, but safety check
        raise ProcessingError(
            "No valid gene sets available for enrichment analysis. "
            "Please provide gene sets via 'gene_sets' parameter or specify a valid 'gene_set_database'."
        )

    # Call appropriate enrichment function based on method
    if params.method == "spatial_enrichmap":
        # Spatial enrichment analysis using EnrichMap
        result_dict = (
            await perform_enrichment_analysis(
                data_id=data_id,
                data_store=data_store,
                gene_sets=gene_sets,
                score_keys=params.score_keys,
                spatial_key=params.spatial_key,
                n_neighbors=params.n_neighbors,
                smoothing=params.smoothing,
                correct_spatial_covariates=params.correct_spatial_covariates,
                batch_key=params.batch_key,
                gene_weights=params.gene_weights,
                species=params.species,
                database=params.gene_set_database,
                context=context,
            )
        ).to_dict()
        if context:
            await context.info(
                "Spatial enrichment analysis complete. Use create_visualization tool with plot_type='spatial_enrichment' to visualize results"
            )
    else:
        # Generic enrichment analysis (GSEA, ORA, ssGSEA, Enrichr)
        from .tools.enrichment import (
            perform_enrichr,
            perform_gsea,
            perform_ora,
            perform_ssgsea,
        )

        if params.method == "pathway_gsea":
            result_dict = (
                await perform_gsea(
                    adata=adata,
                    gene_sets=gene_sets,
                    ranking_key=params.score_keys,
                    permutation_num=params.n_permutations,
                    min_size=params.min_genes,
                    max_size=params.max_genes,
                    species=params.species,
                    database=params.gene_set_database,
                    context=context,
                )
            ).to_dict()
            if context:
                await context.info(
                    "Pathway GSEA analysis complete. Use create_visualization tool with plot_type='pathway_enrichment' to visualize results"
                )
        elif params.method == "pathway_ora":
            result_dict = (
                await perform_ora(
                    adata=adata,
                    gene_sets=gene_sets,
                    pvalue_threshold=params.pvalue_cutoff,
                    min_size=params.min_genes,
                    max_size=params.max_genes,
                    species=params.species,
                    database=params.gene_set_database,
                    context=context,
                )
            ).to_dict()
            if context:
                await context.info(
                    "Pathway ORA analysis complete. Use create_visualization tool with plot_type='pathway_enrichment' to visualize results"
                )
        elif params.method == "pathway_ssgsea":
            result_dict = (
                await perform_ssgsea(
                    adata=adata,
                    gene_sets=gene_sets,
                    min_size=params.min_genes,
                    max_size=params.max_genes,
                    species=params.species,
                    database=params.gene_set_database,
                    context=context,
                )
            ).to_dict()
            if context:
                await context.info(
                    "Pathway ssGSEA analysis complete. Use create_visualization tool with plot_type='pathway_enrichment' to visualize results"
                )
        elif params.method == "pathway_enrichr":
            # For Enrichr, we need a gene list - use HVG or top variable genes
            if "highly_variable" in adata.var:
                gene_list = adata.var_names[adata.var["highly_variable"]].tolist()[:500]
            else:
                # Use top variable genes
                import numpy as np
                from scipy import sparse

                if sparse.issparse(adata.X):
                    # Handle sparse matrix
                    var_scores = np.array(adata.X.toarray().var(axis=0)).flatten()
                else:
                    var_scores = np.array(adata.X.var(axis=0)).flatten()
                top_indices = np.argsort(var_scores)[-500:]
                gene_list = adata.var_names[top_indices].tolist()

            result_dict = (
                await perform_enrichr(
                    gene_list=gene_list,
                    gene_sets=params.gene_set_database,
                    organism=params.species,  # Use explicit species from params
                    context=context,
                )
            ).to_dict()
            if context:
                await context.info(
                    "Pathway Enrichr analysis complete. Use create_visualization tool with plot_type='pathway_enrichment' to visualize results"
                )
        else:
            raise ValueError(f"Unknown enrichment method: {params.method}")

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # === AUTOMATIC TOKEN OPTIMIZATION ===
    # To ensure efficient LLM communication, we automatically limit the number of
    # pathways returned to 50 most significant results. This is transparent to users
    # and prevents excessive token usage while preserving all scientifically important
    # findings. Complete results are stored in adata.uns for downstream analysis.
    #
    # Rationale:
    # - Full enrichment results are preserved in adata for visualization and further analysis
    # - Token efficiency: Returning 50 pathways instead of 500+ saves ~90% of tokens
    # - Scientific validity: Top 50 pathways capture the most biologically significant findings
    # - User experience: Focused results are easier to interpret than massive lists
    MAX_RESULTS_FOR_LLM = 50  # Automatic token optimization (transparent to user)

    # Get top significant gene sets
    adjusted_pvals = result_dict.get("adjusted_pvalues", {})
    pvals = result_dict.get("pvalues", {})

    # Get significant gene sets
    significant_sets = (
        {
            k
            for k, v in adjusted_pvals.items()
            if v is not None and v < params.pvalue_cutoff
        }
        if adjusted_pvals
        else set()
    )

    # Get top gene sets
    top_sets = set(result_dict.get("top_gene_sets", []))
    top_depleted = set(result_dict.get("top_depleted_sets", []))

    # Combine significant and top sets, limit to MAX_RESULTS_FOR_LLM
    display_sets = list((significant_sets | top_sets | top_depleted))[
        :MAX_RESULTS_FOR_LLM
    ]

    # Filter large dictionaries to only include display sets
    filtered_enrichment_scores = {
        k: v
        for k, v in result_dict.get("enrichment_scores", {}).items()
        if k in display_sets
    }
    filtered_pvalues = {k: v for k, v in pvals.items() if k in display_sets}
    filtered_adjusted_pvalues = {
        k: v for k, v in adjusted_pvals.items() if k in display_sets
    }
    filtered_gene_set_statistics = {
        k: v
        for k, v in result_dict.get("gene_set_statistics", {}).items()
        if k in display_sets
    }

    # Filter gene set summaries to display sets
    gene_set_summaries = result_dict.get("gene_set_summaries", {})
    filtered_gene_set_summaries = {
        k: v for k, v in gene_set_summaries.items() if k in display_sets
    }

    # Create EnrichmentResult object
    result = EnrichmentResult(
        method=params.method,
        n_gene_sets=result_dict.get("n_gene_sets", 0),
        n_significant=result_dict.get("n_significant", 0),
        enrichment_scores=filtered_enrichment_scores,
        pvalues=filtered_pvalues,
        adjusted_pvalues=filtered_adjusted_pvalues,
        gene_set_statistics=filtered_gene_set_statistics,
        gene_set_summaries=filtered_gene_set_summaries,
        spatial_metrics=result_dict.get("spatial_metrics"),
        spatial_scores_key=result_dict.get("spatial_scores_key"),
        top_gene_sets=result_dict.get("top_gene_sets", []),
        top_depleted_sets=result_dict.get("top_depleted_sets", []),
        parameters_used=params.model_dump(),
        computation_time=time.time() - start_time,
    )

    # Save enrichment result
    await data_manager.save_result(data_id, "enrichment", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "enrichment", result)

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def find_spatial_genes(
    data_id: str,
    params: SpatialVariableGenesParameters = SpatialVariableGenesParameters(),
    context: Context = None,
) -> SpatialVariableGenesResult:
    """Identify spatially variable genes using various methods

    Args:
        data_id: Dataset ID
        params: Spatial variable gene parameters

    Returns:
        Spatial variable genes result

    Notes:
        Available methods:
        - sparkx: SPARK-X non-parametric method (default, best accuracy)
        - spatialde: SpatialDE Gaussian process-based method (statistically rigorous)

        Method selection via params.method parameter.
        Each method has specific parameters - see SpatialVariableGenesParameters model.

        Performance comparison (3000 spots × 20000 genes):
        - SPARK-X: ~2-5 min (best accuracy)
        - SpatialDE: ~15-30 min (best statistical rigor)
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Lazy import spatial genes tool
    from .tools.spatial_genes import identify_spatial_genes

    # Call spatial genes function
    result = await identify_spatial_genes(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save spatial genes result
    await data_manager.save_result(data_id, "spatial_genes", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(
        data_id, "spatial_genes", result
    )

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def register_spatial_data(
    source_id: str,
    target_id: str,
    method: str = "paste",
    landmarks: Optional[List[Dict[str, Any]]] = None,
    context: Context = None,
) -> Dict[str, Any]:
    """Register/align spatial transcriptomics data across sections

    Args:
        source_id: Source dataset ID
        target_id: Target dataset ID to align to
        method: Registration method (paste, stalign)
        landmarks: Additional parameters for registration methods

    Returns:
        Registration result with transformation matrix
    """
    # Import registration function
    from .tools.spatial_registration import register_spatial_slices_mcp

    # Validate datasets
    validate_dataset(source_id)
    validate_dataset(target_id)

    # Get datasets from data manager
    source_info = await data_manager.get_dataset(source_id)
    target_info = await data_manager.get_dataset(target_id)
    data_store = {source_id: source_info, target_id: target_info}

    # Call registration function using standard architecture
    result = await register_spatial_slices_mcp(
        source_id, target_id, data_store, method, context
    )

    # Update datasets in data manager
    data_manager.data_store[source_id] = data_store[source_id]
    data_manager.data_store[target_id] = data_store[target_id]

    # Save registration result
    await data_manager.save_result(source_id, "registration", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(
        source_id, "registration", result
    )

    return result


# Register tool metadata with the adapter
tool_metadata = {
    "load_data": MCPToolMetadata(
        name="load_data",
        title="Load Spatial Data",
        description="Load spatial transcriptomics data from file",
        read_only_hint=True,
        idempotent_hint=True,
        open_world_hint=True,
    ),
    "preprocess_data": MCPToolMetadata(
        name="preprocess_data",
        title="Preprocess Spatial Data",
        description="Preprocess and normalize spatial data",
        read_only_hint=False,
        idempotent_hint=False,
    ),
    "visualize_data": MCPToolMetadata(
        name="visualize_data",
        title="Visualize Spatial Data",
        description="Create visualizations of spatial data",
        read_only_hint=True,
        idempotent_hint=True,
    ),
    "annotate_cell_types": MCPToolMetadata(
        name="annotate_cell_types",
        title="Annotate Cell Types",
        description="Identify cell types in spatial data",
        read_only_hint=False,
        idempotent_hint=False,
        open_world_hint=True,
    ),
    "analyze_spatial_statistics": MCPToolMetadata(
        name="analyze_spatial_statistics",
        title="Spatial Pattern Analysis",
        description="Analyze spatial patterns and correlations",
        read_only_hint=False,
        idempotent_hint=True,
    ),
    "find_markers": MCPToolMetadata(
        name="find_markers",
        title="Find Marker Genes",
        description="Identify differentially expressed genes",
        read_only_hint=True,
        idempotent_hint=True,
    ),
    "analyze_velocity_data": MCPToolMetadata(
        name="analyze_velocity_data",
        title="RNA Velocity Analysis",
        description="Analyze RNA velocity dynamics",
        read_only_hint=False,
        idempotent_hint=False,
    ),
    "analyze_trajectory_data": MCPToolMetadata(
        name="analyze_trajectory_data",
        title="Trajectory Analysis",
        description="Infer cellular trajectories",
        read_only_hint=False,
        idempotent_hint=False,
    ),
    "integrate_samples": MCPToolMetadata(
        name="integrate_samples",
        title="Integrate Multiple Samples",
        description="Integrate multiple spatial datasets",
        read_only_hint=False,
        idempotent_hint=False,
    ),
    "deconvolve_data": MCPToolMetadata(
        name="deconvolve_data",
        title="Spatial Deconvolution",
        description="Deconvolve spatial spots into cell types",
        read_only_hint=False,
        idempotent_hint=False,
        open_world_hint=True,
    ),
    "identify_spatial_domains": MCPToolMetadata(
        name="identify_spatial_domains",
        title="Identify Spatial Domains",
        description="Find spatial domains and niches",
        read_only_hint=False,
        idempotent_hint=False,
    ),
    "analyze_cell_communication": MCPToolMetadata(
        name="analyze_cell_communication",
        title="Cell Communication Analysis",
        description="Analyze cell-cell communication",
        read_only_hint=False,
        idempotent_hint=True,
        open_world_hint=True,
    ),
    "analyze_enrichment": MCPToolMetadata(
        name="analyze_enrichment",
        title="Gene Set Enrichment Analysis",
        description="Perform gene set enrichment analysis",
        read_only_hint=False,
        idempotent_hint=True,
    ),
    "find_spatial_genes": MCPToolMetadata(
        name="find_spatial_genes",
        title="Find Spatial Variable Genes",
        description="Identify spatially variable genes",
        read_only_hint=False,
        idempotent_hint=False,
    ),
    "register_spatial_data": MCPToolMetadata(
        name="register_spatial_data",
        title="Register Spatial Sections",
        description="Align spatial sections",
        read_only_hint=False,
        idempotent_hint=False,
    ),
}

# Update adapter with tool metadata
for name, metadata in tool_metadata.items():
    adapter._tool_metadata[name] = metadata


# ============== Publication Export Tools ==============


@mcp.tool()
@mcp_tool_error_handler()
async def save_data(
    data_id: str,
    output_path: Optional[str] = None,
    context: Context = None,
) -> str:
    """Manually save dataset to disk

    Saves the current state of the dataset including all analysis results
    and metadata to a compressed H5AD file.

    Args:
        data_id: Dataset ID to save
        output_path: Optional custom save path. If not provided, saves to:
                    - CHATSPATIAL_DATA_DIR environment variable location, or
                    - .chatspatial_saved/ directory next to original data

    Returns:
        Path where data was saved

    Examples:
        # Save to default location
        save_data("data1")

        # Save to custom location
        save_data("data1", output_path="/path/to/save/my_analysis.h5ad")

    Note:
        Saved files include all preprocessing, analysis results, and metadata.
        Use CHATSPATIAL_DATA_DIR environment variable for centralized storage.
    """
    from .utils.adata_persistence import save_adata

    # Validate dataset exists
    validate_dataset(data_id)

    if context:
        await context.info(f"Saving dataset '{data_id}'...")

    # Get dataset info
    dataset_info = data_manager.data_store[data_id]
    adata = dataset_info["adata"]
    original_path = dataset_info.get("path", "")

    try:
        if output_path:
            # User specified custom path
            from pathlib import Path

            # Resolve to absolute path to avoid confusion about save location
            save_path = Path(output_path).resolve()
            save_path.parent.mkdir(parents=True, exist_ok=True)
            adata.write_h5ad(save_path, compression="gzip", compression_opts=4)
        else:
            # Use default location
            save_path = save_adata(data_id, adata, original_path)

        # Always return absolute path so user knows exact location
        absolute_path = save_path.resolve()

        if context:
            await context.info(f"Dataset saved to: {absolute_path}")

        return f"Dataset '{data_id}' saved to: {absolute_path}"

    except Exception as e:
        error_msg = f"Failed to save dataset: {str(e)}"
        if context:
            await context.error(error_msg)
        raise


def main():
    """Run the MCP server"""
    import argparse

    parser = argparse.ArgumentParser(description="ChatSpatial MCP Server")
    parser.add_argument(
        "--transport",
        choices=["stdio", "sse"],
        default="stdio",
        help="Transport protocol to use (default: stdio)",
    )

    args = parser.parse_args()

    print(
        f"Starting ChatSpatial server with {args.transport} transport...",
        file=sys.stderr,
    )
    mcp.run(transport=args.transport)


if __name__ == "__main__":
    main()
