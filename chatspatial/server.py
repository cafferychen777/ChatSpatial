"""
Main server implementation for ChatSpatial.
"""

from typing import Dict, Any, List, Optional, Union
import warnings

# Suppress warnings to speed up startup
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from mcp.server.fastmcp import FastMCP, Context
from mcp.server.fastmcp.utilities.types import Image
from .utils.image_utils import fig_to_image, create_placeholder_image
from .mcp_improvements import (
    TOOL_ANNOTATIONS, SPATIAL_PROMPTS, ErrorType, format_mcp_error,
    get_resource_list, read_resource_content, prompt_to_tool_params,
    get_all_tools_with_annotations
)
from .utils.tool_error_handling import mcp_tool_error_handler
from .utils.pydantic_error_handler import mcp_pydantic_error_handler
from .utils.mcp_parameter_handler import (
    manual_parameter_validation,
    validate_analysis_params,
    validate_visualization_params,
    validate_spatial_analysis_params
)

from .models.data import (
    SpatialDataset,
    AnalysisParameters,
    VisualizationParameters,
    AnnotationParameters,
    SpatialAnalysisParameters,
    RNAVelocityParameters,
    TrajectoryParameters,
    IntegrationParameters,
    DeconvolutionParameters,
    SpatialDomainParameters,
    SpatialVariableGenesParameters,
    CellCommunicationParameters
)
from .models.analysis import (
    PreprocessingResult,
    DifferentialExpressionResult,
    AnnotationResult,
    SpatialAnalysisResult,
    RNAVelocityResult,
    TrajectoryResult,
    IntegrationResult,
    DeconvolutionResult,
    SpatialDomainResult,
    SpatialVariableGenesResult,
    CellCommunicationResult
)
from .tools.annotation import annotate_cell_types
from .tools.spatial_analysis import analyze_spatial_patterns
from .tools.differential import differential_expression
from .tools.trajectory import analyze_rna_velocity
from .tools.deconvolution import deconvolve_spatial_data
from .tools.spatial_genes import identify_spatial_genes
from .utils.data_loader import load_spatial_data

# Create MCP server
mcp = FastMCP("ChatSpatial")

# Store for loaded datasets
data_store: Dict[str, Any] = {}

# Global storage for visualization resources
if not hasattr(mcp, '_visualization_resources'):
    mcp._visualization_resources = {}


def validate_dataset(data_id: str) -> None:
    """Validate that a dataset exists in the data store

    Args:
        data_id: Dataset ID

    Raises:
        ValueError: If the dataset is not found
    """
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")


@mcp.tool()
@mcp_tool_error_handler()
async def load_data(
    data_path: str,
    data_type: str = "auto",
    name: Optional[str] = None,
    context: Context = None
) -> SpatialDataset:
    """Load spatial transcriptomics data

    Args:
        data_path: Path to the data file or directory
        data_type: Type of spatial data (auto, 10x_visium, slide_seq, merfish, seqfish, other, h5ad).
                  If 'auto', will try to determine the type from the file extension or directory structure.
        name: Optional name for the dataset

    Returns:
        Dataset information
    """
    if context:
        await context.info(f"Loading data from {data_path} (type: {data_type})")

    # Load data
    dataset_info = await load_spatial_data(data_path, data_type, name)

    if context:
        await context.info(f"Successfully loaded {dataset_info['type']} data with {dataset_info['n_cells']} cells and {dataset_info['n_genes']} genes")

    # Generate unique ID
    data_id = f"data_{len(data_store) + 1}"

    # Store data
    data_store[data_id] = dataset_info

    # Return dataset information
    return SpatialDataset(
        id=data_id,
        name=dataset_info["name"],
        data_type=data_type,
        description=f"Loaded from {data_path}"
    )


@mcp.tool()
@mcp_tool_error_handler()
@manual_parameter_validation(
    ("params", validate_analysis_params)
)
async def preprocess_data(
    data_id: str,
    params: Any = None,
    context: Context = None
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
        - none: No normalization
        - scvi: Use scVI for normalization and dimensionality reduction
        
        When use_scvi_preprocessing=True, scVI will be used for advanced preprocessing
        including denoising and batch effect correction.
    """
    # Import to avoid name conflict
    from .tools.preprocessing import preprocess_data as preprocess_func

    # Validate dataset
    validate_dataset(data_id)

    # Call preprocessing function
    result = await preprocess_func(data_id, data_store, params, context)

    return result


@mcp.tool()
@mcp_tool_error_handler()
@manual_parameter_validation(
    ("params", validate_visualization_params)
)
async def visualize_data(
    data_id: str,
    params: Any = None,
    context: Context = None
):
    """Visualize spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Visualization parameters including:
            - plot_type: Type of visualization (spatial, heatmap, violin, umap,
                        spatial_domains, cell_communication, deconvolution, trajectory,
                        spatial_analysis, multi_gene, lr_pairs, gene_correlation,
                        gaston_isodepth, gaston_domains, gaston_genes)
            - feature: Gene or feature to visualize (for spatial plots)
            - colormap: Color scheme for visualization
            - figure_size: Size of the output figure

    Returns:
        Visualization image
    """
    # Import to avoid name conflict
    from .tools.visualization import visualize_data as visualize_func

    # Validate dataset
    validate_dataset(data_id)

    # Handle different parameter formats for backward compatibility
    if isinstance(params, str):
        # Handle simple string format like "gene:CCL21"
        if params.startswith("gene:"):
            feature = params.split(":", 1)[1]
            params = VisualizationParameters(feature=feature, plot_type="spatial")
        else:
            # Treat as feature name
            params = VisualizationParameters(feature=params, plot_type="spatial")
    elif isinstance(params, dict):
        # Handle dictionary format
        params = VisualizationParameters(**params)
    elif params is None:
        # Handle None
        params = VisualizationParameters()

    # Call visualization function
    image = await visualize_func(data_id, data_store, params, context)

    # Save image as MCP Resource and return reference
    if image is not None:
        import time
        from pathlib import Path

        # Create visualization resources directory
        resources_dir = Path("visualization_resources")
        resources_dir.mkdir(exist_ok=True)

        # Generate unique resource URI
        timestamp = int(time.time())
        resource_id = f"viz_{data_id}_{params.plot_type}_{timestamp}"
        resource_uri = f"visualization://{resource_id}"
        resource_file = resources_dir / f"{resource_id}.png"

        # Save image to file
        with open(resource_file, 'wb') as f:
            f.write(image.data)

        # Store resource info globally for MCP resource handler
        if not hasattr(mcp, '_visualization_resources'):
            mcp._visualization_resources = {}

        mcp._visualization_resources[resource_uri] = {
            'file_path': str(resource_file),
            'name': f"{params.plot_type} - {getattr(params, 'feature', 'N/A')}",
            'description': f"Visualization of {data_id}",
            'mime_type': 'image/png',
            'size': len(image.data)
        }

        file_size_kb = len(image.data) / 1024

        if context:
            await context.info(f"Image saved as resource: {resource_uri} ({file_size_kb:.1f}KB)")

        return f"""✅ 可视化已成功生成！

📊 **可视化信息**:
- 类型: {params.plot_type}
- 特征: {getattr(params, 'feature', 'N/A')}
- 数据集: {data_id}
- 图像大小: {file_size_kb:.1f} KB

🖼️ **图像资源**: `{resource_uri}`

💡 **提示**: 图像已保存为 MCP 资源。在 Claude Desktop 中，你可以通过资源面板访问和查看此图像。资源 URI: {resource_uri}"""

    else:
        # Return error message if no image was generated
        return "❌ 可视化生成失败，请检查数据和参数设置。"


@mcp.tool()
@mcp_tool_error_handler()
async def annotate_cells(
    data_id: str,
    params: AnnotationParameters = AnnotationParameters(),
    context: Context = None
) -> AnnotationResult:
    """Annotate cell types in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Annotation parameters

    Returns:
        Annotation result with cell type information and optional visualization

    Notes:
        Available methods:
        - marker_genes: Use known marker genes for annotation
        - tangram: Maps single-cell data to spatial data (requires reference_data_id)
        - scanvi: Semi-supervised annotation using scANVI from scvi-tools (requires reference_data_id)
        - cellassign: Probabilistic cell type assignment using marker genes (from scvi-tools)
        - correlation, supervised, popv, gptcelltype, scrgcl: Other available methods
        
        For methods requiring reference data (tangram, scanvi), reference_data_id must point to a loaded single-cell dataset.
    """
    # Validate dataset
    validate_dataset(data_id)

    # Validate reference data for methods that require it
    if params.method in ["tangram", "scanvi"] and params.reference_data_id:
        if params.reference_data_id not in data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")
        if context:
            await context.info(f"Using reference dataset {params.reference_data_id} for {params.method} method")
    elif params.method in ["tangram", "scanvi"] and not params.reference_data_id:
        raise ValueError(f"{params.method} method requires a reference_data_id")

    # Call annotation function
    result = await annotate_cell_types(data_id, data_store, params, context)

    # Log results
    if context:
        await context.info(f"Annotation completed with {len(result.cell_types)} cell types identified")
        if params.method == "tangram" and result.tangram_mapping_score is not None:
            await context.info(f"Tangram mapping score: {result.tangram_mapping_score:.4f}")

    return result


@mcp.tool()
@mcp_tool_error_handler()
@mcp_pydantic_error_handler()
async def analyze_spatial_data(
    data_id: str,
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Context = None
) -> SpatialAnalysisResult:
    """Perform spatial analysis on transcriptomics data and return the results.

    This tool runs the specified spatial analysis (e.g., neighborhood, Moran's I)
    and stores the results in the AnnData object. It returns a summary of the analysis.
    To visualize the results, use the 'visualize_data' tool with
    plot_type='spatial_analysis' and the corresponding 'analysis_sub_type'.

    Args:
        data_id: Dataset ID
        params: Spatial analysis parameters

    Returns:
        A SpatialAnalysisResult object containing statistics and metadata.
    """
    # Validate dataset
    validate_dataset(data_id)

    # Call spatial analysis unified function (without return_type parameter)
    result = await analyze_spatial_patterns(data_id, data_store, params, context)

    # Log and return the result object
    if context:
        await context.info(f"Spatial analysis '{params.analysis_type}' completed.")
        await context.info(f"Results stored in adata.uns['{result.statistics.get('analysis_key_in_adata', 'N/A')}']")
        await context.info("To visualize the results, use the 'visualize_data' tool with:")
        await context.info(f"  plot_type='spatial_analysis', analysis_sub_type='{params.analysis_type}'")

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def find_markers(
    data_id: str,
    group_key: str,
    group1: str,
    group2: str,
    n_top_genes: int = 50,
    method: str = "wilcoxon",
    context: Context = None
) -> DifferentialExpressionResult:
    """Find marker genes between groups

    Args:
        data_id: Dataset ID
        group_key: Key in adata.obs for grouping cells
        group1: First group for comparison
        group2: Second group for comparison
        n_top_genes: Number of top differentially expressed genes to return
        method: Statistical method for DE analysis

    Returns:
        Differential expression analysis result
    """
    # Validate dataset
    validate_dataset(data_id)

    # Call differential expression function
    return await differential_expression(
        data_id, data_store, group_key, group1, group2, n_top_genes, method, context
    )


# Add RNA velocity analysis tool
@mcp.tool()
@mcp_tool_error_handler()
async def analyze_velocity_data(
    data_id: str,
    params: RNAVelocityParameters = RNAVelocityParameters(),
    context: Context = None
) -> RNAVelocityResult:
    """Analyze RNA velocity in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: RNA velocity analysis parameters

    Returns:
        RNA velocity analysis result
    """
    # Validate dataset
    validate_dataset(data_id)

    # Call RNA velocity analysis function
    result = await analyze_rna_velocity(data_id, data_store, params, context)

    # Log results
    if context:
        if result.velocity_computed:
            await context.info(f"RNA velocity analysis completed using {params.mode} mode.")
            await context.info("To visualize results, use the 'visualize_data' tool with plot_type='spatial' and feature='velocity'.")
        else:
            await context.warning("Could not compute RNA velocity.")

    return result


# Add trajectory analysis tool
@mcp.tool()
@mcp_tool_error_handler()
async def analyze_trajectory_data(
    data_id: str,
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Context = None
) -> TrajectoryResult:
    """Analyze trajectory and return results

    This function performs trajectory analysis and returns computation results.
    For visualization, use the 'visualize_data' tool with plot_type='trajectory'.

    Args:
        data_id: Dataset ID
        params: Trajectory analysis parameters
        context: MCP context

    Returns:
        Trajectory analysis result
        
    Notes:
        Available methods:
        - cellrank: Uses CellRank for trajectory inference (default)
        - palantir: Uses Palantir for trajectory inference
        - velovi: Uses VeloVI from scvi-tools for velocity-based trajectory inference
        
        VeloVI provides probabilistic velocity modeling and is particularly useful
        for spatial transcriptomics data with velocity information.
    """
    # Validate dataset
    validate_dataset(data_id)

    # Import the trajectory analysis function
    from .tools.trajectory import analyze_trajectory

    # Call trajectory analysis function
    result = await analyze_trajectory(data_id, data_store, params, context)

    # Log results
    if context:
        await context.info(f"Trajectory analysis with method '{params.method}' completed.")
        await context.info(f"Results stored in adata.obs['{result.pseudotime_key}'].")
        await context.info("To visualize results, use the 'visualize_data' tool with plot_type='trajectory'.")

    return result


# Add tool for integrating multiple samples
@mcp.tool()
@mcp_tool_error_handler()
async def integrate_samples(
    data_ids: List[str],
    params: IntegrationParameters = IntegrationParameters(),
    context: Context = None
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples

    Args:
        data_ids: List of dataset IDs to integrate
        params: Integration parameters

    Returns:
        Integration result with visualizations
        
    Notes:
        Available methods:
        - harmony: Uses Harmony for batch effect correction
        - bbknn: Uses BBKNN for batch-aware k-NN graph construction
        - scanorama: Uses Scanorama for integration
        - mnn: Uses mutual nearest neighbors for integration
        - scvi: Uses scVI from scvi-tools for probabilistic integration
        - multivi: Uses MultiVI from scvi-tools for multi-modal integration
        - totalvi: Uses TotalVI from scvi-tools for protein+RNA integration
    """
    if context:
        await context.info(f"Integrating {len(data_ids)} samples using {params.method} method")

    # Validate all dataset IDs
    for data_id in data_ids:
        validate_dataset(data_id)

    # Call integration function
    from .tools.integration import integrate_samples as integrate_func
    return await integrate_func(data_ids, data_store, params, context)


# Add tool for spatial deconvolution
@mcp.tool()
@mcp_tool_error_handler()
async def deconvolve_data(
    data_id: str,
    params: DeconvolutionParameters = DeconvolutionParameters(),
    context: Context = None
) -> DeconvolutionResult:
    """Deconvolve spatial transcriptomics data to estimate cell type proportions

    Args:
        data_id: Dataset ID
        params: Deconvolution parameters

    Returns:
        Deconvolution result with cell type proportions and visualization
        
    Notes:
        Available methods:
        - cell2location: Uses Cell2location for spatial deconvolution (requires reference_data_id)
        - spotiphy: Uses Spotiphy for spatial deconvolution (requires reference_data_id)
        - rctd: Uses RCTD (Robust Cell Type Decomposition) (requires reference_data_id)
        - destvi: Uses DestVI from scvi-tools for spatial deconvolution (requires reference_data_id)
        - stereoscope: Uses Stereoscope from scvi-tools for spatial deconvolution (requires reference_data_id)
        
        All methods require a reference single-cell dataset specified by reference_data_id.
    """
    try:
        if context:
            await context.info(f"Deconvolving spatial data using {params.method} method")
            await context.info(f"Parameters: {params}")

        # Validate dataset
        validate_dataset(data_id)

        # Validate reference data if provided
        if params.reference_data_id and params.reference_data_id not in data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")

        # Call deconvolution function
        result = await deconvolve_spatial_data(data_id, data_store, params, context)

        # Automatically visualize the results if visualization_params are provided
        if result.visualization_params:
            if context:
                await context.info("Automatically visualizing deconvolution results")

            # Create visualization parameters
            vis_params = VisualizationParameters(**result.visualization_params)

            # Call visualization function
            try:
                # This will create a multi-panel figure with top cell types
                from .tools.visualization import visualize_data as visualize_func
                await visualize_func(data_id, data_store, vis_params, context)
            except Exception as e:
                if context:
                    await context.warning(f"Failed to automatically visualize deconvolution results: {str(e)}")

        return result
    except Exception as e:
        error_msg = f"Error in deconvolution: {str(e)}"
        if context:
            await context.warning(error_msg)
            await context.info("If you're having issues with parameter formatting, try using this format:")
            await context.info('deconvolve_data(data_id="data_1", params={"method": "cell2location", "reference_data_id": "data_2", "cell_type_key": "CellType"})')
            await context.info("Available methods: cell2location, spotiphy, rctd, destvi, stereoscope")
        raise


# Add tool for spatial domain identification
@mcp.tool()
@mcp_tool_error_handler()
async def identify_spatial_domains(
    data_id: str,
    params: SpatialDomainParameters = SpatialDomainParameters(),
    context: Context = None
) -> SpatialDomainResult:
    """Identify spatial domains in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Spatial domain identification parameters

    Returns:
        Spatial domain identification result with domain labels and statistics
    """
    if context:
        await context.info(f"Identifying spatial domains using {params.method} method")

    # Validate dataset
    validate_dataset(data_id)

    # Call spatial domain identification function
    from .tools.spatial_domains import identify_spatial_domains as identify_domains_func
    result = await identify_domains_func(data_id, data_store, params, context)

    if context:
        await context.info(f"Successfully identified {result.n_domains} spatial domains")
        await context.info(f"Domain labels stored in adata.obs['{result.domain_key}']")
        if result.refined_domain_key:
            await context.info(f"Refined domain labels stored in adata.obs['{result.refined_domain_key}']")

    return result


# Add tool for cell communication analysis
@mcp.tool()
@mcp_tool_error_handler()
async def analyze_cell_communication(
    data_id: str,
    params: CellCommunicationParameters = CellCommunicationParameters(),
    context: Context = None
) -> CellCommunicationResult:
    """Analyze cell-cell communication in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Cell communication analysis parameters

    Returns:
        Cell communication analysis result with LR pairs and communication networks
    """
    if context:
        await context.info(f"Analyzing cell communication using {params.method} method")

    # Validate dataset
    validate_dataset(data_id)

    # Call cell communication analysis function
    from .tools.cell_communication import analyze_cell_communication as analyze_comm_func
    result = await analyze_comm_func(data_id, data_store, params, context)

    if context:
        await context.info(f"Successfully analyzed {result.n_significant_pairs} significant LR pairs")
        if result.global_results_key:
            await context.info(f"Global results stored in adata.uns['{result.global_results_key}']")
        if result.local_analysis_performed and result.local_results_key:
            await context.info(f"Local results stored in adata.uns['{result.local_results_key}']")
        if result.top_lr_pairs:
            await context.info(f"Top LR pair: {result.top_lr_pairs[0]}")

    return result


# Add tool for enrichment analysis using EnrichMap
@mcp.tool()
@mcp_tool_error_handler()
async def analyze_enrichment(
    data_id: str,
    gene_sets: Union[List[str], Dict[str, List[str]]],
    score_keys: Optional[Union[str, List[str]]] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = True,
    batch_key: Optional[str] = None,
    context: Context = None
) -> Dict[str, Any]:
    """Perform spatially-aware gene set enrichment analysis using EnrichMap
    
    EnrichMap computes and visualizes enrichment scores of gene sets in spatial transcriptomics
    datasets with spatial smoothing and covariate correction.
    
    Args:
        data_id: Dataset ID
        gene_sets: Either a single gene list or a dictionary of gene sets where keys are 
                  signature names and values are lists of genes
        score_keys: Names for the gene signatures if gene_sets is a list
        spatial_key: Key in adata.obsm containing spatial coordinates (default: "spatial")
        n_neighbors: Number of nearest spatial neighbors for smoothing (default: 6)
        smoothing: Whether to perform spatial smoothing (default: True)
        correct_spatial_covariates: Whether to correct for spatial covariates using GAM (default: True)
        batch_key: Column in adata.obs for batch-wise normalization
    
    Returns:
        Enrichment analysis result with scores, gene contributions, and statistics
        
    Examples:
        # Single gene set
        analyze_enrichment(data_id="data_1", gene_sets=["CD3D", "CD3E", "CD8A"], score_keys="T_cell")
        
        # Multiple gene sets
        analyze_enrichment(data_id="data_1", gene_sets={
            "T_cell": ["CD3D", "CD3E", "CD8A"],
            "B_cell": ["CD19", "MS4A1", "CD79A"]
        })
    """
    if context:
        await context.info("Performing spatially-aware enrichment analysis using EnrichMap")
    
    # Validate dataset
    validate_dataset(data_id)
    
    # Import and call enrichment analysis function
    from .tools.enrichment_analysis import perform_enrichment_analysis
    result = await perform_enrichment_analysis(
        data_id=data_id,
        data_store=data_store,
        gene_sets=gene_sets,
        score_keys=score_keys,
        spatial_key=spatial_key,
        n_neighbors=n_neighbors,
        smoothing=smoothing,
        correct_spatial_covariates=correct_spatial_covariates,
        batch_key=batch_key,
        context=context
    )
    
    if context:
        await context.info(f"Successfully computed enrichment scores for {len(result['signatures'])} signatures")
        for sig in result['signatures']:
            stats = result['summary_stats'][sig]
            await context.info(f"  {sig}: mean={stats['mean']:.3f}, std={stats['std']:.3f}, n_genes={stats['n_genes']}")
    
    return result


# Add tool for spatial variable genes identification using GASTON
@mcp.tool()
@mcp_tool_error_handler()
async def find_spatial_genes(
    data_id: str,
    params: SpatialVariableGenesParameters = SpatialVariableGenesParameters(),
    context: Context = None
) -> SpatialVariableGenesResult:
    """Identify spatial variable genes using GASTON method

    GASTON (Generative Adversarial Spatial Transcriptomics Optimization Network) learns
    a topographic map of tissue slices by modeling gene expression through isodepth coordinates.

    Key features:
    - Learns 1D isodepth coordinate that varies smoothly across tissue
    - Identifies spatial domains and continuous/discontinuous gene expression patterns
    - Uses interpretable deep learning with neural networks
    - Supports GLM-PCA or Pearson residuals preprocessing

    Args:
        data_id: Dataset ID
        params: GASTON analysis parameters including:
            - preprocessing_method: "glmpca" (recommended) or "pearson_residuals"
            - spatial_hidden_layers: Architecture for spatial embedding network
            - expression_hidden_layers: Architecture for expression function network
            - epochs: Number of training epochs
            - learning_rate: Learning rate for optimization
            - use_positional_encoding: Whether to use positional encoding
            - n_domains: Number of spatial domains to identify (default: 5)
            - num_bins: Number of bins for isodepth binning (default: 70)
            - continuous_quantile: Quantile threshold for continuous genes (default: 0.9)
            - discontinuous_quantile: Quantile threshold for discontinuous genes (default: 0.9)
            - umi_threshold: Minimum UMI count threshold for genes (default: 500)

    Returns:
        Spatial variable genes identification result with spatial domains and gene classification.
        Use visualize_data with plot_type="gaston_isodepth", "gaston_domains", or "gaston_genes"
        to visualize the results.
    """
    if context:
        await context.info(f"Identifying spatial variable genes using GASTON method")

    # Validate dataset
    validate_dataset(data_id)

    # Call GASTON spatial variable genes identification function
    result = await identify_spatial_genes(data_id, data_store, params, context)

    if context:
        await context.info(f"Successfully completed GASTON analysis")
        await context.info(f"Identified {result.n_spatial_domains} spatial domains")
        await context.info(f"Found {result.n_continuous_genes} genes with continuous gradients")
        await context.info(f"Found {result.n_discontinuous_genes} genes with discontinuities")
        await context.info(f"Final training loss: {result.final_loss:.6f}")
        await context.info(f"Model R²: {result.model_performance.get('r2', 'N/A')}")
        await context.info(f"Results stored with keys: {result.isodepth_key}, {result.spatial_domains_key}")
        await context.info("Use visualize_data tool with plot_type='gaston_isodepth', 'gaston_domains', or 'gaston_genes' to visualize results")

    return result


# Add resource for dataset information
@mcp.resource("dataset://{data_id}")
def get_dataset_info(data_id: str) -> Dict[str, Any]:
    """Get information about a dataset"""
    # Validate dataset
    validate_dataset(data_id)

    return data_store[data_id]


# The list_resources handler is automatically provided by FastMCP
# We just need to register resources using mcp.add_resource()


# # Implement read_resource handler
# @mcp.read_resource
# async def read_resource(uri: str) -> str:
#     """Read resource content"""
#     try:
#         return read_resource_content(uri, data_store)
#     except Exception as e:
#         error = format_mcp_error(
#             ErrorType.INTERNAL_ERROR,
#             f"Failed to read resource: {str(e)}",
#             {"uri": uri}
#         )
#         raise ValueError(error)


# # # Implement list_prompts handler
# # @mcp.list_prompts
# # async def list_prompts() -> List[Dict[str, Any]]:
# #     """List available prompts for spatial analysis"""
# #     return [
# #         {
# #             "name": p.name,
# #             "description": p.description,
# #             "arguments": p.arguments
# #         }
# #         for p in SPATIAL_PROMPTS
# #     ]
# 
# 
# # # Implement get_prompt handler
# # @mcp.get_prompt()
# # async def get_prompt(name: str, arguments: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
#     """Get a specific prompt with filled arguments"""
#     # Find the prompt
#     prompt = next((p for p in SPATIAL_PROMPTS if p.name == name), None)
#     if not prompt:
#         error = format_mcp_error(
#             ErrorType.INVALID_REQUEST,
#             f"Prompt '{name}' not found"
#         )
#         raise ValueError(error)
#     
#     # Convert to tool parameters
#     try:
#         if arguments:
#             tool_params = prompt_to_tool_params(name, arguments, data_store)
#             # Format as a helpful message
#             tool_name = tool_params["tool"]
#             params = tool_params["params"]
#             
#             # Create a formatted prompt message
#             messages = [
#                 {
#                     "role": "user",
#                     "content": f"Please run the {tool_name} tool with these parameters: {params}"
#                 }
#             ]
#         else:
#             # Return the prompt template
#             messages = [
#                 {
#                     "role": "user", 
#                     "content": f"Please help me with: {prompt.description}"
#                 }
#             ]
#         
#         return {
#             "description": prompt.description,
#             "messages": messages
#         }
#     except Exception as e:
#         error = format_mcp_error(
#             ErrorType.INTERNAL_ERROR,
#             f"Failed to process prompt: {str(e)}"
#         )
#         raise ValueError(error)
# 
# 
# # # Override list_tools to include annotations
# # @mcp.list_tools
# # async def list_tools() -> List[Dict[str, Any]]:
#     """List all available tools with annotations"""
#     # Get the default tools list from FastMCP
#     tools = []
#     
#     # Get registered tools from mcp instance
#     for tool_name, tool_info in mcp._tool_handlers.items():
#         # Get annotations if available
#         annotations = TOOL_ANNOTATIONS.get(tool_name, {})
        
#         # Create tool info with annotations
#         tool_data = {
#             "name": tool_name,
#             "description": tool_info.description or "",
#             "inputSchema": tool_info.parameters_schema
#         }
#         
#         # Add annotations if available
#         if annotations:
#             tool_data["annotations"] = {
#                 "title": annotations.get("title", tool_name),
#                 "readOnlyHint": annotations.get("readOnlyHint", False),
#                 "destructiveHint": annotations.get("destructiveHint", False),
#                 "idempotentHint": annotations.get("idempotentHint", False),
#                 "openWorldHint": annotations.get("openWorldHint", False)
#             }
#         
#         tools.append(tool_data)
#     
#     return tools


# MCP Resources handlers
@mcp.resource()
async def list_visualization_resources():
    """List available visualization resources"""
    resources = []

    if hasattr(mcp, '_visualization_resources'):
        for uri, info in mcp._visualization_resources.items():
            resources.append({
                "uri": uri,
                "name": info['name'],
                "description": info['description'],
                "mimeType": info['mime_type'],
                "size": info['size']
            })

    return resources

@mcp.resource()
async def read_visualization_resource(uri: str):
    """Read a visualization resource"""
    if not hasattr(mcp, '_visualization_resources'):
        raise ValueError(f"Resource not found: {uri}")

    if uri not in mcp._visualization_resources:
        raise ValueError(f"Resource not found: {uri}")

    resource_info = mcp._visualization_resources[uri]
    file_path = resource_info['file_path']

    try:
        with open(file_path, 'rb') as f:
            image_data = f.read()

        import base64
        base64_data = base64.b64encode(image_data).decode('utf-8')

        return {
            "uri": uri,
            "mimeType": resource_info['mime_type'],
            "blob": base64_data
        }
    except FileNotFoundError:
        raise ValueError(f"Resource file not found: {file_path}")


def main():
    """Run the MCP server"""
    import argparse
    
    parser = argparse.ArgumentParser(description="ChatSpatial MCP Server")
    parser.add_argument(
        "--transport",
        choices=["stdio", "sse"],
        default="stdio",
        help="Transport protocol to use (default: stdio)"
    )
    
    args = parser.parse_args()
    
    print(f"Starting ChatSpatial server with {args.transport} transport...")
    mcp.run(transport=args.transport)


if __name__ == "__main__":
    main()