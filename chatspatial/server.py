"""
Main server implementation for ChatSpatial using the Spatial MCP Adapter.
"""

from typing import Dict, Any, List, Optional, Union
import warnings
import asyncio
import logging

# Suppress warnings to speed up startup
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image

from .spatial_mcp_adapter import (
    create_spatial_mcp_server,
    SpatialMCPAdapter,
    DefaultSpatialDataManager,
    MCPToolMetadata
)

from .utils.tool_error_handling import mcp_tool_error_handler
from .utils.pydantic_error_handler import mcp_pydantic_error_handler
from .utils.error_handling import ProcessingError
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
    CellCommunicationParameters,
    EnrichmentParameters
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
    CellCommunicationResult,
    EnrichmentResult
)
from .tools.annotation import annotate_cell_types
from .tools.spatial_analysis import analyze_spatial_patterns
from .tools.differential import differential_expression
from .tools.trajectory import analyze_rna_velocity
from .tools.deconvolution import deconvolve_spatial_data
from .tools.spatial_genes import identify_spatial_genes
from .utils.data_loader import load_spatial_data

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

    # Load data using data manager
    data_id = await data_manager.load_dataset(data_path, data_type, name)
    dataset_info = await data_manager.get_dataset(data_id)

    if context:
        await context.info(f"Successfully loaded {dataset_info['type']} data with {dataset_info['n_cells']} cells and {dataset_info['n_genes']} genes")

    # Create resource for the dataset
    await adapter.resource_manager.create_dataset_resource(data_id, dataset_info)

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
@mcp_tool_error_handler()  # ⚠️ CRITICAL: This decorator has special Image handling - see /docs/CRITICAL_IMAGE_DISPLAY_BUG.md
@manual_parameter_validation(
    ("params", validate_visualization_params)
)
async def visualize_data(
    data_id: str,
    params: Any = None,
    context: Context = None
) -> Image:  # ⚠️ CRITICAL: Must return Image object, NOT ImageContent or dict!
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

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

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

    # Create visualization resource and return the image
    if image is not None:
        import time
        viz_id = f"{data_id}_{params.plot_type}_{int(time.time())}"
        
        metadata = {
            "data_id": data_id,
            "plot_type": params.plot_type,
            "feature": getattr(params, 'feature', 'N/A'),
            "timestamp": int(time.time()),
            "name": f"{params.plot_type} - {getattr(params, 'feature', 'N/A')}",
            "description": f"Visualization of {data_id}"
        }
        
        await adapter.resource_manager.create_visualization_resource(
            viz_id, image.data, metadata
        )

        file_size_kb = len(image.data) / 1024

        if context:
            await context.info(f"Image saved as resource: spatial://visualizations/{viz_id} ({file_size_kb:.1f}KB)")
            await context.info(f"Visualization type: {params.plot_type}, feature: {getattr(params, 'feature', 'N/A')}")

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

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Validate reference data for methods that require it
    if params.method in ["tangram", "scanvi"] and params.reference_data_id:
        if params.reference_data_id not in data_manager.data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")
        ref_info = await data_manager.get_dataset(params.reference_data_id)
        data_store[params.reference_data_id] = ref_info

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
@manual_parameter_validation(
    ("params", validate_spatial_analysis_params)
)
async def analyze_spatial_data(
    data_id: str,
    params: Any = None,
    context: Context = None
) -> SpatialAnalysisResult:
    """Analyze spatial patterns and relationships in the data

    Args:
        data_id: Dataset ID
        params: Analysis parameters

    Returns:
        Spatial analysis result with statistics and optional visualization

    Notes:
        Available analysis types:
        - spatial_autocorrelation: Moran's I for spatial patterns
        - nearest_neighbor: Analyze nearest neighbor relationships
        - co_expression: Spatial co-expression analysis
        - spatial_enrichment: Enrichment of features in spatial regions
        - interaction_analysis: Cell-cell interaction based on proximity
        - ligand_receptor: Ligand-receptor interaction analysis
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Call spatial analysis function
    result = await analyze_spatial_patterns(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save spatial analysis result
    await data_manager.save_result(data_id, "spatial_analysis", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "spatial_analysis", result)

    # Create visualization if available
    if params and getattr(params, 'visualize', True):
        vis_image = await adapter.create_visualization_from_result(
            data_id, "spatial_analysis", result, context
        )
        if vis_image:
            result.visualization = vis_image

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def find_markers(
    data_id: str,
    group_key: str,
    group1: Optional[str] = None,
    group2: Optional[str] = None,
    method: str = "wilcoxon",
    n_genes: int = 25,
    context: Context = None
) -> DifferentialExpressionResult:
    """Find differentially expressed genes between groups

    Args:
        data_id: Dataset ID
        group_key: Column name defining groups
        group1: First group (if None, compare against all others)
        group2: Second group (if None, compare group1 against all others)
        method: Statistical test method
        n_genes: Number of top genes to return

    Returns:
        Differential expression result with top marker genes
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Call differential expression function
    result = await differential_expression(
        data_id=data_id,
        data_store=data_store,
        group_key=group_key,
        group1=group1,
        group2=group2,
        method=method,
        n_top_genes=n_genes,
        context=context
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
async def analyze_velocity_data(
    data_id: str,
    params: RNAVelocityParameters = RNAVelocityParameters(),
    context: Context = None
) -> RNAVelocityResult:
    """Analyze RNA velocity to understand cellular dynamics

    Args:
        data_id: Dataset ID
        params: RNA velocity parameters

    Returns:
        RNA velocity analysis result
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

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
    context: Context = None
) -> TrajectoryResult:
    """Infer cellular trajectories and pseudotime

    Args:
        data_id: Dataset ID
        params: Trajectory analysis parameters

    Returns:
        Trajectory analysis result

    Notes:
        Available methods:
        - paga: Partition-based graph abstraction
        - palantir: Probabilistic trajectory inference
        - cellrank: RNA velocity-based trajectory inference
        - dpt: Diffusion pseudotime
        - gaston: GASTON spatial trajectory analysis
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
    context: Context = None
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples

    Args:
        data_ids: List of dataset IDs to integrate
        params: Integration parameters

    Returns:
        Integration result with integrated dataset ID

    Notes:
        Available methods:
        - harmony: Fast batch correction
        - scvi: Deep learning-based integration
        - combat: ComBat batch correction
        - mnn: Mutual nearest neighbors
        - PASTE: Probabilistic Alignment of Spatial Transcriptomics Experiments
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
    integrated_id = result.integrated_data_id
    if integrated_id and integrated_id in data_store:
        data_manager.data_store[integrated_id] = data_store[integrated_id]
        
        # Create resource for integrated dataset
        await adapter.resource_manager.create_dataset_resource(
            integrated_id, data_store[integrated_id]
        )

    # Save integration result
    await data_manager.save_result(integrated_id, "integration", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(integrated_id, "integration", result)

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def deconvolve_data(
    data_id: str,
    params: DeconvolutionParameters = DeconvolutionParameters(),
    context: Context = None
) -> DeconvolutionResult:
    """Deconvolve spatial spots to estimate cell type proportions

    Args:
        data_id: Dataset ID
        params: Deconvolution parameters

    Returns:
        Deconvolution result with cell type proportions

    Notes:
        Available methods:
        - spotlight: SPOTlight deconvolution
        - cell2location: Probabilistic mapping
        - rctd: Robust Cell Type Decomposition
        - stereoscope: Probabilistic model
        - destvi: DestVI from scvi-tools
        - tangram: Tangram spatial mapping
        - mrvi: MRVI multi-resolution analysis
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

    # Call deconvolution function
    result = await deconvolve_spatial_data(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save deconvolution result
    await data_manager.save_result(data_id, "deconvolution", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "deconvolution", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def identify_spatial_domains(
    data_id: str,
    params: SpatialDomainParameters = SpatialDomainParameters(),
    context: Context = None
) -> SpatialDomainResult:
    """Identify spatial domains and tissue architecture

    Args:
        data_id: Dataset ID
        params: Spatial domain parameters

    Returns:
        Spatial domain result with identified domains

    Notes:
        Available methods:
        - stlearn: stLearn spatial domains
        - spagcn: SpaGCN graph convolutional network
        - sedr: SEDR deep learning method
        - bayesspace: BayesSpace Bayesian method
        - banksy: BANKSY spatial domains
        - stagate: STAGATE spatial domains
        - gaston: GASTON spatial domains
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
async def analyze_cell_communication(
    data_id: str,
    params: CellCommunicationParameters = CellCommunicationParameters(),
    context: Context = None
) -> CellCommunicationResult:
    """Analyze cell-cell communication patterns

    Args:
        data_id: Dataset ID
        params: Cell communication parameters

    Returns:
        Cell communication analysis result

    Notes:
        Available methods:
        - liana: LIANA framework with multiple methods
        - cellphonedb: CellPhoneDB statistical framework
        - cellchat: CellChat network analysis
        - nichenet: NicheNet ligand-target analysis
        - connectome: Connectome weighted networks
        - cytotalk: CytoTalk crosstalk analysis
        - squidpy: Squidpy permutation test
    """
    # Import cell communication function
    from .tools.cell_communication import analyze_cell_communication as analyze_comm_func

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
    await adapter.resource_manager.create_result_resource(data_id, "communication", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def analyze_enrichment(
    data_id: str,
    params: EnrichmentParameters = EnrichmentParameters(),
    context: Context = None
) -> EnrichmentResult:
    """Perform gene set enrichment analysis

    Args:
        data_id: Dataset ID
        params: Enrichment analysis parameters

    Returns:
        Enrichment analysis result

    Notes:
        Available methods:
        - gsea: Gene Set Enrichment Analysis
        - ora: Over-representation analysis
        - enrichr: Enrichr web service
        - enrichmap: EnrichMap spatial enrichment
        
        Available gene sets:
        - GO_Molecular_Function, GO_Biological_Process, GO_Cellular_Component
        - KEGG, Reactome, WikiPathways
        - MSigDB collections (H, C1-C8)
        - Custom gene sets via gene_sets parameter
    """
    # Import enrichment analysis function
    from .tools.enrichment_analysis import perform_enrichment_analysis
    import time

    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Start timing
    start_time = time.time()

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
            from .tools.gene_set_loader import load_gene_sets
        except ImportError:
            raise ProcessingError("gseapy package is required for gene set loading. Install with: pip install gseapy")
        try:
            # Determine species from data if available
            species = adata.uns.get('species', 'human')
            
            gene_sets = await load_gene_sets(
                database=params.gene_set_database,
                species=species,
                min_genes=params.min_genes,
                max_genes=params.max_genes if hasattr(params, 'max_genes') else 500,
                context=context
            )
            
            if context:
                await context.info(f"Loaded {len(gene_sets)} gene sets from {params.gene_set_database}")
                
        except Exception as e:
            # Fallback: use highly variable genes as a default gene set
            if context:
                await context.info(f"Failed to load gene sets from {params.gene_set_database}: {e}")
                await context.info("Using highly variable genes as fallback")
            
            # Get highly variable genes if available
            if 'highly_variable' in adata.var.columns:
                hvg_genes = adata.var_names[adata.var.highly_variable].tolist()
                if len(hvg_genes) >= params.min_genes:
                    gene_sets = {"highly_variable_genes": hvg_genes[:200]}  # Use top 200 HVGs
                else:
                    # Use all available genes as last resort
                    gene_sets = {"all_genes": adata.var_names.tolist()[:500]}
            else:
                # Use top genes by mean expression
                import numpy as np
                gene_means = np.array(adata.X.mean(axis=0)).flatten()
                top_gene_indices = np.argsort(gene_means)[-500:]
                gene_sets = {"top_expressed_genes": adata.var_names[top_gene_indices].tolist()}
    
    # If still no gene sets, use default strategy
    if gene_sets is None:
        if context:
            await context.info("No gene sets provided or loaded. Using default gene sets based on data.")
        
        # Use highly variable genes or top expressed genes
        if 'highly_variable' in adata.var.columns:
            hvg_genes = adata.var_names[adata.var.highly_variable].tolist()
            gene_sets = {"highly_variable_genes": hvg_genes[:200]}
        else:
            import numpy as np
            gene_means = np.array(adata.X.mean(axis=0)).flatten()
            top_gene_indices = np.argsort(gene_means)[-200:]
            gene_sets = {"top_expressed_genes": adata.var_names[top_gene_indices].tolist()}
    
    if gene_sets is None or len(gene_sets) == 0:
        raise ValueError("No valid gene sets available for enrichment analysis")

    # Call appropriate enrichment function based on method
    if params.method == "enrichmap":
        # Spatial enrichment analysis using EnrichMap
        result_dict = await perform_enrichment_analysis(
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
            context=context
        )
    else:
        # Generic enrichment analysis (GSEA, ORA, ssGSEA, Enrichr)
        from .tools.generic_enrichment import perform_gsea, perform_ora, perform_ssgsea, perform_enrichr
        
        if params.method == "gsea":
            result_dict = await perform_gsea(
                adata=adata,
                gene_sets=gene_sets,
                ranking_key=params.score_keys,
                permutation_num=params.n_permutations,
                min_size=params.min_genes,
                max_size=params.max_genes,
                context=context
            )
        elif params.method == "ora":
            result_dict = await perform_ora(
                adata=adata,
                gene_sets=gene_sets,
                pvalue_threshold=params.pvalue_cutoff,
                min_size=params.min_genes,
                max_size=params.max_genes,
                context=context
            )
        elif params.method == "ssgsea":
            result_dict = await perform_ssgsea(
                adata=adata,
                gene_sets=gene_sets,
                min_size=params.min_genes,
                max_size=params.max_genes,
                context=context
            )
        elif params.method == "enrichr":
            # For Enrichr, we need a gene list
            if hasattr(params, 'query_genes') and params.query_genes:
                gene_list = params.query_genes
            else:
                # Use highly variable genes or DEGs
                if 'highly_variable' in adata.var:
                    gene_list = adata.var_names[adata.var['highly_variable']].tolist()[:500]
                else:
                    # Use top variable genes
                    var_scores = np.array(adata.X.var(axis=0)).flatten()
                    top_indices = np.argsort(var_scores)[-500:]
                    gene_list = adata.var_names[top_indices].tolist()
            
            result_dict = await perform_enrichr(
                gene_list=gene_list,
                gene_sets=params.gene_set_database,
                organism=adata.uns.get('species', 'human'),
                context=context
            )
        else:
            raise ValueError(f"Unknown enrichment method: {params.method}")

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Create EnrichmentResult object
    result = EnrichmentResult(
        method=params.method,
        n_gene_sets=result_dict.get('n_gene_sets', 0),
        n_significant=result_dict.get('n_significant', 0),
        enrichment_scores=result_dict.get('enrichment_scores', {}),
        pvalues=result_dict.get('pvalues', {}),
        adjusted_pvalues=result_dict.get('adjusted_pvalues', {}),
        gene_set_statistics=result_dict.get('gene_set_statistics', {}),
        spatial_metrics=result_dict.get('spatial_metrics'),
        spatial_scores_key=result_dict.get('spatial_scores_key'),
        gene_sets_used=result_dict.get('gene_sets_used', {}),
        genes_found=result_dict.get('genes_found', {}),
        top_gene_sets=result_dict.get('top_gene_sets', []),
        top_depleted_sets=result_dict.get('top_depleted_sets', []),
        visualization=result_dict.get('visualization'),
        spatial_visualization=result_dict.get('spatial_visualization'),
        parameters_used=params.model_dump(),
        computation_time=time.time() - start_time
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
    context: Context = None
) -> SpatialVariableGenesResult:
    """Identify spatially variable genes

    Args:
        data_id: Dataset ID
        params: Spatial variable gene parameters

    Returns:
        Spatial variable genes result

    Notes:
        Available methods:
        - squidpy: Moran's I and Geary's C
        - sepal: SEPAL neural network
        - somde: SOMDE Gaussian process
        - spark: SPARK statistical model
        - spatialde: SpatialDE Gaussian process
        - hotspot: Hotspot local autocorrelation
    """
    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Call spatial genes function
    result = await identify_spatial_genes(data_id, data_store, params, context)

    # Update dataset in data manager
    data_manager.data_store[data_id] = data_store[data_id]

    # Save spatial genes result
    await data_manager.save_result(data_id, "spatial_genes", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(data_id, "spatial_genes", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def register_spatial_data(
    source_id: str,
    target_id: str,
    method: str = "paste",
    landmarks: Optional[List[Dict[str, Any]]] = None,
    context: Context = None
) -> Dict[str, Any]:
    """Register/align spatial transcriptomics data across sections

    Args:
        source_id: Source dataset ID
        target_id: Target dataset ID to align to
        method: Registration method (paste, manual)
        landmarks: Manual landmarks for alignment (if method='manual')

    Returns:
        Registration result with transformation matrix
    """
    # Import registration function
    from .tools.spatial_registration import register_spatial_slices

    # Validate datasets
    validate_dataset(source_id)
    validate_dataset(target_id)

    # Get datasets from data manager
    source_info = await data_manager.get_dataset(source_id)
    target_info = await data_manager.get_dataset(target_id)
    data_store = {
        source_id: source_info,
        target_id: target_info
    }

    # Call registration function
    result = await register_spatial_slices(
        source_id, target_id, data_store, method, landmarks, context
    )

    # Update datasets in data manager
    data_manager.data_store[source_id] = data_store[source_id]
    data_manager.data_store[target_id] = data_store[target_id]

    # Save registration result
    await data_manager.save_result(source_id, "registration", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(source_id, "registration", result)

    return result


@mcp.tool()
@mcp_tool_error_handler()
async def calculate_spatial_statistics(
    data_id: str,
    feature: str,
    statistic: str = "morans_i",
    n_neighbors: int = 6,
    context: Context = None
) -> Dict[str, Any]:
    """Calculate spatial statistics for features

    Args:
        data_id: Dataset ID
        feature: Feature/gene to analyze
        statistic: Type of statistic (morans_i, gearys_c, local_morans)
        n_neighbors: Number of neighbors for spatial graph

    Returns:
        Spatial statistics result
    """
    # Import spatial statistics function
    from .tools.spatial_statistics import calculate_spatial_stats

    # Validate dataset
    validate_dataset(data_id)

    # Get dataset from data manager
    dataset_info = await data_manager.get_dataset(data_id)
    data_store = {data_id: dataset_info}

    # Call spatial statistics function
    result = await calculate_spatial_stats(
        data_id, data_store, feature, statistic, n_neighbors, context
    )

    # Save statistics result
    await data_manager.save_result(data_id, f"spatial_stats_{feature}", result)

    # Create result resource
    await adapter.resource_manager.create_result_resource(
        data_id, f"statistics_{feature}", result
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
        open_world_hint=True
    ),
    "preprocess_data": MCPToolMetadata(
        name="preprocess_data",
        title="Preprocess Spatial Data",
        description="Preprocess and normalize spatial data",
        read_only_hint=False,
        idempotent_hint=False
    ),
    "visualize_data": MCPToolMetadata(
        name="visualize_data",
        title="Visualize Spatial Data",
        description="Create visualizations of spatial data",
        read_only_hint=True,
        idempotent_hint=True
    ),
    "annotate_cells": MCPToolMetadata(
        name="annotate_cells",
        title="Annotate Cell Types",
        description="Identify cell types in spatial data",
        read_only_hint=False,
        idempotent_hint=False,
        open_world_hint=True
    ),
    "analyze_spatial_data": MCPToolMetadata(
        name="analyze_spatial_data",
        title="Spatial Pattern Analysis",
        description="Analyze spatial patterns and correlations",
        read_only_hint=False,
        idempotent_hint=True
    ),
    "find_markers": MCPToolMetadata(
        name="find_markers",
        title="Find Marker Genes",
        description="Identify differentially expressed genes",
        read_only_hint=True,
        idempotent_hint=True
    ),
    "analyze_velocity_data": MCPToolMetadata(
        name="analyze_velocity_data",
        title="RNA Velocity Analysis",
        description="Analyze RNA velocity dynamics",
        read_only_hint=False,
        idempotent_hint=False
    ),
    "analyze_trajectory_data": MCPToolMetadata(
        name="analyze_trajectory_data",
        title="Trajectory Analysis",
        description="Infer cellular trajectories",
        read_only_hint=False,
        idempotent_hint=False
    ),
    "integrate_samples": MCPToolMetadata(
        name="integrate_samples",
        title="Integrate Multiple Samples",
        description="Integrate multiple spatial datasets",
        read_only_hint=False,
        idempotent_hint=False
    ),
    "deconvolve_data": MCPToolMetadata(
        name="deconvolve_data",
        title="Spatial Deconvolution",
        description="Deconvolve spatial spots into cell types",
        read_only_hint=False,
        idempotent_hint=False,
        open_world_hint=True
    ),
    "identify_spatial_domains": MCPToolMetadata(
        name="identify_spatial_domains",
        title="Identify Spatial Domains",
        description="Find spatial domains and niches",
        read_only_hint=False,
        idempotent_hint=False
    ),
    "analyze_cell_communication": MCPToolMetadata(
        name="analyze_cell_communication",
        title="Cell Communication Analysis",
        description="Analyze cell-cell communication",
        read_only_hint=False,
        idempotent_hint=True,
        open_world_hint=True
    ),
    "analyze_enrichment": MCPToolMetadata(
        name="analyze_enrichment",
        title="Gene Set Enrichment Analysis",
        description="Perform gene set enrichment analysis",
        read_only_hint=False,
        idempotent_hint=True
    ),
    "find_spatial_genes": MCPToolMetadata(
        name="find_spatial_genes",
        title="Find Spatial Variable Genes",
        description="Identify spatially variable genes",
        read_only_hint=False,
        idempotent_hint=False
    ),
    "register_spatial_data": MCPToolMetadata(
        name="register_spatial_data",
        title="Register Spatial Sections",
        description="Align spatial sections",
        read_only_hint=False,
        idempotent_hint=False
    ),
    "calculate_spatial_statistics": MCPToolMetadata(
        name="calculate_spatial_statistics",
        title="Calculate Spatial Statistics",
        description="Calculate spatial statistics for features",
        read_only_hint=True,
        idempotent_hint=True
    )
}

# Update adapter with tool metadata
for name, metadata in tool_metadata.items():
    adapter._tool_metadata[name] = metadata


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