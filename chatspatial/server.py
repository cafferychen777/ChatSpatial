"""
Main server implementation for ChatSpatial using the Spatial MCP Adapter.
"""

import os
import sys
import warnings
from typing import Any, Literal, Optional, cast

# Suppress warnings to speed up startup
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# CRITICAL: Disable progress bars to prevent stdout pollution
# This protects against accidental stdout usage if server is imported directly
os.environ["TQDM_DISABLE"] = "1"

# Suppress scanpy/squidpy verbosity and enable multi-threading
try:
    import scanpy as sc

    sc.settings.verbosity = 0
    sc.settings.n_jobs = -1  # Use all CPU cores for parallel computation
except ImportError:
    pass

from mcp.server.fastmcp import Context  # noqa: E402

from .models.analysis import AnnotationResult  # noqa: E402
from .models.analysis import CellCommunicationResult  # noqa: E402
from .models.analysis import CNVResult  # noqa: E402
from .models.analysis import ConditionComparisonResult  # noqa: E402
from .models.analysis import DeconvolutionResult  # noqa: E402
from .models.analysis import DifferentialExpressionResult  # noqa: E402
from .models.analysis import EnrichmentResult  # noqa: E402
from .models.analysis import IntegrationResult  # noqa: E402
from .models.analysis import PreprocessingResult  # noqa: E402
from .models.analysis import RNAVelocityResult  # noqa: E402
from .models.analysis import SpatialDomainResult  # noqa: E402
from .models.analysis import SpatialStatisticsResult  # noqa: E402
from .models.analysis import SpatialVariableGenesResult  # noqa: E402
from .models.analysis import TrajectoryResult  # noqa: E402
from .models.data import AnnotationParameters  # noqa: E402
from .models.data import CellCommunicationParameters  # noqa: E402
from .models.data import CNVParameters  # noqa: E402
from .models.data import ColumnInfo  # noqa: E402
from .models.data import ConditionComparisonParameters  # noqa: E402
from .models.data import DeconvolutionParameters  # noqa: E402
from .models.data import DifferentialExpressionParameters  # noqa: E402
from .models.data import EnrichmentParameters  # noqa: E402
from .models.data import IntegrationParameters  # noqa: E402
from .models.data import PreprocessingParameters  # noqa: E402
from .models.data import RNAVelocityParameters  # noqa: E402
from .models.data import SpatialDataset  # noqa: E402
from .models.data import SpatialDomainParameters  # noqa: E402
from .models.data import SpatialStatisticsParameters  # noqa: E402
from .models.data import SpatialVariableGenesParameters  # noqa: E402
from .models.data import TrajectoryParameters  # noqa: E402
from .models.data import VisualizationParameters  # noqa: E402
from .spatial_mcp_adapter import ToolContext  # noqa: E402
from .spatial_mcp_adapter import create_spatial_mcp_server  # noqa: E402
from .spatial_mcp_adapter import get_tool_annotations  # noqa: E402
from .utils.exceptions import DataNotFoundError  # noqa: E402
from .utils.mcp_utils import mcp_tool_error_handler  # noqa: E402

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
    if not data_manager.dataset_exists(data_id):
        raise DataNotFoundError(f"Dataset {data_id} not found")


@mcp.tool(annotations=get_tool_annotations("load_data"))
@mcp_tool_error_handler()
async def load_data(
    data_path: str,
    data_type: str,
    name: Optional[str] = None,
    context: Optional[Context] = None,
) -> SpatialDataset:
    """Load spatial transcriptomics data with comprehensive metadata profile.

    Args:
        data_path: Path to data file or directory
        data_type: 'visium', 'xenium', 'slide_seq', 'merfish', 'seqfish', or 'generic'
        name: Optional dataset name

    Returns:
        SpatialDataset with cell/gene counts and metadata profiles
    """
    # Create ToolContext for consistent logging
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    await ctx.info(f"Loading data from {data_path} (type: {data_type})")

    # Load data using data manager
    data_id = await data_manager.load_dataset(data_path, data_type, name)
    dataset_info = await data_manager.get_dataset(data_id)

    await ctx.info(
        f"Successfully loaded {dataset_info['type']} data with "
        f"{dataset_info['n_cells']} cells and {dataset_info['n_genes']} genes"
    )

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


@mcp.tool(annotations=get_tool_annotations("preprocess_data"))
@mcp_tool_error_handler()
async def preprocess_data(
    data_id: str,
    params: PreprocessingParameters = PreprocessingParameters(),
    context: Optional[Context] = None,
) -> PreprocessingResult:
    """Preprocess spatial transcriptomics data.

    Args:
        data_id: Dataset ID
        params: Preprocessing parameters

    Returns:
        PreprocessingResult with HVGs, PCA, clustering, and spatial neighbors
    """
    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import (avoid name conflict with MCP tool)
    from .tools.preprocessing import preprocess_data as preprocess_func

    # Call preprocessing function
    result = await preprocess_func(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save preprocessing result
    await data_manager.save_result(data_id, "preprocessing", result)

    return result


@mcp.tool(annotations=get_tool_annotations("compute_embeddings"))
@mcp_tool_error_handler()
async def compute_embeddings(
    data_id: str,
    compute_pca: bool = True,
    compute_neighbors: bool = True,
    compute_umap: bool = True,
    compute_clustering: bool = True,
    compute_diffmap: bool = False,
    compute_spatial_neighbors: bool = True,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    clustering_resolution: float = 1.0,
    clustering_method: str = "leiden",
    force: bool = False,
    context: Optional[Context] = None,
) -> dict[str, Any]:
    """Compute dimensionality reduction, clustering, and neighbor graphs.

    Args:
        data_id: Dataset ID
        compute_*: Boolean flags for each computation type
        force: Force recomputation even if results exist

    Returns:
        Summary of computed embeddings
    """
    # Validate dataset
    validate_dataset(data_id)

    # Lazy import
    from .tools.embeddings import EmbeddingParameters
    from .tools.embeddings import compute_embeddings as compute_embeddings_func

    # Create parameters - cast clustering_method to Literal type
    clustering_method_literal = cast(Literal["leiden", "louvain"], clustering_method)
    params = EmbeddingParameters(
        compute_pca=compute_pca,
        compute_neighbors=compute_neighbors,
        compute_umap=compute_umap,
        compute_clustering=compute_clustering,
        compute_diffmap=compute_diffmap,
        compute_spatial_neighbors=compute_spatial_neighbors,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        clustering_resolution=clustering_resolution,
        clustering_method=clustering_method_literal,
        force=force,
    )

    # Create ToolContext
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Call function
    result = await compute_embeddings_func(data_id, ctx, params)

    return result.model_dump()


@mcp.tool(annotations=get_tool_annotations("visualize_data"))
@mcp_tool_error_handler()
async def visualize_data(
    data_id: str,
    params: VisualizationParameters = VisualizationParameters(),
    context: Optional[Context] = None,
) -> str:
    """Visualize spatial transcriptomics data.

    Args:
        data_id: Dataset ID
        params: Visualization parameters

    Plot types (11 unified types):
        - feature: Spatial/UMAP feature visualization (use basis='spatial'|'umap')
        - expression: Aggregated expression (subtype='heatmap'|'violin'|'dotplot'|'correlation')
        - deconvolution: Cell type proportions (subtype='spatial_multi'|'pie'|'dominant'|'diversity'|'umap'|'imputation')
        - communication: Cell-cell communication (subtype='heatmap'|'chord'|'dotplot'|'tileplot'|'circle_plot')
        - interaction: Spatial ligand-receptor pairs
        - trajectory: Pseudotime and fate analysis (subtype='pseudotime'|'circular'|'fate_map'|'gene_trends'|'fate_heatmap'|'palantir')
        - velocity: RNA velocity visualization (subtype='stream'|'phase'|'proportions'|'heatmap'|'paga')
        - statistics: Spatial statistics (subtype='neighborhood'|'co_occurrence'|'ripley'|'moran'|'centrality'|'getis_ord')
        - enrichment: Pathway enrichment (subtype='barplot'|'dotplot')
        - cnv: Copy number variation (subtype='heatmap'|'spatial')
        - integration: Batch integration quality (subtype='batch'|'cluster'|'highlight')

    Export options (in params):
        - output_path: Custom save path (default: ./visualizations/)
        - output_format: png, pdf, svg, eps, tiff (default: png)
        - dpi: Resolution (default: 300)

    Returns:
        Path to saved visualization file
    """
    from .tools.visualization import visualize_data as visualize_func

    validate_dataset(data_id)

    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    result = await visualize_func(data_id, ctx, params)

    if result:
        return result
    else:
        return "Visualization generation failed, please check the data and parameter settings."


@mcp.tool(annotations=get_tool_annotations("annotate_cell_types"))
@mcp_tool_error_handler()
async def annotate_cell_types(
    data_id: str,
    params: AnnotationParameters = AnnotationParameters(),
    context: Optional[Context] = None,
) -> AnnotationResult:
    """Annotate cell types in spatial transcriptomics data.

    Args:
        data_id: Dataset ID
        params: Annotation parameters

    Key requirements:
        - Reference methods (tangram, scanvi): reference_data_id must be preprocessed first
        - cell_type_key: Auto-detected if None

    Returns:
        AnnotationResult with cell type assignments and confidence scores
    """
    # Validate dataset
    validate_dataset(data_id)

    # Validate reference data for methods that require it
    if (
        params.method in ["tangram", "scanvi", "singler"]
        and params.reference_data_id
        and not data_manager.dataset_exists(params.reference_data_id)
    ):
        raise DataNotFoundError(
            f"Reference dataset {params.reference_data_id} not found"
        )

    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import annotation tool (avoids slow startup)
    from .tools.annotation import annotate_cell_types

    # Call annotation function with ToolContext
    result = await annotate_cell_types(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save annotation result
    await data_manager.save_result(data_id, "annotation", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool(annotations=get_tool_annotations("analyze_spatial_statistics"))
@mcp_tool_error_handler()
async def analyze_spatial_statistics(
    data_id: str,
    params: SpatialStatisticsParameters = SpatialStatisticsParameters(),
    context: Optional[Context] = None,
) -> SpatialStatisticsResult:
    """Analyze spatial statistics and autocorrelation patterns.

    Args:
        data_id: Dataset ID
        params: Analysis parameters (analysis_type, cluster_key, genes)

    Analysis types:
        - Gene-based: moran, local_moran, geary, getis_ord, bivariate_moran
        - Group-based (requires cluster_key): neighborhood, co_occurrence, ripley
        - Categorical: join_count (binary), local_join_count (multi-category)
        - Network: centrality, network_properties

    Returns:
        SpatialStatisticsResult with statistics and p-values
    """
    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import spatial_statistics (squidpy is slow to import)
    from .tools.spatial_statistics import (
        analyze_spatial_statistics as _analyze_spatial_statistics,
    )

    # Call spatial statistics analysis function with ToolContext
    result = await _analyze_spatial_statistics(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save spatial statistics result
    await data_manager.save_result(data_id, "spatial_statistics", result)

    # Note: Visualization should be created separately using create_visualization tool
    # This maintains clean separation between analysis and visualization

    return result


@mcp.tool(annotations=get_tool_annotations("find_markers"))
@mcp_tool_error_handler()
async def find_markers(
    data_id: str,
    group_key: str,
    group1: Optional[str] = None,
    group2: Optional[str] = None,
    method: str = "wilcoxon",
    n_top_genes: int = 25,  # Number of top differentially expressed genes to return
    pseudocount: float = 1.0,  # Pseudocount for log2 fold change calculation
    min_cells: int = 3,  # Minimum cells per group for statistical testing
    sample_key: Optional[str] = None,  # Sample key for pseudobulk (pydeseq2)
    context: Optional[Context] = None,
) -> DifferentialExpressionResult:
    """Find differentially expressed genes between groups.

    Args:
        data_id: Dataset ID
        group_key: Column name defining groups
        group1: First group (if None, compare each group vs rest)
        group2: Second group (if None, compare group1 vs all others)
        method: Statistical test ('wilcoxon', 't-test', 'pydeseq2')
        sample_key: Required for 'pydeseq2' pseudobulk method

    Returns:
        DifferentialExpressionResult with top marker genes
    """
    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Create params object for unified signature pattern
    params = DifferentialExpressionParameters(
        group_key=group_key,
        group1=group1,
        group2=group2,
        method=method,  # type: ignore[arg-type]
        n_top_genes=n_top_genes,
        pseudocount=pseudocount,
        min_cells=min_cells,
        sample_key=sample_key,
    )

    # Lazy import differential expression tool
    from .tools.differential import differential_expression

    # Call differential expression function with unified (data_id, ctx, params) signature
    result = await differential_expression(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save differential expression result
    await data_manager.save_result(data_id, "differential_expression", result)

    return result


@mcp.tool(annotations=get_tool_annotations("compare_conditions"))
@mcp_tool_error_handler()
async def compare_conditions(
    data_id: str,
    condition_key: str,
    condition1: str,
    condition2: str,
    sample_key: str,
    cell_type_key: Optional[str] = None,
    method: str = "pseudobulk",
    n_top_genes: int = 50,
    min_cells_per_sample: int = 10,
    min_samples_per_condition: int = 2,
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 0.0,
    context: Optional[Context] = None,
) -> ConditionComparisonResult:
    """Compare experimental conditions using pseudobulk differential expression (DESeq2).

    Args:
        data_id: Dataset ID
        condition_key: Column with experimental conditions (e.g., 'treatment')
        condition1: Experimental group
        condition2: Control group
        sample_key: Column identifying biological replicates (e.g., 'patient_id')
        cell_type_key: Optional - if provided, analysis is stratified by cell type

    Returns:
        ConditionComparisonResult with differential expression results
    """
    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Create params object
    params = ConditionComparisonParameters(
        condition_key=condition_key,
        condition1=condition1,
        condition2=condition2,
        sample_key=sample_key,
        cell_type_key=cell_type_key,
        method=method,  # type: ignore[arg-type]
        n_top_genes=n_top_genes,
        min_cells_per_sample=min_cells_per_sample,
        min_samples_per_condition=min_samples_per_condition,
        padj_threshold=padj_threshold,
        log2fc_threshold=log2fc_threshold,
    )

    # Lazy import
    from .tools.condition_comparison import compare_conditions as _compare_conditions

    # Run analysis
    result = await _compare_conditions(data_id, ctx, params)

    # Save result
    await data_manager.save_result(data_id, "condition_comparison", result)

    return result


@mcp.tool(annotations=get_tool_annotations("analyze_cnv"))
@mcp_tool_error_handler()
async def analyze_cnv(
    data_id: str,
    reference_key: str,
    reference_categories: list[str],
    method: str = "infercnvpy",
    window_size: int = 100,
    step: int = 10,
    exclude_chromosomes: Optional[list[str]] = None,
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
    context: Optional[Context] = None,
) -> CNVResult:
    """Analyze copy number variations (CNVs) in spatial transcriptomics data.

    Args:
        data_id: Dataset identifier
        reference_key: Column with cell type labels
        reference_categories: Cell types to use as normal reference (e.g., ['T cells', 'B cells'])
        method: 'infercnvpy' (default, expression-based) or 'numbat' (allele-based, requires R)

    Returns:
        CNVResult with CNV scores and optional clone assignments
    """
    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

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

    # Call CNV inference function with ToolContext
    result = await infer_cnv(data_id=data_id, ctx=ctx, params=params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save CNV result
    await data_manager.save_result(data_id, "cnv_analysis", result)

    return result


@mcp.tool(annotations=get_tool_annotations("analyze_velocity_data"))
@mcp_tool_error_handler()
async def analyze_velocity_data(
    data_id: str,
    params: RNAVelocityParameters = RNAVelocityParameters(),
    context: Optional[Context] = None,
) -> RNAVelocityResult:
    """Analyze RNA velocity to understand cellular dynamics.

    Args:
        data_id: Dataset ID (must have 'spliced' and 'unspliced' layers)
        params: method='scvelo' (modes: deterministic/stochastic/dynamical) or 'velovi'

    Returns:
        RNAVelocityResult with velocity vectors and latent time
    """
    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import velocity analysis tool
    from .tools.velocity import analyze_rna_velocity

    # Call RNA velocity function with ToolContext
    result = await analyze_rna_velocity(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save velocity result
    await data_manager.save_result(data_id, "rna_velocity", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool(annotations=get_tool_annotations("analyze_trajectory_data"))
@mcp_tool_error_handler()
async def analyze_trajectory_data(
    data_id: str,
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Optional[Context] = None,
) -> TrajectoryResult:
    """Infer cellular trajectories and pseudotime.

    Args:
        data_id: Dataset ID
        params: method='cellrank' (requires velocity), 'palantir', or 'dpt'

    Returns:
        TrajectoryResult with pseudotime and fate probabilities
    """
    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import trajectory function
    from .tools.trajectory import analyze_trajectory

    # Call trajectory function
    result = await analyze_trajectory(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save trajectory result
    await data_manager.save_result(data_id, "trajectory", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool(annotations=get_tool_annotations("integrate_samples"))
@mcp_tool_error_handler()
async def integrate_samples(
    data_ids: list[str],
    params: IntegrationParameters = IntegrationParameters(),
    context: Optional[Context] = None,
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples.

    Args:
        data_ids: List of dataset IDs to integrate
        params: method='harmony' (default), 'bbknn', 'scanorama', or 'scvi'

    Returns:
        IntegrationResult with integrated dataset ID
    """
    # Validate all datasets first
    for data_id in data_ids:
        validate_dataset(data_id)

    # Create ToolContext for clean data access
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import to avoid slow startup
    from .tools.integration import integrate_samples as integrate_func

    # Call integration function with ToolContext
    # Note: integrate_func uses ctx.add_dataset() to store the integrated dataset
    result = await integrate_func(data_ids, ctx, params)

    # Save integration result
    integrated_id = result.data_id
    await data_manager.save_result(integrated_id, "integration", result)

    return result


@mcp.tool(annotations=get_tool_annotations("deconvolve_data"))
@mcp_tool_error_handler()
async def deconvolve_data(
    data_id: str,
    params: DeconvolutionParameters,  # No default - LLM must provide parameters
    context: Optional[Context] = None,
) -> DeconvolutionResult:
    """Deconvolve spatial spots to estimate cell type proportions.

    Args:
        data_id: Dataset ID
        params: Required - method, cell_type_key, reference_data_id

    Methods:
        - flashdeconv: Fast sketch-based method (default, recommended)
        - cell2location: Deep learning, accurate but slow (requires scvi-tools)
        - rctd: R-based, modes: doublet (high-res), full (Visium), multi
        - destvi, stereoscope, tangram: Alternative deep learning methods
        - spotlight, card: R-based methods (card supports spatial imputation)

    Returns:
        DeconvolutionResult with cell type proportions per spot
    """
    # Validate dataset
    validate_dataset(data_id)

    # Validate reference data if provided
    if params.reference_data_id and not data_manager.dataset_exists(
        params.reference_data_id
    ):
        raise DataNotFoundError(
            f"Reference dataset {params.reference_data_id} not found"
        )

    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import deconvolution tool
    from .tools.deconvolution import deconvolve_spatial_data

    # Call deconvolution function with ToolContext
    result = await deconvolve_spatial_data(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save deconvolution result
    await data_manager.save_result(data_id, "deconvolution", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool(annotations=get_tool_annotations("identify_spatial_domains"))
@mcp_tool_error_handler()
async def identify_spatial_domains(
    data_id: str,
    params: SpatialDomainParameters = SpatialDomainParameters(),
    context: Optional[Context] = None,
) -> SpatialDomainResult:
    """Identify spatial domains and tissue architecture.

    Args:
        data_id: Dataset ID
        params: method='spagcn' (default, uses histology), 'leiden', 'louvain', 'stagate', 'graphst'

    Returns:
        SpatialDomainResult with domain_key for visualization
    """
    # Validate dataset first
    validate_dataset(data_id)

    # Create ToolContext for clean data access
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import to avoid slow startup
    from .tools.spatial_domains import identify_spatial_domains as identify_domains_func

    # Call spatial domains function with ToolContext
    result = await identify_domains_func(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save spatial domains result
    await data_manager.save_result(data_id, "spatial_domains", result)

    return result


@mcp.tool(annotations=get_tool_annotations("analyze_cell_communication"))
@mcp_tool_error_handler()
async def analyze_cell_communication(
    data_id: str,
    params: CellCommunicationParameters,  # No default - LLM must provide parameters
    context: Optional[Context] = None,
) -> CellCommunicationResult:
    """Analyze cell-cell communication patterns.

    Args:
        data_id: Dataset ID
        params: Required - species, cell_type_key, and method

    Methods:
        - liana: Multi-method consensus (default). Use liana_resource for database selection
        - cellphonedb: Statistical permutation-based analysis
        - cellchat_r: R-based CellChat (requires rpy2)
        - fastccc: Fast C++ implementation (human only)

    Species configuration:
        - human: liana_resource="consensus" (default)
        - mouse: liana_resource="mouseconsensus"

    Returns:
        CellCommunicationResult with significant ligand-receptor interactions
    """
    # Validate dataset first
    validate_dataset(data_id)

    # Create ToolContext for clean data access
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import to avoid slow startup
    from .tools.cell_communication import (
        analyze_cell_communication as analyze_comm_func,
    )

    # Call cell communication function with ToolContext
    result = await analyze_comm_func(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save communication result
    await data_manager.save_result(data_id, "cell_communication", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool(annotations=get_tool_annotations("analyze_enrichment"))
@mcp_tool_error_handler()
async def analyze_enrichment(
    data_id: str,
    params: Optional[EnrichmentParameters] = None,
    context: Optional[Context] = None,
) -> EnrichmentResult:
    """Perform gene set enrichment analysis.

    Args:
        data_id: Dataset ID
        params: Required - species must be specified ('human' or 'mouse')

    Methods:
        - pathway_ora: Over-representation analysis (default)
        - pathway_gsea: Gene Set Enrichment Analysis
        - spatial_enrichmap: Spatial enrichment mapping
        - pathway_enrichr, pathway_ssgsea: Alternative methods

    Databases (gene_set_database):
        - KEGG_Pathways (recommended), Reactome_Pathways, MSigDB_Hallmark
        - GO_Biological_Process (default), GO_Molecular_Function, GO_Cellular_Component

    Returns:
        EnrichmentResult with enriched pathways and statistics
    """
    from .tools.enrichment import analyze_enrichment as analyze_enrichment_func

    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Use default parameters if not provided (species is required by analyze_enrichment_func)
    if params is None:
        raise ValueError(
            "EnrichmentParameters is required. Please specify at least 'species' parameter."
        )

    # Call enrichment analysis (all business logic is in tools/enrichment.py)
    result = await analyze_enrichment_func(data_id, ctx, params)

    # Save result
    await data_manager.save_result(data_id, "enrichment", result)

    return result


@mcp.tool(annotations=get_tool_annotations("find_spatial_genes"))
@mcp_tool_error_handler()
async def find_spatial_genes(
    data_id: str,
    params: SpatialVariableGenesParameters = SpatialVariableGenesParameters(),
    context: Optional[Context] = None,
) -> SpatialVariableGenesResult:
    """Identify spatially variable genes.

    Args:
        data_id: Dataset ID
        params: method='sparkx' (default, fast) or 'spatialde' (Gaussian process)

    Returns:
        SpatialVariableGenesResult with ranked genes and statistics
    """
    # Validate dataset
    validate_dataset(data_id)

    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import spatial genes tool
    from .tools.spatial_genes import identify_spatial_genes

    # Call spatial genes function with ToolContext
    result = await identify_spatial_genes(data_id, ctx, params)

    # Note: No writeback needed - adata modifications are in-place on the same object

    # Save spatial genes result
    await data_manager.save_result(data_id, "spatial_genes", result)

    # Visualization should be done separately via visualization tools

    return result


@mcp.tool(annotations=get_tool_annotations("register_spatial_data"))
@mcp_tool_error_handler()
async def register_spatial_data(
    source_id: str,
    target_id: str,
    method: str = "paste",
    context: Optional[Context] = None,
) -> dict[str, Any]:
    """Register/align spatial transcriptomics data across sections

    Args:
        source_id: Source dataset ID
        target_id: Target dataset ID to align to
        method: Registration method (paste, stalign)

    Returns:
        Registration result with transformation matrix
    """
    # Validate datasets first
    validate_dataset(source_id)
    validate_dataset(target_id)

    # Create ToolContext for unified data access
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import to avoid slow startup
    from .tools.spatial_registration import register_spatial_slices_mcp

    # Call registration function using ToolContext
    # Note: registration modifies adata in-place, changes reflected via reference
    result = await register_spatial_slices_mcp(source_id, target_id, ctx, method)

    # Save registration result
    await data_manager.save_result(source_id, "registration", result)

    return result


# ============== Data Export/Reload Tools ==============


@mcp.tool(annotations=get_tool_annotations("export_data"))
@mcp_tool_error_handler()
async def export_data(
    data_id: str,
    path: Optional[str] = None,
    context: Optional[Context] = None,
) -> str:
    """Export dataset to disk for external script access.

    Args:
        data_id: Dataset ID to export
        path: Custom path (default: ~/.chatspatial/active/{data_id}.h5ad)

    Returns:
        Absolute path where data was exported
    """
    from pathlib import Path as PathLib

    from .utils.persistence import export_adata

    # Validate dataset exists
    validate_dataset(data_id)

    if context:
        await context.info(f"Exporting dataset '{data_id}'...")

    # Get dataset info
    dataset_info = await data_manager.get_dataset(data_id)
    adata = dataset_info["adata"]

    try:
        export_path = export_adata(data_id, adata, PathLib(path) if path else None)
        absolute_path = export_path.resolve()

        if context:
            await context.info(f"Dataset exported to: {absolute_path}")

        return f"Dataset '{data_id}' exported to: {absolute_path}"

    except Exception as e:
        error_msg = f"Failed to export dataset: {e}"
        if context:
            await context.error(error_msg)
        raise


@mcp.tool(annotations=get_tool_annotations("reload_data"))
@mcp_tool_error_handler()
async def reload_data(
    data_id: str,
    path: Optional[str] = None,
    context: Optional[Context] = None,
) -> str:
    """Reload dataset from disk after external script modifications.

    Args:
        data_id: Dataset ID to reload (must exist in MCP memory)
        path: Custom path (default: ~/.chatspatial/active/{data_id}.h5ad)

    Returns:
        Summary of reloaded dataset
    """
    from pathlib import Path as PathLib

    from .utils.persistence import load_adata_from_active

    # Validate dataset exists in memory
    validate_dataset(data_id)

    if context:
        await context.info(f"Reloading dataset '{data_id}'...")

    try:
        adata = load_adata_from_active(data_id, PathLib(path) if path else None)

        # Update the in-memory dataset
        await data_manager.update_adata(data_id, adata)

        if context:
            await context.info(f"Dataset '{data_id}' reloaded successfully")

        return (
            f"Dataset '{data_id}' reloaded: "
            f"{adata.n_obs} cells × {adata.n_vars} genes"
        )

    except FileNotFoundError as e:
        error_msg = str(e)
        if context:
            await context.error(error_msg)
        raise
    except Exception as e:
        error_msg = f"Failed to reload dataset: {e}"
        if context:
            await context.error(error_msg)
        raise


# CLI entry point is in __main__.py (single source of truth)
# Use: python -m chatspatial server [options]
