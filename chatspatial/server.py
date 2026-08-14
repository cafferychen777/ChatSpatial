"""
Main server implementation for ChatSpatial using the Spatial MCP Adapter.
"""

from typing import Annotated, Optional, TypeVar

from mcp.server.mcpserver import Context
from mcp.types import ToolAnnotations
from pydantic import Field

# Initialize runtime configuration at the application boundary. Importing the
# configuration module alone remains side-effect free for library consumers.
# Runtime initialization configures:
# - Environment variables (TQDM_DISABLE, DASK_*)
# - Warning filters
# - Scanpy settings
from . import __version__, config

config.init_runtime()

from .models.analysis import AnnotationResult  # noqa: E402
from .models.analysis import CellCommunicationResult  # noqa: E402
from .models.analysis import CNVResult  # noqa: E402
from .models.analysis import ConditionComparisonResult  # noqa: E402
from .models.analysis import DeconvolutionResult  # noqa: E402
from .models.analysis import DifferentialExpressionResult  # noqa: E402
from .models.analysis import EmbeddingResult  # noqa: E402
from .models.analysis import EnrichmentResult  # noqa: E402
from .models.analysis import IntegrationResult  # noqa: E402
from .models.analysis import PreprocessingResult  # noqa: E402
from .models.analysis import RNAVelocityResult  # noqa: E402
from .models.analysis import SpatialDomainResult  # noqa: E402
from .models.analysis import SpatialRegistrationResult  # noqa: E402
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
from .models.data import RegistrationParameters  # noqa: E402
from .models.data import RNAVelocityParameters  # noqa: E402
from .models.data import SpatialDataset  # noqa: E402
from .models.data import SpatialDomainParameters  # noqa: E402
from .models.data import SpatialPlatform  # noqa: E402
from .models.data import SpatialStatisticsParameters  # noqa: E402
from .models.data import SpatialVariableGenesParameters  # noqa: E402
from .models.data import TrajectoryParameters  # noqa: E402
from .models.data import VisualizationParameters  # noqa: E402
from .spatial_mcp_adapter import ToolContext  # noqa: E402
from .spatial_mcp_adapter import create_spatial_mcp_server  # noqa: E402
from .spatial_mcp_adapter import serialize_dataset_access  # noqa: E402
from .tools.embeddings import EmbeddingParameters  # noqa: E402
from .utils.mcp_utils import (  # noqa: E402
    mcp_tool_error_handler,
    suppress_output_async,
)

# Create MCP server and adapter
mcp, adapter = create_spatial_mcp_server(
    "ChatSpatial",
    server_version=__version__,
)

# Get data manager from adapter
data_manager = adapter.data_manager

P = TypeVar("P")
R = TypeVar("R")


def _resolve_params(params: Optional[P], default_factory: type[P]) -> P:
    """Resolve optional params to a concrete model instance."""
    return params if params is not None else default_factory()


def _finalize_result(ctx: ToolContext, result: R) -> R:
    """Attach request-scoped warnings to a completed public result."""
    return ctx.finalize(result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=False,
        idempotent_hint=False,
        open_world_hint=True,
    )
)
@mcp_tool_error_handler()
async def load_data(
    data_path: str,
    data_type: SpatialPlatform,
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

    # Stage the dataset first. It becomes visible only after the complete public
    # result has passed validation below.
    data_id, dataset_info = await data_manager.stage_dataset(data_path, data_type, name)

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

    result = SpatialDataset(
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
    await ctx.info(
        f"Successfully loaded {dataset_info['type']} data with "
        f"{dataset_info['n_cells']} cells and {dataset_info['n_genes']} genes"
    )
    await data_manager.publish_dataset(data_id, dataset_info)
    return result


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def preprocess_data(
    data_id: str,
    params: Optional[PreprocessingParameters] = None,
    context: Optional[Context] = None,
) -> PreprocessingResult:
    """Run QC, filtering, normalization, and highly variable gene selection.

    This tool does not compute PCA, UMAP, clustering, or neighbor graphs. Run
    compute_embeddings() afterward when downstream tools require those artifacts.

    Args:
        data_id: Dataset ID
        params: Preprocessing parameters (all have sensible defaults)
    """
    # Create ToolContext
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import (avoid name conflict with MCP tool)
    from .tools.preprocessing import preprocess_data as preprocess_func

    # Resolve optional wrapper input to concrete params for tool-level contract
    resolved_params = _resolve_params(params, PreprocessingParameters)

    # Call preprocessing function
    result = await preprocess_func(data_id, ctx, resolved_params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def compute_embeddings(
    data_id: str,
    params: Optional[EmbeddingParameters] = None,
    context: Optional[Context] = None,
) -> EmbeddingResult:
    """Compute dimensionality reduction (PCA, UMAP), clustering, and neighbor graphs.

    Args:
        data_id: Dataset ID
        params: Embedding parameters (PCA, UMAP, clustering, etc.)
    """
    from .tools.embeddings import compute_embeddings as compute_embeddings_func

    resolved = _resolve_params(params, EmbeddingParameters)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)
    async with suppress_output_async():
        result = await compute_embeddings_func(data_id, ctx, resolved)
    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=True,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def visualize_data(
    data_id: str,
    params: Optional[VisualizationParameters] = None,
    context: Optional[Context] = None,
) -> str:
    """Visualize spatial transcriptomics data. Set plot_type and subtype in params; see VisualizationParameters schema for all options.

    Args:
        data_id: Dataset ID
        params: Visualization parameters (plot_type, subtype, genes, output_format, dpi, etc.)
    """
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    resolved_params = _resolve_params(params, VisualizationParameters)
    async with suppress_output_async():
        from .tools.visualization.main import visualize_data as visualize_func

        result = await visualize_func(data_id, ctx, resolved_params)

    if not result:
        result = (
            "Visualization generation failed, please check the data and parameter "
            "settings."
        )
    return ctx.finalize(result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=True,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id", "params.reference_data_id")
async def annotate_cell_types(
    data_id: str,
    params: Optional[AnnotationParameters] = None,
    context: Optional[Context] = None,
) -> AnnotationResult:
    """Annotate cell types in spatial transcriptomics data.

    Args:
        data_id: Dataset ID
        params: Annotation parameters (method, reference_data_id, cell_type_key, etc.)

    Note: Reference methods (tangram, scanvi) require reference_data_id to be preprocessed first.
    """
    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import annotation tool (avoids slow startup)
    from .tools.annotation import annotate_cell_types

    resolved_params = _resolve_params(params, AnnotationParameters)

    # Call annotation function with ToolContext
    async with suppress_output_async():
        result = await annotate_cell_types(data_id, ctx, resolved_params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def analyze_spatial_statistics(
    data_id: str,
    params: Optional[SpatialStatisticsParameters] = None,
    context: Optional[Context] = None,
) -> SpatialStatisticsResult:
    """Analyze spatial statistics and autocorrelation patterns.

    Args:
        data_id: Dataset ID
        params: Analysis parameters (analysis_type, cluster_key, genes). See SpatialStatisticsParameters for all types.
    """
    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import spatial_statistics (squidpy is slow to import)
    from .tools.spatial_statistics import (
        analyze_spatial_statistics as _analyze_spatial_statistics,
    )

    resolved_params = _resolve_params(params, SpatialStatisticsParameters)

    # Call spatial statistics analysis function with ToolContext
    result = await _analyze_spatial_statistics(data_id, ctx, resolved_params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def find_markers(
    data_id: str,
    params: DifferentialExpressionParameters,
    context: Optional[Context] = None,
) -> DifferentialExpressionResult:
    """Find differentially expressed genes between groups.

    Args:
        data_id: Dataset ID
        params: Required - group_key and optional method, group1/group2, n_top_genes, etc.
    """
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    from .tools.differential import differential_expression

    result = await differential_expression(data_id, ctx, params)
    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def compare_conditions(
    data_id: str,
    params: ConditionComparisonParameters,
    context: Optional[Context] = None,
) -> ConditionComparisonResult:
    """Compare experimental conditions using pseudobulk differential expression (DESeq2).

    Args:
        data_id: Dataset ID
        params: Required - condition_key, condition1, condition2, sample_key, etc.
    """
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    from .tools.condition_comparison import compare_conditions as _compare_conditions

    result = await _compare_conditions(data_id, ctx, params)
    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def analyze_cnv(
    data_id: str,
    params: CNVParameters,
    context: Optional[Context] = None,
) -> CNVResult:
    """Analyze copy number variations (CNVs) in spatial transcriptomics data.

    Args:
        data_id: Dataset identifier
        params: Required - reference_key, reference_categories, and optional method/thresholds.
    """
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    from .tools.cnv_analysis import infer_cnv

    result = await infer_cnv(data_id=data_id, ctx=ctx, params=params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def analyze_velocity_data(
    data_id: str,
    params: Optional[RNAVelocityParameters] = None,
    context: Optional[Context] = None,
) -> RNAVelocityResult:
    """Analyze RNA velocity to understand cellular dynamics. Requires 'spliced' and 'unspliced' layers.

    Args:
        data_id: Dataset ID
        params: Velocity parameters (method, scvelo_mode, etc.)
    """
    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import velocity analysis tool
    from .tools.velocity import analyze_rna_velocity

    resolved_params = _resolve_params(params, RNAVelocityParameters)

    # Call RNA velocity function with ToolContext
    result = await analyze_rna_velocity(data_id, ctx, resolved_params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def analyze_trajectory_data(
    data_id: str,
    params: Optional[TrajectoryParameters] = None,
    context: Optional[Context] = None,
) -> TrajectoryResult:
    """Infer cellular trajectories and pseudotime ordering.

    Args:
        data_id: Dataset ID
        params: Trajectory parameters (method, root_cell, spatial_weight, etc.)
    """
    # Create ToolContext
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import trajectory function
    from .tools.trajectory import analyze_trajectory

    resolved_params = _resolve_params(params, TrajectoryParameters)

    # Call trajectory function
    result = await analyze_trajectory(data_id, ctx, resolved_params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_ids")
async def integrate_samples(
    data_ids: Annotated[
        list[str],
        Field(min_length=2, description="Two or more distinct dataset IDs."),
    ],
    params: Optional[IntegrationParameters] = None,
    context: Optional[Context] = None,
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples into a unified dataset.

    Args:
        data_ids: List of dataset IDs to integrate
        params: Integration parameters (method, batch_key, n_pcs, etc.)
    """
    # Create ToolContext for clean data access
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import to avoid slow startup
    from .tools.integration import integrate_samples as integrate_func

    resolved_params = _resolve_params(params, IntegrationParameters)

    # Call integration function with ToolContext
    # Note: integrate_func uses ctx.add_dataset() to store the integrated dataset
    result = await integrate_func(data_ids, ctx, resolved_params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=True,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id", "params.reference_data_id")
async def deconvolve_data(
    data_id: str,
    params: DeconvolutionParameters,  # No default - LLM must provide parameters
    context: Optional[Context] = None,
) -> DeconvolutionResult:
    """Deconvolve spatial spots to estimate cell type proportions.

    Args:
        data_id: Dataset ID
        params: Required - method, cell_type_key, reference_data_id. See DeconvolutionParameters for all methods and options.
    """
    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import deconvolution tool
    from .tools.deconvolution import deconvolve_spatial_data

    # Call deconvolution function with ToolContext
    result = await deconvolve_spatial_data(data_id, ctx, params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def identify_spatial_domains(
    data_id: str,
    params: Optional[SpatialDomainParameters] = None,
    context: Optional[Context] = None,
) -> SpatialDomainResult:
    """Identify spatial domains and tissue architecture.

    Args:
        data_id: Dataset ID
        params: Spatial domain parameters (method, n_domains, resolution, etc.)
    """
    # Create ToolContext for clean data access
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import to avoid slow startup
    from .tools.spatial_domains import identify_spatial_domains as identify_domains_func

    resolved_params = _resolve_params(params, SpatialDomainParameters)

    # Call spatial domains function with ToolContext
    result = await identify_domains_func(data_id, ctx, resolved_params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=True,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def analyze_cell_communication(
    data_id: str,
    params: CellCommunicationParameters,  # No default - LLM must provide parameters
    context: Optional[Context] = None,
) -> CellCommunicationResult:
    """Analyze cell-cell communication and ligand-receptor interaction patterns.

    Args:
        data_id: Dataset ID
        params: Required - species, cell_type_key, and method. For mouse with liana, set liana_resource='mouseconsensus'.
    """
    # Create ToolContext for clean data access
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import to avoid slow startup
    from .tools.cell_communication import (
        analyze_cell_communication as analyze_comm_func,
    )

    # Call cell communication function with ToolContext
    async with suppress_output_async():
        result = await analyze_comm_func(data_id, ctx, params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=True,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def analyze_enrichment(
    data_id: str,
    params: EnrichmentParameters,
    context: Optional[Context] = None,
) -> EnrichmentResult:
    """Perform gene set enrichment analysis.

    Args:
        data_id: Dataset ID
        params: Required - species must be specified. See EnrichmentParameters for methods and gene_set_database options.
    """
    from .tools.enrichment import analyze_enrichment as analyze_enrichment_func

    # Create ToolContext
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Call enrichment analysis (all business logic is in tools/enrichment.py)
    result = await analyze_enrichment_func(data_id, ctx, params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
async def find_spatial_genes(
    data_id: str,
    params: Optional[SpatialVariableGenesParameters] = None,
    context: Optional[Context] = None,
) -> SpatialVariableGenesResult:
    """Identify spatially variable genes.

    Args:
        data_id: Dataset ID
        params: Spatial variable gene parameters (method, n_top_genes, etc.)
    """
    # Create ToolContext for clean data access (no redundant dict wrapping)
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import spatial genes tool
    from .tools.spatial_genes import identify_spatial_genes

    resolved_params = _resolve_params(params, SpatialVariableGenesParameters)

    # Call spatial genes function with ToolContext
    result = await identify_spatial_genes(data_id, ctx, resolved_params)

    return _finalize_result(ctx, result)


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=False,
        open_world_hint=False,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "source_id", "target_id")
async def register_spatial_data(
    source_id: str,
    target_id: str,
    params: Optional[RegistrationParameters] = None,
    context: Optional[Context] = None,
) -> SpatialRegistrationResult:
    """Register/align spatial transcriptomics data across sections

    Args:
        source_id: Source dataset ID
        target_id: Target dataset ID to align to
        params: Registration parameters (method, alignment settings, etc.)

    Returns:
        Registration result with method, dataset IDs, spot counts, and registered spatial key
    """
    # Create ToolContext for unified data access
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)

    # Lazy import to avoid slow startup
    from .tools.spatial_registration import register_spatial_slices_mcp

    resolved_params = _resolve_params(params, RegistrationParameters)

    # Call registration function using ToolContext
    # Note: registration modifies adata in-place, changes reflected via reference
    result = await register_spatial_slices_mcp(
        source_id, target_id, ctx, resolved_params
    )

    return _finalize_result(ctx, result)


# ============== Data Export/Reload Tools ==============


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=True,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
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

    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)
    await ctx.info(f"Exporting dataset '{data_id}'...")

    # Get dataset info
    dataset_info = await data_manager.get_dataset(data_id)
    adata = dataset_info["adata"]

    export_path = export_adata(data_id, adata, PathLib(path) if path else None)
    absolute_path = export_path.resolve()
    await ctx.info(f"Dataset exported to: {absolute_path}")
    return ctx.finalize(f"Dataset '{data_id}' exported to: {absolute_path}")


@mcp.tool(
    annotations=ToolAnnotations(
        read_only_hint=False,
        destructive_hint=True,
        idempotent_hint=True,
        open_world_hint=True,
    )
)
@mcp_tool_error_handler()
@serialize_dataset_access(data_manager, "data_id")
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

    from .utils.persistence import (
        check_reload_identity,
        get_active_path,
        load_adata_from_active,
    )

    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)
    await ctx.info(f"Reloading dataset '{data_id}'...")

    load_path = PathLib(path).expanduser() if path else get_active_path(data_id)
    adata = load_adata_from_active(data_id, load_path)
    dataset_info = await data_manager.get_dataset(data_id)

    identity_warning = check_reload_identity(
        dataset_info["adata"], adata, data_id, load_path
    )
    if identity_warning:
        await ctx.warning(identity_warning)

    await data_manager.update_adata(data_id, adata)
    await ctx.info(f"Dataset '{data_id}' reloaded successfully")

    return ctx.finalize(
        f"Dataset '{data_id}' reloaded: {adata.n_obs} cells × {adata.n_vars} genes"
    )


# CLI entry point is in __main__.py (single source of truth)
# Use: python -m chatspatial server [options]
