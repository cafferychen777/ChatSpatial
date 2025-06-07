"""
Main server implementation for ChatSpatial.
"""

from typing import Dict, Any, List, Optional
import warnings

# Suppress warnings to speed up startup
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from mcp.server.fastmcp import FastMCP, Context
from mcp.server.fastmcp.utilities.types import Image
from .utils.image_utils import fig_to_image, create_placeholder_image

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
from .tools.preprocessing import preprocess_data
from .tools.visualization import visualize_data
from .tools.annotation import annotate_cell_types
from .tools.spatial_analysis import analyze_spatial_unified
from .tools.differential import differential_expression
from .tools.trajectory import analyze_rna_velocity, analyze_trajectory
from .tools.integration import integrate_samples
from .tools.deconvolution import deconvolve_spatial_data
from .tools.spatial_domains import identify_spatial_domains

from .tools.cell_communication import analyze_cell_communication
from .tools.gaston_spatial_genes import identify_spatial_variable_genes_gaston
from .utils.data_loader import load_spatial_data

# Create MCP server
mcp = FastMCP("ChatSpatial")

# Store for loaded datasets
data_store: Dict[str, Any] = {}


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
async def preprocess(
    data_id: str,
    params: AnalysisParameters = AnalysisParameters(),
    context: Context = None
) -> PreprocessingResult:
    """Preprocess spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Preprocessing parameters

    Returns:
        Preprocessing result
    """
    # Validate dataset
    validate_dataset(data_id)

    # Call preprocessing function
    result = await preprocess_data(data_id, data_store, params, context)

    return result


@mcp.tool()
async def visualize(
    data_id: str,
    params: VisualizationParameters = VisualizationParameters(),
    context: Context = None
):
    """Visualize spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Visualization parameters

    Returns:
        Visualization image
    """
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
    return await visualize_data(data_id, data_store, params, context)


@mcp.tool()
async def annotate(
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
        For Tangram method, a reference_data_id must be provided pointing to a single-cell dataset.
        The Tangram method maps single-cell data to spatial data to infer cell types in space.
    """
    # Validate dataset
    validate_dataset(data_id)

    # Validate reference data for Tangram method
    if params.method == "tangram" and params.reference_data_id:
        if params.reference_data_id not in data_store:
            raise ValueError(f"Reference dataset {params.reference_data_id} not found")
        if context:
            await context.info(f"Using reference dataset {params.reference_data_id} for Tangram mapping")

    # Call annotation function
    result = await annotate_cell_types(data_id, data_store, params, context)

    # Log results
    if context:
        await context.info(f"Annotation completed with {len(result.cell_types)} cell types identified")
        if params.method == "tangram" and result.tangram_mapping_score is not None:
            await context.info(f"Tangram mapping score: {result.tangram_mapping_score:.4f}")

    return result


@mcp.tool()
async def analyze_spatial(
    data_id: str,
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Context = None
) -> Image:
    """Perform spatial analysis on transcriptomics data

    Args:
        data_id: Dataset ID
        params: Spatial analysis parameters

    Returns:
        Spatial analysis visualization image
    """
    # Validate dataset
    validate_dataset(data_id)

    # Call spatial analysis unified function with return_type="image"
    return await analyze_spatial_unified(data_id, data_store, params, context, return_type="image")


@mcp.tool()
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


@mcp.tool()
async def list_datasets(context: Context = None) -> List[SpatialDataset]:
    """List all loaded datasets

    Returns:
        List of dataset information
    """
    if context:
        await context.info("Listing all datasets")

    datasets = []
    for data_id, info in data_store.items():
        datasets.append(
            SpatialDataset(
                id=data_id,
                name=info["name"],
                data_type=info["type"],
                description=f"Loaded from {info.get('path', 'unknown')}"
            )
        )

    return datasets


# Add RNA velocity analysis tool
@mcp.tool()
async def analyze_velocity(
    data_id: str,
    params: RNAVelocityParameters = RNAVelocityParameters(),
    context: Context = None
) -> Image:
    """Analyze RNA velocity in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: RNA velocity analysis parameters

    Returns:
        RNA velocity visualization image
    """
    # Validate dataset
    validate_dataset(data_id)

    # Call RNA velocity analysis function
    result = await analyze_rna_velocity(data_id, data_store, params, context)

    # Return the visualization image directly
    if result.visualization is None:
        # Create a simple placeholder image if no visualization is available
        return create_placeholder_image("No velocity visualization available")

    return result.visualization


# Add trajectory analysis tool with visualization type parameter
@mcp.tool()
async def analyze_trajectory_visualization(
    data_id: str,
    params: TrajectoryParameters = TrajectoryParameters(),
    visualization_type: str = "pseudotime",  # "pseudotime" or "velocity"
    context: Context = None
) -> Image:
    """Analyze trajectory and visualize results

    This function combines the functionality of analyze_spatial_trajectory and
    visualize_trajectory_velocity into a single function with a visualization_type parameter.

    Args:
        data_id: Dataset ID
        params: Trajectory analysis parameters
        visualization_type: Type of visualization to return ("pseudotime" or "velocity")
        context: MCP context

    Returns:
        Trajectory visualization image (pseudotime or velocity)
    """
    # Validate dataset
    validate_dataset(data_id)

    # Validate visualization type
    if visualization_type not in ["pseudotime", "velocity"]:
        raise ValueError(f"Invalid visualization_type: {visualization_type}. Must be 'pseudotime' or 'velocity'")

    # Call trajectory analysis function
    result = await analyze_trajectory(data_id, data_store, params, context)

    # Return the requested visualization
    if visualization_type == "velocity":
        if result.velocity_visualization is None:
            # Create a simple placeholder image if no visualization is available
            return create_placeholder_image("No velocity stream visualization available")
        return result.velocity_visualization
    else:  # pseudotime
        return result.pseudotime_visualization


# For backward compatibility
@mcp.tool()
async def analyze_spatial_trajectory(
    data_id: str,
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Context = None
) -> Image:
    """Analyze trajectory and cell state transitions in spatial transcriptomics data

    This function is maintained for backward compatibility.
    Consider using analyze_trajectory_visualization with visualization_type="pseudotime" instead.

    Args:
        data_id: Dataset ID
        params: Trajectory analysis parameters

    Returns:
        Pseudotime visualization image
    """
    return await analyze_trajectory_visualization(data_id, params, "pseudotime", context)


# For backward compatibility
@mcp.tool()
async def visualize_trajectory_velocity(
    data_id: str,
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Context = None
) -> Image:
    """Visualize velocity streams from trajectory analysis

    This function is maintained for backward compatibility.
    Consider using analyze_trajectory_visualization with visualization_type="velocity" instead.

    Args:
        data_id: Dataset ID
        params: Trajectory analysis parameters

    Returns:
        Velocity stream visualization image
    """
    return await analyze_trajectory_visualization(data_id, params, "velocity", context)


# Add tool to get spatial analysis statistics
@mcp.tool()
async def get_spatial_analysis_stats(
    data_id: str,
    params: SpatialAnalysisParameters = SpatialAnalysisParameters(),
    context: Context = None
) -> SpatialAnalysisResult:
    """Get statistics from spatial analysis without image

    Args:
        data_id: Dataset ID
        params: Spatial analysis parameters

    Returns:
        Spatial analysis result with statistics but no image
    """
    # Validate dataset
    validate_dataset(data_id)

    # Don't include image to save bandwidth
    params.include_image = False

    # Call spatial analysis unified function with return_type="result"
    return await analyze_spatial_unified(data_id, data_store, params, context, return_type="result")


# Add tool for integrating multiple samples
@mcp.tool()
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
    """
    if context:
        await context.info(f"Integrating {len(data_ids)} samples using {params.method} method")

    # Validate all dataset IDs
    for data_id in data_ids:
        validate_dataset(data_id)

    # Call integration function
    return await integrate_samples(data_ids, data_store, params, context)


# Add tool for spatial deconvolution
@mcp.tool()
async def deconvolve(
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
                await visualize_data(data_id, data_store, vis_params, context)
            except Exception as e:
                if context:
                    await context.warning(f"Failed to automatically visualize deconvolution results: {str(e)}")

        return result
    except Exception as e:
        error_msg = f"Error in deconvolution: {str(e)}"
        if context:
            await context.warning(error_msg)
            await context.info("If you're having issues with parameter formatting, try using this format:")
            await context.info('deconvolve(data_id="data_1", params={"method": "nnls", "reference_data_id": "data_2", "cell_type_key": "CellType"})')
        raise

# Add a helper tool for deconvolution with simpler parameter format
@mcp.tool()
async def deconvolve_data(
    data_id: str,
    method: str = "nnls",
    reference_data_id: str = None,
    cell_type_key: str = "cell_type",
    n_top_genes: int = 2000,
    context: Context = None
) -> DeconvolutionResult:
    """Deconvolve spatial data with simpler parameter format

    This is a helper function that provides a simpler interface to the deconvolve tool.

    Args:
        data_id: Dataset ID
        method: Deconvolution method (nnls, cell2location, spotiphy)
        reference_data_id: Reference dataset ID
        cell_type_key: Key in reference data for cell type information
        n_top_genes: Number of top genes to use

    Returns:
        Deconvolution result with cell type proportions and visualization
    """
    if context:
        await context.info(f"Deconvolving spatial data using {method} method")
        await context.info(f"Reference dataset: {reference_data_id}, Cell type key: {cell_type_key}")

    # Create parameters object
    params = DeconvolutionParameters(
        method=method,
        reference_data_id=reference_data_id,
        cell_type_key=cell_type_key,
        n_top_genes=n_top_genes
    )

    # Call the main deconvolve function
    return await deconvolve(data_id, params, context)


# Add a specialized tool for visualizing UMAP with cell types from deconvolution
@mcp.tool()
async def visualize_umap_with_celltypes(
    data_id: str,
    deconv_key: str = "deconvolution_spotiphy",
    context: Context = None
) -> Image:
    """Visualize UMAP with cell types from deconvolution results

    This specialized tool creates a UMAP visualization colored by the dominant cell type
    from deconvolution results. It's particularly useful for visualizing cell type
    distributions in the UMAP embedding.

    Args:
        data_id: Dataset ID
        deconv_key: Key for deconvolution results in adata.obsm

    Returns:
        UMAP visualization image colored by cell types
    """
    # Validate dataset
    validate_dataset(data_id)

    if context:
        await context.info(f"Creating UMAP visualization with cell types from {deconv_key}")

    # Get the AnnData object
    adata = data_store[data_id]["adata"]

    # Check if deconvolution results exist
    if deconv_key not in adata.obsm:
        raise ValueError(f"Deconvolution results '{deconv_key}' not found in dataset {data_id}")

    # Check if cell_type column already exists
    if 'cell_type' not in adata.obs:
        if context:
            await context.info("Creating cell_type column from deconvolution results")

        # Get cell types from uns
        cell_types_key = f"{deconv_key}_cell_types"
        if cell_types_key not in adata.uns:
            raise ValueError(f"Cell types not found for {deconv_key}")

        # Get deconvolution results as DataFrame
        from .tools.visualization import get_deconvolution_dataframe
        deconv_df = get_deconvolution_dataframe(adata, deconv_key)

        if deconv_df is None:
            raise ValueError(f"Could not get deconvolution dataframe for {deconv_key}")

        # Determine the dominant cell type for each spot
        dominant_cell_types = []
        for i in range(deconv_df.shape[0]):
            row = deconv_df.iloc[i]
            max_idx = row.argmax()
            dominant_cell_types.append(deconv_df.columns[max_idx])

        # Add to adata.obs
        adata.obs['cell_type'] = dominant_cell_types

        # Make it categorical
        adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

        if context:
            await context.info(f"Created cell_type annotation with {len(adata.uns[cell_types_key])} cell types")

    # Check if UMAP has been computed
    if 'X_umap' not in adata.obsm:
        if context:
            await context.info("Computing UMAP...")

        # Compute neighbors and UMAP
        import scanpy as sc
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

    # Create UMAP visualization
    import matplotlib.pyplot as plt
    import scanpy as sc
    from .utils.image_utils import fig_to_image

    # Create figure
    fig = sc.pl.umap(adata, color='cell_type', show=False, return_fig=True)

    # Convert to image
    return fig_to_image(fig, format="png")

# Note: The original visualize_umap function has been removed as it was redundant.
# Users can use the visualize function with plot_type="umap" instead.

# Add tool to get annotation visualization
@mcp.tool()
async def get_annotation_visualization(
    data_id: str,
    context: Context = None
) -> Image:
    """Get visualization of cell type annotation if available

    This is particularly useful for Tangram annotations which include spatial visualization
    of cell type distributions.

    Args:
        data_id: Dataset ID

    Returns:
        Visualization image of cell type annotation
    """
    # Validate dataset
    validate_dataset(data_id)

    # Check if dataset has cell type annotation
    adata = data_store[data_id]["adata"]

    # Check for cell type annotation or deconvolution results
    has_cell_type = "cell_type" in adata.obs
    has_deconvolution = any(key.startswith("deconvolution_") for key in adata.obsm.keys())

    if not has_cell_type and not has_deconvolution:
        raise ValueError(f"Dataset {data_id} does not have cell type annotations or deconvolution results. Run annotate or deconvolve first.")

    # If we have deconvolution results but no cell type annotation, create it
    if not has_cell_type and has_deconvolution:
        if context:
            await context.info("Creating cell type annotation from deconvolution results")

        # Find the first deconvolution key
        deconv_key = next(key for key in adata.obsm.keys() if key.startswith("deconvolution_"))

        # Get cell types from uns
        cell_types_key = f"{deconv_key}_cell_types"
        if cell_types_key in adata.uns:
            # Get deconvolution results as DataFrame
            from .tools.visualization import get_deconvolution_dataframe
            deconv_df = get_deconvolution_dataframe(adata, deconv_key)

            if deconv_df is not None:
                # Determine the dominant cell type for each spot
                dominant_cell_types = []
                for i in range(deconv_df.shape[0]):
                    row = deconv_df.iloc[i]
                    max_idx = row.argmax()
                    dominant_cell_types.append(deconv_df.columns[max_idx])

                # Add to adata.obs
                adata.obs['cell_type'] = dominant_cell_types

                # Make it categorical
                adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

                if context:
                    await context.info(f"Created cell type annotation with {len(adata.uns[cell_types_key])} cell types")
            else:
                raise ValueError(f"Could not create cell type annotation from deconvolution results. Please run annotate first.")

    # Check if there's a tangram_ct_pred in obsm
    if "tangram_ct_pred" in adata.obsm:
        if context:
            await context.info("Found Tangram cell type predictions. Creating visualization...")

        import matplotlib.pyplot as plt
        import scanpy as sc
        from .utils.image_utils import fig_to_image

        # Get cell types
        cell_types = list(adata.obsm["tangram_ct_pred"].columns)

        # Create a multi-panel figure with cell type distributions
        n_cols = min(3, len(cell_types))
        n_rows = (len(cell_types) + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
        axes = axes.flatten() if n_rows * n_cols > 1 else [axes]

        # Plot each cell type
        for i, cell_type in enumerate(cell_types):
            if i < len(axes):
                ax = axes[i]
                # Plot spatial distribution of cell type probability
                sc.pl.spatial(adata, color=cell_type, ax=ax, show=False,
                             title=f"{cell_type}", size=10, color_map='viridis')

        # Hide empty axes
        for i in range(len(cell_types), len(axes)):
            axes[i].axis('off')

        plt.tight_layout()

        # Convert figure to image
        image = fig_to_image(fig, format="png")
        plt.close(fig)

        return image
    else:
        # Create a simple visualization of cell type annotation
        if context:
            await context.info("Creating basic cell type annotation visualization...")

        import matplotlib.pyplot as plt
        import scanpy as sc
        from .utils.image_utils import fig_to_image

        fig, ax = plt.subplots(figsize=(8, 8))
        sc.pl.spatial(adata, color="cell_type", ax=ax, show=False,
                     title="Cell Type Annotation", size=10)

        # Convert figure to image
        image = fig_to_image(fig, format="png")
        plt.close(fig)

        return image

# Add tool for spatial domain identification
@mcp.tool()
async def identify_domains(
    data_id: str,
    params: SpatialDomainParameters = SpatialDomainParameters(),
    context: Context = None
) -> SpatialDomainResult:
    """Identify spatial domains in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Spatial domain identification parameters

    Returns:
        Spatial domain identification result with domain labels and visualization
    """
    if context:
        await context.info(f"Identifying spatial domains using {params.method} method")

    # Validate dataset
    validate_dataset(data_id)

    # Call spatial domain identification function
    result = await identify_spatial_domains(data_id, data_store, params, context)

    if context:
        await context.info(f"Successfully identified {result.n_domains} spatial domains")
        await context.info(f"Domain labels stored in adata.obs['{result.domain_key}']")
        if result.refined_domain_key:
            await context.info(f"Refined domain labels stored in adata.obs['{result.refined_domain_key}']")

    return result


# Add helper tool for spatial domain identification with simpler parameters
@mcp.tool()
async def identify_spatial_domains_simple(
    data_id: str,
    method: str = "stagate",
    n_domains: int = 7,
    context: Context = None
) -> SpatialDomainResult:
    """Identify spatial domains with simplified parameters

    This is a helper function that provides a simpler interface to the identify_domains tool.

    Args:
        data_id: Dataset ID
        method: Method to use (stagate, spagcn, leiden, louvain)
        n_domains: Number of spatial domains to identify

    Returns:
        Spatial domain identification result
    """
    # Create parameters object
    params = SpatialDomainParameters(
        method=method,
        n_domains=n_domains
    )

    # Call the main function
    return await identify_domains(data_id, params, context)





# Add tool for cell communication analysis
@mcp.tool()
async def analyze_communication(
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
    result = await analyze_cell_communication(data_id, data_store, params, context)

    if context:
        await context.info(f"Successfully analyzed {result.n_significant_pairs} significant LR pairs")
        if result.global_results_key:
            await context.info(f"Global results stored in adata.uns['{result.global_results_key}']")
        if result.local_analysis_performed and result.local_results_key:
            await context.info(f"Local results stored in adata.uns['{result.local_results_key}']")
        if result.top_lr_pairs:
            await context.info(f"Top LR pair: {result.top_lr_pairs[0]}")

    return result


# Add tool for spatial variable genes identification using GASTON
@mcp.tool()
async def find_spatial_genes_gaston(
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

    Returns:
        Spatial variable genes identification result with isodepth map, spatial domains,
        and gene classification into continuous/discontinuous patterns
    """
    if context:
        await context.info(f"Identifying spatial variable genes using GASTON method")

    # Validate dataset
    validate_dataset(data_id)

    # Call GASTON spatial variable genes identification function
    result = await identify_spatial_variable_genes_gaston(data_id, data_store, params, context)

    if context:
        await context.info(f"Successfully completed GASTON analysis")
        await context.info(f"Identified {result.n_spatial_domains} spatial domains")
        await context.info(f"Found {result.n_continuous_genes} genes with continuous gradients")
        await context.info(f"Found {result.n_discontinuous_genes} genes with discontinuities")
        await context.info(f"Final training loss: {result.final_loss:.6f}")
        await context.info(f"Model R²: {result.model_performance.get('r2', 'N/A')}")
        await context.info(f"Results stored with keys: {result.isodepth_key}, {result.spatial_domains_key}")

    return result


# Add helper tool for cell communication analysis with simpler parameters
@mcp.tool()
async def analyze_communication_simple(
    data_id: str,
    method: str = "liana",
    species: str = "human",
    database: str = "cellchat",
    context: Context = None
) -> CellCommunicationResult:
    """Analyze cell communication with simplified parameters

    This is a helper function that provides a simpler interface to the analyze_communication tool.

    Args:
        data_id: Dataset ID
        method: Method to use (only liana is supported)
        species: Species for ligand-receptor database (human, mouse, zebrafish)
        database: Ligand-receptor database (LIANA+ uses its own resource system)

    Returns:
        Cell communication analysis result
    """
    # Create parameters object
    params = CellCommunicationParameters(
        method=method,
        species=species,
        database=database
    )

    # Call the main function
    return await analyze_communication(data_id, params, context)


# Add resource for dataset information
@mcp.resource("dataset://{data_id}")
def get_dataset_info(data_id: str) -> Dict[str, Any]:
    """Get information about a dataset"""
    # Validate dataset
    validate_dataset(data_id)

    return data_store[data_id]


def main():
    """Run the MCP server"""
    print("Starting ChatSpatial server...")
    mcp.run(transport='stdio')


if __name__ == "__main__":
    main()