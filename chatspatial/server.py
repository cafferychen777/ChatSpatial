"""
Main server implementation for ChatSpatial.
"""

from typing import Dict, Any, List, Optional
import matplotlib.pyplot as plt
import numpy as np

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
    DeconvolutionParameters
)
from .models.analysis import (
    PreprocessingResult,
    DifferentialExpressionResult,
    AnnotationResult,
    SpatialAnalysisResult,
    RNAVelocityResult,
    TrajectoryResult,
    IntegrationResult,
    DeconvolutionResult
)
from .tools.preprocessing import preprocess_data
from .tools.visualization import visualize_data
from .tools.annotation import annotate_cell_types
from .tools.spatial_analysis import analyze_spatial as perform_spatial_analysis
from .tools.spatial_analysis import analyze_spatial_with_image
from .tools.spatial_analysis import analyze_spatial_unified
from .tools.differential import differential_expression
from .tools.trajectory import analyze_rna_velocity, analyze_trajectory
from .tools.integration import integrate_samples as perform_integration
from .tools.deconvolution import deconvolve_spatial_data
from .utils.data_loader import load_spatial_data

# Create MCP server
mcp = FastMCP("ChatSpatial")

# Store for loaded datasets
data_store: Dict[str, Any] = {}


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
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    # Call preprocessing function
    result = await preprocess_data(data_id, data_store, params, context)

    return result


@mcp.tool()
async def visualize(
    data_id: str,
    params: VisualizationParameters = VisualizationParameters(),
    context: Context = None
) -> Image:
    """Visualize spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Visualization parameters

    Returns:
        Visualization image
    """
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

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
        Annotation result
    """
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    # Call annotation function
    return await annotate_cell_types(data_id, data_store, params, context)


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
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    # Call unified spatial analysis function with return_type="image"
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
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

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
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    # Call RNA velocity analysis function
    result = await analyze_rna_velocity(data_id, data_store, params, context)

    # Return the visualization image directly
    if result.visualization is None:
        # Create a simple placeholder image if no visualization is available
        return create_placeholder_image("No velocity visualization available")

    return result.visualization


# Add trajectory analysis tool
@mcp.tool()
async def analyze_spatial_trajectory(
    data_id: str,
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Context = None
) -> Image:
    """Analyze trajectory and cell state transitions in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        params: Trajectory analysis parameters

    Returns:
        Pseudotime visualization image
    """
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    # Call trajectory analysis function
    result = await analyze_trajectory(data_id, data_store, params, context)

    # Return the pseudotime visualization image directly
    return result.pseudotime_visualization


# Add tool to visualize trajectory velocity
@mcp.tool()
async def visualize_trajectory_velocity(
    data_id: str,
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Context = None
) -> Image:
    """Visualize velocity streams from trajectory analysis

    Args:
        data_id: Dataset ID
        params: Trajectory analysis parameters

    Returns:
        Velocity stream visualization image
    """
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    # Call trajectory analysis function
    result = await analyze_trajectory(data_id, data_store, params, context)

    # Return the velocity visualization image directly
    if result.velocity_visualization is None:
        # Create a simple placeholder image if no visualization is available
        return create_placeholder_image("No velocity stream visualization available")

    return result.velocity_visualization


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
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    # Don't include image to save bandwidth
    params.include_image = False

    # Call unified spatial analysis function with return_type="result"
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
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found")

    # Call integration function
    return await perform_integration(data_ids, data_store, params, context)


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
    if context:
        await context.info(f"Deconvolving spatial data using {params.method} method")

    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    # Validate reference data if provided
    if params.reference_data_id and params.reference_data_id not in data_store:
        raise ValueError(f"Reference dataset {params.reference_data_id} not found")

    # Call deconvolution function
    return await deconvolve_spatial_data(data_id, data_store, params, context)


# Initialize the MCP server
mcp = FastMCP("chatspatial")

# Register all tools
mcp.tool()(load_data)
mcp.tool()(preprocess)
mcp.tool()(visualize)
mcp.tool()(annotate)
mcp.tool()(analyze_spatial)
mcp.tool()(find_markers)
mcp.tool()(list_datasets)
mcp.tool()(analyze_velocity)
mcp.tool()(analyze_spatial_trajectory)
mcp.tool()(visualize_trajectory_velocity)
mcp.tool()(get_spatial_analysis_stats)
mcp.tool()(integrate_samples)
mcp.tool()(deconvolve)

# Add resource for dataset information
@mcp.resource("dataset://{data_id}")
def get_dataset_info(data_id: str) -> Dict[str, Any]:
    """Get information about a dataset"""
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")

    return data_store[data_id]

# Global data store for datasets
data_store = {}


def main():
    """Run the MCP server"""
    print("Starting ChatSpatial server...")
    mcp.run(transport='stdio')


if __name__ == "__main__":
    main()