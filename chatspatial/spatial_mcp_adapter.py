"""
Spatial MCP Adapter for ChatSpatial

This module provides a clean abstraction layer between MCP protocol requirements
and ChatSpatial's spatial analysis functionality.
"""

from typing import Dict, Any, List, Optional, Union, Callable
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import json
import asyncio
import logging
from pathlib import Path

from mcp.server.fastmcp import FastMCP, Context
from mcp.server.fastmcp.utilities.types import Image

from .utils.tool_error_handling import mcp_tool_error_handler
from .utils.mcp_parameter_handler import manual_parameter_validation
from .utils.image_utils import fig_to_image, create_placeholder_image

from .models.data import (
    SpatialDataset,
    AnalysisParameters,
    VisualizationParameters,
    AnnotationParameters,
    SpatialAnalysisParameters,
    IntegrationParameters,
    DeconvolutionParameters,
    SpatialDomainParameters,
    CellCommunicationParameters,
    SpatialVariableGenesParameters
)

from .models.analysis import (
    PreprocessingResult,
    DifferentialExpressionResult,
    AnnotationResult,
    SpatialAnalysisResult,
    IntegrationResult,
    DeconvolutionResult,
    SpatialDomainResult,
    CellCommunicationResult,
    SpatialVariableGenesResult
)

# Import MCP improvements
from .mcp.errors import ErrorType, format_mcp_error
from .mcp.resources import Resource, ResourceManager
from .mcp.prompts import Prompt, PromptArgument, PromptManager
from .mcp.annotations import get_tool_annotation, get_all_tools_with_annotations

logger = logging.getLogger(__name__)


@dataclass
class MCPResource:
    """MCP Resource representation"""
    uri: str
    name: str
    mime_type: str
    description: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    content_provider: Optional[Callable[[], Union[str, bytes]]] = None


@dataclass
class MCPPrompt:
    """MCP Prompt representation"""
    name: str
    description: str
    arguments: List[Dict[str, Any]] = field(default_factory=list)
    handler: Optional[Callable[[Dict[str, Any]], Dict[str, Any]]] = None


@dataclass
class MCPToolMetadata:
    """Enhanced tool metadata including MCP annotations"""
    name: str
    title: str
    description: str
    read_only_hint: bool = False
    destructive_hint: bool = False
    idempotent_hint: bool = True
    open_world_hint: bool = False
    parameters_schema: Optional[Dict[str, Any]] = None


class SpatialDataManager(ABC):
    """Abstract interface for spatial data management"""
    
    @abstractmethod
    async def load_dataset(self, path: str, data_type: str, name: Optional[str] = None) -> str:
        """Load a spatial dataset and return its ID"""
        pass
    
    @abstractmethod
    async def get_dataset(self, data_id: str) -> Any:
        """Get a dataset by ID"""
        pass
    
    @abstractmethod
    async def list_datasets(self) -> List[Dict[str, Any]]:
        """List all loaded datasets"""
        pass
    
    @abstractmethod
    async def save_result(self, data_id: str, result_type: str, result: Any) -> None:
        """Save analysis results"""
        pass
    
    @abstractmethod
    async def get_result(self, data_id: str, result_type: str) -> Any:
        """Get analysis results"""
        pass


class SpatialResourceManager:
    """Manages MCP resources for spatial data"""
    
    def __init__(self, data_manager: SpatialDataManager):
        self.data_manager = data_manager
        self._resources: Dict[str, MCPResource] = {}
        self._visualization_cache: Dict[str, Any] = {}
        
    async def create_dataset_resource(self, data_id: str, dataset_info: Dict[str, Any]) -> MCPResource:
        """Create a resource for a dataset"""
        resource = MCPResource(
            uri=f"spatial://datasets/{data_id}",
            name=dataset_info.get("name", f"Dataset {data_id}"),
            mime_type="application/x-anndata",
            description=f"Spatial dataset with {dataset_info.get('n_cells', 0)} cells",
            metadata={
                "n_cells": dataset_info.get("n_cells", 0),
                "n_genes": dataset_info.get("n_genes", 0),
                "data_type": dataset_info.get("type", "unknown"),
                "has_spatial": dataset_info.get("has_spatial", False)
            },
            content_provider=lambda: json.dumps(dataset_info, indent=2)
        )
        self._resources[resource.uri] = resource
        return resource
    
    async def create_result_resource(self, data_id: str, result_type: str, result: Any) -> MCPResource:
        """Create a resource for analysis results"""
        resource = MCPResource(
            uri=f"spatial://results/{data_id}/{result_type}",
            name=f"{result_type.title()} results for {data_id}",
            mime_type="application/json",
            description=f"Analysis results: {result_type}",
            content_provider=lambda: json.dumps(self._serialize_result(result), indent=2)
        )
        self._resources[resource.uri] = resource
        return resource
    
    async def create_visualization_resource(self, viz_id: str, image_data: bytes, metadata: Dict[str, Any]) -> MCPResource:
        """Create a resource for visualization"""
        resource = MCPResource(
            uri=f"spatial://visualizations/{viz_id}",
            name=metadata.get("name", f"Visualization {viz_id}"),
            mime_type="image/png",
            description=metadata.get("description", "Spatial visualization"),
            metadata=metadata,
            content_provider=lambda: image_data
        )
        self._resources[resource.uri] = resource
        self._visualization_cache[viz_id] = image_data
        return resource
    
    async def get_resource(self, uri: str) -> Optional[MCPResource]:
        """Get a resource by URI"""
        return self._resources.get(uri)
    
    async def list_resources(self) -> List[MCPResource]:
        """List all available resources"""
        return list(self._resources.values())
    
    async def read_resource_content(self, uri: str) -> Union[str, bytes]:
        """Read resource content"""
        resource = await self.get_resource(uri)
        if not resource:
            raise ValueError(f"Resource not found: {uri}")
        
        if resource.content_provider:
            return resource.content_provider()
        
        return json.dumps({"error": "No content available"})
    
    def _serialize_result(self, result: Any) -> Dict[str, Any]:
        """Serialize analysis results for JSON output"""
        if hasattr(result, 'dict'):
            return result.dict()
        elif hasattr(result, '__dict__'):
            return {k: v for k, v in result.__dict__.items() if not k.startswith('_')}
        else:
            return {"result": str(result)}


class SpatialPromptManager:
    """Manages MCP prompts for spatial analysis"""
    
    def __init__(self, data_manager: SpatialDataManager):
        self.data_manager = data_manager
        self._prompts: Dict[str, MCPPrompt] = {}
        self._initialize_prompts()
    
    def _initialize_prompts(self):
        """Initialize spatial analysis prompts"""
        self._prompts = {
            "analyze-spatial-expression": MCPPrompt(
                name="analyze-spatial-expression",
                description="Analyze spatial gene expression patterns",
                arguments=[
                    {"name": "genes", "description": "Genes to analyze", "required": True},
                    {"name": "method", "description": "Analysis method", "required": False}
                ]
            ),
            "find-cell-types": MCPPrompt(
                name="find-cell-types",
                description="Identify cell types in spatial data",
                arguments=[
                    {"name": "method", "description": "Cell type identification method", "required": False},
                    {"name": "reference_data", "description": "Reference dataset path", "required": False}
                ]
            ),
            "compare-conditions": MCPPrompt(
                name="compare-conditions",
                description="Compare spatial patterns between conditions",
                arguments=[
                    {"name": "condition_key", "description": "Column defining conditions", "required": True},
                    {"name": "groups", "description": "Groups to compare", "required": True}
                ]
            ),
            "generate-visualization": MCPPrompt(
                name="generate-visualization",
                description="Generate spatial visualization",
                arguments=[
                    {"name": "plot_type", "description": "Type of visualization", "required": True},
                    {"name": "feature", "description": "Feature(s) to visualize (single gene or list of genes)", "required": False},
                    {"name": "save_path", "description": "Path to save figure", "required": False}
                ]
            ),
            "quality-control": MCPPrompt(
                name="quality-control",
                description="Perform quality control on spatial data",
                arguments=[
                    {"name": "metrics", "description": "QC metrics to compute", "required": False},
                    {"name": "thresholds", "description": "QC thresholds", "required": False}
                ]
            ),
            "batch-correction": MCPPrompt(
                name="batch-correction",
                description="Correct batch effects in integrated data",
                arguments=[
                    {"name": "batch_key", "description": "Column defining batches", "required": True},
                    {"name": "method", "description": "Batch correction method", "required": False}
                ]
            ),
            "spatial-clustering": MCPPrompt(
                name="spatial-clustering",
                description="Perform spatial clustering analysis",
                arguments=[
                    {"name": "method", "description": "Clustering method", "required": False},
                    {"name": "n_clusters", "description": "Number of clusters", "required": False},
                    {"name": "resolution", "description": "Clustering resolution", "required": False}
                ]
            ),
            "cellular-communication": MCPPrompt(
                name="cellular-communication",
                description="Analyze cell-cell communication",
                arguments=[
                    {"name": "method", "description": "Communication analysis method", "required": False},
                    {"name": "lr_database", "description": "Ligand-receptor database", "required": False}
                ]
            ),
            "gene-enrichment": MCPPrompt(
                name="gene-enrichment",
                description="Perform gene set enrichment analysis",
                arguments=[
                    {"name": "gene_sets", "description": "Gene sets to test", "required": True},
                    {"name": "method", "description": "Enrichment method", "required": False}
                ]
            ),
            "spatial-deconvolution": MCPPrompt(
                name="spatial-deconvolution",
                description="Deconvolve spatial spots into cell types",
                arguments=[
                    {"name": "reference_data", "description": "Single-cell reference data", "required": True},
                    {"name": "method", "description": "Deconvolution method", "required": False}
                ]
            )
        }
    
    async def get_prompt(self, name: str) -> Optional[MCPPrompt]:
        """Get a prompt by name"""
        return self._prompts.get(name)
    
    async def list_prompts(self) -> List[MCPPrompt]:
        """List all available prompts"""
        return list(self._prompts.values())


class SpatialMCPAdapter:
    """Main adapter class that bridges MCP and spatial analysis functionality"""
    
    def __init__(self, mcp_server: FastMCP, data_manager: SpatialDataManager):
        self.mcp = mcp_server
        self.data_manager = data_manager
        self.resource_manager = SpatialResourceManager(data_manager)
        self.prompt_manager = SpatialPromptManager(data_manager)
        self._tool_metadata: Dict[str, MCPToolMetadata] = {}
        self._initialize_tools()
    
    def _initialize_tools(self):
        """Initialize tool metadata"""
        self._tool_metadata = {
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
                idempotent_hint=False,
                open_world_hint=False
            ),
            "visualize_data": MCPToolMetadata(
                name="visualize_data",
                title="Visualize Spatial Data",
                description="Create visualizations of spatial data",
                read_only_hint=True,
                idempotent_hint=True,
                open_world_hint=False
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
                idempotent_hint=True,
                open_world_hint=False
            ),
            "find_markers": MCPToolMetadata(
                name="find_markers",
                title="Find Marker Genes",
                description="Identify differentially expressed genes",
                read_only_hint=True,
                idempotent_hint=True,
                open_world_hint=False
            ),
            "integrate_samples": MCPToolMetadata(
                name="integrate_samples",
                title="Integrate Multiple Samples",
                description="Integrate multiple spatial datasets",
                read_only_hint=False,
                idempotent_hint=False,
                open_world_hint=False
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
                idempotent_hint=False,
                open_world_hint=False
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
                idempotent_hint=True,
                open_world_hint=False
            ),
            "find_spatial_genes": MCPToolMetadata(
                name="find_spatial_genes",
                title="Find Spatial Variable Genes",
                description="Identify spatially variable genes",
                read_only_hint=False,
                idempotent_hint=True,
                open_world_hint=False
            )
        }
    
    def get_tool_metadata(self, tool_name: str) -> Optional[MCPToolMetadata]:
        """Get metadata for a tool"""
        return self._tool_metadata.get(tool_name)
    
    def get_all_tools_metadata(self) -> List[MCPToolMetadata]:
        """Get metadata for all tools"""
        return list(self._tool_metadata.values())
    
    async def handle_resource_list(self) -> List[Dict[str, Any]]:
        """Handle MCP resource list request"""
        resources = await self.resource_manager.list_resources()
        return [
            {
                "uri": r.uri,
                "name": r.name,
                "mimeType": r.mime_type,
                "description": r.description,
                "metadata": r.metadata
            }
            for r in resources
        ]
    
    async def handle_resource_read(self, uri: str) -> Dict[str, Any]:
        """Handle MCP resource read request"""
        content = await self.resource_manager.read_resource_content(uri)
        
        if isinstance(content, bytes):
            # Binary content (like images)
            import base64
            return {
                "contents": [{
                    "uri": uri,
                    "mimeType": "image/png",
                    "blob": base64.b64encode(content).decode('utf-8')
                }]
            }
        else:
            # Text content
            return {
                "contents": [{
                    "uri": uri,
                    "mimeType": "application/json",
                    "text": content
                }]
            }
    
    async def handle_prompt_list(self) -> List[Dict[str, Any]]:
        """Handle MCP prompt list request"""
        prompts = await self.prompt_manager.list_prompts()
        return [
            {
                "name": p.name,
                "description": p.description,
                "arguments": p.arguments
            }
            for p in prompts
        ]
    
    async def handle_prompt_execute(self, name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Handle MCP prompt execution request"""
        result = await self.prompt_manager.execute_prompt(name, arguments)
        return result
    
    def register_tool(self, tool_func: Callable, metadata: MCPToolMetadata) -> None:
        """Register a tool with the MCP server"""
        # Wrap the tool with MCP decorators
        decorated_tool = mcp_tool_error_handler()(tool_func)
        
        # Register with FastMCP
        self.mcp.tool()(decorated_tool)
        
        # Store metadata
        self._tool_metadata[metadata.name] = metadata
    
    async def create_visualization_from_result(
        self,
        data_id: str,
        plot_type: str,
        result: Any,
        context: Optional[Context] = None
    ) -> Optional[Image]:
        """Create visualization from analysis result"""
        try:
            # Import visualization function
            from .tools.visualization import visualize_data
            
            # Create visualization parameters
            params = VisualizationParameters(plot_type=plot_type)
            
            # Get dataset
            dataset_info = await self.data_manager.get_dataset(data_id)
            
            # Call visualization
            image = await visualize_data(data_id, {"data_id": dataset_info}, params, context)
            
            if image:
                # Create resource for the visualization
                import time
                viz_id = f"{data_id}_{plot_type}_{int(time.time())}"
                
                metadata = {
                    "data_id": data_id,
                    "plot_type": plot_type,
                    "timestamp": int(time.time()),
                    "name": f"{plot_type} visualization",
                    "description": f"Visualization of {plot_type} for dataset {data_id}"
                }
                
                await self.resource_manager.create_visualization_resource(
                    viz_id, image.data, metadata
                )
                
                if context:
                    await context.info(f"Created visualization resource: spatial://visualizations/{viz_id}")
            
            return image
            
        except Exception as e:
            logger.error(f"Error creating visualization: {e}")
            if context:
                await context.error(f"Failed to create visualization: {str(e)}")
            return None


class DefaultSpatialDataManager(SpatialDataManager):
    """Default implementation of SpatialDataManager"""
    
    def __init__(self):
        self.data_store: Dict[str, Any] = {}
        self._next_id = 1
    
    async def load_dataset(self, path: str, data_type: str, name: Optional[str] = None) -> str:
        """Load a spatial dataset and return its ID"""
        from .utils.data_loader import load_spatial_data
        
        # Load data
        dataset_info = await load_spatial_data(path, data_type, name)
        
        # Generate ID
        data_id = f"data_{self._next_id}"
        self._next_id += 1
        
        # Store data
        self.data_store[data_id] = dataset_info
        
        return data_id
    
    async def get_dataset(self, data_id: str) -> Any:
        """Get a dataset by ID"""
        if data_id not in self.data_store:
            raise ValueError(f"Dataset {data_id} not found")
        return self.data_store[data_id]
    
    async def list_datasets(self) -> List[Dict[str, Any]]:
        """List all loaded datasets"""
        return [
            {
                "id": data_id,
                "name": info.get("name", f"Dataset {data_id}"),
                "type": info.get("type", "unknown"),
                "n_cells": info.get("n_cells", 0),
                "n_genes": info.get("n_genes", 0)
            }
            for data_id, info in self.data_store.items()
        ]
    
    async def save_result(self, data_id: str, result_type: str, result: Any) -> None:
        """Save analysis results"""
        if data_id not in self.data_store:
            raise ValueError(f"Dataset {data_id} not found")
        
        if "results" not in self.data_store[data_id]:
            self.data_store[data_id]["results"] = {}
        
        self.data_store[data_id]["results"][result_type] = result
    
    async def get_result(self, data_id: str, result_type: str) -> Any:
        """Get analysis results"""
        if data_id not in self.data_store:
            raise ValueError(f"Dataset {data_id} not found")
        
        results = self.data_store[data_id].get("results", {})
        if result_type not in results:
            raise ValueError(f"No {result_type} results found for dataset {data_id}")
        
        return results[result_type]


def create_spatial_mcp_server(
    server_name: str = "ChatSpatial",
    data_manager: Optional[SpatialDataManager] = None
) -> tuple[FastMCP, SpatialMCPAdapter]:
    """
    Create and configure a spatial MCP server with adapter
    
    Args:
        server_name: Name of the MCP server
        data_manager: Optional custom data manager (uses default if None)
    
    Returns:
        Tuple of (FastMCP server instance, SpatialMCPAdapter instance)
    """
    # Create MCP server
    mcp = FastMCP(server_name)
    
    # Create data manager if not provided
    if data_manager is None:
        data_manager = DefaultSpatialDataManager()
    
    # Create adapter
    adapter = SpatialMCPAdapter(mcp, data_manager)
    
    # Configure resource handlers
    # Note: Uncomment these when FastMCP supports resource/prompt decorators
    # @mcp.list_resources
    # async def handle_list_resources():
    #     return await adapter.handle_resource_list()
    # 
    # @mcp.read_resource
    # async def handle_read_resource(uri: str):
    #     return await adapter.handle_resource_read(uri)
    # 
    # # Configure prompt handlers
    # @mcp.list_prompts
    # async def handle_list_prompts():
    #     return await adapter.handle_prompt_list()
    # 
    # @mcp.get_prompt
    # async def handle_get_prompt(name: str, arguments: Dict[str, Any]):
    #     return await adapter.handle_prompt_execute(name, arguments)
    
    return mcp, adapter