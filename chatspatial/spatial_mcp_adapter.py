"""
Spatial MCP Adapter for ChatSpatial

This module provides a clean abstraction layer between MCP protocol requirements
and ChatSpatial's spatial analysis functionality.
"""

import json
import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Union

from mcp.server.fastmcp import Context, FastMCP
from mcp.types import EmbeddedResource, ImageContent

# Import MCP improvements
from .models.data import VisualizationParameters
from .utils.tool_error_handling import mcp_tool_error_handler

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


class SpatialResourceManager:
    """Manages MCP resources for spatial data"""

    def __init__(self, data_manager: "DefaultSpatialDataManager"):
        self.data_manager = data_manager
        self._resources: Dict[str, MCPResource] = {}
        self._visualization_cache: Dict[str, Any] = {}

    async def create_dataset_resource(
        self, data_id: str, dataset_info: Dict[str, Any]
    ) -> MCPResource:
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
                "has_spatial": dataset_info.get("has_spatial", False),
            },
            content_provider=lambda: json.dumps(dataset_info, indent=2),
        )
        self._resources[resource.uri] = resource
        return resource

    async def create_result_resource(
        self, data_id: str, result_type: str, result: Any
    ) -> MCPResource:
        """Create a resource for analysis results"""
        resource = MCPResource(
            uri=f"spatial://results/{data_id}/{result_type}",
            name=f"{result_type.title()} results for {data_id}",
            mime_type="application/json",
            description=f"Analysis results: {result_type}",
            content_provider=lambda: json.dumps(
                self._serialize_result(result), indent=2
            ),
        )
        self._resources[resource.uri] = resource
        return resource

    async def create_visualization_resource(
        self, viz_id: str, image_data: bytes, metadata: Dict[str, Any]
    ) -> MCPResource:
        """Create a resource for visualization

        Note: This function only creates a resource for MCP protocol.
        The actual visualization caching is handled in server.py using cache_key
        (without timestamp) for easier retrieval by save/export functions.
        """
        resource = MCPResource(
            uri=f"spatial://visualizations/{viz_id}",
            name=metadata.get("name", f"Visualization {viz_id}"),
            mime_type="image/png",
            description=metadata.get("description", "Spatial visualization"),
            metadata=metadata,
            content_provider=lambda: image_data,
        )
        self._resources[resource.uri] = resource
        # DO NOT cache with viz_id - caching is done in server.py with cache_key
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
        """Serialize analysis results for JSON output with size control"""
        if hasattr(result, "dict"):
            result_dict = result.dict()

            # Size control for CellCommunicationResult to prevent token overflow
            if hasattr(result, "method") and "liana" in getattr(result, "method", ""):
                return self._safe_serialize_communication_result(result_dict)

            return result_dict
        elif hasattr(result, "__dict__"):
            return {k: v for k, v in result.__dict__.items() if not k.startswith("_")}
        else:
            return {"result": str(result)}

    def _safe_serialize_communication_result(
        self, result_dict: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Size-controlled serialization for cell communication results"""
        # ULTRATHINK: Prevent token overflow by applying truncation rules to large fields
        safe_dict = result_dict.copy()

        # Rule 1: Limit top_lr_pairs to prevent overflow from large L-R pair lists
        if "top_lr_pairs" in safe_dict and len(safe_dict["top_lr_pairs"]) > 10:
            safe_dict["top_lr_pairs"] = safe_dict["top_lr_pairs"][:10]
            safe_dict["top_lr_pairs_truncated"] = True

        # Rule 2: Filter statistics to remove large objects while keeping basic metrics
        if "statistics" in safe_dict and isinstance(safe_dict["statistics"], dict):
            stats = safe_dict["statistics"]
            # Keep only simple key-value pairs, exclude complex objects
            safe_stats = {
                k: v
                for k, v in stats.items()
                if not isinstance(v, (list, dict)) or len(str(v)) < 1000
            }
            safe_dict["statistics"] = safe_stats

        # Rule 3: Add size control marker for debugging
        safe_dict["_serialization_controlled"] = True

        return safe_dict


class SpatialMCPAdapter:
    """Main adapter class that bridges MCP and spatial analysis functionality"""

    def __init__(self, mcp_server: FastMCP, data_manager: "DefaultSpatialDataManager"):
        self.mcp = mcp_server
        self.data_manager = data_manager
        self.resource_manager = SpatialResourceManager(data_manager)
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
                open_world_hint=True,
            ),
            "preprocess_data": MCPToolMetadata(
                name="preprocess_data",
                title="Preprocess Spatial Data",
                description="Preprocess and normalize spatial data",
                read_only_hint=False,
                idempotent_hint=False,
                open_world_hint=False,
            ),
            "visualize_data": MCPToolMetadata(
                name="visualize_data",
                title="Visualize Spatial Data",
                description="Create visualizations of spatial data",
                read_only_hint=True,
                idempotent_hint=True,
                open_world_hint=False,
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
                open_world_hint=False,
            ),
            "find_markers": MCPToolMetadata(
                name="find_markers",
                title="Find Marker Genes",
                description="Identify differentially expressed genes",
                read_only_hint=True,
                idempotent_hint=True,
                open_world_hint=False,
            ),
            "integrate_samples": MCPToolMetadata(
                name="integrate_samples",
                title="Integrate Multiple Samples",
                description="Integrate multiple spatial datasets",
                read_only_hint=False,
                idempotent_hint=False,
                open_world_hint=False,
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
                open_world_hint=False,
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
                open_world_hint=False,
            ),
            "find_spatial_genes": MCPToolMetadata(
                name="find_spatial_genes",
                title="Find Spatial Variable Genes",
                description="Identify spatially variable genes",
                read_only_hint=False,
                idempotent_hint=True,
                open_world_hint=False,
            ),
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
                "metadata": r.metadata,
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
                "contents": [
                    {
                        "uri": uri,
                        "mimeType": "image/png",
                        "blob": base64.b64encode(content).decode("utf-8"),
                    }
                ]
            }
        else:
            # Text content
            return {
                "contents": [
                    {"uri": uri, "mimeType": "application/json", "text": content}
                ]
            }

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
        context: Optional[Context] = None,
    ) -> Optional[Union[ImageContent, tuple[ImageContent, EmbeddedResource]]]:
        """Create visualization from analysis result"""
        try:
            # Import visualization function
            from .tools.visualization import visualize_data

            # Create visualization parameters
            params = VisualizationParameters(plot_type=plot_type)

            # Get dataset
            dataset_info = await self.data_manager.get_dataset(data_id)

            # Call visualization
            image = await visualize_data(
                data_id, {"data_id": dataset_info}, params, context
            )

            if image:
                # Create resource for the visualization
                import time

                # Use consistent cache key (no timestamp for easier lookup)
                cache_key = f"{data_id}_{plot_type}"
                viz_id = f"{cache_key}_{int(time.time())}"  # Resource ID with timestamp

                metadata = {
                    "data_id": data_id,
                    "plot_type": plot_type,
                    "timestamp": int(time.time()),
                    "name": f"{plot_type} visualization",
                    "description": f"Visualization of {plot_type} for dataset {data_id}",
                }

                # Decode base64 string to bytes before caching
                import base64

                image_bytes = base64.b64decode(image.data)

                # Store with consistent cache_key (for save_visualization lookup)
                await self.resource_manager.create_visualization_resource(
                    cache_key, image_bytes, metadata
                )

                if context:
                    await context.info(
                        f"Created visualization resource: spatial://visualizations/{viz_id}"
                    )

            return image

        except Exception as e:
            logger.error(f"Error creating visualization: {e}")
            if context:
                await context.error(f"Failed to create visualization: {str(e)}")
            return None


class DefaultSpatialDataManager:
    """Default implementation of spatial data management"""

    def __init__(self):
        self.data_store: Dict[str, Any] = {}
        self._next_id = 1

    async def load_dataset(
        self, path: str, data_type: str, name: Optional[str] = None
    ) -> str:
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
                "n_genes": info.get("n_genes", 0),
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
    data_manager: Optional[DefaultSpatialDataManager] = None,
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
    # Note: Uncomment these when FastMCP supports resource decorators
    # @mcp.list_resources
    # async def handle_list_resources():
    #     return await adapter.handle_resource_list()
    #
    # @mcp.read_resource
    # async def handle_read_resource(uri: str):
    #     return await adapter.handle_resource_read(uri)

    return mcp, adapter
