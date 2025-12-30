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
from mcp.types import EmbeddedResource, ImageContent, ToolAnnotations

# Import MCP improvements
from .models.data import VisualizationParameters
from .utils.exceptions import DataNotFoundError, ParameterError

logger = logging.getLogger(__name__)


# =============================================================================
# TOOL ANNOTATIONS - Single Source of Truth
# =============================================================================
# These annotations are passed to FastMCP's @mcp.tool() decorator to inform
# LLM clients about tool behavior characteristics.
#
# Annotation meanings (from MCP spec):
# - readOnlyHint: Tool only reads data, doesn't modify state
# - idempotentHint: Repeated calls with same args have no additional effect
# - openWorldHint: Tool may interact with external entities (network, files)
# =============================================================================

TOOL_ANNOTATIONS: Dict[str, ToolAnnotations] = {
    # Data I/O tools
    "load_data": ToolAnnotations(
        readOnlyHint=True,  # Reads from filesystem, doesn't modify data
        idempotentHint=True,  # Loading same file yields same result
        openWorldHint=True,  # Accesses filesystem
    ),
    "save_data": ToolAnnotations(
        readOnlyHint=False,  # Writes to filesystem
        idempotentHint=True,  # Saving same data to same path is idempotent
        openWorldHint=True,  # Accesses filesystem
    ),
    # Preprocessing - modifies data in-place
    "preprocess_data": ToolAnnotations(
        readOnlyHint=False,  # Modifies adata in-place
        idempotentHint=False,  # Re-running changes state
    ),
    # Visualization - read-only analysis
    "visualize_data": ToolAnnotations(
        readOnlyHint=True,  # Only reads data to generate plots
        idempotentHint=True,  # Same params yield same plot
    ),
    "save_visualization": ToolAnnotations(
        readOnlyHint=False,  # Writes to filesystem
        idempotentHint=True,  # Saving same viz is idempotent
        openWorldHint=True,  # Accesses filesystem
    ),
    "export_all_visualizations": ToolAnnotations(
        readOnlyHint=False,
        idempotentHint=True,
        openWorldHint=True,
    ),
    "clear_visualization_cache": ToolAnnotations(
        readOnlyHint=False,  # Clears cache
        idempotentHint=True,  # Clearing empty cache is idempotent
    ),
    # Analysis tools - modify adata by adding results
    "annotate_cell_types": ToolAnnotations(
        readOnlyHint=False,  # Adds cell type annotations to adata.obs
        idempotentHint=False,  # Re-running may yield different results
        openWorldHint=True,  # May use external references
    ),
    "analyze_spatial_statistics": ToolAnnotations(
        readOnlyHint=False,  # Adds statistics to adata.uns
        idempotentHint=True,  # Same params yield same statistics
    ),
    "find_markers": ToolAnnotations(
        readOnlyHint=True,  # Computes markers without modifying adata
        idempotentHint=True,  # Deterministic computation
    ),
    "analyze_velocity_data": ToolAnnotations(
        readOnlyHint=False,  # Adds velocity to adata
        idempotentHint=False,  # Stochastic methods
    ),
    "analyze_trajectory_data": ToolAnnotations(
        readOnlyHint=False,  # Adds trajectory info to adata
        idempotentHint=False,  # May have stochastic elements
    ),
    "integrate_samples": ToolAnnotations(
        readOnlyHint=False,  # Creates new integrated dataset
        idempotentHint=False,  # Creates new dataset each time
    ),
    "deconvolve_data": ToolAnnotations(
        readOnlyHint=False,  # Adds deconvolution results to adata
        idempotentHint=False,  # Deep learning methods are stochastic
        openWorldHint=True,  # May use external references
    ),
    "identify_spatial_domains": ToolAnnotations(
        readOnlyHint=False,  # Adds domain labels to adata.obs
        idempotentHint=False,  # Clustering can vary
    ),
    "analyze_cell_communication": ToolAnnotations(
        readOnlyHint=False,  # Adds communication results to adata.uns
        idempotentHint=True,  # Deterministic given same inputs
        openWorldHint=True,  # Uses LR databases
    ),
    "analyze_enrichment": ToolAnnotations(
        readOnlyHint=False,  # Adds enrichment scores to adata
        idempotentHint=True,  # Deterministic
    ),
    "find_spatial_genes": ToolAnnotations(
        readOnlyHint=False,  # Adds spatial gene info to adata.var
        idempotentHint=True,  # Deterministic methods
    ),
    "analyze_cnv": ToolAnnotations(
        readOnlyHint=False,  # Adds CNV results to adata
        idempotentHint=True,  # Deterministic
    ),
    "register_spatial_data": ToolAnnotations(
        readOnlyHint=False,  # Modifies spatial coordinates
        idempotentHint=False,  # Registration can vary
    ),
}


def get_tool_annotations(tool_name: str) -> ToolAnnotations:
    """Get annotations for a tool by name.

    Args:
        tool_name: Name of the tool (e.g., 'load_data', 'preprocess_data')

    Returns:
        ToolAnnotations object for the tool. Returns conservative defaults
        if tool is not in registry.

    Usage:
        @mcp.tool(annotations=get_tool_annotations("load_data"))
        async def load_data(...): ...
    """
    return TOOL_ANNOTATIONS.get(
        tool_name,
        # Conservative defaults: assume tool modifies state and is not idempotent
        ToolAnnotations(readOnlyHint=False, idempotentHint=False),
    )


@dataclass
class MCPResource:
    """MCP Resource representation"""

    uri: str
    name: str
    mime_type: str
    description: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    content_provider: Optional[Callable[[], Union[str, bytes]]] = None


class SpatialResourceManager:
    """Manages MCP resources for spatial data"""

    def __init__(self, data_manager: "DefaultSpatialDataManager"):
        self.data_manager = data_manager
        self._resources: Dict[str, MCPResource] = {}

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
            raise DataNotFoundError(f"Resource not found: {uri}")

        if resource.content_provider:
            return resource.content_provider()

        raise DataNotFoundError(f"Resource has no content provider: {uri}")

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
        # Session-level cache for rendered visualizations
        # Stored at adapter level (not ResourceManager) for clear ownership
        self.visualization_cache: Dict[str, Any] = {}

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
            params = VisualizationParameters(plot_type=plot_type)  # type: ignore[call-arg]

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
            raise DataNotFoundError(f"Dataset {data_id} not found")
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
            raise DataNotFoundError(f"Dataset {data_id} not found")

        if "results" not in self.data_store[data_id]:
            self.data_store[data_id]["results"] = {}

        self.data_store[data_id]["results"][result_type] = result

    async def get_result(self, data_id: str, result_type: str) -> Any:
        """Get analysis results"""
        if data_id not in self.data_store:
            raise DataNotFoundError(f"Dataset {data_id} not found")

        results = self.data_store[data_id].get("results", {})
        if result_type not in results:
            raise DataNotFoundError(
                f"No {result_type} results found for dataset {data_id}"
            )

        return results[result_type]

    def dataset_exists(self, data_id: str) -> bool:
        """Check if a dataset exists.

        Args:
            data_id: Dataset identifier

        Returns:
            True if the dataset exists, False otherwise
        """
        return data_id in self.data_store

    async def update_adata(self, data_id: str, adata: Any) -> None:
        """Update the adata object for an existing dataset.

        Use this when preprocessing creates a new adata object (e.g., copy,
        subsample, or format conversion).

        Args:
            data_id: Dataset identifier
            adata: New AnnData object to store

        Raises:
            DataNotFoundError: If dataset not found
        """
        if data_id not in self.data_store:
            raise DataNotFoundError(f"Dataset {data_id} not found")
        self.data_store[data_id]["adata"] = adata

    async def create_dataset(
        self,
        data_id: str,
        adata: Any,
        name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Create a new dataset with specified ID.

        Use this when creating derived datasets (e.g., integration results,
        subset data).

        Args:
            data_id: Unique identifier for the new dataset
            adata: AnnData object to store
            name: Optional display name for the dataset
            metadata: Optional additional metadata dict

        Raises:
            ParameterError: If dataset with same ID already exists
        """
        if data_id in self.data_store:
            raise ParameterError(
                f"Dataset {data_id} already exists. Use update_adata() to update."
            )
        dataset_info: Dict[str, Any] = {"adata": adata}
        if name:
            dataset_info["name"] = name
        if metadata:
            dataset_info.update(metadata)
        self.data_store[data_id] = dataset_info


@dataclass
class ToolContext:
    """Unified context for ChatSpatial tool execution.

    This class provides a clean interface for tools to access data and logging
    without the redundant data_store dict wrapping pattern.

    Design Rationale:
    - Python dict assignment is reference, not copy. The old pattern of wrapping
      dataset_info in a temp dict and "writing back" was completely unnecessary.
    - Tools should access adata directly via get_adata(), not through dict wrapping.
    - Logging methods fall back gracefully when MCP context is unavailable.

    Logging Strategy:
    - User-visible messages: await ctx.info(), await ctx.warning(), await ctx.error()
      These appear in Claude's conversation and provide user-friendly progress updates.
    - Developer debugging: ctx.debug()
      This writes to Python logger for debugging, not visible to users.

    Usage:
        async def my_tool(data_id: str, ctx: ToolContext, params: Params) -> Result:
            adata = await ctx.get_adata(data_id)
            await ctx.info(f"Processing {adata.n_obs} cells")  # User sees this
            ctx.debug(f"Internal state: {some_detail}")  # Developer log only
            # ... analysis logic ...
            return result
    """

    _data_manager: "DefaultSpatialDataManager"
    _mcp_context: Optional[Context] = None
    _visualization_cache: Optional[Dict[str, Any]] = None
    _logger: Optional[logging.Logger] = field(default=None, repr=False)

    def __post_init__(self) -> None:
        """Initialize the logger for debug messages."""
        if self._logger is None:
            self._logger = logging.getLogger("chatspatial.tools")

    def debug(self, msg: str) -> None:
        """Log debug message for developers (not visible to users).

        Use this for detailed technical information that helps with debugging
        but would be noise for end users. These messages go to Python logger.

        Args:
            msg: Debug message to log
        """
        if self._logger:
            self._logger.debug(msg)

    def log_config(self, title: str, config: Dict[str, Any]) -> None:
        """Log configuration details for developers.

        Convenience method for logging parameter configurations in a
        structured format. Goes to Python logger, not user-visible.

        Args:
            title: Configuration section title
            config: Dictionary of configuration key-value pairs
        """
        if self._logger:
            self._logger.info("=" * 50)
            self._logger.info(f"{title}:")
            for key, value in config.items():
                self._logger.info(f"  {key}: {value}")
            self._logger.info("=" * 50)

    async def get_adata(self, data_id: str) -> Any:
        """Get AnnData object directly by ID.

        This is the primary data access method for tools. Returns the AnnData
        object directly without intermediate dict wrapping.

        Args:
            data_id: Dataset identifier

        Returns:
            AnnData object for the dataset

        Raises:
            ValueError: If dataset not found
        """
        dataset_info = await self._data_manager.get_dataset(data_id)
        return dataset_info["adata"]

    async def get_dataset_info(self, data_id: str) -> Dict[str, Any]:
        """Get full dataset info dict when metadata is needed.

        Use this only when you need access to metadata beyond adata,
        such as 'name', 'type', 'source_path', etc.
        """
        return await self._data_manager.get_dataset(data_id)

    async def set_adata(self, data_id: str, adata: Any) -> None:
        """Update the AnnData object for a dataset.

        Use this when preprocessing creates a new adata object (e.g., copy,
        subsample, or format conversion). This updates the reference in the
        data manager's store.

        Args:
            data_id: Dataset identifier
            adata: New AnnData object to store

        Raises:
            ValueError: If dataset not found
        """
        await self._data_manager.update_adata(data_id, adata)

    async def add_dataset(
        self,
        data_id: str,
        adata: Any,
        name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Add a new dataset to the data store.

        Use this when creating new datasets (e.g., integration results,
        subset data, or derived datasets).

        Args:
            data_id: Unique identifier for the new dataset
            adata: AnnData object to store
            name: Optional display name for the dataset
            metadata: Optional additional metadata dict

        Raises:
            ValueError: If dataset with same ID already exists
        """
        await self._data_manager.create_dataset(data_id, adata, name, metadata)

    async def info(self, msg: str) -> None:
        """Log info message to MCP context if available."""
        if self._mcp_context:
            await self._mcp_context.info(msg)

    async def warning(self, msg: str) -> None:
        """Log warning message to MCP context if available."""
        if self._mcp_context:
            await self._mcp_context.warning(msg)

    async def error(self, msg: str) -> None:
        """Log error message to MCP context if available."""
        if self._mcp_context:
            await self._mcp_context.error(msg)

    def get_visualization_cache(self) -> Dict[str, Any]:
        """Get the visualization cache dict.

        Returns:
            The visualization cache dictionary. Returns empty dict if not set.
        """
        if self._visualization_cache is None:
            return {}
        return self._visualization_cache

    def set_visualization(self, key: str, value: Any) -> None:
        """Store a visualization in the cache.

        Args:
            key: Cache key for the visualization
            value: Visualization data (bytes, dict, or other)
        """
        if self._visualization_cache is not None:
            self._visualization_cache[key] = value

    def get_visualization(self, key: str) -> Optional[Any]:
        """Get a visualization from the cache.

        Args:
            key: Cache key for the visualization

        Returns:
            Visualization data if found, None otherwise
        """
        if self._visualization_cache is None:
            return None
        return self._visualization_cache.get(key)

    def clear_visualizations(self, prefix: Optional[str] = None) -> int:
        """Clear visualizations from the cache.

        Args:
            prefix: Optional prefix to filter which keys to clear.
                   If None, clears all visualizations.

        Returns:
            Number of visualizations cleared
        """
        if self._visualization_cache is None:
            return 0

        if prefix is None:
            count = len(self._visualization_cache)
            self._visualization_cache.clear()
            return count

        keys_to_remove = [k for k in self._visualization_cache if k.startswith(prefix)]
        for key in keys_to_remove:
            del self._visualization_cache[key]
        return len(keys_to_remove)


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
    # Server instructions for LLM guidance on tool usage
    instructions = """ChatSpatial provides spatial transcriptomics analysis through 60+ integrated methods across 12 analytical categories.

CORE WORKFLOW PATTERN:
1. Always start with load_data() to import spatial transcriptomics data
2. Run preprocess_data() before most analytical tools (required for clustering, spatial analysis, etc.)
3. Use visualize_data() to inspect results after each analysis step

CRITICAL OPERATIONAL CONSTRAINTS:
- Preprocessing creates filtered gene sets for efficiency but preserves raw data in adata.raw
- Cell communication analysis automatically uses adata.raw when available for comprehensive gene coverage
- Species-specific parameters are critical: set species="mouse" or "human" and use appropriate resources (e.g., liana_resource="mouseconsensus" for mouse)
- Reference data for annotation methods (tangram, scanvi) must be PREPROCESSED before use

PLATFORM-SPECIFIC GUIDANCE:
- Spot-based platforms (Visium, Slide-seq): Deconvolution is recommended to infer cell type compositions
- Single-cell platforms (MERFISH, Xenium, CosMx): Skip deconvolution - native single-cell resolution provided
- Visium with histology images: Use SpaGCN for spatial domain identification
- High-resolution data without images: Use STAGATE or GraphST

TOOL RELATIONSHIPS:
- Spatial domain identification → Enables spatial statistics (neighborhood enrichment, co-occurrence)
- Cell type annotation → Required for cell communication analysis
- Deconvolution results → Can be used for downstream spatial statistics
- Integration → Recommended before cross-sample comparative analyses

PARAMETER GUIDANCE:
All tools include comprehensive parameter documentation in their schemas. Refer to tool descriptions for default values, platform-specific optimizations, and method-specific requirements.

For multi-step analyses, preserve data_id across operations to maintain analysis continuity."""

    # Create MCP server with instructions
    mcp = FastMCP(server_name, instructions=instructions)

    # Create data manager if not provided
    if data_manager is None:
        data_manager = DefaultSpatialDataManager()

    # Create adapter
    adapter = SpatialMCPAdapter(mcp, data_manager)

    return mcp, adapter
