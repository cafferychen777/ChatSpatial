"""
Spatial MCP Adapter for ChatSpatial

This module provides a clean abstraction layer between MCP protocol requirements
and ChatSpatial's spatial analysis functionality.

Design Principles:
- Single source of truth: One registry for visualization state
- Store params, not bytes: Images can be regenerated on demand
- LRU eviction: Automatic memory management
"""

import logging
import os
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from mcp.server.fastmcp import Context, FastMCP
from mcp.types import ToolAnnotations

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


# =============================================================================
# VISUALIZATION REGISTRY - Single Source of Truth
# =============================================================================
# First Principles Design:
# - Store params, not bytes: Images can be regenerated on demand from params
# - LRU eviction: Automatic memory management prevents unbounded growth
# - Unified cache key: {data_id}_{plot_type}[_{subtype}]
# =============================================================================


@dataclass
class VisualizationEntry:
    """A single visualization entry in the registry.

    Stores the parameters needed to regenerate a visualization, not the image bytes.
    This follows the first principle: data + params = reproducible output.
    """

    params: VisualizationParameters
    file_path: Optional[str] = None  # Path to saved PNG (if large image)
    timestamp: float = field(default_factory=time.time)


class VisualizationRegistry:
    """Single source of truth for visualization state.

    Design Principles:
    - Store params, not bytes: Images regenerated on demand
    - LRU eviction: Oldest entries removed when max_entries exceeded
    - Unified interface: One place for all visualization state
    """

    def __init__(self, max_entries: int = 100):
        """Initialize the registry with LRU capacity.

        Args:
            max_entries: Maximum number of entries before LRU eviction
        """
        self._entries: OrderedDict[str, VisualizationEntry] = OrderedDict()
        self._max_entries = max_entries

    def store(
        self,
        key: str,
        params: VisualizationParameters,
        file_path: Optional[str] = None,
    ) -> None:
        """Store a visualization entry.

        Args:
            key: Cache key format: {data_id}_{plot_type}[_{subtype}]
            params: Visualization parameters for regeneration
            file_path: Optional path to saved image file
        """
        # Move to end if exists (LRU behavior)
        if key in self._entries:
            self._entries.move_to_end(key)

        self._entries[key] = VisualizationEntry(
            params=params,
            file_path=file_path,
            timestamp=time.time(),
        )

        # LRU eviction
        while len(self._entries) > self._max_entries:
            oldest_key, oldest_entry = self._entries.popitem(last=False)
            # Clean up file if exists
            if oldest_entry.file_path and os.path.exists(oldest_entry.file_path):
                try:
                    os.remove(oldest_entry.file_path)
                except OSError:
                    pass

    def get(self, key: str) -> Optional[VisualizationEntry]:
        """Get a visualization entry.

        Args:
            key: Cache key to look up

        Returns:
            VisualizationEntry if found, None otherwise
        """
        if key in self._entries:
            self._entries.move_to_end(key)  # LRU touch
            return self._entries[key]
        return None

    def exists(self, key: str) -> bool:
        """Check if a visualization exists in registry."""
        return key in self._entries

    def list_for_dataset(self, data_id: str) -> List[str]:
        """List all visualization keys for a dataset.

        Args:
            data_id: Dataset identifier

        Returns:
            List of cache keys matching the dataset
        """
        prefix = f"{data_id}_"
        return [k for k in self._entries.keys() if k.startswith(prefix)]

    def clear(self, prefix: Optional[str] = None) -> int:
        """Clear visualizations from registry.

        Args:
            prefix: Optional prefix to filter which keys to clear.
                   If None, clears all visualizations.

        Returns:
            Number of entries cleared
        """
        if prefix is None:
            count = len(self._entries)
            # Clean up all files
            for entry in self._entries.values():
                if entry.file_path and os.path.exists(entry.file_path):
                    try:
                        os.remove(entry.file_path)
                    except OSError:
                        pass
            self._entries.clear()
            return count

        keys_to_remove = [k for k in self._entries if k.startswith(prefix)]
        for key in keys_to_remove:
            entry = self._entries.pop(key)
            if entry.file_path and os.path.exists(entry.file_path):
                try:
                    os.remove(entry.file_path)
                except OSError:
                    pass
        return len(keys_to_remove)

    def keys(self) -> List[str]:
        """Return all keys in the registry."""
        return list(self._entries.keys())


class SpatialMCPAdapter:
    """Main adapter class that bridges MCP and spatial analysis functionality.

    Simplified design: Only manages data and visualization registry.
    Removed dead code (ResourceManager was never registered to MCP server).
    """

    def __init__(self, mcp_server: FastMCP, data_manager: "DefaultSpatialDataManager"):
        self.mcp = mcp_server
        self.data_manager = data_manager
        self.visualization_registry = VisualizationRegistry(max_entries=100)


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
    _visualization_registry: Optional["VisualizationRegistry"] = None
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

    def get_visualization_registry(self) -> Optional["VisualizationRegistry"]:
        """Get the visualization registry.

        Returns:
            VisualizationRegistry instance if set, None otherwise
        """
        return self._visualization_registry

    def store_visualization(
        self,
        key: str,
        params: VisualizationParameters,
        file_path: Optional[str] = None,
    ) -> None:
        """Store a visualization in the registry.

        Args:
            key: Cache key format: {data_id}_{plot_type}[_{subtype}]
            params: Visualization parameters for regeneration
            file_path: Optional path to saved image file
        """
        if self._visualization_registry is not None:
            self._visualization_registry.store(key, params, file_path)

    def get_visualization(self, key: str) -> Optional[VisualizationEntry]:
        """Get a visualization entry from the registry.

        Args:
            key: Cache key for the visualization

        Returns:
            VisualizationEntry if found, None otherwise
        """
        if self._visualization_registry is None:
            return None
        return self._visualization_registry.get(key)

    def visualization_exists(self, key: str) -> bool:
        """Check if a visualization exists in registry.

        Args:
            key: Cache key to check

        Returns:
            True if exists, False otherwise
        """
        if self._visualization_registry is None:
            return False
        return self._visualization_registry.exists(key)

    def list_visualizations(self, data_id: str) -> List[str]:
        """List all visualization keys for a dataset.

        Args:
            data_id: Dataset identifier

        Returns:
            List of cache keys matching the dataset
        """
        if self._visualization_registry is None:
            return []
        return self._visualization_registry.list_for_dataset(data_id)

    def clear_visualizations(self, prefix: Optional[str] = None) -> int:
        """Clear visualizations from the registry.

        Args:
            prefix: Optional prefix to filter which keys to clear.
                   If None, clears all visualizations.

        Returns:
            Number of visualizations cleared
        """
        if self._visualization_registry is None:
            return 0
        return self._visualization_registry.clear(prefix)


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
    instructions = """ChatSpatial provides spatial transcriptomics analysis through 60+ integrated methods across 15 analytical categories.

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
