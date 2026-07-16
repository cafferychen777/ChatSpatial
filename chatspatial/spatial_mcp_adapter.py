"""
Spatial MCP Adapter for ChatSpatial

This module provides a clean abstraction layer between MCP protocol requirements
and ChatSpatial's spatial analysis functionality.
"""

import itertools
import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Optional

from mcp.server.fastmcp import Context, FastMCP

from .utils.adata_utils import get_spatial_key, has_tissue_image
from .utils.exceptions import DataNotFoundError

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .models.data import SpatialPlatform


class SpatialMCPAdapter:
    """Main adapter class that bridges MCP and spatial analysis functionality."""

    def __init__(
        self, mcp_server: FastMCP, data_manager: "DefaultSpatialDataManager"
    ) -> None:
        self.mcp = mcp_server
        self.data_manager = data_manager


class DefaultSpatialDataManager:
    """In-memory spatial data management with async interface.

    Design Note:
        Methods are async for interface consistency and future extensibility
        (e.g., remote storage, database backends), even though current
        implementation is synchronous. This is intentional - async overhead
        is negligible and changing the interface later would break 20+ call sites.
    """

    def __init__(self) -> None:
        self.data_store: dict[str, Any] = {}
        self._id_counter = itertools.count(1)

    @staticmethod
    def _extract_adata_metadata(adata: Any) -> dict[str, Any]:
        """Extract lightweight metadata from an AnnData object.

        Pure function: reads adata attributes, returns a new dict.
        Does not modify the input.
        """
        meta: dict[str, Any] = {}

        n_obs = getattr(adata, "n_obs", None)
        n_vars = getattr(adata, "n_vars", None)
        meta["n_cells"] = int(n_obs) if isinstance(n_obs, int) else 0
        meta["n_genes"] = int(n_vars) if isinstance(n_vars, int) else 0

        obsm = getattr(adata, "obsm", None)
        if obsm is not None:
            meta["obsm_keys"] = list(obsm.keys())
            spatial_key = get_spatial_key(adata)
            has_spatial = spatial_key is not None and len(obsm[spatial_key]) > 0
            meta["spatial_coordinates_available"] = bool(has_spatial)
        else:
            meta["obsm_keys"] = []
            meta["spatial_coordinates_available"] = False

        uns = getattr(adata, "uns", None)
        if uns is not None:
            meta["uns_keys"] = list(uns.keys())
            meta["tissue_image_available"] = has_tissue_image(adata)
        else:
            meta["uns_keys"] = []
            meta["tissue_image_available"] = False

        return meta

    def _materialize_metadata(self, dataset_info: dict[str, Any], data_id: str) -> None:
        """Write adata-derived metadata into dataset_info.

        Call this on write paths only (load, create, update_adata).
        """
        dataset_info.setdefault("name", f"Dataset {data_id}")
        dataset_info.setdefault("type", "unknown")
        adata = dataset_info.get("adata")
        if adata is not None:
            dataset_info.update(self._extract_adata_metadata(adata))

    async def load_dataset(
        self, path: str, data_type: "SpatialPlatform", name: Optional[str] = None
    ) -> str:
        """Load a spatial dataset and return its ID"""
        from .utils.data_loader import load_spatial_data

        dataset_info = await load_spatial_data(path, data_type, name)

        # Generate ID — single expression, structurally atomic
        data_id = f"data_{next(self._id_counter)}"

        # Store data
        self.data_store[data_id] = dataset_info
        self._materialize_metadata(dataset_info, data_id)

        return data_id

    async def get_dataset(self, data_id: str) -> Any:
        """Get a dataset by ID."""
        if data_id not in self.data_store:
            raise DataNotFoundError(f"Dataset {data_id} not found")
        return self.data_store[data_id]

    async def list_datasets(self) -> list[dict[str, Any]]:
        """List all loaded datasets with current cell/gene counts."""
        listed: list[dict[str, Any]] = []
        for data_id, info in self.data_store.items():
            adata = info.get("adata")
            n_obs = getattr(adata, "n_obs", None)
            n_vars = getattr(adata, "n_vars", None)
            listed.append(
                {
                    "id": data_id,
                    "name": info.get("name", f"Dataset {data_id}"),
                    "type": info.get("type", "unknown"),
                    "n_cells": n_obs if isinstance(n_obs, int) else 0,
                    "n_genes": n_vars if isinstance(n_vars, int) else 0,
                }
            )
        return listed

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
        self._materialize_metadata(self.data_store[data_id], data_id)

    async def create_dataset(
        self,
        adata: Any,
        prefix: str = "data",
        name: Optional[str] = None,
        metadata: Optional[dict[str, Any]] = None,
    ) -> str:
        """Create a new dataset with an auto-generated unique ID.

        ID generation follows the same ``{prefix}_{counter}`` pattern as
        ``load_dataset`` to keep identifiers short, unique, and collision-free.

        Args:
            adata: AnnData object to store
            prefix: ID prefix (e.g. ``"integrated"``, ``"subset"``)
            name: Optional display name for the dataset
            metadata: Optional additional metadata dict (reserved keys
                ``adata``, ``name``, ``results`` are silently dropped)

        Returns:
            The generated dataset ID
        """
        data_id = f"{prefix}_{next(self._id_counter)}"
        dataset_info: dict[str, Any] = {
            "adata": adata,
            "name": name or f"Dataset {data_id}",
            "type": "unknown",
        }
        if metadata:
            _RESERVED_KEYS = {"adata", "name", "results"}
            safe = {k: v for k, v in metadata.items() if k not in _RESERVED_KEYS}
            dataset_info.update(safe)
        self._materialize_metadata(dataset_info, data_id)
        self.data_store[data_id] = dataset_info
        return data_id


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

    def log_config(self, title: str, config: dict[str, Any]) -> None:
        """Log configuration details for developers.

        Convenience method for logging parameter configurations in a
        structured format. Goes to Python logger, not user-visible.

        Args:
            title: Configuration section title
            config: Dictionary of configuration key-value pairs
        """
        if self._logger:
            self._logger.debug("=" * 50)
            self._logger.debug(f"{title}:")
            for key, value in config.items():
                self._logger.debug(f"  {key}: {value}")
            self._logger.debug("=" * 50)

    async def get_adata(self, data_id: str) -> Any:
        """Get AnnData object directly by ID.

        This is the primary data access method for tools. Returns the AnnData
        object directly without intermediate dict wrapping.

        Args:
            data_id: Dataset identifier

        Returns:
            AnnData object for the dataset

        Raises:
            DataNotFoundError: If dataset not found
        """
        dataset_info = await self._data_manager.get_dataset(data_id)
        return dataset_info["adata"]

    async def get_dataset_info(self, data_id: str) -> dict[str, Any]:
        """Get full dataset info dict when metadata is needed.

        Use this only when you need access to metadata beyond adata,
        such as 'name', 'type', 'source_path', etc.
        """
        return await self._data_manager.get_dataset(data_id)

    async def set_adata(self, data_id: str, adata: Any) -> None:
        """Update the AnnData object for a dataset.

        Use this to publish a completed candidate after a mutating tool has
        finished computation, validation, metadata, and export. This replaces
        the reference in the data manager's store in one ownership handoff.

        Args:
            data_id: Dataset identifier
            adata: New AnnData object to store

        Raises:
            DataNotFoundError: If dataset not found
        """
        await self._data_manager.update_adata(data_id, adata)

    async def add_dataset(
        self,
        adata: Any,
        prefix: str = "data",
        name: Optional[str] = None,
        metadata: Optional[dict[str, Any]] = None,
    ) -> str:
        """Add a new dataset to the data store.

        Use this when creating new datasets (e.g., integration results,
        subset data, or derived datasets).

        Args:
            adata: AnnData object to store
            prefix: ID prefix (e.g. ``"integrated"``, ``"subset"``)
            name: Optional display name for the dataset
            metadata: Optional additional metadata dict

        Returns:
            The generated dataset ID
        """
        return await self._data_manager.create_dataset(adata, prefix, name, metadata)

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
    instructions = """ChatSpatial provides spatial transcriptomics analysis through 65 methods across 15 analytical categories.

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
