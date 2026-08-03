"""
Spatial MCP Adapter for ChatSpatial

This module provides a clean abstraction layer between MCP protocol requirements
and ChatSpatial's spatial analysis functionality.
"""

import itertools
import logging
from collections.abc import AsyncIterator, Callable, Sequence
from contextlib import AsyncExitStack, asynccontextmanager
from dataclasses import dataclass, field
from functools import wraps
from inspect import signature
from typing import TYPE_CHECKING, Any, Optional, ParamSpec, TypeVar

import anyio
from mcp.server import CacheHint, MCPServer
from mcp.server.mcpserver import Context
from pydantic import BaseModel

from . import __version__
from .utils.adata_utils import get_spatial_key, has_tissue_image
from .utils.exceptions import DataNotFoundError, ParameterError

P = ParamSpec("P")
R = TypeVar("R")

if TYPE_CHECKING:
    from .models.data import SpatialPlatform


class SpatialMCPAdapter:
    """Main adapter class that bridges MCP and spatial analysis functionality."""

    def __init__(
        self, mcp_server: MCPServer, data_manager: "DefaultSpatialDataManager"
    ) -> None:
        self.mcp = mcp_server
        self.data_manager = data_manager


class StrictMCPServer(MCPServer):
    """MCP server whose generated top-level JSON models reject extras.

    MCP SDK 2.0.0 generates permissive argument models for function signatures,
    and for scalar output wrappers, even when nested Pydantic models are strict.
    Tightening both models at registration time keeps discovery and runtime
    validation aligned and prevents unknown values from being silently discarded.
    """

    def add_tool(
        self,
        fn: Callable[..., Any],
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().add_tool(fn, *args, **kwargs)

        explicit_name = kwargs.get("name")
        if explicit_name is None and args:
            explicit_name = args[0]
        tool_name = explicit_name or fn.__name__
        tool = self._tool_manager.get_tool(tool_name)
        if tool is None:  # pragma: no cover - defensive SDK compatibility guard
            raise RuntimeError(f"Registered MCP tool {tool_name!r} is unavailable")

        argument_model = tool.fn_metadata.arg_model
        argument_model.model_config["extra"] = "forbid"
        argument_model.model_rebuild(force=True)
        tool.parameters = argument_model.model_json_schema(by_alias=True)

        output_model = tool.fn_metadata.output_model
        if output_model is not None:
            output_model.model_config["extra"] = "forbid"
            output_model.model_rebuild(force=True)
            tool.fn_metadata.output_schema = output_model.model_json_schema()


@dataclass
class _DatasetLockEntry:
    """A dataset lock and the active calls that still reference it."""

    lock: anyio.Lock = field(default_factory=anyio.Lock)
    references: int = 0


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
        self._dataset_locks: dict[str, _DatasetLockEntry] = {}

    def reset(self) -> None:
        """Reset all process-local datasets and synchronization state."""
        self.data_store.clear()
        self._dataset_locks.clear()
        self._id_counter = itertools.count(1)

    @asynccontextmanager
    async def lock_datasets(self, data_ids: Sequence[str]) -> AsyncIterator[None]:
        """Serialize access to datasets in a deadlock-safe order.

        MCP transports may execute tool calls concurrently. AnnData objects are
        mutable, so both readers and writers must hold the same per-dataset lock
        for the complete tool call. Sorting unique IDs gives multi-dataset tools
        a stable acquisition order. Lock entries are reference-counted and
        removed after the last active or waiting call releases them.
        """
        ordered_ids = sorted(set(data_ids))
        entries: list[tuple[str, _DatasetLockEntry]] = []
        try:
            for data_id in ordered_ids:
                entry = self._dataset_locks.get(data_id)
                if entry is None:
                    entry = _DatasetLockEntry()
                    self._dataset_locks[data_id] = entry
                entry.references += 1
                entries.append((data_id, entry))

            async with AsyncExitStack() as stack:
                for _, entry in entries:
                    await stack.enter_async_context(entry.lock)
                yield
        finally:
            for data_id, entry in entries:
                entry.references -= 1
                if entry.references == 0 and self._dataset_locks.get(data_id) is entry:
                    del self._dataset_locks[data_id]

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
        """Load, validate metadata for, and publish a spatial dataset."""
        data_id, dataset_info = await self.stage_dataset(path, data_type, name)
        await self.publish_dataset(data_id, dataset_info)
        return data_id

    async def stage_dataset(
        self, path: str, data_type: "SpatialPlatform", name: Optional[str] = None
    ) -> tuple[str, dict[str, Any]]:
        """Load and normalize a dataset without making it externally visible.

        Tool boundaries can validate their public response against the staged
        metadata before calling :meth:`publish_dataset`. This keeps failed loads
        from leaving an unreachable dataset behind in process-local state.
        """
        from .utils.data_loader import load_spatial_data

        dataset_info = await load_spatial_data(path, data_type, name)
        data_id = self.reserve_dataset_id("data")
        prepared = dict(dataset_info)
        self._materialize_metadata(prepared, data_id)
        return data_id, prepared

    async def publish_dataset(self, data_id: str, dataset_info: dict[str, Any]) -> None:
        """Atomically publish a fully prepared dataset under a reserved ID."""
        if data_id in self.data_store:
            raise ParameterError(f"Dataset {data_id} already exists")

        prepared = dict(dataset_info)
        self._materialize_metadata(prepared, data_id)
        self.data_store[data_id] = prepared

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
        await self.update_adatas({data_id: adata})

    async def update_adatas(self, updates: dict[str, Any]) -> None:
        """Atomically publish multiple completed AnnData candidates.

        Every target is validated and all derived metadata is materialized in
        detached dataset-info dictionaries before the managed store changes.
        This prevents bilateral tools from publishing only one side when the
        other target or its metadata is invalid.
        """
        missing = [data_id for data_id in updates if data_id not in self.data_store]
        if missing:
            if len(missing) == 1:
                raise DataNotFoundError(f"Dataset {missing[0]} not found")
            joined = ", ".join(sorted(missing))
            raise DataNotFoundError(f"Dataset(s) not found: {joined}")

        prepared: dict[str, dict[str, Any]] = {}
        for data_id, adata in updates.items():
            dataset_info = dict(self.data_store[data_id])
            dataset_info["adata"] = adata
            self._materialize_metadata(dataset_info, data_id)
            prepared[data_id] = dataset_info

        self.data_store.update(prepared)

    async def create_dataset(
        self,
        adata: Any,
        prefix: str = "data",
        name: Optional[str] = None,
        metadata: Optional[dict[str, Any]] = None,
        data_id: Optional[str] = None,
    ) -> str:
        """Create a new dataset with a generated or previously reserved ID.

        Generated IDs follow the same ``{prefix}_{counter}`` pattern as
        ``load_dataset`` to keep identifiers short, unique, and collision-free.

        Args:
            adata: AnnData object to store
            prefix: ID prefix (e.g. ``"integrated"``, ``"subset"``)
            name: Optional display name for the dataset
            metadata: Optional additional metadata dict (reserved keys
                ``adata`` and ``name`` are silently dropped)
            data_id: Previously reserved ID. None allocates a new ID.

        Returns:
            The generated dataset ID
        """
        data_id = data_id or self.reserve_dataset_id(prefix)
        if data_id in self.data_store:
            raise ParameterError(f"Dataset {data_id} already exists")
        dataset_info: dict[str, Any] = {
            "adata": adata,
            "name": name or f"Dataset {data_id}",
            "type": "unknown",
        }
        if metadata:
            _RESERVED_KEYS = {"adata", "name"}
            safe = {k: v for k, v in metadata.items() if k not in _RESERVED_KEYS}
            dataset_info.update(safe)
        await self.publish_dataset(data_id, dataset_info)
        return data_id

    def reserve_dataset_id(self, prefix: str = "data") -> str:
        """Allocate a unique process-local ID before publishing a dataset."""
        while True:
            data_id = f"{prefix}_{next(self._id_counter)}"
            if data_id not in self.data_store:
                return data_id


@dataclass
class ToolContext:
    """Unified context for ChatSpatial tool execution.

    This class provides a clean interface for tools to access data and status
    without the redundant data_store dict wrapping pattern.

    Design Rationale:
    - Python dict assignment is reference, not copy. The old pattern of wrapping
      dataset_info in a temp dict and "writing back" was completely unnecessary.
    - Tools should access adata directly via get_adata(), not through dict wrapping.
    - Status methods fall back gracefully when MCP context is unavailable.

    Status strategy:
    - Progress messages use MCP progress notifications when the client requests them.
    - All messages also use standard Python logging for stderr/OpenTelemetry capture.
    - Non-fatal warnings are retained and attached to the structured tool result.

    Usage:
        async def my_tool(data_id: str, ctx: ToolContext, params: Params) -> Result:
            adata = await ctx.get_adata(data_id)
            await ctx.info(f"Processing {adata.n_obs} cells")  # Progress status
            ctx.debug(f"Internal state: {some_detail}")  # Developer log only
            # ... analysis logic ...
            return result
    """

    _data_manager: "DefaultSpatialDataManager"
    _mcp_context: Optional[Context] = None
    _logger: Optional[logging.Logger] = field(default=None, repr=False)
    _warnings: list[str] = field(default_factory=list, init=False, repr=False)
    _progress_step: int = field(default=0, init=False, repr=False)

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

    async def set_adatas(self, updates: dict[str, Any]) -> None:
        """Atomically publish candidates for a multi-dataset tool call."""
        await self._data_manager.update_adatas(updates)

    async def add_dataset(
        self,
        adata: Any,
        prefix: str = "data",
        name: Optional[str] = None,
        metadata: Optional[dict[str, Any]] = None,
        data_id: Optional[str] = None,
    ) -> str:
        """Add a new dataset to the data store.

        Use this when creating new datasets (e.g., integration results,
        subset data, or derived datasets).

        Args:
            adata: AnnData object to store
            prefix: ID prefix (e.g. ``"integrated"``, ``"subset"``)
            name: Optional display name for the dataset
            metadata: Optional additional metadata dict
            data_id: Previously reserved ID. None allocates a new ID.

        Returns:
            The generated dataset ID
        """
        return await self._data_manager.create_dataset(
            adata,
            prefix,
            name,
            metadata,
            data_id,
        )

    def reserve_dataset_id(self, prefix: str = "data") -> str:
        """Reserve an ID for a candidate that is not ready to publish yet."""
        return self._data_manager.reserve_dataset_id(prefix)

    async def info(self, msg: str) -> None:
        """Record an informational status and report request progress."""
        if self._logger:
            self._logger.info(msg)
        await self._report_status(msg)

    async def warning(self, msg: str) -> None:
        """Record a non-fatal warning and report request progress."""
        if msg not in self._warnings:
            self._warnings.append(msg)
        if self._logger:
            self._logger.warning(msg)
        await self._report_status(f"Warning: {msg}")

    async def error(self, msg: str) -> None:
        """Record an error status before the tool raises its exception."""
        if self._logger:
            self._logger.error(msg)
        await self._report_status(f"Error: {msg}")

    async def _report_status(self, msg: str) -> None:
        """Best-effort monotonic progress reporting.

        Progress is auxiliary transport state, not part of the scientific
        operation. A client that cannot accept a progress notification must not
        turn an otherwise valid analysis into a failed or partially committed
        tool call. Cancellation still propagates because cancellation exceptions
        do not derive from ``Exception``.
        """
        if self._mcp_context is None:
            return
        self._progress_step += 1
        try:
            await self._mcp_context.report_progress(
                float(self._progress_step),
                message=msg,
            )
        except Exception:
            if self._logger:
                self._logger.warning(
                    "Unable to deliver MCP progress update",
                    exc_info=True,
                )

    @property
    def warnings(self) -> tuple[str, ...]:
        """Return immutable non-fatal warnings collected during the tool call."""
        return tuple(self._warnings)

    def finalize(self, result: R) -> R:
        """Attach collected warnings to a tool result without mutating it."""
        if not self._warnings:
            return result

        warnings = list(self._warnings)
        if (
            isinstance(result, BaseModel)
            and "warnings" in result.__class__.model_fields
        ):
            existing = list(getattr(result, "warnings", []))
            merged = list(dict.fromkeys([*existing, *warnings]))
            return result.model_copy(update={"warnings": merged})
        if isinstance(result, dict):
            finalized = dict(result)
            existing = finalized.get("warnings", [])
            if not isinstance(existing, list):
                existing = [str(existing)]
            finalized["warnings"] = list(dict.fromkeys([*existing, *warnings]))
            return finalized  # type: ignore[return-value]
        if isinstance(result, str):
            warning_text = "\n".join(f"- {warning}" for warning in warnings)
            return f"{result}\n\nWarnings:\n{warning_text}"  # type: ignore[return-value]
        return result


def serialize_dataset_access(
    data_manager: DefaultSpatialDataManager,
    *parameter_paths: str,
) -> Callable[[Callable[P, Any]], Callable[P, Any]]:
    """Serialize a tool call across dataset IDs resolved from parameter paths."""

    def decorator(func: Callable[P, Any]) -> Callable[P, Any]:
        func_signature = signature(func)

        @wraps(func)
        async def wrapper(*args: P.args, **kwargs: P.kwargs) -> Any:
            bound = func_signature.bind(*args, **kwargs)
            bound.apply_defaults()
            data_ids: list[str] = []
            for parameter_path in parameter_paths:
                root_name, *attributes = parameter_path.split(".")
                value = bound.arguments[root_name]
                for attribute in attributes:
                    if value is None:
                        break
                    value = getattr(value, attribute)
                if value is None:
                    continue
                if isinstance(value, str):
                    data_ids.append(value)
                elif isinstance(value, Sequence) and all(
                    isinstance(item, str) for item in value
                ):
                    data_ids.extend(value)
                else:
                    raise TypeError(
                        f"{parameter_path} must be a dataset ID or a sequence of dataset IDs"
                    )
            async with data_manager.lock_datasets(data_ids):
                return await func(*args, **kwargs)

        return wrapper

    return decorator


def create_spatial_mcp_server(
    server_name: str = "ChatSpatial",
    data_manager: Optional[DefaultSpatialDataManager] = None,
    *,
    server_version: str = __version__,
) -> tuple[MCPServer, SpatialMCPAdapter]:
    """
    Create and configure a spatial MCP server with adapter

    Args:
        server_name: Name of the MCP server
        data_manager: Optional custom data manager (uses default if None)
        server_version: ChatSpatial package version reported to MCP clients

    Returns:
        Tuple of (MCPServer instance, SpatialMCPAdapter instance)
    """
    # Server instructions for LLM guidance on tool usage
    instructions = """ChatSpatial provides spatial transcriptomics analysis through 66 methods across 15 analytical categories.

CORE WORKFLOW PATTERN:
1. Always start with load_data() to import spatial transcriptomics data
2. Run preprocess_data() for QC, normalization, and highly variable gene selection
3. Run compute_embeddings() before tools that require PCA, UMAP, clusters, or neighbor graphs
4. Run the desired biological or spatial analysis tools
5. Use visualize_data() to inspect results after each analysis step

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

    # Tool definitions are process-global and identical for every caller, so a
    # public cache hint is safe and improves client prompt-cache stability.
    mcp = StrictMCPServer(
        server_name,
        title="ChatSpatial",
        description=(
            "Schema-enforced spatial transcriptomics analysis with 20 tools "
            "covering 66 methods across 15 analytical categories."
        ),
        instructions=instructions,
        website_url="https://github.com/cafferychen777/ChatSpatial",
        version=server_version,
        cache_hints={
            "tools/list": CacheHint(ttl_ms=3_600_000, scope="public"),
        },
    )

    # Create data manager if not provided
    if data_manager is None:
        data_manager = DefaultSpatialDataManager()

    # Create adapter
    adapter = SpatialMCPAdapter(mcp, data_manager)

    return mcp, adapter
