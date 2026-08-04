"""Virtual spatial transcriptomics from pre-cut H&E tiles.

This module turns a directory of pre-cut 224x224 H&E tiles plus a coordinate
manifest into a spatial AnnData whose ``X`` holds DeepSpot-M predicted
log1p-CPM expression. It is a data entry path, not a backend for an existing
analysis tool: nothing here reads an existing dataset.

Scope is deliberately narrow. One slide and one tile directory per invocation,
tiles that the caller has already cut, and no persistence. Whole-slide reading
and tiling stay outside ChatSpatial.

Manifest v1
-----------
A CSV with the required columns ``tile_id``, ``tile_path``, ``slide_id``,
``x_px``, ``y_px``, ``mpp_x`` and ``mpp_y``. ``slide_id``, ``mpp_x`` and
``mpp_y`` are manifest-wide values repeated on every row; the repetition is
intentional so the manifest stays a portable single file, and the values must
be consistent across rows. Because the v1 schema fixes the coordinate
convention, no companion metadata file is needed.

``x_px`` and ``y_px`` are the upper-left corner of the tile in the native
level-0 slide frame: top-left origin, x increasing rightward, y increasing
downward. ``obsm["spatial"]`` holds tile centers, ``(x_px + width / 2,
y_px + height / 2)``, which is the convention scanpy and squidpy expect.

``grid_row`` and ``grid_col`` are optional, must appear together, and are
published only when a single affine mapping relates them to the pixel
coordinates. Irregular manifests stay pixel-only. These are generic lattice
indices, not Visium ``array_row``/``array_col``.

Licensing
---------
DeepSpot-M code is PolyForm Noncommercial 1.0.0 and its published weights are
CC-BY-NC-SA-4.0 and gated on Hugging Face. The package is an optional extra
that ChatSpatial imports lazily, so installing ChatSpatial attaches neither
license to it. Output is for research use only.
"""

from __future__ import annotations

import asyncio
import concurrent.futures
import csv
import hashlib
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, Sequence

import numpy as np
import pandas as pd

from ..models.analysis import VirtualExpressionResult
from ..models.data import HistologyExpressionParameters
from ..utils.dependency_manager import require
from ..utils.device_utils import get_device
from ..utils.exceptions import (
    ChatSpatialError,
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)
from ..utils.provenance import (
    PROVENANCE_PREDICTED,
    UNITS_LOG1P_CPM,
    set_expression_provenance,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..spatial_mcp_adapter import ToolContext

# =============================================================================
# Contract constants
# =============================================================================

#: Version of the CSV manifest schema this module accepts.
MANIFEST_SCHEMA_VERSION = 1

#: Columns every v1 manifest must provide.
REQUIRED_MANIFEST_COLUMNS: tuple[str, ...] = (
    "tile_id",
    "tile_path",
    "slide_id",
    "x_px",
    "y_px",
    "mpp_x",
    "mpp_y",
)

#: Optional lattice columns. Both or neither.
LATTICE_COLUMNS: tuple[str, ...] = ("grid_row", "grid_col")

#: DeepSpot-M consumes 224x224 RGB tiles at roughly 20x (about 0.5 mpp).
TILE_WIDTH_PX = 224
TILE_HEIGHT_PX = 224

#: Written to ``uns`` so downstream readers do not have to guess the frame.
COORDINATE_CONVENTION = (
    "tile upper-left corner in the native level-0 slide frame, top-left "
    "origin, x increases rightward, y increases downward"
)

#: Identifier used in the generic expression metadata block.
PRODUCER = "deepspotm"

#: Licensing and use metadata, published on every generated dataset.
CODE_LICENSE = "PolyForm-Noncommercial-1.0.0"
WEIGHTS_LICENSE = "CC-BY-NC-SA-4.0"
RESEARCH_USE_NOTICE = "Research use only. Not for clinical or diagnostic use."
REFERENCE_DOI = "https://doi.org/10.64898/2026.06.19.26356060"

#: Default checkpoint reference resolved when the caller pins nothing.
DEFAULT_REVISION = "main"

#: Files that define one local DeepSpot-M checkpoint.
LOCAL_CHECKPOINT_FILES: tuple[str, ...] = (
    "config.json",
    "model.safetensors",
    "tokens.csv",
)

#: How many offending rows an error message lists before summarizing.
_MAX_REPORTED = 5


def _summarize(items: Sequence[str]) -> str:
    """Render a bounded, deterministic sample of offending manifest entries."""
    listed = list(items)
    shown = ", ".join(listed[:_MAX_REPORTED])
    if len(listed) > _MAX_REPORTED:
        return f"{shown} (and {len(listed) - _MAX_REPORTED} more)"
    return shown


# =============================================================================
# Manifest
# =============================================================================


@dataclass(frozen=True)
class TileManifest:
    """A validated v1 manifest bound to one slide and one tile directory."""

    manifest_path: Path
    tile_directory: Path
    digest: str
    tile_ids: list[str]
    tile_paths: list[str]
    resolved_paths: list[Path]
    x_px: np.ndarray
    y_px: np.ndarray
    slide_id: str
    mpp_x: float
    mpp_y: float
    grid_row: Optional[np.ndarray] = None
    grid_col: Optional[np.ndarray] = None
    row_stride_px: Optional[int] = None
    col_stride_px: Optional[int] = None

    @property
    def n_tiles(self) -> int:
        return len(self.tile_ids)

    @property
    def has_lattice(self) -> bool:
        return self.grid_row is not None and self.grid_col is not None

    def centers(self) -> np.ndarray:
        """Return tile-center coordinates in the level-0 slide frame."""
        centers = np.empty((self.n_tiles, 2), dtype=np.float64)
        centers[:, 0] = self.x_px + TILE_WIDTH_PX / 2.0
        centers[:, 1] = self.y_px + TILE_HEIGHT_PX / 2.0
        return centers


def _digest_file(path: Path) -> str:
    """Return the SHA-256 digest of a file, so a manifest is identifiable."""
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def _read_header(manifest_path: Path) -> list[str]:
    """Return the raw header row, before pandas can rename duplicates."""
    with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
        try:
            header = next(csv.reader(handle))
        except StopIteration as exc:
            raise DataError(
                f"Tile manifest '{manifest_path}' is empty. A v1 manifest needs "
                f"a header row with the columns "
                f"{', '.join(REQUIRED_MANIFEST_COLUMNS)} and at least one tile row."
            ) from exc
    return [name.strip() for name in header]


def _validate_header(manifest_path: Path, header: Sequence[str]) -> None:
    """Reject missing required columns, duplicates, and half-declared lattices."""
    duplicates = sorted({name for name in header if list(header).count(name) > 1})
    if duplicates:
        raise DataError(
            f"Tile manifest '{manifest_path}' declares duplicate columns: "
            f"{', '.join(duplicates)}. Each column may appear only once."
        )

    missing = [name for name in REQUIRED_MANIFEST_COLUMNS if name not in header]
    if missing:
        raise DataError(
            f"Tile manifest '{manifest_path}' is missing required v1 columns: "
            f"{', '.join(missing)}. A v1 manifest requires "
            f"{', '.join(REQUIRED_MANIFEST_COLUMNS)}, with slide_id, mpp_x and "
            f"mpp_y repeated identically on every row."
        )

    present_lattice = [name for name in LATTICE_COLUMNS if name in header]
    if len(present_lattice) == 1:
        raise DataError(
            f"Tile manifest '{manifest_path}' declares '{present_lattice[0]}' "
            f"without the matching column. grid_row and grid_col are optional, "
            f"but they must appear together. Omit both for an irregular "
            f"manifest; the pixel coordinates remain the source of truth."
        )


def _read_frame(manifest_path: Path) -> pd.DataFrame:
    """Read the manifest as raw text so every conversion is explicit."""
    try:
        frame = pd.read_csv(
            manifest_path,
            dtype=str,
            keep_default_na=False,
            na_filter=False,
            encoding="utf-8-sig",
        )
    except pd.errors.EmptyDataError as exc:
        raise DataError(
            f"Tile manifest '{manifest_path}' is empty. A v1 manifest needs a "
            f"header row and at least one tile row."
        ) from exc
    except (OSError, UnicodeDecodeError, pd.errors.ParserError) as exc:
        raise DataError(
            f"Tile manifest '{manifest_path}' could not be read as UTF-8 CSV: {exc}"
        ) from exc

    if frame.empty:
        raise DataError(
            f"Tile manifest '{manifest_path}' has a header but no tile rows. "
            f"At least two tiles are required because a single observation "
            f"cannot support spatial-neighborhood analysis."
        )
    return frame


def _column_values(frame: pd.DataFrame, column: str) -> list[str]:
    """Return a stripped string column, preserving manifest row order."""
    return [str(value).strip() for value in frame[column].tolist()]


def _require_non_empty(manifest_path: Path, column: str, values: Sequence[str]) -> None:
    blank = [str(index + 2) for index, value in enumerate(values) if not value]
    if blank:
        raise DataError(
            f"Tile manifest '{manifest_path}' has blank '{column}' values at "
            f"CSV line(s) {_summarize(blank)}. Every row must provide a value."
        )


def _single_value(manifest_path: Path, column: str, values: Sequence[str]) -> str:
    """Return the one manifest-wide value of a column, or reject the manifest."""
    distinct = sorted(set(values))
    if len(distinct) != 1:
        raise DataError(
            f"Tile manifest '{manifest_path}' mixes {len(distinct)} '{column}' "
            f"values: {_summarize(distinct)}. This column is manifest-wide and "
            f"must repeat the same value on every row. Split the tiles into one "
            f"manifest per value and run this tool once per slide."
        )
    return distinct[0]


#: Trailing guidance appended to coordinate and lattice parse errors.
_PIXEL_GUIDANCE = (
    "Pixel coordinates are the tile's upper-left corner in the level-0 slide "
    "frame, measured from the top-left origin, so they are whole numbers and "
    "never negative."
)
_LATTICE_GUIDANCE = (
    "Lattice indices are non-negative whole numbers on a regular axis-aligned "
    "grid. Omit grid_row and grid_col to keep the manifest pixel-only."
)


def _parse_integer_column(
    manifest_path: Path,
    column: str,
    values: Sequence[str],
    guidance: str = _PIXEL_GUIDANCE,
) -> np.ndarray:
    """Parse a column as finite, non-negative whole numbers."""
    parsed = np.empty(len(values), dtype=np.int64)
    for index, raw in enumerate(values):
        line = index + 2
        try:
            number = float(raw)
        except ValueError as exc:
            raise DataError(
                f"Tile manifest '{manifest_path}' has a non-numeric '{column}' "
                f"value {raw!r} at CSV line {line}. {guidance}"
            ) from exc
        if not np.isfinite(number):
            raise DataError(
                f"Tile manifest '{manifest_path}' has a non-finite '{column}' "
                f"value {raw!r} at CSV line {line}. {guidance}"
            )
        if number < 0:
            raise DataError(
                f"Tile manifest '{manifest_path}' has a negative '{column}' "
                f"value {raw!r} at CSV line {line}. {guidance}"
            )
        if number != int(number):
            raise DataError(
                f"Tile manifest '{manifest_path}' has a fractional '{column}' "
                f"value {raw!r} at CSV line {line}. {guidance}"
            )
        parsed[index] = int(number)
    return parsed


def _parse_manifest_mpp(
    manifest_path: Path, column: str, values: Sequence[str]
) -> float:
    """Parse and enforce one numeric microns-per-pixel value."""
    parsed: list[float] = []
    for index, raw in enumerate(values):
        try:
            number = float(raw)
        except ValueError as exc:
            raise DataError(
                f"Tile manifest '{manifest_path}' has a non-numeric '{column}' "
                f"value {raw!r} at CSV line {index + 2}. Declare the slide's "
                f"microns per pixel at level 0."
            ) from exc
        if not np.isfinite(number) or number <= 0:
            raise DataError(
                f"Tile manifest '{manifest_path}' has an invalid '{column}' value "
                f"{raw!r} at CSV line {index + 2}. Microns per pixel must be "
                f"finite and greater than zero."
            )
        parsed.append(number)

    distinct = sorted(set(parsed))
    if len(distinct) != 1:
        rendered = [str(value) for value in distinct]
        raise DataError(
            f"Tile manifest '{manifest_path}' mixes {len(distinct)} numeric "
            f"'{column}' values: {_summarize(rendered)}. This column is "
            f"manifest-wide and must repeat the same value on every row."
        )
    return distinct[0]


def _reject_duplicate_identifiers(manifest_path: Path, tile_ids: Sequence[str]) -> None:
    seen: dict[str, int] = {}
    duplicates: list[str] = []
    for index, tile_id in enumerate(tile_ids):
        if tile_id in seen:
            duplicates.append(
                f"{tile_id!r} (CSV lines {seen[tile_id] + 2} and {index + 2})"
            )
        else:
            seen[tile_id] = index
    if duplicates:
        raise DataError(
            f"Tile manifest '{manifest_path}' repeats tile identifiers: "
            f"{_summarize(duplicates)}. tile_id must be unique because it becomes "
            f"the observation name."
        )


def _reject_duplicate_coordinates(
    manifest_path: Path, x_px: np.ndarray, y_px: np.ndarray
) -> None:
    seen: dict[tuple[int, int], int] = {}
    duplicates: list[str] = []
    for index, (x, y) in enumerate(zip(x_px.tolist(), y_px.tolist())):
        key = (int(x), int(y))
        if key in seen:
            duplicates.append(
                f"({key[0]}, {key[1]}) at CSV lines {seen[key] + 2} and {index + 2}"
            )
        else:
            seen[key] = index
    if duplicates:
        raise DataError(
            f"Tile manifest '{manifest_path}' places two tiles at the same "
            f"upper-left corner: {_summarize(duplicates)}. Each tile needs a "
            f"distinct position in the level-0 slide frame."
        )


def _resolve_tile_paths(
    manifest_path: Path, tile_directory: Path, relative_paths: Sequence[str]
) -> list[Path]:
    """Resolve each tile_path inside the declared directory, rejecting escapes."""
    absolute = [value for value in relative_paths if Path(value).is_absolute()]
    if absolute:
        raise DataError(
            f"Tile manifest '{manifest_path}' contains absolute tile_path values: "
            f"{_summarize(absolute)}. tile_path resolves relative to the declared "
            f"tile directory '{tile_directory}'."
        )

    resolved: list[Path] = []
    escaping: list[str] = []
    for value in relative_paths:
        candidate = (tile_directory / value).resolve()
        if not candidate.is_relative_to(tile_directory):
            escaping.append(value)
        resolved.append(candidate)

    if escaping:
        raise DataError(
            f"Tile manifest '{manifest_path}' contains tile_path values that "
            f"resolve outside the declared tile directory '{tile_directory}': "
            f"{_summarize(escaping)}. Tiles must live inside that directory."
        )

    missing = [
        str(path.relative_to(tile_directory)) for path in resolved if not path.is_file()
    ]
    if missing:
        raise DataNotFoundError(
            f"Tile manifest '{manifest_path}' references {len(missing)} file(s) "
            f"that do not exist under '{tile_directory}': {_summarize(missing)}."
        )

    seen: dict[Path, int] = {}
    duplicates: list[str] = []
    for index, path in enumerate(resolved):
        if path in seen:
            duplicates.append(
                f"{relative_paths[index]!r} (CSV lines {seen[path] + 2} and "
                f"{index + 2})"
            )
        else:
            seen[path] = index
    if duplicates:
        raise DataError(
            f"Tile manifest '{manifest_path}' references the same decoded tile "
            f"more than once: {_summarize(duplicates)}. Each observation must "
            f"come from a distinct tile file."
        )
    return resolved


def _solve_axis_affine(
    manifest_path: Path,
    grid_name: str,
    pixel_name: str,
    grid: np.ndarray,
    pixel: np.ndarray,
) -> Optional[int]:
    """Return the stride of one exact affine mapping ``pixel = stride * grid + offset``.

    Returns ``None`` when the axis is degenerate, that is when the manifest
    describes a single grid line and the stride is therefore unobservable.
    """
    distinct = np.unique(grid)
    if distinct.size == 1:
        if np.unique(pixel).size != 1:
            raise DataError(
                f"Tile manifest '{manifest_path}' assigns one '{grid_name}' value "
                f"to several different '{pixel_name}' values. A lattice index and "
                f"a pixel coordinate must agree."
            )
        return None

    reference = int(np.argmin(grid))
    grid_reference = int(grid[reference])
    pixel_reference = int(pixel[reference])

    deltas_grid = grid - grid_reference
    candidates = np.flatnonzero(deltas_grid)
    nearest = candidates[int(np.argmin(np.abs(deltas_grid[candidates])))]
    delta_grid = int(deltas_grid[nearest])
    delta_pixel = int(pixel[nearest]) - pixel_reference

    if delta_pixel % delta_grid != 0:
        raise DataError(
            f"Tile manifest '{manifest_path}' has no whole-pixel stride between "
            f"'{grid_name}' and '{pixel_name}'. Declare grid indices only for a "
            f"regular axis-aligned lattice, or omit grid_row and grid_col and "
            f"keep the manifest pixel-only."
        )

    stride = delta_pixel // delta_grid
    if stride <= 0:
        raise DataError(
            f"Tile manifest '{manifest_path}' maps increasing '{grid_name}' onto "
            f"non-increasing '{pixel_name}'. A lattice index must increase with "
            f"the pixel coordinate."
        )

    offset = pixel_reference - stride * grid_reference
    if not np.array_equal(pixel, stride * grid + offset):
        raise DataError(
            f"Tile manifest '{manifest_path}' does not follow one affine mapping "
            f"between '{grid_name}' and '{pixel_name}'. The first consistent "
            f"stride is {stride} px with offset {offset} px, which does not hold "
            f"for every row. Declare grid indices only for a regular "
            f"axis-aligned lattice, or omit grid_row and grid_col and keep the "
            f"manifest pixel-only."
        )
    return stride


def _parse_lattice(
    manifest_path: Path,
    frame: pd.DataFrame,
    x_px: np.ndarray,
    y_px: np.ndarray,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[int], Optional[int]]:
    """Validate optional lattice columns against the pixel coordinates."""
    if not all(name in frame.columns for name in LATTICE_COLUMNS):
        return None, None, None, None

    grid_row = _parse_integer_column(
        manifest_path, "grid_row", _column_values(frame, "grid_row"), _LATTICE_GUIDANCE
    )
    grid_col = _parse_integer_column(
        manifest_path, "grid_col", _column_values(frame, "grid_col"), _LATTICE_GUIDANCE
    )

    pairs = set(zip(grid_row.tolist(), grid_col.tolist()))
    if len(pairs) != grid_row.size:
        raise DataError(
            f"Tile manifest '{manifest_path}' repeats (grid_row, grid_col) pairs. "
            f"Each lattice position may hold only one tile."
        )

    row_stride = _solve_axis_affine(manifest_path, "grid_row", "y_px", grid_row, y_px)
    col_stride = _solve_axis_affine(manifest_path, "grid_col", "x_px", grid_col, x_px)
    return grid_row, grid_col, row_stride, col_stride


def read_tile_manifest(manifest_path: str, tile_directory: str) -> TileManifest:
    """Read and fully validate a v1 tile manifest.

    Every check runs before any model is loaded, so an invalid manifest fails
    fast and cheaply. Manifest row order is preserved throughout.

    Args:
        manifest_path: Path to the CSV manifest.
        tile_directory: Directory that every ``tile_path`` resolves against.

    Returns:
        A validated :class:`TileManifest`.

    Raises:
        ParameterError: If the manifest or tile directory path is unusable.
        DataError: If the manifest violates the v1 schema.
        DataNotFoundError: If a referenced tile file is missing.
    """
    manifest = Path(manifest_path).expanduser()
    directory = Path(tile_directory).expanduser()

    if not manifest.is_file():
        raise ParameterError(
            f"Tile manifest '{manifest}' does not exist or is not a file. Provide "
            f"the path to a v1 CSV manifest with the columns "
            f"{', '.join(REQUIRED_MANIFEST_COLUMNS)}."
        )
    if not directory.is_dir():
        raise ParameterError(
            f"Tile directory '{directory}' does not exist or is not a directory. "
            f"Every manifest tile_path resolves relative to this directory."
        )

    manifest = manifest.resolve()
    directory = directory.resolve()

    header = _read_header(manifest)
    _validate_header(manifest, header)
    frame = _read_frame(manifest)

    tile_ids = _column_values(frame, "tile_id")
    _require_non_empty(manifest, "tile_id", tile_ids)
    _reject_duplicate_identifiers(manifest, tile_ids)

    tile_paths = _column_values(frame, "tile_path")
    _require_non_empty(manifest, "tile_path", tile_paths)

    # The manifest-wide columns repeat on every row so the manifest stays a
    # portable single file. Collapse each to its one value, or reject.
    wide: dict[str, str] = {}
    for column in ("slide_id",):
        values = _column_values(frame, column)
        _require_non_empty(manifest, column, values)
        wide[column] = _single_value(manifest, column, values)

    slide_id = wide["slide_id"]
    mpp_values = {
        column: _column_values(frame, column) for column in ("mpp_x", "mpp_y")
    }
    for column, values in mpp_values.items():
        _require_non_empty(manifest, column, values)
    mpp_x = _parse_manifest_mpp(manifest, "mpp_x", mpp_values["mpp_x"])
    mpp_y = _parse_manifest_mpp(manifest, "mpp_y", mpp_values["mpp_y"])

    x_px = _parse_integer_column(manifest, "x_px", _column_values(frame, "x_px"))
    y_px = _parse_integer_column(manifest, "y_px", _column_values(frame, "y_px"))
    _reject_duplicate_coordinates(manifest, x_px, y_px)

    grid_row, grid_col, row_stride, col_stride = _parse_lattice(
        manifest, frame, x_px, y_px
    )
    resolved_paths = _resolve_tile_paths(manifest, directory, tile_paths)
    if len(frame.index) < 2:
        raise DataError(
            f"Tile manifest '{manifest}' contains only one tile row. Manifest "
            f"v1 requires at least two pre-cut tiles; single-tile prediction is "
            f"outside this tool's spatial-analysis boundary."
        )

    return TileManifest(
        manifest_path=manifest,
        tile_directory=directory,
        digest=_digest_file(manifest),
        tile_ids=tile_ids,
        tile_paths=tile_paths,
        resolved_paths=resolved_paths,
        x_px=x_px,
        y_px=y_px,
        slide_id=slide_id,
        mpp_x=mpp_x,
        mpp_y=mpp_y,
        grid_row=grid_row,
        grid_col=grid_col,
        row_stride_px=row_stride,
        col_stride_px=col_stride,
    )


# =============================================================================
# Tiles
# =============================================================================


def _describe_tile(image: Any) -> str:
    """Describe a decoded tile for an error message."""
    mode = getattr(image, "mode", "unknown")
    size = getattr(image, "size", ("unknown", "unknown"))
    return f"{mode} {size[0]}x{size[1]}"


def _check_tile_image(image: Any, relative_path: str) -> None:
    """Reject a tile whose decoded geometry or color mode is not the contract.

    Resampling or color conversion here would change the model input without
    telling the user, so an out-of-contract tile is an error rather than
    something to fix silently.
    """
    if image.mode != "RGB":
        raise DataError(
            f"Tile '{relative_path}' decodes as {_describe_tile(image)}, but "
            f"DeepSpot-M consumes RGB tiles. Convert the tiles to RGB before "
            f"running this tool. ChatSpatial does not convert color modes "
            f"silently."
        )
    if image.size != (TILE_WIDTH_PX, TILE_HEIGHT_PX):
        raise DataError(
            f"Tile '{relative_path}' decodes as {_describe_tile(image)}, but "
            f"DeepSpot-M consumes {TILE_WIDTH_PX}x{TILE_HEIGHT_PX} tiles cut at "
            f"roughly 20x (about 0.5 microns per pixel). ChatSpatial neither "
            f"resamples tiles nor infers magnification, so re-cut the tiles at "
            f"the required geometry."
        )


def validate_tile_geometry(manifest: TileManifest) -> None:
    """Validate every tile's color mode and size before any model is loaded."""
    from PIL import Image

    for relative_path, path in zip(manifest.tile_paths, manifest.resolved_paths):
        try:
            with Image.open(path) as image:
                _check_tile_image(image, relative_path)
        except ChatSpatialError:
            raise
        except Exception as exc:
            raise DataError(
                f"Tile '{relative_path}' could not be decoded as an image: {exc}"
            ) from exc


def _open_tile(path: Path, relative_path: str) -> Any:
    """Decode one tile and re-check it against the geometry contract."""
    from PIL import Image

    with Image.open(path) as image:
        _check_tile_image(image, relative_path)
        return image.copy()


# =============================================================================
# Inference
# =============================================================================


def _resolve_gene_symbols(model: Any, requested: Optional[list[str]]) -> list[str]:
    """Resolve the gene panel to predict, preserving the requested order."""
    if requested:
        return list(requested)

    panel = getattr(model, "gene_names", None)
    if panel is None:
        raise ProcessingError(
            "The loaded DeepSpot-M model does not expose 'gene_names', so the "
            "default gene panel cannot be resolved. Pass an explicit list of "
            "HGNC symbols in params.genes."
        )
    symbols = [str(symbol).strip() for symbol in panel]
    if not symbols or any(not symbol for symbol in symbols):
        raise ProcessingError(
            "The loaded DeepSpot-M model exposes an empty or invalid default "
            "gene panel. Pass an explicit non-empty list of HGNC symbols in "
            "params.genes."
        )
    if len(set(symbols)) != len(symbols):
        raise ProcessingError(
            "The loaded DeepSpot-M model exposes duplicate symbols in its "
            "default gene panel. Pass an explicit unique list of HGNC symbols "
            "in params.genes."
        )
    return symbols


def _local_checkpoint_revision(directory: Path) -> str:
    """Return a content identity for the files that define a local checkpoint."""
    missing = [
        name for name in LOCAL_CHECKPOINT_FILES if not (directory / name).is_file()
    ]
    if missing:
        raise DataNotFoundError(
            f"Local DeepSpot-M checkpoint '{directory}' is missing: "
            f"{', '.join(missing)}. A checkpoint directory must contain "
            f"{', '.join(LOCAL_CHECKPOINT_FILES)}."
        )

    hasher = hashlib.sha256()
    for name in LOCAL_CHECKPOINT_FILES:
        encoded_name = name.encode("utf-8")
        hasher.update(len(encoded_name).to_bytes(4, "big"))
        hasher.update(encoded_name)
        with (directory / name).open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                hasher.update(block)
    return f"local-sha256:{hasher.hexdigest()}"


def _resolve_checkpoint_revision(params: HistologyExpressionParameters) -> str:
    """Resolve a mutable model reference to an immutable checkpoint identity."""
    require("deepspotm", feature="virtual spatial transcriptomics from H&E tiles")
    repository = Path(params.model_repository).expanduser()
    if repository.is_dir():
        if params.model_revision is not None:
            raise ParameterError(
                "params.model_revision cannot be used with a local DeepSpot-M "
                "checkpoint directory. The local checkpoint is identified by "
                "a content digest automatically."
            )
        return _local_checkpoint_revision(repository.resolve())

    huggingface_hub = require(
        "huggingface_hub", feature="resolving the DeepSpot-M checkpoint revision"
    )
    revision = params.model_revision or DEFAULT_REVISION
    try:
        info = huggingface_hub.HfApi().model_info(
            params.model_repository,
            revision=revision,
        )
    except Exception as exc:
        raise DependencyError(
            f"Could not resolve DeepSpot-M repository "
            f"'{params.model_repository}' at revision '{revision}': {exc}. "
            f"Confirm the repository and revision, accept the gated model terms "
            f"at https://huggingface.co/{params.model_repository}, and authenticate "
            f"with 'hf auth login'."
        ) from exc

    resolved = getattr(info, "sha", None)
    if not isinstance(resolved, str) or not resolved.strip():
        raise ProcessingError(
            f"Hugging Face did not return an immutable commit SHA for "
            f"'{params.model_repository}' at revision '{revision}'."
        )
    return resolved.strip()


def _load_model(
    params: HistologyExpressionParameters,
    device: str,
    checkpoint_revision: str,
) -> tuple[Any, Any]:
    """Load the published DeepSpot-M checkpoint. Blocking; run off the loop."""
    deepspotm = require(
        "deepspotm", feature="virtual spatial transcriptomics from H&E tiles"
    )
    estimator = getattr(deepspotm, "DeepSpotM", None)
    if estimator is None:
        raise DependencyError(
            "The installed deepspotm package does not expose the DeepSpotM "
            "model. Install a release that provides the from_pretrained API: "
            "pip install 'chatspatial[deepspotm]'"
        )

    torch = require("torch", feature="DeepSpot-M inference")
    return estimator.from_pretrained(
        params.model_repository,
        source=params.gene_embedding_source,
        device=torch.device(device),
        revision=(
            None
            if Path(params.model_repository).expanduser().is_dir()
            else checkpoint_revision
        ),
    )


def _predict_expression(
    model: Any,
    image_processor: Any,
    manifest: TileManifest,
    gene_symbols: list[str],
    batch_size: int,
    device: str,
    cancel_event: Optional[threading.Event] = None,
) -> np.ndarray:
    """Score every tile in manifest order. Blocking; run off the loop."""
    torch = require("torch", feature="DeepSpot-M inference")

    blocks: list[np.ndarray] = []
    for start in range(0, manifest.n_tiles, batch_size):
        if cancel_event is not None and cancel_event.is_set():
            raise ProcessingError("DeepSpot-M prediction was cancelled.")
        stop = min(start + batch_size, manifest.n_tiles)
        tiles = [
            _open_tile(
                manifest.resolved_paths[index],
                manifest.tile_paths[index],
            )
            for index in range(start, stop)
        ]
        pixel_values = torch.stack([image_processor(tile) for tile in tiles]).to(device)
        with torch.no_grad():
            predictions = model.predict_genes(pixel_values, gene_symbols)
        blocks.append(np.asarray(predictions.float().cpu().numpy(), dtype=np.float32))

    return np.concatenate(blocks, axis=0)


def _validate_predictions(
    expression: np.ndarray, n_tiles: int, n_genes: int
) -> np.ndarray:
    """Reject a prediction matrix that does not match the requested work."""
    try:
        matrix = np.asarray(expression, dtype=np.float32)
    except (TypeError, ValueError) as exc:
        raise ProcessingError(
            "DeepSpot-M returned predictions that cannot be represented as a "
            "numeric float32 matrix."
        ) from exc
    if matrix.shape != (n_tiles, n_genes):
        raise ProcessingError(
            f"DeepSpot-M returned a {matrix.shape} matrix for {n_tiles} tiles "
            f"and {n_genes} genes."
        )
    if not np.all(np.isfinite(matrix)):
        raise ProcessingError(
            "DeepSpot-M returned non-finite predicted expression values."
        )
    return matrix


# =============================================================================
# AnnData assembly
# =============================================================================


def build_virtual_anndata(
    manifest: TileManifest,
    expression: np.ndarray,
    gene_symbols: Sequence[str],
    *,
    checkpoint_revision: str,
    gene_embedding_source: str,
    model_repository: str,
    package_version: str,
) -> Any:
    """Assemble the spatial AnnData for a scored manifest.

    Predicted log1p-CPM goes into ``X`` once. The dataset carries no measured
    expression matrix that would compete for that slot, and duplicating a full
    transcriptome into a layer would double its memory footprint.
    """
    import anndata as ad

    obs = pd.DataFrame(
        {
            "slide_id": pd.Categorical([manifest.slide_id] * manifest.n_tiles),
            "tile_path": manifest.tile_paths,
            "x_px": manifest.x_px,
            "y_px": manifest.y_px,
        },
        index=pd.Index(manifest.tile_ids, name="tile_id"),
    )
    if manifest.has_lattice:
        obs["grid_row"] = manifest.grid_row
        obs["grid_col"] = manifest.grid_col

    adata = ad.AnnData(
        X=expression,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(list(gene_symbols), name="gene_symbol")),
    )
    adata.obsm["spatial"] = manifest.centers()

    set_expression_provenance(
        adata,
        provenance=PROVENANCE_PREDICTED,
        units=UNITS_LOG1P_CPM,
        producer=PRODUCER,
    )

    tile_geometry: dict[str, Any] = {
        "tile_width_px": TILE_WIDTH_PX,
        "tile_height_px": TILE_HEIGHT_PX,
        "coordinate_convention": COORDINATE_CONVENTION,
        "spatial_key": "spatial",
        "spatial_coordinates": "tile centers, (x_px + width / 2, y_px + height / 2)",
        "lattice_columns_published": bool(manifest.has_lattice),
    }
    # A single-row or single-column manifest has no observable stride on that
    # axis. Omit the key rather than writing a null that no h5ad can hold.
    if manifest.row_stride_px is not None:
        tile_geometry["grid_row_stride_px"] = manifest.row_stride_px
    if manifest.col_stride_px is not None:
        tile_geometry["grid_col_stride_px"] = manifest.col_stride_px

    adata.uns[PRODUCER] = {
        "model_repository": model_repository,
        "checkpoint_revision": checkpoint_revision,
        "gene_embedding_source": gene_embedding_source,
        "tile_geometry": tile_geometry,
        "mpp_x": manifest.mpp_x,
        "mpp_y": manifest.mpp_y,
        "slide_id": manifest.slide_id,
        "manifest": {
            "schema_version": MANIFEST_SCHEMA_VERSION,
            "format": "csv",
            "filename": manifest.manifest_path.name,
            "sha256": manifest.digest,
            "n_tiles": manifest.n_tiles,
        },
        "package_version": package_version,
        "code_license": CODE_LICENSE,
        "weights_license": WEIGHTS_LICENSE,
        "research_use_only": RESEARCH_USE_NOTICE,
        "reference": REFERENCE_DOI,
    }
    return adata


# =============================================================================
# Tool entry point
# =============================================================================


def _package_version() -> str:
    """Report the installed deepspotm version without importing metadata eagerly."""
    from importlib import metadata

    try:
        return metadata.version("deepspotm")
    except metadata.PackageNotFoundError:  # pragma: no cover - defensive
        return "unknown"


def _model_repository_identity(repository: str) -> str:
    """Return a portable repository identity without persisting host paths."""
    return "local-checkpoint" if Path(repository).expanduser().is_dir() else repository


async def predict_spatial_expression_from_histology(
    manifest_path: str,
    tile_directory: str,
    ctx: "ToolContext",
    params: Optional[HistologyExpressionParameters] = None,
) -> VirtualExpressionResult:
    """Generate virtual spatial transcriptomics from pre-cut H&E tiles.

    Args:
        manifest_path: Path to the v1 CSV coordinate manifest.
        tile_directory: Directory that every manifest ``tile_path`` resolves against.
        ctx: Tool context for data access and progress reporting.
        params: Model, gene panel, and execution parameters.

    Returns:
        A :class:`VirtualExpressionResult` naming the registered dataset.
    """
    if params is None:
        params = HistologyExpressionParameters()

    # Manifest parsing and tile validation both touch the filesystem once per
    # tile, so they stay off the event loop too. Validating everything up front
    # keeps an invalid manifest from paying for a checkpoint download.
    manifest = await asyncio.to_thread(
        read_tile_manifest, manifest_path, tile_directory
    )
    await ctx.info(
        f"Validated manifest for slide '{manifest.slide_id}': "
        f"{manifest.n_tiles} tiles at {manifest.mpp_x} x {manifest.mpp_y} "
        f"microns per pixel"
    )

    await asyncio.to_thread(validate_tile_geometry, manifest)

    device = get_device(prefer_gpu=params.use_gpu, allow_mps=False)
    if params.use_gpu and device == "cpu":
        await ctx.warning("GPU requested but not available, using CPU")

    timeout_seconds = params.timeout if params.timeout is not None else 600
    requested_revision = params.model_revision or DEFAULT_REVISION

    await ctx.info(
        f"Resolving and loading {params.model_repository} at revision "
        f"'{requested_revision}' "
        f"(gene embedding source: {params.gene_embedding_source}) on {device}"
    )

    # Checkpoint loading and inference are blocking and can take minutes. Run
    # them on a worker thread so the MCP event loop keeps serving requests and
    # can return at the configured deadline.
    loop = asyncio.get_running_loop()
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    cancel_event = threading.Event()
    inference_completed = False
    try:

        def _run() -> tuple[np.ndarray, list[str], str]:
            checkpoint_revision = _resolve_checkpoint_revision(params)
            model, image_processor = _load_model(
                params,
                device,
                checkpoint_revision,
            )
            gene_symbols = _resolve_gene_symbols(model, params.genes)
            expression = _predict_expression(
                model,
                image_processor,
                manifest,
                gene_symbols,
                params.batch_size,
                device,
                cancel_event,
            )
            return expression, gene_symbols, checkpoint_revision

        expression, gene_symbols, checkpoint_revision = await asyncio.wait_for(
            loop.run_in_executor(executor, _run),
            timeout=timeout_seconds,
        )
        inference_completed = True
    except asyncio.TimeoutError as exc:
        raise ProcessingError(
            f"DeepSpot-M timed out after {timeout_seconds} seconds while scoring "
            f"{manifest.n_tiles} tiles. Raise params.timeout, lower "
            f"params.batch_size, or split the manifest."
        ) from exc
    except KeyError as exc:
        raise ParameterError(
            f"DeepSpot-M does not know the requested gene symbol(s): {exc}. Pass "
            f"HGNC symbols from the model's panel."
        ) from exc
    except ChatSpatialError:
        raise
    except Exception as exc:
        if type(exc).__name__ in {
            "GatedRepoError",
            "RepositoryNotFoundError",
            "HfHubHTTPError",
        }:
            raise DependencyError(
                f"DeepSpot-M checkpoint access failed: {exc}. Accept the model "
                f"terms at https://huggingface.co/{params.model_repository} and "
                f"authenticate with 'hf auth login'."
            ) from exc
        raise ProcessingError(f"DeepSpot-M prediction failed: {exc}") from exc
    finally:
        if not inference_completed:
            cancel_event.set()
        # A context manager waits for the worker on exit and therefore defeats
        # wait_for(). Do not block the MCP response after a timeout; completed
        # work receives an orderly shutdown, while active inference observes
        # the cancellation flag before starting another batch.
        executor.shutdown(wait=inference_completed, cancel_futures=True)

    matrix = _validate_predictions(expression, manifest.n_tiles, len(gene_symbols))

    model_repository = _model_repository_identity(params.model_repository)
    adata = build_virtual_anndata(
        manifest,
        matrix,
        gene_symbols,
        checkpoint_revision=checkpoint_revision,
        gene_embedding_source=params.gene_embedding_source,
        model_repository=model_repository,
        package_version=_package_version(),
    )

    # Reserve first, validate the public result, then publish once. The dataset
    # becomes reachable only after a complete and valid AnnData exists.
    data_id = ctx.reserve_dataset_id("histology")
    result = VirtualExpressionResult(
        data_id=data_id,
        slide_id=manifest.slide_id,
        n_tiles=manifest.n_tiles,
        n_genes=len(gene_symbols),
        model_repository=model_repository,
        checkpoint_revision=checkpoint_revision,
        gene_embedding_source=params.gene_embedding_source,
        expression_provenance=PROVENANCE_PREDICTED,
        expression_units=UNITS_LOG1P_CPM,
        tile_width_px=TILE_WIDTH_PX,
        tile_height_px=TILE_HEIGHT_PX,
        mpp_x=manifest.mpp_x,
        mpp_y=manifest.mpp_y,
        manifest_sha256=manifest.digest,
        lattice_columns_published=manifest.has_lattice,
    )
    await ctx.add_dataset(
        adata,
        prefix="histology",
        name=params.name or f"Virtual spatial transcriptomics: {manifest.slide_id}",
        metadata={"type": "virtual_histology"},
        data_id=data_id,
    )
    await ctx.info(
        f"Registered virtual spatial transcriptomics dataset '{data_id}' with "
        f"{manifest.n_tiles} tiles and {len(gene_symbols)} genes. Use "
        f"export_data to write it to disk."
    )
    return result
