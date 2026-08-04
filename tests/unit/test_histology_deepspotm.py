"""Unit contracts for DeepSpot-M virtual spatial transcriptomics.

Everything here is mocked. No test installs deepspotm, imports torch, reaches
Hugging Face, or touches the network in any other way.
"""

from __future__ import annotations

import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from PIL import Image

from chatspatial.tools import histology
from chatspatial.utils import dependency_manager as dm
from chatspatial.utils.exceptions import (
    DataCompatibilityError,
    DataError,
    DataNotFoundError,
    DependencyError,
    ParameterError,
    ProcessingError,
)
from chatspatial.utils.provenance import (
    is_predicted_expression,
    set_expression_provenance,
)

MANIFEST_HEADER = "tile_id,tile_path,slide_id,x_px,y_px,mpp_x,mpp_y"
LATTICE_HEADER = MANIFEST_HEADER + ",grid_row,grid_col"


# =============================================================================
# Fixtures
# =============================================================================


def _write_tile(
    path: Path,
    *,
    size: tuple[int, int] = (224, 224),
    mode: str = "RGB",
    level: int = 32,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    color = level if mode == "L" else (level, level, level)
    Image.new(mode, size, color=color).save(path)


@pytest.fixture
def tile_dir(tmp_path: Path) -> Path:
    """Four distinguishable tiles, so tile order is observable in the output."""
    directory = tmp_path / "tiles"
    directory.mkdir()
    for level, name in enumerate(("a.png", "b.png", "c.png", "d.png"), start=1):
        _write_tile(directory / name, level=level * 10)
    return directory


def _manifest(tmp_path: Path, rows: list[str], header: str = MANIFEST_HEADER) -> str:
    path = tmp_path / "manifest.csv"
    path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
    return str(path)


@pytest.fixture
def simple_manifest(tmp_path: Path, tile_dir: Path) -> str:
    return _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,0,0,0.5,0.5",
            "t1,b.png,SLIDE-1,224,0,0.5,0.5",
            "t2,c.png,SLIDE-1,0,224,0.5,0.5",
            "t3,d.png,SLIDE-1,224,224,0.5,0.5",
        ],
    )


class _RecordingContext:
    """Minimal ToolContext stand-in that records dataset publication."""

    def __init__(self) -> None:
        self.messages: list[str] = []
        self.published: dict[str, object] = {}
        self.reserved: list[str] = []

    async def info(self, msg: str) -> None:
        self.messages.append(msg)

    async def warning(self, msg: str) -> None:
        self.messages.append(f"warning: {msg}")

    def reserve_dataset_id(self, prefix: str = "data") -> str:
        data_id = f"{prefix}_{len(self.reserved) + 1}"
        self.reserved.append(data_id)
        return data_id

    async def add_dataset(
        self, adata, prefix="data", name=None, metadata=None, data_id=None
    ):
        assert data_id in self.reserved, "dataset must publish under a reserved id"
        self.published[data_id] = adata
        return data_id


# =============================================================================
# Manifest schema
# =============================================================================


@pytest.mark.unit
def test_valid_manifest_is_parsed_in_row_order(simple_manifest: str, tile_dir: Path):
    manifest = histology.read_tile_manifest(simple_manifest, str(tile_dir))

    assert manifest.tile_ids == ["t0", "t1", "t2", "t3"]
    assert manifest.slide_id == "SLIDE-1"
    assert manifest.mpp_x == 0.5
    assert manifest.mpp_y == 0.5
    assert manifest.n_tiles == 4
    assert manifest.has_lattice is False
    assert np.array_equal(manifest.x_px, np.array([0, 224, 0, 224]))
    assert np.array_equal(manifest.y_px, np.array([0, 0, 224, 224]))


@pytest.mark.unit
def test_tile_centers_use_the_upper_left_corner_plus_half_the_tile(
    simple_manifest: str, tile_dir: Path
):
    manifest = histology.read_tile_manifest(simple_manifest, str(tile_dir))

    expected = np.array(
        [[112.0, 112.0], [336.0, 112.0], [112.0, 336.0], [336.0, 336.0]]
    )
    assert np.allclose(manifest.centers(), expected)


@pytest.mark.unit
def test_missing_manifest_file_is_rejected(tile_dir: Path, tmp_path: Path):
    with pytest.raises(ParameterError, match="does not exist"):
        histology.read_tile_manifest(str(tmp_path / "absent.csv"), str(tile_dir))


@pytest.mark.unit
def test_missing_tile_directory_is_rejected(simple_manifest: str, tmp_path: Path):
    with pytest.raises(ParameterError, match="Tile directory"):
        histology.read_tile_manifest(simple_manifest, str(tmp_path / "absent"))


@pytest.mark.unit
def test_empty_manifest_is_rejected(tmp_path: Path, tile_dir: Path):
    path = tmp_path / "manifest.csv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(DataError, match="empty"):
        histology.read_tile_manifest(str(path), str(tile_dir))


@pytest.mark.unit
def test_header_only_manifest_is_rejected(tmp_path: Path, tile_dir: Path):
    with pytest.raises(DataError, match="no tile rows"):
        histology.read_tile_manifest(_manifest(tmp_path, []), str(tile_dir))


@pytest.mark.unit
def test_missing_required_columns_are_named(tmp_path: Path, tile_dir: Path):
    path = tmp_path / "manifest.csv"
    path.write_text("tile_id,tile_path\nt0,a.png\n", encoding="utf-8")
    with pytest.raises(DataError, match="missing required v1 columns"):
        histology.read_tile_manifest(str(path), str(tile_dir))


@pytest.mark.unit
def test_duplicate_columns_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5,1"],
        header=MANIFEST_HEADER + ",x_px",
    )
    with pytest.raises(DataError, match="duplicate columns"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_duplicate_tile_identifiers_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5", "t0,b.png,SLIDE-1,224,0,0.5,0.5"],
    )
    with pytest.raises(DataError, match="repeats tile identifiers"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_duplicate_coordinates_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5", "t1,b.png,SLIDE-1,0,0,0.5,0.5"],
    )
    with pytest.raises(DataError, match="same upper-left corner"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_mixed_slide_identifiers_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5", "t1,b.png,SLIDE-2,224,0,0.5,0.5"],
    )
    with pytest.raises(DataError, match="mixes 2 'slide_id' values"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_inconsistent_mpp_is_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5", "t1,b.png,SLIDE-1,224,0,0.25,0.5"],
    )
    with pytest.raises(DataError, match="mixes 2 numeric 'mpp_x' values"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_equivalent_mpp_spellings_are_consistent(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.50", "t1,b.png,SLIDE-1,224,0,0.50,0.5"],
    )
    manifest = histology.read_tile_manifest(path, str(tile_dir))
    assert manifest.mpp_x == 0.5
    assert manifest.mpp_y == 0.5


@pytest.mark.unit
def test_single_tile_manifest_is_outside_the_spatial_boundary(
    tmp_path: Path, tile_dir: Path
):
    path = _manifest(tmp_path, ["t0,a.png,SLIDE-1,0,0,0.5,0.5"])
    with pytest.raises(DataError, match="at least two pre-cut tiles"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_duplicate_resolved_tile_paths_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5", "t1,./a.png,SLIDE-1,224,0,0.5,0.5"],
    )
    with pytest.raises(DataError, match="same decoded tile more than once"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
@pytest.mark.parametrize("bad", ["nan", "inf", "-inf"])
def test_non_finite_coordinates_are_rejected(tmp_path: Path, tile_dir: Path, bad: str):
    path = _manifest(tmp_path, [f"t0,a.png,SLIDE-1,{bad},0,0.5,0.5"])
    with pytest.raises(DataError, match="non-finite|negative"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_negative_coordinates_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(tmp_path, ["t0,a.png,SLIDE-1,-1,0,0.5,0.5"])
    with pytest.raises(DataError, match="negative 'x_px'"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_fractional_coordinates_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(tmp_path, ["t0,a.png,SLIDE-1,0.5,0,0.5,0.5"])
    with pytest.raises(DataError, match="fractional 'x_px'"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_non_numeric_coordinates_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(tmp_path, ["t0,a.png,SLIDE-1,left,0,0.5,0.5"])
    with pytest.raises(DataError, match="non-numeric 'x_px'"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
@pytest.mark.parametrize("value", ["0", "-0.5", "not-a-number"])
def test_invalid_mpp_is_rejected(tmp_path: Path, tile_dir: Path, value: str):
    path = _manifest(tmp_path, [f"t0,a.png,SLIDE-1,0,0,{value},0.5"])
    with pytest.raises(DataError, match="'mpp_x'"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_blank_tile_identifier_is_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(tmp_path, [",a.png,SLIDE-1,0,0,0.5,0.5"])
    with pytest.raises(DataError, match="blank 'tile_id'"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_missing_tile_file_is_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(tmp_path, ["t0,absent.png,SLIDE-1,0,0,0.5,0.5"])
    with pytest.raises(DataNotFoundError, match="do not exist"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_tile_path_escaping_the_tile_directory_is_rejected(
    tmp_path: Path, tile_dir: Path
):
    outside = tmp_path / "outside.png"
    _write_tile(outside)
    path = _manifest(tmp_path, ["t0,../outside.png,SLIDE-1,0,0,0.5,0.5"])
    with pytest.raises(DataError, match="resolve outside the declared tile directory"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_absolute_tile_path_is_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(tmp_path, [f"t0,{tile_dir / 'a.png'},SLIDE-1,0,0,0.5,0.5"])
    with pytest.raises(DataError, match="absolute tile_path"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_manifest_digest_identifies_the_file(simple_manifest: str, tile_dir: Path):
    first = histology.read_tile_manifest(simple_manifest, str(tile_dir))
    Path(simple_manifest).write_text(
        Path(simple_manifest).read_text(encoding="utf-8").replace("t3", "t9"),
        encoding="utf-8",
    )
    second = histology.read_tile_manifest(simple_manifest, str(tile_dir))

    assert len(first.digest) == 64
    assert first.digest != second.digest


# =============================================================================
# Optional lattice columns
# =============================================================================


@pytest.mark.unit
def test_consistent_lattice_is_published_with_its_stride(
    tmp_path: Path, tile_dir: Path
):
    path = _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,0,0,0.5,0.5,0,0",
            "t1,b.png,SLIDE-1,224,0,0.5,0.5,0,1",
            "t2,c.png,SLIDE-1,0,224,0.5,0.5,1,0",
            "t3,d.png,SLIDE-1,224,224,0.5,0.5,1,1",
        ],
        header=LATTICE_HEADER,
    )
    manifest = histology.read_tile_manifest(path, str(tile_dir))

    assert manifest.has_lattice is True
    assert manifest.row_stride_px == 224
    assert manifest.col_stride_px == 224


@pytest.mark.unit
def test_lattice_offset_does_not_have_to_be_zero(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,1000,500,0.5,0.5,2,5",
            "t1,b.png,SLIDE-1,1224,500,0.5,0.5,2,6",
            "t2,c.png,SLIDE-1,1000,724,0.5,0.5,3,5",
        ],
        header=LATTICE_HEADER,
    )
    manifest = histology.read_tile_manifest(path, str(tile_dir))

    assert manifest.col_stride_px == 224
    assert manifest.row_stride_px == 224


@pytest.mark.unit
def test_lattice_without_its_partner_column_is_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5,0"],
        header=MANIFEST_HEADER + ",grid_row",
    )
    with pytest.raises(DataError, match="must appear together"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_irregular_lattice_is_rejected_rather_than_published(
    tmp_path: Path, tile_dir: Path
):
    path = _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,0,0,0.5,0.5,0,0",
            "t1,b.png,SLIDE-1,224,0,0.5,0.5,0,1",
            "t2,c.png,SLIDE-1,500,0,0.5,0.5,0,2",
        ],
        header=LATTICE_HEADER,
    )
    with pytest.raises(DataError, match="one affine mapping"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_fractional_lattice_stride_is_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,0,0,0.5,0.5,0,0",
            "t1,b.png,SLIDE-1,225,0,0.5,0.5,0,2",
        ],
        header=LATTICE_HEADER,
    )
    with pytest.raises(DataError, match="whole-pixel stride"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_duplicate_lattice_positions_are_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,0,0,0.5,0.5,0,0",
            "t1,b.png,SLIDE-1,224,0,0.5,0.5,0,0",
        ],
        header=LATTICE_HEADER,
    )
    with pytest.raises(DataError, match="repeats \\(grid_row, grid_col\\) pairs"):
        histology.read_tile_manifest(path, str(tile_dir))


@pytest.mark.unit
def test_single_grid_line_stays_degenerate_without_failing(
    tmp_path: Path, tile_dir: Path
):
    path = _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,0,0,0.5,0.5,0,0",
            "t1,b.png,SLIDE-1,224,0,0.5,0.5,0,1",
        ],
        header=LATTICE_HEADER,
    )
    manifest = histology.read_tile_manifest(path, str(tile_dir))

    assert manifest.has_lattice is True
    assert manifest.row_stride_px is None
    assert manifest.col_stride_px == 224


@pytest.mark.unit
def test_lattice_index_contradicting_pixels_is_rejected(tmp_path: Path, tile_dir: Path):
    path = _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,0,0,0.5,0.5,0,0",
            "t1,b.png,SLIDE-1,224,0,0.5,0.5,0,1",
            "t2,c.png,SLIDE-1,0,224,0.5,0.5,0,2",
        ],
        header=LATTICE_HEADER,
    )
    with pytest.raises(DataError, match="several different 'y_px' values"):
        histology.read_tile_manifest(path, str(tile_dir))


# =============================================================================
# Tile geometry
# =============================================================================


@pytest.mark.unit
def test_tiles_matching_the_contract_validate(simple_manifest: str, tile_dir: Path):
    manifest = histology.read_tile_manifest(simple_manifest, str(tile_dir))
    histology.validate_tile_geometry(manifest)


@pytest.mark.unit
def test_wrong_tile_size_is_rejected_without_resampling(tmp_path: Path, tile_dir: Path):
    _write_tile(tile_dir / "a.png", size=(256, 256))
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5", "t1,b.png,SLIDE-1,224,0,0.5,0.5"],
    )
    manifest = histology.read_tile_manifest(path, str(tile_dir))

    with pytest.raises(DataError, match="neither resamples tiles nor infers"):
        histology.validate_tile_geometry(manifest)


@pytest.mark.unit
def test_non_rgb_tile_is_rejected_without_conversion(tmp_path: Path, tile_dir: Path):
    _write_tile(tile_dir / "a.png", mode="L")
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5", "t1,b.png,SLIDE-1,224,0,0.5,0.5"],
    )
    manifest = histology.read_tile_manifest(path, str(tile_dir))

    with pytest.raises(DataError, match="does not convert color modes"):
        histology.validate_tile_geometry(manifest)


@pytest.mark.unit
def test_undecodable_tile_is_rejected(tmp_path: Path, tile_dir: Path):
    (tile_dir / "a.png").write_text("not an image", encoding="utf-8")
    path = _manifest(
        tmp_path,
        ["t0,a.png,SLIDE-1,0,0,0.5,0.5", "t1,b.png,SLIDE-1,224,0,0.5,0.5"],
    )
    manifest = histology.read_tile_manifest(path, str(tile_dir))

    with pytest.raises(DataError, match="could not be decoded"):
        histology.validate_tile_geometry(manifest)


# =============================================================================
# Prediction and AnnData assembly
# =============================================================================


@pytest.mark.unit
@pytest.mark.asyncio
async def test_prediction_registers_a_dataset_with_the_agreed_layout(
    simple_manifest: str, tile_dir: Path, deepspotm_stack: dict[str, object]
):
    from chatspatial.models.data import HistologyExpressionParameters

    ctx = _RecordingContext()
    result = await histology.predict_spatial_expression_from_histology(
        simple_manifest,
        str(tile_dir),
        ctx,  # type: ignore[arg-type]
        HistologyExpressionParameters(
            genes=["EPCAM", "CD3D"], batch_size=3, use_gpu=False
        ),
    )

    assert result.data_id in ctx.published
    adata = ctx.published[result.data_id]

    # Expression lives in X once, never duplicated into a layer.
    assert adata.X.shape == (4, 2)
    assert adata.X.dtype == np.float32
    # AnnData exposes X itself under the None key, so only named layers count.
    assert [key for key in adata.layers if key is not None] == []

    # Manifest order is preserved and centers follow the agreed formula.
    assert list(adata.obs_names) == ["t0", "t1", "t2", "t3"]
    assert list(adata.var_names) == ["EPCAM", "CD3D"]
    assert np.allclose(
        adata.obsm["spatial"],
        np.array([[112.0, 112.0], [336.0, 112.0], [112.0, 336.0], [336.0, 336.0]]),
    )
    assert list(adata.obs["x_px"]) == [0, 224, 0, 224]
    assert "grid_row" not in adata.obs

    # Row i of X comes from tile i of the manifest, in manifest order.
    assert np.allclose(adata.X[:, 0], np.array([10.0, 20.0, 30.0, 40.0]))
    assert np.allclose(adata.X[:, 1], 2 * adata.X[:, 0])

    # Batching respects params.batch_size without dropping or reordering tiles.
    assert deepspotm_stack["model"].batches == [3, 1]  # type: ignore[union-attr]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_generic_and_backend_provenance_blocks_are_published(
    simple_manifest: str, tile_dir: Path, deepspotm_stack: dict[str, object]
):
    from chatspatial.models.data import HistologyExpressionParameters

    ctx = _RecordingContext()
    result = await histology.predict_spatial_expression_from_histology(
        simple_manifest,
        str(tile_dir),
        ctx,  # type: ignore[arg-type]
        HistologyExpressionParameters(model_revision="v1.0.0", use_gpu=False),
    )
    adata = ctx.published[result.data_id]

    generic = adata.uns["chatspatial"]["expression"]
    assert generic == {
        "schema_version": 1,
        "provenance": "predicted",
        "units": "log1p_cpm",
        "producer": "deepspotm",
    }
    assert is_predicted_expression(adata)

    backend = adata.uns["deepspotm"]
    assert backend["model_repository"] == "ratschlab/DeepSpotM"
    assert backend["checkpoint_revision"] == "v1.0.0"
    assert backend["gene_embedding_source"] == "scgpt"
    assert backend["mpp_x"] == 0.5
    assert backend["mpp_y"] == 0.5
    assert backend["slide_id"] == "SLIDE-1"
    assert backend["tile_geometry"]["tile_width_px"] == 224
    assert backend["tile_geometry"]["tile_height_px"] == 224
    assert "level-0" in backend["tile_geometry"]["coordinate_convention"]
    assert backend["manifest"]["schema_version"] == 1
    assert backend["manifest"]["filename"] == "manifest.csv"
    assert len(backend["manifest"]["sha256"]) == 64
    assert "path" not in backend["manifest"]
    assert "tile_directory" not in backend["manifest"]
    assert backend["code_license"] == "PolyForm-Noncommercial-1.0.0"
    assert backend["weights_license"] == "CC-BY-NC-SA-4.0"
    assert "Research use only" in backend["research_use_only"]
    assert "10.64898" in backend["reference"]

    # The full model panel is used when no gene list is requested.
    assert result.n_genes == 3
    assert result.checkpoint_revision == "v1.0.0"
    assert result.lattice_columns_published is False


@pytest.mark.unit
@pytest.mark.asyncio
async def test_validated_lattice_columns_reach_obs(
    tmp_path: Path, tile_dir: Path, deepspotm_stack: dict[str, object]
):
    from chatspatial.models.data import HistologyExpressionParameters

    path = _manifest(
        tmp_path,
        [
            "t0,a.png,SLIDE-1,0,0,0.5,0.5,0,0",
            "t1,b.png,SLIDE-1,224,0,0.5,0.5,0,1",
            "t2,c.png,SLIDE-1,0,224,0.5,0.5,1,0",
        ],
        header=LATTICE_HEADER,
    )
    ctx = _RecordingContext()
    result = await histology.predict_spatial_expression_from_histology(
        path,
        str(tile_dir),
        ctx,  # type: ignore[arg-type]
        HistologyExpressionParameters(use_gpu=False),
    )
    adata = ctx.published[result.data_id]

    assert result.lattice_columns_published is True
    assert list(adata.obs["grid_row"]) == [0, 0, 1]
    assert list(adata.obs["grid_col"]) == [0, 1, 0]
    # Generic lattice indices only. Visium keys are never invented.
    assert "array_row" not in adata.obs
    assert "array_col" not in adata.obs


@pytest.mark.unit
@pytest.mark.asyncio
async def test_model_receives_the_resolved_repository_and_revision(
    simple_manifest: str, tile_dir: Path, deepspotm_stack: dict[str, object]
):
    from chatspatial.models.data import HistologyExpressionParameters

    ctx = _RecordingContext()
    await histology.predict_spatial_expression_from_histology(
        simple_manifest,
        str(tile_dir),
        ctx,  # type: ignore[arg-type]
        HistologyExpressionParameters(
            gene_embedding_source="evo2", model_revision="abc123", use_gpu=False
        ),
    )

    assert deepspotm_stack["repo"] == "ratschlab/DeepSpotM"
    assert deepspotm_stack["source"] == "evo2"
    assert deepspotm_stack["revision"] == "abc123"


@pytest.mark.unit
@pytest.mark.asyncio
async def test_default_revision_is_resolved_and_pinned(
    simple_manifest: str, tile_dir: Path, deepspotm_stack: dict[str, object]
):
    from chatspatial.models.data import HistologyExpressionParameters

    ctx = _RecordingContext()
    result = await histology.predict_spatial_expression_from_histology(
        simple_manifest,
        str(tile_dir),
        ctx,  # type: ignore[arg-type]
        HistologyExpressionParameters(use_gpu=False),
    )

    assert deepspotm_stack["revision"] == "resolved-main-sha"
    assert result.checkpoint_revision == "resolved-main-sha"


@pytest.mark.unit
def test_remote_revision_resolves_to_hugging_face_commit(
    monkeypatch: pytest.MonkeyPatch,
):
    from chatspatial.models.data import HistologyExpressionParameters

    calls: list[tuple[str, str]] = []

    class FakeApi:
        def model_info(self, repository: str, *, revision: str):
            calls.append((repository, revision))
            return SimpleNamespace(sha="0123456789abcdef")

    def fake_require(name: str, **_kwargs):
        if name == "huggingface_hub":
            return SimpleNamespace(HfApi=FakeApi)
        return object()

    monkeypatch.setattr(histology, "require", fake_require)
    params = HistologyExpressionParameters(model_revision="release-v1")

    assert histology._resolve_checkpoint_revision(params) == "0123456789abcdef"
    assert calls == [("ratschlab/DeepSpotM", "release-v1")]


@pytest.mark.unit
def test_local_checkpoint_revision_is_content_addressed(tmp_path: Path):
    checkpoint = tmp_path / "checkpoint"
    checkpoint.mkdir()
    for name in histology.LOCAL_CHECKPOINT_FILES:
        (checkpoint / name).write_bytes(f"first:{name}".encode())

    first = histology._local_checkpoint_revision(checkpoint)
    (checkpoint / "tokens.csv").write_text("second", encoding="utf-8")
    second = histology._local_checkpoint_revision(checkpoint)

    assert first.startswith("local-sha256:")
    assert len(first.removeprefix("local-sha256:")) == 64
    assert second != first


@pytest.mark.unit
@pytest.mark.asyncio
async def test_timeout_returns_without_waiting_for_the_worker(
    simple_manifest: str,
    tile_dir: Path,
    deepspotm_stack: dict[str, object],
    monkeypatch: pytest.MonkeyPatch,
):
    from chatspatial.models.data import HistologyExpressionParameters

    del deepspotm_stack

    def slow_resolution(_params):
        time.sleep(1.5)
        raise ProcessingError("worker stopped after the public timeout")

    monkeypatch.setattr(histology, "_resolve_checkpoint_revision", slow_resolution)
    started = time.monotonic()
    with pytest.raises(ProcessingError, match="timed out after 1 seconds"):
        await histology.predict_spatial_expression_from_histology(
            simple_manifest,
            str(tile_dir),
            _RecordingContext(),  # type: ignore[arg-type]
            HistologyExpressionParameters(use_gpu=False, timeout=1),
        )
    assert time.monotonic() - started < 1.4


@pytest.mark.unit
@pytest.mark.asyncio
async def test_unknown_gene_symbols_produce_an_actionable_error(
    simple_manifest: str, tile_dir: Path, deepspotm_stack: dict[str, object]
):
    from chatspatial.models.data import HistologyExpressionParameters

    ctx = _RecordingContext()
    with pytest.raises(ParameterError, match="does not know the requested gene"):
        await histology.predict_spatial_expression_from_histology(
            simple_manifest,
            str(tile_dir),
            ctx,  # type: ignore[arg-type]
            HistologyExpressionParameters(genes=["NOT_A_GENE"], use_gpu=False),
        )
    assert ctx.published == {}


@pytest.mark.unit
@pytest.mark.asyncio
async def test_missing_package_reports_the_dedicated_extra(
    simple_manifest: str, tile_dir: Path, monkeypatch: pytest.MonkeyPatch
):
    from chatspatial.models.data import HistologyExpressionParameters

    monkeypatch.setattr(dm, "_try_import", lambda module_name: None)
    monkeypatch.setattr(dm, "_check_spec", lambda module_name: False)

    ctx = _RecordingContext()
    with pytest.raises(DependencyError) as excinfo:
        await histology.predict_spatial_expression_from_histology(
            simple_manifest,
            str(tile_dir),
            ctx,  # type: ignore[arg-type]
            HistologyExpressionParameters(use_gpu=False),
        )

    message = str(excinfo.value)
    assert "chatspatial[deepspotm]" in message
    assert "Non-commercial" in message
    assert "huggingface-cli login" in message


@pytest.mark.unit
def test_prediction_shape_mismatch_is_rejected():
    with pytest.raises(ProcessingError, match="matrix for 4 tiles"):
        histology._validate_predictions(np.zeros((3, 2), dtype=np.float32), 4, 2)


@pytest.mark.unit
def test_non_finite_predictions_are_rejected():
    matrix = np.zeros((2, 2), dtype=np.float32)
    matrix[0, 0] = np.nan
    with pytest.raises(ProcessingError, match="non-finite"):
        histology._validate_predictions(matrix, 2, 2)


# =============================================================================
# Shared provenance guards
# =============================================================================


def _predicted_adata():
    import anndata as ad

    adata = ad.AnnData(X=np.ones((4, 3), dtype=np.float32))
    set_expression_provenance(
        adata, provenance="predicted", units="log1p_cpm", producer="deepspotm"
    )
    return adata


@pytest.mark.unit
def test_integer_count_access_rejects_predicted_expression_before_heuristics():
    from chatspatial.utils.adata_utils import get_raw_data_source

    with pytest.raises(DataCompatibilityError, match="Raw count access"):
        get_raw_data_source(_predicted_adata(), require_integer_counts=True)


@pytest.mark.unit
def test_normalized_expression_access_accepts_predicted_expression():
    from chatspatial.utils.adata_utils import get_raw_data_source

    result = get_raw_data_source(_predicted_adata(), require_integer_counts=False)
    assert result.source == "current"
    assert result.is_integer_counts is False


@pytest.mark.unit
def test_counts_layer_creation_rejects_predicted_expression():
    from chatspatial.utils.adata_utils import ensure_counts_layer

    with pytest.raises(DataCompatibilityError, match="counts layer"):
        ensure_counts_layer(_predicted_adata())


@pytest.mark.unit
def test_measured_datasets_are_untouched_by_the_guard():
    import anndata as ad

    from chatspatial.utils.adata_utils import get_raw_data_source

    adata = ad.AnnData(X=np.ones((4, 3), dtype=np.float32))
    set_expression_provenance(
        adata, provenance="measured", units="counts", producer="load_data"
    )
    assert get_raw_data_source(adata).source == "current"


@pytest.mark.unit
def test_datasets_without_metadata_are_untouched_by_the_guard():
    import anndata as ad

    from chatspatial.utils.adata_utils import get_raw_data_source

    adata = ad.AnnData(X=np.ones((4, 3), dtype=np.float32))
    assert get_raw_data_source(adata).source == "current"


@pytest.mark.unit
def test_unknown_provenance_values_are_rejected_at_write_time():
    import anndata as ad

    adata = ad.AnnData(X=np.ones((2, 2), dtype=np.float32))
    with pytest.raises(ParameterError, match="Unknown expression provenance"):
        set_expression_provenance(
            adata, provenance="guessed", units="log1p_cpm", producer="x"
        )
