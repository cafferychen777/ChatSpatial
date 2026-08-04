"""MCP contract tests for the virtual spatial transcriptomics tool.

These exercise the registered tool through the MCP surface and the shared
provenance guards through the public server entry points. Nothing here installs
deepspotm or reaches the network.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from mcp.client import Client
from PIL import Image

from chatspatial.models.analysis import VirtualExpressionResult
from chatspatial.models.data import HistologyExpressionParameters
from chatspatial.server import (
    compute_embeddings,
    data_manager,
    export_data,
    mcp,
    predict_spatial_expression_from_histology,
    preprocess_data,
    reload_data,
)
from chatspatial.tools.embeddings import EmbeddingParameters
from chatspatial.utils.exceptions import DataCompatibilityError, ParameterError
from chatspatial.utils.provenance import set_expression_provenance

TOOL_NAME = "predict_spatial_expression_from_histology"


@pytest.fixture
def tile_workspace(tmp_path: Path) -> tuple[str, str]:
    """A two-tile slide plus its v1 manifest."""
    tiles = tmp_path / "tiles"
    tiles.mkdir()
    for level, name in ((10, "a.png"), (20, "b.png")):
        Image.new("RGB", (224, 224), color=(level, level, level)).save(tiles / name)

    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "tile_id,tile_path,slide_id,x_px,y_px,mpp_x,mpp_y\n"
        "t0,a.png,SLIDE-1,0,0,0.5,0.5\n"
        "t1,b.png,SLIDE-1,224,0,0.5,0.5\n",
        encoding="utf-8",
    )
    return str(manifest), str(tiles)


# =============================================================================
# Tool registration
# =============================================================================


@pytest.mark.integration
@pytest.mark.asyncio
async def test_tool_is_registered_with_closed_schemas_and_safety_hints():
    tools = {tool.name: tool for tool in await mcp.list_tools()}

    assert TOOL_NAME in tools
    tool = tools[TOOL_NAME]

    assert tool.input_schema["additionalProperties"] is False
    assert tool.output_schema is not None
    assert tool.output_schema["additionalProperties"] is False

    assert tool.annotations is not None
    assert tool.annotations.read_only_hint is False
    assert tool.annotations.destructive_hint is False
    assert tool.annotations.idempotent_hint is False
    assert tool.annotations.open_world_hint is True


@pytest.mark.integration
@pytest.mark.asyncio
async def test_input_schema_requires_the_manifest_and_tile_directory():
    tools = {tool.name: tool for tool in await mcp.list_tools()}
    schema = tools[TOOL_NAME].input_schema

    assert set(schema["required"]) == {"manifest_path", "tile_directory"}
    assert "params" in schema["properties"]


@pytest.mark.integration
@pytest.mark.asyncio
async def test_output_schema_publishes_provenance_and_model_identity():
    tools = {tool.name: tool for tool in await mcp.list_tools()}
    schema = tools[TOOL_NAME].output_schema

    assert {
        "data_id",
        "slide_id",
        "n_tiles",
        "n_genes",
        "model_repository",
        "checkpoint_revision",
        "gene_embedding_source",
        "expression_provenance",
        "expression_units",
        "tile_width_px",
        "tile_height_px",
        "mpp_x",
        "mpp_y",
        "manifest_sha256",
        "lattice_columns_published",
    } <= set(schema["properties"])

    assert schema["properties"]["expression_provenance"]["const"] == "predicted"
    assert schema["properties"]["expression_units"]["const"] == "log1p_cpm"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_parameter_schema_exposes_the_published_gene_embedding_sources():
    tools = {tool.name: tool for tool in await mcp.list_tools()}
    schema = tools[TOOL_NAME].input_schema

    definitions = schema.get("$defs", {})
    params_schema = definitions["HistologyExpressionParameters"]
    source = params_schema["properties"]["gene_embedding_source"]

    assert source["default"] == "scgpt"
    assert set(source["enum"]) == {"evo2", "orthrus", "prott5", "scgpt", "apertus"}
    assert params_schema["additionalProperties"] is False


@pytest.mark.integration
def test_parameters_reject_unknown_and_invalid_values():
    with pytest.raises(Exception):
        HistologyExpressionParameters(unknown_option=True)
    with pytest.raises(Exception):
        HistologyExpressionParameters(gene_embedding_source="not-a-source")
    with pytest.raises(Exception):
        HistologyExpressionParameters(batch_size=0)
    with pytest.raises(Exception):
        HistologyExpressionParameters(genes=["EPCAM", "EPCAM"])


# =============================================================================
# Tool behavior through the server boundary
# =============================================================================


@pytest.mark.integration
@pytest.mark.asyncio
async def test_invalid_manifest_fails_before_any_model_is_loaded(
    tmp_path: Path, reset_data_manager
):
    tiles = tmp_path / "tiles"
    tiles.mkdir()

    with pytest.raises(ParameterError, match="does not exist"):
        await predict_spatial_expression_from_histology(
            str(tmp_path / "absent.csv"), str(tiles)
        )

    assert await data_manager.list_datasets() == []


@pytest.mark.integration
@pytest.mark.asyncio
async def test_successful_prediction_registers_the_dataset_and_returns_a_data_id(
    tile_workspace: tuple[str, str],
    deepspotm_stack: dict[str, object],
    reset_data_manager,
):
    manifest_path, tile_directory = tile_workspace

    result = await predict_spatial_expression_from_histology(
        manifest_path,
        tile_directory,
        params=HistologyExpressionParameters(genes=["EPCAM"], use_gpu=False),
    )

    assert isinstance(result, VirtualExpressionResult)
    assert result.n_tiles == 2
    assert result.n_genes == 1
    assert result.slide_id == "SLIDE-1"
    assert result.expression_provenance == "predicted"
    assert result.expression_units == "log1p_cpm"

    listed = await data_manager.list_datasets()
    assert [entry["id"] for entry in listed] == [result.data_id]

    dataset = await data_manager.get_dataset(result.data_id)
    adata = dataset["adata"]
    assert adata.shape == (2, 1)
    assert dataset["spatial_coordinates_available"] is True
    assert adata.uns["chatspatial"]["expression"]["provenance"] == "predicted"


@pytest.mark.integration
@pytest.mark.asyncio
async def test_prediction_runs_through_the_mcp_transport(
    tile_workspace: tuple[str, str],
    deepspotm_stack: dict[str, object],
    reset_data_manager,
):
    del deepspotm_stack
    manifest_path, tile_directory = tile_workspace

    async with Client(mcp) as client:
        response = await client.call_tool(
            TOOL_NAME,
            {
                "manifest_path": manifest_path,
                "tile_directory": tile_directory,
                "params": {
                    "genes": ["EPCAM", "CD3D"],
                    "use_gpu": False,
                },
            },
        )

    assert response.is_error is False
    assert response.structured_content is not None
    assert response.structured_content["n_tiles"] == 2
    assert response.structured_content["n_genes"] == 2
    data_id = response.structured_content["data_id"]
    dataset = await data_manager.get_dataset(data_id)
    assert dataset["type"] == "virtual_histology"
    assert dataset["adata"].shape == (2, 2)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_generated_dataset_survives_the_export_and_reload_round_trip(
    tile_workspace: tuple[str, str],
    deepspotm_stack: dict[str, object],
    reset_data_manager,
    tmp_path: Path,
):
    """Persistence stays with export_data, so its output has to be complete."""
    manifest_path, tile_directory = tile_workspace

    result = await predict_spatial_expression_from_histology(
        manifest_path,
        tile_directory,
        params=HistologyExpressionParameters(use_gpu=False),
    )
    target = tmp_path / "virtual.h5ad"
    await export_data(result.data_id, str(target))
    await reload_data(result.data_id, str(target))

    adata = (await data_manager.get_dataset(result.data_id))["adata"]
    assert adata.uns["chatspatial"]["expression"]["provenance"] == "predicted"
    assert adata.uns["deepspotm"]["slide_id"] == "SLIDE-1"
    assert adata.uns["deepspotm"]["weights_license"] == "CC-BY-NC-SA-4.0"
    assert adata.obsm["spatial"].shape == (2, 2)


@pytest.mark.integration
@pytest.mark.asyncio
async def test_a_failed_prediction_leaves_no_dataset_behind(
    tile_workspace: tuple[str, str],
    deepspotm_stack: dict[str, object],
    reset_data_manager,
):
    manifest_path, tile_directory = tile_workspace

    with pytest.raises(ParameterError):
        await predict_spatial_expression_from_histology(
            manifest_path,
            tile_directory,
            params=HistologyExpressionParameters(genes=["NOPE"], use_gpu=False),
        )

    assert await data_manager.list_datasets() == []


# =============================================================================
# Shared provenance guard through the public tools
# =============================================================================


def _predicted_dataset(n_obs: int = 24, n_vars: int = 12):
    rng = np.random.default_rng(0)
    adata = ad.AnnData(X=np.abs(rng.normal(size=(n_obs, n_vars))).astype(np.float32))
    adata.obs_names = [f"tile_{i}" for i in range(n_obs)]
    adata.var_names = [f"GENE{i}" for i in range(n_vars)]
    adata.obsm["spatial"] = rng.uniform(0, 1000, size=(n_obs, 2))
    set_expression_provenance(
        adata, provenance="predicted", units="log1p_cpm", producer="deepspotm"
    )
    return adata


@pytest.mark.integration
@pytest.mark.asyncio
async def test_preprocess_data_rejects_predicted_expression_with_guidance(
    reset_data_manager,
):
    data_id = await data_manager.create_dataset(
        _predicted_dataset(n_obs=8, n_vars=5), prefix="histology"
    )

    with pytest.raises(DataCompatibilityError) as excinfo:
        await preprocess_data(data_id)

    message = str(excinfo.value)
    assert "preprocess_data" in message
    assert "compute_embeddings" in message


@pytest.mark.integration
@pytest.mark.asyncio
async def test_compute_embeddings_accepts_predicted_expression(reset_data_manager):
    """The documented next step after prediction has to actually work."""
    data_id = await data_manager.create_dataset(
        _predicted_dataset(), prefix="histology"
    )

    result = await compute_embeddings(
        data_id, params=EmbeddingParameters(n_pcs=5, compute_umap=False)
    )

    assert result.data_id == data_id
    adata = (await data_manager.get_dataset(data_id))["adata"]
    assert "X_pca" in adata.obsm
    # The guard is provenance-aware, not destructive: the metadata survives.
    assert adata.uns["chatspatial"]["expression"]["provenance"] == "predicted"
