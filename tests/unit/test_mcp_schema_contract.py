"""MCP JSON Schema compatibility contracts."""

from __future__ import annotations

import pytest


def _collect_schema_array_issues(
    schema: object, path: tuple[str, ...] = ()
) -> list[str]:
    issues: list[str] = []
    if isinstance(schema, dict):
        if schema.get("type") == "array":
            if "items" not in schema:
                issues.append(f"{'/'.join(path)} is an array without items")
            if "prefixItems" in schema:
                issues.append(f"{'/'.join(path)} uses prefixItems")
        for key, value in schema.items():
            issues.extend(_collect_schema_array_issues(value, (*path, str(key))))
    elif isinstance(schema, list):
        for index, value in enumerate(schema):
            issues.extend(_collect_schema_array_issues(value, (*path, str(index))))
    return issues


@pytest.mark.unit
@pytest.mark.asyncio
async def test_exported_mcp_tool_input_schemas_are_anthropic_compatible():
    from chatspatial.server import mcp

    issues: list[str] = []
    for tool in await mcp.list_tools():
        issues.extend(
            _collect_schema_array_issues(
                tool.input_schema, path=(tool.name, "input_schema")
            )
        )

    assert not issues


@pytest.mark.unit
@pytest.mark.asyncio
async def test_all_top_level_tool_json_schemas_are_closed():
    from chatspatial.server import mcp

    for tool in await mcp.list_tools():
        assert tool.input_schema["additionalProperties"] is False, tool.name
        assert tool.output_schema["additionalProperties"] is False, tool.name


@pytest.mark.unit
@pytest.mark.asyncio
async def test_known_tuple_risk_tools_are_registered_and_schema_safe():
    from chatspatial.server import mcp

    tools = await mcp.list_tools()
    tool_names = {tool.name for tool in tools}
    risky_tools = {
        "visualize_data",
        "analyze_spatial_statistics",
        "analyze_trajectory_data",
        "register_spatial_data",
    }

    assert risky_tools <= tool_names

    tools_by_name = {tool.name: tool for tool in tools}
    for tool_name in risky_tools:
        issues = _collect_schema_array_issues(
            tools_by_name[tool_name].input_schema,
            path=(tool_name, "input_schema"),
        )
        assert not issues


@pytest.mark.unit
@pytest.mark.asyncio
async def test_analyze_enrichment_params_required_in_mcp_schema():
    from chatspatial.server import mcp

    tools_by_name = {tool.name: tool for tool in await mcp.list_tools()}
    schema = tools_by_name["analyze_enrichment"].input_schema

    assert "params" in schema["required"]

    params_schema = schema["properties"]["params"]
    assert "default" not in params_schema
    assert params_schema.get("type") != "null"
    assert not any(
        option.get("type") == "null"
        for option in params_schema.get("anyOf", [])
        if isinstance(option, dict)
    )


@pytest.mark.unit
@pytest.mark.asyncio
async def test_load_data_schema_enumerates_supported_spatial_platforms():
    from chatspatial.server import mcp

    tools_by_name = {tool.name: tool for tool in await mcp.list_tools()}
    data_type_schema = tools_by_name["load_data"].input_schema["properties"][
        "data_type"
    ]

    assert data_type_schema["enum"] == [
        "visium",
        "xenium",
        "slide_seq",
        "merfish",
        "seqfish",
        "generic",
    ]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_all_tools_publish_output_schemas_and_explicit_safety_hints():
    from chatspatial.server import mcp

    tools = await mcp.list_tools()

    assert len(tools) == 21
    assert all(tool.output_schema is not None for tool in tools)
    assert all(tool.annotations is not None for tool in tools)
    assert all(tool.annotations.read_only_hint is not None for tool in tools)
    assert all(tool.annotations.destructive_hint is not None for tool in tools)
    assert all(tool.annotations.idempotent_hint is not None for tool in tools)
    assert all(tool.annotations.open_world_hint is not None for tool in tools)


@pytest.mark.unit
@pytest.mark.asyncio
async def test_embedding_and_registration_outputs_have_closed_schemas():
    from chatspatial.server import mcp

    tools = {tool.name: tool for tool in await mcp.list_tools()}

    embedding_schema = tools["compute_embeddings"].output_schema
    registration_schema = tools["register_spatial_data"].output_schema

    assert embedding_schema["additionalProperties"] is False
    assert registration_schema["additionalProperties"] is False
    assert set(embedding_schema["required"]) == {"data_id", "computed", "skipped"}
    assert set(registration_schema["required"]) == {
        "method",
        "source_id",
        "target_id",
        "n_source_spots",
        "n_target_spots",
        "spatial_key_registered",
    }
    assert embedding_schema["properties"]["computed"]["items"]["type"] == "string"
    assert registration_schema["properties"]["method"]["enum"] == [
        "paste",
        "stalign",
    ]


@pytest.mark.unit
@pytest.mark.asyncio
async def test_internal_result_fields_are_absent_from_mcp_output_schemas():
    from chatspatial.server import mcp

    tools = {tool.name: tool for tool in await mcp.list_tools()}
    excluded_by_tool = {
        "find_markers": {"statistics"},
        "analyze_spatial_statistics": {"statistics"},
        "deconvolve_data": {"statistics"},
        "identify_spatial_domains": {"statistics"},
        "analyze_cell_communication": {"statistics"},
        "analyze_cnv": {"statistics"},
        "analyze_enrichment": {
            "enrichment_scores",
            "pvalues",
            "adjusted_pvalues",
            "gene_set_statistics",
            "spatial_metrics",
        },
        "find_spatial_genes": {
            "gene_statistics",
            "p_values",
            "q_values",
            "spatialde_results",
            "sparkx_results",
            "flashs_results",
        },
    }

    for tool_name, internal_fields in excluded_by_tool.items():
        public_fields = set(tools[tool_name].output_schema["properties"])
        assert public_fields.isdisjoint(internal_fields), tool_name


@pytest.mark.unit
@pytest.mark.asyncio
async def test_all_public_parameter_models_reject_unknown_fields():
    from chatspatial.server import mcp

    unchecked_models: list[str] = []
    for tool in await mcp.list_tools():
        for model_name, schema in tool.input_schema.get("$defs", {}).items():
            if (
                model_name.endswith("Parameters")
                and schema.get("additionalProperties") is not False
            ):
                unchecked_models.append(f"{tool.name}.{model_name}")

    assert not unchecked_models


@pytest.mark.unit
@pytest.mark.asyncio
async def test_preprocessing_contract_routes_embedding_work_to_explicit_tool():
    from chatspatial.server import mcp

    tools_by_name = {tool.name: tool for tool in await mcp.list_tools()}
    preprocessing = tools_by_name["preprocess_data"]

    assert "does not compute PCA, UMAP, clustering, or neighbor graphs" in (
        preprocessing.description or ""
    )
    assert "compute_embeddings" in (preprocessing.description or "")
    assert "compute_embeddings" in (mcp.instructions or "")
    assert mcp.instructions.index("preprocess_data") < mcp.instructions.index(
        "compute_embeddings"
    )
    assert (
        "compute_embeddings"
        in preprocessing.output_schema["properties"]["clusters"]["description"]
    )

    preprocessing_fields = preprocessing.input_schema["$defs"][
        "PreprocessingParameters"
    ]["properties"]
    embedding_only_fields = {
        "n_pcs",
        "n_neighbors",
        "clustering_resolution",
        "clustering_key",
    }
    assert embedding_only_fields.isdisjoint(preprocessing_fields)

    embedding_fields = tools_by_name["compute_embeddings"].input_schema["$defs"][
        "EmbeddingParameters"
    ]["properties"]
    assert embedding_only_fields <= embedding_fields.keys()


@pytest.mark.unit
@pytest.mark.asyncio
async def test_integration_schema_requires_multiple_dataset_ids():
    from chatspatial.server import mcp

    tools_by_name = {tool.name: tool for tool in await mcp.list_tools()}
    data_ids = tools_by_name["integrate_samples"].input_schema["properties"]["data_ids"]

    assert data_ids["minItems"] == 2
    assert "distinct" in data_ids["description"]
