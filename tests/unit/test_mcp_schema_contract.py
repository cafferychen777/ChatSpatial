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
def test_exported_mcp_tool_input_schemas_are_anthropic_compatible():
    from chatspatial.server import mcp

    issues: list[str] = []
    for tool in mcp._tool_manager.list_tools():
        issues.extend(
            _collect_schema_array_issues(
                tool.parameters, path=(tool.name, "parameters")
            )
        )

    assert not issues


@pytest.mark.unit
def test_known_tuple_risk_tools_are_registered_and_schema_safe():
    from chatspatial.server import mcp

    tool_names = {tool.name for tool in mcp._tool_manager.list_tools()}
    risky_tools = {
        "visualize_data",
        "analyze_spatial_statistics",
        "analyze_trajectory_data",
        "register_spatial_data",
    }

    assert risky_tools <= tool_names

    tools_by_name = {tool.name: tool for tool in mcp._tool_manager.list_tools()}
    for tool_name in risky_tools:
        issues = _collect_schema_array_issues(
            tools_by_name[tool_name].parameters, path=(tool_name, "parameters")
        )
        assert not issues


@pytest.mark.unit
def test_analyze_enrichment_params_required_in_mcp_schema():
    from chatspatial.server import mcp

    tools_by_name = {tool.name: tool for tool in mcp._tool_manager.list_tools()}
    schema = tools_by_name["analyze_enrichment"].parameters

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
def test_load_data_schema_enumerates_supported_spatial_platforms():
    from chatspatial.server import mcp

    tools_by_name = {tool.name: tool for tool in mcp._tool_manager.list_tools()}
    data_type_schema = tools_by_name["load_data"].parameters["properties"]["data_type"]

    assert data_type_schema["enum"] == [
        "visium",
        "xenium",
        "slide_seq",
        "merfish",
        "seqfish",
        "generic",
    ]
