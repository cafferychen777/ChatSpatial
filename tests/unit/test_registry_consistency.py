"""Structural consistency tests: verify registries match parameter model Literals.

These tests ensure that runtime dispatch registries stay in sync with the
Pydantic Literal types that define valid values in the MCP JSON schema.
If a method is added to a registry but not to the Literal (or vice versa),
these tests fail immediately.
"""

from __future__ import annotations

import json
import tomllib
from pathlib import Path

import pytest
from packaging.specifiers import SpecifierSet


@pytest.mark.unit
def test_runtime_version_matches_pyproject():
    """Runtime package version should match pyproject.toml."""
    import chatspatial

    pyproject_path = Path(__file__).resolve().parents[2] / "pyproject.toml"
    with pyproject_path.open("rb") as f:
        expected_version = tomllib.load(f)["project"]["version"]

    assert chatspatial.__version__ == expected_version


@pytest.mark.unit
def test_mcp_registry_versions_match_pyproject():
    """MCP Registry metadata versions should match pyproject.toml."""
    repo_root = Path(__file__).resolve().parents[2]
    with (repo_root / "pyproject.toml").open("rb") as f:
        expected_version = tomllib.load(f)["project"]["version"]
    with (repo_root / "server.json").open("r", encoding="utf-8") as f:
        server = json.load(f)

    assert server["version"] == expected_version
    assert all(package["version"] == expected_version for package in server["packages"])


@pytest.mark.unit
def test_mcp_registry_pypi_package_is_uvx_runnable():
    """Registry clients should be able to derive the zero-environment command."""
    repo_root = Path(__file__).resolve().parents[2]
    with (repo_root / "server.json").open("r", encoding="utf-8") as f:
        server = json.load(f)

    pypi_packages = [
        package for package in server["packages"] if package["registryType"] == "pypi"
    ]
    assert len(pypi_packages) == 1

    package = pypi_packages[0]
    assert package["runtimeHint"] == "uvx"
    assert package["packageArguments"] == [{"type": "positional", "value": "server"}]
    assert package["transport"] == {"type": "stdio"}


@pytest.mark.unit
def test_docs_build_python_satisfies_package_requires_python():
    """Read the Docs build Python must satisfy the package support policy."""
    repo_root = Path(__file__).resolve().parents[2]

    with (repo_root / "pyproject.toml").open("rb") as f:
        requires_python = tomllib.load(f)["project"]["requires-python"]
    with (repo_root / ".readthedocs.yaml").open("r", encoding="utf-8") as f:
        rtd_config = f.read()

    declared_python = next(
        line.split('"')[1]
        for line in rtd_config.splitlines()
        if "python:" in line and '"' in line
    )

    assert f"{declared_python}.0" in SpecifierSet(requires_python)


@pytest.mark.unit
def test_deconvolution_registry_matches_literal():
    """METHOD_REGISTRY keys must equal DeconvolutionParameters.method Literal."""
    from chatspatial.models.data import DeconvolutionParameters
    from chatspatial.tools.deconvolution import METHOD_REGISTRY

    literal_methods = set(
        DeconvolutionParameters.model_fields["method"].annotation.__args__
    )
    assert set(METHOD_REGISTRY) == literal_methods


@pytest.mark.unit
def test_spatial_statistics_registry_matches_literal():
    """_ANALYSIS_REGISTRY keys must equal SpatialStatisticsParameters.analysis_type Literal."""
    from chatspatial.models.data import SpatialStatisticsParameters
    from chatspatial.tools.spatial_statistics import _ANALYSIS_REGISTRY

    literal_types = set(
        SpatialStatisticsParameters.model_fields["analysis_type"].annotation.__args__
    )
    assert set(_ANALYSIS_REGISTRY) == literal_types


@pytest.mark.unit
@pytest.mark.asyncio
async def test_every_tool_has_annotations():
    """Every registered MCP tool must have ToolAnnotations (inline, not a dict)."""
    from chatspatial.server import mcp

    tools = await mcp.list_tools()
    missing = [t.name for t in tools if t.annotations is None]
    assert not missing, f"Tools without annotations: {missing}"
