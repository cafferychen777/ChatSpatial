"""Contracts for optional dependency families and shared constraints."""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from packaging.requirements import Requirement

from chatspatial.utils.dependency_manager import DEPENDENCY_REGISTRY

REPO_ROOT = Path(__file__).resolve().parents[2]


def _project_metadata() -> dict[str, object]:
    with (REPO_ROOT / "pyproject.toml").open("rb") as file:
        return tomllib.load(file)["project"]


def _requirements(items: list[str]) -> dict[str, list[Requirement]]:
    parsed: dict[str, list[Requirement]] = {}
    for item in items:
        requirement = Requirement(item)
        parsed.setdefault(requirement.name.lower(), []).append(requirement)
    return parsed


@pytest.mark.unit
def test_dependency_families_keep_conflicting_backends_separate() -> None:
    """PyPI extras must remain independently resolvable."""
    project = _project_metadata()
    extras = project["optional-dependencies"]

    core = _requirements(project["dependencies"])
    full = _requirements(extras["full"])
    trajectory = _requirements(extras["trajectory"])
    communication = _requirements(extras["cell-communication"])

    assert "jinja2" not in core
    assert "cellrank" not in full
    assert "pygam" not in full
    assert {"cellrank", "palantir", "pygam", "scipy"} <= trajectory.keys()
    assert communication["fastccc"][0].specifier.contains("1.0.1")
    assert communication["ktplotspy"][0].specifier.contains("0.3.5")


@pytest.mark.unit
def test_spatial_domain_extra_pins_a_coherent_graph_stack() -> None:
    """Louvain and Leiden must resolve against the same igraph generation."""
    project = _project_metadata()
    requirements = _requirements(project["optional-dependencies"]["spatial-domains"])

    expected_versions = {
        "spagcn": "1.2.7",
        "louvain": "0.8.2",
        "igraph": "0.11.9",
        "python-igraph": "0.11.9",
        "leidenalg": "0.10.2",
        "pybanksy": "1.3.5",
    }
    for name, version in expected_versions.items():
        assert requirements[name][0].specifier.contains(version)

    assert not requirements["igraph"][0].specifier.contains("1.0.0")
    assert not requirements["leidenalg"][0].specifier.contains("0.11.0")


@pytest.mark.unit
def test_aestetik_extra_is_bounded_and_excludes_python_314() -> None:
    """AESTETIK must stay isolated from full and unsupported Python releases."""
    project = _project_metadata()
    extras = project["optional-dependencies"]
    requirement = _requirements(extras["aestetik"])["aestetik"][0]

    assert requirement.specifier.contains("0.3.1")
    assert not requirement.specifier.contains("0.4.0")
    assert requirement.marker is not None
    assert requirement.marker.evaluate({"python_version": "3.13"})
    assert not requirement.marker.evaluate({"python_version": "3.14"})
    assert "aestetik" not in _requirements(extras["full"])


@pytest.mark.unit
def test_deepspotm_extra_is_bounded_and_excluded_from_full() -> None:
    """A non-commercial backend must stay an explicit, isolated opt-in."""
    project = _project_metadata()
    extras = project["optional-dependencies"]
    requirement = _requirements(extras["deepspotm"])["deepspotm"][0]

    assert requirement.specifier.contains("1.0.0")
    assert not requirement.specifier.contains("2.0.0")
    assert requirement.marker is not None
    assert requirement.marker.evaluate({"python_version": "3.13"})
    assert not requirement.marker.evaluate({"python_version": "3.14"})

    assert "deepspotm" not in _requirements(extras["full"])
    assert "deepspotm" not in _requirements(project["dependencies"])


@pytest.mark.unit
def test_deepspotm_guidance_states_the_license_and_weight_gate() -> None:
    """Runtime guidance must name the non-commercial terms and the gate."""
    info = DEPENDENCY_REGISTRY["deepspotm"]

    assert info.install_cmd == "pip install 'chatspatial[deepspotm]'"
    assert "Non-commercial" in info.description
    assert "PolyForm Noncommercial 1.0.0" in info.description
    assert "CC-BY-NC-SA-4.0" in info.description
    assert "huggingface-cli login" in info.description


@pytest.mark.unit
def test_shared_environment_constraints_cover_conflict_boundaries() -> None:
    """The combined repository environment must pin every conflict boundary."""
    constraint_lines = {
        line.strip()
        for line in (REPO_ROOT / "constraints/shared-py312.txt")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip() and not line.startswith("#")
    }

    expected_prefixes = {
        "cellrank==",
        "fastccc==",
        "igraph==",
        "ktplotspy==",
        "leidenalg==",
        "louvain==",
        "pybanksy==",
        "pygpcca @ git+https://github.com/msmdev/pyGPCCA.git@",
        "python-igraph==",
        "scipy==",
    }
    for prefix in expected_prefixes:
        assert any(line.startswith(prefix) for line in constraint_lines)


@pytest.mark.unit
def test_install_guidance_uses_compatible_optional_families() -> None:
    """Runtime errors should direct users to curated extras, not raw packages."""
    expected_extras = {
        "cellrank": "trajectory",
        "palantir": "trajectory",
        "pygam": "trajectory",
        "SpaGCN": "spatial-domains",
        "banksy": "spatial-domains",
        "louvain": "spatial-domains",
        "aestetik": "aestetik",
        "deepspotm": "deepspotm",
    }

    for dependency, extra in expected_extras.items():
        assert DEPENDENCY_REGISTRY[dependency].install_cmd == (
            f"pip install 'chatspatial[{extra}]'"
        )


@pytest.mark.unit
def test_dev_extra_includes_request_type_stubs() -> None:
    """The declared mypy environment must include third-party request stubs."""
    project = _project_metadata()
    dev = _requirements(project["optional-dependencies"]["dev"])

    assert "types-requests" in dev
