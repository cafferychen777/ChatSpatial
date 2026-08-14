"""Contracts for optional dependency families and shared constraints."""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from packaging.requirements import Requirement

from chatspatial.models.data import CellCommunicationParameters
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
    fastccc = _requirements(extras["fastccc"])
    r_backends = _requirements(extras["r-backends"])

    assert "jinja2" not in core
    assert "cellrank" not in full
    assert "pygam" not in full
    assert full["spatialde-modern"][0].specifier.contains("1.1.3.post1")
    assert "spatialde" not in full
    assert {"cellrank", "palantir", "pygam", "scipy"} <= trajectory.keys()
    assert "fastccc" not in full
    assert "fastccc" not in communication
    assert "rpy2" not in full
    assert "anndata2ri" not in full
    assert {"rpy2", "anndata2ri"} <= r_backends.keys()
    assert r_backends["rpy2"][0].specifier.contains("3.6.7")
    liana = communication["liana"][0]
    assert liana.specifier.contains("1.4.0")
    assert liana.marker is not None
    assert liana.marker.evaluate({"python_version": "3.13"})
    assert not liana.marker.evaluate({"python_version": "3.14"})
    assert communication["ktplotspy"][0].specifier.contains("0.3.5")
    requirement = fastccc["fastccc"][0]
    assert requirement.specifier.contains("1.0.1")
    assert requirement.marker is not None
    assert requirement.marker.evaluate({"python_version": "3.12"})
    assert not requirement.marker.evaluate({"python_version": "3.13"})


@pytest.mark.unit
def test_default_flashs_backend_is_a_core_dependency() -> None:
    """A clean default install must support the default spatial gene method."""
    project = _project_metadata()
    core = _requirements(project["dependencies"])

    assert "flashs" in core
    requirement = core["flashs"][0]
    assert requirement.specifier.contains("0.1.1")
    assert not requirement.specifier.contains("0.2.0")


@pytest.mark.unit
def test_direct_core_imports_have_explicit_dependencies() -> None:
    """Core code must not rely on optional or transitive installations."""
    project = _project_metadata()
    core = _requirements(project["dependencies"])

    assert {"statsmodels", "typing-extensions", "zarr"} <= core.keys()


@pytest.mark.unit
def test_unused_system_stacks_are_not_advertised_as_extras() -> None:
    """Extras should correspond to code paths owned by ChatSpatial."""
    project = _project_metadata()
    extras = project["optional-dependencies"]

    assert "hpc" not in extras
    assert all(
        "pysal" not in _requirements(requirements) for requirements in extras.values()
    )


@pytest.mark.unit
def test_spatial_domain_extra_uses_maintained_spagcn_graph_stack() -> None:
    """SpaGCN must use Leiden without pulling the obsolete Louvain extension."""
    project = _project_metadata()
    requirements = _requirements(project["optional-dependencies"]["spatial-domains"])

    assert requirements["spagcn-modern"][0].specifier.contains("1.2.7.post2")
    banksy = requirements["pybanksy"][0]
    assert banksy.specifier.contains("1.3.5")
    assert banksy.marker is not None
    assert banksy.marker.evaluate({"python_version": "3.13"})
    assert not banksy.marker.evaluate({"python_version": "3.14"})
    assert "python-igraph" not in requirements
    assert "igraph" not in requirements
    assert "leidenalg" not in requirements
    assert "louvain" not in requirements


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
def test_install_guidance_uses_compatible_optional_families() -> None:
    """Runtime errors should direct users to curated extras, not raw packages."""
    expected_extras = {
        "cellrank": "trajectory",
        "palantir": "trajectory",
        "pygam": "trajectory",
        "SpaGCN": "spatial-domains",
        "banksy": "spatial-domains",
        "aestetik": "aestetik",
        "fastccc": "fastccc",
        "rpy2": "r-backends",
        "anndata2ri": "r-backends",
    }

    for dependency, extra in expected_extras.items():
        assert DEPENDENCY_REGISTRY[dependency].install_cmd == (
            f"pip install 'chatspatial[{extra}]'"
        )

    assert DEPENDENCY_REGISTRY["louvain"].install_cmd == "pip install louvain"


@pytest.mark.unit
def test_default_communication_backend_is_in_the_composable_extra() -> None:
    """Default parameters must not select the isolated FastCCC backend."""
    params = CellCommunicationParameters(species="human", cell_type_key="cell_type")

    assert params.method == "liana"


@pytest.mark.unit
def test_dev_extra_includes_request_type_stubs() -> None:
    """The declared mypy environment must include third-party request stubs."""
    project = _project_metadata()
    dev = _requirements(project["optional-dependencies"]["dev"])

    assert "types-requests" in dev
