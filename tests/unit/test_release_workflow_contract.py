"""Structural contracts for the release workflow dependency graph."""

from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CI_WORKFLOW = REPOSITORY_ROOT / ".github" / "workflows" / "ci.yml"
RELEASE_WORKFLOW = REPOSITORY_ROOT / ".github" / "workflows" / "release.yml"
DOCKER_WORKFLOW = REPOSITORY_ROOT / ".github" / "workflows" / "docker.yml"


def test_release_workflow_uses_one_explicit_quality_gate() -> None:
    """Release artifacts must be built inside the tagged release run."""
    ci = CI_WORKFLOW.read_text(encoding="utf-8")
    release = RELEASE_WORKFLOW.read_text(encoding="utf-8")

    assert "workflow_call:" in ci
    assert "build_release_distributions:" in ci
    assert "uses: ./.github/workflows/ci.yml" in release
    assert "build_release_distributions: true" in release
    assert "needs: provenance" in release
    assert "needs: quality" in release


def test_release_workflow_has_no_cross_run_discovery() -> None:
    """Release ordering must use job dependencies, never eventual API state."""
    release = RELEASE_WORKFLOW.read_text(encoding="utf-8")

    assert "workflow_dispatch:" not in release
    forbidden_fragments = (
        "gh run list",
        "ci_run_id",
        "run-id:",
        "github-token:",
        "actions: read",
        "workflow_run",
    )
    for fragment in forbidden_fragments:
        assert fragment not in release


def test_release_workflow_owns_docker_publication() -> None:
    """A release must build one image only after its public GitHub Release."""
    release = RELEASE_WORKFLOW.read_text(encoding="utf-8")
    docker = DOCKER_WORKFLOW.read_text(encoding="utf-8")

    assert "uses: ./.github/workflows/docker.yml" in release
    assert "needs: publish-github" in release
    assert "workflow_call:" in docker
    assert "\n  push:\n" not in docker
    assert "github.ref_type == 'tag'" in docker
    assert "group: docker-${{ github.ref }}" in docker
    assert "cancel-in-progress: false" in docker
