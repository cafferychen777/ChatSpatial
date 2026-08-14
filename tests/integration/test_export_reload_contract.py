"""Integration contract tests for export_data and reload_data tools."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import anndata as ad
import pytest

from chatspatial.server import data_manager, export_data, reload_data
from chatspatial.utils.exceptions import DataError, DataNotFoundError
from tests.fixtures.helpers import load_generic_dataset


def _extract_export_path(msg: str) -> Path:
    # "Dataset 'data_1' exported to: /path/to/file.h5ad"
    return Path(msg.split(" exported to: ", 1)[1].strip())


class DummyContext:
    """Minimal MCP context mock for progress-notification contracts."""

    def __init__(self):
        self.progress_updates: list[tuple[float, float | None, str | None]] = []

    async def report_progress(
        self,
        progress: float,
        total: float | None = None,
        message: str | None = None,
    ) -> None:
        self.progress_updates.append((progress, total, message))


@pytest.mark.integration
@pytest.mark.asyncio
async def test_export_data_writes_file_and_returns_path(
    spatial_dataset_path, tmp_path, reset_data_manager
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="export_contract")

    target = tmp_path / "exports" / "dataset_export.h5ad"
    msg = await export_data(dataset.id, path=str(target))

    out = _extract_export_path(msg)
    assert out.exists()
    assert out == target.resolve()


@pytest.mark.integration
@pytest.mark.asyncio
async def test_reload_data_applies_external_changes(
    spatial_dataset_path, tmp_path, reset_data_manager
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="reload_contract")

    target = tmp_path / "exports" / "dataset_reload.h5ad"
    await export_data(dataset.id, path=str(target))

    # Simulate external script modification: subset genes and write back
    adata = ad.read_h5ad(target)
    adata_mod = adata[:, :7].copy()
    adata_mod.write_h5ad(target)

    msg = await reload_data(dataset.id, path=str(target))
    assert "reloaded" in msg
    assert "× 7 genes" in msg

    stored = await data_manager.get_dataset(dataset.id)
    assert stored["adata"].n_vars == 7


@pytest.mark.integration
@pytest.mark.asyncio
async def test_reload_data_emits_progress_on_success(
    spatial_dataset_path, tmp_path, reset_data_manager
):
    dataset = await load_generic_dataset(
        spatial_dataset_path, name="reload_ctx_success"
    )

    target = tmp_path / "exports" / "dataset_reload_ctx.h5ad"
    await export_data(dataset.id, path=str(target))

    adata = ad.read_h5ad(target)
    adata_mod = adata[:, :5].copy()
    adata_mod.write_h5ad(target)

    ctx = DummyContext()
    msg = await reload_data(dataset.id, path=str(target), context=ctx)

    assert "reloaded" in msg
    assert ctx.progress_updates == [
        (1.0, None, f"Reloading dataset '{dataset.id}'..."),
        (2.0, None, f"Dataset '{dataset.id}' reloaded successfully"),
    ]


@pytest.mark.integration
@pytest.mark.asyncio
async def test_reload_data_missing_file_raises_file_not_found(
    spatial_dataset_path, tmp_path, reset_data_manager
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="reload_missing")

    missing = tmp_path / "does_not_exist.h5ad"
    with pytest.raises(FileNotFoundError, match="Data file not found"):
        await reload_data(dataset.id, path=str(missing))


@pytest.mark.integration
@pytest.mark.asyncio
async def test_export_data_missing_dataset_raises_not_found(reset_data_manager):
    with pytest.raises(DataNotFoundError, match="Dataset absent_data not found"):
        await export_data("absent_data")


@pytest.mark.integration
@pytest.mark.asyncio
async def test_export_data_emits_progress(
    spatial_dataset_path, tmp_path, reset_data_manager
):
    dataset = await load_generic_dataset(spatial_dataset_path, name="export_ctx")
    target = tmp_path / "exports" / "ctx_export.h5ad"
    ctx = DummyContext()

    msg = await export_data(dataset.id, path=str(target), context=ctx)

    assert "exported to:" in msg
    assert ctx.progress_updates[0] == (
        1.0,
        None,
        f"Exporting dataset '{dataset.id}'...",
    )
    assert ctx.progress_updates[1][0:2] == (2.0, None)
    assert ctx.progress_updates[1][2] is not None
    assert ctx.progress_updates[1][2].startswith("Dataset exported to:")


@pytest.mark.integration
@pytest.mark.asyncio
async def test_export_data_does_not_emit_duplicate_error_status(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    class BoomError(RuntimeError):
        pass

    async def fake_get_dataset(data_id: str):
        class FakeAdata:
            n_obs = 1
            n_vars = 1

        return {"adata": FakeAdata()}

    def fake_export_adata(data_id, adata, path):
        raise BoomError("boom")

    monkeypatch.setattr(data_manager, "get_dataset", fake_get_dataset)
    monkeypatch.setitem(
        sys.modules,
        "chatspatial.utils.persistence",
        SimpleNamespace(export_adata=fake_export_adata),
    )

    ctx = DummyContext()
    with pytest.raises(BoomError, match="boom"):
        await export_data("d_boom", context=ctx)

    assert ctx.progress_updates == [(1.0, None, "Exporting dataset 'd_boom'...")]


@pytest.mark.integration
@pytest.mark.asyncio
async def test_reload_data_does_not_emit_duplicate_missing_file_error_status(
    spatial_dataset_path, tmp_path, reset_data_manager
):
    dataset = await load_generic_dataset(
        spatial_dataset_path, name="reload_ctx_missing"
    )
    ctx = DummyContext()
    missing = tmp_path / "missing_reload_file.h5ad"

    with pytest.raises(FileNotFoundError, match="Data file not found"):
        await reload_data(dataset.id, path=str(missing), context=ctx)

    assert ctx.progress_updates == [
        (1.0, None, f"Reloading dataset '{dataset.id}'..."),
    ]


@pytest.mark.integration
@pytest.mark.asyncio
async def test_reload_data_does_not_emit_duplicate_generic_error_status(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    def fake_load_adata_from_active(data_id, path):
        raise RuntimeError("unexpected reload failure")

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.utils.persistence",
        SimpleNamespace(
            load_adata_from_active=fake_load_adata_from_active,
            check_reload_identity=lambda *_args, **_kwargs: None,
            get_active_path=lambda data_id: Path(f"/tmp/{data_id}.h5ad"),
        ),
    )

    ctx = DummyContext()
    with pytest.raises(RuntimeError, match="unexpected reload failure"):
        await reload_data("d_fail", path="/tmp/not_used.h5ad", context=ctx)

    assert ctx.progress_updates == [(1.0, None, "Reloading dataset 'd_fail'...")]


@pytest.mark.integration
@pytest.mark.asyncio
async def test_reload_data_refuses_file_holding_a_different_dataset(
    spatial_dataset_path, tmp_path, reset_data_manager
):
    """Dataset IDs restart every session, so the active file may be someone else's.

    Reloading it would silently swap the live dataset for unrelated data.
    """
    dataset = await load_generic_dataset(spatial_dataset_path, name="reload_identity")

    target = tmp_path / "exports" / "stale.h5ad"
    await export_data(dataset.id, path=str(target))

    stale = ad.read_h5ad(target)
    stale.obs_names = [f"other_{i}" for i in range(stale.n_obs)]
    stale.write_h5ad(target)

    with pytest.raises(DataError, match="different dataset"):
        await reload_data(dataset.id, path=str(target))

    stored = await data_manager.get_dataset(dataset.id)
    assert list(stored["adata"].obs_names)[0] != "other_0"
