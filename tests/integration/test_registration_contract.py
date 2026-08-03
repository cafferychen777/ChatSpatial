"""Integration contract tests for spatial registration tool entry."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest
from pydantic import ValidationError

from chatspatial.models.analysis import SpatialRegistrationResult
from chatspatial.models.data import RegistrationParameters
from chatspatial.server import register_spatial_data
from tests.fixtures.helpers import load_generic_dataset


@pytest.mark.integration
@pytest.mark.asyncio
async def test_register_spatial_data_invalid_method_raises_validation_error(
    spatial_dataset_path, reset_data_manager
):
    source = await load_generic_dataset(spatial_dataset_path, name="reg_src")
    target = await load_generic_dataset(spatial_dataset_path, name="reg_tgt")

    with pytest.raises(ValidationError, match="Input should be 'paste' or 'stalign'"):
        await register_spatial_data(
            source.id, target.id, params=RegistrationParameters(method="invalid_method")
        )


@pytest.mark.integration
@pytest.mark.asyncio
async def test_register_spatial_data_success_returns_registration_result(
    reset_data_manager, monkeypatch: pytest.MonkeyPatch
):
    async def fake_register(source_id, target_id, ctx, params=None):
        return SpatialRegistrationResult(
            source_id=source_id,
            target_id=target_id,
            method=params.method if params else "paste",
            n_source_spots=10,
            n_target_spots=12,
            spatial_key_registered="spatial_registered",
        )

    monkeypatch.setitem(
        sys.modules,
        "chatspatial.tools.spatial_registration",
        SimpleNamespace(register_spatial_slices_mcp=fake_register),
    )

    result = await register_spatial_data(
        "source_1", "target_1", params=RegistrationParameters(method="paste")
    )

    assert result.source_id == "source_1"
    assert result.target_id == "target_1"
    assert result.method == "paste"
