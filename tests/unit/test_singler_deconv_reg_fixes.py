"""Regression tests for deconvolution count validation and registration fixes."""

from __future__ import annotations

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# _prepare_counts with require_counts=True must raise DataError when only
# normalized (non-integer) data is available.
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_prepare_counts_require_counts_raises_for_normalized_data(
    minimal_spatial_adata,
):
    """_prepare_counts(require_counts=True) should raise DataError when the
    dataset contains only normalized (float) expression data."""
    from chatspatial.tools.deconvolution.base import _prepare_counts
    from chatspatial.utils.exceptions import DataError

    adata = minimal_spatial_adata.copy()
    # Replace X with clearly normalized (non-integer) data
    adata.X = (
        np.random.default_rng(42).uniform(0, 5, size=adata.shape).astype(np.float32)
    )
    # Remove any counts layer that might exist
    if "counts" in adata.layers:
        del adata.layers["counts"]

    class FakeCtx:
        async def warning(self, msg: str) -> None:
            pass

    ctx = FakeCtx()

    with pytest.raises(DataError, match="No raw integer counts found"):
        await _prepare_counts(
            adata,
            "Spatial",
            ctx,
            require_int_dtype=False,
            require_counts=True,
        )


@pytest.mark.asyncio
async def test_prepare_counts_require_counts_succeeds_for_integer_data(
    minimal_spatial_adata,
):
    """_prepare_counts(require_counts=True) should succeed when integer
    counts are available."""
    from chatspatial.tools.deconvolution.base import _prepare_counts

    adata = minimal_spatial_adata.copy()
    # Fixture X is already Poisson-distributed integers
    adata.X = np.round(adata.X).astype(np.float32)

    class FakeCtx:
        async def warning(self, msg: str) -> None:
            pass

    ctx = FakeCtx()

    result = await _prepare_counts(
        adata,
        "Spatial",
        ctx,
        require_int_dtype=False,
        require_counts=True,
    )
    assert result.shape == adata.shape


# ---------------------------------------------------------------------------
# Registration coordinate-system regressions
# ---------------------------------------------------------------------------


def test_stalign_image_transform_round_trip_preserves_xy_order():
    """Image row/column conversion must be reversible for asymmetric axes."""
    from chatspatial.tools import spatial_registration as reg

    coords = np.array(
        [[100.0, 1000.0], [140.0, 1150.0], [180.0, 1300.0]],
        dtype=np.float32,
    )
    transform = reg._ImageCoordinateTransform.from_coords(coords, (21, 41))

    image_coords = transform.to_image(coords)
    assert image_coords.shape == coords.shape
    assert not np.allclose(image_coords[:, 0], image_coords[:, 1])
    np.testing.assert_allclose(transform.from_image(image_coords), coords, atol=1e-5)


def test_stalign_target_image_transform_restores_target_coordinate_system():
    """Transformed image points must be decoded with the target slice scale."""
    from chatspatial.tools import spatial_registration as reg

    source = np.array([[0.0, 0.0], [10.0, 20.0]], dtype=np.float32)
    target = np.array([[100.0, 1000.0], [300.0, 1400.0]], dtype=np.float32)
    source_transform = reg._ImageCoordinateTransform.from_coords(source, (31, 51))
    target_transform = reg._ImageCoordinateTransform.from_coords(target, (31, 51))

    source_pixels = source_transform.to_image(source)
    restored = target_transform.from_image(source_pixels)
    np.testing.assert_allclose(restored, target, atol=1e-5)
