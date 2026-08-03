"""Unit contracts for compatibility patch utilities."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from chatspatial.utils import compat


def test_numpy2_compat_assert_array_equal_works_with_legacy_keywords():
    with compat.numpy2_compat():
        np.testing.assert_array_equal(x=np.array([1, 2]), y=np.array([1, 2]))


def test_numpy2_compat_wrapper_preserves_modern_positional_arguments():
    received = []

    def original(*args, **kwargs):
        received.append((args, kwargs))

    wrapper = compat._make_numpy2_compat_wrapper(original)
    actual = np.array([1])
    desired = np.array([1])

    wrapper(actual, desired, "detail", False, strict=True)

    assert received == [((actual, desired, "detail", False), {"strict": True})]


def test_numpy2_compat_overlapping_scopes_restore_after_last_owner():
    original = np.testing.assert_array_equal
    first_cleanup = compat._patch_numpy_testing()
    wrapper = np.testing.assert_array_equal
    second_cleanup = compat._patch_numpy_testing()

    first_cleanup()
    assert np.testing.assert_array_equal is wrapper
    np.testing.assert_array_equal(x=np.array([1]), y=np.array([1]))

    second_cleanup()
    assert np.testing.assert_array_equal is original


def test_numpy2_compat_unpatch_leaves_no_residual_state():
    """After unpatch, np.testing.assert_array_equal has no compat marker."""
    original = np.testing.assert_array_equal
    with compat.numpy2_compat():
        assert getattr(
            np.testing.assert_array_equal, compat._NUMPY2_COMPAT_MARKER, False
        )
    # After context exit, original is restored with no leaked attributes
    assert np.testing.assert_array_equal is original
    assert not getattr(
        np.testing.assert_array_equal, compat._NUMPY2_COMPAT_MARKER, False
    )


def test_derivative_compat_first_and_second_order():
    def f(x):
        return x**2

    assert compat._derivative_compat(f, 3.0, dx=1e-4, n=1, order=3) == pytest.approx(
        6.0, rel=1e-3
    )
    assert compat._derivative_compat(f, 3.0, dx=1e-4, n=2, order=3) == pytest.approx(
        2.0, rel=1e-3
    )


def test_derivative_compat_validates_order():
    def f(x):
        return x

    with pytest.raises(ValueError, match="must be at least"):
        compat._derivative_compat(f, 1.0, n=2, order=2)
    with pytest.raises(ValueError, match="must be odd"):
        compat._derivative_compat(f, 1.0, n=1, order=4)


def test_patch_scipy_sparse_matrix_A_restores_alias():
    compat.patch_scipy_sparse_matrix_A()
    arr = csr_matrix([[1, 0], [0, 2]])
    assert hasattr(arr, "A")
    np.testing.assert_array_equal(arr.A, arr.toarray())


def test_cellrank_compat_decorator_runs_cleanup_on_error(
    monkeypatch: pytest.MonkeyPatch,
):
    state = {"cleaned": False}

    def _fake_ensure():
        return lambda: state.__setitem__("cleaned", True)

    monkeypatch.setattr(compat, "ensure_cellrank_compat", _fake_ensure)

    @compat.cellrank_compat
    def _boom():
        raise RuntimeError("x")

    with pytest.raises(RuntimeError):
        _boom()
    assert state["cleaned"]


def test_ensure_spatialde_compat_calls_both_patches(monkeypatch: pytest.MonkeyPatch):
    called = {"misc": False, "aliases": False}
    monkeypatch.setattr(
        compat, "patch_scipy_misc_derivative", lambda: called.__setitem__("misc", True)
    )
    monkeypatch.setattr(
        compat, "patch_scipy_numpy_aliases", lambda: called.__setitem__("aliases", True)
    )

    compat.ensure_spatialde_compat()
    assert called["misc"] and called["aliases"]


def test_get_compatibility_info_and_derivative_status_have_expected_keys():
    info = compat.get_compatibility_info()
    status = compat.check_scipy_derivative_status()

    assert "known_issues" in info
    assert "patches_available" in info
    assert "has_derivative" in status
    assert "needs_patch" in status
