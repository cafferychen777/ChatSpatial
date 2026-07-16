"""Unit tests for shared statistical post-processing."""

from types import ModuleType

import numpy as np
import pytest

from chatspatial.utils import statistics


def test_adjust_pvalues_preserves_invalid_positions(
    monkeypatch: pytest.MonkeyPatch,
):
    captured: dict[str, object] = {}
    multitest = ModuleType("statsmodels.stats.multitest")

    def _multipletests(values, **kwargs):
        captured["values"] = np.asarray(values)
        captured.update(kwargs)
        return np.array([True, False]), np.array([0.02, 0.4]), 0.0, 0.0

    multitest.multipletests = _multipletests
    monkeypatch.setattr(
        statistics,
        "require_module",
        lambda *_args, **_kwargs: multitest,
    )

    reject, adjusted = statistics.adjust_pvalues(
        [0.01, np.nan, 0.2, np.inf],
        alpha=0.05,
        method="fdr_bh",
        feature="test correction",
    )

    np.testing.assert_array_equal(reject, [True, False, False, False])
    np.testing.assert_allclose(adjusted[[0, 2]], [0.02, 0.4])
    assert np.isnan(adjusted[[1, 3]]).all()
    np.testing.assert_allclose(captured["values"], [0.01, 0.2])
    assert captured["method"] == "fdr_bh"
    assert captured["is_sorted"] is False
    assert captured["returnsorted"] is False


def test_adjust_pvalues_skips_dependency_when_no_finite_values(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(
        statistics,
        "require_module",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("empty correction family must not load statsmodels")
        ),
    )

    reject, adjusted = statistics.adjust_pvalues([np.nan, np.inf])

    assert not reject.any()
    assert np.isnan(adjusted).all()


def test_adjust_pvalues_rejects_non_vector_input():
    with pytest.raises(ValueError, match="one-dimensional"):
        statistics.adjust_pvalues([[0.1, 0.2]])
