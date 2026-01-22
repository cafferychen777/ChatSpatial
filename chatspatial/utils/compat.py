"""
Compatibility layer for handling version mismatches in scientific Python dependencies.

This module provides patches and workarounds for known incompatibilities between
different versions of NumPy, SciPy, and dependent packages (CellRank, SpatialDE, etc.).

Design Principles:
1. Minimal Intervention - Only patch when absolutely necessary
2. Centralized Management - All compatibility fixes in one place
3. Traceable - Document the issue source, affected versions, and fix
4. Reversible - Patches should be easy to remove when upstream fixes are released

Sections:
- NumPy 2.x Compatibility (assert_array_equal parameter names)
- SciPy Compatibility (deprecated scipy.misc.derivative, removed csr_matrix.A)
- Package-specific utilities (CellRank, SpatialDE, SpaGCN)
- Diagnostic functions
"""

import functools
from contextlib import contextmanager
from typing import Any, Callable, Generator

import numpy as np
from scipy import misc as scipy_misc

# =============================================================================
# NumPy 2.x Compatibility
# =============================================================================

# Issue: NumPy 2.x changed np.testing.assert_array_equal parameter names
#        from (x, y, ...) to (actual, desired, ...)
#
# Affected packages:
#   - CellRank 2.0.7 uses: np.testing.assert_array_equal(x=..., y=...)
#   - This fails with NumPy 2.x: TypeError: got an unexpected keyword argument 'x'
#
# Fix status:
#   - Fixed in CellRank main branch (uses positional args)
#   - Not yet released to PyPI as of 2025-01
#
# Reference:
#   - NumPy 2.0 Migration Guide
#   - https://github.com/theislab/cellrank/blob/main/src/cellrank/kernels/_velocity_kernel.py


def _numpy2_compat_assert_array_equal(
    *args: Any,
    x: Any = None,
    y: Any = None,
    actual: Any = None,
    desired: Any = None,
    **kwargs: Any,
) -> None:
    """Wrapper for np.testing.assert_array_equal accepting both old and new parameter names.

    NumPy 2.x changed parameter names from (x, y) to (actual, desired).
    This wrapper accepts both conventions for backward compatibility.
    """
    # Resolve arrays from either naming convention
    if args:
        # Positional arguments take precedence
        arr_actual = args[0]
        arr_desired = args[1] if len(args) > 1 else desired or y
    else:
        # Named arguments: prefer new names, fall back to old
        arr_actual = actual if actual is not None else x
        arr_desired = desired if desired is not None else y

    if arr_actual is None or arr_desired is None:
        raise ValueError(
            "assert_array_equal requires two arrays. "
            "Use positional args or (actual, desired) / (x, y) keyword args."
        )

    # Call original function with positional arguments (works in all NumPy versions)
    # Use getattr for dynamic attribute access (__wrapped__ is set at runtime)
    original_func = getattr(np.testing.assert_array_equal, "__wrapped__", None)
    if original_func is not None:
        original_func(arr_actual, arr_desired, **kwargs)
    else:
        # Fallback if not patched (shouldn't happen in normal usage)
        np.testing.assert_array_equal(arr_actual, arr_desired, **kwargs)


def _is_numpy2() -> bool:
    """Check if NumPy major version is 2 or higher."""
    major_version = int(np.__version__.split(".")[0])
    return major_version >= 2


def _patch_numpy_testing() -> Callable[[], None]:
    """Patch np.testing.assert_array_equal for NumPy 2.x compatibility.

    Returns:
        Unpatch function to restore original behavior.
    """
    original_func = np.testing.assert_array_equal

    # Only patch if not already patched
    if hasattr(original_func, "__wrapped__"):
        return lambda: None  # No-op unpatch

    # Store original function using setattr for dynamic attribute setting
    setattr(_numpy2_compat_assert_array_equal, "__wrapped__", original_func)

    # Apply patch
    np.testing.assert_array_equal = _numpy2_compat_assert_array_equal

    def unpatch() -> None:
        np.testing.assert_array_equal = original_func

    return unpatch


@contextmanager
def numpy2_compat() -> Generator[None, None, None]:
    """Context manager for NumPy 2.x compatibility patches.

    Use this when calling code that uses old NumPy parameter conventions.

    Example:
        with numpy2_compat():
            import cellrank as cr
            vk = cr.kernels.VelocityKernel(adata)
            vk.compute_transition_matrix()
    """
    unpatch = _patch_numpy_testing()
    try:
        yield
    finally:
        unpatch()


# =============================================================================
# SciPy Compatibility
# =============================================================================

# Issue: scipy.misc.derivative was removed in scipy 1.14.0
#
# Affected packages:
#   - SpatialDE 1.1.3 imports from scipy.misc.derivative
#   - This fails with scipy >= 1.14: ImportError
#
# Timeline:
#   - scipy 1.10.0: deprecated scipy.misc.derivative
#   - scipy 1.14.0: removed scipy.misc.derivative
#
# Reference:
#   - https://docs.scipy.org/doc/scipy/release/1.14.0-notes.html


def _derivative_compat(
    func: Callable[..., float],
    x0: float,
    dx: float = 1.0,
    n: int = 1,
    args: tuple = (),
    order: int = 3,
) -> float:
    """Compute the nth derivative of a function at a point.

    This is a compatibility implementation of the deprecated scipy.misc.derivative,
    using central difference formulas.

    For SpatialDE, only n=1 and n=2 are used with default dx=1.0 and order=3.
    """
    if order < n + 1:
        raise ValueError(
            f"'order' ({order}) must be at least the derivative order 'n' ({n}) plus 1"
        )
    if order % 2 == 0:
        raise ValueError(f"'order' ({order}) must be odd")

    # Central difference formulas (most common cases used by SpatialDE)
    if n == 1 and order >= 3:
        # f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
        return (func(x0 + dx, *args) - func(x0 - dx, *args)) / (2 * dx)
    elif n == 2 and order >= 3:
        # f''(x) ≈ (f(x+h) - 2f(x) + f(x-h)) / h²
        return (func(x0 + dx, *args) - 2 * func(x0, *args) + func(x0 - dx, *args)) / (
            dx**2
        )

    # For higher derivatives, use Richardson extrapolation
    return _derivative_richardson(func, x0, dx, n, args, order)


def _derivative_richardson(
    func: Callable[..., float],
    x0: float,
    dx: float,
    n: int,
    args: tuple,
    order: int,
) -> float:
    """Compute derivative using Richardson extrapolation for higher orders."""
    half = order // 2
    points = np.arange(-half, half + 1) * dx + x0
    f_values = np.array([func(p, *args) for p in points])

    result = f_values.copy()
    for _ in range(n):
        result = np.diff(result) / dx

    return result[len(result) // 2]


def patch_scipy_misc_derivative() -> None:
    """Patch scipy.misc to include the deprecated derivative function.

    This enables packages like SpatialDE that import `from scipy.misc import derivative`
    to work with scipy >= 1.14.

    This patch is idempotent - calling it multiple times has no additional effect.

    Example:
        >>> patch_scipy_misc_derivative()
        >>> import SpatialDE  # Now works with scipy >= 1.14
    """
    if not hasattr(scipy_misc, "derivative"):
        scipy_misc.derivative = _derivative_compat


# Issue: scipy.sparse.csr_matrix.A was removed in scipy 1.13.0
#
# Affected packages:
#   - SpaGCN 1.2.7 uses csr_matrix.A in multiple places (SpaGCN.py, util.py)
#   - This fails with scipy >= 1.13: AttributeError: 'csr_matrix' has no attribute 'A'
#
# Timeline:
#   - scipy 1.11.0: deprecated csr_matrix.A (use .toarray() instead)
#   - scipy 1.13.0: removed csr_matrix.A
#
# Reference:
#   - https://github.com/scipy/scipy/issues/21049
#   - https://scipy.github.io/devdocs/reference/sparse.migration_to_sparray.html


def _is_scipy_sparse_patched() -> bool:
    """Check if scipy.sparse matrices already have the .A property."""
    from scipy.sparse import csr_matrix

    return hasattr(csr_matrix, "A") and isinstance(
        getattr(type(csr_matrix([0])), "A", None), property
    )


def patch_scipy_sparse_matrix_A() -> None:
    """Patch scipy.sparse matrix classes to restore the deprecated .A property.

    The .A property was a shorthand for .toarray() that was removed in scipy 1.13.
    This patch restores it for packages like SpaGCN that depend on this behavior.

    This patch is idempotent - calling it multiple times has no additional effect.

    Example:
        >>> patch_scipy_sparse_matrix_A()
        >>> import SpaGCN  # Now works with scipy >= 1.13
    """
    if _is_scipy_sparse_patched():
        return

    from scipy.sparse import (
        bsr_matrix,
        coo_matrix,
        csc_matrix,
        csr_matrix,
        dia_matrix,
        dok_matrix,
        lil_matrix,
    )

    # Define the .A property as an alias to .toarray()
    def _toarray_property(self):  # type: ignore[no-untyped-def]
        return self.toarray()

    # Patch all sparse matrix classes that previously had .A
    for sparse_class in [
        csr_matrix,
        csc_matrix,
        coo_matrix,
        bsr_matrix,
        dia_matrix,
        dok_matrix,
        lil_matrix,
    ]:
        if not hasattr(sparse_class, "A"):
            sparse_class.A = property(_toarray_property)


# =============================================================================
# Package-specific Compatibility Utilities
# =============================================================================


def ensure_cellrank_compat() -> Callable[[], None]:
    """Apply all necessary patches for CellRank compatibility.

    Call this before importing or using CellRank with NumPy 2.x.

    Returns:
        Cleanup function to restore original behavior.

    Example:
        cleanup = ensure_cellrank_compat()
        try:
            import cellrank as cr
            vk = cr.kernels.VelocityKernel(adata)
            vk.compute_transition_matrix()
        finally:
            cleanup()
    """
    cleanups: list[Callable[[], None]] = []

    if _is_numpy2():
        cleanups.append(_patch_numpy_testing())

    def cleanup_all() -> None:
        for cleanup in cleanups:
            cleanup()

    return cleanup_all


def cellrank_compat(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to apply CellRank compatibility patches.

    Example:
        @cellrank_compat
        def run_cellrank_analysis(adata):
            import cellrank as cr
            vk = cr.kernels.VelocityKernel(adata)
            vk.compute_transition_matrix()
            return vk
    """

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        cleanup = ensure_cellrank_compat()
        try:
            return func(*args, **kwargs)
        finally:
            cleanup()

    return wrapper


def ensure_spatialde_compat() -> None:
    """Apply all necessary patches for SpatialDE compatibility.

    Call this before importing SpatialDE with scipy >= 1.14.

    Example:
        ensure_spatialde_compat()
        import SpatialDE  # Now works
    """
    patch_scipy_misc_derivative()


def ensure_spagcn_compat() -> None:
    """Apply all necessary patches for SpaGCN compatibility.

    Call this before importing SpaGCN with scipy >= 1.13.

    SpaGCN uses csr_matrix.A which was removed in scipy 1.13.
    This patch restores the .A property as an alias to .toarray().

    Example:
        ensure_spagcn_compat()
        import SpaGCN  # Now works with scipy >= 1.13
    """
    patch_scipy_sparse_matrix_A()


# =============================================================================
# Diagnostic Functions
# =============================================================================

KNOWN_ISSUES = {
    "cellrank_numpy2": {
        "description": "CellRank 2.0.7 uses np.testing.assert_array_equal(x=, y=) "
        "which fails with NumPy 2.x (parameter names changed to actual, desired)",
        "affected_versions": {"cellrank": "<=2.0.7", "numpy": ">=2.0.0"},
        "status": "Fixed in CellRank main branch, awaiting PyPI release",
        "workaround": "Use ensure_cellrank_compat() or @cellrank_compat decorator",
    },
    "spatialde_scipy": {
        "description": "SpatialDE 1.1.3 imports scipy.misc.derivative which was "
        "removed in scipy 1.14.0",
        "affected_versions": {"spatialde": "<=1.1.3", "scipy": ">=1.14.0"},
        "status": "Upstream not fixed",
        "workaround": "Use ensure_spatialde_compat() before importing SpatialDE",
    },
    "spagcn_scipy": {
        "description": "SpaGCN 1.2.7 uses csr_matrix.A which was removed in scipy 1.13.0",
        "affected_versions": {"spagcn": "<=1.2.7", "scipy": ">=1.13.0"},
        "status": "Upstream not fixed",
        "workaround": "Use ensure_spagcn_compat() before importing SpaGCN",
    },
}


def get_compatibility_info() -> dict[str, Any]:
    """Get information about current compatibility status.

    Returns:
        Dictionary with version info, compatibility status, and applied patches.
    """
    import scipy

    return {
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "is_numpy2": _is_numpy2(),
        "scipy_has_derivative": hasattr(scipy_misc, "derivative"),
        "known_issues": KNOWN_ISSUES,
        "patches_available": [
            "numpy2_compat",
            "cellrank_compat",
            "ensure_cellrank_compat",
            "ensure_spatialde_compat",
            "ensure_spagcn_compat",
            "patch_scipy_misc_derivative",
            "patch_scipy_sparse_matrix_A",
        ],
    }


def check_scipy_derivative_status() -> dict[str, Any]:
    """Check the status of scipy.misc.derivative and return diagnostic info."""
    import scipy

    has_derivative = hasattr(scipy_misc, "derivative")
    return {
        "scipy_version": scipy.__version__,
        "has_derivative": has_derivative,
        "needs_patch": not has_derivative,
    }
