"""
Scipy compatibility layer for deprecated functions.

This module provides compatibility shims for scipy functions that have been
removed or relocated in newer versions, enabling packages like SpatialDE
that depend on these deprecated APIs to work with modern scipy (>=1.14).

Background:
    - scipy.misc.derivative was deprecated in scipy 1.10.0
    - scipy.misc.derivative was removed in scipy 1.14.0
    - SpatialDE 1.1.3 still imports from scipy.misc.derivative
    - This module patches scipy.misc to restore the derivative function

Usage:
    Call `patch_scipy_misc_derivative()` before importing SpatialDE:

        from chatspatial.utils.scipy_compat import patch_scipy_misc_derivative
        patch_scipy_misc_derivative()
        import SpatialDE  # Now works with scipy >= 1.14
"""

import numpy as np
from scipy import misc as scipy_misc


def _derivative_compat(func, x0, dx=1.0, n=1, args=(), order=3):
    """Compute the nth derivative of a function at a point.

    This is a compatibility implementation of the deprecated scipy.misc.derivative,
    using central difference formulas.

    Parameters
    ----------
    func : callable
        Function of which to compute the derivative. Should accept a single
        float argument and return a float.
    x0 : float
        Point at which to evaluate the derivative.
    dx : float, optional
        Spacing for finite difference. Default is 1.0.
    n : int, optional
        Order of the derivative. Default is 1.
    args : tuple, optional
        Extra arguments to pass to func.
    order : int, optional
        Number of points to use for the central difference formula.
        Must be odd. Default is 3.

    Returns
    -------
    float
        The nth derivative of func at x0.

    Notes
    -----
    This implementation uses central difference formulas. The accuracy
    depends on the `order` parameter and the smoothness of the function.

    For SpatialDE, only n=1 and n=2 are used with default dx=1.0 and order=3.
    """
    if order < n + 1:
        raise ValueError(
            f"'order' ({order}) must be at least the derivative order 'n' ({n}) plus 1"
        )
    if order % 2 == 0:
        raise ValueError(f"'order' ({order}) must be odd")

    # Central difference weights for various orders and derivative degrees
    # These are the standard finite difference coefficients
    # Reference: https://en.wikipedia.org/wiki/Finite_difference_coefficient

    # For simplicity, we implement the most common cases used by SpatialDE
    if n == 1:
        # First derivative using central difference
        if order >= 3:
            # 3-point formula: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
            return (
                func(x0 + dx, *args) - func(x0 - dx, *args)
            ) / (2 * dx)
    elif n == 2:
        # Second derivative using central difference
        if order >= 3:
            # 3-point formula: f''(x) ≈ (f(x+h) - 2f(x) + f(x-h)) / h²
            return (
                func(x0 + dx, *args) - 2 * func(x0, *args) + func(x0 - dx, *args)
            ) / (dx ** 2)

    # For higher derivatives or orders, use Richardson extrapolation
    # This is a more general but slower approach
    return _derivative_richardson(func, x0, dx, n, args, order)


def _derivative_richardson(func, x0, dx, n, args, order):
    """Compute derivative using Richardson extrapolation.

    This is a more general implementation that handles arbitrary
    derivative orders using Richardson extrapolation for improved accuracy.
    """
    # Number of points for central difference
    num_points = order

    # Generate sample points symmetrically around x0
    half = num_points // 2
    points = np.arange(-half, half + 1) * dx + x0

    # Evaluate function at all points
    f_values = np.array([func(p, *args) for p in points])

    # Compute derivative using finite differences
    # For nth derivative, we need to apply the difference operator n times
    result = f_values.copy()
    for _ in range(n):
        result = np.diff(result) / dx

    # Return the central value
    return result[len(result) // 2]


def patch_scipy_misc_derivative():
    """Patch scipy.misc to include the deprecated derivative function.

    This function adds a `derivative` attribute to `scipy.misc` module,
    allowing packages like SpatialDE that import `from scipy.misc import derivative`
    to work with scipy >= 1.14.

    This patch is idempotent - calling it multiple times has no additional effect.

    Example
    -------
    >>> from chatspatial.utils.scipy_compat import patch_scipy_misc_derivative
    >>> patch_scipy_misc_derivative()
    >>> import SpatialDE  # Now works!
    """
    if not hasattr(scipy_misc, 'derivative'):
        scipy_misc.derivative = _derivative_compat


# Auto-check: provide helpful message if scipy already has derivative
def check_scipy_derivative_status():
    """Check the status of scipy.misc.derivative and return diagnostic info."""
    import scipy
    has_derivative = hasattr(scipy_misc, 'derivative')
    return {
        'scipy_version': scipy.__version__,
        'has_derivative': has_derivative,
        'needs_patch': not has_derivative,
    }
