"""
Device utilities for compute backend selection.

This module provides lazy-loaded device detection and selection functions
for GPU/CPU computation. Follows the same design principles as compute.py.

Design Principles:
1. Lazy Loading: torch is only imported when needed
2. Pure Functions: No side effects, callers decide how to handle results
3. String Returns: Callers convert to torch.device if needed
4. Composable: Basic building blocks for various use cases

Usage:
    # Simple device selection
    device = get_device(use_gpu=True)

    # With warning when GPU unavailable
    device = get_device(use_gpu=params.use_gpu)
    if params.use_gpu and device == "cpu":
        await ctx.warning("GPU requested but not available")

    # Convert to torch.device when needed
    import torch
    device = torch.device(get_device(use_gpu=True, allow_mps=True))
"""

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


# =============================================================================
# Availability Checks (has_* pattern)
# =============================================================================


def cuda_available() -> bool:
    """Check if CUDA GPU is available.

    Lazy imports torch to avoid loading it when not needed.

    Returns:
        True if CUDA is available, False otherwise
    """
    try:
        import torch

        return torch.cuda.is_available()
    except ImportError:
        return False


def mps_available() -> bool:
    """Check if Apple Silicon MPS is available.

    Lazy imports torch to avoid loading it when not needed.

    Returns:
        True if MPS is available, False otherwise
    """
    try:
        import torch

        return hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
    except ImportError:
        return False


# =============================================================================
# Device Selection (core function)
# =============================================================================


def get_device(
    prefer_gpu: bool = False,
    allow_mps: bool = False,
) -> str:
    """Select compute device based on preference and availability.

    This is THE single source of truth for device selection across ChatSpatial.
    Returns a device string that can be used directly or converted to torch.device.

    Args:
        prefer_gpu: If True, try to use GPU (CUDA first, then MPS if allowed)
        allow_mps: If True, allow MPS as fallback when CUDA unavailable

    Returns:
        Device string: "cuda:0", "mps", or "cpu"

    Examples:
        # Basic usage
        device = get_device(use_gpu=True)  # "cuda:0" or "cpu"

        # With MPS support (Apple Silicon)
        device = get_device(use_gpu=True, allow_mps=True)  # "cuda:0", "mps", or "cpu"

        # Convert to torch.device
        import torch
        device = torch.device(get_device(prefer_gpu=True))

        # With warning when requested but unavailable
        device = get_device(params.use_gpu)
        if params.use_gpu and device == "cpu":
            await ctx.warning("GPU requested but not available - using CPU")
    """
    if prefer_gpu:
        if cuda_available():
            return "cuda:0"
        if allow_mps and mps_available():
            return "mps"
    return "cpu"


# =============================================================================
# Async Helper with Context Warning
# =============================================================================


async def resolve_device_async(
    prefer_gpu: bool,
    ctx: "ToolContext",
    allow_mps: bool = False,
    warn_on_fallback: bool = True,
) -> str:
    """Select device with optional warning when GPU unavailable.

    Convenience function for async tools that want automatic warning.

    Args:
        prefer_gpu: If True, try to use GPU
        ctx: ToolContext for logging warnings
        allow_mps: If True, allow MPS as fallback
        warn_on_fallback: If True, warn when requested GPU is unavailable

    Returns:
        Device string: "cuda:0", "mps", or "cpu"
    """
    device = get_device(prefer_gpu=prefer_gpu, allow_mps=allow_mps)

    if warn_on_fallback and prefer_gpu and device == "cpu":
        await ctx.warning("GPU requested but not available - using CPU")

    return device


# =============================================================================
# Specialized Backend Functions
# =============================================================================


def get_ot_backend(use_gpu: bool = False) -> Any:
    """Get optimal transport backend for PASTE alignment.

    Args:
        use_gpu: If True, try to use TorchBackend with CUDA

    Returns:
        POT backend (TorchBackend if CUDA available and requested, else NumpyBackend)
    """
    import ot

    if use_gpu and cuda_available():
        return ot.backend.TorchBackend()
    return ot.backend.NumpyBackend()
