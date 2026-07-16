"""
Device utilities for compute backend selection.

This module provides lazy-loaded device detection and selection functions
for GPU/CPU computation. Follows the same design principles as compute.py.

Design Principles:
1. Lazy Loading: torch is only imported when needed
2. Pure Functions: No side effects, callers decide how to handle results
3. String Returns: Callers convert to torch.device if needed
4. Composable: Basic building blocks for various use cases
5. Single Source of Truth: All GPU/MPS configuration in one place

Usage:
    # Simple device selection (direct PyTorch)
    device = get_device(prefer_gpu=True)

    # With warning when GPU unavailable
    device = get_device(prefer_gpu=params.use_gpu)
    if params.use_gpu and device == "cpu":
        await ctx.warning("GPU requested but not available")

    # Convert to torch.device when needed
    import torch
    device = torch.device(get_device(prefer_gpu=True))

    # For scvi-tools / PyTorch Lightning methods
    accelerator = get_accelerator(prefer_gpu=True)  # "gpu" or "cpu"
"""

import os
from typing import TYPE_CHECKING, Any

from .dependency_manager import get, require

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


# =============================================================================
# MPS Configuration (Apple Silicon)
# =============================================================================

# Track if MPS has been configured (idempotent)
_mps_configured = False


def _configure_mps() -> None:
    """Configure MPS backend for optimal performance.

    This function configures PyTorch MPS settings for Apple Silicon Macs.
    It is called automatically when MPS is selected via get_device().

    Key configuration:
    - PYTORCH_MPS_HIGH_WATERMARK_RATIO=0.0: Disables memory limit enforcement.
      Mac's Unified Memory allows GPU/CPU to share system RAM dynamically.
      The default limit (often ~36GB) is too conservative for large models
      like STAGATE/GraphST that need flexible memory allocation.

    This is idempotent - safe to call multiple times.
    """
    global _mps_configured
    if _mps_configured:
        return

    # Disable MPS memory limit - Mac Unified Memory handles this dynamically
    # Without this, large models fail with "MPS backend out of memory"
    # even when system has available RAM
    os.environ.setdefault("PYTORCH_MPS_HIGH_WATERMARK_RATIO", "0.0")

    _mps_configured = True


# =============================================================================
# Availability Checks (has_* pattern)
# =============================================================================


def cuda_available() -> bool:
    """Check if CUDA GPU is available.

    Lazy imports torch to avoid loading it when not needed.

    Returns:
        True if CUDA is available, False otherwise
    """
    torch = get("torch")
    return torch is not None and torch.cuda.is_available()


def mps_available() -> bool:
    """Check if Apple Silicon MPS is available.

    Lazy imports torch to avoid loading it when not needed.

    Returns:
        True if MPS is available, False otherwise
    """
    torch = get("torch")
    return (
        torch is not None
        and hasattr(torch.backends, "mps")
        and torch.backends.mps.is_available()
    )


# =============================================================================
# Device Selection (core function)
# =============================================================================


def get_device(
    prefer_gpu: bool = False,
    allow_mps: bool = True,
) -> str:
    """Select compute device based on preference and availability.

    This is THE single source of truth for direct PyTorch device selection.
    Returns a device string that can be used directly or converted to
    torch.device.

    For scvi-tools / PyTorch Lightning methods, use get_accelerator() instead.

    Args:
        prefer_gpu: If True, try to use GPU (CUDA first, then MPS)
        allow_mps: If True, allow MPS as fallback when CUDA unavailable.
            Most PyTorch methods support MPS; set False to opt out.

    Returns:
        Device string: "cuda:0", "mps", or "cpu"

    Examples:
        # Basic usage
        device = get_device(prefer_gpu=True)  # "cuda:0", "mps", or "cpu"

        # Convert to torch.device
        import torch
        device = torch.device(get_device(prefer_gpu=True))

        # With warning when requested but unavailable
        device = get_device(prefer_gpu=params.use_gpu)
        if params.use_gpu and device == "cpu":
            await ctx.warning("GPU requested but not available - using CPU")
    """
    if prefer_gpu:
        if cuda_available():
            return "cuda:0"
        if allow_mps and mps_available():
            # Configure MPS before first use (idempotent)
            _configure_mps()
            return "mps"
    return "cpu"


def get_accelerator(prefer_gpu: bool = False) -> str:
    """Select PyTorch Lightning accelerator for scvi-tools methods.

    Lightning's accelerator="gpu" auto-detects CUDA vs MPS.
    Use this for Cell2location, DestVI, Stereoscope, scVI, VeloVI.

    For direct PyTorch methods, use get_device() instead.

    Args:
        prefer_gpu: If True, return "gpu" when any GPU is available

    Returns:
        Accelerator string: "gpu" or "cpu"
    """
    if prefer_gpu:
        if cuda_available():
            return "gpu"
        if mps_available():
            _configure_mps()
            return "gpu"
    return "cpu"


# =============================================================================
# Async Helper with Context Warning
# =============================================================================


async def resolve_device_async(
    prefer_gpu: bool,
    ctx: "ToolContext",
    allow_mps: bool = True,
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
    ot = require("ot", feature="optimal transport backend selection")

    if use_gpu and cuda_available():
        return ot.backend.TorchBackend()
    return ot.backend.NumpyBackend()
