"""
Configuration module for ChatSpatial - Single Source of Truth (SSOT).

This module centralizes all runtime configuration:
- Environment variables (TQDM, Dask, etc.)
- Warning filters
- Scanpy settings
- Path constants and utilities

Design principles:
1. Single Source of Truth: All config defined here, imported elsewhere
2. Idempotent: Safe to call init_runtime() multiple times
3. No side effects on cwd: Never change working directory
4. Package directory protection: Outputs never fall into package dir
"""

import os
import sys
import warnings
from pathlib import Path
from typing import Optional

# =============================================================================
# Path Constants (computed once at module load)
# =============================================================================

# Package root directory (where this config.py lives)
# This is always correct regardless of cwd
PACKAGE_ROOT: Path = Path(__file__).parent.resolve()

# User's home directory
HOME_DIR: Path = Path.home()

# Default output directory (outside package, in user's home)
DEFAULT_OUTPUT_DIR: Path = HOME_DIR / "chatspatial_outputs"


# =============================================================================
# Path Utilities
# =============================================================================


def is_inside_package_dir(path: Optional[Path] = None) -> bool:
    """Check if a path is inside the package directory.

    This is used to prevent outputs from polluting the source code directory.

    Args:
        path: Path to check. If None, checks current working directory.

    Returns:
        True if the path is inside the package directory.

    Examples:
        >>> is_inside_package_dir(Path("/Users/me/chatspatial/chatspatial"))
        True
        >>> is_inside_package_dir(Path("/Users/me/data"))
        False
    """
    if path is None:
        path = Path.cwd()

    try:
        # Resolve to absolute path for accurate comparison
        resolved_path = path.resolve()
        # Check if path is inside package directory
        resolved_path.relative_to(PACKAGE_ROOT)
        return True
    except ValueError:
        # relative_to raises ValueError if path is not relative to PACKAGE_ROOT
        return False


def get_default_output_dir() -> Path:
    """Get the default output directory for generated files.

    Priority:
    1. CHATSPATIAL_OUTPUT_DIR environment variable (if set)
    2. Current working directory (if not inside package and writable)
    3. ~/chatspatial_outputs (user's home directory)
    4. /tmp/chatspatial/outputs (fallback)

    Returns:
        Absolute path to a writable output directory.
    """
    # Priority 1: Environment variable
    env_dir = os.environ.get("CHATSPATIAL_OUTPUT_DIR")
    if env_dir:
        env_path = Path(env_dir)
        if _is_writable_dir(env_path):
            return env_path.resolve()

    # Priority 2: Current working directory (if safe)
    cwd = Path.cwd()
    if not is_inside_package_dir(cwd) and _is_writable_dir(cwd):
        return cwd

    # Priority 3: User's home directory
    if _is_writable_dir(DEFAULT_OUTPUT_DIR, create=True):
        return DEFAULT_OUTPUT_DIR

    # Priority 4: Temp directory (last resort)
    tmp_dir = Path("/tmp/chatspatial/outputs")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    return tmp_dir


def _is_writable_dir(path: Path, create: bool = False) -> bool:
    """Check if a directory is writable.

    Args:
        path: Directory path to check.
        create: If True, attempt to create the directory if missing.

    Returns:
        True if the directory exists and is writable.
    """
    try:
        if create:
            path.mkdir(parents=True, exist_ok=True)

        if not path.exists():
            return False

        # Test write permission
        test_file = path / ".write_test"
        test_file.touch()
        test_file.unlink()
        return True
    except (OSError, PermissionError):
        return False


# =============================================================================
# Warning Configuration
# =============================================================================


def _configure_warnings() -> None:
    """Configure warning filters for known issues.

    This function is idempotent - safe to call multiple times.
    """
    # Broad filters for startup speed (matching original behavior)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    # Specific filters for known library issues
    warnings.filterwarnings(
        "ignore",
        message="The legacy Dask DataFrame implementation is deprecated",
        category=FutureWarning,
    )

    warnings.filterwarnings(
        "ignore",
        message="functools.partial will be a method descriptor",
        category=FutureWarning,
    )

    warnings.filterwarnings(
        "ignore",
        message="nopython is set for njit and is ignored",
        category=RuntimeWarning,
    )

    warnings.filterwarnings(
        "ignore",
        message="Importing read_text from `anndata` is deprecated",
        category=FutureWarning,
    )


# =============================================================================
# Environment Variable Configuration
# =============================================================================


def _configure_environment() -> None:
    """Configure environment variables for dependencies.

    This function is idempotent - safe to call multiple times.
    """
    # CRITICAL: Disable progress bars to prevent stdout pollution
    # MCP uses JSON-RPC over stdio, any non-JSON output breaks communication
    os.environ["TQDM_DISABLE"] = "1"

    # Configure Dask to use new DataFrame implementation
    os.environ.setdefault("DASK_DATAFRAME__QUERY_PLANNING", "True")


# =============================================================================
# Library Configuration
# =============================================================================


def _configure_libraries() -> None:
    """Configure scientific libraries (scanpy, dask, etc.).

    This function is idempotent - safe to call multiple times.

    Dependency classification (aligned with pyproject.toml):
    - REQUIRED: scanpy, anndata, squidpy - direct import, fail fast if missing
    - OPTIONAL: dask - try/except, graceful degradation
    """
    # Configure Dask (OPTIONAL - not in pyproject.toml dependencies)
    # Used by some spatial data formats but not core functionality
    try:
        import dask

        dask.config.set({"dataframe.query-planning": True})
    except ImportError:
        pass  # Dask not installed, graceful degradation

    # Configure Scanpy (REQUIRED - in pyproject.toml dependencies)
    # This is a core dependency - if missing, pip install failed
    import scanpy as sc

    sc.settings.verbosity = 0  # Suppress output
    sc.settings.n_jobs = -1  # Use all CPU cores for parallel computation


# =============================================================================
# Runtime Initialization (Main Entry Point)
# =============================================================================

# Track initialization state to ensure idempotence
_initialized: bool = False


def init_runtime(verbose: bool = False) -> None:
    """Initialize ChatSpatial runtime environment.

    This is the main entry point for configuration. It should be called
    early in the application lifecycle (e.g., in __main__.py or server.py).

    This function is idempotent - calling it multiple times is safe and
    has no additional effect after the first call.

    Args:
        verbose: If True, print initialization info to stderr.

    Example:
        >>> from chatspatial.config import init_runtime
        >>> init_runtime()  # Call once at startup
    """
    global _initialized

    if _initialized:
        return

    # Configure in order of dependency
    _configure_environment()  # Env vars first (may affect library imports)
    _configure_warnings()  # Warnings before library imports
    _configure_libraries()  # Libraries last

    _initialized = True

    if verbose:
        print(
            f"ChatSpatial initialized:\n"
            f"  Package root: {PACKAGE_ROOT}\n"
            f"  Working directory: {Path.cwd()}\n"
            f"  Default output: {get_default_output_dir()}",
            file=sys.stderr,
        )


# =============================================================================
# Auto-initialization on module import
# =============================================================================

# Initialize when this module is imported
# This ensures configuration is applied as early as possible
init_runtime()
