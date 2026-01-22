"""Path handling utilities for ChatSpatial MCP server.

This module provides robust path handling that works correctly regardless of
the current working directory when the MCP server is launched.

Key features:
- Resolves relative paths against CWD (respects user configuration)
- Automatic fallback to /tmp for permission issues
- Write permission validation before returning paths
"""

import warnings
from pathlib import Path

# Project root directory (based on utils module location)
# This is always correct regardless of cwd
_PROJECT_ROOT = Path(__file__).parent.parent.resolve()


def get_project_root() -> Path:
    """Get ChatSpatial project root directory.

    Returns:
        Absolute path to project root, regardless of current working directory.

    Example:
        >>> root = get_project_root()
        >>> print(root)
        /path/to/chatspatial/chatspatial
    """
    return _PROJECT_ROOT


def get_safe_output_path(
    output_dir: str,
    fallback_to_tmp: bool = True,
    create_if_missing: bool = True,
) -> Path:
    """Get safe, writable output directory path.

    This function handles path resolution robustly:
    - Relative paths are resolved against CURRENT WORKING DIRECTORY (respects user config)
    - Tests write permission before returning
    - Falls back to /tmp/chatspatial/outputs if original path not writable

    Args:
        output_dir: User-provided output directory (relative or absolute)
        fallback_to_tmp: If True, fallback to /tmp if output_dir not writable
        create_if_missing: If True, create directory if it doesn't exist

    Returns:
        Absolute path to writable output directory

    Raises:
        PermissionError: If no writable path can be found (when fallback disabled)

    Examples:
        >>> # Relative path (resolved against cwd)
        >>> path = get_safe_output_path("./outputs")
        >>> # Returns: <cwd>/outputs

        >>> # Absolute path
        >>> path = get_safe_output_path("/tmp/my_outputs")
        >>> # Returns: /tmp/my_outputs

        >>> # Read-only path with fallback
        >>> path = get_safe_output_path("/outputs")
        >>> # Returns: /tmp/chatspatial/outputs (with warning)
    """
    # Convert to Path object
    user_path = Path(output_dir)

    # If absolute path, use directly; otherwise resolve against CWD
    if user_path.is_absolute():
        target_path = user_path
    else:
        # For relative paths, resolve against CWD (respects user configuration!)
        # This follows standard Unix/Python conventions
        target_path = Path.cwd() / user_path

    # Try to create/verify the directory
    try:
        if create_if_missing:
            target_path.mkdir(parents=True, exist_ok=True)

        # Test write permission by creating a temporary test file
        test_file = target_path / ".write_test"
        test_file.touch()
        test_file.unlink()

        return target_path

    except (OSError, PermissionError) as e:
        # If fallback enabled, try temp directory
        if fallback_to_tmp:
            warnings.warn(
                f"Cannot write to {target_path}: {e}. "
                f"Falling back to /tmp/chatspatial/outputs",
                UserWarning,
                stacklevel=2,
            )

            fallback_path = Path("/tmp/chatspatial/outputs")
            fallback_path.mkdir(parents=True, exist_ok=True)
            return fallback_path
        else:
            raise PermissionError(
                f"Cannot write to output directory: {target_path}. " f"Error: {e}"
            ) from e
