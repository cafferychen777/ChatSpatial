"""Path handling utilities for ChatSpatial MCP server.

This module provides robust path handling that works correctly regardless of
the current working directory when the MCP server is launched.

Key features:
- Resolves relative paths against project root (not cwd)
- Automatic fallback to /tmp for permission issues
- Security checks to prevent writing outside safe directories
"""

import os
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


def is_safe_output_path(path: Path) -> bool:
    """Check if output path is safe (within project or /tmp).

    This provides security by ensuring files are only written to:
    1. Project directory or its subdirectories
    2. /tmp/chatspatial directory

    Args:
        path: Path to check

    Returns:
        True if path is safe for output, False otherwise

    Examples:
        >>> # Safe paths
        >>> is_safe_output_path(Path("/Users/.../chatspatial/outputs"))
        True
        >>> is_safe_output_path(Path("/tmp/chatspatial/outputs"))
        True

        >>> # Unsafe paths
        >>> is_safe_output_path(Path("/etc/outputs"))
        False
        >>> is_safe_output_path(Path("/Users/other_user/outputs"))
        False
    """
    # Normalize path to absolute
    abs_path = path.resolve() if not path.is_absolute() else path

    project_root = get_project_root()
    tmp_root = Path("/tmp/chatspatial")

    # Allow paths within project directory
    try:
        abs_path.relative_to(project_root)
        return True
    except ValueError:
        pass

    # Allow paths within /tmp/chatspatial
    try:
        abs_path.relative_to(tmp_root)
        return True
    except ValueError:
        pass

    # Disallow all other paths
    return False


def get_output_dir_from_config(default: str = "./outputs") -> str:
    """Get output directory from environment variable or configuration.

    Priority order:
    1. CHATSPATIAL_OUTPUT_DIR environment variable (highest priority)
    2. Default value (usually "./outputs")

    This allows users to configure the output directory via:
    - Claude Desktop config: env.CHATSPATIAL_OUTPUT_DIR
    - Shell environment: export CHATSPATIAL_OUTPUT_DIR=/path/to/outputs

    Args:
        default: Default output directory if no configuration found

    Returns:
        Output directory path (relative or absolute)

    Examples:
        >>> # With environment variable set
        >>> os.environ["CHATSPATIAL_OUTPUT_DIR"] = "/tmp/results"
        >>> get_output_dir_from_config()
        '/tmp/results'

        >>> # Without environment variable
        >>> get_output_dir_from_config()
        './outputs'

    Configuration example for Claude Desktop (claude_desktop_config.json):
    ```json
    {
      "mcpServers": {
        "chatspatial": {
          "command": "python",
          "args": ["-m", "chatspatial"],
          "cwd": "/Users/username/my_project",
          "env": {
            "CHATSPATIAL_OUTPUT_DIR": "./results"
          }
        }
      }
    }
    ```
    """
    # Check environment variable (highest priority)
    if env_dir := os.getenv("CHATSPATIAL_OUTPUT_DIR"):
        return env_dir

    # Return default
    return default


def resolve_output_path(
    output_dir: str,
    default_dir: str = "./outputs",
) -> Path:
    """Resolve output directory path safely.

    Resolves relative paths against current working directory,
    following standard Unix/Python conventions.

    Args:
        output_dir: User-provided output directory
        default_dir: Default directory if output_dir is None

    Returns:
        Resolved absolute path

    Examples:
        >>> resolve_output_path("./outputs")
        PosixPath('<cwd>/outputs')

        >>> resolve_output_path("/tmp/test")
        PosixPath('/tmp/test')
    """
    path = Path(output_dir or default_dir)

    if path.is_absolute():
        return path
    else:
        # Resolve against CWD (respects user configuration)
        return Path.cwd() / path
