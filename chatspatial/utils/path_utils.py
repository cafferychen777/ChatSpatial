"""Path handling utilities for ChatSpatial MCP server.

This module provides robust path handling that works correctly regardless of
the current working directory when the MCP server is launched.

Key features:
- Resolves relative paths intelligently (respects user configuration)
- Prevents outputs from falling into package directory
- Automatic fallback to safe directories for permission issues
- Write permission validation before returning paths

Design principle: Outputs should NEVER pollute the package source directory.
"""

import os
import tempfile
import warnings
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

from ..config import (
    PACKAGE_ROOT,
    get_default_output_dir,
    get_temp_output_dir,
    is_inside_package_dir,
    is_writable_dir,
)
from .exceptions import ParameterError


def validate_path_component(value: str, *, name: str) -> str:
    """Validate an identifier before using it as one filesystem component.

    Paths derived from dataset or analysis identifiers must never escape their
    managed root. Reject separators explicitly so the contract is identical on
    POSIX and Windows, regardless of the host running the validation.
    """
    if (
        not isinstance(value, str)
        or not value
        or value in {".", ".."}
        or "/" in value
        or "\\" in value
        or "\x00" in value
    ):
        raise ParameterError(f"{name} must be a non-empty filesystem-safe identifier")
    return value


@contextmanager
def atomic_output_path(destination: Path) -> Iterator[Path]:
    """Yield a same-directory staging path and atomically publish on success.

    Keeping the staging file beside the destination makes ``os.replace`` an
    atomic rename on the destination filesystem. A failed writer leaves any
    previous destination untouched and removes its incomplete staging file.
    """
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    suffix = f".tmp{destination.suffix}" if destination.suffix else ".tmp"
    descriptor, staging_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=suffix,
    )
    os.close(descriptor)
    staging_path = Path(staging_name)

    try:
        yield staging_path
        os.replace(staging_path, destination)
    finally:
        staging_path.unlink(missing_ok=True)


def get_safe_output_path(
    output_dir: str,
    fallback_to_tmp: bool = True,
    create_if_missing: bool = True,
) -> Path:
    """Get safe, writable output directory path.

    This function handles path resolution robustly with package protection:

    1. For relative paths (e.g., "./visualizations"):
       - If cwd is OUTSIDE package directory: resolve against cwd
       - If cwd is INSIDE package directory: resolve against safe default dir

    2. For absolute paths:
       - If path is INSIDE package directory: reject and use safe default
       - Otherwise: use the absolute path directly

    3. Permission handling:
       - Tests write permission before returning
       - Falls back to the platform temporary output directory if needed

    Args:
        output_dir: User-provided output directory (relative or absolute)
        fallback_to_tmp: If True, use the platform temporary directory when
            output_dir is not writable
        create_if_missing: If True, create directory if it doesn't exist

    Returns:
        Absolute path to writable output directory (guaranteed outside package)

    Raises:
        PermissionError: If no writable path can be found (when fallback disabled)

    Examples:
        >>> # Relative path with cwd outside package
        >>> # cwd = /home/user/data
        >>> path = get_safe_output_path("./outputs")
        >>> # Returns: /home/user/data/outputs

        >>> # Relative path with cwd inside package (PROTECTED)
        >>> # cwd = /path/to/chatspatial/chatspatial
        >>> path = get_safe_output_path("./outputs")
        >>> # Returns: ~/chatspatial_outputs/outputs (not package dir!)

        >>> # Absolute path
        >>> path = get_safe_output_path("/tmp/my_outputs")
        >>> # Returns: /tmp/my_outputs
    """
    user_path = Path(output_dir)
    cwd = Path.cwd()

    # Determine the base directory for path resolution
    if user_path.is_absolute():
        # Absolute path - check if it's inside package
        if is_inside_package_dir(user_path):
            warnings.warn(
                f"Output path {user_path} is inside package directory. "
                f"Using safe default: {get_default_output_dir()}",
                UserWarning,
                stacklevel=2,
            )
            # Try to preserve path structure relative to package root
            try:
                relative_part = user_path.relative_to(PACKAGE_ROOT)
                target_path = get_default_output_dir() / relative_part
            except ValueError:
                # Shouldn't happen since is_inside_package_dir returned True
                target_path = get_default_output_dir() / user_path.name
        else:
            target_path = user_path
    else:
        # Relative path - check if cwd is inside package
        if is_inside_package_dir(cwd):
            # CWD is inside package - use safe default as base
            # Preserve full relative path structure to avoid collisions
            # e.g., "./a/plots" -> ~/chatspatial_outputs/a/plots
            #       "./b/plots" -> ~/chatspatial_outputs/b/plots (no collision)
            base_dir = get_default_output_dir()
            target_path = base_dir / user_path
        else:
            # CWD is safe - resolve against cwd (standard behavior)
            target_path = cwd / user_path

    # Ensure the resolved path is not inside package (double check)
    if is_inside_package_dir(target_path):
        safe_base = get_default_output_dir()
        try:
            relative_part = target_path.relative_to(PACKAGE_ROOT)
            target_path = safe_base / relative_part
        except ValueError:
            target_path = safe_base / target_path.name

    # Try to create/verify the directory
    try:
        if not is_writable_dir(target_path, create=create_if_missing):
            raise PermissionError("directory permission probe failed")

        return target_path

    except (OSError, PermissionError) as e:
        # If fallback enabled, try temp directory
        if fallback_to_tmp:
            fallback_path = get_temp_output_dir()
            warnings.warn(
                f"Cannot write to {target_path}: {e}. Falling back to {fallback_path}",
                UserWarning,
                stacklevel=2,
            )
            return fallback_path
        else:
            raise PermissionError(
                f"Cannot write to output directory: {target_path}. Error: {e}"
            ) from e
