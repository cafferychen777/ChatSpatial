"""
Entry point for ChatSpatial.

This module provides the command-line interface for starting the
ChatSpatial server using either stdio or SSE transport.
"""

import os
import sys
import traceback
import warnings
from pathlib import Path

import click

# Suppress warnings to speed up startup
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# CRITICAL: Disable all progress bars to prevent stdout pollution in MCP protocol
# MCP uses JSON-RPC over stdio, any non-JSON output breaks communication
os.environ["TQDM_DISABLE"] = "1"  # Disable tqdm globally

# Configure scientific libraries to suppress output
try:
    import scanpy as sc

    sc.settings.verbosity = 0  # Suppress scanpy output
except ImportError:
    pass  # scanpy may not be installed yet

# IMPORTANT: Intelligent working directory handling
# Only change cwd when it's clearly problematic, otherwise respect user configuration
PROJECT_ROOT = Path(__file__).parent.resolve()
user_cwd = Path.cwd()

# Identify problematic working directories that should be changed
problematic_cwds = [
    Path("/"),  # Root directory
    Path("/tmp"),  # Temp directory
    Path("/var"),  # System directory
    Path("/usr"),  # System directory
    Path("/etc"),  # System directory
]

# Check if current cwd is problematic
is_problematic = (
    # Check exact match
    user_cwd in problematic_cwds
    or
    # Check if cwd is a parent of problematic directories
    any(user_cwd == p.parent for p in problematic_cwds)
    or
    # Check if directory doesn't exist
    not user_cwd.exists()
    or
    # Check if it's a temporary npx directory (common MCP issue)
    "_npx" in str(user_cwd)
    or ".npm" in str(user_cwd)
)

if is_problematic:
    print(
        f"WARNING:Working directory appears problematic: {user_cwd}\n"
        f"   Changing to project root: {PROJECT_ROOT}\n"
        f"   (This ensures file operations work correctly)",
        file=sys.stderr,
    )
    os.chdir(PROJECT_ROOT)
else:
    print(
        f"Using configured working directory: {user_cwd}\n"
        f"  (Project root: {PROJECT_ROOT})",
        file=sys.stderr,
    )
    # Keep user's configured cwd - don't change it!

from .server import mcp  # noqa: E402


@click.group()
def cli():
    """ChatSpatial - AI-powered spatial transcriptomics analysis"""
    pass


@cli.command()
@click.option("--port", default=8000, help="Port to listen on for SSE transport")
@click.option(
    "--transport",
    type=click.Choice(["stdio", "sse"]),
    default="stdio",
    help="Transport type (stdio or sse)",
)
@click.option(
    "--host",
    default="127.0.0.1",  # nosec B104 - Default to localhost for security
    help="Host to bind to for SSE transport",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
    default="INFO",
    help="Logging level",
)
def server(port: int, transport: str, host: str, log_level: str):
    """Start the ChatSpatial server.

    This command starts the ChatSpatial server using either stdio or SSE transport.
    For stdio transport, the server communicates through standard input/output.
    For SSE transport, the server starts an HTTP server on the specified host and port.
    """
    try:
        # Configure server settings
        print(
            f"Starting ChatSpatial server with {transport} transport...",
            file=sys.stderr,
        )

        # Set server settings
        mcp.settings.host = host
        mcp.settings.port = port
        mcp.settings.log_level = log_level

        # Run the server with the specified transport
        # This is the recommended way to run a FastMCP server
        mcp.run(transport=transport)

    except Exception as e:
        print(f"Error starting MCP server: {str(e)}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


def main():
    """Main entry point for ChatSpatial CLI"""
    cli()


if __name__ == "__main__":
    main()
