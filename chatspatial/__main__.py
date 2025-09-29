"""
Entry point for ChatSpatial.

This module provides the command-line interface for starting the
ChatSpatial server using either stdio or SSE transport.
"""

import os
import sys
import traceback
import warnings

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

from .server import mcp


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
