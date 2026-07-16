"""
Entry point for ChatSpatial.

This module provides the command-line interface for starting the
ChatSpatial server using either stdio or SSE transport.
"""

import sys
import traceback
from typing import Literal, cast

import click

from . import __version__, config
from .server import mcp

# The server module owns runtime initialization. Keep the configuration handle
# here only for the explicit verbose status requested by the CLI.


@click.group()
@click.version_option(__version__, prog_name="ChatSpatial")
def cli():
    """ChatSpatial - AI-powered spatial transcriptomics analysis"""


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
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Print initialization info",
)
def server(port: int, transport: str, host: str, log_level: str, verbose: bool):
    """Start the ChatSpatial server.

    This command starts the ChatSpatial server using either stdio or SSE transport.
    For stdio transport, the server communicates through standard input/output.
    For SSE transport, the server starts an HTTP server on the specified host and port.
    """
    try:
        if verbose:
            # Re-initialize with verbose output
            config.init_runtime(verbose=True)

        print(
            f"Starting ChatSpatial server with {transport} transport...",
            file=sys.stderr,
        )

        # Set server settings
        mcp.settings.host = host
        mcp.settings.port = port
        mcp.settings.log_level = cast(
            Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], log_level
        )

        # Run the server with the specified transport
        mcp.run(transport=cast(Literal["stdio", "sse", "streamable-http"], transport))

    except Exception as e:
        print(f"Error starting MCP server: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


def main():
    """Main entry point for ChatSpatial CLI"""
    cli()


if __name__ == "__main__":  # pragma: no cover
    main()
