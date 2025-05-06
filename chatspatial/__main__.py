"""
Entry point for ChatSpatial.

This module provides the command-line interface for starting the
ChatSpatial server using either stdio or SSE transport.
"""

import click
import sys
import traceback
from .server import mcp

@click.command()
@click.option("--port", default=8000, help="Port to listen on for SSE transport")
@click.option(
    "--transport",
    type=click.Choice(["stdio", "sse"]),
    default="stdio",
    help="Transport type (stdio or sse)"
)
@click.option(
    "--host",
    default="0.0.0.0",
    help="Host to bind to for SSE transport"
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
    default="INFO",
    help="Logging level"
)
def main(port: int, transport: str, host: str, log_level: str):
    """Start the ChatSpatial server.

    This command starts the ChatSpatial server using either stdio or SSE transport.
    For stdio transport, the server communicates through standard input/output.
    For SSE transport, the server starts an HTTP server on the specified host and port.
    """
    try:
        # Configure server settings
        print(f"Starting ChatSpatial server with {transport} transport...", file=sys.stderr)

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

if __name__ == "__main__":
    main()
