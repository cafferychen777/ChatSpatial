"""
Entry point for ChatSpatial.

This module provides the command-line interface for starting the
ChatSpatial server using stdio or HTTP transports.
"""

import logging
import sys
import traceback
from typing import Literal, cast

import click
from mcp.server.transport_security import TransportSecuritySettings

from . import __version__, config
from .server import mcp

_LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1"}
_LOOPBACK_ALLOWED_HOSTS = ("127.0.0.1:*", "localhost:*", "[::1]:*")
_LOOPBACK_ALLOWED_ORIGINS = (
    "http://127.0.0.1:*",
    "http://localhost:*",
    "http://[::1]:*",
)

# The server module owns runtime initialization. Keep the configuration handle
# here only for the explicit verbose status requested by the CLI.


@click.group()
@click.version_option(__version__, prog_name="ChatSpatial")
def cli():
    """ChatSpatial - AI-powered spatial transcriptomics analysis"""


@cli.command()
@click.option("--port", default=8000, help="Port to listen on for HTTP transport")
@click.option(
    "--transport",
    type=click.Choice(["stdio", "streamable-http", "sse"]),
    default="stdio",
    help="Transport type; streamable-http is preferred for HTTP deployments",
)
@click.option(
    "--host",
    default="127.0.0.1",  # nosec B104 - Default to localhost for security
    help="Host to bind to for HTTP transport",
)
@click.option(
    "--allowed-host",
    "allowed_hosts",
    multiple=True,
    help="Allowed HTTP Host header, repeatable; required for non-loopback binding",
)
@click.option(
    "--allowed-origin",
    "allowed_origins",
    multiple=True,
    help="Allowed HTTP Origin header, repeatable",
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
def server(
    port: int,
    transport: str,
    host: str,
    allowed_hosts: tuple[str, ...],
    allowed_origins: tuple[str, ...],
    log_level: str,
    verbose: bool,
):
    """Start the ChatSpatial server.

    STDIO is the secure default for local clients. Streamable HTTP implements
    the MCP 2026-07-28 stateless protocol and also serves legacy MCP clients.
    SSE remains available only for compatibility with older deployments.
    """
    try:
        if verbose:
            # Re-initialize with verbose output
            config.init_runtime(verbose=True)

        print(
            f"Starting ChatSpatial server with {transport} transport...",
            file=sys.stderr,
        )

        resolved_log_level = cast(
            Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            log_level,
        )
        mcp.settings.log_level = resolved_log_level
        logging.getLogger().setLevel(resolved_log_level)

        if transport == "stdio":
            mcp.run(transport="stdio")
            return

        if transport == "sse":
            click.echo(
                "Warning: SSE transport is deprecated; migrate to streamable-http.",
                err=True,
            )

        transport_security = _resolve_transport_security(
            host,
            allowed_hosts,
            allowed_origins,
        )
        if transport == "streamable-http":
            mcp.run(
                transport="streamable-http",
                host=host,
                port=port,
                streamable_http_path="/mcp",
                stateless_http=True,
                transport_security=transport_security,
            )
        else:
            mcp.run(
                transport="sse",
                host=host,
                port=port,
                transport_security=transport_security,
            )

    except click.ClickException:
        raise
    except Exception as e:
        print(f"Error starting MCP server: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


def _resolve_transport_security(
    host: str,
    allowed_hosts: tuple[str, ...],
    allowed_origins: tuple[str, ...],
) -> TransportSecuritySettings | None:
    """Build explicit DNS-rebinding protection for non-loopback HTTP binds."""
    if host in _LOOPBACK_HOSTS:
        if not (allowed_hosts or allowed_origins):
            # MCPServer supplies the same strict loopback allowlist in this case.
            return None
        return TransportSecuritySettings(
            enable_dns_rebinding_protection=True,
            allowed_hosts=list(allowed_hosts or _LOOPBACK_ALLOWED_HOSTS),
            allowed_origins=list(allowed_origins or _LOOPBACK_ALLOWED_ORIGINS),
        )
    if not allowed_hosts:
        raise click.UsageError(
            "Non-loopback HTTP binding requires at least one --allowed-host value."
        )
    return TransportSecuritySettings(
        enable_dns_rebinding_protection=True,
        allowed_hosts=list(allowed_hosts),
        allowed_origins=list(allowed_origins),
    )


def main():
    """Main entry point for ChatSpatial CLI"""
    cli()


if __name__ == "__main__":  # pragma: no cover
    main()
