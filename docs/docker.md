# Docker / GHCR

Use the Docker image when you want to run ChatSpatial without resolving the scientific Python dependency stack on your host machine.

This is the canonical Docker runtime guide. Other pages may remind you to use
container paths such as `/data/sample.h5ad`, but Docker image tags, mounts,
`--rm -i`, `/outputs`, HTTP mode, and Docker-specific failures are maintained
here.

## Image

```text
ghcr.io/cafferychen777/chatspatial:v1.3.0
```

This guide uses the versioned `v1.3.0` tag so commands are reproducible. Use
`:latest` only when you explicitly want the newest published image.

## Pull and verify the image

```bash
docker pull ghcr.io/cafferychen777/chatspatial:v1.3.0
docker run --rm ghcr.io/cafferychen777/chatspatial:v1.3.0 --version
docker run --rm ghcr.io/cafferychen777/chatspatial:v1.3.0 server --help
```

## Run as an MCP server

Use this command shape in MCP clients that support Docker-backed stdio servers:

```bash
docker run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:v1.3.0 server --transport stdio
```

Use `--rm -i`, not `-it`, for MCP stdio. A TTY can corrupt the JSON-RPC stream used by MCP clients.

## Mount data and outputs

ChatSpatial runs inside the container. It cannot see host files unless they are mounted into the container.

| Host path | Docker mount | Path to use in prompts |
|---|---|---|
| `/Users/alice/spatial-data` | `-v /Users/alice/spatial-data:/data:ro` | `/data/sample.h5ad` |
| `/home/alice/chatspatial-outputs` | `-v /home/alice/chatspatial-outputs:/outputs` | generated files appear under `/outputs` |

If your host data directory is mounted as `/data`, prompt ChatSpatial with `/data/sample.h5ad`, not the original host path.

Generated visualizations are written under `/outputs` by default.
If you call `export_data()` without an explicit path, it uses
`~/.chatspatial/active/{data_id}.h5ad` inside the container; pass a path under
`/outputs` when the exported file should persist after the container exits.

## MCP client examples

### Claude Code

```bash
claude mcp add chatspatial-docker docker -- \
  run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:v1.3.0 server --transport stdio
```

### Claude Desktop

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "docker",
      "args": [
        "run",
        "--rm",
        "-i",
        "-v",
        "/absolute/path/to/your/data:/data:ro",
        "-v",
        "/absolute/path/to/outputs:/outputs",
        "ghcr.io/cafferychen777/chatspatial:v1.3.0",
        "server",
        "--transport",
        "stdio"
      ]
    }
  }
}
```

Restart your MCP client after changing configuration.

## Streamable HTTP server

```bash
docker run --rm -p 127.0.0.1:8000:8000 \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:v1.3.0 server \
  --transport streamable-http \
  --host 0.0.0.0 \
  --port 8000 \
  --allowed-host 'localhost:*' \
  --allowed-host '127.0.0.1:*' \
  --allowed-origin 'http://localhost:*' \
  --allowed-origin 'http://127.0.0.1:*'
```

Then connect your MCP client to `http://localhost:8000/mcp`.

The port mapping above binds only to the host loopback interface. ChatSpatial's
HTTP endpoint does not add application authentication, so do not expose it to an
untrusted network. For remote access, place it behind an authenticated reverse
proxy or a private network and replace the Host/Origin allowlists with the exact
public values.

MCP `2026-07-28` is stateless at the protocol layer, but loaded AnnData objects
remain process-local application state. Run one ChatSpatial worker per endpoint;
`data_id` handles remain valid for that process lifetime. Export important
datasets to `/outputs` before restarting or replacing the container.

## Build locally from source

From the repository root:

```bash
docker build -t chatspatial:local .
```

Then replace `ghcr.io/cafferychen777/chatspatial:v1.3.0` with
`chatspatial:local` in the commands above.

## Common Docker issues

| Symptom | Fix |
|---|---|
| `docker: command not found` | Install Docker Desktop or Docker Engine, then restart your MCP client. |
| Pull fails | Check the image name and network access with `docker pull ghcr.io/cafferychen777/chatspatial:v1.3.0`. |
| MCP tools do not appear | Use `--rm -i`, not `-it`, and restart the client. |
| Dataset not found | Mount the host data directory and use the container path in prompts, for example `/data/sample.h5ad`. |
| Permission denied on outputs | Confirm the host output directory exists and Docker has permission to write to it. |
| Works in terminal but not in the client | Use absolute host paths in `-v` mounts and restart the client. |

## Next steps

- [Configuration Guide](advanced/configuration.md) — exact MCP client syntax
- [Quick Start](quickstart.md) — first analysis prompt
- [Troubleshooting](advanced/troubleshooting.md) — symptom-based fixes
