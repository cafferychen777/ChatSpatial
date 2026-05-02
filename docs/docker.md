# Docker / GHCR

Use the Docker image when you want to run ChatSpatial without resolving the scientific Python dependency stack on your host machine.

## Image

```text
ghcr.io/cafferychen777/chatspatial:latest
```

`latest` is convenient for first-time setup. For reproducible analyses, prefer a versioned image tag once one is available.

## Pull and verify the image

```bash
docker pull ghcr.io/cafferychen777/chatspatial:latest
docker run --rm ghcr.io/cafferychen777/chatspatial:latest --version
docker run --rm ghcr.io/cafferychen777/chatspatial:latest server --help
```

## Run as an MCP server

Use this command shape in MCP clients that support Docker-backed stdio servers:

```bash
docker run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:latest server --transport stdio
```

Use `--rm -i`, not `-it`, for MCP stdio. A TTY can corrupt the JSON-RPC stream used by MCP clients.

## Mount data and outputs

ChatSpatial runs inside the container. It cannot see host files unless they are mounted into the container.

| Host path | Docker mount | Path to use in prompts |
|---|---|---|
| `/Users/alice/spatial-data` | `-v /Users/alice/spatial-data:/data:ro` | `/data/sample.h5ad` |
| `/home/alice/chatspatial-outputs` | `-v /home/alice/chatspatial-outputs:/outputs` | generated files appear under `/outputs` |

If your host data directory is mounted as `/data`, prompt ChatSpatial with `/data/sample.h5ad`, not the original host path.

Generated visualizations and exported files are written under `/outputs` by default.

## MCP client examples

### Claude Code

```bash
claude mcp add chatspatial-docker docker -- \
  run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:latest server --transport stdio
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
        "ghcr.io/cafferychen777/chatspatial:latest",
        "server",
        "--transport",
        "stdio"
      ]
    }
  }
}
```

Restart your MCP client after changing configuration.

## SSE server

```bash
docker run --rm -p 8000:8000 \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:latest server --transport sse --host 0.0.0.0 --port 8000
```

Then connect your MCP client to `http://localhost:8000/sse`.

## Common Docker issues

| Symptom | Fix |
|---|---|
| `docker: command not found` | Install Docker Desktop or Docker Engine, then restart your MCP client. |
| Pull fails | Check the image name and network access with `docker pull ghcr.io/cafferychen777/chatspatial:latest`. |
| MCP tools do not appear | Use `--rm -i`, not `-it`, and restart the client. |
| Dataset not found | Mount the host data directory and use the container path in prompts, for example `/data/sample.h5ad`. |
| Permission denied on outputs | Confirm the host output directory exists and Docker has permission to write to it. |
| Works in terminal but not in the client | Use absolute host paths in `-v` mounts and restart the client. |

## Next steps

- [Configuration Guide](advanced/configuration.md) — exact MCP client syntax
- [Quick Start](quickstart.md) — first analysis prompt
- [Troubleshooting](advanced/troubleshooting.md) — symptom-based fixes
