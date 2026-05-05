# Docker / GHCR

Use the Docker image when you want to run ChatSpatial without resolving the scientific Python dependency stack on your host machine.

## Image

```text
ghcr.io/cafferychen777/chatspatial:v1.2.6
```

This guide uses the versioned `v1.2.6` tag so commands are reproducible. Use
`:latest` only when you explicitly want the newest published image.

## Pull and verify the image

```bash
docker pull ghcr.io/cafferychen777/chatspatial:v1.2.6
docker run --rm ghcr.io/cafferychen777/chatspatial:v1.2.6 --version
docker run --rm ghcr.io/cafferychen777/chatspatial:v1.2.6 server --help
```

## Run as an MCP server

Use this command shape in MCP clients that support Docker-backed stdio servers:

```bash
docker run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:v1.2.6 server --transport stdio
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
  ghcr.io/cafferychen777/chatspatial:v1.2.6 server --transport stdio
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
        "ghcr.io/cafferychen777/chatspatial:v1.2.6",
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
  ghcr.io/cafferychen777/chatspatial:v1.2.6 server --transport sse --host 0.0.0.0 --port 8000
```

Then connect your MCP client to `http://localhost:8000/sse`.

## Build locally from source

From the repository root:

```bash
docker build -t chatspatial:local .
```

Then replace `ghcr.io/cafferychen777/chatspatial:v1.2.6` with
`chatspatial:local` in the commands above.

## Common Docker issues

| Symptom | Fix |
|---|---|
| `docker: command not found` | Install Docker Desktop or Docker Engine, then restart your MCP client. |
| Pull fails | Check the image name and network access with `docker pull ghcr.io/cafferychen777/chatspatial:v1.2.6`. |
| MCP tools do not appear | Use `--rm -i`, not `-it`, and restart the client. |
| Dataset not found | Mount the host data directory and use the container path in prompts, for example `/data/sample.h5ad`. |
| Permission denied on outputs | Confirm the host output directory exists and Docker has permission to write to it. |
| Works in terminal but not in the client | Use absolute host paths in `-v` mounts and restart the client. |

## Next steps

- [Configuration Guide](advanced/configuration.md) — exact MCP client syntax
- [Quick Start](quickstart.md) — first analysis prompt
- [Troubleshooting](advanced/troubleshooting.md) — symptom-based fixes
