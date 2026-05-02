# Configuration Guide

This page is the canonical reference for **exact MCP client configuration syntax**.

- To install ChatSpatial, see [Installation](../installation.md).
- To run your first workflow after setup, see [Quick Start](../quickstart.md).
- If configuration fails, see [Troubleshooting](troubleshooting.md).

---

## Configuration Workflow

1. Choose a runtime: Python environment or Docker/GHCR image.
2. For Python, activate the environment and run `which python`.
3. Use the Python command or Docker command shape in your MCP client config.
4. Restart the client after configuration changes.
5. Verify the server can start.

Canonical Python command shape:

```text
/absolute/path/to/python -m chatspatial server
```

Canonical Docker command shape:

```bash
docker run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:latest server --transport stdio
```

Use `--rm -i`, not `-it`, for MCP stdio. If you mount host data to `/data`, prompts must use container paths such as `/data/sample.h5ad`.

---

## Claude Code (Recommended)

```bash
source venv/bin/activate
which python
claude mcp add chatspatial /path/to/venv/bin/python -- -m chatspatial server
claude mcp list
```

**Notes:**
- `--` separates the Python path from module arguments
- use the absolute Python path from `which python`
- use `--scope user` if you want the server available across projects

### Docker-backed Claude Code server

```bash
claude mcp add chatspatial-docker docker -- \
  run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:latest server --transport stdio
```

Use `/data/...` paths in prompts when using this Docker-backed server.

---

## Codex

Codex stores MCP configuration in `~/.codex/config.toml`.

### Add via CLI

```bash
source venv/bin/activate
which python
codex mcp add chatspatial -- /path/to/venv/bin/python -m chatspatial server
```

### Or edit config directly

```toml
[mcp_servers.chatspatial]
command = "/path/to/venv/bin/python"
args = ["-m", "chatspatial", "server"]

[mcp_servers.chatspatial.env]
CHATSPATIAL_DATA_DIR = "/path/to/data"
```

### Advanced options

```toml
[mcp_servers.chatspatial]
command = "/path/to/venv/bin/python"
args = ["-m", "chatspatial", "server"]
startup_timeout_sec = 30
tool_timeout_sec = 120
enabled = true
```

---

## OpenCode

OpenCode stores MCP configuration in:

- global: `~/.config/opencode/opencode.json`
- project: `opencode.json`

Project config takes precedence when both exist.

### Add via CLI

```bash
opencode mcp add
opencode mcp list
```

### Or edit config directly

```json
{
  "$schema": "https://opencode.ai/config.json",
  "mcp": {
    "chatspatial": {
      "type": "local",
      "command": ["/path/to/venv/bin/python", "-m", "chatspatial", "server"],
      "enabled": true,
      "environment": {
        "CHATSPATIAL_DATA_DIR": "/path/to/data"
      }
    }
  }
}
```

**Notes:**
- `command` is an array: `[executable, ...args]`
- use the **absolute** Python path from `which python`
- prefer project-level config for repo-specific settings

---

## Claude Desktop

Edit the Claude Desktop config file:

- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`
- Linux: `~/.config/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```

Docker-backed example:

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

Restart Claude Desktop after saving the file.

---

## Other MCP Clients

ChatSpatial works with any MCP-compatible client.

Minimum requirement:
- configure the executable as your environment’s Python
- pass `-m chatspatial server` as arguments

Use the same absolute Python path pattern shown above.

---

## Environment Variables

### Data storage

```bash
export CHATSPATIAL_DATA_DIR="/path/to/your/spatial/data"
```

When `export_data()` is called without an explicit `path`, ChatSpatial saves to this directory.

Default behavior: `.chatspatial_saved/` next to the original data file.

---

## Verify Configuration

```bash
which python
python -c "import chatspatial; print(f'ChatSpatial {chatspatial.__version__} ready')"
python -m chatspatial server --help
```

If these checks fail, use [Troubleshooting](troubleshooting.md).

---

## Next Steps

- [Docker / GHCR](../docker.md) — container-backed runtime setup
- [Quick Start](../quickstart.md) — first successful analysis
- [Troubleshooting](troubleshooting.md) — fix configuration or runtime issues
- [Methods Reference](methods-reference.md) — exact parameters and defaults
