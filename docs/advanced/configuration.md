# Configuration Guide

This page is the canonical reference for **exact MCP client configuration syntax**.
ChatSpatial works with any MCP-compatible client; Claude Code, Codex, OpenCode,
and Claude Desktop are examples.

- To install ChatSpatial, see [Installation](../installation.md).
- To run your first workflow after setup, see [Quick Start](../quickstart.md).
- If configuration fails, see [Troubleshooting](troubleshooting.md).

---

## Configuration Workflow

1. Choose a runtime: `uvx`, a persistent Python environment, or Docker/GHCR.
2. Prefer `uvx` unless you need a customized persistent environment.
3. Use the corresponding command shape in your MCP client config.
4. Restart the client after configuration changes.
5. Verify the server can start.

Canonical `uvx` command shape:

```text
uvx --from chatspatial chatspatial server
```

To expose every composable Python method family, put the `full` extra on the
`--from` package specification:

```text
uvx --from 'chatspatial[full]' chatspatial server
```

Use `chatspatial[full]==1.5.0` when the MCP environment must be reproducible.
R bridges, AESTETIK, and rctd-py are separate extras; see the installation
guide before adding them.

Canonical persistent-Python command shape:

```text
/absolute/path/to/python -m chatspatial server
```

Canonical Docker command shape:

```bash
docker run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:v1.5.0 server --transport stdio
```

Use `--rm -i`, not `-it`, for MCP stdio. If you mount host data to `/data`, prompts must use container paths such as `/data/sample.h5ad`.

---

## Path Model

ChatSpatial has three separate path concepts:

| Concept | What it means | Default |
|---------|---------------|---------|
| **Input data path** | Path passed to `load_data`; it must be visible to the ChatSpatial runtime | No search path; use an absolute path |
| **Active export/reload path** | Default path used by `export_data()` and `reload_data()` when `path` is omitted | `~/.chatspatial/active/{data_id}.h5ad` inside the runtime |
| **Visualization/output path** | Where generated figures and explicit output files are written | `CHATSPATIAL_OUTPUT_DIR` if set, otherwise a safe writable directory |

For local `uvx` and Python runtimes, input paths are normal host paths such as
`/Users/alice/spatial/sample.h5ad`.

For Docker, the runtime sees container paths. If you mount
`/Users/alice/spatial-data` as `/data`, prompts must use `/data/sample.h5ad`.
The Docker image sets `CHATSPATIAL_OUTPUT_DIR=/outputs`, so mount a writable
host directory to `/outputs` when you want generated files to persist.

Do not rely on a global data search directory. Keep data locations explicit in
prompts or tool calls.

---

## Claude Code

```bash
claude mcp add --scope user chatspatial -- \
  uvx --from chatspatial chatspatial server
claude mcp list
```

**Notes:**
- `--` separates Claude Code options from the server command
- `uvx` creates and caches an isolated ChatSpatial environment
- remove `--scope user` or choose another scope for project-specific setup

For a persistent environment, replace the command after `--` with
`/path/to/venv/bin/python -m chatspatial server`.

### Docker-backed Claude Code server

```bash
claude mcp add chatspatial-docker docker -- \
  run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:v1.5.0 server --transport stdio
```

Use `/data/...` paths in prompts when using this Docker-backed server.

---

## Codex

Codex stores MCP configuration in `~/.codex/config.toml`.

### Add via CLI

```bash
codex mcp add chatspatial -- uvx --from chatspatial chatspatial server
```

### Or edit config directly

```toml
[mcp_servers.chatspatial]
command = "uvx"
args = ["--from", "chatspatial", "chatspatial", "server"]
```

The ChatGPT desktop app, Codex CLI, and Codex IDE extension share this Codex
configuration. Restart the active client after changing it. For a persistent
environment, use its absolute Python path as `command` and
`["-m", "chatspatial", "server"]` as `args`.

In the ChatGPT desktop app, you can configure the same command under
**Settings → MCP servers → Add server**: choose **STDIO**, enter `uvx` as the
command, and enter `--from chatspatial chatspatial server` as the arguments.
See the [official OpenAI MCP documentation](https://learn.chatgpt.com/docs/extend/mcp)
for the current UI and shared-configuration behavior.

### Advanced options

```toml
[mcp_servers.chatspatial]
command = "uvx"
args = ["--from", "chatspatial", "chatspatial", "server"]
startup_timeout_sec = 30
tool_timeout_sec = 120
enabled = true

[mcp_servers.chatspatial.env]
CHATSPATIAL_OUTPUT_DIR = "/absolute/path/to/chatspatial-outputs"
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
      "command": ["uvx", "--from", "chatspatial", "chatspatial", "server"],
      "enabled": true,
      "environment": {
        "CHATSPATIAL_OUTPUT_DIR": "/absolute/path/to/chatspatial-outputs"
      }
    }
  }
}
```

**Notes:**
- `command` is an array: `[executable, ...args]`
- confirm `uvx --version` works before restarting OpenCode
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
      "command": "uvx",
      "args": ["--from", "chatspatial", "chatspatial", "server"]
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
        "ghcr.io/cafferychen777/chatspatial:v1.5.0",
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
- configure the executable as `uvx`
- pass `--from chatspatial chatspatial server` as arguments

Use a persistent Python command only when you need a custom environment.

---

## Streamable HTTP

STDIO is preferred for local desktop and coding clients. Use Streamable HTTP
only when a client cannot launch ChatSpatial as a subprocess or when the server
runs in a separately managed environment.

Loopback-only server:

```bash
python -m chatspatial server \
  --transport streamable-http \
  --host 127.0.0.1 \
  --port 8000
```

Connect to `http://127.0.0.1:8000/mcp`. The MCP SDK supplies a strict loopback
Host/Origin allowlist automatically.

Binding to a non-loopback interface requires explicit Host values:

```bash
python -m chatspatial server \
  --transport streamable-http \
  --host 0.0.0.0 \
  --port 8000 \
  --allowed-host 'chatspatial.example.org:*' \
  --allowed-origin 'https://chatspatial.example.org'
```

The built-in HTTP endpoint has no application authentication. Put non-local
deployments behind an authenticated reverse proxy or private network, terminate
TLS there, and use exact Host/Origin allowlists. The legacy `sse` transport is
retained for compatibility but is deprecated.

The MCP protocol is stateless, while ChatSpatial's loaded AnnData objects are
intentionally process-local. Use one worker per endpoint. A `data_id` is valid
for the lifetime of that process; call `export_data()` before a restart and
`reload_data()` afterward when the analysis must continue.

---

## Output Configuration

Set `CHATSPATIAL_OUTPUT_DIR` when you want generated figures and default output
files in a predictable location:

```bash
export CHATSPATIAL_OUTPUT_DIR="/absolute/path/to/chatspatial-outputs"
```

When `export_data()` or `reload_data()` is called without an explicit `path`,
ChatSpatial uses `~/.chatspatial/active/{data_id}.h5ad` inside the runtime. In
Docker, pass an explicit `/outputs/...` path if the exported file should persist
after the container exits.

---

## Verify Configuration

```bash
uvx --version
uvx --from chatspatial chatspatial --version
```

For a persistent environment, also verify its absolute Python path and run
`python -m chatspatial server --help`. If these checks fail, use
[Troubleshooting](troubleshooting.md).

---

## Next Steps

- [Docker / GHCR](../docker.md) — container-backed runtime setup
- [Quick Start](../quickstart.md) — first successful analysis
- [Troubleshooting](troubleshooting.md) — fix configuration or runtime issues
- [Methods Reference](methods-reference.md) — exact parameters and defaults
