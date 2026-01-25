# Configuration Guide

Configure ChatSpatial MCP server for your environment.

---

## MCP Client Configuration

### Claude Code (Recommended)

```bash
# Find your virtual environment Python path
source chatspatial_env/bin/activate
which python
# Copy the output path

# Add ChatSpatial MCP server
claude mcp add chatspatial /path/to/chatspatial_env/bin/python -- -m chatspatial server

# Verify connection
claude mcp list
# Should show: chatspatial: ... - Connected
```

**Key points:**
- The `--` separates the Python path from the module arguments
- Always use absolute path from `which python`
- Use `--scope user` to make ChatSpatial available across all projects

---

### Codex (CLI and IDE Extension)

Codex stores MCP configuration in `~/.codex/config.toml`. The CLI and IDE extension share this configuration.

**Add via CLI:**

```bash
# Find your virtual environment Python path
source chatspatial_env/bin/activate
which python
# Copy the output path

# Add ChatSpatial MCP server
codex mcp add chatspatial -- /path/to/chatspatial_env/bin/python -m chatspatial server

# Verify in Codex TUI
/mcp
```

**Or edit `~/.codex/config.toml` directly:**

```toml
[mcp_servers.chatspatial]
command = "/path/to/chatspatial_env/bin/python"
args = ["-m", "chatspatial", "server"]

# Optional: Environment variables
[mcp_servers.chatspatial.env]
CHATSPATIAL_DATA_DIR = "/path/to/data"
```

**Advanced options:**

```toml
[mcp_servers.chatspatial]
command = "/path/to/chatspatial_env/bin/python"
args = ["-m", "chatspatial", "server"]
startup_timeout_sec = 30    # Default: 10
tool_timeout_sec = 120      # Default: 60
enabled = true              # Set to false to disable without deleting
```

**Key points:**
- Use `[mcp_servers.chatspatial]` (underscore, not hyphen)
- Configuration is shared between CLI and IDE extension
- Use `/mcp` in Codex TUI to verify connection

---

### Claude Desktop

Edit Claude Desktop configuration file:

**Location:**
- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%/Claude/claude_desktop_config.json`
- Linux: `~/.config/Claude/claude_desktop_config.json`

**Configuration:**

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```

**Example with actual path:**

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/Users/yourname/Projects/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```

**Important:** Restart Claude Desktop after configuration changes.

---

### Other MCP Clients

For other MCP-compatible clients:

1. **Find Python path:** Activate virtual environment and run `which python`
2. **Use in configuration:** Replace `"python"` with the full path

ChatSpatial follows the standard MCP protocol and works with any MCP-compatible client.

---

## Environment Variables (Optional)

Configure ChatSpatial behavior using environment variables:

### Data Storage

```bash
# Set custom data directory for saved datasets
export CHATSPATIAL_DATA_DIR="/path/to/your/spatial/data"
```

**Usage:** When you use `export_data()` without specifying `path`, datasets are saved to this directory.

**Default:** `.chatspatial_saved/` next to original data file

---

## Troubleshooting Configuration

### Common Issues

| Problem | Solution |
|---------|----------|
| "python not found" | Use full path to virtual environment Python |
| "module not found" | Ensure virtual environment is activated before adding server |
| Claude can't connect | Check JSON syntax and restart Claude Desktop |
| Server not showing up | Verify Python path is correct with `which python` |

### Verify Configuration

```bash
# Make sure you're in the virtual environment
which python
# Should show virtual environment path, not system Python

# Test ChatSpatial import
python -c "import chatspatial; print(f'ChatSpatial {chatspatial.__version__} ready')"

# Test MCP server
python -m chatspatial server --help
# Should display server options
```

---

## Next Steps

- [Quick Start](../quickstart.md) - Start analyzing data
- [Troubleshooting](troubleshooting.md) - Solve common problems
- [Methods Reference](methods-reference.md) - Explore available tools
