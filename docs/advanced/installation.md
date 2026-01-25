# Installation Guide

Detailed setup for all platforms and configurations.

**Quick start?** See [Quick Start](../quickstart.md) for 5-minute setup.

**Full guide?** See [INSTALLATION.md on GitHub](https://github.com/cafferychen777/ChatSpatial/blob/main/INSTALLATION.md) for comprehensive instructions including R dependencies and troubleshooting.

---

## Install Options

```bash
# Standard (recommended for most users)
pip install chatspatial

# Full (all methods including deep learning)
pip install chatspatial[full]
```

| Option | Features | Install Time |
|--------|----------|--------------|
| Standard | Core analysis, basic deconvolution, cell communication | 3-5 min |
| Full | + Deep learning, Cell2location, CellRank, SpaGCN, STAGATE | 6-10 min |

---

## Requirements

- **Python 3.11+** (3.12 recommended, 3.13 supported)
- **8GB+ RAM** (16GB+ for large datasets)
- **macOS, Linux, or Windows**

```bash
python --version  # Should be 3.11+
```

---

## Virtual Environment (Recommended)

```bash
# Create
python3 -m venv chatspatial_env
source chatspatial_env/bin/activate  # macOS/Linux
# chatspatial_env\Scripts\activate   # Windows

# Install
pip install chatspatial[full]
```

---

## Configure MCP

### Claude Code (Recommended)

```bash
claude mcp add chatspatial /path/to/chatspatial_env/bin/python -- -m chatspatial server
```

### Codex (CLI or IDE Extension)

MCP config is shared in `~/.codex/config.toml`.

**Option A: Add via CLI**

```bash
codex mcp add chatspatial -- /path/to/chatspatial_env/bin/python -m chatspatial server
```

**Option B: Edit config.toml**

```toml
[mcp_servers.chatspatial]
command = "/path/to/chatspatial_env/bin/python"
args = ["-m", "chatspatial", "server"]
```

Use `/mcp` in Codex TUI to verify the server is active.

### Claude Desktop

Edit config file (location varies by OS):
- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%\Claude\claude_desktop_config.json`
- Linux: `~/.config/Claude/claude_desktop_config.json`

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

**Restart your client after configuration.**

---

## Platform Notes

### Windows

**Not available:** SingleR, PETSc acceleration

**Alternatives:**
- Cell type annotation: Tangram, scANVI, CellAssign, mLLMCelltype
- CellRank: Works without PETSc

### R-based Methods (Optional)

For RCTD, SPOTlight, CARD, scType, SingleR, CellChat:

1. Install R 4.4+ from [CRAN](https://cran.r-project.org/)
2. Run our installer: `Rscript install_r_dependencies.R`

See [full installation guide](https://github.com/cafferychen777/ChatSpatial/blob/main/INSTALLATION.md#r-dependencies-for-r-based-methods) for details.

---

## Verify

```bash
python -m chatspatial server --help
python -c "import chatspatial; print('OK')"
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Import errors | `pip install --upgrade pip chatspatial[full]` |
| Claude not seeing tools | Check Python path, restart Claude |
| R methods failing | Install R and required packages |
| "mcp>=1.17.0 not found" | Python version too old, need 3.11+ |

See [Troubleshooting Guide](troubleshooting.md) for more.
