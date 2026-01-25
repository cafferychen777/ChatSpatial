# Installation Guide

Detailed setup for all platforms and configurations.

**Quick start?** See [Quick Start](../quickstart.md) for 5-minute setup.

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

### Claude Desktop

Edit `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS):

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```

> Windows: `%APPDATA%\Claude\claude_desktop_config.json`
> Linux: `~/.config/Claude/claude_desktop_config.json`

### Claude Code

```bash
claude mcp add chatspatial python -- -m chatspatial server
```

### Codex (CLI or IDE Extension)

MCP config is shared in `~/.codex/config.toml`.

**Option A: Add via CLI**

```bash
codex mcp add chatspatial -- python -m chatspatial server
```

**Option B: Edit config.toml**

```toml
[mcp_servers.chatspatial]
command = "python"
args = ["-m", "chatspatial", "server"]
```

**Virtual environment note:** To use a specific Python environment:

```bash
codex mcp add chatspatial -- /path/to/venv/bin/python -m chatspatial server
```

Or in config.toml:

```toml
[mcp_servers.chatspatial]
command = "/path/to/venv/bin/python"
args = ["-m", "chatspatial", "server"]
```

Use `/mcp` in Codex TUI to verify the server is active.

**Restart your client after configuration.**

---

## Platform Notes

### Windows

**Not available:** SingleR, PETSc acceleration

**Alternatives:**
- Cell type annotation: Tangram, scANVI, CellAssign, mLLMCelltype
- CellRank: Works without PETSc

### R-based Methods

Optional. Requires R 4.4+ installation:

| Method | R Package | Install Command |
|--------|-----------|-----------------|
| RCTD | spacexr | `install.packages("spacexr")` |
| SPOTlight | SPOTlight | `BiocManager::install("SPOTlight")` |
| CARD | CARD | `devtools::install_github("YingMa0107/CARD")` |
| scType | scType | `devtools::install_github("IanevskiAleksandr/scType")` |
| SingleR | SingleR | `BiocManager::install("SingleR")` |
| CellChat | CellChat | `devtools::install_github("jinworks/CellChat")` |
| SCTransform | sctransform | `install.packages("sctransform")` |

```bash
# Install R from https://cran.r-project.org/
# Install BiocManager for Bioconductor packages
R -e 'install.packages("BiocManager")'
```

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

See [Troubleshooting Guide](troubleshooting.md) for more.
