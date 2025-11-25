# Installation Guide

Detailed installation instructions and platform-specific guidance.

**Quick start?** See [Quick Start](../quickstart.md) for 5-minute setup.

---

## Installation Options

| Option | Command | Features | Time |
|--------|---------|----------|------|
| **Full (Recommended)** | `pip install -e ".[full]"` | 100% features | 6-10 min |
| Standard | `pip install -e .` | 80% features | 3-5 min |

### Standard Installation
- Core spatial analysis (Moran's I, Getis-Ord)
- Basic deconvolution with scvi-tools
- RNA velocity with scVelo
- Cell communication (LIANA, CellPhoneDB)
- Batch integration (Harmony, BBKNN)
- Gene enrichment analysis

### Full Installation (Recommended)
Everything in Standard, plus:
- Deep learning methods (PyTorch)
- Advanced deconvolution (Cell2location, CARD)
- Advanced trajectory (CellRank, Palantir)
- Spatial domains (SpaGCN, STAGATE, GraphST)
- Spatial variable genes (SPARK-X, SpatialDE)
- CNV analysis (inferCNVpy, Numbat)

---

## Requirements

- **Python 3.10-3.13** (MCP requires 3.10+)
- **5-10 GB disk space** (for full installation)
- **macOS, Linux, or Windows**

### Check Python Version

```bash
python --version
# Should show: Python 3.10.x, 3.11.x, 3.12.x, or 3.13.x
```

**Python 3.13 Support:** ChatSpatial is fully compatible with Python 3.13. Some FutureWarning messages from dependencies may appear but do not affect functionality.

---

## Platform-Specific Limitations

### Windows Limitations

**Not Available on Windows:**
- **SingleR** cell type annotation - C++ compilation issues
- **PETSc/SLEPc** acceleration for CellRank

**Windows Alternatives:**
- Cell type annotation: Use Tangram, scANVI, CellAssign, or mllmcelltype instead of SingleR
- CellRank: Works without PETSc (uses 'brandts' method for small-medium datasets)

**All other features work on Windows** including:
- R-based methods (RCTD, SPOTlight, Numbat) via rpy2
- All deconvolution methods (Cell2location, DestVI, Stereoscope, CARD)
- All trajectory and spatial analysis methods

---

## Detailed Installation Steps

### 1. Create Virtual Environment

**Using venv (Python built-in):**
```bash
git clone https://github.com/cafferychen777/ChatSpatial.git
cd chatspatial

python3 -m venv chatspatial_env
source chatspatial_env/bin/activate  # macOS/Linux
# chatspatial_env\Scripts\activate  # Windows
```

**Using conda:**
```bash
conda create -n chatspatial python=3.10
conda activate chatspatial
```

**Using uv (fast and modern):**
```bash
uv venv
source .venv/bin/activate
```

**Why virtual environment?** Isolates dependencies, prevents conflicts, and makes MCP configuration cleaner.

### 2. Install ChatSpatial

```bash
# Full installation (recommended)
pip install -e ".[full]"

# Standard installation (faster, fewer features)
pip install -e .
```

---

## Configure MCP Client

**Important:** Use your virtual environment's Python path.

### Find Virtual Environment Python Path

```bash
# After activating your virtual environment
which python
# Output: /Users/yourname/Projects/chatspatial_env/bin/python
```

### Claude Desktop

Edit configuration file:
- **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
- **Linux**: `~/.config/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/your/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```

**Restart Claude Desktop** after configuration changes.

### Claude Code

```bash
# Find Python path
source chatspatial_env/bin/activate
which python  # Copy this path

# Add ChatSpatial MCP server
claude mcp add chatspatial -- /Users/yourname/chatspatial_env/bin/python -m chatspatial server

# Verify
claude mcp list
# Should show: chatspatial: ... - ✓ Connected
```

**Key points:**
- The `--` separates Claude CLI options from server command
- Always use **absolute path** from `which python`
- Use `--scope user` for global availability: `claude mcp add --scope user chatspatial -- ...`

---

## Verify Installation

### Step 1: Check Virtual Environment

```bash
# Ensure you're in virtual environment
which python
# Should show: /path/to/chatspatial_env/bin/python (NOT system Python)

# Check Python version
python --version
# Should show: Python 3.10.x, 3.11.x, 3.12.x, or 3.13.x
```

### Step 2: Test ChatSpatial

```bash
# Test command-line interface
python -m chatspatial server --help
# Should display server options

# Test Python import
python -c "import chatspatial; print(f'✅ ChatSpatial {chatspatial.__version__} ready')"
```

### Step 3: Test Core Dependencies

```bash
python -c "
import scanpy as sc
import squidpy as sq
import numpy as np
print('✅ Core packages working:')
print(f'  - scanpy: {sc.__version__}')
print(f'  - squidpy: {sq.__version__}')
print(f'  - numpy: {np.__version__}')
"
```

**Expected:** All commands complete without errors. FutureWarning messages (especially with Python 3.13) are normal.

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Import errors | `pip install --upgrade pip` then reinstall |
| Package conflicts | `pip install --force-reinstall -e ".[full]"` |
| Claude Desktop doesn't see server | Verify command path points to virtual environment Python |
| Claude Code "command not found" | Install CLI: `npm install -g @anthropic-ai/claude-code` |
| Claude Code connection error | Use absolute path; verify with `claude mcp list` |
| Wrong Python version | Recreate virtual environment with correct Python |

---

## Getting Help

- **GitHub Issues**: [Report problems](https://github.com/cafferychen777/ChatSpatial/issues)
- **Troubleshooting Guide**: [Advanced troubleshooting](troubleshooting.md)
- **FAQ**: [Common questions](faq.md)

---

**Next:** [Quick Start](../quickstart.md) to run your first analysis
