
title: Installation
---

# Installation Guide

## Table of contents

## Quick Start

### 1. Create Virtual Environment (Strongly Recommended)

```bash
# Clone repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd chatspatial

# Create virtual environment (choose one method)
# Option A: Using venv (Python built-in)
python3 -m venv chatspatial_env
source chatspatial_env/bin/activate  # On macOS/Linux
# chatspatial_env\Scripts\activate    # On Windows

# Option B: Using conda
conda create -n chatspatial python=3.10
conda activate chatspatial

# Option C: Using uv (fast and modern)
uv venv
source .venv/bin/activate
```

### 2. Install ChatSpatial

```bash
# Recommended: Install with all features (in activated virtual environment)
pip install -e ".[full]"
```

💡 **Why virtual environment?** Isolates dependencies, prevents conflicts, and makes configuration with MCP clients cleaner.

💡 **Why [full]?** Enables all 16+ analysis methods. **Installation time:** 6-10 minutes depending on your internet speed and system (includes compiling PETSc/SLEPc on first install).

## Installation Options

| Option | Command | Features | Time |
|--------|---------|----------|------|
| **Full (Recommended)** | `pip install -e ".[full]"` | 100% features | 6-10 min |
| Standard | `pip install -e .` | 80% features | 3-5 min |

### Standard Installation (Default)
- ✅ Core spatial analysis (Moran's I, Getis-Ord)
- ✅ Basic deconvolution with scvi-tools
- ✅ RNA velocity with scVelo
- ✅ Cell communication (LIANA, CellPhoneDB)
- ✅ Batch integration (Harmony, BBKNN)
- ✅ Gene enrichment analysis

### Full Installation (Recommended)
Everything in Standard, plus:
- ✅ Deep learning methods (PyTorch)
- ✅ Advanced deconvolution (Cell2location, CARD)
- ✅ Advanced trajectory (CellRank, Palantir)
- ✅ Spatial domains (SpaGCN, STAGATE, GraphST)
- ✅ Spatial variable genes (SPARK-X, SpatialDE)
- ✅ CNV analysis (inferCNVpy, Numbat)

## Requirements

- Python 3.10-3.13 (MCP requires Python 3.10+)
- 5-10 GB disk space (for full installation)
- macOS, Linux, or Windows

### Check Your Python Version

```bash
# Check your Python version
python --version
# or
python3 --version

# Should output: Python 3.10.x, 3.11.x, 3.12.x, or 3.13.x
```

💡 **Python 3.13 Support:** ChatSpatial is fully compatible with Python 3.13. Some minor FutureWarning messages from dependencies may appear but do not affect functionality.

### Platform-Specific Limitations

**Windows Users:** SingleR and PETSc/SLEPc acceleration are not available on Windows due to C++ compilation limitations.

#### Windows Limitations

**❌ Not Available on Windows:**
- **SingleR** cell type annotation - C++ dependencies (mattress, knncolle) fail to compile due to MinGW compiler limitations
- **PETSc/SLEPc** acceleration for CellRank - Requires Cygwin Python (not native Windows Python)

**✅ Windows Alternatives:**
- Cell type annotation: Use Tangram, scANVI, CellAssign, or mllmcelltype instead of SingleR
- CellRank: Works without PETSc (automatically uses 'brandts' method for small-medium datasets)

**✅ All other features work on Windows** including:
- R-based methods (RCTD, SPOTlight, Numbat) via rpy2
- All deconvolution methods (Cell2location, DestVI, Stereoscope, CARD)
- All trajectory and spatial analysis methods

**Technical Note:** GitHub Actions CI installs R 4.4.1 on all platforms, enabling rpy2 compilation on Windows. However, SingleR's C++ dependencies still fail due to compiler limitations.

## Configure Your MCP Client

**Important:** Use your virtual environment's Python path for reliable operation.

### Find Your Virtual Environment Python Path

```bash
# After activating your virtual environment
which python
# This will output something like:
# /Users/yourname/Projects/chatspatial_env/bin/python
```

### Option A: Claude Desktop (GUI Application)

Edit Claude Desktop configuration file:
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

**Real Example:**
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

### Option B: Claude Code (Terminal/IDE)

Install Claude Code CLI and add ChatSpatial:

```bash
# Install Claude Code CLI (if not installed)
npm install -g @anthropic-ai/claude-code

# Step 1: Find your virtual environment Python path
# (Make sure your virtual environment is activated first)
which python
# Copy the output path (e.g., /Users/yourname/chatspatial_env/bin/python)

# Step 2: Add ChatSpatial MCP server using the FULL path
claude mcp add chatspatial -- /Users/yourname/chatspatial_env/bin/python -m chatspatial server

# Step 3: Verify the server is connected
claude mcp list
# You should see "chatspatial: ... - ✓ Connected"
```

**Real Example:**
```bash
# After running "which python" in activated environment
# Output: /Users/alice/Projects/chatspatial_env/bin/python

# Add the server with that exact path:
claude mcp add chatspatial -- /Users/alice/Projects/chatspatial_env/bin/python -m chatspatial server
```

💡 **Key points:**
- The `--` separates Claude CLI options from the server command
- Always use the **absolute path** from `which python` (not relative paths)
- Use `--scope user` to make ChatSpatial available across all projects: `claude mcp add --scope user chatspatial -- ...`

⚠️ **Never use system Python:** Always use the Python from your virtual environment to avoid dependency conflicts. Test with `which python` to ensure you're using the right one.

## Verify Installation

### Step 1: Check Virtual Environment

```bash
# First, ensure you're in the virtual environment
which python
# Should output path to your virtual environment, not system Python
# Example: /Users/yourname/chatspatial_env/bin/python

# Verify Python version
python --version
# Should show Python 3.10.x, 3.11.x, 3.12.x, or 3.13.x
```

### Step 2: Test ChatSpatial Installation

```bash
# Test the command-line interface
python -m chatspatial server --help
# Should display server options without errors

# Test Python import
python -c "import chatspatial; print(f'✅ ChatSpatial {chatspatial.__version__} installed successfully')"
```

### Step 3: Test Core Dependencies

```bash
# Quick functionality test
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

**Expected output:** All commands should complete without errors. You may see some FutureWarning messages (especially with Python 3.13), which are normal and don't affect functionality.

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Import errors | `pip install --upgrade pip` |
| Package conflicts | `pip install --force-reinstall -e ".[full]"` |
| Claude Desktop doesn't see server | Check that command path points to virtual environment Python |
| Claude Code "command not found" | Install CLI: `npm install -g @anthropic-ai/claude-code` |
| Claude Code connection error | Use absolute path to Python; run `claude mcp list` to verify |
| Wrong Python version | Recreate virtual environment with correct Python version |

## Getting Help

- **GitHub Issues**: [Report problems](https://github.com/cafferychen777/ChatSpatial/issues)
- **Documentation**: Check docstrings with `help(function_name)`

---

**Next:** [Quick Start Guide](quick-start.md) to run your first analysis