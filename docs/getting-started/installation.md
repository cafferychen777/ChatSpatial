---
layout: default
title: Installation
parent: Getting Started
nav_order: 1
---

# Installation Guide
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

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

{: .highlight }
💡 **Why virtual environment?** Isolates dependencies, prevents conflicts, and makes configuration with MCP clients cleaner.

{: .highlight }
💡 **Why [full]?** Enables all 16+ analysis methods. Takes ~13 minutes but provides the complete ChatSpatial experience.

## Installation Options

| Option | Command | Features | Time |
|--------|---------|----------|------|
| **Full (Recommended)** | `pip install -e ".[full]"` | 100% features | ~13 min |
| Standard | `pip install -e .` | 80% features | ~6 min |

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

- Python 3.10-3.12 (MCP requires Python 3.10+)
- 5-10 GB disk space (for full installation)
- macOS, Linux, or Windows

## Configure Your MCP Client

**Important:** Use your virtual environment's Python path for reliable operation.

### Find Your Virtual Environment Python Path

```bash
# After activating your virtual environment
which python
# This will output something like:
# /Users/apple/Research/SpatialTrans_MCP/chatspatial_env/bin/python
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
      "command": "/Users/apple/Research/SpatialTrans_MCP/chatspatial_env/bin/python",
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

# Add ChatSpatial MCP server
claude mcp add chatspatial /path/to/your/chatspatial_env/bin/python -- -m chatspatial server

# Verify installation
claude mcp list
```

{: .note }
💡 The `--` separates Claude's options from the server command. Use `--scope user` to make ChatSpatial available across all projects.

{: .warning }
⚠️ **Never use system Python:** Always use the Python from your virtual environment to avoid dependency conflicts.

## Verify Installation

```bash
# First, ensure you're in the virtual environment
which python
# Should output path to your virtual environment, not system Python

# Check installation
python -m chatspatial --help

# Test import
python -c "import chatspatial; print('✅ Installation successful')"
```

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

**Next:** [Quick Start Guide](quick-start.html) to run your first analysis