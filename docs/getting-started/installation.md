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
💡 **Why virtual environment?** Isolates dependencies, prevents conflicts, and makes configuration with Claude Desktop cleaner.

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
- ✅ Advanced deconvolution (Cell2location)
- ✅ Advanced trajectory (CellRank, Palantir)
- ✅ Spatial domains (SpaGCN, STAGATE)
- ✅ Spatial variable genes (GASTON, SpatialDE)

## Requirements

- Python 3.8-3.12
- 5-10 GB disk space (for full installation)
- macOS, Linux, or Windows

## Configure with Claude Desktop

**Important:** Use your virtual environment's Python path for reliable operation.

### Find Your Virtual Environment Python Path

```bash
# After activating your virtual environment
which python
# This will output something like:
# /Users/apple/Research/SpatialTrans_MCP/chatspatial_env/bin/python
```

### Add to Claude Desktop Configuration

Edit `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/your/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial"]
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
      "args": ["-m", "chatspatial"]
    }
  }
}
```

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
| Claude doesn't see server | Check that command path points to virtual environment Python |
| "command not found" error | Use absolute path to Python in virtual environment |
| Wrong Python version | Recreate virtual environment with correct Python version |

## Getting Help

- **GitHub Issues**: [Report problems](https://github.com/cafferychen777/ChatSpatial/issues)
- **Documentation**: Check docstrings with `help(function_name)`

---

**Next:** [Quick Start Guide](quick-start.html) to run your first analysis