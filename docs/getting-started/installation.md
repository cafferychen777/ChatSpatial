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

```bash
# Clone repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd chatspatial

# Recommended: Install with all features
pip install -e ".[full]"
```

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

Add to `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "python",
      "args": ["-m", "chatspatial"]
    }
  }
}
```

## Verify Installation

```bash
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
| Claude doesn't see server | Restart Claude Desktop |

## Getting Help

- **GitHub Issues**: [Report problems](https://github.com/cafferychen777/ChatSpatial/issues)
- **Documentation**: Check docstrings with `help(function_name)`

---

**Next:** [Quick Start Guide](quick-start.html) to run your first analysis