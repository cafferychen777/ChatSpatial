# ChatSpatial Installation

## Quick Start (2 minutes)

### Step 1: Install ChatSpatial

```bash
# Clone repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd chatspatial

# Recommended: Install with all features
pip install -e ".[full]"
```

> 💡 **Why [full]?** Enables all 16+ analysis methods including Cell2location, CellRank, SpaGCN, and more. Takes ~13 minutes but worth it for the complete experience.

### Step 2: Configure Claude Desktop

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

### Step 3: Restart Claude Desktop

That's it! Start analyzing your spatial data with natural language.

## Installation Options

| Option | Install Command | Features | Install Time |
|--------|----------------|----------|--------------|
| **Full (Recommended)** | `pip install -e ".[full]"` | 100% features | ~13 minutes |
| Standard | `pip install -e .` | 80% features | ~6 minutes |

### What's included in each option?

**Standard Installation** (Default):
- ✅ Core spatial analysis (Moran's I, Getis-Ord, etc.)
- ✅ Basic deconvolution with scvi-tools
- ✅ RNA velocity with scVelo
- ✅ Cell communication (LIANA, CellPhoneDB)
- ✅ Batch integration (Harmony, BBKNN)
- ✅ Gene enrichment analysis

**Full Installation** (Recommended):
- ✅ Everything in Standard, plus:
- ✅ Deep learning methods (PyTorch)
- ✅ Advanced deconvolution (Cell2location)
- ✅ Advanced trajectory (CellRank, Palantir)
- ✅ Spatial domains (SpaGCN, STAGATE)
- ✅ Spatial variable genes (GASTON, SpatialDE)
- ✅ R-based methods (RCTD, Spotlight)

## Requirements

- Python 3.8-3.12
- macOS, Linux, or Windows
- 5-10 GB disk space (for full installation)

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Import errors | Update pip: `pip install --upgrade pip` |
| Package conflicts | `pip install --force-reinstall -e ".[full]"` |
| Claude doesn't see server | Restart Claude Desktop after configuration |

## Verify Installation

```bash
# Check if ChatSpatial is installed
python -m chatspatial --help

# Test in Python
python -c "import chatspatial; print('✅ Installation successful')"
```

## Getting Help

- **Issues**: [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
- **Documentation**: Check docstrings with `help(function_name)`
- **Community**: Ask questions in Discussions