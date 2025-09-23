# ChatSpatial Installation

## Quick Start (2 minutes)

### Step 1: Create Virtual Environment (Strongly Recommended)

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

# Option C: Using uv (fast, modern)
uv venv
source .venv/bin/activate
```

### Step 2: Install ChatSpatial

```bash
# Recommended: Install with all features
pip install -e ".[full]"
```

> 💡 **Why [full]?** Enables all 16+ analysis methods including Cell2location, CellRank, SpaGCN, and more. Takes ~13 minutes but worth it for the complete experience.

### Step 3: Configure Claude Desktop with Virtual Environment

**Important:** Use your virtual environment's Python path in the configuration.

Add to `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/your/chatspatial/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial"]
    }
  }
}
```

**Example with actual path:**
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

> ⚠️ **Critical:** Replace `/path/to/your/chatspatial/` with your actual installation directory!

### Step 4: Restart Claude Desktop

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
| Claude doesn't see server | Check that command path points to virtual environment Python |
| "python not found" error | Use absolute path to virtual environment Python in config |
| Virtual environment issues | Recreate environment: `rm -rf chatspatial_env && python3 -m venv chatspatial_env` |

## Verify Installation

```bash
# Make sure you're in the virtual environment
which python  # Should show path to your virtual environment

# Check if ChatSpatial is installed
python -m chatspatial --help

# Test in Python
python -c "import chatspatial; print('✅ Installation successful')"
```

## Getting Help

- **Issues**: [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
- **Documentation**: Check docstrings with `help(function_name)`
- **Community**: Ask questions in Discussions