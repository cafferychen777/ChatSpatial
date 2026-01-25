# ChatSpatial Installation

## Step 1: Create Virtual Environment

```bash
# Python 3.11+ required
python3 -m venv venv
source venv/bin/activate  # macOS/Linux
# venv\Scripts\activate   # Windows

# Or use conda
conda create -n chatspatial python=3.12
conda activate chatspatial
```

## Step 2: Install

```bash
pip install --upgrade pip
pip install chatspatial[full]
```

| Option | Command | Features | Time |
|--------|---------|----------|------|
| **Full** | `pip install chatspatial[full]` | All 60+ methods | ~13 min |
| Standard | `pip install chatspatial` | Core methods | ~6 min |

## Step 3: Configure MCP Client

```bash
# Get your Python path
which python  # e.g., /Users/you/venv/bin/python

# Add to Claude Code (recommended)
claude mcp add chatspatial /path/to/venv/bin/python -- -m chatspatial server
```

**Other clients:** See [Configuration Guide](docs/advanced/configuration.md) for Codex, Claude Desktop, and advanced options.

## Step 4: Verify

```bash
python -c "import chatspatial; print('Ready')"
```

Restart your client, then see [Quick Start](docs/quickstart.md) to begin analyzing.

---

## Requirements

- **Python 3.11+** (3.12 recommended)
- **8GB+ RAM** (16GB+ for large datasets)
- **macOS, Linux, or Windows**

---

## Platform Notes

### Windows

**Not available:** SingleR, PETSc

**Use instead:** Tangram, scANVI, CellAssign for annotation; CellRank works without PETSc.

### Python Version Error

If you see `mcp>=1.17.0 not found`:

```bash
rm -rf venv
python3.12 -m venv venv
source venv/bin/activate
pip install chatspatial[full]
```

---

## Optional Dependencies

### R Methods

For RCTD, SPOTlight, CARD, scType, Numbat:

```bash
# Install R 4.4+ from https://cran.r-project.org/
Rscript install_r_dependencies.R
```

### STAGATE

```bash
git clone https://github.com/QIFEIDKN/STAGATE_pyG.git
cd STAGATE_pyG && python setup.py install
```

---

## Help

- [Configuration Guide](docs/advanced/configuration.md) — MCP client setup
- [Troubleshooting](docs/advanced/troubleshooting.md) — Common issues
- [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
