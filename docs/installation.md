# Installation

This page is the docs-site version of the installation guide. It is the canonical installation page for the documentation site.

- For exact MCP client syntax, see [Configuration Guide](advanced/configuration.md).
- For your first workflow after setup, see [Quick Start](quickstart.md).
- For installation failures, see [Troubleshooting](advanced/troubleshooting.md).

---

## Requirements

- **Python 3.11-3.13** (3.12 recommended)
- **8GB+ RAM** (16GB+ for large datasets)
- **macOS, Linux, or Windows**

---

## Step 1: Create an Environment

```bash
# venv
python3 -m venv venv
source venv/bin/activate  # macOS/Linux
# venv\Scripts\activate   # Windows

# or conda
conda create -n chatspatial python=3.12
conda activate chatspatial
```

---

## Step 2: Install ChatSpatial

**Recommended: use `uv` for dependency resolution**

```bash
# Install uv if needed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install the MCP server and core analysis stack
uv pip install chatspatial
```

> **Why `uv`?** ChatSpatial depends on a large scientific Python stack. Standard `pip` can fail on deep dependency resolution; `uv` is more reliable for this environment.

### Install options

| Option | Command | Use when |
|--------|---------|----------|
| **Standard** | `uv pip install chatspatial` | You want the MCP server, data loading, preprocessing, embeddings, visualization, and core analysis |
| Method extras | `uv pip install 'chatspatial[cell-communication,velocity]'` | You need specific advanced method families |
| Full | `uv pip install 'chatspatial[full]'` | You want the broadest method coverage on a workstation and accept a large install |

<details>
<summary>Alternative: pip</summary>

```bash
pip install --upgrade pip
pip install chatspatial
```

If you hit `resolution-too-deep`, switch to `uv`.
</details>

### Optional method families

Install only the method families you plan to use:

```bash
uv pip install 'chatspatial[cell-communication]'  # LIANA+, CellPhoneDB, FastCCC
uv pip install 'chatspatial[velocity]'            # scVelo
uv pip install 'chatspatial[deep-learning]'       # scVI, scANVI, VeloVI, DestVI backend
uv pip install 'chatspatial[integration]'         # Harmony, BBKNN, Scanorama
uv pip install 'chatspatial[deconvolution]'       # FlashDeconv, Cell2location
uv pip install 'chatspatial[spatial-stats]'       # PySAL/ESDA extensions
```

ChatSpatial tools fail with targeted installation guidance if you call a method
whose optional dependency is not installed.

---

## Step 3: Register ChatSpatial in Your MCP Client

1. Activate your environment.
2. Get the **absolute** Python path:

```bash
which python
```

3. Register ChatSpatial using this command shape:

```text
/absolute/path/to/python -m chatspatial server
```

4. Restart your client after configuration changes.

For exact client-specific syntax, use the [Configuration Guide](advanced/configuration.md).

---

## Step 4: Verify the Installation

```bash
python -c "import chatspatial; print(f'ChatSpatial {chatspatial.__version__} ready')"
python -m chatspatial server --help
```

If both commands work, continue to [Quick Start](quickstart.md).

---

## Platform Notes

### Windows

**Not available:** SingleR, PETSc

**Use instead:** Tangram, scANVI, CellAssign for annotation; CellRank works without PETSc.

### If Python or MCP dependencies fail to resolve

```bash
rm -rf venv
python3.12 -m venv venv
source venv/bin/activate
uv pip install chatspatial[full]
```

---

## Optional Dependencies

### R-based methods

For RCTD, SPOTlight, CARD, CellChat, SPARK-X, scType, Numbat, and SCTransform:

```bash
# Install R 4.4+
Rscript install_r_dependencies.R
```

### STAGATE

```bash
git clone https://github.com/QIFEIDKN/STAGATE_pyG.git
cd STAGATE_pyG && python setup.py install
```

---

## Next Steps

- [Configuration Guide](advanced/configuration.md) — exact client setup
- [Quick Start](quickstart.md) — first successful analysis
- [Troubleshooting](advanced/troubleshooting.md) — fix install or runtime issues
