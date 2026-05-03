# ChatSpatial Installation

This page covers **environment setup and package installation**.

- For exact MCP client configuration, see the [Configuration Guide](docs/advanced/configuration.md).
- For your first analysis after setup, see [Quick Start](docs/quickstart.md).
- For installation failures, see [Troubleshooting](docs/advanced/troubleshooting.md).

---

## Requirements

- **Python 3.11-3.13** (3.12 recommended)
- **8GB+ RAM** (16GB+ for large datasets)
- **macOS, Linux, or Windows**

---


## Docker Installation

If local installation fails or you want the most reproducible setup, use the published Docker image:

```bash
# Pull and verify the image
docker pull ghcr.io/cafferychen777/chatspatial:latest
docker run --rm ghcr.io/cafferychen777/chatspatial:latest --version
docker run --rm ghcr.io/cafferychen777/chatspatial:latest server --help
```

Use this command shape in MCP clients that support Docker-backed stdio servers:

```bash
docker run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:latest server --transport stdio
```

Mounted data are available under `/data`, and generated outputs are written under `/outputs` by default. In prompts, use container paths such as `/data/sample.h5ad`, not the original host path.

See [DOCKER.md](DOCKER.md) for MCP client examples, SSE mode, volume details, and local source builds.

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
| Full | `uv pip install 'chatspatial[full]'` | You want the broadest method coverage on a workstation and accept a large install. **Requires R ≥ 4.4 on PATH** (see [R-based methods](#r-based-methods) below) |

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

For exact client-specific syntax, use the [Configuration Guide](docs/advanced/configuration.md).

---

## Step 4: Verify the Installation

```bash
python -c "import chatspatial; print(f'ChatSpatial {chatspatial.__version__} ready')"
python -m chatspatial server --help
```

If both commands work, continue to [Quick Start](docs/quickstart.md).

---

## Platform Notes

### macOS (Intel / x86_64)

Some dependencies in `chatspatial[full]` do not publish pre-built wheels for Intel Macs:

- **gseapy** — requires the Rust toolchain to compile from source
- **llvmlite** (via numba) — requires LLVM to compile from source

Install the prerequisites first:

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"

# Install LLVM (for llvmlite)
brew install llvm
export LLVM_CONFIG="$(brew --prefix llvm)/bin/llvm-config"

# Then install ChatSpatial with all optional Python methods
uv pip install chatspatial[full]
```

> Apple Silicon Macs (M1/M2/M3/M4) have pre-built wheels for all dependencies and do not require these steps.

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

The `[full]` extra includes `rpy2`, which **requires R to be available on your `PATH` at install time** (it compiles against your R installation). On HPC systems where R is provided via modules, run `module load R` (or equivalent) before installing. If you only need Python-based methods, the standard install (`uv pip install chatspatial`) does not require R.

Once R is available, install the R packages used by ChatSpatial:

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

- [Docker Guide](DOCKER.md) — Run ChatSpatial from GHCR without local Python dependency resolution
- [Configuration Guide](docs/advanced/configuration.md) — Exact client setup
- [Quick Start](docs/quickstart.md) — First successful analysis
- [Troubleshooting](docs/advanced/troubleshooting.md) — Fix install or runtime issues
