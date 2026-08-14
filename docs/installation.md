# Installation

Use `uvx` for the shortest setup: it runs ChatSpatial in an isolated,
automatically managed environment. Use a persistent environment when you need
optional method families or tighter dependency control. If you want a
containerized runtime, use the [Docker / GHCR](docker.md) guide instead.

- For exact MCP client syntax, see [Configuration Guide](advanced/configuration.md).
- For your first workflow after setup, see [Quick Start](quickstart.md).
- For installation failures, see [Troubleshooting](advanced/troubleshooting.md).

---

## Requirements

- **uv** for the recommended zero-environment setup
- **Python 3.11-3.14** (3.12 recommended) for persistent environments
- **MCP Python SDK 2.x** (installed automatically with ChatSpatial)
- **8GB+ RAM** (16GB+ for large datasets)
- **macOS, Linux, or Windows**
- **Docker** only if you choose the container runtime

---

## Choose a Runtime

| Runtime | Use when | Guide |
|---------|----------|-------|
| **uvx (recommended)** | You want the shortest setup with an isolated, cached environment | Continue below |
| **Persistent Python environment** | You need optional methods or direct control of packages | [Persistent installation](#persistent-python-installation) |
| **Docker / GHCR** | You want the most reproducible runtime or local dependency resolution fails | [Docker / GHCR](docker.md) |

---

## Recommended: Run with uvx

### Step 1: Install uv

```bash
# macOS and Linux
curl -LsSf https://astral.sh/uv/install.sh | sh
```

On Windows, use the installer documented by the
[uv project](https://docs.astral.sh/uv/getting-started/installation/).

### Step 2: Add ChatSpatial to your MCP client

Codex:

```bash
codex mcp add chatspatial -- uvx --from chatspatial chatspatial server
```

Claude Code:

```bash
claude mcp add --scope user chatspatial -- \
  uvx --from chatspatial chatspatial server
```

For other clients, use `uvx` as the command and the following arguments:

```text
--from chatspatial chatspatial server
```

The first launch downloads the core scientific stack into an isolated cache.
Later launches reuse that cache. No virtual environment or Python executable
path needs to be managed manually.

### Step 3: Verify

```bash
uvx --from chatspatial chatspatial --version
```

Restart the MCP client, then confirm that ChatSpatial exposes its tools. In
Codex, use `/mcp`; in Claude Code, run `claude mcp list`.

### Pin a release

Use an exact version when a reproducible runtime matters:

```bash
uvx --from 'chatspatial==1.3.8' chatspatial server
```

Without a pin, `uvx` resolves the current PyPI release and reuses its cached
environment. Pin the package for a manuscript or long-running analysis.

### Run optional method families with uvx

Extras can be selected in the `--from` package specification:

```bash
uvx --from 'chatspatial[cell-communication,velocity]' chatspatial server
```

Because MCP clients store a single command, update that command when you change
extras. For long-lived, heavily customized stacks, use a persistent environment
instead.

---

(persistent-python-installation)=
## Persistent Python Installation

### Step 1: Create an environment

```bash
# venv
python3.12 -m venv venv
source venv/bin/activate  # macOS/Linux
# venv\Scripts\activate   # Windows

# or conda
conda create -n chatspatial python=3.12
conda activate chatspatial
```

### Step 2: Install ChatSpatial

```bash
uv pip install chatspatial
```

ChatSpatial depends on a large scientific Python stack. `uv` generally resolves
it faster and more reliably than `pip`.

### Install options

| Option | Command | Use when |
|--------|---------|----------|
| **Standard** | `uv pip install chatspatial` | You want the MCP server, data loading, preprocessing, embeddings, visualization, and core analysis |
| Method extras | `uv pip install 'chatspatial[cell-communication,velocity]'` | You need specific advanced method families |
| Full | `uv pip install 'chatspatial[full]'` | You want the broadest mutually compatible Python method set on a workstation |

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
uv pip install 'chatspatial[cell-communication]'  # LIANA+ and CellPhoneDB
uv pip install 'chatspatial[fastccc]'             # FastCCC (Python 3.11-3.12)
uv pip install 'chatspatial[velocity]'            # scVelo
uv pip install 'chatspatial[trajectory]'          # CellRank (3.12+), Palantir
uv pip install 'chatspatial[deep-learning]'       # scVI, scANVI, VeloVI, DestVI backend
uv pip install 'chatspatial[integration]'         # Harmony, BBKNN, Scanorama
uv pip install 'chatspatial[deconvolution]'       # FlashDeconv, Cell2location
uv pip install 'chatspatial[annotation]'          # Tangram, SingleR, mLLMCellType
uv pip install 'chatspatial[enrichment]'          # GSEA and enrichment maps
uv pip install 'chatspatial[cnv]'                 # infercnvpy
uv pip install 'chatspatial[differential]'        # PyDESeq2
uv pip install 'chatspatial[registration]'        # PASTE and STalign
uv pip install 'chatspatial[spatial-genes]'       # SpatialDE
uv pip install 'chatspatial[rctd-python]'         # PyTorch RCTD backend (rctd-py)
uv pip install 'chatspatial[r-backends]'          # Python bridges for R-based methods
uv pip install 'chatspatial[spatial-stats]'       # PySAL/ESDA extensions
uv pip install 'chatspatial[spatial-domains]'     # GraphST, STAGATE, SpaGCN, BANKSY
```

CellRank is installed on Python 3.12 and newer; Python 3.11 receives Palantir
without the obsolete CellRank 2.0 compatibility patch. LIANA and BANKSY
currently support Python through 3.13. Their extras remain
installable on Python 3.14, but those individual backends are omitted until
their upstream packages add 3.14 support. SpaGCN supports Python 3.14 through
the maintained `spagcn-modern` distribution. GraphST, STAGATE, SpatialDE,
PASTE, and STalign are installed from focused maintained PyPI distributions;
no Git URL, local wheel, or source build command is required.

ChatSpatial tools fail with targeted installation guidance if you call a method
whose optional dependency is not installed.

The `rctd-python` extra is intentionally separate from `deconvolution` and
`full`. Select it with `rctd_backend="python"`; the default remains the
spacexr R backend. On first use, rctd-py downloads an approximately 400 MB
likelihood table into `~/.cache/rctd`, so the first run needs network access
and additional disk space. ChatSpatial downloads this cache with the bundled
certificate authority and publishes it atomically so failed or concurrent
downloads do not leave a partial cache file.

FastCCC is isolated because it supports Python 3.11-3.12 and requires
Jinja2 >= 3.1.6, while the PyPI release of pyGPCCA pulled by CellRank pins
Jinja2 3.0.3. Install `fastccc` and `trajectory` in separate environments.
The `cell-communication`, `full`, and `trajectory` extras can otherwise be
combined normally.

### Shared repository environment

The workspace environment combines development and the mutually compatible
optional method families directly from the project metadata:

```bash
python -m pip install \
  -e '.[full,dev]'
```

---

### Step 3: Connect the environment to an MCP client

After installation, register the environment's Python executable in your MCP
client. The command shape is:

```text
/absolute/path/to/python -m chatspatial server
```

Use the [Configuration Guide](advanced/configuration.md) for exact client syntax,
absolute-path rules, Docker-backed client examples, and the runtime path model.

---

### Step 4: Verify the installation

```bash
python -c "import chatspatial; print(f'ChatSpatial {chatspatial.__version__} ready')"
python -m chatspatial server --help
```

If both commands work, continue to [Quick Start](quickstart.md).

---

## Platform Notes

### macOS (Intel / x86_64)

Some dependencies in `chatspatial[full]` do not publish pre-built wheels for
Intel Macs:

- **gseapy** requires the Rust toolchain to compile from source
- **llvmlite** (via numba) requires LLVM to compile from source

Install those prerequisites before the full optional stack:

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"

# Install LLVM for llvmlite
brew install llvm
export LLVM_CONFIG="$(brew --prefix llvm)/bin/llvm-config"

# Then install ChatSpatial with all optional Python methods
uv pip install 'chatspatial[full]'
```

Apple Silicon Macs (M1/M2/M3/M4) have pre-built wheels for all dependencies and
do not require these steps.

### Windows

**Not available:** SingleR, PETSc

**Use instead:** Tangram, scANVI, CellAssign for annotation; CellRank works without PETSc.

### If Python or MCP dependencies fail to resolve

```bash
rm -rf venv
python3.12 -m venv venv
source venv/bin/activate
uv pip install 'chatspatial[full]'
```

---

## Optional Dependencies

### R-based methods

The `[r-backends]` extra includes `rpy2`, which **requires R 4.5 or newer to be
available on your `PATH` at install time** because it links against that R installation.
R bridges are deliberately excluded from `[full]`: installing the Python bridge
does not install the R implementations used by RCTD, SPOTlight, CellChat,
SCTransform, or other R-backed methods. On HPC systems where R is provided via
modules, run `module load R` (or equivalent) first.

```bash
uv pip install 'chatspatial[r-backends]'
```

Once R is available, install the R packages used by ChatSpatial:

```bash
# Install R 4.5+
Rscript install_r_dependencies.R
```

## Next Steps

- [Docker / GHCR](docker.md) — run ChatSpatial without local Python dependency resolution
- [Configuration Guide](advanced/configuration.md) — exact client setup
- [Quick Start](quickstart.md) — first successful analysis
- [Troubleshooting](advanced/troubleshooting.md) — fix install or runtime issues
