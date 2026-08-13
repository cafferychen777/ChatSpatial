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
uvx --from 'chatspatial==1.3.5' chatspatial server
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
| Full | `uv pip install 'chatspatial[full]'` | You want the broadest mutually compatible PyPI method set on a workstation. **Requires R ≥ 4.4 on PATH** (see below) |

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
uv pip install 'chatspatial[trajectory]'          # CellRank, Palantir
uv pip install 'chatspatial[deep-learning]'       # scVI, scANVI, VeloVI, DestVI backend
uv pip install 'chatspatial[integration]'         # Harmony, BBKNN, Scanorama
uv pip install 'chatspatial[deconvolution]'       # FlashDeconv, Cell2location
uv pip install 'chatspatial[rctd-python]'         # PyTorch RCTD backend (rctd-py)
uv pip install 'chatspatial[spatial-stats]'       # PySAL/ESDA extensions
uv pip install 'chatspatial[spatial-domains]'     # SpaGCN and BANKSY
```

ChatSpatial tools fail with targeted installation guidance if you call a method
whose optional dependency is not installed.

The `rctd-python` extra is intentionally separate from `deconvolution` and
`full`. Select it with `rctd_backend="python"`; the default remains the
spacexr R backend. On first use, rctd-py downloads an approximately 400 MB
likelihood table into `~/.cache/rctd`, so the first run needs network access
and additional disk space. ChatSpatial downloads this cache with the bundled
certificate authority and publishes it atomically so failed or concurrent
downloads do not leave a partial cache file.

FastCCC and the PyPI release of pyGPCCA currently require incompatible Jinja2
versions. Keep `cell-communication`/`full` and `trajectory` in separate runtime
environments for standard PyPI installations. Repository developers who need
both use the pinned shared-environment constraints documented below.

### Shared repository environment

The workspace environment intentionally combines development, documentation,
and every optional method family. Install it with the repository constraint
set so upstream metadata cannot leave only part of a dependency family upgraded:

```bash
# Replace the stale PyPI metadata even when pyGPCCA 1.0.4 is already installed.
python -m pip install --force-reinstall --no-deps \
  -c constraints/shared-py312.txt pygpcca

python -m pip install -c constraints/shared-py312.txt \
  -e '.[full,trajectory,spatial-domains,dev]'
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

- [Docker / GHCR](docker.md) — run ChatSpatial without local Python dependency resolution
- [Configuration Guide](advanced/configuration.md) — exact client setup
- [Quick Start](quickstart.md) — first successful analysis
- [Troubleshooting](advanced/troubleshooting.md) — fix install or runtime issues
