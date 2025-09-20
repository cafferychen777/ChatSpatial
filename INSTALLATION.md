# ChatSpatial Installation Guide

ChatSpatial provides AI-powered spatial transcriptomics analysis through a simple, reliable installation process optimized for compatibility and performance.

## Quick Start

```bash
# Core installation (recommended for most users)
pip install -e .

# Advanced features (deep learning, specialized methods)
pip install -e ".[advanced]"

# Development setup
pip install -e ".[dev]"
```

## Installation Options

### Core Installation
**Compatible with Python 3.8-3.12** (3.10+ recommended)

Essential spatial transcriptomics analysis capabilities:
- Data loading and preprocessing (10X Visium, Slide-seq, etc.)
- Clustering and dimensionality reduction (Leiden, UMAP)
- Spatial analysis (Squidpy integration)
- Visualization (matplotlib, seaborn, scanpy)
- MCP server functionality for AI integration

```bash
pip install -e .
```

### Advanced Installation
**Compatible with Python 3.8-3.12** (3.10+ recommended)

Adds cutting-edge spatial analysis methods:
- Deep learning deconvolution (Cell2location, scvi-tools)
- RNA velocity analysis (scVelo, CellRank)
- Advanced integration methods (Harmony, BBKNN)
- Spatial domain identification (SpaGCN, STAGATE)
- Cell communication analysis (LIANA, CellPhoneDB)
- Spatial variable gene detection (GASTON, SpatialDE)

```bash
pip install -e ".[advanced]"
```

### Development Installation
For contributors and advanced users:

```bash
pip install -e ".[dev]"
```

Includes testing and development tools (pytest, black, mypy).

## Python Version Compatibility

| Version | Core | Advanced | Development |
|---------|------|----------|-------------|
| Python 3.8 | ✅ | ✅ | ✅ |
| Python 3.9 | ✅ | ✅ | ✅ |
| Python 3.10 | ✅ | ✅ | ✅ |
| Python 3.11 | ✅ | ✅ | ✅ |
| Python 3.12 | ✅ | ✅ | ✅ |

## Key Dependencies

### Core Dependencies
- **numpy** 2.x (latest stable, significant performance improvements)
- **pandas** 2.x (latest features and performance)
- **scanpy** 1.x (single-cell analysis)
- **squidpy** 1.x (spatial analysis)
- **fastapi** (MCP server framework)

### Advanced Dependencies
- **torch** 2.x (deep learning framework)
- **scvi-tools** (probabilistic models)
- **cellrank** (trajectory analysis)
- **liana** (cell communication)

## Installation Examples

### Basic Setup
```bash
# Clone repository
git clone https://github.com/your-org/chatspatial.git
cd chatspatial

# Install core features
pip install -e .

# Test installation
python -c "import chatspatial; print('✅ Installation successful')"
```

### GPU-Accelerated Setup
```bash
# Install PyTorch with CUDA support first
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# Install ChatSpatial with advanced features
pip install -e ".[advanced]"
```

### Conda Environment
```bash
# Create clean environment
conda create -n chatspatial python=3.11
conda activate chatspatial

# Install ChatSpatial
pip install -e ".[advanced]"
```

## R Integration (Optional)

For methods requiring R (e.g., RCTD deconvolution):

```bash
# Install R (if needed)
# macOS: brew install r
# Ubuntu: sudo apt install r-base r-base-dev

# Install R interface
pip install rpy2>=3.4.0
```

## Testing Installation

```bash
# Test core functionality
python -c "import chatspatial; print('Core installation OK')"

# Test MCP server
python -m chatspatial --help

# Test advanced features (if installed)
python -c "import scvi; print('Advanced features OK')"
```

## Common Installation Tips

### Virtual Environment (Recommended)
```bash
python -m venv chatspatial_env
source chatspatial_env/bin/activate  # Linux/macOS
# or: chatspatial_env\Scripts\activate  # Windows
pip install -e ".[advanced]"
```

### Upgrading Dependencies
```bash
pip install --upgrade -e ".[advanced]"
```

### Clean Installation
```bash
pip cache purge
pip install -e ".[advanced]" --force-reinstall
```

## Docker Installation

```dockerfile
FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential git && \
    rm -rf /var/lib/apt/lists/*

# Install ChatSpatial
COPY . /app/chatspatial
WORKDIR /app/chatspatial
RUN pip install -e ".[advanced]"

# Verify installation
RUN python -c "import chatspatial; print('Docker installation successful')"
```

## Troubleshooting

### Installation Fails
1. Update pip: `pip install --upgrade pip`
2. Check Python version: `python --version`
3. Try clean installation: `pip cache purge`

### Import Errors
1. Verify installation: `pip list | grep chatspatial`
2. Check Python path: `python -c "import sys; print(sys.path)"`
3. Reinstall: `pip install --force-reinstall -e .`

### Performance Issues
1. Ensure numpy 2.x: `python -c "import numpy; print(numpy.__version__)"`
2. For GPU: Verify CUDA installation with `nvidia-smi`
3. For CPU: Consider installing Intel MKL: `pip install mkl`

## Getting Help

1. **Check installation**: `python -c "import chatspatial; print('OK')"`
2. **Verify dependencies**: `pip list`
3. **GitHub Issues**: Report installation problems with system details
4. **Documentation**: Check function-specific requirements in docstrings

The streamlined installation system ensures reliable dependency resolution and optimal performance across different environments and Python versions.