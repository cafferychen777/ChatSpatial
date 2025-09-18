# Installation Guide

Complete installation guide for ChatSpatial MCP server.

## System Requirements

### Minimum Requirements
- **Operating System**: macOS, Linux, or Windows 10+
- **Python**: 3.10 or higher (3.11 recommended)
- **Memory**: 8GB RAM minimum, 16GB recommended
- **Storage**: 5GB free space for installation and data
- **Internet**: Required for package downloads and MCP communication

### Recommended Setup
- **Python**: 3.11 (best compatibility)
- **Memory**: 16GB+ RAM for large datasets
- **Storage**: SSD with 20GB+ free space
- **GPU**: Optional, for accelerated deep learning methods

## Installation Methods

### Method 1: Quick Install (Recommended)

```bash
# Create conda environment
conda create -n chatspatial python=3.11
conda activate chatspatial

# Clone repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# Install ChatSpatial
pip install -e .

# Verify installation
chatspatial --help
```

### Method 2: Development Install

```bash
# Clone with development dependencies
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# Create environment with dev dependencies
conda env create -f environment.yml
conda activate chatspatial

# Install in development mode
pip install -e ".[dev]"

# Run tests to verify
pytest tests/
```

### Method 3: Docker Install

```bash
# Pull Docker image
docker pull cafferychen777/chatspatial:latest

# Run container
docker run -it --rm \
  -v $(pwd)/data:/app/data \
  -p 8000:8000 \
  cafferychen777/chatspatial:latest
```

## Dependencies

### Core Dependencies
- **scanpy**: Single-cell analysis
- **squidpy**: Spatial transcriptomics
- **anndata**: Data structures
- **pandas**: Data manipulation
- **numpy**: Numerical computing
- **scipy**: Scientific computing

### Optional Dependencies
- **torch**: Deep learning methods
- **tensorflow**: Alternative deep learning
- **R**: For scType and other R-based methods
- **CUDA**: GPU acceleration

### Install Optional Dependencies

```bash
# Deep learning support
pip install torch torchvision torchaudio

# R integration
conda install -c conda-forge r-base r-essentials

# GPU support (NVIDIA)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

## Platform-Specific Instructions

### macOS

```bash
# Install Homebrew (if not installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install conda
brew install miniconda

# Follow standard installation
conda create -n chatspatial python=3.11
conda activate chatspatial
pip install -e .
```

### Linux (Ubuntu/Debian)

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install system dependencies
sudo apt install -y python3-dev python3-pip git build-essential

# Install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Follow standard installation
conda create -n chatspatial python=3.11
conda activate chatspatial
pip install -e .
```

### Windows

```powershell
# Install Miniconda from https://docs.conda.io/en/latest/miniconda.html
# Open Anaconda Prompt and run:

conda create -n chatspatial python=3.11
conda activate chatspatial

# Clone repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# Install
pip install -e .
```

## Verification

### Test Installation

```bash
# Activate environment
conda activate chatspatial

# Test CLI
chatspatial --version
chatspatial --help

# Test Python import
python -c "import chatspatial; print('Installation successful!')"

# Test MCP server
python -m chatspatial --test
```

### Run Example Analysis

```bash
# Download test data
python -c "
import chatspatial
from chatspatial.data import download_example_data
download_example_data()
"

# Run basic test
python -c "
from chatspatial.tools.data_management import load_data
result = load_data('data/example_visium.h5ad', 'test')
print(f'Loaded {result.n_obs} spots, {result.n_vars} genes')
"
```

## Troubleshooting

### Common Issues

#### 1. Import Errors
```bash
# Solution: Reinstall with dependencies
pip uninstall chatspatial
pip install -e ".[all]"
```

#### 2. Memory Errors
```bash
# Solution: Increase swap space or use smaller datasets
# For large datasets, consider subsampling:
python -c "
import scanpy as sc
adata = sc.read_h5ad('large_data.h5ad')
adata = adata[::10, :].copy()  # Subsample every 10th spot
adata.write('smaller_data.h5ad')
"
```

#### 3. R Integration Issues
```bash
# Install R packages
conda install -c conda-forge r-base r-essentials
R -e "install.packages(c('Seurat', 'dplyr', 'HGNChelper'))"
```

#### 4. GPU Issues
```bash
# Check CUDA availability
python -c "
import torch
print(f'CUDA available: {torch.cuda.is_available()}')
print(f'CUDA version: {torch.version.cuda}')
"
```

### Getting Help

If you encounter issues:

1. **Check logs**: Look for error messages in terminal output
2. **Update packages**: Ensure all dependencies are up-to-date
3. **GitHub Issues**: Report bugs at [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
4. **Discussions**: Ask questions at [GitHub Discussions](https://github.com/cafferychen777/ChatSpatial/discussions)

## Next Steps

After successful installation:

1. **Configure MCP**: Set up Claude Desktop integration
2. **Download Data**: Get example datasets for testing
3. **Run Tutorial**: Follow the basic spatial analysis tutorial
4. **Explore API**: Check out the complete API reference

See [Quick Start](./quick-start.md) for the next steps!
