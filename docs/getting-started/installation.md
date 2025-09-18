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

ChatSpatial uses a tiered dependency system designed for maximum compatibility and user control. This guide explains how to install ChatSpatial with the right level of dependencies for your needs.

## Quick Start

```bash
# Minimal installation (recommended for most users)
pip install -e .

# Advanced features (deep learning, deconvolution)  
pip install -e ".[advanced]"

# All features (experimental, requires R)
pip install -e ".[experimental]"
```

## Installation Tiers

### Tier 1: Minimal (Core Features)

{: .highlight }
**Compatible with Python 3.8-3.12**

This provides essential spatial transcriptomics analysis capabilities:
- Data loading and preprocessing
- Basic visualization (matplotlib, seaborn)
- Clustering (Leiden, k-means)
- Dimensionality reduction (PCA, UMAP)
- Spatial analysis (Squidpy integration)
- MCP server functionality

```bash
pip install -e .
# OR 
pip install -r requirements-minimal.txt
```

**What you get:**
- ✅ Load spatial data (10X Visium, Slide-seq, etc.)
- ✅ Quality control and filtering
- ✅ Clustering and UMAP embedding
- ✅ Basic spatial statistics
- ✅ Standard visualizations
- ✅ MCP server for AI integration

### Tier 2: Advanced Features

{: .highlight }
**Compatible with Python 3.8-3.11** (deep learning limitations)

Adds cutting-edge spatial analysis methods:
- Deep learning deconvolution (Cell2location)
- Advanced integration methods
- RNA velocity analysis
- Enhanced visualizations

```bash
pip install -e ".[advanced]"
# OR
pip install -r requirements-advanced.txt
```

**Additional features:**
- ✅ Cell2location spatial deconvolution
- ✅ scvi-tools ecosystem (scANVI, DestVI)
- ✅ RNA velocity with scVelo
- ✅ Advanced trajectory analysis (CellRank)
- ✅ Batch integration (Harmony, BBKNN, Scanorama)
- ✅ Spatial variable gene detection
- ✅ Interactive visualizations (Plotly)

### Tier 3: Experimental

{: .warning }
**Compatible with Python 3.8-3.11** (R interface limitations)

Cutting-edge and experimental features:
- R-based methods (RCTD deconvolution)
- High-performance computing
- Emerging spatial methods

```bash
pip install -e ".[experimental]"
# OR  
pip install -r requirements-experimental.txt
```

**Experimental features:**
- ✅ RCTD deconvolution (requires R)
- ✅ High-performance linear algebra (PETSc/SLEPc)
- ✅ Emerging spatial frameworks
- ⚠️ May have compatibility issues

## Python Version Compatibility

| Feature | Python 3.8 | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 |
|---------|-------------|-------------|--------------|--------------|--------------|
| **Minimal** | ✅ | ✅ | ✅ | ✅ | ✅ |
| **Advanced** | ✅ | ✅ | ✅ | ✅ | ⚠️¹ |
| **Experimental** | ✅ | ✅ | ✅ | ⚠️² | ❌³ |

¹ Some deep learning dependencies may lag Python 3.12 support  
² R interface (rpy2) has limited Python 3.11 support  
³ Most experimental features not yet compatible with Python 3.12

## Common Installation Issues

### Issue: TensorFlow 1.x Conflicts

{: .note }
**Problem:** Some legacy tools require TensorFlow 1.x, which is incompatible with modern Python.

**Solution:** ChatSpatial no longer supports TensorFlow 1.x. Use PyTorch-based alternatives:
- Instead of STAGATE → Use SpaGCN
- Instead of old scVI → Use scvi-tools 1.0+

### Issue: PyTorch Version Conflicts

{: .note }
**Problem:** Different tools require different PyTorch versions.

**Solution:** ChatSpatial standardizes on PyTorch 2.x:

```bash
pip install torch>=2.0.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

For GPU support:
```bash
pip install torch>=2.0.0 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

### Issue: R Dependencies

{: .note }
**Problem:** RCTD requires R installation and rpy2.

**Solution:** Install R first, then Python interface:

```bash
# macOS
brew install r

# Ubuntu/Debian  
sudo apt install r-base r-base-dev

# Then install Python interface
pip install rpy2>=3.4.0
```

## Testing Your Installation

```bash
# Test dependency detection
python test_dependencies.py

# Test specific installation tier
python test_dependencies.py --level minimal
python test_dependencies.py --level advanced  
python test_dependencies.py --level experimental

# Verbose output for debugging
python test_dependencies.py --verbose
```

## Docker Installation

For reproducible environments:

```dockerfile
FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install ChatSpatial
COPY . /app/chatspatial
WORKDIR /app/chatspatial
RUN pip install -e ".[advanced]"

# Test installation
RUN python test_dependencies.py --level advanced
```

## Troubleshooting

### Get Dependency Report

```python
from chatspatial.utils.dependency_manager import dependency_manager
dependency_manager.print_dependency_report()
```

### Check Available Methods

```python  
from chatspatial.utils.dependency_manager import get_available_methods

# What deconvolution methods can I use?
deconv_deps = {
    "cell2location": ["cell2location", "torch"],
    "rctd": ["rpy2"], 
    "destvi": ["scvi-tools", "torch"]
}
print("Available:", get_available_methods(deconv_deps))
```

### Force Reinstall Dependencies

```bash
pip install --force-reinstall --no-cache-dir -e ".[advanced]"
```

## Getting Help

1. **Check installation:** `python test_dependencies.py --verbose`
2. **View dependency report:** Run the dependency manager as shown above
3. **Open an issue:** Include the output of the dependency test
4. **Join discussions:** ChatSpatial community forums

{: .highlight }
The tiered installation system ensures you get maximum functionality while maintaining compatibility with your specific environment and Python version.

---

**Next:** [Quick Start Guide](quick-start.html) to run your first analysis