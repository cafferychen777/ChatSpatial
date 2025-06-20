# Spotiphy Installation Guide

## Overview

Spotiphy is a spatial transcriptomics deconvolution method that uses Bayesian inference for cell type proportion estimation. Due to its complex dependencies, it requires special installation steps.

**Note**: Spotiphy is available on PyPI (https://pypi.org/project/spotiphy/) but has strict dependency requirements that may conflict with newer Python versions.

## Python Version Requirements

**IMPORTANT**: Spotiphy requires Python 3.8-3.11. It is NOT compatible with Python 3.12+ due to TensorFlow 2.12.0 dependency.

## Installation Steps

### 1. Check Python Version

```bash
python --version
```

Ensure you have Python 3.8, 3.9, 3.10, or 3.11.

### 2. Create a Virtual Environment (Recommended)

```bash
# Using conda
conda create -n chatspatial python=3.10
conda activate chatspatial

# Or using venv
python3.10 -m venv chatspatial_env
source chatspatial_env/bin/activate  # On Windows: chatspatial_env\Scripts\activate
```

### 3. Install ChatSpatial Base Dependencies

```bash
cd /path/to/chatspatial
pip install -e .
```

### 4. Install Spotiphy Dependencies

First, install the core dependencies:

```bash
pip install torch pyro-ppl opencv-python tqdm
```

### 5. Install Spotiphy

You have two options:

#### Option A: Install from PyPI (Recommended)
```bash
pip install spotiphy
```

#### Option B: Install from Local Source
```bash
# From the ChatSpatial directory
pip install -e third_party/Spotiphy/
```

Note: Installation may take several minutes as it needs to compile scikit-learn 1.2.2 and install TensorFlow. The PyPI version is identical to our local copy.

### 6. Verify Installation

```python
from chatspatial.tools.deconvolution import is_spotiphy_available

is_available, error_msg = is_spotiphy_available()
if is_available:
    print("✓ Spotiphy is successfully installed!")
else:
    print(f"✗ Spotiphy installation issue: {error_msg}")
```

## Troubleshooting

### Common Issues

1. **Python Version Error**: 
   - Error: "No matching distribution found for tensorflow"
   - Solution: Use Python 3.8-3.11

2. **scikit-learn Compilation Error**:
   - Error: "error: Microsoft Visual C++ 14.0 or greater is required" (Windows)
   - Solution: Install Visual Studio Build Tools
   
   - Error: "clang: error: unsupported option '-fopenmp'" (macOS)
   - Solution: Install libomp: `brew install libomp`

3. **Memory Error During Installation**:
   - Solution: Close other applications and retry, or use a machine with more RAM

4. **GPU Support**:
   - For CUDA support: Install PyTorch with CUDA
   - For Apple Silicon: PyTorch should automatically use MPS

### Alternative Installation (Pre-built Wheels)

If compilation fails, try installing pre-built wheels:

```bash
# Install dependencies with specific versions
pip install numpy==1.22.4
pip install scipy==1.9.1
pip install scikit-learn==1.2.2
pip install pandas==2.0.0
pip install scanpy==1.9.3
pip install squidpy==1.2.3
pip install tensorflow==2.12.0
pip install stardist==0.8.3

# Then install Spotiphy
pip install -e third_party/Spotiphy/ --no-deps
```

## Using Spotiphy in ChatSpatial

Once installed, you can use Spotiphy through ChatSpatial's deconvolution tool:

```python
# In your ChatSpatial prompts
"Please deconvolve my spatial data using Spotiphy method with the reference dataset"
```

Or programmatically:

```python
from chatspatial.tools.deconvolution import deconvolve_spatial_data
from chatspatial.models.data import DeconvolutionParameters

params = DeconvolutionParameters(
    method="spotiphy",
    reference_data_id="my_reference",
    cell_type_key="cell_type",
    n_epochs=8000,
    use_gpu=True  # If available
)

result = await deconvolve_spatial_data(
    data_id="my_spatial_data",
    data_store=data_store,
    params=params
)
```

## System Requirements

- **RAM**: Minimum 16GB, 32GB recommended for large datasets
- **Storage**: ~5GB for dependencies
- **GPU**: Optional but recommended for faster computation
  - NVIDIA GPU with CUDA support
  - Apple Silicon (M1/M2) with MPS support

## Additional Notes

- Spotiphy uses probabilistic inference which can be computationally intensive
- Runtime scales with dataset size and number of epochs
- Consider using fewer epochs (e.g., 1000-2000) for initial testing