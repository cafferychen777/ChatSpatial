
title: Troubleshooting
description: Common issues and solutions for ChatSpatial
---

# ChatSpatial Troubleshooting Guide

Welcome to the ChatSpatial troubleshooting guide! This comprehensive resource will help you resolve common issues and get back to analyzing spatial transcriptomics data quickly.

## Quick Problem Finder

**Having issues right now?** Jump to the relevant section:

- 🚀 [Installation Problems](#installation-and-dependency-issues)
- 💾 [Memory and Performance](#memory-and-performance-problems)  
- 🔬 [Analysis Not Working](#analysis-specific-issues)
- 🔗 [Cell Communication Errors](#cell-communication-analysis)
- 📁 [Data Loading Issues](#data-format-and-loading-problems)
- ⚠️ [Error Messages](#common-error-messages-with-solutions)
- 💬 [Agent Conversation Issues](#conversation-troubleshooting)

## Installation and Dependency Issues

### Python Environment Problems

#### Problem: ImportError or ModuleNotFoundError

**Common error messages:**
```bash
ImportError: No module named 'chatspatial'
ModuleNotFoundError: No module named 'scanpy'
ImportError: cannot import name 'X' from 'Y'
```

**Step-by-step solution:**
1. **Check your environment:**
   ```bash
   which python
   python --version  # Should be 3.10+ (required for MCP)
   ```

2. **Verify ChatSpatial installation:**
   ```bash
   pip list | grep chatspatial
   # OR
   python -c "import chatspatial; print('SUCCESS')"
   ```

3. **Check dependency status:**
   ```bash
   pip list | grep -E "(numpy|pandas|scanpy|squidpy)"
   # Or check if advanced features are available
   python -c "import scvi; print('Advanced features OK')" 2>/dev/null || echo "Advanced features not installed"
   ```

4. **Reinstall if needed:**
   ```bash
   pip uninstall chatspatial
   pip cache purge
   pip install -e ".[advanced]"  # For full features
   ```

> **💡 Prevention Tips:**
> - Always use conda environments: `conda create -n chatspatial python=3.10`
> - Install with dependency groups: `pip install -e ".[advanced]"` for most users
> - Regularly update: `pip install --upgrade -e ".[advanced]"`

#### Problem: Conda vs Pip Conflicts
**Common symptoms:**
- Package versions don't match requirements
- Mix of conda and pip installations causing conflicts
- `pip check` shows dependency conflicts

**Solution:**
1. **Create clean environment:**
   ```bash
   conda create -n chatspatial-clean python=3.10
   conda activate chatspatial-clean
   ```

2. **Install core scientific packages with conda:**
   ```bash
   conda install numpy pandas scipy scikit-learn matplotlib seaborn
   conda install -c conda-forge scanpy anndata
   ```

3. **Install ChatSpatial with pip:**
   ```bash
   pip install -e ".[advanced]"
   ```

> **💡 Prevention Tip:** Stick to one package manager per environment when possible.

### R Dependencies (for Advanced Features)

#### Problem: rpy2 Installation Fails

**Error messages:**
```bash
ERROR: Failed building wheel for rpy2
fatal error: 'R.h' file not found
```

**Solution:**
1. **Install R first:**
   ```bash
   # macOS
   brew install r
   
   # Ubuntu/Debian
   sudo apt install r-base r-base-dev
   
   # Windows (use WSL or R installer)
   ```

2. **Set R environment variables:**
   ```bash
   export R_HOME=$(R RHOME)
   export R_USER=$(whoami)
   ```

3. **Install rpy2:**
   ```bash
   pip install rpy2>=3.4.0
   ```

> **💡 Prevention Tip:** Only install R dependencies if you need RCTD deconvolution.

### macOS Library Loading Issues

#### Problem: llvmlite/numba libc++.1.dylib Not Found (Apple Silicon)

**Error messages:**
```bash
OSError: dlopen(libllvmlite.dylib): Library not loaded: @rpath/libc++.1.dylib
Referenced from: .../llvmlite/binding/libllvmlite.dylib
Reason: tried: '/opt/homebrew/lib/libc++.1.dylib' (no such file)
```

**Common symptoms:**
- Squidpy neighborhood analysis fails with multiprocessing error
- Numba works in main process but fails in subprocesses
- Any parallel computation using numba crashes

**Solution 1: Fix Library Links (Recommended for quick fix)**
```bash
# Find the exact llvmlite location
find ~/your_env -name "libllvmlite.dylib"

# Fix the library references (requires sudo)
sudo install_name_tool -change "@rpath/libc++.1.dylib" "/usr/lib/libc++.1.dylib" /path/to/llvmlite/binding/libllvmlite.dylib
sudo install_name_tool -change "@rpath/libz.1.dylib" "/usr/lib/libz.1.dylib" /path/to/llvmlite/binding/libllvmlite.dylib

# Verify the fix
otool -L /path/to/llvmlite/binding/libllvmlite.dylib | head -5
# Should show /usr/lib/ paths instead of @rpath
```

**Solution 2: Use Conda/Miniforge (Recommended for long-term)**
```bash
# Install Miniforge for Apple Silicon
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh

# Create new environment
conda create -n spatial python=3.10
conda activate spatial

# Install from conda-forge
conda install -c conda-forge numba llvmlite squidpy
```

**Solution 3: Temporary Workaround**
```bash
# Set library path before running
export DYLD_LIBRARY_PATH=/usr/lib:$DYLD_LIBRARY_PATH

# Or use single-threaded mode in squidpy
# Pass n_jobs=1 to neighborhood analysis functions
```

> **💡 Prevention Tips:**
> - On Apple Silicon Macs, prefer conda/miniforge over pip for scientific packages
> - Use ARM64-native conda distributions
> - Avoid mixing x86_64 and ARM64 packages
> 
> **⚠️ Impact:** This issue affects all squidpy spatial analysis functions that use parallel processing (neighborhood enrichment, co-occurrence analysis, etc.)

### MCP Connection Issues

#### Problem: MCP Server Won't Start

**Error messages:**
```bash
Failed to start MCP server
Connection refused
Server process exited with code 1
```

**Step-by-step diagnosis:**
1. **Test server manually:**
   ```bash
   python -m chatspatial --help
   # Should show help message without errors
   ```

2. **Check Python path in MCP config:**
   ```bash
   which python
   # Copy this EXACT path to your MCP configuration
   ```

3. **Verify MCP configuration format:**
   ```json
   {
     "mcpServers": {
       "chatspatial": {
         "command": "/full/path/to/python",
         "args": ["-m", "chatspatial"],
         "env": {}
       }
     }
   }
   ```

4. **Test with MCP Inspector:**
   ```bash
   npx @modelcontextprotocol/inspector python -m chatspatial
   ```

> **⚠️ When to Seek Help:** If MCP Inspector also fails, report the exact error message.

## Memory and Performance Problems

### Large Dataset Handling

#### Problem: Out of Memory Errors

**Error messages:**
```bash
MemoryError: Unable to allocate X GB
RuntimeError: CUDA out of memory
Process killed (signal 9)
```

**Immediate solutions:**
1. **Check available memory:**
   ```python
   import psutil
   print(f"Available RAM: {psutil.virtual_memory().available / (1024**3):.1f} GB")
   ```

2. **Subsample your data:**
   ```python
   # Keep every 10th cell for initial exploration
   adata_subset = adata[::10].copy()
   
   # Or randomly sample
   import numpy as np
   n_sample = 5000
   indices = np.random.choice(adata.n_obs, n_sample, replace=False)
   adata_subset = adata[indices].copy()
   ```

3. **Use sparse matrices:**
   ```python
   import scipy.sparse as sp
   if not sp.issparse(adata.X):
       adata.X = sp.csr_matrix(adata.X)
   ```

4. **Process in chunks:**
   ```text
   "Process my dataset in batches of 1000 cells"
   "Run analysis on spatial domains separately"
   ```

> **🔧 Long-term Strategies:**
> - **Increase swap space** (Linux/macOS): Add virtual memory
> - **Use high-memory nodes** on computing clusters
> - **Consider cloud computing** with more RAM (32GB+ for large datasets)

#### Problem: Slow Performance
**Symptoms:**
- Analysis takes hours instead of minutes
- LLM agent times out
- Progress bars move very slowly

**Optimization strategies:**

1. **Use Cherry Studio instead of Claude Desktop:**
   - Configure timeout to 3600 seconds
   - Better for computationally intensive tasks
   - More stable for long-running processes

2. **Enable multi-threading:**
   ```python
   import scanpy as sc
   sc.settings.n_jobs = 4  # Adjust based on your CPU cores
   ```

3. **Optimize specific methods:**
   ```text
   "Use Leiden clustering for fastest spatial domain identification"
   "Use SPARK-X for faster spatial gene detection on large datasets"
   "Use basic visualization instead of interactive plots"
   ```

> **💡 Prevention Tips:**
> - Start with small subsets for method testing
> - Use progress monitoring: "Show progress updates every 100 cells"
> - Monitor system resources during analysis

### GPU-Related Issues

#### Problem: CUDA/GPU Not Available

**Error messages:**
```bash
RuntimeError: No CUDA GPUs are available
torch.cuda.is_available() returns False
```

**Solutions:**
1. **Check GPU setup:**
   ```bash
   nvidia-smi  # Should show GPU information
   python -c "import torch; print(torch.cuda.is_available())"
   ```

2. **Install CUDA-compatible PyTorch:**
   ```bash
   pip install torch>=2.0.0 --index-url https://download.pytorch.org/whl/cu118
   ```

3. **Fall back to CPU:**
   ```text
   "Use CPU-only mode for analysis"
   "Disable GPU acceleration for this analysis"
   ```

> **ℹ️ When GPU Isn't Needed:** Most ChatSpatial analyses work fine on CPU. Only deep learning methods (Cell2location, scVI) benefit significantly from GPU.

## Analysis-Specific Issues

### Clustering and Spatial Domains

#### Problem: No Clusters Found or Too Many Clusters
**Symptoms:**
- All cells assigned to one cluster
- Hundreds of tiny clusters
- Spatial domains don't make biological sense

**Diagnostic approach:**
```text
"Show me the clustering parameters used"
"Visualize quality control metrics"
"Check if data needs better preprocessing"
```

**Solutions:**
1. **Adjust resolution parameters:**
   ```text
   "Try Leiden clustering with resolution 0.5"
   "Use SpaGCN with higher alpha parameter"
   "Increase k for STAGATE spatial domains"
   ```

2. **Improve preprocessing:**
   ```text
   "Filter genes present in at least 10 cells"
   "Remove cells with very low gene counts"
   "Apply log normalization before clustering"
   ```

3. **Try different methods:**
   ```text
   "Compare SpaGCN vs STAGATE results"
   "Try Leiden instead of Louvain clustering"
   ```

#### Problem: Cell Type Annotation Fails

**Error messages:**
```bash
No significant markers found
Reference dataset incompatible
Annotation method failed
```

**Step-by-step troubleshooting:**
1. **Check marker genes:**
   ```text
   "Show available marker genes for my tissue type"
   "List genes present in my dataset"
   "Use mouse markers instead of human"
   ```

2. **Try different annotation methods:**
   ```text
   "Try scType instead of marker-based annotation"
   "Use Tangram with reference data"
   "Try CellAssign with custom markers"
   ```

3. **Verify data compatibility:**
   ```text
   "Check if gene symbols match reference data"
   "Convert gene names to match annotation database"
   ```

### Visualization Problems

#### Problem: Plots Not Showing or Corrupted

**Error messages:**
```bash
Image could not be displayed
Visualization failed
Empty plot generated
```

**Solutions:**
1. **Check data availability:**
   ```text
   "Show me data summary before plotting"
   "Check if spatial coordinates exist"
   "Verify embedding was computed"
   ```

2. **Simplify visualization:**
   ```text
   "Create basic scatter plot instead of complex visualization"
   "Use default parameters for plotting"
   "Show first 1000 cells only"
   ```

3. **Try alternative plot types:**
   ```text
   "Create violin plot instead of spatial plot"
   "Use heatmap instead of UMAP"
   "Generate summary statistics table"
   ```

> **💡 Prevention Tip:** Always start with simple plots before creating complex visualizations.

### Cell Communication Analysis

#### Problem: CellPhoneDB KeyError with Complex Names

**Error messages:**
```bash
KeyError: 'ICAM3_integrin_aDb2_complex'
KeyError: 'GP1BA_integrin_aMb2_complex'
KeyError: 'BMP8A_ACVR_1A2B_receptor'
```

**Known Issue:** CellPhoneDB v5.0.1 has database inconsistencies with certain protein complexes.

**Solutions:**
1. **Use LIANA instead (recommended):**
   ```text
   "Use LIANA method for cell communication analysis"
   "Try LIANA with consensus database"
   "Use cell_type_handling='create_from_column' if needed"
   ```

2. **Enable gene filtering (partial fix):**
   ```text
   "Enable gene filtering with moderate strategy"
   "Use cellphonedb_use_microenvironments=true with filtering"
   ```

3. **Workarounds:**
   - Remove problematic genes manually (ICAM3, ITGAD, GP1BA, BMP8A)
   - Use older CellPhoneDB version (v4.x)
   - Switch to alternative databases

> **⚠️ Note:** This is a known bug in CellPhoneDB v5.0.1, not a ChatSpatial issue. The development team is aware of these database inconsistencies.

#### Problem: CellPhoneDB 'significant_means' KeyError

**Error message:**
```bash
KeyError: 'significant_means'
```

**Causes:**
- Insufficient ligand-receptor gene coverage
- Species mismatch (mouse genes vs human database)
- Dataset too sparse (< 5000 genes)

**Solutions:**
1. **Check gene coverage:**
   ```text
   "Check how many L-R genes are in my dataset"
   "Use data_source='raw' for more genes"
   ```

2. **Species issues:**
   ```text
   "Ensure species parameter matches your data"
   "Convert mouse gene names to human format if needed"
   ```

3. **Use alternatives:**
   ```text
   "Try LIANA which handles sparse data better"
   "Use CellChat via LIANA for mouse data"
   ```

### RNA Velocity and Trajectory Analysis

#### Problem: Velocity Analysis Fails

**Error messages:**
```bash
Missing spliced/unspliced layers
Velocity vectors could not be computed
Insufficient velocity genes
```

**Solutions:**
1. **Check velocity data:**
   ```text
   "Show me what layers are available in my data"
   "Check if RNA velocity preprocessing was done"
   "Verify spliced and unspliced counts exist"
   ```

2. **Use alternative approaches:**
   ```text
   "Try pseudotime analysis with Palantir instead"
   "Use DPT for trajectory inference"
   "Focus on spatial patterns without velocity"
   ```

## Data Format and Loading Problems

### Path Issues (Most Common!)

#### Problem: File Not Found - Wrong Path Type

> ⚠️ **CRITICAL REQUIREMENT**: ChatSpatial's MCP server **requires absolute paths**. Relative paths will always fail!

**Error messages:**
```bash
File not found
FileNotFoundError: [Errno 2] No such file or directory: 'data.h5ad'
Cannot locate file: ./folder/data.h5ad
```

**Most common mistake:**
```text
❌ Wrong: "Load data.h5ad"
❌ Wrong: "Load ./Downloads/data.h5ad"  
❌ Wrong: "Load ~/Downloads/data.h5ad"
✅ Correct: "Load /Users/yourname/Downloads/data.h5ad"
```

**Quick fix:**
1. **Find the absolute path:**
   ```bash
   cd /your/download/folder
   pwd  # Shows: /Users/yourname/Downloads
   ls *.h5ad  # Shows: card_spatial.h5ad
   ```

2. **Use in ChatSpatial:**
   ```text
   "Load /Users/yourname/Downloads/card_spatial.h5ad"
   ```

### File Format Issues

#### Problem: Data Won't Load

**Error messages:**
```bash
Unsupported file format
Corrupted data file
KeyError: spatial coordinates not found
```

**File format checklist:**
1. **Verify file exists and is accessible:**
   ```bash
   ls -la /Users/yourname/Downloads/data.h5ad
   file /Users/yourname/Downloads/data.h5ad  # Check file type
   ```

2. **Check file format support:**
   
   **Supported formats:**
   - `.h5ad` (recommended)
   - `.h5` (10X format) 
   - `.mtx + .tsv` (Matrix Market)
   - `.csv/.tsv` (text files)
   - `.xlsx` (Excel, limited)

3. **Test with ChatSpatial (remember: absolute paths only!):**
   ```text
   "Load data from /Users/yourname/Downloads/card_spatial.h5ad"
   "Check what file formats are supported"
   "Try loading with different format options"
   ```

#### Problem: Missing Spatial Information

**Error messages:**
```bash
No spatial coordinates found
KeyError: 'spatial'
Spatial analysis requires coordinates
```

**Solutions:**
1. **Check for spatial data:**
   ```text
   "Show me what's available in adata.obsm"
   "List all coordinate systems in the data"
   "Check if X_spatial or spatial key exists"
   ```

2. **Load coordinates separately:**
   ```text
   "Load coordinates from tissue_positions_list.csv"
   "Add spatial coordinates from separate file"
   "Use pixel coordinates as spatial data"
   ```

3. **Use non-spatial analysis:**
   ```text
   "Perform standard single-cell analysis instead"
   "Focus on cell type identification without spatial context"
   ```

### 10X Visium Specific Issues

#### Problem: Spatial Folder Structure Missing

**Expected structure:**
```bash
project/
├── filtered_feature_bc_matrix.h5
└── spatial/
    ├── tissue_positions_list.csv
    ├── scalefactors_json.json
    ├── tissue_hires_image.png
    └── tissue_lowres_image.png
```

**Solutions:**
1. **Check folder structure:**
   ```bash
   find /path/to/data -name "*.csv" -o -name "*.json" -o -name "*.png"
   ```

2. **Load components separately:**
   ```text
   "Load count matrix from filtered_feature_bc_matrix.h5"
   "Add spatial coordinates from tissue_positions_list.csv"
   "Include image data from spatial folder"
   ```

## Common Error Messages with Solutions

### "Dataset not found in data store"

**Meaning:** ChatSpatial can't find your loaded dataset.

**Solutions:**
1. `"Show me what datasets are currently loaded"`
2. `"Load my data again from /full/path"`
3. `"List available datasets in memory"`

### "Missing required dependencies"

**Meaning:** Optional package needed for specific analysis is not installed.

**Solutions:**
1. **Check what's needed:**
   ```bash
   pip list | grep -E "(chatspatial|numpy|pandas|scanpy|squidpy)"
   ```
2. **Install missing package:**
   ```bash
   pip install package_name
   ```
3. **Use alternative method:**
   ```text
   "Try different method that doesn't require X"
   ```

### "Memory allocation failed"

**Meaning:** Not enough RAM for the analysis.

**Quick fixes:**
1. `"Subsample data to 5000 cells"`
2. `"Use sparse matrices"`
3. `"Process data in smaller chunks"`

### "CUDA out of memory"

**Meaning:** GPU doesn't have enough memory.

**Solutions:**
1. `"Use CPU instead of GPU"`
2. `"Reduce batch size for deep learning"`
3. `"Free GPU memory before analysis"`

### "No significant results found"

**Meaning:** Statistical analysis didn't find meaningful patterns.

**Troubleshooting:**
1. `"Check data quality metrics"`
2. `"Try less stringent parameters"`
3. `"Use different statistical method"`

### "Visualization failed to render"

**Meaning:** Plot generation encountered an error.

**Solutions:**
1. `"Create simple scatter plot first"`
2. `"Check if data has required information"`
3. `"Use basic plotting parameters"`

## Conversation Troubleshooting

### Agent Doesn't Understand Your Request

#### Problem: Generic or Unhelpful Responses
**Symptoms:**
- Agent suggests irrelevant analyses
- Asks for clarification on clear requests
- Provides generic instead of specific help

**Better communication strategies:**

1. **Be specific about your data:**
   ```text
   ❌ "Analyze my data"
   ✅ "Load 10X Visium data from /path/data.h5ad and identify spatial domains using SpaGCN"
   ```

2. **Specify the biological context:**
   ```text
   ❌ "Find cell types"
   ✅ "Annotate cell types in mouse brain tissue using known cortical markers"
   ```

3. **Include parameters when needed:**
   ```text
   ❌ "Do clustering"
   ✅ "Perform Leiden clustering with resolution 0.8 and visualize on UMAP"
   ```

#### Problem: Agent Suggests Wrong Analysis Method
**Common mix-ups:**
- Suggests single-cell methods for spatial data
- Recommends incompatible analysis combinations
- Proposes computationally expensive methods unnecessarily

**Clarification approaches:**
```text
"I have spatial transcriptomics data, not single-cell RNA-seq"
"My dataset is large (50K cells), suggest fast methods"
"Focus on spatial-specific analyses only"
"I need methods compatible with Visium data"
```

### Parameter and Method Issues

#### Problem: Agent Chooses Inappropriate Parameters
**Symptoms:**
- Analysis fails due to incompatible settings
- Results don't make biological sense
- Methods take too long to complete

**Parameter guidance strategies:**
```text
"Use standard parameters for mouse brain Visium analysis"
"Choose parameters suitable for large datasets"
"Optimize for speed over precision"
"Use conservative statistical thresholds"
```

#### Problem: Method Compatibility Issues
**Common conflicts:**
- Trying to use methods requiring specific data types
- Combining incompatible analysis steps
- Using experimental methods with production data

**Compatibility checking:**
```text
"Check if my data is compatible with Cell2location"
"What methods work best with my data type?"
"Suggest analysis pipeline for 10X Visium data"
"List methods that work without GPU"
```

### Progress and Status Monitoring

#### Problem: Long Analysis with No Updates
**Solutions:**
1. **Request progress updates:**
   ```text
   "Show progress updates every minute"
   "Report when each analysis step completes"
   "Provide estimated time remaining"
   ```

2. **Use Cherry Studio for better monitoring:**
   - Set timeout to 3600 seconds
   - Better progress reporting
   - Can handle interruptions gracefully

3. **Break into smaller steps:**
   ```text
   "First load and preprocess data, then tell me when ready"
   "Do quality control first, then we'll proceed with analysis"
   ```

## When to Seek Help

### Gather Information First
Before seeking help, collect this information:
1. **System details:**
   ```bash
   python --version
   pip list | grep chatspatial
   uname -a  # Linux/macOS
   ```

2. **Error details:**
   - Full error message
   - What you were trying to do
   - Steps that led to the error

3. **Data information:**
   - File format and size
   - Data source (10X, Slide-seq, etc.)
   - Any preprocessing done

### Where to Get Help

1. **GitHub Issues:** For bugs and feature requests
   - Include system info and error messages
   - Provide minimal reproducible example
   - Check existing issues first

2. **Community Forums:** For usage questions
   - Describe your biological question
   - Include what you've already tried
   - Share relevant code snippets

3. **Documentation:** For method-specific help
   - Check tutorial sections
   - Review API documentation
   - Look at example workflows

### Emergency Workarounds

When you need results quickly:

1. **Use minimal installation:** `pip install -e .` for basic features
2. **Process small subsets:** Test with 1000 cells first
3. **Use simple methods:** Basic clustering and visualization
4. **Try different LLM client:** Cherry Studio vs Claude Desktop
5. **Run analysis step-by-step:** Break complex workflows into pieces

> **📝 Remember:** ChatSpatial is designed to be robust and provide helpful error messages. Most issues can be resolved by following the specific guidance above or trying alternative approaches suggested by the LLM agent.

## Quick Reference Card

### Installation Commands
```bash
# Basic installation
pip install -e .

# Advanced features  
pip install -e ".[advanced]"

# Check dependencies
pip list | grep -E "(chatspatial|numpy|pandas|scanpy|squidpy)"

# Test MCP server
python -m chatspatial --help
npx @modelcontextprotocol/inspector python -m chatspatial
```

### Common Agent Requests
```text
# Data loading
"Load 10X Visium data from /full/path/to/data.h5ad"

# Quick analysis
"Preprocess data with standard QC, then identify spatial domains"

# Troubleshooting
"Show me data summary and available analysis options"
"Check system resources and suggest optimization"

# Method alternatives
"Suggest fast methods for large dataset analysis"
"What methods work without GPU?"
```

This troubleshooting guide covers the most common issues you'll encounter with ChatSpatial. Keep it bookmarked for quick reference during your spatial transcriptomics analyses!