# Troubleshooting Guide

Quick solutions to common problems.

---

## Quick Problem Finder

Jump to:
- [MCP Connection Issues](#mcp-connection-issues)
- [Installation Problems](#installation-problems)
- [Data Loading Errors](#data-loading-errors)
- [Analysis Failures](#analysis-failures)
- [Memory Issues](#memory-issues)

---

## MCP Connection Issues

### Problem: ChatSpatial not showing up in Claude

**Symptoms:**
- No ChatSpatial tools appear in Claude
- MCP server not connecting

**Solution:**

1. **Check Python path is correct:**
   ```bash
   # Activate virtual environment
   source chatspatial_env/bin/activate
   which python
   # Copy this path to your MCP config
   ```

2. **Verify configuration file:**
   - macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
   - Windows: `%APPDATA%\Claude\claude_desktop_config.json`
   - Linux: `~/.config/Claude/claude_desktop_config.json`

3. **Restart Claude Desktop** after any configuration changes

4. **Test ChatSpatial installation:**
   ```bash
   python -m chatspatial server --help
   # Should display server options
   ```

---

## Installation Problems

### Problem: ImportError or ModuleNotFoundError

**Common errors:**
```
ImportError: No module named 'chatspatial'
ModuleNotFoundError: No module named 'scanpy'
```

**Solution:**

1. **Check you're in the virtual environment:**
   ```bash
   which python
   # Should show path to virtual environment
   ```

2. **Reinstall ChatSpatial:**
   ```bash
   pip install --upgrade pip
   pip install -e ".[full]"
   ```

3. **Verify installation:**
   ```bash
   python -c "import chatspatial; print('OK')"
   ```

### Problem: R-based methods failing (RCTD, SPOTlight, etc.)

**Error:**
```
rpy2 not available or R packages not installed
```

**Solution:**

1. **Install R:** [Download from CRAN](https://cran.r-project.org/)

2. **Install R packages:**
   ```r
   # In R console
   install.packages("spacexr")  # For RCTD
   install.packages("SPOTlight")  # For SPOTlight
   ```

3. **Verify R installation:**
   ```bash
   python -c "import rpy2; print('R available')"
   ```

**Note:** R-based methods are optional. Use Python alternatives if R installation is problematic.

---

## Data Loading Errors

### Problem: "Dataset not found" or path errors

**Most common issue!** ChatSpatial requires **absolute paths**.

**Wrong:**
```
~/data/sample.h5ad
./data/sample.h5ad
data/sample.h5ad
```

**Correct:**
```
/Users/yourname/data/sample.h5ad
/home/yourname/data/sample.h5ad
C:/Users/yourname/data/sample.h5ad
```

**Solution:**
```bash
# Get absolute path
pwd
# Use full path: /Users/yourname/current/directory/data.h5ad
```

### Problem: File format not recognized

**Error:**
```
ValueError: Could not read file format
```

**Solution:**

1. **Check file format:**
   ```bash
   file yourdata.h5ad
   # Should show HDF5 data
   ```

2. **For Visium data:** Point to the directory, not individual files
   ```
   /path/to/visium_output/  (directory containing spatial/ and filtered_feature_bc_matrix.h5)
   ```

3. **For H5AD files:** Ensure they're valid AnnData objects
   ```python
   import scanpy as sc
   adata = sc.read_h5ad("yourdata.h5ad")  # Test in Python first
   ```

---

## Analysis Failures

### Problem: "Run preprocessing first"

**Error:**
```
Analysis requires preprocessed data
```

**Solution:**

Ask Claude: "Preprocess the data first"

Preprocessing computes necessary information (PCA, neighbors, clustering) required by most analysis tools.

### Problem: "No significant results found"

**Common in:**
- Differential expression
- Spatial variable genes
- Cell communication

**Solutions:**

1. **Check your data quality:**
   - Sufficient cells/spots (>500 recommended)
   - Adequate gene coverage (>1000 genes)
   - Proper normalization

2. **Adjust parameters:**
   - Lower significance thresholds
   - Increase number of neighbors
   - Try different methods

3. **Verify biological expectations:**
   - Some datasets may genuinely lack strong signals
   - Check positive controls

### Problem: Cell communication analysis fails

**Error:**
```
Too few features from resource found in data
```

**Solution:**

1. **Check species parameter:**
   ```
   For mouse data: species="mouse", liana_resource="mouseconsensus"
   For human data: species="human", liana_resource="consensus"
   ```

2. **Ensure preprocessing was done** (creates adata.raw for comprehensive gene coverage)

3. **Verify cell type annotations exist** in your data

---

## Memory Issues

### Problem: "Memory allocation failed" or system freezes

**Symptoms:**
- Analysis crashes mid-run
- System becomes unresponsive
- "MemoryError" or "killed" messages

**Solutions:**

1. **Reduce dataset size for testing:**
   - Subsample spots/cells
   - Use fewer genes (top HVGs)

2. **Use memory-efficient parameters:**
   - Smaller batch sizes
   - Fewer epochs for deep learning
   - Chunk processing

3. **Monitor memory usage:**
   ```bash
   # macOS/Linux
   top
   # Press 'q' to quit
   ```

4. **For large datasets:**
   - Use server with more RAM
   - Consider cloud computing (AWS, GCP)
   - Use methods optimized for large data

### Problem: "CUDA out of memory" (GPU)

**Error:**
```
RuntimeError: CUDA out of memory
```

**Solutions:**

1. **Reduce batch size** in method parameters

2. **Use CPU instead:**
   Set `use_gpu=False` in parameters

3. **Free GPU memory:**
   ```python
   import torch
   torch.cuda.empty_cache()
   ```

---

## Common Error Messages

### "Missing required dependencies"

**Solution:** Install the required feature group:
```bash
pip install -e ".[full]"  # All features
```

### "Server not responding"

**Solution:**
1. Check MCP logs for errors
2. Restart Claude Desktop
3. Verify Python path in config

### "JSON decode error"

**Solution:** Check MCP configuration file syntax - missing commas, brackets, quotes

---

## When to Seek Help

If problems persist after trying these solutions:

1. **Check GitHub Issues:** [Existing issues](https://github.com/cafferychen777/ChatSpatial/issues)

2. **Report new issue** with:
   - Error message (full traceback)
   - Steps to reproduce
   - System info (OS, Python version)
   - Dataset info (size, format)

3. **Join discussions:** [GitHub Discussions](https://github.com/cafferychen777/ChatSpatial/discussions)

---

## Quick Reference

### Installation
```bash
# Create environment
python3 -m venv chatspatial_env
source chatspatial_env/bin/activate

# Install
pip install -e ".[full]"

# Verify
python -c "import chatspatial; print('OK')"
```

### MCP Configuration (Claude Desktop)
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/full/path/to/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```

### Common Analysis Flow
```
1. Load data (absolute path!)
2. Preprocess data
3. Run analysis
4. Visualize results
```

---

**Still stuck?** [Open a GitHub issue](https://github.com/cafferychen777/ChatSpatial/issues/new)
