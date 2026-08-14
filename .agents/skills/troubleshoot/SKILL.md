---
name: troubleshoot
description: |
  Diagnose and resolve common issues in spatial transcriptomics analysis.
  Use when analysis fails, produces unexpected results, or user encounters errors.
  Triggers: "error", "failed", "not working", "issue", "problem", "help",
  "unexpected result", "debug", "fix", "wrong output".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob, Bash(python *)
---

# Troubleshooting Guide

## Overview

This skill helps diagnose and resolve common issues in ChatSpatial analysis.

## Quick Diagnosis Tree

```
Q: What type of problem?
│
├─ Analysis fails with error
│   └─ See: Error Messages section
│
├─ Analysis runs but results look wrong
│   └─ See: Result Validation section
│
├─ Performance issues (slow/memory)
│   └─ See: Performance section
│
└─ Dependency/installation problems
    └─ See: Environment section
```

## Error Messages

### Data Loading Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `FileNotFoundError` | Path incorrect | Verify file path exists |
| `KeyError: 'spatial'` | No spatial coordinates | Check data format, use correct loader |
| `ValueError: sparse matrix` | Format mismatch | Convert to compatible format |
| `MemoryError` | File too large | Use chunked loading or subsample |

**Diagnostic steps:**
```
1. Verify file exists: ls -la <path>
2. Check file size: du -h <path>
3. Try reading header only to test format
```

### Preprocessing Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `No highly variable genes` | Filtering too strict | Lower min_mean/max_mean thresholds |
| `Negative values after log` | Data not normalized | Use raw counts first |
| `Empty AnnData after QC` | QC too strict | Relax min_genes/min_cells |

### Method-Specific Errors

#### Deconvolution

| Error | Cause | Solution |
|-------|-------|----------|
| `No shared genes` | Reference mismatch | Check gene naming (symbols vs IDs) |
| `RCTD failed` | R/rpy2 issue | Verify R installation, check rpy2 |
| `Cell2location GPU error` | CUDA issue | Use CPU mode, check GPU availability |
| `Reference cell types missing` | Bad annotation column | Verify cell_type_key exists |

#### Cell Communication

| Error | Cause | Solution |
|-------|-------|----------|
| `CellChat R error` | R dependency | Install CellChat in R first |
| `No interactions found` | Stringent filtering | Lower p-value threshold |
| `Species database not found` | Wrong species config | Use "human" or "mouse" |

#### Trajectory/Velocity

| Error | Cause | Solution |
|-------|-------|----------|
| `No spliced/unspliced` | Missing layers | Process with velocyto first |
| `Root cell not found` | Bad root specification | Use valid cell barcode or cluster |
| `Disconnected graph` | Poor connectivity | Increase n_neighbors |

## Result Validation

### Unexpected Clustering Results

**Symptoms:**
- Too few/many clusters
- Clusters don't match histology
- Known cell types not separated

**Solutions:**
1. Adjust resolution parameter (higher = more clusters)
2. Try different clustering methods
3. Verify preprocessing quality
4. Check for batch effects

### Poor Deconvolution Results

**Symptoms:**
- All spots have same composition
- Proportions don't sum to 1
- Expected cell types missing

**Diagnostic checklist:**
- [ ] Reference and spatial data normalized consistently?
- [ ] Sufficient gene overlap (>2000 genes)?
- [ ] Reference cell types well-defined?
- [ ] Spatial data quality passed QC?

### Meaningless Communication Results

**Symptoms:**
- Random interactions without biological sense
- Same interactions everywhere
- Missing expected pathways

**Solutions:**
1. Verify cell type annotations are correct
2. Check species setting matches data
3. Try different communication database
4. Increase minimum expression threshold

## Performance Issues

### Memory Problems

**Symptoms:**
- `MemoryError`
- System becomes unresponsive
- Kernel crashes

**Solutions by data size:**

| Data Size | Strategy |
|-----------|----------|
| <50k cells | Should work, check other processes |
| 50k-100k | Close other applications, increase swap |
| 100k-500k | Use chunked processing, subsample |
| >500k | Requires HPC or cloud resources |

**Specific techniques:**
```
1. Subsample for exploration:
   sc.pp.subsample(adata, n_obs=50000)

2. Use sparse matrices:
   adata.X = scipy.sparse.csr_matrix(adata.X)

3. Process in chunks:
   Use batch processing for large datasets
```

### Slow Execution

**Common causes and solutions:**

| Bottleneck | Solution |
|------------|----------|
| Neighbor computation | Reduce n_neighbors, use approximate methods |
| Clustering | Use faster method (leiden > louvain) |
| Deconvolution | Use FlashDeconv instead of Cell2location |
| Visualization | Reduce point size, subsample for plotting |

## Environment Issues

### Python Dependencies

**Check installation:**
```bash
python -c "import scanpy; print(scanpy.__version__)"
python -c "import squidpy; print(squidpy.__version__)"
```

**Common fixes:**
```bash
pip install --upgrade scanpy squidpy
pip install scvi-tools  # For deep learning methods
```

### R Dependencies (rpy2)

**Check R availability:**
```python
import rpy2.robjects as ro
print(ro.r('R.version.string'))
```

**Install R packages:**
```r
# In R console
install.packages("Seurat")
BiocManager::install("RCTD")
devtools::install_github("sqjin/CellChat")
```

**Common rpy2 issues:**
```
1. R_HOME not set:
   export R_HOME=/path/to/R

2. Library path issues:
   Check .libPaths() in R

3. Version mismatch:
   Ensure rpy2 matches R version
```

### GPU/CUDA Issues

**Check GPU availability:**
```python
import torch
print(torch.cuda.is_available())
print(torch.cuda.device_count())
```

**Common fixes:**
- Update CUDA drivers
- Install correct torch version for CUDA
- Use CPU fallback: `use_gpu=False`

## Data Format Issues

### AnnData Structure Problems

**Verify structure:**
```python
print(adata)
print(adata.obs.columns.tolist())
print(adata.var.columns.tolist())
print(list(adata.obsm.keys()))
print(list(adata.layers.keys()))
```

**Common fixes:**

| Issue | Solution |
|-------|----------|
| Missing raw counts | `adata.raw = adata.copy()` before normalization |
| Wrong X matrix | Check if log-transformed when needed |
| Missing obsm keys | Re-run PCA/UMAP computation |

### Coordinate System Issues

**Spatial coordinates check:**
```python
if 'spatial' in adata.obsm:
    print(adata.obsm['spatial'].shape)
    print(adata.obsm['spatial'][:5])
```

**Common problems:**
- Coordinates in wrong units
- Flipped/rotated coordinates
- Missing spatial key

## Getting Help

### Information to Provide

When reporting issues, include:
1. **Error message**: Full traceback
2. **Data description**: n_cells, n_genes, platform
3. **Command used**: Exact parameters
4. **Environment**: Python version, package versions

### Quick Diagnostic Commands

```python
# Dataset summary
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Spatial: {'spatial' in adata.obsm}")
print(f"Layers: {list(adata.layers.keys())}")
print(f"Obs columns: {adata.obs.columns.tolist()[:10]}")
```
