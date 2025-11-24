
title: Data Formats
description: Supported data formats and preparation guide
---

# Data Formats Guide

Comprehensive guide to supported data formats and how to prepare your spatial transcriptomics data for ChatSpatial.

## Supported Formats

### Primary Formats

| Format | Extension | Description | Recommended Use |
|--------|-----------|-------------|-----------------|
| **AnnData** | `.h5ad` | HDF5-based format with spatial coordinates | **Primary and ONLY directly supported format** - Best performance |
| **10x H5** | `.h5` | 10x Genomics HDF5 format | 10x Visium data (auto-converts to AnnData) |
| **CSV** | `.csv` | Comma-separated values | **Requires manual conversion to .h5ad first** (see examples below) |

**Note**: All spatial data types (slide_seq, merfish, seqfish, other) ultimately load `.h5ad` files using scanpy's `read_h5ad()`. The `data_type` parameter only affects metadata labeling, not the underlying file format.

### Platform-Specific Data Types

| Platform | Required Format | Auto-Detection | Implementation |
|----------|----------------|----------------|----------------|
| **10x Visium** | Space Ranger directory or `.h5` file | ✅ | Native support via `sc.read_visium()` / `sc.read_10x_h5()` |
| **MERFISH** | `.h5ad` (pre-converted) | ✅ | Loads via `sc.read_h5ad()` |
| **Slide-seq** | `.h5ad` (pre-converted) | ✅ | Loads via `sc.read_h5ad()` |
| **seqFISH+** | `.h5ad` (pre-converted) | ✅ | Loads via `sc.read_h5ad()` |

**Important**: MERFISH, Slide-seq, and seqFISH+ data must be pre-converted to `.h5ad` format using scanpy/squidpy before loading into ChatSpatial. The `data_type` parameter only labels the technology type in metadata.

## Data Requirements

### Essential Components

1. **Gene Expression Matrix**
   - Genes × Spots/Cells
   - Raw or normalized counts
   - Sparse or dense format

2. **Spatial Coordinates**
   - X, Y positions for each spot/cell
   - Optional Z coordinate for 3D data
   - Pixel or physical coordinates

3. **Metadata** (Optional)
   - Cell type annotations
   - Batch information
   - Quality metrics
   - Experimental conditions

### Data Structure

```python
# AnnData structure
adata.X          # Gene expression matrix (n_obs × n_vars)
adata.obs        # Spot/cell metadata (n_obs × metadata)
adata.var        # Gene metadata (n_vars × gene_info)
adata.obsm       # Multi-dimensional annotations
  ['spatial']    # Spatial coordinates (n_obs × 2 or 3)
  ['X_pca']      # PCA coordinates (optional)
  ['X_umap']     # UMAP coordinates (optional)
adata.uns        # Unstructured annotations
  ['spatial']    # Spatial metadata (images, scale factors)
```

## Format-Specific Guides

### AnnData (.h5ad) - Recommended

**Advantages:**
- Native format for scanpy/squidpy
- Efficient storage and loading
- Preserves all metadata
- Best performance in ChatSpatial

**Example creation:**
```python
import scanpy as sc
import pandas as pd
import numpy as np

# Load expression data
expression = pd.read_csv('expression_matrix.csv', index_col=0)
coordinates = pd.read_csv('spatial_coordinates.csv', index_col=0)

# Create AnnData object
adata = sc.AnnData(X=expression.values)
adata.obs_names = expression.index
adata.var_names = expression.columns

# Add spatial coordinates
adata.obsm['spatial'] = coordinates[['x', 'y']].values

# Add metadata
adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()

# Save
adata.write('spatial_data.h5ad')
```

### CSV Format

**Structure:**
```
data/
├── expression_matrix.csv    # Genes × Spots
├── spatial_coordinates.csv  # Spot coordinates
└── metadata.csv            # Optional metadata
```

**Expression Matrix:**
```csv
,spot_1,spot_2,spot_3,...
gene_1,10,5,20,...
gene_2,0,15,8,...
gene_3,25,0,12,...
```

**Spatial Coordinates:**
```csv
spot_id,x,y
spot_1,100.5,200.3
spot_2,105.2,198.7
spot_3,110.1,205.9
```

**Loading CSV data:**
```python
# First convert CSV to AnnData format (recommended)
import pandas as pd
import scanpy as sc

# Load CSV files
expr = pd.read_csv("expression_matrix.csv", index_col=0)
coords = pd.read_csv("spatial_coordinates.csv", index_col=0)

# Create AnnData object
adata = sc.AnnData(X=expr.T)  # Transpose for spots × genes
adata.obsm['spatial'] = coords[['x', 'y']].values
adata.write('converted_data.h5ad')

# Then load with ChatSpatial
result = load_data(
    data_path="converted_data.h5ad",
    data_type="auto",  # or "other" for generic h5ad files
    name="my_data"
)
```

### 10x Visium Format

**Directory structure:**
```
visium_data/
├── filtered_feature_bc_matrix.h5
├── spatial/
│   ├── tissue_positions_list.csv
│   ├── scalefactors_json.json
│   ├── tissue_hires_image.png
│   └── tissue_lowres_image.png
└── analysis/
    └── clustering/
```

**Loading Visium data:**
```python
import scanpy as sc

# Load 10x Visium data
adata = sc.read_visium('path/to/visium_data')

# Make variable names unique
adata.var_names_unique()

# Save as h5ad for ChatSpatial
adata.write('visium_data.h5ad')
```

### HDF5 Format

**Custom HDF5 structure:**
```python
import h5py
import numpy as np

# Create HDF5 file
with h5py.File('spatial_data.h5', 'w') as f:
    # Expression data
    f.create_dataset('X', data=expression_matrix)
    f.create_dataset('obs_names', data=spot_names)
    f.create_dataset('var_names', data=gene_names)
    
    # Spatial coordinates
    f.create_dataset('spatial/coordinates', data=coordinates)
    
    # Metadata
    f.create_dataset('obs/total_counts', data=total_counts)
```

## Data Preparation

### Quality Control

```python
import scanpy as sc

# Load data
adata = sc.read_h5ad('raw_data.h5ad')

# Basic QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Filter by QC metrics
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# Save filtered data
adata.write('filtered_data.h5ad')
```

### Coordinate Normalization

```python
# Normalize spatial coordinates
coords = adata.obsm['spatial'].copy()

# Center coordinates
coords = coords - coords.mean(axis=0)

# Scale to unit variance
coords = coords / coords.std(axis=0)

# Update coordinates
adata.obsm['spatial'] = coords
```

### Data Conversion

#### CSV to AnnData
```python
def csv_to_anndata(expression_file, coordinates_file):
    import pandas as pd
    import scanpy as sc
    
    # Load data
    expr = pd.read_csv(expression_file, index_col=0)
    coords = pd.read_csv(coordinates_file, index_col=0)
    
    # Create AnnData
    adata = sc.AnnData(X=expr.T)  # Transpose for spots × genes
    adata.obsm['spatial'] = coords[['x', 'y']].values
    
    return adata
```

#### Seurat to AnnData
```r
# In R
library(Seurat)
library(SeuratDisk)

# Convert Seurat object to h5Seurat
SaveH5Seurat(seurat_obj, filename = "data.h5Seurat")
Convert("data.h5Seurat", dest = "h5ad")
```

```python
# In Python
import scanpy as sc
adata = sc.read_h5ad('data.h5ad')
```

## Validation

### Data Validation Checklist

```python
def validate_spatial_data(adata):
    """Validate spatial transcriptomics data."""
    
    checks = []
    
    # Check basic structure
    checks.append(("Has expression data", adata.X is not None))
    checks.append(("Has spatial coordinates", 'spatial' in adata.obsm))
    checks.append(("Coordinates shape", adata.obsm['spatial'].shape[1] >= 2))
    
    # Check data quality
    checks.append(("No empty spots", (adata.X.sum(axis=1) > 0).all()))
    checks.append(("No empty genes", (adata.X.sum(axis=0) > 0).all()))
    
    # Check coordinate validity
    coords = adata.obsm['spatial']
    checks.append(("Valid coordinates", not np.isnan(coords).any()))
    checks.append(("Positive coordinates", (coords >= 0).all()))
    
    # Print results
    for check, passed in checks:
        status = "✅" if passed else "❌"
        print(f"{status} {check}")
    
    return all(passed for _, passed in checks)

# Validate your data
validate_spatial_data(adata)
```

### Common Issues and Solutions

#### 1. Missing Spatial Coordinates
```python
# Solution: Add dummy coordinates
n_spots = adata.n_obs
adata.obsm['spatial'] = np.random.rand(n_spots, 2) * 1000
```

#### 2. Inconsistent Gene Names
```python
# Solution: Make gene names unique
adata.var_names_unique()
```

#### 3. Wrong Matrix Orientation
```python
# Solution: Transpose if needed
if adata.n_obs < adata.n_vars:
    adata = adata.T
```

#### 4. Large File Sizes
```python
# Solution: Use sparse matrices
import scipy.sparse as sp
adata.X = sp.csr_matrix(adata.X)
```

## Best Practices

### File Organization
```
project/
├── raw_data/
│   ├── sample1.h5ad
│   ├── sample2.h5ad
│   └── metadata.csv
├── processed_data/
│   ├── sample1_processed.h5ad
│   └── sample2_processed.h5ad
└── results/
    ├── analysis_results.h5ad
    └── plots/
```

### Naming Conventions
- Use descriptive file names: `mouse_brain_visium_sample1.h5ad`
- Include date and version: `data_2024_01_15_v1.h5ad`
- Separate raw and processed data
- Document processing steps

### Performance Tips
- Use `.h5ad` format for best performance
- Store large datasets on SSD
- Use sparse matrices for count data
- Compress files when storing long-term

## Next Steps

After preparing your data:

1. **Load into ChatSpatial**: Use the `load_data` tool
2. **Run Quality Control**: Check data quality metrics
3. **Start Analysis**: Follow the basic tutorial
4. **Explore Tools**: Try different analysis methods

See [Getting Started](../getting-started/README.md) for analysis workflows!