# SingleR Integration Documentation

## Overview

The `marker_genes` annotation method has been enhanced with SingleR integration, providing more accurate cell type annotation using reference-based correlation analysis instead of simple marker gene overlap.

## Key Improvements

### Old Method (marker_genes with scanpy)
- **Approach**: Simple overlap counting between cluster markers and reference markers
- **Resolution**: Cluster-level annotation
- **Accuracy**: Limited by marker gene selection
- **Statistical Method**: Jaccard index or overlap count

### New Method (SingleR)
- **Approach**: Spearman correlation with reference expression profiles
- **Resolution**: Single-cell level annotation
- **Accuracy**: Higher, uses full transcriptome information
- **Statistical Method**: Correlation analysis with fine-tuning

## Installation

```bash
# Install SingleR and dependencies
pip install singler singlecellexperiment celldex

# Or install with advanced features
pip install "chatspatial[advanced]"
```

## Usage

### Basic Usage with Default Reference

```python
# The marker_genes method now automatically uses SingleR if available
params = {
    "method": "marker_genes",
    "marker_genes": {
        "T_cells": ["CD3D", "CD3E", "CD3G"],
        "B_cells": ["CD19", "MS4A1", "CD79A"],
        # ... more markers
    }
}
```

### Using Pre-built References

```python
params = {
    "method": "marker_genes",
    "singler_reference": "blueprint_encode",  # or "dice", "hpca", etc.
    "num_threads": 4
}
```

### Using Custom Reference Dataset

```python
# First load reference data
reference_params = {
    "data_path": "path/to/reference.h5ad",
    "data_type": "h5ad",
    "name": "MyReference"
}

# Then use for annotation
annotation_params = {
    "method": "marker_genes",
    "reference_data_id": "reference_data_id",
    "num_threads": 4
}
```

### Integrated Annotation with Multiple References

```python
params = {
    "method": "marker_genes",
    "singler_reference": "blueprint_encode",
    "singler_integrated": true,
    "num_threads": 6
}
```

## Available Reference Datasets

Through `celldex`, the following pre-built references are available:

1. **blueprint_encode**: Immune and blood cells
2. **dice**: Immune cells from DICE project
3. **hpca**: Human Primary Cell Atlas
4. **monaco_immune**: Monaco immune dataset
5. **immgen**: Mouse immune cells
6. **mouse_rnaseq**: Mouse tissues

## Fallback Mechanism

If SingleR is not installed or fails, the system automatically falls back to the original scanpy marker gene overlap method, ensuring backward compatibility.

## Performance Considerations

- **Memory**: SingleR requires more memory than simple marker overlap
- **Speed**: SingleR is slower but more accurate
- **Parallelization**: Use `num_threads` parameter for faster processing

## Comparison Results

| Method | Accuracy | Speed | Memory | Resolution |
|--------|----------|-------|--------|------------|
| Old marker_genes | Medium | Fast | Low | Cluster-level |
| SingleR | High | Moderate | High | Cell-level |

## Technical Details

### Algorithm Steps

1. **Data Conversion**: AnnData → SingleCellExperiment
2. **Reference Preparation**: Load or create reference profiles
3. **Correlation Analysis**: Compute Spearman correlation
4. **Fine-tuning**: Refine assignments for ambiguous cells
5. **Result Integration**: Add annotations back to AnnData

### Output Format

The method returns:
- `cell_types`: List of unique cell types
- `counts`: Dictionary of cell type counts
- `confidence_scores`: Dictionary of confidence scores per type
- Results stored in `adata.obs['singler_celltype']`

## Troubleshooting

### SingleR Not Available
```
Error: SingleR is not installed
Solution: pip install singler singlecellexperiment
```

### Reference Not Found
```
Error: Could not find labels in reference
Solution: Ensure reference data has cell type annotations
```

### Memory Issues
```
Error: Out of memory
Solution: Reduce num_threads or use smaller reference
```

## Citation

If using SingleR, please cite:
```
Aran et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology.
```