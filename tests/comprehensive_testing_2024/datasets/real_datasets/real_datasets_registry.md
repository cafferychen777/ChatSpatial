# ChatSpatial Real Datasets Registry

**Validation Date**: 2025-08-24T03:24:47.340188

## Summary Statistics

- **Total Datasets**: 35
- **Successfully Validated**: 22
- **Total Cells**: 71,887
- **Total Genes**: 101,613
- **Total Size**: 2581.5 MB
- **Spatial Datasets**: 13

## Quality Distribution

- **Medium**: 2 datasets
- **High**: 17 datasets
- **Low**: 3 datasets

## Compatibility Results

- **Tested**: 10 datasets
- **Load Success**: 10
- **Spatial Compatible**: 1
- **Preprocessing Compatible**: 10
- **Analysis Compatible**: 10

## Dataset Details

### Core Datasets

| Dataset | Cells | Genes | Size (MB) | Spatial | Quality |
|---------|-------|-------|-----------|---------|----------|
| slideseq_MOp_1217.h5ad | 9,852 | 24,518 | 923.2 | ✗ | Medium |

### Harmony Integration

| Dataset | Cells | Genes | Size (MB) | Spatial | Quality |
|---------|-------|-------|-----------|---------|----------|
| filtered_feature_bc_matrix.h5ad | 1,600 | 2,000 | 24.8 | ✗ | High |
| jurkat_293t_mixture_simulated.h5ad | 4,600 | 2,000 | 70.9 | ✗ | High |
| mixture_simulated.h5ad | 1,600 | 2,000 | 24.8 | ✗ | High |
| pure_293t_simulated.h5ad | 1,500 | 2,000 | 23.3 | ✗ | High |
| pure_jurkat_simulated.h5ad | 1,500 | 2,000 | 23.3 | ✗ | High |
| quick_demo_293t.h5ad | 500 | 800 | 3.2 | ✗ | High |
| quick_demo_combined.h5ad | 1,000 | 800 | 6.3 | ✗ | High |
| quick_demo_jurkat.h5ad | 500 | 800 | 3.2 | ✗ | High |

### Other

| Dataset | Cells | Genes | Size (MB) | Spatial | Quality |
|---------|-------|-------|-----------|---------|----------|
| hdst_squamous_carcinoma.h5ad | 3,850 | 1,530 | 45.3 | ✓ | High |
| osmfish_somatosensory.h5ad | 993 | 46 | 0.4 | ✓ | Low |
| pixelseq_data.h5ad | 4,409 | 3,535 | 119.7 | ✓ | High |
| slideseq_cerebellum.h5ad | 15,000 | 10,000 | 574.6 | ✓ | High |
| slideseq_v2_hippocampus.h5ad | 8,000 | 12,000 | 368.1 | ✓ | High |
| slideseq_v2_olfactory.h5ad | 6,000 | 11,000 | 253.3 | ✓ | High |
| squidpy_imc.h5ad | 4,668 | 34 | 1.5 | ✓ | Low |
| st_human_breast_cancer.h5ad | 600 | 8,000 | 19.1 | ✓ | High |
| st_mouse_brain.h5ad | 1,000 | 15,000 | 58.6 | ✓ | High |
| stereoseq_mouse_embryo.h5ad | 2,115 | 1,850 | 30.2 | ✓ | High |

### Test Datasets

| Dataset | Cells | Genes | Size (MB) | Spatial | Quality |
|---------|-------|-------|-----------|---------|----------|
| perf_test_100_500.h5ad | 100 | 500 | 0.5 | ✓ | Medium |
| synthetic_merfish.h5ad | 2,000 | 200 | 3.4 | ✓ | Low |
| synthetic_visium.h5ad | 500 | 1,000 | 4.0 | ✓ | High |

## Usage Instructions

```python
import scanpy as sc
from pathlib import Path

# Load a dataset
dataset_path = Path('datasets/real_datasets/core/ST_mouse_brain.h5ad')
adata = sc.read_h5ad(dataset_path)

# Check spatial coordinates
if 'spatial' in adata.obsm:
    print(f'Spatial coords: {adata.obsm["spatial"].shape}')
```
