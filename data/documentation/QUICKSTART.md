# ChatSpatial 数据集快速开始

## 🚀 立即开始

```python
import scanpy as sc

# 1. 快速演示
adata = sc.read_h5ad('data/demo_datasets/visium_demo.h5ad')
print(f'Demo data: {adata.n_obs} cells, {adata.n_vars} genes')

# 2. 空间分析
adata = sc.read_h5ad('data/spatial_datasets/squidpy_merfish.h5ad')
print(f'Spatial data: {adata.n_obs} spots, {adata.n_vars} genes')
print('Spatial coordinates:', adata.obsm['spatial'].shape)

# 3. 大规模数据
adata = sc.read_h5ad('data/spatial_datasets/squidpy_slideseqv2.h5ad')
print(f'Large dataset: {adata.n_obs} cells')
```

## 📋 数据集选择指南

| 用途 | 推荐数据集 | 位置 |
|------|------------|------|
| 快速测试 | empty_velocity_layers.h5ad | spatial_datasets/ |
| 快速测试 | st_mouse_brain_backup_9MB.h5ad | spatial_datasets/ |
| 快速测试 | pancreas_subset_for_cellrank.h5ad | spatial_datasets/ |
| 空间分析 | squidpy_merfish.h5ad | spatial_datasets/ |
| 空间分析 | squidpy_slideseqv2.h5ad | spatial_datasets/ |
| 空间分析 | squidpy_seqfish.h5ad | spatial_datasets/ |
| 性能测试 | squidpy_slideseqv2.h5ad | spatial_datasets/ |
| 性能测试 | slideseq_cerebellum.h5ad | spatial_datasets/ |
| 性能测试 | squidpy_seqfish.h5ad | spatial_datasets/ |
