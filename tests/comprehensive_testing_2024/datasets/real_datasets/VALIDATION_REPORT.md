# 真实空间转录组数据集验证报告

验证时间: 2025-08-24 03:26:47

## 验证概览

- **总数据集**: 16
- **验证通过**: 16 ✅
- **有问题**: 0 ⚠️  
- **加载失败**: 0 ❌
- **成功率**: 100.0%

## 数据集详情

| 数据集 | 状态 | Spots | Genes | 技术 | 稀疏性 | 文件大小 | 演示数据 |
|--------|------|-------|-------|------|--------|----------|----------|
| hdst_squamous_carcinoma | ✅ valid | 3,850 | 1,530 | Unknown | 0.001 | 45.3 MB | 否 |
| osmfish_somatosensory | ✅ valid | 993 | 46 | Unknown | 0.638 | 0.4 MB | 否 |
| pixelseq_data | ✅ valid | 4,409 | 3,535 | Unknown | 0.016 | 119.7 MB | 否 |
| seqfish_demo | ✅ valid | 1,000 | 351 | Unknown | 0.726 | 1.6 MB | 否 |
| slideseq_cerebellum | ✅ valid | 15,000 | 10,000 | Slide-seq | 0.924 | 574.6 MB | 是 |
| slideseq_v2_hippocampus | ✅ valid | 8,000 | 12,000 | Slide-seq V2 | 0.924 | 368.1 MB | 是 |
| slideseq_v2_olfactory | ✅ valid | 6,000 | 11,000 | Slide-seq V2 | 0.924 | 253.3 MB | 是 |
| squidpy_imc | ✅ valid | 4,668 | 34 | Unknown | 0.015 | 1.5 MB | 否 |
| squidpy_merfish | ✅ valid | 73,655 | 161 | Unknown | 0.607 | 49.2 MB | 否 |
| squidpy_mibitof | ✅ valid | 3,309 | 36 | Unknown | 0.000 | 19.4 MB | 否 |
| squidpy_seqfish | ✅ valid | 19,416 | 351 | Unknown | 0.723 | 30.7 MB | 否 |
| squidpy_slideseqv2 | ✅ valid | 41,786 | 4,000 | Unknown | 0.974 | 251.3 MB | 否 |
| st_human_breast_cancer | ✅ valid | 600 | 8,000 | Spatial Transcriptomics | 0.654 | 19.1 MB | 是 |
| st_mouse_brain | ✅ valid | 1,000 | 15,000 | Spatial Transcriptomics | 0.654 | 58.6 MB | 是 |
| stereoseq_mouse_embryo | ✅ valid | 2,115 | 1,850 | Unknown | 0.000 | 30.2 MB | 否 |
| visium_demo | ✅ valid | 1,000 | 18,078 | Unknown | 0.685 | 164.4 MB | 否 |


## 技术分类统计

### 按技术类型分组
- **Slide-seq**: 1 个数据集
- **Slide-seq V2**: 2 个数据集
- **Spatial Transcriptomics**: 2 个数据集
- **Unknown**: 11 个数据集


### 数据规模分布

#### Spots数量分布
- **0-1,000**: 2 个数据集
- **1,000-5,000**: 8 个数据集
- **5,000-20,000**: 4 个数据集
- **20,000+**: 2 个数据集


## 质量评估

### 高质量数据集（推荐用于分析）
适合正式分析的验证通过数据集：

- **hdst_squamous_carcinoma**: 3,850 spots, 1,530 genes, Unknown
- **osmfish_somatosensory**: 993 spots, 46 genes, Unknown
- **pixelseq_data**: 4,409 spots, 3,535 genes, Unknown
- **seqfish_demo**: 1,000 spots, 351 genes, Unknown
- **slideseq_cerebellum**: 15,000 spots, 10,000 genes, Slide-seq
- **slideseq_v2_hippocampus**: 8,000 spots, 12,000 genes, Slide-seq V2
- **slideseq_v2_olfactory**: 6,000 spots, 11,000 genes, Slide-seq V2
- **squidpy_imc**: 4,668 spots, 34 genes, Unknown
- **squidpy_merfish**: 73,655 spots, 161 genes, Unknown
- **squidpy_mibitof**: 3,309 spots, 36 genes, Unknown
- **squidpy_seqfish**: 19,416 spots, 351 genes, Unknown
- **squidpy_slideseqv2**: 41,786 spots, 4,000 genes, Unknown
- **st_human_breast_cancer**: 600 spots, 8,000 genes, Spatial Transcriptomics
- **st_mouse_brain**: 1,000 spots, 15,000 genes, Spatial Transcriptomics
- **stereoseq_mouse_embryo**: 2,115 spots, 1,850 genes, Unknown
- **visium_demo**: 1,000 spots, 18,078 genes, Unknown


## 使用建议

### 数据加载
```python
import scanpy as sc

# 加载验证通过的数据集
adata = sc.read_h5ad('datasets/real_datasets/slideseq_v2_hippocampus.h5ad')

# 检查数据完整性
print(f"Spots: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Has spatial: {'spatial' in adata.obsm}")
```

### 数据过滤建议
对于有问题的数据集，建议进行以下过滤：
```python
# 过滤低质量spots
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 检查空间坐标
if 'spatial' in adata.obsm:
    coords = adata.obsm['spatial']
    # 移除坐标为NaN的spots
    valid_coords = ~np.isnan(coords).any(axis=1)
    adata = adata[valid_coords, :]
```

---
*此报告由 DatasetValidator 自动生成*
