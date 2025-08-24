# 真实空间转录组数据集使用指南

## 快速开始

### 1. 加载数据集

```python
import scanpy as sc
import pandas as pd

# 加载seqFISH数据集
adata = sc.read_h5ad('squidpy_seqfish.h5ad')

# 查看基本信息
print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")
print(f"空间坐标: {'spatial' in adata.obsm}")

# 查看空间坐标
if 'spatial' in adata.obsm:
    print(f"空间坐标形状: {adata.obsm['spatial'].shape}")
    print(f"坐标范围: x({adata.obsm['spatial'][:, 0].min():.1f}, {adata.obsm['spatial'][:, 0].max():.1f}), y({adata.obsm['spatial'][:, 1].min():.1f}, {adata.obsm['spatial'][:, 1].max():.1f})")
```

### 2. 可视化空间数据

```python
import matplotlib.pyplot as plt
import squidpy as sq

# 空间散点图
sq.pl.spatial_scatter(adata, color='total_counts')

# 如果有细胞类型信息
if 'cell_type' in adata.obs.columns:
    sq.pl.spatial_scatter(adata, color='cell_type')
```

## 推荐数据集

### 技术演示用途

**seqFISH+ (小规模)**
```python
# 1000个细胞的演示数据，适合快速测试
adata = sc.read_h5ad('seqfish_demo.h5ad')
```

**Visium (中等规模)**
```python
# 1000个spots，高基因数，适合算法测试
adata = sc.read_h5ad('visium_demo.h5ad')
```

### 大规模数据测试

**MERFISH (高细胞数)**
```python
# 73,655个细胞，适合性能测试
adata = sc.read_h5ad('squidpy_merfish.h5ad')
```

**Slide-seqV2 (高基因数)**
```python
# 41,786个细胞，4000个基因
adata = sc.read_h5ad('squidpy_slideseqv2.h5ad')
```

## 数据集选择指南

### 按技术类型

| 技术 | 数据集 | 细胞数 | 基因数 | 特点 |
|------|--------|--------|--------|------|
| seqFISH+ | `squidpy_seqfish.h5ad` | 19,416 | 351 | 单细胞分辨率 |
| MERFISH | `squidpy_merfish.h5ad` | 73,655 | 161 | 高通量 |
| Visium | `squidpy_visium.h5ad` | 684 | 18,078 | 高基因覆盖 |
| Slide-seqV2 | `squidpy_slideseqv2.h5ad` | 41,786 | 4,000 | 高细胞数 |

### 按组织类型

| 组织 | 数据集 | 物种 | 技术 |
|------|--------|------|------|
| 大脑皮层 | `ST_mouse_brain.h5ad` | 小鼠 | ST |
| 乳腺癌 | `st_human_breast_cancer.h5ad` | 人类 | ST |
| 海马 | `slideseq_v2_hippocampus.h5ad` | 小鼠 | Slide-seqV2 |
| 小脑 | `slideseq_cerebellum.h5ad` | 小鼠 | Slide-seq |

### 按数据大小

| 用途 | 推荐数据集 | 大小 | 加载时间 |
|------|------------|------|----------|
| 快速测试 | `seqfish_demo.h5ad` | 1.6 MB | <1秒 |
| 算法开发 | `squidpy_slideseqv2.h5ad` | 251 MB | ~5秒 |
| 性能测试 | `slideseq_cerebellum.h5ad` | 575 MB | ~10秒 |

## 常见分析流程

### 1. 质量控制

```python
# 计算质量指标
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# 可视化质量指标
sq.pl.spatial_scatter(adata, color='total_counts')
sq.pl.spatial_scatter(adata, color='n_genes_by_counts')
```

### 2. 数据预处理

```python
# 基本预处理
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

### 3. 空间分析

```python
# 构建空间邻接图
sq.gr.spatial_neighbors(adata)

# 计算空间自相关
sq.gr.spatial_autocorr(adata)

# 空间可变基因
sq.gr.spatial_autocorr(
    adata,
    mode='moran'
)
```

## 故障排除

### 常见问题

1. **文件读取失败**
   ```python
   # 检查文件是否存在
   from pathlib import Path
   if not Path('squidpy_seqfish.h5ad').exists():
       print("文件不存在，请检查路径")
   ```

2. **内存不足**
   ```python
   # 使用backed模式读取大文件
   adata = sc.read_h5ad('large_dataset.h5ad', backed='r')
   ```

3. **缺少空间坐标**
   ```python
   # 检查空间信息
   if 'spatial' not in adata.obsm:
       print("该数据集不包含空间坐标")
   ```

### 性能优化

```python
# 对于大数据集，可以先进行子采样
sc.pp.subsample(adata, n_obs=5000, random_state=42)

# 或者只加载部分基因
adata_subset = adata[:, :1000].copy()
```

## 元数据信息

每个数据集都包含丰富的元数据，可以通过以下方式查看：

```python
# 查看观察值(细胞)的元数据
print("细胞元数据列:")
print(adata.obs.columns.tolist())

# 查看基因元数据
print("基因元数据列:")
print(adata.var.columns.tolist())

# 查看附加信息
print("附加信息:")
print(list(adata.uns.keys()))

# 查看多维数据
print("多维观察数据:")
print(list(adata.obsm.keys()))
```

## 引用要求

使用这些数据集发表研究成果时，请引用相应的原始论文。详细的引用信息可以在以下文件中找到：
- `datasets_metadata.json` - 包含每个数据集的citation字段
- `REAL_DATASETS_REPORT.md` - 包含主要技术论文的完整引用

## 技术支持

如有问题或需要帮助：
1. 查看`validation_report.json`了解数据集验证状态
2. 查看`REAL_DATASETS_REPORT.md`了解详细技术信息
3. 通过GitHub issue报告问题

---
**更新时间**: 2025-08-24  
**版本**: 1.0