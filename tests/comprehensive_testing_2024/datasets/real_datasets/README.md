# 真实空间转录组数据集

本目录包含用于comprehensive testing的真实空间转录组数据集。

## 数据集来源

### 1. Squidpy内置数据集
- 经过验证的标准数据集
- 格式标准化，可直接使用
- 包含多种空间转录组技术

### 2. 10X Genomics公开数据集  
- Visium技术的官方数据
- 高质量的人类和小鼠组织数据
- 包含完整的表达和空间信息

### 3. 已发表论文数据
- 来自高影响因子期刊的原始数据
- 代表各种空间转录组技术的最新进展
- 用于验证方法在真实数据上的性能

## 数据集列表

详细的数据集信息请查看 `datasets_metadata.json` 文件。

## 使用说明

所有数据集都保存为 `.h5ad` 格式，可直接用scanpy加载：

```python
import scanpy as sc

# 加载数据集
adata = sc.read_h5ad('squidpy_seqfish.h5ad')

# 查看基本信息
print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")
print(f"是否有空间坐标: {'spatial' in adata.obsm}")
```

## 引用信息

使用这些数据集时请引用相应的原始论文，详见metadata中的citation字段。

## 更新时间

数据集最后更新时间: 2025-08-24 03:23:45
