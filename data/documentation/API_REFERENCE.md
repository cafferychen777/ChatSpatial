# ChatSpatial 数据集 API 参考

## 数据加载 API

```python
# 基本加载
import scanpy as sc
adata = sc.read_h5ad('data/spatial_datasets/dataset.h5ad')

# 批量加载
import json
with open('data/metadata/comprehensive_catalog.json') as f:
    catalog = json.load(f)

# 按类别加载
spatial_datasets = catalog['datasets_by_category']['spatial_datasets']
for dataset in spatial_datasets:
    adata = sc.read_h5ad(f'data/{dataset["path"]}')
```

## 数据集查询 API

```python
# 查询空间数据集
spatial_datasets = [d for d in catalog['spatial_datasets'] if d['has_spatial']]

# 查询大规模数据集
large_datasets = [d for d in catalog['datasets_by_size']['large']]

# 查询特定技术
visium_datasets = catalog['datasets_by_technology']['visium']
```
