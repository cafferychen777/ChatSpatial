# 空间统计分析使用指南

## 概述

ChatSpatial 提供了两种方式来计算空间统计：

1. **analyze_spatial_patterns** - 用于多种空间分析类型（包括 Moran's I）
2. **calculate_spatial_statistics** - 专门用于高级空间统计（Geary's C、Local Moran's I）

## ⚠️ **重要更新**

**Moran's I 分析已统一到 `analyze_spatial_patterns()` 工具中**。请使用 `analysis_type="moran"` 而不是 `calculate_spatial_statistics()`。

## 1. analyze_spatial_patterns

用于执行各种空间分析，**包括 Moran's I 统计**。

### Moran's I 分析参数格式

```python
params = {
    "analysis_type": "moran",  # 使用这里进行 Moran's I 分析
    "morans_i_gene": "GENE_NAME",  # 指定要分析的基因
    "n_neighbors": 15
}
```

### 支持的 analysis_type 值

- `"neighborhood"` - 邻域分析
- `"co_occurrence"` - 共现分析
- `"ripley"` - Ripley's K 函数
- `"moran"` - **Moran's I 空间自相关** ← 使用这个
- `"centrality"` - 中心性分析
- `"getis_ord"` - Getis-Ord Gi* 热点分析

## 2. calculate_spatial_statistics

专门用于计算**高级空间统计**（不包括 Moran's I）。

### ⚠️ **注意**: 此工具不再支持 Moran's I

```python
# ❌ 错误 - 不再支持
result = calculate_spatial_statistics(
    statistic="morans_i"  # 这会报错
)

# ✅ 正确 - 使用 analyze_spatial_patterns
result = analyze_spatial_patterns(
    analysis_type="moran",
    morans_i_gene="Gad1"
)
```

### 使用示例

```python
# 计算单个基因的 Geary's C
result = calculate_spatial_statistics(
    data_id="data_1",
    feature="Gad1",
    statistic="gearys_c",
    n_neighbors=6
)
```

### 支持的统计类型

- `"gearys_c"` - 全局 Geary's C
- `"local_morans"` - 局部 Moran's I (LISA)

**注意**: Moran's I 已移至 `analyze_spatial_patterns()`

### 返回结果格式

对于 Geary's C：
```python
{
    "feature": "Gad1",
    "statistic": "gearys_c",
    "value": 0.456,          # 统计值
    "expected": -0.001,      # 期望值
    "variance": 0.002,       # 方差
    "z_score": 10.23,        # Z分数
    "p_value": 0.0001,       # P值
    "q_value": 0.001,        # 校正后的P值
    "n_neighbors": 6
}
```

对于局部统计（local_morans）：
```python
{
    "feature": "Gad1",
    "statistic": "local_morans",
    "local_values": [...],    # 每个点的局部I值
    "p_values": [...],        # 每个点的P值
    "clusters": [...],        # 聚类分配
    "n_neighbors": 6
}
```

## 批量分析多个基因

如果需要分析多个基因，可以循环调用：

```python
genes = ["Gad1", "Camk2a", "Slc17a7", "Pvalb", "Sst"]
results = []

for gene in genes:
    try:
        result = calculate_spatial_statistics(
            data_id="data_1",
            feature=gene,
            statistic="morans_i",
            n_neighbors=6
        )
        results.append(result)
    except Exception as e:
        print(f"Error analyzing {gene}: {e}")
```

## 依赖项

这些功能需要安装 `pysal` 库：
```bash
pip install pysal
```

## 注意事项

1. 基因名称必须存在于数据集中
2. 空间坐标必须存储在 `adata.obsm['spatial']`
3. n_neighbors 参数影响空间权重矩阵的构建
4. 较小的 p_value 表示显著的空间模式

---

*创建日期：2025年1月*