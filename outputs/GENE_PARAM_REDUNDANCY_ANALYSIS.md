# 基因数量参数冗余性分析报告

## 问题发现

用户报告：Moran's I 请求 n_top_genes=100，但只分析了 20 个基因

## 根本原因

### 1. 参数设计不一致

```python
# 5个 Parameters 类中有 4 个定义了 n_top_genes
RNAVelocityParameters.n_top_genes = 2000
DeconvolutionParameters.n_top_genes = 2000  
SpatialVariableGenesParameters.n_top_genes = Optional[int]
EnrichmentParameters: min_genes, max_genes (不同用途)

# ❌ 唯独 SpatialAnalysisParameters 没有！
SpatialAnalysisParameters:
  - moran_n_genes = 20 (max=100)
  - getis_ord_n_genes = 20 (max=100)
  - ❌ 没有 n_top_genes
```

### 2. 用户传入参数被忽略

```python
# 用户调用
analyze_spatial_statistics(
    data_id="data_1",
    params={"analysis_type": "moran", "n_top_genes": 100}  # ← 被忽略
)

# Pydantic 处理
class SpatialAnalysisParameters(BaseModel):
    moran_n_genes: int = 20  # ← 使用默认值
    # n_top_genes 未定义 → 被丢弃
```

### 3. Workaround 代码暴露设计问题

```python
# spatial_statistics.py:465
n_genes_to_analyze = getattr(params, 'n_top_genes', None) or params.moran_n_genes
# ↑ 这个 getattr 尝试获取不存在的属性 → 设计缺陷的证据
```

## 参数冗余度分析

### 当前状态

| 方法 | 专用参数 | 默认值 | 上限 | 实际作用 |
|------|----------|--------|------|----------|
| Moran's I | `moran_n_genes` | 20 | 100 | ✅ 有效使用 |
| Geary's C | `moran_n_genes` | 20 | 100 | ✅ 复用 Moran 参数 |
| Getis-Ord | `getis_ord_n_genes` | 20 | 100 | ✅ 有效使用 |
| Local Moran | - | - | - | ❌ 未使用任何参数 |
| Bivariate Moran | - | - | - | ✅ 使用 gene_pairs |
| 其他方法 | - | - | - | ❌ 不分析基因 |

**结论**：只有 **3个方法** 需要基因数量参数，却定义了 **2个专用参数**

### 冗余度评估

```
专用参数使用率 = 3个使用方法 / 12个总方法 = 25%

实际需要专用参数的理由：
  - Moran's I: 默认20，计算开销大
  - Getis-Ord: 默认20，计算开销大  
  - Local Moran: 可能需要不同默认值？(但当前未实现)

是否需要分开？
  ❌ 不需要！两者都是空间自相关分析，默认值相同，上限相同
```

## 代码实证

### 使用情况统计

```bash
$ grep -c "moran_n_genes" chatspatial/tools/spatial_statistics.py
4  # 使用4次

$ grep -c "getis_ord_n_genes" chatspatial/tools/spatial_statistics.py  
1  # 使用1次

$ grep -c "n_top_genes" chatspatial/tools/spatial_statistics.py
1  # 尝试读取1次 (失败的 workaround)
```

### 使用模式分析

```python
# Moran's I (line 465)
n_genes_to_analyze = getattr(params, 'n_top_genes', None) or params.moran_n_genes
# ↑ 尝试使用通用参数，回退到专用参数

# Geary's C (line 532)
genes = adata.var_names[adata.var["highly_variable"]][:params.moran_n_genes].tolist()
# ↑ 复用 moran_n_genes

# Getis-Ord (line 669)
genes = adata.var_names[adata.var["highly_variable"]][:params.getis_ord_n_genes].tolist()
# ↑ 唯一使用 getis_ord_n_genes 的地方

# Local Moran (line 1030+)
# ❌ 未使用任何基因数量参数！直接分析指定基因列表
```

**发现**：
- `moran_n_genes` 被 2 个方法共享 (Moran, Geary)
- `getis_ord_n_genes` 只被 1 个方法使用
- 没有方法真正需要独立的专用参数

## 推荐方案

### 方案 A: 添加 n_top_genes + 保留向后兼容 (推荐 ⭐)

```python
class SpatialAnalysisParameters(BaseModel):
    # 通用参数 (优先级最高)
    n_top_genes: Optional[Annotated[int, Field(gt=0, le=500)]] = None
    
    # 方法专用参数 (向后兼容，n_top_genes=None 时使用)
    moran_n_genes: int = 20
    getis_ord_n_genes: int = 20
    
# 使用逻辑 (spatial_statistics.py)
n_genes = params.n_top_genes if params.n_top_genes is not None else params.moran_n_genes
```

**优点**：
✅ 与其他 Parameters 类一致  
✅ 用户体验友好 (统一参数名)  
✅ 向后兼容 (已有代码不受影响)  
✅ 灵活性：可以全局设置或单独设置  

**缺点**：
⚠️ 增加一个参数字段  

---

### 方案 B: 统一为 n_top_genes (激进)

```python
class SpatialAnalysisParameters(BaseModel):
    # 统一参数
    n_top_genes: Annotated[int, Field(gt=0, le=500)] = 20
    
    # ❌ 移除: moran_n_genes, getis_ord_n_genes
```

**优点**：
✅ 最简洁  
✅ 减少参数数量  
✅ 消除冗余  

**缺点**：
❌ **破坏性变更** (影响已有用户)  
❌ 失去方法级默认值灵活性  

---

### 方案 C: 重命名统一参数名 (不推荐)

```python
class SpatialAnalysisParameters(BaseModel):
    # 统一为 n_genes (与 genes 对应)
    genes: Optional[List[str]] = None
    n_genes: Optional[int] = None
```

**优点**：
✅ 语义清晰  

**缺点**：
❌ 与其他类不一致 (都叫 n_top_genes)  
❌ 容易混淆 (n_genes vs n_hvgs vs n_top_genes)  

## 最终建议

采用 **方案 A**，理由：

1. **修复当前 Bug**：解决 n_top_genes=100 被忽略的问题
2. **API 一致性**：与其他 5 个 Parameters 类保持一致
3. **向后兼容**：不破坏现有代码
4. **用户体验**：提供直观的参数名
5. **最小改动**：只需添加 1 个字段 + 修改 3 处使用逻辑

### 实现步骤

1. **修改 data.py** (添加 n_top_genes)
2. **修改 spatial_statistics.py** (更新 3 个方法的参数读取逻辑)
3. **添加单元测试** (验证参数优先级)
4. **更新文档** (说明参数优先级)

### 参数优先级

```python
# 优先级: n_top_genes > 方法专用参数 > 默认值
if params.n_top_genes is not None:
    n = params.n_top_genes
elif params.analysis_type == "moran":
    n = params.moran_n_genes
elif params.analysis_type == "getis_ord":
    n = params.getis_ord_n_genes
else:
    n = 20  # fallback
```

## 长期重构建议

未来版本可以考虑：
1. 废弃 `moran_n_genes` 和 `getis_ord_n_genes` (添加 deprecation 警告)
2. 统一所有基因选择参数为 `n_top_genes`
3. 提供版本迁移指南

但现在 **不应该破坏向后兼容性**。
