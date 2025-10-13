# Deconvolution "0 genes in common" 问题修复总结

## 问题根源

Deconvolution方法(Cell2location, DestVI, RCTD等)报错"Only 0 genes in common"，原因是：

1. **数据验证顺序错误**（已修复）
   - 原逻辑：先验证当前var_names → 再prepare counts → 用空的common_genes subset
   - 问题：验证在prepare之前，检查的是HVG子集(可能0重叠)

2. **_prepare_anndata_for_counts恢复数据不完整**（已修复）
   - 原逻辑：只提取当前var_names对应的counts，不恢复完整基因集
   - 问题：即使.raw有2000个基因，只返回当前的150个HVG

3. **preprocessing.py的.raw保存机制问题**（已修复 - 根本原因！）
   - 原代码：`adata.raw = adata` 创建了一个VIEW，不是frozen copy
   - 问题：后续normalization修改adata.X时，也修改了.raw.X
   - 结果：.raw虽然保留所有基因，但数据是normalized而非counts

## 修复内容

### 1. deconvolution.py - 移除提前验证

**修改的函数** (6个):
- `run_cell2location_deconvolution` (line 584-630)
- `run_rctd_deconvolution` (line 882-928)
- `run_destvi_deconvolution` (line 1800-1828)
- `run_stereoscope_deconvolution` (line 1994-2023)
- `run_spotlight_deconvolution` (line 2193-2232)
- `run_tangram_deconvolution` (line 2472-2501)

**新流程**:
```python
# 1. 验证cell type key存在
if cell_type_key not in reference_adata.obs:
    raise ValueError(...)

# 2. 先prepare (恢复完整基因集)
ref = _prepare_anndata_for_counts(reference_adata.copy(), "Reference")
sp = _prepare_anndata_for_counts(spatial_adata.copy(), "Spatial")

# 3. 在prepare后找common genes
common_genes = list(set(ref.var_names) & set(sp.var_names))

# 4. 此时才验证common genes是否足够
if len(common_genes) < min_common_genes:
    raise ValueError(...)

# 5. Subset到common genes
ref = ref[:, common_genes].copy()
sp = sp[:, common_genes].copy()
```

### 2. deconvolution.py - 修改_prepare_anndata_for_counts

**位置**: `chatspatial/tools/deconvolution.py` line 215-250

**新优先级逻辑**:
1. **检查.raw是否包含counts** (sample检查，避免全量读取)
   - 如果是counts → 返回完整.raw.to_adata() (所有基因)
   - 如果是normalized → 跳过，尝试layers['counts']
2. **检查layers['counts']**
   - 如果存在 → 使用(但只有当前var_names)
3. **使用当前.X**
   - Fallback

**关键改进**:
```python
# Step 1: 先验证.raw是否是counts
if adata_copy.raw is not None:
    raw_adata = adata_copy.raw.to_adata()

    # Sample检查避免全量读取
    sample_X = raw_X[:100, :100]
    has_decimals = not np.allclose(sample_X, np.round(sample_X), atol=1e-6)
    has_negatives = sample_X.min() < 0

    if not has_decimals and not has_negatives:
        # Raw是counts，返回完整基因集
        adata_copy = raw_adata  # 所有基因!
        data_source = "raw"
```

### 3. preprocessing.py - 正确保存.raw counts (根本修复!)

**位置**: `chatspatial/tools/preprocessing.py` line 313-322

**问题代码**:
```python
# 旧代码 (line 314)
adata.raw = adata  # ← 创建VIEW，后续normalization会修改它!
```

**修复代码**:
```python
# 新代码 (line 315-322)
import anndata as ad_module

adata.raw = ad_module.AnnData(
    X=adata.X.copy(),  # ← 真正的copy!
    var=adata.var.copy(),
    obs=adata.obs.copy(),
    uns=adata.uns.copy()
)
```

**为什么这个修复至关重要**:
- 旧代码：`adata.raw = adata` 创建的是VIEW，与主对象共享底层数据
- 新代码：创建独立的AnnData对象，完全隔离
- 结果：normalization(line 403-417)只修改adata.X，不影响adata.raw.X
- 最终：HVG subsetting(line 596)后，.raw保留完整基因集的counts!

## 数据流验证

使用preprocessing.py处理后的数据结构:

```
1. 初始状态 (line 313):
   adata.X: 1000 genes, counts

2. 保存raw (line 315-322):
   adata.raw: 1000 genes, counts (独立copy)
   adata.layers['counts']: 1000 genes, counts (copy)

3. 归一化后 (line 403-417):
   adata.X: 1000 genes, normalized
   adata.raw.X: 1000 genes, counts (未受影响!)
   adata.layers['counts']: 1000 genes, counts (未受影响)

4. HVG subsetting后 (line 596):
   adata.X: 200 genes (HVG), normalized
   adata.raw.X: 1000 genes, counts (完整基因集!)
   adata.layers['counts']: 200 genes (HVG), counts (被subset)
```

## Deconvolution流程(修复后)

```
用户数据 (两个数据集都经过preprocessing):
- ref: 150 HVG genes (normalized)
  - ref.raw: 2000 genes (counts)
- spatial: 150 HVG genes (normalized)
  - spatial.raw: 13948 genes (counts)

当前var_names共同基因: 可能是0 (不同的HVG)

Deconvolution调用:
1. _prepare_anndata_for_counts(ref)
   → 检查ref.raw是counts → 返回2000 genes的counts

2. _prepare_anndata_for_counts(spatial)
   → 检查spatial.raw是counts → 返回13948 genes的counts

3. 找common genes
   → set(2000) & set(13948) = 1888 genes (充足!)

4. Subset到common genes
   → ref_final: 1888 genes
   → spatial_final: 1888 genes

5. 运行deconvolution
   → 成功!
```

## 测试验证

### 已通过的测试脚本:

1. **test_preprocessing_fix.py** - 验证.raw正确保存
   - ✓ .raw在normalization后仍然是counts
   - ✓ HVG subsetting后保留完整基因集

2. **test_prepare_with_proper_data.py** - 验证_prepare行为
   - ✓ .raw是counts时返回完整基因集
   - ✓ .raw是normalized时fallback到layers['counts']

### 待测试 (需要重启MCP):

使用真实MCP workflow测试:
1. Load两个数据集
2. Preprocess (会创建正确的.raw)
3. 模拟极端HVG不重叠场景
4. 运行deconvolution (Cell2location, RCTD, DestVI等)
5. 验证成功运行

## 影响范围

### 改动的文件:
1. `chatspatial/tools/deconvolution.py`
   - 修改了6个deconvolution函数
   - 修改了_prepare_anndata_for_counts

2. `chatspatial/tools/preprocessing.py`
   - 修改了.raw保存逻辑 (line 315-322)

### 向后兼容性:
- ✓ 对已有数据：如果.raw是normalized会fallback到layers['counts']
- ✓ 对新数据：preprocessing会创建正确的.raw
- ✓ 不影响其他分析工具

## 下一步

**请重启MCP使修改生效，然后测试完整的deconvolution流程！**

测试建议:
1. 加载两个空间转录组数据集
2. 对两个数据集使用不同的HVG参数(确保基因子集不重叠)
3. 运行deconvolution (如Cell2location)
4. 验证能找到足够的common genes并成功运行
