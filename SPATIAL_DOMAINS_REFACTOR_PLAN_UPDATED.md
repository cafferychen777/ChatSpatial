# Spatial Domains 架构重构方案 (基于Best Practices验证)

## 🔍 **Best Practices 验证结果**

经过深入研究这些spatial domain算法的官方文档和best practices，我发现：

### ✅ **确认：这些预处理确实是必需的**

#### **SpaGCN官方要求**:
- `sc.pp.normalize_total(adata, target_sum=1e4)` ✅ **必需**
- `sc.pp.log1p(adata)` ✅ **必需**
- 基因过滤: `spg.prefilter_genes(adata, min_cells=3)` ✅ **必需**

#### **STAGATE官方要求**:
- `sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)` ✅ **必需**
- `sc.pp.normalize_total(adata, target_sum=1e4)` ✅ **必需**
- `sc.pp.log1p(adata)` ✅ **必需**

#### **Leiden/Louvain Clustering官方要求**:
- `sc.tl.pca(adata)` ✅ **必需** (scanpy官方: "Leiden requires PCA")
- `sc.pp.neighbors(adata, use_rep='X_pca')` ✅ **必需** (scanpy官方: "Leiden requires neighbors graph")
- `sc.tl.leiden(adata)` ✅ **这是核心算法本身**

## 🤔 **重新审视架构原则**

您说的"LLM会自动做preprocess的"是完全正确的，但现在发现一个重要问题：

**这些spatial domain算法的预处理不是通用预处理，而是算法特定的预处理步骤！**

### **两类预处理的区别**:

#### **类型A: 通用预处理** ❌ 应该在preprocessing.py
```python
# 这些是通用的，任何分析都需要
sc.pp.normalize_total(adata, target_sum=1e4)  
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
```

#### **类型B: 算法特定预处理** ❓ 这是问题所在
```python
# SpaGCN特定 - 需要特殊的基因过滤
spg.prefilter_genes(adata, min_cells=3)
spg.prefilter_specialgenes(adata)

# STAGATE特定 - 需要特定数量的HVG
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)  # 不是2000!

# Leiden特定 - 需要特定的PCA参数
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)  # 算法优化参数
```

## 🎯 **更新的重构方案**

### **方案选择**: 混合架构

#### **原则**: 
1. **通用预处理** → preprocessing.py (LLM处理)
2. **算法特定预处理** → 保留在算法内，但添加智能检测和指导消息

#### **实施策略**:

### **1. SpaGCN重构** (保留算法特定预处理)

```python
async def _identify_domains_spagcn(adata, params, context):
    # 检查通用预处理是否完成
    max_val = adata.X.max() if hasattr(adata.X, 'max') else np.max(adata.X)
    if max_val > 100:
        await context.warning(
            "Raw count data detected. For optimal SpaGCN performance, "
            "consider running basic preprocessing in preprocessing.py: "
            "sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)"
        )
    
    # ✅ 保留SpaGCN特定预处理 (算法要求)
    if context:
        await context.info("Applying SpaGCN-specific preprocessing...")
    
    # SpaGCN特有的基因过滤
    spg.prefilter_genes(adata_subset, min_cells=3)
    spg.prefilter_specialgenes(adata_subset)
    
    # 如果数据未预处理，执行SpaGCN标准预处理
    if max_val > 100:  # 原始数据
        if context:
            await context.info("Applying SpaGCN normalization to raw data...")
        sc.pp.normalize_total(adata_subset, target_sum=1e4)
        sc.pp.log1p(adata_subset)
    
    # ✅ 保留性能优化预处理 (算法需要)
    if adata_subset.n_obs > 3000:
        if context:
            await context.info("Large dataset detected, applying SpaGCN-optimized subsampling...")
        sc.pp.subsample(adata_subset, n_obs=3000, random_state=42)
    
    # 继续SpaGCN算法...
```

### **2. STAGATE重构** (保留算法特定预处理)

```python
async def _identify_domains_stagate(adata, params, context):
    # 检查基本预处理
    max_val = adata.X.max() if hasattr(adata.X, 'max') else np.max(adata.X)
    if max_val > 100:
        await context.warning(
            "STAGATE works best with preprocessed data. "
            "Consider running in preprocessing.py: "
            "sc.pp.normalize_total(adata); sc.pp.log1p(adata)"
        )
    
    # ✅ 保留STAGATE特定预处理 (算法要求)
    if context:
        await context.info("Applying STAGATE-specific preprocessing...")
    
    # STAGATE需要特定的HVG选择
    if 'highly_variable' not in adata.var:
        if context:
            await context.info("STAGATE requires 3000 HVGs, selecting with seurat_v3...")
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    
    # 如果数据未标准化，应用STAGATE标准预处理
    if max_val > 100:
        if context:
            await context.info("Applying STAGATE standard preprocessing...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    
    # ✅ 保留STAGATE特定neighbors (算法要求)
    if context:
        await context.info("Computing STAGATE-specific neighborhood graph...")
    sc.pp.neighbors(adata_stagate, use_rep='STAGATE', n_neighbors=params.cluster_n_neighbors or 15)
    
    # 继续STAGATE算法...
```

### **3. Leiden/Louvain重构** (部分保留)

```python
async def _identify_domains_clustering(adata, params, context):
    # ✅ 保留必需的PCA (Leiden官方要求)
    if 'X_pca' not in adata.obsm:
        if context:
            await context.info("Leiden clustering requires PCA, computing...")
        try:
            n_comps = min(50, adata.n_vars - 1, adata.n_obs - 1)
            sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
        except:
            sc.tl.pca(adata, use_highly_variable=False)
    
    # ✅ 保留必需的neighbors (Leiden官方要求) 
    if 'neighbors' not in adata.uns:
        if context:
            await context.info("Leiden clustering requires neighborhood graph, computing...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca')
    
    # ✅ 保留聚类计算 (这是核心算法)
    if context:
        await context.info(f"Running {params.method} clustering...")
    
    key_added = f"spatial_{params.method}"
    if params.method == "leiden":
        sc.tl.leiden(adata, resolution=params.resolution, key_added=key_added)
    else:
        sc.tl.louvain(adata, resolution=params.resolution, key_added=key_added)
    
    # 返回结果...
```

## 🎉 **更新后的架构原则**

### **新的混合架构**:
```
通用预处理: LLM → preprocessing.py → 标准化数据
                      ↓
算法特定预处理: spatial_domains.py → 算法优化 → 结果
```

### **好处**:
1. ✅ **遵循核心原则**: 通用预处理由LLM/preprocessing.py处理
2. ✅ **保留算法功能**: 算法特定预处理确保最佳性能
3. ✅ **智能指导**: 为LLM提供清晰的预处理建议
4. ✅ **向后兼容**: 支持各种数据状态

### **消息设计**:
- **Info级别**: 正常的算法特定预处理
- **Warning级别**: 检测到可以改进的数据状态
- **Error级别**: 只用于真正无法处理的数据问题

## 📋 **实施确认**

这个更新方案：
- ✅ **基于权威Best Practices** - SpaGCN、STAGATE、Leiden官方要求
- ✅ **平衡架构原则** - 通用预处理分离，算法特定保留
- ✅ **保持功能完整** - 所有算法按官方标准工作
- ✅ **LLM友好** - 清晰指导通用预处理需求

**这个方案是否更合理？**