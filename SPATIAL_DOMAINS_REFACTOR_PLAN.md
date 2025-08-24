# Spatial Domains 架构重构方案

## 🚨 **发现的严重架构违规**

在 `spatial_domains.py` 中发现了 **15 个严重的预处理违规**，违背了核心架构原则："LLM会自动做preprocess的，它们可以做preprocessing.py的"

## 📊 **违规分类统计**

### 1. **标准化和Log变换违规** (4处)
```python
# 违规代码位置: 行148-149, 161-162
sc.pp.normalize_total(adata_subset, target_sum=1e4)
sc.pp.log1p(adata_subset)
```
**违规场景**: SpaGCN算法中对原始数据和raw数据进行标准化

### 2. **下采样违规** (1处)
```python
# 违规代码位置: 行195
sc.pp.subsample(adata_subset, n_obs=3000, random_state=42)
```
**违规场景**: SpaGCN性能优化时对大数据集下采样

### 3. **HVG选择违规** (2处)
```python
# 违规代码位置: 行207, 212
sc.pp.highly_variable_genes(adata_subset, n_top_genes=max_genes)
```
**违规场景**: SpaGCN内存优化时选择高变异基因

### 4. **PCA计算违规** (2处)
```python
# 违规代码位置: 行508, 513
sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
sc.tl.pca(adata, use_highly_variable=False)
```
**违规场景**: Leiden/Louvain聚类前的降维

### 5. **邻域图计算违规** (3处)
```python
# 违规代码位置: 行525, 531, 728
sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca')
sc.pp.neighbors(adata_stagate, use_rep='STAGATE', n_neighbors=params.cluster_n_neighbors or 15)
```
**违规场景**: 聚类算法前的邻域图构建

### 6. **聚类计算违规** (3处)
```python
# 违规代码位置: 行594, 597, 602, 754
sc.tl.leiden(adata, resolution=params.resolution, key_added=key_added)
sc.tl.louvain(adata, resolution=params.resolution, key_added=key_added)
sc.tl.leiden(adata_stagate, resolution=params.cluster_resolution or 1.0)
```
**违规场景**: 多种聚类算法的实际执行

## 🔧 **重构方案设计**

### **核心原则**
1. **严格遵循架构**: 移除所有预处理，期望接收预处理好的数据
2. **智能错误提示**: 为每种算法提供具体的预处理要求指导
3. **LLM友好消息**: 明确告诉LLM每个算法需要什么样的预处理

### **方案1: SpaGCN重构**

**当前违规**: 执行标准化、log变换、下采样、HVG选择

**重构后**:
```python
async def _identify_domains_spagcn(adata, params, context):
    # 1. 验证数据预处理状态
    max_val = adata.X.max() if hasattr(adata.X, 'max') else np.max(adata.X)
    min_val = adata.X.min() if hasattr(adata.X, 'min') else np.min(adata.X)
    
    # 检查是否为原始计数数据
    if min_val >= 0 and max_val > 100:
        raise ValueError(
            "SpaGCN requires preprocessed data but raw counts detected. "
            "Required preprocessing: "
            "1) Count normalization: sc.pp.normalize_total(adata, target_sum=1e4) "
            "2) Log transformation: sc.pp.log1p(adata) "
            "3) Optional: HVG selection for large datasets "
            "Please use preprocessing.py with SpaGCN-specific parameters."
        )
    
    # 检查scaled数据问题
    if min_val < -1:  # 可能是scaled数据
        if not (hasattr(adata, 'raw') and adata.raw is not None):
            raise ValueError(
                "SpaGCN detected scaled data but raw data unavailable. "
                "SpaGCN requires normalized (not scaled) data. "
                "Please ensure preprocessing preserves raw data or run SpaGCN before scaling."
            )
        # 使用raw数据但验证其预处理状态
        raw_max = adata.raw.X.max() if hasattr(adata.raw.X, 'max') else np.max(adata.raw.X)
        if raw_max > 100:
            raise ValueError(
                "Raw data appears unnormalized. SpaGCN requires: "
                "sc.pp.normalize_total(adata.raw, target_sum=1e4); sc.pp.log1p(adata.raw)"
            )
        # 使用预处理好的raw数据
        adata_subset = adata.raw.to_adata()
    
    # 检查数据集大小
    if adata.n_obs > 3000:
        await context.warning(
            f"Large dataset ({adata.n_obs} spots) may cause SpaGCN performance issues. "
            "Consider subsampling in preprocessing.py: sc.pp.subsample(adata, n_obs=3000)"
        )
    
    # 检查基因数量
    if adata.n_vars > 2000:
        if 'highly_variable' not in adata.var:
            raise ValueError(
                f"Large gene set ({adata.n_vars} genes) requires HVG selection. "
                "Please run HVG selection in preprocessing.py: "
                "sc.pp.highly_variable_genes(adata, n_top_genes=2000)"
            )
        # 使用预选的HVG
        adata_subset = adata[:, adata.var.highly_variable].copy()
    
    # 继续SpaGCN核心算法...
```

### **方案2: Leiden/Louvain重构**

**当前违规**: PCA计算、邻域图构建、聚类执行

**重构后**:
```python
async def _identify_domains_clustering(adata, params, context):
    # 验证PCA
    if 'X_pca' not in adata.obsm:
        raise ValueError(
            f"{params.method} clustering requires PCA. "
            "Please run PCA in preprocessing.py: "
            "sc.tl.pca(adata, n_comps=50)"
        )
    
    # 验证邻域图
    if 'neighbors' not in adata.uns:
        raise ValueError(
            f"{params.method} clustering requires neighborhood graph. "
            "Please compute neighbors in preprocessing.py: "
            "sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')"
        )
    
    # 验证是否已有聚类结果
    cluster_key = f"spatial_{params.method}"
    if cluster_key not in adata.obs:
        raise ValueError(
            f"{params.method} clustering not found. "
            "Please run clustering in preprocessing.py: "
            f"sc.tl.{params.method}(adata, resolution={params.resolution}, key_added='{cluster_key}')"
        )
    
    # 使用预计算的聚类结果
    domain_labels = adata.obs[cluster_key].astype(str)
    # 返回结果...
```

### **方案3: STAGATE重构**

**当前违规**: 邻域图构建、聚类执行

**重构后**:
```python
async def _identify_domains_stagate(adata, params, context):
    # 验证STAGATE embeddings
    if 'STAGATE' not in adata.obsm:
        raise ValueError(
            "STAGATE embeddings not found. "
            "Please run STAGATE preprocessing: "
            "1) Install STAGATE: pip install stagate "
            "2) Run in preprocessing.py: STAGATE.train_STAGATE(adata) "
            "3) Generate embeddings: STAGATE.cal_STAGATE(adata)"
        )
    
    # 验证聚类结果
    if 'leiden' not in adata.obs:
        raise ValueError(
            "STAGATE clustering not found. "
            "Please run clustering on STAGATE embeddings in preprocessing.py: "
            "sc.pp.neighbors(adata, use_rep='STAGATE'); "
            "sc.tl.leiden(adata, resolution=1.0)"
        )
    
    # 使用预计算结果
    domain_labels = adata.obs['leiden'].astype(str)
```

## 🎯 **LLM指导消息设计**

每个错误消息都包含：
1. **问题识别**: 明确说明缺失什么
2. **解决方案**: 具体的代码示例
3. **预处理位置**: 明确指向preprocessing.py
4. **算法特殊要求**: 每个算法的特定需求

### **消息模板**:
```python
raise ValueError(
    f"{algorithm} requires {requirement} but {current_state} detected. "
    f"Required preprocessing in preprocessing.py: "
    f"{specific_code_example} "
    f"Algorithm-specific note: {special_requirements}"
)
```

## 📋 **实施计划**

### **Phase 1: SpaGCN重构** (最复杂)
- 移除4处标准化违规
- 移除下采样违规  
- 移除HVG违规
- 添加智能数据验证

### **Phase 2: Clustering重构** 
- 移除PCA违规
- 移除neighbors违规
- 移除leiden/louvain违规
- 期望预计算结果

### **Phase 3: STAGATE/BANKSY重构**
- 移除neighbors违规
- 移除聚类违规
- 期望预计算embeddings

### **Phase 4: 测试验证**
- 创建测试用例验证错误消息
- 确保LLM能理解预处理要求
- 验证算法功能完整性

## ⚠️ **风险评估**

1. **算法破坏风险**: 某些算法可能有硬编码的预处理依赖
2. **用户体验**: 错误消息必须足够清晰指导用户
3. **向后兼容**: 需要考虑现有使用者
4. **性能影响**: 验证逻辑不应显著影响性能

## 🎉 **预期收益**

1. **架构一致性**: 完全符合"LLM自动预处理"原则
2. **责任分离**: preprocessing.py专门处理预处理
3. **可维护性**: 算法逻辑和预处理逻辑分离
4. **LLM友好**: 清晰的错误消息指导自动预处理

## 📝 **实施确认**

请确认此方案是否符合架构要求：
- ✅ 严格移除所有预处理逻辑
- ✅ 添加智能验证和错误指导  
- ✅ 为LLM提供清晰的预处理要求
- ✅ 保持算法功能完整性

**方案确认后开始实施重构。**