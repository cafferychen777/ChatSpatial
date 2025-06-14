# Spatial Domains 模块问题诊断和修复

## 发现的主要问题

### 1. 依赖环境问题 ⚠️ **主要原因**
- **问题**: squidpy 导入失败，因为 dask.dataframe 的 legacy implementation 不再支持
- **影响**: 导致整个模块无法导入，这是"莫名其妙跑不出来"的根本原因
- **修复**: 添加了 squidpy 的导入保护和手动空间图构建作为后备方案

### 2. SpaGCN 超时问题 
- **问题**: SpaGCN 的 `detect_spatial_domains_ez_mode` 可能会无响应挂起，特别是在以下情况：
  - 基因数量过多 (>3000)
  - 参数设置不当
  - 数据包含 NaN/无穷值
- **修复**: 
  - 添加了5分钟超时保护
  - 使用线程池执行避免阻塞
  - 限制基因数量到3000以内
  - 清理数据中的问题值

### 3. 数据验证缺失
- **问题**: 缺少对输入数据的充分验证
- **修复**: 添加了多层验证：
  - 空间坐标验证
  - 数据形状匹配检查
  - NaN/无穷值处理
  - 退化空间坐标检测

### 4. 聚类方法健壮性
- **问题**: PCA和邻居图计算可能因数据大小不当而失败
- **修复**: 
  - 动态调整PCA组件数量
  - 根据数据大小调整邻居数量
  - 添加后备参数

### 5. 空间图构建依赖
- **问题**: 完全依赖 squidpy 进行空间图构建
- **修复**: 实现了手动空间图构建作为后备

## 具体修复内容

### 导入保护
```python
# 导入 squidpy 时使用 try-except
try:
    import squidpy as sq
    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
```

### SpaGCN 超时保护
```python
# 使用线程池和超时
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
    future = loop.run_in_executor(executor, spagcn_function)
    domain_labels = await asyncio.wait_for(future, timeout=300.0)
```

### 数据预处理增强
```python
# 限制基因数量
if adata_subset.n_vars > 3000:
    # 选择前3000个高变基因
    
# 清理问题数据
if np.any(np.isnan(adata_subset.X)) or np.any(np.isinf(adata_subset.X)):
    adata_subset.X = np.nan_to_num(adata_subset.X, nan=0.0, posinf=0.0, neginf=0.0)
```

### 手动空间图构建
```python
# 不依赖 squidpy 的空间图构建
coords = adata.obsm['spatial']
nbrs = NearestNeighbors(n_neighbors=6).fit(coords)
spatial_distances, spatial_indices = nbrs.kneighbors(coords)
# 构建稀疏连接矩阵...
```

## 建议的使用方式

1. **检查环境**: 确保依赖包版本兼容
2. **数据预处理**: 在调用前清理数据
3. **参数调优**: 
   - 对大数据集使用较少的 `n_domains`
   - 调整 SpaGCN 参数 (s, b, p)
   - 考虑使用 `use_highly_variable=True`

## 仍存在的限制

1. **环境依赖**: 如果 squidpy 完全不可用，空间图功能会降级
2. **性能**: 大数据集仍可能较慢
3. **SpaGCN**: 某些极端参数组合仍可能导致问题

## 测试建议

在使用修复后的代码时，建议：
1. 先用小数据集测试
2. 监控执行时间
3. 检查日志输出中的警告信息
4. 如果遇到超时，尝试调整参数或使用不同方法

## 文件位置

- 主文件: `/Users/apple/Research/SpatialTrans_MCP/chatspatial/chatspatial/tools/spatial_domains.py`
- 测试文件: `/Users/apple/Research/SpatialTrans_MCP/chatspatial/chatspatial/test_spatial_domains.py`
- 调试脚本: `/Users/apple/Research/SpatialTrans_MCP/chatspatial/debug_spatial_domains.py`

**修复的核心目标**: 让模块在各种环境和数据条件下都能稳定运行，避免神秘的挂起和崩溃。