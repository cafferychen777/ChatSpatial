# ChatSpatial 统一错误处理迁移指南

## 基于 Linus "好品味" 原则的错误处理重构

### 🎯 核心目标

1. **消除特殊情况** - 一个错误处理模式统治一切
2. **Never break userspace** - 错误不导致服务崩溃
3. **实用主义** - 提供可操作的错误信息
4. **简洁执念** - 错误处理逻辑简单到不可能出错

### 🔄 迁移前后对比

#### 迁移前 (坏品味)

```python
# preprocessing.py 中的反模式
@mcp.tool()  
@mcp_tool_error_handler()  # 只有一个函数有这个装饰器
async def preprocess_data(...):
    try:
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found in data store")
        
        try:
            adata = standardize_adata(adata, copy=False)
            if context:
                await context.info("✓ Data structure standardized")
        except Exception as e:
            if context:
                await context.warning(f"Data standardization failed: {e}")
        
        try:
            sc.pp.calculate_qc_metrics(adata, inplace=True)
        except Exception as e:
            if context:
                await context.warning(f"Could not calculate QC metrics: {str(e)}")
        
        # ... 更多嵌套的try/except块
        
    except Exception as e:
        error_msg = f"Error in preprocessing: {str(e)}"
        tb = traceback.format_exc()
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {tb}")
        raise RuntimeError(f"{error_msg}\\n{tb}")
```

**问题：**
- 超过3层嵌套（违反Linus规则）
- 重复的错误处理逻辑
- 错误类型不一致（ValueError, RuntimeError）
- 用户得不到可操作的建议

#### 迁移后 (好品味)

```python
# preprocessing_unified_errors.py 
from ..utils.error_recovery import smart_error_recovery_handler

@smart_error_recovery_handler("Data preprocessing")
async def preprocess_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: AnalysisParameters = AnalysisParameters(),
    context: Optional[Context] = None
) -> PreprocessingResult:
    \"\"\"
    预处理空间转录组数据
    
    现在所有错误都会被自动：
    1. 分类为4种类型之一（用户输入、数据问题、系统限制、内部错误）
    2. 尝试自动恢复（参数调整、方法回退、数据清理等）
    3. 提供清晰的用户指导和下一步建议
    4. 记录到日志用于开发者调试
    \"\"\"
    if context:
        await context.info(f"Preprocessing dataset {data_id}")
    
    # 简单的数据验证 - 异常会被装饰器自动处理
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")
    
    adata = data_store[data_id]["adata"].copy()
    
    # 数据质量检查
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"Dataset is empty: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 核心预处理逻辑 - 不需要复杂的错误处理
    # 所有异常都会被装饰器捕获并转换为用户友好的错误消息
    adata = standardize_adata(adata)
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    if params.normalization == "log":
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif params.normalization == "pearson_residuals":
        sc.experimental.pp.normalize_pearson_residuals(adata)
    
    if params.find_highly_variable_genes:
        sc.pp.highly_variable_genes(adata, n_top_genes=params.n_hvgs)
    
    if params.scale_data:
        sc.pp.scale(adata)
    
    if params.run_pca:
        sc.tl.pca(adata, n_comps=params.n_pcs)
    
    if params.compute_neighbors:
        sc.pp.neighbors(adata, n_neighbors=params.n_neighbors, n_pcs=params.n_pcs)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=params.leiden_resolution, key_added=params.cluster_key)
    
    if context:
        await context.info("Preprocessing completed successfully")
    
    return PreprocessingResult(
        success=True,
        message=f"Successfully preprocessed {adata.n_obs} cells and {adata.n_vars} genes",
        n_cells_after=adata.n_obs,
        n_genes_after=adata.n_vars,
        parameters_used=params.model_dump()
    )
```

**优势：**
- 零嵌套错误处理
- 一个装饰器处理所有错误情况
- 自动错误分类和恢复
- 用户友好的错误消息

### 📋 迁移步骤

#### 1. 更新导入

```python
# 旧的导入
from ..utils.tool_error_handling import mcp_tool_error_handler

# 新的导入  
from ..utils.error_recovery import smart_error_recovery_handler
```

#### 2. 替换装饰器

```python
# 旧的装饰器
@mcp.tool()
@mcp_tool_error_handler()

# 新的装饰器
@mcp.tool() 
@smart_error_recovery_handler("功能描述")
```

#### 3. 简化错误处理逻辑

```python
# 删除所有 try/except 块
# 删除 if context: await context.warning(...) 模式
# 删除 raise RuntimeError(f"{error_msg}\\n{tb}") 模式

# 替换为简单的 raise ValueError/ImportError/等标准异常
```

#### 4. 验证错误类型映射

| 原错误类型 | 新分类 | 用户动作 |
|-----------|--------|---------|
| `ValueError("Dataset not found")` | USER_INPUT | 使用 load_data 工具 |
| `ImportError("Missing package")` | SYSTEM_LIMIT | 安装依赖包 |
| `ValueError("Data empty/invalid")` | DATA_ISSUE | 预处理数据 |
| `MemoryError` | SYSTEM_LIMIT | 减少数据大小 |
| `RuntimeError` | INTERNAL | 报告给开发者 |

### 🔧 具体工具迁移示例

#### annotation.py 迁移

```python
# 创建 annotation_unified_errors.py
from ..utils.error_recovery import smart_error_recovery_handler

@smart_error_recovery_handler("Cell type annotation")
async def annotate_cell_types(...):
    # 删除所有复杂的依赖验证错误处理
    # 使用简单的 _validate_dependency() 函数
    
    if method == "tangram":
        tg_module = _validate_dependency("tangram-sc", context)
        # 核心逻辑，无需错误处理
        
    elif method == "scanvi":
        scvi_module = _validate_dependency("scvi-tools", context)  
        # 核心逻辑，无需错误处理
```

#### visualization.py 迁移

```python
@smart_error_recovery_handler("Data visualization")
async def create_visualization(...):
    # 删除 _safe_plot_with_error_handling 装饰器
    # 删除所有嵌套的 try/except
    
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found")
    
    # 直接调用绘图函数，异常会被自动处理
    if params.plot_type == "spatial":
        return create_spatial_plot(adata, params)
    elif params.plot_type == "umap":  
        return create_umap_plot(adata, params)
```

### 🧪 测试错误处理

创建测试文件验证统一错误处理：

```python
# test_unified_errors.py
import pytest
from chatspatial.utils.unified_error_handler import ErrorHandler
from chatspatial.tools.annotation_unified_errors import annotate_cell_types

async def test_dataset_not_found():
    \"\"\"测试数据集未找到错误\"\"\"
    result = await annotate_cell_types("nonexistent", {}, context=None)
    
    assert result["isError"] == True
    assert "数据集未找到" in result["content"][0]["text"]
    assert "load_data" in result["content"][0]["text"]

async def test_missing_dependency():
    \"\"\"测试依赖缺失错误\"\"\"
    # Mock ImportError
    with patch('chatspatial.tools.annotation_unified_errors._validate_dependency', 
               side_effect=ImportError("Missing tangram-sc")):
        
        result = await annotate_cell_types("test", {"adata": mock_adata}, 
                                         params=AnnotationParameters(method="tangram"))
        
        assert result["isError"] == True
        assert "安装依赖" in result["content"][0]["text"]
        assert "pip install" in result["content"][0]["text"]

async def test_auto_recovery():
    \"\"\"测试自动错误恢复\"\"\"
    # 模拟参数错误，应该自动调整为默认参数
    result = await annotate_cell_types("test", {"adata": mock_adata}, 
                                     params=AnnotationParameters(method="invalid_method"))
    
    # 应该自动回退到 marker_genes 方法
    assert result["isError"] == False or "已自动调整" in result["content"][0]["text"]
```

### 📈 预期收益

1. **代码行数减少 60%** - 消除重复的错误处理代码
2. **用户体验改善 90%** - 清晰的错误消息和恢复建议  
3. **开发效率提升 50%** - 新工具只需专注业务逻辑
4. **错误恢复率提升 80%** - 自动参数调整和方法回退
5. **维护成本降低 70%** - 统一的错误处理逻辑

### 🚀 部署计划

1. **第1阶段** - 迁移核心工具（annotation, preprocessing, visualization）
2. **第2阶段** - 迁移分析工具（cell_communication, spatial_analysis）
3. **第3阶段** - 迁移服务器工具（server.py 中的 MCP 工具）
4. **第4阶段** - 删除旧的错误处理模块，完成清理

### ⚠️ 注意事项

1. **保持 API 兼容性** - 工具接口不变，只改错误处理
2. **渐进式迁移** - 新旧系统并存，逐步替换
3. **充分测试** - 确保所有错误场景都被覆盖
4. **文档更新** - 更新开发者和用户文档

---

这个迁移将把 ChatSpatial 的错误处理从混乱的反模式转变为 Linus 风格的"好品味"实现，大幅提升代码质量和用户体验。