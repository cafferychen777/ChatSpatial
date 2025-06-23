# MCP 工具错误处理迁移指南

根据 MCP 规范，工具错误应该在结果对象中报告，而不是作为协议级错误抛出。这允许 LLM 看到错误并可能采取纠正措施。

## 为什么要改进？

当前实现的问题：
- 工具抛出异常时，LLM 无法看到具体错误信息
- 错误处理不符合 MCP 规范
- 缺少 `isError` 标志，客户端无法正确识别错误

## 迁移方案

### 方案 1：使用装饰器（推荐）

最简单的方式是使用 `@mcp_tool_error_handler` 装饰器：

```python
from chatspatial.utils.tool_error_handling import mcp_tool_error_handler

@mcp.tool()
@mcp_tool_error_handler()  # 添加这个装饰器
async def preprocess_data(
    data_id: str,
    params: AnalysisParameters = AnalysisParameters(),
    context: Context = None
):
    # 原有代码保持不变
    # 装饰器会自动处理异常并返回正确格式
    ...
```

### 方案 2：手动改造

对于需要更精细控制的情况：

```python
from chatspatial.utils.tool_error_handling import (
    create_success_result, 
    create_error_result,
    dataset_not_found_error
)

@mcp.tool()
async def preprocess_data(
    data_id: str,
    params: AnalysisParameters = AnalysisParameters(),
    context: Context = None
) -> Dict[str, Any]:  # 注意返回类型改为 Dict
    # 验证数据集
    if data_id not in data_store:
        return dataset_not_found_error(data_id).to_dict()
    
    try:
        # 原有的处理逻辑
        adata = data_store[data_id]["adata"]
        
        # ... 处理过程 ...
        
        result = PreprocessingResult(
            data_id=data_id,
            n_cells=adata.n_obs,
            n_genes=adata.n_vars,
            # ...
        )
        
        # 成功时返回
        return create_success_result(result).to_dict()
        
    except ValueError as e:
        # 参数错误
        return create_error_result(e, include_traceback=False).to_dict()
        
    except Exception as e:
        # 其他错误
        if context:
            await context.warning(f"Processing failed: {str(e)}")
        return create_error_result(e, include_traceback=True).to_dict()
```

## 具体迁移步骤

### 1. 更新 server.py

在 server.py 中导入错误处理工具：

```python
from .utils.tool_error_handling import mcp_tool_error_handler
```

### 2. 更新每个工具

为每个工具添加装饰器：

```python
# 之前
@mcp.tool()
@mcp_error_handler  # 旧的错误处理
async def load_data(...):
    ...

# 之后
@mcp.tool()
@mcp_tool_error_handler()  # 新的 MCP 兼容错误处理
async def load_data(...):
    ...
```

### 3. 处理特殊情况

某些工具可能需要特殊处理：

```python
@mcp.tool()
async def visualize_data(
    data_id: str,
    params: VisualizationParameters = VisualizationParameters(),
    context: Context = None
) -> Dict[str, Any]:
    if data_id not in data_store:
        return dataset_not_found_error(data_id).to_dict()
    
    try:
        # 生成图像
        fig = create_figure(...)
        
        # 返回图像结果
        return {
            "content": [{
                "type": "image",
                "data": fig_to_base64(fig),
                "mimeType": "image/png"
            }],
            "isError": False
        }
        
    except Exception as e:
        # 可视化错误时返回文本说明
        return create_error_result(e).to_dict()
```

## 返回格式示例

### 成功返回：

```json
{
    "content": [{
        "type": "text",
        "text": "{\"data_id\": \"data_1\", \"n_cells\": 5000, \"n_genes\": 20000}"
    }],
    "isError": false
}
```

### 错误返回：

```json
{
    "content": [{
        "type": "text",
        "text": "Error: Dataset 'data_99' not found. Please load a dataset first using the 'load_data' tool."
    }],
    "isError": true
}
```

### 带追踪信息的错误：

```json
{
    "content": [
        {
            "type": "text",
            "text": "Error: Matrix multiplication failed: shapes (100,50) and (60,30) not aligned"
        },
        {
            "type": "text",
            "text": "Traceback:\n  File \"preprocessing.py\", line 123, in preprocess_data\n    ..."
        }
    ],
    "isError": true
}
```

## 测试迁移

创建测试脚本验证错误处理：

```python
# test_error_handling.py
async def test_tool_errors():
    # 测试数据集不存在
    result = await preprocess_data("invalid_id")
    assert result["isError"] == True
    assert "not found" in result["content"][0]["text"]
    
    # 测试处理错误
    result = await preprocess_data("data_1", 
        params=AnalysisParameters(n_pcs=99999))  # 无效参数
    assert result["isError"] == True
    
    print("✅ 错误处理测试通过")
```

## 注意事项

1. **不要混用**：要么使用装饰器，要么手动处理，不要同时使用
2. **保持一致性**：所有工具都应该使用相同的错误处理方式
3. **日志记录**：仍然可以使用 `context.warning()` 记录错误，但主要错误信息应该在返回值中
4. **向后兼容**：装饰器会自动处理旧格式的返回值

## 迁移优先级

1. **高优先级**：经常出错的工具
   - `load_data` - 文件不存在
   - `preprocess_data` - 参数错误
   - `visualize_data` - 绘图错误

2. **中优先级**：计算密集型工具
   - `analyze_spatial_data`
   - `deconvolve_data`
   - `find_spatial_genes`

3. **低优先级**：简单查询工具
   - `find_markers`
   - 其他稳定工具

## 完成标志

迁移完成后，所有工具应该：
- ✅ 永不抛出未捕获的异常
- ✅ 错误时返回包含 `isError: true` 的结果
- ✅ 提供清晰的错误信息给 LLM
- ✅ 保留必要的追踪信息用于调试