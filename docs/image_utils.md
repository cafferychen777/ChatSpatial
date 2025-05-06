# 图像处理模块

空间转录组MCP提供了一个标准化的图像处理模块，用于处理和转换可视化结果。该模块位于`spatial_transcriptomics_mcp/utils/image_utils.py`，提供了一系列函数，确保所有可视化函数都返回标准化的Image对象。

## 主要功能

### fig_to_image

将matplotlib图形转换为Image对象，适用于MCP框架。

```python
from spatial_transcriptomics_mcp.utils.image_utils import fig_to_image
import matplotlib.pyplot as plt

# 创建matplotlib图形
fig, ax = plt.subplots()
ax.plot([1, 2, 3, 4], [1, 4, 9, 16])
ax.set_title("Sample Plot")

# 转换为Image对象
image = fig_to_image(fig, dpi=100, format='png', max_size_kb=500, close_fig=True)

# 现在可以直接返回image对象给MCP框架
```

参数说明：
- `fig`: Matplotlib图形对象
- `dpi`: 图像分辨率（点/英寸），默认为100
- `format`: 图像格式，默认为'png'
- `max_size_kb`: 图像最大大小（KB），默认为500
- `close_fig`: 是否在转换后关闭图形，默认为True

### create_placeholder_image

创建带有消息的占位图像，当无法生成实际可视化时使用。

```python
from spatial_transcriptomics_mcp.utils.image_utils import create_placeholder_image

# 创建占位图像
image = create_placeholder_image(
    message="No visualization available", 
    figsize=(6, 6), 
    format='png'
)

# 直接返回image对象给MCP框架
```

参数说明：
- `message`: 要显示在占位图像中的消息，默认为"No visualization available"
- `figsize`: 图像大小（英寸），默认为(6, 6)
- `format`: 图像格式，默认为'png'

### fig_to_base64

将matplotlib图形转换为base64编码的字符串（向后兼容）。

```python
from spatial_transcriptomics_mcp.utils.image_utils import fig_to_base64
import matplotlib.pyplot as plt

# 创建matplotlib图形
fig, ax = plt.subplots()
ax.plot([1, 2, 3, 4], [1, 4, 9, 16])
ax.set_title("Sample Plot")

# 转换为base64字符串
base64_str = fig_to_base64(fig, dpi=100, format='png', max_size_mb=5, close_fig=True)
```

参数说明：
- `fig`: Matplotlib图形对象
- `dpi`: 图像分辨率（点/英寸），默认为100
- `format`: 图像格式，默认为'png'
- `max_size_mb`: 图像最大大小（MB），默认为5
- `close_fig`: 是否在转换后关闭图形，默认为True

## 使用建议

1. **优先使用fig_to_image**：所有新代码应该使用`fig_to_image`函数，而不是`fig_to_base64`，以确保返回标准化的Image对象。

2. **处理大图像**：对于复杂的可视化，可能会生成较大的图像。`fig_to_image`函数会自动尝试压缩图像以满足大小限制。

3. **错误处理**：当图像生成失败时，可以使用`create_placeholder_image`函数创建占位图像，而不是返回None或抛出异常。

4. **在MCP工具中的使用**：

```python
from mcp.server.fastmcp import Context
from mcp.server.fastmcp.utilities.types import Image
from spatial_transcriptomics_mcp.utils.image_utils import fig_to_image

async def my_visualization_tool(data_id: str, data_store: dict, params: dict, context: Context) -> Image:
    """示例可视化工具函数"""
    # 获取数据
    data = data_store.get(data_id)
    if data is None:
        await context.warning(f"Dataset {data_id} not found")
        return create_placeholder_image("Dataset not found")
    
    # 创建可视化
    fig, ax = plt.subplots()
    # ... 绘图代码 ...
    
    # 返回Image对象
    return fig_to_image(fig)
```

## 测试

图像处理模块包含完整的单元测试，可以通过以下命令运行：

```bash
python -m tests.test_image_utils
```

集成测试也可用于测试图像处理模块与MCP工具的集成：

```bash
python -m tests.test_visualization_integration
```
