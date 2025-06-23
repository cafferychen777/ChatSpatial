# ChatSpatial MCP 改进计划

基于 MCP 官方文档分析，我们的 ChatSpatial 可以在以下方面进行完善：

## 1. 🎯 核心功能增强

### 1.1 完整实现 MCP 核心概念
- [ ] **Resources**: 目前未实现资源系统
  - 添加空间数据资源暴露（如：`spatial://data/adata_file.h5ad`）
  - 支持资源订阅和更新通知
  - 实现资源模板，支持动态资源 URI

- [ ] **Prompts**: 目前未实现提示系统
  - 创建空间分析相关的提示模板
  - 支持参数化提示（如：分析类型、数据路径等）
  - 实现提示发现和动态更新

- [x] **Tools**: 已实现基础工具系统
  - 需要添加工具注解（annotations）来提供更好的 UX
  - 实现工具发现更新通知
  - 改进错误处理，使用 `isError` 标志

### 1.2 工具注解（Tool Annotations）
根据文档，应该为每个工具添加以下注解：
```typescript
{
  title: string;           // 人类可读的标题
  readOnlyHint: boolean;   // 是否只读操作
  destructiveHint: boolean;// 是否破坏性操作
  idempotentHint: boolean; // 是否幂等操作
  openWorldHint: boolean;  // 是否与外部交互
}
```

## 2. 🔄 传输层改进

### 2.1 支持多种传输方式
当前只支持 stdio，应该添加：
- [ ] **Streamable HTTP Transport**: 支持 Web 集成
  - 实现 HTTP POST 和 SSE 流
  - 支持会话管理
  - 实现断线重连

### 2.2 安全性增强
- [ ] 实现 Origin 验证（防止 DNS 重绑定攻击）
- [ ] 绑定到 localhost（127.0.0.1）而非 0.0.0.0
- [ ] 添加认证机制
- [ ] 实现速率限制

## 3. 📊 空间分析特定功能

### 3.1 Resources 实现建议
```python
# 空间数据资源
- spatial://datasets/{dataset_id}  # 数据集
- spatial://plots/{plot_id}        # 图表
- spatial://results/{analysis_id}  # 分析结果
- spatial://logs/{session_id}      # 分析日志
```

### 3.2 Prompts 实现建议
```python
# 空间分析提示模板
- "analyze-spatial-expression"     # 空间表达分析
- "find-cell-types"               # 细胞类型识别
- "compare-conditions"            # 条件比较
- "generate-visualization"        # 可视化生成
```

### 3.3 改进现有工具
为每个工具添加适当的注解，例如：
```python
"load_data": {
    "annotations": {
        "title": "Load Spatial Data",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": True  # 读取外部文件
    }
}
```

## 4. 🛠️ 开发者体验改进

### 4.1 日志和调试
- [ ] 实现结构化日志（使用 `send_log_message`）
- [ ] 添加性能指标记录
- [ ] 实现请求追踪（request IDs）

### 4.2 错误处理
- [ ] 统一错误格式
- [ ] 实现优雅的错误恢复
- [ ] 添加重试机制

### 4.3 文档和示例
- [ ] 创建详细的 API 文档
- [ ] 提供更多使用示例
- [ ] 添加集成测试

## 5. 🔌 客户端兼容性

### 5.1 确保与主流客户端兼容
根据文档，重点支持：
- Claude Desktop（已支持）
- VS Code GitHub Copilot
- Continue
- Cursor

### 5.2 功能支持矩阵
确保我们的实现支持：
- ✅ Tools（已支持）
- ❌ Resources（待实现）
- ❌ Prompts（待实现）
- ❌ Discovery（待实现）
- ❌ Sampling（可选）
- ❌ Roots（可选）

## 6. 🚀 实施优先级

### Phase 1（高优先级）
1. 实现 Resources 系统
2. 添加工具注解
3. 改进错误处理

### Phase 2（中优先级）
1. 实现 Prompts 系统
2. 添加 Streamable HTTP 传输
3. 增强安全性

### Phase 3（低优先级）
1. 实现 Sampling（如果需要）
2. 支持 Roots
3. 高级功能（如多服务器协调）

## 7. 📝 具体实施步骤

### 7.1 Resources 系统实现
```python
# 在 server.py 中添加
@app.list_resources()
async def list_resources() -> list[Resource]:
    """列出可用的空间数据资源"""
    return [
        Resource(
            uri="spatial://datasets/current",
            name="Current Dataset",
            mimeType="application/x-anndata",
            description="Currently loaded spatial dataset"
        )
    ]

@app.read_resource()
async def read_resource(uri: str) -> str:
    """读取资源内容"""
    if uri == "spatial://datasets/current":
        # 返回当前数据集的摘要信息
        return get_dataset_summary()
```

### 7.2 Prompts 系统实现
```python
@app.list_prompts()
async def list_prompts() -> list[Prompt]:
    """列出可用的分析提示"""
    return [
        Prompt(
            name="analyze-spatial-expression",
            description="Analyze spatial gene expression patterns",
            arguments=[
                PromptArgument(
                    name="genes",
                    description="Genes to analyze",
                    required=True
                ),
                PromptArgument(
                    name="method",
                    description="Analysis method",
                    required=False
                )
            ]
        )
    ]
```

### 7.3 工具注解添加
```python
TOOL_ANNOTATIONS = {
    "load_data": {
        "title": "Load Spatial Data",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": True
    },
    "spatial_analysis": {
        "title": "Run Spatial Analysis",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": False,
        "openWorldHint": False
    }
}
```

## 8. 🎯 预期成果

完成这些改进后，ChatSpatial 将：
1. 完全符合 MCP 规范
2. 提供更好的开发者体验
3. 支持更多客户端
4. 具有更强的安全性和可靠性
5. 成为空间转录组学分析的标准 MCP 实现