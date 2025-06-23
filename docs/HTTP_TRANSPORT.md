# HTTP Transport for ChatSpatial MCP Server

ChatSpatial 现在支持通过 HTTP 传输协议提供 MCP 服务，这使得 Web 应用和其他 HTTP 客户端可以访问空间转录组分析功能。

## 功能特性

### 1. 多种传输方式
- **stdio**: 标准输入/输出（默认）
- **sse**: Server-Sent Events（FastMCP 内置）
- **http**: 完整的 HTTP/REST API（新增）

### 2. HTTP 服务器特性
- RESTful API 端点
- Server-Sent Events (SSE) 流式传输
- 会话管理
- CORS 支持（默认仅限 localhost）
- 安全中间件（Origin 验证、安全头）
- 速率限制（每分钟 100 请求）

### 3. 安全特性
- 默认绑定到 127.0.0.1（localhost）
- Origin 验证防止 DNS 重绑定攻击
- 安全响应头（X-Content-Type-Options, X-Frame-Options 等）
- 会话隔离

## 快速开始

### 1. 启动 HTTP 服务器

```bash
# 使用默认设置（localhost:8000）
chatspatial-http

# 自定义端口
chatspatial-http --port 8080

# 开发模式（自动重载）
chatspatial-http --reload

# 允许外部连接（不推荐）
chatspatial-http --host 0.0.0.0 --allow-external
```

### 2. API 端点

#### 基础信息
```bash
GET http://localhost:8000/
```

#### 健康检查
```bash
GET http://localhost:8000/health
```

#### 会话管理
```bash
# 创建会话
POST http://localhost:8000/sessions

# 获取会话信息
GET http://localhost:8000/sessions/{session_id}

# 删除会话
DELETE http://localhost:8000/sessions/{session_id}
```

#### RPC 调用
```bash
POST http://localhost:8000/rpc
Content-Type: application/json

{
  "jsonrpc": "2.0",
  "method": "tools/list",
  "id": "1"
}
```

#### SSE 流
```bash
GET http://localhost:8000/sse
```

## Python 客户端示例

```python
import asyncio
from chatspatial.examples.http_client_example import ChatSpatialHTTPClient

async def main():
    # 创建客户端
    client = ChatSpatialHTTPClient("http://localhost:8000")
    
    # 创建会话
    session_id = await client.create_session()
    print(f"Session: {session_id}")
    
    # 列出工具
    tools = await client.list_tools()
    for tool in tools:
        print(f"Tool: {tool['name']}")
    
    # 调用工具
    result = await client.call_tool(
        "load_data",
        data_path="/path/to/data.h5ad",
        data_type="10x_visium"
    )
    print(f"Result: {result}")

asyncio.run(main())
```

## JavaScript/TypeScript 客户端示例

```javascript
// 创建会话
const sessionResp = await fetch('http://localhost:8000/sessions', {
  method: 'POST'
});
const { session_id } = await sessionResp.json();

// RPC 调用
const rpcResp = await fetch('http://localhost:8000/rpc', {
  method: 'POST',
  headers: {
    'Content-Type': 'application/json',
    'X-Session-Id': session_id
  },
  body: JSON.stringify({
    jsonrpc: '2.0',
    method: 'tools/call',
    params: {
      name: 'load_data',
      arguments: {
        data_path: '/path/to/data.h5ad',
        data_type: '10x_visium'
      }
    },
    id: '1'
  })
});

const result = await rpcResp.json();
console.log(result);
```

## SSE 流式连接

```javascript
const eventSource = new EventSource('http://localhost:8000/sse');

eventSource.onmessage = (event) => {
  const data = JSON.parse(event.data);
  console.log('Received:', data);
};

eventSource.onerror = (error) => {
  console.error('SSE Error:', error);
};
```

## 安全建议

1. **生产环境**：始终使用 HTTPS 和适当的认证机制
2. **防火墙**：确保服务器端口仅对授权客户端开放
3. **速率限制**：根据需要调整速率限制设置
4. **会话管理**：定期清理过期会话

## 与 MCP 客户端集成

HTTP 传输使得 ChatSpatial 可以与更多类型的客户端集成：

1. **Web 应用**：通过 JavaScript/TypeScript 直接访问
2. **移动应用**：通过 HTTP API 访问
3. **其他语言**：任何支持 HTTP 的语言都可以访问

## 故障排除

### 端口已被占用
```bash
# 使用其他端口
chatspatial-http --port 8081
```

### CORS 错误
确保客户端运行在 localhost，或者配置适当的 CORS 设置。

### 连接被拒绝
检查防火墙设置，确保端口未被阻止。

## 性能优化

1. **使用会话**：避免每次请求都重新加载数据
2. **批量操作**：尽可能合并多个操作
3. **SSE 流**：对于长时间运行的操作，使用 SSE 获取进度更新

## 总结

HTTP 传输支持使 ChatSpatial 更加灵活和易于集成。它保持了与 MCP 规范的兼容性，同时提供了现代 Web 应用所需的功能。