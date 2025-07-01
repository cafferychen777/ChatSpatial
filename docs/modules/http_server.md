# HTTP Server Module Documentation

## Overview

The `http_server.py` module provides an HTTP/SSE (Server-Sent Events) transport layer for the ChatSpatial MCP server. It wraps the core MCP server functionality in a FastAPI application, enabling web-based clients to interact with the spatial transcriptomics analysis tools through HTTP endpoints.

## Architecture

### Core Components

1. **FastAPI Application**
   - Modern async web framework
   - Automatic API documentation
   - Request/response validation
   - CORS middleware for cross-origin requests

2. **Session Management**
   - UUID-based session identification
   - Session-specific data stores
   - Activity tracking
   - Session lifecycle management

3. **Security Features**
   - CORS restricted to localhost by default
   - Origin validation middleware
   - Security headers (X-Frame-Options, X-XSS-Protection, etc.)
   - Rate limiting (100 requests/minute per IP)

4. **Transport Protocols**
   - HTTP POST for RPC calls
   - Server-Sent Events (SSE) for streaming
   - RESTful endpoints for session management

## API Endpoints

### General Endpoints

#### `GET /`
Returns server information and available endpoints.

**Response:**
```json
{
  "name": "ChatSpatial MCP Server",
  "version": "1.0.0",
  "protocol": "MCP",
  "transports": ["stdio", "sse", "http"],
  "endpoints": {
    "rpc": "/rpc",
    "sse": "/sse",
    "health": "/health",
    "sessions": "/sessions"
  }
}
```

#### `GET /health`
Health check endpoint for monitoring.

**Response:**
```json
{
  "status": "healthy",
  "sessions": 2,
  "datasets": 5
}
```

### Session Management

#### `POST /sessions`
Creates a new session with isolated data store.

**Response:**
```json
{
  "session_id": "uuid-string"
}
```
**Headers:**
- `X-Session-Id`: The created session ID

#### `GET /sessions/{session_id}`
Retrieves session information.

**Response:**
```json
{
  "id": "session-uuid",
  "created_at": 1234567890.0,
  "last_activity": 1234567900.0,
  "datasets": ["data_1", "data_2"]
}
```

#### `DELETE /sessions/{session_id}`
Deletes a session and its associated data.

**Response:**
```json
{
  "message": "Session deleted"
}
```

### MCP Protocol Endpoints

#### `POST /rpc`
Main RPC endpoint for MCP protocol communication.

**Request Body:**
```json
{
  "jsonrpc": "2.0",
  "method": "tools/call",
  "params": {
    "name": "load_data",
    "arguments": {
      "data_path": "/path/to/data.h5ad"
    }
  },
  "id": "request-id"
}
```

**Response:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    // Tool-specific result
  },
  "id": "request-id"
}
```

**Supported Methods:**
- `initialize`: Initialize MCP connection
- `tools/list`: List available tools
- `tools/call`: Execute a tool
- `resources/list`: List available resources
- `resources/read`: Read a resource
- `prompts/list`: List available prompts
- `prompts/get`: Get a specific prompt

#### `GET /sse`
Server-Sent Events endpoint for streaming responses.

**Event Types:**
- `connection`: Initial connection confirmation
- `heartbeat`: Keep-alive messages (every 30 seconds)
- Tool execution progress events (when implemented)

## Data Models

### MCPRequest
```python
class MCPRequest(BaseModel):
    jsonrpc: str = "2.0"
    method: str
    params: Optional[Dict[str, Any]] = None
    id: Optional[str] = None
```

### MCPResponse
```python
class MCPResponse(BaseModel):
    jsonrpc: str = "2.0"
    result: Optional[Any] = None
    error: Optional[Dict[str, Any]] = None
    id: Optional[str] = None
```

## Security Features

### CORS Configuration
- Default: Only localhost origins allowed
- Configurable via `--allow-external` flag
- Credentials supported for session management

### Security Middleware
1. **Origin Validation**
   - Checks Origin header for non-GET requests
   - Rejects requests from non-localhost origins
   - Returns 403 Forbidden for invalid origins

2. **Security Headers**
   - `X-Content-Type-Options: nosniff`
   - `X-Frame-Options: DENY`
   - `X-XSS-Protection: 1; mode=block`

### Rate Limiting
- Simple in-memory rate limiter
- 100 requests per minute per IP address
- Returns 429 Too Many Requests when exceeded
- 1-minute sliding window

## Session Management

### Session Structure
```python
{
    "id": "session-uuid",
    "created_at": timestamp,
    "last_activity": timestamp,
    "data_store": {}  # Isolated data storage
}
```

### Session Isolation
- Each session has its own data store
- Tool calls use session-specific storage
- Prevents data conflicts between users
- Automatic activity tracking

## Command-Line Interface

### Usage
```bash
python -m chatspatial.http_server [options]
```

### Options
- `--host`: Host to bind to (default: 127.0.0.1)
- `--port`: Port to bind to (default: 8000)
- `--reload`: Enable auto-reload for development
- `--allow-external`: Allow external connections (security warning)

### Examples
```bash
# Development mode with auto-reload
python -m chatspatial.http_server --reload

# Custom port
python -m chatspatial.http_server --port 8080

# Allow external connections (not recommended)
python -m chatspatial.http_server --host 0.0.0.0 --allow-external
```

## Implementation Details

### Request Flow
1. Client sends HTTP request to endpoint
2. Security middleware validates origin
3. Rate limiting middleware checks request count
4. Request routed to appropriate handler
5. For tool calls, session-specific data store is used
6. Response sent with security headers

### Tool Call Handling
```python
async def handle_tool_call(params: Dict[str, Any], session_id: Optional[str]) -> Any:
    # 1. Extract tool name and arguments
    # 2. Validate tool exists
    # 3. Switch to session-specific data store
    # 4. Execute tool
    # 5. Restore global data store
    # 6. Return result
```

### Error Handling
- MCP-compliant error responses
- Error codes follow JSON-RPC 2.0 specification
- Detailed error logging for debugging
- Client-friendly error messages

## Integration Examples

### JavaScript Client
```javascript
// Create session
const response = await fetch('http://localhost:8000/sessions', {
  method: 'POST'
});
const { session_id } = await response.json();

// Call tool
const rpcResponse = await fetch('http://localhost:8000/rpc', {
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
        data_path: '/path/to/data.h5ad'
      }
    },
    id: '1'
  })
});
```

### Python Client
```python
import requests

# Create session
session_resp = requests.post('http://localhost:8000/sessions')
session_id = session_resp.json()['session_id']

# Call tool
rpc_payload = {
    'jsonrpc': '2.0',
    'method': 'tools/call',
    'params': {
        'name': 'load_data',
        'arguments': {
            'data_path': '/path/to/data.h5ad'
        }
    },
    'id': '1'
}

headers = {'X-Session-Id': session_id}
result = requests.post('http://localhost:8000/rpc', 
                      json=rpc_payload, 
                      headers=headers)
```

## Best Practices

1. **Security**
   - Always run on localhost for local analysis
   - Use HTTPS proxy for production deployments
   - Implement authentication for multi-user scenarios

2. **Session Management**
   - Create new session for each analysis workflow
   - Clean up sessions when done
   - Monitor session count for resource management

3. **Error Handling**
   - Check response status codes
   - Handle rate limiting gracefully
   - Implement retry logic for transient failures

4. **Performance**
   - Use SSE for long-running operations
   - Batch multiple tool calls when possible
   - Monitor server health endpoint

## Future Enhancements

1. **WebSocket Support**: For bidirectional streaming
2. **Authentication**: OAuth2/JWT integration
3. **Metrics**: Prometheus endpoints for monitoring
4. **Clustering**: Redis-based session storage for scaling
5. **Progress Tracking**: Real-time updates via SSE