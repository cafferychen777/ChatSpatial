# What is Model Context Protocol (MCP)?

**Model Context Protocol (MCP)** is an open standard that connects AI applications to the systems where your data and tools live. Think of it as a universal adapter, like USB-C, but for AI. Just as USB-C provides a standard way to connect devices to peripherals, MCP offers a standardized, secure, and simple way to connect AI assistants like Claude to real-world applications, databases, and services.

This enables AI to move beyond its training data, giving it the context it needs to read files from your computer, search through internal knowledge bases, or update tasks in project management tools.

## The Problem MCP Solves

Before MCP, AI assistants were powerful but limited to information you manually provided. They couldn't:

- Access real-time information
- Interact with external tools and services (APIs, databases)
- Perform actions beyond text generation
- Connect to specialized software securely

Every integration required bespoke, custom development, making it difficult and time-consuming to scale AI capabilities. MCP solves this with a single, open protocol, creating a growing ecosystem of interoperable AI applications and tools.

## How MCP Works

MCP creates a standardized bridge between an **AI application (Host)** and external **Servers**. The Host (e.g., Claude Desktop, VS Code) runs **Clients**, with each client maintaining a dedicated, one-to-one connection to a server.

```mermaid
graph TB
    subgraph "AI Application (Host Process)"
        H[Host] --> C1[Client 1]
        H --> C2[Client 2]
    end

    subgraph "Local Machine"
        S1[Server 1<br/>Filesystem]
        R1[("Local Files")]
        C1 --> S1
        S1 <--> R1
    end

    subgraph "Internet"
        S2[Server 2<br/>External API (e.g., GitHub)]
        R2[("Remote API")]
        C2 --> S2
        S2 <--> R2
    end

    style H fill:#e1f5fe
    style C1,C2 fill:#f3e5f5
    style S1,S2 fill:#fff3e0
    style R1,R2 fill:#fce4ec
```

### Key Architectural Layers

1. **MCP Host**: The AI application (like Claude or VS Code) that manages one or more MCP clients. It handles the user interface, security policies, and AI model integration.
2. **MCP Client**: A component within the host that connects to a single MCP server, handling protocol negotiation and routing messages.
3. **MCP Server**: A program that provides context and capabilities by exposing tools, resources, and prompts. Servers can run locally (stdio transport) or remotely (Streamable HTTP transport).
4. **Transport Layer**: Defines the communication channel:
   - **`stdio`**: For local processes, exchanging JSON-RPC messages via stdin/stdout
   - **`Streamable HTTP`**: For remote services, using POST/GET requests with SSE support
5. **Data Layer**: Defines the JSON-RPC 2.0 based protocol for exchanging messages

## Core Concepts (The Building Blocks)

Servers provide functionality through three core primitives, each with different control:

| Primitive | Who Controls It | Purpose | Real-World Example |
| :--- | :--- | :--- | :--- |
| **Tools** | **Model**-controlled | Enables AI to perform actions | Search flights, send messages, create files |
| **Resources** | **Application**-controlled | Provides data for context | Documents, calendars, database schemas |
| **Prompts** | **User**-controlled | Reusable interaction templates | `/plan-vacation`, `/summarize-meetings` |

### 1. Tools (AI Actions)
Functions that AI assistants can call to perform actions. Each tool has a schema-defined interface. The model requests tool execution, but users must provide explicit approval.

```json
{
  "name": "searchFlights",
  "description": "Search for available flights",
  "inputSchema": {
    "type": "object",
    "properties": {
      "origin": { "type": "string", "description": "Departure city" },
      "destination": { "type": "string", "description": "Arrival city" },
      "date": { "type": "string", "format": "date" }
    },
    "required": ["origin", "destination", "date"]
  }
}
```

### 2. Resources (Context Data)
Data sources that AI assistants can read to gain context. Resources are identified by URI and can be anything from local files to API endpoints. The application decides how to retrieve and use this data.

```json
{
  "uri": "file:///Documents/Travel/passport.pdf",
  "name": "passport.pdf",
  "mimeType": "application/pdf"
}
```

### 3. Prompts (Interaction Templates)
Reusable, user-controlled templates that structure common tasks. They accept arguments to create consistent workflows. Users typically invoke them via slash commands or command palette.

```json
{
  "name": "plan-vacation",
  "title": "Plan a vacation",
  "description": "Guide through the vacation planning process",
  "arguments": [
    { "name": "destination", "type": "string", "required": true },
    { "name": "duration", "type": "number", "description": "days" }
  ]
}
```

## Protocol Implementation Details

### Communication Format
- All messages use **JSON-RPC 2.0** format
- UTF-8 encoded, newline-delimited in stdio mode
- Support for requests, responses, and notifications

### Lifecycle & Handshake
Connection establishment follows a strict sequence:
1. Client sends `initialize` request with capabilities
2. Server responds with its capabilities
3. Client sends `initialized` notification
4. Connection is ready for operation

### Capability Negotiation
During initialization, both sides declare supported features:
- Tools, resources, prompts support
- Sampling capabilities (server requesting LLM inference)
- Root directory awareness
- Authorization mechanisms

## Benefits of MCP

### For Users
- **Natural Language Interface**: Interact with complex tools conversationally
- **Access to Your Context**: AI can securely access your documents, data, and tools
- **Real-time Results**: Live data and immediate analysis
- **Secure by Design**: Full control with explicit permission for every action

### For Developers
- **Standardized API**: Build once, work with all MCP-compatible AI clients
- **Reduced Complexity**: Focus on features instead of custom connectors
- **Growing Ecosystem**: Leverage open-source servers from Anthropic and community
- **Future-proof**: Open standard ensures compatibility as AI evolves

## MCP vs. Traditional Approaches

| Aspect | Traditional Tools | MCP-Enabled Tools |
| :--- | :--- | :--- |
| **Interface** | Command line, GUI | Natural language |
| **Learning Curve** | Steep, tool-specific | Minimal, conversational |
| **Integration** | Manual scripting | Automatic discovery |
| **Error Handling** | Manual debugging | AI-assisted troubleshooting |
| **Workflow** | Linear, rigid | Flexible, composable |

## Real-World Example: Spatial Transcriptomics

### Traditional Workflow (CLI)
```bash
# 1. Find relevant files
grep -r "Q3-roadmap" ~/Documents/project-alpha

# 2. Open and read each file
cat ~/Documents/project-alpha/meeting-notes.md

# 3. Load spatial data
import scanpy as sc
adata = sc.read_h5ad("data.h5ad")

# 4. Preprocess
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# 5. Analyze
import squidpy as sq
sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(adata)
```

### MCP-Powered Workflow with ChatSpatial
```
👤 User: "Load my Visium data and identify spatial domains"

🤖 AI: I'll analyze your spatial transcriptomics data:
1. Loading your Visium dataset
2. Performing preprocessing 
3. Identifying spatial domains using SpaGCN
4. Creating visualizations

[AI automatically executes using ChatSpatial MCP server]
```

## MCP Ecosystem

<div align="center">
  <strong>70+</strong> Compatible Clients &nbsp;&nbsp;&nbsp;&nbsp;
  <strong>1000+</strong> Available Servers &nbsp;&nbsp;&nbsp;&nbsp;
  <strong>9</strong> Official SDKs
</div>

### Compatible AI Clients
- **Claude Desktop & Claude.ai**: Native MCP support
- **VS Code**: Integration via GitHub Copilot
- **LM Studio**: Connect local models to tools
- **Cursor, Warp, Zed**: AI-native editors with MCP support

### Popular MCP Servers
- **File Systems**: Secure access to local and cloud files
- **Databases**: PostgreSQL, SQLite, MongoDB
- **APIs**: GitHub, Slack, Google Drive, Sentry
- **Development Tools**: Git, Docker, Kubernetes
- **Scientific Tools**: **ChatSpatial** for spatial transcriptomics

## Getting Started

### 1. Install an MCP Client
Download [Claude Desktop](https://claude.ai/download) or [VS Code](https://code.visualstudio.com/).

### 2. Configure MCP Servers
Configure your client's config file:

**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
**Windows**: `%APPDATA%\Claude\claude_desktop_config.json`

Example configuration for ChatSpatial:
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "python",
      "args": ["-m", "chatspatial"],
      "env": {}
    },
    "filesystem": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-filesystem",
        "/Users/your_username/Desktop"
      ]
    }
  }
}
```

### 3. Start Using Tools
Restart your client. You'll see an indicator that MCP tools are available. Start asking your AI assistant to perform tasks naturally.

## Security and Authorization

Security is a core principle of MCP:

- **Explicit User Consent**: No tool execution without approval
- **Granular Permissions**: Configure which tools AI can use
- **Sandboxing**: Isolated server connections
- **Secure Authorization**: OAuth 2.1 for remote servers
- **Local-First Option**: Keep sensitive data on your machine with stdio transport

## Key Code Patterns (Python SDK Example)

### Server Implementation - Defining Tools
```python
from mcp.server.fastmcp import FastMCP

# Initialize server
mcp = FastMCP("weather_server")

@mcp.tool()
async def get_forecast(latitude: float, longitude: float) -> str:
    """
    Get weather forecast for a location.
    
    Args:
        latitude: Latitude of the location.
        longitude: Longitude of the location.
    """
    # Implementation
    forecast_data = await call_weather_api(latitude, longitude)
    return f"Forecast: {forecast_data}"

if __name__ == "__main__":
    mcp.run(transport='stdio')
```

The `@mcp.tool()` decorator automatically converts functions into MCP tools, parsing function names, types, and docstrings into tool schemas.

## The Future of MCP (Roadmap)

MCP is rapidly evolving with focus on:

- **Agents**: Support for long-running, asynchronous operations
- **Authentication & Security**: Enterprise-grade authorization and SSO
- **Validation & Tooling**: Reference implementations and compliance testing
- **Registry**: Centralized server discovery and distribution
- **Multimodality**: Support for video and streaming interactive experiences

## Learn More

### 📚 Official Resources
- **🌐 Official Website**: [modelcontextprotocol.io](https://modelcontextprotocol.io) - Complete documentation
- **📜 Specification**: [spec.modelcontextprotocol.io](https://spec.modelcontextprotocol.io) - Technical specification
- **💻 GitHub**: [github.com/modelcontextprotocol](https://github.com/modelcontextprotocol) - Open source repositories
- **💬 Discord**: [Join Discord](https://discord.gg/6CSzBmMkjX) - Community discussion

### 🎥 Video Resources
- **📺 Simple MCP Demo**: [YouTube Tutorial](https://www.youtube.com/watch?v=sfCBCyNyw7U) - Clear explanation with examples
- **🎬 Claude Code Introduction**: [Official Video](https://www.youtube.com/watch?v=AJpK3YTTKZ4) - Anthropic announcement
- **📖 Claude Code Tutorial**: [Setup Guide](https://www.youtube.com/watch?v=SUysp3sJHbA) - Installation and usage

### 🛠️ Developer Resources
- **📦 MCP Servers Collection**: [modelcontextprotocol/servers](https://github.com/modelcontextprotocol/servers) - Reference implementations
- **🔨 MCP SDKs**: TypeScript, Python, Go, Kotlin, Swift, Java, C#, Ruby, Rust
- **🧬 ChatSpatial Example**: See MCP enabling spatial transcriptomics analysis

---

Ready to unlock the full potential of your AI assistant? Try ChatSpatial with Claude Desktop and experience how natural language transforms spatial transcriptomics workflows!