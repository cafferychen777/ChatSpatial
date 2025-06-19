---
name: MCP Integration Issue
about: Report issues related to Model Context Protocol integration
title: '[MCP] '
labels: 'mcp'
assignees: ''

---

**MCP Client Information**
- Client: [e.g. Claude Desktop, other MCP client]
- Client version: [e.g. 1.0.0]
- Transport method: [e.g. stdio, SSE]

**Issue Description**
A clear and concise description of the MCP integration issue.

**Server Configuration**
Please share your MCP server configuration (remove any sensitive information):
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "...",
      "args": [...],
      "env": {...}
    }
  }
}
```

**Error Logs**
If available, please include relevant error logs from:
- MCP client logs
- ChatSpatial server logs
- Terminal output

**Steps to Reproduce**
1. Start ChatSpatial MCP server with '...'
2. Connect from MCP client
3. Execute command '...'
4. See error

**Expected Behavior**
What should have happened?

**Additional Context**
Any other information that might help diagnose the issue.
