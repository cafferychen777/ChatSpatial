# Deployment Guide

This section helps you choose and implement the right deployment strategy for ChatSpatial based on your specific needs and resources.

## Why Two Deployment Options?

ChatSpatial supports two distinct deployment scenarios because different users have different needs:

### 1. Local Deployment (Individual Use)

**What it is:** Install ChatSpatial directly on your personal computer and use it with Claude Desktop, Continue, or other MCP-compatible clients.

**Best for:**

- Individual researchers with sufficient local computing resources
- Quick prototyping and small dataset analysis
- Users who already have Claude Desktop or similar MCP clients
- Privacy-sensitive data that cannot leave local machines

**How it works:** ChatSpatial runs as an MCP server on your machine, directly communicating with your AI assistant through the MCP protocol.

**Get Started:** [Installation Guide](../getting-started/installation.md)

### 2. Server Deployment (Team Collaboration)

**What it is:** Deploy ChatSpatial on a centralized server with a web interface (Open WebUI), allowing multiple team members to access it through their browsers.

**Best for:**

- Research teams sharing computational resources
- Large datasets requiring 128GB+ RAM
- Groups without individual high-performance computers
- Centralized data management and collaboration

**How it works:** ChatSpatial runs on a server, accessed through Open WebUI. Since Open WebUI does not natively support MCP, we use an mcpo proxy to bridge the protocols.

**Get Started:** [Server Deployment Guide](./server-deployment.md)

## Quick Decision Guide

Ask yourself these questions to choose the right deployment:

1. **Do you have Claude Desktop or another MCP client installed?**
   - Yes → Consider local deployment
   - No → Consider server deployment

2. **Will multiple people need access?**
   - Yes → Server deployment is better
   - No → Local deployment is sufficient

3. **Do you have a powerful personal computer (32GB+ RAM)?**
   - Yes → Local deployment will work
   - No → Server deployment recommended

4. **Is your data sensitive and cannot leave your machine?**
   - Yes → Local deployment required
   - No → Either option works

## Comparison Table

| Aspect | Local Deployment | Server Deployment |
|--------|------------------|-------------------|
| **Users** | Single user | Multiple users |
| **Access** | MCP client (Claude Desktop) | Web browser |
| **Setup Complexity** | Simple (pip install) | Moderate (Docker, nginx) |
| **Resource Sharing** | No | Yes |
| **Minimum RAM** | 16GB | 64GB |
| **Data Location** | Local machine | Centralized server |
| **Cost** | LLM API only | LLM API + Server |
| **Best For** | Individual analysis | Team collaboration |

## Technical Architecture Differences

### Local Deployment
```
Claude Desktop/MCP Client
        ↓
    ChatSpatial MCP Server
        ↓
    Local Computing Resources
```

### Server Deployment
```
Team Members (Browsers)
        ↓
    Open WebUI (Web Interface)
        ↓
    mcpo Proxy (Protocol Bridge)
        ↓
    ChatSpatial MCP Server
        ↓
    Server Computing Resources
```

## Next Steps

- **For individual use with Claude Desktop** → [Installation Guide](../getting-started/installation.md)
- **For team collaboration via web browser** → [Server Deployment Guide](./server-deployment.md)
- **Need help deciding?** → [Contact Support](https://github.com/cafferychen777/ChatSpatial/issues)