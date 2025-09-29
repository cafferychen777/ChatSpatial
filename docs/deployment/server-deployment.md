---
layout: default
title: Server Deployment
parent: Deployment
nav_order: 1
description: "Deploy ChatSpatial for team access"
---

# Server Deployment Guide
{: .fs-7 }

Deploy ChatSpatial on a server for team-wide access to spatial transcriptomics analysis.

## MCP Server Deployment

ChatSpatial implements the Model Context Protocol (MCP) with two transport modes:

### Local Mode (stdio)
For Claude Desktop and local clients:

```bash
# Install ChatSpatial
pip install -e ".[full]"

# Run in stdio mode (default)
python -m chatspatial
```

### Remote Mode (SSE)
For remote access and HTTP-based clients:

```bash
# SSE mode with HTTP/WebSocket support
python -m chatspatial server --transport sse --port 8000 --host 0.0.0.0
```

**Note**: Use with mcp-remote proxy for Claude Desktop remote connections:
```bash
# Client-side proxy setup
npx mcp-remote http://your-server:8000/mcp/
```

## Requirements

- Python 3.10-3.12 (MCP requires Python 3.10+)
- 64 GB RAM recommended for large datasets
- Fast SSD storage
- Network access for team

## Configuration

Create `config.yaml` for server settings:

```yaml
server:
  host: 0.0.0.0
  port: 8000
  max_connections: 10
  
data:
  cache_dir: /data/chatspatial_cache
  max_dataset_size: 10GB
```

## Security Notes

- Use VPN or internal network only
- Add authentication if exposing to internet
- Regular backups of analysis results

## Getting Help

- GitHub Issues for deployment questions
- Check logs in `/var/log/chatspatial/`

---

**Next:** [Configuration Guide](configuration.html)