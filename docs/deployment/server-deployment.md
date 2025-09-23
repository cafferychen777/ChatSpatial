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

## Simple HTTP Server

Run ChatSpatial as an HTTP server for browser-based access:

```bash
# Install ChatSpatial
pip install -e ".[full]"

# Run HTTP server
python -m chatspatial.http_server --port 8000
```

Access at `http://your-server:8000`

## MCP Server Mode

For integration with Claude Desktop or other MCP clients:

```bash
# Standard stdio mode (for Claude Desktop)
python -m chatspatial

# SSE mode (for HTTP clients)
python -m chatspatial server --transport sse --port 8000
```

## Requirements

- Python 3.8-3.12
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