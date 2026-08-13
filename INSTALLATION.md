# Installation

This repository-level page is intentionally short. The canonical installation
guide lives in [docs/installation.md](docs/installation.md), so zero-environment
`uvx` setup, persistent environments, platform notes, and optional dependency
guidance stay in one place.

For most users, install `uv` once and register the server directly:

```bash
codex mcp add chatspatial -- uvx --from chatspatial chatspatial server
```

Use the docs in this order:

1. [Installation](docs/installation.md) - choose uvx, a persistent Python environment, or Docker
2. [Configuration](docs/advanced/configuration.md) - register ChatSpatial in an MCP-compatible client
3. [Quick Start](docs/quickstart.md) - run the first analysis prompt
4. [Troubleshooting](docs/advanced/troubleshooting.md) - fix install or runtime issues

For container-specific details, use [Docker / GHCR](docs/docker.md).
