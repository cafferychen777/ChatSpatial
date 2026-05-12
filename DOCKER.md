# Docker / GHCR

This repository-level page is intentionally short. The canonical Docker guide
lives in [docs/docker.md](docs/docker.md), so image tags, volume mounts, MCP
client examples, SSE mode, local builds, and troubleshooting stay in one place.

Key rules:

- Use a versioned image such as `ghcr.io/cafferychen777/chatspatial:v1.2.9` for reproducible analyses.
- Mount host data into the container, usually as `/data`, and prompt ChatSpatial with the container path.
- Mount a writable output directory as `/outputs`; the Docker image sets `CHATSPATIAL_OUTPUT_DIR=/outputs`.
- Use `--rm -i`, not `-it`, for MCP stdio servers.

Full guide: [docs/docker.md](docs/docker.md).
