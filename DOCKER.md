# Docker user guide

This guide provides the fastest reproducible way to run ChatSpatial without solving the scientific Python dependency stack on your host machine.

## Pull the published image

```bash
docker pull ghcr.io/cafferychen777/chatspatial:latest
```

## Verify the image

```bash
docker run --rm ghcr.io/cafferychen777/chatspatial:latest --version
docker run --rm ghcr.io/cafferychen777/chatspatial:latest server --help
```

## Run the MCP server with stdio transport

Use this command in an MCP client configuration:

```bash
docker run --rm -i \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:latest server --transport stdio
```

Inside the container, mounted data are available under `/data` and generated outputs are written under `/outputs` by default.

## Run the SSE server

```bash
docker run --rm -p 8000:8000 \
  -v /absolute/path/to/your/data:/data:ro \
  -v /absolute/path/to/outputs:/outputs \
  ghcr.io/cafferychen777/chatspatial:latest server --transport sse --host 0.0.0.0 --port 8000
```

Then connect your MCP client to `http://localhost:8000/sse`.

## Build locally from source

From the repository root:

```bash
docker build -t chatspatial:latest .
```

Then replace `ghcr.io/cafferychen777/chatspatial:latest` with `chatspatial:latest` in the commands above.

## Notes

- The Docker image installs ChatSpatial in a controlled Linux environment so users do not need local compilers for packages such as `gseapy` or `llvmlite`.
- The default container working directory is `/workspace`.
- `CHATSPATIAL_OUTPUT_DIR` defaults to `/outputs`.
