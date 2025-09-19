---
layout: default
title: Server Deployment
parent: Deployment
nav_order: 1
description: "Deploy ChatSpatial as a production server for team access"
---

# Server Deployment Guide
{: .fs-7 }

!!! success "Updated Deployment Options (2025-08-28)"
    **Based on extensive testing, we recommend the following deployment approaches:**
    
    ### ✅ **Recommended: Native Deployment with AnythingLLM**
    - Simple 5-minute setup without Docker complexity
    - Native MCP support via stdio protocol
    - Proven to work with ChatSpatial's 16 tools
    - See [AnythingLLM Deployment](#anythingllm-deployment) section below
    
    ### ⚠️ **Not Recommended: Open WebUI** 
    - Tool integration is currently broken (Open WebUI bug)
    - Even simple tools fail to invoke properly
    - Retained below for historical reference only
    
    ### 🚀 **Alternative: Remote Development**
    - Use Cursor/Windsurf with SSH connection to server
    - Zero configuration for end users
    - Built-in LLM support without API keys

Deploy ChatSpatial on a centralized server to enable team-wide spatial transcriptomics analysis through various interfaces.

## Overview

This guide describes how to deploy ChatSpatial with Open WebUI, providing a browser-based interface for spatial transcriptomics analysis. The deployment uses high-performance language models (LLMs) through commercial APIs to ensure accurate analysis results.

## Architecture

The server deployment consists of three main components:

1. **Open WebUI**: Web interface for user interaction
2. **mcpo Proxy**: Bridge between Open WebUI and MCP servers
3. **ChatSpatial MCP**: Spatial transcriptomics analysis tools

### Deployment Architectures

#### Option A: Full Docker Deployment (Production)
```
Users (Web Browser)
        ↓
    Open WebUI (Docker)
        ├── Commercial LLM APIs (Claude/GPT-5)
        └── mcpo + ChatSpatial (Docker) → Server Resources
```

#### Option B: Hybrid Deployment (Development)
```
Users (Web Browser)
        ↓
    Open WebUI (Docker)
        ├── Commercial LLM APIs (Claude/GPT-5)
        └── mcpo + ChatSpatial (Host) → Server Resources
```

!!! note "MCP Integration"
    Open WebUI does not natively support the Model Context Protocol (MCP). The mcpo proxy server is required to bridge MCP's stdio-based communication with Open WebUI's REST API requirements.

!!! important "Deployment Recommendations"
    - **Production**: Use Docker for everything (Option A) - provides best isolation, security, and portability
    - **Development**: Run mcpo and ChatSpatial on host (Option B) - easier debugging and faster iteration
    - **Both ChatSpatial and mcpo support Python 3.11**, simplifying deployment requirements

## Requirements

### Hardware Requirements

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| CPU | 16 cores | 32 cores | For parallel processing |
| RAM | 64 GB | 128-256 GB | Large datasets require more memory |
| Storage | 500 GB SSD | 1 TB NVMe SSD | Fast I/O is critical |
| GPU | Optional | NVIDIA with CUDA | For deep learning models |
| Network | 1 Gbps | 10 Gbps | For data transfer |

### Software Requirements

- Ubuntu 22.04 LTS or Rocky Linux 9
- Docker and Docker Compose v2 (included with Docker Desktop)
- Python 3.11+ (required for mcpo proxy)
- Conda or Miniconda

### API Requirements

You will need API keys from at least one commercial LLM provider:

- **Anthropic Claude API**: Excellent reasoning and code understanding
- **OpenAI GPT-5 API**: Strong analysis capabilities
- **Alternative**: DeepSeek or other compatible APIs

## Installation Steps

### Step 1: System Preparation

Update the system and install required packages:

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install dependencies
sudo apt install -y curl wget git build-essential \
    python3-pip python3-venv nginx certbot \
    python3-certbot-nginx htop tmux

# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER

# Docker Compose v2 is included with Docker
# Verify it is installed
docker compose version

# Log out and back in for group changes
```

### Step 2: Install Python 3.11+

The mcpo proxy requires Python 3.11 or newer. Install via Conda:

```bash
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

# Initialize
~/miniconda3/bin/conda init bash
source ~/.bashrc

# Create Python 3.11 environment
conda create -n spatial-platform python=3.11 -y
conda activate spatial-platform
```

### Step 3: Create Project Directory

```bash
# Option A: System-wide installation (requires sudo)
sudo mkdir -p /opt/spatial-platform
cd /opt/spatial-platform
sudo mkdir -p data/{uploads,processed,results} configs logs scripts

# Set permissions (adjust based on your OS)
# Linux: use $USER:$USER
# macOS: use $USER:staff  
sudo chown -R $USER:$USER /opt/spatial-platform  # Linux
# sudo chown -R $USER:staff /opt/spatial-platform  # macOS

# Option B: User directory (recommended for testing)
mkdir -p ~/spatial-platform
cd ~/spatial-platform
mkdir -p data/{uploads,processed,results} configs logs scripts
```

### Step 4: Configure LLM APIs

Create environment configuration:

```bash
cat > /opt/spatial-platform/.env << 'EOF'
# Claude API (if using Anthropic)
ANTHROPIC_API_KEY=your-anthropic-api-key

# OpenAI API (if using GPT-5)
OPENAI_API_KEY=your-openai-api-key
OPENAI_API_BASE=https://api.openai.com/v1

# Security
WEBUI_SECRET_KEY=generate-a-secure-key-here
MCP_API_KEY=generate-another-secure-key

# Configuration
DEFAULT_MODEL=claude-sonnet-4-20250514
WEBUI_AUTH=True
ENABLE_SIGNUP=False
EOF
```

### Step 5: Install ChatSpatial

```bash
# Stay in the Python 3.11 environment
conda activate spatial-platform

# Clone and install ChatSpatial
cd /opt/spatial-platform
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial
pip install -e .

# Verify installation
chatspatial --version
```

### Step 6: Deploy ChatSpatial and mcpo Proxy

The mcpo proxy converts MCP protocol to OpenAPI format for Open WebUI. Choose your deployment strategy based on your needs:

#### Option A: Everything in Docker (Recommended for Production)

This approach provides the best isolation and is easiest to manage:

```bash
# Create a Dockerfile for ChatSpatial with mcpo
cat > /opt/spatial-platform/Dockerfile.chatspatial << 'EOF'
FROM python:3.11-slim

# Install dependencies
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

# Install ChatSpatial and mcpo
WORKDIR /app
COPY ./ChatSpatial /app/ChatSpatial
RUN pip install -e /app/ChatSpatial
RUN pip install mcpo

# Expose port
EXPOSE 8000

# Start mcpo with ChatSpatial
CMD ["mcpo", "--host", "0.0.0.0", "--port", "8000", "--api-key", "${MCP_API_KEY}", "--", "python", "-m", "chatspatial"]
EOF

# Build the image
docker build -f Dockerfile.chatspatial -t chatspatial-mcpo:latest .
```

#### Option B: Host Deployment (Simpler for Development)

Run both ChatSpatial and mcpo directly on the host:

```bash
# Install mcpo in the Python 3.11 environment
conda activate spatial-platform
pip install mcpo

# Create a startup script
cat > /opt/spatial-platform/scripts/start-mcpo.sh << 'EOF'
#!/bin/bash
source ~/miniconda3/bin/activate spatial-platform
cd /opt/spatial-platform

# Start mcpo with ChatSpatial
mcpo --host 0.0.0.0 --port 8000 \
     --api-key ${MCP_API_KEY} \
     -- python -m chatspatial
EOF

chmod +x /opt/spatial-platform/scripts/start-mcpo.sh

# Run mcpo in the background
nohup /opt/spatial-platform/scripts/start-mcpo.sh > /opt/spatial-platform/logs/mcpo.log 2>&1 &
```

#### Option C: Docker-outside-of-Docker (Advanced)

For managing multiple MCP services with better isolation:

```bash
# Create configuration for DooD approach
cat > /opt/spatial-platform/configs/dood-config.json << 'EOF'
{
  "mcpServers": {
    "chatspatial": {
      "command": "docker",
      "args": [
        "run", "--rm", "-i",
        "--network", "host",
        "chatspatial:latest",
        "python", "-m", "chatspatial"
      ]
    }
  }
}
EOF
```

!!! tip "Which option to choose?"
    - **Production servers**: Use Option A (Everything in Docker) for best isolation and easy management
    - **Development/Testing**: Use Option B (Host deployment) for easier debugging
    - **Complex environments**: Use Option C (DooD) when managing multiple MCP services

!!! note "Testing Status"
    Option B (Host deployment) has been partially tested. Option A (Docker deployment) is based on best practices from Open WebUI documentation. We recommend testing in your specific environment before production deployment.

### Step 7: Deploy Complete Stack with Docker Compose

!!! note "Image Download Size"
    The Open WebUI Docker image is approximately 1.2GB in size. The initial download may take 10-30 minutes depending on your internet connection speed. This is normal and only needs to be done once.

Choose the docker-compose configuration based on your deployment option from Step 6:

#### For Option A: Everything in Docker

```bash
cat > /opt/spatial-platform/docker-compose.yml << 'EOF'
services:
  chatspatial-mcpo:
    image: chatspatial-mcpo:latest
    container_name: chatspatial-mcpo
    restart: always
    ports:
      - "8000:8000"
    environment:
      - MCP_API_KEY=${MCP_API_KEY}
    volumes:
      - ./data/spatial:/data

  open-webui:
    image: ghcr.io/open-webui/open-webui:main
    container_name: open-webui
    restart: always
    ports:
      - "3000:8080"
    volumes:
      - ./data/webui:/app/backend/data
    environment:
      - OPENAI_API_KEY=${OPENAI_API_KEY}
      - OPENAI_API_BASE=${OPENAI_API_BASE}
      - WEBUI_SECRET_KEY=${WEBUI_SECRET_KEY}
      - WEBUI_AUTH=${WEBUI_AUTH}
      - ENABLE_SIGNUP=${ENABLE_SIGNUP}
      - WEBUI_NAME=Spatial Analysis Platform
    depends_on:
      - chatspatial-mcpo

volumes:
  data:
EOF
```

#### For Option B: Host Deployment

```bash
cat > /opt/spatial-platform/docker-compose.yml << 'EOF'
services:
  open-webui:
    image: ghcr.io/open-webui/open-webui:main
    container_name: open-webui
    restart: always
    ports:
      - "3000:8080"
    volumes:
      - ./data:/app/backend/data
    environment:
      - OPENAI_API_KEY=${OPENAI_API_KEY}
      - OPENAI_API_BASE=${OPENAI_API_BASE}
      - WEBUI_SECRET_KEY=${WEBUI_SECRET_KEY}
      - WEBUI_AUTH=${WEBUI_AUTH}
      - ENABLE_SIGNUP=${ENABLE_SIGNUP}
      - WEBUI_NAME=Spatial Analysis Platform
    extra_hosts:
      - "host.docker.internal:host-gateway"

volumes:
  data:
EOF
```

Start the services:

```bash
# Start services
docker compose up -d

# Wait for initialization
sleep 60

# Check if services are running
docker compose ps

# For Option B, also verify mcpo is running on host
ps aux | grep mcpo
```

### Step 8: Configure Web Server

Setup nginx as reverse proxy with SSL:

```bash
# Create nginx configuration
sudo tee /etc/nginx/sites-available/spatial-platform << 'EOF'
server {
    listen 80;
    server_name your-domain.com;
    return 301 https://$server_name$request_uri;
}

server {
    listen 443 ssl http2;
    server_name your-domain.com;

    # SSL certificates (will be configured by certbot)
    ssl_certificate /etc/letsencrypt/live/your-domain.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/your-domain.com/privkey.pem;
    
    # Security headers
    add_header X-Frame-Options "SAMEORIGIN" always;
    add_header X-Content-Type-Options "nosniff" always;
    add_header X-XSS-Protection "1; mode=block" always;

    # Open WebUI
    location / {
        proxy_pass http://localhost:3000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 86400;
    }

    # MCP API
    location /api/mcp/ {
        proxy_pass http://localhost:8000/;
        proxy_set_header X-API-Key $http_x_api_key;
        proxy_set_header Host $host;
    }
}
EOF

# Enable site
sudo ln -s /etc/nginx/sites-available/spatial-platform /etc/nginx/sites-enabled/

# Get SSL certificate
sudo certbot --nginx -d your-domain.com

# Reload nginx
sudo systemctl reload nginx
```

## Configuration

### Open WebUI Setup

1. **Verify Service is Running**:

   ```bash
   # Check container health
   docker compose ps
   
   # Test API endpoint
   curl http://localhost:3000/api/config
   ```

2. **Access Web Interface**:
   - Navigate to `https://your-domain.com` (or `http://localhost:3000` for testing)
   - First user to register becomes administrator
   - If auth is disabled for testing, you can access directly

3. **Configure MCP Connection** (via Settings → Connections):
   - Add OpenAPI Server:
     - Name: ChatSpatial
     - Base URL:
       - For Option A (Docker): `http://chatspatial-mcpo:8000`
       - For Option B (Host): `http://host.docker.internal:8000`
     - API Key: Your MCP_API_KEY from `.env`
   - Test the connection by clicking "Test Connection"
   - If successful, ChatSpatial tools will be available to the LLM

### Adding Claude Support

Since Open WebUI does not natively support Claude API, add it through Functions:

1. Go to Settings → Functions
2. Import the Anthropic function from: [Open WebUI Pipelines](https://github.com/open-webui/pipelines)
3. Configure with your Anthropic API key

### User Management

Configure authentication methods:

```bash
# Update .env for LDAP/OAuth
cat >> /opt/spatial-platform/.env << 'EOF'
# Optional: OAuth configuration
OAUTH_PROVIDER=google
OAUTH_CLIENT_ID=your-client-id
OAUTH_CLIENT_SECRET=your-client-secret
EOF

# Restart Open WebUI
cd /opt/spatial-platform
docker compose restart open-webui
```

## Usage

### For End Users

1. **Access the Platform**: Navigate to `https://your-domain.com`
2. **Login**: Use your credentials
3. **Start Analysis**: Click "New Chat"
4. **Upload Data**: Use the file upload feature
5. **Interact**: Describe your analysis needs in natural language

Example queries:

- "Load the mouse brain Visium dataset and show spatial domains"
- "Identify spatially variable genes in my data"
- "Perform trajectory analysis on the uploaded dataset"

### For Administrators

Monitor system health:

```bash
# Check service status
docker compose ps

# View logs
docker compose logs -f open-webui
docker compose logs -f mcpo

# Resource usage
docker stats
```

## Maintenance

### Backup Strategy

Create automated backups:

```bash
#!/bin/bash
# Save as /opt/spatial-platform/scripts/backup.sh

BACKUP_DIR="/opt/backups/spatial-platform"
DATE=$(date +%Y%m%d_%H%M%S)

# Create backup
mkdir -p "$BACKUP_DIR"
tar czf "$BACKUP_DIR/backup_$DATE.tar.gz" \
    /opt/spatial-platform/data \
    /opt/spatial-platform/configs

# Keep only last 7 days
find "$BACKUP_DIR" -type f -mtime +7 -delete
```

Schedule with cron:

```bash
# Add to crontab
0 2 * * * /opt/spatial-platform/scripts/backup.sh
```

### Updates

Update components regularly:

```bash
# Update Docker images
docker compose pull
docker compose up -d

# Update ChatSpatial
cd /opt/spatial-platform/ChatSpatial
git pull
pip install -e . --upgrade
```

## Testing the Integration

### Verify mcpo is Working

Test that mcpo correctly exposes ChatSpatial tools:

```bash
# Check mcpo is running
ps aux | grep mcpo

# Test the API endpoint (replace with your actual API key)
curl -H "X-API-Key: test-mcp-key-12345" http://localhost:8000/tools

# You should see a JSON response listing ChatSpatial tools like:
# load_data, preprocessing, spatial_analysis, etc.
```

### Test in Open WebUI

1. Access Open WebUI at `http://localhost:3000`
2. Go to Settings → Connections → OpenAPI
3. Add the mcpo server:
   - URL: `http://localhost:8000`
   - API Key: Your MCP_API_KEY
4. Click "Test Connection"
5. Start a new chat and try: "List available spatial transcriptomics analysis tools"

## Troubleshooting

### Common Issues

#### mcpo Cannot Start

Check Python version and environment:

```bash
conda activate spatial-platform
python --version  # Should be 3.11+
which mcpo  # Should be in the conda environment

# If mcpo is not installed:
pip install mcpo
```

#### Open WebUI Cannot Connect to mcpo

Verify mcpo is running and accessible:

```bash
# Check if mcpo process is running
ps aux | grep mcpo

# Test direct connection
curl -H "X-API-Key: your-key" http://localhost:8000/tools

# Check mcpo logs
tail -f /opt/spatial-platform/logs/mcpo.log
```

#### High Memory Usage

Monitor and limit container resources:

```bash
# Add to docker-compose.yml under each service
deploy:
  resources:
    limits:
      memory: 8G
```

#### API Rate Limits

Implement usage monitoring:

```bash
# Check API usage in logs
docker compose logs open-webui | grep -i "api"
```

### Getting Help

- ChatSpatial Issues: [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
- Open WebUI Documentation: [docs.openwebui.com](https://docs.openwebui.com/)
- MCP Protocol: [modelcontextprotocol.io](https://modelcontextprotocol.io/)

## Security Considerations

### Essential Security Measures

1. **API Key Management**
   - Store keys in environment variables
   - Rotate keys regularly
   - Monitor usage for anomalies

2. **Network Security**
   - Use HTTPS with valid certificates
   - Configure firewall rules
   - Limit access to necessary ports

3. **Data Protection**
   - Encrypt sensitive data at rest
   - Implement access controls
   - Regular security audits

### Firewall Configuration

```bash
# Configure UFW firewall
sudo ufw allow 22/tcp   # SSH
sudo ufw allow 80/tcp   # HTTP
sudo ufw allow 443/tcp  # HTTPS
sudo ufw enable
```

## Alternative Solution 1: AnythingLLM with Native MCP Support

!!! success "Tested and Recommended (2025-08-28)"
    AnythingLLM has been successfully tested with ChatSpatial. The native deployment (non-Docker) provides the simplest and most reliable solution for laboratory servers.

### Overview

AnythingLLM is an all-in-one AI application that includes full MCP compatibility. Unlike Open WebUI, it has a built-in MCP Hypervisor that can manage MCP servers directly, potentially avoiding the tool integration issues we encountered with Open WebUI.

### Architecture

```
Users (Web Browser)
        ↓
    AnythingLLM (Native Installation)
        ├── LLM Provider (Claude/GPT/Ollama)
        └── MCP Hypervisor (stdio)
             └── ChatSpatial MCP Server
```

### Key Features

- **Native MCP Support**: Built-in MCP Hypervisor with stdio protocol
- **Multi-user Support**: Web-based interface for team access
- **Multiple LLM Providers**: Claude, GPT, Ollama, and more
- **Agent System**: AI agents with direct ChatSpatial tool access
- **Simple Setup**: 5-minute installation without containers

### Deployment Approach

!!! success "Native Installation Required"
    AnythingLLM must be installed directly on the host system (not in Docker) for MCP to work properly.
    The stdio protocol used by MCP requires direct process communication, which doesn't work across container boundaries.

3. **Enable in Agent Settings**:
   - Go to Admin → Agents
   - Click on MCP Servers
   - Enable ChatSpatial server
   - Tools will be available as `@@mcp_chatspatial`

### Usage After Native Deployment

1. **Access Web Interface**: http://your-server:3001
2. **Configure LLM Provider**: Settings → LLM (Claude/GPT/Ollama)
3. **Enable MCP**: Settings → Agent Skills → MCP Hypervisor
4. **Test ChatSpatial**: "List all spatial transcriptomics tools"
5. **Start Analysis**: Upload data and use natural language commands

### Why Native Deployment Works Best

- **Direct stdio Communication**: No container boundaries to cross
- **Simple Path Configuration**: Absolute paths just work
- **Zero Network Overhead**: Local process communication
- **Easy Debugging**: Standard Linux tools for troubleshooting
- **Proven Reliability**: Tested and verified with ChatSpatial

### Testing Results

!!! success "Native Deployment Verified (2025-08-28)"
    Testing confirms native deployment is the best approach:
    - ✅ **Native deployment**: Works perfectly with stdio MCP protocol
    - ✅ **All 16 ChatSpatial tools**: Fully accessible and functional
    - ✅ **Setup time**: Only 5 minutes
    - ❌ **Docker deployment**: Does NOT work due to stdio limitations

### Native Deployment Guide (Recommended)

!!! tip "Best Practice for Lab Servers"
    Based on our testing (2025-08-28), native deployment without Docker is the simplest and most reliable approach.
    Complete installation takes about 30 minutes.

#### Prerequisites

- Linux server (Ubuntu 20.04+ or CentOS 8+) or macOS
- Python 3.10+ (3.11 recommended)
- Node.js 18+ (for building AnythingLLM from source)
- 8GB+ RAM (ChatSpatial with dependencies uses ~1.5GB)
- 20GB+ disk space (dependencies are large)
- Internet connection for package downloads

#### Quick Start Script

```bash
#!/bin/bash
# ChatSpatial + AnythingLLM Lab Server Deployment
set -e

DEPLOY_DIR="/opt/spatial-platform"
PYTHON_VERSION="python3.11"

echo "Creating deployment directory..."
mkdir -p $DEPLOY_DIR
cd $DEPLOY_DIR

echo "Setting up Python environment..."
$PYTHON_VERSION -m venv venv
source venv/bin/activate

echo "Installing ChatSpatial (5-10 minutes)..."
pip install git+https://github.com/cafferychen777/ChatSpatial.git

echo "Creating MCP configuration..."
mkdir -p storage/plugins
cat > storage/plugins/anythingllm_mcp_servers.json << 'EOF'
{
  "mcpServers": {
    "chatspatial": {
      "command": "/opt/spatial-platform/venv/bin/python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
EOF

echo "Deployment complete! Now install AnythingLLM (see below)"
```

#### Detailed Step-by-Step Installation

##### 1. Install ChatSpatial

```bash
# Create working directory
sudo mkdir -p /opt/spatial-platform
cd /opt/spatial-platform
sudo chown -R $USER:$USER /opt/spatial-platform

# Create Python virtual environment
python3.11 -m venv venv
source venv/bin/activate

# Upgrade pip for better dependency resolution
pip install --upgrade pip

# Install ChatSpatial from source (not available on PyPI yet)
# Option 1: Install from local source (for development)
pip install -e /path/to/chatspatial

# Option 2: Install from GitHub (recommended)
pip install git+https://github.com/cafferychen777/ChatSpatial.git

# Verify installation (should show version number)
/opt/spatial-platform/venv/bin/python -c "from chatspatial import __version__; print(__version__)"
```

!!! note "Installation Time"
    ChatSpatial installation takes 5-10 minutes due to numerous scientific computing dependencies (100+ packages, ~1.5GB total)

##### 2. Install AnythingLLM

**Option A: Pre-built Binary (If Available)**

```bash
# For Linux x64
wget https://s3.us-west-1.amazonaws.com/public.useanything.com/latest/anythingllm-linux-x64.tar.gz
tar -xzf anythingllm-linux-x64.tar.gz
chmod +x anythingllm-linux-x64

# Create storage directory
mkdir -p /opt/spatial-platform/storage
```

**Option B: Build from Source (Recommended if binary unavailable)**

```bash
# Install Node.js dependencies
npm install -g yarn

# Clone AnythingLLM source
git clone https://github.com/Mintplex-Labs/anything-llm.git anythingllm-source
cd anythingllm-source

# Setup environment files
yarn setup:envs

# Configure storage path
echo "STORAGE_DIR=/opt/spatial-platform/storage" >> server/.env.development

# Install dependencies
cd server && yarn install
cd ../frontend && yarn install
cd ..

# Initialize database
yarn prisma:generate
yarn prisma:migrate

# Start services (in separate terminals)
# Terminal 1:
cd server && yarn dev  # Backend on port 3001

# Terminal 2:
cd frontend && yarn dev  # Frontend on port 3000
```

##### 3. Configure MCP Servers

```bash
# For pre-built binary installation:
mkdir -p /opt/spatial-platform/storage/plugins
cat > /opt/spatial-platform/storage/plugins/anythingllm_mcp_servers.json << 'EOF'
{
  "mcpServers": {
    "chatspatial": {
      "command": "/opt/spatial-platform/venv/bin/python",
      "args": ["-m", "chatspatial"],
      "cwd": "/opt/spatial-platform",
      "transport": "stdio"
    }
  }
}
EOF

# For source code installation:
mkdir -p /opt/spatial-platform/anythingllm-source/server/storage/plugins
cat > /opt/spatial-platform/anythingllm-source/server/storage/plugins/anythingllm_mcp_servers.json << 'EOF'
{
  "mcpServers": {
    "chatspatial": {
      "command": "/opt/spatial-platform/venv/bin/python",
      "args": ["-m", "chatspatial"],
      "cwd": "/opt/spatial-platform",
      "transport": "stdio"
    }
  }
}
EOF
```

!!! important "Critical Configuration Points"
    - **Correct file name**: Must be `anythingllm_mcp_servers.json` (NOT `mcp_servers.json`)
    - **Correct location**: Must be in `storage/plugins/` directory
    - **Must use absolute paths**: Relative paths will not work
    - **Working directory (cwd)**: Must specify where ChatSpatial should run from
    - **Transport type**: Explicitly set to `"stdio"` for clarity
    - **Python environment**: Use virtual environment's Python binary directly
    - **No 'server' argument**: Just use `["-m", "chatspatial"]` without additional arguments
    - **Development vs Production**: Source builds use `server/storage/plugins/`, binaries use `STORAGE_DIR/plugins/`

##### 4. Create Systemd Service

```bash
# Create service file
sudo tee /etc/systemd/system/anythingllm.service > /dev/null << 'EOF'
[Unit]
Description=AnythingLLM Service
After=network.target

[Service]
Type=simple
WorkingDirectory=/opt/spatial-platform
Environment="STORAGE_DIR=/opt/spatial-platform/storage"
ExecStart=/opt/spatial-platform/anythingllm-linux-x64
Restart=always
User=root

[Install]
WantedBy=multi-user.target
EOF

# Enable and start service
sudo systemctl daemon-reload
sudo systemctl enable anythingllm
sudo systemctl start anythingllm
```

##### 5. Configure Nginx (Optional)

For HTTPS access:

```bash
# Install Nginx
sudo apt install nginx

# Configure reverse proxy
sudo tee /etc/nginx/sites-available/anythingllm > /dev/null << 'EOF'
server {
    listen 80;
    server_name spatial.your-lab.edu;
    
    location / {
        proxy_pass http://localhost:3001;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
    }
}
EOF

# Enable site
sudo ln -s /etc/nginx/sites-available/anythingllm /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl reload nginx
```

##### 6. Access and Configure

1. **Access Web Interface**: 
   - Built from source: http://your-server:3000 (frontend)
   - Pre-built binary: http://your-server:3001

2. **Initial Setup**:
   - First user to register becomes administrator
   - Set a secure password

3. **Configure LLM Provider** (Settings → LLM):
   - Anthropic Claude (recommended)
   - OpenAI GPT
   - Local Ollama
   - Custom endpoints

4. **Enable MCP Hypervisor** (Settings → Agent Skills):
   - Toggle "MCP Hypervisor" to ON
   - ChatSpatial tools will be available as `@@mcp_chatspatial`

5. **Test ChatSpatial Integration**:
   - Start a new chat
   - Type: "List all available ChatSpatial tools"
   - You should see all 16 tools listed

6. **Verify Tool Functionality**:
   ```
   "Load example spatial transcriptomics data"
   "Preprocess the data and show quality metrics"
   "Identify spatial domains using SpaGCN"
   ```

#### Verification

```bash
# Check service status (if using systemd)
sudo systemctl status anythingllm

# Test ChatSpatial directly
/opt/spatial-platform/venv/bin/python -c "from chatspatial import __version__; print(f'ChatSpatial version: {__version__}')"

# Verify MCP configuration
cat /opt/spatial-platform/storage/plugins/anythingllm_mcp_servers.json  # For binary
# OR
cat /opt/spatial-platform/anythingllm-source/server/storage/plugins/anythingllm_mcp_servers.json  # For source

# Check if AnythingLLM can access ChatSpatial
curl http://localhost:3001/api/system/mcp-servers

# View logs
# For systemd:
journalctl -u anythingllm -f
# For development:
tail -f /opt/spatial-platform/anythingllm-source/server/logs/*.log
```

#### Production Deployment Tips

1. **Use a Process Manager**:
   ```bash
   # Install PM2
   npm install -g pm2
   
   # Start AnythingLLM with PM2
   pm2 start /opt/spatial-platform/anythingllm-linux-x64 --name anythingllm
   pm2 save
   pm2 startup
   ```

2. **Set Resource Limits**:
   ```bash
   # In systemd service file
   LimitNOFILE=65536
   LimitMEMLOCK=infinity
   ```

3. **Enable Logging**:
   ```bash
   # Configure detailed logging
   export DEBUG=mcp:*
   export LOG_LEVEL=debug
   ```

#### Troubleshooting

| Issue | Solution |
|-------|----------|
| MCP tools not appearing | 1. Check file name is `anythingllm_mcp_servers.json`<br>2. Verify it's in `storage/plugins/` directory<br>3. Ensure absolute paths in the config<br>4. Add 'server' argument to ChatSpatial args<br>5. Refresh MCP servers page in web UI |
| Service won't start | 1. Check permissions: `sudo chown -R $USER /opt/spatial-platform`<br>2. Verify Python path exists<br>3. Check port 3001 is not in use |
| Cannot access web UI | 1. Check firewall: `sudo ufw allow 3001`<br>2. Verify service is running: `ps aux | grep anythingllm`<br>3. Try accessing locally first: `curl http://localhost:3001` |
| ChatSpatial not found | 1. Activate venv: `source /opt/spatial-platform/venv/bin/activate`<br>2. Reinstall: `pip install git+https://github.com/cafferychen777/ChatSpatial.git`<br>3. Verify: `python -m chatspatial --version` |
| Tools fail to execute | 1. Check ChatSpatial logs<br>2. Ensure sufficient memory (8GB+ recommended)<br>3. Verify all dependencies installed correctly |
| Database errors | Run migrations: `cd anythingllm-source && yarn prisma:migrate` |

## Alternative Solution 2: Remote Development with Cursor/Windsurf

!!! success "Recommended Working Solution"
    This approach uses SSH remote development with AI-powered IDEs to access server resources while maintaining a local-like development experience. This method is currently the most reliable way to use ChatSpatial with server resources.

### Overview

Instead of deploying a web interface, this approach uses modern AI-powered IDEs (Cursor or Windsurf) to connect to your server via SSH. Once connected, you can install and use ChatSpatial's MCP server as if it were running locally, but with access to the server's computational resources.

### Architecture

```
Developer Machine (Cursor/Windsurf)
         ↓ (SSH)
    Server with ChatSpatial
         ├── MCP Server (stdio)
         ├── Python Environment
         └── Server GPU/CPU Resources
```

### Setup Instructions

#### Step 1: Server Preparation

Connect to your server and install ChatSpatial:

```bash
# Install Python 3.11 if not available
sudo apt update
sudo apt install python3.11 python3.11-venv

# Create virtual environment
python3.11 -m venv ~/chatspatial-env

# Install ChatSpatial
~/chatspatial-env/bin/pip install chatspatial
```

#### Step 2: Configure SSH Access

Ensure your SSH configuration is optimal for remote development:

```bash
# On your local machine, edit ~/.ssh/config
Host research-server
    HostName your.server.address
    User your_username
    ForwardAgent yes
    ServerAliveInterval 60
    ServerAliveCountMax 3
```

#### Step 3: Connect with Cursor/Windsurf

1. **Install Cursor or Windsurf** on your local machine
2. **Open Remote Connection**:
   - Cursor: Cmd/Ctrl+Shift+P → "Remote-SSH: Connect to Host"
   - Windsurf: Similar SSH remote development feature
3. **Select your server** from the SSH hosts list

#### Step 4: Configure MCP

Once connected to the server through your IDE, configure MCP:

**For Cursor:**
Create `~/.cursor/config.json` on the server:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "~/chatspatial-env/bin/python",
      "args": ["-m", "chatspatial"]
    }
  }
}
```

**For Windsurf:**
Create `~/.windsurf/config.json` with similar configuration.

!!! info "No API Key Required"
    Cursor and Windsurf have their own built-in LLM integration. You don't need to configure API keys on the server - the IDE handles all LLM interactions using its own credentials.

### Advantages of This Approach

1. **Full Server Resources**: Direct access to server's GPU, CPU, and storage
2. **Native MCP Support**: Works exactly like local MCP installation
3. **No Web Interface Issues**: Avoids all Open WebUI compatibility problems
4. **Familiar Development Environment**: Use your preferred AI-powered IDE
5. **Real-time Collaboration**: Can work on the same server with teammates
6. **Integrated Terminal**: Run computational tasks directly on the server
7. **No API Key Management**: IDE handles LLM credentials automatically
8. **Cost Efficiency**: Team uses their IDE's included LLM quota

### Usage Example

After setup, you can use ChatSpatial normally within your IDE:

```python
# In Cursor/Windsurf chat:
"Load the spatial transcriptomics data at /data/visium/sample1"
"Preprocess the data and identify spatial domains"
"Find spatially variable genes using GASTON"
```

The AI assistant will use the ChatSpatial MCP tools running on the server, utilizing the server's resources for computation while providing results in your local IDE.

### Team Deployment

For team usage:

1. Each team member gets their own user account on the server

2. Each user installs their preferred IDE (Cursor/Windsurf) locally

3. All team members connect to the same server via SSH

4. Data and results are shared on the server's filesystem

### Security Considerations

- Use SSH keys instead of passwords
- Configure firewall to only allow SSH from trusted IPs
- Set up separate user accounts for team members
- Use file permissions to control data access

## Cost Optimization

### API Usage Management

Monitor and control LLM API costs:

1. **Set User Quotas**: Limit API calls per user
2. **Cache Responses**: Store common analysis results
3. **Use Appropriate Models**: Choose cost-effective models for simple tasks
4. **Monitor Usage**: Track API consumption regularly

### Estimated Costs

| Component | Monthly Cost |
|-----------|--------------|
| Server hosting | $100-500 |
| Claude API | $500-2000* |
| GPT-5 API | $500-2000* |
| Storage | $50-200 |
| Backup | $20-100 |

*Depends on usage volume

## Performance Tuning

### System Optimization

```bash
# Optimize kernel parameters
sudo tee -a /etc/sysctl.conf << 'EOF'
# Network optimizations
net.core.rmem_max = 134217728
net.core.wmem_max = 134217728

# Memory optimizations
vm.swappiness = 10
vm.dirty_ratio = 15
EOF

sudo sysctl -p
```

### Container Optimization

Configure Docker for better performance:

```json
{
  "storage-driver": "overlay2",
  "log-driver": "json-file",
  "log-opts": {
    "max-size": "10m",
    "max-file": "3"
  }
}
```

## Next Steps

- Review [API Reference](../reference/api/README.md) for tool details
- Explore [Tutorials](../tutorials/README.md) for analysis workflows
- Check [Troubleshooting Guide](../reference/troubleshooting/common_issues.md) for solutions
