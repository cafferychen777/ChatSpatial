# ChatSpatial Installation

## Quick Start (2 minutes)

### Step 1: Create Virtual Environment (Strongly Recommended)

```bash
# Clone repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd chatspatial

# Create virtual environment (choose one method)
# Option A: Using venv (Python built-in)
python3 -m venv chatspatial_env
source chatspatial_env/bin/activate  # On macOS/Linux
# chatspatial_env\Scripts\activate    # On Windows

# Option B: Using conda
conda create -n chatspatial python=3.10
conda activate chatspatial

# Option C: Using uv (fast, modern)
uv venv
source .venv/bin/activate
```

### Step 2: Install ChatSpatial

```bash
# Recommended: Install with all features
pip install -e ".[full]"
```

> 💡 **Why [full]?** Enables all 16+ analysis methods including Cell2location, CellRank, SpaGCN, and more. Takes ~13 minutes but worth it for the complete experience.

### Step 3: Configure Your MCP Client

Choose your client below:

<a name="claude-desktop"></a>
#### Option A: Claude Desktop (GUI Application)

**Important:** Use your virtual environment's Python path in the configuration.

Add to Claude Desktop config file:
- **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
- **Linux**: `~/.config/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/your/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial"]
    }
  }
}
```

**Example with actual path:**
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/Users/apple/Research/SpatialTrans_MCP/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial"]
    }
  }
}
```

<a name="claude-code"></a>
#### Option B: Claude Code (Terminal/IDE)

**Step 1: Install Claude Code CLI**
```bash
# Install via npm (recommended)
npm install -g @anthropic-ai/claude-code

# Verify installation
claude --version
```

**Step 2: Add ChatSpatial MCP Server**
```bash
# Get your virtual environment Python path
which python  # Should show: /path/to/chatspatial_env/bin/python

# Add ChatSpatial (local scope - current project only)
claude mcp add chatspatial /path/to/chatspatial_env/bin/python -- -m chatspatial

# OR add globally (all projects)
claude mcp add --scope user chatspatial /path/to/chatspatial_env/bin/python -- -m chatspatial
```

**Step 3: Verify Configuration**
```bash
# List all MCP servers
claude mcp list

# Check ChatSpatial details
claude mcp get chatspatial
```

> 💡 **Note:** The `--` separates Claude's options from the server command. Everything after `--` is the command to run ChatSpatial.

### Step 4: Restart Your Client

- **Claude Desktop**: Restart the application
- **Claude Code**: The server is ready immediately (use `/mcp` in Claude Code to check status)

That's it! Start analyzing your spatial data with natural language.

## Installation Options

| Option | Install Command | Features | Install Time |
|--------|----------------|----------|--------------|
| **Full (Recommended)** | `pip install -e ".[full]"` | 100% features | ~13 minutes |
| Standard | `pip install -e .` | 80% features | ~6 minutes |

### What's included in each option?

**Standard Installation** (Default):
- ✅ Core spatial analysis (Moran's I, Getis-Ord, etc.)
- ✅ Basic deconvolution with scvi-tools
- ✅ RNA velocity with scVelo
- ✅ Cell communication (LIANA, CellPhoneDB)
- ✅ Batch integration (Harmony, BBKNN)
- ✅ Gene enrichment analysis

**Full Installation** (Recommended):
- ✅ Everything in Standard, plus:
- ✅ Deep learning methods (PyTorch)
- ✅ Advanced deconvolution (Cell2location)
- ✅ Advanced trajectory (CellRank, Palantir)
- ✅ Spatial domains (SpaGCN, STAGATE)
- ✅ Spatial variable genes (GASTON, SpatialDE)
- ✅ R-based methods (RCTD, Spotlight)

## Requirements

- Python 3.8-3.12
- macOS, Linux, or Windows
- 5-10 GB disk space (for full installation)

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Import errors | Update pip: `pip install --upgrade pip` |
| Package conflicts | `pip install --force-reinstall -e ".[full]"` |
| Claude Desktop doesn't see server | Check that command path points to virtual environment Python |
| Claude Code connection error | Ensure absolute path to Python; run `claude mcp list` to verify |
| "python not found" error | Use absolute path to virtual environment Python in config |
| Virtual environment issues | Recreate environment: `rm -rf chatspatial_env && python3 -m venv chatspatial_env` |

## Verify Installation

```bash
# Make sure you're in the virtual environment
which python  # Should show path to your virtual environment

# Check if ChatSpatial is installed
python -m chatspatial --help

# Test in Python
python -c "import chatspatial; print('✅ Installation successful')"
```

## Getting Help

- **Issues**: [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
- **Documentation**: Check docstrings with `help(function_name)`
- **Community**: Ask questions in Discussions

## Advanced Configuration

### Claude Code Scopes

Claude Code supports three configuration scopes:

| Scope | Command Flag | Visibility | Use Case |
|-------|-------------|------------|----------|
| **Local** | (default) | Current project only | Personal projects, experiments |
| **Project** | `--scope project` | Team via `.mcp.json` | Shared team configuration |
| **User** | `--scope user` | All your projects | Personal tools across projects |

**Examples:**
```bash
# Local scope (default)
claude mcp add chatspatial /path/to/venv/bin/python -- -m chatspatial

# Project scope (creates .mcp.json for team sharing)
claude mcp add --scope project chatspatial /path/to/venv/bin/python -- -m chatspatial

# User scope (available in all projects)
claude mcp add --scope user chatspatial /path/to/venv/bin/python -- -m chatspatial
```

### Managing MCP Servers

```bash
# List all servers
claude mcp list

# Remove a server
claude mcp remove chatspatial

# Check server details
claude mcp get chatspatial

# Check status in Claude Code
/mcp  # Type this in Claude Code chat
```

### Environment Variables

```bash
# Add with custom environment variables
claude mcp add chatspatial \
  --env CHATSPATIAL_DATA_DIR=/path/to/data \
  --env CHATSPATIAL_CACHE_DIR=/path/to/cache \
  /path/to/venv/bin/python -- -m chatspatial
```

### Windows-Specific Notes

For Windows users not using WSL:
```cmd
# Native Windows path example
claude mcp add chatspatial ^
  C:\Users\user\chatspatial_env\Scripts\python.exe ^
  -- -m chatspatial
```