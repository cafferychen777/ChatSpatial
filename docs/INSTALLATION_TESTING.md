# ChatSpatial Installation Testing

## GitHub Action Workflow

**File**: `.github/workflows/test-installation.yml`

This workflow automatically tests ChatSpatial's MCP installation on every PR and push to main.

---

## What Gets Tested

### 1. Basic Installation (`pip install -e .`)
- **Platforms**: Python 3.10, 3.11, 3.12
- **Tests**:
  - ✅ Package installation succeeds
  - ✅ Core imports work (`import chatspatial`)
  - ✅ MCP server imports successfully
  - ✅ CLI commands are accessible

**Use Case**: Users who want minimal dependencies for basic spatial analysis

### 2. Full Installation (`pip install -e ".[full]"`)
- **Platforms**: Python 3.10, 3.11
- **Tests**:
  - ✅ All dependencies install correctly (PyTorch, scvi-tools, etc.)
  - ✅ Advanced features available
  - ✅ MCP server starts with full capabilities
  - ✅ Key packages verified (torch, scvi, tangram, etc.)

**Use Case**: Power users wanting all 60+ analysis methods

### 3. Dev Installation (`pip install -e ".[dev]"`)
- **Platform**: Python 3.10
- **Tests**:
  - ✅ Development tools installed (pytest, black, mypy, ruff)
  - ✅ Code formatting tools work
  - ✅ Testing framework ready

**Use Case**: Contributors developing ChatSpatial

### 4. MCP Protocol Communication
- **Tests**:
  - ✅ MCP server initializes correctly
  - ✅ Tools are registered (load_data, preprocess_data, visualize_data, etc.)
  - ✅ CLI interface functional
  - ✅ Server can be invoked

**Use Case**: Ensure MCP protocol integration works

---

## Workflow Triggers

The installation test runs on:

### Pull Requests
```yaml
paths:
  - 'pyproject.toml'          # Dependency changes
  - 'chatspatial/**'          # Code changes
  - '.github/workflows/test-installation.yml'
```

### Push to Main
- Same path filters as PRs
- Validates production deployments

### Manual Dispatch
- Can be triggered manually from GitHub Actions tab
- Useful for testing before releases

---

## Test Results

### Success Criteria

**Critical Tests** (Must Pass):
- ✅ Basic installation on all Python versions
- ✅ MCP protocol communication test

**Important Tests** (Should Pass):
- ✅ Full installation (may fail on dependency issues)
- ✅ Dev installation

### Viewing Results

1. **GitHub UI**:
   - Go to repository → Actions tab
   - Click on "Test MCP Installation" workflow
   - View detailed logs for each job

2. **PR Status Checks**:
   - Installation test appears as status check on PRs
   - Must pass before merge (if configured)

3. **Summary Output**:
   ```
   📊 Installation Test Results:
   ==============================
   ✅ Basic installation: PASSED
   ✅ Full installation: PASSED
   ✅ Dev installation: PASSED
   ✅ MCP protocol: PASSED
   ==============================
   ✅ All critical installation tests passed!
   ```

---

## Common Failure Scenarios

### 1. Dependency Conflict
**Symptom**: `pip install -e ".[full]"` fails with version conflicts

**Fix**:
```bash
# Update pyproject.toml with compatible versions
# Test locally first:
pip install -e ".[full]" --dry-run
```

### 2. Import Error
**Symptom**: `import chatspatial` fails in test

**Fix**: Check for missing dependencies in `pyproject.toml` core requirements

### 3. MCP Protocol Failure
**Symptom**: MCP tools not registered

**Fix**:
- Verify `server.py` has all tools decorated with `@mcp.tool()`
- Check for import errors in tool modules

### 4. CLI Not Accessible
**Symptom**: `python -m chatspatial --help` fails

**Fix**: Ensure `__main__.py` is properly configured

---

## Local Testing

Replicate GitHub Actions tests locally:

### Test Basic Installation
```bash
# Clean environment
python3.10 -m venv test_env
source test_env/bin/activate

# Install and test
pip install --upgrade pip
pip install -e .

# Verify
python -c "import chatspatial; print(chatspatial.__version__)"
python -m chatspatial server --help
```

### Test Full Installation
```bash
# Clean environment
python3.10 -m venv test_full_env
source test_full_env/bin/activate

# Install and test (takes ~5-10 min)
pip install --upgrade pip
pip install -e ".[full]"

# Verify advanced features
python -c "import torch; import scvi; import chatspatial"
```

### Test MCP Protocol
```bash
# Run the test script
cat << 'EOF' > test_mcp_local.py
from chatspatial.server import mcp

# Check tools registered
tools = mcp._tools if hasattr(mcp, '_tools') else {}
print(f"MCP server has {len(tools)} tools registered")

# List some tools
for tool_name in list(tools.keys())[:5]:
    print(f"  - {tool_name}")
EOF

python test_mcp_local.py
```

---

## Maintenance

### Adding New Tests

Edit `.github/workflows/test-installation.yml`:

```yaml
- name: Test new feature
  run: |
    python -c "import chatspatial.new_module"
    echo "✅ New feature works"
```

### Updating Python Versions

When dropping/adding Python version support:

1. Update `matrix.python-version` in workflow
2. Update README.md requirements
3. Update pyproject.toml `requires-python`

### Optimizing Test Speed

**Current optimizations**:
- Pip caching per Python version
- Fail-fast disabled (run all tests even if one fails)
- Parallel job execution
- Basic install: ~2-3 minutes
- Full install: ~8-15 minutes

**Future optimizations**:
- Pre-built Docker images
- Artifact caching for compiled dependencies
- Conditional full install (only on dependency changes)

---

## CI/CD Integration

### Branch Protection

Recommend requiring this check before merge:

```
Settings → Branches → Branch protection rules → main
☑ Require status checks to pass before merging
  ☑ Test MCP Installation / installation-summary
```

### Release Process

Before creating a release:

1. ✅ All installation tests pass on main
2. ✅ Manual test on fresh environment
3. ✅ Test with Claude Desktop integration
4. ✅ Create GitHub release

---

## Badge for README

Add this badge to show installation test status:

```markdown
[![Installation Tests](https://github.com/cafferychen777/ChatSpatial/workflows/Test%20MCP%20Installation/badge.svg)](https://github.com/cafferychen777/ChatSpatial/actions/workflows/test-installation.yml)
```

---

**Last Updated**: 2025-10-21
**Status**: ✅ Active and Monitoring
