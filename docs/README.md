# ChatSpatial Documentation

Welcome to the ChatSpatial documentation! This directory contains all user guides, technical documentation, and reference materials for the ChatSpatial MCP server.

## 🚀 Getting Started

### New Users
1. **[INSTALLATION.md](INSTALLATION.md)** - Complete installation guide
2. **[QUICK_START.md](QUICK_START.md)** - Get up and running in minutes
3. **[API_REFERENCE.md](API_REFERENCE.md)** - Complete API documentation for all 16 MCP tools
4. **[user_guides/COMPREHENSIVE_USER_GUIDE.md](user_guides/COMPREHENSIVE_USER_GUIDE.md)** - Detailed usage guide

### Developers
1. **[CRITICAL_IMAGE_DISPLAY_BUG.md](CRITICAL_IMAGE_DISPLAY_BUG.md)** - ⚠️ **MUST READ** before modifying visualization code
2. **[technical_docs/](technical_docs/)** - MCP protocol and technical specifications
3. **[modules/](modules/)** - API documentation for tools, models, and utilities

## 📁 Documentation Structure

### Core Documentation
- **[INSTALLATION.md](INSTALLATION.md)** - Installation and setup instructions
- **[QUICK_START.md](QUICK_START.md)** - Quick start guide and basic workflows  
- **[API_REFERENCE.md](API_REFERENCE.md)** - Complete API documentation for all MCP tools
- **[DEPENDENCY_WARNINGS.md](DEPENDENCY_WARNINGS.md)** - Known issues and dependency warnings

### User Guides
- **[user_guides/](user_guides/)** - Comprehensive user documentation
  - **COMPREHENSIVE_USER_GUIDE.md** - Complete user manual
  - **ERROR_HANDLING_GUIDE.md** - Error resolution strategies
  - **LIANA_COMPREHENSIVE_GUIDE.md** - Cell communication analysis guide
  - **SPATIAL_STATISTICS_USAGE.md** - Spatial statistics workflows
  - **FIND_MARKERS_USAGE.md** - Marker gene identification guide
  - And more...

### Technical Documentation
- **[technical_docs/](technical_docs/)** - Technical specifications and guides
  - **MCP_SERVER_SPECIFICATION.md** - MCP protocol implementation details
  - **MCP_TOOLS_QUICK_REFERENCE.md** - Tool usage reference
  - **HTTP_TRANSPORT.md** - HTTP transport configuration
  - **MCP_ERROR_HANDLING_GUIDE.md** - Error handling implementation
  - **ENRICHMAP_INTEGRATION.md** - EnrichMap integration guide

### API Documentation
- **[modules/](modules/)** - Module-specific documentation
  - **[tools/](modules/tools/)** - Individual tool documentation
  - **[models/](modules/models/)** - Data model documentation
  - **[utils/](modules/utils/)** - Utility function documentation

### Reference Materials
- **[reference/](reference/)** - Reference files and lookup tables

### Archive
- **[archive/](archive/)** - Historical documentation and development reports
  - See **[archive/README.md](archive/README.md)** for navigation guide

## 🚨 Critical Information

### For Developers
**⚠️ MUST READ: [CRITICAL_IMAGE_DISPLAY_BUG.md](CRITICAL_IMAGE_DISPLAY_BUG.md)**

This document describes a critical bug that took 2 weeks to identify and fix. Read this before modifying:
- Image handling in error decorators
- The `visualize_data` function  
- The `mcp_tool_error_handler` decorator
- Any FastMCP return value processing

### For Users
- Start with **[QUICK_START.md](QUICK_START.md)** for immediate usage
- Check **[DEPENDENCY_WARNINGS.md](DEPENDENCY_WARNINGS.md)** for known issues
- Use **[user_guides/ERROR_HANDLING_GUIDE.md](user_guides/ERROR_HANDLING_GUIDE.md)** for troubleshooting

## 🔍 Quick Navigation

| Need | Go To |
|------|-------|
| Install ChatSpatial | [INSTALLATION.md](INSTALLATION.md) |
| Learn basic usage | [QUICK_START.md](QUICK_START.md) |
| Complete API reference | [API_REFERENCE.md](API_REFERENCE.md) |
| Complete user guide | [user_guides/COMPREHENSIVE_USER_GUIDE.md](user_guides/COMPREHENSIVE_USER_GUIDE.md) |
| Fix errors | [user_guides/ERROR_HANDLING_GUIDE.md](user_guides/ERROR_HANDLING_GUIDE.md) |
| Tool reference | [technical_docs/MCP_TOOLS_QUICK_REFERENCE.md](technical_docs/MCP_TOOLS_QUICK_REFERENCE.md) |
| Module documentation | [modules/](modules/) |
| Development history | [archive/](archive/) |

## 📝 Contributing

When adding documentation:
1. **User-facing docs**: Add to [user_guides/](user_guides/)
2. **Technical specs**: Add to [technical_docs/](technical_docs/)
3. **API docs**: Add to [modules/](modules/)
4. **Historical records**: Archive in [archive/](archive/)

---

**Need help?** Check the [user_guides/ERROR_HANDLING_GUIDE.md](user_guides/ERROR_HANDLING_GUIDE.md) or create an issue on GitHub.