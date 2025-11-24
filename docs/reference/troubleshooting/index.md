# Troubleshooting Guide

Solutions to common issues and frequently asked questions.

## Quick Links

- [Common Issues](common_issues.md) - Solutions to frequently encountered problems
- [FAQ](faq.md) - Frequently asked questions and answers

## Getting Help

If you can't find a solution here:

1. **Check the documentation** - Search this site for relevant information
2. **Review error messages** - They often contain helpful diagnostic information
3. **Search GitHub issues** - Someone may have encountered the same problem
4. **Ask the community** - Post in [GitHub Discussions](https://github.com/cafferychen777/ChatSpatial/discussions)
5. **Report bugs** - Create an [issue on GitHub](https://github.com/cafferychen777/ChatSpatial/issues)

## Common Problem Categories

### Installation Issues

- Python version compatibility
- Missing dependencies
- Virtual environment setup
- MCP server configuration

See: [Common Issues - Installation](common_issues.md#installation)

### Data Loading Problems

- Unsupported file formats
- Corrupted or incomplete data
- Memory issues with large datasets
- Missing spatial coordinates

See: [Common Issues - Data Loading](common_issues.md#data-loading)

### Analysis Errors

- Parameter validation failures
- Method-specific errors
- Insufficient data quality
- Computational resource limitations

See: [Common Issues - Analysis](common_issues.md#analysis)

### MCP Connection Issues

- Claude Desktop not connecting
- Tool calls failing
- Server crashes or timeouts
- Communication protocol errors

See: [Common Issues - MCP](common_issues.md#mcp-connection)

## Best Practices for Troubleshooting

1. **Read error messages carefully** - They usually indicate what went wrong
2. **Check prerequisites** - Ensure data is properly preprocessed
3. **Verify parameters** - Use default values when unsure
4. **Start simple** - Test with small datasets first
5. **Keep logs** - Save error messages for reference

## Reporting Issues

When reporting a bug, please include:

- **ChatSpatial version** - Run `pip show chatspatial`
- **Python version** - Run `python --version`
- **Operating system** - macOS, Linux, Windows
- **Error message** - Complete error text
- **Steps to reproduce** - What you did to trigger the error
- **Expected behavior** - What should have happened
- **Data information** - Dataset type, size, format (if applicable)

## Additional Resources

- [Installation Guide](../../getting-started/installation.md)
- [Quick Reference](../quick-reference/index.md)
- [API Documentation](../api/index.md)
- [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues)
