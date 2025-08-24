# ChatSpatial Project Structure

## Root Directory Overview

```
chatspatial/
├── README.md                    # Main project documentation
├── CHANGELOG.md                 # Version history and changes
├── CONTRIBUTING.md              # Contribution guidelines
├── LICENSE                      # Project license
├── SECURITY.md                  # Security policy
├── CLAUDE.md                    # Claude AI configuration and instructions
├── GEMINI.md                    # Gemini AI configuration
├── llms-full.txt                # LLM reference and configuration data
├── pyproject.toml               # Python project configuration
├── setup.cfg                    # Setup configuration
├── pytest.ini                  # Testing configuration
├── requirements-dev.txt         # Development dependencies
├── requirements-full.txt        # Full production dependencies
├── chatspatial/                 # Main package source code
├── tests/                       # Test suite
├── scripts/                     # Utility scripts and tools
├── docs/                        # Documentation
├── examples/                    # Usage examples
├── data/                        # All datasets organized by purpose
│   ├── core/                    # Core example datasets
│   ├── harmony/                 # Harmony integration test data
│   └── test/                    # Testing datasets
├── archived_tests/              # Historical test scripts (organized)
├── archived_docs/               # Historical documentation
├── paper/                       # Research paper and figures
└── third_party/                 # External dependencies and integrations
```

## Main Package Structure (`chatspatial/`)

```
chatspatial/
├── __init__.py                  # Package initialization
├── __main__.py                  # CLI entry point
├── config.py                    # Configuration management
├── server.py                    # Main MCP server
├── http_server.py               # HTTP transport server
├── spatial_mcp_adapter.py       # Spatial data MCP adapter
├── mcp/                         # MCP-specific modules
│   ├── annotations.py           # MCP annotations
│   ├── errors.py                # MCP error handling
│   ├── prompts.py               # MCP prompts
│   └── resources.py             # MCP resources
├── models/                      # Data models
│   ├── analysis.py              # Analysis result models
│   └── data.py                  # Data parameter models
├── tools/                       # Analysis tools
│   ├── annotation.py            # Cell type annotation
│   ├── cell_communication.py    # Cell-cell communication
│   ├── deconvolution.py         # Spot deconvolution
│   ├── differential.py          # Differential expression
│   ├── integration.py           # Data integration
│   ├── pathway_enrichment.py    # Pathway analysis
│   ├── preprocessing.py         # Data preprocessing
│   ├── spatial_analysis.py      # Spatial analysis
│   ├── spatial_domains.py       # Spatial domain identification
│   ├── spatial_genes.py         # Spatially variable genes
│   ├── trajectory.py            # Trajectory analysis
│   └── visualization.py         # Data visualization
└── utils/                       # Utility functions
    ├── data_loader.py           # Data loading utilities
    ├── error_handling.py        # Data validation and processing errors
    ├── tool_error_handling.py   # MCP protocol error handling
    ├── image_utils.py           # Image processing utilities
    └── mcp_parameter_handler.py # MCP parameter handling
```

## Organized Archive Directories

### `archived_tests/` - Historical Test Scripts
```
archived_tests/
├── scanvi/                      # scANVI method tests
│   └── test_scanvi_comprehensive.py
├── dependency_validation/       # Dependency validation tests
│   └── test_dependency_validation.py
├── tangram/                     # Tangram method tests
│   ├── test_tangram_fix.py
│   └── debug_tangram.py
└── cellassign/                  # CellAssign method tests
    └── test_cellassign_comprehensive.py
```

### `archived_docs/` - Historical Documentation
```
archived_docs/
├── reports/                     # Historical test and analysis reports
│   ├── SCANVI_TEST_REPORT.md
│   ├── TANGRAM_IMPROVEMENT_PLAN.md
│   ├── TANGRAM_TESTING_REPORT.md
│   └── cellassign_test_results.json
└── historical/                  # Historical reference files
    ├── CLAUDE.md                # Historical Claude configuration
    ├── GEMINI.md                # Historical Gemini documentation
    └── llms-full.txt            # Historical LLM reference
```

## Key Directories

### `/tests/` - Active Test Suite
- **Unit tests**: Individual component testing
- **Integration tests**: Multi-component workflow testing
- **E2E tests**: End-to-end MCP integration testing
- **Stress tests**: Performance and reliability testing

### `/docs/` - Active Documentation
- **API documentation**: Tool and method documentation
- **User guides**: Usage instructions and tutorials
- **Development guides**: Contributing and development setup
- **Reports**: Current analysis and status reports

### `/scripts/` - Utility Scripts
- **Development tools**: Testing and debugging scripts
- **Setup scripts**: Environment and dependency setup
- **Data preparation**: Dataset creation and formatting
- **Maintenance**: Project cleanup and organization

### `/examples/` - Usage Examples
- **Basic usage**: Simple tool demonstrations
- **Advanced workflows**: Complex analysis pipelines
- **MCP integration**: Server setup and client examples
- **Visualization demos**: Plotting and figure generation

## File Organization Principles

1. **Clean Root**: Essential project files and active AI configurations in root directory
2. **Logical Grouping**: Related files organized in appropriate subdirectories  
3. **Clear Naming**: Descriptive file and directory names
4. **Archive Separation**: Historical files separated from active codebase
5. **Type-based Organization**: Tests, docs, examples in dedicated directories
6. **Active AI Configs**: AI configuration files (CLAUDE.md, GEMINI.md, llms-full.txt) remain in root for easy access

## Getting Started

1. **Installation**: See `README.md` for setup instructions
2. **Development**: See `CONTRIBUTING.md` for development guidelines
3. **Testing**: Use `pytest` from root directory
4. **Examples**: Check `/examples/` for usage demonstrations
5. **Documentation**: Browse `/docs/` for detailed guides

## Maintenance Notes

- Archived directories contain historical files for reference
- Active development should focus on main package and `/tests/`
- Regular cleanup moves temporary files to appropriate archives
- Documentation is kept current in `/docs/` directory