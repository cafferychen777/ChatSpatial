# Contributing to ChatSpatial

Thank you for your interest in contributing to ChatSpatial! This document provides guidelines for contributing to the project.

## Getting Started

### Prerequisites
- Python 3.8 or higher (Python 3.10-3.11 recommended for optimal compatibility)
- Git
- Basic understanding of spatial transcriptomics
- Familiarity with Model Context Protocol (MCP) is helpful
- For visualization testing: Understanding of matplotlib and image processing
- For advanced methods: Knowledge of PyTorch, R (for R-based methods)

### Development Setup

1. **Fork and clone the repository**
   ```bash
   git clone https://github.com/your-username/ChatSpatial.git
   cd ChatSpatial
   ```

2. **Create a virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install development dependencies**
   ```bash
   # Install ChatSpatial with all optional dependencies
   pip install -e .[all]
   
   # Install development tools
   pip install pytest pytest-cov black isort mypy flake8
   
   # Install MCP development tools
   pip install mcp[dev]
   
   # Optional: Install R dependencies for R-based methods
   pip install rpy2  # Requires R to be installed separately
   ```

4. **Run tests to ensure everything works**
   ```bash
   # Run basic test suite
   pytest tests/
   
   # Run with coverage reporting
   pytest tests/ --cov=chatspatial --cov-report=html
   
   # Run specific test categories
   pytest tests/unit/          # Unit tests
   pytest tests/integration/   # Integration tests
   pytest tests/e2e/          # End-to-end MCP tests
   
   # Test MCP server directly
   python -m chatspatial --help
   ```

## Project Architecture

### Core Components

ChatSpatial follows a modular MCP server architecture:

```
chatspatial/
├── server.py                    # Main MCP server with tool definitions
├── spatial_mcp_adapter.py       # Spatial data adapter for MCP
├── tools/                       # Analysis tool implementations
│   ├── annotation.py            # Cell type annotation methods
│   ├── cell_communication.py    # LIANA+, CellPhoneDB, CellChat
│   ├── spatial_genes.py         # GASTON, SpatialDE, SPARK
│   ├── spatial_domains.py       # SpaGCN, STAGATE, BANKSY
│   ├── visualization.py         # All plot types and image generation
│   └── ...
├── models/                      # Pydantic data models
│   ├── data.py                  # Parameter models for each tool
│   └── analysis.py              # Result models for analysis outputs
├── utils/                       # Utility functions
│   ├── error_handling.py        # Data validation and business logic errors
│   ├── tool_error_handling.py   # MCP protocol error handling
│   ├── data_loader.py           # Data loading utilities
│   └── image_utils.py           # Image processing and conversion
└── mcp/                         # MCP-specific modules
    ├── errors.py                # MCP error definitions
    ├── resources.py             # MCP resource management
    └── ...
```

### Key Design Principles

1. **Separation of Analysis and Visualization**: Each analysis tool focuses on computation, visualization is handled separately
2. **Two-Layer Error Handling**: Data validation errors vs. MCP protocol errors
3. **Pydantic Parameter Validation**: All tool parameters are validated using Pydantic models
4. **Image Object Handling**: Critical handling of matplotlib images for MCP display
5. **Resource Management**: Automatic cleanup and memory management for large datasets

## Development Workflow

### Branching Strategy
- `main`: Stable release branch (production-ready)
- `develop`: Development branch for new features (optional - direct to main is acceptable)
- Feature branches: `feature/your-feature-name`
- Bug fixes: `bugfix/issue-description`
- Documentation: `docs/update-description`

### Making Changes

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**
   - Follow the existing code style
   - Add tests for new functionality
   - Update documentation as needed

3. **Run tests and linting**
   ```bash
   # Run tests
   pytest tests/
   
   # Format code
   black .
   isort .
   
   # Check linting
   flake8 .
   mypy chatspatial
   ```

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "feat: add new spatial analysis method"
   ```

5. **Push and create a Pull Request**
   ```bash
   git push origin feature/your-feature-name
   ```

## Code Style Guidelines

### Python Code Style
- Follow PEP 8
- Use Black for code formatting
- Use isort for import sorting
- Maximum line length: 88 characters (Black default)
- Use type hints where appropriate

### Commit Message Format
Follow conventional commits:
- `feat:` for new features
- `fix:` for bug fixes
- `docs:` for documentation changes
- `test:` for adding tests
- `refactor:` for code refactoring
- `perf:` for performance improvements

### Documentation
- Update README.md for new features
- Add docstrings to all public functions and classes
- Include examples in docstrings
- Update API documentation

## Testing Guidelines

### Test Structure

ChatSpatial uses a comprehensive testing strategy:

```
tests/
├── unit/                        # Unit tests for individual components
│   ├── algorithms/             # Algorithm-specific tests
│   ├── models/                 # Pydantic model tests
│   └── utils/                  # Utility function tests
├── integration/                # Integration tests for workflows
│   ├── test_enrichment_analysis.py
│   └── test_trajectory_analysis.py
├── e2e/                        # End-to-end MCP integration tests
│   └── test_mcp_integration.py
├── stress_tests/               # Performance and reliability tests
│   └── test_comprehensive_stress.py
├── visualization_tests/        # Image generation and display tests
│   └── test_visualization_comprehensive.py
└── fixtures/                   # Shared test fixtures
    └── __init__.py
```

### Test Categories and When to Use

1. **Unit Tests**: Individual function/method testing
   ```python
   # Test single function behavior
   async def test_validate_adata():
       # Test data validation logic
   ```

2. **Integration Tests**: Multi-component workflows
   ```python
   # Test complete analysis workflows
   async def test_gaston_spatial_genes_workflow():
       # Test load → preprocess → analyze → visualize
   ```

3. **E2E Tests**: Full MCP protocol testing
   ```python
   # Test MCP tool calls end-to-end
   async def test_mcp_spatial_analysis_complete():
       # Test through MCP protocol
   ```

4. **Stress Tests**: Large dataset and performance testing
   ```python
   # Test with realistic dataset sizes
   async def test_large_dataset_performance():
       # Test memory usage and timeout handling
   ```

5. **Visualization Tests**: Image generation testing
   ```python
   # Test image creation and format
   async def test_spatial_plot_generation():
       # Test matplotlib → MCP Image object conversion
   ```

### Critical Testing Requirements

#### **Image Object Testing (CRITICAL)**

**⚠️ NEVER modify image handling code without extensive testing**

The image display functionality has a critical bug that took 2 weeks to resolve. When testing visualization:

```python
# CRITICAL: Test image object handling
async def test_image_object_handling():
    """Test that Image objects are handled correctly for MCP display."""
    
    # Test visualization tool
    result = await visualize_data("test_data", params)
    
    # CRITICAL: Verify Image object is returned correctly
    assert isinstance(result, Image)  # Must be raw Image object
    assert result.data  # Must contain base64 data
    assert result.mimeType == "image/png"  # Must be PNG format
    
    # CRITICAL: Test that wrapper functions don't break Image objects
    # See docs/CRITICAL_IMAGE_DISPLAY_BUG.md for details
```

#### **MCP Error Handling Testing**

Test both error handling layers:

```python
# Test data validation layer
async def test_data_validation_errors():
    with pytest.raises(ProcessingError):
        await analyze_with_invalid_data()

# Test MCP protocol layer  
async def test_mcp_error_formatting():
    result = await mcp_tool_with_error()
    assert result["isError"] == True
    assert "Error:" in result["content"][0]["text"]
```

### Test Data Management

#### **Synthetic Data Creation**

```python
# Create realistic synthetic data for testing
def create_synthetic_spatial_data(n_spots=200, n_genes=100):
    """Create synthetic spatial transcriptomics data."""
    
    # Create spatial coordinates
    spatial_coords = generate_spatial_pattern()
    
    # Create expression data with spatial patterns
    X = create_spatial_expression(n_spots, n_genes, spatial_coords)
    
    # Create AnnData with proper structure
    adata = sc.AnnData(X)
    adata.obsm['spatial'] = spatial_coords
    adata.obs['cell_type'] = assign_cell_types(spatial_coords)
    
    return adata
```

#### **Test Data Guidelines**

- **Small datasets**: <1000 spots, <500 genes for unit tests
- **Medium datasets**: 1000-3000 spots for integration tests  
- **Large datasets**: >3000 spots only for stress tests
- **Don't commit large files**: Use `.gitignore` for test data >10MB
- **Reproducible data**: Set random seeds for consistent test results

### Performance Testing

#### **Memory Usage Testing**

```python
import psutil
import os

async def test_memory_usage():
    """Test memory usage during analysis."""
    
    # Get initial memory
    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info().rss
    
    # Run analysis
    result = await analyze_large_dataset()
    
    # Check memory increase is reasonable
    final_memory = process.memory_info().rss
    memory_increase = (final_memory - initial_memory) / 1024 / 1024  # MB
    
    assert memory_increase < 500  # Less than 500MB increase
```

#### **Timeout Testing**

```python
import asyncio

async def test_analysis_timeout():
    """Test that long-running analysis completes within reasonable time."""
    
    start_time = time.time()
    
    # Run analysis with timeout
    try:
        result = await asyncio.wait_for(
            analyze_complex_data(), 
            timeout=300  # 5 minutes max
        )
        assert result.success
    except asyncio.TimeoutError:
        pytest.fail("Analysis took too long (>5 minutes)")
    
    end_time = time.time()
    assert end_time - start_time < 300
```

### Test Fixtures and Utilities

Create reusable test fixtures:

```python
# tests/fixtures/__init__.py
import pytest

@pytest.fixture
def sample_spatial_data():
    """Standard spatial transcriptomics test data."""
    return create_synthetic_spatial_data(n_spots=100, n_genes=50)

@pytest.fixture  
def mock_context():
    """Mock MCP context for testing."""
    class MockContext:
        async def info(self, msg): pass
        async def warning(self, msg): pass
    return MockContext()

@pytest.fixture
def data_store(sample_spatial_data):
    """Standard data store for testing."""
    return {
        "test_data": {
            "adata": sample_spatial_data,
            "name": "Test Dataset",
            "type": "synthetic"
        }
    }
```

### Test Automation and CI

#### **Pre-commit Testing**

```bash
# Test before committing
python -m pytest tests/unit/ -x          # Stop on first failure
python -m pytest tests/integration/ -v   # Verbose integration tests
python -m pytest tests/e2e/ --timeout=300  # E2E with timeout
```

#### **Coverage Requirements**

```bash
# Generate coverage report
pytest tests/ --cov=chatspatial --cov-report=html --cov-report=term

# Minimum coverage targets:
# - Unit tests: >90% coverage
# - Integration tests: >80% coverage
# - Overall: >85% coverage
```

#### **Test Performance Monitoring**

```bash
# Run with performance profiling
pytest tests/ --durations=10  # Show 10 slowest tests

# Memory profiling
pytest tests/ --memory-profiler  # If installed
```

### Debugging Test Failures

#### **Common Test Failure Patterns**

1. **Image Display Issues**: Check Image object handling
2. **Memory Errors**: Reduce test dataset size or increase system memory
3. **Timeout Errors**: Use Cherry Studio configuration for testing
4. **Import Errors**: Check optional dependencies are installed
5. **Data Format Errors**: Verify AnnData structure matches expectations

#### **Test Debugging Tools**

```python
# Add debug logging to tests
import logging
logging.basicConfig(level=logging.DEBUG)

# Use pytest debugging
pytest tests/test_your_tool.py -v -s --pdb  # Drop to debugger on failure

# Capture stdout/stderr
pytest tests/test_your_tool.py -s  # Don't capture output
```

## Adding New Features

### Adding New MCP Tools

Follow this pattern for adding new analysis tools:

#### 1. Create Parameter Model (`chatspatial/models/data.py`)

```python
from pydantic import BaseModel, Field
from typing import Optional, List

class YourAnalysisParameters(BaseModel):
    """Parameters for your new analysis method."""
    
    method: str = Field("default_method", description="Analysis method to use")
    parameter1: int = Field(10, description="Description of parameter1")
    parameter2: Optional[float] = Field(None, description="Optional parameter")
    genes: Optional[List[str]] = Field(None, description="Genes to analyze")
```

#### 2. Create Result Model (`chatspatial/models/analysis.py`)

```python
from pydantic import BaseModel
from typing import Optional, Dict, Any

class YourAnalysisResult(BaseModel):
    """Result from your analysis method."""
    
    method: str
    success: bool
    n_features_analyzed: int
    results_key: Optional[str] = None  # Key in adata.uns where results are stored
    statistics: Optional[Dict[str, Any]] = None
    message: str
```

#### 3. Implement Tool Logic (`chatspatial/tools/your_tool.py`)

```python
from chatspatial.utils.error_handling import validate_adata, ProcessingError
from chatspatial.utils.tool_error_handling import mcp_tool_error_handler
from chatspatial.models.data import YourAnalysisParameters
from chatspatial.models.analysis import YourAnalysisResult

async def your_analysis_function(
    data_id: str,
    data_store: dict,
    params: YourAnalysisParameters,
    context
) -> YourAnalysisResult:
    """
    Implement your analysis method.
    
    Args:
        data_id: Dataset identifier
        data_store: Data storage dictionary
        params: Analysis parameters
        context: MCP context for logging
        
    Returns:
        YourAnalysisResult: Analysis results
    """
    
    # Validate input data
    adata = data_store[data_id]["adata"]
    validate_adata(adata, {
        'obs': ['required_column'],  # Add required columns
        'obsm': ['spatial']          # Spatial coordinates usually required
    }, check_spatial=True)
    
    # Implement your analysis logic
    try:
        # Your analysis code here
        results = run_your_analysis(adata, params)
        
        # Store results in adata
        adata.uns['your_analysis_results'] = results
        
        return YourAnalysisResult(
            method=params.method,
            success=True,
            n_features_analyzed=len(results),
            results_key='your_analysis_results',
            message=f"Analysis completed successfully with {len(results)} features"
        )
        
    except Exception as e:
        raise ProcessingError(f"Analysis failed: {str(e)}")
```

#### 4. Register MCP Tool (`chatspatial/server.py`)

```python
from chatspatial.tools.your_tool import your_analysis_function

@mcp.tool()
@mcp_tool_error_handler()
async def your_analysis_tool(
    data_id: str,
    params: Optional[Dict[str, Any]] = None,
    context: Context = None
) -> YourAnalysisResult:
    """
    Brief description of your analysis tool.
    
    Args:
        data_id: Dataset identifier from load_data
        params: Analysis parameters (optional)
        
    Returns:
        Analysis results
    """
    
    # Parse parameters
    analysis_params = YourAnalysisParameters(**(params or {}))
    
    # Run analysis
    return await your_analysis_function(
        data_id, data_store, analysis_params, context
    )
```

#### 5. Add Tests

Create comprehensive tests in `tests/`:

```python
# tests/test_your_tool.py
import pytest
from chatspatial.tools.your_tool import your_analysis_function
from chatspatial.models.data import YourAnalysisParameters

@pytest.fixture
def sample_data():
    # Create or load sample data for testing
    pass

async def test_your_analysis_basic(sample_data):
    """Test basic functionality."""
    params = YourAnalysisParameters(method="basic")
    result = await your_analysis_function("test_data", sample_data, params, None)
    
    assert result.success
    assert result.n_features_analyzed > 0
    assert result.results_key in sample_data["test_data"]["adata"].uns

async def test_your_analysis_error_handling(sample_data):
    """Test error handling."""
    # Test with invalid parameters
    params = YourAnalysisParameters(parameter1=-1)  # Invalid value
    
    with pytest.raises(ProcessingError):
        await your_analysis_function("test_data", sample_data, params, None)
```

### Adding New Visualization Types

To add support for visualizing your analysis results:

#### 1. Add Plot Type to Visualization Tool (`chatspatial/tools/visualization.py`)

```python
def create_your_analysis_plot(adata, params):
    """Create visualization for your analysis results."""
    
    # Get your results from adata.uns
    results = adata.uns.get(params.results_key, None)
    if results is None:
        raise ValueError(f"Results not found: {params.results_key}")
    
    # Create your plot
    fig, ax = plt.subplots(figsize=params.figure_size)
    
    # Your plotting logic here
    
    return fig
```

#### 2. Register Plot Type in the Main Visualization Function

Add your plot type to the main `visualize_data` function's plot type mapping.

### Spatial Analysis Methods

For new spatial analysis methods (like adding support for a new spatial variable genes method):

1. **Study Existing Patterns**: Look at `spatial_genes.py` for GASTON, SpatialDE, SPARK implementations
2. **Follow Error Handling**: Use both data validation and MCP error handling layers
3. **Parameter Validation**: Create comprehensive Pydantic models
4. **Result Standardization**: Follow existing result model patterns
5. **Testing**: Add unit tests, integration tests, and real data tests

### Method Integration Checklist

When adding a new method:

- [ ] Parameter model created with proper validation
- [ ] Result model follows existing patterns
- [ ] Implementation handles all error cases
- [ ] MCP tool registered with proper decorators
- [ ] Comprehensive tests added
- [ ] Documentation updated
- [ ] Example usage provided
- [ ] Integration with visualization tool (if applicable)
- [ ] Memory usage optimized for large datasets
- [ ] Optional dependencies properly handled

### Dependencies
- Add new dependencies to `pyproject.toml`
- Use optional dependencies for specialized features
- Document installation requirements

## Documentation

### API Documentation
- Use Google-style docstrings
- Include parameter types and descriptions
- Provide usage examples
- Document return values

### User Documentation
- Update README.md for new features
- Add examples to `examples/` directory
- Create tutorials for complex workflows

## Issue Reporting

### Bug Reports
- Use the bug report template
- Include minimal reproducible example
- Provide environment information
- Include error messages and stack traces

### Feature Requests
- Use the feature request template
- Describe the use case clearly
- Suggest implementation approach
- Consider backwards compatibility

## Review Process

### Pull Request Requirements
- All tests must pass
- Code coverage should not decrease
- Documentation must be updated
- Follow the PR template
- Get approval from maintainers

### Review Criteria
- Code quality and style
- Test coverage and quality
- Documentation completeness
- Backwards compatibility
- Performance impact

## Community Guidelines

### Code of Conduct
- Be respectful and inclusive
- Focus on constructive feedback
- Help newcomers learn
- Acknowledge contributions
- Provide clear, reproducible examples when reporting issues
- Respect the time constraints of maintainers

### Communication Channels

1. **GitHub Issues**: Bug reports, feature requests, method proposals
2. **GitHub Discussions**: General questions, usage help, method comparisons
3. **Pull Request Reviews**: Code-specific feedback and suggestions
4. **Documentation**: Update guides and examples to help future contributors

### When Contributing

#### **Before Starting Work**

1. **Check existing issues** to avoid duplication
2. **Discuss major changes** in GitHub issues first
3. **Review the USER_GUIDE.md** to understand current capabilities
4. **Check method compatibility** with existing analysis workflows

#### **Quality Standards**

- **Code must be production-ready**: All methods should be thoroughly tested
- **Follow existing patterns**: Study current implementations before adding new ones
- **Document everything**: Include comprehensive docstrings and usage examples
- **Test comprehensively**: Unit tests, integration tests, and real-data validation

### Common Pitfalls and How to Avoid Them

#### **1. Image Display Bugs**
**Problem**: Modifying image handling breaks visualization in AI assistants
**Solution**: 
- Never modify image object handling in `tool_error_handling.py`
- Always test image display end-to-end
- Read `/docs/CRITICAL_IMAGE_DISPLAY_BUG.md` before touching visualization code

#### **2. Memory Issues with Large Datasets**
**Problem**: New methods consume too much memory
**Solution**:
- Test with datasets >3000 spots during development
- Implement subsampling options for large datasets  
- Use memory profiling during testing
- Provide clear memory requirements in documentation

#### **3. Optional Dependency Hell**
**Problem**: New methods introduce complex dependency conflicts
**Solution**:
- Use optional dependencies in `pyproject.toml`
- Provide clear installation instructions
- Test installation in clean environments
- Handle import errors gracefully with informative messages

#### **4. Breaking Existing Workflows**
**Problem**: Changes break existing analysis patterns
**Solution**:
- Maintain backward compatibility for parameter models
- Test integration with existing workflows
- Follow the two-step analysis pattern (analyze → visualize)
- Update examples and documentation when changing APIs

#### **5. Inadequate Error Handling**
**Problem**: New tools don't follow the two-layer error handling pattern
**Solution**:
- Use `validate_adata()` for data validation (business logic layer)
- Use `@mcp_tool_error_handler()` for MCP tools (protocol layer)
- Provide informative error messages
- Test error conditions comprehensively

### Method Integration Best Practices

#### **When Adding New Analysis Methods**

1. **Research Integration**: Study how similar tools (GASTON, LIANA+, SpaGCN) are integrated
2. **Parameter Design**: Create comprehensive Pydantic models with validation
3. **Result Standardization**: Follow existing result model patterns
4. **Visualization Support**: Ensure results can be visualized using existing plot types
5. **Documentation**: Provide usage examples and method comparison guidance

#### **Performance Considerations**

- **Fast methods** (<2 min): Can be used freely in interactive workflows
- **Moderate methods** (2-10 min): Should provide progress feedback
- **Slow methods** (>10 min): Must support parameter optimization for speed
- **Memory-intensive methods**: Must handle large datasets gracefully

#### **Testing New Methods**

1. **Synthetic data testing**: Create controlled test cases
2. **Real data validation**: Test with actual spatial transcriptomics datasets
3. **Comparison testing**: Compare results with reference implementations
4. **Performance testing**: Measure memory usage and execution time
5. **Integration testing**: Test full workflows including visualization

### Contributing Workflow Checklist

Before submitting a pull request:

- [ ] **Code follows existing patterns** and uses established error handling
- [ ] **All tests pass** including unit, integration, and visualization tests
- [ ] **Documentation is complete** with usage examples and parameter descriptions
- [ ] **Method comparison guide updated** if adding new analysis capability
- [ ] **Memory usage tested** with realistic dataset sizes
- [ ] **Optional dependencies handled** correctly with clear installation instructions
- [ ] **Backward compatibility maintained** for existing parameter models
- [ ] **Image handling tested** if visualization is involved
- [ ] **Performance benchmarked** and optimization options provided
- [ ] **Real data validation** completed for new analysis methods

### Common Review Comments and How to Address Them

#### **"This breaks the image display"**
- **Cause**: Modified Image object handling incorrectly
- **Fix**: Revert to returning raw Image objects, never wrap them

#### **"Memory usage is too high"**
- **Cause**: New method doesn't handle large datasets efficiently
- **Fix**: Add subsampling options, implement memory-efficient algorithms

#### **"Error handling doesn't follow project patterns"**
- **Cause**: Not using the two-layer error handling system
- **Fix**: Use `validate_adata()` and `@mcp_tool_error_handler()` appropriately

#### **"Tests are insufficient"**
- **Cause**: Missing critical test scenarios
- **Fix**: Add unit tests, integration tests, error condition tests, and performance tests

#### **"Documentation is incomplete"**
- **Cause**: Missing usage examples or parameter descriptions
- **Fix**: Add comprehensive docstrings, usage examples, and update user guides

## Release Process

### Version Numbering
- Follow semantic versioning (SemVer)
- Major: Breaking changes
- Minor: New features, backwards compatible
- Patch: Bug fixes, backwards compatible

### Release Checklist
- Update version in `pyproject.toml`
- Update CHANGELOG.md
- Create release notes
- Tag release in Git
- Publish to PyPI (when ready)

## Getting Help

- Check existing issues and documentation
- Ask questions in GitHub discussions
- Contact maintainers for complex issues

Thank you for contributing to ChatSpatial!
