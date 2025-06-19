# Contributing to ChatSpatial

Thank you for your interest in contributing to ChatSpatial! This document provides guidelines for contributing to the project.

## Getting Started

### Prerequisites
- Python 3.8 or higher
- Git
- Basic understanding of spatial transcriptomics
- Familiarity with Model Context Protocol (MCP) is helpful

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
   pip install -e .[all]
   pip install pytest pytest-cov black isort mypy flake8
   ```

4. **Run tests to ensure everything works**
   ```bash
   pytest tests/
   ```

## Development Workflow

### Branching Strategy
- `main`: Stable release branch
- `develop`: Development branch for new features
- Feature branches: `feature/your-feature-name`
- Bug fixes: `bugfix/issue-description`

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
- Unit tests in `tests/`
- Integration tests in `tests/integration/`
- Use pytest for all tests
- Aim for >80% code coverage

### Writing Tests
```python
import pytest
from chatspatial.tools.your_module import your_function

def test_your_function():
    """Test description."""
    # Arrange
    input_data = ...
    
    # Act
    result = your_function(input_data)
    
    # Assert
    assert result == expected_result
```

### Test Data
- Use small, synthetic datasets for tests
- Store test data in `tests/data/`
- Don't commit large data files

## Adding New Features

### Spatial Analysis Methods
1. Create module in `chatspatial/tools/`
2. Follow existing patterns for parameter validation
3. Return standardized result objects
4. Add comprehensive tests
5. Update documentation

### MCP Tools
1. Add tool definition in `chatspatial/server.py`
2. Implement tool handler function
3. Add input/output models in `chatspatial/models/`
4. Test MCP integration

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

### Communication
- Use GitHub issues for bug reports and feature requests
- Use GitHub discussions for questions and ideas
- Be clear and concise in communications

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
