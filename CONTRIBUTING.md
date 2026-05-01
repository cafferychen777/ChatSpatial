# Contributing to ChatSpatial

Contributions are welcome — bug reports, new analysis methods, documentation improvements, and feature requests.

## Getting Started

```bash
# Fork and clone
git clone https://github.com/YOUR_USERNAME/ChatSpatial.git
cd ChatSpatial

# Create environment and install
python3 -m venv venv && source venv/bin/activate
pip install -e ".[dev]"

# Verify
pytest tests/unit/ -x
```

**Prerequisites**: Python 3.11-3.13, Git. For R-based methods (RCTD, CellChat, SPARK-X, etc.): R 4.4+ and rpy2.

## Project Structure

```
chatspatial/
├── server.py                 # MCP tool definitions (entry point)
├── spatial_mcp_adapter.py    # ToolContext and data manager
├── config.py                 # Runtime configuration
├── tools/                    # Analysis implementations
│   ├── spatial_genes.py      # SpatialDE, SPARK-X, FlashS
│   ├── spatial_domains.py    # SpaGCN, STAGATE, GraphST, BANKSY, Leiden
│   ├── cell_communication.py # FastCCC, LIANA, CellPhoneDB, CellChat (`cellchat_r`)
│   ├── deconvolution/        # FlashDeconv, Cell2location, RCTD, etc.
│   ├── visualization/        # 11 plot types
│   └── ...
├── models/
│   ├── data.py               # Pydantic parameter models
│   └── analysis.py           # Pydantic result models
└── utils/
    ├── mcp_utils.py           # @mcp_tool_error_handler decorator
    ├── exceptions.py          # Custom exception classes
    ├── adata_utils.py         # AnnData validation helpers
    └── dependency_manager.py  # Optional dependency checking
```

## Adding a New Analysis Method

This is the most common contribution. Follow the existing pattern:

### 1. Parameter model (`models/data.py`)

```python
class YourMethodParameters(BaseModel):
    method: Literal["method_a", "method_b"] = Field(
        default="method_a",
        description="Which algorithm to use.",
    )
    n_top_genes: Optional[int] = Field(
        default=None, description="Number of top genes to return."
    )
```

### 2. Result model (`models/analysis.py`)

```python
class YourMethodResult(BaseModel):
    data_id: str
    method: str
    n_genes_analyzed: int
    results_key: Optional[str] = None
```

### 3. Tool implementation (`tools/your_tool.py`)

```python
from ..utils.exceptions import DataError, ProcessingError
from ..utils.dependency_manager import require

async def your_method(
    data_id: str,
    ctx: "ToolContext",
    params: YourMethodParameters,
) -> YourMethodResult:
    """Implement your analysis."""
    require("optional_package")  # Checks at runtime, clear error if missing

    adata = await ctx.get_adata(data_id)
    # ... analysis logic ...
    return YourMethodResult(...)
```

### 4. Register in `server.py`

```python
@mcp.tool()
@mcp_tool_error_handler()
async def your_tool(
    data_id: str,
    params: Optional[YourMethodParameters] = None,
    context: Optional[Context] = None,
) -> YourMethodResult:
    """Brief description for LLM tool selection."""
    ctx = ToolContext(_data_manager=data_manager, _mcp_context=context)
    p = _resolve_params(params, YourMethodParameters)
    return await your_method(data_id, ctx, p)
```

### 5. Add tests

```python
# tests/unit/test_your_tool.py
@pytest.mark.asyncio
async def test_your_method_basic(minimal_spatial_adata, monkeypatch):
    # Mock external dependencies, test logic
    ...
```

### Checklist

- [ ] Parameter model with Pydantic validation
- [ ] Result model following existing patterns
- [ ] Implementation using `ToolContext` (not raw data_store dict)
- [ ] Optional dependencies handled via `require()`
- [ ] MCP tool registered with `@mcp_tool_error_handler()`
- [ ] Unit tests with mocked dependencies
- [ ] Docstrings on public functions

## Code Style

```bash
# Format and lint
black chatspatial/
isort chatspatial/
ruff check chatspatial/ --fix

# Type check
mypy chatspatial/
```

- Max line length: 88 (Black default)
- Type hints on all public functions
- Imports: stdlib, third-party, local (isort handles this)

## Testing

```bash
pytest tests/unit/           # Fast, no external deps
pytest tests/integration/    # Multi-component workflows
pytest tests/e2e/            # Full MCP tool calls

# Pre-PR quality gate
make test-gates
```

- Unit tests: mock external packages, test logic in isolation
- Integration tests: test tool dispatch and result storage
- Keep test data small (<1000 spots, <500 genes)
- Set random seeds for reproducibility

## Submitting Changes

1. Create a branch: `git checkout -b feature/your-feature`
2. Make changes, run tests and linting
3. Commit with clear messages: `feat: add X method for Y analysis`
4. Open a PR against `main`

### Commit style

```
feat: add new spatial analysis method
fix: handle edge case in deconvolution
docs: update methods reference
test: add integration test for trajectory
```

## Reporting Issues

- **Bugs**: include a minimal reproducible example, error traceback, and `pip show chatspatial` output
- **Feature requests**: describe the use case and suggest which tool category it fits

## Questions?

Open a [GitHub Discussion](https://github.com/cafferychen777/ChatSpatial/discussions) or check the [docs](https://docs.cafferyang.com/).
