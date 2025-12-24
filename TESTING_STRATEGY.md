# ChatSpatial 自动化测试体系设计方案

## 📋 执行摘要

**目标**: 建立 80%+ 代码覆盖率的自动化测试体系，满足学术发表要求

**预计时间**: 2-4 周完整实现

**关键挑战**:
- 异步 MCP 工具函数（20个 @mcp.tool 装饰器）
- 大型外部依赖（scVI, Tangram, R packages）
- AnnData 对象依赖（需要真实或mock数据）
- 可视化测试（matplotlib输出）

---

## 🏗️ 代码架构分析

### 核心组件

```
chatspatial/
├── server.py              # 20个 MCP tools (异步函数)
├── models/
│   ├── data.py           # 28个 Pydantic 参数模型
│   └── analysis.py       # 结果模型
├── tools/                # 16个分析模块
│   ├── preprocessing.py  # 必须优先测试
│   ├── differential.py   # 相对简单，适合早期测试
│   ├── annotation.py     # 复杂，需要mock
│   ├── deconvolution.py  # 最复杂，需要mock
│   └── ...
├── utils/                # 11个工具模块
│   ├── data_loader.py    # 核心，必须测试
│   ├── data_adapter.py   # 核心，必须测试
│   ├── error_handling.py # 独立，容易测试
│   └── ...
└── spatial_mcp_adapter.py # MCP 适配器
```

### 依赖关系层级

```
Level 1 (无外部依赖，优先测试):
  - models/data.py         # Pydantic 模型验证
  - models/analysis.py     # 结果模型
  - utils/error_handling.py # 错误类和装饰器
  - utils/validation.py    # 验证函数

Level 2 (基础依赖，scanpy/squidpy):
  - utils/data_loader.py   # 数据加载
  - utils/data_adapter.py  # 数据标准化
  - utils/data_validator.py

Level 3 (核心工具，需要 AnnData):
  - tools/preprocessing.py # 预处理
  - tools/differential.py  # 差异表达
  - tools/visualization.py # 可视化

Level 4 (高级工具，需要 mock):
  - tools/annotation.py    # 细胞注释
  - tools/deconvolution.py # 去卷积
  - tools/spatial_statistics.py

Level 5 (集成测试):
  - server.py             # MCP 工具端到端测试
```

---

## 📁 测试目录结构

```
tests/
├── conftest.py                     # pytest 配置和全局 fixtures
├── fixtures/                       # 测试数据和工具
│   ├── __init__.py
│   ├── mock_adata.py              # 生成mock AnnData对象
│   ├── mock_models.py             # Mock外部模型（scVI, Tangram等）
│   ├── sample_data/               # 小型测试数据集 (<10MB)
│   │   ├── visium_sample.h5ad     # 100 spots × 100 genes
│   │   ├── merfish_sample.h5ad    # 小型MERFISH数据
│   │   └── reference.h5ad         # 参考数据集
│   └── expected_outputs/          # 预期输出用于回归测试
│       ├── preprocessing_output.json
│       └── differential_output.json
│
├── unit/                          # 单元测试 (目标: 70%覆盖率)
│   ├── test_models/
│   │   ├── test_data_models.py   # Pydantic参数模型验证
│   │   └── test_analysis_models.py
│   ├── test_utils/
│   │   ├── test_error_handling.py
│   │   ├── test_data_loader.py
│   │   ├── test_data_adapter.py
│   │   ├── test_data_validator.py
│   │   ├── test_validation.py
│   │   └── test_image_utils.py
│   └── test_tools/
│       ├── test_preprocessing.py
│       ├── test_differential.py
│       ├── test_visualization.py
│       ├── test_annotation_mock.py     # 使用mock
│       ├── test_deconvolution_mock.py  # 使用mock
│       ├── test_spatial_statistics.py
│       └── ...
│
├── integration/                   # 集成测试 (目标: 10个关键工作流)
│   ├── test_basic_workflow.py    # 加载 → 预处理 → 聚类
│   ├── test_annotation_workflow.py
│   ├── test_spatial_workflow.py  # 空间统计 → 空间域 → 可视化
│   ├── test_communication_workflow.py
│   └── test_end_to_end.py        # 完整分析流程
│
├── mcp/                          # MCP 协议测试
│   ├── test_mcp_tools.py         # 测试20个 @mcp.tool 函数
│   ├── test_mcp_adapter.py       # 测试 SpatialResourceManager
│   └── test_mcp_error_handling.py
│
├── performance/                   # 性能测试 (可选)
│   ├── test_memory_usage.py      # 内存泄漏检测
│   └── test_speed_benchmarks.py  # 性能回归
│
└── data/                         # 测试专用数据 (gitignored)
    └── .gitkeep
```

---

## 🧪 测试策略详解

### 1. 单元测试策略

#### 1.1 Pydantic 模型测试 (优先级: P0)

**测试文件**: `tests/unit/test_models/test_data_models.py`

**测试内容**:
- ✅ 有效参数验证通过
- ✅ 无效参数被拒绝（Field constraints）
- ✅ 默认值正确
- ✅ 模型序列化/反序列化

**示例**:
```python
import pytest
from pydantic import ValidationError
from chatspatial.models.data import PreprocessingParameters

def test_preprocessing_params_valid():
    """测试有效的预处理参数"""
    params = PreprocessingParameters(
        normalization="log",
        n_hvgs=2000,
        n_pcs=30
    )
    assert params.normalization == "log"
    assert params.n_hvgs == 2000

def test_preprocessing_params_invalid_normalization():
    """测试无效的归一化方法"""
    with pytest.raises(ValidationError):
        PreprocessingParameters(normalization="invalid_method")

def test_preprocessing_params_out_of_range():
    """测试超出范围的参数"""
    with pytest.raises(ValidationError):
        PreprocessingParameters(n_hvgs=10000)  # 超过5000上限
```

**覆盖目标**: 28个模型类 × 平均5个测试 = 140个测试用例

---

#### 1.2 工具函数测试 (优先级: P0)

**测试文件**: `tests/unit/test_utils/test_error_handling.py`

**测试内容**:
```python
import pytest
from chatspatial.utils.error_handling import (
    validate_adata,
    DataNotFoundError,
    suppress_output
)
from tests.fixtures.mock_adata import create_mock_adata

def test_validate_adata_valid():
    """测试验证有效的 AnnData 对象"""
    adata = create_mock_adata(n_obs=100, n_vars=50)
    adata.obs['leiden'] = ['0', '1'] * 50

    # 不应该抛出异常
    validate_adata(
        adata,
        required_keys={'obs': ['leiden']}
    )

def test_validate_adata_missing_key():
    """测试缺失必需键时抛出异常"""
    adata = create_mock_adata()

    with pytest.raises(DataNotFoundError, match="Missing required keys"):
        validate_adata(
            adata,
            required_keys={'obs': ['non_existent_column']}
        )

def test_suppress_output_context():
    """测试输出抑制上下文管理器"""
    import warnings

    with suppress_output():
        warnings.warn("This warning should be suppressed")
        print("This print should be suppressed")
    # 验证没有输出
```

---

#### 1.3 数据加载器测试 (优先级: P0)

**测试文件**: `tests/unit/test_utils/test_data_loader.py`

**测试内容**:
```python
import pytest
import os
from chatspatial.utils.data_loader import load_spatial_data
from tests.fixtures.sample_data import get_sample_visium_path

@pytest.mark.asyncio
async def test_load_visium_h5ad():
    """测试加载 h5ad 文件"""
    data_path = get_sample_visium_path()
    result = await load_spatial_data(data_path, data_type="auto")

    assert "adata" in result
    assert result["adata"].n_obs > 0
    assert result["adata"].n_vars > 0
    assert "spatial" in result["adata"].obsm

@pytest.mark.asyncio
async def test_load_nonexistent_file():
    """测试加载不存在的文件"""
    with pytest.raises(FileNotFoundError):
        await load_spatial_data("/nonexistent/path.h5ad")

@pytest.mark.asyncio
async def test_auto_detect_file_type():
    """测试自动检测文件类型"""
    # H5AD文件
    result = await load_spatial_data("test.h5ad", data_type="auto")
    # 验证正确识别为 h5ad
```

---

#### 1.4 预处理工具测试 (优先级: P1)

**测试文件**: `tests/unit/test_tools/test_preprocessing.py`

**测试内容**:
```python
import pytest
from chatspatial.tools.preprocessing import preprocess_data
from chatspatial.models.data import PreprocessingParameters
from tests.fixtures.mock_adata import create_mock_adata

@pytest.mark.asyncio
async def test_preprocess_log_normalization():
    """测试 log 归一化"""
    adata = create_mock_adata(n_obs=100, n_vars=200)
    data_store = {"test_id": {"adata": adata}}

    params = PreprocessingParameters(
        normalization="log",
        filter_genes_min_cells=3,
        filter_cells_min_genes=30
    )

    result = await preprocess_data("test_id", data_store, params)

    assert result.success
    assert "normalized" in result.message.lower()
    assert data_store["test_id"]["adata"].n_obs > 0

@pytest.mark.asyncio
async def test_preprocess_filters_empty_genes():
    """测试过滤空基因"""
    adata = create_mock_adata_with_zeros(n_obs=100, n_vars=200)
    # 添加一些全0的基因

    data_store = {"test_id": {"adata": adata}}
    result = await preprocess_data("test_id", data_store)

    # 验证基因数减少
    assert data_store["test_id"]["adata"].n_vars < 200

@pytest.mark.asyncio
async def test_preprocess_missing_dataset():
    """测试数据集不存在时的错误处理"""
    with pytest.raises(ValueError, match="not found"):
        await preprocess_data("nonexistent_id", {})
```

**覆盖目标**:
- 5种归一化方法
- 过滤逻辑
- QC指标计算
- 错误处理

---

### 2. 集成测试策略

#### 2.1 基础工作流测试 (优先级: P1)

**测试文件**: `tests/integration/test_basic_workflow.py`

**测试内容**:
```python
import pytest
from chatspatial.utils.data_loader import load_spatial_data
from chatspatial.tools.preprocessing import preprocess_data
from chatspatial.tools.differential import differential_expression
from chatspatial.models.data import PreprocessingParameters

@pytest.mark.asyncio
@pytest.mark.integration
async def test_complete_basic_workflow():
    """测试完整的基础分析工作流"""

    # 1. 加载数据
    data_path = "tests/fixtures/sample_data/visium_sample.h5ad"
    result = await load_spatial_data(data_path)
    adata = result["adata"]

    data_store = {"workflow_test": {"adata": adata}}

    # 2. 预处理
    prep_params = PreprocessingParameters(
        normalization="log",
        n_hvgs=500,
        n_pcs=20
    )
    prep_result = await preprocess_data(
        "workflow_test",
        data_store,
        prep_params
    )
    assert prep_result.success

    # 3. 聚类 (通过预处理自动完成)
    adata = data_store["workflow_test"]["adata"]
    assert "leiden" in adata.obs.columns

    # 4. 差异表达
    de_result = await differential_expression(
        "workflow_test",
        data_store,
        group_key="leiden",
        n_top_genes=20
    )
    assert len(de_result.top_genes) > 0

    # 验证数据流完整性
    assert adata.n_obs > 0
    assert "X_pca" in adata.obsm
```

---

### 3. MCP 协议测试

#### 3.1 MCP Tool 测试 (优先级: P2)

**测试文件**: `tests/mcp/test_mcp_tools.py`

**测试内容**:
```python
import pytest
from chatspatial import server
from mcp.server.fastmcp import Context

@pytest.mark.asyncio
async def test_mcp_load_data_tool():
    """测试 load_data MCP 工具"""
    mock_context = MockContext()

    result = await server.load_data(
        data_path="tests/fixtures/sample_data/visium_sample.h5ad",
        data_type="auto",
        context=mock_context
    )

    assert result.n_cells > 0
    assert result.n_genes > 0
    assert result.spatial_coordinates_available

@pytest.mark.asyncio
async def test_mcp_preprocess_data_tool():
    """测试 preprocess_data MCP 工具"""
    # 先加载数据
    load_result = await server.load_data(...)

    # 预处理
    prep_result = await server.preprocess_data(
        data_id=load_result.id,
        params=PreprocessingParameters()
    )

    assert prep_result.success

@pytest.mark.asyncio
async def test_mcp_error_handling():
    """测试 MCP 错误处理装饰器"""
    with pytest.raises(ValueError):
        await server.preprocess_data(
            data_id="nonexistent",
            params=PreprocessingParameters()
        )
```

**覆盖目标**: 20个 MCP 工具 × 3个测试场景 = 60个测试

---

### 4. Mock 策略（高级工具）

#### 4.1 细胞注释 Mock 测试

**测试文件**: `tests/unit/test_tools/test_annotation_mock.py`

**策略**: Mock外部模型（Tangram, scANVI等）

```python
import pytest
from unittest.mock import patch, MagicMock
from chatspatial.tools.annotation import annotate_cell_types
from chatspatial.models.data import AnnotationParameters

@pytest.mark.asyncio
@patch('chatspatial.tools.annotation.tangram')
async def test_tangram_annotation_mocked(mock_tangram):
    """使用 mock 测试 Tangram 注释"""
    # Mock Tangram 的返回值
    mock_result = MagicMock()
    mock_result.X = np.random.rand(100, 10)
    mock_tangram.mapping_utils.map_cells_to_space.return_value = mock_result

    data_store = {
        "test_id": {"adata": create_mock_adata()},
        "ref_id": {"adata": create_mock_reference()}
    }

    params = AnnotationParameters(
        method="tangram",
        reference_id="ref_id"
    )

    result = await annotate_cell_types("test_id", data_store, params)

    # 验证调用
    mock_tangram.mapping_utils.map_cells_to_space.assert_called_once()
    assert result.success

@pytest.mark.asyncio
@patch('chatspatial.tools.annotation.scvi')
async def test_scanvi_annotation_mocked(mock_scvi):
    """使用 mock 测试 scANVI 注释"""
    # Mock scVI 模型
    mock_model = MagicMock()
    mock_model.predict.return_value = np.array(['CellType1'] * 100)
    mock_scvi.model.SCANVI.load.return_value = mock_model

    # 测试逻辑...
```

---

## 🛠️ 测试工具和 Fixtures

### conftest.py 配置

```python
"""
pytest 配置和全局 fixtures
"""
import pytest
import numpy as np
import scanpy as sc
from anndata import AnnData

# ========== Pytest 配置 ==========

def pytest_configure(config):
    """配置 pytest markers"""
    config.addinivalue_line("markers", "unit: 单元测试")
    config.addinivalue_line("markers", "integration: 集成测试")
    config.addinivalue_line("markers", "slow: 慢速测试 (>5s)")
    config.addinivalue_line("markers", "requires_r: 需要 R 环境")
    config.addinivalue_line("markers", "requires_gpu: 需要 GPU")

# ========== 全局 Fixtures ==========

@pytest.fixture
def mock_adata():
    """创建基础 mock AnnData 对象"""
    from tests.fixtures.mock_adata import create_mock_adata
    return create_mock_adata(n_obs=100, n_vars=200)

@pytest.fixture
def mock_adata_with_spatial():
    """创建带空间坐标的 mock AnnData"""
    from tests.fixtures.mock_adata import create_mock_adata
    return create_mock_adata(
        n_obs=100,
        n_vars=200,
        add_spatial=True
    )

@pytest.fixture
def data_store(mock_adata):
    """创建测试用的 data_store"""
    return {"test_dataset": {"adata": mock_adata}}

@pytest.fixture
def sample_visium_path():
    """返回样本 Visium 数据路径"""
    return "tests/fixtures/sample_data/visium_sample.h5ad"

@pytest.fixture(scope="session")
def create_sample_data():
    """Session级别的样本数据创建 (只运行一次)"""
    from tests.fixtures.mock_adata import generate_sample_datasets
    generate_sample_datasets()

# ========== Mock Classes ==========

class MockContext:
    """Mock MCP Context for testing"""
    def __init__(self):
        self.messages = []

    async def info(self, message: str):
        self.messages.append(("info", message))

    async def warning(self, message: str):
        self.messages.append(("warning", message))

    async def error(self, message: str):
        self.messages.append(("error", message))

@pytest.fixture
def mock_context():
    """创建 mock MCP context"""
    return MockContext()
```

---

### fixtures/mock_adata.py

```python
"""
生成 mock AnnData 对象的工具函数
"""
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

def create_mock_adata(
    n_obs: int = 100,
    n_vars: int = 200,
    add_spatial: bool = False,
    add_clusters: bool = False,
    sparse_matrix: bool = True,
    random_seed: int = 42
) -> AnnData:
    """
    创建用于测试的 mock AnnData 对象

    Args:
        n_obs: 细胞数量
        n_vars: 基因数量
        add_spatial: 是否添加空间坐标
        add_clusters: 是否添加聚类标签
        sparse_matrix: 是否使用稀疏矩阵
        random_seed: 随机种子

    Returns:
        Mock AnnData 对象
    """
    np.random.seed(random_seed)

    # 生成表达矩阵
    if sparse_matrix:
        # 稀疏矩阵 (90% 零值)
        X = sparse.random(n_obs, n_vars, density=0.1, format='csr')
        X.data = np.abs(X.data) * 10  # 正值表达量
    else:
        X = np.random.rand(n_obs, n_vars) * 10

    # 创建 obs (细胞元数据)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])

    # 创建 var (基因元数据)
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)])

    adata = AnnData(X=X, obs=obs, var=var)

    # 添加空间坐标
    if add_spatial:
        spatial_coords = np.random.rand(n_obs, 2) * 100
        adata.obsm['spatial'] = spatial_coords

    # 添加聚类标签
    if add_clusters:
        n_clusters = min(5, n_obs // 20)
        adata.obs['leiden'] = [f"{i % n_clusters}" for i in range(n_obs)]

    return adata

def create_mock_reference_adata(
    n_obs: int = 500,
    n_vars: int = 200,
    cell_types: list = None
) -> AnnData:
    """创建用于细胞注释的参考数据集"""
    if cell_types is None:
        cell_types = ['CellTypeA', 'CellTypeB', 'CellTypeC']

    adata = create_mock_adata(n_obs=n_obs, n_vars=n_vars)

    # 添加细胞类型标签
    adata.obs['cell_type'] = np.random.choice(
        cell_types,
        size=n_obs
    )

    return adata

def generate_sample_datasets():
    """生成所有测试用的样本数据集 (运行一次)"""
    import os

    output_dir = "tests/fixtures/sample_data"
    os.makedirs(output_dir, exist_ok=True)

    # 1. Visium 样本
    visium = create_mock_adata(
        n_obs=100,
        n_vars=100,
        add_spatial=True,
        add_clusters=True
    )
    visium.write_h5ad(f"{output_dir}/visium_sample.h5ad")

    # 2. 参考数据集
    reference = create_mock_reference_adata(n_obs=500, n_vars=100)
    reference.write_h5ad(f"{output_dir}/reference.h5ad")

    print(f"✅ 生成测试数据集到 {output_dir}/")
```

---

## 🎯 测试优先级和时间线

### Week 1: 基础设施 (P0)

**目标**: 建立测试框架，达到30%覆盖率

1. **Day 1-2**: 设置测试环境
   - [ ] 创建 `tests/` 目录结构
   - [ ] 编写 `conftest.py` 和 fixtures
   - [ ] 生成样本数据集
   - [ ] 配置 pytest.ini

2. **Day 3-4**: 模型和工具测试
   - [ ] `test_data_models.py` - 28个模型类
   - [ ] `test_error_handling.py` - 错误处理
   - [ ] `test_validation.py` - 验证函数

3. **Day 5-7**: 核心工具测试
   - [ ] `test_data_loader.py`
   - [ ] `test_data_adapter.py`
   - [ ] `test_preprocessing.py` (基础场景)

**里程碑**: 30% 覆盖率，CI/CD 集成

---

### Week 2: 核心功能 (P1)

**目标**: 覆盖主要工具模块，达到60%覆盖率

1. **Day 1-2**: 分析工具测试
   - [ ] `test_differential.py`
   - [ ] `test_visualization.py` (基础图表)
   - [ ] `test_spatial_statistics.py` (部分方法)

2. **Day 3-4**: Mock 高级工具
   - [ ] `test_annotation_mock.py`
   - [ ] `test_deconvolution_mock.py`
   - [ ] `test_cell_communication_mock.py`

3. **Day 5-7**: 集成测试
   - [ ] `test_basic_workflow.py`
   - [ ] `test_spatial_workflow.py`
   - [ ] `test_annotation_workflow.py`

**里程碑**: 60% 覆盖率，主要工作流可测试

---

### Week 3: MCP 和高级功能 (P2)

**目标**: 达到80%覆盖率

1. **Day 1-3**: MCP 协议测试
   - [ ] `test_mcp_tools.py` - 20个工具
   - [ ] `test_mcp_adapter.py`
   - [ ] `test_mcp_error_handling.py`

2. **Day 4-5**: 剩余工具模块
   - [ ] `test_integration.py` (批次校正)
   - [ ] `test_trajectory.py` (mock)
   - [ ] `test_cnv_analysis.py` (mock)

3. **Day 6-7**: 完善和修复
   - [ ] 修复所有失败测试
   - [ ] 补充边界情况测试
   - [ ] 代码覆盖率分析

**里程碑**: 80% 覆盖率，所有核心功能覆盖

---

### Week 4: CI/CD 和文档 (P3)

**目标**: 完整的测试基础设施

1. **Day 1-2**: CI/CD 配置
   - [ ] 更新 GitHub Actions 运行测试
   - [ ] 配置 codecov 覆盖率报告
   - [ ] 添加测试徽章到 README

2. **Day 3-4**: 性能和压力测试
   - [ ] `test_memory_usage.py`
   - [ ] `test_speed_benchmarks.py`
   - [ ] 大数据集测试

3. **Day 5-7**: 文档和清理
   - [ ] 编写测试文档
   - [ ] 添加测试示例到 CONTRIBUTING.md
   - [ ] 代码审查和重构

**里程碑**: 完整测试体系，可发表

---

## 📊 覆盖率目标

```
目标覆盖率: ≥ 80%

分模块目标:
├── models/           95% (Pydantic 模型容易测试)
├── utils/            90% (核心工具)
├── tools/
│   ├── preprocessing.py    85%
│   ├── differential.py     85%
│   ├── visualization.py    70% (部分图表类型)
│   ├── spatial_statistics.py  75%
│   ├── annotation.py       60% (mock 主要方法)
│   ├── deconvolution.py    60% (mock 主要方法)
│   └── 其他                 65%
└── server.py         75% (MCP 工具)

总体: 80%+
```

---

## 🚀 CI/CD 集成

### 更新 `.github/workflows/test.yml`

```yaml
name: Tests

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main ]

jobs:
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
        python-version: ['3.10', '3.11']

    steps:
    - uses: actions/checkout@v3

    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: ${{ matrix.python-version }}

    - name: Cache dependencies
      uses: actions/cache@v3
      with:
        path: ~/.cache/pip
        key: ${{ runner.os }}-pip-${{ hashFiles('**/pyproject.toml') }}

    - name: Install dependencies
      run: |
        pip install -e ".[dev]"
        pip install pytest-cov pytest-asyncio

    - name: Generate test data
      run: |
        python -c "from tests.fixtures.mock_adata import generate_sample_datasets; generate_sample_datasets()"

    - name: Run unit tests
      run: |
        pytest tests/unit -v --cov=chatspatial --cov-report=xml

    - name: Run integration tests
      run: |
        pytest tests/integration -v --cov=chatspatial --cov-append --cov-report=xml

    - name: Upload coverage to Codecov
      uses: codecov/codecov-action@v3
      with:
        file: ./coverage.xml
        flags: unittests
        name: codecov-${{ matrix.os }}-py${{ matrix.python-version }}
```

---

## 🔧 pytest.ini 配置

```ini
[tool:pytest]
minversion = 6.0
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*

# Async support
asyncio_mode = auto
asyncio_default_fixture_loop_scope = function

# Markers
markers =
    unit: Unit tests
    integration: Integration tests
    slow: Slow tests (> 5 seconds)
    requires_r: Tests requiring R environment
    requires_gpu: Tests requiring GPU

# Coverage
addopts =
    -v
    --strict-markers
    --tb=short
    --cov-report=term-missing
    --cov-report=html

# Warnings
filterwarnings =
    ignore::DeprecationWarning
    ignore::PendingDeprecationWarning
    ignore:.*jinja2.*:UserWarning
    ignore:.*numpy.*version.*:UserWarning
```

---

## 📝 测试最佳实践

### 1. 命名约定
- 测试文件: `test_<module_name>.py`
- 测试类: `Test<ClassName>`
- 测试函数: `test_<what_is_tested>`
- Fixtures: `<descriptive_name>` (无 test_ 前缀)

### 2. 测试结构 (AAA Pattern)
```python
def test_function():
    # Arrange - 准备测试数据
    adata = create_mock_adata()
    params = PreprocessingParameters()

    # Act - 执行被测试的功能
    result = preprocess_data(adata, params)

    # Assert - 验证结果
    assert result.success
    assert adata.n_obs > 0
```

### 3. 参数化测试
```python
@pytest.mark.parametrize("normalization,expected", [
    ("log", "log1p"),
    ("sct", "sctransform"),
    ("none", "none"),
])
def test_normalization_methods(normalization, expected):
    params = PreprocessingParameters(normalization=normalization)
    # 测试逻辑...
```

### 4. Mock 外部依赖
```python
@patch('chatspatial.tools.annotation.scvi')
def test_with_mock(mock_scvi):
    mock_scvi.model.SCANVI.return_value = mock_model
    # 测试逻辑...
```

### 5. 异步测试
```python
@pytest.mark.asyncio
async def test_async_function():
    result = await async_function()
    assert result is not None
```

---

## 🎓 实施建议

### 优先级排序

**立即开始 (Week 1)**:
1. ✅ 创建测试目录结构
2. ✅ 编写 conftest.py 和基础 fixtures
3. ✅ 测试 Pydantic 模型（最简单）
4. ✅ 测试工具函数（error_handling, validation）

**核心功能 (Week 2)**:
5. ✅ 测试 data_loader 和 data_adapter
6. ✅ 测试 preprocessing (最重要)
7. ✅ 测试 differential expression
8. ✅ 第一个集成测试

**扩展覆盖 (Week 3-4)**:
9. ✅ Mock 高级工具（annotation, deconvolution）
10. ✅ MCP 协议测试
11. ✅ CI/CD 集成
12. ✅ 文档完善

### 快速开始命令

```bash
# 1. 创建测试目录
mkdir -p tests/{unit/{test_models,test_utils,test_tools},integration,mcp,fixtures/sample_data}

# 2. 安装测试依赖
pip install pytest pytest-asyncio pytest-cov pytest-mock

# 3. 生成样本数据
python -c "from tests.fixtures.mock_adata import generate_sample_datasets; generate_sample_datasets()"

# 4. 运行测试
pytest tests/unit -v

# 5. 检查覆盖率
pytest --cov=chatspatial --cov-report=html
open htmlcov/index.html  # 查看覆盖率报告
```

---

## 📈 成功指标

### 定量指标
- ✅ 代码覆盖率 ≥ 80%
- ✅ 所有 20 个 MCP 工具有测试
- ✅ 至少 5 个集成测试工作流
- ✅ 测试运行时间 < 5 分钟 (unit tests)
- ✅ CI/CD 自动运行测试

### 定性指标
- ✅ 测试文档完整
- ✅ 所有核心功能有测试
- ✅ 错误处理被测试
- ✅ 边界情况被覆盖
- ✅ Mock 策略清晰

---

## 🔗 参考资源

1. **pytest 文档**: https://docs.pytest.org/
2. **pytest-asyncio**: https://pytest-asyncio.readthedocs.io/
3. **pytest-cov**: https://pytest-cov.readthedocs.io/
4. **unittest.mock**: https://docs.python.org/3/library/unittest.mock.html
5. **Codecov**: https://about.codecov.io/

---

## 📞 下一步

需要我帮你：
1. ✅ 创建测试目录结构？
2. ✅ 编写第一批测试用例（models, utils）？
3. ✅ 实现 mock_adata.py 和 conftest.py？
4. ✅ 更新 GitHub Actions？

让我知道从哪里开始！
