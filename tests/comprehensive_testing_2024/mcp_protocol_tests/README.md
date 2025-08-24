# ChatSpatial MCP Protocol Tests

这是ChatSpatial MCP服务器协议层的全面测试套件。按照Linus Torvalds的"好品味"原则设计 - 简单、直接、可靠。

## 测试架构

### 核心原则
1. **数据结构优先** - 测试MCP协议的核心数据结构（工具、参数、响应）
2. **消除特殊情况** - 所有工具使用统一的测试模式  
3. **简洁实现** - 直接测试MCP消息收发，避免过度抽象
4. **零破坏性** - 测试不修改生产代码

## 测试模块

### 1. 服务器启动测试 (`test_server_startup.py`)
- **目的**: 验证MCP服务器在不同传输模式下的启动能力
- **覆盖**: stdio传输、HTTP传输、启动性能、内存使用
- **性能要求**: 
  - 平均启动时间 < 2秒
  - 最大启动时间 < 5秒  
  - 内存增长 < 100MB

### 2. 工具注册测试 (`test_tool_registration.py`)
- **目的**: 验证15个MCP工具的正确注册和元数据格式
- **覆盖**: 工具发现、元数据完整性、JSON Schema兼容性
- **验证工具**:
  ```
  load_data, preprocess_data, visualize_data, annotate_cells,
  analyze_spatial_data, find_markers, analyze_velocity_data,
  analyze_trajectory_data, integrate_samples, deconvolve_data,
  identify_spatial_domains, analyze_cell_communication,
  analyze_enrichment, find_spatial_genes, register_spatial_data,
  calculate_spatial_statistics
  ```

### 3. 参数验证测试 (`test_parameter_validation.py`)
- **目的**: 测试参数验证、类型转换、边界检查
- **覆盖**: Pydantic模型验证、默认值处理、错误消息质量
- **性能要求**: 参数验证 < 0.001秒/次

### 4. 错误响应测试 (`test_error_responses.py`)
- **目的**: 验证MCP协议错误处理和JSON-RPC 2.0兼容性
- **覆盖**: 标准错误代码、自定义错误、错误恢复、并发错误处理
- **错误代码范围**:
  - 标准JSON-RPC: -32603至-32600
  - ChatSpatial自定义: -32099至-32000

### 5. HTTP传输测试 (`test_http_transport.py`)
- **目的**: 验证HTTP传输层、REST API、并发处理
- **覆盖**: 基本端点、会话管理、CORS配置、并发处理
- **性能要求**: 
  - 平均请求时间 < 0.1秒
  - 吞吐量 > 10请求/秒
  - 并发成功率 > 80%

## 快速开始

### 运行所有测试
```bash
cd /Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/mcp_protocol_tests/

# 基本运行
python run_protocol_tests.py

# 详细输出
python run_protocol_tests.py --verbose

# 包含性能测试
python run_protocol_tests.py --performance --verbose
```

### 运行特定测试模块
```bash
# 只测试服务器启动
python run_protocol_tests.py --modules test_server_startup

# 测试工具注册和参数验证
python run_protocol_tests.py --modules test_tool_registration test_parameter_validation
```

### 使用pytest直接运行
```bash
# 运行单个测试文件
pytest test_server_startup.py -v

# 运行所有测试
pytest -v

# 生成覆盖率报告
pytest --cov=chatspatial --cov-report=html
```

## 测试报告

测试运行后会生成以下报告：

- **JSON报告**: `reports/protocol_test_report_YYYYMMDD_HHMMSS.json`
- **HTML报告**: `reports/protocol_test_report_YYYYMMDD_HHMMSS.html` 
- **性能基准**: `reports/performance_benchmark_YYYYMMDD_HHMMSS.json`
- **最新报告**: `reports/latest_report.html` (符号链接)

### 报告内容
- 每个测试模块的详细结果
- 性能指标和基准
- 错误消息和堆栈跟踪
- 系统信息和环境详情

## 性能基准

### 当前基准指标
| 测试类别 | 指标 | 要求 |
|---------|------|------|
| 服务器启动 | 平均时间 | < 2秒 |
| 工具发现 | 平均时间 | < 0.1秒 |
| 参数验证 | 平均时间 | < 0.001秒 |
| HTTP请求 | 平均响应 | < 0.1秒 |
| 并发处理 | 成功率 | > 80% |

### 性能测试标记
使用`@pytest.mark.performance`标记的测试：
```bash
# 只运行性能测试
pytest -m performance

# 跳过性能测试  
pytest -m "not performance"
```

## 依赖要求

### 核心依赖
- `pytest >= 7.0.0`
- `httpx >= 0.24.0`
- `fastapi >= 0.68.0`
- `pydantic >= 1.8.0`

### 可选依赖
- `pytest-cov` - 覆盖率报告
- `pytest-json-report` - JSON格式报告
- `pytest-html` - HTML格式报告

### 安装依赖
```bash
pip install pytest pytest-cov pytest-json-report pytest-html httpx
```

## 故障排除

### 常见问题

1. **模块导入失败**
   ```bash
   # 确保在正确的目录
   cd /Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/mcp_protocol_tests/
   
   # 检查Python路径
   python -c "import sys; print('\\n'.join(sys.path))"
   ```

2. **HTTP测试失败**
   ```bash
   # 检查端口是否被占用
   netstat -an | grep 28765
   
   # 或使用不同端口
   export TEST_PORT=38765
   ```

3. **性能测试不稳定**
   - 在资源受限环境中调整性能阈值
   - 使用`--performance`标志包含性能测试

### 调试模式
```bash
# 启用详细日志
python run_protocol_tests.py --verbose

# 使用pytest调试
pytest --pdb test_server_startup.py::TestServerStartup::test_stdio_server_creation

# 查看详细错误
pytest -vvv --tb=long
```

## 贡献指南

### 添加新测试
1. 遵循现有的命名约定：`test_<category>.py`
2. 使用类组织相关测试：`TestCategoryName`
3. 添加适当的文档字符串和类型提示
4. 更新`__init__.py`中的导出列表

### 测试质量标准
- 每个测试应该有清晰的目的和断言
- 使用有意义的测试数据和边界条件  
- 包含性能要求和质量指标
- 提供有用的错误消息

### 性能测试规则
- 标记为`@pytest.mark.performance`
- 设置合理的性能阈值
- 包含系统信息上下文
- 测量并报告关键指标

## License

This test suite is part of ChatSpatial and follows the same license terms.