# ChatSpatial MCP 测试快速启动指南

## 🚀 快速开始

### 一键运行所有测试
```bash
cd /Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024

# 运行所有可用测试模块
python run_all_tests.py

# 静默模式运行
python run_all_tests.py --quiet

# 只运行特定模块
python run_all_tests.py --modules mcp_protocol tool_functionality
```

### 分模块运行测试

#### 1. MCP协议层测试
```bash
cd mcp_protocol_tests/
python run_protocol_tests.py
```

#### 2. 核心工具功能测试
```bash
cd tool_functionality_tests/
# 快速测试
python run_comprehensive_tests.py --type quick

# 完整测试
python run_comprehensive_tests.py --type all --verbose

# 性能基准测试
python run_comprehensive_tests.py --type performance
```

#### 3. 集成工作流测试
```bash
cd integration_workflow_tests/
python run_all_workflow_tests.py
```

#### 4. 错误处理和兼容性测试
```bash
# 在根目录运行
python run_comprehensive_error_and_compatibility_tests.py
```

## 📊 测试数据集

### 测试数据位置
```
datasets/
├── small_synthetic.h5ad          # 100x500 - 快速功能测试
├── medium_synthetic.h5ad         # 1000x2000 - 标准测试
├── large_synthetic.h5ad          # 5000x3000 - 性能测试
├── benchmark_5kx5k.h5ad          # 5000x5000 - 压力测试
├── squidpy_visium.h5ad           # 真实Visium数据
├── squidpy_seqfish.h5ad          # 真实seqFISH数据
├── empty_dataset.h5ad            # 边界条件测试
├── single_cell.h5ad              # 单细胞边界测试
├── high_sparsity.h5ad            # 高稀疏度测试
├── no_spatial.h5ad               # 无空间坐标测试
└── [其他特殊测试数据]
```

### 数据集重新生成
```bash
# 重新生成所有测试数据集
python download_test_datasets.py
```

## 📝 查看测试报告

### 报告位置
测试完成后，报告会保存在 `reports/` 目录：
```
reports/
├── comprehensive_test_report_[timestamp].json    # 详细JSON报告
├── comprehensive_test_report_[timestamp].md      # 可读性报告
├── protocol_test_report.html                     # MCP协议测试HTML报告
├── tool_functionality_report.html                # 工具功能测试报告
└── [其他模块报告]
```

### 实时查看测试进度
```bash
# 运行测试时会显示实时进度
python run_all_tests.py

# 输出示例:
# 2024-12-19 10:30:15 - INFO - 开始运行: MCP协议层测试
# 2024-12-19 10:30:45 - INFO - MCP协议层测试: ✅ 成功 (30.2秒)
# 2024-12-19 10:30:45 - INFO - 开始运行: 核心工具功能测试
# ...
```

## ⚙️ 高级选项

### 自定义测试执行
```bash
# 指定输出文件名
python run_all_tests.py --output my_test_report

# 只运行必要的测试（跳过性能测试）
python run_all_tests.py --modules mcp_protocol tool_functionality error_compatibility

# 结合其他选项
python run_all_tests.py --modules mcp_protocol --output protocol_only --quiet
```

### 环境检查
```bash
# 验证测试环境
cd mcp_protocol_tests/
python validate_setup.py

cd ../tool_functionality_tests/
python -c "import pytest; print('pytest available')"
```

## 🐛 故障排除

### 常见问题

#### 1. 模块导入错误
```bash
# 确保在正确的Python环境中
which python
python -c "import chatspatial; print('✅ chatspatial可用')"
```

#### 2. 缺少依赖
```bash
# 安装测试依赖
pip install pytest pytest-asyncio pytest-html memory-profiler
```

#### 3. 测试超时
```bash
# 单独运行超时的测试模块，增加详细输出
cd tool_functionality_tests/
python run_comprehensive_tests.py --type basic --verbose
```

#### 4. 数据集问题
```bash
# 检查数据集完整性
python -c "
import scanpy as sc
import pandas as pd
try:
    adata = sc.read_h5ad('datasets/small_synthetic.h5ad')
    print(f'✅ 数据集可用: {adata.n_obs} cells, {adata.n_vars} genes')
except Exception as e:
    print(f'❌ 数据集问题: {e}')
"
```

### 性能调优

#### 减少测试时间
```bash
# 只运行快速测试
python run_all_tests.py --modules mcp_protocol error_compatibility

# 使用小数据集
cd tool_functionality_tests/
python run_comprehensive_tests.py --type quick
```

#### 详细诊断
```bash
# 启用详细日志
python run_all_tests.py --modules tool_functionality 2>&1 | tee test_log.txt

# 检查特定错误
grep "ERROR\|FAILED" test_log.txt
```

## 📈 持续集成

### GitHub Actions示例
```yaml
name: ChatSpatial Tests
on: [push, pull_request]

jobs:
  comprehensive-tests:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    - name: Setup Python
      uses: actions/setup-python@v2
      with:
        python-version: '3.9'
    - name: Install dependencies
      run: |
        pip install -r requirements-full.txt
        pip install pytest pytest-html
    - name: Run comprehensive tests
      run: |
        cd tests/comprehensive_testing_2024
        python run_all_tests.py --quiet
    - name: Upload test reports
      uses: actions/upload-artifact@v2
      with:
        name: test-reports
        path: tests/comprehensive_testing_2024/reports/
```

### 本地CI脚本
```bash
#!/bin/bash
# local_ci.sh - 本地持续集成脚本

set -e  # 遇到错误时退出

echo "🧪 ChatSpatial本地CI测试开始"

# 检查环境
python --version
pip list | grep -E "(scanpy|pandas|numpy)"

# 运行测试
cd tests/comprehensive_testing_2024
python run_all_tests.py

# 检查结果
if [ $? -eq 0 ]; then
    echo "✅ 所有测试通过"
else
    echo "❌ 测试失败"
    exit 1
fi
```

---

## 📞 获取帮助

1. **查看详细文档**: `README.md` 在各个测试目录
2. **检查实现报告**: `COMPREHENSIVE_TESTING_REPORT.md`
3. **查看测试结果**: `reports/` 目录中的生成报告
4. **调试模式**: 使用 `--verbose` 参数获取详细输出