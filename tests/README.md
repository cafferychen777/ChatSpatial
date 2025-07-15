# ChatSpatial 测试套件

本目录包含ChatSpatial的所有测试脚本，按功能分类组织。

## 目录结构

```
tests/
├── README.md                          # 本文档
├── visualization_tests/               # 可视化功能测试
│   ├── test_visualization_comprehensive.py
│   ├── test_visualization_scenarios.py
│   ├── test_direct_visualization.py
│   ├── test_final_visualization.py
│   ├── test_mcp_server_visualization.py
│   └── test_image_extraction.py
├── claude_tests/                      # Claude前端交互测试
│   ├── test_claude_conversation_simulation.py    # 原始对话模拟
│   ├── test_claude_conversation_improved.py      # 改进版（有一些失败）
│   ├── test_claude_conversation_fixed.py         # 修复版（100%成功）
│   ├── test_claude_frontend_complete.py
│   ├── test_claude_frontend_simulation.py
│   └── test_claude_functionality_simple.py
├── stress_tests/                      # 压力测试和真实场景测试
│   ├── test_comprehensive_stress.py
│   ├── test_real_analysis_session.py
│   └── test_real_world_scenarios.py
├── test_all_features_claude.py        # 全功能测试
└── test_fixes_verification.py         # 修复验证测试
```

## 运行测试

### 1. 运行所有功能测试（推荐）
```bash
python tests/test_all_features_claude.py
```

### 2. 运行Claude对话模拟测试（最完整）
```bash
python tests/claude_tests/test_claude_conversation_fixed.py
```

### 3. 运行可视化测试
```bash
python tests/visualization_tests/test_visualization_scenarios.py
```

### 4. 运行压力测试
```bash
python tests/stress_tests/test_comprehensive_stress.py
```

## 测试结果

最新的测试结果：
- **成功率**: 100% (所有测试通过)
- **测试报告位置**: `docs/test_reports/FINAL_TEST_REPORT.md`

## 依赖要求

运行测试前，请确保安装了所有必需的依赖：

```bash
# 基础依赖
pip install -r requirements.txt

# 完整依赖（包括高级功能）
pip install -r requirements-full.txt
```

## 注意事项

1. 某些测试会生成临时数据和可视化文件
2. 测试运行时间约2-5分钟（取决于测试范围）
3. 如果遇到依赖问题，请参考`requirements-full.txt`