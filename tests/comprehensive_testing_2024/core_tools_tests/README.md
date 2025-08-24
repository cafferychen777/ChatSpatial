# ChatSpatial 高级分析工具测试套件

本测试套件按照 **Linus Torvalds 的软件质量哲学**设计，专门验证 ChatSpatial 的四个高级分析工具模块。

## 测试目标

### 核心工具模块
1. **trajectory.py** - 轨迹分析 (CellRank, Palantir, DPT)
2. **pathway_enrichment.py** - 通路富集分析 (GSEA, ORA, ssGSEA)
3. **spatial_enrichment.py** - 空间富集分析 (EnrichMap集成)
4. **spatial_registration.py** - 空间配准 (PASTE算法)

### 测试数据集
- **轨迹分析**: `data/spatial_datasets/squidpy_slideseqv2.h5ad` (42K细胞)
- **富集分析**: `data/spatial_datasets/squidpy_merfish.h5ad` (73K细胞)
- **快速测试**: `data/demo_datasets/visium_demo.h5ad` (1K细胞)
- **配准测试**: `data/benchmark_datasets/benchmark_5kx5k.h5ad` (5K×5K)

## 快速开始

### 基本用法

```bash
# 快速测试（推荐用于日常验证）
python run_advanced_analysis_tests.py --quick

# 完整测试（包括性能基准）
python run_advanced_analysis_tests.py --full

# 仅检查依赖状态
python run_advanced_analysis_tests.py --report-only
```

## 依赖要求

### 核心依赖 (Critical)
```bash
pip install scanpy pandas numpy scipy matplotlib seaborn scikit-learn
```

### 轨迹分析
```bash
pip install scvelo>=0.2.5 cellrank>=2.0.0 palantir-sc>=1.0.0
```

### 通路富集分析
```bash
pip install gseapy>=1.0.0
```

### 空间富集分析
```bash
pip install pygam>=0.8.0 adjustText>=0.7.0
```

### 空间配准
```bash
pip install POT>=0.8.0 paste-bio>=1.0.0
```

## 输出文件

- **`advanced_analysis_test_results.json`** - 完整测试结果
- **`advanced_analysis_test_summary.md`** - 测试摘要
- **`dependency_report.json`** - 依赖状态详细报告
- **`dependency_report.md`** - 依赖状态摘要

## 测试架构

遵循 Linus 的 "Good Taste" 原则：
- **数据结构简洁**: 统一的 data_store 字典
- **消除特殊情况**: 标准化测试流程
- **Never break userspace**: 全面的fallback测试

## 设计哲学 (Linus式原则)

### 1. "好品味" (Good Taste)
- 数据结构简洁明了
- 消除特殊情况和边界条件
- 函数职责单一

### 2. "Never Break Userspace"
- 全面的Fallback测试机制
- 向后兼容验证
- 错误恢复和优雅降级

### 3. 实用主义优先
- 测试真实场景，不是假想威胁
- 性能基准测试
- 生物学结果验证

### 4. 简洁执念
- 三层测试架构：依赖检查 → 功能测试 → 性能评估
- 清晰的错误信息
- 最小复杂度设计

## 关键测试验证点

### 轨迹分析 (trajectory.py)
- ✅ CellRank集成测试
- ✅ RNA velocity计算
- ✅ 伪时间分析
- ✅ 轨迹可视化
- ✅ 分支点识别

### 通路富集分析 (pathway_enrichment.py)
- ✅ GSEA分析准确性
- ✅ 通路数据库集成
- ✅ 富集结果可视化
- ✅ 统计显著性验证
- ✅ 多重检验校正

### 空间富集分析 (spatial_enrichment.py)
- ✅ 空间模式富集
- ✅ 局部富集检测
- ✅ 空间自相关分析
- ✅ 富集热图生成

### 空间配准 (spatial_registration.py)
- ✅ PASTE算法集成
- ✅ 配准精度评估
- ✅ 多切片配准
- ✅ 配准结果可视化

## 使用示例

```bash
# 运行完整测试套件
cd /Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/core_tools_tests
python run_advanced_analysis_tests.py --full --output-dir ./test_results

# 检查测试结果
ls ./test_results/
# 输出:
# - advanced_analysis_test_results.json
# - advanced_analysis_test_summary.md
# - dependency_report.json
# - dependency_report.md
# - performance_benchmark.json
```

## 故障排除

### 常见问题

1. **依赖缺失**:
   ```bash
   python run_advanced_analysis_tests.py --report-only
   # 查看缺失的依赖并按建议安装
   ```

2. **内存不足**:
   ```bash
   python run_advanced_analysis_tests.py --quick
   # 使用快速模式，减少内存使用
   ```

3. **第三方工具问题**:
   - CellRank安装: `pip install cellrank>=2.0.0`
   - PASTE安装: `pip install paste-bio`
   - EnrichMap: 检查 `third_party/EnrichMap` 目录

### 查看详细日志
```bash
tail -f advanced_analysis_test.log
```

## 性能基准

- **快速测试**: 30-120秒
- **完整测试**: 5-20分钟
- **内存需求**: 2-8GB RAM (取决于数据集大小)

## 贡献

遵循Linus的代码质量标准：
1. Good Taste - 简洁明了的代码
2. 实用主义 - 解决真实问题
3. Never break userspace - 保持向后兼容
4. 充分测试 - 每个功能都要有测试

---

**"理论和实践有时会冲突。理论输。每一次。"** - Linus Torvalds