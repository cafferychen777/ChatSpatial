# ChatSpatial Harmony Integration Demo

快速演示ChatSpatial的Harmony批次效应校正功能。

**数据目录位置**: `data/harmony/` (已从 `harmony_datasets/` 重构)

## 🚀 Quick Start (2分钟演示)

```bash
# 1. 创建快速演示数据集
python create_quick_demo.py

# 2. 运行快速演示  
python quick_integration.py
```

## 📋 Claude Desktop使用方法

### 基础集成
```python
from chatspatial.tools.integration import integrate_multiple_samples

# 加载数据集
adata1 = sc.read_h5ad("dataset1.h5ad")
adata2 = sc.read_h5ad("dataset2.h5ad")

# 运行Harmony集成
combined = integrate_multiple_samples(
    [adata1, adata2], 
    batch_key='batch',
    method='harmony'
)
```

### 推荐测试数据集

| 数据集 | 细胞数 | 预计时间 | 适用场景 |
|--------|--------|----------|----------|
| `quick_demo_*.h5ad` | 1000 | 2分钟 | 快速测试 |
| `pure_*.h5ad` | 3000 | 5分钟 | 标准测试 |  
| `jurkat_293t_mixture.h5ad` | 4600 | 10分钟 | 完整测试 |

## 🔧 数据集说明

- **quick_demo_293t.h5ad**: 500个293T细胞，800个基因
- **quick_demo_jurkat.h5ad**: 500个Jurkat细胞，800个基因
- **quick_demo_combined.h5ad**: 预合并的1000细胞数据集

## ⚡ 性能优化建议

1. **小规模测试**: 先用`quick_demo_*.h5ad`验证功能
2. **参数调整**: 减少`n_pcs`参数可提升速度
3. **数据预处理**: 预先过滤低质量细胞和基因

## 🎯 Claude Desktop示例对话

```text
请帮我测试ChatSpatial的Harmony集成：

1. 使用quick_demo_293t.h5ad和quick_demo_jurkat.h5ad
2. 运行Harmony批次校正
3. 生成UMAP可视化对比集成前后效果
4. 报告集成质量指标
```

## 📁 文件说明

- `create_quick_demo.py` - 创建演示数据集
- `quick_integration.py` - 快速演示脚本 (非生产用)
- `test_chatspatial_integration.py` - 完整功能测试
- `quick_validation.py` - 安全修复验证

## ✅ 验证清单

- [ ] 数据加载正常
- [ ] Harmony集成完成
- [ ] 生成X_harmony embedding
- [ ] 生成X_umap可视化  
- [ ] 批次混合效果良好
- [ ] 细胞类型分离清晰