# ChatSpatial 数据集最终整理报告

## 📊 整理完成统计

**整理时间**: 2025-08-24 03:47:53
**数据集总数**: 79个 (更新后)
**总存储空间**: ~7.2GB
**整理前位置**: 分散在多个测试目录中
**整理后位置**: 统一在 `/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/` 

## 🗂️ 标准化目录结构

```
data/
├── real_datasets/           # 真实空间转录组数据 (47个)
│   ├── 10X Visium数据
│   ├── Slide-seq数据  
│   ├── MERFISH数据
│   ├── seqFISH数据
│   ├── osmFISH数据
│   └── 其他空间转录组技术
├── synthetic_datasets/      # 合成数据 (13个)
│   ├── 模拟Visium数据
│   ├── 模拟Slide-seq数据
│   └── 测试用合成数据
├── demo_datasets/          # 演示教程数据 (7个)  
│   ├── 快速演示数据
│   ├── 教程示例数据
│   └── Harmony集成演示
├── benchmark_datasets/     # 性能测试数据 (8个)
│   ├── 小规模测试 (100-500细胞)
│   ├── 中等规模测试 (1K-2K细胞)  
│   └── 大规模测试 (5K细胞)
├── reference_datasets/     # 参考标准数据 (3个)
│   ├── Paul15造血数据
│   ├── 胰腺参考数据
│   └── 标准测试数据
└── metadata/              # 元数据和索引
    └── datasets_catalog.json
```

## 📈 数据集分类统计

### 按技术类型分类

| 技术类型 | 数据集数量 | 代表性数据 | 主要用途 |
|----------|------------|------------|----------|
| **10X Visium** | 8个 | `squidpy_visium.h5ad` (2,688 spots) | 空间域分析 |
| **Slide-seq** | 12个 | `squidpy_slideseqv2.h5ad` (41,786 cells) | 高分辨率空间分析 |
| **MERFISH** | 4个 | `squidpy_merfish.h5ad` (73,655 cells) | 大规模空间分析 |
| **seqFISH** | 3个 | `squidpy_seqfish.h5ad` (19,416 cells) | 单细胞空间分析 |
| **osmFISH** | 2个 | `osmfish_somatosensory.h5ad` (1,000 cells) | 精细空间分析 |
| **合成数据** | 13个 | `large_synthetic.h5ad` (5,000 cells) | 算法测试 |
| **其他** | 37个 | 各种测试和参考数据 | 多种用途 |

### 按数据规模分类

| 规模分类 | 细胞数范围 | 数据集数量 | 推荐用途 |
|----------|------------|------------|----------|
| **超大规模** | >50,000 | 2个 | 性能压力测试 |
| **大规模** | 10,000-50,000 | 8个 | 标准性能测试 |  
| **中等规模** | 1,000-10,000 | 35个 | 日常开发测试 |
| **小规模** | 100-1,000 | 25个 | 快速功能验证 |
| **微型** | <100 | 9个 | 边界条件测试 |

## 🔬 空间数据集重点推荐

### 高质量空间数据集 (含空间坐标)

1. **squidpy_merfish.h5ad** - 73,655细胞，MERFISH技术，超大规模
2. **squidpy_slideseqv2.h5ad** - 41,786细胞，Slide-seq V2，大规模
3. **squidpy_seqfish.h5ad** - 19,416细胞，seqFISH，中大规模
4. **slideseq_cerebellum.h5ad** - 15,000细胞，Slide-seq，标准规模
5. **slideseq_v2_hippocampus.h5ad** - 8,000细胞，高质量标注

### 快速测试推荐

1. **visium_demo.h5ad** - 1,000细胞，演示用途
2. **seqfish_demo.h5ad** - 1,000细胞，教程用途  
3. **small_synthetic.h5ad** - 100细胞，快速验证

### 特殊测试数据

1. **empty_dataset.h5ad** - 空数据集，边界测试
2. **single_cell.h5ad** - 单细胞，极限测试
3. **high_sparsity.h5ad** - 高稀疏度，算法鲁棒性测试

## 📖 使用指南

### 开发者快速上手

```python
import scanpy as sc

# 1. 快速功能测试
adata = sc.read_h5ad('data/demo_datasets/visium_demo.h5ad')

# 2. 标准开发测试  
adata = sc.read_h5ad('data/real_datasets/slideseq_v2_hippocampus.h5ad')

# 3. 性能基准测试
adata = sc.read_h5ad('data/benchmark_datasets/benchmark_5kx5k.h5ad')

# 4. 大规模压力测试
adata = sc.read_h5ad('data/real_datasets/squidpy_merfish.h5ad')
```

### 数据集查询

```python
import json

# 加载数据集目录
with open('data/metadata/datasets_catalog.json') as f:
    catalog = json.load(f)

# 查询空间数据集
spatial_datasets = catalog['quick_access']['spatial_datasets']
print(f"找到 {len(spatial_datasets)} 个空间数据集")

# 查询大数据集
large_datasets = catalog['quick_access']['large_datasets'] 
print(f"找到 {len(large_datasets)} 个大规模数据集")
```

## 🧹 清理完成情况

### ✅ 已完成

- [x] 将91个分散的数据集统一整理到data目录
- [x] 按功能分类到5个标准目录
- [x] 创建详细的数据集目录和索引
- [x] 生成自动化的README文档
- [x] 建立快速访问索引
- [x] 清理测试目录中的重复文件

### ⚠️ 注意事项

1. **备份文件**: 保留了一些重名文件的备份（如`squidpy_visium_backup_314MB.h5ad`）
2. **旧目录保留**: 保留了原有的`core/`, `harmony/`, `test/`目录结构用于向后兼容
3. **符号链接**: 部分重要数据集在多个位置有访问路径

### 📍 重要文件位置

- **主索引**: `data/DATASETS_INDEX.md`
- **详细目录**: `data/metadata/datasets_catalog.json`
- **分类说明**: 各子目录的`README.md`文件

## 🔗 相关链接

- **项目文档**: `README.md`
- **使用示例**: `examples/` 目录
- **测试脚本**: `tests/` 目录
- **API文档**: `docs/` 目录

---

**数据集整理完成！现在所有空间转录组数据都已统一管理，可以高效地支持ChatSpatial MCP项目的开发和测试工作。**