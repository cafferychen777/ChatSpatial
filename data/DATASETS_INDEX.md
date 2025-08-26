# ChatSpatial 数据集索引

**更新时间**: 2025-08-24 (大规模数据整合完成)  
**总数据集**: 90+ (包含大规模数据集)  
**总大小**: 约 27GB (含大规模数据)
**大规模细胞总数**: **10,224,779 细胞** (1020万+)

## 🚀 **大规模数据集 (NEW!)** - 推荐用于serious development

**详细报告**: [BIG_DATASETS_FINAL_REPORT.md](BIG_DATASETS_FINAL_REPORT.md)

### 🏆 超大规模数据 (>1M 细胞)
- **MERFISH Allen Institute**: 4,334,174 细胞 × 2 (8.7M 细胞总计, 14.2GB)
- **synthetic_1M_cells**: 1,000,000 细胞 (现代空间技术模拟, 0.62GB) ⭐⭐⭐
- **synthetic_500K_cells**: 500,000 细胞 × 2 (两个版本, 1.5GB总计) ⭐⭐

### 🎯 大规模真实数据 (10K-100K 细胞)  
- **merfish_hypothalamus.h5ad**: 73,655 细胞 (真实MERFISH + 空间坐标) ⭐⭐
- **MOP_sn_tutorial.h5ad**: 26,431 细胞 (全转录组, 1.9GB)
- **seqfish_embryo.h5ad**: 19,416 细胞 (发育时间序列)
- **多个Visium数据集**: 9,425 细胞总计 (5个数据集)

## 🌟 标准测试数据集 (快速开发用)

**位置**: `real_datasets/`  
**特点**: 100%可靠，经过验证，遵循Linus的"简单实用"原则

### 🎯 快速测试 (< 5K cells)
- **synthetic_small_spatial.h5ad**: 1,000 cells (合成空间, 2.1MB) ⭐
- **pbmc68k_reduced.h5ad**: 700 cells (单细胞, 4.3MB) ⭐  
- **pbmc3k_reference.h5ad**: 2,700 cells (单细胞, 20.5MB)
- **pbmc3k_processed.h5ad**: 2,638 cells (单细胞, 38.4MB)
- **visium_hne_adata.h5ad**: 2,688 cells (Visium, 314.1MB)

### ⚖️ 标准测试 (5K-50K cells)
- **synthetic_medium_spatial.h5ad**: 5,000 cells (合成空间, 10.1MB) ⭐
- **seqfish.h5ad**: 19,416 cells (seqFISH, 30.7MB) ⭐
- **paul15_reference.h5ad**: 2,730 cells (轨迹推断, 36.2MB)

### 🚀 性能测试 (> 50K cells)  
- **synthetic_large_spatial.h5ad**: 20,000 cells (合成空间, 154.8MB) ⭐
- **slideseqv2.h5ad**: 41,786 cells (Slide-seq, 251.3MB) ⭐
- **merfish.h5ad**: 73,655 cells (MERFISH, 49.2MB) ⭐

## 📋 快速选择指南

### 按技术平台
- **Visium**: `visium_hne_adata.h5ad`
- **MERFISH**: `merfish.h5ad` 
- **seqFISH**: `seqfish.h5ad`
- **Slide-seq**: `slideseqv2.h5ad`
- **单细胞**: `pbmc3k_processed.h5ad`, `paul15_reference.h5ad`
- **合成数据**: `synthetic_small/medium/large_spatial.h5ad`

### 按使用场景
- **工具调试**: `synthetic_small_spatial.h5ad` (最快)
- **功能验证**: `synthetic_medium_spatial.h5ad` 
- **方法测试**: 选择对应技术的真实数据
- **性能测试**: `synthetic_large_spatial.h5ad`

## 传统数据集 (兼容性保留)

### 🔬 空间数据集 (spatial_datasets/)
- **squidpy_merfish.h5ad**: 73,655 cells (MERFISH)
- **squidpy_slideseqv2.h5ad**: 41,786 cells (Slide-seq)
- **squidpy_seqfish.h5ad**: 19,416 cells (seqFISH)
- **slideseq_cerebellum.h5ad**: 15,000 cells
- **osmfish_somatosensory.h5ad**: 中等规模
- **squidpy_visium.h5ad**: Visium参考

### 🎯 演示数据集 (demo_datasets/)
- **visium_demo.h5ad**: 1,000 cells
- **seqfish_demo.h5ad**: 1,000 cells
- **jurkat_293t_mixture_simulated.h5ad**: 混合细胞类型

## 详细分类

- **[Real Datasets](real_datasets/README.md)**: 44 datasets
- **[Demo Datasets](demo_datasets/README.md)**: 6 datasets
- **[Reference Datasets](reference_datasets/README.md)**: 3 datasets
- **[Synthetic Datasets](synthetic_datasets/README.md)**: 13 datasets
- **[Benchmark Datasets](benchmark_datasets/README.md)**: 5 datasets
