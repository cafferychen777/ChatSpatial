#!/usr/bin/env python3
"""
处理和验证真实空间转录组数据集的脚本
展示如何处理不同格式的Slide-seq和ST数据
"""

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io
import warnings
warnings.filterwarnings('ignore')

# 设置scanpy
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80, facecolor='white')

BASE_DIR = Path(__file__).parent / 'datasets' / 'real_datasets'

class SpatialDataProcessor:
    """空间转录组数据处理器"""
    
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        self.processed_dir = self.data_dir / 'processed'
        self.processed_dir.mkdir(exist_ok=True)
        
    def load_and_validate_dataset(self, dataset_name):
        """加载和验证数据集"""
        print(f"\n{'='*60}")
        print(f"处理数据集: {dataset_name}")
        print(f"{'='*60}")
        
        dataset_path = self.data_dir / f"{dataset_name}.h5ad"
        
        if not dataset_path.exists():
            print(f"❌ 数据集文件不存在: {dataset_path}")
            return None
        
        # 加载数据
        print("📁 加载数据集...")
        adata = sc.read_h5ad(dataset_path)
        
        # 基本信息
        print(f"📊 数据集信息:")
        print(f"   - Spots: {adata.n_obs:,}")
        print(f"   - Genes: {adata.n_vars:,}") 
        print(f"   - 数据类型: {type(adata.X)}")
        print(f"   - 稀疏性: {(adata.X == 0).sum() / (adata.X.shape[0] * adata.X.shape[1]):.3f}")
        
        # 空间信息
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            print(f"   - 空间坐标: {coords.shape}")
            print(f"   - X范围: {coords[:, 0].min():.1f} - {coords[:, 0].max():.1f}")
            print(f"   - Y范围: {coords[:, 1].min():.1f} - {coords[:, 1].max():.1f}")
        else:
            print("   - ⚠️ 缺少空间坐标信息")
        
        # 技术信息
        if 'spatial' in adata.uns:
            spatial_info = adata.uns['spatial']
            if 'technology' in spatial_info:
                tech = spatial_info['technology']
                print(f"   - 技术: {tech.get('technology', 'Unknown')}")
            if 'dataset_info' in spatial_info:
                info = spatial_info['dataset_info']
                print(f"   - 数据源: {info.get('source', 'Unknown')}")
                print(f"   - 演示数据: {info.get('is_demo', False)}")
        
        return adata
    
    def process_slideseq_format(self, count_matrix_path, coordinates_path, output_name):
        """
        处理标准Slide-seq格式数据
        
        Parameters:
        -----------
        count_matrix_path : str
            计数矩阵文件路径 (.csv, .mtx, .h5)
        coordinates_path : str  
            坐标文件路径 (.csv)
        output_name : str
            输出文件名
        """
        print(f"\n处理Slide-seq格式数据: {output_name}")
        
        # 读取计数矩阵
        if count_matrix_path.suffix == '.csv':
            print("  读取CSV格式计数矩阵...")
            expr_df = pd.read_csv(count_matrix_path, index_col=0)
            X = expr_df.values.T  # 转置：genes x cells -> cells x genes
            gene_names = expr_df.index.tolist()
            cell_names = expr_df.columns.tolist()
            
        elif count_matrix_path.suffix == '.mtx':
            print("  读取MTX格式计数矩阵...")
            X = scipy.io.mmread(count_matrix_path).T.tocsr()
            # 需要额外的基因和细胞名称文件
            gene_names = [f'Gene_{i}' for i in range(X.shape[1])]
            cell_names = [f'Cell_{i}' for i in range(X.shape[0])]
            
        elif count_matrix_path.suffix == '.h5':
            print("  读取H5格式计数矩阵...")
            import h5py
            with h5py.File(count_matrix_path, 'r') as f:
                # 具体结构取决于10X或其他格式
                pass
        
        # 读取坐标信息
        print("  读取坐标信息...")
        coords_df = pd.read_csv(coordinates_path)
        
        # 常见列名变种
        x_col = None
        y_col = None
        for col in coords_df.columns:
            if col.lower() in ['x', 'xcoord', 'x_coord', 'x_coordinate']:
                x_col = col
            elif col.lower() in ['y', 'ycoord', 'y_coord', 'y_coordinate']:
                y_col = col
        
        if x_col is None or y_col is None:
            raise ValueError("无法识别坐标列")
        
        # 匹配细胞/beads
        if len(coords_df) != len(cell_names):
            print(f"  ⚠️ 坐标数量({len(coords_df)})与细胞数量({len(cell_names)})不匹配")
            # 尝试通过索引匹配
            if 'barcode' in coords_df.columns:
                coords_df = coords_df.set_index('barcode').reindex(cell_names)
        
        # 创建AnnData对象
        print("  创建AnnData对象...")
        adata = sc.AnnData(
            X=X,
            obs=pd.DataFrame(index=cell_names),
            var=pd.DataFrame({'gene_symbol': gene_names}, index=gene_names)
        )
        
        adata.obsm['spatial'] = coords_df[[x_col, y_col]].values
        
        # 添加技术信息
        adata.uns['spatial'] = {
            'technology': {
                'technology': 'Slide-seq',
                'bead_diameter': 10,
                'resolution': 'subcellular'
            }
        }
        
        # 保存
        output_path = self.processed_dir / f"{output_name}.h5ad"
        adata.write(output_path)
        print(f"  ✅ 保存至: {output_path}")
        
        return adata
    
    def process_st_format(self, expression_file, coordinate_file, output_name):
        """
        处理标准ST格式数据
        
        Parameters:
        -----------
        expression_file : str
            表达文件路径 (.tsv, .csv)
        coordinate_file : str
            坐标文件路径 (.tsv, .csv)
        output_name : str
            输出文件名
        """
        print(f"\n处理ST格式数据: {output_name}")
        
        # 读取表达数据
        print("  读取表达矩阵...")
        if expression_file.suffix in ['.tsv', '.txt']:
            expr_df = pd.read_csv(expression_file, sep='\t', index_col=0)
        else:
            expr_df = pd.read_csv(expression_file, index_col=0)
        
        # 读取坐标数据
        print("  读取坐标信息...")
        if coordinate_file.suffix in ['.tsv', '.txt']:
            coords_df = pd.read_csv(coordinate_file, sep='\t', index_col=0)
        else:
            coords_df = pd.read_csv(coordinate_file, index_col=0)
        
        # ST数据通常是 spots x genes
        X = expr_df.values
        spot_names = expr_df.index.tolist()
        gene_names = expr_df.columns.tolist()
        
        # 匹配坐标
        coords_df = coords_df.reindex(spot_names)
        
        # 创建AnnData
        print("  创建AnnData对象...")
        adata = sc.AnnData(
            X=X,
            obs=pd.DataFrame(index=spot_names),
            var=pd.DataFrame({'gene_symbol': gene_names}, index=gene_names)
        )
        
        # 添加空间坐标
        coord_cols = coords_df.columns[:2]  # 通常前两列是x, y
        adata.obsm['spatial'] = coords_df[coord_cols].values
        
        # 添加技术信息
        adata.uns['spatial'] = {
            'technology': {
                'technology': 'Spatial Transcriptomics',
                'spot_diameter': 100,
                'resolution': 'multi-cellular'
            }
        }
        
        # 保存
        output_path = self.processed_dir / f"{output_name}.h5ad"
        adata.write(output_path)
        print(f"  ✅ 保存至: {output_path}")
        
        return adata
    
    def validate_data_integrity(self, adata, dataset_name):
        """验证数据完整性"""
        print(f"\n🔍 验证数据完整性: {dataset_name}")
        
        issues = []
        
        # 基本结构检查
        if adata.n_obs == 0:
            issues.append("❌ 无观测数据(spots)")
        if adata.n_vars == 0:
            issues.append("❌ 无基因数据")
        
        # 空间坐标检查
        if 'spatial' not in adata.obsm:
            issues.append("❌ 缺少空间坐标")
        else:
            coords = adata.obsm['spatial']
            if coords.shape[1] != 2:
                issues.append(f"❌ 空间坐标维度错误: {coords.shape[1]}")
            if np.any(np.isnan(coords)):
                issues.append("⚠️ 空间坐标包含NaN值")
        
        # 表达数据检查
        if np.any(np.isnan(adata.X)):
            issues.append("⚠️ 表达矩阵包含NaN值")
        if np.any(adata.X < 0):
            issues.append("⚠️ 表达矩阵包含负值")
        
        # 统计信息
        sparsity = (adata.X == 0).sum() / (adata.X.shape[0] * adata.X.shape[1])
        total_counts = adata.X.sum()
        
        print(f"  📊 统计信息:")
        print(f"     - 稀疏性: {sparsity:.3f}")
        print(f"     - 总计数: {total_counts:,.0f}")
        print(f"     - 平均每spot计数: {total_counts / adata.n_obs:.1f}")
        print(f"     - 平均每基因计数: {total_counts / adata.n_vars:.1f}")
        
        if issues:
            print(f"  ⚠️ 发现问题:")
            for issue in issues:
                print(f"     {issue}")
        else:
            print(f"  ✅ 数据完整性验证通过")
        
        return len([i for i in issues if i.startswith("❌")]) == 0
    
    def create_quality_report(self, dataset_name):
        """创建数据质量报告"""
        dataset_path = self.data_dir / f"{dataset_name}.h5ad"
        
        if not dataset_path.exists():
            return
        
        adata = sc.read_h5ad(dataset_path)
        
        # 计算质量指标
        adata.var['n_cells'] = (adata.X > 0).sum(axis=0).A1
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
        adata.obs['total_counts'] = adata.X.sum(axis=1).A1
        
        # 创建质量图表
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'数据质量报告: {dataset_name}', fontsize=16)
        
        # 1. 每spot基因数分布
        axes[0, 0].hist(adata.obs['n_genes'], bins=50, alpha=0.7)
        axes[0, 0].set_xlabel('Genes per spot')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Gene count distribution')
        
        # 2. 每spot总计数分布
        axes[0, 1].hist(adata.obs['total_counts'], bins=50, alpha=0.7)
        axes[0, 1].set_xlabel('Total counts per spot') 
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('UMI count distribution')
        
        # 3. 每基因检出细胞数分布
        axes[0, 2].hist(adata.var['n_cells'], bins=50, alpha=0.7)
        axes[0, 2].set_xlabel('Spots per gene')
        axes[0, 2].set_ylabel('Frequency')
        axes[0, 2].set_title('Gene detection distribution')
        
        # 4. 空间分布
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            scatter = axes[1, 0].scatter(coords[:, 0], coords[:, 1], 
                                       c=adata.obs['total_counts'], s=1, alpha=0.6)
            axes[1, 0].set_xlabel('X coordinate')
            axes[1, 0].set_ylabel('Y coordinate')
            axes[1, 0].set_title('Spatial UMI distribution')
            plt.colorbar(scatter, ax=axes[1, 0])
        
        # 5. 基因数空间分布
        if 'spatial' in adata.obsm:
            scatter = axes[1, 1].scatter(coords[:, 0], coords[:, 1],
                                       c=adata.obs['n_genes'], s=1, alpha=0.6)
            axes[1, 1].set_xlabel('X coordinate')
            axes[1, 1].set_ylabel('Y coordinate')
            axes[1, 1].set_title('Spatial gene count distribution')
            plt.colorbar(scatter, ax=axes[1, 1])
        
        # 6. 技术统计
        axes[1, 2].axis('off')
        stats_text = f"""
数据集统计:
• Spots: {adata.n_obs:,}
• Genes: {adata.n_vars:,}
• 稀疏性: {(adata.X == 0).sum() / adata.X.size:.3f}
• 中位基因数/spot: {np.median(adata.obs['n_genes']):.0f}
• 中位UMI数/spot: {np.median(adata.obs['total_counts']):.0f}
• 检出基因数: {(adata.var['n_cells'] > 0).sum():,}
        """
        axes[1, 2].text(0.1, 0.5, stats_text, fontsize=11, verticalalignment='center')
        axes[1, 2].set_title('Dataset Statistics')
        
        plt.tight_layout()
        
        # 保存图表
        plot_path = self.processed_dir / f"{dataset_name}_quality_report.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  📊 质量报告已保存: {plot_path}")
    
    def process_all_datasets(self):
        """处理目录中的所有数据集"""
        print("开始处理所有真实数据集...")
        
        datasets = []
        for h5ad_file in self.data_dir.glob('*.h5ad'):
            dataset_name = h5ad_file.stem
            datasets.append(dataset_name)
        
        print(f"发现 {len(datasets)} 个数据集:")
        for ds in datasets:
            print(f"  - {ds}")
        
        # 处理每个数据集
        results = {}
        for dataset_name in datasets:
            try:
                adata = self.load_and_validate_dataset(dataset_name)
                if adata is not None:
                    valid = self.validate_data_integrity(adata, dataset_name)
                    self.create_quality_report(dataset_name)
                    results[dataset_name] = {'status': 'success', 'valid': valid}
                else:
                    results[dataset_name] = {'status': 'failed', 'valid': False}
                    
            except Exception as e:
                print(f"❌ 处理 {dataset_name} 时出错: {e}")
                results[dataset_name] = {'status': 'error', 'error': str(e)}
        
        # 生成总结报告
        self._generate_processing_summary(results)
        
        return results
    
    def _generate_processing_summary(self, results):
        """生成处理总结报告"""
        summary_path = self.processed_dir / 'PROCESSING_SUMMARY.md'
        
        success_count = sum(1 for r in results.values() if r.get('status') == 'success')
        valid_count = sum(1 for r in results.values() if r.get('valid') == True)
        
        summary_content = f"""# 真实数据集处理总结报告

处理时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## 处理结果概览

- **总数据集**: {len(results)}
- **成功处理**: {success_count}
- **数据完整**: {valid_count}
- **成功率**: {success_count/len(results)*100:.1f}%

## 详细结果

"""
        
        for dataset_name, result in results.items():
            status = result.get('status', 'unknown')
            valid = result.get('valid', False)
            
            if status == 'success':
                status_icon = '✅'
                valid_icon = '✅' if valid else '⚠️'
            elif status == 'failed':
                status_icon = '❌'
                valid_icon = '❌'
            else:
                status_icon = '❌'
                valid_icon = '❌'
            
            summary_content += f"### {dataset_name}\n"
            summary_content += f"- **处理状态**: {status_icon} {status}\n"
            summary_content += f"- **数据完整性**: {valid_icon} {'通过' if valid else '有问题'}\n"
            
            if 'error' in result:
                summary_content += f"- **错误信息**: {result['error']}\n"
            
            summary_content += "\n"
        
        summary_content += """
## 使用建议

### 高质量数据集
推荐用于正式分析的数据集（处理成功且数据完整）

### 问题数据集  
需要进一步清理或修复的数据集

### 质量报告
每个数据集的详细质量报告图表已保存在 `processed/` 目录中

---
*此报告由 SpatialDataProcessor 自动生成*
"""
        
        with open(summary_path, 'w', encoding='utf-8') as f:
            f.write(summary_content)
        
        print(f"\n📋 处理总结报告已保存: {summary_path}")


def main():
    """主函数"""
    print("=" * 80)
    print("真实空间转录组数据集处理器")
    print("=" * 80)
    
    processor = SpatialDataProcessor(BASE_DIR)
    results = processor.process_all_datasets()
    
    print("\n" + "=" * 80)
    print("处理完成！")
    print(f"处理结果保存在: {processor.processed_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()