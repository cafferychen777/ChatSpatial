#!/usr/bin/env python3
"""
深度整理 ChatSpatial data 目录
- 处理重复文件
- 清理旧目录结构  
- 规范命名
- 重新分类数据集
- 创建完整索引
"""

import os
import shutil
from pathlib import Path
import json
import scanpy as sc
import pandas as pd
import warnings
from collections import defaultdict
import hashlib
warnings.filterwarnings('ignore')

class DataDirectoryCleaner:
    """深度清理data目录"""
    
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        self.duplicates_found = []
        self.moved_files = []
        self.removed_files = []
        self.errors = []
        
        # 标准目录结构
        self.standard_dirs = {
            'spatial_datasets': '空间转录组数据集',
            'single_cell_datasets': '单细胞RNA-seq数据集', 
            'synthetic_datasets': '合成和模拟数据集',
            'demo_datasets': '演示和教程数据集',
            'benchmark_datasets': '性能基准测试数据集',
            'reference_datasets': '参考和标准数据集',
            'metadata': '元数据和索引文件',
            'documentation': '文档和说明',
            'scripts': '数据处理脚本',
            'archive': '归档和备份文件'
        }
    
    def find_duplicate_files(self):
        """查找重复文件"""
        print("🔍 查找重复文件...")
        
        file_hashes = defaultdict(list)
        file_sizes = defaultdict(list)
        
        # 收集所有h5ad文件
        for h5ad_file in self.data_dir.rglob("*.h5ad"):
            if h5ad_file.is_file():
                # 计算文件大小
                size = h5ad_file.stat().st_size
                file_sizes[size].append(h5ad_file)
        
        # 对相同大小的文件计算哈希
        for size, files in file_sizes.items():
            if len(files) > 1:
                for file_path in files:
                    try:
                        with open(file_path, 'rb') as f:
                            # 只读取前1MB来计算哈希，加速处理
                            content = f.read(1024 * 1024)
                            hash_value = hashlib.md5(content).hexdigest()
                            file_hashes[hash_value].append(file_path)
                    except Exception as e:
                        self.errors.append(f"计算哈希失败 {file_path}: {e}")
        
        # 找出真正的重复文件
        for hash_value, files in file_hashes.items():
            if len(files) > 1:
                self.duplicates_found.append({
                    'hash': hash_value,
                    'files': files,
                    'size_mb': files[0].stat().st_size / (1024*1024)
                })
                print(f"  发现重复: {[f.name for f in files]}")
        
        print(f"找到 {len(self.duplicates_found)} 组重复文件")
        
    def handle_duplicates(self):
        """处理重复文件"""
        print("🧹 处理重复文件...")
        
        for dup_group in self.duplicates_found:
            files = dup_group['files']
            
            # 选择保留的文件（优先选择在正确分类目录中的）
            keep_file = None
            remove_files = []
            
            # 优先级：标准目录 > 非backup > 文件名简洁
            for file_path in files:
                is_in_standard_dir = any(std_dir in str(file_path) for std_dir in self.standard_dirs.keys())
                is_backup = 'backup' in file_path.name.lower()
                
                if keep_file is None:
                    keep_file = file_path
                else:
                    # 当前文件更好的条件
                    current_better = (
                        is_in_standard_dir and not any(std_dir in str(keep_file) for std_dir in self.standard_dirs.keys()) or
                        not is_backup and 'backup' in keep_file.name.lower() or
                        len(file_path.name) < len(keep_file.name)
                    )
                    
                    if current_better:
                        remove_files.append(keep_file)
                        keep_file = file_path
                    else:
                        remove_files.append(file_path)
            
            # 移动重复文件到archive
            archive_dir = self.data_dir / 'archive' / 'duplicates'
            archive_dir.mkdir(parents=True, exist_ok=True)
            
            for remove_file in remove_files:
                archive_path = archive_dir / remove_file.name
                counter = 1
                while archive_path.exists():
                    stem = remove_file.stem
                    suffix = remove_file.suffix
                    archive_path = archive_dir / f"{stem}_dup{counter}{suffix}"
                    counter += 1
                
                try:
                    shutil.move(remove_file, archive_path)
                    self.removed_files.append(str(remove_file))
                    print(f"  归档重复文件: {remove_file.name} -> archive/duplicates/")
                except Exception as e:
                    self.errors.append(f"归档失败 {remove_file}: {e}")
    
    def categorize_dataset_smart(self, file_path):
        """智能分类数据集"""
        filename = file_path.name.lower()
        
        try:
            # 尝试读取数据集获取更多信息
            adata = sc.read_h5ad(file_path)
            has_spatial = 'spatial' in adata.obsm
            n_cells = adata.n_obs
            n_genes = adata.n_vars
            
        except:
            has_spatial = False
            n_cells = 0
            n_genes = 0
        
        # 智能分类规则
        if any(x in filename for x in ['benchmark', 'perf', 'stress']) or (100 <= n_cells <= 5000 and 500 <= n_genes <= 5000):
            return 'benchmark_datasets'
        elif any(x in filename for x in ['demo', 'tutorial', 'quick', 'example']) or 'mop_sn_tutorial' in filename:
            return 'demo_datasets'
        elif any(x in filename for x in ['synthetic', 'simulated', 'artificial', 'generated']):
            return 'synthetic_datasets'
        elif any(x in filename for x in ['paul15', 'reference', 'standard']) or n_cells > 0 and n_cells < 1000 and 'test' in filename:
            return 'reference_datasets'
        elif has_spatial or any(x in filename for x in ['visium', 'slideseq', 'merfish', 'seqfish', 'spatial', 'st_']):
            return 'spatial_datasets'
        elif n_cells > 0 and not has_spatial:
            return 'single_cell_datasets'
        else:
            return 'spatial_datasets'  # 默认归类为空间数据
    
    def reorganize_datasets(self):
        """重新组织数据集"""
        print("📁 重新组织数据集分类...")
        
        # 创建标准目录
        for dir_name, description in self.standard_dirs.items():
            dir_path = self.data_dir / dir_name
            dir_path.mkdir(exist_ok=True)
        
        # 收集所有需要分类的文件
        files_to_process = []
        
        # 从旧的分类目录收集文件
        old_dirs = ['real_datasets', 'benchmark_datasets', 'demo_datasets', 'synthetic_datasets', 'reference_datasets']
        
        for old_dir in old_dirs:
            old_path = self.data_dir / old_dir
            if old_path.exists():
                for h5ad_file in old_path.glob("*.h5ad"):
                    files_to_process.append(h5ad_file)
        
        # 处理根目录和其他位置的文件
        for h5ad_file in self.data_dir.rglob("*.h5ad"):
            if h5ad_file.parent.name not in self.standard_dirs and h5ad_file not in files_to_process:
                files_to_process.append(h5ad_file)
        
        # 重新分类每个文件
        moved_count = 0
        for file_path in files_to_process:
            if not file_path.exists():
                continue
                
            new_category = self.categorize_dataset_smart(file_path)
            target_dir = self.data_dir / new_category
            target_path = target_dir / file_path.name
            
            # 检查是否需要移动
            if file_path.parent != target_dir:
                # 处理同名文件
                if target_path.exists():
                    counter = 1
                    stem = file_path.stem
                    suffix = file_path.suffix
                    while target_path.exists():
                        target_path = target_dir / f"{stem}_v{counter}{suffix}"
                        counter += 1
                
                try:
                    shutil.move(file_path, target_path)
                    self.moved_files.append(f"{file_path.name} -> {new_category}/")
                    moved_count += 1
                    print(f"  移动: {file_path.name} -> {new_category}/")
                except Exception as e:
                    self.errors.append(f"移动失败 {file_path}: {e}")
        
        print(f"重新分类移动了 {moved_count} 个文件")
    
    def clean_old_directories(self):
        """清理旧的目录结构"""
        print("🗂️ 清理旧目录结构...")
        
        old_dirs_to_clean = ['real_datasets', 'core', 'paul15', 'test', 'test_datasets']
        archive_dir = self.data_dir / 'archive'
        archive_dir.mkdir(exist_ok=True)
        
        for old_dir_name in old_dirs_to_clean:
            old_dir = self.data_dir / old_dir_name
            if old_dir.exists():
                # 检查目录是否为空或只含非h5ad文件
                h5ad_files = list(old_dir.rglob("*.h5ad"))
                
                if not h5ad_files:
                    # 目录没有h5ad文件，可以移动到archive
                    archive_path = archive_dir / f"old_{old_dir_name}"
                    if archive_path.exists():
                        shutil.rmtree(archive_path)
                    
                    try:
                        shutil.move(old_dir, archive_path)
                        print(f"  归档旧目录: {old_dir_name} -> archive/")
                    except Exception as e:
                        self.errors.append(f"归档目录失败 {old_dir}: {e}")
                else:
                    print(f"  保留 {old_dir_name} (仍有 {len(h5ad_files)} 个数据文件)")
    
    def organize_harmony_data(self):
        """整理harmony相关数据"""
        print("🎵 整理Harmony数据...")
        
        harmony_dir = self.data_dir / 'harmony'
        if not harmony_dir.exists():
            return
        
        # 移动Python脚本到scripts目录
        scripts_dir = self.data_dir / 'scripts' / 'harmony'
        scripts_dir.mkdir(parents=True, exist_ok=True)
        
        for py_file in harmony_dir.glob("*.py"):
            target_path = scripts_dir / py_file.name
            if not target_path.exists():
                try:
                    shutil.move(py_file, target_path)
                    print(f"  移动脚本: {py_file.name} -> scripts/harmony/")
                except Exception as e:
                    self.errors.append(f"移动脚本失败 {py_file}: {e}")
        
        # 移动图片到documentation
        docs_dir = self.data_dir / 'documentation' / 'harmony_figures'
        figures_dir = harmony_dir / 'figures'
        if figures_dir.exists():
            docs_dir.mkdir(parents=True, exist_ok=True)
            try:
                for fig_file in figures_dir.iterdir():
                    target_path = docs_dir / fig_file.name
                    if not target_path.exists():
                        shutil.move(fig_file, target_path)
                print(f"  移动图片: figures/ -> documentation/harmony_figures/")
            except Exception as e:
                self.errors.append(f"移动图片失败: {e}")
        
        # 保留README和重要的h5ad文件
        for h5ad_file in harmony_dir.rglob("*.h5ad"):
            # 判断是否为演示数据
            if 'demo' in h5ad_file.name.lower() or 'quick' in h5ad_file.name.lower():
                target_dir = self.data_dir / 'demo_datasets'
                target_path = target_dir / h5ad_file.name
                
                if not target_path.exists():
                    try:
                        shutil.move(h5ad_file, target_path)
                        print(f"  移动数据: {h5ad_file.name} -> demo_datasets/")
                    except Exception as e:
                        self.errors.append(f"移动数据失败 {h5ad_file}: {e}")
    
    def create_comprehensive_index(self):
        """创建全面的数据集索引"""
        print("📊 创建综合数据集索引...")
        
        catalog = {
            'metadata': {
                'created_at': pd.Timestamp.now().isoformat(),
                'total_datasets': 0,
                'total_size_gb': 0,
                'directories': {}
            },
            'datasets_by_category': {},
            'datasets_by_technology': defaultdict(list),
            'datasets_by_size': {
                'small': [],      # < 1000 cells
                'medium': [],     # 1000-10000 cells  
                'large': [],      # 10000-50000 cells
                'xl': []          # > 50000 cells
            },
            'spatial_datasets': [],
            'quick_access': {
                'recommended_spatial': [],
                'recommended_demo': [],
                'recommended_benchmark': []
            }
        }
        
        total_size_bytes = 0
        total_datasets = 0
        
        # 扫描每个标准目录
        for dir_name in self.standard_dirs.keys():
            dir_path = self.data_dir / dir_name
            if not dir_path.exists():
                continue
            
            datasets = []
            
            for h5ad_file in dir_path.glob("*.h5ad"):
                try:
                    # 分析数据集
                    adata = sc.read_h5ad(h5ad_file)
                    file_size = h5ad_file.stat().st_size
                    total_size_bytes += file_size
                    total_datasets += 1
                    
                    dataset_info = {
                        'filename': h5ad_file.name,
                        'path': f"{dir_name}/{h5ad_file.name}",
                        'n_cells': adata.n_obs,
                        'n_genes': adata.n_vars,
                        'has_spatial': 'spatial' in adata.obsm,
                        'spatial_dims': adata.obsm['spatial'].shape if 'spatial' in adata.obsm else None,
                        'has_clusters': len(adata.obs.select_dtypes(include=['category', 'object']).columns) > 0,
                        'has_leiden': 'leiden' in adata.obs.columns,
                        'sparsity': float((adata.X == 0).sum() / adata.X.size) if hasattr(adata.X, 'size') else 0.0,
                        'file_size_mb': float(file_size / (1024 * 1024)),
                        'category': dir_name
                    }
                    
                    # 推断技术类型
                    dataset_info['technology'] = self.infer_technology(h5ad_file.name, adata)
                    
                    datasets.append(dataset_info)
                    
                    # 分类到大小类别
                    n_cells = adata.n_obs
                    if n_cells < 1000:
                        catalog['datasets_by_size']['small'].append(dataset_info)
                    elif n_cells < 10000:
                        catalog['datasets_by_size']['medium'].append(dataset_info)
                    elif n_cells < 50000:
                        catalog['datasets_by_size']['large'].append(dataset_info)
                    else:
                        catalog['datasets_by_size']['xl'].append(dataset_info)
                    
                    # 技术分类
                    catalog['datasets_by_technology'][dataset_info['technology']].append(dataset_info)
                    
                    # 空间数据集
                    if dataset_info['has_spatial']:
                        catalog['spatial_datasets'].append(dataset_info)
                    
                except Exception as e:
                    self.errors.append(f"分析数据集失败 {h5ad_file}: {e}")
                    continue
            
            catalog['datasets_by_category'][dir_name] = datasets
            catalog['metadata']['directories'][dir_name] = {
                'count': len(datasets),
                'description': self.standard_dirs[dir_name]
            }
        
        # 更新元数据
        catalog['metadata']['total_datasets'] = total_datasets
        catalog['metadata']['total_size_gb'] = total_size_bytes / (1024**3)
        
        # 生成推荐列表
        spatial_sorted = sorted(catalog['spatial_datasets'], key=lambda x: x['n_cells'], reverse=True)
        catalog['quick_access']['recommended_spatial'] = spatial_sorted[:10]
        
        demo_datasets = catalog['datasets_by_category'].get('demo_datasets', [])
        catalog['quick_access']['recommended_demo'] = sorted(demo_datasets, key=lambda x: x['n_cells'])[:5]
        
        benchmark_datasets = catalog['datasets_by_category'].get('benchmark_datasets', [])
        catalog['quick_access']['recommended_benchmark'] = sorted(benchmark_datasets, key=lambda x: x['n_cells'])
        
        # 保存索引文件
        metadata_dir = self.data_dir / 'metadata'
        metadata_dir.mkdir(exist_ok=True)
        
        with open(metadata_dir / 'comprehensive_catalog.json', 'w') as f:
            json.dump(catalog, f, indent=2)
        
        print(f"创建了包含 {total_datasets} 个数据集的综合索引")
        return catalog
    
    def infer_technology(self, filename, adata):
        """推断空间转录组技术"""
        filename = filename.lower()
        
        if 'visium' in filename:
            return 'visium'
        elif 'slideseq' in filename:
            return 'slide-seq'
        elif 'merfish' in filename:
            return 'merfish'
        elif 'seqfish' in filename:
            return 'seqfish'
        elif 'osmfish' in filename:
            return 'osmfish'
        elif 'st_' in filename or 'spatial' in filename:
            return 'spatial_transcriptomics'
        elif 'synthetic' in filename or 'simulated' in filename:
            return 'synthetic'
        elif 'spatial' in adata.obsm:
            # 根据数据特征推断
            n_spots = adata.n_obs
            if n_spots > 30000:
                return 'slide-seq'
            elif n_spots > 2000:
                return 'visium'
            else:
                return 'unknown_spatial'
        else:
            return 'single_cell_rna_seq'
    
    def create_documentation(self, catalog):
        """创建完整文档"""
        print("📚 创建文档...")
        
        docs_dir = self.data_dir / 'documentation'
        docs_dir.mkdir(exist_ok=True)
        
        # 主索引文档
        self.create_main_index(catalog)
        
        # 快速开始指南
        self.create_quickstart_guide(catalog)
        
        # API参考
        self.create_api_reference(catalog)
        
        # 更新各目录的README
        self.update_directory_readmes(catalog)
    
    def create_main_index(self, catalog):
        """创建主索引文档"""
        docs_dir = self.data_dir / 'documentation'
        
        with open(docs_dir / 'MAIN_INDEX.md', 'w') as f:
            f.write(f"# ChatSpatial 数据集完整索引\n\n")
            f.write(f"**更新时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**数据集总数**: {catalog['metadata']['total_datasets']}\n")
            f.write(f"**总存储空间**: {catalog['metadata']['total_size_gb']:.1f} GB\n\n")
            
            f.write("## 📊 按类别统计\n\n")
            for dir_name, info in catalog['metadata']['directories'].items():
                f.write(f"- **{dir_name.replace('_', ' ').title()}**: {info['count']} 个数据集\n")
            
            f.write(f"\n## 🔬 按技术分类\n\n")
            for tech, datasets in catalog['datasets_by_technology'].items():
                f.write(f"- **{tech.replace('_', ' ').title()}**: {len(datasets)} 个数据集\n")
            
            f.write(f"\n## 📏 按规模分类\n\n")
            for size_cat, datasets in catalog['datasets_by_size'].items():
                f.write(f"- **{size_cat.title()}**: {len(datasets)} 个数据集\n")
            
            f.write(f"\n## 🎯 推荐数据集\n\n")
            
            f.write("### 空间分析推荐\n")
            for i, dataset in enumerate(catalog['quick_access']['recommended_spatial'][:5], 1):
                f.write(f"{i}. **{dataset['filename']}** - {dataset['n_cells']:,} cells ({dataset['technology']})\n")
            
            f.write(f"\n### 演示教程推荐\n")
            for i, dataset in enumerate(catalog['quick_access']['recommended_demo'][:3], 1):
                f.write(f"{i}. **{dataset['filename']}** - {dataset['n_cells']:,} cells\n")
            
            f.write(f"\n## 📖 详细目录\n\n")
            for dir_name in catalog['metadata']['directories'].keys():
                f.write(f"- [{dir_name.replace('_', ' ').title()}](../{dir_name}/README.md)\n")
    
    def create_quickstart_guide(self, catalog):
        """创建快速开始指南"""
        docs_dir = self.data_dir / 'documentation'
        
        with open(docs_dir / 'QUICKSTART.md', 'w') as f:
            f.write("# ChatSpatial 数据集快速开始\n\n")
            
            f.write("## 🚀 立即开始\n\n")
            f.write("```python\nimport scanpy as sc\n\n")
            
            # 推荐几个典型使用案例
            if catalog['quick_access']['recommended_demo']:
                demo = catalog['quick_access']['recommended_demo'][0]
                f.write(f"# 1. 快速演示\n")
                f.write(f"adata = sc.read_h5ad('data/demo_datasets/{demo['filename']}')\n")
                f.write(f"print(f'Demo data: {{adata.n_obs}} cells, {{adata.n_vars}} genes')\n\n")
            
            if catalog['quick_access']['recommended_spatial']:
                spatial = catalog['quick_access']['recommended_spatial'][0]
                f.write(f"# 2. 空间分析\n")
                f.write(f"adata = sc.read_h5ad('data/spatial_datasets/{spatial['filename']}')\n")
                f.write(f"print(f'Spatial data: {{adata.n_obs}} spots, {{adata.n_vars}} genes')\n")
                f.write(f"print('Spatial coordinates:', adata.obsm['spatial'].shape)\n\n")
            
            if catalog['datasets_by_size']['large']:
                large = catalog['datasets_by_size']['large'][0]
                f.write(f"# 3. 大规模数据\n")
                f.write(f"adata = sc.read_h5ad('data/{large['category']}/{large['filename']}')\n")
                f.write(f"print(f'Large dataset: {{adata.n_obs}} cells')\n")
            
            f.write("```\n\n")
            
            f.write("## 📋 数据集选择指南\n\n")
            f.write("| 用途 | 推荐数据集 | 位置 |\n")
            f.write("|------|------------|------|\n")
            
            for purpose, datasets, location in [
                ('快速测试', catalog['datasets_by_size']['small'][:3], 'demo_datasets/'),
                ('空间分析', catalog['quick_access']['recommended_spatial'][:3], 'spatial_datasets/'),
                ('性能测试', catalog['datasets_by_size']['large'][:3], 'benchmark_datasets/'),
            ]:
                for dataset in datasets:
                    f.write(f"| {purpose} | {dataset['filename']} | {dataset['category']}/ |\n")
    
    def create_api_reference(self, catalog):
        """创建API参考文档"""
        docs_dir = self.data_dir / 'documentation'
        
        with open(docs_dir / 'API_REFERENCE.md', 'w') as f:
            f.write("# ChatSpatial 数据集 API 参考\n\n")
            
            f.write("## 数据加载 API\n\n")
            f.write("```python\n")
            f.write("# 基本加载\n")
            f.write("import scanpy as sc\n")
            f.write("adata = sc.read_h5ad('data/spatial_datasets/dataset.h5ad')\n\n")
            
            f.write("# 批量加载\n")
            f.write("import json\n")
            f.write("with open('data/metadata/comprehensive_catalog.json') as f:\n")
            f.write("    catalog = json.load(f)\n\n")
            
            f.write("# 按类别加载\n")
            f.write("spatial_datasets = catalog['datasets_by_category']['spatial_datasets']\n")
            f.write("for dataset in spatial_datasets:\n")
            f.write("    adata = sc.read_h5ad(f'data/{dataset[\"path\"]}')\n")
            f.write("```\n\n")
            
            f.write("## 数据集查询 API\n\n")
            f.write("```python\n")
            f.write("# 查询空间数据集\n")
            f.write("spatial_datasets = [d for d in catalog['spatial_datasets'] if d['has_spatial']]\n\n")
            
            f.write("# 查询大规模数据集\n")
            f.write("large_datasets = [d for d in catalog['datasets_by_size']['large']]\n\n")
            
            f.write("# 查询特定技术\n")
            f.write("visium_datasets = catalog['datasets_by_technology']['visium']\n")
            f.write("```\n")
    
    def update_directory_readmes(self, catalog):
        """更新各目录的README文件"""
        for dir_name, datasets in catalog['datasets_by_category'].items():
            if not datasets:
                continue
                
            readme_path = self.data_dir / dir_name / 'README.md'
            
            with open(readme_path, 'w') as f:
                f.write(f"# {dir_name.replace('_', ' ').title()}\n\n")
                f.write(f"{self.standard_dirs.get(dir_name, '数据集目录')}\n\n")
                
                # 统计信息
                total_cells = sum(d['n_cells'] for d in datasets)
                total_genes = sum(d['n_genes'] for d in datasets)  
                total_size = sum(d['file_size_mb'] for d in datasets)
                spatial_count = sum(1 for d in datasets if d['has_spatial'])
                
                f.write(f"## 📊 统计信息\n\n")
                f.write(f"- **数据集数量**: {len(datasets)}\n")
                f.write(f"- **总细胞数**: {total_cells:,}\n")
                f.write(f"- **总基因数**: {total_genes:,}\n")
                f.write(f"- **总文件大小**: {total_size:.1f} MB\n")
                f.write(f"- **空间数据集**: {spatial_count}/{len(datasets)}\n\n")
                
                # 数据集列表
                f.write(f"## 📋 数据集列表\n\n")
                f.write("| 数据集 | 细胞数 | 基因数 | 大小(MB) | 技术 | 空间 | 聚类 |\n")
                f.write("|--------|--------|--------|----------|------|------|------|\n")
                
                for dataset in sorted(datasets, key=lambda x: x['n_cells'], reverse=True):
                    spatial_icon = "✅" if dataset['has_spatial'] else "❌"
                    cluster_icon = "✅" if dataset['has_leiden'] else "❌"
                    
                    f.write(f"| {dataset['filename']} | {dataset['n_cells']:,} | {dataset['n_genes']:,} | {dataset['file_size_mb']:.1f} | {dataset['technology']} | {spatial_icon} | {cluster_icon} |\n")
    
    def run_deep_cleanup(self):
        """运行深度清理"""
        print("🧹 开始深度清理 ChatSpatial data 目录")
        print("=" * 60)
        
        try:
            # 1. 查找和处理重复文件
            self.find_duplicate_files()
            self.handle_duplicates()
            
            # 2. 重新组织数据集分类
            self.reorganize_datasets()
            
            # 3. 整理特殊目录
            self.organize_harmony_data()
            
            # 4. 清理旧目录结构
            self.clean_old_directories()
            
            # 5. 创建comprehensive索引
            catalog = self.create_comprehensive_index()
            
            # 6. 创建完整文档
            self.create_documentation(catalog)
            
            # 7. 生成最终报告
            self.generate_final_report(catalog)
            
        except Exception as e:
            self.errors.append(f"深度清理过程出错: {e}")
            raise
    
    def generate_final_report(self, catalog):
        """生成最终清理报告"""
        print("\n" + "=" * 60)
        print("深度清理完成报告")
        print("=" * 60)
        
        print(f"📊 数据集统计:")
        print(f"  - 总数据集: {catalog['metadata']['total_datasets']}")
        print(f"  - 总大小: {catalog['metadata']['total_size_gb']:.1f} GB")
        print(f"  - 空间数据集: {len(catalog['spatial_datasets'])}")
        
        print(f"\n🔄 处理统计:")
        print(f"  - 发现重复文件组: {len(self.duplicates_found)}")
        print(f"  - 移动文件: {len(self.moved_files)}")
        print(f"  - 归档文件: {len(self.removed_files)}")
        
        if self.errors:
            print(f"\n⚠️  错误 ({len(self.errors)}):")
            for error in self.errors[:5]:  # 只显示前5个错误
                print(f"  - {error}")
        
        print(f"\n📍 新的目录结构:")
        for dir_name, info in catalog['metadata']['directories'].items():
            if info['count'] > 0:
                print(f"  - {dir_name}: {info['count']} 个数据集")
        
        print(f"\n📚 文档位置:")
        print(f"  - 主索引: data/documentation/MAIN_INDEX.md")
        print(f"  - 快速开始: data/documentation/QUICKSTART.md")
        print(f"  - API参考: data/documentation/API_REFERENCE.md")
        print(f"  - 详细目录: data/metadata/comprehensive_catalog.json")


def main():
    """主函数"""
    data_dir = "/Users/apple/Research/SpatialTrans_MCP/chatspatial/data"
    
    cleaner = DataDirectoryCleaner(data_dir)
    cleaner.run_deep_cleanup()


if __name__ == "__main__":
    main()