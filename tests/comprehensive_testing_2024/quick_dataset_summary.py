#!/usr/bin/env python3
"""
快速数据集摘要 - Linus原则：简单高效
"""

import os
from pathlib import Path
import time

def quick_scan(base_dir: Path):
    """
    快速扫描h5ad文件，只统计基本信息
    """
    h5ad_files = list(base_dir.rglob("*.h5ad"))
    
    total_files = len(h5ad_files)
    total_size = sum(f.stat().st_size for f in h5ad_files) / (1024 * 1024)  # MB
    
    # 按目录分类
    categories = {}
    for f in h5ad_files:
        rel_path = f.relative_to(base_dir)
        category = str(rel_path.parent) if rel_path.parent != Path(".") else "root"
        
        if category not in categories:
            categories[category] = {'count': 0, 'size_mb': 0, 'files': []}
        
        categories[category]['count'] += 1
        categories[category]['size_mb'] += f.stat().st_size / (1024 * 1024)
        categories[category]['files'].append(f.name)
    
    return {
        'total_files': total_files,
        'total_size_mb': round(total_size, 1),
        'categories': categories,
        'scan_time': time.strftime('%Y-%m-%d %H:%M:%S')
    }

def main():
    print("快速数据集摘要")
    print("=" * 30)
    
    base_dir = Path("/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/datasets/real_datasets")
    
    if not base_dir.exists():
        print(f"目录不存在: {base_dir}")
        return
    
    results = quick_scan(base_dir)
    
    print(f"扫描完成时间: {results['scan_time']}")
    print(f"总文件数: {results['total_files']}")
    print(f"总大小: {results['total_size_mb']} MB ({results['total_size_mb']/1024:.1f} GB)")
    print(f"平均文件大小: {results['total_size_mb']/results['total_files']:.1f} MB")
    
    print("\n按类别统计:")
    print("-" * 50)
    print(f"{'类别':<20} {'文件数':<8} {'大小(MB)':<12} {'示例文件'}")
    print("-" * 50)
    
    for category, info in sorted(results['categories'].items()):
        example_file = info['files'][0] if info['files'] else ""
        if len(example_file) > 30:
            example_file = example_file[:27] + "..."
        
        print(f"{category:<20} {info['count']:<8} {info['size_mb']:<12.1f} {example_file}")
    
    # 保存简单摘要
    summary_file = base_dir / "QUICK_SUMMARY.md"
    with open(summary_file, 'w') as f:
        f.write("# 数据集快速摘要\n\n")
        f.write(f"扫描时间: {results['scan_time']}\n\n")
        f.write(f"- **总文件数**: {results['total_files']}\n")
        f.write(f"- **总大小**: {results['total_size_mb']} MB ({results['total_size_mb']/1024:.1f} GB)\n")
        f.write(f"- **平均大小**: {results['total_size_mb']/results['total_files']:.1f} MB\n\n")
        
        f.write("## 按类别统计\n\n")
        for category, info in sorted(results['categories'].items()):
            f.write(f"### {category}\n")
            f.write(f"- 文件数: {info['count']}\n")
            f.write(f"- 大小: {info['size_mb']:.1f} MB\n")
            f.write(f"- 文件列表:\n")
            for file_name in sorted(info['files']):
                f.write(f"  - {file_name}\n")
            f.write("\n")
    
    print(f"\n摘要保存到: {summary_file}")

if __name__ == "__main__":
    main()