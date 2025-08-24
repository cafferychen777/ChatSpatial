#!/usr/bin/env python3
"""
Final Real Datasets Integration Report Generator
Linus-style: One script to rule them all. Simple, comprehensive, no bullshit.
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime
import sys

def generate_comprehensive_report(datasets_dir: Path):
    """Generate the final comprehensive report."""
    
    datasets_dir = Path(datasets_dir)
    
    print("ChatSpatial Real Datasets - Final Integration Report")
    print("=" * 60)
    print(f"Report Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Base Directory: {datasets_dir}")
    
    # Load validation results
    validation_file = datasets_dir / 'real_datasets_summary.json'
    compatibility_file = datasets_dir / 'compatibility_test_results.json'
    registry_file = datasets_dir / 'real_datasets_registry.md'
    
    validation_data = {}
    compatibility_data = {}
    
    if validation_file.exists():
        with open(validation_file, 'r') as f:
            validation_data = json.load(f)
    
    if compatibility_file.exists():
        with open(compatibility_file, 'r') as f:
            compatibility_data = json.load(f)
    
    print("\n" + "=" * 60)
    print("EXECUTIVE SUMMARY")
    print("=" * 60)
    
    # Executive summary
    summary_stats = validation_data.get('summary_stats', {})
    
    print(f"✓ DATASETS DISCOVERED: {summary_stats.get('total_datasets', 0)}")
    print(f"✓ SUCCESSFULLY VALIDATED: {summary_stats.get('successful_validations', 0)}")
    print(f"✓ TOTAL CELLS: {summary_stats.get('total_cells', 0):,}")
    print(f"✓ TOTAL GENES: {summary_stats.get('total_genes', 0):,}")
    print(f"✓ TOTAL DATA SIZE: {summary_stats.get('total_size_mb', 0):.1f} MB")
    print(f"✓ SPATIAL DATASETS: {summary_stats.get('spatial_datasets', 0)}")
    
    # Quality distribution
    quality_dist = summary_stats.get('quality_distribution', {})
    if quality_dist:
        print(f"\nQUALITY BREAKDOWN:")
        for tier, count in sorted(quality_dist.items()):
            print(f"  {tier.upper()}: {count} datasets")
    
    # Compatibility results
    compat_tested = compatibility_data.get('tested_datasets', 0)
    compat_success = compatibility_data.get('successful_loads', 0)
    if compat_tested > 0:
        success_rate = compat_success / compat_tested * 100
        print(f"\nCOMPATIBILITY TESTING:")
        print(f"  Datasets tested: {compat_tested}")
        print(f"  Load success rate: {success_rate:.1f}% ({compat_success}/{compat_tested})")
        print(f"  Preprocessing compatible: {compatibility_data.get('preprocessing_compatible', 0)}")
        print(f"  Spatial analysis compatible: {compatibility_data.get('spatial_compatible', 0)}")
    
    print("\n" + "=" * 60)
    print("DATASET CATEGORIES")
    print("=" * 60)
    
    # Analyze dataset distribution by category
    datasets = validation_data.get('datasets', [])
    successful_datasets = [d for d in datasets if d['validation_status'] == 'success']
    
    # Categorize datasets
    categories = {
        'Core Research Datasets': [],
        'Harmony Integration': [],
        'Test/Synthetic': [],
        'Technology-Specific': []
    }
    
    for dataset in successful_datasets:
        dataset_path = dataset['file_path']
        dataset_name = dataset['file_name']
        
        if '/core/' in dataset_path:
            categories['Core Research Datasets'].append(dataset)
        elif '/harmony/' in dataset_path:
            categories['Harmony Integration'].append(dataset)
        elif '/test/' in dataset_path or 'synthetic' in dataset_name.lower():
            categories['Test/Synthetic'].append(dataset)
        else:
            categories['Technology-Specific'].append(dataset)
    
    # Print category summaries
    for category, datasets_list in categories.items():
        if datasets_list:
            print(f"\n{category.upper()}: {len(datasets_list)} datasets")
            
            # Calculate category stats
            total_cells = sum(d.get('n_cells', 0) for d in datasets_list)
            total_size = sum(d.get('file_size_mb', 0) for d in datasets_list)
            spatial_count = sum(1 for d in datasets_list if d.get('has_spatial', False))
            
            print(f"  Cells: {total_cells:,}")
            print(f"  Size: {total_size:.1f} MB") 
            print(f"  Spatial: {spatial_count}/{len(datasets_list)}")
            
            # List top 3 largest datasets in category
            sorted_datasets = sorted(datasets_list, key=lambda x: x.get('n_cells', 0), reverse=True)
            print(f"  Top datasets:")
            for i, dataset in enumerate(sorted_datasets[:3]):
                spatial_icon = "🗺️" if dataset.get('has_spatial', False) else "📊"
                print(f"    {i+1}. {dataset['file_name']:<30} {spatial_icon} "
                      f"{dataset.get('n_cells', 0):,} cells")
    
    print("\n" + "=" * 60)
    print("TECHNOLOGY COVERAGE")
    print("=" * 60)
    
    # Identify different technologies
    technologies = {
        'Visium': [],
        'Slide-seq': [], 
        'MERFISH': [],
        'seqFISH': [],
        'osmFISH': [],
        'STARmap': [],
        'Xenium': [],
        'Generic ST': []
    }
    
    for dataset in successful_datasets:
        name = dataset['file_name'].lower()
        
        if 'visium' in name:
            technologies['Visium'].append(dataset)
        elif 'slideseq' in name:
            technologies['Slide-seq'].append(dataset)
        elif 'merfish' in name:
            technologies['MERFISH'].append(dataset)
        elif 'seqfish' in name:
            technologies['seqFISH'].append(dataset)
        elif 'osmfish' in name:
            technologies['osmFISH'].append(dataset)
        elif 'starmap' in name:
            technologies['STARmap'].append(dataset)
        elif 'xenium' in name:
            technologies['Xenium'].append(dataset)
        elif dataset.get('has_spatial', False):
            technologies['Generic ST'].append(dataset)
    
    # Print technology coverage
    for tech, datasets_list in technologies.items():
        if datasets_list:
            print(f"{tech}: {len(datasets_list)} datasets")
    
    print("\n" + "=" * 60)
    print("DATA QUALITY ASSESSMENT")
    print("=" * 60)
    
    # Quality metrics analysis
    high_quality = [d for d in successful_datasets if d.get('quality_tier') == 'high']
    medium_quality = [d for d in successful_datasets if d.get('quality_tier') == 'medium']
    low_quality = [d for d in successful_datasets if d.get('quality_tier') == 'low']
    
    print(f"HIGH QUALITY ({len(high_quality)} datasets):")
    print(f"  Typical characteristics: >500 genes/cell, >1000 UMIs/cell")
    print(f"  Best for: Production analysis, method development")
    
    print(f"\nMEDIUM QUALITY ({len(medium_quality)} datasets):")
    print(f"  Typical characteristics: 200-500 genes/cell, 500-1000 UMIs/cell")
    print(f"  Best for: Testing workflows, benchmarking")
    
    print(f"\nLOW QUALITY ({len(low_quality)} datasets):")
    print(f"  Typical characteristics: <200 genes/cell, <500 UMIs/cell")
    print(f"  Best for: Stress testing, edge case validation")
    
    print("\n" + "=" * 60)
    print("USAGE RECOMMENDATIONS")
    print("=" * 60)
    
    print("FOR DEVELOPERS:")
    print("  • Start testing with harmony/quick_demo_*.h5ad (fast, reliable)")
    print("  • Use core datasets for comprehensive testing")
    print("  • Test spatial features with slideseq or visium data")
    
    print("\nFOR USERS:")
    print("  • Small analysis: Use datasets with <5K cells")
    print("  • Production workflows: Focus on 'high' quality datasets")
    print("  • Spatial analysis: Prioritize datasets with spatial coordinates")
    
    print("\nFOR BENCHMARKING:")
    print("  • Speed tests: Use perf_test_100_500.h5ad")
    print("  • Memory tests: Use slideseq_MOp_1217.h5ad (large)")
    print("  • Method comparison: Use multiple datasets from same category")
    
    print("\n" + "=" * 60)
    print("KNOWN ISSUES & LIMITATIONS")
    print("=" * 60)
    
    # Failed datasets
    failed_datasets = [d for d in datasets if d['validation_status'] == 'failed']
    if failed_datasets:
        print(f"FAILED VALIDATIONS ({len(failed_datasets)} datasets):")
        for dataset in failed_datasets[:5]:  # Show first 5
            print(f"  • {dataset['file_name']}: {dataset.get('error', 'Unknown error')[:60]}...")
    
    # Compatibility issues
    issues = compatibility_data.get('issues_found', [])
    if issues:
        print(f"\nCOMPATIBILITY ISSUES ({len(issues)} found):")
        for issue in issues[:5]:  # Show first 5
            print(f"  • {issue['dataset']}: {issue['issue']}")
    
    print("\n" + "=" * 60)
    print("FILES GENERATED")
    print("=" * 60)
    
    generated_files = [
        ('real_datasets_summary.json', 'Machine-readable validation results'),
        ('real_datasets_registry.md', 'Human-readable dataset catalog'),
        ('compatibility_test_results.json', 'ChatSpatial compatibility results'),
        ('download_log.txt', 'Processing history log'),
        ('validate_real_datasets.py', 'Dataset validation system'),
        ('test_compatibility_simple.py', 'Simple compatibility tester')
    ]
    
    for filename, description in generated_files:
        file_path = datasets_dir / filename
        status = "✓" if file_path.exists() else "✗"
        size = f"({file_path.stat().st_size / 1024:.1f} KB)" if file_path.exists() else ""
        print(f"{status} {filename:<35} {size:<10} - {description}")
    
    print("\n" + "=" * 60)
    print("QUICK START COMMANDS")
    print("=" * 60)
    
    print("# Load a high-quality spatial dataset")
    print("import scanpy as sc")
    print("adata = sc.read_h5ad('datasets/real_datasets/core/slideseq_MOp_1217.h5ad')")
    print()
    print("# Quick test with harmony data") 
    print("adata = sc.read_h5ad('datasets/real_datasets/harmony/quick_demo_combined.h5ad')")
    print()
    print("# Validate dataset structure")
    print("python datasets/real_datasets/validate_real_datasets.py")
    
    print("\n" + "=" * 60)
    print("INTEGRATION COMPLETE")
    print("=" * 60)
    
    success_rate = summary_stats.get('successful_validations', 0) / summary_stats.get('total_datasets', 1) * 100
    
    if success_rate > 80:
        print("🎉 EXCELLENT: >80% datasets validated successfully")
    elif success_rate > 60:
        print("✅ GOOD: >60% datasets validated successfully")  
    else:
        print("⚠️ ISSUES: <60% validation success rate")
    
    print(f"\nFinal Status: {summary_stats.get('successful_validations', 0)}/{summary_stats.get('total_datasets', 0)} datasets ready for use")
    print("Real datasets integration system is operational.")


def main():
    """Generate the comprehensive report."""
    
    datasets_dir = Path(__file__).parent
    
    try:
        generate_comprehensive_report(datasets_dir)
        print("\n✓ Final report generation completed!")
        
    except Exception as e:
        print(f"\n✗ Report generation failed: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()