#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple test runner for all comprehensive tools testing.
Runs basic imports and generates reports.
"""

import sys
import os
from pathlib import Path
import json
import warnings
import traceback

# Add project root to path
project_root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(project_root))

warnings.filterwarnings('ignore')

def test_basic_imports():
    """Test basic imports for all modules"""
    print("Testing basic imports for all modules...")
    
    results = {}
    
    modules_to_test = [
        ('spatial_genes', 'chatspatial.tools.spatial_genes'),
        ('spatial_registration', 'chatspatial.tools.spatial_registration'),
        ('spatial_domains', 'chatspatial.tools.spatial_domains'),
        ('integration', 'chatspatial.tools.integration'),
        ('pathway_enrichment', 'chatspatial.tools.pathway_enrichment'),
        ('preprocessing', 'chatspatial.tools.preprocessing'),
        ('spatial_analysis', 'chatspatial.tools.spatial_analysis'),
        ('spatial_enrichment', 'chatspatial.tools.spatial_enrichment'),
        ('gene_set_loader', 'chatspatial.tools.gene_set_loader'),
        ('differential', 'chatspatial.tools.differential')
    ]
    
    for module_name, import_path in modules_to_test:
        try:
            __import__(import_path)
            results[module_name] = "SUCCESS"
            print("  SUCCESS {}: Import successful".format(module_name))
        except Exception as e:
            results[module_name] = "FAILED: {}".format(str(e))
            print("  FAILED {}: {}".format(module_name, str(e)))
    
    return results

def check_data_availability():
    """Check available test datasets"""
    print("\nChecking test datasets...")
    
    data_dir = project_root / "data" / "demo_datasets"
    datasets = []
    
    if data_dir.exists():
        for file_path in data_dir.glob("*.h5ad"):
            try:
                # Just check file size, don't load
                size_mb = file_path.stat().st_size / (1024 * 1024)
                datasets.append({
                    'name': file_path.name,
                    'size_mb': round(size_mb, 2)
                })
                print("  SUCCESS {}: {:.1f} MB".format(file_path.name, size_mb))
            except Exception as e:
                print("  WARNING {}: Cannot access".format(file_path.name))
    else:
        print("  ERROR Data directory not found: {}".format(data_dir))
    
    return datasets

def check_dependencies():
    """Check key dependencies"""
    print("\n🔍 Checking key dependencies...")
    
    dependencies = {
        'scanpy': 'scanpy',
        'pandas': 'pandas', 
        'numpy': 'numpy',
        'scipy': 'scipy',
        'sklearn': 'sklearn',
        'matplotlib': 'matplotlib'
    }
    
    results = {}
    
    for dep_name, import_name in dependencies.items():
        try:
            module = __import__(import_name)
            version = getattr(module, '__version__', 'unknown')
            results[dep_name] = f"✅ {version}"
            print(f"  ✅ {dep_name}: {version}")
        except ImportError:
            results[dep_name] = "❌ Not available"
            print(f"  ❌ {dep_name}: Not available")
    
    return results

def generate_summary_report():
    """Generate comprehensive summary report"""
    print("\n📄 Generating comprehensive testing summary...")
    
    # Run all checks
    import_results = test_basic_imports()
    datasets = check_data_availability()
    dependencies = check_dependencies()
    
    # Count successes
    successful_imports = len([r for r in import_results.values() if r.startswith('✅')])
    total_modules = len(import_results)
    
    successful_deps = len([r for r in dependencies.values() if r.startswith('✅')])
    total_deps = len(dependencies)
    
    summary = {
        "test_timestamp": str(Path(__file__).stat().st_mtime),
        "summary": {
            "total_modules_tested": total_modules,
            "successful_imports": successful_imports,
            "import_success_rate": f"{successful_imports}/{total_modules} ({successful_imports/total_modules*100:.1f}%)",
            "dependencies_available": f"{successful_deps}/{total_deps}",
            "datasets_available": len(datasets)
        },
        "module_results": import_results,
        "dependency_results": dependencies,
        "available_datasets": datasets,
        "test_scripts_created": [
            "test_spatial_genes_comprehensive.py",
            "test_spatial_registration_comprehensive.py", 
            "test_spatial_domains_comprehensive.py",
            "test_integration_comprehensive.py",
            "test_pathway_enrichment_comprehensive.py",
            "test_preprocessing_comprehensive.py",
            "test_spatial_analysis_comprehensive.py",
            "test_spatial_enrichment_comprehensive.py",
            "test_gene_set_loader_comprehensive.py",
            "test_differential_comprehensive.py"
        ],
        "test_coverage": {
            "spatial_genes": "1 main function (identify_spatial_genes)",
            "spatial_registration": "3 functions (register, evaluate, mcp_register)",
            "spatial_domains": "1 main function (identify_spatial_domains)", 
            "integration": "6 functions (validation, integration, alignment, trajectory, contrastive_vi, mcp_interface)",
            "pathway_enrichment": "5 functions (gsea, ora, ssgsea, enrichr, availability_check)",
            "preprocessing": "2 functions (preprocess_data, preprocess_with_resolvi)",
            "spatial_analysis": "2 functions (analyze_spatial_patterns, analyze_spatial_with_scviva)",
            "spatial_enrichment": "4 functions (enrichment_analysis, spatial_metrics, gene_correlation, availability_check)",
            "gene_set_loader": "1 main function (load_gene_sets)",
            "differential": "1 main function (differential_expression)"
        }
    }
    
    # Save summary report
    report_path = Path(__file__).parent / "simple_test_summary.json"
    with open(report_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n📄 Summary report saved to: {report_path}")
    return summary

def main():
    """Main test execution"""
    print("🚀 ChatSpatial Comprehensive Tools Testing - Simple Runner")
    print("=" * 65)
    
    try:
        summary = generate_summary_report()
        
        print("\n" + "=" * 65)
        print("🏁 SIMPLE TESTING COMPLETE")
        print("=" * 65)
        
        print(f"✅ Module imports: {summary['summary']['import_success_rate']}")
        print(f"📊 Dependencies: {summary['summary']['dependencies_available']} available") 
        print(f"💾 Datasets: {summary['summary']['datasets_available']} available")
        print(f"📝 Test scripts: {len(summary['test_scripts_created'])} created")
        
        if summary['summary']['successful_imports'] == summary['summary']['total_modules_tested']:
            print("\n🎉 ALL MODULES IMPORT SUCCESSFULLY!")
            print("✅ Ready for comprehensive testing")
        else:
            failed_modules = [name for name, result in summary['module_results'].items() 
                            if not result.startswith('✅')]
            print(f"\n⚠️  Some modules need attention: {failed_modules}")
        
    except Exception as e:
        print(f"❌ Testing failed with error: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    main()