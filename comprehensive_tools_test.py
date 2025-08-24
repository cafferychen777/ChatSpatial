#!/usr/bin/env python3
"""
Comprehensive testing suite for ChatSpatial tools.
Tests all major functionality with real datasets.

Based on Linus principles:
1. Good taste - eliminate special cases
2. Practical testing of real problems
3. Simple data structures, complex behavior

Author: Generated with Claude Code
"""

import os
import sys
import time
import traceback
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

import numpy as np
import pandas as pd
import scanpy as sc
import anndata

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# Add tools to path
sys.path.append('/Users/apple/Research/SpatialTrans_MCP/chatspatial')

class ToolsTestSuite:
    """
    Comprehensive test suite for all ChatSpatial tools.
    
    Follows Linus principle: "Bad programmers worry about the code. 
    Good programmers worry about data structures."
    """
    
    def __init__(self):
        self.data_dir = Path('/Users/apple/Research/SpatialTrans_MCP/chatspatial/data')
        self.results = {}
        self.test_datasets = self._get_test_datasets()
        
    def _get_test_datasets(self) -> Dict[str, Path]:
        """Get carefully selected test datasets for different scenarios."""
        return {
            # Small demo datasets for fast testing
            'small_demo': self.data_dir / 'demo_datasets/quick_demo_jurkat.h5ad',
            'medium_demo': self.data_dir / 'demo_datasets/quick_demo_combined.h5ad',
            
            # Spatial datasets with coordinates
            'spatial_small': self.data_dir / 'spatial_datasets/seqfish_developmental.h5ad',
            'spatial_medium': self.data_dir / 'spatial_datasets/slideseq_cerebellum.h5ad',
            
            # Benchmark datasets for performance
            'benchmark_small': self.data_dir / 'benchmark_datasets/benchmark_500x1k.h5ad',
            'benchmark_medium': self.data_dir / 'benchmark_datasets/benchmark_2kx3k.h5ad',
            
            # Real datasets
            'real_data': self.data_dir / 'real_datasets/mop_sn_tutorial.h5ad'
        }
    
    def _load_dataset(self, name: str) -> Optional[anndata.AnnData]:
        """Load dataset with error handling."""
        try:
            path = self.test_datasets[name]
            if not path.exists():
                print(f"❌ Dataset not found: {path}")
                return None
            
            adata = sc.read_h5ad(path)
            print(f"✅ Loaded {name}: {adata.n_obs} cells, {adata.n_vars} genes")
            return adata
            
        except Exception as e:
            print(f"❌ Failed to load {name}: {e}")
            return None
    
    def _time_function(self, func, *args, **kwargs) -> Tuple[Any, float]:
        """Time function execution."""
        start = time.time()
        try:
            result = func(*args, **kwargs)
            elapsed = time.time() - start
            return result, elapsed
        except Exception as e:
            elapsed = time.time() - start
            raise Exception(f"Function failed after {elapsed:.2f}s: {e}")
    
    def test_preprocessing(self) -> Dict[str, Any]:
        """Test preprocessing functionality."""
        print("\n🧪 Testing Preprocessing Tools...")
        
        try:
            from chatspatial.tools.preprocessing import (
                quality_control, normalize_and_scale, 
                filter_genes_and_cells, calculate_qc_metrics
            )
            
            results = {}
            
            # Test with multiple datasets
            for name in ['small_demo', 'medium_demo']:
                adata = self._load_dataset(name)
                if adata is None:
                    continue
                
                print(f"  Testing preprocessing on {name}...")
                
                # Quality control
                try:
                    qc_result, qc_time = self._time_function(quality_control, adata.copy())
                    results[f'{name}_qc'] = {'success': True, 'time': qc_time}
                    print(f"    ✅ QC completed in {qc_time:.2f}s")
                except Exception as e:
                    results[f'{name}_qc'] = {'success': False, 'error': str(e)}
                    print(f"    ❌ QC failed: {e}")
                
                # Normalization
                try:
                    norm_result, norm_time = self._time_function(normalize_and_scale, adata.copy())
                    results[f'{name}_norm'] = {'success': True, 'time': norm_time}
                    print(f"    ✅ Normalization completed in {norm_time:.2f}s")
                except Exception as e:
                    results[f'{name}_norm'] = {'success': False, 'error': str(e)}
                    print(f"    ❌ Normalization failed: {e}")
            
            return results
            
        except ImportError as e:
            return {'error': f'Import failed: {e}'}
    
    def test_spatial_analysis(self) -> Dict[str, Any]:
        """Test spatial analysis tools."""
        print("\n🗺️  Testing Spatial Analysis Tools...")
        
        try:
            from chatspatial.tools.spatial_analysis import (
                calculate_spatial_autocorrelation, 
                find_spatial_patterns,
                spatial_neighbor_analysis
            )
            
            results = {}
            
            # Test with spatial datasets
            for name in ['spatial_small', 'spatial_medium']:
                adata = self._load_dataset(name)
                if adata is None:
                    continue
                
                # Check if spatial coordinates exist
                if 'spatial' not in adata.obsm and 'X_spatial' not in adata.obsm:
                    print(f"    ⚠️  No spatial coordinates in {name}, skipping...")
                    continue
                
                print(f"  Testing spatial analysis on {name}...")
                
                # Spatial autocorrelation
                try:
                    auto_result, auto_time = self._time_function(
                        calculate_spatial_autocorrelation, adata.copy()
                    )
                    results[f'{name}_autocorr'] = {'success': True, 'time': auto_time}
                    print(f"    ✅ Spatial autocorrelation completed in {auto_time:.2f}s")
                except Exception as e:
                    results[f'{name}_autocorr'] = {'success': False, 'error': str(e)}
                    print(f"    ❌ Spatial autocorrelation failed: {e}")
            
            return results
            
        except ImportError as e:
            return {'error': f'Import failed: {e}'}
    
    def test_spatial_domains(self) -> Dict[str, Any]:
        """Test spatial domain detection."""
        print("\n🎯 Testing Spatial Domain Tools...")
        
        try:
            from chatspatial.tools.spatial_domains import (
                detect_spatial_domains_spagcn,
                detect_spatial_domains_scanpy,
                spatial_clustering
            )
            
            results = {}
            
            # Test with spatial datasets
            for name in ['spatial_small']:  # Start with small dataset
                adata = self._load_dataset(name)
                if adata is None:
                    continue
                
                print(f"  Testing spatial domains on {name}...")
                
                # SpaGCN clustering
                try:
                    spagcn_result, spagcn_time = self._time_function(
                        detect_spatial_domains_spagcn, adata.copy()
                    )
                    results[f'{name}_spagcn'] = {'success': True, 'time': spagcn_time}
                    print(f"    ✅ SpaGCN completed in {spagcn_time:.2f}s")
                except Exception as e:
                    results[f'{name}_spagcn'] = {'success': False, 'error': str(e)}
                    print(f"    ❌ SpaGCN failed: {e}")
                
                # Scanpy clustering
                try:
                    scanpy_result, scanpy_time = self._time_function(
                        detect_spatial_domains_scanpy, adata.copy()
                    )
                    results[f'{name}_scanpy_domains'] = {'success': True, 'time': scanpy_time}
                    print(f"    ✅ Scanpy domains completed in {scanpy_time:.2f}s")
                except Exception as e:
                    results[f'{name}_scanpy_domains'] = {'success': False, 'error': str(e)}
                    print(f"    ❌ Scanpy domains failed: {e}")
            
            return results
            
        except ImportError as e:
            return {'error': f'Import failed: {e}'}
    
    def test_visualization(self) -> Dict[str, Any]:
        """Test visualization tools."""
        print("\n📊 Testing Visualization Tools...")
        
        try:
            from chatspatial.tools.visualization import (
                plot_spatial_expression,
                plot_spatial_domains,
                create_spatial_plot
            )
            
            results = {}
            
            # Test with spatial datasets
            for name in ['small_demo']:
                adata = self._load_dataset(name)
                if adata is None:
                    continue
                
                print(f"  Testing visualization on {name}...")
                
                # Basic spatial plot
                try:
                    plot_result, plot_time = self._time_function(
                        create_spatial_plot, adata.copy()
                    )
                    results[f'{name}_plot'] = {'success': True, 'time': plot_time}
                    print(f"    ✅ Spatial plot completed in {plot_time:.2f}s")
                except Exception as e:
                    results[f'{name}_plot'] = {'success': False, 'error': str(e)}
                    print(f"    ❌ Spatial plot failed: {e}")
            
            return results
            
        except ImportError as e:
            return {'error': f'Import failed: {e}'}
    
    def test_trajectory_analysis(self) -> Dict[str, Any]:
        """Test trajectory analysis tools."""
        print("\n🛤️  Testing Trajectory Analysis Tools...")
        
        try:
            from chatspatial.tools.trajectory import (
                trajectory_analysis,
                pseudotime_calculation,
                trajectory_inference
            )
            
            results = {}
            
            # Test with suitable datasets
            for name in ['medium_demo']:
                adata = self._load_dataset(name)
                if adata is None:
                    continue
                
                print(f"  Testing trajectory analysis on {name}...")
                
                # Basic trajectory analysis
                try:
                    traj_result, traj_time = self._time_function(
                        trajectory_analysis, adata.copy()
                    )
                    results[f'{name}_trajectory'] = {'success': True, 'time': traj_time}
                    print(f"    ✅ Trajectory analysis completed in {traj_time:.2f}s")
                except Exception as e:
                    results[f'{name}_trajectory'] = {'success': False, 'error': str(e)}
                    print(f"    ❌ Trajectory analysis failed: {e}")
            
            return results
            
        except ImportError as e:
            return {'error': f'Import failed: {e}'}
    
    def test_enrichment_analysis(self) -> Dict[str, Any]:
        """Test enrichment analysis tools."""
        print("\n🧬 Testing Enrichment Analysis Tools...")
        
        try:
            from chatspatial.tools.pathway_enrichment import (
                perform_enrichment_analysis,
                gene_set_enrichment,
                pathway_analysis
            )
            
            results = {}
            
            # Test with gene expression datasets
            for name in ['small_demo']:
                adata = self._load_dataset(name)
                if adata is None:
                    continue
                
                print(f"  Testing enrichment analysis on {name}...")
                
                # Basic enrichment analysis
                try:
                    enrich_result, enrich_time = self._time_function(
                        perform_enrichment_analysis, adata.copy()
                    )
                    results[f'{name}_enrichment'] = {'success': True, 'time': enrich_time}
                    print(f"    ✅ Enrichment analysis completed in {enrich_time:.2f}s")
                except Exception as e:
                    results[f'{name}_enrichment'] = {'success': False, 'error': str(e)}
                    print(f"    ❌ Enrichment analysis failed: {e}")
            
            return results
            
        except ImportError as e:
            return {'error': f'Import failed: {e}'}
    
    def run_comprehensive_test(self) -> Dict[str, Any]:
        """Run all tests and generate comprehensive report."""
        print("🚀 Starting Comprehensive ChatSpatial Tools Testing")
        print("=" * 60)
        
        start_time = time.time()
        
        # Run all test modules
        test_modules = [
            ('preprocessing', self.test_preprocessing),
            ('spatial_analysis', self.test_spatial_analysis), 
            ('spatial_domains', self.test_spatial_domains),
            ('visualization', self.test_visualization),
            ('trajectory', self.test_trajectory_analysis),
            ('enrichment', self.test_enrichment_analysis)
        ]
        
        all_results = {}
        
        for module_name, test_func in test_modules:
            try:
                print(f"\n📍 Testing {module_name}...")
                module_results = test_func()
                all_results[module_name] = module_results
                
                # Count successes
                if isinstance(module_results, dict) and 'error' not in module_results:
                    successes = sum(1 for r in module_results.values() 
                                  if isinstance(r, dict) and r.get('success', False))
                    total = len([r for r in module_results.values() 
                               if isinstance(r, dict)])
                    print(f"    📊 {module_name}: {successes}/{total} tests passed")
                
            except Exception as e:
                print(f"    ❌ {module_name} testing failed: {e}")
                all_results[module_name] = {'error': str(e)}
        
        total_time = time.time() - start_time
        
        # Generate summary report
        self._generate_test_report(all_results, total_time)
        
        return all_results
    
    def _generate_test_report(self, results: Dict[str, Any], total_time: float):
        """Generate comprehensive test report."""
        print("\n" + "=" * 60)
        print("📋 COMPREHENSIVE TEST REPORT")
        print("=" * 60)
        
        total_tests = 0
        passed_tests = 0
        
        for module_name, module_results in results.items():
            print(f"\n🔧 {module_name.upper()}")
            print("-" * 30)
            
            if isinstance(module_results, dict) and 'error' not in module_results:
                for test_name, test_result in module_results.items():
                    if isinstance(test_result, dict):
                        total_tests += 1
                        if test_result.get('success', False):
                            passed_tests += 1
                            time_str = f"{test_result['time']:.2f}s" if 'time' in test_result else 'N/A'
                            print(f"  ✅ {test_name}: {time_str}")
                        else:
                            error_msg = test_result.get('error', 'Unknown error')[:50]
                            print(f"  ❌ {test_name}: {error_msg}")
            else:
                error_msg = module_results.get('error', 'Module failed')[:50]
                print(f"  ❌ Module error: {error_msg}")
        
        # Summary statistics
        print(f"\n📊 SUMMARY")
        print("-" * 30)
        print(f"Total tests run: {total_tests}")
        print(f"Tests passed: {passed_tests}")
        print(f"Tests failed: {total_tests - passed_tests}")
        print(f"Success rate: {passed_tests/total_tests*100:.1f}%" if total_tests > 0 else "No tests run")
        print(f"Total time: {total_time:.2f}s")
        
        # Save detailed report to file
        report_path = 'tools_test_report.md'
        self._save_detailed_report(results, total_time, report_path)
        print(f"\n📄 Detailed report saved to: {report_path}")
    
    def _save_detailed_report(self, results: Dict[str, Any], total_time: float, filepath: str):
        """Save detailed markdown report."""
        with open(filepath, 'w') as f:
            f.write("# ChatSpatial Tools Comprehensive Test Report\n\n")
            f.write(f"**Date**: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Total Duration**: {total_time:.2f} seconds\n\n")
            
            # Test datasets used
            f.write("## Test Datasets\n\n")
            for name, path in self.test_datasets.items():
                exists = "✅" if path.exists() else "❌"
                f.write(f"- **{name}**: {exists} `{path.name}`\n")
            
            # Detailed results
            f.write("\n## Test Results by Module\n\n")
            
            for module_name, module_results in results.items():
                f.write(f"### {module_name.title()}\n\n")
                
                if isinstance(module_results, dict) and 'error' not in module_results:
                    for test_name, test_result in module_results.items():
                        if isinstance(test_result, dict):
                            if test_result.get('success', False):
                                time_info = f" ({test_result['time']:.2f}s)" if 'time' in test_result else ""
                                f.write(f"- ✅ **{test_name}**{time_info}\n")
                            else:
                                error_msg = test_result.get('error', 'Unknown error')
                                f.write(f"- ❌ **{test_name}**: {error_msg}\n")
                else:
                    error_msg = module_results.get('error', 'Module failed')
                    f.write(f"- ❌ Module import/setup failed: {error_msg}\n")
                
                f.write("\n")
            
            f.write("---\n")
            f.write("*Generated by ChatSpatial Comprehensive Testing Suite*\n")


if __name__ == "__main__":
    # Linus principle: "Talk is cheap. Show me the code."
    tester = ToolsTestSuite()
    results = tester.run_comprehensive_test()