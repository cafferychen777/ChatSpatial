"""
Performance Benchmark for ChatSpatial Core Tools

Tests memory usage, execution time, and scalability with different dataset sizes.
Focus on measurable performance metrics rather than functional correctness.
"""

import asyncio
import sys
import time
import psutil
import gc
from pathlib import Path
from typing import Dict, Any, List, Tuple
import traceback

# Import paths
script_dir = Path(__file__).parent  
sys.path.insert(0, str(script_dir.parent / "tool_functionality_tests"))

import numpy as np
import pandas as pd
import scanpy as sc
from test_base import create_synthetic_adata, MockContext


class PerformanceBenchmark:
    """Benchmark performance of core tools with different dataset sizes"""
    
    def __init__(self):
        self.results = []
        self.process = psutil.Process()
        
    def get_memory_usage(self) -> float:
        """Get current memory usage in MB"""
        return self.process.memory_info().rss / 1024 / 1024
        
    def measure_performance(self, func_name: str):
        """Decorator to measure function performance"""
        def decorator(func):
            async def wrapper(*args, **kwargs):
                # Measure before
                start_time = time.time()
                start_memory = self.get_memory_usage()
                gc.collect()  # Clean up before measurement
                
                try:
                    result = await func(*args, **kwargs)
                    success = True
                    error = None
                except Exception as e:
                    result = None
                    success = False
                    error = str(e)
                
                # Measure after
                end_time = time.time()
                end_memory = self.get_memory_usage()
                
                # Record performance data
                perf_data = {
                    "function": func_name,
                    "execution_time": end_time - start_time,
                    "memory_before": start_memory,
                    "memory_after": end_memory,
                    "memory_delta": end_memory - start_memory,
                    "success": success,
                    "error": error
                }
                
                if hasattr(self, 'current_dataset_info'):
                    perf_data.update(self.current_dataset_info)
                
                self.results.append(perf_data)
                return result
                
            return wrapper
        return decorator
        
    async def benchmark_data_loading_and_preprocessing(self):
        """Benchmark data loading and basic preprocessing"""
        dataset_sizes = [
            (50, 100, "tiny"),
            (200, 500, "small"), 
            (500, 1000, "medium"),
            (1000, 2000, "large")
        ]
        
        for n_cells, n_genes, size_label in dataset_sizes:
            self.current_dataset_info = {
                "dataset_size": size_label,
                "n_cells": n_cells,
                "n_genes": n_genes
            }
            
            print(f"Benchmarking {size_label} dataset ({n_cells} cells, {n_genes} genes)")
            
            # Test data creation
            @self.measure_performance(f"create_synthetic_data_{size_label}")
            async def create_data():
                return create_synthetic_adata(n_cells, n_genes, spatial_pattern="structured")
            
            adata = await create_data()
            if adata is None:
                continue
                
            # Test preprocessing steps
            @self.measure_performance(f"preprocessing_{size_label}")
            async def preprocess_data():
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)
                sc.pp.highly_variable_genes(adata, n_top_genes=min(1000, n_genes//2))
                return True
                
            await preprocess_data()
            
            # Test clustering
            @self.measure_performance(f"clustering_{size_label}")
            async def perform_clustering():
                sc.pp.pca(adata, n_comps=min(50, min(n_cells-1, n_genes-1)))
                sc.pp.neighbors(adata, n_neighbors=min(15, n_cells//4))
                sc.tl.leiden(adata, resolution=0.3)
                return True
                
            await perform_clustering()
            
            # Cleanup
            del adata
            gc.collect()
            
    async def benchmark_tool_imports(self):
        """Benchmark import time for different tools"""
        tools_to_import = [
            ("cell_communication", "chatspatial.tools.cell_communication", "analyze_cell_communication"),
            ("annotation", "chatspatial.tools.annotation", "annotate_cell_types"), 
            ("deconvolution", "chatspatial.tools.deconvolution", "deconvolve_spatial_data"),
            ("visualization", "chatspatial.tools.visualization", "visualize_data"),
            ("spatial_analysis", "chatspatial.tools.spatial_analysis", "analyze_spatial_patterns")
        ]
        
        for tool_name, module_path, function_name in tools_to_import:
            @self.measure_performance(f"import_{tool_name}")
            async def import_tool():
                module = __import__(module_path, fromlist=[function_name])
                func = getattr(module, function_name)
                return func is not None
                
            await import_tool()
            
    async def benchmark_parameter_creation(self):
        """Benchmark parameter object creation"""
        from chatspatial.models.data import (
            CellCommunicationParameters, 
            AnnotationParameters,
            DeconvolutionParameters,
            SpatialAnalysisParameters
        )
        
        parameter_classes = [
            ("CellCommunicationParameters", CellCommunicationParameters),
            ("AnnotationParameters", AnnotationParameters),
            ("DeconvolutionParameters", DeconvolutionParameters), 
            ("SpatialAnalysisParameters", SpatialAnalysisParameters)
        ]
        
        for param_name, param_class in parameter_classes:
            @self.measure_performance(f"create_params_{param_name}")
            async def create_params():
                # Create with defaults
                params = param_class()
                # Test basic attribute access
                method = getattr(params, 'method', None)
                return method is not None
                
            await create_params()
            
    async def benchmark_dependency_checks(self):
        """Benchmark dependency availability checks"""
        dependencies = [
            "liana", "scvi", "tangram", "cellphonedb", 
            "squidpy", "scanpy", "pandas", "numpy"
        ]
        
        for dep in dependencies:
            @self.measure_performance(f"check_dependency_{dep}")
            async def check_dep():
                try:
                    __import__(dep)
                    return True
                except ImportError:
                    return False
                    
            await check_dep()
            
    async def run_all_benchmarks(self):
        """Run complete benchmark suite"""
        print("🚀 ChatSpatial Performance Benchmark Suite")
        print("=" * 60)
        
        benchmark_suites = [
            ("Tool Imports", self.benchmark_tool_imports),
            ("Parameter Creation", self.benchmark_parameter_creation),
            ("Dependency Checks", self.benchmark_dependency_checks),
            ("Data Processing", self.benchmark_data_loading_and_preprocessing)
        ]
        
        for suite_name, suite_func in benchmark_suites:
            print(f"\n📊 Running: {suite_name}")
            try:
                await suite_func()
                print(f"✅ {suite_name} completed")
            except Exception as e:
                print(f"❌ {suite_name} failed: {str(e)}")
                
    def generate_performance_report(self) -> str:
        """Generate detailed performance report"""
        if not self.results:
            return "No performance data collected."
            
        df = pd.DataFrame(self.results)
        
        report = []
        report.append("# ChatSpatial Core Tools - Performance Benchmark Report")
        report.append(f"Generated: {pd.Timestamp.now()}")
        report.append(f"Total measurements: {len(df)}")
        report.append("")
        
        # Overall statistics
        successful_tests = df[df['success'] == True]
        report.append("## Overall Performance Summary")
        report.append(f"- **Successful tests**: {len(successful_tests)}/{len(df)} ({len(successful_tests)/len(df)*100:.1f}%)")
        if len(successful_tests) > 0:
            report.append(f"- **Average execution time**: {successful_tests['execution_time'].mean():.3f}s")
            report.append(f"- **Average memory usage**: {successful_tests['memory_after'].mean():.1f}MB")
            report.append(f"- **Average memory delta**: {successful_tests['memory_delta'].mean():.1f}MB")
        report.append("")
        
        # Performance by category
        categories = ["import", "create", "check", "preprocess", "cluster", "data"]
        for category in categories:
            cat_data = df[df['function'].str.contains(category, case=False, na=False)]
            if len(cat_data) > 0:
                cat_success = cat_data[cat_data['success'] == True]
                report.append(f"## {category.title()} Performance")
                report.append(f"- Tests: {len(cat_success)}/{len(cat_data)} successful")
                
                if len(cat_success) > 0:
                    report.append(f"- Avg time: {cat_success['execution_time'].mean():.3f}s")
                    report.append(f"- Max time: {cat_success['execution_time'].max():.3f}s")
                    report.append(f"- Avg memory: {cat_success['memory_delta'].mean():.1f}MB")
                    
                    # Show individual results
                    for _, row in cat_success.iterrows():
                        report.append(f"  - `{row['function']}`: {row['execution_time']:.3f}s, {row['memory_delta']:.1f}MB")
                        
                if len(cat_data) > len(cat_success):
                    failed = cat_data[cat_data['success'] == False]
                    report.append("  **Failures**:")
                    for _, row in failed.iterrows():
                        report.append(f"  - `{row['function']}`: {row['error']}")
                report.append("")
        
        # Dataset scalability analysis
        dataset_results = df[df['dataset_size'].notna() if 'dataset_size' in df.columns else df.iloc[:0]]
        if len(dataset_results) > 0:
            report.append("## Dataset Scalability")
            
            for size in ['tiny', 'small', 'medium', 'large']:
                size_data = dataset_results[dataset_results['dataset_size'] == size]
                if len(size_data) > 0:
                    successful_size = size_data[size_data['success'] == True]
                    if len(successful_size) > 0:
                        report.append(f"### {size.title()} Dataset")
                        sample_row = successful_size.iloc[0]
                        report.append(f"- Size: {sample_row.get('n_cells', 'unknown')} cells × {sample_row.get('n_genes', 'unknown')} genes")
                        report.append(f"- Total time: {successful_size['execution_time'].sum():.2f}s")
                        report.append(f"- Peak memory: {successful_size['memory_after'].max():.1f}MB")
                        report.append(f"- Success rate: {len(successful_size)}/{len(size_data)} operations")
                        report.append("")
        
        # Recommendations
        report.append("## Performance Recommendations")
        
        # Memory analysis
        high_memory = successful_tests[successful_tests['memory_delta'] > 100]
        if len(high_memory) > 0:
            report.append("### High Memory Usage Detected")
            for _, row in high_memory.iterrows():
                report.append(f"- `{row['function']}`: {row['memory_delta']:.1f}MB - Consider optimization")
        
        # Timing analysis  
        slow_operations = successful_tests[successful_tests['execution_time'] > 5.0]
        if len(slow_operations) > 0:
            report.append("### Slow Operations Detected") 
            for _, row in slow_operations.iterrows():
                report.append(f"- `{row['function']}`: {row['execution_time']:.2f}s - Consider caching or optimization")
        
        if len(high_memory) == 0 and len(slow_operations) == 0:
            report.append("🎉 All operations performed within acceptable performance thresholds!")
        
        return "\n".join(report)


async def main():
    """Run performance benchmark"""
    
    # Suppress warnings for cleaner output
    import warnings
    warnings.filterwarnings('ignore')
    
    benchmark = PerformanceBenchmark()
    
    try:
        await benchmark.run_all_benchmarks()
        
        # Generate report
        report = benchmark.generate_performance_report()
        
        # Save report
        report_path = Path("performance_benchmark_report.md")
        report_path.write_text(report)
        
        print("\n" + "=" * 60)
        print("📊 PERFORMANCE BENCHMARK COMPLETE")
        print(f"📄 Report saved: {report_path.absolute()}")
        
        # Summary to console
        df = pd.DataFrame(benchmark.results)
        successful_tests = len(df[df['success'] == True]) if len(df) > 0 else 0
        total_tests = len(df)
        
        if total_tests > 0:
            success_rate = successful_tests / total_tests * 100
            print(f"📈 Success Rate: {success_rate:.1f}% ({successful_tests}/{total_tests})")
            
            if successful_tests > 0:
                successful_df = df[df['success'] == True]
                avg_time = successful_df['execution_time'].mean()
                avg_memory = successful_df['memory_delta'].mean()
                print(f"⏱️  Average Time: {avg_time:.3f}s")
                print(f"💾 Average Memory: {avg_memory:.1f}MB")
        
        return 0 if successful_tests > 0 else 1
        
    except Exception as e:
        print(f"💥 Benchmark failed: {str(e)}")
        traceback.print_exc()
        return 2


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)