"""
Base performance testing framework for ChatSpatial.

Following Linus's philosophy:
- "Good taste is eliminating special cases"
- "Never break userspace"
- Simple, practical, no bullshit
"""
import time
import psutil
import gc
import os
import json
import threading
from typing import Dict, List, Any, Optional, Callable
from dataclasses import dataclass, asdict
from contextlib import contextmanager
from pathlib import Path
import pandas as pd
import scanpy as sc
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class PerformanceMetrics:
    """Simple data structure for performance metrics. No special cases."""
    test_name: str
    dataset_size: str  # e.g., "1000x2000"
    execution_time_seconds: float
    peak_memory_mb: float
    cpu_percent_avg: float
    memory_before_mb: float
    memory_after_mb: float
    memory_leaked_mb: float
    success: bool
    error_message: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class MemoryMonitor:
    """Memory monitoring utility. Simple and effective."""
    
    def __init__(self):
        self.process = psutil.Process(os.getpid())
        self.peak_memory = 0
        self.cpu_samples = []
        self.monitoring = False
        self.monitor_thread = None
    
    def start_monitoring(self):
        """Start memory and CPU monitoring in background thread."""
        self.monitoring = True
        self.peak_memory = 0
        self.cpu_samples.clear()
        self.monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self.monitor_thread.start()
    
    def stop_monitoring(self):
        """Stop monitoring and return stats."""
        self.monitoring = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=1.0)
        
        avg_cpu = sum(self.cpu_samples) / len(self.cpu_samples) if self.cpu_samples else 0
        return self.peak_memory, avg_cpu
    
    def _monitor_loop(self):
        """Background monitoring loop."""
        while self.monitoring:
            try:
                mem_info = self.process.memory_info()
                current_memory = mem_info.rss / 1024 / 1024  # MB
                self.peak_memory = max(self.peak_memory, current_memory)
                
                cpu_percent = self.process.cpu_percent()
                self.cpu_samples.append(cpu_percent)
                
                time.sleep(0.1)  # Sample every 100ms
            except:
                break


class BasePerformanceTest:
    """Base class for all performance tests. Eliminate special cases through uniform interface."""
    
    def __init__(self, datasets_dir: str):
        self.datasets_dir = Path(datasets_dir)
        self.results = []
        self.monitor = MemoryMonitor()
        
        # Ensure datasets directory exists
        if not self.datasets_dir.exists():
            raise FileNotFoundError(f"Datasets directory not found: {datasets_dir}")
    
    def get_available_datasets(self) -> List[Path]:
        """Get all available test datasets."""
        datasets = list(self.datasets_dir.glob("*.h5ad"))
        if not datasets:
            raise FileNotFoundError(f"No .h5ad files found in {self.datasets_dir}")
        return sorted(datasets)
    
    def load_dataset_info(self) -> pd.DataFrame:
        """Load dataset summary information."""
        summary_path = self.datasets_dir / "datasets_summary.csv"
        if summary_path.exists():
            return pd.read_csv(summary_path)
        return None
    
    @contextmanager
    def performance_context(self, test_name: str, dataset_path: Path):
        """Context manager for performance measurement. One interface for all tests."""
        # Force garbage collection before test
        gc.collect()
        
        # Get initial memory
        initial_memory = self.monitor.process.memory_info().rss / 1024 / 1024
        
        # Start monitoring
        self.monitor.start_monitoring()
        start_time = time.time()
        
        error_message = None
        success = True
        
        try:
            yield
        except Exception as e:
            error_message = str(e)
            success = False
            logger.error(f"Test {test_name} failed: {e}")
        finally:
            # Stop monitoring
            end_time = time.time()
            execution_time = end_time - start_time
            peak_memory, avg_cpu = self.monitor.stop_monitoring()
            
            # Force garbage collection and get final memory
            gc.collect()
            final_memory = self.monitor.process.memory_info().rss / 1024 / 1024
            memory_leaked = max(0, final_memory - initial_memory)
            
            # Get dataset size info
            dataset_size = self._get_dataset_size_string(dataset_path)
            
            # Create metrics
            metrics = PerformanceMetrics(
                test_name=test_name,
                dataset_size=dataset_size,
                execution_time_seconds=execution_time,
                peak_memory_mb=peak_memory,
                cpu_percent_avg=avg_cpu,
                memory_before_mb=initial_memory,
                memory_after_mb=final_memory,
                memory_leaked_mb=memory_leaked,
                success=success,
                error_message=error_message
            )
            
            self.results.append(metrics)
            logger.info(f"Test completed: {test_name} on {dataset_size} - "
                       f"Time: {execution_time:.2f}s, Peak Memory: {peak_memory:.1f}MB")
    
    def _get_dataset_size_string(self, dataset_path: Path) -> str:
        """Get dataset size as string (e.g., '1000x2000')."""
        try:
            adata = sc.read(dataset_path)
            return f"{adata.n_obs}x{adata.n_vars}"
        except:
            return dataset_path.stem
    
    def run_single_test(self, test_func: Callable, dataset_path: Path, **kwargs):
        """Run a single performance test."""
        test_name = f"{test_func.__name__}_{dataset_path.stem}"
        
        with self.performance_context(test_name, dataset_path):
            # Load data
            adata = sc.read(dataset_path)
            logger.info(f"Running {test_name} on dataset {adata.n_obs}x{adata.n_vars}")
            
            # Execute test function
            test_func(adata, **kwargs)
    
    def run_stress_test(self, test_func: Callable, iterations: int = 10, **kwargs):
        """Run stress test with multiple iterations."""
        datasets = self.get_available_datasets()
        
        for dataset_path in datasets:
            for i in range(iterations):
                test_name = f"{test_func.__name__}_stress_iter_{i}_{dataset_path.stem}"
                
                with self.performance_context(test_name, dataset_path):
                    adata = sc.read(dataset_path)
                    test_func(adata, **kwargs)
    
    def run_memory_leak_test(self, test_func: Callable, iterations: int = 50):
        """Test for memory leaks by running function multiple times."""
        datasets = [d for d in self.get_available_datasets() if "small" in d.name or "500x1k" in d.name]
        
        if not datasets:
            logger.warning("No small datasets found for memory leak testing")
            return
        
        dataset_path = datasets[0]  # Use smallest dataset
        
        for i in range(iterations):
            test_name = f"{test_func.__name__}_memleak_iter_{i}"
            
            with self.performance_context(test_name, dataset_path):
                adata = sc.read(dataset_path)
                test_func(adata)
                # Explicitly delete adata to test cleanup
                del adata
    
    def save_results(self, output_path: str):
        """Save performance results to JSON and CSV."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save as JSON for detailed analysis
        json_data = [result.to_dict() for result in self.results]
        with open(output_path.with_suffix('.json'), 'w') as f:
            json.dump(json_data, f, indent=2)
        
        # Save as CSV for easy viewing
        if self.results:
            df = pd.DataFrame([result.to_dict() for result in self.results])
            df.to_csv(output_path.with_suffix('.csv'), index=False)
        
        logger.info(f"Results saved to {output_path}")
    
    def generate_performance_report(self) -> str:
        """Generate a simple performance report. No fancy formatting, just facts."""
        if not self.results:
            return "No performance results to report."
        
        successful_tests = [r for r in self.results if r.success]
        failed_tests = [r for r in self.results if not r.success]
        
        # Basic statistics
        avg_time = sum(r.execution_time_seconds for r in successful_tests) / len(successful_tests) if successful_tests else 0
        avg_memory = sum(r.peak_memory_mb for r in successful_tests) / len(successful_tests) if successful_tests else 0
        total_memory_leaked = sum(r.memory_leaked_mb for r in self.results)
        
        report = f"""
Performance Test Report
======================

Summary:
- Total tests: {len(self.results)}
- Successful: {len(successful_tests)}
- Failed: {len(failed_tests)}
- Average execution time: {avg_time:.2f} seconds
- Average peak memory: {avg_memory:.1f} MB
- Total memory leaked: {total_memory_leaked:.1f} MB

"""
        
        if failed_tests:
            report += "\nFailed Tests:\n"
            for test in failed_tests:
                report += f"- {test.test_name}: {test.error_message}\n"
        
        # Performance warnings (Linus style: direct and practical)
        if total_memory_leaked > 100:
            report += f"\nWARNING: Memory leak detected! {total_memory_leaked:.1f} MB leaked.\n"
        
        if avg_time > 30:
            report += f"\nWARNING: Tests are slow. Average time {avg_time:.1f}s exceeds 30s threshold.\n"
        
        slowest_test = max(successful_tests, key=lambda x: x.execution_time_seconds) if successful_tests else None
        if slowest_test and slowest_test.execution_time_seconds > 60:
            report += f"\nSlowest test: {slowest_test.test_name} took {slowest_test.execution_time_seconds:.1f}s\n"
        
        return report


def create_performance_baseline(datasets_dir: str, output_dir: str):
    """Create performance baseline for future regression testing."""
    baseline_tests = BasePerformanceTest(datasets_dir)
    
    def dummy_load_test(adata):
        """Baseline test: just load and access basic properties."""
        _ = adata.n_obs
        _ = adata.n_vars
        _ = adata.X.shape
    
    # Run baseline tests on all datasets
    for dataset_path in baseline_tests.get_available_datasets():
        baseline_tests.run_single_test(dummy_load_test, dataset_path)
    
    # Save baseline
    output_path = Path(output_dir) / "performance_baseline"
    baseline_tests.save_results(str(output_path))
    
    return baseline_tests.generate_performance_report()