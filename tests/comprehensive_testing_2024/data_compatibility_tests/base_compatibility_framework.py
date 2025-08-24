"""
Base Compatibility Testing Framework

Simple, systematic approach to compatibility testing.
No over-engineering - just test what matters in production.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Any, Optional, Callable, Tuple
import platform
import sys
import os
import time
import asyncio
import importlib.util

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../../'))


class CompatibilityCategory(Enum):
    """Compatibility test categories"""
    DATA_FORMAT = "data_format" 
    PARAMETER_VERSION = "parameter_version"
    PLATFORM = "platform"
    DEPENDENCY_VERSION = "dependency_version"


@dataclass
class CompatibilityResult:
    """Result of a compatibility test"""
    test_name: str
    category: CompatibilityCategory
    compatible: bool
    version_tested: str = ""
    platform_tested: str = ""
    error_message: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    execution_time: float = 0.0
    additional_info: Dict[str, Any] = field(default_factory=dict)


@dataclass
class CompatibilityTestCase:
    """
    A compatibility test case.
    Simple data structure - no complex inheritance.
    """
    name: str
    category: CompatibilityCategory
    description: str
    test_function: Callable
    version_info: Dict[str, str] = field(default_factory=dict)
    platform_requirements: List[str] = field(default_factory=list)
    setup_function: Optional[Callable] = None
    cleanup_function: Optional[Callable] = None
    
    
class CompatibilityMatrix:
    """
    Compatibility matrix for tracking test results across versions/platforms.
    
    Simple structure: matrix[category][version][platform] = result
    """
    
    def __init__(self):
        self.matrix: Dict[str, Dict[str, Dict[str, CompatibilityResult]]] = {}
        self.metadata = {
            "python_version": sys.version,
            "platform": platform.platform(),
            "architecture": platform.machine(),
            "test_timestamp": None
        }
    
    def add_result(self, result: CompatibilityResult):
        """Add a test result to the matrix"""
        category = result.category.value
        version = result.version_tested or "default"
        platform_info = result.platform_tested or self.metadata["platform"]
        
        if category not in self.matrix:
            self.matrix[category] = {}
        if version not in self.matrix[category]:
            self.matrix[category][version] = {}
        
        self.matrix[category][version][platform_info] = result
    
    def get_compatibility_summary(self) -> Dict[str, Any]:
        """Generate compatibility summary"""
        summary = {
            "total_tests": 0,
            "compatible": 0, 
            "incompatible": 0,
            "by_category": {},
            "critical_issues": [],
            "warnings_count": 0
        }
        
        for category, versions in self.matrix.items():
            summary["by_category"][category] = {
                "total": 0,
                "compatible": 0,
                "incompatible": 0
            }
            
            for version, platforms in versions.items():
                for platform_info, result in platforms.items():
                    summary["total_tests"] += 1
                    summary["by_category"][category]["total"] += 1
                    
                    if result.compatible:
                        summary["compatible"] += 1
                        summary["by_category"][category]["compatible"] += 1
                    else:
                        summary["incompatible"] += 1
                        summary["by_category"][category]["incompatible"] += 1
                        
                        # Track critical issues
                        summary["critical_issues"].append({
                            "test": result.test_name,
                            "category": category,
                            "version": version,
                            "platform": platform_info,
                            "error": result.error_message
                        })
                    
                    summary["warnings_count"] += len(result.warnings)
        
        return summary


class CompatibilityTestRunner:
    """
    Test runner for compatibility tests.
    
    Keep it simple: run tests, collect results, generate reports.
    """
    
    def __init__(self):
        self.matrix = CompatibilityMatrix()
        self.current_platform = platform.platform()
        self.current_python = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    
    def get_dependency_version(self, module_name: str) -> Optional[str]:
        """Get version of a dependency if available"""
        try:
            module = __import__(module_name)
            return getattr(module, '__version__', 'unknown')
        except ImportError:
            return None
    
    def check_platform_requirements(self, requirements: List[str]) -> Tuple[bool, List[str]]:
        """Check if current platform meets requirements"""
        issues = []
        current = platform.system().lower()
        
        for req in requirements:
            req_lower = req.lower()
            if req_lower not in ['any', 'all']:
                if req_lower not in current and current not in req_lower:
                    issues.append(f"Platform requirement '{req}' not met (current: {current})")
        
        return len(issues) == 0, issues
    
    async def run_test_case(self, test_case: CompatibilityTestCase) -> CompatibilityResult:
        """Run a single compatibility test case"""
        start_time = time.time()
        
        # Check platform requirements first
        platform_ok, platform_issues = self.check_platform_requirements(test_case.platform_requirements)
        if not platform_ok:
            return CompatibilityResult(
                test_name=test_case.name,
                category=test_case.category,
                compatible=False,
                platform_tested=self.current_platform,
                error_message="; ".join(platform_issues),
                execution_time=time.time() - start_time
            )
        
        try:
            # Setup if needed
            if test_case.setup_function:
                if asyncio.iscoroutinefunction(test_case.setup_function):
                    await test_case.setup_function()
                else:
                    test_case.setup_function()
            
            # Run test
            test_result = None
            warnings = []
            
            if asyncio.iscoroutinefunction(test_case.test_function):
                test_result = await test_case.test_function()
            else:
                test_result = test_case.test_function()
            
            # Handle test result
            if isinstance(test_result, dict):
                compatible = test_result.get('compatible', True)
                warnings = test_result.get('warnings', [])
                additional_info = test_result.get('additional_info', {})
                error_msg = test_result.get('error_message')
            else:
                # Assume boolean result
                compatible = bool(test_result)
                additional_info = {}
                error_msg = None
            
            result = CompatibilityResult(
                test_name=test_case.name,
                category=test_case.category,
                compatible=compatible,
                platform_tested=self.current_platform,
                version_tested=test_case.version_info.get('version', self.current_python),
                warnings=warnings,
                error_message=error_msg,
                execution_time=time.time() - start_time,
                additional_info=additional_info
            )
            
        except Exception as e:
            result = CompatibilityResult(
                test_name=test_case.name,
                category=test_case.category,
                compatible=False,
                platform_tested=self.current_platform,
                version_tested=test_case.version_info.get('version', self.current_python),
                error_message=str(e),
                execution_time=time.time() - start_time,
                additional_info={"exception_type": type(e).__name__}
            )
        
        # Cleanup if needed
        if test_case.cleanup_function:
            try:
                if asyncio.iscoroutinefunction(test_case.cleanup_function):
                    await test_case.cleanup_function()
                else:
                    test_case.cleanup_function()
            except Exception as cleanup_error:
                result.warnings.append(f"Cleanup warning: {cleanup_error}")
        
        return result
    
    async def run_test_suite(self, test_cases: List[CompatibilityTestCase]) -> List[CompatibilityResult]:
        """Run a suite of compatibility test cases"""
        results = []
        
        print(f"Running compatibility tests on {self.current_platform} with Python {self.current_python}")
        
        for test_case in test_cases:
            print(f"Testing: {test_case.name}")
            result = await self.run_test_case(test_case)
            results.append(result)
            self.matrix.add_result(result)
            
            # Print immediate feedback
            status = "✓" if result.compatible else "✗"
            print(f"  {status} {result.test_name} ({result.execution_time:.2f}s)")
            if not result.compatible:
                print(f"    Issue: {result.error_message}")
            if result.warnings:
                for warning in result.warnings:
                    print(f"    Warning: {warning}")
        
        return results
    
    def generate_compatibility_report(self) -> Dict[str, Any]:
        """Generate comprehensive compatibility report"""
        summary = self.matrix.get_compatibility_summary()
        
        # Add environment information
        environment = {
            "python_version": self.current_python,
            "platform": self.current_platform,
            "architecture": platform.machine(),
            "dependencies": {}
        }
        
        # Check key dependencies
        key_deps = ['scanpy', 'pandas', 'numpy', 'scipy', 'anndata']
        for dep in key_deps:
            version = self.get_dependency_version(dep)
            environment["dependencies"][dep] = version or "not_available"
        
        return {
            "summary": summary,
            "environment": environment,
            "detailed_results": self.matrix.matrix,
            "recommendations": self._generate_recommendations(summary)
        }
    
    def _generate_recommendations(self, summary: Dict[str, Any]) -> List[str]:
        """Generate actionable recommendations based on test results"""
        recommendations = []
        
        if summary["incompatible"] > 0:
            recommendations.append(
                f"Found {summary['incompatible']} compatibility issues that need attention"
            )
        
        critical_count = len(summary.get("critical_issues", []))
        if critical_count > 0:
            recommendations.append(
                f"Address {critical_count} critical compatibility issues before deployment"
            )
        
        if summary.get("warnings_count", 0) > 0:
            recommendations.append(
                f"Review {summary['warnings_count']} compatibility warnings for potential issues"
            )
        
        # Category-specific recommendations
        for category, stats in summary.get("by_category", {}).items():
            if stats["incompatible"] > 0:
                recommendations.append(
                    f"Investigate {category} compatibility issues: {stats['incompatible']} failures"
                )
        
        if not recommendations:
            recommendations.append("All compatibility tests passed successfully")
        
        return recommendations


def get_test_data_info() -> Dict[str, Any]:
    """Get information about available test data"""
    test_data_dir = os.path.join(os.path.dirname(__file__), "../datasets")
    
    data_info = {
        "test_data_directory": test_data_dir,
        "available_datasets": [],
        "missing_datasets": []
    }
    
    expected_datasets = [
        "empty_dataset.h5ad",
        "single_cell.h5ad", 
        "high_sparsity.h5ad",
        "no_spatial.h5ad",
        "small_synthetic.h5ad",
        "large_synthetic.h5ad"
    ]
    
    for dataset in expected_datasets:
        path = os.path.join(test_data_dir, dataset)
        if os.path.exists(path):
            data_info["available_datasets"].append({
                "name": dataset,
                "path": path,
                "size_mb": os.path.getsize(path) / (1024*1024)
            })
        else:
            data_info["missing_datasets"].append(dataset)
    
    return data_info