"""
Base Error Testing Framework

Simple, direct implementation following Linus principles:
- No over-engineering
- Clear data structures  
- Eliminate special cases through good design
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Any, Optional, Callable, Union
import traceback
import time
import sys
import os
import asyncio

# Add project root to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../../'))

from chatspatial.utils.error_handling import (
    SpatialMCPError, 
    DataNotFoundError, 
    InvalidParameterError, 
    ProcessingError,
    DataCompatibilityError
)


class ErrorCategory(Enum):
    """Error categories for systematic testing"""
    EXCEPTIONAL_INPUT = "exceptional_input"
    BOUNDARY_CONDITION = "boundary_condition" 
    RESOURCE_LIMITATION = "resource_limitation"
    DEPENDENCY_FAILURE = "dependency_failure"
    NETWORK_EXCEPTION = "network_exception"


class ErrorSeverity(Enum):
    """Error severity levels"""
    CRITICAL = "critical"  # Should never happen in production
    HIGH = "high"         # Major functionality broken
    MEDIUM = "medium"     # Degraded experience
    LOW = "low"          # Minor issues


@dataclass
class TestResult:
    """Test result with clear success/failure indication"""
    test_name: str
    category: ErrorCategory
    severity: ErrorSeverity
    passed: bool
    error_message: Optional[str] = None
    execution_time: float = 0.0
    expected_exception: Optional[type] = None
    actual_exception: Optional[type] = None
    additional_info: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.additional_info is None:
            self.additional_info = {}


@dataclass  
class ErrorTestCase:
    """
    Simple test case data structure.
    No fancy inheritance - just data and a test function.
    """
    name: str
    category: ErrorCategory
    severity: ErrorSeverity
    description: str
    test_function: Callable
    expected_exception: Optional[type] = None
    setup_data: Dict[str, Any] = None
    cleanup_function: Optional[Callable] = None
    
    def __post_init__(self):
        if self.setup_data is None:
            self.setup_data = {}


class ErrorTestRunner:
    """
    Test runner that executes error tests systematically.
    
    Design principle: Keep it simple. No complex state machines or 
    fancy abstractions. Just run tests and collect results.
    """
    
    def __init__(self):
        self.results: List[TestResult] = []
        self.setup_errors: List[str] = []
    
    async def run_test_case(self, test_case: ErrorTestCase) -> TestResult:
        """Run a single test case and return result"""
        start_time = time.time()
        
        try:
            # Setup if needed
            if test_case.setup_data:
                # Test functions should handle their own setup
                pass
                
            # Execute test
            if asyncio.iscoroutinefunction(test_case.test_function):
                await test_case.test_function(**test_case.setup_data)
            else:
                test_case.test_function(**test_case.setup_data)
            
            # If we get here and expected an exception, test failed
            if test_case.expected_exception:
                result = TestResult(
                    test_name=test_case.name,
                    category=test_case.category,
                    severity=test_case.severity,
                    passed=False,
                    error_message=f"Expected {test_case.expected_exception.__name__} but none was raised",
                    execution_time=time.time() - start_time
                )
            else:
                # Test passed normally
                result = TestResult(
                    test_name=test_case.name,
                    category=test_case.category, 
                    severity=test_case.severity,
                    passed=True,
                    execution_time=time.time() - start_time
                )
                
        except Exception as e:
            # Check if this was the expected exception
            if test_case.expected_exception and isinstance(e, test_case.expected_exception):
                result = TestResult(
                    test_name=test_case.name,
                    category=test_case.category,
                    severity=test_case.severity, 
                    passed=True,
                    execution_time=time.time() - start_time,
                    expected_exception=test_case.expected_exception,
                    actual_exception=type(e)
                )
            else:
                # Unexpected exception
                result = TestResult(
                    test_name=test_case.name,
                    category=test_case.category,
                    severity=test_case.severity,
                    passed=False,
                    error_message=str(e),
                    execution_time=time.time() - start_time,
                    expected_exception=test_case.expected_exception,
                    actual_exception=type(e),
                    additional_info={"traceback": traceback.format_exc()}
                )
        
        # Cleanup if needed
        if test_case.cleanup_function:
            try:
                if asyncio.iscoroutinefunction(test_case.cleanup_function):
                    await test_case.cleanup_function()
                else:
                    test_case.cleanup_function()
            except Exception as cleanup_error:
                result.additional_info["cleanup_error"] = str(cleanup_error)
        
        return result
    
    async def run_test_suite(self, test_cases: List[ErrorTestCase]) -> List[TestResult]:
        """Run a suite of test cases"""
        results = []
        
        for test_case in test_cases:
            print(f"Running: {test_case.name}")
            result = await self.run_test_case(test_case)
            results.append(result)
            
            # Print immediate feedback
            status = "✓" if result.passed else "✗"
            print(f"  {status} {result.test_name} ({result.execution_time:.2f}s)")
            if not result.passed:
                print(f"    Error: {result.error_message}")
        
        self.results.extend(results)
        return results
    
    def generate_summary_report(self) -> Dict[str, Any]:
        """Generate a summary report of all test results"""
        if not self.results:
            return {"error": "No test results available"}
        
        total_tests = len(self.results)
        passed_tests = sum(1 for r in self.results if r.passed)
        failed_tests = total_tests - passed_tests
        
        # Group by category
        by_category = {}
        for result in self.results:
            cat = result.category.value
            if cat not in by_category:
                by_category[cat] = {"total": 0, "passed": 0, "failed": 0}
            by_category[cat]["total"] += 1
            if result.passed:
                by_category[cat]["passed"] += 1
            else:
                by_category[cat]["failed"] += 1
        
        # Group by severity
        by_severity = {}
        for result in self.results:
            sev = result.severity.value
            if sev not in by_severity:
                by_severity[sev] = {"total": 0, "passed": 0, "failed": 0}
            by_severity[sev]["total"] += 1
            if result.passed:
                by_severity[sev]["passed"] += 1
            else:
                by_severity[sev]["failed"] += 1
        
        # Failed tests details
        failed_details = []
        for result in self.results:
            if not result.passed:
                failed_details.append({
                    "name": result.test_name,
                    "category": result.category.value,
                    "severity": result.severity.value,
                    "error": result.error_message,
                    "expected_exception": result.expected_exception.__name__ if result.expected_exception else None,
                    "actual_exception": result.actual_exception.__name__ if result.actual_exception else None
                })
        
        return {
            "summary": {
                "total_tests": total_tests,
                "passed": passed_tests,
                "failed": failed_tests,
                "success_rate": f"{(passed_tests/total_tests)*100:.1f}%"
            },
            "by_category": by_category,
            "by_severity": by_severity, 
            "failed_tests": failed_details,
            "setup_errors": self.setup_errors
        }


def create_test_with_data(test_data_path: str = None) -> Dict[str, Any]:
    """
    Load test data for error handling tests.
    Simple approach: return paths to existing test datasets.
    """
    if test_data_path is None:
        # Use default test data directory
        test_data_path = os.path.join(
            os.path.dirname(__file__), 
            "../datasets"
        )
    
    test_datasets = {
        "empty_dataset": os.path.join(test_data_path, "empty_dataset.h5ad"),
        "single_cell": os.path.join(test_data_path, "single_cell.h5ad"), 
        "high_sparsity": os.path.join(test_data_path, "high_sparsity.h5ad"),
        "no_spatial": os.path.join(test_data_path, "no_spatial.h5ad"),
        "small_synthetic": os.path.join(test_data_path, "small_synthetic.h5ad"),
        "large_synthetic": os.path.join(test_data_path, "large_synthetic.h5ad"),
    }
    
    return test_datasets