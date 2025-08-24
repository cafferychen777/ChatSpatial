"""
Comprehensive test runner for ChatSpatial tool functionality tests.

This script runs all tool functionality tests and generates a comprehensive report
covering correctness, performance, and compatibility across all tool modules.
"""

import pytest
import sys
import json
import time
from pathlib import Path
from typing import Dict, Any, List
import argparse

def run_test_suite(test_type: str = "all", verbose: bool = True) -> Dict[str, Any]:
    """
    Run comprehensive test suite
    
    Args:
        test_type: Type of tests to run ("all", "basic", "performance", "integration")
        verbose: Whether to show verbose output
    
    Returns:
        Test results summary
    """
    test_dir = Path(__file__).parent
    
    # Configure pytest arguments
    pytest_args = [str(test_dir)]
    
    if verbose:
        pytest_args.extend(["-v", "-s"])
    
    # Add markers based on test type
    if test_type == "basic":
        pytest_args.extend(["-m", "not slow and not integration"])
    elif test_type == "performance":
        pytest_args.extend(["-m", "slow"])
    elif test_type == "integration":
        pytest_args.extend(["-m", "integration"])
    elif test_type == "quick":
        pytest_args.extend(["-m", "not slow and not integration", "-x"])  # Stop on first failure
    
    # Add test result output
    report_file = test_dir / f"test_results_{test_type}_{int(time.time())}.json"
    pytest_args.extend(["--json-report", f"--json-report-file={report_file}"])
    
    print(f"Running ChatSpatial {test_type} test suite...")
    print(f"Test directory: {test_dir}")
    print(f"Pytest args: {' '.join(pytest_args)}")
    print("-" * 60)
    
    # Run tests
    start_time = time.time()
    exit_code = pytest.main(pytest_args)
    execution_time = time.time() - start_time
    
    # Parse results
    results = {
        "test_type": test_type,
        "exit_code": exit_code,
        "execution_time": execution_time,
        "success": exit_code == 0
    }
    
    # Load detailed results if available
    if report_file.exists():
        try:
            with open(report_file, 'r') as f:
                detailed_results = json.load(f)
            results["detailed_results"] = detailed_results
        except Exception as e:
            print(f"Warning: Could not load detailed results: {e}")
    
    return results


def generate_summary_report(results: Dict[str, Any]) -> str:
    """Generate a human-readable summary report"""
    
    report_lines = [
        "=" * 80,
        "ChatSpatial Tool Functionality Test Report",
        "=" * 80,
        f"Test Type: {results['test_type']}",
        f"Execution Time: {results['execution_time']:.2f} seconds",
        f"Success: {'✅ PASSED' if results['success'] else '❌ FAILED'}",
        "",
    ]
    
    if "detailed_results" in results:
        detailed = results["detailed_results"]
        summary = detailed.get("summary", {})
        
        report_lines.extend([
            "Test Summary:",
            f"  Total tests: {summary.get('total', 'N/A')}",
            f"  Passed: {summary.get('passed', 'N/A')}",
            f"  Failed: {summary.get('failed', 'N/A')}",
            f"  Errors: {summary.get('error', 'N/A')}",
            f"  Skipped: {summary.get('skipped', 'N/A')}",
            "",
        ])
        
        # Test breakdown by module
        if "tests" in detailed:
            tests = detailed["tests"]
            modules = {}
            
            for test in tests:
                module_name = test.get("nodeid", "").split("::")[0]
                if module_name:
                    module_name = Path(module_name).stem
                    if module_name not in modules:
                        modules[module_name] = {"passed": 0, "failed": 0, "error": 0, "skipped": 0}
                    
                    outcome = test.get("outcome", "unknown")
                    if outcome in modules[module_name]:
                        modules[module_name][outcome] += 1
            
            if modules:
                report_lines.extend([
                    "Results by Module:",
                ])
                
                for module, counts in modules.items():
                    total = sum(counts.values())
                    passed = counts.get("passed", 0)
                    success_rate = (passed / total * 100) if total > 0 else 0
                    
                    status = "✅" if counts.get("failed", 0) == 0 and counts.get("error", 0) == 0 else "❌"
                    
                    report_lines.append(
                        f"  {status} {module}: {passed}/{total} ({success_rate:.1f}% success)"
                    )
                
                report_lines.append("")
        
        # Failed tests details
        if "tests" in detailed:
            failed_tests = [test for test in detailed["tests"] 
                          if test.get("outcome") in ["failed", "error"]]
            
            if failed_tests:
                report_lines.extend([
                    "Failed Tests:",
                ])
                
                for test in failed_tests[:10]:  # Show first 10 failures
                    test_name = test.get("nodeid", "Unknown test")
                    failure_info = ""
                    
                    if "call" in test and "longrepr" in test["call"]:
                        failure_info = test["call"]["longrepr"][:200] + "..."
                    
                    report_lines.extend([
                        f"  ❌ {test_name}",
                        f"     {failure_info}",
                        ""
                    ])
                
                if len(failed_tests) > 10:
                    report_lines.append(f"  ... and {len(failed_tests) - 10} more failures")
                
                report_lines.append("")
    
    # Performance insights (if available)
    if results['test_type'] in ['all', 'performance']:
        perf_file = Path(__file__).parent / "performance_report.json"
        if perf_file.exists():
            try:
                with open(perf_file, 'r') as f:
                    perf_data = json.load(f)
                
                report_lines.extend([
                    "Performance Summary:",
                    f"  Total performance tests: {perf_data.get('summary', {}).get('total_tests', 'N/A')}",
                    f"  Total execution time: {perf_data.get('summary', {}).get('total_execution_time', 'N/A'):.2f}s",
                    f"  Total memory usage: {perf_data.get('summary', {}).get('total_memory_delta', 'N/A'):.1f}MB",
                    ""
                ])
            except Exception as e:
                report_lines.extend([
                    "Performance Summary: Could not load performance data",
                    f"Error: {e}",
                    ""
                ])
    
    # Recommendations
    report_lines.extend([
        "Recommendations:",
    ])
    
    if results['success']:
        report_lines.extend([
            "  ✅ All tests passed! The ChatSpatial tool functionality is working correctly.",
            "  📊 Consider running performance tests if not already done: python run_comprehensive_tests.py --type performance",
            "  🔄 Run tests regularly during development to catch regressions early.",
        ])
    else:
        report_lines.extend([
            "  ❌ Some tests failed. Please review the failure details above.",
            "  🔧 Fix failing tests before deploying to production.",
            "  📝 Update tests if functionality has changed intentionally.",
            "  🚀 Consider running a quick test subset during development: python run_comprehensive_tests.py --type quick",
        ])
    
    report_lines.extend([
        "",
        "=" * 80,
        "End of Report",
        "=" * 80
    ])
    
    return "\n".join(report_lines)


def main():
    """Main test runner function"""
    parser = argparse.ArgumentParser(description="Run ChatSpatial comprehensive tests")
    parser.add_argument(
        "--type", 
        choices=["all", "basic", "performance", "integration", "quick"],
        default="basic",
        help="Type of tests to run (default: basic)"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    parser.add_argument(
        "--report-only",
        action="store_true", 
        help="Only generate report from existing results"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output file for the report"
    )
    
    args = parser.parse_args()
    
    if not args.report_only:
        # Run tests
        results = run_test_suite(args.type, args.verbose)
    else:
        # Load existing results
        results = {"test_type": args.type, "success": False, "execution_time": 0}
    
    # Generate report
    report = generate_summary_report(results)
    
    # Output report
    if args.output:
        with open(args.output, 'w') as f:
            f.write(report)
        print(f"Report written to: {args.output}")
    else:
        print("\n" + report)
    
    # Exit with appropriate code
    sys.exit(0 if results['success'] else 1)


if __name__ == "__main__":
    main()