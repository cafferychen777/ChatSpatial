#!/usr/bin/env python3
"""
ChatSpatial MCP Server - Test Suite Runner
Comprehensive test execution with reporting and CI integration
"""

import os
import sys
import subprocess
import argparse
import json
import time
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict


@dataclass
class TestResult:
    """Test result data structure"""
    suite: str
    passed: int
    failed: int
    skipped: int
    duration: float
    coverage: Optional[float] = None
    details: Optional[Dict] = None


class TestRunner:
    """Main test runner class"""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.test_dir = project_root / "tests"
        self.results: List[TestResult] = []
    
    def run_unit_tests(self, verbose: bool = False) -> TestResult:
        """Run unit tests"""
        print("🧪 Running Unit Tests...")
        
        cmd = [
            sys.executable, "-m", "pytest",
            str(self.test_dir / "unit"),
            "-v" if verbose else "-q",
            "--tb=short",
            "--cov=chatspatial",
            "--cov-report=term-missing",
            "--cov-report=json",
            "-m", "unit"
        ]
        
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.project_root)
        duration = time.time() - start_time
        
        # Parse pytest output
        passed, failed, skipped = self._parse_pytest_output(result.stdout)
        
        # Get coverage
        coverage = self._get_coverage_percentage()
        
        test_result = TestResult(
            suite="unit",
            passed=passed,
            failed=failed,
            skipped=skipped,
            duration=duration,
            coverage=coverage,
            details={"stdout": result.stdout, "stderr": result.stderr, "returncode": result.returncode}
        )
        
        self.results.append(test_result)
        self._print_test_summary("Unit Tests", test_result)
        
        return test_result
    
    def run_tool_tests(self, verbose: bool = False) -> TestResult:
        """Run tool-level tests"""
        print("🔧 Running Tool-Level Tests...")
        
        cmd = [
            sys.executable, "-m", "pytest",
            str(self.test_dir / "tools"),
            "-v" if verbose else "-q",
            "--tb=short",
            "-m", "tool"
        ]
        
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.project_root)
        duration = time.time() - start_time
        
        passed, failed, skipped = self._parse_pytest_output(result.stdout)
        
        test_result = TestResult(
            suite="tool",
            passed=passed,
            failed=failed,
            skipped=skipped,
            duration=duration,
            details={"stdout": result.stdout, "stderr": result.stderr, "returncode": result.returncode}
        )
        
        self.results.append(test_result)
        self._print_test_summary("Tool-Level Tests", test_result)
        
        return test_result
    
    def run_workflow_tests(self, verbose: bool = False) -> TestResult:
        """Run workflow tests"""
        print("🔄 Running Workflow Tests...")
        
        cmd = [
            sys.executable, "-m", "pytest",
            str(self.test_dir / "workflows"),
            "-v" if verbose else "-q",
            "--tb=short",
            "-m", "workflow"
        ]
        
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.project_root)
        duration = time.time() - start_time
        
        passed, failed, skipped = self._parse_pytest_output(result.stdout)
        
        test_result = TestResult(
            suite="workflow",
            passed=passed,
            failed=failed,
            skipped=skipped,
            duration=duration,
            details={"stdout": result.stdout, "stderr": result.stderr, "returncode": result.returncode}
        )
        
        self.results.append(test_result)
        self._print_test_summary("Workflow Tests", test_result)
        
        return test_result
    
    def run_e2e_tests(self, verbose: bool = False) -> TestResult:
        """Run end-to-end tests"""
        print("🌐 Running End-to-End Tests...")
        
        cmd = [
            sys.executable, "-m", "pytest",
            str(self.test_dir / "e2e"),
            "-v" if verbose else "-q",
            "--tb=short",
            "-m", "e2e"
        ]
        
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.project_root)
        duration = time.time() - start_time
        
        passed, failed, skipped = self._parse_pytest_output(result.stdout)
        
        test_result = TestResult(
            suite="e2e",
            passed=passed,
            failed=failed,
            skipped=skipped,
            duration=duration,
            details={"stdout": result.stdout, "stderr": result.stderr, "returncode": result.returncode}
        )
        
        self.results.append(test_result)
        self._print_test_summary("End-to-End Tests", test_result)
        
        return test_result
    
    def run_performance_tests(self, verbose: bool = False) -> TestResult:
        """Run performance tests"""
        print("⚡ Running Performance Tests...")
        
        cmd = [
            sys.executable, "-m", "pytest",
            str(self.test_dir),
            "-v" if verbose else "-q",
            "--tb=short",
            "-m", "slow"
        ]
        
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.project_root)
        duration = time.time() - start_time
        
        passed, failed, skipped = self._parse_pytest_output(result.stdout)
        
        test_result = TestResult(
            suite="performance",
            passed=passed,
            failed=failed,
            skipped=skipped,
            duration=duration,
            details={"stdout": result.stdout, "stderr": result.stderr, "returncode": result.returncode}
        )
        
        self.results.append(test_result)
        self._print_test_summary("Performance Tests", test_result)
        
        return test_result
    
    def run_all_tests(self, verbose: bool = False, include_slow: bool = False) -> List[TestResult]:
        """Run all test suites"""
        print("🚀 Running Complete Test Suite...")
        print("=" * 60)
        
        # Run core test suites
        self.run_unit_tests(verbose)
        self.run_tool_tests(verbose)
        self.run_workflow_tests(verbose)
        
        # E2E tests (may be flaky in CI)
        try:
            self.run_e2e_tests(verbose)
        except Exception as e:
            print(f"⚠️  E2E tests failed to run: {e}")
        
        # Performance tests (optional)
        if include_slow:
            try:
                self.run_performance_tests(verbose)
            except Exception as e:
                print(f"⚠️  Performance tests failed to run: {e}")
        
        return self.results
    
    def generate_report(self, output_file: Optional[Path] = None) -> Dict:
        """Generate comprehensive test report"""
        total_passed = sum(r.passed for r in self.results)
        total_failed = sum(r.failed for r in self.results)
        total_skipped = sum(r.skipped for r in self.results)
        total_duration = sum(r.duration for r in self.results)
        
        # Calculate overall success rate
        total_tests = total_passed + total_failed
        success_rate = (total_passed / total_tests * 100) if total_tests > 0 else 0
        
        report = {
            "timestamp": time.time(),
            "summary": {
                "total_passed": total_passed,
                "total_failed": total_failed,
                "total_skipped": total_skipped,
                "total_duration": total_duration,
                "success_rate": success_rate
            },
            "results": [asdict(r) for r in self.results],
            "environment": {
                "python_version": sys.version,
                "platform": sys.platform,
                "project_root": str(self.project_root)
            }
        }
        
        # Save report if output file specified
        if output_file:
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, 'w') as f:
                json.dump(report, f, indent=2)
            print(f"📄 Test report saved to: {output_file}")
        
        return report
    
    def print_final_summary(self):
        """Print final test summary"""
        print("\n" + "=" * 60)
        print("📊 FINAL TEST SUMMARY")
        print("=" * 60)
        
        total_passed = sum(r.passed for r in self.results)
        total_failed = sum(r.failed for r in self.results)
        total_skipped = sum(r.skipped for r in self.results)
        total_duration = sum(r.duration for r in self.results)
        
        print(f"Total Tests: {total_passed + total_failed + total_skipped}")
        print(f"✅ Passed: {total_passed}")
        print(f"❌ Failed: {total_failed}")
        print(f"⏭️  Skipped: {total_skipped}")
        print(f"⏱️  Duration: {total_duration:.2f}s")
        
        if total_passed + total_failed > 0:
            success_rate = total_passed / (total_passed + total_failed) * 100
            print(f"📈 Success Rate: {success_rate:.1f}%")
        
        # Coverage summary
        unit_result = next((r for r in self.results if r.suite == "unit"), None)
        if unit_result and unit_result.coverage:
            print(f"📊 Code Coverage: {unit_result.coverage:.1f}%")
        
        print("=" * 60)
        
        # Final verdict
        if total_failed == 0:
            print("🎉 ALL TESTS PASSED! 🎉")
            return True
        else:
            print(f"❌ {total_failed} TEST(S) FAILED")
            return False
    
    def _parse_pytest_output(self, output: str) -> tuple:
        """Parse pytest output to extract test counts"""
        passed = failed = skipped = 0
        
        lines = output.split('\n')
        for line in lines:
            if 'passed' in line or 'failed' in line or 'skipped' in line:
                # Look for pytest summary line
                if '=====' in line or '-----' in line:
                    continue
                    
                # Parse individual counts
                if 'passed' in line:
                    try:
                        passed = int(line.split()[0])
                    except (ValueError, IndexError):
                        pass
                
                # Try to parse full summary line like "5 passed, 2 failed, 1 skipped"
                parts = line.split(',')
                for part in parts:
                    part = part.strip()
                    if 'passed' in part:
                        try:
                            passed = int(part.split()[0])
                        except (ValueError, IndexError):
                            pass
                    elif 'failed' in part:
                        try:
                            failed = int(part.split()[0])
                        except (ValueError, IndexError):
                            pass
                    elif 'skipped' in part:
                        try:
                            skipped = int(part.split()[0])
                        except (ValueError, IndexError):
                            pass
        
        return passed, failed, skipped
    
    def _get_coverage_percentage(self) -> Optional[float]:
        """Get coverage percentage from coverage.json"""
        try:
            coverage_file = self.project_root / "coverage.json"
            if coverage_file.exists():
                with open(coverage_file) as f:
                    data = json.load(f)
                    return data.get("totals", {}).get("percent_covered")
        except Exception:
            pass
        return None
    
    def _print_test_summary(self, suite_name: str, result: TestResult):
        """Print summary for a test suite"""
        total = result.passed + result.failed + result.skipped
        status = "✅ PASSED" if result.failed == 0 else "❌ FAILED"
        
        print(f"  {status} - {suite_name}")
        print(f"    Tests: {total} (✅ {result.passed}, ❌ {result.failed}, ⏭️ {result.skipped})")
        print(f"    Duration: {result.duration:.2f}s")
        if result.coverage:
            print(f"    Coverage: {result.coverage:.1f}%")
        print()


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(description="ChatSpatial Test Suite Runner")
    
    parser.add_argument(
        "--suite", 
        choices=["unit", "tool", "workflow", "e2e", "performance", "all"],
        default="all",
        help="Test suite to run"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    
    parser.add_argument(
        "--include-slow",
        action="store_true", 
        help="Include slow/performance tests"
    )
    
    parser.add_argument(
        "--report",
        type=Path,
        help="Save test report to file"
    )
    
    parser.add_argument(
        "--ci",
        action="store_true",
        help="CI mode - exit with error code on test failure"
    )
    
    args = parser.parse_args()
    
    # Setup test runner
    project_root = Path(__file__).parent
    runner = TestRunner(project_root)
    
    # Run specified test suite
    if args.suite == "unit":
        runner.run_unit_tests(args.verbose)
    elif args.suite == "tool":
        runner.run_tool_tests(args.verbose)
    elif args.suite == "workflow":
        runner.run_workflow_tests(args.verbose)
    elif args.suite == "e2e":
        runner.run_e2e_tests(args.verbose)
    elif args.suite == "performance":
        runner.run_performance_tests(args.verbose)
    else:  # all
        runner.run_all_tests(args.verbose, args.include_slow)
    
    # Generate report
    if args.report:
        runner.generate_report(args.report)
    
    # Print final summary and exit
    success = runner.print_final_summary()
    
    if args.ci and not success:
        sys.exit(1)
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())