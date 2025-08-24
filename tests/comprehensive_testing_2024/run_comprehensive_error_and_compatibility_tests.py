"""
Comprehensive Error Handling and Data Compatibility Test Runner

Main script to run all error handling and compatibility tests.
Linus principle: "Test everything, report clearly, fix systematically."

This runs all test suites and generates comprehensive reports.
"""

import os
import sys
import asyncio
import json
import time
from datetime import datetime
from typing import Dict, List, Any

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../'))

# Import error handling tests
from error_handling_tests import (
    ExceptionalInputTests,
    BoundaryConditionTests,
    ResourceLimitationTests,
    DependencyFailureTests,
    NetworkExceptionTests,
    ErrorTestRunner
)

# Import compatibility tests
from data_compatibility_tests import (
    DataFormatCompatibilityTests,
    ParameterCompatibilityTests,
    PlatformCompatibilityTests,
    DependencyCompatibilityTests,
    CompatibilityTestRunner
)


class ComprehensiveTestSuite:
    """
    Master test suite coordinator.
    
    Runs all tests and generates unified reports.
    """
    
    def __init__(self, output_dir: str = None):
        if output_dir is None:
            output_dir = os.path.join(os.path.dirname(__file__), "reports")
        
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.error_runner = ErrorTestRunner()
        self.compatibility_runner = CompatibilityTestRunner()
        
        self.test_start_time = None
        self.test_end_time = None
    
    async def run_all_error_tests(self) -> Dict[str, Any]:
        """Run all error handling tests"""
        print("=" * 80)
        print("RUNNING ERROR HANDLING TESTS")
        print("=" * 80)
        
        all_error_results = []
        
        # Create test instances
        test_classes = [
            ExceptionalInputTests(),
            BoundaryConditionTests(),
            ResourceLimitationTests(),
            DependencyFailureTests(),
            NetworkExceptionTests()
        ]
        
        for test_class in test_classes:
            class_name = test_class.__class__.__name__
            print(f"\\n--- Running {class_name} ---")
            
            test_cases = test_class.get_test_cases()
            results = await self.error_runner.run_test_suite(test_cases)
            all_error_results.extend(results)
        
        # Generate error handling report
        error_report = self.error_runner.generate_summary_report()
        
        return {
            "error_handling_results": all_error_results,
            "error_handling_summary": error_report
        }
    
    async def run_all_compatibility_tests(self) -> Dict[str, Any]:
        """Run all compatibility tests"""
        print("\\n" + "=" * 80)
        print("RUNNING COMPATIBILITY TESTS")
        print("=" * 80)
        
        all_compatibility_results = []
        
        # Create test instances
        test_classes = [
            DataFormatCompatibilityTests(),
            ParameterCompatibilityTests(),
            PlatformCompatibilityTests(),
            DependencyCompatibilityTests()
        ]
        
        for test_class in test_classes:
            class_name = test_class.__class__.__name__
            print(f"\\n--- Running {class_name} ---")
            
            test_cases = test_class.get_test_cases()
            results = await self.compatibility_runner.run_test_suite(test_cases)
            all_compatibility_results.extend(results)
        
        # Generate compatibility report
        compatibility_report = self.compatibility_runner.generate_compatibility_report()
        
        return {
            "compatibility_results": all_compatibility_results,
            "compatibility_summary": compatibility_report
        }
    
    async def run_comprehensive_tests(self) -> Dict[str, Any]:
        """Run all tests and generate comprehensive report"""
        self.test_start_time = datetime.now()
        print(f"Starting comprehensive tests at {self.test_start_time}")
        
        # Run error handling tests
        error_results = await self.run_all_error_tests()
        
        # Run compatibility tests  
        compatibility_results = await self.run_all_compatibility_tests()
        
        self.test_end_time = datetime.now()
        test_duration = (self.test_end_time - self.test_start_time).total_seconds()
        
        # Generate unified report
        comprehensive_report = self.generate_comprehensive_report(
            error_results, compatibility_results, test_duration
        )
        
        # Save reports
        self.save_reports(comprehensive_report, error_results, compatibility_results)
        
        return comprehensive_report
    
    def generate_comprehensive_report(self, error_results: Dict[str, Any], 
                                    compatibility_results: Dict[str, Any],
                                    test_duration: float) -> Dict[str, Any]:
        """Generate comprehensive unified report"""
        
        # Extract key metrics
        error_summary = error_results["error_handling_summary"]["summary"]
        compatibility_summary = compatibility_results["compatibility_summary"]["summary"]
        
        # Calculate overall scores
        total_error_tests = error_summary["total_tests"]
        total_compatibility_tests = compatibility_summary["total_tests"]
        total_tests = total_error_tests + total_compatibility_tests
        
        total_passed = error_summary["passed"] + compatibility_summary["compatible"]
        total_failed = error_summary["failed"] + compatibility_summary["incompatible"]
        
        overall_success_rate = (total_passed / total_tests * 100) if total_tests > 0 else 0
        
        # Create error classification matrix
        error_matrix = self.create_error_classification_matrix(error_results)
        
        # Create compatibility matrix
        compatibility_matrix = self.create_compatibility_matrix(compatibility_results)
        
        # Generate recommendations
        recommendations = self.generate_recommendations(error_results, compatibility_results)
        
        return {
            "test_execution_info": {
                "start_time": self.test_start_time.isoformat(),
                "end_time": self.test_end_time.isoformat(),
                "duration_seconds": test_duration,
                "total_tests": total_tests
            },
            "overall_summary": {
                "total_tests": total_tests,
                "total_passed": total_passed,
                "total_failed": total_failed,
                "overall_success_rate": f"{overall_success_rate:.1f}%",
                "error_handling_success_rate": error_summary["success_rate"],
                "compatibility_success_rate": f"{(compatibility_summary['compatible']/(compatibility_summary['compatible']+compatibility_summary['incompatible'])*100):.1f}%" if (compatibility_summary['compatible']+compatibility_summary['incompatible']) > 0 else "0.0%"
            },
            "error_classification_matrix": error_matrix,
            "compatibility_matrix": compatibility_matrix,
            "critical_issues": self.identify_critical_issues(error_results, compatibility_results),
            "recommendations": recommendations,
            "detailed_results": {
                "error_handling": error_results,
                "compatibility": compatibility_results
            }
        }
    
    def create_error_classification_matrix(self, error_results: Dict[str, Any]) -> Dict[str, Any]:
        """Create error classification matrix"""
        matrix = {
            "by_category": {},
            "by_severity": {},
            "failure_patterns": []
        }
        
        # Group by category and severity
        for result in error_results["error_handling_results"]:
            category = result.category.value
            severity = result.severity.value
            
            # By category
            if category not in matrix["by_category"]:
                matrix["by_category"][category] = {"total": 0, "passed": 0, "failed": 0}
            matrix["by_category"][category]["total"] += 1
            if result.passed:
                matrix["by_category"][category]["passed"] += 1
            else:
                matrix["by_category"][category]["failed"] += 1
            
            # By severity
            if severity not in matrix["by_severity"]:
                matrix["by_severity"][severity] = {"total": 0, "passed": 0, "failed": 0}
            matrix["by_severity"][severity]["total"] += 1
            if result.passed:
                matrix["by_severity"][severity]["passed"] += 1
            else:
                matrix["by_severity"][severity]["failed"] += 1
                
                # Track failure patterns
                matrix["failure_patterns"].append({
                    "test": result.test_name,
                    "category": category,
                    "severity": severity,
                    "error": result.error_message
                })
        
        return matrix
    
    def create_compatibility_matrix(self, compatibility_results: Dict[str, Any]) -> Dict[str, Any]:
        """Create compatibility matrix"""
        matrix = {
            "by_category": {},
            "platform_compatibility": {},
            "dependency_versions": {},
            "incompatibility_issues": []
        }
        
        # Process compatibility results
        for result in compatibility_results["compatibility_results"]:
            category = result.category.value
            
            # By category
            if category not in matrix["by_category"]:
                matrix["by_category"][category] = {"total": 0, "compatible": 0, "incompatible": 0}
            matrix["by_category"][category]["total"] += 1
            if result.compatible:
                matrix["by_category"][category]["compatible"] += 1
            else:
                matrix["by_category"][category]["incompatible"] += 1
                matrix["incompatibility_issues"].append({
                    "test": result.test_name,
                    "category": category,
                    "platform": result.platform_tested,
                    "version": result.version_tested,
                    "error": result.error_message
                })
            
            # Platform info
            if result.platform_tested:
                platform = result.platform_tested
                if platform not in matrix["platform_compatibility"]:
                    matrix["platform_compatibility"][platform] = {"compatible": 0, "incompatible": 0}
                if result.compatible:
                    matrix["platform_compatibility"][platform]["compatible"] += 1
                else:
                    matrix["platform_compatibility"][platform]["incompatible"] += 1
            
            # Version info
            if result.version_tested and result.version_tested != "default":
                version_key = f"{category}_{result.version_tested}"
                matrix["dependency_versions"][version_key] = result.compatible
        
        return matrix
    
    def identify_critical_issues(self, error_results: Dict[str, Any], 
                                compatibility_results: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Identify critical issues that need immediate attention"""
        critical_issues = []
        
        # Critical error handling failures
        for result in error_results["error_handling_results"]:
            if not result.passed and result.severity.value == "critical":
                critical_issues.append({
                    "type": "error_handling",
                    "severity": "critical",
                    "test": result.test_name,
                    "category": result.category.value,
                    "description": result.error_message,
                    "impact": "System may crash or behave unpredictably"
                })
        
        # Critical compatibility issues
        for result in compatibility_results["compatibility_results"]:
            if not result.compatible and result.category.value in ["data_format", "dependency_version"]:
                critical_issues.append({
                    "type": "compatibility",
                    "severity": "high",
                    "test": result.test_name,
                    "category": result.category.value,
                    "platform": result.platform_tested,
                    "version": result.version_tested,
                    "description": result.error_message,
                    "impact": "May not work in production environments"
                })
        
        return critical_issues
    
    def generate_recommendations(self, error_results: Dict[str, Any], 
                               compatibility_results: Dict[str, Any]) -> List[str]:
        """Generate actionable recommendations"""
        recommendations = []
        
        error_summary = error_results["error_handling_summary"]["summary"]
        compatibility_summary = compatibility_results["compatibility_summary"]["summary"]
        
        # Error handling recommendations
        if error_summary["failed"] > 0:
            recommendations.append(
                f"🔧 Fix {error_summary['failed']} error handling issues before production deployment"
            )
        
        # High severity failures
        high_severity_failures = sum(1 for r in error_results["error_handling_results"] 
                                   if not r.passed and r.severity.value in ["critical", "high"])
        if high_severity_failures > 0:
            recommendations.append(
                f"⚠️  Address {high_severity_failures} high/critical severity error handling issues immediately"
            )
        
        # Compatibility recommendations
        if compatibility_summary["incompatible"] > 0:
            recommendations.append(
                f"🔗 Resolve {compatibility_summary['incompatible']} compatibility issues"
            )
        
        # Dependency-specific recommendations
        dep_issues = sum(1 for r in compatibility_results["compatibility_results"] 
                        if not r.compatible and r.category.value == "dependency_version")
        if dep_issues > 0:
            recommendations.append(
                f"📦 Update or fix {dep_issues} dependency version compatibility issues"
            )
        
        # Platform-specific recommendations
        platform_issues = sum(1 for r in compatibility_results["compatibility_results"] 
                             if not r.compatible and r.category.value == "platform")
        if platform_issues > 0:
            recommendations.append(
                f"🖥️  Test and fix {platform_issues} platform-specific issues"
            )
        
        # Overall recommendations
        if len(recommendations) == 0:
            recommendations.append("✅ All tests passed! System appears robust and compatible.")
        else:
            recommendations.append(
                "🧪 Run tests regularly during development to catch regressions early"
            )
            recommendations.append(
                "📋 Consider adding the failed test cases to your CI/CD pipeline"
            )
        
        return recommendations
    
    def save_reports(self, comprehensive_report: Dict[str, Any], 
                    error_results: Dict[str, Any], 
                    compatibility_results: Dict[str, Any]):
        """Save all reports to files"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save comprehensive report
        comprehensive_file = os.path.join(
            self.output_dir, 
            f"comprehensive_test_report_{timestamp}.json"
        )
        with open(comprehensive_file, 'w') as f:
            json.dump(comprehensive_report, f, indent=2, default=str)
        
        # Save detailed error results
        error_file = os.path.join(
            self.output_dir,
            f"error_handling_detailed_{timestamp}.json"
        )
        with open(error_file, 'w') as f:
            json.dump(error_results, f, indent=2, default=str)
        
        # Save detailed compatibility results
        compatibility_file = os.path.join(
            self.output_dir,
            f"compatibility_detailed_{timestamp}.json"
        )
        with open(compatibility_file, 'w') as f:
            json.dump(compatibility_results, f, indent=2, default=str)
        
        # Generate human-readable summary
        summary_file = os.path.join(
            self.output_dir,
            f"test_summary_{timestamp}.md"
        )
        self.generate_markdown_summary(comprehensive_report, summary_file)
        
        print(f"\\n📄 Reports saved to {self.output_dir}")
        print(f"   - Comprehensive: {os.path.basename(comprehensive_file)}")
        print(f"   - Error Details: {os.path.basename(error_file)}")
        print(f"   - Compatibility: {os.path.basename(compatibility_file)}")
        print(f"   - Summary: {os.path.basename(summary_file)}")
    
    def generate_markdown_summary(self, report: Dict[str, Any], output_file: str):
        """Generate human-readable markdown summary"""
        with open(output_file, 'w') as f:
            f.write(f"# ChatSpatial Error Handling & Compatibility Test Report\\n\\n")
            f.write(f"Generated: {report['test_execution_info']['start_time']}\\n")
            f.write(f"Duration: {report['test_execution_info']['duration_seconds']:.1f} seconds\\n\\n")
            
            # Overall Summary
            f.write("## 📊 Overall Summary\\n\\n")
            summary = report['overall_summary']
            f.write(f"- **Total Tests**: {summary['total_tests']}\\n")
            f.write(f"- **Passed**: {summary['total_passed']}\\n")
            f.write(f"- **Failed**: {summary['total_failed']}\\n")
            f.write(f"- **Success Rate**: {summary['overall_success_rate']}\\n\\n")
            
            # Critical Issues
            if report['critical_issues']:
                f.write("## 🚨 Critical Issues\\n\\n")
                for issue in report['critical_issues']:
                    f.write(f"- **{issue['test']}** ({issue['type']})\\n")
                    f.write(f"  - Category: {issue['category']}\\n")
                    f.write(f"  - Impact: {issue['impact']}\\n")
                    f.write(f"  - Details: {issue['description']}\\n\\n")
            
            # Recommendations
            f.write("## 💡 Recommendations\\n\\n")
            for rec in report['recommendations']:
                f.write(f"- {rec}\\n")
            
            f.write("\\n\\n---\\n")
            f.write("*Report generated by ChatSpatial Comprehensive Test Suite*\\n")


async def main():
    """Main function to run comprehensive tests"""
    print("🧪 ChatSpatial Comprehensive Error Handling & Compatibility Tests")
    print("="*80)
    
    # Create test suite
    suite = ComprehensiveTestSuite()
    
    try:
        # Run all tests
        report = await suite.run_comprehensive_tests()
        
        # Print summary
        print("\\n" + "="*80)
        print("🎯 FINAL SUMMARY")
        print("="*80)
        
        overall = report['overall_summary']
        print(f"Total Tests: {overall['total_tests']}")
        print(f"Success Rate: {overall['overall_success_rate']}")
        print(f"Critical Issues: {len(report['critical_issues'])}")
        
        if report['critical_issues']:
            print("\\n⚠️  Critical Issues Found:")
            for issue in report['critical_issues'][:5]:  # Show first 5
                print(f"  - {issue['test']}: {issue['description']}")
        
        print("\\n🎉 Test execution completed successfully!")
        
        return 0
        
    except Exception as e:
        print(f"\\n❌ Test execution failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit_code = asyncio.run(main())