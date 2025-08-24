"""
Comprehensive workflow test runner.

This script runs all integration workflow tests and generates a comprehensive report
of ChatSpatial's end-to-end functionality and performance.
"""

import pytest
import sys
import os
import time
import json
import logging
from pathlib import Path
from typing import Dict, List, Any
import pandas as pd

# Add the chatspatial path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('workflow_tests.log')
    ]
)

logger = logging.getLogger(__name__)


class WorkflowTestRunner:
    """Runs and analyzes integration workflow tests."""
    
    def __init__(self):
        self.test_dir = Path(__file__).parent
        self.results_dir = self.test_dir.parent / "reports"
        self.results_dir.mkdir(exist_ok=True)
        
        self.test_modules = [
            "test_complete_analysis_workflow.py",
            "test_multi_tool_chaining.py", 
            "test_batch_processing_workflow.py",
            "test_real_user_scenarios.py",
            "test_data_flow_validation.py"
        ]
        
        self.results = {}
        self.start_time = None
        self.end_time = None
        
    def run_single_test_module(self, module_name: str) -> Dict[str, Any]:
        """Run a single test module and capture results."""
        logger.info(f"Running test module: {module_name}")
        
        module_path = self.test_dir / module_name
        if not module_path.exists():
            return {
                'module': module_name,
                'status': 'error',
                'error': f"Module {module_name} not found",
                'duration': 0
            }
        
        start_time = time.time()
        
        # Run pytest with detailed output
        exit_code = pytest.main([
            str(module_path),
            "-v",  # Verbose output
            "-s",  # Don't capture output
            "--tb=short",  # Short traceback format
            f"--junit-xml={self.results_dir}/{module_name}.xml"  # JUnit XML output
        ])
        
        duration = time.time() - start_time
        
        result = {
            'module': module_name,
            'status': 'passed' if exit_code == 0 else 'failed',
            'exit_code': exit_code,
            'duration': duration
        }
        
        if exit_code != 0:
            result['error'] = f"Module failed with exit code {exit_code}"
        
        logger.info(f"Completed {module_name} - Status: {result['status']}, Duration: {duration:.2f}s")
        
        return result
    
    def run_all_tests(self) -> Dict[str, Any]:
        """Run all workflow test modules."""
        logger.info("Starting comprehensive workflow testing suite")
        self.start_time = time.time()
        
        overall_results = {
            'start_time': self.start_time,
            'modules': [],
            'summary': {}
        }
        
        for module in self.test_modules:
            try:
                result = self.run_single_test_module(module)
                overall_results['modules'].append(result)
                self.results[module] = result
                
            except Exception as e:
                logger.error(f"Unexpected error running {module}: {str(e)}")
                error_result = {
                    'module': module,
                    'status': 'error',
                    'error': str(e),
                    'duration': 0
                }
                overall_results['modules'].append(error_result)
                self.results[module] = error_result
        
        self.end_time = time.time()
        overall_results['end_time'] = self.end_time
        overall_results['total_duration'] = self.end_time - self.start_time
        
        # Generate summary
        overall_results['summary'] = self._generate_summary()
        
        return overall_results
    
    def _generate_summary(self) -> Dict[str, Any]:
        """Generate summary statistics from test results."""
        if not self.results:
            return {}
        
        total_modules = len(self.results)
        passed_modules = sum(1 for r in self.results.values() if r['status'] == 'passed')
        failed_modules = sum(1 for r in self.results.values() if r['status'] == 'failed')
        error_modules = sum(1 for r in self.results.values() if r['status'] == 'error')
        
        total_duration = sum(r['duration'] for r in self.results.values())
        avg_duration = total_duration / total_modules if total_modules > 0 else 0
        
        success_rate = passed_modules / total_modules if total_modules > 0 else 0
        
        summary = {
            'total_modules': total_modules,
            'passed_modules': passed_modules,
            'failed_modules': failed_modules,
            'error_modules': error_modules,
            'success_rate': success_rate,
            'total_duration': total_duration,
            'avg_duration_per_module': avg_duration,
            'status': 'PASSED' if success_rate >= 0.8 else 'FAILED'
        }
        
        return summary
    
    def generate_detailed_report(self, results: Dict[str, Any]) -> str:
        """Generate detailed HTML report."""
        report_html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>ChatSpatial Integration Workflow Test Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
        .summary {{ background-color: #e8f5e8; padding: 15px; margin: 20px 0; border-radius: 5px; }}
        .summary.failed {{ background-color: #ffe8e8; }}
        .module {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
        .module.passed {{ border-left: 5px solid #4CAF50; }}
        .module.failed {{ border-left: 5px solid #f44336; }}
        .module.error {{ border-left: 5px solid #ff9800; }}
        .metric {{ display: inline-block; margin: 10px 20px 10px 0; }}
        .metric-value {{ font-weight: bold; font-size: 1.2em; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 10px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>ChatSpatial Integration Workflow Test Report</h1>
        <p>Generated on: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(results['start_time']))}</p>
        <p>Total Duration: {results['total_duration']:.2f} seconds</p>
    </div>
    
    <div class="summary {'passed' if results['summary']['status'] == 'PASSED' else 'failed'}">
        <h2>Summary - {results['summary']['status']}</h2>
        <div class="metric">
            <div class="metric-value">{results['summary']['success_rate']:.1%}</div>
            <div>Success Rate</div>
        </div>
        <div class="metric">
            <div class="metric-value">{results['summary']['passed_modules']}</div>
            <div>Passed Modules</div>
        </div>
        <div class="metric">
            <div class="metric-value">{results['summary']['failed_modules']}</div>
            <div>Failed Modules</div>
        </div>
        <div class="metric">
            <div class="metric-value">{results['summary']['error_modules']}</div>
            <div>Error Modules</div>
        </div>
        <div class="metric">
            <div class="metric-value">{results['summary']['avg_duration_per_module']:.1f}s</div>
            <div>Avg Duration</div>
        </div>
    </div>
    
    <h2>Test Module Results</h2>
    <table>
        <tr>
            <th>Module</th>
            <th>Status</th>
            <th>Duration (s)</th>
            <th>Details</th>
        </tr>
"""
        
        for module_result in results['modules']:
            status_class = module_result['status']
            error_info = module_result.get('error', '')
            
            report_html += f"""
        <tr>
            <td>{module_result['module']}</td>
            <td class="{status_class}">{module_result['status'].upper()}</td>
            <td>{module_result['duration']:.2f}</td>
            <td>{error_info}</td>
        </tr>
"""
        
        report_html += """
    </table>
    
    <h2>Test Descriptions</h2>
    <div class="module">
        <h3>test_complete_analysis_workflow.py</h3>
        <p>Tests complete end-to-end spatial transcriptomics analysis workflows from data loading through visualization. Validates basic Visium workflows, seqFISH workflows with fallbacks, synthetic data pipelines, and error recovery mechanisms.</p>
    </div>
    
    <div class="module">
        <h3>test_multi_tool_chaining.py</h3>
        <p>Tests complex workflows that chain multiple ChatSpatial tools together. Validates preprocessing→spatial analysis chains, differential→visualization chains, and complex spatial domain analysis workflows.</p>
    </div>
    
    <div class="module">
        <h3>test_batch_processing_workflow.py</h3>
        <p>Tests batch processing capabilities across multiple datasets. Validates sequential and parallel processing, spatial batch workflows, and error handling across diverse dataset types.</p>
    </div>
    
    <div class="module">
        <h3>test_real_user_scenarios.py</h3>
        <p>Simulates real-world user scenarios including new user exploration, experienced user comprehensive analysis, comparative analysis, troubleshooting sessions, and mixed usage patterns.</p>
    </div>
    
    <div class="module">
        <h3>test_data_flow_validation.py</h3>
        <p>Tests data integrity and consistency across workflow steps. Validates preprocessing data flow, spatial analysis consistency, cross-tool compatibility, and data preservation through complex workflows.</p>
    </div>
    
    <h2>System Information</h2>
    <p>Python Version: """ + sys.version + """</p>
    <p>Test Directory: """ + str(self.test_dir) + """</p>
    
</body>
</html>
"""
        
        return report_html
    
    def save_results(self, results: Dict[str, Any]) -> None:
        """Save test results in multiple formats."""
        timestamp = time.strftime('%Y%m%d_%H%M%S', time.localtime(results['start_time']))
        
        # Save JSON results
        json_path = self.results_dir / f"workflow_test_results_{timestamp}.json"
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        # Save CSV summary
        csv_path = self.results_dir / f"workflow_test_summary_{timestamp}.csv"
        modules_df = pd.DataFrame(results['modules'])
        modules_df.to_csv(csv_path, index=False)
        
        # Save HTML report
        html_path = self.results_dir / f"workflow_test_report_{timestamp}.html"
        html_report = self.generate_detailed_report(results)
        with open(html_path, 'w') as f:
            f.write(html_report)
        
        # Save latest results (overwrites previous)
        latest_json = self.results_dir / "latest_workflow_results.json"
        latest_html = self.results_dir / "latest_workflow_report.html"
        
        with open(latest_json, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        with open(latest_html, 'w') as f:
            f.write(html_report)
        
        logger.info(f"Results saved to:")
        logger.info(f"  JSON: {json_path}")
        logger.info(f"  CSV: {csv_path}")
        logger.info(f"  HTML: {html_path}")
        logger.info(f"  Latest HTML: {latest_html}")


def main():
    """Main function to run all workflow tests."""
    runner = WorkflowTestRunner()
    
    try:
        logger.info("=" * 80)
        logger.info("ChatSpatial Integration Workflow Test Suite")
        logger.info("=" * 80)
        
        results = runner.run_all_tests()
        
        logger.info("=" * 80)
        logger.info("TEST SUITE COMPLETED")
        logger.info("=" * 80)
        
        summary = results['summary']
        logger.info(f"Overall Status: {summary['status']}")
        logger.info(f"Success Rate: {summary['success_rate']:.1%}")
        logger.info(f"Modules Passed: {summary['passed_modules']}/{summary['total_modules']}")
        logger.info(f"Total Duration: {summary['total_duration']:.2f} seconds")
        
        # Save results
        runner.save_results(results)
        
        # Exit with appropriate code
        sys.exit(0 if summary['status'] == 'PASSED' else 1)
        
    except Exception as e:
        logger.error(f"Test suite failed with unexpected error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()