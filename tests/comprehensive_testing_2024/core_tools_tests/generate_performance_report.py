#!/usr/bin/env python3
"""
Integration Performance Report Generator

Generates comprehensive reports from integration test results,
including visualizations and recommendations.

Author: ChatSpatial Development Team
Created: 2024-08-24
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, Any, List
import argparse
from datetime import datetime

class PerformanceReportGenerator:
    """Generates comprehensive performance reports from test results"""
    
    def __init__(self, results_file: Path):
        self.results_file = Path(results_file)
        self.results = self.load_results()
        self.report_dir = self.results_file.parent
        
    def load_results(self) -> Dict[str, Any]:
        """Load test results from JSON file"""
        try:
            with open(self.results_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            raise Exception(f"Failed to load results file: {e}")
    
    def generate_summary_statistics(self) -> Dict[str, Any]:
        """Generate summary statistics"""
        integration_results = self.results['integration_results']
        perf_evals = self.results['performance_evaluations']
        
        # Basic statistics
        total_tests = len(integration_results)
        successful_tests = sum(1 for r in integration_results if r['success'])
        failed_tests = total_tests - successful_tests
        
        # Performance statistics
        performance_passes = sum(1 for p in perf_evals if p.get('overall_performance') == 'PASS')
        time_passes = sum(1 for p in perf_evals if p.get('time_performance') == 'PASS')
        memory_passes = sum(1 for p in perf_evals if p.get('memory_performance') == 'PASS')
        
        # Timing statistics
        times = [r['total_time'] for r in integration_results if r['success']]
        memories = [r['peak_memory_mb'] for r in integration_results if r['success']]
        
        stats = {
            'test_counts': {
                'total': total_tests,
                'successful': successful_tests,
                'failed': failed_tests,
                'success_rate': (successful_tests / total_tests * 100) if total_tests > 0 else 0
            },
            'performance_counts': {
                'overall_passes': performance_passes,
                'time_passes': time_passes,
                'memory_passes': memory_passes,
                'overall_pass_rate': (performance_passes / total_tests * 100) if total_tests > 0 else 0
            },
            'timing_stats': {
                'mean_time': sum(times) / len(times) if times else 0,
                'max_time': max(times) if times else 0,
                'min_time': min(times) if times else 0
            },
            'memory_stats': {
                'mean_memory': sum(memories) / len(memories) if memories else 0,
                'max_memory': max(memories) if memories else 0,
                'min_memory': min(memories) if memories else 0
            }
        }
        
        return stats
    
    def create_performance_visualizations(self):
        """Create performance visualization charts"""
        plt.style.use('default')
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Extract data
        integration_results = self.results['integration_results']
        perf_evals = self.results['performance_evaluations']
        
        datasets = [r['dataset_name'] for r in integration_results]
        times = [r['total_time'] for r in integration_results]
        memories = [r['peak_memory_mb'] for r in integration_results]
        successes = [r['success'] for r in integration_results]
        
        # 1. Execution Time Comparison
        colors = ['green' if s else 'red' for s in successes]
        ax1.bar(datasets, times, color=colors, alpha=0.7)
        ax1.set_title('Execution Time by Dataset', fontsize=14, fontweight='bold')
        ax1.set_ylabel('Time (seconds)')
        ax1.tick_params(axis='x', rotation=45)
        
        # Add benchmark lines
        benchmarks = self.results.get('benchmarks_used', {}).get('time_benchmarks', {})
        for i, dataset in enumerate(datasets):
            if dataset in benchmarks:
                benchmark_time = benchmarks[dataset]['total']
                ax1.axhline(y=benchmark_time, xmin=i/len(datasets), xmax=(i+1)/len(datasets), 
                           color='blue', linestyle='--', alpha=0.8)
        
        # 2. Memory Usage Comparison
        ax2.bar(datasets, memories, color=colors, alpha=0.7)
        ax2.set_title('Peak Memory Usage by Dataset', fontsize=14, fontweight='bold')
        ax2.set_ylabel('Memory (MB)')
        ax2.tick_params(axis='x', rotation=45)
        
        # Add benchmark lines
        memory_benchmarks = self.results.get('benchmarks_used', {}).get('memory_benchmarks', {})
        for i, dataset in enumerate(datasets):
            if dataset in memory_benchmarks:
                benchmark_memory = memory_benchmarks[dataset]
                ax2.axhline(y=benchmark_memory, xmin=i/len(datasets), xmax=(i+1)/len(datasets),
                           color='blue', linestyle='--', alpha=0.8)
        
        # 3. Performance vs Benchmark Ratio
        time_ratios = []
        memory_ratios = []
        for perf_eval in perf_evals:
            if 'actual_total_time' in perf_eval and 'total_time_benchmark' in perf_eval:
                time_ratio = perf_eval['actual_total_time'] / perf_eval['total_time_benchmark']
                time_ratios.append(time_ratio)
            
            if 'actual_peak_memory_mb' in perf_eval and 'memory_benchmark_mb' in perf_eval:
                memory_ratio = perf_eval['actual_peak_memory_mb'] / perf_eval['memory_benchmark_mb']
                memory_ratios.append(memory_ratio)
        
        if time_ratios:
            x_pos = range(len(time_ratios))
            bars = ax3.bar(x_pos, time_ratios, color=['green' if r <= 1 else 'red' for r in time_ratios], alpha=0.7)
            ax3.axhline(y=1, color='blue', linestyle='--', alpha=0.8, label='Benchmark')
            ax3.set_title('Time Performance vs Benchmark (Ratio)', fontsize=14, fontweight='bold')
            ax3.set_ylabel('Actual/Benchmark Ratio')
            ax3.set_xticks(x_pos)
            ax3.set_xticklabels([p['dataset'] for p in perf_evals], rotation=45)
            ax3.legend()
        
        # 4. Success Rate by Step
        step_names = ['preprocessing_workflow', 'spatial_analysis_workflow', 'clustering_workflow', 'visualization_workflow']
        step_success_rates = []
        
        for step_name in step_names:
            successes = 0
            total = 0
            for result in integration_results:
                for step in result['step_metrics']:
                    if step['test_name'] == step_name:
                        total += 1
                        if step['success']:
                            successes += 1
            
            success_rate = (successes / total * 100) if total > 0 else 0
            step_success_rates.append(success_rate)
        
        bars = ax4.bar(range(len(step_names)), step_success_rates, 
                      color=['green' if sr >= 85 else 'orange' if sr >= 70 else 'red' for sr in step_success_rates],
                      alpha=0.7)
        ax4.set_title('Success Rate by Workflow Step', fontsize=14, fontweight='bold')
        ax4.set_ylabel('Success Rate (%)')
        ax4.set_xticks(range(len(step_names)))
        ax4.set_xticklabels([s.replace('_workflow', '').title() for s in step_names], rotation=45)
        ax4.set_ylim(0, 100)
        
        # Add value labels on bars
        for bar, rate in zip(bars, step_success_rates):
            ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                    f'{rate:.1f}%', ha='center', va='bottom')
        
        plt.tight_layout()
        
        # Save chart
        chart_file = self.report_dir / f"performance_charts_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(chart_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"📊 Performance charts saved to: {chart_file}")
        return chart_file
    
    def generate_detailed_report(self) -> str:
        """Generate detailed markdown report"""
        stats = self.generate_summary_statistics()
        
        report_lines = [
            "# ChatSpatial Core Tools Integration Performance Report",
            "",
            f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Test Results File:** {self.results_file.name}",
            "",
            "## Executive Summary",
            "",
            f"- **Total Tests:** {stats['test_counts']['total']}",
            f"- **Success Rate:** {stats['test_counts']['success_rate']:.1f}% ({stats['test_counts']['successful']}/{stats['test_counts']['total']})",
            f"- **Performance Pass Rate:** {stats['performance_counts']['overall_pass_rate']:.1f}%",
            "",
            "### Key Findings",
            ""
        ]
        
        # Add findings based on results
        if stats['test_counts']['success_rate'] >= 95:
            report_lines.append("✅ **Excellent**: All core functionalities are working reliably")
        elif stats['test_counts']['success_rate'] >= 85:
            report_lines.append("✅ **Good**: Core functionalities are mostly stable with minor issues")  
        elif stats['test_counts']['success_rate'] >= 70:
            report_lines.append("⚠️ **Needs Attention**: Some core functionalities have significant issues")
        else:
            report_lines.append("❌ **Critical**: Major stability issues detected")
        
        if stats['performance_counts']['overall_pass_rate'] >= 80:
            report_lines.append("⚡ **Performance**: Meeting or exceeding benchmark expectations")
        elif stats['performance_counts']['overall_pass_rate'] >= 60:
            report_lines.append("⚡ **Performance**: Acceptable performance with room for optimization")
        else:
            report_lines.append("🐌 **Performance**: Significant performance issues detected")
        
        # Detailed Results Section
        report_lines.extend([
            "",
            "## Detailed Test Results",
            "",
            "### Dataset Performance Summary",
            "",
            "| Dataset | Success | Time (s) | Memory (MB) | Time Status | Memory Status | Overall |",
            "|---------|---------|----------|-------------|-------------|---------------|---------|"
        ])
        
        integration_results = self.results['integration_results']
        perf_evals = self.results['performance_evaluations']
        
        for result, perf_eval in zip(integration_results, perf_evals):
            status = "✅" if result['success'] else "❌"
            dataset = result['dataset_name']
            time_s = f"{result['total_time']:.1f}"
            memory_mb = f"{result['peak_memory_mb']:.0f}"
            
            time_status = perf_eval.get('time_performance', 'N/A')
            memory_status = perf_eval.get('memory_performance', 'N/A')
            overall_status = perf_eval.get('overall_performance', 'N/A')
            
            report_lines.append(f"| {dataset} | {status} | {time_s} | {memory_mb} | {time_status} | {memory_status} | {overall_status} |")
        
        # Step-wise Analysis
        report_lines.extend([
            "",
            "### Workflow Step Analysis",
            ""
        ])
        
        step_analysis = {}
        for result in integration_results:
            for step in result['step_metrics']:
                step_name = step['test_name']
                if step_name not in step_analysis:
                    step_analysis[step_name] = {'total': 0, 'successful': 0, 'times': [], 'errors': []}
                
                step_analysis[step_name]['total'] += 1
                if step['success']:
                    step_analysis[step_name]['successful'] += 1
                    step_analysis[step_name]['times'].append(step['execution_time'])
                else:
                    if step.get('error_message'):
                        step_analysis[step_name]['errors'].append(step['error_message'])
        
        for step_name, data in step_analysis.items():
            success_rate = (data['successful'] / data['total'] * 100) if data['total'] > 0 else 0
            avg_time = sum(data['times']) / len(data['times']) if data['times'] else 0
            
            report_lines.extend([
                f"#### {step_name.replace('_', ' ').title()}",
                f"- **Success Rate:** {success_rate:.1f}% ({data['successful']}/{data['total']})",
                f"- **Average Time:** {avg_time:.2f}s",
                ""
            ])
            
            if data['errors']:
                report_lines.append("**Common Errors:**")
                for error in set(data['errors'][:3]):  # Show up to 3 unique errors
                    report_lines.append(f"- {error}")
                report_lines.append("")
        
        # Recommendations
        report_lines.extend([
            "## Recommendations",
            ""
        ])
        
        if stats['test_counts']['failed'] > 0:
            report_lines.extend([
                "### Stability Improvements",
                "- Investigate failed test cases to identify root causes",
                "- Implement additional error handling for edge cases",
                "- Consider more robust data validation",
                ""
            ])
        
        if stats['performance_counts']['overall_pass_rate'] < 80:
            report_lines.extend([
                "### Performance Optimizations",
                "- Profile memory usage during large dataset processing",
                "- Implement data streaming for memory-intensive operations",
                "- Consider parallel processing for independent operations",
                ""
            ])
        
        # Technical Details
        report_lines.extend([
            "## Technical Details",
            "",
            "### Test Environment",
            f"- **Test Framework:** ChatSpatial Integration Test Suite v1.0",
            f"- **Test Datasets:** {len(self.results.get('benchmarks_used', {}).get('time_benchmarks', {}))} datasets",
            f"- **Performance Benchmarks:** Time and memory thresholds defined",
            "",
            "### Benchmark Thresholds",
            ""
        ])
        
        benchmarks = self.results.get('benchmarks_used', {})
        if 'time_benchmarks' in benchmarks:
            report_lines.append("**Time Benchmarks (seconds):**")
            for dataset, times in benchmarks['time_benchmarks'].items():
                report_lines.append(f"- {dataset}: {times.get('total', 'N/A')}s total")
            report_lines.append("")
        
        if 'memory_benchmarks' in benchmarks:
            report_lines.append("**Memory Benchmarks (MB):**")
            for dataset, memory in benchmarks['memory_benchmarks'].items():
                report_lines.append(f"- {dataset}: {memory}MB peak")
            report_lines.append("")
        
        # Failed Tests Details
        failed_results = [r for r in integration_results if not r['success']]
        if failed_results:
            report_lines.extend([
                "## Failed Tests Analysis",
                ""
            ])
            
            for result in failed_results:
                report_lines.extend([
                    f"### {result['dataset_name']}",
                    f"**Error:** {result.get('error_message', 'Unknown error')}",
                    ""
                ])
                
                # Show step details
                for step in result['step_metrics']:
                    status = "✅" if step['success'] else "❌"
                    report_lines.append(f"- {step['test_name']}: {status}")
                    if not step['success'] and step.get('error_message'):
                        report_lines.append(f"  - Error: {step['error_message']}")
                
                report_lines.append("")
        
        report_content = "\n".join(report_lines)
        
        # Save report
        report_file = self.report_dir / f"integration_performance_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
        with open(report_file, 'w') as f:
            f.write(report_content)
        
        print(f"📄 Detailed report saved to: {report_file}")
        return str(report_file)
    
    def generate_ci_summary(self) -> str:
        """Generate CI/CD-friendly summary"""
        stats = self.generate_summary_statistics()
        
        # Simple pass/fail determination
        overall_pass = (stats['test_counts']['success_rate'] >= 85.0 and 
                       stats['performance_counts']['overall_pass_rate'] >= 60.0)
        
        summary = {
            "overall_status": "PASS" if overall_pass else "FAIL",
            "success_rate": stats['test_counts']['success_rate'],
            "performance_rate": stats['performance_counts']['overall_pass_rate'],
            "total_tests": stats['test_counts']['total'],
            "failed_tests": stats['test_counts']['failed']
        }
        
        # Save CI summary
        ci_file = self.report_dir / "ci_summary.json"
        with open(ci_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"🤖 CI summary saved to: {ci_file}")
        return str(ci_file)
    
    def generate_all_reports(self):
        """Generate all report types"""
        print("📊 Generating comprehensive performance reports...")
        
        # Generate visualizations
        chart_file = self.create_performance_visualizations()
        
        # Generate detailed report
        report_file = self.generate_detailed_report()
        
        # Generate CI summary
        ci_file = self.generate_ci_summary()
        
        print("\n✅ All reports generated successfully!")
        print(f"📊 Charts: {chart_file}")
        print(f"📄 Report: {report_file}")
        print(f"🤖 CI Summary: {ci_file}")
        
        return {
            'charts': chart_file,
            'report': report_file, 
            'ci_summary': ci_file
        }


def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(description='Generate performance reports from integration test results')
    parser.add_argument('results_file', help='Path to JSON results file')
    parser.add_argument('--charts-only', action='store_true', help='Generate only performance charts')
    parser.add_argument('--report-only', action='store_true', help='Generate only detailed report')
    
    args = parser.parse_args()
    
    try:
        generator = PerformanceReportGenerator(args.results_file)
        
        if args.charts_only:
            generator.create_performance_visualizations()
        elif args.report_only:
            generator.generate_detailed_report()
        else:
            generator.generate_all_reports()
            
    except Exception as e:
        print(f"❌ Error generating reports: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())