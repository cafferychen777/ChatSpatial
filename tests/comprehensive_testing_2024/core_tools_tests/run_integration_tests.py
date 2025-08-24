#!/usr/bin/env python3
"""
ChatSpatial Integration Tests Runner

Simple runner script for core integration tests with automatic report generation.

Usage:
    python run_integration_tests.py [--quick] [--report-only]
    
Author: ChatSpatial Development Team
Created: 2024-08-24
"""

import sys
import os
import subprocess
import argparse
from pathlib import Path

def run_command(cmd, description=""):
    """Run a command and return success status"""
    print(f"🔄 {description}")
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(f"✅ {description} - SUCCESS")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} - FAILED")
        print(f"   Error: {e.stderr}")
        return False

def find_latest_results_file(test_dir):
    """Find the most recent test results file"""
    results_files = list(test_dir.glob("integration_test_results_*.json"))
    if not results_files:
        return None
    return max(results_files, key=lambda x: x.stat().st_mtime)

def main():
    parser = argparse.ArgumentParser(description='Run ChatSpatial integration tests')
    parser.add_argument('--quick', action='store_true', 
                       help='Run quick tests only (demo dataset)')
    parser.add_argument('--report-only', action='store_true',
                       help='Generate reports from existing results only')
    parser.add_argument('--no-reports', action='store_true',
                       help='Skip report generation')
    
    args = parser.parse_args()
    
    # Setup paths
    test_dir = Path('/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/core_tools_tests')
    test_script = test_dir / 'test_integration_performance.py'
    report_script = test_dir / 'generate_performance_report.py'
    
    print("🚀 ChatSpatial Integration Tests Runner")
    print("=" * 50)
    
    # Check if test scripts exist
    if not test_script.exists():
        print(f"❌ Test script not found: {test_script}")
        return 1
    
    success = True
    
    # Run tests (unless report-only mode)
    if not args.report_only:
        print("\n📋 PHASE 1: Running Integration Tests")
        print("-" * 40)
        
        # Prepare test environment
        os.chdir(test_dir)
        
        # Set environment variables
        env_cmd = f"export PYTHONPATH=/Users/apple/Research/SpatialTrans_MCP/chatspatial:$PYTHONPATH"
        
        if args.quick:
            # Modify test script to run only quick demo
            test_cmd = f"{env_cmd} && python {test_script} --quick"
            success = run_command(test_cmd, "Running quick integration test")
        else:
            # Run full test suite
            test_cmd = f"{env_cmd} && python {test_script}"
            success = run_command(test_cmd, "Running full integration test suite")
    
    # Generate reports (unless disabled)
    if not args.no_reports:
        print(f"\n📊 PHASE 2: Generating Reports")
        print("-" * 40)
        
        # Find latest results file
        results_file = find_latest_results_file(test_dir)
        
        if not results_file:
            print("❌ No test results file found. Run tests first.")
            return 1
        
        print(f"📄 Using results file: {results_file.name}")
        
        # Generate reports
        if report_script.exists():
            report_cmd = f"cd {test_dir} && python {report_script} {results_file}"
            report_success = run_command(report_cmd, "Generating performance reports")
            success = success and report_success
        else:
            print(f"⚠️  Report generator not found: {report_script}")
    
    # Final summary
    print(f"\n{'=' * 50}")
    if success:
        print("🎉 ALL OPERATIONS COMPLETED SUCCESSFULLY!")
        
        # Show output files
        print(f"\n📂 Output files in: {test_dir}")
        
        # List recent files
        recent_files = []
        for pattern in ["*.json", "*.csv", "*.md", "*.png"]:
            recent_files.extend(test_dir.glob(pattern))
        
        recent_files.sort(key=lambda x: x.stat().st_mtime, reverse=True)
        
        print("📄 Recent output files:")
        for f in recent_files[:10]:  # Show last 10 files
            size_mb = f.stat().st_size / 1024 / 1024
            print(f"   {f.name} ({size_mb:.1f}MB)")
        
        return 0
    else:
        print("❌ SOME OPERATIONS FAILED - CHECK LOGS ABOVE")
        return 1

if __name__ == "__main__":
    exit(main())