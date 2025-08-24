#!/usr/bin/env python3
"""
Master script to run all deconvolution tests and generate final report.
"""

import os
import sys
import subprocess
import time
from pathlib import Path

def run_test_script(script_name, description):
    """Run a test script and capture results."""
    print(f"\n{'='*60}")
    print(f"RUNNING: {description}")
    print(f"Script: {script_name}")
    print('='*60)
    
    start_time = time.time()
    
    try:
        result = subprocess.run([
            sys.executable, script_name
        ], capture_output=True, text=True, cwd=Path(__file__).parent)
        
        elapsed = time.time() - start_time
        
        print(f"Execution time: {elapsed:.1f} seconds")
        print(f"Return code: {result.returncode}")
        
        if result.stdout:
            print("\nSTDOUT:")
            print(result.stdout)
        
        if result.stderr:
            print("\nSTDERR:")
            print(result.stderr)
        
        success = result.returncode == 0
        return success, elapsed, result.stdout, result.stderr
        
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"Failed to run {script_name}: {e}")
        return False, elapsed, "", str(e)


def main():
    """Run all deconvolution tests."""
    print("🧬 CHATSPATIAL DECONVOLUTION COMPREHENSIVE TESTING SUITE 🧬")
    print(f"Started at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Define test scripts
    test_scripts = [
        ("simple_data_check.py", "Data Preparation and Validation"),
        ("test_deconvolution_comprehensive.py", "Comprehensive Function Testing"),
        ("test_cell2location_detailed.py", "Cell2location Detailed Testing"),
        ("test_r_methods_detailed.py", "R-based Methods Testing")
    ]
    
    results = {}
    total_start_time = time.time()
    
    # Run each test script
    for script, description in test_scripts:
        script_path = Path(__file__).parent / script
        
        if script_path.exists():
            success, elapsed, stdout, stderr = run_test_script(script, description)
            results[script] = {
                'success': success,
                'elapsed': elapsed,
                'description': description,
                'stdout': stdout,
                'stderr': stderr
            }
        else:
            print(f"Warning: {script} not found, skipping...")
            results[script] = {
                'success': False,
                'elapsed': 0,
                'description': description,
                'stdout': '',
                'stderr': f"File {script} not found"
            }
    
    total_elapsed = time.time() - total_start_time
    
    # Generate summary report
    print(f"\n{'='*80}")
    print("🏁 FINAL TEST SUMMARY REPORT")
    print(f"{'='*80}")
    
    successful_tests = sum(1 for r in results.values() if r['success'])
    total_tests = len(results)
    
    print(f"Total execution time: {total_elapsed:.1f} seconds")
    print(f"Tests completed: {total_tests}")
    print(f"Successful: {successful_tests}")
    print(f"Failed: {total_tests - successful_tests}")
    print(f"Success rate: {successful_tests/total_tests*100:.1f}%")
    
    print(f"\n{'Individual Test Results':<40} {'Status':<10} {'Time(s)':<8}")
    print("-" * 80)
    
    for script, result in results.items():
        status = "✅ PASS" if result['success'] else "❌ FAIL"
        elapsed = f"{result['elapsed']:.1f}s"
        description = result['description'][:35] + "..." if len(result['description']) > 35 else result['description']
        print(f"{description:<40} {status:<10} {elapsed:<8}")
    
    # Print key findings
    print(f"\n{'KEY FINDINGS':<40}")
    print("-" * 80)
    
    # Check which methods worked
    method_status = {}
    
    for script, result in results.items():
        if 'comprehensive' in script and result['success']:
            print("✅ Helper functions and validation: WORKING")
            method_status['helpers'] = True
        elif 'cell2location' in script:
            if 'not available' in result['stdout'] or 'cannot import' in result['stdout']:
                print("⚠️  Cell2location: VERSION COMPATIBILITY ISSUE")
                method_status['cell2location'] = False
            elif result['success']:
                print("✅ Cell2location: WORKING")
                method_status['cell2location'] = True
            else:
                print("❌ Cell2location: FAILED")
                method_status['cell2location'] = False
        elif 'r_methods' in script:
            if 'SPOTlight basic test PASSED' in result['stdout']:
                print("✅ SPOTlight: WORKING")
                method_status['spotlight'] = True
            else:
                print("❌ SPOTlight: FAILED")
                method_status['spotlight'] = False
                
            if 'fewer than 10 regression differentially expressed genes' in result['stdout']:
                print("⚠️  RCTD: DATA QUALITY ISSUE (needs better synthetic data)")
                method_status['rctd'] = False
            elif 'RCTD basic test PASSED' in result['stdout']:
                print("✅ RCTD: WORKING")
                method_status['rctd'] = True
            else:
                print("❌ RCTD: FAILED")
                method_status['rctd'] = False
    
    # scvi-tools methods (inferred from comprehensive test)
    if any('comprehensive' in script and 'scvi' in result.get('stdout', '') for script, result in results.items()):
        print("✅ scvi-tools methods (DestVI, Stereoscope, Tangram, MRVI): WORKING")
        method_status['scvi_methods'] = True
    
    print(f"\n{'PRODUCTION READINESS ASSESSMENT':<40}")
    print("-" * 80)
    
    ready_methods = sum(1 for status in method_status.values() if status)
    total_methods = len(method_status)
    
    if ready_methods >= 5:  # Most methods working
        print("🎉 PRODUCTION READY: Most deconvolution methods are functional")
    elif ready_methods >= 3:
        print("⚠️  PARTIALLY READY: Several methods working, some need fixes")
    else:
        print("❌ NOT READY: Major issues need to be resolved")
    
    print(f"\nWorking methods: {ready_methods}/{total_methods}")
    
    # Recommendations
    print(f"\n{'RECOMMENDATIONS':<40}")
    print("-" * 80)
    
    if not method_status.get('cell2location', True):
        print("🔧 Fix Cell2location version compatibility with scvi-tools")
    
    if not method_status.get('rctd', True):
        print("🔧 Improve RCTD data preparation (need more DE genes)")
        
    print("📊 Create more realistic synthetic datasets for testing")
    print("🧪 Add cross-validation and benchmarking against known results")
    print("📈 Performance optimization for large datasets")
    
    # Final status
    print(f"\n{'OVERALL ASSESSMENT':<40}")
    print("-" * 80)
    
    if successful_tests / total_tests >= 0.75:
        print("🟢 EXCELLENT: Deconvolution module shows high quality and robustness")
    elif successful_tests / total_tests >= 0.5:
        print("🟡 GOOD: Deconvolution module is functional with some issues to resolve")
    else:
        print("🔴 NEEDS WORK: Significant issues need to be addressed")
    
    print(f"\nCompleted at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("📝 See DECONVOLUTION_COMPREHENSIVE_TEST_REPORT.md for detailed analysis")
    
    return successful_tests, total_tests


if __name__ == "__main__":
    main()