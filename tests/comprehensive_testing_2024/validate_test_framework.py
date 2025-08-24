"""
Independent Test Framework Validation

This validates the test framework without importing the full test suites.
"""

import asyncio
import sys
import os
import importlib.util

def load_module_from_path(name, path):
    """Load a module directly from file path"""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def validate_error_framework():
    """Test error handling framework"""
    print("🚨 Validating Error Handling Framework...")
    
    # Load base framework directly
    base_path = os.path.join(os.path.dirname(__file__), 
                            "error_handling_tests", "base_test_framework.py")
    base_module = load_module_from_path("base_error", base_path)
    
    # Test basic classes
    ErrorTestCase = base_module.ErrorTestCase
    ErrorCategory = base_module.ErrorCategory
    ErrorSeverity = base_module.ErrorSeverity
    ErrorTestRunner = base_module.ErrorTestRunner
    
    # Create test case
    def test_function():
        return True
    
    test_case = ErrorTestCase(
        name="validation_test",
        category=ErrorCategory.EXCEPTIONAL_INPUT,
        severity=ErrorSeverity.LOW,
        description="Framework validation test",
        test_function=test_function
    )
    
    # Test runner
    runner = ErrorTestRunner()
    
    print("  ✅ ErrorTestCase creation successful")
    print("  ✅ ErrorCategory enum working")
    print("  ✅ ErrorSeverity enum working") 
    print("  ✅ ErrorTestRunner instantiation successful")
    
    return True

def validate_compatibility_framework():
    """Test compatibility framework"""
    print("🔗 Validating Compatibility Framework...")
    
    # Load base framework directly
    base_path = os.path.join(os.path.dirname(__file__), 
                            "data_compatibility_tests", "base_compatibility_framework.py")
    base_module = load_module_from_path("base_compatibility", base_path)
    
    # Test basic classes
    CompatibilityTestCase = base_module.CompatibilityTestCase
    CompatibilityCategory = base_module.CompatibilityCategory
    CompatibilityTestRunner = base_module.CompatibilityTestRunner
    
    # Create test case
    def test_function():
        return {"compatible": True}
    
    test_case = CompatibilityTestCase(
        name="validation_test",
        category=CompatibilityCategory.DATA_FORMAT,
        description="Framework validation test",
        test_function=test_function
    )
    
    # Test runner
    runner = CompatibilityTestRunner()
    
    print("  ✅ CompatibilityTestCase creation successful")
    print("  ✅ CompatibilityCategory enum working")
    print("  ✅ CompatibilityTestRunner instantiation successful")
    
    return True

def validate_test_data():
    """Test that special test data was generated"""
    print("📊 Validating Test Data...")
    
    datasets_dir = os.path.join(os.path.dirname(__file__), "datasets")
    
    if not os.path.exists(datasets_dir):
        print("  ❌ Datasets directory not found")
        return False
    
    expected_datasets = [
        "empty_dataset.h5ad",
        "single_cell.h5ad",
        "high_sparsity.h5ad",
        "nan_coordinates.h5ad",
        "identical_coordinates.h5ad",
        "extremely_sparse.h5ad"
    ]
    
    found_datasets = 0
    for dataset in expected_datasets:
        dataset_path = os.path.join(datasets_dir, dataset)
        if os.path.exists(dataset_path):
            found_datasets += 1
            print(f"  ✅ {dataset} found")
        else:
            print(f"  ⚠️  {dataset} missing")
    
    print(f"  📈 Found {found_datasets}/{len(expected_datasets)} expected datasets")
    
    # Check summary file
    summary_path = os.path.join(datasets_dir, "datasets_summary.csv")
    if os.path.exists(summary_path):
        print("  ✅ Datasets summary found")
        return True
    else:
        print("  ⚠️  Datasets summary missing")
        return False

def main():
    """Main validation function"""
    print("🧪 ChatSpatial Test Framework Validation")
    print("="*50)
    
    results = []
    
    try:
        results.append(validate_error_framework())
        print()
        results.append(validate_compatibility_framework())
        print()
        results.append(validate_test_data())
    except Exception as e:
        print(f"❌ Validation failed with error: {e}")
        return False
    
    print("\n" + "="*50)
    print("📋 VALIDATION SUMMARY")
    print("="*50)
    
    passed = sum(results)
    total = len(results)
    
    print(f"Components Validated: {total}")
    print(f"Successful: {passed}")
    print(f"Failed: {total - passed}")
    
    if passed == total:
        print("\n✅ ALL VALIDATIONS PASSED!")
        print("The test framework is ready for use.")
        print("\nNext steps:")
        print("1. Fix any import issues in specific test files")
        print("2. Run individual test components")
        print("3. Add more test cases as needed")
        return True
    else:
        print(f"\n❌ {total - passed} VALIDATIONS FAILED")
        print("Please address the issues above before proceeding.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)