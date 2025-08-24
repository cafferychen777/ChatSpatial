#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple test for the smart dependency management system
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, '.')

def test_core_imports():
    """Test that our core modules can be imported"""
    print("Testing core module imports...")
    
    try:
        from chatspatial.utils.smart_import import smart_import
        print("✓ smart_import imported successfully")
        
        # Test basic functionality
        result = smart_import("numpy")
        if result.success:
            print("✓ numpy import test passed")
        else:
            print("✗ numpy import failed: {}".format(result.error))
            
    except Exception as e:
        print("✗ Failed to import smart_import: {}".format(e))
        return False
    
    try:
        from chatspatial.utils.dependency_fallbacks import list_available_fallbacks
        fallbacks = list_available_fallbacks()
        print("✓ dependency_fallbacks imported, found {} fallbacks".format(len(fallbacks)))
    except Exception as e:
        print("✗ Failed to import dependency_fallbacks: {}".format(e))
        return False
        
    return True

def test_basic_functionality():
    """Test basic functionality without complex dependencies"""
    print("\nTesting basic functionality...")
    
    from chatspatial.utils.smart_import import smart_import, get_smart_importer
    
    # Test with a package that should exist
    numpy_result = smart_import("numpy", required=False)
    print("numpy available: {}".format(numpy_result.success))
    
    # Test with a package that probably doesn't exist
    fake_result = smart_import("definitely_not_a_real_package", required=False)
    print("fake package correctly rejected: {}".format(not fake_result.success))
    
    # Test the importer instance
    importer = get_smart_importer()
    print("Smart importer created successfully")
    
    return True

def main():
    """Run simple tests"""
    print("ChatSpatial Smart Dependency System - Simple Test")
    print("=" * 60)
    
    success = True
    
    success &= test_core_imports()
    success &= test_basic_functionality()
    
    print("\n" + "=" * 60)
    if success:
        print("✓ All tests passed! The smart dependency system is working.")
        print("\nNow you can use:")
        print("  - smart_import() for intelligent package importing")
        print("  - Dependency fallbacks for graceful degradation")
        print("  - CLI commands for dependency management")
        return 0
    else:
        print("✗ Some tests failed. Check the errors above.")
        return 1

if __name__ == "__main__":
    exit(main())