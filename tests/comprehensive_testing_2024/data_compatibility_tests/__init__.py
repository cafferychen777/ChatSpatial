"""
ChatSpatial Data Compatibility Test Suite

Systematic testing of data format, parameter, platform, and dependency compatibility.
Linus principle: "Compatibility is a feature, not an afterthought."

Test Categories:
1. Data Format Compatibility - h5ad versions, CSV, HDF5, different encodings
2. Parameter Version Compatibility - backward compatibility of tool parameters  
3. Platform Compatibility - macOS, Linux, Windows behavior differences
4. Dependency Version Compatibility - scanpy, pandas, numpy version combinations

Design Philosophy:
- Test real-world scenarios, not theoretical edge cases
- Clear pass/fail criteria
- Actionable compatibility reports
"""

from .base_compatibility_framework import (
    CompatibilityTestCase,
    CompatibilityCategory, 
    CompatibilityResult,
    CompatibilityTestRunner,
    CompatibilityMatrix
)

from .test_data_format_compatibility import DataFormatCompatibilityTests
from .test_parameter_compatibility import ParameterCompatibilityTests
from .test_platform_compatibility import PlatformCompatibilityTests  
from .test_dependency_compatibility import DependencyCompatibilityTests

__all__ = [
    'CompatibilityTestCase',
    'CompatibilityCategory',
    'CompatibilityResult', 
    'CompatibilityTestRunner',
    'CompatibilityMatrix',
    'DataFormatCompatibilityTests',
    'ParameterCompatibilityTests',
    'PlatformCompatibilityTests',
    'DependencyCompatibilityTests'
]