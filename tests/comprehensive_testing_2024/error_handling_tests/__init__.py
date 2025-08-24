"""
ChatSpatial Error Handling Test Suite

Comprehensive testing framework for error handling and failure scenarios.
Follows Linus principle: "Good code has no special cases."

Test Categories:
1. Exceptional Input Handling - Bad data, malformed parameters
2. Boundary Conditions - Empty data, single cell, extreme sizes
3. Resource Limitations - Memory, disk, computation limits
4. Dependency Failures - Third-party library failures with graceful fallbacks
5. Network Exceptions - HTTP mode connectivity issues
"""

from .base_test_framework import (
    ErrorTestCase,
    ErrorCategory,
    ErrorSeverity,
    TestResult,
    ErrorTestRunner
)

from .test_exceptional_inputs import ExceptionalInputTests
from .test_boundary_conditions import BoundaryConditionTests  
from .test_resource_limitations import ResourceLimitationTests
from .test_dependency_failures import DependencyFailureTests
from .test_network_exceptions import NetworkExceptionTests

__all__ = [
    'ErrorTestCase',
    'ErrorCategory', 
    'ErrorSeverity',
    'TestResult',
    'ErrorTestRunner',
    'ExceptionalInputTests',
    'BoundaryConditionTests',
    'ResourceLimitationTests', 
    'DependencyFailureTests',
    'NetworkExceptionTests'
]