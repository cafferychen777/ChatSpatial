"""
ChatSpatial 高级分析工具测试包

本包包含 ChatSpatial 高级分析工具的综合测试套件，
按照 Linus Torvalds 的软件质量哲学设计。

测试模块：
- test_advanced_analysis_tools.py: 主测试框架
- dependency_validator.py: 依赖验证器
- run_advanced_analysis_tests.py: 测试运行器

设计原则：
1. "好品味" - 简洁的数据结构和清晰的逻辑
2. "Never break userspace" - 全面的fallback机制
3. 实用主义 - 测试真实使用场景
4. 简洁执念 - 消除不必要的复杂性
"""

__version__ = "1.0.0"
__author__ = "ChatSpatial Development Team"

# 导出主要测试组件
from .test_advanced_analysis_tools import AdvancedAnalysisTestFramework
from .dependency_validator import DependencyValidator
from .run_advanced_analysis_tests import AdvancedAnalysisTestRunner

__all__ = [
    'AdvancedAnalysisTestFramework',
    'DependencyValidator', 
    'AdvancedAnalysisTestRunner'
]