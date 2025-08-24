#!/usr/bin/env python3
"""
ChatSpatial MCP 协议测试设置验证脚本

快速验证测试环境是否正确配置，所有依赖是否可用。
运行此脚本以确保测试套件可以正常执行。

作者: Linus 风格的环境验证
"""

import sys
import importlib
from pathlib import Path
from typing import List, Tuple

def check_python_version() -> Tuple[bool, str]:
    """检查Python版本"""
    if sys.version_info < (3, 8):
        return False, f"Python 3.8+ required, found {sys.version_info.major}.{sys.version_info.minor}"
    return True, f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

def check_required_packages() -> List[Tuple[str, bool, str]]:
    """检查必需的包"""
    required_packages = [
        'pytest',
        'httpx', 
        'fastapi',
        'pydantic',
        'requests',
        'psutil',
        'pathlib'
    ]
    
    results = []
    for package in required_packages:
        try:
            importlib.import_module(package)
            results.append((package, True, "OK"))
        except ImportError as e:
            results.append((package, False, str(e)))
    
    return results

def check_test_files() -> List[Tuple[str, bool, str]]:
    """检查测试文件存在性"""
    test_files = [
        "test_server_startup.py",
        "test_tool_registration.py",
        "test_parameter_validation.py", 
        "test_error_responses.py",
        "test_http_transport.py",
        "run_protocol_tests.py"
    ]
    
    results = []
    base_path = Path(__file__).parent
    
    for test_file in test_files:
        file_path = base_path / test_file
        if file_path.exists():
            results.append((test_file, True, f"Found: {file_path}"))
        else:
            results.append((test_file, False, f"Missing: {file_path}"))
    
    return results

def check_chatspatial_imports() -> List[Tuple[str, bool, str]]:
    """检查ChatSpatial模块导入"""
    # 添加项目路径
    project_root = Path(__file__).parents[4]
    sys.path.insert(0, str(project_root))
    
    chatspatial_modules = [
        'chatspatial.server',
        'chatspatial.spatial_mcp_adapter',
        'chatspatial.http_server',
        'chatspatial.mcp.errors',
        'chatspatial.utils.mcp_parameter_handler'
    ]
    
    results = []
    for module in chatspatial_modules:
        try:
            importlib.import_module(module)
            results.append((module, True, "Imported successfully"))
        except ImportError as e:
            results.append((module, False, f"Import failed: {str(e)[:100]}"))
    
    return results

def check_test_runner() -> Tuple[bool, str]:
    """检查测试运行器"""
    try:
        from run_protocol_tests import ProtocolTestRunner
        runner = ProtocolTestRunner()
        return True, f"Test runner ready with {len(runner.test_modules)} modules"
    except Exception as e:
        return False, f"Test runner failed: {str(e)}"

def main():
    """主验证函数"""
    print("=" * 70)
    print("ChatSpatial MCP Protocol Test Environment Validation")
    print("=" * 70)
    
    all_good = True
    
    # 检查Python版本
    print("\n1. Python Version Check:")
    py_ok, py_msg = check_python_version()
    print(f"   {'✅' if py_ok else '❌'} {py_msg}")
    if not py_ok:
        all_good = False
    
    # 检查必需包
    print("\n2. Required Packages:")
    package_results = check_required_packages()
    for package, ok, msg in package_results:
        print(f"   {'✅' if ok else '❌'} {package}: {msg}")
        if not ok:
            all_good = False
    
    # 检查测试文件
    print("\n3. Test Files:")
    file_results = check_test_files()
    for test_file, ok, msg in file_results:
        print(f"   {'✅' if ok else '❌'} {test_file}")
        if not ok:
            all_good = False
    
    # 检查ChatSpatial导入
    print("\n4. ChatSpatial Module Imports:")
    import_results = check_chatspatial_imports()
    for module, ok, msg in import_results:
        print(f"   {'✅' if ok else '❌'} {module}")
        if not ok:
            print(f"      {msg}")
            all_good = False
    
    # 检查测试运行器
    print("\n5. Test Runner:")
    runner_ok, runner_msg = check_test_runner()
    print(f"   {'✅' if runner_ok else '❌'} {runner_msg}")
    if not runner_ok:
        all_good = False
    
    # 总结
    print("\n" + "=" * 70)
    if all_good:
        print("🎉 ALL CHECKS PASSED! Test environment is ready.")
        print("\nYou can now run the protocol tests:")
        print("  python run_protocol_tests.py --verbose")
        print("  python run_protocol_tests.py --modules test_server_startup")
        return 0
    else:
        print("❌ SOME CHECKS FAILED! Please fix the issues above.")
        print("\nCommon fixes:")
        print("  - Install missing packages: pip install pytest httpx fastapi pydantic")
        print("  - Ensure you're in the correct directory")
        print("  - Check that ChatSpatial is properly installed")
        return 1

if __name__ == "__main__":
    sys.exit(main())