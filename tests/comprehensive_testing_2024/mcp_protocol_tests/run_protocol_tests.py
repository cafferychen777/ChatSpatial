"""
ChatSpatial MCP 协议测试运行器

本脚本提供完整的MCP协议层测试套件执行器，包括性能基准、详细报告生成和测试结果汇总。
这是Linus风格的测试运行器 - 简单、直接、无废话。

功能：
1. 运行所有5个测试模块
2. 生成性能基准报告
3. 汇总测试结果和指标
4. 保存详细的JSON和HTML报告
5. 提供CI/CD友好的退出代码

使用方式：
    python run_protocol_tests.py [--verbose] [--performance] [--report-dir DIR]

作者: Linus 风格的测试自动化
"""

import argparse
import json
import sys
import time
import subprocess
import importlib.util
from pathlib import Path
from typing import Dict, Any, List, Optional
from datetime import datetime
import traceback
import os

# 测试配置
PROTOCOL_TEST_MODULES = [
    "test_server_startup",
    "test_tool_registration", 
    "test_parameter_validation",
    "test_error_responses",
    "test_http_transport"
]

class ProtocolTestRunner:
    """MCP协议测试运行器"""
    
    def __init__(self, verbose: bool = False, report_dir: Optional[Path] = None):
        self.verbose = verbose
        self.report_dir = report_dir or Path(__file__).parent
        self.report_dir.mkdir(exist_ok=True)
        
        self.test_results = {}
        self.performance_metrics = {}
        self.start_time = time.time()
        self.test_modules = PROTOCOL_TEST_MODULES  # 默认使用所有模块
        
        # 确保报告目录存在
        self.reports_dir = self.report_dir / "reports"
        self.reports_dir.mkdir(exist_ok=True)
    
    def log(self, message: str, level: str = "INFO"):
        """记录日志消息"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        prefix = f"[{timestamp}] {level}:"
        
        if self.verbose or level in ["ERROR", "WARNING"]:
            print(f"{prefix} {message}")
    
    def run_single_test_module(self, module_name: str) -> Dict[str, Any]:
        """运行单个测试模块"""
        self.log(f"Running {module_name}...")
        
        module_start = time.time()
        
        try:
            # 使用pytest运行特定模块
            test_file = Path(__file__).parent / f"{module_name}.py"
            
            if not test_file.exists():
                return {
                    "status": "FAILED",
                    "error": f"Test file {test_file} not found",
                    "execution_time": 0,
                    "tests_run": 0,
                    "tests_passed": 0,
                    "tests_failed": 1
                }
            
            # 运行pytest并捕获输出
            cmd = [
                sys.executable, "-m", "pytest", 
                str(test_file),
                "-v",
                "--tb=short",
                "--json-report",
                f"--json-report-file={self.reports_dir / f'{module_name}_report.json'}"
            ]
            
            if not self.verbose:
                cmd.extend(["--quiet"])
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5分钟超时
            )
            
            module_time = time.time() - module_start
            
            # 解析pytest结果
            json_report_path = self.reports_dir / f"{module_name}_report.json"
            if json_report_path.exists():
                with open(json_report_path, 'r') as f:
                    pytest_report = json.load(f)
                
                # 从pytest报告提取信息
                summary = pytest_report.get("summary", {})
                tests_run = summary.get("total", 0)
                tests_passed = summary.get("passed", 0)
                tests_failed = summary.get("failed", 0)
                tests_error = summary.get("error", 0)
                
                status = "PASSED" if tests_failed == 0 and tests_error == 0 else "FAILED"
            else:
                # Fallback: 根据退出代码判断
                status = "PASSED" if result.returncode == 0 else "FAILED"
                tests_run = 1  # 至少运行了某些测试
                tests_passed = 1 if status == "PASSED" else 0
                tests_failed = 0 if status == "PASSED" else 1
            
            module_result = {
                "status": status,
                "execution_time": module_time,
                "tests_run": tests_run,
                "tests_passed": tests_passed,
                "tests_failed": tests_failed,
                "stdout": result.stdout[-1000:] if result.stdout else "",  # 最后1000字符
                "stderr": result.stderr[-1000:] if result.stderr else ""
            }
            
            if result.returncode != 0 and not json_report_path.exists():
                module_result["error"] = f"pytest exited with code {result.returncode}"
            
            self.log(f"Completed {module_name}: {status} ({tests_passed}/{tests_run} passed, {module_time:.2f}s)")
            
            return module_result
            
        except subprocess.TimeoutExpired:
            return {
                "status": "TIMEOUT",
                "error": f"Test module {module_name} timed out after 300s",
                "execution_time": time.time() - module_start,
                "tests_run": 0,
                "tests_passed": 0,
                "tests_failed": 1
            }
        except Exception as e:
            return {
                "status": "ERROR",
                "error": f"Exception running {module_name}: {str(e)}",
                "execution_time": time.time() - module_start,
                "tests_run": 0,
                "tests_passed": 0,
                "tests_failed": 1
            }
    
    def run_all_tests(self) -> Dict[str, Any]:
        """运行所有协议测试"""
        self.log("Starting ChatSpatial MCP Protocol Test Suite...")
        self.log(f"Test modules: {', '.join(self.test_modules)}")
        
        overall_results = {
            "start_time": datetime.now().isoformat(),
            "modules": {},
            "summary": {
                "total_modules": len(self.test_modules),
                "modules_passed": 0,
                "modules_failed": 0,
                "total_tests": 0,
                "total_passed": 0,
                "total_failed": 0,
                "total_execution_time": 0
            }
        }
        
        # 运行每个测试模块
        for module_name in self.test_modules:
            module_result = self.run_single_test_module(module_name)
            overall_results["modules"][module_name] = module_result
            
            # 更新汇总统计
            if module_result["status"] == "PASSED":
                overall_results["summary"]["modules_passed"] += 1
            else:
                overall_results["summary"]["modules_failed"] += 1
            
            overall_results["summary"]["total_tests"] += module_result.get("tests_run", 0)
            overall_results["summary"]["total_passed"] += module_result.get("tests_passed", 0)
            overall_results["summary"]["total_failed"] += module_result.get("tests_failed", 0)
            overall_results["summary"]["total_execution_time"] += module_result.get("execution_time", 0)
        
        overall_results["end_time"] = datetime.now().isoformat()
        overall_results["total_wall_time"] = time.time() - self.start_time
        
        return overall_results
    
    def generate_performance_benchmark(self) -> Dict[str, Any]:
        """生成性能基准报告"""
        self.log("Generating performance benchmark...")
        
        benchmark = {
            "timestamp": datetime.now().isoformat(),
            "system_info": self.get_system_info(),
            "benchmarks": {}
        }
        
        # 从各个测试模块收集性能数据
        perf_files = list(self.reports_dir.glob("*performance*.json"))
        
        for perf_file in perf_files:
            try:
                with open(perf_file, 'r') as f:
                    perf_data = json.load(f)
                    benchmark["benchmarks"][perf_file.stem] = perf_data
            except Exception as e:
                self.log(f"Failed to load performance data from {perf_file}: {e}", "WARNING")
        
        # 添加一些基本性能指标
        benchmark["benchmarks"]["test_suite_overall"] = {
            "total_execution_time": time.time() - self.start_time,
            "modules_tested": len(self.test_modules),
            "avg_module_time": (time.time() - self.start_time) / len(self.test_modules) if self.test_modules else 0
        }
        
        return benchmark
    
    def get_system_info(self) -> Dict[str, Any]:
        """获取系统信息"""
        import platform
        import psutil
        
        return {
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "processor": platform.processor(),
            "cpu_count": psutil.cpu_count(),
            "memory_gb": round(psutil.virtual_memory().total / (1024**3), 2),
            "hostname": platform.node()
        }
    
    def generate_html_report(self, results: Dict[str, Any]) -> str:
        """生成HTML测试报告"""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>ChatSpatial MCP Protocol Test Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background: #f0f0f0; padding: 20px; border-radius: 5px; }}
        .summary {{ background: #e8f5e8; padding: 15px; margin: 20px 0; border-radius: 5px; }}
        .module {{ background: #f9f9f9; padding: 15px; margin: 10px 0; border-radius: 5px; }}
        .passed {{ color: #008000; font-weight: bold; }}
        .failed {{ color: #cc0000; font-weight: bold; }}
        .timeout {{ color: #ff8800; font-weight: bold; }}
        .error {{ color: #cc0000; font-weight: bold; }}
        pre {{ background: #f0f0f0; padding: 10px; overflow-x: auto; font-size: 12px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>ChatSpatial MCP Protocol Test Report</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Total Execution Time:</strong> {results.get('total_wall_time', 0):.2f} seconds</p>
    </div>
    
    <div class="summary">
        <h2>Test Summary</h2>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            <tr><td>Total Modules</td><td>{results['summary']['total_modules']}</td></tr>
            <tr><td>Modules Passed</td><td class="passed">{results['summary']['modules_passed']}</td></tr>
            <tr><td>Modules Failed</td><td class="failed">{results['summary']['modules_failed']}</td></tr>
            <tr><td>Total Tests</td><td>{results['summary']['total_tests']}</td></tr>
            <tr><td>Tests Passed</td><td class="passed">{results['summary']['total_passed']}</td></tr>
            <tr><td>Tests Failed</td><td class="failed">{results['summary']['total_failed']}</td></tr>
            <tr><td>Success Rate</td><td>{(results['summary']['total_passed'] / max(results['summary']['total_tests'], 1) * 100):.1f}%</td></tr>
        </table>
    </div>
    
    <h2>Module Results</h2>
"""
        
        # 为每个模块添加详细结果
        for module_name, module_result in results["modules"].items():
            status_class = module_result["status"].lower()
            html_content += f"""
    <div class="module">
        <h3>{module_name} - <span class="{status_class}">{module_result["status"]}</span></h3>
        <p><strong>Execution Time:</strong> {module_result.get('execution_time', 0):.2f} seconds</p>
        <p><strong>Tests:</strong> {module_result.get('tests_passed', 0)}/{module_result.get('tests_run', 0)} passed</p>
        
        {f'<p><strong>Error:</strong> <code>{module_result["error"]}</code></p>' if "error" in module_result else ''}
        
        {f'<details><summary>Standard Output</summary><pre>{module_result["stdout"]}</pre></details>' if module_result.get("stdout") else ''}
        {f'<details><summary>Standard Error</summary><pre>{module_result["stderr"]}</pre></details>' if module_result.get("stderr") else ''}
    </div>
"""
        
        html_content += """
</body>
</html>"""
        
        return html_content
    
    def save_reports(self, results: Dict[str, Any], benchmark: Dict[str, Any]):
        """保存测试报告"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 保存JSON报告
        json_report_path = self.reports_dir / f"protocol_test_report_{timestamp}.json"
        with open(json_report_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        # 保存性能基准
        benchmark_path = self.reports_dir / f"performance_benchmark_{timestamp}.json"
        with open(benchmark_path, 'w') as f:
            json.dump(benchmark, f, indent=2, default=str)
        
        # 保存HTML报告
        html_report_path = self.reports_dir / f"protocol_test_report_{timestamp}.html"
        with open(html_report_path, 'w') as f:
            f.write(self.generate_html_report(results))
        
        # 保存最新报告链接
        latest_json = self.reports_dir / "latest_report.json"
        latest_html = self.reports_dir / "latest_report.html"
        latest_benchmark = self.reports_dir / "latest_benchmark.json"
        
        # 创建符号链接或复制文件
        try:
            if latest_json.exists():
                latest_json.unlink()
            if latest_html.exists():
                latest_html.unlink()  
            if latest_benchmark.exists():
                latest_benchmark.unlink()
            
            # 在Windows上，使用复制；在Unix上，使用符号链接
            if os.name == 'nt':
                import shutil
                shutil.copy2(json_report_path, latest_json)
                shutil.copy2(html_report_path, latest_html)
                shutil.copy2(benchmark_path, latest_benchmark)
            else:
                latest_json.symlink_to(json_report_path.name)
                latest_html.symlink_to(html_report_path.name)
                latest_benchmark.symlink_to(benchmark_path.name)
        except Exception as e:
            self.log(f"Failed to create latest report links: {e}", "WARNING")
        
        self.log(f"Reports saved:")
        self.log(f"  JSON: {json_report_path}")
        self.log(f"  HTML: {html_report_path}")
        self.log(f"  Benchmark: {benchmark_path}")
        
        return {
            "json_report": json_report_path,
            "html_report": html_report_path,
            "benchmark_report": benchmark_path
        }
    
    def print_summary(self, results: Dict[str, Any]):
        """打印测试摘要"""
        print("\\n" + "="*70)
        print("CHATSPATIAL MCP PROTOCOL TEST RESULTS")
        print("="*70)
        
        summary = results["summary"]
        
        print(f"Modules: {summary['modules_passed']}/{summary['total_modules']} passed")
        print(f"Tests:   {summary['total_passed']}/{summary['total_tests']} passed")
        print(f"Time:    {summary['total_execution_time']:.2f} seconds")
        
        success_rate = (summary['total_passed'] / max(summary['total_tests'], 1)) * 100
        print(f"Success: {success_rate:.1f}%")
        
        print("\\nModule Results:")
        for module_name, module_result in results["modules"].items():
            status_symbol = "✅" if module_result["status"] == "PASSED" else "❌"
            print(f"  {status_symbol} {module_name}: {module_result['status']} "
                  f"({module_result.get('tests_passed', 0)}/{module_result.get('tests_run', 0)} tests, "
                  f"{module_result.get('execution_time', 0):.2f}s)")
            
            if "error" in module_result:
                print(f"     Error: {module_result['error'][:100]}...")
        
        if summary['modules_failed'] == 0:
            print("\\n🎉 ALL PROTOCOL TESTS PASSED! MCP协议层健康。")
        else:
            print(f"\\n⚠️  {summary['modules_failed']} modules failed. 检查详细报告。")
        
        print("="*70)


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="Run ChatSpatial MCP Protocol Tests")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--performance", "-p", action="store_true", help="Include performance tests")
    parser.add_argument("--report-dir", type=Path, help="Report output directory")
    parser.add_argument("--modules", nargs="+", choices=PROTOCOL_TEST_MODULES, 
                       help="Specific modules to run")
    
    args = parser.parse_args()
    
    # 创建测试运行器
    runner = ProtocolTestRunner(
        verbose=args.verbose,
        report_dir=args.report_dir
    )
    
    try:
        # 如果指定了特定模块，只运行那些模块
        modules_to_run = args.modules if args.modules else PROTOCOL_TEST_MODULES
        if args.modules:
            runner.log(f"Running specific modules: {args.modules}")
        
        # 更新运行器使用的模块列表
        runner.test_modules = modules_to_run
        
        # 运行测试
        results = runner.run_all_tests()
        
        # 生成性能基准
        benchmark = runner.generate_performance_benchmark()
        
        # 保存报告
        report_paths = runner.save_reports(results, benchmark)
        
        # 打印摘要
        runner.print_summary(results)
        
        # 返回适当的退出代码
        if results["summary"]["modules_failed"] == 0:
            runner.log("All protocol tests passed successfully!")
            sys.exit(0)
        else:
            runner.log(f"{results['summary']['modules_failed']} modules failed", "ERROR")
            sys.exit(1)
            
    except KeyboardInterrupt:
        runner.log("Test execution interrupted by user", "WARNING")
        sys.exit(130)
    except Exception as e:
        runner.log(f"Test runner failed with exception: {str(e)}", "ERROR")
        if args.verbose:
            runner.log(traceback.format_exc(), "ERROR")
        sys.exit(1)


if __name__ == "__main__":
    main()