#!/usr/bin/env python3
"""
ChatSpatial MCP 全面测试执行器

这个脚本会运行所有5个测试模块的完整测试套件。
"""

import os
import sys
import subprocess
import time
import json
import argparse
from pathlib import Path
from typing import Dict, List, Any
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ComprehensiveTestRunner:
    """全面测试运行器"""
    
    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self.results = {}
        self.total_start_time = time.time()
        
        # 测试模块配置
        self.test_modules = {
            'mcp_protocol': {
                'directory': 'mcp_protocol_tests',
                'script': 'run_protocol_tests.py',
                'description': 'MCP协议层测试',
                'timeout': 120
            },
            'tool_functionality': {
                'directory': 'tool_functionality_tests', 
                'script': 'run_comprehensive_tests.py',
                'args': ['--type', 'all', '--verbose'],
                'description': '核心工具功能测试',
                'timeout': 600
            },
            'integration_workflow': {
                'directory': 'integration_workflow_tests',
                'script': 'run_all_workflow_tests.py', 
                'description': '集成工作流测试',
                'timeout': 900
            },
            'performance_stress': {
                'directory': 'performance_stress_tests',
                'script': 'run_performance_tests.py',
                'description': '性能压力测试', 
                'timeout': 1200,
                'optional': True  # 这个模块可能因为Agent超时而不存在
            },
            'error_compatibility': {
                'directory': '.',  # 在根目录
                'script': 'run_comprehensive_error_and_compatibility_tests.py',
                'description': '错误处理和兼容性测试',
                'timeout': 300
            }
        }
    
    def check_module_availability(self) -> Dict[str, bool]:
        """检查各测试模块的可用性"""
        availability = {}
        
        for module_name, config in self.test_modules.items():
            script_path = self.base_dir / config['directory'] / config['script']
            
            if config.get('optional', False):
                # 对于可选模块，检查是否存在
                availability[module_name] = script_path.exists()
                if not script_path.exists():
                    logger.warning(f"可选测试模块 {module_name} 不可用: {script_path}")
            else:
                # 对于必需模块，必须存在
                if not script_path.exists():
                    logger.error(f"必需测试模块 {module_name} 缺失: {script_path}")
                    availability[module_name] = False
                else:
                    availability[module_name] = True
        
        return availability
    
    def run_single_test_module(self, module_name: str, config: Dict[str, Any]) -> Dict[str, Any]:
        """运行单个测试模块"""
        logger.info(f"开始运行: {config['description']}")
        
        # 准备命令
        script_path = self.base_dir / config['directory'] / config['script']
        cmd = [sys.executable, str(script_path)]
        
        # 添加额外参数
        if 'args' in config:
            cmd.extend(config['args'])
        
        # 运行测试
        start_time = time.time()
        
        try:
            result = subprocess.run(
                cmd,
                cwd=self.base_dir / config['directory'],
                capture_output=True,
                text=True,
                timeout=config.get('timeout', 300)
            )
            
            execution_time = time.time() - start_time
            
            return {
                'module': module_name,
                'description': config['description'],
                'success': result.returncode == 0,
                'execution_time': execution_time,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'return_code': result.returncode,
                'command': ' '.join(cmd)
            }
            
        except subprocess.TimeoutExpired:
            execution_time = time.time() - start_time
            logger.error(f"{module_name} 超时 (>{config.get('timeout', 300)}秒)")
            
            return {
                'module': module_name,
                'description': config['description'],
                'success': False,
                'execution_time': execution_time,
                'error': 'TIMEOUT',
                'timeout': config.get('timeout', 300),
                'command': ' '.join(cmd)
            }
            
        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"{module_name} 执行异常: {e}")
            
            return {
                'module': module_name,
                'description': config['description'], 
                'success': False,
                'execution_time': execution_time,
                'error': str(e),
                'command': ' '.join(cmd)
            }
    
    def run_all_tests(self, selected_modules: List[str] = None) -> Dict[str, Any]:
        """运行所有测试"""
        logger.info("=== ChatSpatial MCP 全面测试开始 ===")
        
        # 检查模块可用性
        availability = self.check_module_availability()
        
        # 确定要运行的模块
        if selected_modules:
            modules_to_run = {k: v for k, v in self.test_modules.items() 
                            if k in selected_modules and availability.get(k, False)}
        else:
            modules_to_run = {k: v for k, v in self.test_modules.items() 
                            if availability.get(k, False)}
        
        logger.info(f"将运行 {len(modules_to_run)} 个测试模块:")
        for module_name, config in modules_to_run.items():
            logger.info(f"  - {config['description']}")
        
        # 运行测试
        results = []
        for module_name, config in modules_to_run.items():
            result = self.run_single_test_module(module_name, config)
            results.append(result)
            
            # 实时反馈
            status = "✅ 成功" if result['success'] else "❌ 失败"
            logger.info(f"{config['description']}: {status} ({result['execution_time']:.1f}秒)")
        
        # 生成总结果
        total_time = time.time() - self.total_start_time
        successful_tests = sum(1 for r in results if r['success'])
        total_tests = len(results)
        
        summary = {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_execution_time': total_time,
            'total_modules': total_tests,
            'successful_modules': successful_tests,
            'success_rate': successful_tests / total_tests if total_tests > 0 else 0,
            'overall_success': successful_tests == total_tests,
            'module_results': results,
            'availability': availability
        }
        
        return summary
    
    def generate_report(self, results: Dict[str, Any], output_file: str = None):
        """生成测试报告"""
        if output_file is None:
            output_file = f"comprehensive_test_report_{int(time.time())}.json"
        
        output_path = self.base_dir / 'reports' / output_file
        output_path.parent.mkdir(exist_ok=True)
        
        # 保存JSON报告
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        
        # 生成Markdown报告
        md_report = self.generate_markdown_report(results)
        md_path = output_path.with_suffix('.md')
        
        with open(md_path, 'w', encoding='utf-8') as f:
            f.write(md_report)
        
        logger.info(f"报告已生成:")
        logger.info(f"  JSON: {output_path}")
        logger.info(f"  Markdown: {md_path}")
        
        return output_path
    
    def generate_markdown_report(self, results: Dict[str, Any]) -> str:
        """生成Markdown格式的报告"""
        
        md = f"""# ChatSpatial MCP 全面测试报告

## 测试总结

- **执行时间**: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}
- **总耗时**: {results['total_execution_time']:.1f} 秒
- **测试模块**: {results['successful_modules']}/{results['total_modules']} 成功
- **成功率**: {results['success_rate']:.1%}
- **总体状态**: {'✅ 通过' if results['overall_success'] else '❌ 失败'}

## 模块测试结果

"""
        
        for result in results['module_results']:
            status = "✅ 成功" if result['success'] else "❌ 失败"
            md += f"### {result['description']}\n\n"
            md += f"- **状态**: {status}\n"
            md += f"- **执行时间**: {result['execution_time']:.1f} 秒\n"
            md += f"- **命令**: `{result['command']}`\n"
            
            if not result['success']:
                if 'error' in result:
                    md += f"- **错误**: {result['error']}\n"
                if 'stderr' in result and result['stderr']:
                    md += f"- **错误输出**:\n```\n{result['stderr'][:500]}...\n```\n"
            
            md += "\n"
        
        # 模块可用性
        md += "## 模块可用性\n\n"
        for module, available in results['availability'].items():
            status = "✅ 可用" if available else "❌ 不可用"
            config = self.test_modules.get(module, {})
            optional = " (可选)" if config.get('optional', False) else ""
            md += f"- **{config.get('description', module)}**: {status}{optional}\n"
        
        return md
    
    def print_summary(self, results: Dict[str, Any]):
        """打印测试摘要"""
        print("\n" + "="*60)
        print("ChatSpatial MCP 全面测试完成")
        print("="*60)
        print(f"总耗时: {results['total_execution_time']:.1f} 秒")
        print(f"成功率: {results['success_rate']:.1%} ({results['successful_modules']}/{results['total_modules']})")
        print(f"总体状态: {'✅ 通过' if results['overall_success'] else '❌ 失败'}")
        print()
        
        print("各模块结果:")
        for result in results['module_results']:
            status = "✅" if result['success'] else "❌"
            print(f"  {status} {result['description']} ({result['execution_time']:.1f}s)")
        
        print("\n模块可用性:")
        for module, available in results['availability'].items():
            config = self.test_modules.get(module, {})
            status = "✅" if available else "❌"
            optional = " (可选)" if config.get('optional', False) else ""
            print(f"  {status} {config.get('description', module)}{optional}")

def main():
    parser = argparse.ArgumentParser(description='ChatSpatial MCP 全面测试执行器')
    parser.add_argument('--modules', nargs='+', 
                       choices=['mcp_protocol', 'tool_functionality', 'integration_workflow', 
                               'performance_stress', 'error_compatibility'],
                       help='选择要运行的测试模块')
    parser.add_argument('--output', help='输出报告文件名')
    parser.add_argument('--quiet', action='store_true', help='静默模式')
    
    args = parser.parse_args()
    
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    
    # 初始化测试运行器
    base_dir = Path(__file__).parent
    runner = ComprehensiveTestRunner(base_dir)
    
    try:
        # 运行测试
        results = runner.run_all_tests(selected_modules=args.modules)
        
        # 生成报告
        report_path = runner.generate_report(results, args.output)
        
        # 打印摘要
        if not args.quiet:
            runner.print_summary(results)
        
        # 返回适当的退出码
        sys.exit(0 if results['overall_success'] else 1)
        
    except KeyboardInterrupt:
        logger.error("测试被用户中断")
        sys.exit(130)
    except Exception as e:
        logger.error(f"测试执行异常: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()