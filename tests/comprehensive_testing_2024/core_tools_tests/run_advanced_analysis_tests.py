#!/usr/bin/env python3
"""
ChatSpatial 高级分析工具测试运行器

这个脚本按照Linus的哲学运行完整的测试套件：
1. "解决真实问题" - 测试实际使用场景
2. "消除特殊情况" - 统一的测试框架
3. "Never break userspace" - 验证所有fallback机制
4. "好品味的代码" - 简洁有效的测试逻辑

用法:
    python run_advanced_analysis_tests.py [options]

选项:
    --quick: 快速测试（跳过耗时测试）
    --full: 完整测试（包括性能基准）
    --report-only: 仅生成依赖报告
    --output-dir: 报告输出目录
"""

import sys
import time
import argparse
import logging
import asyncio
import json
from pathlib import Path
from typing import Dict, Any, Optional, List

# 设置项目路径
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

# 导入测试组件
from test_advanced_analysis_tools import AdvancedAnalysisTestFramework
from dependency_validator import DependencyValidator

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('advanced_analysis_test.log')
    ]
)

logger = logging.getLogger(__name__)


class AdvancedAnalysisTestRunner:
    """
    高级分析工具测试运行器
    
    Linus式设计原则：
    - 单一职责：协调测试运行
    - 失败快速：遇到致命错误立即停止
    - 清晰反馈：实时进度和结果报告
    - 实用主义：专注解决实际问题
    """
    
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.start_time = None
        self.results = {}
        
    def run_dependency_check(self) -> Dict[str, Any]:
        """
        运行依赖检查
        
        返回依赖报告，如果关键依赖缺失则警告但不停止
        （遵循实用主义：尽可能多地测试可用功能）
        """
        logger.info("🔍 Step 1: Dependency validation...")
        
        validator = DependencyValidator()
        dependency_report = validator.generate_dependency_report()
        
        # 保存依赖报告
        json_path = self.output_dir / 'dependency_report.json'
        validator.save_report_json(json_path)
        
        md_path = self.output_dir / 'dependency_report.md'
        validator.save_report_markdown(md_path)
        
        # 分析依赖状态
        summary = dependency_report['summary']
        critical_missing = summary['critical_missing']
        
        if critical_missing:
            logger.warning(f"⚠️  Critical dependencies missing: {', '.join(critical_missing)}")
            logger.info("Will test available functionality with fallback mechanisms")
        else:
            logger.info("✅ All critical dependencies available")
        
        logger.info(f"Dependency check completed: {summary['functionality_rate']:.1%} functional")
        
        return dependency_report
    
    async def run_advanced_analysis_tests(self, quick_mode: bool = False) -> Dict[str, Any]:
        """
        运行高级分析工具测试
        
        参数:
            quick_mode: 如果为True，跳过耗时的测试
        """
        logger.info("🧪 Step 2: Advanced analysis tools testing...")
        
        framework = AdvancedAnalysisTestFramework()
        
        # 设置测试数据
        if not framework.setup_test_data():
            raise RuntimeError("Failed to setup test data - cannot continue")
        
        logger.info(f"Test data ready: {len(framework.data_store)} datasets")
        
        # 检查依赖状态
        framework.check_dependencies()
        
        # 运行测试模块
        test_results = {}
        
        try:
            # 轨迹分析测试
            logger.info("Testing trajectory analysis...")
            trajectory_results = await framework.test_trajectory_analysis()
            test_results['trajectory_analysis'] = trajectory_results
            
        except Exception as e:
            logger.error(f"Trajectory analysis testing failed: {e}")
            test_results['trajectory_analysis'] = {'error': str(e)}
        
        try:
            # 通路富集分析测试
            logger.info("Testing pathway enrichment...")
            pathway_results = await framework.test_pathway_enrichment()
            test_results['pathway_enrichment'] = pathway_results
            
        except Exception as e:
            logger.error(f"Pathway enrichment testing failed: {e}")
            test_results['pathway_enrichment'] = {'error': str(e)}
        
        try:
            # 空间富集分析测试
            logger.info("Testing spatial enrichment...")
            spatial_enrichment_results = await framework.test_spatial_enrichment()
            test_results['spatial_enrichment'] = spatial_enrichment_results
            
        except Exception as e:
            logger.error(f"Spatial enrichment testing failed: {e}")
            test_results['spatial_enrichment'] = {'error': str(e)}
        
        try:
            # 空间配准测试
            logger.info("Testing spatial registration...")
            registration_results = await framework.test_spatial_registration()
            test_results['spatial_registration'] = registration_results
            
        except Exception as e:
            logger.error(f"Spatial registration testing failed: {e}")
            test_results['spatial_registration'] = {'error': str(e)}
        
        # 生成性能基准
        performance_benchmarks = framework.generate_performance_benchmarks()
        
        return {
            'test_results': test_results,
            'performance_benchmarks': performance_benchmarks,
            'dependency_status': framework.dependency_status
        }
    
    def generate_comprehensive_report(
        self, 
        dependency_report: Dict[str, Any], 
        analysis_results: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        生成综合测试报告
        
        整合依赖检查和功能测试结果
        """
        logger.info("📊 Step 3: Generating comprehensive report...")
        
        total_time = time.time() - self.start_time if self.start_time else 0
        
        # 统计测试成功率
        test_results = analysis_results.get('test_results', {})
        success_count = 0
        total_count = 0
        
        for tool_name, tool_results in test_results.items():
            for test_name, test_result in tool_results.items():
                if isinstance(test_result, dict) and 'success' in test_result:
                    total_count += 1
                    if test_result.get('success', False):
                        success_count += 1
        
        success_rate = success_count / total_count if total_count > 0 else 0
        
        # 分析性能数据
        performance_data = analysis_results.get('performance_benchmarks', {})
        performance_summary = {}
        
        if 'summary' in performance_data:
            performance_summary = {
                'total_operations': performance_data['summary'].get('n_operations', 0),
                'total_test_time': performance_data['summary'].get('total_test_time', 0),
                'mean_operation_time': performance_data['summary'].get('mean_operation_time', 0),
                'slowest_operation': performance_data['summary'].get('slowest_operation_time', 0)
            }
        
        # 生成最终建议
        recommendations = self._generate_final_recommendations(
            dependency_report, analysis_results, success_rate
        )
        
        comprehensive_report = {
            'metadata': {
                'test_date': time.strftime('%Y-%m-%d %H:%M:%S'),
                'total_runtime': total_time,
                'python_version': sys.version.split()[0],
                'test_framework_version': '1.0.0'
            },
            'executive_summary': {
                'overall_status': self._get_overall_status(dependency_report, success_rate),
                'dependency_health': dependency_report['summary']['functionality_rate'],
                'test_success_rate': success_rate,
                'critical_issues': self._identify_critical_issues(dependency_report, analysis_results),
                'ready_for_production': success_rate > 0.8 and dependency_report['summary']['functionality_rate'] > 0.7
            },
            'dependency_analysis': dependency_report,
            'functional_testing': analysis_results,
            'performance_analysis': performance_summary,
            'recommendations': recommendations,
            'detailed_results': test_results
        }
        
        return comprehensive_report
    
    def _get_overall_status(self, dependency_report: Dict[str, Any], success_rate: float) -> str:
        """确定总体状态"""
        dep_rate = dependency_report['summary']['functionality_rate']
        
        if success_rate > 0.9 and dep_rate > 0.9:
            return "EXCELLENT"
        elif success_rate > 0.7 and dep_rate > 0.7:
            return "GOOD"
        elif success_rate > 0.5 and dep_rate > 0.5:
            return "FAIR"
        elif success_rate > 0.3 or dep_rate > 0.3:
            return "POOR"
        else:
            return "CRITICAL"
    
    def _identify_critical_issues(
        self, 
        dependency_report: Dict[str, Any], 
        analysis_results: Dict[str, Any]
    ) -> List[str]:
        """识别关键问题"""
        issues = []
        
        # 关键依赖缺失
        critical_missing = dependency_report['summary'].get('critical_missing', [])
        if critical_missing:
            issues.append(f"Critical dependencies missing: {', '.join(critical_missing)}")
        
        # 功能测试失败
        test_results = analysis_results.get('test_results', {})
        failed_tools = []
        
        for tool_name, tool_results in test_results.items():
            if 'error' in tool_results:
                failed_tools.append(tool_name)
            else:
                # 检查各个测试的失败情况
                failed_tests = []
                for test_name, test_result in tool_results.items():
                    if isinstance(test_result, dict) and not test_result.get('success', True):
                        failed_tests.append(test_name)
                
                if failed_tests:
                    issues.append(f"{tool_name} failed tests: {', '.join(failed_tests)}")
        
        if failed_tools:
            issues.append(f"Complete tool failures: {', '.join(failed_tools)}")
        
        return issues
    
    def _generate_final_recommendations(
        self,
        dependency_report: Dict[str, Any],
        analysis_results: Dict[str, Any],
        success_rate: float
    ) -> List[str]:
        """生成最终建议"""
        recommendations = []
        
        # 基于依赖状态的建议
        critical_missing = dependency_report['summary'].get('critical_missing', [])
        if critical_missing:
            rec_text = "IMMEDIATE: Install critical dependencies: "
            install_commands = []
            
            for dep in critical_missing:
                detailed_results = dependency_report.get('detailed_results', {})
                if dep in detailed_results and detailed_results[dep].get('installation_command'):
                    install_commands.append(detailed_results[dep]['installation_command'])
            
            if install_commands:
                rec_text += "; ".join(install_commands)
            else:
                rec_text += ", ".join(critical_missing)
            
            recommendations.append(rec_text)
        
        # 基于测试成功率的建议
        if success_rate < 0.5:
            recommendations.append("URGENT: Major system issues detected - review all installations")
        elif success_rate < 0.8:
            recommendations.append("MEDIUM: Some tools non-functional - consider installing optional dependencies")
        
        # 基于性能的建议
        performance_data = analysis_results.get('performance_benchmarks', {})
        if 'summary' in performance_data:
            total_time = performance_data['summary'].get('total_test_time', 0)
            if total_time > 300:  # 5分钟
                recommendations.append("PERFORMANCE: Consider GPU acceleration for large-scale analysis")
            elif total_time > 120:  # 2分钟
                recommendations.append("INFO: Performance acceptable for interactive use")
        
        # 具体工具建议
        test_results = analysis_results.get('test_results', {})
        
        # 轨迹分析建议
        if 'trajectory_analysis' in test_results:
            traj_results = test_results['trajectory_analysis']
            if 'trajectory_inference' in traj_results:
                inference_results = traj_results['trajectory_inference']
                working_methods = [method for method, result in inference_results.items() 
                                 if isinstance(result, dict) and result.get('success', False)]
                
                if not working_methods:
                    recommendations.append("WARNING: No trajectory inference methods working")
                elif len(working_methods) == 1 and working_methods[0] == 'dpt':
                    recommendations.append("INFO: Only basic DPT available - consider installing CellRank or Palantir")
        
        # 如果所有测试都很好
        if success_rate > 0.9 and dependency_report['summary']['functionality_rate'] > 0.9:
            recommendations.append("SUCCESS: All systems functional - ready for production use")
        
        return recommendations
    
    def save_reports(self, comprehensive_report: Dict[str, Any]) -> None:
        """保存所有报告文件"""
        logger.info("💾 Step 4: Saving reports...")
        
        # JSON格式详细报告
        json_path = self.output_dir / 'advanced_analysis_test_results.json'
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(comprehensive_report, f, indent=2, ensure_ascii=False, default=str)
        
        logger.info(f"Detailed JSON report saved: {json_path}")
        
        # Markdown格式摘要报告
        md_path = self.output_dir / 'advanced_analysis_test_summary.md'
        self._save_markdown_summary(comprehensive_report, md_path)
        
        logger.info(f"Summary Markdown report saved: {md_path}")
        
        # 如果有性能数据，保存性能报告
        if comprehensive_report.get('performance_analysis'):
            perf_path = self.output_dir / 'performance_benchmark.json'
            with open(perf_path, 'w', encoding='utf-8') as f:
                json.dump(comprehensive_report['performance_analysis'], f, indent=2)
            
            logger.info(f"Performance benchmark saved: {perf_path}")
    
    def _save_markdown_summary(self, report: Dict[str, Any], filepath: Path) -> None:
        """保存Markdown格式摘要"""
        summary = report['executive_summary']
        metadata = report['metadata']
        
        md_content = f"""# ChatSpatial 高级分析工具测试摘要

**测试日期**: {metadata['test_date']}  
**运行时间**: {metadata['total_runtime']:.2f}秒  
**Python版本**: {metadata['python_version']}  

## 总体状态: {summary['overall_status']}

- **依赖健康度**: {summary['dependency_health']:.1%}
- **测试成功率**: {summary['test_success_rate']:.1%}
- **生产就绪**: {'✅ 是' if summary['ready_for_production'] else '❌ 否'}

## 关键发现

"""
        
        if summary['critical_issues']:
            md_content += "### 🔴 关键问题\n\n"
            for issue in summary['critical_issues']:
                md_content += f"- {issue}\n"
            md_content += "\n"
        
        # 工具状态摘要
        md_content += "### 工具模块状态\n\n"
        test_results = report.get('detailed_results', {})
        
        for tool_name, tool_results in test_results.items():
            if 'error' in tool_results:
                status_icon = "❌"
                status_text = f"失败: {tool_results['error']}"
            else:
                # 计算这个工具的成功率
                success_count = 0
                total_count = 0
                
                for test_name, test_result in tool_results.items():
                    if isinstance(test_result, dict) and 'success' in test_result:
                        total_count += 1
                        if test_result.get('success', False):
                            success_count += 1
                
                if total_count > 0:
                    success_rate = success_count / total_count
                    if success_rate >= 1.0:
                        status_icon = "✅"
                        status_text = "完全正常"
                    elif success_rate >= 0.5:
                        status_icon = "⚠️"
                        status_text = f"部分正常 ({success_count}/{total_count})"
                    else:
                        status_icon = "❌"
                        status_text = f"大部分失败 ({success_count}/{total_count})"
                else:
                    status_icon = "❓"
                    status_text = "未测试"
            
            tool_display_name = tool_name.replace('_', ' ').title()
            md_content += f"- {status_icon} **{tool_display_name}**: {status_text}\n"
        
        # 性能摘要
        if 'performance_analysis' in report and report['performance_analysis']:
            perf = report['performance_analysis']
            md_content += f"\n### 性能概览\n\n"
            md_content += f"- **总操作数**: {perf.get('total_operations', 0)}\n"
            md_content += f"- **测试总时间**: {perf.get('total_test_time', 0):.2f}秒\n"
            md_content += f"- **平均操作时间**: {perf.get('mean_operation_time', 0):.3f}秒\n"
            md_content += f"- **最慢操作时间**: {perf.get('slowest_operation', 0):.3f}秒\n"
        
        # 建议
        if report.get('recommendations'):
            md_content += "\n## 建议\n\n"
            for i, rec in enumerate(report['recommendations'], 1):
                priority = "🔴" if rec.startswith(('IMMEDIATE', 'URGENT', 'CRITICAL')) else \
                          "🟡" if rec.startswith(('MEDIUM', 'WARNING')) else \
                          "ℹ️"
                md_content += f"{i}. {priority} {rec}\n"
        
        md_content += f"\n---\n\n**详细报告**: `advanced_analysis_test_results.json`\n"
        md_content += f"**依赖报告**: `dependency_report.md`\n"
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(md_content)
    
    async def run_complete_test_suite(self, quick_mode: bool = False) -> Dict[str, Any]:
        """
        运行完整测试套件
        
        这是主要的协调函数，按照Linus的原则：
        1. 失败时快速停止
        2. 清晰的进度反馈
        3. 全面的错误记录
        """
        self.start_time = time.time()
        
        logger.info("🚀 Starting ChatSpatial Advanced Analysis Tools Test Suite")
        logger.info(f"Output directory: {self.output_dir}")
        
        try:
            # 步骤1：依赖检查
            dependency_report = self.run_dependency_check()
            
            # 步骤2：功能测试
            analysis_results = await self.run_advanced_analysis_tests(quick_mode)
            
            # 步骤3：生成综合报告
            comprehensive_report = self.generate_comprehensive_report(
                dependency_report, analysis_results
            )
            
            # 步骤4：保存报告
            self.save_reports(comprehensive_report)
            
            # 最终状态
            total_time = time.time() - self.start_time
            status = comprehensive_report['executive_summary']['overall_status']
            
            logger.info(f"🎯 Test suite completed in {total_time:.2f}s")
            logger.info(f"Overall status: {status}")
            
            if comprehensive_report['executive_summary']['ready_for_production']:
                logger.info("✅ System ready for production use")
            else:
                logger.warning("⚠️ System needs attention before production use")
            
            return comprehensive_report
            
        except Exception as e:
            logger.error(f"💥 Test suite failed: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            
            # 保存错误报告
            error_report = {
                'error': str(e),
                'traceback': traceback.format_exc(),
                'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime': time.time() - self.start_time if self.start_time else 0
            }
            
            error_path = self.output_dir / 'test_suite_error.json'
            with open(error_path, 'w') as f:
                json.dump(error_report, f, indent=2)
            
            raise
    
    def print_summary(self, report: Dict[str, Any]) -> None:
        """打印测试摘要到控制台"""
        print("\n" + "="*80)
        print("ChatSpatial 高级分析工具测试摘要")
        print("="*80)
        
        summary = report['executive_summary']
        metadata = report['metadata']
        
        # 总体状态
        status_colors = {
            'EXCELLENT': '🟢',
            'GOOD': '🟢',
            'FAIR': '🟡',
            'POOR': '🟠',
            'CRITICAL': '🔴'
        }
        
        status_icon = status_colors.get(summary['overall_status'], '❓')
        print(f"总体状态: {status_icon} {summary['overall_status']}")
        print(f"测试时间: {metadata['total_runtime']:.2f}秒")
        print(f"依赖健康度: {summary['dependency_health']:.1%}")
        print(f"测试成功率: {summary['test_success_rate']:.1%}")
        print(f"生产就绪: {'✅' if summary['ready_for_production'] else '❌'}")
        
        # 关键问题
        if summary['critical_issues']:
            print(f"\n🔴 关键问题:")
            for issue in summary['critical_issues']:
                print(f"  • {issue}")
        
        # 建议
        if report.get('recommendations'):
            print(f"\n📋 主要建议:")
            for i, rec in enumerate(report['recommendations'][:3], 1):  # 只显示前3个
                print(f"  {i}. {rec}")
            
            if len(report['recommendations']) > 3:
                print(f"  ... 还有 {len(report['recommendations']) - 3} 条建议（见详细报告）")
        
        print(f"\n📊 详细报告已保存到: {self.output_dir}")
        print("="*80)


async def main():
    """主程序入口"""
    parser = argparse.ArgumentParser(
        description="ChatSpatial Advanced Analysis Tools Test Suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_advanced_analysis_tests.py --full
  python run_advanced_analysis_tests.py --quick --output-dir ./test_results
  python run_advanced_analysis_tests.py --report-only
        """
    )
    
    parser.add_argument('--quick', action='store_true', 
                       help='Run quick tests only (skip time-consuming benchmarks)')
    parser.add_argument('--full', action='store_true',
                       help='Run full test suite including performance benchmarks')
    parser.add_argument('--report-only', action='store_true',
                       help='Generate dependency report only (no functional testing)')
    parser.add_argument('--output-dir', type=Path, 
                       default=Path('./advanced_analysis_test_results'),
                       help='Output directory for test reports (default: ./advanced_analysis_test_results)')
    
    args = parser.parse_args()
    
    # 创建测试运行器
    runner = AdvancedAnalysisTestRunner(args.output_dir)
    
    try:
        if args.report_only:
            # 仅生成依赖报告
            logger.info("Running dependency check only...")
            dependency_report = runner.run_dependency_check()
            
            print("\n" + "="*60)
            print("DEPENDENCY CHECK SUMMARY")
            print("="*60)
            
            summary = dependency_report['summary']
            print(f"Dependencies: {summary['functional_dependencies']}/{summary['total_dependencies']} functional")
            
            if summary['critical_missing']:
                print(f"Critical missing: {', '.join(summary['critical_missing'])}")
            
            print(f"Report saved to: {args.output_dir}")
            print("="*60)
            
        else:
            # 运行完整测试
            quick_mode = args.quick or not args.full
            
            if quick_mode:
                logger.info("Running in quick mode (skipping intensive benchmarks)")
            else:
                logger.info("Running full test suite including performance benchmarks")
            
            comprehensive_report = await runner.run_complete_test_suite(quick_mode)
            runner.print_summary(comprehensive_report)
            
            # 返回适当的退出代码
            if comprehensive_report['executive_summary']['ready_for_production']:
                return 0
            elif comprehensive_report['executive_summary']['overall_status'] in ['POOR', 'CRITICAL']:
                return 2
            else:
                return 1
    
    except KeyboardInterrupt:
        logger.info("Test interrupted by user")
        return 130
    
    except Exception as e:
        logger.error(f"Test suite failed: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(asyncio.run(main()))