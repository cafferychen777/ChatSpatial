"""
ChatSpatial MCP 协议测试包

本包包含ChatSpatial MCP服务器协议层的全面测试套件。

测试模块：
- test_server_startup.py: 服务器启动和传输模式测试
- test_tool_registration.py: 工具注册和元数据验证  
- test_parameter_validation.py: 参数验证和边界检查
- test_error_responses.py: 错误处理和JSON-RPC兼容性
- test_http_transport.py: HTTP传输和并发处理
- run_protocol_tests.py: 测试运行器和报告生成器

作者: Linus 风格的务实测试架构
"""

__version__ = "1.0.0"
__author__ = "ChatSpatial Team"

# 测试常量
EXPECTED_TOOL_COUNT = 16  # 当前ChatSpatial的15个工具（可能会增加）
MCP_PROTOCOL_VERSION = "2024-11-05"
JSON_RPC_VERSION = "2.0"

# 测试配置
TEST_CONFIG = {
    "server_startup": {
        "max_startup_time": 5.0,
        "performance_iterations": 5,
        "memory_limit_mb": 100
    },
    "tool_registration": {
        "expected_tools": 15,
        "discovery_timeout": 0.1,
        "metadata_completeness": True
    },
    "parameter_validation": {
        "performance_threshold": 0.001,
        "error_quality_min": 70.0
    },
    "error_responses": {
        "message_quality_min": 70.0,
        "recovery_rate_min": 50.0,
        "performance_overhead_max": 50.0
    },
    "http_transport": {
        "request_timeout": 0.1,
        "throughput_min": 10.0,
        "concurrent_success_min": 80.0
    }
}

# 导出测试类
from .test_server_startup import TestServerStartup
from .test_tool_registration import TestToolRegistration
from .test_parameter_validation import TestParameterValidation
from .test_error_responses import TestErrorResponses
from .test_http_transport import TestHTTPTransport

__all__ = [
    "TestServerStartup",
    "TestToolRegistration", 
    "TestParameterValidation",
    "TestErrorResponses",
    "TestHTTPTransport",
    "TEST_CONFIG"
]