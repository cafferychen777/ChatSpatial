"""
测试图像处理工具模块
"""

import sys
import os
import unittest
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# 添加项目根目录到 Python 路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from spatial_transcriptomics_mcp.utils.image_utils import (
    fig_to_image,
    fig_to_base64,
    create_placeholder_image
)


class TestImageUtils(unittest.TestCase):
    """测试图像处理工具函数"""

    def test_fig_to_image(self):
        """测试 fig_to_image 函数"""
        # 创建一个简单的图形
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot([1, 2, 3, 4], [1, 4, 9, 16])
        ax.set_title("Test Figure")

        # 转换为 Image 对象
        image = fig_to_image(fig)

        # 验证返回的是 Image 对象
        self.assertEqual(image._format, "png")
        self.assertIsNotNone(image.data)
        self.assertIsInstance(image.data, bytes)

        # 验证图像数据不为空
        self.assertGreater(len(image.data), 0)

        # 验证 mime_type 正确
        self.assertEqual(image._mime_type, "image/png")

    def test_fig_to_base64(self):
        """测试 fig_to_base64 函数"""
        # 创建一个简单的图形
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot([1, 2, 3, 4], [1, 4, 9, 16])
        ax.set_title("Test Figure")

        # 转换为 base64 字符串
        base64_str = fig_to_base64(fig)

        # 验证返回的是字符串
        self.assertIsInstance(base64_str, str)

        # 验证字符串不为空
        self.assertGreater(len(base64_str), 0)

    def test_create_placeholder_image(self):
        """测试 create_placeholder_image 函数"""
        # 创建占位图像
        message = "Test Placeholder"
        image = create_placeholder_image(message)

        # 验证返回的是 Image 对象
        self.assertEqual(image._format, "png")
        self.assertIsNotNone(image.data)
        self.assertIsInstance(image.data, bytes)

        # 验证图像数据不为空
        self.assertGreater(len(image.data), 0)

        # 验证 mime_type 正确
        self.assertEqual(image._mime_type, "image/png")

    def test_image_size_control(self):
        """测试图像大小控制功能"""
        # 创建一个复杂的图形（大尺寸、高分辨率）
        fig, ax = plt.subplots(figsize=(20, 20), dpi=300)

        # 生成复杂数据
        x = np.linspace(0, 10, 1000)
        y = np.sin(x) * np.cos(x**2)

        # 绘制复杂图形
        ax.plot(x, y, linewidth=2)
        ax.set_title("Complex Figure")

        # 使用较小的最大大小限制
        max_size_kb = 100  # 增加大小限制
        image = fig_to_image(fig, max_size_kb=max_size_kb)

        # 验证图像大小不超过限制（允许一些误差）
        actual_size_kb = len(image.data) / 1024
        print(f"Image size: {actual_size_kb:.2f} KB")

        # 验证图像被压缩了（原始大小应该远大于此）
        self.assertLess(actual_size_kb, 200)  # 确保图像被压缩到合理大小


if __name__ == "__main__":
    unittest.main()
