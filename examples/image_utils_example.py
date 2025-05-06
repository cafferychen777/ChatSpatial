"""
图像处理模块使用示例
"""

import matplotlib.pyplot as plt
import numpy as np
from mcp.server.fastmcp.utilities.types import Image

# 导入图像处理模块
from spatial_transcriptomics_mcp.utils.image_utils import (
    fig_to_image,
    fig_to_base64,
    create_placeholder_image
)


def create_sample_figure():
    """创建一个示例图形"""
    # 创建数据
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x, y, 'b-', linewidth=2)
    ax.set_title("Sample Sine Wave")
    ax.set_xlabel("X")
    ax.set_ylabel("sin(x)")
    ax.grid(True)
    
    return fig


def example_fig_to_image():
    """fig_to_image 函数示例"""
    fig = create_sample_figure()
    
    # 转换为 Image 对象
    image = fig_to_image(fig, dpi=100, format='png')
    
    print(f"Image object created: format={image._format}, size={len(image.data)/1024:.2f} KB")
    return image


def example_create_placeholder_image():
    """create_placeholder_image 函数示例"""
    # 创建占位图像
    image = create_placeholder_image(
        message="This is a placeholder image",
        figsize=(6, 6),
        format='png'
    )
    
    print(f"Placeholder image created: format={image._format}, size={len(image.data)/1024:.2f} KB")
    return image


def example_fig_to_base64():
    """fig_to_base64 函数示例（向后兼容）"""
    fig = create_sample_figure()
    
    # 转换为 base64 字符串
    base64_str = fig_to_base64(fig, dpi=100, format='png')
    
    print(f"Base64 string created: length={len(base64_str)} characters")
    return base64_str


def main():
    """运行所有示例"""
    print("Running fig_to_image example...")
    example_fig_to_image()
    
    print("\nRunning create_placeholder_image example...")
    example_create_placeholder_image()
    
    print("\nRunning fig_to_base64 example...")
    example_fig_to_base64()
    
    print("\nAll examples completed successfully!")


if __name__ == "__main__":
    main()
