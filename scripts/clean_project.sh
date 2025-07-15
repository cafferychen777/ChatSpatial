#!/bin/bash
# ChatSpatial 项目清理脚本

echo "ChatSpatial 项目清理工具"
echo "========================"

# 清理Python缓存
echo "清理Python缓存..."
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find . -type f -name "*.pyc" -delete 2>/dev/null
find . -type f -name "*.pyo" -delete 2>/dev/null

# 清理临时文件
echo "清理临时文件..."
find . -type f -name ".DS_Store" -delete 2>/dev/null
find . -type f -name "*.log" -delete 2>/dev/null
find . -type f -name "*.tmp" -delete 2>/dev/null

# 清理测试生成的文件（可选）
read -p "是否清理测试生成的可视化文件? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "清理可视化文件..."
    rm -rf visualization_resources/*.png 2>/dev/null
    rm -rf visualization_resources/*.html 2>/dev/null
fi

# 清理egg-info（可选）
read -p "是否清理egg-info目录? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "清理egg-info..."
    rm -rf *.egg-info 2>/dev/null
fi

echo "清理完成！"

# 显示目录大小
echo -e "\n目录大小统计："
du -sh . 2>/dev/null
du -sh data/ 2>/dev/null
du -sh third_party/ 2>/dev/null
du -sh tests/ 2>/dev/null
du -sh docs/ 2>/dev/null