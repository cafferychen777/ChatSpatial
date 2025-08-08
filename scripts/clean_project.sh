#!/bin/bash
# ChatSpatial project cleanup script

echo "ChatSpatial Project Cleanup Tool"
echo "================================"

# Clean Python cache
echo "Cleaning Python cache..."
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find . -type f -name "*.pyc" -delete 2>/dev/null
find . -type f -name "*.pyo" -delete 2>/dev/null

# Clean temporary files
echo "Cleaning temporary files..."
find . -type f -name ".DS_Store" -delete 2>/dev/null
find . -type f -name "*.log" -delete 2>/dev/null
find . -type f -name "*.tmp" -delete 2>/dev/null

# Clean test-generated files (optional)
read -p "Do you want to clean test-generated visualization files? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Cleaning visualization files..."
    rm -rf visualization_resources/*.png 2>/dev/null
    rm -rf visualization_resources/*.html 2>/dev/null
fi

# Clean egg-info (optional)
read -p "Do you want to clean egg-info directories? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Cleaning egg-info..."
    rm -rf *.egg-info 2>/dev/null
fi

echo "Cleanup completed!"

# Show directory sizes
echo -e "\nDirectory size statistics:"
du -sh . 2>/dev/null
du -sh data/ 2>/dev/null
du -sh third_party/ 2>/dev/null
du -sh tests/ 2>/dev/null
du -sh docs/ 2>/dev/null