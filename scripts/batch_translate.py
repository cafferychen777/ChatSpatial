#!/usr/bin/env python3
"""
Batch translation script for Chinese to English in Python files
"""

import os
import re
import sys
from pathlib import Path

# Translation dictionary for common Chinese terms
TRANSLATIONS = {
    # Comments and docstrings
    "Test": "Test",
    "Test用例": "Test case", 
    "Testdata": "Test data",
    "Testresult": "Test result",
    "Testfunction": "Test function",
    "Testmethod": "Test method",
    "Testclass": "Test class",
    "单元Test": "Unit test",
    "集成Test": "Integration test",
    "功能Test": "Functional test",
    
    # Data and analysis
    "data": "data",
    "analysis": "analysis", 
    "visualization": "visualization",
    "processing": "processing",
    "预processing": "preprocessing",
    "后processing": "postprocessing",
    "statistics": "statistics",
    "statisticsinformation": "statistics",
    "result": "result",
    "output": "output",
    "input": "input",
    
    # Parameters and validation
    "parameter": "parameter",
    "parametervalidation": "parameter validation",
    "validation": "validation",
    "validation": "validation",
    "check": "check",
    "verification": "verification",
    
    # Errors and exceptions
    "error": "error",
    "exception": "exception",
    "failure": "failure",
    "success": "success",
    "warning": "warning",
    "information": "information",
    "message": "message",
    
    # Functions and methods
    "function": "function",
    "method": "method",
    "class": "class",
    "module": "module",
    "tool": "tool",
    "实用tool": "utility",
    "帮助function": "helper function",
    
    # File operations
    "file": "file",
    "path": "path",
    "directory": "directory",
    "file夹": "folder",
    "load": "load",
    "save": "save",
    "read": "read",
    "write": "write",
    
    # Types and formats
    "class型": "type",
    "format": "format",
    "object": "object",
    "instance": "instance",
    "attribute": "attribute",
    "field": "field",
    
    # Common phrases
    "create": "create",
    "generate": "generate",
    "build": "build",
    "initialize": "initialize",
    "setup": "setup",
    "configuration": "configuration",
    "option": "option",
    "setting": "setting",
    
    # Status and states
    "completed": "completed",
    "start": "start",
    "end": "end",
    "run": "run",
    "execute": "execute",
    "call": "call",
    "return": "return",
    
    # Spatial analysis specific
    "spatial": "spatial",
    "spatialanalysis": "spatial analysis",
    "spatialdata": "spatial data",
    "spatialcoordinate": "spatial coordinates",
    "neighborhood": "neighborhood",
    "聚class": "clustering",
    "gene": "gene",
    "cell": "cell",
    "tissue": "tissue",
    "sample": "sample",
    
    # Common sentence patterns
    "添加项目根directory到 Python path": "Add project root directory to Python path",
    "Mock MCP context": "Mock MCP context",
    "createTestdata": "Create test data",
    "create一个简单的": "Create a simple",
    "用于Test": "for testing",
    "Testcompleted": "Testing completed",
    "Testsuccess": "Test successful",
    "Testfailure": "Test failed",
}

def translate_file(file_path):
    """Translate Chinese text in a single file"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        original_content = content
        
        # Apply translations
        for chinese, english in TRANSLATIONS.items():
            # Replace exact matches (word boundaries)
            content = re.sub(rf'\b{re.escape(chinese)}\b', english, content)
            # Also replace in comments and strings
            content = content.replace(chinese, english)
        
        # Only write if content changed
        if content != original_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
            return True
        return False
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False

def main():
    """Main function to batch translate files"""
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    else:
        target_dir = "."
    
    # Find all Python files
    python_files = []
    for root, dirs, files in os.walk(target_dir):
        # Skip certain directories
        if any(skip in root for skip in ['third_party', 'node_modules', 'drawio-mcp-server', '.git']):
            continue
            
        for file in files:
            if file.endswith('.py'):
                python_files.append(os.path.join(root, file))
    
    print(f"Found {len(python_files)} Python files to process...")
    
    translated_count = 0
    for file_path in python_files:
        if translate_file(file_path):
            print(f"Translated: {file_path}")
            translated_count += 1
    
    print(f"\nTranslation completed! {translated_count} files were modified.")

if __name__ == "__main__":
    main()
