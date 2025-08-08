#!/usr/bin/env python3
"""
Enhanced translation script for remaining Chinese text
"""

import os
import re
import sys
from pathlib import Path

# Extended translation dictionary
TRANSLATIONS = {
    # Basic terms that might have been missed
    "image": "image",
    "image": "image", 
    "chart": "chart",
    "plot": "plot",
    "display": "display",
    "show": "show",
    "size": "size",
    "dimension": "dimension",
    "width": "width", 
    "height": "height",
    "color": "color",
    "title": "title",
    "label": "label",
    "coordinate": "coordinate",
    
    # More specific terms
    "saved to": "saved to",
    "save to": "save to",
    "saved": "saved",
    "save": "save",
    "load data": "load data",
    "data loading": "data loading",
    "process data": "process data",
    "data processing": "data processing",
    
    # Status messages
    "processing": "processing",
    "processing": "processing",
    "loading": "loading",
    "loading": "loading",
    "正在save": "saving",
    "save中": "saving",
    "analyzing": "analyzing",
    "analyzing": "analyzing",
    
    # Results and outputs
    "analysis result": "analysis result",
    "processing result": "processing result",
    "output result": "output result",
    "final result": "final result",
    "intermediate result": "intermediate result",
    
    # File operations
    "file path": "file path",
    "filename": "filename",
    "file format": "file format",
    "file type": "file type",
    "文件size": "file size",
    
    # More analysis terms
    "spatial genes": "spatial genes",
    "differential genes": "differential genes",
    "highly variable genes": "highly variable genes",
    "marker genes": "marker genes",
    "expression level": "expression level",
    "expression": "expression",
    "基因expression": "gene expression",
    "cell type": "cell type",
    "cell annotation": "cell annotation",
    
    # Visualization terms
    "visualization result": "visualization result",
    "可视化image": "visualization image",
    "scatter plot": "scatter plot",
    "heatmap": "heatmap",
    "violin plot": "violin plot",
    "box plot": "box plot",
    "bar plot": "bar plot",
    "histogram": "histogram",
    
    # Common phrases in comments
    "if": "if",
    "else": "else",
    "when": "when",
    "then": "then",
    "because": "because",
    "so": "so",
    "but": "but",
    "however": "however",
    "therefore": "therefore",
    "additionally": "additionally",
    "furthermore": "furthermore",
    "finally": "finally",
    "first": "first",
    "second": "second",
    "again": "again",
    "meanwhile": "meanwhile",
    "similarly": "similarly",
    "conversely": "conversely",
    "for example": "for example",
    "such as": "such as",
    "including": "including",
    "except": "except",
    "unless": "unless",
    "until": "until",
    "although": "although",
    "despite": "despite",
    "regardless": "regardless",
    "no matter": "no matter",
    
    # More technical terms
    "algorithm": "algorithm",
    "model": "model",
    "training": "training",
    "prediction": "prediction",
    "classification": "classification",
    "regression": "regression",
    "optimization": "optimization",
    "parameter tuning": "parameter tuning",
    "hyperparameter": "hyperparameter",
    "cross validation": "cross validation",
    "validation set": "validation set",
    "test set": "test set",
    "training集": "training set",
    
    # Quality and performance
    "quality": "quality",
    "performance": "performance",
    "efficiency": "efficiency",
    "accuracy": "accuracy",
    "precision": "precision",
    "recall": "recall",
    "speed": "speed",
    "memory": "memory",
    "storage": "storage",
    "computation": "computation",
    "calculation": "calculation",
    
    # Common sentence patterns
    "if有": "if there is",
    "if没有": "if there is no",
    "if存在": "if exists",
    "if不存在": "if does not exist",
    "whether": "whether",
    "whether can": "whether can",
    "can": "can",
    "cannot": "cannot",
    "should": "should",
    "不should": "should not",
    "must": "must",
    "need not": "need not",
    "need": "need",
    "不need": "do not need",
    "require": "require",
    "recommend": "recommend",
    "recommend": "recommend",
    "support": "support",
    "不support": "not support",
    "allow": "allow",
    "不allow": "not allow",
    "prohibit": "prohibit",
    "enable": "enable",
    "disable": "disable",
    "turn on": "turn on",
    "turn off": "turn off",
}

def translate_file(file_path):
    """Translate Chinese text in a single file"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        original_content = content
        changes_made = 0
        
        # Apply translations
        for chinese, english in TRANSLATIONS.items():
            old_content = content
            # Replace exact matches
            content = content.replace(chinese, english)
            if content != old_content:
                changes_made += 1
        
        # Only write if content changed
        if content != original_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
            return changes_made
        return 0
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return 0

def main():
    """Main function to translate remaining files"""
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    else:
        target_dir = "."
    
    # Find files with Chinese content
    files_with_chinese = []
    for root, dirs, files in os.walk(target_dir):
        # Skip certain directories
        if any(skip in root for skip in ['third_party', 'node_modules', 'drawio-mcp-server', '.git']):
            continue
            
        for file in files:
            if file.endswith('.py'):
                file_path = os.path.join(root, file)
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                        if re.search(r'[\u4e00-\u9fff]', content):
                            files_with_chinese.append(file_path)
                except:
                    continue
    
    print(f"Found {len(files_with_chinese)} files with Chinese content...")
    
    total_changes = 0
    translated_files = 0
    
    for file_path in files_with_chinese:
        changes = translate_file(file_path)
        if changes > 0:
            print(f"Translated {changes} terms in: {file_path}")
            total_changes += changes
            translated_files += 1
    
    print(f"\nTranslation completed!")
    print(f"Files modified: {translated_files}")
    print(f"Total translations: {total_changes}")

if __name__ == "__main__":
    main()
