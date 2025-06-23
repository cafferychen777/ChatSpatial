#!/usr/bin/env python3
"""
Script to migrate all tools in ChatSpatial to use MCP-compliant error handling.

This script:
1. Backs up the original server.py
2. Updates imports to use the new error handling utilities
3. Replaces @mcp_error_handler with @mcp_tool_error_handler()
4. Shows a diff of changes before applying
"""

import os
import sys
import shutil
import re
from pathlib import Path
from datetime import datetime
import difflib

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Files to migrate
FILES_TO_MIGRATE = [
    "chatspatial/server.py",
    "chatspatial/tools/preprocessing.py",
    "chatspatial/tools/visualization.py",
    "chatspatial/tools/annotation.py",
    "chatspatial/tools/spatial_analysis.py",
    "chatspatial/tools/differential.py",
    "chatspatial/tools/trajectory.py",
    "chatspatial/tools/deconvolution.py",
    "chatspatial/tools/spatial_genes.py",
    "chatspatial/tools/integration.py",
    "chatspatial/tools/communication.py",
    "chatspatial/tools/enrichment.py",
    "chatspatial/tools/spatial_domains.py"
]


def backup_file(filepath: str) -> str:
    """Create a backup of the file"""
    backup_dir = Path("backups/error_handling_migration")
    backup_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = Path(filepath).name
    backup_path = backup_dir / f"{filename}.{timestamp}.bak"
    
    shutil.copy2(filepath, backup_path)
    return str(backup_path)


def migrate_server_py(content: str) -> str:
    """Migrate server.py to use new error handling"""
    
    # 1. Update imports
    old_import_pattern = r'from \.utils\.error_handling import mcp_error_handler, sync_mcp_error_handler'
    new_import = 'from .utils.tool_error_handling import mcp_tool_error_handler'
    
    if old_import_pattern in content:
        content = content.replace(
            'from .utils.error_handling import mcp_error_handler, sync_mcp_error_handler',
            'from .utils.tool_error_handling import mcp_tool_error_handler'
        )
    else:
        # Add import after other utils imports
        content = re.sub(
            r'(from \.utils\.image_utils.*\n)',
            r'\1from .utils.tool_error_handling import mcp_tool_error_handler\n',
            content
        )
    
    # 2. Replace decorator usage
    content = re.sub(
        r'@mcp_error_handler\s*\n',
        '@mcp_tool_error_handler()\n',
        content
    )
    
    # 3. Update sync_mcp_error_handler references if any
    content = re.sub(
        r'@sync_mcp_error_handler\s*\n',
        '@mcp_tool_error_handler()\n',
        content
    )
    
    # 4. Update return types for tools that need Dict[str, Any]
    # This is more complex and requires parsing, so we'll add comments for manual review
    
    return content


def migrate_tool_file(content: str, filepath: str) -> str:
    """Migrate a tool file to use new error handling"""
    
    # 1. Add import if not present
    if 'from ..utils.tool_error_handling import' not in content:
        # Find a good place to add the import
        import_match = re.search(r'(from \.\.models.*\n)', content)
        if import_match:
            content = content[:import_match.end()] + \
                      'from ..utils.tool_error_handling import mcp_tool_error_handler\n' + \
                      content[import_match.end():]
        else:
            # Add after other imports
            import_match = re.search(r'(import.*\n\n)', content)
            if import_match:
                content = content[:import_match.end()] + \
                          'from ..utils.tool_error_handling import mcp_tool_error_handler\n\n' + \
                          content[import_match.end():]
    
    # 2. Add decorator to main function
    # Find the main processing function (usually has data_id parameter)
    function_pattern = r'(async def \w+\s*\([^)]*data_id:[^)]+\)[^:]*:)'
    
    def add_decorator(match):
        func_def = match.group(0)
        # Check if already has decorator
        pos = match.start()
        before = content[:pos]
        lines_before = before.split('\n')
        
        # Check last few lines for existing decorator
        for i in range(min(3, len(lines_before))):
            if '@mcp_tool_error_handler' in lines_before[-(i+1)]:
                return func_def  # Already has decorator
        
        # Add decorator
        indent = re.match(r'^(\s*)', func_def).group(1)
        return f'{indent}@mcp_tool_error_handler()\n{func_def}'
    
    # Apply decorator to functions
    content = re.sub(function_pattern, add_decorator, content)
    
    return content


def show_diff(original: str, modified: str, filename: str):
    """Show diff between original and modified content"""
    print(f"\n{'='*60}")
    print(f"Changes for {filename}:")
    print('='*60)
    
    diff = list(difflib.unified_diff(
        original.splitlines(keepends=True),
        modified.splitlines(keepends=True),
        fromfile=f"{filename} (original)",
        tofile=f"{filename} (modified)",
        n=3
    ))
    
    if diff:
        for line in diff:
            if line.startswith('+') and not line.startswith('+++'):
                print(f"\033[92m{line}\033[0m", end='')  # Green
            elif line.startswith('-') and not line.startswith('---'):
                print(f"\033[91m{line}\033[0m", end='')  # Red
            else:
                print(line, end='')
    else:
        print("No changes needed")


def main():
    """Main migration function"""
    print("MCP Tool Error Handling Migration Script")
    print("="*60)
    
    # Check if we're in the right directory
    if not os.path.exists("chatspatial/server.py"):
        print("Error: Please run this script from the ChatSpatial root directory")
        sys.exit(1)
    
    changes = []
    
    for filepath in FILES_TO_MIGRATE:
        if not os.path.exists(filepath):
            print(f"Skipping {filepath} (not found)")
            continue
        
        print(f"\nProcessing {filepath}...")
        
        # Read original content
        with open(filepath, 'r') as f:
            original = f.read()
        
        # Migrate based on file type
        if filepath.endswith("server.py"):
            modified = migrate_server_py(original)
        else:
            modified = migrate_tool_file(original, filepath)
        
        # Show diff
        show_diff(original, modified, filepath)
        
        if original != modified:
            changes.append((filepath, original, modified))
    
    if not changes:
        print("\nNo changes needed!")
        return
    
    # Ask for confirmation
    print(f"\n{'='*60}")
    print(f"Ready to apply changes to {len(changes)} files")
    print("This will:")
    print("1. Create backups of all modified files")
    print("2. Apply the MCP-compliant error handling")
    print("3. Update imports and decorators")
    
    response = input("\nProceed with migration? (yes/no): ")
    
    if response.lower() != 'yes':
        print("Migration cancelled")
        return
    
    # Apply changes
    for filepath, original, modified in changes:
        # Create backup
        backup_path = backup_file(filepath)
        print(f"Backed up {filepath} to {backup_path}")
        
        # Write modified content
        with open(filepath, 'w') as f:
            f.write(modified)
        print(f"Updated {filepath}")
    
    print("\n" + "="*60)
    print("Migration completed successfully!")
    print("\nNext steps:")
    print("1. Run the test script to verify error handling:")
    print("   python scripts/tests/test_mcp_error_handling.py")
    print("2. Test each tool to ensure proper error handling")
    print("3. Update any custom error handling in individual tools")
    print("\nNote: Some tools may need manual adjustments for:")
    print("- Return type changes (to Dict[str, Any])")
    print("- Custom error messages using convenience functions")
    print("- Complex error scenarios requiring manual handling")


if __name__ == "__main__":
    main()