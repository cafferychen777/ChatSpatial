#!/usr/bin/env python3
"""
Script to organize and clean up test files
"""

import os
import shutil

# Root directory
ROOT_DIR = '/Users/apple/Research/SpatialTrans_MCP/chatspatial'

# Test files to organize
test_files = {
    # Trajectory tests - consolidated into test_trajectory_analysis.py
    'trajectory': [
        'test_trajectory_simple.py',
        'test_cellrank_complete.py',
        'test_cellrank_debug.py',
        'test_cellrank_detailed.py',
        'test_cellrank_direct.py',
        'test_cellrank_fixed.py',
        'test_cellrank_without_petsc.py',
        'test_mcp_trajectory.py',
        'test_velocity_trajectory.py'
    ],
    # Enrichment tests - consolidated into test_enrichment_analysis.py
    'enrichment': [
        'test_gsea_quick.py',
        'test_generic_enrichment.py',
        'test_ssgsea_debug.py'
    ]
}

# Archive directory for old tests
ARCHIVE_DIR = os.path.join(ROOT_DIR, 'tests', 'archive', 'root_tests')

def main():
    """Organize test files"""
    print("Organizing test files...")
    
    # Create archive directory
    os.makedirs(ARCHIVE_DIR, exist_ok=True)
    
    # Move test files to archive
    for category, files in test_files.items():
        category_dir = os.path.join(ARCHIVE_DIR, category)
        os.makedirs(category_dir, exist_ok=True)
        
        for test_file in files:
            src = os.path.join(ROOT_DIR, test_file)
            if os.path.exists(src):
                dst = os.path.join(category_dir, test_file)
                print(f"Moving {test_file} to archive/{category}/")
                shutil.move(src, dst)
    
    # List remaining test files in root
    remaining_tests = [f for f in os.listdir(ROOT_DIR) 
                      if f.startswith('test_') and f.endswith('.py')]
    
    if remaining_tests:
        print("\nRemaining test files in root:")
        for f in remaining_tests:
            print(f"  - {f}")
    else:
        print("\nAll test files have been organized!")
    
    print(f"\nOld test files archived to: {ARCHIVE_DIR}")
    print("New consolidated tests are in: tests/integration/")

if __name__ == "__main__":
    main()