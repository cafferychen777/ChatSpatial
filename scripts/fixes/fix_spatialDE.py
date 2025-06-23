#!/usr/bin/env python
"""
Fix SpatialDE compatibility issues with newer scipy/pandas versions.
"""

import os
import sys
import shutil


def fix_spatialDE():
    """Fix the derivative function call in SpatialDE."""
    
    try:
        import SpatialDE
        spatialDE_path = SpatialDE.__file__
        spatialDE_dir = os.path.dirname(spatialDE_path)
        base_file = os.path.join(spatialDE_dir, 'base.py')
        
        print(f"Found SpatialDE at: {base_file}")
        
        # Backup original file
        backup_file = base_file + '.bak'
        if not os.path.exists(backup_file):
            shutil.copy(base_file, backup_file)
            print(f"Created backup: {backup_file}")
        
        # Read the file
        with open(base_file, 'r') as f:
            content = f.read()
        
        # Fix 1: Remove the n=2 parameter from derivative call
        fixes_applied = 0
        
        old_line1 = "s2_logdelta = 1. / (derivative(LL_obj, np.log(max_delta), n=2) ** 2)"
        new_line1 = "s2_logdelta = 1. / (derivative(LL_obj, np.log(max_delta)) ** 2)"
        
        if old_line1 in content:
            content = content.replace(old_line1, new_line1)
            fixes_applied += 1
            
        # Fix 2: Remove the n=1 parameter from derivative call in line 231
        old_line2 = "s2_FSV = derivative(FSV, np.log(max_delta), n=1) ** 2 * s2_logdelta"
        new_line2 = "s2_FSV = derivative(FSV, np.log(max_delta)) ** 2 * s2_logdelta"
        
        if old_line2 in content:
            content = content.replace(old_line2, new_line2)
            fixes_applied += 1
            
        if fixes_applied > 0:
            print(f"✓ Fixed {fixes_applied} derivative function calls")
        else:
            print("⚠️  derivative fixes not needed or already applied")
        
        # Write the fixed content back
        with open(base_file, 'w') as f:
            f.write(content)
        
        print("✅ SpatialDE fixes applied successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Error fixing SpatialDE: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_spatialDE():
    """Test SpatialDE after fixes."""
    print("\n🧪 Testing SpatialDE...")
    
    try:
        import numpy as np
        import pandas as pd
        import SpatialDE
        
        # Create simple test data
        n_spots = 30
        n_genes = 10
        
        # Coordinates
        coords = np.random.rand(n_spots, 2) * 10
        
        # Expression data
        counts = pd.DataFrame(
            np.random.poisson(5, size=(n_spots, n_genes)),
            columns=[f'gene_{i}' for i in range(n_genes)]
        )
        
        # Normalize
        total_counts = counts.sum(axis=1)
        norm_counts = counts.div(total_counts, axis=0) * np.median(total_counts)
        log_counts = np.log1p(norm_counts)
        
        print("✓ Test data created")
        
        # Run SpatialDE
        results = SpatialDE.run(coords, log_counts)
        
        print("✅ SpatialDE test successful!")
        print(f"   Results shape: {results.shape}")
        print(f"   Columns: {list(results.columns)}")
        
        return True
        
    except Exception as e:
        print(f"❌ SpatialDE test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main function."""
    print("🔧 SpatialDE Compatibility Fixer")
    print("=" * 50)
    
    # Apply fixes
    success = fix_spatialDE()
    
    if success:
        # Test the fixes
        test_success = test_spatialDE()
        
        if test_success:
            print("\n✅ All fixes applied and tested successfully!")
        else:
            print("\n⚠️  Fixes applied but test failed. May need additional adjustments.")
    else:
        print("\n❌ Failed to apply fixes.")


if __name__ == "__main__":
    main()