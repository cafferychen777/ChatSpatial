#!/usr/bin/env python
"""
Fix SpatialDE compatibility issues with newer scipy/pandas versions.
"""

import os
import sys
import shutil


def fix_spatialDE():
    """Fix SpatialDE compatibility issues with newer scipy/pandas versions."""
    
    try:
        import SpatialDE
        spatialDE_path = SpatialDE.__file__
        spatialDE_dir = os.path.dirname(spatialDE_path)
        
        # Files to fix
        files_to_fix = [
            ('base.py', os.path.join(spatialDE_dir, 'base.py')),
            ('util.py', os.path.join(spatialDE_dir, 'util.py'))
        ]
        
        total_fixes = 0
        
        for filename, filepath in files_to_fix:
            if not os.path.exists(filepath):
                print(f"⚠️  File not found: {filepath}")
                continue
                
            print(f"\n🔧 Fixing {filename}...")
            
            # Backup original file
            backup_file = filepath + '.bak'
            if not os.path.exists(backup_file):
                shutil.copy(filepath, backup_file)
                print(f"   Created backup: {backup_file}")
            
            # Read the file
            with open(filepath, 'r') as f:
                content = f.read()
            
            original_content = content
            fixes_applied = 0
            
            # Fix 1: Replace scipy functions with numpy equivalents
            scipy_fixes = [
                ('sp.arange', 'np.arange'),
                ('sp.argsort', 'np.argsort'),
                ('sp.array', 'np.array'),
                ('sp.zeros', 'np.zeros'),
                ('sp.ones', 'np.ones'),
                ('sp.sort', 'np.sort')
            ]
            
            for old_func, new_func in scipy_fixes:
                if 'import scipy as sp' in content and old_func in content:
                    content = content.replace(old_func, new_func)
                    fixes_applied += 1
                    print(f"   ✓ Fixed {old_func} -> {new_func}")
                
            # Fix 2: Add numpy import if needed
            if 'np.arange' in content and 'import numpy as np' not in content:
                # Add numpy import after scipy import
                if 'import scipy as sp' in content:
                    content = content.replace('import scipy as sp', 'import scipy as sp\nimport numpy as np')
                    fixes_applied += 1
                    print(f"   ✓ Added numpy import")
            
            # Fix 3: Remove the n=2 parameter from derivative call
            old_line1 = "s2_logdelta = 1. / (derivative(LL_obj, np.log(max_delta), n=2) ** 2)"
            new_line1 = "s2_logdelta = 1. / (derivative(LL_obj, np.log(max_delta)) ** 2)"
            
            if old_line1 in content:
                content = content.replace(old_line1, new_line1)
                fixes_applied += 1
                print(f"   ✓ Fixed derivative function call (n=2)")
                
            # Fix 4: Remove the n=1 parameter from derivative call
            old_line2 = "s2_FSV = derivative(FSV, np.log(max_delta), n=1) ** 2 * s2_logdelta"
            new_line2 = "s2_FSV = derivative(FSV, np.log(max_delta)) ** 2 * s2_logdelta"
            
            if old_line2 in content:
                content = content.replace(old_line2, new_line2)
                fixes_applied += 1
                print(f"   ✓ Fixed derivative function call (n=1)")
                
            # Fix 5: pandas Series.ravel deprecation - more specific fixes
            pandas_fixes = [
                ('pv.ravel()', 'pv.to_numpy().ravel()'),
                ('.ravel()  # flattens', '.to_numpy().ravel()  # flattens'),
                ('pvals.ravel()', 'pvals.to_numpy().ravel()'),
            ]
            
            for old_pattern, new_pattern in pandas_fixes:
                if old_pattern in content:
                    content = content.replace(old_pattern, new_pattern)
                    fixes_applied += 1
                    print(f"   ✓ Fixed pandas ravel: {old_pattern[:20]}...")
            
            # Write the fixed content back if changes were made
            if content != original_content:
                with open(filepath, 'w') as f:
                    f.write(content)
                print(f"   ✅ Applied {fixes_applied} fixes to {filename}")
            else:
                print(f"   ⚠️  No fixes needed or already applied in {filename}")
                
            total_fixes += fixes_applied
        
        print(f"\n✅ Total fixes applied: {total_fixes}")
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