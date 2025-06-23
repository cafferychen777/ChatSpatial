#!/usr/bin/env python
"""
Fix PASTE compatibility with newer POT versions.
"""

import os
import sys
import shutil


def fix_paste_pot_compatibility():
    """Fix the line_search signature issue in PASTE."""
    
    try:
        import paste
        paste_path = paste.__file__
        paste_dir = os.path.dirname(paste_path)
        paste_file = os.path.join(paste_dir, 'PASTE.py')
        
        print(f"Found PASTE at: {paste_file}")
        
        # Backup original file
        backup_file = paste_file + '.bak'
        if not os.path.exists(backup_file):
            shutil.copy(paste_file, backup_file)
            print(f"Created backup: {backup_file}")
        
        # Read the file
        with open(paste_file, 'r') as f:
            content = f.read()
        
        # Fix 1: Update line_search function signatures
        # The new POT version expects 6 arguments but the old signature had 5
        old_line_search = """def line_search(cost, G, deltaG, Mi, cost_G, **kwargs):
            return ot.optim.line_search_armijo(cost, G, deltaG, Mi, cost_G, nx=nx, **kwargs)"""
        
        new_line_search = """def line_search(cost, G, deltaG, Mi, cost_G, alpha_min=0., alpha_max=1., **kwargs):
            return ot.optim.line_search_armijo(cost, G, deltaG, Mi, cost_G, alpha_min=alpha_min, alpha_max=alpha_max, **kwargs)"""
        
        if "def line_search(cost, G, deltaG, Mi, cost_G, **kwargs):" in content:
            content = content.replace(
                "def line_search(cost, G, deltaG, Mi, cost_G, **kwargs):",
                "def line_search(cost, G, deltaG, Mi, cost_G, alpha_min=0., alpha_max=1., **kwargs):"
            )
            print("✓ Fixed line_search function signatures")
        
        # Fix 2: Update the ot.optim.cg call to match new API
        # Check POT version to apply appropriate fix
        import ot
        pot_version = tuple(map(int, ot.__version__.split('.')[:2]))
        
        if pot_version >= (0, 9):
            # For POT >= 0.9, the API has changed
            # Replace the cg function call
            old_cg_call = "ot.optim.cg(p, q, (1 - alpha) * M, alpha, f, df, G0, line_search"
            
            # Check if we need to fix the cg call
            if old_cg_call in content:
                # Modern POT uses a different API
                print("✓ Detected POT >= 0.9, applying compatibility fix")
                
                # Add a wrapper function
                wrapper = '''
# Compatibility wrapper for POT >= 0.9
def cg_compat(p, q, M, alpha, f, df, G0, line_search, **kwargs):
    """Compatibility wrapper for different POT versions."""
    import ot
    try:
        # Try new API first (POT >= 0.9)
        from ot.optim import generic_conditional_gradient
        def loss_fun(G):
            return f(G)
        
        def gradient(G):
            return df(G)
        
        # Adjust line_search wrapper for new API
        def line_search_wrapper(cost, G, deltaG, Mi, cost_G, **ls_kwargs):
            # Handle both old and new signatures
            try:
                return line_search(cost, G, deltaG, Mi, cost_G, **ls_kwargs)
            except TypeError:
                # If it fails, try with alpha parameters
                return line_search(cost, G, deltaG, Mi, cost_G, 0., 1., **ls_kwargs)
        
        return generic_conditional_gradient(
            p, q, M, loss_fun, gradient, G0, 
            linesearch=line_search_wrapper,
            **kwargs
        )
    except (ImportError, AttributeError, TypeError):
        # Fall back to old API
        return ot.optim.cg(p, q, M, alpha, f, df, G0, line_search, **kwargs)

'''
                # Insert the wrapper before the my_fused_gromov_wasserstein function
                if "def my_fused_gromov_wasserstein" in content:
                    insert_pos = content.find("def my_fused_gromov_wasserstein")
                    content = content[:insert_pos] + wrapper + content[insert_pos:]
                    
                    # Replace the cg calls with our wrapper
                    content = content.replace(
                        "res, log = ot.optim.cg(p, q, (1 - alpha) * M, alpha, f, df, G0, line_search, log=True",
                        "res, log = cg_compat(p, q, (1 - alpha) * M, alpha, f, df, G0, line_search, log=True"
                    )
                    content = content.replace(
                        "return ot.optim.cg(p, q, (1 - alpha) * M, alpha, f, df, G0, line_search",
                        "return cg_compat(p, q, (1 - alpha) * M, alpha, f, df, G0, line_search"
                    )
                    print("✓ Added compatibility wrapper for ot.optim.cg")
        
        # Write the fixed content back
        with open(paste_file, 'w') as f:
            f.write(content)
        
        print("✅ PASTE compatibility fixes applied successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Error fixing PASTE: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_paste_after_fix():
    """Test PASTE after applying fixes."""
    print("\n🧪 Testing PASTE after fixes...")
    
    try:
        import numpy as np
        import anndata as ad
        import paste
        
        # Create simple test data
        n_spots = 30
        n_genes = 50
        
        adata1 = ad.AnnData(X=np.random.rand(n_spots, n_genes))
        adata1.obsm['spatial'] = np.random.rand(n_spots, 2) * 10
        adata1.var_names = [f"gene_{i}" for i in range(n_genes)]
        
        adata2 = ad.AnnData(X=np.random.rand(n_spots, n_genes))
        adata2.obsm['spatial'] = np.random.rand(n_spots, 2) * 10 + 2
        adata2.var_names = [f"gene_{i}" for i in range(n_genes)]
        
        print("✓ Test data created")
        
        # Try alignment with minimal parameters
        pi = paste.pairwise_align(
            adata1, adata2, 
            alpha=0.1,
            numItermax=100,  # Reduce iterations for testing
            verbose=False
        )
        
        print("✅ PASTE alignment successful!")
        print(f"   Alignment matrix shape: {pi.shape}")
        print(f"   Sum of transport plan: {pi.sum():.4f}")
        
        return True
        
    except Exception as e:
        print(f"❌ PASTE test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main function."""
    print("🔧 PASTE Compatibility Fixer")
    print("=" * 50)
    
    # Apply fixes
    success = fix_paste_pot_compatibility()
    
    if success:
        # Test the fixes
        test_success = test_paste_after_fix()
        
        if test_success:
            print("\n✅ All fixes applied and tested successfully!")
        else:
            print("\n⚠️  Fixes applied but test failed. May need additional adjustments.")
    else:
        print("\n❌ Failed to apply fixes.")


if __name__ == "__main__":
    main()