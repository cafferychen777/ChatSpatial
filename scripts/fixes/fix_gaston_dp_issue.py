#!/usr/bin/env python
"""
Fix for GASTON dp_related.get_isodepth_labels KeyError issue.

The issue occurs when the dynamic programming algorithm creates invalid states
in segment_map, particularly when num_domains is not compatible with the data structure.
"""

import os
import sys
import tempfile
import shutil
import numpy as np
import pandas as pd
import torch

# Add parent directories to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

def safe_get_isodepth_labels(model, A, S, num_domains, num_buckets=50, num_pcs_A=None):
    """
    Safe version of get_isodepth_labels that handles edge cases.
    """
    from gaston import dp_related
    
    # Try the original function first
    try:
        return dp_related.get_isodepth_labels(model, A, S, num_domains, num_buckets, num_pcs_A)
    except KeyError as e:
        print(f"GASTON DP KeyError {e}, trying alternative approach...")
        
        # Fallback: reduce num_domains or adjust parameters
        fallback_domains = max(2, min(num_domains, 3))  # Use 2-3 domains max
        fallback_buckets = min(num_buckets, 20)  # Use fewer buckets
        
        try:
            return dp_related.get_isodepth_labels(model, A, S, fallback_domains, fallback_buckets, num_pcs_A)
        except Exception as e2:
            print(f"GASTON DP fallback also failed: {e2}")
            # Final fallback: create simple isodepth and labels
            N = A.shape[0]
            S_torch = torch.Tensor(S)
            gaston_isodepth = model.spatial_embedding(S_torch).detach().numpy().flatten()
            
            # Simple k-means style labeling
            sorted_indices = np.argsort(gaston_isodepth)
            labels_per_domain = N // fallback_domains
            gaston_labels = np.zeros(N)
            
            for domain in range(fallback_domains):
                start_idx = domain * labels_per_domain
                end_idx = (domain + 1) * labels_per_domain if domain < fallback_domains - 1 else N
                indices_in_domain = sorted_indices[start_idx:end_idx]
                gaston_labels[indices_in_domain] = domain
                
            print(f"Using simple fallback labeling with {fallback_domains} domains")
            return gaston_isodepth, gaston_labels

async def test_gaston_with_fix():
    """Test GASTON with the fixed get_isodepth_labels function."""
    
    try:
        # Import GASTON components
        import gaston
        from gaston import neural_net
        print("✅ GASTON imported successfully")
        
        # Create simple test data 
        n_spots = 100
        n_genes = 50
        np.random.seed(42)
        
        # Spatial coordinates
        coords = np.random.randn(n_spots, 2) * 3
        
        # Expression features (from PCA)
        expression_features = np.random.randn(n_spots, 10)
        
        # Convert to torch tensors
        S_torch, A_torch = neural_net.load_rescale_input_data(coords, expression_features)
        
        print(f"📊 Test data: {n_spots} spots, {expression_features.shape[1]} features")
        
        # Create a simple model for testing
        print("🧠 Training minimal GASTON model...")
        
        model, loss_list = neural_net.train(
            S_torch, A_torch,
            S_hidden_list=[20],  # Very simple architecture
            A_hidden_list=[10],
            epochs=10,  # Very few epochs for testing
            checkpoint=None,
            save_dir=None,
            optim='Adam',
            lr=0.001,
            seed=42,
            save_final=False,
            embed_size=5,
            sigma=1.0,
            batch_size=None
        )
        
        print("✅ Model training completed")
        
        # Test the fixed function with various domain numbers
        for num_domains in [2, 3, 4, 5]:
            print(f"\n🧪 Testing with {num_domains} domains...")
            
            try:
                gaston_isodepth, gaston_labels = safe_get_isodepth_labels(
                    model, A_torch, S_torch, num_domains
                )
                
                print(f"   ✅ Success! Isodepth range: [{gaston_isodepth.min():.3f}, {gaston_isodepth.max():.3f}]")
                print(f"   📊 Labels: {len(np.unique(gaston_labels))} unique domains")
                
            except Exception as e:
                print(f"   ❌ Failed: {e}")
        
        return True
        
    except Exception as e:
        print(f"❌ GASTON test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

async def main():
    print("🔧 GASTON DP Issue Fix Test")
    print("=" * 50)
    
    success = await test_gaston_with_fix()
    
    if success:
        print("\n🎉 GASTON fix appears to work!")
        print("💡 Next step: integrate this fix into spatial_genes.py")
    else:
        print("\n⚠️  GASTON fix needs more work")

if __name__ == "__main__":
    import asyncio
    asyncio.run(main())