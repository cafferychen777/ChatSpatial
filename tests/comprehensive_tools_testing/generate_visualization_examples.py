#!/usr/bin/env python3
"""
Generate example visualizations to demonstrate ChatSpatial's capabilities.

This script creates a gallery of visualization examples using different 
plot types and spatial data formats to showcase the visualization module's
comprehensive functionality.
"""

import asyncio
import sys
from pathlib import Path
import matplotlib.pyplot as plt

# Add chatspatial to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import visualize_data
from chatspatial.utils.data_loader import load_spatial_data
import scanpy as sc

plt.switch_backend('Agg')  # Use non-interactive backend

class VisualizationGalleryGenerator:
    """Generate a gallery of visualization examples"""
    
    def __init__(self):
        self.output_dir = Path(__file__).parent / "visualization_examples"
        self.output_dir.mkdir(exist_ok=True)
        
    async def generate_gallery(self):
        """Generate complete visualization gallery"""
        print("🎨 Generating ChatSpatial Visualization Gallery")
        print("=" * 60)
        
        # Load test datasets
        datasets = await self._load_datasets()
        
        # Generate different types of visualizations
        examples = [
            # Basic spatial visualizations
            ("spatial_expression", "spatial", {"feature": "Gene_0001"}, 
             "Spatial gene expression visualization"),
            ("spatial_clusters", "spatial", {"feature": "leiden"}, 
             "Spatial clustering visualization"),
            
            # UMAP visualizations  
            ("umap_clusters", "umap", {"feature": "leiden"}, 
             "UMAP clustering visualization"),
            ("umap_expression", "umap", {"feature": "Gene_0001"}, 
             "UMAP gene expression visualization"),
            
            # Heatmaps
            ("expression_heatmap", "heatmap", {}, 
             "Gene expression heatmap"),
            
            # Multi-gene visualizations
            ("multi_gene_spatial", "multi_gene", 
             {"feature": ["Gene_0001", "Gene_0002", "Gene_0003", "Gene_0004"]}, 
             "Multi-gene spatial visualization"),
            
            # Violin plots
            ("gene_violin", "violin", {"feature": "Gene_0001"}, 
             "Gene expression violin plot"),
            
            # Deconvolution results
            ("deconvolution_results", "deconvolution", {}, 
             "Cell type deconvolution results"),
        ]
        
        # Generate examples for each dataset
        for dataset_name, dataset in datasets.items():
            print(f"\n📊 Generating examples for {dataset_name}")
            print("-" * 40)
            
            for example_name, plot_type, params, description in examples:
                await self._generate_example(
                    f"{dataset_name}_{example_name}",
                    dataset_name, 
                    dataset,
                    plot_type,
                    params,
                    description
                )
        
        # Generate summary
        await self._generate_summary()
        
        print(f"\n✅ Gallery generated successfully!")
        print(f"📁 Output directory: {self.output_dir}")
        
    async def _load_datasets(self):
        """Load available test datasets"""
        datasets = {}
        
        # Try to load existing datasets
        data_paths = [
            ("enhanced_visium", "data/test_datasets/enhanced_visium.h5ad"),
            ("synthetic_visium", "data/test/synthetic_visium.h5ad"),
            ("high_res_merfish", "data/test_datasets/high_res_merfish.h5ad"),
        ]
        
        base_path = Path(__file__).parent.parent.parent
        
        for name, rel_path in data_paths:
            full_path = base_path / rel_path
            if full_path.exists():
                try:
                    print(f"📂 Loading {name}...")
                    adata = sc.read_h5ad(full_path)
                    datasets[name] = {"adata": adata}
                    print(f"   ✅ Loaded: {adata.shape}")
                except Exception as e:
                    print(f"   ❌ Failed to load {name}: {e}")
        
        if not datasets:
            print("⚠️  No datasets found, creating minimal synthetic data")
            datasets = await self._create_minimal_dataset()
            
        return datasets
    
    async def _create_minimal_dataset(self):
        """Create minimal synthetic dataset for examples"""
        import numpy as np
        import pandas as pd
        
        # Create simple synthetic data
        n_cells, n_genes = 200, 100
        coords = np.random.uniform(0, 50, size=(n_cells, 2))
        X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
        
        adata = sc.AnnData(X=X)
        adata.obsm['spatial'] = coords
        adata.var_names = [f'Gene_{i:04d}' for i in range(n_genes)]
        adata.var_names_make_unique()
        
        # Add basic preprocessing
        adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
        sc.tl.umap(adata)
        
        # Add mock deconvolution
        cell_types = ['CellType_A', 'CellType_B', 'CellType_C']
        deconv_matrix = np.random.dirichlet(np.ones(len(cell_types)), size=n_cells)
        adata.obsm['deconvolution_mock'] = deconv_matrix
        adata.uns['deconvolution_mock_cell_types'] = cell_types
        
        return {"synthetic_example": {"adata": adata}}
    
    async def _generate_example(self, filename, dataset_name, dataset, 
                              plot_type, params, description):
        """Generate a single visualization example"""
        try:
            # Create parameters
            viz_params = VisualizationParameters(plot_type=plot_type, **params)
            
            # Generate visualization
            result = await visualize_data(
                dataset_name, 
                {dataset_name: dataset}, 
                viz_params
            )
            
            # Save the image
            if hasattr(result, '_data') or hasattr(result, 'data'):
                data = getattr(result, '_data', getattr(result, 'data'))
                output_path = self.output_dir / f"{filename}.png"
                
                with open(output_path, 'wb') as f:
                    f.write(data)
                
                print(f"   ✅ {filename}: {description}")
                
            else:
                print(f"   ⚠️  {filename}: Unexpected result type")
                
        except Exception as e:
            print(f"   ❌ {filename}: Failed - {e}")
    
    async def _generate_summary(self):
        """Generate a summary HTML file"""
        html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>ChatSpatial Visualization Gallery</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .gallery { display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }
        .example { border: 1px solid #ddd; padding: 15px; border-radius: 8px; }
        .example img { max-width: 100%; height: auto; }
        .example h3 { margin-top: 0; color: #2c3e50; }
        .dataset-section { margin-bottom: 40px; }
        .dataset-title { color: #e74c3c; border-bottom: 2px solid #e74c3c; padding-bottom: 5px; }
    </style>
</head>
<body>
    <h1>ChatSpatial Visualization Gallery</h1>
    <p>This gallery demonstrates the comprehensive visualization capabilities of ChatSpatial across different spatial transcriptomics technologies and plot types.</p>
"""
        
        # List all generated images
        image_files = list(self.output_dir.glob("*.png"))
        
        # Group by dataset
        datasets = {}
        for img_file in image_files:
            parts = img_file.stem.split('_', 1)
            if len(parts) >= 2:
                dataset = parts[0]
                example = parts[1]
                if dataset not in datasets:
                    datasets[dataset] = []
                datasets[dataset].append((img_file.name, example))
        
        for dataset, examples in datasets.items():
            html_content += f'''
    <div class="dataset-section">
        <h2 class="dataset-title">{dataset.replace('_', ' ').title()}</h2>
        <div class="gallery">
'''
            for img_name, example_name in examples:
                title = example_name.replace('_', ' ').title()
                html_content += f'''
            <div class="example">
                <h3>{title}</h3>
                <img src="{img_name}" alt="{title}">
            </div>
'''
            html_content += '''
        </div>
    </div>
'''
        
        html_content += '''
    <footer style="margin-top: 60px; padding-top: 20px; border-top: 1px solid #ddd; color: #666;">
        <p>Generated by ChatSpatial Comprehensive Visualization Test Suite</p>
        <p>All visualizations created using the chatspatial.tools.visualization module</p>
    </footer>
</body>
</html>
'''
        
        # Save HTML file
        html_path = self.output_dir / "gallery.html"
        with open(html_path, 'w') as f:
            f.write(html_content)
        
        print(f"📄 Gallery summary saved: {html_path}")

async def main():
    """Main function"""
    generator = VisualizationGalleryGenerator()
    await generator.generate_gallery()

if __name__ == "__main__":
    asyncio.run(main())