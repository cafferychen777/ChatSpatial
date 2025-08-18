#!/usr/bin/env python3
"""
Comprehensive sc-type testing across 10+ tissue types and different scales
"""

import asyncio
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
import traceback
import time
import psutil
import gc
from typing import Dict, List, Tuple

# Tissue-specific marker genes (based on sc-type database and literature)
TISSUE_MARKERS = {
    "Brain": {
        "Neurons": ["RBFOX3", "MAP2", "TUBB3", "SYP", "SNAP25", "GRIN1", "ENO2"],
        "Astrocytes": ["GFAP", "AQP4", "S100B", "ALDH1L1", "GJA1", "SLC1A3"],
        "Oligodendrocytes": ["MBP", "MOG", "PLP1", "OLIG2", "MAG", "MOBP"],
        "Microglia": ["AIF1", "CX3CR1", "P2RY12", "TMEM119", "TREM2"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "FLT1"]
    },
    "Heart": {
        "Cardiomyocytes": ["MYH6", "MYH7", "TNNT2", "TNNI3", "ACTC1", "MYL2"],
        "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "THY1"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "KDR"],
        "Smooth_muscle": ["ACTA2", "MYH11", "TAGLN", "CNN1"],
        "Immune": ["PTPRC", "CD68", "CD3D", "LYZ"]
    },
    "Liver": {
        "Hepatocytes": ["ALB", "CYP3A4", "CYP2E1", "APOE", "FGB", "SERPINA1"],
        "Cholangiocytes": ["KRT19", "SOX9", "CFTR", "EPCAM"],
        "Stellate_cells": ["ACTA2", "COL1A1", "PDGFRA", "RBP1"],
        "Kupffer_cells": ["CD68", "CD163", "MARCO", "LYZ"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "FCN3"]
    },
    "Lung": {
        "Pneumocytes_I": ["AGER", "PDPN", "CAV1", "RTKN2"],
        "Pneumocytes_II": ["SFTPC", "SFTPA1", "SFTPB", "LAMP3"],
        "Alveolar_macrophages": ["CD68", "MARCO", "CHIT1", "FABP4"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "CA4"],
        "Fibroblasts": ["COL1A1", "COL3A1", "PDGFRA", "LUM"],
        "Ciliated_epithelial": ["FOXJ1", "DNAI1", "CCDC39", "PIFO"]
    },
    "Kidney": {
        "Podocytes": ["NPHS1", "NPHS2", "WT1", "SYNPO"],
        "Proximal_tubule": ["LRP2", "SLC34A1", "CUBN", "SLC5A12"],
        "Distal_tubule": ["NCCT", "TRPV5", "CALBINDIN1", "SLC8A1"],
        "Collecting_duct": ["AQP2", "ATP6V1G3", "ATP6V0D2", "SCNN1A"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "PLVAP"],
        "Immune": ["PTPRC", "CD68", "CD3D", "LYZ"]
    },
    "Pancreas": {
        "Beta_cells": ["INS", "IAPP", "MAFA", "NKX6-1", "PCSK1", "G6PC2"],
        "Alpha_cells": ["GCG", "ARX", "TTR", "GC", "LOXL4"],
        "Delta_cells": ["SST", "HHEX", "LEPR", "RGS16"],
        "PP_cells": ["PPY", "SSTR2", "ETV1"],
        "Acinar_cells": ["CPA1", "PRSS1", "AMY2A", "CELA3B"],
        "Ductal_cells": ["KRT19", "SOX9", "CFTR", "CA2"]
    },
    "Muscle": {
        "Skeletal_myocytes": ["MYH1", "MYH2", "MYH4", "ACTN3", "CKM", "MYLPF"],
        "Satellite_cells": ["PAX7", "MYF5", "MYOD1", "VCAM1"],
        "Fibroblasts": ["COL1A1", "COL3A1", "PDGFRA", "TCF21"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "KDR"],
        "Smooth_muscle": ["ACTA2", "MYH11", "TAGLN", "MYLK"]
    },
    "Skin": {
        "Keratinocytes": ["KRT1", "KRT10", "KRT14", "KRT5", "TP63", "IVL"],
        "Fibroblasts": ["COL1A1", "COL3A1", "ELN", "DCN"],
        "Melanocytes": ["MLANA", "TYR", "TYRP1", "DCT", "PMEL"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "PLVAP"],
        "Immune": ["PTPRC", "CD68", "CD3D", "CD207"],
        "Sebaceous": ["PPARG", "FASN", "ELOVL3", "PLIN2"]
    },
    "Intestine": {
        "Enterocytes": ["ALPI", "SI", "FABP2", "APOA1", "RBP2"],
        "Goblet_cells": ["MUC2", "AGR2", "SPDEF", "FCGBP"],
        "Paneth_cells": ["LYZ", "DEFA5", "DEFA6", "PLA2G2A"],
        "Enteroendocrine": ["CHGA", "CHGB", "GCG", "GIP", "CCK"],
        "Stem_cells": ["LGR5", "ASCL2", "OLFM4", "SMOC2"],
        "Immune": ["PTPRC", "CD68", "CD3D", "FOXP3"]
    },
    "Prostate": {
        "Luminal_epithelial": ["KRT8", "KRT18", "AR", "NKX3-1", "KLK3"],
        "Basal_epithelial": ["KRT5", "KRT14", "TP63", "MYC", "GSTP1"],
        "Stromal_fibroblasts": ["COL1A1", "ACTA2", "PDGFRA", "FAP"],
        "Smooth_muscle": ["ACTA2", "MYH11", "TAGLN", "CNN1"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "KDR"],
        "Immune": ["PTPRC", "CD68", "CD3D", "CD19"]
    },
    "Mammary_gland": {
        "Luminal_mature": ["KRT8", "KRT18", "LALBA", "CSN3", "GLYCAM1"],
        "Luminal_progenitor": ["KRT8", "KRT19", "ELF5", "MSX1"],
        "Basal_myoepithelial": ["KRT14", "KRT5", "ACTA2", "TP63", "OXTR"],
        "Stromal_fibroblasts": ["COL1A1", "DCN", "LUM", "PDGFRA"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "PLVAP"],
        "Immune": ["PTPRC", "CD68", "CD3D", "CD19"]
    },
    "Immune_system": {
        "T_cells": ["CD3D", "CD3E", "CD3G", "CD8A", "CD4", "IL7R"],
        "B_cells": ["CD19", "MS4A1", "CD79A", "CD79B", "IGHM"],
        "NK_cells": ["KLRD1", "NCR1", "GNLY", "PRF1", "GZMB"],
        "Macrophages": ["CD68", "CD14", "CSF1R", "MARCO", "MSR1"],
        "Dendritic_cells": ["CD1C", "CLEC9A", "XCR1", "IRF8", "ZBTB46"],
        "Monocytes": ["CD14", "FCGR3A", "LYZ", "S100A8", "S100A9"]
    }
}

def create_tissue_specific_data(tissue_name: str, n_cells: int, n_genes: int, seed: int = 42) -> ad.AnnData:
    """Create tissue-specific synthetic data with realistic marker expression"""
    np.random.seed(seed)
    
    if tissue_name not in TISSUE_MARKERS:
        raise ValueError(f"Unknown tissue: {tissue_name}. Available: {list(TISSUE_MARKERS.keys())}")
    
    markers = TISSUE_MARKERS[tissue_name]
    cell_types = list(markers.keys())
    
    # Create base expression matrix
    X = np.random.negative_binomial(n=2, p=0.3, size=(n_cells, n_genes)).astype(float)
    
    # Collect all marker genes
    all_markers = []
    for cell_type, marker_list in markers.items():
        all_markers.extend(marker_list)
    
    # Remove duplicates and create gene names
    unique_markers = list(set(all_markers))
    other_genes = [f'{tissue_name.upper()}_GENE_{i:04d}' for i in range(len(unique_markers), n_genes)]
    gene_names = unique_markers + other_genes
    
    # Create gene name to index mapping
    gene_to_idx = {gene: i for i, gene in enumerate(gene_names)}
    
    # Assign cells to types roughly equally
    cells_per_type = n_cells // len(cell_types)
    cell_type_assignments = []
    
    for i, cell_type in enumerate(cell_types):
        start_idx = i * cells_per_type
        end_idx = start_idx + cells_per_type if i < len(cell_types) - 1 else n_cells
        cell_type_assignments.extend([cell_type] * (end_idx - start_idx))
    
    # Add tissue-specific marker expression patterns
    for i, assigned_type in enumerate(cell_type_assignments):
        if assigned_type in markers:
            for marker_gene in markers[assigned_type]:
                if marker_gene in gene_to_idx:
                    gene_idx = gene_to_idx[marker_gene]
                    # Add higher expression for this marker in cells of this type
                    additional_expr = np.random.negative_binomial(n=8, p=0.2)
                    X[i, gene_idx] += additional_expr
    
    # Create spatial coordinates (simulate spatial arrangement)
    x_coords = np.random.uniform(0, 50, n_cells)
    y_coords = np.random.uniform(0, 50, n_cells)
    
    # Create AnnData object
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({
            'cell_id': [f"{tissue_name}_{i:05d}" for i in range(n_cells)],
            'true_celltype': cell_type_assignments,
            'tissue': [tissue_name] * n_cells
        }),
        var=pd.DataFrame({
            'gene_id': gene_names,
            'highly_variable': [gene in unique_markers for gene in gene_names]
        })
    )
    
    adata.obs_names = [f"{tissue_name}_{i:05d}" for i in range(n_cells)]
    adata.var_names = gene_names
    adata.obsm['spatial'] = np.column_stack([x_coords, y_coords])
    
    # Standard preprocessing
    adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['n_genes_by_counts'] = np.array((adata.X > 0).sum(axis=1)).flatten()
    
    # Normalization and scaling
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.scale(adata)
    adata.layers['scaled'] = adata.X.copy()
    adata.X = adata.raw.X  # Restore counts
    
    return adata

def extract_mcp_result(mcp_response):
    """Extract result from MCP response"""
    if isinstance(mcp_response, dict):
        if 'isError' in mcp_response and mcp_response['isError']:
            raise Exception(f"MCP Error: {mcp_response.get('content', 'Unknown error')}")
        
        if 'content' in mcp_response and mcp_response['content']:
            content = mcp_response['content'][0]
            if content['type'] == 'text':
                import json
                return json.loads(content['text'])
    return mcp_response

async def test_tissue_dataset(tissue_name: str, n_cells: int, n_genes: int) -> Dict:
    """Test sc-type on a specific tissue dataset"""
    from chatspatial.server import load_data, annotate_cells
    
    # Create tissue data
    adata = create_tissue_specific_data(tissue_name, n_cells, n_genes)
    test_path = Path(f"sctype_test_{tissue_name.lower()}_{n_cells}cells.h5ad")
    adata.write(test_path)
    
    results = {
        'tissue': tissue_name,
        'n_cells': n_cells,
        'n_genes': n_genes,
        'success': False,
        'cell_types_found': 0,
        'cells_annotated': 0,
        'processing_time': 0,
        'memory_usage_mb': 0,
        'error': None,
        'true_cell_types': len(set(adata.obs['true_celltype'])),
        'marker_genes_used': len([g for g in adata.var_names if g in sum(TISSUE_MARKERS[tissue_name].values(), [])])
    }
    
    try:
        start_time = time.time()
        start_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        # Load data
        load_result = await load_data(str(test_path), "h5ad", f"test_{tissue_name}")
        dataset = extract_mcp_result(load_result)
        
        # Test sc-type with appropriate tissue type
        tissue_mapping = {
            "Brain": "Brain",
            "Heart": "Heart", 
            "Liver": "Liver",
            "Lung": "Lung",
            "Kidney": "Kidney",
            "Pancreas": "Pancreas",
            "Muscle": "Muscle",
            "Skin": "Skin",
            "Intestine": "Intestine", 
            "Prostate": "Prostate",
            "Mammary_gland": "Mammary gland",
            "Immune_system": "Immune system"
        }
        
        sctype_tissue = tissue_mapping.get(tissue_name, "Immune system")  # Default fallback
        
        annotation_result = await annotate_cells(
            dataset['id'],
            {
                "method": "sctype",
                "sctype_tissue": sctype_tissue,
                "sctype_scaled": True,
                "sctype_use_cache": True
            }
        )
        
        annotation_data = extract_mcp_result(annotation_result)
        
        end_time = time.time()
        end_memory = psutil.Process().memory_info().rss / 1024 / 1024
        
        results.update({
            'success': True,
            'cell_types_found': len(annotation_data.get('counts', {})),
            'cells_annotated': sum(annotation_data.get('counts', {}).values()),
            'processing_time': end_time - start_time,
            'memory_usage_mb': end_memory - start_memory,
            'annotation_counts': annotation_data.get('counts', {}),
            'confidence_scores': annotation_data.get('confidence_scores', {})
        })
        
    except Exception as e:
        results['error'] = str(e)
        
    finally:
        # Cleanup
        if test_path.exists():
            test_path.unlink()
        # Force garbage collection
        gc.collect()
    
    return results

async def run_comprehensive_tissue_tests():
    """Run comprehensive tests across all tissue types and scales"""
    print("🧪 COMPREHENSIVE TISSUE-SPECIFIC SC-TYPE TESTING")
    print("=" * 70)
    
    # Check availability
    try:
        from chatspatial.tools.annotation import SCTYPE_AVAILABLE
        if not SCTYPE_AVAILABLE:
            print("❌ SC-type not available, skipping tests")
            return
    except Exception as e:
        print(f"❌ Failed to check SC-type availability: {e}")
        return
    
    # Test configurations: (tissue, cells, genes)
    test_configs = [
        # Small scale tests (all tissues)
        ("Brain", 100, 50),
        ("Heart", 100, 50), 
        ("Liver", 100, 50),
        ("Lung", 100, 50),
        ("Kidney", 100, 50),
        ("Pancreas", 100, 50),
        ("Muscle", 100, 50),
        ("Skin", 100, 50),
        ("Intestine", 100, 50),
        ("Prostate", 100, 50),
        ("Mammary_gland", 100, 50),
        ("Immune_system", 100, 50),
        
        # Medium scale tests (selected tissues)
        ("Brain", 500, 200),
        ("Heart", 500, 200),
        ("Liver", 500, 200),
        ("Lung", 500, 200),
        ("Immune_system", 500, 200),
        
        # Large scale tests (key tissues)
        ("Brain", 2000, 1000),
        ("Immune_system", 2000, 1000),
        ("Liver", 2000, 1000),
        
        # Very large scale tests (stress test)
        ("Brain", 5000, 2000),
        ("Immune_system", 5000, 2000),
    ]
    
    print(f"🎯 Running {len(test_configs)} test configurations...")
    print(f"📊 Testing {len(set(t[0] for t in test_configs))} different tissue types")
    print(f"📈 Scale range: {min(t[1] for t in test_configs)} - {max(t[1] for t in test_configs)} cells")
    print()
    
    results = []
    successful_tests = 0
    total_processing_time = 0
    
    for i, (tissue, n_cells, n_genes) in enumerate(test_configs):
        print(f"[{i+1:2d}/{len(test_configs)}] Testing {tissue} ({n_cells} cells, {n_genes} genes)...")
        
        try:
            result = await test_tissue_dataset(tissue, n_cells, n_genes)
            results.append(result)
            
            if result['success']:
                successful_tests += 1
                total_processing_time += result['processing_time']
                status = "✅"
                detail = f"{result['cell_types_found']} types, {result['processing_time']:.2f}s"
            else:
                status = "❌"
                detail = f"Error: {result['error'][:50]}..."
                
            print(f"    {status} {detail}")
            
        except Exception as e:
            print(f"    ❌ Test failed: {e}")
            results.append({
                'tissue': tissue,
                'n_cells': n_cells, 
                'n_genes': n_genes,
                'success': False,
                'error': str(e)
            })
    
    # Generate comprehensive report
    print("\n" + "=" * 70)
    print("COMPREHENSIVE TEST REPORT")
    print("=" * 70)
    
    print(f"📊 Overall Success Rate: {successful_tests}/{len(test_configs)} ({successful_tests/len(test_configs)*100:.1f}%)")
    print(f"⏱️  Total Processing Time: {total_processing_time:.2f}s")
    print(f"⚡ Average Time per Test: {total_processing_time/max(successful_tests,1):.2f}s")
    
    # Tissue-wise summary
    tissue_stats = {}
    for result in results:
        tissue = result['tissue']
        if tissue not in tissue_stats:
            tissue_stats[tissue] = {'success': 0, 'total': 0, 'avg_time': 0, 'times': []}
        
        tissue_stats[tissue]['total'] += 1
        if result['success']:
            tissue_stats[tissue]['success'] += 1
            if 'processing_time' in result:
                tissue_stats[tissue]['times'].append(result['processing_time'])
    
    print(f"\n📋 TISSUE-WISE RESULTS ({len(tissue_stats)} tissues tested):")
    for tissue, stats in sorted(tissue_stats.items()):
        success_rate = stats['success'] / stats['total'] * 100
        avg_time = sum(stats['times']) / len(stats['times']) if stats['times'] else 0
        print(f"  {tissue:15s}: {stats['success']:2d}/{stats['total']:2d} ({success_rate:5.1f}%) - Avg: {avg_time:5.2f}s")
    
    # Scale analysis  
    scale_stats = {}
    for result in results:
        if result['success']:
            scale = "Small" if result['n_cells'] <= 200 else "Medium" if result['n_cells'] <= 1000 else "Large" if result['n_cells'] <= 3000 else "Very Large"
            if scale not in scale_stats:
                scale_stats[scale] = {'count': 0, 'times': [], 'memory': []}
            scale_stats[scale]['count'] += 1
            if 'processing_time' in result:
                scale_stats[scale]['times'].append(result['processing_time'])
            if 'memory_usage_mb' in result:
                scale_stats[scale]['memory'].append(result['memory_usage_mb'])
    
    print(f"\n📈 SCALE ANALYSIS:")
    for scale, stats in scale_stats.items():
        avg_time = sum(stats['times']) / len(stats['times']) if stats['times'] else 0
        avg_memory = sum(stats['memory']) / len(stats['memory']) if stats['memory'] else 0
        print(f"  {scale:10s}: {stats['count']:2d} tests - Avg time: {avg_time:5.2f}s - Avg memory: {avg_memory:+6.1f}MB")
    
    # Performance insights
    successful_results = [r for r in results if r['success']]
    if successful_results:
        fastest = min(successful_results, key=lambda x: x.get('processing_time', float('inf')))
        slowest = max(successful_results, key=lambda x: x.get('processing_time', 0))
        
        print(f"\n⚡ PERFORMANCE INSIGHTS:")
        print(f"  Fastest: {fastest['tissue']} ({fastest['n_cells']} cells) - {fastest['processing_time']:.2f}s")
        print(f"  Slowest: {slowest['tissue']} ({slowest['n_cells']} cells) - {slowest['processing_time']:.2f}s")
    
    # Error analysis
    failed_results = [r for r in results if not r['success']]
    if failed_results:
        print(f"\n❌ FAILED TESTS ({len(failed_results)}):")
        for result in failed_results:
            print(f"  {result['tissue']} ({result['n_cells']} cells): {result.get('error', 'Unknown error')[:60]}...")
    
    # Final assessment
    print(f"\n🎯 FINAL ASSESSMENT:")
    if successful_tests / len(test_configs) >= 0.9:
        print("🎉 EXCELLENT: SC-type integration is highly robust across tissues and scales!")
    elif successful_tests / len(test_configs) >= 0.8:
        print("✅ GOOD: SC-type integration works well with minor issues in some cases.")
    elif successful_tests / len(test_configs) >= 0.7:
        print("⚠️  ACCEPTABLE: SC-type integration works but needs improvement in some areas.")
    else:
        print("❌ NEEDS WORK: SC-type integration has significant issues that need addressing.")
    
    print("=" * 70)
    
    return results

if __name__ == "__main__":
    asyncio.run(run_comprehensive_tissue_tests())