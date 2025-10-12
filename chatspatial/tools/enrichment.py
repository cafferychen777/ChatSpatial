"""
Enrichment analysis tools for spatial transcriptomics data.

This module provides both standard and spatially-aware enrichment analysis methods:
- Standard methods: GSEA, ORA, ssGSEA, Enrichr (via gseapy)
- Spatial methods: EnrichMap-based spatial enrichment analysis
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from mcp.server.fastmcp import Context
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..utils.error_handling import ProcessingError
from ..utils.metadata_storage import store_analysis_metadata

logger = logging.getLogger(__name__)


# ============================================================================
# GENE FORMAT CONVERSION UTILITIES
# ============================================================================

def _convert_gene_format_for_matching(
    pathway_genes: List[str], 
    dataset_genes: set, 
    species: str
) -> Tuple[List[str], Dict[str, str]]:
    """
    Rule-based gene format conversion to match dataset format.
    
    Handles common gene format variations between pathway databases and datasets:
    - Uppercase (GENE) vs Title case (Gene) vs lowercase (gene)
    - Species-specific formatting rules
    - Special prefixes like Gm/GM/gm for mouse genes
    
    Args:
        pathway_genes: Gene names from pathway database (usually uppercase from gseapy)
        dataset_genes: Available gene names in dataset  
        species: Species specified by user ("mouse" or "human")
    
    Returns:
        (dataset_format_genes, conversion_map)
        dataset_format_genes: Gene names in dataset format that can be found
        conversion_map: Maps dataset_format -> original_pathway_format
    """
    dataset_format_genes = []
    conversion_map = {}
    
    for gene in pathway_genes:
        # Try direct match first
        if gene in dataset_genes:
            dataset_format_genes.append(gene)
            conversion_map[gene] = gene
            continue
            
        # Apply multiple format conversion rules
        format_variations = []
        
        if species == "mouse":
            # Mouse-specific format rules (order matters for efficiency)
            # Rule 1: Title case (most common): Cd5l, Gbp2b
            if len(gene) > 1:
                format_variations.append(gene[0].upper() + gene[1:].lower())
            # Rule 2: All lowercase: cd5l, gbp2b  
            format_variations.append(gene.lower())
            # Rule 3: All uppercase: CD5L, GBP2B
            format_variations.append(gene.upper())
            # Rule 4: Capitalize first letter only
            format_variations.append(gene.capitalize())
            
            # Special rule for Gm-prefixed genes (common in mouse)
            if gene.upper().startswith('GM'):
                format_variations.extend([
                    'gm' + gene[2:].lower(),  # gm42418
                    'Gm' + gene[2:].lower(),  # Gm42418  
                    'GM' + gene[2:].upper(),  # GM42418
                ])
                
        elif species == "human":
            # Human-specific format rules
            # Rule 1: All uppercase (most common): HES1, FABP4
            format_variations.append(gene.upper())
            # Rule 2: All lowercase: hes1, fabp4
            format_variations.append(gene.lower())
            # Rule 3: Capitalize first letter
            format_variations.append(gene.capitalize())
        
        # Remove duplicates while preserving order
        seen = set()
        unique_variations = []
        for variation in format_variations:
            if variation not in seen and variation != gene:  # Skip if same as original
                seen.add(variation)
                unique_variations.append(variation)
        
        # Try each format variation against dataset
        for variant in unique_variations:
            if variant in dataset_genes:
                dataset_format_genes.append(variant)  # Use dataset's actual format
                conversion_map[variant] = gene
                break  # Stop after first match
    
    return dataset_format_genes, conversion_map


# ============================================================================
# MCP PROTOCOL COMPLIANCE: STANDARDIZED DATA STRUCTURES
# ============================================================================

@dataclass
class EnrichmentInternalResult:
    """Internal standardized result for enrichment analysis.
    
    Optimized for MCP protocol compliance with token efficiency.
    Focuses on statistical value and essential summaries rather than complete gene lists.
    
    Attributes:
        method: Analysis method name
        n_gene_sets: Number of gene sets analyzed
        n_significant: Number of significant gene sets
        enrichment_scores: Enrichment scores dictionary
        pvalues: Raw p-values dictionary
        adjusted_pvalues: Adjusted p-values dictionary
        gene_set_statistics: Gene set statistical information
        gene_set_summaries: Compact summaries with counts and sample genes
        top_gene_sets: Top enriched gene sets
        top_depleted_sets: Top depleted gene sets
        spatial_metrics: Spatial analysis metrics (optional)
        spatial_scores_key: Spatial scores key (optional)
        method_specific_data: Method-specific extension data
    """
    # Basic information
    method: str
    n_gene_sets: int
    n_significant: int
    
    # Core statistical results
    enrichment_scores: Dict[str, float]
    pvalues: Dict[str, Optional[float]]
    adjusted_pvalues: Dict[str, Optional[float]]
    gene_set_statistics: Dict[str, Dict[str, Any]]
    
    # Gene set summaries (token-optimized)
    gene_set_summaries: Dict[str, Dict[str, Any]]  # Contains count + sample genes only
    
    # Ranking results
    top_gene_sets: List[str]
    top_depleted_sets: List[str]
    
    # Optional spatial analysis results
    spatial_metrics: Optional[Dict[str, Any]] = None
    spatial_scores_key: Optional[str] = None
    
    # Method-specific extension fields
    method_specific_data: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format optimized for MCP token limits.
        
        Returns:
            Token-optimized dictionary focusing on statistical value and essential summaries.
            Prioritizes enrichment statistics over complete gene lists.
        """
        # Optimize gene_set_statistics by removing complete gene lists
        optimized_gene_set_statistics = {}
        if self.gene_set_statistics:
            for pathway, stats in self.gene_set_statistics.items():
                if isinstance(stats, dict):
                    # Create optimized version excluding complete gene lists
                    optimized_stats = {k: v for k, v in stats.items() if k != "genes"}
                    # Add gene count info if genes were present
                    if "genes" in stats and isinstance(stats["genes"], list):
                        optimized_stats["gene_count"] = len(stats["genes"])
                        optimized_stats["sample_genes"] = stats["genes"][:3]  # Only first 3 genes
                    optimized_gene_set_statistics[pathway] = optimized_stats
                else:
                    # Non-dict statistics, pass through as-is
                    optimized_gene_set_statistics[pathway] = stats

        base_dict = {
            "method": self.method,
            "n_gene_sets": self.n_gene_sets,
            "n_significant": self.n_significant,
            "enrichment_scores": self.enrichment_scores,
            "pvalues": self.pvalues,
            "adjusted_pvalues": self.adjusted_pvalues,
            "gene_set_statistics": optimized_gene_set_statistics,  # Optimized version
            "gene_set_summaries": self.gene_set_summaries,  # Optimized field
            "top_gene_sets": self.top_gene_sets,
            "top_depleted_sets": self.top_depleted_sets,
        }
        
        # Add spatial analysis results
        if self.spatial_metrics is not None:
            base_dict["spatial_metrics"] = self.spatial_metrics
        if self.spatial_scores_key is not None:
            base_dict["spatial_scores_key"] = self.spatial_scores_key
            
        # Add method-specific data
        base_dict.update(self.method_specific_data)
        
        return base_dict


@dataclass
class SpatialMetricResult:
    """Spatial metrics calculation result.
    
    Standardized return format for _compute_spatial_metric function.
    """
    data_id: str
    score_key: str
    metrics: Dict[str, Dict[str, float]]
    parameters: Dict[str, Any]


@dataclass
class ClusterCorrelationResult:
    """Cluster correlation analysis result.
    
    Standardized return format for _correlate_with_clusters function.
    """
    data_id: str
    signature_name: str
    cluster_key: str
    correlation_method: str
    cluster_correlations: Dict[str, Any]


@dataclass
class SSGSEAResult:
    """Result from single-sample Gene Set Enrichment Analysis."""
    
    method: str
    n_gene_sets: int
    n_samples: int
    enrichment_scores: Dict[str, List[float]]  # gene_set -> list of sample scores
    gene_set_summaries: Dict[str, Dict[str, Any]]  # Compact summaries with counts and sample genes
    summary_stats: Dict[str, Dict[str, float]]  # gene_set -> {mean, std, etc.}
    score_matrix: Optional[np.ndarray] = None  # samples x gene_sets matrix
    method_specific_data: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary optimized for MCP token limits."""
        return {
            "method": self.method,
            "n_gene_sets": self.n_gene_sets,
            "n_samples": self.n_samples,
            "enrichment_scores": self.enrichment_scores,
            "gene_set_summaries": self.gene_set_summaries,
            "summary_stats": self.summary_stats,
            "score_matrix": self.score_matrix.tolist() if self.score_matrix is not None else None,
            **self.method_specific_data
        }


@dataclass  
class EnrichrResult:
    """Result from Enrichr web service analysis."""
    
    method: str
    n_gene_sets: int
    n_significant: int
    input_gene_summary: Dict[str, Any]  # Summary of input genes: count + sample
    enrichment_scores: Dict[str, float]  # gene_set -> combined score
    pvalues: Dict[str, float]
    adjusted_pvalues: Dict[str, float] 
    odds_ratios: Dict[str, float]
    gene_set_statistics: Dict[str, Dict[str, Any]]
    top_gene_sets: List[str]
    libraries_used: List[str]
    method_specific_data: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary optimized for MCP token limits."""
        return {
            "method": self.method,
            "n_gene_sets": self.n_gene_sets,
            "n_significant": self.n_significant,
            "input_gene_summary": self.input_gene_summary,
            "enrichment_scores": self.enrichment_scores,
            "pvalues": self.pvalues,
            "adjusted_pvalues": self.adjusted_pvalues,
            "odds_ratios": self.odds_ratios,
            "gene_set_statistics": self.gene_set_statistics,
            "gene_set_summaries": {},  # Enrichr doesn't have local gene sets
            "top_gene_sets": self.top_gene_sets,
            "libraries_used": self.libraries_used,
            **self.method_specific_data
        }


# ============================================================================
# INSIGHTS-BASED DATA STRUCTURES (New Value-Driven Approach)
# ============================================================================

@dataclass
class PathwayResult:
    """Individual pathway finding with essential information only."""
    name: str                    # Clean pathway name
    category: str               # "Immune", "Metabolism", etc.
    significance: str           # "Very High (p<0.001)"
    effect_size: str            # "2.3x enriched"
    key_genes: List[str]        # 2-3 representative genes
    description: str            # Brief biological context

@dataclass
class CategorySummary:
    """Grouped pathway insights."""
    count: int
    significance_level: str
    top_themes: List[str]

@dataclass
class EnrichmentInsights:
    """Value-driven enrichment analysis result optimized for insights and conversation."""
    
    # === PRIMARY INSIGHTS ===
    summary: str                                    # Key findings in 1-2 sentences
    significance_overview: str                      # Overall significance assessment
    
    # === TOP FINDINGS ===
    top_pathways: List[PathwayResult]              # 10-20 most significant only
    pathway_categories: Dict[str, CategorySummary]  # Grouped insights
    
    # === ACTIONABLE INFORMATION ===
    key_genes: List[str]                           # 5-10 most important genes
    suggested_analyses: List[str]                  # Follow-up suggestions
    interpretation: str                            # What this means
    
    # === TECHNICAL METADATA ===
    method: str
    total_pathways_tested: int
    significant_count: int
    stats_summary: Dict[str, Any]                  # Minimal reproducibility info
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format for server.py compatibility while maintaining insights focus."""
        
        # Create compatible data structure that server.py expects but with insights content
        enrichment_scores = {}
        pvalues = {}
        adjusted_pvalues = {}
        gene_set_statistics = {}
        
        # Transform insights back to server-expected format
        for pathway in self.top_pathways:
            # Use cleaned pathway name as key
            key = pathway.name.replace(" ", "_").lower()
            
            # Create mock enrichment scores (server expects these)
            enrichment_scores[key] = 2.0  # Default enriched value
            pvalues[key] = None  # Spatial analysis doesn't have p-values
            adjusted_pvalues[key] = None
            
            # Create minimal gene_set_statistics (token-optimized)
            gene_set_statistics[key] = {
                "category": pathway.category,
                "effect_size": pathway.effect_size,
                "significance": pathway.significance,
                "description": pathway.description,
                "key_genes": pathway.key_genes  # Only 2 genes each
                # No "genes" field with complete gene lists!
            }
        
        # Return server.py compatible format with insights content
        return {
            # === CORE INSIGHTS (Primary Information) ===
            "summary": self.summary,
            "significance_overview": self.significance_overview,
            "interpretation": self.interpretation,
            "key_genes": self.key_genes,
            "suggested_analyses": self.suggested_analyses,
            "pathway_categories": {
                cat: {
                    "count": summary.count,
                    "significance_level": summary.significance_level,
                    "top_themes": summary.top_themes
                }
                for cat, summary in self.pathway_categories.items()
            },
            "top_pathways": [
                {
                    "name": p.name,
                    "category": p.category,
                    "significance": p.significance,
                    "effect_size": p.effect_size,
                    "key_genes": p.key_genes,
                    "description": p.description
                }
                for p in self.top_pathways
            ],
            
            # === SERVER COMPATIBILITY FIELDS ===
            "method": self.method,
            "n_gene_sets": self.total_pathways_tested,
            "n_significant": self.significant_count,
            "enrichment_scores": enrichment_scores,
            "pvalues": pvalues,
            "adjusted_pvalues": adjusted_pvalues,
            "gene_set_statistics": gene_set_statistics,  # Optimized, no full gene lists
            "gene_set_summaries": {},  # Insights are in pathway_categories
            "top_gene_sets": [p.name.replace(" ", "_").lower() for p in self.top_pathways[:10]],
            "top_depleted_sets": [],  # Spatial enrichment doesn't have depleted sets
            "spatial_metrics": None,
            "spatial_scores_key": None,
            
            # === TECHNICAL METADATA ===
            "stats_summary": self.stats_summary
        }


# ============================================================================
# INSIGHTS GENERATION HELPER FUNCTIONS
# ============================================================================

def clean_pathway_name(go_term: str) -> str:
    """Convert GO:0006955_immune_response to 'Immune Response'."""
    if "_" in go_term:
        # Split on underscore and clean up
        parts = go_term.split("_", 1)
        if len(parts) > 1:
            name_part = parts[1]
            # Replace underscores with spaces and title case
            return name_part.replace("_", " ").title()
    return go_term.replace("_", " ").title()

def categorize_pathway(pathway_name: str) -> str:
    """Map pathway to high-level category."""
    pathway_lower = pathway_name.lower()
    
    # Define category mappings
    categories = {
        "Immune": ["immune", "defense", "inflammation", "interferon", "cytokine", "antigen"],
        "Metabolism": ["metabolic", "glycol", "lipid", "glucose", "energy", "atp", "respiratory"],
        "Cell Cycle": ["cycle", "division", "mitosis", "meiosis", "proliferation", "growth"],
        "Signaling": ["signal", "transduction", "response", "pathway", "cascade"],
        "Transport": ["transport", "localization", "trafficking", "secretion", "endocytosis"],
        "Development": ["development", "differentiation", "morphogenesis", "embryo"],
        "Cell Death": ["death", "apoptosis", "necrosis", "autophagy"],
        "DNA/RNA": ["transcription", "translation", "dna", "rna", "splicing", "repair"],
        "Protein": ["protein", "folding", "modification", "ubiquitin", "proteolysis"],
        "Structure": ["cytoskeleton", "organization", "assembly", "maintenance"]
    }
    
    # Find best category match
    for category, keywords in categories.items():
        if any(keyword in pathway_lower for keyword in keywords):
            return category
    
    return "Other"


def format_effect_size(score: float, baseline: float = 1.0) -> str:
    """Convert score to interpretable effect size."""
    if score > baseline:
        fold_change = score / baseline
        if fold_change >= 3.0:
            return f"{fold_change:.1f}x enriched (very strong)"
        elif fold_change >= 2.0:
            return f"{fold_change:.1f}x enriched (strong)"
        elif fold_change >= 1.5:
            return f"{fold_change:.1f}x enriched (moderate)"
        else:
            return f"{fold_change:.1f}x enriched (mild)"
    else:
        return f"{score:.2f} (below baseline)"

def generate_insights_summary(pathways: List[PathwayResult]) -> str:
    """Generate 1-2 sentence summary of key findings."""
    if not pathways:
        return "No significant pathways found."
    
    # Count categories
    categories = {}
    for pathway in pathways:
        categories[pathway.category] = categories.get(pathway.category, 0) + 1
    
    # Get top categories
    top_categories = sorted(categories.items(), key=lambda x: x[1], reverse=True)[:3]
    
    if len(top_categories) == 1:
        return f"Strong enrichment in {top_categories[0][0].lower()} pathways with {len(pathways)} significant findings."
    elif len(top_categories) == 2:
        return f"Significant enrichment in {top_categories[0][0].lower()} and {top_categories[1][0].lower()} pathways."
    else:
        top_names = [cat[0].lower() for cat in top_categories]
        return f"Enrichment across multiple processes including {', '.join(top_names[:-1])}, and {top_names[-1]}."

def suggest_followup_analyses(pathways: List[PathwayResult]) -> List[str]:
    """Suggest relevant follow-up analyses based on findings."""
    suggestions = []
    
    # Get categories present
    categories = set(pathway.category for pathway in pathways)
    
    # Category-specific suggestions
    if "Immune" in categories:
        suggestions.append("Perform cell communication analysis to identify immune signaling patterns")
        suggestions.append("Validate key inflammatory markers with targeted gene expression analysis")
    
    if "Metabolism" in categories:
        suggestions.append("Analyze metabolic pathway activity across spatial regions")
        suggestions.append("Check for metabolic heterogeneity between cell types")
    
    if "Cell Cycle" in categories:
        suggestions.append("Perform trajectory analysis to understand proliferation dynamics")
        suggestions.append("Identify proliferative zones with spatial domain analysis")
    
    if "Signaling" in categories:
        suggestions.append("Map signaling pathway activity across tissue regions")
        suggestions.append("Investigate pathway crosstalk with network analysis")
    
    # General suggestions
    suggestions.append("Visualize spatial distribution of top enriched pathways")
    suggestions.append("Validate findings with complementary experimental methods")
    
    # Return top 4-6 suggestions
    return suggestions[:6]

def generate_pathway_interpretation(pathways: List[PathwayResult]) -> str:
    """Generate biological interpretation from pathway patterns."""
    if not pathways:
        return "No significant pathways found for interpretation."
    
    # Analyze categories and their significance
    categories = {}
    for pathway in pathways:
        if pathway.category not in categories:
            categories[pathway.category] = []
        categories[pathway.category].append(pathway)
    
    # Generate interpretation based on patterns
    interpretations = []
    
    if "Immune" in categories and len(categories["Immune"]) >= 2:
        interpretations.append("Strong immune response activation suggests inflammatory processes or tissue response to stimuli")
    
    if "Metabolism" in categories and "Cell Cycle" in categories:
        interpretations.append("Coordinated metabolic and proliferative activity indicates active tissue remodeling or repair")
    
    if "Signaling" in categories and len(categories["Signaling"]) >= 2:
        interpretations.append("Multiple signaling pathway activation suggests complex cellular communication and coordination")
    
    if len(categories) >= 4:
        interpretations.append("Broad pathway activation across multiple biological processes indicates significant physiological changes")
    
    # Default interpretation
    if not interpretations:
        top_category = max(categories.keys(), key=lambda x: len(categories[x]))
        interpretations.append(f"Primary enrichment in {top_category.lower()} processes suggests focused biological activity in this domain")
    
    return ". ".join(interpretations) + "."


def map_gene_set_database_to_enrichr_library(database_name: str, species: str) -> str:
    """Map user-friendly database names to actual Enrichr library names.
    
    Args:
        database_name: User-friendly database name from MCP interface
        species: Species ('human', 'mouse', or 'zebrafish')
        
    Returns:
        Actual Enrichr library name
        
    Raises:
        ValueError: If database_name is not supported
    """
    mapping = {
        "GO_Biological_Process": "GO_Biological_Process_2025",
        "GO_Molecular_Function": "GO_Molecular_Function_2025", 
        "GO_Cellular_Component": "GO_Cellular_Component_2025",
        "KEGG_Pathways": "KEGG_2021_Human" if species.lower() == "human" else "KEGG_2019_Mouse",
        "Reactome_Pathways": "Reactome_Pathways_2024",
        "MSigDB_Hallmark": "MSigDB_Hallmark_2020",
        "Cell_Type_Markers": "CellMarker_Augmented_2021"
    }
    
    if database_name not in mapping:
        available_options = list(mapping.keys())
        raise ValueError(
            f"Unknown gene set database: {database_name}. "
            f"Available options: {available_options}"
        )
    
    return mapping[database_name]


# ============================================================================
# Helper Functions
# ============================================================================


def create_gene_set_summaries(
    gene_sets: Dict[str, List[str]], 
    sample_size: int = 3
) -> Dict[str, Dict[str, Any]]:
    """Create compact gene set summaries for token optimization.
    
    Args:
        gene_sets: Dictionary of gene sets with full gene lists
        sample_size: Number of sample genes to include per set
        
    Returns:
        Dictionary with compact summaries containing count and sample genes only
    """
    summaries = {}
    for pathway_name, genes in gene_sets.items():
        summaries[pathway_name] = {
            "gene_count": len(genes),
            "sample_genes": genes[:sample_size],  # First N genes as samples
            "total_available": len(genes)
        }
    return summaries


def create_input_gene_summary(
    input_genes: List[str],
    genes_found: List[str] = None,
    sample_size: int = 3
) -> Dict[str, Any]:
    """Create compact input gene summary for Enrichr results.
    
    Args:
        input_genes: List of input genes
        genes_found: List of genes found in database (optional)
        sample_size: Number of sample genes to include
        
    Returns:
        Dictionary with compact summary of input genes
    """
    return {
        "total_input": len(input_genes),
        "genes_found_count": len(genes_found) if genes_found else len(input_genes),
        "sample_input_genes": input_genes[:sample_size],
        "sample_found_genes": genes_found[:sample_size] if genes_found else input_genes[:sample_size]
    }


# ============================================================================
# Standard Enrichment Analysis Functions (Non-spatial)
# ============================================================================


def is_gseapy_available() -> Tuple[bool, str]:
    """Check if gseapy is available."""
    try:
        import gseapy

        return True, ""
    except ImportError:
        return False, "gseapy not installed. Install with: pip install gseapy"


async def perform_gsea(
    adata,
    gene_sets: Dict[str, List[str]],
    ranking_key: Optional[str] = None,
    method: str = "signal_to_noise",
    permutation_num: int = 1000,
    min_size: int = 10,
    max_size: int = 500,
    species: Optional[str] = None,
    database: Optional[str] = None,
    context=None,
) -> EnrichmentInternalResult:
    """
    Perform Gene Set Enrichment Analysis (GSEA).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : Dict[str, List[str]]
        Gene sets to test
    ranking_key : Optional[str]
        Key in adata.var for pre-computed ranking. If None, compute from expression
    method : str
        Method for ranking genes if ranking_key is None
    permutation_num : int
        Number of permutations
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size
    species : Optional[str]
        Species for the analysis (e.g., 'mouse', 'human')
    database : Optional[str]
        Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process')
    context : Optional
        MCP context

    Returns
    -------
    Dict containing enrichment results
    """
    is_available, error_msg = is_gseapy_available()
    if not is_available:
        raise ImportError(error_msg)

    import gseapy as gp

    if context:
        await context.info("Running GSEA analysis...")

    # Prepare ranking
    if ranking_key and ranking_key in adata.var:
        # Use pre-computed ranking
        ranking = adata.var[ranking_key].to_dict()
    else:
        # Compute ranking from expression data
        if "log1p" in adata.uns:
            X = adata.X
        else:
            X = adata.raw.X if adata.raw else adata.X

        # Simple fold change if we have groups
        if "condition" in adata.obs or "group" in adata.obs:
            group_key = "condition" if "condition" in adata.obs else "group"
            groups = adata.obs[group_key].unique()
            if len(groups) == 2:
                # Binary comparison
                group1_mask = adata.obs[group_key] == groups[0]
                group2_mask = adata.obs[group_key] == groups[1]

                mean1 = np.array(X[group1_mask, :].mean(axis=0)).flatten()
                mean2 = np.array(X[group2_mask, :].mean(axis=0)).flatten()

                # Compute fold change
                fc = np.log2((mean2 + 1) / (mean1 + 1))
                ranking = dict(zip(adata.var_names, fc))
            else:
                # Use variance as ranking for multi-group
                # Handle sparse matrices
                if hasattr(X, "todense"):
                    var_scores = np.array(X.todense().var(axis=0)).flatten()
                else:
                    var_scores = np.array(X.var(axis=0)).flatten()
                ranking = dict(zip(adata.var_names, var_scores))
        else:
            # Use variance as default ranking
            # Handle sparse matrices
            if hasattr(X, "todense"):
                var_scores = np.array(X.todense().var(axis=0)).flatten()
            else:
                var_scores = np.array(X.var(axis=0)).flatten()
            ranking = dict(zip(adata.var_names, var_scores))

    # Run GSEA preranked
    try:
        # Convert ranking dict to DataFrame for gseapy
        ranking_df = pd.DataFrame.from_dict(ranking, orient="index", columns=["score"])
        ranking_df.index.name = "gene"
        ranking_df = ranking_df.sort_values("score", ascending=False)

        res = gp.prerank(
            rnk=ranking_df,  # Pass DataFrame instead of dict
            gene_sets=gene_sets,
            processes=1,
            permutation_num=permutation_num,
            min_size=min_size,
            max_size=max_size,
            seed=42,
            verbose=False,
            no_plot=True,
            outdir=None,
        )

        # Extract results
        results_df = res.res2d

        # Prepare output
        enrichment_scores = {}
        pvalues = {}
        adjusted_pvalues = {}
        gene_set_statistics = {}

        for idx, row in results_df.iterrows():
            term = row["Term"]
            enrichment_scores[term] = row["ES"]
            pvalues[term] = row["NOM p-val"]
            adjusted_pvalues[term] = row["FDR q-val"]
            gene_set_statistics[term] = {
                "es": row["ES"],
                "nes": row["NES"],
                "pval": row["NOM p-val"],
                "fdr": row["FDR q-val"],
                "size": row.get(
                    "Matched_size", row.get("Gene %", 0)
                ),  # Different versions use different column names
                "lead_genes": (
                    row.get("Lead_genes", "").split(";")[:10]
                    if "Lead_genes" in row
                    else []
                ),
            }

        # Get top enriched and depleted
        results_df_sorted = results_df.sort_values("NES", ascending=False)
        top_enriched = (
            results_df_sorted[results_df_sorted["NES"] > 0].head(10)["Term"].tolist()
        )
        top_depleted = (
            results_df_sorted[results_df_sorted["NES"] < 0].head(10)["Term"].tolist()
        )

        # Save results to adata.uns for visualization
        # Store full results DataFrame for visualization
        adata.uns["gsea_results"] = results_df

        # Store gene set membership for validation
        adata.uns["enrichment_gene_sets"] = gene_sets

        # Store metadata for scientific provenance tracking
        store_analysis_metadata(
            adata,
            analysis_name="enrichment_gsea",
            method="gsea",
            parameters={
                "permutation_num": permutation_num,
                "ranking_method": method,
                "min_size": min_size,
                "max_size": max_size,
                "ranking_key": ranking_key,
            },
            results_keys={"uns": ["gsea_results", "enrichment_gene_sets"]},
            statistics={
                "n_gene_sets": len(gene_sets),
                "n_significant": len(results_df[results_df["FDR q-val"] < 0.05]),
            },
            species=species,
            database=database,
        )

        # Inform user about visualization options
        if context:
            await context.info(
                "GSEA analysis complete. Use create_visualization tool with plot_type='pathway_enrichment' to visualize results"
            )

        return EnrichmentInternalResult(
            method="gsea",
            n_gene_sets=len(gene_sets),
            n_significant=len(results_df[results_df["FDR q-val"] < 0.05]),
            enrichment_scores=enrichment_scores,
            pvalues=pvalues,
            adjusted_pvalues=adjusted_pvalues,
            gene_set_statistics=gene_set_statistics,
            gene_set_summaries=create_gene_set_summaries(gene_sets),
            top_gene_sets=top_enriched,
            top_depleted_sets=top_depleted,
            method_specific_data={
                "gene_sets_size_summary": {k: len(v) for k, v in gene_sets.items()},
                "permutation_num": permutation_num,
                "ranking_method": method,
                "min_size": min_size,
                "max_size": max_size,
            }
        )

    except Exception as e:
        logger.error(f"GSEA failed: {e}")
        raise


async def perform_ora(
    adata,
    gene_sets: Dict[str, List[str]],
    gene_list: Optional[List[str]] = None,
    pvalue_threshold: float = 0.05,
    logfc_threshold: float = 1.0,
    min_size: int = 10,
    max_size: int = 500,
    species: Optional[str] = None,
    database: Optional[str] = None,
    context=None,
) -> EnrichmentInternalResult:
    """
    Perform Over-Representation Analysis (ORA).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : Dict[str, List[str]]
        Gene sets to test
    gene_list : Optional[List[str]]
        List of genes to test. If None, use DEGs
    pvalue_threshold : float
        P-value threshold for selecting DEGs
    logfc_threshold : float
        Log fold change threshold for selecting DEGs
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size
    species : Optional[str]
        Species for the analysis (e.g., 'mouse', 'human')
    database : Optional[str]
        Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process')
    context : Optional
        MCP context

    Returns
    -------
    Dict containing enrichment results
    """
    if context:
        await context.info("Running Over-Representation Analysis...")

    # Get gene list if not provided
    if gene_list is None:
        # Try to get DEGs from adata
        if "rank_genes_groups" in adata.uns:
            # Get DEGs
            result = adata.uns["rank_genes_groups"]
            names = result["names"]

            # Check if pvals exist (not all rank_genes_groups have pvals)
            pvals = None
            if "pvals_adj" in result:
                pvals = result["pvals_adj"]
            elif "pvals" in result:
                pvals = result["pvals"]

            # Check if logfcs exist
            logfcs = result.get("logfoldchanges", None)

            # Get first group's DEGs
            degs = []
            for i in range(len(names[0])):
                # If pvals and logfcs exist, filter by thresholds
                if pvals is not None and logfcs is not None:
                    if (
                        pvals[0][i] < pvalue_threshold
                        and abs(logfcs[0][i]) > logfc_threshold
                    ):
                        degs.append(names[0][i])
                else:
                    # If no pvals/logfcs, just use top N genes (e.g., top 100)
                    if i < 100:
                        degs.append(names[0][i])

            gene_list = degs

            if context:
                if pvals is not None:
                    await context.info(f"Using {len(gene_list)} DEGs for ORA (filtered by p-value and log2FC)")
                else:
                    await context.info(f"Using top {len(gene_list)} ranked genes for ORA (no p-values available)")
        else:
            # Use highly variable genes
            if "highly_variable" in adata.var:
                gene_list = adata.var_names[adata.var["highly_variable"]].tolist()
            else:
                # Use top variable genes
                # Handle sparse matrices
                if hasattr(adata.X, "todense"):
                    var_scores = np.array(adata.X.todense().var(axis=0)).flatten()
                else:
                    var_scores = np.array(adata.X.var(axis=0)).flatten()
                top_indices = np.argsort(var_scores)[-500:]
                gene_list = adata.var_names[top_indices].tolist()

    # Background genes
    background_genes = set(adata.var_names)
    query_genes = set(gene_list) & background_genes

    # Perform hypergeometric test for each gene set
    enrichment_scores = {}
    pvalues = {}
    gene_set_statistics = {}

    for gs_name, gs_genes in gene_sets.items():
        gs_genes_set = set(gs_genes) & background_genes

        if len(gs_genes_set) < min_size or len(gs_genes_set) > max_size:
            continue

        # Hypergeometric test
        # a: genes in both query and gene set
        # b: genes in query but not in gene set
        # c: genes in gene set but not in query
        # d: genes in neither

        a = len(query_genes & gs_genes_set)
        b = len(query_genes - gs_genes_set)
        c = len(gs_genes_set - query_genes)
        d = len(background_genes - query_genes - gs_genes_set)

        # Fisher's exact test
        odds_ratio, p_value = stats.fisher_exact(
            [[a, b], [c, d]], alternative="greater"
        )

        enrichment_scores[gs_name] = odds_ratio
        pvalues[gs_name] = p_value

        gene_set_statistics[gs_name] = {
            "odds_ratio": odds_ratio,
            "pval": p_value,
            "overlap": a,
            "query_size": len(query_genes),
            "gs_size": len(gs_genes_set),
            "overlapping_genes": list(query_genes & gs_genes_set)[:20],  # Top 20
        }

    # Multiple testing correction
    if pvalues:
        pval_array = np.array(list(pvalues.values()))
        _, adjusted_pvals, _, _ = multipletests(pval_array, method="fdr_bh")
        adjusted_pvalues = dict(zip(pvalues.keys(), adjusted_pvals))
    else:
        adjusted_pvalues = {}

    # Get top results
    sorted_by_pval = sorted(pvalues.items(), key=lambda x: x[1])
    top_gene_sets = [x[0] for x in sorted_by_pval[:10]]

    # Save results to adata.uns for visualization
    # Create DataFrame for visualization compatibility
    import pandas as pd

    ora_df = pd.DataFrame(
        {
            "pathway": list(enrichment_scores.keys()),
            "odds_ratio": list(enrichment_scores.values()),
            "pvalue": [pvalues.get(k, 1.0) for k in enrichment_scores.keys()],
            "adjusted_pvalue": [
                adjusted_pvalues.get(k, 1.0) for k in enrichment_scores.keys()
            ],
        }
    )
    ora_df["NES"] = ora_df["odds_ratio"]  # Use odds_ratio as score for visualization
    ora_df = ora_df.sort_values("pvalue")

    adata.uns["ora_results"] = ora_df
    adata.uns["gsea_results"] = (
        ora_df  # Also save as gsea_results for visualization compatibility
    )

    # Store gene set membership for validation
    adata.uns["enrichment_gene_sets"] = gene_sets

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name="enrichment_ora",
        method="ora",
        parameters={
            "pvalue_threshold": pvalue_threshold,
            "logfc_threshold": logfc_threshold,
            "min_size": min_size,
            "max_size": max_size,
            "n_query_genes": len(query_genes),
        },
        results_keys={"uns": ["ora_results", "enrichment_gene_sets"]},
        statistics={
            "n_gene_sets": len(gene_sets),
            "n_significant": sum(
                1 for p in adjusted_pvalues.values() if p is not None and p < 0.05
            ),
            "n_query_genes": len(query_genes),
        },
        species=species,
        database=database,
    )

    # Inform user about visualization options
    if context:
        await context.info(
            "ORA analysis complete. Use create_visualization tool with plot_type='pathway_enrichment' to visualize results"
        )

    return EnrichmentInternalResult(
        method="ora",
        n_gene_sets=len(gene_sets),
        n_significant=sum(1 for p in adjusted_pvalues.values() if p is not None and p < 0.05),
        enrichment_scores=enrichment_scores,
        pvalues=pvalues,
        adjusted_pvalues=adjusted_pvalues,
        gene_set_statistics=gene_set_statistics,
        gene_set_summaries=create_gene_set_summaries(gene_sets),
        top_gene_sets=top_gene_sets,
        top_depleted_sets=[],  # ORA does not produce depleted gene sets
        method_specific_data={
            "query_genes": list(query_genes),
            "gene_sets_size_summary": {k: len(v) for k, v in gene_sets.items()},
            "pvalue_threshold": pvalue_threshold,
            "logfc_threshold": logfc_threshold,
            "min_size": min_size,
            "max_size": max_size,
        }
    )


async def perform_ssgsea(
    adata,
    gene_sets: Dict[str, List[str]],
    min_size: int = 10,
    max_size: int = 500,
    species: Optional[str] = None,
    database: Optional[str] = None,
    context=None,
) -> SSGSEAResult:
    """
    Perform single-sample Gene Set Enrichment Analysis (ssGSEA).

    This calculates enrichment scores for each sample independently.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene_sets : Dict[str, List[str]]
        Gene sets to test
    min_size : int
        Minimum gene set size
    max_size : int
        Maximum gene set size
    species : Optional[str]
        Species for the analysis (e.g., 'mouse', 'human')
    database : Optional[str]
        Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process')
    context : Optional
        MCP context

    Returns
    -------
    Dict containing enrichment results
    """
    is_available, error_msg = is_gseapy_available()
    if not is_available:
        raise ImportError(error_msg)

    import gseapy as gp

    if context:
        await context.info("Running ssGSEA analysis...")

    # Prepare expression data
    if hasattr(adata.X, "todense"):
        expr_df = pd.DataFrame(
            adata.X.todense().T, index=adata.var_names, columns=adata.obs_names
        )
    else:
        expr_df = pd.DataFrame(
            adata.X.T, index=adata.var_names, columns=adata.obs_names
        )

    # Run ssGSEA
    try:
        res = gp.ssgsea(
            data=expr_df,
            gene_sets=gene_sets,
            min_size=min_size,
            max_size=max_size,
            permutation_num=0,  # No permutation for ssGSEA
            no_plot=True,
            processes=1,
            seed=42,
        )

        # Extract results - ssGSEA stores enrichment scores in res.results
        if hasattr(res, "results") and isinstance(res.results, dict):
            # res.results is a dict where keys are sample names and values are DataFrames
            # We need to reorganize this into gene sets x samples format
            all_samples = list(res.results.keys())
            all_gene_sets = set()

            # Get all gene sets
            for sample_df in res.results.values():
                if isinstance(sample_df, pd.DataFrame) and "Term" in sample_df.columns:
                    all_gene_sets.update(sample_df["Term"].values)

            all_gene_sets = list(all_gene_sets)

            # Create scores matrix
            scores_matrix = pd.DataFrame(
                index=all_gene_sets, columns=all_samples, dtype=float
            )

            # Fill in scores
            for sample, df in res.results.items():
                if (
                    isinstance(df, pd.DataFrame)
                    and "Term" in df.columns
                    and "ES" in df.columns
                ):
                    for _, row in df.iterrows():
                        if row["Term"] in scores_matrix.index:
                            scores_matrix.loc[row["Term"], sample] = row["ES"]

            scores_df = scores_matrix.fillna(0)  # Fill missing values with 0
        else:
            # ssGSEA results format not recognized - fail honestly instead of returning empty results
            error_msg = (
                "ssGSEA results format not recognized. Expected 'res' attribute with enrichment scores matrix. "
                "This may indicate an issue with the gseapy installation, gene set database, or input data format. "
                "Please check your gene sets and ensure gseapy is properly configured."
            )
            logger.error(error_msg)
            raise ProcessingError(error_msg)

        # Calculate statistics across samples
        enrichment_scores = {}
        gene_set_statistics = {}

        if not scores_df.empty:
            for gs_name in scores_df.index:
                scores = scores_df.loc[gs_name].values
                enrichment_scores[gs_name] = float(np.mean(scores))

                gene_set_statistics[gs_name] = {
                    "mean_score": float(np.mean(scores)),
                    "std_score": float(np.std(scores)),
                    "min_score": float(np.min(scores)),
                    "max_score": float(np.max(scores)),
                    "size": len(gene_sets.get(gs_name, [])),
                }

            # Add scores to adata
            for gs_name in scores_df.index:
                adata.obs[f"ssgsea_{gs_name}"] = scores_df.loc[gs_name].values

            # Store gene set membership for validation
            adata.uns["enrichment_gene_sets"] = gene_sets

            # Store metadata for scientific provenance tracking
            obs_keys = [f"ssgsea_{gs_name}" for gs_name in scores_df.index]
            store_analysis_metadata(
                adata,
                analysis_name="enrichment_ssgsea",
                method="ssgsea",
                parameters={
                    "min_size": min_size,
                    "max_size": max_size,
                },
                results_keys={"obs": obs_keys, "uns": ["enrichment_gene_sets"]},
                statistics={
                    "n_gene_sets": len(gene_sets),
                    "n_samples": adata.n_obs,
                },
                species=species,
                database=database,
            )

        # Get top gene sets by mean enrichment
        sorted_by_mean = sorted(
            enrichment_scores.items(), key=lambda x: x[1], reverse=True
        )
        top_gene_sets = [x[0] for x in sorted_by_mean[:10]]

        # Convert enrichment scores to SSGSEAResult format (gene_set -> list of scores)
        enrichment_scores_lists = {}
        for gs_name in scores_df.index:
            enrichment_scores_lists[gs_name] = scores_df.loc[gs_name].values.tolist()
        
        # Convert gene_set_statistics to summary_stats format
        summary_stats = {}
        for gs_name, stats in gene_set_statistics.items():
            summary_stats[gs_name] = {
                "mean": stats["mean_score"],
                "std": stats["std_score"],
                "min": stats["min_score"],
                "max": stats["max_score"],
                "size": stats["size"]
            }

        return SSGSEAResult(
            method="ssgsea",
            n_gene_sets=len(gene_sets),
            n_samples=adata.n_obs,
            enrichment_scores=enrichment_scores_lists,
            gene_set_summaries=create_gene_set_summaries(gene_sets),
            summary_stats=summary_stats,
            score_matrix=scores_df.values if scores_df is not None else None,
            method_specific_data={
                "pvalues": {},  # ssGSEA doesn't provide p-values
                "adjusted_pvalues": {},
                "top_gene_sets": top_gene_sets,
                "top_depleted_sets": [],
                "scores_added_to_obs": True,
                "n_significant": len(gene_sets),  # All gene sets get scores in ssGSEA
            }
        )

    except Exception as e:
        logger.error(f"ssGSEA failed: {e}")
        raise


async def perform_enrichr(
    gene_list: List[str],
    gene_sets: Optional[str] = None,
    organism: str = "human",
    context=None,
) -> EnrichrResult:
    """
    Perform enrichment analysis using Enrichr web service.

    Parameters
    ----------
    gene_list : List[str]
        List of genes to analyze
    gene_sets : Optional[str]
        Enrichr library name. If None, use default libraries
    organism : str
        Organism ('human' or 'mouse')
    context : Optional
        MCP context

    Returns
    -------
    Dict containing enrichment results
    """
    is_available, error_msg = is_gseapy_available()
    if not is_available:
        raise ImportError(error_msg)

    import gseapy as gp

    if context:
        await context.info("Running Enrichr analysis...")

    # Default gene set libraries
    if gene_sets is None:
        gene_sets = [
            "GO_Biological_Process_2023",
            "GO_Molecular_Function_2023",
            "GO_Cellular_Component_2023",
            "KEGG_2021_Human" if organism == "human" else "KEGG_2019_Mouse",
            "Reactome_2022",
            "MSigDB_Hallmark_2020",
        ]
    elif isinstance(gene_sets, str):
        # Map user-friendly database name to actual Enrichr library name
        enrichr_library = map_gene_set_database_to_enrichr_library(gene_sets, organism)
        gene_sets = [enrichr_library]

    # Run Enrichr
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism.capitalize(),
            outdir=None,
            cutoff=0.05,
        )

        # Get results - enr.results is already a DataFrame
        all_results = enr.results

        # Prepare output
        enrichment_scores = {}
        pvalues = {}
        adjusted_pvalues = {}
        gene_set_statistics = {}

        for idx, row in all_results.iterrows():
            term = row["Term"]
            enrichment_scores[term] = row["Combined Score"]
            pvalues[term] = row["P-value"]
            adjusted_pvalues[term] = row["Adjusted P-value"]

            gene_set_statistics[term] = {
                "combined_score": row["Combined Score"],
                "pval": row["P-value"],
                "adjusted_pval": row["Adjusted P-value"],
                "z_score": row.get("Z-score", np.nan),
                "overlap": row["Overlap"],
                "genes": (
                    row["Genes"].split(";") if isinstance(row["Genes"], str) else []
                ),
            }

        # Get top results
        all_results_sorted = all_results.sort_values("Combined Score", ascending=False)
        top_gene_sets = all_results_sorted.head(10)["Term"].tolist()

        # Extract odds ratios from results
        odds_ratios = {}
        genes_found_in_results = []
        for idx, row in all_results.iterrows():
            term = row["Term"]
            odds_ratios[term] = row.get("Odds Ratio", 1.0)
            # Collect genes found in any enriched term
            if isinstance(row["Genes"], str):
                genes_found_in_results.extend(row["Genes"].split(";"))
        
        # Remove duplicates and get unique genes found
        genes_found_unique = list(set(genes_found_in_results))
        
        return EnrichrResult(
            method="enrichr", 
            n_gene_sets=len(all_results),
            n_significant=len(all_results[all_results["Adjusted P-value"] < 0.05]),
            input_gene_summary=create_input_gene_summary(gene_list, genes_found_unique),
            enrichment_scores=enrichment_scores,
            pvalues=pvalues,
            adjusted_pvalues=adjusted_pvalues,
            odds_ratios=odds_ratios,
            gene_set_statistics=gene_set_statistics,
            top_gene_sets=top_gene_sets,
            libraries_used=gene_sets if isinstance(gene_sets, list) else [gene_sets] if gene_sets else [],
            method_specific_data={
                "top_depleted_sets": [],
                "n_input_genes": len(gene_list),
                "organism": organism,
            }
        )

    except Exception as e:
        logger.error(f"Enrichr failed: {e}")
        raise


# ============================================================================
# Spatial Enrichment Analysis Functions (EnrichMap-based)
# ============================================================================


def is_enrichmap_available() -> Tuple[bool, str]:
    """Check if EnrichMap is available and all dependencies are met."""
    try:
        import enrichmap as em

        # Check for required dependencies
        # Map package names to their import names
        module_mapping = {
            "scanpy": "scanpy",
            "squidpy": "squidpy",
            "scipy": "scipy",
            "scikit-learn": "sklearn",
            "statsmodels": "statsmodels",
            "pygam": "pygam",
            "scikit-gstat": "skgstat",
            "adjustText": "adjustText",
            "splot": "splot",
        }

        missing = []
        for package, module in module_mapping.items():
            try:
                __import__(module)
            except ImportError:
                missing.append(package)

        if missing:
            return False, f"Missing EnrichMap dependencies: {', '.join(missing)}"

        return True, ""
    except ImportError:
        return False, "EnrichMap not installed. Install with: pip install enrichmap"


async def perform_spatial_enrichment(
    data_id: str,
    data_store: Dict[str, Any],
    gene_sets: Union[List[str], Dict[str, List[str]]],
    score_keys: Optional[Union[str, List[str]]] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    correct_spatial_covariates: bool = True,
    batch_key: Optional[str] = None,
    gene_weights: Optional[Dict[str, Dict[str, float]]] = None,
    species: str = "unknown",
    database: Optional[str] = None,
    context: Optional[Context] = None,
) -> EnrichmentInsights:
    """
    Perform spatially-aware gene set enrichment analysis using EnrichMap.

    Parameters
    ----------
    data_id : str
        Identifier for the spatial data in the data store
    data_store : Dict[str, Any]
        Dictionary containing the data
    gene_sets : Union[List[str], Dict[str, List[str]]]
        Either a single gene list or a dictionary of gene sets where keys are
        signature names and values are lists of genes
    score_keys : Optional[Union[str, List[str]]]
        Names for the gene signatures if gene_sets is a list. Ignored if gene_sets
        is already a dictionary
    spatial_key : str
        Key in adata.obsm containing spatial coordinates (default: "spatial")
    n_neighbors : int
        Number of nearest spatial neighbors for smoothing (default: 6)
    smoothing : bool
        Whether to perform spatial smoothing (default: True)
    correct_spatial_covariates : bool
        Whether to correct for spatial covariates using GAM (default: True)
    batch_key : Optional[str]
        Column in adata.obs for batch-wise normalization
    gene_weights : Optional[Dict[str, Dict[str, float]]]
        Pre-computed gene weights for each signature
    species : str
        Species for the analysis (e.g., 'mouse', 'human')
    database : Optional[str]
        Gene set database used (e.g., 'KEGG_Pathways', 'GO_Biological_Process')
    context : Optional[Context]
        Execution context

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - data_id: ID of the data with enrichment scores
        - signatures: List of computed signatures
        - score_columns: List of column names containing scores
        - gene_contributions: Dictionary of gene contributions per signature
        - summary_stats: Summary statistics for each signature
    """
    # Check if EnrichMap is available
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")

    # Import EnrichMap
    import enrichmap as em

    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

    # Validate spatial coordinates
    if spatial_key not in adata.obsm:
        raise ProcessingError(
            f"Spatial coordinates '{spatial_key}' not found in adata.obsm"
        )

    # Convert single gene list to dictionary format
    if isinstance(gene_sets, list):
        if score_keys is None:
            score_keys = "enrichmap_signature"
        gene_sets = {score_keys: gene_sets}

    # Validate gene sets with format conversion
    available_genes = set(adata.var_names)
    validated_gene_sets = {}

    # Debug: log total available genes
    logger.info(f"Total available genes in dataset: {len(available_genes)}")
    logger.info(f"First 10 genes: {list(available_genes)[:10]}")

    for sig_name, genes in gene_sets.items():
        # Try direct matching first
        common_genes = [gene for gene in genes if gene in available_genes]
        
        # If few matches and we know the species, try format conversion
        if len(common_genes) < len(genes) * 0.5 and species != "unknown":
            dataset_format_genes, conversion_map = _convert_gene_format_for_matching(
                genes, available_genes, species
            )
            
            if len(dataset_format_genes) > len(common_genes):
                # Format conversion helped, use dataset format genes for EnrichMap
                common_genes = dataset_format_genes
                if context:
                    await context.info(
                        f"Applied gene format conversion for '{sig_name}': "
                        f"{len(dataset_format_genes)}/{len(genes)} genes matched after conversion"
                    )
        
        logger.info(
            f"Checking signature '{sig_name}': requested {genes[:3]}... found {len(common_genes)}/{len(genes)}"
        )
        if len(common_genes) < 2:
            logger.warning(
                f"Signature '{sig_name}' has {len(common_genes)} genes in the dataset. Skipping."
            )
            continue
        validated_gene_sets[sig_name] = common_genes
        logger.info(
            f"Signature '{sig_name}': {len(common_genes)}/{len(genes)} genes found"
        )

    if not validated_gene_sets:
        # Collect diagnostic information for better error reporting
        sample_dataset_genes = list(available_genes)[:10]
        sample_requested_genes = []
        for sig_name, genes in gene_sets.items():
            sample_requested_genes.extend(genes[:3])
        sample_requested_genes = list(set(sample_requested_genes))[:10]
        
        error_msg = (
            f"No valid gene signatures found with at least 2 genes.\n\n"
            f"Diagnostic information:\n"
            f"• Dataset contains {len(available_genes)} genes\n"
            f"• Sample dataset genes: {sample_dataset_genes}\n"
            f"• Sample requested genes: {sample_requested_genes}\n"
            f"• Gene signatures provided: {list(gene_sets.keys())}\n\n"
            f"Common causes and solutions:\n"
            f"1. Species mismatch (human vs mouse gene naming)\n"
            f"   - Check if your data species matches the gene set database\n"
            f"2. Gene name format differences (e.g., 'CD3D' vs 'Cd3d')\n"
            f"   - Ensure gene names match your dataset's format\n"
            f"3. Outdated or incompatible gene sets\n"
            f"   - Try a different gene_set_database option\n"
            f"4. Custom gene sets with incorrect gene symbols\n"
            f"   - Verify gene symbols are current and correct"
        )
        raise ProcessingError(error_msg)

    # Run EnrichMap scoring - process each gene set individually
    failed_signatures = []
    successful_signatures = []
    
    for sig_name, genes in validated_gene_sets.items():
        try:
            if context:
                await context.info(f"Processing gene set '{sig_name}' with {len(genes)} genes")
            
            em.tl.score(
                adata=adata,
                gene_set=genes,  # Fixed: use gene_set (correct API parameter name)
                score_key=sig_name,  # Fixed: provide explicit score_key
                spatial_key=spatial_key,
                n_neighbors=n_neighbors,
                smoothing=smoothing,
                correct_spatial_covariates=correct_spatial_covariates,
                batch_key=batch_key,
            )
            successful_signatures.append(sig_name)
            
        except Exception as e:
            logger.error(f"EnrichMap failed for '{sig_name}': {e}")
            failed_signatures.append((sig_name, str(e)))
    
    # Check if any signatures were processed successfully
    if not successful_signatures:
        error_details = "; ".join([f"{name}: {error}" for name, error in failed_signatures])
        raise ProcessingError(
            f"All EnrichMap scoring failed. This may indicate:\n"
            f"1. EnrichMap package installation issues\n"
            f"2. Incompatible gene names or data format\n"
            f"3. Insufficient spatial information\n"
            f"Details: {error_details}"
        )
    
    # Update validated_gene_sets to only include successful ones
    validated_gene_sets = {sig: validated_gene_sets[sig] for sig in successful_signatures}
    
    if context and failed_signatures:
        await context.warning(f"Failed to process {len(failed_signatures)} gene sets: {[name for name, _ in failed_signatures]}")

    # Collect results
    score_columns = [f"{sig}_score" for sig in validated_gene_sets.keys()]

    # Calculate summary statistics
    summary_stats = {}
    for sig_name in validated_gene_sets.keys():
        score_col = f"{sig_name}_score"
        scores = adata.obs[score_col]

        summary_stats[sig_name] = {
            "mean": float(scores.mean()),
            "std": float(scores.std()),
            "min": float(scores.min()),
            "max": float(scores.max()),
            "median": float(scores.median()),
            "q25": float(scores.quantile(0.25)),
            "q75": float(scores.quantile(0.75)),
            "n_genes": len(validated_gene_sets[sig_name]),
        }

    # Get gene contributions
    gene_contributions = {}
    if "gene_contributions" in adata.uns:
        gene_contributions = {
            sig: {gene: float(contrib.mean()) for gene, contrib in contribs.items()}
            for sig, contribs in adata.uns["gene_contributions"].items()
        }

    # Store gene set membership for validation
    adata.uns["enrichment_gene_sets"] = validated_gene_sets

    # Store metadata for scientific provenance tracking
    store_analysis_metadata(
        adata,
        analysis_name="enrichment_spatial",
        method="spatial_enrichmap",
        parameters={
            "spatial_key": spatial_key,
            "n_neighbors": n_neighbors,
            "smoothing": smoothing,
            "correct_spatial_covariates": correct_spatial_covariates,
            "batch_key": batch_key,
        },
        results_keys={
            "obs": score_columns,
            "uns": ["gene_contributions", "enrichment_gene_sets"],
        },
        statistics={
            "n_gene_sets": len(validated_gene_sets),
            "n_successful_signatures": len(successful_signatures),
            "n_failed_signatures": len(failed_signatures),
        },
        species=species,
        database=database,
    )

    # Inform user about visualization options
    if context:
        await context.info(
            "Spatial enrichment analysis complete. Use create_visualization tool with plot_type='enrichment' to visualize results"
        )

    # === INSIGHTS-BASED RESULT CREATION ===
    # Note: Skip creating data for all pathways - only create for top findings

    # === TRANSFORM TO INSIGHTS-BASED RESULT ===
    
    # Sort signatures by their max scores and limit to top findings
    sorted_sigs = sorted(
        summary_stats.keys(), key=lambda x: summary_stats[x]["max"], reverse=True
    )
    
    # Create PathwayResult objects for top pathways only
    top_pathways = []
    max_pathways = min(15, len(sorted_sigs))  # Limit to top 15 pathways
    
    for sig_name in sorted_sigs[:max_pathways]:
        stats = summary_stats[sig_name]
        genes = validated_gene_sets.get(sig_name, [])
        
        pathway_result = PathwayResult(
            name=clean_pathway_name(sig_name),
            category=categorize_pathway(sig_name),
            significance="Spatial Enrichment Score",  # Spatial analysis doesn't have p-values
            effect_size=format_effect_size(stats["max"], baseline=1.0),
            key_genes=genes[:2],  # Only keep 2 representative genes
            description=f"Enrichment score range: {stats['min']:.2f} to {stats['max']:.2f}"
        )
        top_pathways.append(pathway_result)
    
    # Group pathways by category
    pathway_categories = {}
    for pathway in top_pathways:
        if pathway.category not in pathway_categories:
            pathway_categories[pathway.category] = CategorySummary(
                count=0,
                significance_level="Variable",
                top_themes=[]
            )
        pathway_categories[pathway.category].count += 1
    
    # Set significance levels based on count
    for category, summary in pathway_categories.items():
        if summary.count >= 4:
            summary.significance_level = "High"
        elif summary.count >= 2:
            summary.significance_level = "Moderate"
        else:
            summary.significance_level = "Low"
        
        # Extract themes from pathway names in this category
        category_pathways = [p for p in top_pathways if p.category == category]
        summary.top_themes = [p.name for p in category_pathways[:3]]
    
    # Collect key genes from top pathways
    all_key_genes = []
    for pathway in top_pathways[:5]:  # Only from top 5 pathways
        all_key_genes.extend(pathway.key_genes)
    key_genes = list(dict.fromkeys(all_key_genes))[:8]  # Remove duplicates, keep top 8
    
    # Generate insights
    summary = generate_insights_summary(top_pathways)
    interpretation = generate_pathway_interpretation(top_pathways)
    suggestions = suggest_followup_analyses(top_pathways)
    
    # Determine overall significance
    if len(top_pathways) >= 10:
        significance_overview = f"High significance with {len(top_pathways)} enriched pathways across {len(pathway_categories)} categories"
    elif len(top_pathways) >= 5:
        significance_overview = f"Moderate significance with {len(top_pathways)} enriched pathways"
    else:
        significance_overview = f"Limited significance with {len(top_pathways)} enriched pathways"
    
    return EnrichmentInsights(
        summary=summary,
        significance_overview=significance_overview,
        top_pathways=top_pathways,
        pathway_categories=pathway_categories,
        key_genes=key_genes,
        suggested_analyses=suggestions,
        interpretation=interpretation,
        method="spatial_enrichmap",
        total_pathways_tested=len(validated_gene_sets),
        significant_count=len(top_pathways),
        stats_summary={
            "max_enrichment_score": max(stats["max"] for stats in summary_stats.values()) if summary_stats else 0,
            "mean_enrichment_score": np.mean([stats["max"] for stats in summary_stats.values()]) if summary_stats else 0,
            "spatial_parameters": {
                "n_neighbors": n_neighbors,
                "smoothing": smoothing,
                "correct_spatial_covariates": correct_spatial_covariates
            }
        }
    )


async def compute_spatial_metrics(
    data_id: str,
    data_store: Dict[str, Any],
    score_key: str,
    metrics: Optional[List[str]] = None,
    n_neighbors: int = 6,
    n_perms: int = 999,
    context: Optional[Context] = None,
) -> SpatialMetricResult:
    """
    Compute spatial metrics for enrichment scores.

    Parameters
    ----------
    data_id : str
        Identifier for the spatial data
    data_store : Dict[str, Any]
        Dictionary containing the data
    score_key : str
        Column name containing the enrichment scores
    metrics : Optional[List[str]]
        List of metrics to compute. Options: ['morans_i', 'getis_ord', 'variance']
        If None, computes all metrics
    n_neighbors : int
        Number of spatial neighbors (default: 6)
    n_perms : int
        Number of permutations for significance testing (default: 999)
    context : Optional[Context]
        Execution context

    Returns
    -------
    Dict[str, Any]
        Dictionary containing computed spatial metrics
    """
    # Check if EnrichMap is available
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")

    import enrichmap as em

    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

    # Validate score column
    if score_key not in adata.obs.columns:
        raise ProcessingError(f"Score column '{score_key}' not found in adata.obs")

    # Default metrics
    if metrics is None:
        metrics = ["morans_i", "getis_ord", "variance"]

    try:
        # Compute spatial metrics
        result = em.tl.compute_spatial_metrics(
            adata=adata,
            score_keys=[score_key],
            metrics=metrics,
            n_neighs=n_neighbors,
            n_perms=n_perms,
        )

        # Extract results for the single score key
        metric_results = {}
        for metric in metrics:
            if metric in result:
                metric_results[metric] = {
                    "value": float(result[metric][score_key]),
                    "p_value": float(
                        result.get(f"{metric}_pval", {}).get(score_key, np.nan)
                    ),
                }

        return SpatialMetricResult(
            data_id=data_id,
            score_key=score_key,
            metrics=metric_results,
            parameters={"n_neighbors": n_neighbors, "n_perms": n_perms},
        )

    except Exception as e:
        raise ProcessingError(f"Spatial metrics computation failed: {str(e)}")


async def cluster_gene_correlation(
    data_id: str,
    data_store: Dict[str, Any],
    signature_name: str,
    cluster_key: str = "leiden",
    correlation_method: str = "pearson",
    context: Optional[Context] = None,
) -> ClusterCorrelationResult:
    """
    Compute correlation between gene expression and enrichment scores per cluster.

    Parameters
    ----------
    data_id : str
        Identifier for the spatial data
    data_store : Dict[str, Any]
        Dictionary containing the data
    signature_name : str
        Name of the signature to analyze
    cluster_key : str
        Column in adata.obs containing cluster labels (default: "leiden")
    correlation_method : str
        Correlation method: 'pearson' or 'spearman' (default: "pearson")
    context : Optional[Context]
        Execution context

    Returns
    -------
    Dict[str, Any]
        Dictionary containing correlation results per cluster
    """
    # Check if EnrichMap is available
    is_available, error_msg = is_enrichmap_available()
    if not is_available:
        raise ProcessingError(f"EnrichMap is not available: {error_msg}")

    import enrichmap as em

    # Get data
    if data_id not in data_store:
        raise ProcessingError(f"Data '{data_id}' not found in data store")

    adata = data_store[data_id]["adata"]

    # Validate inputs
    score_key = f"{signature_name}_score"
    if score_key not in adata.obs.columns:
        raise ProcessingError(
            f"Score column '{score_key}' not found. Run enrichment analysis first."
        )

    if cluster_key not in adata.obs.columns:
        raise ProcessingError(f"Cluster column '{cluster_key}' not found in adata.obs")

    try:
        # Compute cluster-gene correlations
        result = em.tl.cluster_gene_correlation(
            adata=adata,
            signature_name=signature_name,
            cluster_key=cluster_key,
            correlation_method=correlation_method,
        )

        # Convert results to serializable format
        correlation_results = {}
        for cluster, corr_df in result.items():
            correlation_results[str(cluster)] = {
                "genes": corr_df.index.tolist(),
                "correlations": corr_df["correlation"].tolist(),
                "top_positive": corr_df.nlargest(10, "correlation")[
                    ["correlation"]
                ].to_dict(),
                "top_negative": corr_df.nsmallest(10, "correlation")[
                    ["correlation"]
                ].to_dict(),
            }

        return ClusterCorrelationResult(
            data_id=data_id,
            signature_name=signature_name,
            cluster_key=cluster_key,
            correlation_method=correlation_method,
            cluster_correlations=correlation_results,
        )

    except Exception as e:
        raise ProcessingError(f"Cluster gene correlation failed: {str(e)}")


# ============================================================================
# Gene Set Loading Functions
# ============================================================================


class GeneSetLoader:
    """Load gene sets from various databases and sources."""

    def __init__(self, species: str = "human"):
        """
        Initialize gene set loader.

        Parameters
        ----------
        species : str
            Species for gene sets ('human' or 'mouse')
        """
        self.species = species.lower()
        self.organism = "Homo sapiens" if species == "human" else "Mus musculus"

    def load_msigdb(
        self,
        collection: str = "H",
        subcollection: Optional[str] = None,
        min_size: int = 10,
        max_size: int = 500,
    ) -> Dict[str, List[str]]:
        """
        Load gene sets from MSigDB using gseapy.

        Parameters
        ----------
        collection : str
            MSigDB collection name:
            - H: hallmark gene sets
            - C1: positional gene sets
            - C2: curated gene sets (e.g., CGP, CP:KEGG, CP:REACTOME)
            - C3: motif gene sets
            - C4: computational gene sets
            - C5: GO gene sets (CC, BP, MF)
            - C6: oncogenic signatures
            - C7: immunologic signatures
            - C8: cell type signatures
        subcollection : Optional[str]
            Subcollection for specific databases (e.g., 'CP:KEGG', 'GO:BP')
        min_size : int
            Minimum gene set size
        max_size : int
            Maximum gene set size

        Returns
        -------
        Dict[str, List[str]]
            Dictionary of gene sets
        """
        try:
            import gseapy as gp

            # Get available gene sets
            gene_sets_dict = {}

            if collection == "H":
                # Hallmark gene sets
                gene_sets = gp.get_library_name(organism=self.organism)
                if "MSigDB_Hallmark_2020" in gene_sets:
                    gene_sets_dict = gp.get_library(
                        "MSigDB_Hallmark_2020", organism=self.organism
                    )

            elif collection == "C2" and subcollection == "CP:KEGG":
                # KEGG pathways
                if self.species == "human":
                    gene_sets_dict = gp.get_library(
                        "KEGG_2021_Human", organism=self.organism
                    )
                else:
                    gene_sets_dict = gp.get_library(
                        "KEGG_2019_Mouse", organism=self.organism
                    )

            elif collection == "C2" and subcollection == "CP:REACTOME":
                # Reactome pathways
                gene_sets_dict = gp.get_library("Reactome_2022", organism=self.organism)

            elif collection == "C5":
                # GO gene sets
                if subcollection == "GO:BP" or subcollection is None:
                    gene_sets_dict.update(
                        gp.get_library(
                            "GO_Biological_Process_2023", organism=self.organism
                        )
                    )
                if subcollection == "GO:MF" or subcollection is None:
                    gene_sets_dict.update(
                        gp.get_library(
                            "GO_Molecular_Function_2023", organism=self.organism
                        )
                    )
                if subcollection == "GO:CC" or subcollection is None:
                    gene_sets_dict.update(
                        gp.get_library(
                            "GO_Cellular_Component_2023", organism=self.organism
                        )
                    )

            elif collection == "C8":
                # Cell type signatures
                gene_sets_dict = gp.get_library(
                    "CellMarker_Augmented_2021", organism=self.organism
                )

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets_dict.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(
                f"Loaded {len(filtered_sets)} gene sets from MSigDB {collection}"
            )
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load MSigDB gene sets: {e}")
            return {}

    def load_go_terms(
        self, aspect: str = "BP", min_size: int = 10, max_size: int = 500
    ) -> Dict[str, List[str]]:
        """
        Load GO terms using gseapy.

        Parameters
        ----------
        aspect : str
            GO aspect: 'BP' (biological process), 'MF' (molecular function), 'CC' (cellular component)
        min_size : int
            Minimum gene set size
        max_size : int
            Maximum gene set size

        Returns
        -------
        Dict[str, List[str]]
            Dictionary of GO gene sets
        """
        aspect_map = {
            "BP": "GO_Biological_Process_2023",
            "MF": "GO_Molecular_Function_2023",
            "CC": "GO_Cellular_Component_2023",
        }

        if aspect not in aspect_map:
            raise ValueError(f"Invalid GO aspect: {aspect}")

        try:
            import gseapy as gp

            gene_sets = gp.get_library(aspect_map[aspect], organism=self.organism)

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(f"Loaded {len(filtered_sets)} GO {aspect} gene sets")
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load GO gene sets: {e}")
            return {}

    def load_kegg_pathways(
        self, min_size: int = 10, max_size: int = 500
    ) -> Dict[str, List[str]]:
        """Load KEGG pathways."""
        try:
            import gseapy as gp

            if self.species == "human":
                gene_sets = gp.get_library("KEGG_2021_Human", organism=self.organism)
            else:
                gene_sets = gp.get_library("KEGG_2019_Mouse", organism=self.organism)

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(f"Loaded {len(filtered_sets)} KEGG pathways")
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load KEGG pathways: {e}")
            return {}

    def load_reactome_pathways(
        self, min_size: int = 10, max_size: int = 500
    ) -> Dict[str, List[str]]:
        """Load Reactome pathways."""
        try:
            import gseapy as gp

            gene_sets = gp.get_library("Reactome_2022", organism=self.organism)

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(f"Loaded {len(filtered_sets)} Reactome pathways")
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load Reactome pathways: {e}")
            return {}

    def load_cell_markers(
        self, min_size: int = 5, max_size: int = 200
    ) -> Dict[str, List[str]]:
        """Load cell type marker gene sets."""
        try:
            import gseapy as gp

            gene_sets = gp.get_library(
                "CellMarker_Augmented_2021", organism=self.organism
            )

            # Filter by size
            filtered_sets = {}
            for name, genes in gene_sets.items():
                if min_size <= len(genes) <= max_size:
                    filtered_sets[name] = genes

            logger.info(f"Loaded {len(filtered_sets)} cell type marker sets")
            return filtered_sets

        except Exception as e:
            logger.error(f"Failed to load cell markers: {e}")
            return {}


async def load_gene_sets(
    database: str,
    species: str = "human",
    min_genes: int = 10,
    max_genes: int = 500,
    context=None,
) -> Dict[str, List[str]]:
    """
    Load gene sets from specified database.

    Parameters
    ----------
    database : str
        Database name:
        - GO_Biological_Process, GO_Molecular_Function, GO_Cellular_Component
        - KEGG_Pathways
        - Reactome_Pathways
        - MSigDB_Hallmark
        - Cell_Type_Markers
    species : str
        Species ('human' or 'mouse')
    min_genes : int
        Minimum gene set size
    max_genes : int
        Maximum gene set size
    context : Optional
        MCP context for logging

    Returns
    -------
    Dict[str, List[str]]
        Dictionary of gene sets
    """
    loader = GeneSetLoader(species=species)

    database_map = {
        "GO_Biological_Process": lambda: loader.load_go_terms(
            "BP", min_genes, max_genes
        ),
        "GO_Molecular_Function": lambda: loader.load_go_terms(
            "MF", min_genes, max_genes
        ),
        "GO_Cellular_Component": lambda: loader.load_go_terms(
            "CC", min_genes, max_genes
        ),
        "KEGG_Pathways": lambda: loader.load_kegg_pathways(min_genes, max_genes),
        "Reactome_Pathways": lambda: loader.load_reactome_pathways(
            min_genes, max_genes
        ),
        "MSigDB_Hallmark": lambda: loader.load_msigdb("H", None, min_genes, max_genes),
        "Cell_Type_Markers": lambda: loader.load_cell_markers(min_genes, max_genes),
    }

    if database not in database_map:
        raise ValueError(
            f"Unknown database: {database}. Available: {list(database_map.keys())}"
        )

    if context:
        await context.info(f"Loading gene sets from {database} for {species}")

    gene_sets = database_map[database]()

    if context:
        await context.info(f"Loaded {len(gene_sets)} gene sets from {database}")

    return gene_sets
