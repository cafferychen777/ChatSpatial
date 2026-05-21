#!/usr/bin/env python3
"""
Build all Supplementary Tables (S1-S4) as standalone CSV files.

S1: Sample-level effect sizes for OSCC enrichment
S2: Integrated method library
S3: Test scenarios
S4: AI agent comparison

Requires: compute_case_study_stats.py must be run first to generate
  data/Case2_OSCC/sample_level_tests.csv and effect_size_summary.csv
"""

import csv
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_DIR = SCRIPT_DIR.parent
OUT_DIR = PAPER_DIR / "supplementary_tables"
OUT_DIR.mkdir(exist_ok=True)
CASE_DIR = PAPER_DIR / "data" / "Case2_OSCC"


# ═══════════════════════════════════════════════════════════════════
# S2: Integrated method library
# ═══════════════════════════════════════════════════════════════════

METHODS = [
    # (Category, Method, Version, Citation)
    ("Data Loading & Preprocessing", "Scanpy", "1.9.3", "Wolf et al., 2018"),
    ("Data Loading & Preprocessing", "SpatialData", "0.2.0", "Marconato et al., 2024"),
    ("Data Loading & Preprocessing", "SCTransform", "0.4.0", "Hafemeister & Satija, 2019"),
    ("Multi-sample Integration", "Harmony", "0.0.9+", "Korsunsky et al., 2019"),
    ("Multi-sample Integration", "scVI", "1.0.0+", "Lopez et al., 2018"),
    ("Multi-sample Integration", "BBKNN", "1.5.0+", "Polanski et al., 2020"),
    ("Multi-sample Integration", "Scanorama", "1.7.0+ (optional)", "Hie et al., 2019"),
    ("Spatial Registration", "PASTE", "1.3.0", "Zeira et al., 2022"),
    ("Spatial Registration", "STalign", "2.0.0", "Clifton et al., 2023"),
    ("Visualization", "Squidpy", "1.3.0", "Palla et al., 2022"),
    ("Visualization", "Matplotlib / Seaborn", "3.7.1 / 0.12.2", "Hunter, 2007"),
    ("Cell Type Annotation", "Tangram", "1.0.3", "Biancalani et al., 2021"),
    ("Cell Type Annotation", "scANVI", "1.0.3", "Xu et al., 2021"),
    ("Cell Type Annotation", "CellAssign", "0.1.0", "Zhang et al., 2019"),
    ("Cell Type Annotation", "mLLMCellType", "1.0.0", "Yang et al., 2025*"),
    ("Cell Type Annotation", "scType", "via R", "Ianevski et al., 2022"),
    ("Cell Type Annotation", "SingleR", "via R", "Aran et al., 2019"),
    ("Differential Expression", "Wilcoxon", "via scanpy", "Wolf et al., 2018"),
    ("Differential Expression", "t-test", "via scanpy", "Wolf et al., 2018"),
    ("Differential Expression", "logreg", "via scanpy", "Wolf et al., 2018"),
    ("Spatial Domain Identification", "SpaGCN", "1.2.5", "Hu et al., 2021"),
    ("Spatial Domain Identification", "STAGATE", "1.1.0", "Dong et al., 2022"),
    ("Spatial Domain Identification", "GraphST", "1.1.0+", "Long et al., 2023"),
    ("Spatial Domain Identification", "Leiden", "via scanpy", "Traag et al., 2019"),
    ("Spatial Domain Identification", "Louvain", "via scanpy", "Blondel et al., 2008"),
    ("Spatially Variable Genes", "SPARK-X", "1.3.0", "Sun et al., 2020"),
    ("Spatially Variable Genes", "FlashS", "0.1.0", "Yang et al., 2026*"),
    ("Spatially Variable Genes", "SpatialDE", "1.1.3", "Svensson et al., 2018"),
    ("Cell-Cell Communication", "LIANA+", "1.0.0", "Dimitrov et al., 2022"),
    ("Cell-Cell Communication", "CellPhoneDB", "via LIANA", "Efremova et al., 2020"),
    ("Cell-Cell Communication", "CellChat", "native R", "Jin et al., 2021"),
    ("Cell-Cell Communication", "FastCCC", "0.1.0", "Hou et al., 2025*"),
    ("Deconvolution", "FlashDeconv", "0.1.0", "Yang et al., 2025*"),
    ("Deconvolution", "RCTD", "2.2.0", "Cable et al., 2022"),
    ("Deconvolution", "Cell2location", "0.1.3", "Kleshchevnikov et al., 2022"),
    ("Deconvolution", "SPOTlight", "0.1.7", "Elosua-Bayes et al., 2021"),
    ("Deconvolution", "Stereoscope", "0.3.0", "Andersson et al., 2020"),
    ("Deconvolution", "DestVI", "0.1.0", "Lopez et al., 2022"),
    ("Deconvolution", "Tangram (mapping)", "1.0.3", "Biancalani et al., 2021"),
    ("Deconvolution", "CARD", "2.0.0", "Ma & Zhou, 2022"),
    ("CNV Analysis", "infercnvpy", "0.6.1", "Weiler et al., 2025; Patel et al., 2014"),
    ("CNV Analysis", "Numbat", "1.3.2", "Gao et al., 2023"),
    ("Trajectory Inference", "CellRank", "2.0.0", "Lange et al., 2022"),
    ("Trajectory Inference", "Palantir", "1.3.0", "Setty et al., 2019"),
    ("Trajectory Inference", "DPT", "via scanpy", "Haghverdi et al., 2016"),
    ("RNA Velocity", "scVelo", "0.2.5", "Bergen et al., 2020"),
    ("RNA Velocity", "VeloVI", "1.0.0", "Gayoso et al., 2023"),
    ("Enrichment Analysis", "GSEApy", "1.0.6", "Subramanian et al., 2005"),
    ("Enrichment Analysis", "Pathway GSEA", "via GSEApy", "Subramanian et al., 2005"),
    ("Enrichment Analysis", "Pathway ORA", "via GSEApy", "Subramanian et al., 2005"),
    ("Enrichment Analysis", "ssGSEA", "via GSEApy", "Barbie et al., 2009"),
    ("Enrichment Analysis", "EnrichMap", "0.1.10", "Merico et al., 2010"),
    ("Spatial Statistics", "Moran's I", "via squidpy", "Palla et al., 2022"),
    ("Spatial Statistics", "Local Moran's I (LISA)", "via squidpy", "Palla et al., 2022"),
    ("Spatial Statistics", "Geary's C", "via squidpy", "Palla et al., 2022"),
    ("Spatial Statistics", "Getis-Ord Gi*", "via esda/PySAL", "Palla et al., 2022"),
    ("Spatial Statistics", "Neighborhood Enrichment", "via squidpy", "Palla et al., 2022"),
    ("Spatial Statistics", "Co-occurrence", "via squidpy", "Palla et al., 2022"),
    ("Spatial Statistics", "Ripley's K/L", "via squidpy", "Palla et al., 2022"),
    ("Spatial Statistics", "Centrality Scores", "via squidpy", "Palla et al., 2022"),
    ("Spatial Statistics", "Bivariate Moran's I", "4.6.0", "Palla et al., 2022"),
    ("Spatial Statistics", "Join Count Statistics", "via esda/PySAL", "Palla et al., 2022"),
    ("Spatial Statistics", "Local Join Count", "via esda/PySAL", "Palla et al., 2022"),
    ("Spatial Statistics", "Network Properties", "via squidpy", "Palla et al., 2022"),
    ("Spatial Statistics", "Spatial Centrality", "2.6.0", "Palla et al., 2022"),
]


# ═══════════════════════════════════════════════════════════════════
# S3: Test scenarios
# ═══════════════════════════════════════════════════════════════════

SCENARIOS = [
    # (ID, Category, Scenario, Input Data, Command/Goal, Expected Outcome, Result, Notes)
    (1, "Data Handling & Preprocessing", "Standard Visium Input", "Visium, 4992 spots, h5ad", "Load and preprocess the Visium data.", "Data loaded, QC metrics calculated.", "Pass", "Standard workflow."),
    (2, "Data Handling & Preprocessing", "High-res Xenium Input", "Xenium, 150k cells", "Load Xenium data and create AnnData.", "Data loaded, valid AnnData created.", "Pass", "Non-h5ad input."),
    (3, "Data Handling & Preprocessing", "Low-quality Data", "Visium, <500 genes/spot", "Filter low-quality spots.", "Low-quality spots removed.", "Pass", "QC validation."),
    (4, "Data Handling & Preprocessing", "Multi-sample Integration", "3 Visium slides", "Integrate samples with Harmony.", "Batch effects reduced.", "Pass", "Integration test."),
    (5, "Data Handling & Preprocessing", "Missing Spatial Coords", "No spatial coordinates", "Plot spatial MALAT1.", 'Error: "No spatial coordinates."', "Pass", "Error handling."),
    (6, "Core Spatial Analysis", "Basic Clustering", "DLPFC Visium (spatialLIBD)", "Identify spatial domains using SpaGCN.", "Anatomical layers identified.", "Pass", "Core function."),
    (7, "Core Spatial Analysis", "High-density Clustering", "MERFISH, 50k cells", "Use GraphST for fine domains.", "Sub-clusters identified.", "Pass", "Dense data test."),
    (8, "Core Spatial Analysis", "Spatial Variable Genes", "Mouse brain", "Find top 50 SVGs with SPARK-X.", "Known markers ranked high.", "Pass", "SVG detection."),
    (9, "Core Spatial Analysis", "Cell Deconvolution", "Visium + scRNA ref", "Deconvolve cell types using RCTD.", "Cell proportions generated.", "Pass", "Deconvolution."),
    (10, "Core Spatial Analysis", "Cell Annotation", "Xenium + atlas", "Annotate cells with Tangram.", "Major types assigned.", "Pass", "Annotation test."),
    (11, "Core Spatial Analysis", "Trajectory Analysis", "Dev. tissue data", "Infer trajectory using CellRank.", "Biological pathway identified.", "Pass", "Advanced test."),
    (12, "Core Spatial Analysis", "RNA Velocity", "Visium + spliced counts", "Calculate spatial RNA velocity.", "Dynamics visualized.", "Pass", "Complex analysis."),
    (13, "Core Spatial Analysis", "Spatial Statistics", "Mouse brain Visium", "Compute Moran's I for top marker genes.", "Spatial autocorrelation scores generated.", "Pass", "Statistics test."),
    (14, "Core Spatial Analysis", "Pathway Enrichment", "Cancer Visium", "Run GSEA on cluster markers using GO database.", "Enriched pathways identified.", "Pass", "Enrichment test."),
    (15, "Core Spatial Analysis", "CNV Detection", "Tumor Visium + ref", "Infer CNV using infercnvpy with immune cells as reference.", "CNV profiles generated.", "Pass", "CNV test."),
    (16, "Core Spatial Analysis", "Spatial Registration", "Serial sections", "Align 3 consecutive slices with PASTE.", "Slices registered to common space.", "Pass", "Registration test."),
    (17, "Conversational Workflows", "Iterative Analysis", "Tumor microenv.", "Find clusters. Now, markers for cluster 3.", "Context preserved.", "Pass", "Context test."),
    (18, "Conversational Workflows", "Ambiguous Command", "Any dataset", "Analyze the data.", "Clarification requested.", "Pass", "Proactive help."),
    (19, "Conversational Workflows", "Parameter Override", "Visium data", "Run STAGATE with resolution 1.5.", "User params used.", "Pass", "User control."),
    (20, "Conversational Workflows", "Cross-category Flow", "Cancer + ref", "Cluster, then find tumor-immune comm.", "Sequential execution.", "Pass", "Workflow test."),
    (21, "Conversational Workflows", "Error Recovery", "Large dataset", "Run cell communication analysis.", "Alternatives suggested.", "Pass", "Error handling."),
    (22, "Scalability & Stress Tests", "Large Dataset", "100k+ Visium HD", "Run full analysis pipeline.", "Pipeline completed.", "Pass", "Scale."),
    (23, "Scalability & Stress Tests", "Small Dataset", "300-spot STARmap", "Cluster the data.", "Runs despite noise.", "Pass", "Edge case."),
    (24, "Scalability & Stress Tests", "Sparse Genes", ">99.9% sparsity", "Find variable genes.", "Warning issued.", "Pass", "Robust."),
    (25, "Scalability & Stress Tests", "Noisy Reference", "Ambiguous scRNA ref", "Annotate using noisy ref.", "Low confidence scores.", "Pass", "Quality aware."),
    (26, "Scalability & Stress Tests", "Concurrent Users", "2 users, diff. data", "U1: Cluster brain. U2: Annotate tumor.", "No interference.", "Pass", "Isolation test."),
    (27, "Scalability & Stress Tests", "Long Context", "20+ commands", "Result of first clustering?", "History retrieved.", "Pass", "Memory test."),
    (28, "Scalability & Stress Tests", "Invalid Input", "Non-scientific cmd", "Write Python code to download a file from the internet.", "Polite refusal; redirected to spatial analysis scope.", "Pass", "Scope boundary test."),
    (29, "Known Limitations", "Unsupported method", "Visium data", "Run BayesSpace for spatial domains.", "Tool invoked but method not in enum; error returned listing available methods.", "Partial", "Method outside curated set."),
    (30, "Known Limitations", "Multi-step biological reasoning", "Tumor Visium", "Which subclone is most likely to metastasize?", "System identified relevant analyses but could not synthesize a biological conclusion; deferred to researcher.", "Partial", "Requires domain judgment beyond orchestration."),
    (31, "Known Limitations", "Very large Visium HD", "500k+ spots", "Run SpaGCN spatial domains.", "SpaGCN timed out (600s limit); system suggested Leiden as scalable alternative.", "Partial", "Algorithm-level scalability limit, not framework limit."),
]


# ═══════════════════════════════════════════════════════════════════
# S4: AI agent comparison
# ═══════════════════════════════════════════════════════════════════

AI_COMPARISON = [
    # (Feature, ChatSpatial, STAgent, SpatialAgent)
    (
        "Primary Philosophy",
        "Schema-Enforced Orchestrator. Reliably executes user-defined analytical strategies with scientist in full control.",
        "Interactive Analysis Assistant. User-guided exploration via Streamlit UI with literature synthesis and conflict detection.",
        "Code-Generation Agent. Supports both autonomous single-prompt and interactive multi-turn modes.",
    ),
    (
        "Planning & Execution",
        "Schema-Enforced Tool-calling via MCP.",
        "Hybrid tool-calling and code generation; Squidpy-specific RAG for API code.",
        "LLM-generated Python code calling injected tool functions; 17 skill templates guide workflow selection, with LLM-based (default) or embedding-based tool retrieval.",
    ),
    (
        "Reproducibility mechanism",
        "Schema-constrained parameter validation (Pydantic); measured in ablation.",
        "Tool-generated code with RAG assistance; conflict logging tracks claims against literature.",
        "Skill templates guide planning; LLM-generated code with JSONL observation logging.",
    ),
    (
        "Tool Scope",
        "Cross-Ecosystem (Python & R). 65 methods via 20 MCP tools.",
        "Python-only. ~19 tools (Scanpy/Squidpy core; STAligner, Tangram extensions).",
        "Primarily Python. 72 tools spanning spatial analysis, deconvolution, CCC, trajectory, and literature search; R scripts callable via bash.",
    ),
    (
        "Extensibility",
        "Open standard (MCP). New tools added by defining schema.",
        "LangChain tool registration. Squidpy RAG requires re-indexing for extensions.",
        "LangGraph framework. New tools registered as Python functions injected into REPL namespace.",
    ),
]


# ═══════════════════════════════════════════════════════════════════
# Write CSVs
# ═══════════════════════════════════════════════════════════════════

DISPLAY_NAMES = {
    "-Fibroblast": "CAF (Myofibroblast)",
    "B cell": "B cell",
    "Cancer": "Cancer",
    "Dendritic": "Dendritic",
    "Endothelial": "Endothelial",
    "Fibroblast": "Fibroblast",
    "Macrophage": "Macrophage",
    "Mast": "Mast cell",
    "T cell": "T cell",
    "myocyte": "Myocyte",
}


def write_effect_sizes():
    """S1: Merge sample-level tests and effect size summary into one table."""
    path = OUT_DIR / "Supplementary_Table_1_Effect_Sizes.csv"
    tests_path = CASE_DIR / "sample_level_tests.csv"
    effects_path = CASE_DIR / "effect_size_summary.csv"
    if not tests_path.exists() or not effects_path.exists():
        if path.exists() and path.stat().st_size > 0:
            print(f"  S1: {path.name} (existing committed table; source Case2_OSCC files not present)")
            return
        raise FileNotFoundError(
            "Supplementary Table S1 requires data/Case2_OSCC/sample_level_tests.csv "
            "and data/Case2_OSCC/effect_size_summary.csv. Run "
            "scripts/compute_case_study_stats.py first, or use the committed S1 table."
        )

    tests = pd.read_csv(tests_path)
    effects = pd.read_csv(effects_path)
    merged = tests.merge(effects, on="cell_type", how="left", suffixes=("", "_eff"))

    rows = []
    for _, r in merged.iterrows():
        ct = r["cell_type"]
        q = r["wilcoxon_p_fdr"]
        med = r["median_log2fc"]
        ci_lo, ci_hi = r["median_ci_lower"], r["median_ci_upper"]
        mean_r_all = r["mean_abs_r_all"]
        mean_r_sig = r["mean_abs_r_sig"]
        rows.append({
            "Cell Type": DISPLAY_NAMES.get(ct, ct),
            "Direction": r["direction"],
            "n (patients)": int(r["n_samples"]),
            "Wilcoxon signed-rank q (BH-FDR)": round(q, 4) if pd.notna(q) else "",
            "Median log2FC": round(med, 2) if pd.notna(med) else "",
            "Median log2FC 95% CI lower": round(ci_lo, 2) if pd.notna(ci_lo) else "",
            "Median log2FC 95% CI upper": round(ci_hi, 2) if pd.notna(ci_hi) else "",
            "Sign test (positive / negative)": f"{int(r['sign_n_pos'])}/{int(r['sign_n_neg'])}",
            "Sign test p": round(r["sign_p"], 4) if pd.notna(r["sign_p"]) else "",
            "Mean |r| (all samples)": round(mean_r_all, 3) if pd.notna(mean_r_all) else "",
            "n (all)": int(r["n_all"]),
            "Mean |r| (significant samples)": round(mean_r_sig, 3) if pd.notna(mean_r_sig) else "",
            "n (significant)": int(r["n_significant"]),
            "Consensus (% samples enriched)": f"{r['consensus_pct']:.0f}%",
        })

    out = pd.DataFrame(rows)
    out["_q"] = pd.to_numeric(out["Wilcoxon signed-rank q (BH-FDR)"], errors="coerce")
    out = out.sort_values("_q").drop(columns="_q")

    out.to_csv(path, index=False)
    print(f"  S1: {path.name} ({len(out)} cell types)")


def write_methods():
    path = OUT_DIR / "Supplementary_Table_2_Integrated_Methods.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["Analytical Category", "Method Name", "Version", "Primary Citation"])
        for row in METHODS:
            w.writerow(row)
    print(f"  S2: {path.name} ({len(METHODS)} methods)")


def write_scenarios():
    path = OUT_DIR / "Supplementary_Table_3_Test_Scenarios.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["ID", "Category", "Scenario", "Input Data", "Command / Goal",
                     "Expected Outcome", "Result", "Notes"])
        for row in SCENARIOS:
            w.writerow(row)
    print(f"  S3: {path.name} ({len(SCENARIOS)} scenarios)")


def write_ai_comparison():
    path = OUT_DIR / "Supplementary_Table_4_AI_Agent_Comparison.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["Feature", "ChatSpatial (This Work)", "STAgent", "SpatialAgent"])
        for row in AI_COMPARISON:
            w.writerow(row)
    print(f"  S4: {path.name} ({len(AI_COMPARISON)} rows)")


if __name__ == "__main__":
    print("Building supplementary tables...")
    write_effect_sizes()
    write_methods()
    write_scenarios()
    write_ai_comparison()
    print("Done.")
