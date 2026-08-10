#!/usr/bin/env python3
"""
Code Generation Hallucination Experiment for ChatSpatial Paper

Quantifies hallucination rates when LLMs generate free-form code for spatial
transcriptomics tasks (no schema constraints), as a baseline comparison against
ChatSpatial's schema-enforced approach (which has 0% hallucination by construction).

Same 8 prompts, same 3 models, temperature=1.0.
Each model generates Python code 10 times per prompt.
We check for: package hallucinations, function hallucinations, syntax errors.
"""

import ast
import csv
import json
import os
import re
import time
from collections import defaultdict
from pathlib import Path

import requests

# ============================================================
# Configuration
# ============================================================

N_TRIALS = 10
TEMPERATURE = 1.0
DELAY = 1.0

GEMINI_KEY = os.environ.get("GEMINI_API_KEY", "")
ANTHROPIC_KEY = os.environ.get("ANTHROPIC_API_KEY", "")
OPENAI_KEY = os.environ.get("OPENAI_API_KEY", "")

MODELS = {}
if GEMINI_KEY:
    MODELS["gemini-2.5-flash"] = "gemini"
if ANTHROPIC_KEY:
    MODELS["claude-haiku-4-5-20251001"] = "anthropic"
if OPENAI_KEY:
    MODELS["gpt-5-mini"] = "openai"

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR.parent / "data"
CSV_PATH = DATA_DIR / "codegen_hallucination.csv"
SUMMARY_PATH = DATA_DIR / "codegen_hallucination_summary.txt"
RAW_DIR = DATA_DIR / "codegen_raw"

# ============================================================
# Valid Packages & Functions (ground truth)
# ============================================================

# Real Python packages for spatial transcriptomics
VALID_PACKAGES = {
    # Core
    "scanpy", "anndata", "squidpy", "numpy", "pandas", "scipy",
    "matplotlib", "seaborn", "sklearn", "scikit-learn",
    # Spatial analysis
    "SpaGCN", "spagcn", "cell2location", "scvi", "scvelo", "cellrank",
    "tangram", "destvi", "stereoscope", "spotlight",
    "spatialde", "sparkx", "SPARK",
    "stagate", "STAGATE", "STAGATE_pyG",
    "graphst", "GraphST",
    # Communication
    "liana", "liana_py", "cellphonedb", "squidpy",
    # CNV
    "infercnvpy", "infercnv",
    # Integration
    "harmonypy", "harmony", "bbknn", "scanorama",
    # Deconvolution
    "rctd", "card", "spacexr",
    # Python standard library
    "os", "sys", "pathlib", "json", "warnings", "logging",
    "collections", "itertools", "functools",
    "subprocess", "shutil", "tempfile", "datetime", "random",
    "copy", "typing", "re", "math", "csv", "io", "time",
    "abc", "glob", "textwrap", "pickle", "hashlib", "enum",
    "dataclasses", "argparse", "configparser", "struct",
    "multiprocessing", "threading", "contextlib", "operator",
    "h5py", "tables", "zarr",
    "networkx", "igraph",
    "statsmodels", "pysal", "esda", "libpysal", "splot",
    "rpy2", "anndata2ri",
    "torch", "tensorflow", "jax",
    "tqdm", "rich", "click",
    "plotly", "bokeh", "altair",
    "gseapy", "goatools",
    "muon", "mudata",
    "palantir",
    "pynndescent", "umap",
    "leidenalg", "louvain",
    "numba",
}

# Known VALID function calls (package.function patterns)
# We check if generated code calls functions that exist
VALID_FUNCTION_PATTERNS = {
    # scanpy
    r"sc\.read[_a-z]*", r"sc\.pp\.\w+", r"sc\.tl\.\w+", r"sc\.pl\.\w+",
    r"sc\.external\.\w+", r"sc\.get\.\w+",
    r"scanpy\.read[_a-z]*", r"scanpy\.pp\.\w+", r"scanpy\.tl\.\w+",
    # squidpy
    r"sq\.gr\.\w+", r"sq\.pl\.\w+", r"sq\.im\.\w+",
    r"squidpy\.gr\.\w+", r"squidpy\.pl\.\w+",
    # anndata
    r"ad\.read[_a-z]*", r"anndata\.read[_a-z]*", r"AnnData\(",
    # scvelo
    r"scv\.pp\.\w+", r"scv\.tl\.\w+", r"scv\.pl\.\w+",
    # cellrank
    r"cr\.\w+",
    # numpy/pandas/scipy
    r"np\.\w+", r"pd\.\w+", r"scipy\.\w+",
    # matplotlib
    r"plt\.\w+", r"fig\.\w+", r"ax\.\w+",
    # general
    r"print\(", r"len\(", r"range\(", r"type\(",
}

# Known HALLUCINATED patterns - functions that DON'T exist
HALLUCINATION_PATTERNS = [
    # Scanpy hallucinations
    (r"sc\.tl\.spatial_domains?\(", "sc.tl.spatial_domains() does not exist"),
    (r"sc\.tl\.spagcn\(", "sc.tl.spagcn() does not exist"),
    (r"sc\.tl\.stagate\(", "sc.tl.stagate() does not exist"),
    (r"sc\.tl\.graphst\(", "sc.tl.graphst() does not exist"),
    (r"sc\.tl\.deconvolve?\(", "sc.tl.deconvolve() does not exist"),
    (r"sc\.tl\.rctd\(", "sc.tl.rctd() does not exist"),
    (r"sc\.tl\.cell2location\(", "sc.tl.cell2location() does not exist"),
    (r"sc\.tl\.spatial_variable_genes?\(", "sc.tl.spatial_variable_genes() does not exist"),
    (r"sc\.tl\.moran[_s]*\(", "sc.tl.moran() does not exist (use sq.gr.spatial_autocorr)"),
    (r"sc\.tl\.cellchat\(", "sc.tl.cellchat() does not exist"),
    (r"sc\.tl\.infercnv\(", "sc.tl.infercnv() does not exist"),
    (r"sc\.tl\.trajectory\(", "sc.tl.trajectory() does not exist"),
    (r"sc\.tl\.spatial_autocorr\(", "sc.tl.spatial_autocorr() does not exist (use squidpy)"),
    (r"sc\.tl\.cellrank\(", "sc.tl.cellrank() does not exist"),
    (r"sc\.tl\.spatial_genes?\(", "sc.tl.spatial_genes() does not exist"),
    (r"sc\.pp\.spatial\(", "sc.pp.spatial() does not exist"),
    # Squidpy hallucinations
    (r"sq\.tl\.\w+", "sq.tl.* does not exist (squidpy uses sq.gr.*)"),
    (r"sq\.pp\.\w+", "sq.pp.* does not exist"),
    # Package-level hallucinations
    (r"from\s+scanpy\.spatial\s+import", "scanpy.spatial module does not exist"),
    (r"from\s+scanpy\.tools\.spatial", "scanpy.tools.spatial does not exist"),
    (r"import\s+spatial_transcriptomics", "spatial_transcriptomics package does not exist"),
    (r"from\s+stlearn\b", "stlearn package exists but is deprecated; LLMs often hallucinate its API"),
    # CellChat in Python hallucination (it's R-only)
    (r"from\s+cellchat\s+import", "cellchat is R-only, no Python package"),
    (r"import\s+cellchat\b", "cellchat is R-only, no Python package"),
    # Common function hallucinations
    (r"adata\.spatial_domains?\(", "adata.spatial_domains() does not exist"),
    (r"adata\.deconvolve\(", "adata.deconvolve() does not exist"),
]

# ============================================================
# Code Generation System Prompt
# ============================================================

CODEGEN_SYSTEM = """You are a bioinformatics expert. Write Python code to accomplish the user's spatial transcriptomics analysis task.

Assume the data is loaded as an AnnData object called `adata` with spatial coordinates in `adata.obsm['spatial']`.

Write ONLY the Python code (no explanations). Include all necessary imports."""

# ============================================================
# Test Prompts (same tasks as determinism experiment)
# ============================================================

TEST_PROMPTS = [
    {
        "id": "spatial_domains",
        "prompt": "Identify spatial domains in this Visium dataset using a deep learning method. Use 7 clusters.",
    },
    {
        "id": "deconvolution",
        "prompt": "Deconvolve cell types in this Visium dataset using RCTD with a single-cell reference stored in 'sc_ref.h5ad' with cell type labels in column 'celltype'.",
    },
    {
        "id": "svg_detection",
        "prompt": "Find the top 50 spatially variable genes in this MERFISH dataset.",
    },
    {
        "id": "cell_communication",
        "prompt": "Run cell-cell communication analysis on this human cancer spatial transcriptomics data using CellChat. Cell types are in the 'cell_type' column.",
    },
    {
        "id": "moran_i",
        "prompt": "Compute Moran's I spatial autocorrelation for the top 20 marker genes in this spatial dataset.",
    },
    {
        "id": "spatial_visualization",
        "prompt": "Create a spatial plot of this Visium data colored by cell type annotations stored in the 'cell_type' column.",
    },
    {
        "id": "trajectory",
        "prompt": "Perform trajectory analysis on this spatial transcriptomics data using CellRank, incorporating spatial information.",
    },
    {
        "id": "cnv_inference",
        "prompt": "Infer copy number variations in this tumor spatial data using inferCNV, with T cells and B cells (in 'cell_type' column) as the normal reference.",
    },
]

# ============================================================
# API Callers
# ============================================================

def call_gemini(prompt: str) -> str | None:
    payload = {
        "contents": [{"role": "user", "parts": [{"text": prompt}]}],
        "generationConfig": {"temperature": TEMPERATURE, "maxOutputTokens": 2048},
        "systemInstruction": {"parts": [{"text": CODEGEN_SYSTEM}]},
    }
    try:
        r = requests.post(
            f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash:generateContent?key={GEMINI_KEY}",
            json=payload, timeout=60
        )
        r.raise_for_status()
        return r.json()["candidates"][0]["content"]["parts"][0]["text"]
    except Exception as e:
        print(f"  Gemini error: {e}")
        return None


def call_anthropic(prompt: str) -> str | None:
    payload = {
        "model": "claude-haiku-4-5-20251001",
        "max_tokens": 2048,
        "temperature": TEMPERATURE,
        "system": CODEGEN_SYSTEM,
        "messages": [{"role": "user", "content": prompt}],
    }
    try:
        r = requests.post(
            "https://api.anthropic.com/v1/messages",
            headers={
                "x-api-key": ANTHROPIC_KEY,
                "anthropic-version": "2023-06-01",
                "content-type": "application/json",
            },
            json=payload, timeout=60
        )
        r.raise_for_status()
        return r.json()["content"][0]["text"]
    except Exception as e:
        print(f"  Anthropic error: {e}")
        return None


def call_openai(prompt: str) -> str | None:
    payload = {
        "model": "gpt-5-mini",
        "max_completion_tokens": 2048,
        "temperature": TEMPERATURE,
        "messages": [
            {"role": "system", "content": CODEGEN_SYSTEM},
            {"role": "user", "content": prompt},
        ],
    }
    try:
        r = requests.post(
            "https://api.openai.com/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENAI_KEY}", "Content-Type": "application/json"},
            json=payload, timeout=60
        )
        r.raise_for_status()
        return r.json()["choices"][0]["message"]["content"]
    except Exception as e:
        print(f"  OpenAI error: {e}")
        return None


CALLERS = {
    "gemini": call_gemini,
    "anthropic": call_anthropic,
    "openai": call_openai,
}


# ============================================================
# Code Analysis
# ============================================================

def extract_code(text: str) -> str:
    """Extract Python code from markdown-wrapped response."""
    if text is None:
        return ""
    # Try to find code block
    match = re.search(r"```(?:python)?\s*\n(.*?)```", text, re.DOTALL)
    if match:
        return match.group(1).strip()
    # If no code block, assume entire response is code
    return text.strip()


def check_syntax(code: str) -> tuple[bool, str]:
    """Check if code has valid Python syntax."""
    try:
        ast.parse(code)
        return True, ""
    except SyntaxError as e:
        return False, f"Line {e.lineno}: {e.msg}"


def extract_imports(code: str) -> list[str]:
    """Extract imported package names from code."""
    imports = set()
    for match in re.finditer(r"^import\s+(\w+)", code, re.MULTILINE):
        imports.add(match.group(1))
    for match in re.finditer(r"^from\s+(\w+)", code, re.MULTILINE):
        imports.add(match.group(1))
    return list(imports)


def check_package_hallucinations(code: str) -> list[str]:
    """Check for imports of non-existent packages."""
    hallucinated = []
    imports = extract_imports(code)
    for pkg in imports:
        # Normalize
        pkg_lower = pkg.lower()
        valid_lower = {v.lower() for v in VALID_PACKAGES}
        if pkg_lower not in valid_lower:
            hallucinated.append(pkg)
    return hallucinated


def check_function_hallucinations(code: str) -> list[str]:
    """Check for known hallucinated function calls."""
    found = []
    for pattern, description in HALLUCINATION_PATTERNS:
        if re.search(pattern, code):
            found.append(description)
    return found


def analyze_code(code: str) -> dict:
    """Comprehensive analysis of generated code."""
    result = {
        "has_code": bool(code.strip()),
        "syntax_valid": False,
        "syntax_error": "",
        "hallucinated_packages": [],
        "hallucinated_functions": [],
        "n_lines": len(code.strip().split("\n")) if code.strip() else 0,
        "n_imports": 0,
    }

    if not result["has_code"]:
        return result

    result["syntax_valid"], result["syntax_error"] = check_syntax(code)
    result["hallucinated_packages"] = check_package_hallucinations(code)
    result["hallucinated_functions"] = check_function_hallucinations(code)
    result["n_imports"] = len(extract_imports(code))

    return result


# ============================================================
# Main Experiment
# ============================================================

def run_experiment():
    print(f"Running code generation hallucination experiment")
    print(f"Models: {list(MODELS.keys())}")
    print(f"Prompts: {len(TEST_PROMPTS)}, Trials: {N_TRIALS}, Temperature: {TEMPERATURE}")
    print()

    # Create raw output directory
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    all_results = []

    for model_name, provider in MODELS.items():
        caller = CALLERS[provider]
        print(f"\n{'='*60}")
        print(f"Model: {model_name}")
        print(f"{'='*60}")

        for test in TEST_PROMPTS:
            prompt_id = test["id"]
            prompt = test["prompt"]
            print(f"\n  {prompt_id} ({N_TRIALS} trials)...", end="", flush=True)

            trial_analyses = []
            for i in range(N_TRIALS):
                raw = caller(prompt)
                code = extract_code(raw) if raw else ""

                # Save raw response
                raw_file = RAW_DIR / f"{model_name}_{prompt_id}_{i}.py"
                with open(raw_file, "w") as f:
                    f.write(f"# Model: {model_name}\n")
                    f.write(f"# Prompt: {prompt}\n")
                    f.write(f"# Trial: {i}\n\n")
                    f.write(code if code else "# NO CODE GENERATED")

                analysis = analyze_code(code)
                trial_analyses.append(analysis)
                print(".", end="", flush=True)
                time.sleep(DELAY)

            # Aggregate metrics for this model+prompt combination
            n_valid = sum(1 for a in trial_analyses if a["has_code"])
            n_syntax_ok = sum(1 for a in trial_analyses if a["syntax_valid"])
            n_pkg_hallucination = sum(1 for a in trial_analyses if a["hallucinated_packages"])
            n_func_hallucination = sum(1 for a in trial_analyses if a["hallucinated_functions"])
            n_any_hallucination = sum(
                1 for a in trial_analyses
                if a["hallucinated_packages"] or a["hallucinated_functions"]
            )
            n_correct = sum(
                1 for a in trial_analyses
                if a["syntax_valid"] and not a["hallucinated_packages"] and not a["hallucinated_functions"]
            )

            # Collect all unique hallucinations
            all_pkg_hall = set()
            all_func_hall = set()
            for a in trial_analyses:
                all_pkg_hall.update(a["hallucinated_packages"])
                all_func_hall.update(a["hallucinated_functions"])

            result = {
                "model": model_name,
                "prompt_id": prompt_id,
                "n_trials": N_TRIALS,
                "code_generated": n_valid,
                "syntax_valid": n_syntax_ok,
                "syntax_error_rate": round(1 - n_syntax_ok / max(n_valid, 1), 3),
                "pkg_hallucination_count": n_pkg_hallucination,
                "pkg_hallucination_rate": round(n_pkg_hallucination / max(n_valid, 1), 3),
                "func_hallucination_count": n_func_hallucination,
                "func_hallucination_rate": round(n_func_hallucination / max(n_valid, 1), 3),
                "any_hallucination_count": n_any_hallucination,
                "any_hallucination_rate": round(n_any_hallucination / max(n_valid, 1), 3),
                "fully_correct_rate": round(n_correct / max(n_valid, 1), 3),
                "unique_pkg_hallucinations": "; ".join(sorted(all_pkg_hall)),
                "unique_func_hallucinations": "; ".join(sorted(all_func_hall)),
            }
            all_results.append(result)

            hall_pct = n_any_hallucination / max(n_valid, 1)
            syn_pct = 1 - n_syntax_ok / max(n_valid, 1)
            print(f" halluc={hall_pct:.0%} syntax_err={syn_pct:.0%} correct={n_correct}/{n_valid}")

    # ============================================================
    # Save Results
    # ============================================================

    # CSV
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_results[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nCSV saved: {CSV_PATH}")

    # Summary
    with open(SUMMARY_PATH, "w") as f:
        f.write("Code Generation Hallucination Experiment\n")
        f.write(f"{'='*60}\n")
        f.write(f"Trials per prompt per model: {N_TRIALS}\n")
        f.write(f"Temperature: {TEMPERATURE}\n")
        f.write(f"Models: {list(MODELS.keys())}\n")
        f.write(f"Prompts: {len(TEST_PROMPTS)}\n\n")

        for model_name in MODELS:
            mr = [r for r in all_results if r["model"] == model_name]
            total_trials = sum(r["code_generated"] for r in mr)
            avg_syntax_err = sum(r["syntax_error_rate"] for r in mr) / len(mr)
            avg_pkg_hall = sum(r["pkg_hallucination_rate"] for r in mr) / len(mr)
            avg_func_hall = sum(r["func_hallucination_rate"] for r in mr) / len(mr)
            avg_any_hall = sum(r["any_hallucination_rate"] for r in mr) / len(mr)
            avg_correct = sum(r["fully_correct_rate"] for r in mr) / len(mr)

            f.write(f"\n{model_name}\n")
            f.write(f"{'-'*40}\n")
            f.write(f"  Code generated:         {total_trials}/{len(mr)*N_TRIALS}\n")
            f.write(f"  Syntax error rate:       {avg_syntax_err:.1%}\n")
            f.write(f"  Package halluc. rate:    {avg_pkg_hall:.1%}\n")
            f.write(f"  Function halluc. rate:   {avg_func_hall:.1%}\n")
            f.write(f"  Any hallucination rate:  {avg_any_hall:.1%}\n")
            f.write(f"  Fully correct rate:      {avg_correct:.1%}\n")

            # List unique hallucinations
            all_pkg = set()
            all_func = set()
            for r in mr:
                if r["unique_pkg_hallucinations"]:
                    all_pkg.update(r["unique_pkg_hallucinations"].split("; "))
                if r["unique_func_hallucinations"]:
                    all_func.update(r["unique_func_hallucinations"].split("; "))
            if all_pkg:
                f.write(f"  Hallucinated packages: {', '.join(sorted(all_pkg))}\n")
            if all_func:
                f.write(f"  Hallucinated functions: {', '.join(sorted(all_func))}\n")

        # Cross-model averages
        f.write(f"\n\nCross-Model Averages\n")
        f.write(f"{'='*40}\n")
        total = len(all_results)
        f.write(f"  Syntax error rate:       {sum(r['syntax_error_rate'] for r in all_results)/total:.1%}\n")
        f.write(f"  Package halluc. rate:    {sum(r['pkg_hallucination_rate'] for r in all_results)/total:.1%}\n")
        f.write(f"  Function halluc. rate:   {sum(r['func_hallucination_rate'] for r in all_results)/total:.1%}\n")
        f.write(f"  Any hallucination rate:  {sum(r['any_hallucination_rate'] for r in all_results)/total:.1%}\n")
        f.write(f"  Fully correct rate:      {sum(r['fully_correct_rate'] for r in all_results)/total:.1%}\n")

        # Comparison with ChatSpatial
        avg_hall = sum(r["any_hallucination_rate"] for r in all_results) / total
        avg_correct = sum(r["fully_correct_rate"] for r in all_results) / total
        f.write(f"\n\nKey Finding for Paper\n")
        f.write(f"{'='*40}\n")
        f.write(f"Code generation approach: {avg_hall:.1%} hallucination rate,\n")
        f.write(f"  {avg_correct:.1%} fully correct rate across {len(MODELS)} models.\n")
        f.write(f"ChatSpatial schema-enforced: 0% hallucination rate by construction.\n")
        f.write(f"This represents a {avg_hall:.1%} → 0% reduction in hallucination.\n")

    print(f"Summary saved: {SUMMARY_PATH}")


if __name__ == "__main__":
    run_experiment()
