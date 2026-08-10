#!/usr/bin/env python3
"""
Schema-Enforcement Ablation Study — Part 2: End-to-End Output Experiment

Runs the same three ablation conditions but ACTUALLY EXECUTES the tool calls
on OSCC Sample 1, measuring whether schema enforcement improves output
concordance across repeated trials.

Three tasks:
  1. Spatial domain identification (Leiden) → ARI / NMI
  2. SVG detection (FlashS)                → Jaccard of top-100
  3. Deconvolution (FlashDeconv)           → Pearson r of proportions

Trial matrix: 3 tasks × 3 models × 3 conditions × 10 reps = 270 trials.
"""

import asyncio
import csv
import json
import os
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from ablation_invocation import (
    SCHEMA_VARIANTS,
    call_anthropic,
    call_gemini,
    call_openai,
    parse_json,
    resolve_tool_name,
)
from paths import find_chatspatial_code_dir, load_env_file

from chatspatial.models.data import (
    DeconvolutionParameters,
    SpatialDomainParameters,
    SpatialVariableGenesParameters,
)
from chatspatial.tools.deconvolution import deconvolve_spatial_data
from chatspatial.tools.spatial_domains import identify_spatial_domains
from chatspatial.tools.spatial_genes import identify_spatial_genes

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
CODE_ROOT = find_chatspatial_code_dir(required=True)
DATA_DIR = SCRIPT_DIR.parent / "data" / "ablation" / "e2e"
OUTPUT_DIR = DATA_DIR / "outputs"

RAW_PATH = DATA_DIR / "ablation_e2e_raw.jsonl"
CSV_PATH = DATA_DIR / "ablation_e2e_results.csv"

OSCC_S1_PATH = (
    CODE_ROOT
    / "data"
    / "geo_data"
    / "analyzed_projects"
    / "GSE208253"
    / "processed_h5ad"
    / "s1_processed.h5ad"
)
PURAM_REF_PATH = (
    CODE_ROOT
    / "data"
    / "reference_datasets"
    / "puram_hnscc"
    / "puram_hnscc_reference.h5ad"
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_TRIALS = 10
TEMPERATURE = 1.0
DELAY = 1.0
MAX_RETRIES = 3
RETRY_BACKOFF = 2.0

load_env_file()

GEMINI_KEY = os.environ.get("GEMINI_API_KEY", "")
ANTHROPIC_KEY = os.environ.get("ANTHROPIC_API_KEY", "")
OPENAI_KEY = os.environ.get("OPENAI_API_KEY", "")
OPENAI_API_URL = os.environ.get("OPENAI_API_URL", "")

# OpenAI GPT-5.4 via Anthropic Messages-compatible API gateway
OPENAI_MODEL = "gpt-5.4"

MODELS: dict[str, str] = {}
if GEMINI_KEY:
    MODELS["gemini-2.5-flash"] = "gemini"
if ANTHROPIC_KEY:
    MODELS["claude-haiku-4-5-20251001"] = "anthropic"
if OPENAI_KEY and OPENAI_API_URL:
    MODELS[OPENAI_MODEL] = "openai"

# ---------------------------------------------------------------------------
# Minimal ToolContext for direct execution (bypasses MCP server)
# ---------------------------------------------------------------------------


class AblationCtx:
    """Minimal ToolContext matching the interface in spatial_mcp_adapter.py."""

    def __init__(self, datasets: dict[str, ad.AnnData]):
        self._datasets = datasets
        import logging

        self._logger = logging.getLogger("ablation")

    async def get_adata(self, data_id: str):
        if data_id not in self._datasets:
            raise ValueError(
                f"Dataset '{data_id}' not found. Available: {list(self._datasets.keys())}"
            )
        return self._datasets[data_id]

    async def get_dataset_info(self, data_id: str) -> dict:
        adata = self._datasets[data_id]
        return {"adata": adata, "name": data_id, "type": "visium"}

    async def set_adata(self, data_id: str, adata):
        self._datasets[data_id] = adata

    async def add_dataset(self, adata, prefix="data", name=None, metadata=None) -> str:
        new_id = f"{prefix}_{len(self._datasets)}"
        self._datasets[new_id] = adata
        return new_id

    async def info(self, msg: str):
        pass  # Suppress MCP info messages during experiment

    async def warning(self, msg: str):
        pass

    async def error(self, msg: str):
        pass

    def debug(self, msg: str):
        pass

    def log_config(self, title: str, config: dict):
        pass


# ---------------------------------------------------------------------------
# Schema Variants (identical to Part 1)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# E2E Task Definitions
# ---------------------------------------------------------------------------

E2E_TASKS = [
    {
        "id": "spatial_domains",
        "prompt": "Identify spatial domains in this Visium OSCC dataset. Use graph-based clustering.",
        "expected_tool": "identify_spatial_domains",
        "expected_method": "leiden",
    },
    {
        "id": "svg_detection",
        "prompt": "Find spatially variable genes in this Visium dataset.",
        "expected_tool": "find_spatial_genes",
        "expected_method": "flashs",
    },
    {
        "id": "deconvolution",
        "prompt": "Deconvolve cell types using the loaded reference (reference_data_id='puram_ref', cell_type_key='cell_type').",
        "expected_tool": "deconvolve_data",
        "expected_method": "flashdeconv",
    },
]

# Parameter aliases for Condition C (map free-form names to Pydantic fields)
PARAM_ALIASES = {
    "identify_spatial_domains": {
        "method": [
            "method",
            "clustering_method",
            "algorithm",
            "approach",
            "clustering_algorithm",
        ],
        "n_domains": [
            "n_domains",
            "num_domains",
            "n_clusters",
            "num_clusters",
            "k",
            "n_groups",
        ],
        "resolution": ["resolution", "res", "clustering_resolution"],
    },
    "find_spatial_genes": {
        "method": ["method", "svg_method", "algorithm", "test_method"],
        "n_top_genes": ["n_top_genes", "top_n", "num_genes", "n_genes", "top_genes"],
    },
    "deconvolve_data": {
        "method": ["method", "deconv_method", "algorithm"],
        "reference_data_id": [
            "reference_data_id",
            "ref_id",
            "reference_id",
            "ref_data_id",
            "reference",
        ],
        "cell_type_key": [
            "cell_type_key",
            "celltype_key",
            "cell_type_column",
            "celltype_column",
            "cell_type",
            "celltype",
            "annotation_key",
        ],
    },
}


def remap_params(tool_name: str, raw_params: dict) -> dict:
    """Best-effort mapping from free-form parameter names to canonical names."""
    aliases = PARAM_ALIASES.get(tool_name, {})
    if not aliases:
        return raw_params

    reverse_map = {}
    for canonical, alias_list in aliases.items():
        for alias in alias_list:
            reverse_map[alias.lower()] = canonical

    remapped = {}
    for k, v in raw_params.items():
        canonical = reverse_map.get(k.lower(), k)
        remapped[canonical] = v
    return remapped


CALLERS = {
    "gemini": call_gemini,
    "anthropic": call_anthropic,
    "openai": call_openai,
}

# ---------------------------------------------------------------------------
# Checkpoint / Resume
# ---------------------------------------------------------------------------


def load_completed(path: Path) -> set[str]:
    completed = set()
    if path.exists():
        with open(path) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    completed.add(rec["trial_key"])
                except (json.JSONDecodeError, KeyError):
                    pass
    return completed


def append_raw(path: Path, record: dict):
    with open(path, "a") as f:
        f.write(json.dumps(record, ensure_ascii=False, default=str) + "\n")


# ---------------------------------------------------------------------------
# Tool Execution
# ---------------------------------------------------------------------------

TOOL_PARAM_CLS = {
    "identify_spatial_domains": SpatialDomainParameters,
    "find_spatial_genes": SpatialVariableGenesParameters,
    "deconvolve_data": DeconvolutionParameters,
}

TOOL_FUNCS = {
    "identify_spatial_domains": identify_spatial_domains,
    "find_spatial_genes": identify_spatial_genes,
    "deconvolve_data": deconvolve_spatial_data,
}


async def execute_tool(
    tool_name: str,
    params_dict: dict,
    adata_original: ad.AnnData,
    ref_original: ad.AnnData | None,
) -> tuple[bool, str, dict]:
    """
    Execute a ChatSpatial tool with given parameters on a fresh copy of data.

    Returns: (success, error_message, output_artifacts)
    """
    param_cls = TOOL_PARAM_CLS.get(tool_name)
    tool_func = TOOL_FUNCS.get(tool_name)
    if param_cls is None or tool_func is None:
        return False, f"no_executor_for:{tool_name}", {}

    # Validate and instantiate parameters
    try:
        params = param_cls(**params_dict)
    except Exception as e:
        return False, f"pydantic_error:{str(e)[:300]}", {}

    # Fresh data copy for isolation
    adata = adata_original.copy()
    datasets = {"oscc_s1": adata}
    if ref_original is not None:
        datasets["puram_ref"] = ref_original.copy()

    ctx = AblationCtx(datasets)

    try:
        result = await tool_func("oscc_s1", ctx, params)
    except Exception as e:
        return False, f"execution_error:{type(e).__name__}:{str(e)[:300]}", {}

    # Collect output artifacts depending on tool type
    artifacts = {}
    # Re-fetch adata from ctx (tools may have modified it)
    adata = await ctx.get_adata("oscc_s1")

    if tool_name == "identify_spatial_domains":
        domain_key = result.domain_key
        if domain_key in adata.obs:
            labels = adata.obs[domain_key].values
            artifacts["domain_labels"] = labels.tolist()
            artifacts["n_domains"] = result.n_domains
            artifacts["domain_key"] = domain_key
        else:
            return False, f"missing_output_key:{domain_key}", {}

    elif tool_name == "find_spatial_genes":
        artifacts["spatial_genes"] = result.spatial_genes[:200]  # top 200
        artifacts["n_significant"] = result.n_significant_genes
        artifacts["n_analyzed"] = result.n_genes_analyzed

    elif tool_name == "deconvolve_data":
        prop_key = result.proportions_key
        if prop_key in adata.obsm:
            props = adata.obsm[prop_key]
            if isinstance(props, pd.DataFrame):
                artifacts["proportions"] = props.values.tolist()
                artifacts["cell_types"] = props.columns.tolist()
            else:
                artifacts["proportions"] = np.asarray(props).tolist()
                artifacts["cell_types"] = result.cell_types
            artifacts["n_cell_types"] = result.n_cell_types
        else:
            return False, f"missing_output_key:{prop_key}", {}

    return True, "", artifacts


# ---------------------------------------------------------------------------
# Main Experiment
# ---------------------------------------------------------------------------


async def run_experiment():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for subdirectory in ("spatial_domains", "svg", "deconvolution"):
        (OUTPUT_DIR / subdirectory).mkdir(parents=True, exist_ok=True)

    print("Schema-Enforcement Ablation — Part 2: End-to-End Output")
    print(f"Models: {list(MODELS.keys())}")
    print(f"Tasks: {[t['id'] for t in E2E_TASKS]}")
    total = len(MODELS) * len(SCHEMA_VARIANTS) * len(E2E_TASKS) * N_TRIALS
    print(f"Total trials: {total}")
    print()

    # Load data once
    print("Loading OSCC Sample 1...", end=" ", flush=True)
    adata_original = ad.read_h5ad(OSCC_S1_PATH)
    print(f"{adata_original.shape}")

    print("Loading Puram reference...", end=" ", flush=True)
    ref_original = ad.read_h5ad(PURAM_REF_PATH)
    print(f"{ref_original.shape}")

    # Preprocessing check: ensure neighbors and HVGs exist for Leiden
    import scanpy as sc

    adata_check = adata_original.copy()
    needs_preprocess = "highly_variable" not in adata_check.var.columns
    if needs_preprocess:
        print("Preprocessing OSCC S1 (HVGs + neighbors)...")
        sc.pp.normalize_total(adata_check, target_sum=1e4)
        sc.pp.log1p(adata_check)
        adata_check.raw = adata_check
        sc.pp.highly_variable_genes(adata_check, n_top_genes=2000)
        sc.pp.pca(adata_check)
        sc.pp.neighbors(adata_check)
        adata_original = adata_check
        print("  Done. HVGs and neighbors computed.")
    elif "connectivities" not in adata_original.obsp:
        print("Computing neighbors graph...")
        sc.pp.neighbors(adata_original)

    completed = load_completed(RAW_PATH)
    print(f"Resuming: {len(completed)} trials already completed\n")

    for condition, schema_text in SCHEMA_VARIANTS.items():
        for model_name, provider in MODELS.items():
            caller = CALLERS[provider]

            for task in E2E_TASKS:
                task_id = task["id"]
                prompt = task["prompt"]
                print(f"[{condition}|{model_name}|{task_id}] ", end="", flush=True)

                for rep in range(N_TRIALS):
                    trial_key = f"{condition}|{model_name}|{task_id}|{rep}"
                    if trial_key in completed:
                        print("s", end="", flush=True)
                        continue

                    # Step 1: Get LLM response
                    raw_text = caller(schema_text, prompt)
                    parsed = parse_json(raw_text) if raw_text else None
                    raw_tool = parsed.get("tool_name") if parsed else None
                    canonical_tool = resolve_tool_name(raw_tool)
                    raw_params = parsed.get("parameters", {}) if parsed else {}

                    # Step 2: Try Pydantic validation (Tier 1)
                    tier1_valid = False
                    tier2_valid = False
                    exec_success = False
                    exec_error = ""
                    artifacts = {}

                    if canonical_tool and raw_params is not None:
                        param_cls = TOOL_PARAM_CLS.get(canonical_tool)
                        if param_cls:
                            try:
                                param_cls(**raw_params)
                                tier1_valid = True
                            except Exception:
                                pass

                        # Tier 2: Alias remapping
                        if not tier1_valid:
                            remapped = remap_params(canonical_tool, raw_params)
                            if param_cls:
                                try:
                                    param_cls(**remapped)
                                    tier2_valid = True
                                    raw_params = remapped  # Use remapped for execution
                                except Exception:
                                    pass

                    # Step 3: Execute if validation passed
                    use_params = raw_params
                    if tier1_valid or tier2_valid:
                        exec_success, exec_error, artifacts = await execute_tool(
                            canonical_tool,
                            use_params,
                            adata_original,
                            ref_original,
                        )

                    # Step 4: Save output artifacts
                    artifact_path = ""
                    if exec_success and artifacts:
                        fname = f"{condition}_{model_name}_{task_id}_{rep}"
                        if task_id == "spatial_domains":
                            p = OUTPUT_DIR / "spatial_domains" / f"{fname}.json"
                            with open(p, "w") as f:
                                json.dump(artifacts, f)
                            artifact_path = str(p)
                        elif task_id == "svg_detection":
                            p = OUTPUT_DIR / "svg" / f"{fname}.json"
                            with open(p, "w") as f:
                                json.dump(artifacts, f)
                            artifact_path = str(p)
                        elif task_id == "deconvolution":
                            p = OUTPUT_DIR / "deconvolution" / f"{fname}.json"
                            with open(p, "w") as f:
                                json.dump(artifacts, f)
                            artifact_path = str(p)

                    record = {
                        "trial_key": trial_key,
                        "condition": condition,
                        "model": model_name,
                        "task_id": task_id,
                        "rep": rep,
                        "parse_success": parsed is not None,
                        "raw_tool_name": raw_tool,
                        "canonical_tool": canonical_tool,
                        "tier1_valid": tier1_valid,
                        "tier2_valid": tier2_valid,
                        "exec_success": exec_success,
                        "exec_error": exec_error[:300],
                        "artifact_path": artifact_path,
                        "raw_params": json.dumps(raw_params)[:500],
                    }
                    append_raw(RAW_PATH, record)
                    completed.add(trial_key)

                    status = (
                        "O"
                        if exec_success
                        else (
                            "v"
                            if (tier1_valid or tier2_valid)
                            else ("p" if parsed else "x")
                        )
                    )
                    print(status, end="", flush=True)
                    time.sleep(DELAY)

                print()

    # ---- Summary ----
    print("\n\nGenerating summary...")
    all_records = []
    with open(RAW_PATH) as f:
        for line in f:
            try:
                all_records.append(json.loads(line))
            except json.JSONDecodeError:
                pass

    rows = []
    from collections import defaultdict

    groups = defaultdict(list)
    for r in all_records:
        groups[(r["condition"], r["model"], r["task_id"])].append(r)

    for (cond, model, task_id), recs in sorted(groups.items()):
        n = len(recs)
        rows.append(
            {
                "condition": cond,
                "model": model,
                "task_id": task_id,
                "n_trials": n,
                "parse_rate": round(sum(r["parse_success"] for r in recs) / n, 3),
                "tier1_rate": round(sum(r["tier1_valid"] for r in recs) / n, 3),
                "tier2_rate": round(sum(r["tier2_valid"] for r in recs) / n, 3),
                "any_valid_rate": round(
                    sum(r["tier1_valid"] or r["tier2_valid"] for r in recs) / n, 3
                ),
                "exec_success_rate": round(sum(r["exec_success"] for r in recs) / n, 3),
            }
        )

    if rows:
        with open(CSV_PATH, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys(), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"Results CSV: {CSV_PATH}")

    # Print quick summary
    for cond in SCHEMA_VARIANTS:
        cond_recs = [r for r in all_records if r["condition"] == cond]
        if not cond_recs:
            continue
        n = len(cond_recs)
        print(
            f"\n{cond.upper()}: "
            f"parse={sum(r['parse_success'] for r in cond_recs)/n:.0%} "
            f"valid={sum(r['tier1_valid'] or r['tier2_valid'] for r in cond_recs)/n:.0%} "
            f"exec={sum(r['exec_success'] for r in cond_recs)/n:.0%}"
        )

    print("\nDone! Run ablation_analysis.py to compute concordance metrics.")


if __name__ == "__main__":
    asyncio.run(run_experiment())
