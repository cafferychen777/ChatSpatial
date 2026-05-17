#!/usr/bin/env python3
"""
Case-Study Reproducibility Bridge Experiment

Tests whether the ablation-measured concordance extends to the specific
case-study workflow: OSCC CARD deconvolution (Step 2 from §2.3.1).

Runs the case-study prompt across four LLMs under full-schema and no-schema
conditions on two patient samples (S1, S9), measuring pairwise Pearson
correlation of the resulting proportion matrices.

Trial matrix: 2 samples × 2 conditions × 4 models × 10 reps = 160 trials.

CARD is deterministic given fixed parameters and data, so any output
variation is attributable entirely to schema-condition-induced parameter
variation.
"""

import asyncio
import csv
import json
import os
import sys
import time
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import requests
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
CODE_ROOT = SCRIPT_DIR.parent.parent / "code"
DATA_DIR = SCRIPT_DIR.parent / "data" / "casestudy_reproducibility"
PROP_DIR = DATA_DIR / "proportions"
DATA_DIR.mkdir(parents=True, exist_ok=True)
PROP_DIR.mkdir(parents=True, exist_ok=True)

RAW_PATH = DATA_DIR / "casestudy_repro_raw.jsonl"
METRICS_CSV = DATA_DIR / "casestudy_repro_metrics.csv"
SUMMARY_PATH = DATA_DIR / "casestudy_repro_summary.txt"

SAMPLES = {
    "oscc_s1": CODE_ROOT / "data" / "geo_data" / "analyzed_projects"
    / "GSE208253" / "processed_h5ad" / "s1_processed.h5ad",
    "oscc_s9": CODE_ROOT / "data" / "geo_data" / "analyzed_projects"
    / "GSE208253" / "processed_h5ad" / "s9_processed.h5ad",
}
PURAM_REF_PATH = (
    CODE_ROOT / "data" / "reference_datasets" / "puram_hnscc"
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

# Load .env from workspace root if present
_env_path = Path(__file__).resolve().parents[2] / ".env"
if _env_path.exists():
    with open(_env_path) as _f:
        for _line in _f:
            _line = _line.strip()
            if _line and not _line.startswith("#") and "=" in _line:
                _k, _, _v = _line.partition("=")
                os.environ.setdefault(_k.strip(), _v.strip())

GEMINI_KEY = os.environ.get("GEMINI_API_KEY", "")
ANTHROPIC_KEY = os.environ.get("ANTHROPIC_API_KEY", "")

OPENAI_API_URL = os.environ.get("OPENAI_API_URL", "https://api.openai.com/v1/messages")
OPENAI_MODEL = "gpt-5.4"

# Case-study prompt matching §2.3.1 Step 2
CASESTUDY_PROMPT = (
    "Deconvolve cell types in this OSCC sample using CARD with "
    "the HNSCC single-cell reference from Puram et al. "
    "(reference_data_id='puram_ref', cell_type_key='cell_type')."
)

# ---------------------------------------------------------------------------
# Ensure chatspatial importable
# ---------------------------------------------------------------------------

if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

import anndata as ad
import pandas as pd

from chatspatial.models.data import DeconvolutionParameters
from chatspatial.tools.deconvolution import deconvolve_spatial_data

# Import from ablation scripts
from ablation_invocation import (
    SCHEMA_VARIANTS,
    RESPONSE_INSTRUCTION,
    resolve_tool_name,
    parse_json,
    call_gemini,
    call_openai,
)
from ablation_e2e import (
    AblationCtx,
    remap_params,
    PARAM_ALIASES,
)

# ---------------------------------------------------------------------------
# Schema conditions (subset of ablation)
# ---------------------------------------------------------------------------

CONDITIONS = {
    "full_schema": SCHEMA_VARIANTS["full_schema"],
    "no_schema": SCHEMA_VARIANTS["no_schema"],
}

# ---------------------------------------------------------------------------
# Model configuration (ablation models + case-study model)
# ---------------------------------------------------------------------------


def _call_anthropic_impl(system: str, prompt: str, model: str) -> str | None:
    """Anthropic API caller with configurable model."""
    for attempt in range(MAX_RETRIES):
        try:
            payload = {
                "model": model,
                "max_tokens": 512,
                "temperature": TEMPERATURE,
                "system": system + RESPONSE_INSTRUCTION,
                "messages": [{"role": "user", "content": prompt}],
            }
            r = requests.post(
                "https://api.anthropic.com/v1/messages",
                headers={
                    "x-api-key": ANTHROPIC_KEY,
                    "anthropic-version": "2023-06-01",
                    "content-type": "application/json",
                },
                json=payload,
                timeout=60,
            )
            r.raise_for_status()
            return r.json()["content"][0]["text"]
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF * (attempt + 1))
            else:
                print(f"  Anthropic ({model}) error after {MAX_RETRIES} retries: {e}")
                return None


def call_anthropic_haiku(system: str, prompt: str) -> str | None:
    return _call_anthropic_impl(system, prompt, "claude-haiku-4-5-20251001")


def call_anthropic_sonnet(system: str, prompt: str) -> str | None:
    return _call_anthropic_impl(system, prompt, "claude-sonnet-4-5")


# Build model registry
MODELS: dict[str, str] = {}
CALLERS: dict[str, callable] = {}

if GEMINI_KEY:
    MODELS["gemini-2.5-flash"] = "gemini"
    CALLERS["gemini"] = call_gemini

if ANTHROPIC_KEY:
    MODELS["claude-haiku-4.5"] = "anthropic_haiku"
    CALLERS["anthropic_haiku"] = call_anthropic_haiku
    MODELS["claude-sonnet-4.5"] = "anthropic_sonnet"
    CALLERS["anthropic_sonnet"] = call_anthropic_sonnet

MODELS[OPENAI_MODEL] = "openai"
CALLERS["openai"] = call_openai


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
# Tool execution (reuses ablation pattern)
# ---------------------------------------------------------------------------

async def execute_deconv(
    params_dict: dict,
    adata_original: ad.AnnData,
    ref_original: ad.AnnData,
    data_id: str = "oscc_sample",
) -> tuple[bool, str, dict]:
    """Execute CARD deconvolution on a fresh data copy.

    Returns: (success, error_message, artifacts)
    """
    try:
        params = DeconvolutionParameters(**params_dict)
    except Exception as e:
        return False, f"pydantic_error:{str(e)[:300]}", {}

    adata = adata_original.copy()
    datasets = {data_id: adata, "puram_ref": ref_original.copy()}
    ctx = AblationCtx(datasets)

    try:
        result = await deconvolve_spatial_data(data_id, ctx, params)
    except Exception as e:
        return False, f"execution_error:{type(e).__name__}:{str(e)[:300]}", {}

    adata = await ctx.get_adata(data_id)
    artifacts = {}
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
        artifacts["method"] = params.method
    else:
        return False, f"missing_output_key:{prop_key}", {}

    return True, "", artifacts


# ---------------------------------------------------------------------------
# Concordance metrics
# ---------------------------------------------------------------------------

def pairwise_pearson(prop_matrices: list[np.ndarray]) -> list[float]:
    """Pairwise Pearson r of flattened proportion matrices."""
    correlations = []
    for i, j in combinations(range(len(prop_matrices)), 2):
        flat_i = prop_matrices[i].flatten()
        flat_j = prop_matrices[j].flatten()
        if len(flat_i) < 2 or np.std(flat_i) == 0 or np.std(flat_j) == 0:
            correlations.append(1.0 if np.array_equal(flat_i, flat_j) else 0.0)
        else:
            r, _ = pearsonr(flat_i, flat_j)
            correlations.append(r)
    return correlations


def bootstrap_ci(
    values: list[float],
    n_bootstrap: int = 10000,
    ci: float = 0.95,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Percentile bootstrap for mean. Returns (mean, ci_lower, ci_upper)."""
    if not values:
        return (np.nan, np.nan, np.nan)
    arr = np.array(values, dtype=float)
    if len(arr) == 1:
        return (arr[0], arr[0], arr[0])
    rng = np.random.default_rng(seed)
    boot_means = [
        np.mean(rng.choice(arr, size=len(arr), replace=True))
        for _ in range(n_bootstrap)
    ]
    alpha = (1 - ci) / 2
    return (
        float(np.mean(arr)),
        float(np.percentile(boot_means, alpha * 100)),
        float(np.percentile(boot_means, (1 - alpha) * 100)),
    )


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

async def run_experiment():
    print("=" * 70)
    print("Case-Study Reproducibility Bridge Experiment")
    print("=" * 70)
    print(f"Samples: {list(SAMPLES.keys())}")
    print(f"Conditions: {list(CONDITIONS.keys())}")
    print(f"Models: {list(MODELS.keys())}")
    print(f"Replicates: {N_TRIALS}")
    total = len(SAMPLES) * len(CONDITIONS) * len(MODELS) * N_TRIALS
    print(f"Total trials: {total}")
    print()

    # Load reference once
    print("Loading Puram reference...", end=" ", flush=True)
    ref_original = ad.read_h5ad(PURAM_REF_PATH)
    print(f"{ref_original.shape}")

    # Preload all samples
    import scanpy as sc

    sample_data: dict[str, ad.AnnData] = {}
    for sample_id, sample_path in SAMPLES.items():
        print(f"Loading {sample_id}...", end=" ", flush=True)
        adata = ad.read_h5ad(sample_path)
        print(f"{adata.shape}", end="")

        # Ensure preprocessing (same as ablation_e2e.py)
        if "highly_variable" not in adata.var.columns:
            print(" → preprocessing (HVGs + neighbors)...", end="")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            adata.raw = adata
            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
        elif "connectivities" not in adata.obsp:
            print(" → computing neighbors...", end="")
            sc.pp.neighbors(adata)

        sample_data[sample_id] = adata
        print(" ✓")

    # Resume
    completed = load_completed(RAW_PATH)
    print(f"\nResuming: {len(completed)} trials already completed\n")

    # Experiment loop
    for sample_id, adata_original in sample_data.items():
        for condition, schema_text in CONDITIONS.items():
            for model_name, provider in MODELS.items():
                caller = CALLERS[provider]
                label = f"[{sample_id}|{condition}|{model_name}]"
                print(f"{label} ", end="", flush=True)

                for rep in range(N_TRIALS):
                    trial_key = f"{sample_id}|{condition}|{model_name}|{rep}"
                    if trial_key in completed:
                        print("s", end="", flush=True)
                        continue

                    # Step 1: LLM invocation
                    raw_text = caller(schema_text, CASESTUDY_PROMPT)
                    parsed = parse_json(raw_text) if raw_text else None
                    raw_tool = parsed.get("tool_name") if parsed else None
                    canonical_tool = resolve_tool_name(raw_tool)
                    raw_params = parsed.get("parameters", {}) if parsed else {}

                    # Step 2: Pydantic validation
                    tier1_valid = False
                    tier2_valid = False
                    exec_success = False
                    exec_error = ""
                    artifacts = {}
                    use_params = raw_params

                    if canonical_tool == "deconvolve_data" and raw_params is not None:
                        # Tier 1: direct validation
                        try:
                            DeconvolutionParameters(**raw_params)
                            tier1_valid = True
                        except Exception:
                            pass

                        # Tier 2: alias remapping
                        if not tier1_valid:
                            remapped = remap_params("deconvolve_data", raw_params)
                            try:
                                DeconvolutionParameters(**remapped)
                                tier2_valid = True
                                use_params = remapped
                            except Exception:
                                pass

                    # Step 3: Execute
                    if tier1_valid or tier2_valid:
                        exec_success, exec_error, artifacts = await execute_deconv(
                            use_params, adata_original, ref_original,
                            data_id=sample_id,
                        )

                    # Step 4: Save artifacts
                    artifact_path = ""
                    if exec_success and artifacts:
                        fname = f"{sample_id}_{condition}_{model_name}_{rep}"
                        artifact_path = str(PROP_DIR / f"{fname}.json")
                        with open(artifact_path, "w") as f:
                            json.dump(artifacts, f)

                    # Step 5: Record
                    record = {
                        "trial_key": trial_key,
                        "sample_id": sample_id,
                        "condition": condition,
                        "model": model_name,
                        "rep": rep,
                        "parse_success": parsed is not None,
                        "raw_tool_name": raw_tool,
                        "canonical_tool": canonical_tool,
                        "tier1_valid": tier1_valid,
                        "tier2_valid": tier2_valid,
                        "exec_success": exec_success,
                        "exec_error": exec_error[:300] if exec_error else "",
                        "artifact_path": artifact_path,
                        "method_selected": artifacts.get("method", ""),
                        "raw_params": json.dumps(use_params, default=str)[:500],
                    }
                    append_raw(RAW_PATH, record)

                    status = "." if exec_success else "x"
                    print(status, end="", flush=True)
                    time.sleep(DELAY)

                print()  # newline after reps

    print("\n" + "=" * 70)
    print("Experiment complete. Running analysis...")
    print("=" * 70 + "\n")


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def run_analysis():
    """Compute concordance metrics from completed trials."""
    if not RAW_PATH.exists():
        print("No raw data found. Run experiment first.")
        return

    # Load all records
    records = []
    with open(RAW_PATH) as f:
        for line in f:
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                pass

    # Group by (sample, condition)
    groups = defaultdict(list)
    for rec in records:
        key = (rec["sample_id"], rec["condition"])
        groups[key].append(rec)

    metrics_rows = []
    summary_lines = []
    summary_lines.append("Case-Study Reproducibility Bridge — Summary")
    summary_lines.append("=" * 60)
    summary_lines.append(f"Total trials: {len(records)}")
    summary_lines.append("")

    for (sample_id, condition), cell_records in sorted(groups.items()):
        n_total = len(cell_records)
        n_success = sum(1 for r in cell_records if r["exec_success"])
        n_parse = sum(1 for r in cell_records if r["parse_success"])
        n_valid = sum(1 for r in cell_records if r["tier1_valid"] or r["tier2_valid"])

        summary_lines.append(f"--- {sample_id} | {condition} ---")
        summary_lines.append(f"  Trials: {n_total}  Parse: {n_parse}  Valid: {n_valid}  Exec: {n_success}")

        # Method distribution
        methods = [r["method_selected"] for r in cell_records if r["method_selected"]]
        if methods:
            from collections import Counter
            method_counts = Counter(methods)
            summary_lines.append(f"  Methods: {dict(method_counts)}")

        # Load proportion matrices for successful trials
        prop_matrices = []
        for rec in cell_records:
            if rec["exec_success"] and rec["artifact_path"]:
                try:
                    with open(rec["artifact_path"]) as f:
                        art = json.load(f)
                    prop_matrices.append(np.array(art["proportions"]))
                except (FileNotFoundError, json.JSONDecodeError, KeyError):
                    pass

        # Compute concordance
        if len(prop_matrices) >= 2:
            pearson_values = pairwise_pearson(prop_matrices)
            mean_r, ci_lo, ci_hi = bootstrap_ci(pearson_values)
            summary_lines.append(
                f"  Pearson r: {mean_r:.4f} [95% CI: {ci_lo:.4f}–{ci_hi:.4f}]"
                f"  (n_pairs={len(pearson_values)})"
            )
        elif len(prop_matrices) == 1:
            mean_r, ci_lo, ci_hi = np.nan, np.nan, np.nan
            pearson_values = []
            summary_lines.append("  Pearson r: N/A (only 1 successful trial)")
        else:
            mean_r, ci_lo, ci_hi = np.nan, np.nan, np.nan
            pearson_values = []
            summary_lines.append("  Pearson r: N/A (no successful trials)")

        # Per-model breakdown
        by_model = defaultdict(list)
        for rec in cell_records:
            by_model[rec["model"]].append(rec)

        for model_name, model_recs in sorted(by_model.items()):
            m_total = len(model_recs)
            m_success = sum(1 for r in model_recs if r["exec_success"])

            model_props = []
            for rec in model_recs:
                if rec["exec_success"] and rec["artifact_path"]:
                    try:
                        with open(rec["artifact_path"]) as f:
                            art = json.load(f)
                        model_props.append(np.array(art["proportions"]))
                    except (FileNotFoundError, json.JSONDecodeError, KeyError):
                        pass

            if len(model_props) >= 2:
                model_pearson = pairwise_pearson(model_props)
                m_mean, m_lo, m_hi = bootstrap_ci(model_pearson)
                summary_lines.append(
                    f"    {model_name}: {m_success}/{m_total} success, "
                    f"r={m_mean:.4f} [{m_lo:.4f}–{m_hi:.4f}]"
                )
            else:
                summary_lines.append(
                    f"    {model_name}: {m_success}/{m_total} success"
                )

        # Error analysis for failures
        failures = [r for r in cell_records if not r["exec_success"]]
        if failures:
            error_types = defaultdict(int)
            for r in failures:
                if not r["parse_success"]:
                    error_types["parse_failure"] += 1
                elif r["canonical_tool"] != "deconvolve_data":
                    error_types[f"wrong_tool:{r['canonical_tool']}"] += 1
                elif not (r["tier1_valid"] or r["tier2_valid"]):
                    error_types["pydantic_validation_failure"] += 1
                else:
                    err_short = r["exec_error"].split(":")[0] if r["exec_error"] else "unknown"
                    error_types[f"exec:{err_short}"] += 1
            summary_lines.append(f"  Failures: {dict(error_types)}")

        # Record metrics
        metrics_rows.append({
            "sample": sample_id,
            "condition": condition,
            "n_total": n_total,
            "n_success": n_success,
            "success_rate": n_success / n_total if n_total else 0,
            "pearson_mean": mean_r,
            "pearson_ci_lo": ci_lo,
            "pearson_ci_hi": ci_hi,
            "n_pairs": len(pearson_values),
        })

        summary_lines.append("")

    # Write metrics CSV
    if metrics_rows:
        with open(METRICS_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=metrics_rows[0].keys())
            writer.writeheader()
            writer.writerows(metrics_rows)
        print(f"Metrics written to {METRICS_CSV}")

    # Write summary
    summary_text = "\n".join(summary_lines)
    with open(SUMMARY_PATH, "w") as f:
        f.write(summary_text)
    print(f"Summary written to {SUMMARY_PATH}")
    print()
    print(summary_text)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    asyncio.run(run_experiment())
    run_analysis()
