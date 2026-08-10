#!/usr/bin/env python3
"""
End-to-End Benchmark: ChatSpatial vs STAgent vs SpatialAgent

Runs all three systems on the same dataset (human lymph node Visium),
same tasks (spatial domains, spatially variable genes), same LLM
(Claude Sonnet 4), and measures execution success, method consistency,
and output concordance.

Trial matrix: 3 tasks x 3 systems x 5 reps = 45 runs.
"""

import asyncio
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from paths import find_benchmarks_dir, find_chatspatial_code_dir, load_env_file

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent
load_env_file()
CODE_ROOT = find_chatspatial_code_dir(required=True)
BENCH_ROOT = find_benchmarks_dir(required=True)

DATA_DIR = REPO_ROOT / "data" / "e2e_benchmark"
OUTPUT_DIR = DATA_DIR / "outputs"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RAW_PATH = DATA_DIR / "e2e_benchmark_raw.jsonl"

# Benchmark dataset
DATASET_PATH = CODE_ROOT / "data" / "processed_datasets" / "visium" / "human_lymph_node.h5ad"
DATASET_PATH_CCC = CODE_ROOT / "data" / "processed_datasets" / "visium" / "human_lymph_node_cellphonedb.h5ad"

# System venvs
STAGENT_VENV = BENCH_ROOT / "STAgent" / ".venv" / "bin" / "python"
SPATIALAGENT_VENV = BENCH_ROOT / "SpatialAgent" / ".venv" / "bin" / "python"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_REPS = 5
MODEL = "claude-sonnet-4-20250514"
TEMPERATURE = 1.0
TIMEOUT_SECONDS = 900  # 15 minutes per run

# API keys from environment
ANTHROPIC_KEY = os.environ.get("ANTHROPIC_API_KEY", "")

# ---------------------------------------------------------------------------
# Tasks
# ---------------------------------------------------------------------------

TASKS = [
    {
        "id": "spatial_domains",
        "description": "Spatial domain identification via Leiden clustering",
        "prompt_template": (
            "Load the Visium dataset at '{data_path}'. "
            "Preprocess the data (normalize, log-transform, select highly variable genes, "
            "compute PCA). Compute spatial neighbors using squidpy. "
            "Identify spatial domains using Leiden clustering. "
            "Save the resulting AnnData object (with cluster labels in obs['leiden']) "
            "to '{output_path}'."
        ),
        "output_file": "result.h5ad",
        "chatspatial_tool": "identify_spatial_domains",
        "chatspatial_params": {
            "method": "leiden",
            "n_domains": 7,
            "resolution": 0.5,
        },
    },
    {
        "id": "svg_detection",
        "description": "Spatially variable gene detection",
        "prompt_template": (
            "Load the Visium dataset at '{data_path}'. "
            "Preprocess the data (normalize, log-transform, select highly variable genes). "
            "Compute spatial neighbors. "
            "Find the top 100 spatially variable genes. "
            "Save the ranked gene list as a CSV file (columns: gene, score, pvalue) "
            "to '{output_path}'."
        ),
        "output_file": "svg_results.csv",
        "chatspatial_tool": "find_spatial_genes",
        "chatspatial_params": {
            "method": "flashs",
            "n_top_genes": 100,
            "filter_mt_genes": True,
        },
    },
    {
        "id": "cell_communication",
        "description": "Cell-cell communication via ligand-receptor analysis",
        "prompt_template": (
            "Load the Visium dataset at '{data_path}'. "
            "The data is human lymph node tissue with Leiden cluster labels in obs['leiden']. "
            "Analyze cell-cell communication using ligand-receptor analysis. "
            "Use species='human' and cell_type_key='leiden'. "
            "Identify the top 50 significant ligand-receptor interactions. "
            "Save the results as a CSV file with columns: ligand, receptor, score "
            "(or equivalent significance measure) to '{output_path}'."
        ),
        "output_file": "ccc_results.csv",
        "chatspatial_tool": "analyze_cell_communication",
        "chatspatial_params": {
            "method": "fastccc",
            "species": "human",
            "cell_type_key": "leiden",
        },
    },
]

# ---------------------------------------------------------------------------
# Ensure chatspatial importable
# ---------------------------------------------------------------------------

if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))


# ---------------------------------------------------------------------------
# Checkpoint / Resume
# ---------------------------------------------------------------------------

def load_completed():
    """Load set of completed (system, task, rep) tuples."""
    completed = set()
    if RAW_PATH.exists():
        with open(RAW_PATH) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    completed.add((rec["system"], rec["task_id"], rec["rep"]))
                except (json.JSONDecodeError, KeyError):
                    pass
    return completed


def append_raw(record: dict):
    """Append one trial record to the raw JSONL."""
    with open(RAW_PATH, "a") as f:
        f.write(json.dumps(record, default=str) + "\n")


# ---------------------------------------------------------------------------
# ChatSpatial Execution (direct, same process)
# ---------------------------------------------------------------------------

async def run_chatspatial(task: dict, output_dir: Path, data_path: Path) -> dict:
    """
    Run ChatSpatial's full pipeline: LLM selects tool+params via schema,
    then deterministic execution.

    This mirrors the real ChatSpatial workflow: natural-language prompt →
    schema-constrained JSON → Pydantic validation → tool execution.
    """
    import anndata as ad
    import requests
    from chatspatial.models.data import (
        SpatialDomainParameters,
        SpatialVariableGenesParameters,
        SpatialStatisticsParameters,
        CellCommunicationParameters,
    )
    from chatspatial.tools.spatial_domains import identify_spatial_domains
    from chatspatial.tools.spatial_genes import identify_spatial_genes
    from chatspatial.tools.spatial_statistics import analyze_spatial_statistics
    from chatspatial.tools.cell_communication import analyze_cell_communication

    # Import schema and parsing from ablation scripts
    sys.path.insert(0, str(SCRIPT_DIR))
    from ablation_invocation import SCHEMA_FULL, RESPONSE_INSTRUCTION, parse_json, resolve_tool_name

    tool_map = {
        "identify_spatial_domains": (SpatialDomainParameters, identify_spatial_domains),
        "find_spatial_genes": (SpatialVariableGenesParameters, identify_spatial_genes),
        "analyze_spatial_statistics": (SpatialStatisticsParameters, analyze_spatial_statistics),
        "analyze_cell_communication": (CellCommunicationParameters, analyze_cell_communication),
    }

    # Minimal context (same as AblationCtx from ablation_e2e.py)
    class BenchCtx:
        def __init__(self, datasets):
            self._datasets = datasets
        async def get_adata(self, data_id):
            return self._datasets[data_id]
        async def get_dataset_info(self, data_id):
            return {"adata": self._datasets[data_id], "name": data_id, "type": "visium"}
        async def set_adata(self, data_id, adata):
            self._datasets[data_id] = adata
        async def add_dataset(self, adata, prefix="data", name=None, metadata=None):
            new_id = f"{prefix}_{len(self._datasets)}"
            self._datasets[new_id] = adata
            return new_id
        async def info(self, msg): pass
        async def warning(self, msg): pass
        async def error(self, msg): pass
        def debug(self, msg): pass
        def log_config(self, title, config): pass

    result = {
        "success": False,
        "error": "",
        "wall_time": 0.0,
        "method": "",
        "output_files": [],
        "llm_response": "",
        "parsed_params": {},
    }

    t0 = time.time()
    try:
        # Step 1: LLM selects tool and parameters via schema
        user_prompt = task["prompt_template"].format(
            data_path=str(data_path),
            output_path=str(output_dir / task["output_file"]),
        )

        payload = {
            "model": MODEL,
            "max_tokens": 512,
            "temperature": TEMPERATURE,
            "system": SCHEMA_FULL + RESPONSE_INSTRUCTION,
            "messages": [{"role": "user", "content": user_prompt}],
        }
        r = requests.post(
            "https://api.anthropic.com/v1/messages",
            headers={
                "x-api-key": ANTHROPIC_KEY,
                "anthropic-version": "2023-06-01",
                "content-type": "application/json",
            },
            json=payload, timeout=60,
        )
        r.raise_for_status()
        llm_text = r.json()["content"][0]["text"]
        result["llm_response"] = llm_text

        # Step 2: Parse JSON response
        parsed = parse_json(llm_text)
        if parsed is None:
            result["error"] = "json_parse_failed"
            result["wall_time"] = time.time() - t0
            return result

        tool_name = resolve_tool_name(parsed.get("tool_name"))
        params_dict = parsed.get("parameters", {})
        result["parsed_params"] = params_dict
        result["invoked_tool"] = tool_name

        if tool_name not in tool_map:
            result["error"] = f"unknown_tool:{tool_name}"
            result["wall_time"] = time.time() - t0
            return result

        param_cls, tool_func = tool_map[tool_name]

        # Step 3: Pydantic validation
        params = param_cls(**params_dict)
        result["method"] = getattr(params, "method", "")

        # Step 4: Execute tool on fresh data copy
        adata = ad.read_h5ad(data_path)
        adata.var_names_make_unique()
        datasets = {"bench_data": adata}
        ctx = BenchCtx(datasets)

        tool_result = await tool_func("bench_data", ctx, params)
        result["wall_time"] = time.time() - t0

        # Save outputs
        adata_out = await ctx.get_adata("bench_data")

        if tool_name == "identify_spatial_domains":
            out_path = output_dir / "result.h5ad"
            adata_out.write_h5ad(out_path)
            result["output_files"].append(str(out_path))
            result["domain_key"] = tool_result.domain_key
            result["n_domains"] = tool_result.n_domains
            if tool_result.domain_key in adata_out.obs:
                result["labels"] = adata_out.obs[tool_result.domain_key].values.tolist()

        elif tool_name == "find_spatial_genes":
            import pandas as pd
            genes = tool_result.spatial_genes[:100]
            out_path = output_dir / "svg_results.csv"
            pd.DataFrame({"gene": genes}).to_csv(out_path, index=False)
            result["output_files"].append(str(out_path))
            result["spatial_genes"] = genes
            result["n_significant"] = tool_result.n_significant_genes

        elif tool_name == "analyze_spatial_statistics":
            # Moran's I / spatial autocorrelation produces gene-level stats
            import pandas as pd
            genes = tool_result.top_genes[:100] if hasattr(tool_result, "top_genes") else []
            out_path = output_dir / "svg_results.csv"
            pd.DataFrame({"gene": genes}).to_csv(out_path, index=False)
            result["output_files"].append(str(out_path))
            result["spatial_genes"] = genes

        elif tool_name == "analyze_cell_communication":
            import pandas as pd
            lr_pairs = tool_result.top_lr_pairs[:50]
            out_path = output_dir / "ccc_results.csv"
            rows = []
            for pair in lr_pairs:
                parts = pair.split("^", 1) if "^" in pair else pair.split("_", 1)
                ligand = parts[0] if len(parts) > 0 else pair
                receptor = parts[1] if len(parts) > 1 else ""
                rows.append({"ligand": ligand, "receptor": receptor})
            pd.DataFrame(rows).to_csv(out_path, index=False)
            result["output_files"].append(str(out_path))
            result["lr_pairs"] = lr_pairs
            result["n_significant"] = tool_result.n_significant_pairs

        result["success"] = True

    except Exception as e:
        result["wall_time"] = time.time() - t0
        result["error"] = f"{type(e).__name__}: {str(e)[:500]}"
        import traceback as tb
        result["traceback"] = tb.format_exc()[-1000:]

    return result


# ---------------------------------------------------------------------------
# STAgent / SpatialAgent Execution (subprocess into separate venvs)
# ---------------------------------------------------------------------------

def run_external_system(
    system: str,
    venv_python: Path,
    driver_script: Path,
    task: dict,
    output_dir: Path,
    data_path: Path,
) -> dict:
    """Run STAgent or SpatialAgent via subprocess."""
    result_json = output_dir / "driver_result.json"
    out_path = output_dir / task["output_file"]

    prompt = task["prompt_template"].format(
        data_path=str(data_path),
        output_path=str(out_path),
    )

    cmd = [
        str(venv_python),
        str(driver_script),
        "--data_path", str(data_path),
        "--output_dir", str(output_dir),
        "--prompt", prompt,
        "--model", MODEL,
        "--result_json", str(result_json),
    ]
    # SpatialAgent driver accepts --temperature
    if system == "spatialagent":
        cmd.extend(["--temperature", str(TEMPERATURE)])

    env = os.environ.copy()
    env["ANTHROPIC_API_KEY"] = ANTHROPIC_KEY
    env["PYTHONUNBUFFERED"] = "1"
    # STAgent needs this for chromadb/protobuf compatibility
    env["PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION"] = "python"
    # STAgent initializes all providers on startup; set placeholder keys for unused ones
    env.setdefault("GOOGLE_API_KEY", "not-used")
    env.setdefault("OPENAI_API_KEY", "not-used")

    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SECONDS,
            env=env,
            cwd=str(output_dir),
        )
        wall_time = time.time() - t0

        # Read driver result
        if result_json.exists():
            with open(result_json) as f:
                result = json.load(f)
            result["wall_time"] = wall_time
            result["returncode"] = proc.returncode
        else:
            result = {
                "success": False,
                "error": f"Driver did not produce result JSON. "
                         f"returncode={proc.returncode}",
                "wall_time": wall_time,
                "stdout": proc.stdout[-2000:] if proc.stdout else "",
                "stderr": proc.stderr[-2000:] if proc.stderr else "",
            }

        # Capture stdout/stderr for debugging
        if proc.stdout:
            result["stdout_tail"] = proc.stdout[-1000:]
        if proc.stderr:
            result["stderr_tail"] = proc.stderr[-1000:]

    except subprocess.TimeoutExpired:
        result = {
            "success": False,
            "error": f"timeout_after_{TIMEOUT_SECONDS}s",
            "wall_time": time.time() - t0,
        }
    except Exception as e:
        result = {
            "success": False,
            "error": f"{type(e).__name__}: {str(e)[:500]}",
            "wall_time": time.time() - t0,
        }

    return result


# ---------------------------------------------------------------------------
# Main Experiment Loop
# ---------------------------------------------------------------------------

async def run_experiment():
    print("=" * 70)
    print("End-to-End Benchmark: ChatSpatial vs STAgent vs SpatialAgent")
    print("=" * 70)
    print(f"Dataset: {DATASET_PATH}")
    print(f"Model: {MODEL}")
    print(f"Temperature: {TEMPERATURE}")
    print(f"Reps: {N_REPS}")
    print(f"Tasks: {[t['id'] for t in TASKS]}")
    print(f"Total trials: {len(TASKS) * 3 * N_REPS}")
    print()

    if not DATASET_PATH.exists():
        print(f"ERROR: Dataset not found: {DATASET_PATH}")
        sys.exit(1)

    if not ANTHROPIC_KEY:
        print("ERROR: ANTHROPIC_API_KEY not set")
        sys.exit(1)

    completed = load_completed()
    print(f"Completed trials: {len(completed)}")

    systems = [
        ("chatspatial", None, None),
        ("stagent", STAGENT_VENV, SCRIPT_DIR / "e2e_driver_stagent.py"),
        ("spatialagent", SPATIALAGENT_VENV, SCRIPT_DIR / "e2e_driver_spatialagent.py"),
    ]

    total = len(TASKS) * len(systems) * N_REPS
    done = 0

    for task in TASKS:
        for system_name, venv, driver in systems:
            for rep in range(N_REPS):
                trial_key = (system_name, task["id"], rep)
                if trial_key in completed:
                    done += 1
                    continue

                done += 1
                print(f"\n[{done}/{total}] {system_name} / {task['id']} / rep {rep}")

                # Create isolated output directory
                trial_dir = OUTPUT_DIR / system_name / task["id"] / f"rep_{rep}"
                trial_dir.mkdir(parents=True, exist_ok=True)

                dp = DATASET_PATH_CCC if task["id"] == "cell_communication" else DATASET_PATH

                if system_name == "chatspatial":
                    result = await run_chatspatial(task, trial_dir, dp)
                else:
                    if not venv.exists():
                        result = {
                            "success": False,
                            "error": f"venv not found: {venv}",
                            "wall_time": 0.0,
                        }
                    else:
                        result = run_external_system(
                            system_name, venv, driver, task, trial_dir, dp,
                        )

                # Build record
                record = {
                    "system": system_name,
                    "task_id": task["id"],
                    "rep": rep,
                    "model": MODEL,
                    "temperature": TEMPERATURE,
                    "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
                    **result,
                }

                # Remove large fields before JSONL
                for key in ("messages_text", "stdout", "stderr",
                            "stdout_tail", "stderr_tail", "traceback",
                            "labels"):
                    if key in record and isinstance(record.get(key), str) and len(record[key]) > 2000:
                        record[key] = record[key][:2000] + "..."

                append_raw(record)
                status = "OK" if result.get("success") else f"FAIL: {result.get('error', '?')[:80]}"
                print(f"  -> {status} ({result.get('wall_time', 0):.1f}s)")

                # Brief delay between API calls
                if done < total:
                    time.sleep(2.0)

    print(f"\nDone. Results in {RAW_PATH}")


if __name__ == "__main__":
    asyncio.run(run_experiment())
