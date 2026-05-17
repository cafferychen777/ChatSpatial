#!/usr/bin/env python3
"""
DLPFC Ground-Truth Benchmark: ChatSpatial vs STAgent vs SpatialAgent

Runs all three systems on DLPFC Visium samples with expert-annotated
cortical layers as ground truth. Uses an open-ended prompt (no method
or parameter hints) to test end-to-end analytical capability.

Trial matrix: 3 samples x 3 systems x 10 reps = 90 runs.

Ground truth labels are stored separately and never exposed to systems.
ARI/NMI against ground truth is computed by dlpfc_benchmark_analysis.py.
"""

import asyncio
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from dotenv import load_dotenv

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent
PAPER_ROOT = SCRIPT_DIR.parent
CODE_ROOT = PAPER_ROOT.parent / "code"
BENCH_ROOT = PAPER_ROOT.parent / "benchmarks"

DATA_DIR = PAPER_ROOT / "data" / "dlpfc_benchmark"
SAMPLES_DIR = DATA_DIR / "samples"
OUTPUT_DIR = DATA_DIR / "outputs"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RAW_PATH = DATA_DIR / "dlpfc_benchmark_raw.jsonl"

# Load .env from workspace root
load_dotenv(PAPER_ROOT.parent / ".env")

# System venvs (same as e2e_benchmark.py)
STAGENT_VENV = BENCH_ROOT / "STAgent" / ".venv" / "bin" / "python"
SPATIALAGENT_VENV = BENCH_ROOT / "SpatialAgent" / ".venv" / "bin" / "python"

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

N_REPS = 10
MODEL = "claude-sonnet-4-20250514"
TEMPERATURE = 1.0
TIMEOUT_SECONDS = 900  # 15 minutes per run

ANTHROPIC_KEY = os.environ.get("ANTHROPIC_API_KEY", "")

# Target samples
SAMPLE_IDS = ["151673", "151507", "151669"]

# Open-ended prompt — no method, resolution, or n_domains specified
DLPFC_PROMPT = (
    "This is human dorsolateral prefrontal cortex (DLPFC) Visium data. "
    "Identify spatial domains in this dataset."
)

# ---------------------------------------------------------------------------
# Ensure chatspatial importable
# ---------------------------------------------------------------------------

if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))

# ---------------------------------------------------------------------------
# Checkpoint / Resume
# ---------------------------------------------------------------------------

def load_completed():
    """Load set of completed (system, sample_id, rep) tuples."""
    completed = set()
    if RAW_PATH.exists():
        with open(RAW_PATH) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    completed.add((rec["system"], rec["sample_id"], rec["rep"]))
                except (json.JSONDecodeError, KeyError):
                    pass
    return completed


def append_raw(record: dict):
    with open(RAW_PATH, "a") as f:
        f.write(json.dumps(record, default=str) + "\n")


# ---------------------------------------------------------------------------
# ChatSpatial Execution
# ---------------------------------------------------------------------------

async def run_chatspatial(sample_id: str, output_dir: Path, data_path: Path) -> dict:
    """Run ChatSpatial: LLM selects tool+params via schema, then execute."""
    import anndata as ad
    import requests
    from chatspatial.models.data import SpatialDomainParameters
    from chatspatial.tools.spatial_domains import identify_spatial_domains

    sys.path.insert(0, str(SCRIPT_DIR))
    from ablation_invocation import SCHEMA_FULL, RESPONSE_INSTRUCTION, parse_json, resolve_tool_name

    # Minimal context
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
        # Step 1: LLM selects tool and parameters
        payload = {
            "model": MODEL,
            "max_tokens": 512,
            "temperature": TEMPERATURE,
            "system": SCHEMA_FULL + RESPONSE_INSTRUCTION,
            "messages": [{"role": "user", "content": DLPFC_PROMPT}],
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

        if tool_name != "identify_spatial_domains":
            result["error"] = f"wrong_tool:{tool_name}"
            result["wall_time"] = time.time() - t0
            return result

        # Step 3: Pydantic validation
        params = SpatialDomainParameters(**params_dict)
        result["method"] = params.method
        result["n_domains"] = params.n_domains
        result["resolution"] = params.resolution

        # Step 4: Execute on fresh data copy
        adata = ad.read_h5ad(data_path)
        adata.var_names_make_unique()
        datasets = {"bench_data": adata}
        ctx = BenchCtx(datasets)

        tool_result = await identify_spatial_domains("bench_data", ctx, params)
        result["wall_time"] = time.time() - t0

        # Save outputs
        adata_out = await ctx.get_adata("bench_data")
        out_path = output_dir / "result.h5ad"
        adata_out.write_h5ad(out_path)
        result["output_files"].append(str(out_path))
        result["domain_key"] = tool_result.domain_key
        result["n_domains_found"] = tool_result.n_domains

        # Store labels in the record for quick analysis
        if tool_result.domain_key in adata_out.obs:
            labels = adata_out.obs[tool_result.domain_key].values.tolist()
            result["labels"] = labels

        result["success"] = True

    except Exception as e:
        result["wall_time"] = time.time() - t0
        result["error"] = f"{type(e).__name__}: {str(e)[:500]}"
        import traceback as tb
        result["traceback"] = tb.format_exc()[-1000:]

    return result


# ---------------------------------------------------------------------------
# External System Execution (subprocess)
# ---------------------------------------------------------------------------

def run_external_system(
    system: str,
    venv_python: Path,
    driver_script: Path,
    sample_id: str,
    output_dir: Path,
    data_path: Path,
) -> dict:
    """Run STAgent or SpatialAgent via subprocess."""
    import subprocess

    result_json = output_dir / "driver_result.json"
    out_path = output_dir / "result.h5ad"

    prompt = (
        f"Load the Visium dataset at '{data_path}'. "
        f"This is human dorsolateral prefrontal cortex (DLPFC) tissue. "
        f"Identify spatial domains in this dataset. "
        f"Save the resulting AnnData object (with cluster labels) "
        f"to '{out_path}'."
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
    if system == "spatialagent":
        cmd.extend(["--temperature", str(TEMPERATURE)])

    env = os.environ.copy()
    env["ANTHROPIC_API_KEY"] = ANTHROPIC_KEY
    env["PYTHONUNBUFFERED"] = "1"
    env["PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION"] = "python"
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
    print("DLPFC Ground-Truth Benchmark: ChatSpatial vs STAgent vs SpatialAgent")
    print("=" * 70)
    print(f"Samples: {SAMPLE_IDS}")
    print(f"Model: {MODEL}")
    print(f"Temperature: {TEMPERATURE}")
    print(f"Reps: {N_REPS}")
    print(f"Total trials: {len(SAMPLE_IDS) * 3 * N_REPS}")
    print()

    # Verify samples exist
    for sid in SAMPLE_IDS:
        p = SAMPLES_DIR / f"{sid}.h5ad"
        if not p.exists():
            print(f"ERROR: Sample not found: {p}")
            print("Run dlpfc_prepare_data.py first.")
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

    total = len(SAMPLE_IDS) * len(systems) * N_REPS
    done = 0

    for sample_id in SAMPLE_IDS:
        data_path = SAMPLES_DIR / f"{sample_id}.h5ad"

        for system_name, venv, driver in systems:
            for rep in range(N_REPS):
                trial_key = (system_name, sample_id, rep)
                if trial_key in completed:
                    done += 1
                    continue

                done += 1
                print(f"\n[{done}/{total}] {system_name} / {sample_id} / rep {rep}")

                trial_dir = OUTPUT_DIR / system_name / sample_id / f"rep_{rep}"
                trial_dir.mkdir(parents=True, exist_ok=True)

                if system_name == "chatspatial":
                    result = await run_chatspatial(sample_id, trial_dir, data_path)
                else:
                    if not venv.exists():
                        result = {
                            "success": False,
                            "error": f"venv not found: {venv}",
                            "wall_time": 0.0,
                        }
                    else:
                        result = run_external_system(
                            system_name, venv, driver,
                            sample_id, trial_dir, data_path,
                        )

                record = {
                    "system": system_name,
                    "sample_id": sample_id,
                    "rep": rep,
                    "model": MODEL,
                    "temperature": TEMPERATURE,
                    "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
                    **result,
                }

                # Truncate large fields
                for key in ("messages_text", "stdout", "stderr",
                            "stdout_tail", "stderr_tail", "traceback",
                            "labels"):
                    val = record.get(key)
                    if isinstance(val, str) and len(val) > 2000:
                        record[key] = val[:2000] + "..."
                    elif isinstance(val, list) and len(str(val)) > 5000:
                        record[key] = f"[{len(val)} items, truncated]"

                append_raw(record)
                status = "OK" if result.get("success") else f"FAIL: {result.get('error', '?')[:80]}"
                print(f"  -> {status} ({result.get('wall_time', 0):.1f}s)")

                if done < total:
                    time.sleep(2.0)

    print(f"\nDone. Results in {RAW_PATH}")


if __name__ == "__main__":
    asyncio.run(run_experiment())
