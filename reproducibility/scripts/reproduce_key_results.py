#!/usr/bin/env python3
"""One-command, LLM-free reproduction of key manuscript results.

Re-executes the tool calls that the language models recorded during the
manuscript experiments, without calling any LLM, and compares the regenerated
numbers with the reference values committed under ``reproducibility/data``:

  dlpfc   DLPFC ground-truth ARI (slices 151673, 151507, 151669; ChatSpatial
          arm of ``dlpfc_benchmark.py``; manuscript pooled ARI 0.374).
  card    OSCC CARD cross-model concordance (full-schema arm of
          ``casestudy_reproducibility.py``; four models, Pearson r = 1.000),
          plus the schema-validation outcome of every recorded call and the
          deviation from the recorded proportion matrices.
  replay  Generic replay of an additional recorded tool-call log given with
          ``--replay-log`` (JSONL, one call per line; see ``run_replay``).
          Skipped when no log is supplied.

Analysis code is reused from the original experiment scripts; this file only
replays the recorded calls and compares outputs.

Inputs are read from ``--data-dir`` (default ``/data/chatspatial-repro``).
If they are missing and ``--data-url`` (or ``CHATSPATIAL_REPRO_DATA_URL``)
points to the input archive, it is downloaded and unpacked there first.
Expected layout::

    dlpfc/{151673,151507,151669}.h5ad
    oscc/s1_processed.h5ad
    oscc/s9_processed.h5ad
    oscc/puram_hnscc_reference.h5ad
    logs/dlpfc_benchmark_raw.jsonl
    logs/casestudy_repro_raw.jsonl
    recorded/casestudy_proportions/*.json
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import platform
import shutil
import sys
import tarfile
import tempfile
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPRO_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

DEFAULT_DATA_DIR = Path(
    os.environ.get("CHATSPATIAL_REPRO_DATA_DIR", "/data/chatspatial-repro")
)
DEFAULT_OUT_DIR = Path(os.environ.get("CHATSPATIAL_OUTPUT_DIR", "outputs")) / (
    "reproduction"
)
DLPFC_SAMPLES = ["151673", "151507", "151669"]
OSCC_SAMPLES = {"oscc_s1": "s1_processed.h5ad", "oscc_s9": "s9_processed.h5ad"}


def log(msg: str) -> None:
    print(msg, flush=True)


def read_jsonl(path: Path) -> list[dict]:
    with open(path) as handle:
        return [json.loads(line) for line in handle if line.strip()]


def canonical(params) -> str:
    """Stable key for a validated Pydantic parameter object."""
    return json.dumps(params.model_dump(mode="json"), sort_keys=True, default=str)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def ensure_inputs(data_dir: Path, data_url: str | None, steps: list[str]) -> None:
    needed = []
    if "dlpfc" in steps:
        needed += [f"dlpfc/{s}.h5ad" for s in DLPFC_SAMPLES]
        needed += ["logs/dlpfc_benchmark_raw.jsonl"]
    if "card" in steps:
        needed += [f"oscc/{f}" for f in OSCC_SAMPLES.values()]
        needed += ["oscc/puram_hnscc_reference.h5ad", "logs/casestudy_repro_raw.jsonl"]
    missing = [p for p in needed if not (data_dir / p).exists()]
    if not missing:
        return
    if not data_url:
        raise SystemExit(
            f"Missing inputs under {data_dir}: {missing[:3]}... Pass --data-url "
            "(or set CHATSPATIAL_REPRO_DATA_URL) or mount the unpacked archive."
        )
    data_dir.mkdir(parents=True, exist_ok=True)
    log(f"Downloading inputs from {data_url} ...")
    t0 = time.time()
    with tempfile.NamedTemporaryFile(suffix=".tar.gz", dir=data_dir) as tmp:
        with urllib.request.urlopen(data_url) as resp:
            shutil.copyfileobj(resp, tmp)
        tmp.flush()
        with tarfile.open(tmp.name) as tar:
            tar.extractall(data_dir, filter="data")
    log(f"  done in {time.time() - t0:.0f}s")
    still = [p for p in needed if not (data_dir / p).exists()]
    if still:
        raise SystemExit(f"Archive did not contain: {still}")


# ---------------------------------------------------------------------------
# Step 1: DLPFC ARI
# ---------------------------------------------------------------------------


async def run_dlpfc(data_dir: Path, out_dir: Path, every_trial: bool) -> list[dict]:
    import anndata as ad
    import numpy as np
    import pandas as pd

    from ablation_e2e import AblationCtx
    from dlpfc_benchmark_analysis import compute_gt_metrics, pairwise_ari

    from chatspatial.models.data import SpatialDomainParameters
    from chatspatial.tools.spatial_domains import identify_spatial_domains

    # The manuscript ran the original SpaGCN 1.2.7, which uses the sparse
    # ``.A`` alias removed in SciPy 1.14. ChatSpatial runs SpaGCN in a child
    # process, so the alias is restored through sitecustomize on PYTHONPATH.
    compat = str(SCRIPT_DIR / "compat")
    os.environ["PYTHONPATH"] = os.pathsep.join(
        [compat] + [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p]
    )
    sys.path.insert(0, compat)
    import sitecustomize  # noqa: F401  (also patches this process)

    ref = pd.read_csv(REPRO_ROOT / "data/dlpfc_benchmark/dlpfc_benchmark_results.csv")
    ref = ref[ref.system == "chatspatial"]
    ref = ref.set_index(ref.sample_id.astype(str))
    trials = [
        r
        for r in read_jsonl(data_dir / "logs/dlpfc_benchmark_raw.jsonl")
        if r["system"] == "chatspatial" and r["success"]
    ]
    log(f"[dlpfc] {len(trials)} recorded ChatSpatial calls")

    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for t in trials:
        if t.get("invoked_tool") != "identify_spatial_domains":
            raise RuntimeError(f"unexpected tool in recorded call: {t}")
        params = SpatialDomainParameters(**t["parsed_params"])
        key = (str(t["sample_id"]), canonical(params))
        groups[key].append(t)

    labels_by_trial: dict[tuple[str, int], dict] = {}
    for (sample_id, key), members in sorted(groups.items()):
        params = SpatialDomainParameters(**json.loads(key))
        runs = members if every_trial else members[:1]
        log(
            f"[dlpfc] {sample_id}: {len(members)} recorded calls -> 1 validated "
            f"parameter set; executing {len(runs)}x ({params.method})"
        )
        for t in runs:
            adata = ad.read_h5ad(data_dir / f"dlpfc/{sample_id}.h5ad")
            adata.var_names_make_unique()
            ctx = AblationCtx({"bench_data": adata})
            t0 = time.time()
            res = await identify_spatial_domains("bench_data", ctx, params)
            out = await ctx.get_adata("bench_data")
            labels = dict(zip(out.obs_names, out.obs[res.domain_key].astype(str)))
            log(f"          rep {t['rep']}: {time.time() - t0:.1f}s")
            targets = members if not every_trial else [t]
            for m in targets:
                labels_by_trial[(sample_id, m["rep"])] = labels

    rows, pooled = [], []
    per_trial = []
    for sample_id in DLPFC_SAMPLES:
        gt = pd.read_csv(
            REPRO_ROOT / f"data/dlpfc_benchmark/ground_truth/{sample_id}_labels.csv",
            index_col=0,
        )["ground_truth"]
        sample_labels = [
            v for (s, _), v in sorted(labels_by_trial.items()) if s == sample_id
        ]
        aris = [compute_gt_metrics(lab, gt)["ari"] for lab in sample_labels]
        per_trial += [(sample_id, a) for a in aris]
        pooled += aris
        mean_ari = float(np.mean(aris))
        cross = pairwise_ari(sample_labels) if every_trial else []
        rows.append(
            {
                "result": "DLPFC ARI",
                "item": sample_id,
                "reference": float(ref.loc[sample_id, "gt_ari_mean"]),
                "reproduced": mean_ari,
                "n": len(aris),
                "note": (f"min cross-replicate ARI {min(cross):.6f}" if cross else ""),
            }
        )
    pd.DataFrame(per_trial, columns=["sample_id", "ari"]).to_csv(
        out_dir / "dlpfc_per_trial_ari.csv", index=False
    )
    ref_pooled = (
        pd.read_csv(REPRO_ROOT / "data/dlpfc_benchmark/dlpfc_benchmark_per_trial.csv")
        .query("system == 'chatspatial' and success")["ari"]
        .mean()
    )
    rows.append(
        {
            "result": "DLPFC ARI",
            "item": "pooled (3 slices x 10 calls)",
            "reference": float(ref_pooled),
            "reproduced": float(np.mean(pooled)),
            "n": len(pooled),
            "note": "manuscript reports 0.374",
        }
    )
    return rows


# ---------------------------------------------------------------------------
# Step 2: CARD cross-model concordance
# ---------------------------------------------------------------------------


async def run_card(data_dir: Path, out_dir: Path) -> list[dict]:
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from ablation_e2e import remap_params
    from casestudy_reproducibility import execute_deconv, pairwise_pearson

    from chatspatial.models.data import DeconvolutionParameters

    ref_metrics = pd.read_csv(
        REPRO_ROOT / "data/casestudy_reproducibility/casestudy_repro_metrics.csv"
    ).set_index(["sample", "condition"])
    trials = read_jsonl(data_dir / "logs/casestudy_repro_raw.jsonl")
    log(
        f"[card] {len(trials)} recorded calls (2 samples x 2 schema conditions x 4 models x 10)"
    )

    # Replay schema validation exactly as the original experiment did.
    valid_count: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for t in trials:
        cell = (t["sample_id"], t["condition"])
        valid_count[cell][1] += 1
        if t["canonical_tool"] != "deconvolve_data":
            continue
        raw = json.loads(t["raw_params"])
        params = None
        for candidate in (raw, remap_params("deconvolve_data", raw)):
            try:
                params = DeconvolutionParameters(**candidate)
                break
            except Exception:
                continue
        if params is None:
            continue
        valid_count[cell][0] += 1
        groups[(t["sample_id"], canonical(params))].append(t)

    rows = []
    for (sample_id, cond), (n_valid, n_total) in sorted(valid_count.items()):
        ref_n = int(ref_metrics.loc[(sample_id, cond), "n_success"])
        rows.append(
            {
                "result": "CARD schema replay",
                "item": f"{sample_id} {cond}: valid calls",
                "reference": ref_n,
                "reproduced": n_valid,
                "n": n_total,
                "note": "",
            }
        )

    log("[card] loading Puram reference ...")
    ref_adata = ad.read_h5ad(data_dir / "oscc/puram_hnscc_reference.h5ad")
    rec_dir = data_dir / "recorded/casestudy_proportions"
    for sample_id, fname in OSCC_SAMPLES.items():
        adata = ad.read_h5ad(data_dir / "oscc" / fname)
        if "highly_variable" not in adata.var.columns:  # same as original script
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            adata.raw = adata
            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            sc.pp.pca(adata)
            sc.pp.neighbors(adata)
        elif "connectivities" not in adata.obsp:
            sc.pp.neighbors(adata)

        sample_groups = {k: v for k, v in groups.items() if k[0] == sample_id}
        matrices, max_dev, n_compared = [], 0.0, 0
        for (_, key), members in sample_groups.items():
            log(
                f"[card] {sample_id}: {len(members)} valid calls -> 1 parameter "
                f"set; executing CARD once"
            )
            t0 = time.time()
            ok, err, art = await execute_deconv(
                json.loads(key), adata, ref_adata, data_id=sample_id
            )
            if not ok:
                raise RuntimeError(f"CARD failed for {sample_id}: {err}")
            log(f"          {time.time() - t0:.0f}s")
            new = pd.DataFrame(art["proportions"], columns=art["cell_types"])
            new.to_csv(out_dir / f"card_proportions_{sample_id}.csv", index=False)
            for m in members:
                matrices.append(new.to_numpy())
                rec_path = rec_dir / Path(m.get("artifact_path", "")).name
                if rec_path.is_file():
                    rec = json.loads(rec_path.read_text())
                    old = pd.DataFrame(rec["proportions"], columns=rec["cell_types"])
                    diff = new[old.columns].to_numpy() - old.to_numpy()
                    max_dev = max(max_dev, float(np.abs(diff).max()))
                    n_compared += 1
        r = pairwise_pearson(matrices)
        rows.append(
            {
                "result": "CARD cross-model concordance",
                "item": f"{sample_id}: mean pairwise Pearson r (4 models x 10)",
                "reference": float(
                    ref_metrics.loc[(sample_id, "full_schema"), "pearson_mean"]
                ),
                "reproduced": float(np.mean(r)),
                "n": len(r),
                "note": "",
            }
        )
        rows.append(
            {
                "result": "CARD proportions vs recorded",
                "item": f"{sample_id}: max |proportion difference|",
                "reference": 0.0,
                "reproduced": max_dev if n_compared else float("nan"),
                "n": n_compared,
                "note": "" if n_compared else "recorded proportions not supplied",
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Step 3: generic replay of a recorded tool-call log
# ---------------------------------------------------------------------------


async def run_replay(
    data_dir: Path, out_dir: Path, replay_log: Path | None
) -> list[dict]:
    """Replay a recorded log without an LLM and report exact-match rates.

    Each JSONL line: {"call_id", "tool": "identify_spatial_domains" |
    "deconvolve_data", "params": {...}, "input": "<path under data dir>",
    "reference_input": "<path, deconvolve_data only>",
    "recorded_output": "<CSV under data dir: one column of domain labels, or a
    spots x cell-types proportion table>"}.
    """
    if replay_log is None:
        log("[replay] SKIPPED: no --replay-log supplied")
        return [
            {
                "result": "Replay",
                "item": "recorded tool-call log",
                "reference": float("nan"),
                "reproduced": float("nan"),
                "n": 0,
                "note": "skipped: no log supplied",
            }
        ]
    import anndata as ad
    import numpy as np
    import pandas as pd
    from ablation_e2e import AblationCtx

    from chatspatial.models.data import DeconvolutionParameters, SpatialDomainParameters
    from chatspatial.tools.deconvolution import deconvolve_spatial_data
    from chatspatial.tools.spatial_domains import identify_spatial_domains

    calls = read_jsonl(replay_log)
    exact, records = 0, []
    for c in calls:
        adata = ad.read_h5ad(data_dir / c["input"])
        adata.var_names_make_unique()
        datasets = {"d": adata}
        recorded = pd.read_csv(data_dir / c["recorded_output"], index_col=0)
        if c["tool"] == "identify_spatial_domains":
            ctx = AblationCtx(datasets)
            res = await identify_spatial_domains(
                "d", ctx, SpatialDomainParameters(**c["params"])
            )
            new = (await ctx.get_adata("d")).obs[res.domain_key].astype(str)
            same = (
                new.reindex(recorded.index).tolist()
                == recorded.iloc[:, 0].astype(str).tolist()
            )
            dev = 0.0 if same else float("nan")
        elif c["tool"] == "deconvolve_data":
            p = dict(c["params"], reference_data_id="ref")
            datasets["ref"] = ad.read_h5ad(data_dir / c["reference_input"])
            ctx = AblationCtx(datasets)
            res = await deconvolve_spatial_data("d", ctx, DeconvolutionParameters(**p))
            new = pd.DataFrame((await ctx.get_adata("d")).obsm[res.proportions_key])
            new.index = adata.obs_names
            dev = float(
                np.abs(
                    new.loc[recorded.index, recorded.columns].to_numpy()
                    - recorded.to_numpy()
                ).max()
            )
            same = dev == 0.0
        else:
            raise ValueError(f"replay does not support tool {c['tool']!r}")
        exact += bool(same)
        records.append({"call_id": c.get("call_id"), "exact": same, "max_abs_dev": dev})
    pd.DataFrame(records).to_csv(out_dir / "replay_calls.csv", index=False)
    return [
        {
            "result": "Replay",
            "item": "calls reproduced exactly",
            "reference": len(calls),
            "reproduced": exact,
            "n": len(calls),
            "note": "",
        }
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def environment() -> dict:
    import chatspatial

    return {
        "chatspatial": chatspatial.__version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "image_digest": os.environ.get("CHATSPATIAL_IMAGE_DIGEST", "unknown"),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--steps", default="dlpfc,card,replay")
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    ap.add_argument("--data-url", default=os.environ.get("CHATSPATIAL_REPRO_DATA_URL"))
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--replay-log", type=Path, default=None)
    ap.add_argument(
        "--every-trial",
        action="store_true",
        help="DLPFC: execute every recorded call instead of each distinct "
        "validated parameter set once (about 10x slower).",
    )
    ap.add_argument(
        "--tolerance",
        type=float,
        default=1e-3,
        help="Maximum absolute deviation from reference counted as reproduced.",
    )
    args = ap.parse_args(argv)
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    args.out.mkdir(parents=True, exist_ok=True)

    import pandas as pd

    t_start = time.time()
    env = environment()
    log(f"Environment: {json.dumps(env)}")
    ensure_inputs(args.data_dir, args.data_url, steps)

    rows, timings = [], {}
    for step in steps:
        t0 = time.time()
        if step == "dlpfc":
            step_rows = asyncio.run(
                run_dlpfc(args.data_dir, args.out, args.every_trial)
            )
        elif step == "card":
            step_rows = asyncio.run(run_card(args.data_dir, args.out))
        elif step == "replay":
            step_rows = asyncio.run(
                run_replay(args.data_dir, args.out, args.replay_log)
            )
        else:
            raise SystemExit(f"unknown step {step!r}")
        timings[step] = round(time.time() - t0, 1)
        rows += step_rows
        # Report each step as soon as it finishes, so a later failure does
        # not hide results that were already reproduced.
        for r in step_rows:
            log(
                f"  -> {r['item']}: reference {r['reference']}, reproduced {r['reproduced']}"
            )
        pd.DataFrame(rows).to_csv(args.out / "reproduction_partial.csv", index=False)

    table = pd.DataFrame(rows)
    table["abs_deviation"] = (
        pd.to_numeric(table["reproduced"], errors="coerce")
        - pd.to_numeric(table["reference"], errors="coerce")
    ).abs()
    # The "vs recorded" rows already hold a deviation in `reproduced`.
    is_dev = table["result"] == "CARD proportions vs recorded"
    table.loc[is_dev, "abs_deviation"] = table.loc[is_dev, "reproduced"]
    checked = table["abs_deviation"].notna()
    table["reproduced_within_tol"] = table["abs_deviation"] <= args.tolerance
    table.to_csv(args.out / "reproduction_summary.csv", index=False)

    total = round(time.time() - t_start, 1)
    meta = dict(
        env, step_seconds=timings, total_seconds=total, tolerance=args.tolerance
    )
    (args.out / "environment.json").write_text(json.dumps(meta, indent=2))

    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        log("\n" + table.drop(columns=["note"]).to_string(index=False))
    max_dev = table.loc[checked, "abs_deviation"].max()
    ok = bool(table.loc[checked, "reproduced_within_tol"].all())
    log(
        f"\nMax absolute deviation: {max_dev:.3g} | total runtime {total / 60:.1f} min | "
        f"{'REPRODUCED' if ok else 'DEVIATION ABOVE TOLERANCE'}"
    )
    log(f"Outputs: {args.out.resolve()}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
