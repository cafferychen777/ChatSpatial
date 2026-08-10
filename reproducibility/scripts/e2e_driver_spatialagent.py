#!/usr/bin/env python3
"""
SpatialAgent driver — runs a single trial inside the SpatialAgent venv.

Usage:
    .venv/bin/python e2e_driver_spatialagent.py \
        --data_path /path/to/data.h5ad \
        --output_dir /path/to/output \
        --prompt "..." \
        --model claude-sonnet-4-20250514 \
        --result_json /path/to/result.json

Writes a JSON file with:
  {"success": bool, "error": str, "wall_time": float, "n_steps": int,
   "messages_text": str, "output_files": [...]}
"""

import argparse
import json
import os
import sys
import time
import traceback

from paths import find_competitor_dir


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--prompt", required=True)
    parser.add_argument("--model", default="claude-sonnet-4-20250514")
    parser.add_argument("--result_json", required=True)
    parser.add_argument("--recursion_limit", type=int, default=50)
    parser.add_argument("--temperature", type=float, default=1.0)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Add SpatialAgent to path. Override with SPATIALAGENT_ROOT when needed.
    sa_root = str(
        find_competitor_dir("SpatialAgent", "SPATIALAGENT_ROOT", required=True)
    )
    if sa_root not in sys.path:
        sys.path.insert(0, sa_root)

    result = {
        "success": False,
        "error": "",
        "wall_time": 0.0,
        "n_steps": 0,
        "messages_text": "",
        "output_files": [],
    }

    t0 = time.time()
    try:
        from spatialagent.agent import SpatialAgent, make_llm

        llm = make_llm(
            args.model,
            temperature=args.temperature,
            track_cost=False,
        )

        agent = SpatialAgent(
            llm=llm,
            data_path=os.path.dirname(args.data_path),
            save_path=args.output_dir,
            tool_retrieval=True,
            tool_retrieval_method="llm",
            auto_interpret_figures=False,
            act_timeout=600,
        )

        out = agent.run(
            user_query=args.prompt,
            config={"recursion_limit": args.recursion_limit},
        )

        result["wall_time"] = time.time() - t0

        # Extract messages
        out_msgs = out.get("messages", [])
        result["n_steps"] = len(out_msgs)
        result["next_step"] = out.get("next_step", "unknown")

        # Collect message content
        texts = []
        for msg in out_msgs:
            content = getattr(msg, "content", "")
            if isinstance(content, str) and content.strip():
                texts.append(content[:2000])  # Truncate long messages
        result["messages_text"] = "\n---\n".join(texts[-10:])

        # List output files
        output_files = []
        for root, dirs, files in os.walk(args.output_dir):
            for f in files:
                rel = os.path.relpath(os.path.join(root, f), args.output_dir)
                output_files.append(rel)
        result["output_files"] = output_files

        result["success"] = True

    except Exception as e:
        result["wall_time"] = time.time() - t0
        result["error"] = f"{type(e).__name__}: {str(e)[:500]}"
        result["traceback"] = traceback.format_exc()[-1000:]

    with open(args.result_json, "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(
        f"SpatialAgent driver done. success={result['success']}, "
        f"wall_time={result['wall_time']:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
