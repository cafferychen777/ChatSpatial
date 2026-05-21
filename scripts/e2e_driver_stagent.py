#!/usr/bin/env python3
"""
STAgent driver — runs a single trial inside the STAgent venv.

Usage:
    .venv/bin/python e2e_driver_stagent.py \
        --data_path /path/to/data.h5ad \
        --output_dir /path/to/output \
        --prompt "..." \
        --model claude-sonnet-4-20250514 \
        --result_json /path/to/result.json

Writes a JSON file with:
  {"success": bool, "error": str, "wall_time": float, "n_steps": int,
   "messages_text": str, "code_blocks": [...]}
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
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Set environment for STAgent. Override with STAGENT_ROOT when needed.
    stagent_root = str(find_competitor_dir("STAgent", "STAGENT_ROOT", required=True))
    src_dir = os.path.join(stagent_root, "src")

    # Load .env from STAgent src/
    env_file = os.path.join(src_dir, ".env")
    if os.path.exists(env_file):
        from dotenv import load_dotenv
        load_dotenv(env_file)

    # Configure STAgent environment variables
    os.environ["STAGENT_DATA_PATH"] = args.data_path
    os.environ["STAGENT_PLOT_DIR"] = args.output_dir
    os.environ.setdefault("STAGENT_RECURSION_LIMIT", str(args.recursion_limit))

    # Add STAgent src to path
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)

    result = {
        "success": False,
        "error": "",
        "wall_time": 0.0,
        "n_steps": 0,
        "messages_text": "",
        "code_blocks": [],
    }

    t0 = time.time()
    try:
        from langchain_core.messages import HumanMessage
        from graph_unified import invoke_our_graph

        messages = [HumanMessage(content=args.prompt)]
        out = invoke_our_graph(
            messages=messages,
            model_name=args.model,
            session_id=f"e2e_bench_{int(time.time())}",
            emit_conflict_message=False,
        )

        result["wall_time"] = time.time() - t0

        # Extract messages
        out_msgs = out.get("messages", [])
        result["n_steps"] = len(out_msgs)

        # Collect all message content
        texts = []
        code_blocks = []
        for msg in out_msgs:
            content = getattr(msg, "content", "")
            if isinstance(content, str):
                texts.append(content)
                # Extract Python code blocks
                import re
                for match in re.finditer(r"```python\n(.*?)```", content, re.DOTALL):
                    code_blocks.append(match.group(1))
            elif isinstance(content, list):
                for item in content:
                    if isinstance(item, dict) and "text" in item:
                        texts.append(item["text"])

        result["messages_text"] = "\n---\n".join(texts[-10:])  # last 10 messages
        result["code_blocks"] = code_blocks[-5:]  # last 5 code blocks

        # Check if output files were created
        output_files = []
        # Directories that STAgent always creates (not real output)
        infra_names = {"conversation_histories", "research_reports",
                       "driver_result.json", ".DS_Store"}
        real_files = []
        for f in os.listdir(args.output_dir):
            output_files.append(f)
            if f not in infra_names:
                fpath = os.path.join(args.output_dir, f)
                if os.path.isfile(fpath):
                    real_files.append(f)
        result["output_files"] = output_files

        # Determine success: must have real output files AND enough steps
        # to indicate actual analysis was performed
        full_text = result["messages_text"].lower()
        has_api_error = any(kw in full_text for kw in [
            "apiconnectionerror", "api connection error",
            "model provider request failed", "rate limit",
        ])
        if has_api_error and not real_files:
            result["success"] = False
            result["error"] = "API error during execution, no output produced"
        elif len(real_files) == 0 and result["n_steps"] <= 3:
            result["success"] = False
            result["error"] = "No real output files produced (only infrastructure dirs)"
        else:
            result["success"] = True

    except Exception as e:
        result["wall_time"] = time.time() - t0
        result["error"] = f"{type(e).__name__}: {str(e)[:500]}"
        result["traceback"] = traceback.format_exc()[-1000:]

    with open(args.result_json, "w") as f:
        json.dump(result, f, indent=2, default=str)

    print(f"STAgent driver done. success={result['success']}, "
          f"wall_time={result['wall_time']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
