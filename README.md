# ChatSpatial Reproducibility

Scripts and instructions for reproducing the experiments in:

> **ChatSpatial: Schema-Enforced Agentic Orchestration for Reproducible and Cross-Platform Spatial Transcriptomics**
>
> Chen Yang, Xianyang Zhang, Jun Chen

## Requirements

- Python 3.11+
- [ChatSpatial](https://github.com/cafferychen777/ChatSpatial) (`pip install chatspatial`)
- `requests` (`pip install requests`)

API keys (set as environment variables for the LLM experiments):

```bash
export GEMINI_API_KEY="..."
export ANTHROPIC_API_KEY="..."
export OPENAI_API_KEY="..."
```

## Scripts

| Script | Description | Paper Section |
|--------|-------------|---------------|
| `scripts/reproducibility_analysis.py` | Static analysis of schema constraint coverage across all 20 MCP tools | Results: Schema Constraint Analysis |
| `scripts/determinism_experiment.py` | Single-model determinism (Gemini 2.0 Flash, 20 trials per prompt, T=1.0) | Results: Single-Model Determinism |
| `scripts/determinism_multimodel.py` | Cross-model determinism (3 LLMs x 8 prompts x 10 reps, T=1.0) | Results: Cross-Model Determinism |
| `scripts/codegen_hallucination.py` | Code generation baseline comparison (same 3 LLMs, free-form code) | Results: Code-Generation Baseline |

### Running

```bash
# Schema constraint analysis (no API keys needed, requires ChatSpatial installed)
python scripts/reproducibility_analysis.py

# Single-model determinism (requires GEMINI_API_KEY)
python scripts/determinism_experiment.py

# Cross-model determinism (requires all 3 API keys)
python scripts/determinism_multimodel.py

# Code generation baseline (requires all 3 API keys)
python scripts/codegen_hallucination.py
```

Output CSVs are written to `data/`.

## Case Studies

The two case study analyses (OSCC and HGSOC) were performed through natural-language conversations with ChatSpatial via the Model Context Protocol. The exact prompts and workflows are described in the main text of the paper (Results: Case Study sections).

To reproduce, install ChatSpatial, connect it to an MCP-compatible client, and follow the analysis steps described in the paper.

## Datasets

| Dataset | Accession | Reference |
|---------|-----------|-----------|
| OSCC (12 Visium samples) | [GSE208253](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208253) | Arora et al., *Nat. Commun.* 2023 |
| HNSCC scRNA-seq reference | [GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322) | Puram et al., *Cell* 2017 |
| HGSOC (8 Visium samples) | [GSE211956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211956) | Denisenko et al., *Nat. Commun.* 2024 |
| SPOTS benchmark | [GSE198353](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198353) | Kleshchevnikov et al., *Nat. Biotechnol.* 2022 |
| Visium multi-sample | [GSE254652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254652), [GSE243275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243275) | See Supplementary Table S2 |
| MERFISH | [GSE113576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576) | Moffitt et al., *Science* 2018 |
| seqFISH | [GSE133244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133244) | Lohoff et al., *Nat. Biotechnol.* 2022 |
| Slide-seq | [SCP354](https://singlecell.broadinstitute.org/single_cell/study/SCP354) | Rodriques et al., *Science* 2019 |
| Xenium | [10x Genomics](https://www.10xgenomics.com/datasets) | 10x Genomics public dataset |

## License

MIT
