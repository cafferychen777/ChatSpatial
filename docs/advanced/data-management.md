# Data Management

This project treats large datasets as **external assets**, not repository contents.

## Principles

1. Keep Git history lightweight and reviewable.
2. Keep tests reproducible with small synthetic fixtures.
3. Keep real datasets in external storage (local cache, object store, or shared filesystem).

## Dataset Registry

Use `scripts/data/dataset_registry.json` as the source of truth for downloadable datasets.

- `id`: stable dataset identifier
- `url`: download URL
- `filename`: local filename under dataset folder
- `sha256`: optional integrity check
- `extract`: auto-extract archives when true
- `enabled`: include/exclude from fetch runs

## Fetch Datasets

```bash
# list enabled datasets
scripts/data/fetch_datasets.py --list

# fetch selected datasets to default cache (~/.cache/chatspatial/datasets)
scripts/data/fetch_datasets.py --dataset your_dataset_id

# fetch all enabled datasets to a custom location
scripts/data/fetch_datasets.py --dest /data/chatspatial
```

## Register Local Paths (Team/Personal)

For datasets that should not be downloaded publicly:

```bash
scripts/data/register_external_dataset.py my_dataset /absolute/path/to/data.h5ad
```

This creates `data/datasets.local.json` (gitignored).

## Workspace Cleanup

Clean local build/test artifacts without touching source files:

```bash
scripts/maintenance/clean_workspace.sh
```

## Testing Contract

- Default test suite must not depend on `data/`.
- Use fixtures that generate temporary `.h5ad` files under `tmp_path`.
- Heavy dependency checks are marked `@pytest.mark.slow` and excluded by default.
