"""Shared gene-class annotations for QC and downstream gene selection."""

import numpy as np
import pandas as pd

from .exceptions import DataError


def mitochondrial_mask(adata, gene_names=None) -> tuple[np.ndarray, str]:
    """Resolve mitochondrial identity without changing gene identifiers.

    Boolean ``var['mt']`` annotations take precedence over symbol heuristics.
    Raw annotations cover genes removed from the active matrix; active
    annotations override them for shared genes. Matching is by gene identifier.
    """
    names = pd.Index(adata.var_names if gene_names is None else gene_names)
    mask = np.asarray(names.str.upper().str.startswith(("MT-", "MT.")), dtype=bool)
    source = "gene_symbols"
    annotations = []
    if adata.raw is not None:
        annotations.append((adata.raw.var, "raw.var.mt"))
    annotations.append((adata.var, "var.mt"))
    for var, label in annotations:
        if "mt" not in var:
            continue
        mt = var["mt"]
        if not pd.api.types.is_bool_dtype(mt.dtype) or mt.isna().any():
            raise DataError(
                f"{label} must contain boolean values without missing entries."
            )
        aligned = mt.reindex(names)
        present = aligned.notna().to_numpy()
        if present.any():
            mask[present] = aligned[present].to_numpy(dtype=bool)
            source = label if present.all() else "annotations_and_gene_symbols"
    return mask, source
