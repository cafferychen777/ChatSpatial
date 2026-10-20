"""Gene classes that analyses can exclude by name (immunoglobulin, ribosomal, ...).

Shared by analyze_cnv and find_markers so both tools exclude exactly the same
genes for the same class.
"""

import re
from collections.abc import Iterable
from fnmatch import fnmatchcase

import numpy as np
import pandas as pd

# Name rules of the gene classes, matched case-insensitively so that the mouse
# spellings (Ighm, Rpl13, mt-Co1, H2-K1) follow the human ones.
GENE_CLASS_PATTERNS = {
    # V, D, J and constant gene segments of the IGH, IGK and IGL loci (and
    # their pseudogenes), heavy-chain constant genes, plus the protein-coding
    # genes of these loci expressed by B and plasma cells.
    "immunoglobulin": re.compile(
        r"^(?:IG[HKL][VDJC]|IGH[AGMDE]|JCHAIN$|IGLL[15]$)", re.IGNORECASE
    ),
    "ribosomal": re.compile(r"^RP[LS]", re.IGNORECASE),
    "mitochondrial": re.compile(r"^MT-", re.IGNORECASE),
    "hla": re.compile(r"^(?:HLA-|H2-)", re.IGNORECASE),
}
# Genes that match a class prefix but do not belong to the class: a helicase in
# the IGH name space and the ribosomal protein S6 kinases.
_GENE_CLASS_EXCEPTIONS = {
    "immunoglobulin": re.compile(r"^IGHMBP2$", re.IGNORECASE),
    "ribosomal": re.compile(r"^(?:RPS6K|RPS19BP1$)", re.IGNORECASE),
}


def gene_class_mask(gene_names: Iterable[str], gene_class: str) -> np.ndarray:
    """Return a boolean mask of the genes that belong to a gene class."""
    names = pd.Series(np.asarray(list(gene_names), dtype=str))
    mask = names.str.contains(GENE_CLASS_PATTERNS[gene_class])
    exception = _GENE_CLASS_EXCEPTIONS.get(gene_class)
    if exception is not None:
        mask &= ~names.str.contains(exception)
    return mask.to_numpy(dtype=bool)


def excluded_gene_mask(
    gene_names: Iterable[str],
    gene_classes: Iterable[str] = (),
    gene_patterns: Iterable[str] = (),
) -> tuple[np.ndarray, dict[str, int]]:
    """Combine gene classes and case-sensitive wildcard patterns into one mask.

    Returns the mask and the number of genes each class or pattern matches
    (a gene matched by several entries is counted under each). Patterns are
    reported as ``"pattern:<pattern>"``.
    """
    names = np.asarray(list(gene_names), dtype=str)
    excluded = np.zeros(len(names), dtype=bool)
    counts: dict[str, int] = {}
    for gene_class in dict.fromkeys(gene_classes):
        mask = gene_class_mask(names, gene_class)
        counts[gene_class] = int(mask.sum())
        excluded |= mask
    for pattern in dict.fromkeys(gene_patterns):
        mask = np.fromiter(
            (fnmatchcase(name, pattern) for name in names), bool, len(names)
        )
        counts[f"pattern:{pattern}"] = int(mask.sum())
        excluded |= mask
    return excluded, counts
