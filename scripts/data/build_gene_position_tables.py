#!/usr/bin/env python3
"""Build the gene-position tables bundled for ``analyze_cnv``.

``analyze_cnv`` (infercnvpy) needs the chromosome, start and end of every gene.
When ``adata.var`` lacks them, the tool reads the tables written here from
``chatspatial/resources/gene_positions/``.

Source: Ensembl release 110 GTF gene records, the annotation behind the
10x Genomics 2024-A references (GENCODE v44 human, GENCODE vM33 mouse).

- human: https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
- mouse: https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

Ensembl data are released without restriction on use or redistribution
(https://www.ensembl.org/info/about/legal/disclaimer.html).

Each output table is a gzip TSV with columns ``gene_id`` (Ensembl gene ID
without version), ``gene_name`` (empty when Ensembl gives none), ``chromosome``
(``chr``-prefixed), ``start`` and ``end`` (1-based, inclusive), restricted to
the primary assembled chromosomes (1-22/1-19, X, Y) and sorted by position.

Usage:
    python scripts/data/build_gene_position_tables.py --download-dir /path/to/cache
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from urllib.request import urlopen

ENSEMBL_RELEASE = 110
SOURCES = {
    "human": (
        "homo_sapiens/Homo_sapiens.GRCh38",
        [str(i) for i in range(1, 23)] + ["X", "Y"],
        "human_grch38_ensembl110.tsv.gz",
    ),
    "mouse": (
        "mus_musculus/Mus_musculus.GRCm39",
        [str(i) for i in range(1, 20)] + ["X", "Y"],
        "mouse_grcm39_ensembl110.tsv.gz",
    ),
}
OUT_DIR = (
    Path(__file__).resolve().parents[2] / "chatspatial" / "resources" / "gene_positions"
)

_GENE_ID = re.compile(r'gene_id "([^"]+)"')
_GENE_NAME = re.compile(r'gene_name "([^"]+)"')
_GENE_BIOTYPE = re.compile(r'gene_biotype "([^"]+)"')


def _url(species: str) -> str:
    stem = SOURCES[species][0]
    return (
        f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/gtf/"
        f"{stem}.{ENSEMBL_RELEASE}.gtf.gz"
    )


def _download(url: str, out_file: Path) -> None:
    if out_file.exists():
        return
    out_file.parent.mkdir(parents=True, exist_ok=True)
    partial = out_file.with_suffix(out_file.suffix + ".part")
    with urlopen(url) as resp, partial.open("wb") as handle:
        while chunk := resp.read(1024 * 1024):
            handle.write(chunk)
    partial.rename(out_file)


def _parse_genes(gtf: Path, chromosomes: list[str]) -> list[tuple]:
    keep = set(chromosomes)
    order = {c: i for i, c in enumerate(chromosomes)}
    rows = []
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "gene" or fields[0] not in keep:
                continue
            attributes = fields[8]
            gene_id = _GENE_ID.search(attributes).group(1).split(".")[0]
            name = _GENE_NAME.search(attributes)
            biotype = _GENE_BIOTYPE.search(attributes)
            rows.append(
                (
                    order[fields[0]],
                    int(fields[3]),
                    gene_id,
                    name.group(1) if name else "",
                    "chr" + fields[0],
                    int(fields[4]),
                    biotype.group(1) if biotype else "",
                )
            )
    rows.sort()
    return rows


def _report_immunoglobulin_rule(rows: list[tuple]) -> None:
    """Compare the name rule of analyze_cnv's immunoglobulin class with IG_* biotypes."""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    import pandas as pd

    from chatspatial.utils.gene_classes import gene_class_mask

    names = pd.Index([r[3] for r in rows])
    by_name = gene_class_mask(names, "immunoglobulin")
    by_biotype = pd.Series([r[6].startswith("IG_") for r in rows]).to_numpy()
    print(
        f"  immunoglobulin: name rule {int(by_name.sum())}, IG_* biotype "
        f"{int(by_biotype.sum())}, both {int((by_name & by_biotype).sum())}"
    )
    print(f"  name rule only: {sorted(set(names[by_name & ~by_biotype]))}")
    print(f"  biotype only:   {sorted(set(names[by_biotype & ~by_name]))}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--download-dir", type=Path, required=True)
    parser.add_argument("--species", nargs="+", default=list(SOURCES))
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for species in args.species:
        url = _url(species)
        gtf = args.download_dir / url.rsplit("/", 1)[1]
        _download(url, gtf)
        rows = _parse_genes(gtf, SOURCES[species][1])
        out = OUT_DIR / SOURCES[species][2]
        # mtime=0 keeps the gzip bytes identical across rebuilds.
        with (
            open(out, "wb") as raw,
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz,
        ):
            gz.write(b"gene_id\tgene_name\tchromosome\tstart\tend\n")
            for _, start, gene_id, name, chrom, end, _ in rows:
                gz.write(f"{gene_id}\t{name}\t{chrom}\t{start}\t{end}\n".encode())
        print(f"{species}: {len(rows)} genes -> {out} ({out.stat().st_size} bytes)")
        _report_immunoglobulin_rule(rows)


if __name__ == "__main__":
    main()
