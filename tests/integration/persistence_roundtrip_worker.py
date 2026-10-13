"""Create or round-trip a semantic H5AD fixture in a selected Python runtime."""

import json
import sys
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from chatspatial.utils.persistence import export_adata


def fixture():
    adata = ad.AnnData(sparse.csr_matrix(np.arange(16).reshape(4, 4)))
    adata.obs_names = pd.Index(["001", "002", "003", "004"], dtype="string")
    adata.var_names = pd.Index(["01", "02", "03", "04"], dtype="string")
    for table in [adata.obs, adata.var]:
        table["label"] = pd.array(["001", None, "", "β"], dtype="string")
    adata.obs["ordered"] = pd.Categorical(
        ["b", None, "a", "b"],
        categories=pd.Index(["b", "a", "unused"], dtype="string"),
        ordered=True,
    )
    adata.obs["integer"] = pd.array([1, None, 3, 4], dtype="Int64")
    adata.obs["boolean"] = pd.array([True, None, False, True], dtype="boolean")
    adata.uns["nested"] = {"table": pd.DataFrame({"label": ["001", None, "", "β"]})}
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["spatial"] = np.arange(8).reshape(4, 2)
    adata.raw = adata.copy()
    return adata[:, :3].copy()


def snapshot(adata):
    def frame(df):
        return {
            "index": list(map(str, df.index)),
            "values": df.astype(object).where(df.notna(), None).to_dict("list"),
        }

    return {
        "X": adata.X.toarray().tolist(),
        "counts": adata.layers["counts"].toarray().tolist(),
        "spatial": adata.obsm["spatial"].tolist(),
        "obs": frame(adata.obs),
        "var": frame(adata.var),
        "raw_X": adata.raw.X.toarray().tolist(),
        "raw_var": frame(adata.raw.var),
        "uns": frame(adata.uns["nested"]["table"]),
        "categories": adata.obs["ordered"].cat.categories.tolist(),
        "ordered": adata.obs["ordered"].cat.ordered,
    }


if __name__ == "__main__":
    source, destination = sys.argv[1:]
    adata = fixture() if source == "create" else ad.read_h5ad(source)
    before = snapshot(adata)
    original_setting = ad.settings.allow_write_nullable_strings
    export_adata("compatibility", adata, Path(destination))
    assert snapshot(adata) == before, "Export mutated the source"
    assert ad.settings.allow_write_nullable_strings == original_setting
    assert snapshot(ad.read_h5ad(destination)) == before, "Round-trip changed values"
    Path(destination + ".json").write_text(json.dumps(before, ensure_ascii=False))
    print(
        json.dumps(
            {
                "anndata": version("anndata"),
                "pandas": version("pandas"),
                "output": destination,
            }
        )
    )
