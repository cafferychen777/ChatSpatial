"""Shared test helpers to keep behavior checks consistent and DRY."""

from __future__ import annotations

from pathlib import Path

from chatspatial.server import load_data


def extract_saved_path(viz_message: str) -> Path:
    """Parse saved file path from visualize_data return text."""
    prefix = "Visualization saved: "
    for line in viz_message.splitlines():
        if line.startswith(prefix):
            return Path(line[len(prefix) :].strip())
    raise AssertionError(f"Cannot parse visualization output path from: {viz_message}")


async def load_generic_dataset(spatial_dataset_path, *, name: str):
    """Load the shared temporary h5ad fixture via public load_data API."""
    return await load_data(str(spatial_dataset_path), "generic", name=name)


# =============================================================================
# DeepSpot-M test doubles
# =============================================================================
#
# The virtual spatial transcriptomics tool must be exercised without installing
# deepspotm, importing torch, or reaching Hugging Face. These stand-ins cover
# exactly the API surface the integration touches and nothing more.


class FakeTensor:
    """The narrow slice of the tensor API the DeepSpot-M integration uses."""

    def __init__(self, array):
        self.array = array
        self.device = None

    def to(self, device):
        self.device = device
        return self

    def float(self):
        return self

    def cpu(self):
        return self

    def numpy(self):
        return self.array


def fake_torch_module():
    """Build a stand-in ``torch`` module for tile batching and inference."""
    import contextlib
    import types

    import numpy as np

    module = types.ModuleType("torch")
    module.device = lambda spec: f"device({spec})"
    module.stack = lambda items: FakeTensor(np.stack(items))
    module.no_grad = contextlib.nullcontext
    return module


class FakeDeepSpotMModel:
    """Predicts a deterministic value derived from each tile's mean pixel."""

    gene_names = ("EPCAM", "CD3D", "PTPRC")

    def __init__(self):
        self.batches: list[int] = []
        self.requested_genes: list[list[str]] = []

    def predict_genes(self, pixel_values, genes):
        import numpy as np

        unknown = [symbol for symbol in genes if symbol not in self.gene_names]
        if unknown:
            raise KeyError(", ".join(unknown))

        n_tiles = pixel_values.array.shape[0]
        self.batches.append(n_tiles)
        self.requested_genes.append(list(genes))
        means = pixel_values.array.reshape(n_tiles, -1).mean(axis=1)
        return FakeTensor(
            np.outer(means, np.arange(1, len(genes) + 1)).astype(np.float32)
        )


def fake_deepspotm_module(recorder: dict):
    """Build a stand-in ``deepspotm`` module recording its load arguments."""
    import types

    import numpy as np

    module = types.ModuleType("deepspotm")

    class DeepSpotM:
        @classmethod
        def from_pretrained(
            cls, repo_id_or_path, *, source=None, device=None, revision=None
        ):
            recorder["repo"] = repo_id_or_path
            recorder["source"] = source
            recorder["device"] = device
            recorder["revision"] = revision
            model = FakeDeepSpotMModel()
            recorder["model"] = model

            def image_processor(tile):
                assert tile.mode == "RGB"
                assert tile.size == (224, 224)
                return np.asarray(tile, dtype=np.float32).transpose(2, 0, 1)

            return model, image_processor

    module.DeepSpotM = DeepSpotM
    return module
