"""Contracts for the optional rctd-py deconvolution adapter."""

from __future__ import annotations

import io
from types import SimpleNamespace

import numpy as np
import pytest

import chatspatial.tools.deconvolution as deconv_module
from chatspatial.tools.deconvolution import rctd as rctd_router
from chatspatial.tools.deconvolution import rctd_python
from chatspatial.tools.deconvolution.base import PreparedDeconvolutionData
from chatspatial.utils.exceptions import (
    DependencyError,
    ParameterError,
    ProcessingError,
)


class _Ctx:
    async def warning(self, _message: str) -> None:
        return None


def _prepared_data(minimal_spatial_adata) -> PreparedDeconvolutionData:
    spatial = minimal_spatial_adata.copy()
    reference = minimal_spatial_adata.copy()
    reference.obs["cell_type"] = ["A"] * (reference.n_obs // 2) + ["B"] * (
        reference.n_obs - reference.n_obs // 2
    )
    return PreparedDeconvolutionData(
        spatial=spatial,
        reference=reference,
        cell_type_key="cell_type",
        cell_types=["A", "B"],
        common_genes=list(spatial.var_names),
        spatial_coords=np.asarray(spatial.obsm["spatial"]),
        ctx=_Ctx(),
    )


def _fake_module(result, captured: dict[str, object]) -> SimpleNamespace:
    class _Reference:
        def __init__(self, adata, **kwargs) -> None:
            captured["reference"] = adata
            captured["reference_kwargs"] = kwargs
            self.cell_type_names = ["A", "B"]

    def _config(**kwargs):
        captured["config"] = kwargs
        return kwargs

    def _run(spatial, reference, **kwargs):
        captured["spatial"] = spatial
        captured["reference_object"] = reference
        captured["run_kwargs"] = kwargs
        return result

    return SimpleNamespace(
        Reference=_Reference,
        RCTDConfig=_config,
        run_rctd=_run,
        SPOT_CLASS_NAMES=[
            "reject",
            "singlet",
            "doublet_certain",
            "doublet_uncertain",
        ],
        __version__="0.3.7-test",
    )


def test_rctd_python_full_expands_filtered_spots_and_maps_config(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    data = _prepared_data(minimal_spatial_adata)
    mask = np.ones(data.n_spots, dtype=bool)
    mask[1] = False
    result = SimpleNamespace(
        weights=np.tile([0.75, 0.25], (int(mask.sum()), 1)),
        cell_type_names=["A", "B"],
        converged=np.ones(int(mask.sum()), dtype=bool),
        pixel_mask=mask,
    )
    captured: dict[str, object] = {}

    def _require(name: str, **kwargs):
        captured["dependency"] = (name, kwargs)
        return _fake_module(result, captured)

    monkeypatch.setattr(rctd_python, "require", _require)
    monkeypatch.setattr(rctd_python, "_ensure_likelihood_cache", lambda _module: None)
    monkeypatch.setattr(
        rctd_python,
        "version",
        lambda _name: (_ for _ in ()).throw(rctd_python.PackageNotFoundError),
    )
    proportions, stats = rctd_python.deconvolve(
        data,
        mode="full",
        device="cpu",
        batch_size=1024,
        dtype="float32",
        sigma_override=80,
    )

    assert proportions.shape == (data.n_spots, 2)
    assert proportions.iloc[1].tolist() == [0.0, 0.0]
    assert captured["dependency"] == (
        "rctd-py",
        {"feature": "RCTD Python deconvolution"},
    )
    assert captured["reference_kwargs"] == {
        "cell_type_col": "cell_type",
        "cell_min": 25,
        "min_UMI": 5,
    }
    assert captured["config"] == {
        "UMI_min_sigma": 10,
        "MAX_MULTI_TYPES": 4,
        "CONFIDENCE_THRESHOLD": 10.0,
        "DOUBLET_THRESHOLD": 25.0,
        "dtype": "float32",
        "device": "cpu",
    }
    assert captured["run_kwargs"] == {
        "mode": "full",
        "config": captured["config"],
        "batch_size": 1024,
        "sigma_override": 80,
    }
    assert stats["backend"] == "python"
    assert stats["backend_version"] == "0.3.7-test"
    assert stats["n_filtered_spots"] == 1
    outputs = stats["_backend_outputs"]
    assert outputs["obs"]["rctd_status"][1] == "filtered"
    assert outputs["obs"]["rctd_converged"][1] == np.False_


def test_rctd_python_preserves_doublet_annotations(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    data = _prepared_data(minimal_spatial_adata)
    mask = np.ones(data.n_spots, dtype=bool)
    mask[-1] = False
    n_analyzed = int(mask.sum())
    result = SimpleNamespace(
        weights=np.tile([0.6, 0.4], (n_analyzed, 1)),
        weights_doublet=np.tile([0.6, 0.4], (n_analyzed, 1)),
        spot_class=np.full(n_analyzed, 2),
        first_type=np.zeros(n_analyzed, dtype=int),
        second_type=np.ones(n_analyzed, dtype=int),
        first_class=np.ones(n_analyzed, dtype=bool),
        second_class=np.zeros(n_analyzed, dtype=bool),
        min_score=np.full(n_analyzed, 3.0),
        singlet_score=np.full(n_analyzed, 7.0),
        first_class_name=None,
        second_class_name=None,
        cell_type_names=["A", "B"],
        pixel_mask=mask,
    )
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        rctd_python,
        "require",
        lambda *_args, **_kwargs: _fake_module(result, captured),
    )
    monkeypatch.setattr(rctd_python, "_ensure_likelihood_cache", lambda _module: None)

    _proportions, stats = rctd_python.deconvolve(data, mode="doublet", device="cpu")
    outputs = stats["_backend_outputs"]

    assert outputs["obs"]["rctd_spot_class"][0] == "doublet_certain"
    assert outputs["obs"]["rctd_first_type"][0] == "A"
    assert outputs["obs"]["rctd_second_type"][0] == "B"
    assert outputs["obs"]["rctd_spot_class"][-1] == "filtered"
    assert outputs["obsm"]["rctd_doublet_weights"][-1].tolist() == [0.0, 0.0]


def test_rctd_python_rejects_invalid_weight_shape(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    data = _prepared_data(minimal_spatial_adata)
    result = SimpleNamespace(
        weights=np.ones((data.n_spots - 1, 2)),
        cell_type_names=["A", "B"],
        converged=np.ones(data.n_spots - 1, dtype=bool),
        pixel_mask=None,
    )
    monkeypatch.setattr(
        rctd_python,
        "require",
        lambda *_args, **_kwargs: _fake_module(result, {}),
    )
    monkeypatch.setattr(rctd_python, "_ensure_likelihood_cache", lambda _module: None)

    with pytest.raises(ProcessingError, match="invalid weights matrix"):
        rctd_python.deconvolve(data, device="cpu")


def test_rctd_python_preserves_dependency_errors(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    data = _prepared_data(minimal_spatial_adata)
    monkeypatch.setattr(
        rctd_python,
        "require",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("rctd-py unavailable")
        ),
    )

    with pytest.raises(DependencyError, match="rctd-py unavailable"):
        rctd_python.deconvolve(data, device="cpu")


def test_likelihood_cache_download_uses_certifi_and_atomic_file(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    archive = io.BytesIO()
    np.savez(archive, Q_100=np.ones((2, 2)), X_vals=np.ones(2))
    calls: dict[str, object] = {}

    class _Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            self.close()
            return False

    def _urlopen(url, *, timeout, context):
        calls["url"] = url
        calls["timeout"] = timeout
        calls["context"] = context
        return _Response(archive.getvalue())

    tls_context = object()
    monkeypatch.setattr(rctd_python.certifi, "where", lambda: "/tmp/test-ca.pem")

    def _create_default_context(*, cafile):
        calls["cafile"] = cafile
        return tls_context

    monkeypatch.setattr(
        rctd_python.ssl, "create_default_context", _create_default_context
    )
    monkeypatch.setattr(rctd_python.urllib.request, "urlopen", _urlopen)
    monkeypatch.setattr(
        rctd_python.importlib,
        "import_module",
        lambda name: SimpleNamespace(
            _Q_MATRICES_URL="https://example.test/q_matrices.npz"
        ),
    )
    target = tmp_path / "cache" / "q_matrices.npz"

    rctd_python._ensure_likelihood_cache(
        SimpleNamespace(__name__="rctd"),
        cache_path=target,
    )

    assert target.read_bytes() == archive.getvalue()
    assert calls["url"] == "https://example.test/q_matrices.npz"
    assert calls["timeout"] == 60
    assert calls["cafile"] == "/tmp/test-ca.pem"
    assert not list(target.parent.glob("*.part"))


def test_rctd_router_rejects_unknown_backend_before_loading_dependencies(
    minimal_spatial_adata,
) -> None:
    data = _prepared_data(minimal_spatial_adata)

    with pytest.raises(ParameterError, match="Invalid RCTD backend"):
        rctd_router.deconvolve(data, backend="unsupported")


def test_store_and_clear_rctd_backend_outputs(minimal_spatial_adata) -> None:
    adata = minimal_spatial_adata.copy()
    n_obs = adata.n_obs
    outputs = {
        "obs": {
            "rctd_status": np.array(["analyzed"] * n_obs, dtype=object),
            "rctd_converged": np.ones(n_obs, dtype=bool),
        },
        "obsm": {"rctd_doublet_weights": np.ones((n_obs, 2))},
    }
    stats = {
        "backend": "python",
        "backend_package": "rctd-py",
        "backend_version": "0.3.7",
        "mode": "full",
        "sigma_override": None,
    }

    stored = deconv_module._store_rctd_backend_outputs(adata, outputs, stats)

    assert stored == {
        "obs": ["rctd_status", "rctd_converged"],
        "obsm": ["rctd_doublet_weights"],
        "uns": ["rctd_backend"],
    }
    assert adata.uns["rctd_backend"] == {
        "backend": "python",
        "backend_package": "rctd-py",
        "backend_version": "0.3.7",
        "mode": "full",
    }

    deconv_module._clear_rctd_backend_outputs(adata)

    assert "rctd_status" not in adata.obs
    assert "rctd_converged" not in adata.obs
    assert "rctd_doublet_weights" not in adata.obsm
    assert "rctd_backend" not in adata.uns
