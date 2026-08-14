"""Unit tests for spatial_registration routing and MCP wrapper contracts."""

from __future__ import annotations

import sys
import types

import numpy as np
import pytest

from chatspatial.models.data import RegistrationParameters
from chatspatial.tools import spatial_registration as reg
from chatspatial.utils.exceptions import ParameterError, ProcessingError


class DummyCtx:
    def __init__(self, datasets: dict[str, object]):
        self.datasets = datasets
        self.warnings: list[str] = []

    async def get_adata(self, data_id: str):
        return self.datasets[data_id]

    async def set_adatas(self, updates: dict[str, object]) -> None:
        self.datasets.update(updates)

    async def warning(self, message: str) -> None:
        self.warnings.append(message)


class _IdentityImageTransform:
    def to_image(self, coords):
        return np.asarray(coords, dtype=np.float32)

    def from_image(self, coords):
        return np.asarray(coords, dtype=np.float32)


def _fake_image(_coords, _intensity, _image_size, **_kwargs):
    return [np.array([0.0]), np.array([0.0])], "img", _IdentityImageTransform()


def test_validate_spatial_coords_raises_for_missing_spatial(minimal_spatial_adata):
    adata = minimal_spatial_adata.copy()
    del adata.obsm["spatial"]
    with pytest.raises(ParameterError, match="missing spatial coordinates"):
        reg._validate_spatial_coords([adata])


def test_register_slices_requires_at_least_two_slices(minimal_spatial_adata):
    with pytest.raises(ParameterError, match="at least 2 slices"):
        reg.register_slices([minimal_spatial_adata.copy()], RegistrationParameters())


def test_register_slices_dispatches_to_paste_and_stalign(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()
    calls: list[str] = []

    fake_paste_module = object()
    fake_stalign_module = object()
    fake_torch_module = object()

    def _require(name, *_args, **_kwargs):
        assert name in {"paste", "torch"}
        return fake_paste_module if name == "paste" else fake_torch_module

    def _require_module(name, module_name, *_args, **_kwargs):
        assert (name, module_name) == ("stalign", "STalign.STalign")
        return fake_stalign_module

    def _fake_paste(adata_list, params, spatial_key="spatial", *, pst=None):
        assert pst is fake_paste_module
        calls.append(f"paste:{spatial_key}")
        return adata_list

    def _fake_stalign(
        adata_list,
        params,
        spatial_key="spatial",
        *,
        stalign_module=None,
        torch_module=None,
        expression_sources=None,
    ):
        assert stalign_module is fake_stalign_module
        assert torch_module is fake_torch_module
        calls.append(f"stalign:{spatial_key}")
        return adata_list

    monkeypatch.setattr(reg, "require", _require)
    monkeypatch.setattr(reg, "require_module", _require_module)
    monkeypatch.setattr(reg, "_register_paste", _fake_paste)
    monkeypatch.setattr(reg, "_register_stalign", _fake_stalign)

    out1 = reg.register_slices([ad1, ad2], RegistrationParameters(method="paste"))
    out2 = reg.register_slices([ad1, ad2], RegistrationParameters(method="stalign"))
    assert out1[0] is ad1
    assert out2[0] is ad1
    assert calls == ["paste:spatial", "stalign:spatial"]


@pytest.mark.asyncio
async def test_register_spatial_slices_mcp_happy_path_records_metadata(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    src = minimal_spatial_adata.copy()
    tgt = minimal_spatial_adata.copy()
    captured: list[dict[str, object]] = []

    def _fake_register_slices(adata_list, params, *, ctx=None, expression_sources=None):
        assert isinstance(ctx, DummyCtx)
        for i, adata in enumerate(adata_list):
            adata.obsm["spatial_registered"] = adata.obsm["spatial"] + i
        return adata_list

    monkeypatch.setattr(reg, "register_slices", _fake_register_slices)
    monkeypatch.setattr(
        reg,
        "store_analysis_metadata",
        lambda _adata, **kwargs: captured.append(kwargs),
    )
    monkeypatch.setattr(reg, "export_analysis_result", lambda *_args, **_kwargs: [])

    from chatspatial.models.data import RegistrationParameters

    out = await reg.register_spatial_slices_mcp(
        "src",
        "tgt",
        DummyCtx({"src": src, "tgt": tgt}),
        params=RegistrationParameters(method="paste"),
    )
    assert out.registration_completed is True
    assert out.spatial_key_registered == "spatial_registered"
    assert len(captured) == 2
    assert captured[0]["analysis_name"] == "registration_paste"
    assert captured[0]["results_keys"] == {"obsm": ["spatial_registered"]}
    assert "spatial_registered" not in src.obsm
    assert "spatial_registered" not in tgt.obsm


@pytest.mark.asyncio
async def test_register_spatial_slices_mcp_late_failure_preserves_both_sources(
    minimal_spatial_adata,
    monkeypatch: pytest.MonkeyPatch,
):
    src = minimal_spatial_adata.copy()
    tgt = minimal_spatial_adata.copy()
    ctx = DummyCtx({"src": src, "tgt": tgt})

    def _fake_register_slices(adata_list, params, *, ctx=None, expression_sources=None):
        del params, ctx
        for adata in adata_list:
            adata.obsm["spatial_registered"] = adata.obsm["spatial"] + 1
        return adata_list

    export_calls = 0

    def _fail_second_export(*_args, **_kwargs):
        nonlocal export_calls
        export_calls += 1
        if export_calls == 2:
            raise RuntimeError("target export failed")
        return []

    monkeypatch.setattr(reg, "register_slices", _fake_register_slices)
    monkeypatch.setattr(reg, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(reg, "export_analysis_result", _fail_second_export)

    with pytest.raises(ProcessingError, match="target export failed"):
        await reg.register_spatial_slices_mcp(
            "src",
            "tgt",
            ctx,
            params=RegistrationParameters(method="paste"),
        )

    assert ctx.datasets["src"] is src
    assert ctx.datasets["tgt"] is tgt
    assert "spatial_registered" not in src.obsm
    assert "spatial_registered" not in tgt.obsm


@pytest.mark.asyncio
async def test_register_spatial_slices_mcp_rejects_same_dataset(
    minimal_spatial_adata,
):
    ctx = DummyCtx({"same": minimal_spatial_adata.copy()})

    with pytest.raises(ParameterError, match="distinct datasets"):
        await reg.register_spatial_slices_mcp(
            "same",
            "same",
            ctx,
            params=RegistrationParameters(method="paste"),
        )


@pytest.mark.asyncio
async def test_register_spatial_slices_mcp_wraps_runtime_errors(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    src = minimal_spatial_adata.copy()
    tgt = minimal_spatial_adata.copy()
    monkeypatch.setattr(reg, "require", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        reg,
        "register_slices",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
    )

    with pytest.raises(ProcessingError, match="Registration failed: boom"):
        await reg.register_spatial_slices_mcp(
            "src",
            "tgt",
            DummyCtx({"src": src, "tgt": tgt}),
            params=RegistrationParameters(method="paste"),
        )


def test_get_common_genes_handles_duplicate_gene_names(minimal_spatial_adata):
    ad1 = minimal_spatial_adata[:, :4].copy()
    ad2 = minimal_spatial_adata[:, 2:6].copy()

    ad1.var_names = ["g0", "g0", "g2", "g3"]
    ad2.var_names = ["g2", "g3", "g4", "g4"]

    common = reg._get_common_genes([ad1, ad2])

    assert set(common) == {"g2", "g3"}


def test_register_slices_unknown_method_raises_parameter_error(minimal_spatial_adata):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()
    params = RegistrationParameters(method="paste").model_copy(
        update={"method": "unknown"}
    )

    with pytest.raises(ParameterError, match="Unknown method"):
        reg.register_slices([ad1, ad2], params)


def test_register_stalign_rejects_non_pairwise_input(
    minimal_spatial_adata,
):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()
    ad3 = minimal_spatial_adata.copy()

    with pytest.raises(ParameterError, match="only supports pairwise registration"):
        reg._register_stalign([ad1, ad2, ad3], RegistrationParameters(method="stalign"))


def test_register_stalign_invalid_transform_payload_preserves_domain_error(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()

    monkeypatch.setattr(
        reg,
        "_prepare_stalign_image",
        _fake_image,
    )
    monkeypatch.setattr(reg, "get_device", lambda prefer_gpu=False: "cpu")

    fake_torch = types.SimpleNamespace(
        float32="float32",
        tensor=lambda x, dtype=None, device=None: np.asarray(x),
        Tensor=np.ndarray,
    )

    fake_st = types.ModuleType("STalign.STalign")
    fake_st.LDDMM = lambda **_kwargs: {"A": None, "v": None, "xv": None}
    fake_st.transform_points_source_to_target = lambda xv, v, A, points: points

    with pytest.raises(
        ProcessingError, match="STalign did not return valid transformation"
    ):
        reg._register_stalign(
            [ad1, ad2],
            RegistrationParameters(method="stalign"),
            stalign_module=fake_st,
            torch_module=fake_torch,
        )


def test_register_stalign_expression_intensity_uses_dense_sum_path(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()
    ad1.X = np.asarray(ad1.X, dtype=float)
    ad2.X = np.asarray(ad2.X, dtype=float)

    captured: dict[str, np.ndarray] = {}

    def _fake_prepare(coords, intensity, image_size, **_kwargs):
        del coords, image_size
        key = (
            "source_intensity"
            if "source_intensity" not in captured
            else "target_intensity"
        )
        captured[key] = np.asarray(intensity)
        return (
            [np.array([0.0]), np.array([0.0])],
            np.ones((1, 2, 2), dtype=np.float32),
            _IdentityImageTransform(),
        )

    monkeypatch.setattr(reg, "_prepare_stalign_image", _fake_prepare)
    monkeypatch.setattr(reg, "get_device", lambda prefer_gpu=False: "cpu")

    fake_torch = types.SimpleNamespace(
        float32="float32",
        tensor=lambda x, dtype=None, device=None: np.asarray(x, dtype=np.float32),
        Tensor=type("Tensor", (), {}),
    )

    fake_st = types.ModuleType("STalign.STalign")
    fake_st.LDDMM = lambda **_kwargs: {
        "A": np.eye(2),
        "v": np.zeros((1, 2)),
        "xv": np.zeros((1, 2)),
    }
    fake_st.transform_points_source_to_target = lambda xv, v, A, points: points
    out = reg._register_stalign(
        [ad1, ad2],
        RegistrationParameters(method="stalign", stalign_use_expression=True),
        stalign_module=fake_st,
        torch_module=fake_torch,
    )

    assert out[0].obsm["spatial_registered"].shape == ad1.obsm["spatial"].shape
    assert captured["source_intensity"].shape[0] == ad1.n_obs
    assert np.allclose(captured["source_intensity"], ad1.X.sum(axis=1))
    assert np.allclose(captured["target_intensity"], ad2.X.sum(axis=1))


def test_prepare_stalign_image_returns_normalized_tensor(
    monkeypatch: pytest.MonkeyPatch,
):
    fake_torch = types.SimpleNamespace(
        float32="float32",
        arange=lambda stop, dtype=None: np.arange(stop, dtype=np.float32),
        tensor=lambda x, dtype=None: np.asarray(x, dtype=np.float32),
    )

    coords = np.array([[0.0, 0.0], [10.0, 0.0], [5.0, 8.0]], dtype=np.float32)
    intensity = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    coords_before = coords.copy()
    intensity_before = intensity.copy()
    xgrid, image, transform = reg._prepare_stalign_image(
        coords,
        intensity,
        (16, 12),
        torch_module=fake_torch,
    )

    assert len(xgrid) == 2
    assert image.shape == (1, 16, 12)
    assert float(image.max()) <= 1.0
    assert float(image.min()) >= 0.0
    np.testing.assert_allclose(coords, coords_before)
    np.testing.assert_allclose(intensity, intensity_before)
    np.testing.assert_allclose(transform.from_image(transform.to_image(coords)), coords)


def test_register_paste_pairwise_populates_registered_coords(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()

    def _pairwise_align(_slice1, _slice2, **_kwargs):
        return np.eye(_slice1.n_obs)

    def _stack_pairwise(slices, _pis):
        out0 = slices[0].copy()
        out1 = slices[1].copy()
        out0.obsm["spatial"] = out0.obsm["spatial"] + 1.0
        out1.obsm["spatial"] = out1.obsm["spatial"] + 2.0
        return [out0, out1]

    fake_paste = types.ModuleType("paste")
    fake_paste.pairwise_align = _pairwise_align
    fake_paste.stack_slices_pairwise = _stack_pairwise

    fake_scanpy = types.ModuleType("scanpy")
    fake_scanpy.pp = types.SimpleNamespace(
        normalize_total=lambda *_args, **_kwargs: None,
        log1p=lambda *_args, **_kwargs: None,
    )

    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setattr(reg, "get_ot_backend", lambda _use_gpu: "numpy")

    out = reg._register_paste(
        [ad1, ad2],
        RegistrationParameters(method="paste"),
        pst=fake_paste,
    )
    np.testing.assert_allclose(
        out[0].obsm["spatial_registered"], ad1.obsm["spatial"] + 1.0
    )
    np.testing.assert_allclose(
        out[1].obsm["spatial_registered"], ad2.obsm["spatial"] + 2.0
    )
    assert out[0].X is ad1.X
    assert out[1].X is ad2.X


def test_register_paste_multi_slice_uses_center_alignment(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata[:40].copy()
    ad3 = minimal_spatial_adata[:30].copy()

    pairwise_calls: list[dict[str, object]] = []

    def _pairwise_align(ref_slice, moving_slice, **kwargs):
        pairwise_calls.append(
            {"ref": ref_slice.n_obs, "moving": moving_slice.n_obs, **kwargs}
        )
        return np.full(
            (ref_slice.n_obs, moving_slice.n_obs),
            1.0 / (ref_slice.n_obs * moving_slice.n_obs),
        )

    def _center_align(_ref, slices, **_kwargs):
        plans = [
            np.full(
                (_ref.n_obs, slice_data.n_obs),
                1.0 / (_ref.n_obs * slice_data.n_obs),
            )
            for slice_data in slices
        ]
        return _ref, plans

    def _stack_center(center, slices, plans):
        assert all(
            plan.shape == (center.n_obs, slice_data.n_obs)
            for plan, slice_data in zip(plans, slices, strict=True)
        )
        aligned = [slice_data.copy() for slice_data in slices]
        for index, slice_data in enumerate(aligned):
            slice_data.obsm["spatial"] = slice_data.obsm["spatial"] + index
        return center, aligned

    fake_paste = types.ModuleType("paste")
    fake_paste.pairwise_align = _pairwise_align
    fake_paste.center_align = _center_align
    fake_paste.stack_slices_center = _stack_center

    fake_scanpy = types.ModuleType("scanpy")
    fake_scanpy.pp = types.SimpleNamespace(
        normalize_total=lambda *_args, **_kwargs: None,
        log1p=lambda *_args, **_kwargs: None,
    )

    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    monkeypatch.setattr(reg, "get_ot_backend", lambda _use_gpu: "numpy")

    params = RegistrationParameters(method="paste", reference_idx=0, use_gpu=False)
    out = reg._register_paste([ad1, ad2, ad3], params, pst=fake_paste)

    assert len(pairwise_calls) == 2
    assert pairwise_calls[0]["backend"] == "numpy"
    np.testing.assert_allclose(out[0].obsm["spatial_registered"], ad1.obsm["spatial"])
    np.testing.assert_allclose(
        out[1].obsm["spatial_registered"], ad2.obsm["spatial"] + 1
    )
    np.testing.assert_allclose(
        out[2].obsm["spatial_registered"], ad3.obsm["spatial"] + 2
    )
    assert out[1].obsm["spatial_registered"].shape == ad2.obsm["spatial"].shape
    assert out[2].obsm["spatial_registered"].shape == ad3.obsm["spatial"].shape
    assert out[0].X is ad1.X
    assert out[1].X is ad2.X
    assert out[2].X is ad3.X


def test_register_paste_rejects_missing_center_transport_plan(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    slices = [minimal_spatial_adata.copy() for _ in range(3)]

    fake_paste = types.ModuleType("paste")
    fake_paste.pairwise_align = lambda ref, moving, **_kwargs: np.eye(
        moving.n_obs, ref.n_obs
    )
    fake_paste.center_align = lambda ref, _slices, **_kwargs: (
        ref,
        [np.eye(ref.n_obs), np.eye(ref.n_obs)],
    )
    monkeypatch.setattr(reg, "get_ot_backend", lambda _use_gpu: "numpy")

    with pytest.raises(ProcessingError, match="expected 3, received 2"):
        reg._register_paste(
            slices,
            RegistrationParameters(method="paste", reference_idx=0, use_gpu=False),
            pst=fake_paste,
        )


def test_register_stalign_success_with_uniform_intensity(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()

    class FakeTensor(np.ndarray):
        def numpy(self):
            return np.asarray(self)

    def _to_tensor(x, dtype=None, device=None):
        arr = np.asarray(x, dtype=np.float32)
        return arr.view(FakeTensor)

    fake_torch = types.SimpleNamespace(
        float32="float32",
        tensor=_to_tensor,
        Tensor=FakeTensor,
    )

    fake_st = types.ModuleType("STalign.STalign")
    fake_st.LDDMM = lambda **_kwargs: {"A": "A", "v": "v", "xv": "xv"}
    fake_st.transform_points_source_to_target = lambda _xv, _v, _A, points: _to_tensor(
        np.asarray(points) + 3.0
    )

    monkeypatch.setattr(reg, "_prepare_stalign_image", _fake_image)
    monkeypatch.setattr(reg, "get_device", lambda prefer_gpu=False: "cpu")

    params = RegistrationParameters(method="stalign", stalign_use_expression=False)
    out = reg._register_stalign(
        [ad1, ad2],
        params,
        stalign_module=fake_st,
        torch_module=fake_torch,
    )

    np.testing.assert_allclose(
        out[0].obsm["spatial_registered"], ad1.obsm["spatial"] + 3.0
    )
    np.testing.assert_allclose(out[1].obsm["spatial_registered"], ad2.obsm["spatial"])
    assert out[0].X is ad1.X
    assert out[1].X is ad2.X


def test_register_slices_defaults_to_paste_when_params_none(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()

    called = {"method": None}

    fake_paste_module = object()

    def _fake_paste(adata_list, params, spatial_key="spatial", *, pst=None):
        assert pst is fake_paste_module
        called["method"] = params.method
        return adata_list

    monkeypatch.setattr(reg, "require", lambda *_args, **_kwargs: fake_paste_module)
    monkeypatch.setattr(reg, "_register_paste", _fake_paste)
    out = reg.register_slices([ad1, ad2], params=None)

    assert called["method"] == "paste"
    assert len(out) == 2


@pytest.mark.asyncio
async def test_register_spatial_slices_mcp_passes_context_to_dependency_owner(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    src = minimal_spatial_adata.copy()
    tgt = minimal_spatial_adata.copy()
    seen: dict[str, object] = {}

    def _register(adata_list, _params, *, ctx=None, expression_sources=None):
        seen["ctx"] = ctx
        for adata in adata_list:
            adata.obsm["spatial_registered"] = adata.obsm["spatial"].copy()
        return adata_list

    monkeypatch.setattr(
        reg,
        "register_slices",
        _register,
    )
    monkeypatch.setattr(reg, "store_analysis_metadata", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(reg, "export_analysis_result", lambda *_args, **_kwargs: None)

    ctx = DummyCtx({"src": src, "tgt": tgt})
    out = await reg.register_spatial_slices_mcp(
        "src",
        "tgt",
        ctx,
        params=RegistrationParameters(method="stalign"),
    )

    assert out.method == "stalign"
    assert seen["ctx"] is ctx
    assert ctx.datasets["src"] is not src
    assert ctx.datasets["tgt"] is not tgt


# =============================================================================
# Issue 1 regression: normalization consistency across pairwise / multi-slice
# =============================================================================


def test_prepare_paste_slices_normalizes_uniformly(minimal_spatial_adata):
    """_prepare_paste_slices applies normalize_total + log1p to all slices."""
    import scipy.sparse as sp

    a = minimal_spatial_adata.copy()
    b = minimal_spatial_adata.copy()
    genes = list(a.var_names[:5])

    slices = reg._prepare_paste_slices([a, b], genes, "spatial")

    for s in slices:
        x = s.X.toarray() if sp.issparse(s.X) else s.X
        assert np.issubdtype(x.dtype, np.floating)
        assert s.n_vars == len(genes)
        # After log1p the max should be well below raw integer range
        assert x.max() < 15


def test_prepare_paste_slices_copies_alt_key_to_spatial(minimal_spatial_adata):
    """When spatial_key != 'spatial', obsm['spatial'] is created for PASTE."""
    a = minimal_spatial_adata.copy()
    coords = a.obsm.pop("spatial")
    a.obsm["X_spatial"] = coords
    genes = list(a.var_names[:3])

    slices = reg._prepare_paste_slices([a], genes, "X_spatial")
    s = slices[0]

    assert "spatial" in s.obsm
    assert "X_spatial" in s.obsm
    np.testing.assert_array_equal(s.obsm["spatial"], s.obsm["X_spatial"])


def test_register_paste_pairwise_normalizes_before_align(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """Pairwise PASTE branch must pass normalized data to pairwise_align."""
    import scipy.sparse as sp

    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()

    captured = []

    def _pairwise_align(s1, s2, **_kw):
        # Capture expression data that PASTE receives
        x1 = s1.X.toarray() if sp.issparse(s1.X) else s1.X
        captured.append(x1.max())
        return np.eye(s1.n_obs, s2.n_obs) / s2.n_obs

    def _stack(slices, pis):
        return slices

    fake_paste = types.ModuleType("paste")
    fake_paste.pairwise_align = _pairwise_align
    fake_paste.stack_slices_pairwise = _stack
    monkeypatch.setattr(reg, "get_ot_backend", lambda _use_gpu: "numpy")

    reg._register_paste(
        [ad1, ad2],
        RegistrationParameters(method="paste"),
        pst=fake_paste,
    )

    assert len(captured) == 1
    # After normalize_total + log1p, max should be << raw counts
    assert captured[0] < 15, "Pairwise branch should normalize before PASTE"


def test_register_paste_multi_slice_normalizes_before_align(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """Multi-slice PASTE branch must also pass normalized data."""
    import scipy.sparse as sp

    slices = [minimal_spatial_adata.copy() for _ in range(3)]

    captured = []

    def _pairwise_align(s1, s2, **_kw):
        x1 = s1.X.toarray() if sp.issparse(s1.X) else s1.X
        captured.append(x1.max())
        return np.eye(s1.n_obs, s2.n_obs) / s2.n_obs

    def _center_align(_ref, all_slices, **_kw):
        return _ref, [np.eye(s.n_obs) for s in all_slices]

    def _stack_center(center, all_slices, _plans):
        return center, all_slices

    fake_paste = types.ModuleType("paste")
    fake_paste.pairwise_align = _pairwise_align
    fake_paste.center_align = _center_align
    fake_paste.stack_slices_center = _stack_center
    monkeypatch.setattr(reg, "get_ot_backend", lambda _use_gpu: "numpy")

    reg._register_paste(
        slices,
        RegistrationParameters(method="paste", use_gpu=False),
        pst=fake_paste,
    )

    assert len(captured) >= 1
    for val in captured:
        assert val < 15, "Multi-slice branch should normalize before PASTE"


# =============================================================================
# Issue 2 regression: pairwise branch must respect alternative spatial_key
# =============================================================================


def test_register_paste_pairwise_works_with_alt_spatial_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """Pairwise PASTE must succeed when coords are in X_spatial, not spatial."""
    ad1 = minimal_spatial_adata.copy()
    ad2 = minimal_spatial_adata.copy()
    # Move coords to alternative key
    for ad in (ad1, ad2):
        ad.obsm["X_spatial"] = ad.obsm.pop("spatial")

    def _pairwise_align(s1, s2, **_kw):
        assert "spatial" in s1.obsm, "PASTE needs obsm['spatial']"
        assert "spatial" in s2.obsm, "PASTE needs obsm['spatial']"
        return np.eye(s1.n_obs, s2.n_obs) / s2.n_obs

    def _stack(slices, pis):
        return slices

    fake_paste = types.ModuleType("paste")
    fake_paste.pairwise_align = _pairwise_align
    fake_paste.stack_slices_pairwise = _stack
    monkeypatch.setattr(reg, "get_ot_backend", lambda _use_gpu: "numpy")

    result = reg._register_paste(
        [ad1, ad2],
        RegistrationParameters(method="paste"),
        spatial_key="X_spatial",
        pst=fake_paste,
    )

    assert "spatial_registered" in result[0].obsm
    assert "spatial_registered" in result[1].obsm


def test_register_paste_multi_slice_works_with_alt_spatial_key(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """Multi-slice PASTE must succeed when coords are in X_spatial."""
    slices = [minimal_spatial_adata.copy() for _ in range(3)]
    for ad in slices:
        ad.obsm["X_spatial"] = ad.obsm.pop("spatial")

    def _pairwise_align(s1, s2, **_kw):
        assert "spatial" in s1.obsm
        assert "spatial" in s2.obsm
        return np.eye(s1.n_obs, s2.n_obs) / s2.n_obs

    def _center_align(_ref, all_slices, **_kw):
        return _ref, [np.eye(s.n_obs) for s in all_slices]

    def _stack_center(center, all_slices, _plans):
        return center, all_slices

    fake_paste = types.ModuleType("paste")
    fake_paste.pairwise_align = _pairwise_align
    fake_paste.center_align = _center_align
    fake_paste.stack_slices_center = _stack_center
    monkeypatch.setattr(reg, "get_ot_backend", lambda _use_gpu: "numpy")

    result = reg._register_paste(
        slices,
        RegistrationParameters(method="paste", use_gpu=False),
        spatial_key="X_spatial",
        pst=fake_paste,
    )

    for r in result:
        assert "spatial_registered" in r.obsm


# ---------------------------------------------------------------------------
# Regression: registration metadata must include method-specific parameters
# ---------------------------------------------------------------------------


@pytest.mark.asyncio
async def test_register_mcp_metadata_contains_paste_params(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """store_analysis_metadata must record PASTE-specific params for provenance."""
    source = minimal_spatial_adata.copy()
    target = minimal_spatial_adata.copy()
    ctx = DummyCtx({"src": source, "tgt": target})

    def _fake_register(slices, params, *, ctx=None, expression_sources=None):
        assert ctx is not None
        for s in slices:
            s.obsm["spatial_registered"] = s.obsm["spatial"]
        return slices

    monkeypatch.setattr(reg, "register_slices", _fake_register)

    captured_params: list[dict] = []

    def _capture_store(_adata, **kwargs):
        captured_params.append(kwargs.get("parameters", {}))

    monkeypatch.setattr(reg, "store_analysis_metadata", _capture_store)
    monkeypatch.setattr(reg, "export_analysis_result", lambda *_a, **_k: None)

    params = RegistrationParameters(
        method="paste",
        paste_alpha=0.3,
        paste_n_components=20,
        paste_numItermax=100,
        use_gpu=False,
    )
    result = await reg.register_spatial_slices_mcp("src", "tgt", ctx, params)
    assert result.registration_completed is True

    # Both source and target metadata calls must have method-specific params
    assert len(captured_params) == 2
    for p in captured_params:
        assert p["paste_alpha"] == 0.3
        assert p["paste_n_components"] == 20
        assert p["paste_numItermax"] == 100
        assert p["use_gpu"] is False
        assert "stalign_niter" not in p  # no cross-method leakage


@pytest.mark.asyncio
async def test_register_mcp_metadata_contains_stalign_params(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """store_analysis_metadata must record STalign-specific params."""
    source = minimal_spatial_adata.copy()
    target = minimal_spatial_adata.copy()
    ctx = DummyCtx({"src": source, "tgt": target})

    def _fake_register(slices, params, *, ctx=None, expression_sources=None):
        assert ctx is not None
        for s in slices:
            s.obsm["spatial_registered"] = s.obsm["spatial"]
        return slices

    monkeypatch.setattr(reg, "register_slices", _fake_register)

    captured_params: list[dict] = []

    def _capture_store(_adata, **kwargs):
        captured_params.append(kwargs.get("parameters", {}))

    monkeypatch.setattr(reg, "store_analysis_metadata", _capture_store)
    monkeypatch.setattr(reg, "export_analysis_result", lambda *_a, **_k: None)

    params = RegistrationParameters(
        method="stalign",
        stalign_niter=80,
        stalign_a=300.0,
        stalign_use_expression=False,
    )
    result = await reg.register_spatial_slices_mcp("src", "tgt", ctx, params)
    assert result.registration_completed is True

    assert len(captured_params) == 2
    for p in captured_params:
        assert p["stalign_niter"] == 80
        assert p["stalign_a"] == 300.0
        assert p["stalign_use_expression"] is False
        assert "paste_alpha" not in p


@pytest.mark.asyncio
async def test_register_mcp_recovers_expression_for_stalign_from_raw(
    minimal_spatial_adata, monkeypatch: pytest.MonkeyPatch
):
    """STalign rasterizes summed expression, which cannot be negative.

    Variance-stabilizing normalization is ChatSpatial's default, so reading X
    would make expression-based registration unusable right after preprocessing.
    """
    source = minimal_spatial_adata.copy()
    target = minimal_spatial_adata.copy()
    for adata in (source, target):
        adata.raw = adata.copy()
        adata.X = np.full((adata.n_obs, adata.n_vars), -1.5, dtype=np.float32)

    seen: dict[str, object] = {}

    def _fake_register(slices, params, *, ctx=None, expression_sources=None):
        seen["expression_sources"] = expression_sources
        for adata in slices:
            adata.obsm["spatial_registered"] = np.asarray(
                adata.obsm["spatial"], dtype=np.float32
            )
        return slices

    monkeypatch.setattr(reg, "register_slices", _fake_register)
    monkeypatch.setattr(reg, "store_analysis_metadata", lambda *_a, **_k: None)
    monkeypatch.setattr(reg, "export_analysis_result", lambda *_a, **_k: [])

    ctx = DummyCtx({"s": source, "t": target})
    await reg.register_spatial_slices_mcp(
        "s", "t", ctx, RegistrationParameters(method="stalign")
    )

    sources = seen["expression_sources"]
    assert sources is not None
    assert all(float(np.asarray(a.X).min()) >= 0 for a in sources)
    assert any("adata.raw" in warning for warning in ctx.warnings)
