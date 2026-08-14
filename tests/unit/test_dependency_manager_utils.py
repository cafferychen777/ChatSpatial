"""Unit contracts for dependency manager behavior and error messaging."""

from __future__ import annotations

import sys
from contextlib import contextmanager
from types import ModuleType, SimpleNamespace

import pytest

from chatspatial.utils import dependency_manager as dm
from chatspatial.utils.exceptions import DependencyError


@pytest.fixture(autouse=True)
def _clear_dependency_caches():
    """Keep fake dependency modules isolated between tests."""
    dm._try_import.cache_clear()
    dm._check_spec.cache_clear()
    yield
    dm._try_import.cache_clear()
    dm._check_spec.cache_clear()


def _install_fake_rpy2(
    monkeypatch: pytest.MonkeyPatch,
    *,
    missing_packages: set[str] | None = None,
    unloadable_packages: set[str] | None = None,
) -> None:
    """Install lightweight fake rpy2/anndata2ri modules into sys.modules."""

    missing_packages = missing_packages or set()
    unloadable_packages = unloadable_packages or set()

    class _PackageNotInstalledError(ImportError):
        pass

    class _LibraryError(ImportError):
        pass

    class _Lock:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    class _LocalConverter:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return None

        def __exit__(self, exc_type, exc, tb):
            return False

    @contextmanager
    def _local_r_context(*, use_rlock: bool = True):
        del use_rlock
        local_environment: dict[str, object] = {}
        fake_ro.local_context_calls += 1
        try:
            yield local_environment
        finally:
            fake_ro.last_local_environment = local_environment
            fake_ro.local_context_exits += 1

    fake_anndata2ri = ModuleType("anndata2ri")
    fake_ro = ModuleType("rpy2.robjects")
    fake_conversion = ModuleType("rpy2.robjects.conversion")
    fake_packages = ModuleType("rpy2.robjects.packages")
    fake_rinterface_lib = ModuleType("rpy2.rinterface_lib")
    fake_numpy2ri = ModuleType("rpy2.robjects.numpy2ri")
    fake_pandas2ri = ModuleType("rpy2.robjects.pandas2ri")
    fake_scvi = ModuleType("scvi")
    fake_scvi_model = ModuleType("scvi.model")
    fake_scvi_external = ModuleType("scvi.external")

    def _fake_importr(name: str):
        if name in missing_packages:
            raise _PackageNotInstalledError(f"missing {name}")
        if name in unloadable_packages:
            raise _LibraryError(f"cannot load {name}")
        return object()

    fake_conversion.localconverter = _LocalConverter
    fake_ro.conversion = fake_conversion
    fake_ro.default_converter = object()
    fake_ro.globalenv = {}
    fake_ro.local_context = _local_r_context
    fake_ro.local_context_calls = 0
    fake_ro.local_context_exits = 0
    fake_ro.last_local_environment = None
    fake_numpy2ri.converter = object()
    fake_pandas2ri.converter = object()
    fake_ro.numpy2ri = fake_numpy2ri
    fake_ro.pandas2ri = fake_pandas2ri
    fake_ro.r = lambda _expr: "ok"
    fake_packages.importr = _fake_importr
    fake_packages.PackageNotInstalledError = _PackageNotInstalledError
    fake_packages.LibraryError = _LibraryError
    fake_rinterface_lib.openrlib = SimpleNamespace(rlock=_Lock())
    fake_scvi.model = fake_scvi_model
    fake_scvi.external = fake_scvi_external

    monkeypatch.setitem(sys.modules, "anndata2ri", fake_anndata2ri)
    monkeypatch.setitem(sys.modules, "rpy2", ModuleType("rpy2"))
    monkeypatch.setitem(sys.modules, "rpy2.robjects", fake_ro)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.conversion", fake_conversion)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.packages", fake_packages)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.numpy2ri", fake_numpy2ri)
    monkeypatch.setitem(sys.modules, "rpy2.robjects.pandas2ri", fake_pandas2ri)
    monkeypatch.setitem(sys.modules, "rpy2.rinterface_lib", fake_rinterface_lib)
    monkeypatch.setitem(sys.modules, "scvi", fake_scvi)
    monkeypatch.setitem(sys.modules, "scvi.model", fake_scvi_model)
    monkeypatch.setitem(sys.modules, "scvi.external", fake_scvi_external)


def test_get_info_supports_registered_alias_and_unknown_defaults():
    spatialde_by_name = dm._get_info("spatialde")
    spatialde_by_module = dm._get_info("NaiveDE")
    velovi = dm._get_info("velovi")
    unknown = dm._get_info("some-new-package")

    assert spatialde_by_name.module_name == "SpatialDE"
    assert spatialde_by_module.install_cmd == "pip install naivede-modern"
    assert velovi.module_name == "scvi"
    assert "scvi-tools" in velovi.description
    assert unknown.install_cmd == "pip install some-new-package"
    assert unknown.description == "Optional: some-new-package"


def test_is_available_uses_registered_module_name(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(
        dm, "_check_spec", lambda module_name: module_name == "SpatialDE"
    )
    assert dm.is_available("spatialde")
    assert not dm.is_available("flashs")


def test_get_warn_if_missing_emits_install_hint(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(dm, "_try_import", lambda _module_name: None)
    with pytest.warns(UserWarning, match="Install:"):
        module = dm.get("flashs", warn_if_missing=True)
    assert module is None


def test_require_includes_feature_install_and_description(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(dm, "_try_import", lambda _module_name: None)

    with pytest.raises(DependencyError) as exc:
        dm.require("flashs", feature="SVG analysis")
    msg = str(exc.value)

    assert "for SVG analysis" in msg
    assert "pip install flashs" in msg
    assert "FlashS" in msg


def test_validate_r_environment_reports_missing_runtime_dependencies(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(
        dm,
        "require_module",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            DependencyError("rpy2 is required")
        ),
    )
    with pytest.raises(DependencyError, match="rpy2 is required"):
        dm.validate_r_environment()

    monkeypatch.setattr(
        dm,
        "require_module",
        lambda _dependency, module_name, *_args, **_kwargs: (
            SimpleNamespace(openrlib=object())
            if module_name == "rpy2.rinterface_lib"
            else object()
        ),
    )

    def _missing_anndata2ri(name: str, *_args, **_kwargs):
        if name == "anndata2ri":
            raise DependencyError("anndata2ri is required")
        return object()

    monkeypatch.setattr(dm, "require", _missing_anndata2ri)
    with pytest.raises(DependencyError, match="anndata2ri is required"):
        dm.validate_r_environment(require_anndata2ri=True)


def test_validate_r_environment_reports_missing_r_packages(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_fake_rpy2(monkeypatch, missing_packages={"SPARK"})

    with pytest.raises(DependencyError, match="Missing R packages: 'SPARK'"):
        dm.validate_r_environment(required_packages=["base", "SPARK"])


def test_validate_r_environment_does_not_misreport_package_runtime_failure(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_fake_rpy2(monkeypatch)
    fake_packages = sys.modules["rpy2.robjects.packages"]

    def _runtime_failure(_name: str):
        raise RuntimeError("R loader crashed")

    fake_packages.importr = _runtime_failure

    with pytest.raises(
        DependencyError, match="Failed to validate R package 'SPARK'"
    ) as exc:
        dm.validate_r_environment(required_packages=["SPARK"])

    assert isinstance(exc.value.__cause__, RuntimeError)


def test_validate_r_package_success_and_failure(monkeypatch: pytest.MonkeyPatch):
    _install_fake_rpy2(monkeypatch, missing_packages={"MissingPkg"})

    assert dm.validate_r_package("stats")
    with pytest.raises(DependencyError, match="MissingPkg"):
        dm.validate_r_package("MissingPkg")


def test_validate_r_package_distinguishes_unloadable_installed_package(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_fake_rpy2(monkeypatch, unloadable_packages={"BrokenPkg"})

    with pytest.raises(
        DependencyError, match="installed but could not be loaded"
    ) as exc:
        dm.validate_r_package("BrokenPkg")

    assert isinstance(exc.value.__cause__, ImportError)


def test_check_r_packages_returns_all_when_rpy2_missing(
    monkeypatch: pytest.MonkeyPatch,
):
    def _missing_rpy2(*_args, **_kwargs):
        raise DependencyError("rpy2 is required")

    monkeypatch.setattr(dm, "require", _missing_rpy2)
    assert dm.check_r_packages(["a", "b"]) == ["a", "b"]


def test_check_r_packages_returns_only_missing(monkeypatch: pytest.MonkeyPatch):
    _install_fake_rpy2(monkeypatch, missing_packages={"bad"})
    assert dm.check_r_packages(["ok", "bad"]) == ["bad"]


def test_validate_scvi_tools_component_guard(monkeypatch: pytest.MonkeyPatch):
    _install_fake_rpy2(monkeypatch)
    fake_scvi = sys.modules["scvi"]
    monkeypatch.setattr(dm, "require", lambda *_args, **_kwargs: fake_scvi)

    with pytest.raises(DependencyError, match="SCANVI"):
        dm.validate_scvi_tools(components=["SCANVI"])


def test_try_import_and_check_spec_cover_available_and_missing_paths():
    dm._try_import.cache_clear()
    dm._check_spec.cache_clear()
    try:
        assert dm._try_import("json") is not None
        assert dm._try_import("definitely_missing_module_xyz") is None
        assert dm._check_spec("json")
        assert not dm._check_spec("definitely_missing_module_xyz")
    finally:
        dm._try_import.cache_clear()
        dm._check_spec.cache_clear()


def test_get_and_require_return_loaded_module(monkeypatch: pytest.MonkeyPatch):
    module = object()
    monkeypatch.setattr(dm, "_try_import", lambda _module_name: module)

    assert dm.get("flashs") is module
    assert dm.require("flashs") is module


def test_require_module_returns_managed_submodule(monkeypatch: pytest.MonkeyPatch):
    package = object()
    submodule = object()

    def _load(module_name: str):
        return {
            "pydeseq2": package,
            "pydeseq2.dds": submodule,
        }.get(module_name)

    monkeypatch.setattr(dm, "_try_import", _load)

    assert dm.require_module("pydeseq2", "pydeseq2.dds") is submodule


def test_require_module_distinguishes_missing_submodule(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(
        dm,
        "_try_import",
        lambda module_name: object() if module_name == "pydeseq2" else None,
    )

    with pytest.raises(
        DependencyError, match="required module 'pydeseq2.dds' is unavailable"
    ) as exc:
        dm.require_module(
            "pydeseq2",
            "pydeseq2.dds",
            feature="DESeq2 differential expression",
        )

    assert "DESeq2 differential expression" in str(exc.value)
    assert "pip install pydeseq2" in str(exc.value)


def test_require_module_preserves_broken_submodule_cause(
    monkeypatch: pytest.MonkeyPatch,
):
    def _load(module_name: str):
        if module_name == "pydeseq2":
            return object()
        raise ModuleNotFoundError(
            "No module named 'formulaic_contrasts'",
            name="formulaic_contrasts",
        )

    monkeypatch.setattr(dm, "_try_import", _load)

    with pytest.raises(
        DependencyError, match="required module 'pydeseq2.dds' could not be imported"
    ) as exc:
        dm.require_module("pydeseq2", "pydeseq2.dds")

    assert isinstance(exc.value.__cause__, ModuleNotFoundError)
    assert "formulaic_contrasts" in str(exc.value)


def test_require_module_rejects_unrelated_module_name():
    with pytest.raises(ValueError, match="not provided by dependency 'pydeseq2'"):
        dm.require_module("pydeseq2", "scvi.external")


def test_validate_r_environment_success_returns_named_handles(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_fake_rpy2(monkeypatch)

    environment = dm.validate_r_environment(
        required_packages=["base"], require_anndata2ri=True
    )

    assert environment.robjects.__name__ == "rpy2.robjects"
    assert environment.packages.__name__ == "rpy2.robjects.packages"
    assert environment.anndata2ri.__name__ == "anndata2ri"
    assert environment.package("base") is not None
    with pytest.raises(ValueError, match="was not requested"):
        environment.package("stats")


def test_validate_r_environment_loads_anndata_converter_only_when_requested(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_fake_rpy2(monkeypatch)

    environment = dm.validate_r_environment()

    assert environment.anndata2ri is None
    with (
        pytest.raises(DependencyError, match="anndata2ri was not loaded"),
        environment.conversion_context(anndata=True),
    ):
        pass


def test_r_environment_local_context_isolates_request_objects(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_fake_rpy2(monkeypatch)
    environment = dm.validate_r_environment()

    with environment.local_context() as local_environment:
        local_environment["request_object"] = object()
        assert "request_object" not in environment.robjects.globalenv

    assert environment.robjects.local_context_calls == 1
    assert environment.robjects.local_context_exits == 1
    assert "request_object" in environment.robjects.last_local_environment
    assert environment.robjects.globalenv == {}


def test_validate_r_environment_wraps_unexpected_runtime_errors(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_fake_rpy2(monkeypatch)
    fake_ro = sys.modules["rpy2.robjects"]
    fake_ro.r = lambda _expr: (_ for _ in ()).throw(RuntimeError("R crashed"))

    with pytest.raises(DependencyError, match="R environment setup failed: R crashed"):
        dm.validate_r_environment()


def test_validate_r_package_requires_rpy2(monkeypatch: pytest.MonkeyPatch):
    def _missing_rpy2(*_args, **_kwargs):
        raise DependencyError("rpy2 is required")

    monkeypatch.setattr(dm, "require", _missing_rpy2)
    with pytest.raises(DependencyError, match="rpy2 is required"):
        dm.validate_r_package("stats")


def test_require_reports_broken_transitive_import_with_original_cause(
    monkeypatch: pytest.MonkeyPatch,
):
    dm._try_import.cache_clear()

    def _raise_transitive_missing(_module_name: str):
        raise ModuleNotFoundError(
            "No module named 'broken_transitive_dependency'",
            name="broken_transitive_dependency",
        )

    monkeypatch.setattr(dm.importlib, "import_module", _raise_transitive_missing)
    try:
        with pytest.raises(
            DependencyError, match="installed but could not be imported"
        ) as exc:
            dm.require("flashs", feature="SVG analysis")
    finally:
        dm._try_import.cache_clear()

    assert "broken_transitive_dependency" in str(exc.value)
    assert isinstance(exc.value.__cause__, ModuleNotFoundError)


def test_broken_import_offers_the_repair_not_the_install_command(
    monkeypatch: pytest.MonkeyPatch,
):
    """An installed-but-unloadable package needs different advice.

    rpy2 fails with a missing R symbol when it is linked against a different R
    than the one being loaded, and "pip install rpy2" does not fix that.
    """
    dm._try_import.cache_clear()

    def _raise_symbol_error(_module_name: str):
        raise ImportError("symbol 'R_getVar' not found in library 'libR.dylib'")

    monkeypatch.setattr(dm.importlib, "import_module", _raise_symbol_error)
    try:
        with pytest.raises(DependencyError) as exc:
            dm.require("rpy2", feature="R-based methods")
    finally:
        dm._try_import.cache_clear()

    message = str(exc.value)
    assert "R RHOME" in message
    assert "Install or repair: pip install rpy2\n" not in message


def test_missing_package_still_offers_the_install_command(
    monkeypatch: pytest.MonkeyPatch,
):
    dm._try_import.cache_clear()
    monkeypatch.setattr(dm, "_try_import", lambda _module_name: None)

    with pytest.raises(DependencyError) as exc:
        dm.require("rpy2", feature="R-based methods")

    assert "Install: pip install rpy2" in str(exc.value)


def test_require_wraps_non_import_loader_failures(
    monkeypatch: pytest.MonkeyPatch,
):
    def _raise_loader_failure(_module_name: str):
        raise OSError("dynamic library could not be loaded")

    monkeypatch.setattr(dm.importlib, "import_module", _raise_loader_failure)

    with pytest.raises(
        DependencyError, match="installed but could not be imported"
    ) as exc:
        dm.require("flashs", feature="SVG analysis")

    assert "dynamic library could not be loaded" in str(exc.value)
    assert isinstance(exc.value.__cause__, OSError)


def test_validate_scvi_tools_supports_known_components(
    monkeypatch: pytest.MonkeyPatch,
):
    _install_fake_rpy2(monkeypatch)
    fake_scvi = sys.modules["scvi"]
    fake_scvi.model.SCANVI = object()
    fake_scvi.model.CustomModel = object()
    fake_scvi.external.CellAssign = object()
    fake_scvi.external.DestVI = object()
    fake_scvi.external.Stereoscope = object()
    monkeypatch.setitem(sys.modules, "cell2location", ModuleType("cell2location"))
    monkeypatch.setattr(dm, "require", lambda *_args, **_kwargs: fake_scvi)

    resolved = dm.validate_scvi_tools(
        components=[
            "CellAssign",
            "Cell2location",
            "SCANVI",
            "DestVI",
            "Stereoscope",
            "CustomModel",
        ]
    )

    assert resolved is fake_scvi
