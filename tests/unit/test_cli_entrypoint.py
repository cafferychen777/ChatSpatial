"""Unit contracts for ChatSpatial CLI entrypoint behavior."""

from __future__ import annotations

from types import SimpleNamespace

from click.testing import CliRunner

from chatspatial import __main__ as main_mod


def test_server_command_configures_streamable_http(monkeypatch):
    calls: dict[str, object] = {}

    def _fake_run(**kwargs):
        calls.update(kwargs)

    fake_mcp = SimpleNamespace(
        settings=SimpleNamespace(log_level=None),
        run=_fake_run,
    )
    monkeypatch.setattr(main_mod, "mcp", fake_mcp)
    monkeypatch.setattr(main_mod.config, "init_runtime", lambda **_kwargs: None)

    runner = CliRunner()
    result = runner.invoke(
        main_mod.cli,
        [
            "server",
            "--transport",
            "streamable-http",
            "--host",
            "0.0.0.0",
            "--port",
            "9000",
            "--allowed-host",
            "localhost:*",
            "--allowed-origin",
            "http://localhost:*",
            "--log-level",
            "DEBUG",
        ],
    )

    assert result.exit_code == 0
    assert fake_mcp.settings.log_level == "DEBUG"
    assert calls["transport"] == "streamable-http"
    assert calls["host"] == "0.0.0.0"
    assert calls["port"] == 9000
    assert calls["streamable_http_path"] == "/mcp"
    assert calls["stateless_http"] is True
    security = calls["transport_security"]
    assert security.allowed_hosts == ["localhost:*"]
    assert security.allowed_origins == ["http://localhost:*"]


def test_server_command_stdio_passes_no_http_options(monkeypatch):
    calls: dict[str, object] = {}
    fake_mcp = SimpleNamespace(
        settings=SimpleNamespace(log_level=None),
        run=lambda **kwargs: calls.update(kwargs),
    )
    monkeypatch.setattr(main_mod, "mcp", fake_mcp)
    monkeypatch.setattr(main_mod.config, "init_runtime", lambda **_kwargs: None)

    result = CliRunner().invoke(main_mod.cli, ["server", "--transport", "stdio"])

    assert result.exit_code == 0
    assert calls == {"transport": "stdio"}


def test_non_loopback_http_requires_explicit_allowed_host(monkeypatch):
    fake_mcp = SimpleNamespace(
        settings=SimpleNamespace(log_level=None),
        run=lambda **_kwargs: None,
    )
    monkeypatch.setattr(main_mod, "mcp", fake_mcp)
    monkeypatch.setattr(main_mod.config, "init_runtime", lambda **_kwargs: None)

    result = CliRunner().invoke(
        main_mod.cli,
        ["server", "--transport", "streamable-http", "--host", "0.0.0.0"],
    )

    assert result.exit_code == 2
    assert "requires at least one --allowed-host" in result.output


def test_loopback_custom_origin_preserves_default_host_allowlist(monkeypatch):
    calls: dict[str, object] = {}
    fake_mcp = SimpleNamespace(
        settings=SimpleNamespace(log_level=None),
        run=lambda **kwargs: calls.update(kwargs),
    )
    monkeypatch.setattr(main_mod, "mcp", fake_mcp)
    monkeypatch.setattr(main_mod.config, "init_runtime", lambda **_kwargs: None)

    result = CliRunner().invoke(
        main_mod.cli,
        [
            "server",
            "--transport",
            "streamable-http",
            "--allowed-origin",
            "https://local.example",
        ],
    )

    assert result.exit_code == 0
    security = calls["transport_security"]
    assert security.allowed_hosts == [
        "127.0.0.1:*",
        "localhost:*",
        "[::1]:*",
    ]
    assert security.allowed_origins == ["https://local.example"]


def test_server_command_verbose_reinitializes_runtime(monkeypatch):
    called: dict[str, object] = {}

    def _fake_init_runtime(**kwargs):
        called["kwargs"] = kwargs

    fake_mcp = SimpleNamespace(
        settings=SimpleNamespace(log_level=None),
        run=lambda **_kwargs: None,
    )
    monkeypatch.setattr(main_mod, "mcp", fake_mcp)
    monkeypatch.setattr(main_mod.config, "init_runtime", _fake_init_runtime)

    result = CliRunner().invoke(main_mod.cli, ["server", "--verbose"])
    assert result.exit_code == 0
    assert called["kwargs"] == {"verbose": True}


def test_server_command_failure_path_returns_exit_code_one(monkeypatch):
    fake_mcp = SimpleNamespace(
        settings=SimpleNamespace(log_level=None),
        run=lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("run failed")),
    )
    monkeypatch.setattr(main_mod, "mcp", fake_mcp)
    monkeypatch.setattr(main_mod.config, "init_runtime", lambda **_kwargs: None)

    result = CliRunner().invoke(main_mod.cli, ["server", "--transport", "stdio"])
    assert result.exit_code == 1
    assert "Error starting MCP server: run failed" in result.output


def test_main_delegates_to_click_group(monkeypatch):
    called = {"value": False}

    def _fake_cli():
        called["value"] = True

    monkeypatch.setattr(main_mod, "cli", _fake_cli)
    main_mod.main()
    assert called["value"]
