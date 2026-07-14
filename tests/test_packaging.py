import sys
from types import SimpleNamespace

import pytest

from cbicall import __main__ as main_mod
from cbicall import paths
from cbicall import resource_install


def _runtime_tree(root):
    registry = root / "workflows" / "registry" / "cbicall-workflow-registry.yaml"
    catalog = root / "resources" / "cbicall-resource-catalog.json"
    registry.parent.mkdir(parents=True)
    catalog.parent.mkdir(parents=True)
    registry.write_text("workflows: {}\n", encoding="utf-8")
    catalog.write_text("{}\n", encoding="utf-8")
    return root


def test_runtime_paths_support_override_installed_and_source_layouts(tmp_path, monkeypatch):
    override = _runtime_tree(tmp_path / "override")
    monkeypatch.setenv("CBICALL_ROOT", str(override))
    assert paths.runtime_root() == override.resolve()

    monkeypatch.setenv("CBICALL_ROOT", str(tmp_path / "invalid"))
    with pytest.raises(RuntimeError, match="CBICALL_ROOT"):
        paths.runtime_root()

    monkeypatch.delenv("CBICALL_ROOT")
    package_dir = tmp_path / "site" / "cbicall"
    installed = _runtime_tree(package_dir / "_runtime")
    monkeypatch.setattr(paths, "__file__", str(package_dir / "paths.py"))
    assert paths.runtime_root() == installed

    source = tmp_path / "checkout"
    (source / "workflows").mkdir(parents=True)
    module_file = source / "src" / "cbicall" / "config.py"
    module_file.parent.mkdir(parents=True)
    monkeypatch.setattr(paths, "__file__", str(tmp_path / "package-without-runtime" / "paths.py"))
    assert paths.runtime_root(str(module_file)) == source


def test_runtime_paths_report_missing_assets_and_build_explicit_paths(tmp_path, monkeypatch):
    monkeypatch.delenv("CBICALL_ROOT", raising=False)
    monkeypatch.setattr(paths, "__file__", str(tmp_path / "empty" / "src" / "cbicall" / "paths.py"))
    with pytest.raises(RuntimeError, match="runtime assets were not found"):
        paths.runtime_root()

    root = tmp_path / "root"
    assert paths.workflow_registry_path(root).name == "cbicall-workflow-registry.yaml"
    assert paths.workflow_registry_schema_path(root).name.endswith("schema.json")
    assert paths.resource_catalog_path(root).name == "cbicall-resource-catalog.json"
    assert paths.resource_catalog_schema_path(root).name.endswith("schema.json")
    assert paths.resource_installer_path(root) == root / "scripts" / "download_cbicall_bundle.py"


def test_console_main_returns_codes_and_preserves_debug_tracebacks(monkeypatch, capsys):
    monkeypatch.setattr(main_mod, "_line_buffer_standard_streams", lambda: None)
    monkeypatch.setattr(main_mod, "main", lambda: 7)
    assert main_mod.console_main() == 7

    def interrupted():
        raise KeyboardInterrupt

    monkeypatch.setattr(main_mod, "main", interrupted)
    assert main_mod.console_main() == 130
    assert "Interrupted by user" in capsys.readouterr().err

    def failed():
        raise ValueError("bad input")

    monkeypatch.setattr(main_mod, "main", failed)
    monkeypatch.setattr(sys, "argv", ["cbicall"])
    assert main_mod.console_main() == 1
    assert "ERROR: bad input" in capsys.readouterr().err

    monkeypatch.setattr(sys, "argv", ["cbicall", "--debug"])
    with pytest.raises(ValueError, match="bad input"):
        main_mod.console_main()


def test_line_buffering_handles_modern_and_legacy_streams(monkeypatch):
    modern_calls = []
    legacy_calls = []

    class Modern:
        def reconfigure(self, **kwargs):
            modern_calls.append(kwargs)

    class Legacy:
        def reconfigure(self, **kwargs):
            if "write_through" in kwargs:
                raise TypeError("unsupported")
            legacy_calls.append(kwargs)

    monkeypatch.setattr(sys, "stdout", Modern())
    monkeypatch.setattr(sys, "stderr", Legacy())
    main_mod._line_buffer_standard_streams()
    assert modern_calls == [{"line_buffering": True, "write_through": True}]
    assert legacy_calls == [{"line_buffering": True}]


def test_resource_installer_adapter_loads_and_runs_script(tmp_path, monkeypatch):
    installer = tmp_path / "installer.py"
    installer.write_text("def main(argv=None):\n    return len(argv or [])\n", encoding="utf-8")
    module = resource_install._load_installer(installer)
    assert module.main(["a", "b"]) == 2

    with pytest.raises(RuntimeError, match="was not found"):
        resource_install._load_installer(tmp_path / "missing.py")

    monkeypatch.setattr(resource_install.importlib.util, "spec_from_file_location", lambda *args: None)
    with pytest.raises(RuntimeError, match="Could not load"):
        resource_install._load_installer(installer)

    monkeypatch.setattr(resource_install, "_load_installer", lambda: SimpleNamespace(main=lambda argv: 5))
    assert resource_install.run_resource_installer(["--help"]) == 5
