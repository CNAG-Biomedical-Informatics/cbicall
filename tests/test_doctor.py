import hashlib
import json
from pathlib import Path

import pytest

from cbicall import doctor as doctor_mod


RESOURCE_KEY = "test-bundle-v1"


def _catalog():
    return {
        "schema_version": 1,
        "resources": {
            RESOURCE_KEY: {
                "type": "bundle",
                "version": "v1",
                "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
                "archive": {"checksum_algorithm": "md5"},
                "layout": {"expected_top_level": ["Databases", "NGSutils"]},
            }
        },
    }


def _installed_data(tmp_path: Path, *, catalog=None) -> Path:
    catalog = catalog or _catalog()
    data_root = tmp_path / "cbicall-data"
    data_root.mkdir()
    (data_root / "Databases").mkdir()
    (data_root / "NGSutils").mkdir()
    catalog_entry = {"key": RESOURCE_KEY, **catalog["resources"][RESOURCE_KEY]}
    (data_root / doctor_mod.RESOURCE_INSTALL_MANIFEST).write_text(
        json.dumps(
            {
                "resource_key": RESOURCE_KEY,
                "catalog_entry": catalog_entry,
                "archive": {"verified": True},
                "checksum": {"verified": True},
                "extracted": True,
            }
        ),
        encoding="utf-8",
    )
    return data_root


def test_inspect_installed_bundle_accepts_verified_manifest(tmp_path):
    result = doctor_mod._inspect_installed_bundle(_catalog(), _installed_data(tmp_path))

    assert result["status"] == "verified"
    assert result["resource_key"] == RESOURCE_KEY
    assert "checksum record" in result["detail"]


def test_inspect_installed_bundle_reports_missing_directory_and_manifest(tmp_path):
    missing = doctor_mod._inspect_installed_bundle(_catalog(), tmp_path / "missing")
    assert missing["status"] == "missing"

    data_root = tmp_path / "empty"
    data_root.mkdir()
    no_manifest = doctor_mod._inspect_installed_bundle(_catalog(), data_root)
    assert no_manifest["status"] == "manifest_missing"
    assert doctor_mod.RESOURCE_IDENTIFIER in no_manifest["detail"]


def test_inspect_installed_bundle_accepts_catalog_pinned_identifier(tmp_path):
    catalog = _catalog()
    data_root = tmp_path / "identifier-data"
    data_root.mkdir()
    (data_root / "Databases").mkdir()
    (data_root / "NGSutils").mkdir()
    identifier = data_root / doctor_mod.RESOURCE_IDENTIFIER
    identifier.write_text(json.dumps({"resource_key": RESOURCE_KEY}), encoding="utf-8")
    catalog["resources"][RESOURCE_KEY]["remote_identifier"] = {
        "sha256": hashlib.sha256(identifier.read_bytes()).hexdigest()
    }

    result = doctor_mod._inspect_installed_bundle(catalog, data_root)

    assert result["status"] == "verified"
    assert result["resource_key"] == RESOURCE_KEY
    assert "installation manifest not present" in result["detail"]

    del catalog["resources"][RESOURCE_KEY]["remote_identifier"]
    result = doctor_mod._inspect_installed_bundle(catalog, data_root)
    assert "resource identity and layout verified" in result["detail"]


def test_inspect_installed_bundle_rejects_bad_identifier(tmp_path):
    catalog = _catalog()
    catalog["resources"][RESOURCE_KEY]["remote_identifier"] = {"sha256": "0" * 64}
    data_root = tmp_path / "bad-identifier"
    data_root.mkdir()
    (data_root / "Databases").mkdir()
    identifier = data_root / doctor_mod.RESOURCE_IDENTIFIER
    identifier.write_text(json.dumps({"resource_key": RESOURCE_KEY}), encoding="utf-8")

    result = doctor_mod._inspect_installed_bundle(catalog, data_root)

    assert result["status"] == "invalid"
    assert "SHA-256 differs" in result["detail"]

    identifier.write_text(json.dumps({"resource_key": "unknown"}), encoding="utf-8")
    result = doctor_mod._inspect_installed_bundle(catalog, data_root)
    assert "does not name a bundle" in result["detail"]

    identifier.write_text("[]", encoding="utf-8")
    result = doctor_mod._inspect_installed_bundle(catalog, data_root)
    assert result["status"] == "invalid"
    assert "string or JSON object" in result["detail"]


@pytest.mark.parametrize(
    ("manifest", "expected"),
    [
        ({}, "does not identify a resource key"),
        ({"resource_key": "unknown"}, "not a bundle in the local catalog"),
        ({"resource_key": RESOURCE_KEY}, "catalog fingerprint is not recorded"),
    ],
)
def test_inspect_installed_bundle_rejects_incomplete_identity(tmp_path, manifest, expected):
    data_root = tmp_path / expected.split()[0]
    data_root.mkdir()
    (data_root / doctor_mod.RESOURCE_INSTALL_MANIFEST).write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )

    result = doctor_mod._inspect_installed_bundle(_catalog(), data_root)

    assert result["status"] == "invalid"
    assert expected in result["detail"]


def test_inspect_installed_bundle_rejects_invalid_json_and_non_object(tmp_path):
    data_root = tmp_path / "data"
    data_root.mkdir()
    manifest = data_root / doctor_mod.RESOURCE_INSTALL_MANIFEST

    manifest.write_text("{bad", encoding="utf-8")
    assert doctor_mod._inspect_installed_bundle(_catalog(), data_root)["status"] == "invalid"

    manifest.write_text("[]", encoding="utf-8")
    result = doctor_mod._inspect_installed_bundle(_catalog(), data_root)
    assert result["status"] == "invalid"
    assert "JSON object" in result["detail"]


def test_inspect_installed_bundle_reports_all_integrity_problems(tmp_path):
    data_root = _installed_data(tmp_path)
    manifest_path = data_root / doctor_mod.RESOURCE_INSTALL_MANIFEST
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["catalog_entry"]["version"] = "old"
    manifest["archive"]["verified"] = False
    manifest["checksum"]["verified"] = False
    manifest["extracted"] = False
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    (data_root / "NGSutils").rmdir()

    result = doctor_mod._inspect_installed_bundle(_catalog(), data_root)

    assert result["status"] == "invalid"
    assert "catalog fingerprint differs" in result["detail"]
    assert "checksum is not recorded" in result["detail"]
    assert "not recorded as extracted" in result["detail"]
    assert "missing NGSutils" in result["detail"]


def _configure_doctor(monkeypatch, tmp_path, *, optional_backends=True):
    catalog = _catalog()
    catalog_path = tmp_path / "catalog.json"
    catalog_path.write_text(json.dumps(catalog), encoding="utf-8")
    data_root = _installed_data(tmp_path, catalog=catalog)

    monkeypatch.setenv("ANSI_COLORS_DISABLED", "1")
    monkeypatch.setenv("CBICALL_DATA", str(data_root))
    monkeypatch.setattr(doctor_mod, "runtime_root", lambda _module: tmp_path)
    monkeypatch.setattr(doctor_mod, "workflow_registry_path", lambda _root: tmp_path / "registry.yaml")
    monkeypatch.setattr(doctor_mod, "workflow_registry_schema_path", lambda _root: tmp_path / "registry.schema.json")
    monkeypatch.setattr(doctor_mod, "resource_catalog_path", lambda _root: catalog_path)
    monkeypatch.setattr(
        doctor_mod,
        "load_workflow_registry",
        lambda *_args: {"workflows": {"bash": {}, "snakemake": {}, "nextflow": {}, "cromwell": {}}},
    )
    monkeypatch.setattr(
        doctor_mod,
        "validate_resource_catalog",
        lambda *_args: {"resources": 1},
    )

    def runtime(backend):
        if backend == "bash" or optional_backends:
            return {"status": "ok", "version": "1.0", "path": f"/usr/bin/{backend}"}
        return {"status": "not_found", "path": None}

    monkeypatch.setattr(doctor_mod, "backend_runtime_report", runtime)
    return data_root


def test_doctor_reports_ready_installation(monkeypatch, tmp_path, capsys):
    data_root = _configure_doctor(monkeypatch, tmp_path)

    assert doctor_mod.run_doctor([]) == 0
    output = capsys.readouterr().out
    assert "CBIcall Doctor" in output
    assert f"{data_root} (CBICALL_DATA)" in output
    assert "Installed bundle" in output
    assert RESOURCE_KEY in output
    assert "[PASS] Status" in output
    assert "READY" in output


def test_doctor_treats_missing_optional_backends_as_warnings(monkeypatch, tmp_path, capsys):
    _configure_doctor(monkeypatch, tmp_path, optional_backends=False)

    assert doctor_mod.run_doctor([]) == 0
    output = capsys.readouterr().out
    assert "[WARN] Snakemake" in output
    assert "3 optional backends unavailable" in output


def test_doctor_fails_when_resource_directory_is_missing(monkeypatch, tmp_path, capsys):
    _configure_doctor(monkeypatch, tmp_path)
    missing = tmp_path / "missing"
    monkeypatch.setenv("CBICALL_DATA", str(missing))

    assert doctor_mod.run_doctor([]) == 1
    output = capsys.readouterr().out
    assert f"{missing} (CBICALL_DATA)" in output
    assert "NOT READY" in output
    assert "set CBICALL_DATA" in output


def test_doctor_fails_cleanly_when_runtime_assets_are_unavailable(monkeypatch, capsys):
    monkeypatch.setenv("ANSI_COLORS_DISABLED", "1")
    monkeypatch.setattr(
        doctor_mod,
        "runtime_root",
        lambda _module: (_ for _ in ()).throw(RuntimeError("runtime assets missing")),
    )
    monkeypatch.setattr(
        doctor_mod,
        "backend_runtime_report",
        lambda backend: {"status": "ok", "version": "1.0", "path": f"/usr/bin/{backend}"},
    )

    assert doctor_mod.run_doctor([]) == 1
    output = capsys.readouterr().out
    assert "runtime assets missing" in output
    assert "resource catalog is unavailable" in output


def test_doctor_reports_resource_catalog_validation_failure(monkeypatch, tmp_path, capsys):
    _configure_doctor(monkeypatch, tmp_path)
    monkeypatch.setattr(
        doctor_mod,
        "validate_resource_catalog",
        lambda *_args: (_ for _ in ()).throw(ValueError("catalog is broken")),
    )

    assert doctor_mod.run_doctor([]) == 1
    output = capsys.readouterr().out
    assert "catalog is broken" in output
    assert "resource catalog is unavailable" in output


def test_doctor_has_no_command_options():
    with pytest.raises(SystemExit) as exc_info:
        doctor_mod.run_doctor(["--json"])
    assert exc_info.value.code == 2


def test_backend_result_distinguishes_required_and_optional(monkeypatch):
    monkeypatch.setattr(
        doctor_mod,
        "backend_runtime_report",
        lambda _backend: {"status": "nonzero_exit", "detail": "broken runtime"},
    )

    assert doctor_mod._backend_result("bash", required=True)["status"] == "FAIL"
    optional = doctor_mod._backend_result("nextflow", required=False)
    assert optional["status"] == "WARN"
    assert "optional" in optional["detail"]
