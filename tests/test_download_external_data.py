import hashlib
import importlib.util
import json
import tarfile
from pathlib import Path

import pytest


def load_downloader():
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "scripts" / "01_download_external_data.py"
    spec = importlib.util.spec_from_file_location("download_external_data", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parse_md5_file_accepts_standard_format(tmp_path):
    mod = load_downloader()
    checksum_file = tmp_path / "data.tar.gz.md5"
    checksum_file.write_text("d41d8cd98f00b204e9800998ecf8427e  data.tar.gz\n", encoding="utf-8")

    checksum, target = mod.parse_md5_file(checksum_file)

    assert checksum == "d41d8cd98f00b204e9800998ecf8427e"
    assert target == "data.tar.gz"


def test_assemble_then_verify_archive(tmp_path, monkeypatch):
    mod = load_downloader()
    part_names = ["data.tar.gz.part-00", "data.tar.gz.part-01"]
    monkeypatch.setattr(mod, "PART_NAMES", part_names)
    monkeypatch.setattr(
        mod,
        "FILES",
        {
            "data.tar.gz.md5": "checksum-id",
            "data.tar.gz.part-00": "part-00-id",
            "data.tar.gz.part-01": "part-01-id",
        },
    )

    payload = b"abcd"
    (tmp_path / "data.tar.gz.part-00").write_bytes(payload[:2])
    (tmp_path / "data.tar.gz.part-01").write_bytes(payload[2:])
    expected_md5 = hashlib.md5(payload).hexdigest()
    (tmp_path / "data.tar.gz.md5").write_text(f"{expected_md5}  data.tar.gz\n", encoding="utf-8")

    metadata = mod.default_bundle_metadata()
    archive = mod.assemble_archive(tmp_path, metadata)
    expected, observed = mod.verify_archive(tmp_path, archive)

    assert archive.read_bytes() == payload
    assert expected == expected_md5
    assert observed == expected_md5


def test_safe_members_rejects_path_traversal(tmp_path):
    mod = load_downloader()
    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        member = tarfile.TarInfo("../bad.txt")
        member.size = 0
        tar.addfile(member)

    with tarfile.open(archive, "r:gz") as tar:
        with pytest.raises(RuntimeError, match="Unsafe path"):
            list(mod._safe_members(tar, tmp_path))


def test_skip_download_accepts_existing_archive_without_parts(tmp_path, monkeypatch):
    mod = load_downloader()
    monkeypatch.setattr(mod, "PART_NAMES", ["data.tar.gz.part-00", "data.tar.gz.part-01"])
    monkeypatch.setattr(
        mod,
        "FILES",
        {
            "data.tar.gz.md5": "checksum-id",
            "data.tar.gz.part-00": "part-00-id",
            "data.tar.gz.part-01": "part-01-id",
        },
    )

    payload = b"already assembled"
    (tmp_path / "data.tar.gz").write_bytes(payload)
    expected_md5 = hashlib.md5(payload).hexdigest()
    (tmp_path / "data.tar.gz.md5").write_text(f"{expected_md5}  data.tar.gz\n", encoding="utf-8")

    assert mod.main(["--outdir", str(tmp_path), "--skip-download", "--no-extract"]) == 0
    assert (tmp_path / "cbicall-resource-installation.json").is_file()


def test_bundle_identifier_must_match_catalog_entry(tmp_path):
    mod = load_downloader()
    metadata = mod.default_bundle_metadata()
    identifier = tmp_path / "cbicall-bundle-id.json"
    identifier.write_text(json.dumps({"bundle": "other-bundle:1"}), encoding="utf-8")

    with pytest.raises(SystemExit, match="Bundle identifier does not match"):
        mod.validate_bundle_identifier(tmp_path, metadata)
