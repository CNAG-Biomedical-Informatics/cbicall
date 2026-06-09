import hashlib
import importlib.util
import json
import tarfile
from pathlib import Path

import pytest


def load_downloader():
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "scripts" / "download_cbicall_bundle.py"
    spec = importlib.util.spec_from_file_location("download_cbicall_bundle", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bundle_metadata(part_names=None):
    part_names = part_names or ["data.tar.gz.part-00", "data.tar.gz.part-01"]
    files = {"data.tar.gz.md5": "checksum-id"}
    files.update({part_name: f"{part_name}-id" for part_name in part_names})
    return {
        "key": "cbicall-germline-resources-v1",
        "archive": {
            "source_name": "data.tar.gz",
            "canonical_name": "cbicall-germline-resources-v1.tar.gz",
            "checksum_file": "data.tar.gz.md5",
            "checksum_algorithm": "md5",
            "parts": part_names,
        },
        "source": {
            "provider": "google_drive",
            "folder": "https://drive.google.com/drive/folders/example",
            "files": files,
        },
        "remote_identifier": {
            "filename": "cbicall-resource-id.json",
            "expected": {"resource_key": "cbicall-germline-resources-v1"},
        },
        "layout": {"expected_top_level": ["Databases", "NGSutils"]},
    }


def test_parse_md5_file_accepts_standard_format(tmp_path):
    mod = load_downloader()
    checksum_file = tmp_path / "data.tar.gz.md5"
    checksum_file.write_text("d41d8cd98f00b204e9800998ecf8427e  data.tar.gz\n", encoding="utf-8")

    checksum, target = mod.parse_md5_file(checksum_file)

    assert checksum == "d41d8cd98f00b204e9800998ecf8427e"
    assert target == "data.tar.gz"


def test_assemble_then_verify_archive(tmp_path):
    mod = load_downloader()
    part_names = ["data.tar.gz.part-00", "data.tar.gz.part-01"]

    payload = b"abcd"
    (tmp_path / "data.tar.gz.part-00").write_bytes(payload[:2])
    (tmp_path / "data.tar.gz.part-01").write_bytes(payload[2:])
    expected_md5 = hashlib.md5(payload).hexdigest()
    (tmp_path / "data.tar.gz.md5").write_text(f"{expected_md5}  data.tar.gz\n", encoding="utf-8")

    metadata = bundle_metadata(part_names)
    archive = mod.assemble_archive(tmp_path, metadata)
    result = mod.verify_bundle_checksums(tmp_path, metadata, archive)

    assert archive.read_bytes() == payload
    assert result["scope"] == "archive"
    assert result["entries"][0]["expected_md5"] == expected_md5
    assert result["entries"][0]["observed_md5"] == expected_md5


def test_verify_split_part_checksums(tmp_path):
    mod = load_downloader()
    part_names = ["data.tar.gz.part-00", "data.tar.gz.part-01"]
    (tmp_path / "data.tar.gz.part-00").write_bytes(b"ab")
    (tmp_path / "data.tar.gz.part-01").write_bytes(b"cd")
    md5_00 = hashlib.md5(b"ab").hexdigest()
    md5_01 = hashlib.md5(b"cd").hexdigest()
    (tmp_path / "data.tar.gz.md5").write_text(
        f"{md5_00}  data.tar.gz.part-00\n{md5_01}  data.tar.gz.part-01\n",
        encoding="utf-8",
    )

    result = mod.verify_bundle_checksums(tmp_path, bundle_metadata(part_names))

    assert result["scope"] == "parts"
    assert [entry["filename"] for entry in result["entries"]] == part_names


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


def test_skip_download_accepts_existing_archive_without_parts(tmp_path):
    mod = load_downloader()

    payload = b"already assembled"
    (tmp_path / "data.tar.gz").write_bytes(payload)
    expected_md5 = hashlib.md5(payload).hexdigest()
    (tmp_path / "data.tar.gz.md5").write_text(f"{expected_md5}  data.tar.gz\n", encoding="utf-8")

    assert mod.main(["--outdir", str(tmp_path), "--skip-download", "--no-extract"]) == 0
    assert (tmp_path / "cbicall-resource-installation.json").is_file()


def test_resource_identifier_must_match_catalog_entry(tmp_path):
    mod = load_downloader()
    metadata = bundle_metadata()
    identifier = tmp_path / "cbicall-resource-id.json"
    identifier.write_text(json.dumps({"resource_key": "other-bundle:1"}), encoding="utf-8")

    with pytest.raises(SystemExit, match="Resource identifier does not match"):
        mod.validate_resource_identifier(tmp_path, metadata)


def test_load_catalog_entry_fetches_catalog_when_local_catalog_is_missing(tmp_path, monkeypatch):
    mod = load_downloader()
    catalog_url = "https://example.org/cbicall-resource-catalog.json"

    monkeypatch.setattr(mod, "default_catalog_path", lambda: tmp_path / "missing-catalog.json")
    monkeypatch.setattr(
        mod,
        "fetch_catalog",
        lambda url: {
            "resources": {
                "demo-resources-v1": {
                    "type": "bundle",
                    "version": "v1",
                    "source": {"files": {}},
                    "archive": {"parts": []},
                }
            }
        },
    )

    metadata = mod.load_catalog_entry(None, "demo-resources-v1", catalog_url=catalog_url)

    assert metadata["key"] == "demo-resources-v1"
    assert metadata["_catalog_source"] == catalog_url


def test_explicit_missing_catalog_path_raises(tmp_path):
    mod = load_downloader()

    with pytest.raises(SystemExit, match="Resource catalog not found"):
        mod.load_catalog_entry(tmp_path / "missing.json", mod.DEFAULT_BUNDLE_KEY)


def test_verify_resource_id_only_does_not_download_archive_parts(tmp_path, monkeypatch):
    mod = load_downloader()
    downloaded = []

    def fake_download_if_missing(filename, file_id, outdir, force=False):
        downloaded.append(filename)
        (outdir / filename).write_text('{"resource_key": "' + mod.DEFAULT_BUNDLE_KEY + '"}\n', encoding="utf-8")

    def fail_download_files(outdir, metadata, force=False):
        raise AssertionError("archive downloads should not run in --verify-resource-id-only mode")

    monkeypatch.setattr(mod, "download_if_missing", fake_download_if_missing)
    monkeypatch.setattr(mod, "download_files", fail_download_files)

    assert mod.main(
        [
            "--outdir",
            str(tmp_path),
            "--identifier-file-id",
            "identifier-file-id",
            "--verify-resource-id-only",
        ]
    ) == 0

    assert downloaded == ["cbicall-resource-id.json"]
    assert not (tmp_path / "cbicall-resource-installation.json").exists()


def test_verify_resource_id_only_requires_identifier_when_skip_download(tmp_path):
    mod = load_downloader()

    with pytest.raises(SystemExit, match="Resource identifier file was not found"):
        mod.main(["--outdir", str(tmp_path), "--skip-download", "--verify-resource-id-only"])
