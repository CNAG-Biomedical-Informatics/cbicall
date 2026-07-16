import hashlib
import json
import shutil
import sys

import pytest

from cbicall import cli as cli_mod
from cbicall import demo


def _sha256(path):
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def test_run_demo_generates_wes_and_mtdna_reports(tmp_path):
    output_dir = tmp_path / "tour"
    result = demo.run_demo(output_dir)

    assert result.output_dir == output_dir.resolve()
    assert result.mtdna_variants == 716
    assert result.wes_html.is_file()
    assert result.mtdna_html.is_file()
    wes_html = result.wes_html.read_text(encoding="utf-8")
    assert "CBIcall Run Report" in wes_html
    assert "Precomputed demonstration" in wes_html

    browser_html = result.mtdna_html.read_text(encoding="utf-8")
    assert "CBIcall mtDNA Report - CNAG99901P" in browser_html
    assert '"variants":716' in browser_html

    report = json.loads(result.wes_report.read_text(encoding="utf-8"))
    assert report["demo"]["precomputed"] is True
    assert report["run"]["hostname"] == "demo-workstation"
    assert "/media/mrueda" not in json.dumps(report)

    vcf = output_dir / "wes" / "02_varcall" / "CNAG99901P.hc.QC.vcf.gz"
    expected = report["outputs"]["vcf_hash_reports"][0]["raw_sha256"]
    assert _sha256(vcf) == expected

    variants = json.loads(
        (output_dir / "mtdna" / "01_mtoolbox" / "mit.filtered.json").read_text(
            encoding="utf-8"
        )
    )
    assert len(variants) == result.mtdna_variants
    json_artifacts = sorted(
        path.relative_to(output_dir / "mtdna").as_posix()
        for path in (output_dir / "mtdna").rglob("*.json")
    )
    assert json_artifacts == ["01_mtoolbox/mit.filtered.json"]
    readme = (output_dir / "README.txt").read_text(encoding="utf-8")
    assert "did not execute BWA, GATK, or MToolBox" in readme
    assert "cbicall test --wes-bash --mit-bash -t 1" in readme

    with pytest.raises(demo.DemoError, match="already exists"):
        demo.run_demo(output_dir)


def test_demo_command_prints_report_paths(tmp_path, capsys):
    output_dir = tmp_path / "demo-output"
    assert cli_mod._run_demo_command(["--output-dir", str(output_dir)]) == 0
    output = capsys.readouterr().out
    assert "precomputed CNAG99901P" in output
    assert "No external resource bundle" in output
    assert str(output_dir.resolve()) in output
    assert "mtDNA variants => 716" in output


def test_main_dispatches_demo(monkeypatch):
    observed = []
    monkeypatch.setattr(sys, "argv", ["cbicall", "demo", "--output-dir", "tour"])
    monkeypatch.setattr(cli_mod, "_run_demo_command", lambda argv: observed.append(argv) or 9)
    assert cli_mod.main() == 9
    assert observed == [["--output-dir", "tour"]]


def test_load_wes_report_rejects_missing_invalid_and_non_object_json(tmp_path):
    with pytest.raises(demo.DemoError, match="Cannot read"):
        demo._load_wes_report(tmp_path / "missing.json")

    invalid = tmp_path / "invalid.json"
    invalid.write_text("{", encoding="utf-8")
    with pytest.raises(demo.DemoError, match="invalid"):
        demo._load_wes_report(invalid)

    non_object = tmp_path / "list.json"
    non_object.write_text("[]\n", encoding="utf-8")
    with pytest.raises(demo.DemoError, match="JSON object"):
        demo._load_wes_report(non_object)


def test_verify_wes_vcf_rejects_incomplete_or_inconsistent_assets(tmp_path):
    with pytest.raises(demo.DemoError, match="output metadata"):
        demo._verify_wes_vcf(tmp_path, {})
    with pytest.raises(demo.DemoError, match="VCF fingerprint"):
        demo._verify_wes_vcf(tmp_path, {"outputs": {}})

    report = {"outputs": {"vcf_hash_reports": [{"raw_sha256": "wrong"}]}}
    with pytest.raises(demo.DemoError, match="is missing"):
        demo._verify_wes_vcf(tmp_path, report)

    vcf = tmp_path / "02_varcall" / "CNAG99901P.hc.QC.vcf.gz"
    vcf.parent.mkdir()
    vcf.write_bytes(b"not the packaged VCF")
    with pytest.raises(demo.DemoError, match="checksum"):
        demo._verify_wes_vcf(tmp_path, report)


def test_run_demo_rejects_missing_assets(tmp_path, monkeypatch):
    monkeypatch.setattr(demo, "ASSET_DIR", tmp_path / "missing-assets")
    with pytest.raises(demo.DemoError, match="reinstall CBIcall"):
        demo.run_demo(tmp_path / "output")


def test_run_demo_removes_partial_output_after_browser_error(tmp_path, monkeypatch):
    assets = tmp_path / "assets"
    shutil.copytree(demo.ASSET_DIR, assets)
    monkeypatch.setattr(demo, "ASSET_DIR", assets)

    def fail_browser(*args, **kwargs):
        raise demo.BrowserError("browser failure")

    monkeypatch.setattr(demo, "generate_browser_report", fail_browser)
    output_dir = tmp_path / "partial"
    with pytest.raises(demo.DemoError, match="browser failure"):
        demo.run_demo(output_dir)
    assert not output_dir.exists()


def test_run_demo_rejects_disagreeing_mtdna_counts(tmp_path, monkeypatch):
    monkeypatch.setattr(demo, "generate_browser_report", lambda *args, **kwargs: {"variants": 0})
    output_dir = tmp_path / "count-mismatch"
    with pytest.raises(demo.DemoError, match="variant counts differ"):
        demo.run_demo(output_dir)
    assert not output_dir.exists()
