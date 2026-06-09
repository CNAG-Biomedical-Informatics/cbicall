from pathlib import Path

import pytest
import yaml

from cbicall import multiqc


def test_multiqc_helper_fallbacks_and_invalid_values(tmp_path):
    assert multiqc._short_hash(None) is None
    assert multiqc._short_hash("short") == "short"
    assert multiqc._run_label({"_report_alias": " alias "}) == "alias"
    assert multiqc._run_label({"run": {"run_id": "RID"}, "workflow": {"backend": "bash", "pipeline": "wes"}}) == "bash_wes_RID"
    assert multiqc._run_label({}, tmp_path / "run_a" / "run-report.json") == "run_a"
    assert multiqc._run_label({"_report_path": str(tmp_path / "run_b" / "run-report.json")}) == "run_b"
    assert multiqc._run_label({}) == "cbicall_run"

    assert multiqc._provider({}) == "cbicall"
    assert multiqc._provider({"workflow": {"metadata": {"provider": "nf-core"}}}) == "nf-core"
    assert multiqc._provider({"execution_contract": {"provider": "external"}}) == "external"

    assert multiqc._as_float(None) is None
    assert multiqc._as_float("bad") is None
    assert multiqc._as_int(None) is None
    assert multiqc._as_int("bad") is None
    assert multiqc._bytes_to_gib("bad") is None
    assert multiqc._bytes_to_mib("bad") is None
    assert multiqc._elapsed_minutes({"elapsed_seconds": "bad"}) is None
    assert multiqc._first_vcf_hash({"outputs": {"vcf_hash_reports": [{"sha256": "raw"}]}}) == "raw"
    assert multiqc._first_vcf_hash({"outputs": {"vcf_hash_reports": ["bad"]}}) is None
    assert multiqc._first_vcf_records({"outputs": {"vcf_hash_reports": [{"normalized_records": "bad"}]}}) is None


def test_multiqc_native_qc_parsers_ignore_unusable_files(tmp_path):
    stats_dir = tmp_path / "03_stats"
    stats_dir.mkdir()
    (stats_dir / "bad.coverage.txt").write_bytes(b"\xff\xfe")
    (stats_dir / "no_header.coverage.txt").write_text("sample\nvalue\n", encoding="utf-8")
    (stats_dir / "short.coverage.txt").write_text("sampleID\tmode\nonlyonefield\n", encoding="utf-8")
    (stats_dir / "bad.sex.txt").write_bytes(b"\xff\xfe")
    (stats_dir / "note.sex.txt").write_text("no equals here\n", encoding="utf-8")

    payload = multiqc.build_native_sample_qc_payload(tmp_path / "run-report.json")
    assert payload is not None
    assert payload["data"] == {"bad": {}, "note": {}}
    assert multiqc.build_native_sample_qc_payload(tmp_path / "missing" / "run-report.json") is None


def test_multiqc_final_outputs_and_bundle_rewrite(tmp_path):
    report_path = tmp_path / "run-report.json"
    payload = {
        "outputs": {
            "vcf_hash_reports": [],
            "canonical_outputs": [
                "ignore-me",
                {"name": "final_vcf", "type": "vcf", "status": "missing", "pattern": "*.vcf.gz", "matches": []},
            ],
        }
    }
    final_outputs = multiqc.build_final_outputs_payload(report_path, payload)
    assert final_outputs is not None
    assert final_outputs["data"]["final_vcf"]["Status"] == "missing"
    assert final_outputs["data"]["final_vcf"]["Matches"] == 0
    assert multiqc.build_final_outputs_payload(report_path, {"outputs": {"canonical_outputs": ["bad"]}}) is None

    outdir = tmp_path / "cbicall_mqc"
    outdir.mkdir()
    stale = outdir / "stale_mqc.yaml"
    stale.write_text("stale: true\n", encoding="utf-8")
    keep = outdir / "keep.txt"
    keep.write_text("keep\n", encoding="utf-8")
    run_payload = {
        "status": "success",
        "elapsed_seconds": 30,
        "workflow": {"backend": "bash", "pipeline": "wes"},
        "run": {"run_id": "RID", "threads": "2"},
        "outputs": {"file_inventory": {"entries": "1", "total_bytes": 1024}, "vcf_hash_reports": []},
    }
    written = multiqc.write_multiqc_report(report_path, run_payload, outdir)
    assert written == outdir
    assert not stale.exists()
    assert keep.exists()
    assert (outdir / "cbicall_run_general_stats_mqc.yaml").is_file()
    data = yaml.safe_load((outdir / "cbicall_run_general_stats_mqc.yaml").read_text(encoding="utf-8"))
    assert data["data"]["bash_wes_RID"]["threads"] == 2

    compare_payloads = [dict(run_payload, _report_alias="local"), dict(run_payload, _report_alias="cloud")]
    pair_summary = multiqc.build_compare_pair_summary_payload(compare_payloads)
    assert pair_summary["data"]["local vs cloud"]["Overall category"] == "same"
    assert pair_summary["data"]["local vs cloud"]["Overall similarity"] == 1.0
    heatmap = multiqc.build_compare_similarity_heatmap_payload(compare_payloads)
    assert heatmap["plot_type"] == "heatmap"
    assert heatmap["pconfig"]["ycats_samples"] == ["local", "cloud"]
    assert heatmap["data"][0][1] == 1.0

    with pytest.raises(ValueError, match="directory bundle"):
        multiqc.write_compare_multiqc_report([run_payload, run_payload], tmp_path / "old.yaml")
