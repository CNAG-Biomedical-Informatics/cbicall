from cbicall import report_utils as ru


def _report(**overrides):
    base = {
        "_report_path": "/tmp/run_a/run-report.json",
        "framework": {"version": "1.0"},
        "runtime": {"python": {"version": "3.12"}, "backend": {"version": "bash"}},
        "workflow": {
            "key": "bash/wes/single/gatk-4.6/v1",
            "files": [{"role": "entrypoint", "path": "/tmp/wes.sh", "sha256": "a" * 64}, "bad"],
        },
        "execution_contract": {"generated_files": [{"role": "params", "normalized_sha256": "b" * 64}, {}]},
        "outputs": {
            "file_inventory": {"entries": 1, "total_bytes": 1536, "sha256": "manifest"},
            "vcf_hash_reports": [{"file": "/tmp/sample.vcf.gz", "normalized_sha256": "c" * 64, "call_sha256": "e" * 64}, "bad"],
        },
        "resources": {"bundle": {"key": "bundle", "version": "v1", "fingerprint": "d" * 64}},
    }
    base.update(overrides)
    return base


def _report_with_inventory(total_bytes=None, sha256="manifest"):
    report = _report()
    report["outputs"]["file_inventory"] = {"sha256": sha256}
    if total_bytes is not None:
        report["outputs"]["file_inventory"]["total_bytes"] = total_bytes
    return report


def test_report_utils_scalar_helpers():
    assert ru._nested({"a": 1}, "a", "b") is None
    assert ru._nested("bad", "a") is None
    assert ru._same_text("x", "x") == "same"
    assert ru._same_text("x", "y") == "different"
    assert ru._run_label({"_report_alias": " alias ", "_report_path": "/tmp/a/run-report.json"}) == "alias"
    assert ru._run_label({"_report_path": "/tmp/a/run-report.json"}).endswith("run-report.json")
    assert ru._short_hash(None) is None
    assert ru._short_hash("not-a-hex-hash-but-long-enough-to-skip") == "not-a-hex-hash-but-long-enough-to-skip"
    assert ru._short_hash("a" * 64) == "aaaaaaaaaaaa...aaaaaaaa"
    assert ru._format_bytes("bad") == "bad"
    assert ru._format_bytes(-1) == "-1"
    assert ru._format_bytes(9) == "9 B (9 bytes)"
    assert ru._format_bytes(15 * 1024) == "15.0 KiB (15360 bytes)"
    assert ru._format_bytes(150 * 1024) == "150 KiB (153600 bytes)"
    assert ru._format_optional_bytes(None) is None
    assert ru._comparison_value(None) == "(missing)"
    assert ru._comparison_value([None, "x", "/tmp/a/b.txt"]) == "x, /tmp/a/b.txt"


def test_report_utils_maps_statuses_and_sections():
    report = _report()
    assert ru._workflow_file_map(report)["entrypoint"]["sha256"] == "a" * 64
    assert ru._execution_file_map(report)["params"]["normalized_sha256"] == "b" * 64
    assert ru._execution_file_value("missing", report) is None
    assert ru._vcf_hash_map(report)["sample.vcf.gz"]["normalized_sha256"] == "c" * 64
    assert ru._multi_workflow_file_value("entrypoint", report) == "a" * 64
    assert ru._multi_workflow_file_value("missing", report) is None
    assert ru._multi_vcf_hash_value("missing.vcf.gz", report) is None
    assert ru._multi_vcf_call_value("sample.vcf.gz", report) == "e" * 64

    same_hash_other_path = _report(
        workflow={
            "key": "bash/wes/single/gatk-4.6/v1",
            "files": [{"role": "entrypoint", "path": "/different/wes.sh", "sha256": "a" * 64}],
        }
    )
    workflow_row = next(row for section in ru._comparison_specs([report, same_hash_other_path]) if section["section"] == "Workflow Files" for row in section["rows"] if row["label"] == "entrypoint")
    assert ru._row_pair_status(workflow_row, report, same_hash_other_path)[0] == "same"

    assert ru._status_for_values(None, None) == ("unavailable", "not available", "")
    assert ru._status_for_values(None, "x")[0] == "missing"
    assert ru._status_for_values("x", "x") == ("same", "same", "x")
    assert ru._status_for_values("x", "y")[0] == "different"

    same_size = _report_with_inventory(total_bytes=1536)
    note_size = _report_with_inventory(total_bytes=2048)
    diff_size = _report_with_inventory(total_bytes=2048, sha256="other")
    missing_size = _report_with_inventory()
    assert ru._inventory_size_status(report, same_size)[0] == "same"
    assert ru._inventory_size_status(report, note_size)[0] == "note"
    assert ru._inventory_size_status(report, diff_size)[0] == "different"
    assert ru._inventory_size_status(report, missing_size)[0] == "missing"
    assert ru._inventory_size_status(missing_size, missing_size)[0] == "unavailable"

    assert ru._row_pair_status({"kind": "inventory_size"}, report, note_size)[0] == "note"
    assert ru._aggregate_status(["unavailable"]) == "unavailable"
    assert ru._aggregate_status(["same", "note"]) == "note"
    assert ru._aggregate_status(["same", "missing"]) == "missing"
    assert ru._aggregate_status(["same", "different"]) == "different"
    assert ru._aggregate_status(["same", "same"]) == "same"

    sections = ru._comparison_sections_with_overall([report, note_size])
    names = [section["section"] for section in sections]
    assert names[0] == "Overall"
    assert "Final VCF" in names
    assert any(row["label"] == "Inventory size" for row in sections[0]["rows"])
    final_vcf = next(section for section in sections if section["section"] == "Final VCF")
    assert any(row["label"] == "sample.vcf.gz calls" for row in final_vcf["rows"])
    labels = [row["label"] for row in final_vcf["rows"]]
    assert labels.index("sample.vcf.gz calls") < labels.index("sample.vcf.gz strict records")
