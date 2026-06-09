from cbicall.audit_similarity import build_audit_similarity, qualitative_similarity_label


def _report(alias="run", workflow_hash="workflow", vcf_hash="vcf", backend="1.0", workflow_path="/tmp/a/wes.sh"):
    return {
        "_report_alias": alias,
        "framework": {"version": "1.0.1"},
        "runtime": {
            "python": {"version": "3.10.8"},
            "java": {"version": "17.0.10"},
            "configured_java": {"version": "17.0.10"},
            "backend": {"version": backend},
        },
        "workflow": {
            "key": "bash/wes/single/gatk-4.6/v1",
            "registry_version": "v1",
            "entrypoint": workflow_path,
            "fingerprint": workflow_hash,
            "files": [{"role": "entrypoint", "path": workflow_path, "sha256": workflow_hash}],
        },
        "execution_contract": {
            "fingerprint": "contract",
            "normalized_command_sha256": "command",
            "generated_files": [],
        },
        "software_versions": {"sha256": "software"},
        "resources": {"bundle": {"key": "bundle", "version": "v1", "fingerprint": "resource"}},
        "execution_trace": {"tasks": 2, "max_peak_rss": {"bytes": 1000}, "max_peak_vmem": {"bytes": 2000}},
        "outputs": {
            "file_inventory": {"entries": 3, "total_bytes": 4096, "sha256": "inventory"},
            "vcf_hash_reports": [{"file": "/tmp/a/sample.vcf.gz", "normalized_sha256": vcf_hash}],
        },
    }


def _layer(similarity, name):
    return next(layer for layer in similarity["layers"] if layer["name"] == name)


def test_audit_similarity_identical_reports_score_one():
    similarity = build_audit_similarity([_report("local"), _report("cloud")])

    assert similarity["runs"] == ["local", "cloud"]
    assert _layer(similarity, "Overall")["rows"][0][1]["score"] == 1.0
    assert _layer(similarity, "Final VCF")["rows"][0][1]["status"] == "same"


def test_audit_similarity_reports_layer_specific_differences():
    similarity = build_audit_similarity([_report("local"), _report("hpc", backend="2.0", vcf_hash="other")])

    assert _layer(similarity, "Overall")["rows"][0][1]["score"] < 1.0
    assert _layer(similarity, "Framework")["rows"][0][1]["status"] in {"medium", "low"}
    assert _layer(similarity, "Final VCF")["rows"][0][1]["status"] == "low"


def test_audit_similarity_unavailable_when_layer_has_no_information():
    left = _report("left")
    right = _report("right")
    left.pop("execution_trace")
    right.pop("execution_trace")

    similarity = build_audit_similarity([left, right])

    assert _layer(similarity, "Execution")["rows"][0][1]["status"] == "unavailable"


def test_audit_similarity_workflow_file_hash_ignores_path_noise():
    left = _report("left", workflow_hash="a" * 64, workflow_path="/tmp/left/wes.sh")
    right = _report("right", workflow_hash="a" * 64, workflow_path="/tmp/right/wes.sh")

    similarity = build_audit_similarity([left, right])

    assert _layer(similarity, "Workflow Files")["rows"][0][1]["score"] == 1.0


def test_qualitative_similarity_label_categories():
    assert qualitative_similarity_label("same", 1.0) == "same"
    assert qualitative_similarity_label("different", 0.98) == "near"
    assert qualitative_similarity_label("different", 0.85) == "partial"
    assert qualitative_similarity_label("different", 0.60) == "diverged"
    assert qualitative_similarity_label("missing", 0.95) == "missing"
    assert qualitative_similarity_label("unavailable", None) == "n/a"
