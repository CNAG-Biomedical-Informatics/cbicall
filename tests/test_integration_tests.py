import gzip
import hashlib
import json
from types import SimpleNamespace
from pathlib import Path

import pytest

from cbicall import integration_tests as integration_mod
from cbicall.integration_tests import IntegrationTestError, TESTS, validate_contract


def _write_report(run_dir: Path, *, status="success", backend="bash") -> None:
    (run_dir / "run-report.json").write_text(
        json.dumps(
            {
                "status": status,
                "workflow": {"backend": backend, "pipeline": "wes"},
                "resources": {"bundle": {"key": "cbicall-germline-resources-v1"}},
            }
        ),
        encoding="utf-8",
    )


def _normalized_hash(lines):
    payload = ("\n".join(sorted(lines)) + "\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def test_validate_contract_accepts_required_files_json_and_normalized_gzip_hash(tmp_path):
    _write_report(tmp_path)
    (tmp_path / "workflow.log").write_text("ok\n", encoding="utf-8")
    with gzip.open(tmp_path / "sample.vcf.gz", "wt", encoding="utf-8") as handle:
        handle.write("##header\n")
        handle.write("#CHROM\tPOS\n")
        handle.write("2\t200\n")
        handle.write("1\t100\n")

    contract = {
        "required_files": ["run-report.json", "workflow.log", "sample.vcf.gz"],
        "json_expectations": [
            {"path": "status", "equals": "success"},
            {"path": "workflow.backend", "equals": "bash"},
            {"path": "resources.bundle.key", "equals": "cbicall-germline-resources-v1"},
        ],
        "hashes": [
            {
                "name": "sample",
                "type": "normalized_vcf",
                "path": "sample.vcf.gz",
                "pattern": "^#",
                "records": 2,
                "sha256": _normalized_hash(["1\t100", "2\t200"]),
            }
        ],
    }

    validate_contract(tmp_path, contract)


def test_validate_contract_reports_missing_glob(tmp_path):
    _write_report(tmp_path)
    contract = {"required_globs": ["pipeline_info/params_*.json"]}

    with pytest.raises(IntegrationTestError, match="missing glob"):
        validate_contract(tmp_path, contract)


def test_validate_contract_reports_json_mismatch(tmp_path):
    _write_report(tmp_path, status="failed")
    contract = {"json_expectations": [{"path": "status", "equals": "success"}]}

    with pytest.raises(IntegrationTestError, match="JSON expectation failed"):
        validate_contract(tmp_path, contract)


def test_validate_contract_reports_hash_mismatch(tmp_path):
    _write_report(tmp_path)
    (tmp_path / "out.txt").write_text("actual\n", encoding="utf-8")
    contract = {
        "hashes": [
            {
                "type": "normalized_text",
                "path": "out.txt",
                "sha256": _normalized_hash(["expected"]),
            }
        ]
    }

    with pytest.raises(IntegrationTestError, match="Hash mismatch"):
        validate_contract(tmp_path, contract)


def test_validate_contract_accepts_canonical_json_hash(tmp_path):
    _write_report(tmp_path)
    payload = {"b": [2, 1], "a": {"z": "value"}}
    (tmp_path / "out.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    contract = {
        "hashes": [
            {
                "type": "canonical_json",
                "path": "out.json",
                "sha256": hashlib.sha256(canonical).hexdigest(),
            }
        ]
    }

    validate_contract(tmp_path, contract)


def _write_fixture(project_root: Path, name: str, payload: dict) -> None:
    fixture_dir = project_root / "tests" / "fixtures" / "integration"
    fixture_dir.mkdir(parents=True)
    import yaml

    (fixture_dir / name).write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def test_run_integration_tests_executes_and_cleans_contract_run(tmp_path, monkeypatch, capsys):
    workdir = tmp_path / "examples" / "input"
    workdir.mkdir(parents=True)
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "cbicall").write_text("#!/bin/sh\n", encoding="utf-8")
    (workdir / "param.yaml").write_text("mode: single\n", encoding="utf-8")
    _write_fixture(
        tmp_path,
        "native-wes-bash.yaml",
        {
            "workdir": "examples/input",
            "workflow_log": "workflow.log",
            "run": {
                "parameter_file": "param.yaml",
                "base_dir": ".",
                "run_glob": "cbicall_bash_test_*",
            },
            "required_files": ["run-report.json", "run-report.html", "workflow.log"],
            "json_expectations": [{"path": "status", "equals": "success"}],
            "cleanup_paths": ["work"],
        },
    )

    def fake_run(cmd, cwd, stdout, stderr, text, check):
        run_dir = workdir / "cbicall_bash_test_001"
        run_dir.mkdir()
        (run_dir / "run-report.json").write_text(json.dumps({"status": "success"}), encoding="utf-8")
        (run_dir / "run-report.html").write_text("<html></html>\n", encoding="utf-8")
        (run_dir / "workflow.log").write_text("ok\n", encoding="utf-8")
        (run_dir / "work").mkdir()
        return SimpleNamespace(returncode=0, stdout=f"  Report     => {run_dir / 'run-report.json'}\n")

    monkeypatch.setattr(integration_mod.subprocess, "run", fake_run)

    rc = integration_mod.run_integration_tests(
        project_root=tmp_path,
        selected=[TESTS["wes-bash"]],
        threads=1,
        runtime_profile="local",
        keep_external_work=False,
    )

    assert rc == 0
    assert not (workdir / "cbicall_bash_test_001" / "work").exists()
    out = capsys.readouterr().out
    assert "SUCCESS: WES Bash integration contract passed." in out


def test_run_integration_tests_keeps_external_work_when_requested(tmp_path, monkeypatch):
    workdir = tmp_path / "examples" / "input"
    workdir.mkdir(parents=True)
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "cbicall").write_text("#!/bin/sh\n", encoding="utf-8")
    _write_fixture(
        tmp_path,
        "nf-core-demo.yaml",
        {
            "workdir": "examples/input",
            "workflow_log": "nf-core_demo_single.log",
            "run": {
                "parameter_file": "nf-core-demo.yaml",
                "base_dir": ".",
                "run_glob": "cbicall_nextflow_nf-core_demo_single_no-genome_*",
            },
            "required_files": ["run-report.json", "nf-core_demo_single.log"],
            "json_expectations": [{"path": "workflow.metadata.provider", "equals": "nf-core"}],
            "cleanup_paths": ["work"],
        },
    )
    (workdir / "nf-core-demo.yaml").write_text("workflow_provider: nf-core\n", encoding="utf-8")

    def fake_run(cmd, cwd, stdout, stderr, text, check):
        run_dir = workdir / "cbicall_nextflow_nf-core_demo_single_no-genome_001"
        run_dir.mkdir()
        (run_dir / "run-report.json").write_text(
            json.dumps({"workflow": {"metadata": {"provider": "nf-core"}}}),
            encoding="utf-8",
        )
        (run_dir / "nf-core_demo_single.log").write_text("ok\n", encoding="utf-8")
        (run_dir / "work").mkdir()
        return SimpleNamespace(returncode=0, stdout=f"Working directory: {run_dir}\n")

    monkeypatch.setattr(integration_mod.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(integration_mod.subprocess, "run", fake_run)

    rc = integration_mod.run_integration_tests(
        project_root=tmp_path,
        selected=[TESTS["nf-core-demo"]],
        threads=2,
        runtime_profile="local",
        keep_external_work=True,
    )

    assert rc == 0
    assert (workdir / "cbicall_nextflow_nf-core_demo_single_no-genome_001" / "work").is_dir()


def test_run_integration_tests_skips_missing_optional_backend_under_all(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(integration_mod.shutil, "which", lambda name: None)

    rc = integration_mod.run_integration_tests(
        project_root=tmp_path,
        selected=[TESTS["wes-snakemake"]],
        threads=1,
        runtime_profile="local",
        skip_missing_optional=True,
    )

    assert rc == 0
    assert "skipped" in capsys.readouterr().out


def test_run_integration_tests_fails_missing_explicit_backend(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(integration_mod.shutil, "which", lambda name: None)

    rc = integration_mod.run_integration_tests(
        project_root=tmp_path,
        selected=[TESTS["wes-snakemake"]],
        threads=1,
        runtime_profile="local",
        skip_missing_optional=False,
    )

    assert rc == 1
    assert "requires backend executable snakemake" in capsys.readouterr().out


def _write_release_report(run_dir: Path, *, backend="bash", vcf_hash="same-hash", records="6") -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "run-report.json").write_text(
        json.dumps(
            {
                "status": "success",
                "workflow": {
                    "backend": backend,
                    "pipeline": "wes",
                    "mode": "single",
                    "software_stack": "gatk-4.6",
                },
                "run": {"display_genome": "b37"},
                "resources": {"bundle": {"key": "cbicall-germline-resources-v1"}},
                "outputs": {
                    "vcf_hash_reports": [
                        {
                            "file": "CNAG99901P.hc.QC.vcf.gz",
                            "normalized_sha256": vcf_hash,
                            "normalized_records": records,
                        }
                    ]
                },
            }
        ),
        encoding="utf-8",
    )


def test_run_release_equivalence_compares_available_native_backend(tmp_path, monkeypatch, capsys):
    bash_run = tmp_path / "bash_run"
    snakemake_run = tmp_path / "snakemake_run"
    _write_release_report(bash_run, backend="bash")
    _write_release_report(snakemake_run, backend="snakemake")

    def fake_run_one(**kwargs):
        selection = kwargs["selection"]
        if selection.key == "wes-bash":
            return selection.label, "passed", str(bash_run)
        if selection.key == "wes-snakemake":
            return selection.label, "passed", str(snakemake_run)
        raise AssertionError(f"unexpected run: {selection.key}")

    monkeypatch.setattr(integration_mod, "_run_one", fake_run_one)
    monkeypatch.setattr(integration_mod, "load_contract", lambda project_root, selection: {})
    monkeypatch.setattr(
        integration_mod,
        "_backend_is_available",
        lambda selection: (selection.key == "wes-snakemake", "missing test backend"),
    )

    rc = integration_mod.run_release_equivalence_test(
        project_root=tmp_path,
        threads=2,
        runtime_profile="cnag-hpc",
    )

    assert rc == 0
    out = capsys.readouterr().out
    assert "Baseline" in out
    assert "Backend equivalence" in out
    assert "WES Snakemake => same final VCF | same-hash | 6 records" in out
    assert "WES Bash      => passed | same-hash | 6 records" in out
    assert "Compared non-Bash backends: 1" in out
    assert "Status: PASSED" in out


def test_run_release_equivalence_fails_without_non_bash_comparator(tmp_path, monkeypatch, capsys):
    bash_run = tmp_path / "bash_run"
    _write_release_report(bash_run, backend="bash")

    def fake_run_one(**kwargs):
        selection = kwargs["selection"]
        assert selection.key == "wes-bash"
        return selection.label, "passed", str(bash_run)

    monkeypatch.setattr(integration_mod, "_run_one", fake_run_one)
    monkeypatch.setattr(integration_mod, "load_contract", lambda project_root, selection: {})
    monkeypatch.setattr(integration_mod, "_backend_is_available", lambda selection: (False, "missing test backend"))

    rc = integration_mod.run_release_equivalence_test(
        project_root=tmp_path,
        threads=1,
        runtime_profile="local",
    )

    assert rc == 1
    assert "no non-Bash backend was available" in capsys.readouterr().out


def test_run_release_equivalence_fails_on_vcf_hash_mismatch(tmp_path, monkeypatch, capsys):
    bash_run = tmp_path / "bash_run"
    nextflow_run = tmp_path / "nextflow_run"
    _write_release_report(bash_run, backend="bash", vcf_hash="baseline")
    _write_release_report(nextflow_run, backend="nextflow", vcf_hash="different")

    def fake_run_one(**kwargs):
        selection = kwargs["selection"]
        if selection.key == "wes-bash":
            return selection.label, "passed", str(bash_run)
        if selection.key == "wes-nextflow":
            return selection.label, "passed", str(nextflow_run)
        raise AssertionError(f"unexpected run: {selection.key}")

    monkeypatch.setattr(integration_mod, "_run_one", fake_run_one)
    monkeypatch.setattr(integration_mod, "load_contract", lambda project_root, selection: {})
    monkeypatch.setattr(
        integration_mod,
        "_backend_is_available",
        lambda selection: (selection.key == "wes-nextflow", "missing test backend"),
    )

    rc = integration_mod.run_release_equivalence_test(
        project_root=tmp_path,
        threads=1,
        runtime_profile="local",
    )

    assert rc == 1
    out = capsys.readouterr().out
    assert "normalized VCF hash differs from Bash baseline" in out
