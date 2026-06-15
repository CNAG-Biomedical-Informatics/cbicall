import gzip
import hashlib
import json
from types import SimpleNamespace
from pathlib import Path

import pytest
import yaml

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
    (tmp_path / "sample.coverage.txt").write_text("region\tsampleID\nchr22\tCNAG99901P\n", encoding="utf-8")
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
        "text_expectations": [
            {"path": "sample.coverage.txt", "startswith": "region\tsampleID", "contains": "CNAG99901P", "matches": "^(chr)?22\t"},
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


def test_write_setup_files_from_runs_creates_fake_sample_map(tmp_path):
    workdir = tmp_path / "work"
    run_dir = tmp_path / "single-run"
    gvcf = run_dir / "02_varcall" / "sample.hc.g.vcf.gz"
    gvcf.parent.mkdir(parents=True)
    workdir.mkdir()
    gvcf.write_text("gvcf\n", encoding="utf-8")
    (gvcf.parent / "sample.hc.g.vcf.gz.tbi").write_text("index\n", encoding="utf-8")

    contract = {
        "setup_files_from_runs": [
            {
                "run": "wes-bash",
                "source": "02_varcall/sample.hc.g.vcf.gz",
                "source_index": "02_varcall/sample.hc.g.vcf.gz.tbi",
                "sample_map": "cohort-map.tsv",
                "samples": ["FAKE_01", "FAKE_02"],
            }
        ]
    }

    paths = integration_mod._write_setup_files_from_runs(workdir, contract, {"wes-bash": run_dir})

    assert paths == [workdir / "cohort-map.tsv"]
    assert (workdir / "cohort-map.tsv").read_text(encoding="utf-8") == (
        f"FAKE_01\t{gvcf}\n"
        f"FAKE_02\t{gvcf}\n"
    )


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
            "notes": ["This test prints a setup note."],
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
    assert "Note: This test prints a setup note." in out


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



def test_contract_loader_path_helpers_and_selection_errors(tmp_path):
    selection = integration_mod.TestSelection("demo", "Demo", "missing.yaml")
    with pytest.raises(FileNotFoundError):
        integration_mod.load_contract(tmp_path, selection)

    fixture = tmp_path / "tests" / "fixtures" / "integration"
    fixture.mkdir(parents=True)
    (fixture / "bad.yaml").write_text("- not\n- a mapping\n", encoding="utf-8")
    with pytest.raises(IntegrationTestError, match="must be a mapping"):
        integration_mod.load_contract(tmp_path, integration_mod.TestSelection("bad", "Bad", "bad.yaml"))

    assert integration_mod._project_path(tmp_path, "relative") == tmp_path / "relative"
    assert integration_mod._project_path(tmp_path, str(tmp_path / "absolute")) == tmp_path / "absolute"
    workdir = tmp_path / "work"
    assert integration_mod._work_path(workdir, "run") == workdir / "run"
    assert integration_mod._work_path(workdir, str(tmp_path / "run")) == tmp_path / "run"
    assert integration_mod.list_run_dirs(tmp_path / "missing", "cbicall_*") == []
    assert integration_mod._parse_run_dir_from_stdout("nothing useful\n") is None
    assert integration_mod._new_run_dir([], []) is None


def test_validate_contract_suffix_contains_and_hash_error_paths(tmp_path):
    _write_report(tmp_path)
    report = json.loads((tmp_path / "run-report.json").read_text(encoding="utf-8"))
    report["links"] = {"report": "sample-report.html", "message": "contains expected token"}
    (tmp_path / "run-report.json").write_text(json.dumps(report), encoding="utf-8")
    (tmp_path / "out.txt").write_text("A\nB\n", encoding="utf-8")

    validate_contract(
        tmp_path,
        {
            "json_expectations": [
                {"path": "links.report", "endswith": "report.html"},
                {"path": "links.message", "contains": "expected"},
            ],
            "hashes": [{"type": "normalized_text", "path": "out.txt", "records": 2, "sha256": _normalized_hash(["A", "B"])}],
        },
    )

    with pytest.raises(IntegrationTestError, match="Missing JSON path"):
        validate_contract(tmp_path, {"json_expectations": [{"path": "links.missing", "equals": "x"}]})
    with pytest.raises(IntegrationTestError, match="missing file: run-report.json"):
        validate_contract(tmp_path / "no-report", {"json_expectations": [{"path": "status", "equals": "success"}]})
    with pytest.raises(IntegrationTestError, match="Hash target does not exist"):
        validate_contract(tmp_path, {"hashes": [{"path": "missing.txt", "sha256": "x"}]})
    with pytest.raises(IntegrationTestError, match="Record count mismatch"):
        validate_contract(tmp_path, {"hashes": [{"path": "out.txt", "records": 3, "sha256": "x"}]})
    with pytest.raises(IntegrationTestError, match="Unsupported hash type"):
        validate_contract(tmp_path, {"hashes": [{"path": "out.txt", "type": "md5", "sha256": "x"}]})
    with pytest.raises(IntegrationTestError, match="Text expectation target does not exist"):
        validate_contract(tmp_path, {"text_expectations": [{"path": "missing.txt", "contains": "x"}]})
    with pytest.raises(IntegrationTestError, match="Text expectation failed"):
        validate_contract(tmp_path, {"text_expectations": [{"path": "out.txt", "contains": "missing"}]})
    with pytest.raises(IntegrationTestError, match="Text expectation failed"):
        validate_contract(tmp_path, {"text_expectations": [{"path": "out.txt", "startswith": "B"}]})
    with pytest.raises(IntegrationTestError, match="Text expectation failed"):
        validate_contract(tmp_path, {"text_expectations": [{"path": "out.txt", "endswith": "C"}]})
    with pytest.raises(IntegrationTestError, match="Text expectation failed"):
        validate_contract(tmp_path, {"text_expectations": [{"path": "out.txt", "matches": "^missing"}]})


def test_run_one_handles_inline_parameters_and_missing_run_dir(tmp_path, monkeypatch, capsys):
    workdir = tmp_path / "examples" / "input"
    workdir.mkdir(parents=True)
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "cbicall").write_text("#!/bin/sh\n", encoding="utf-8")
    _write_fixture(
        tmp_path,
        "inline.yaml",
        {
            "workdir": "examples/input",
            "workflow_log": "workflow.log",
            "run": {"parameters": {"mode": "single"}, "base_dir": ".", "run_glob": "cbicall_inline_*"},
            "required_files": ["run-report.json", "workflow.log"],
            "cleanup_paths": ["scratch.tmp"],
        },
    )
    selection = integration_mod.TestSelection("inline", "Inline", "inline.yaml")

    def fake_run(cmd, cwd, stdout, stderr, text, check):
        param_file = Path(cmd[cmd.index("-p") + 1])
        assert param_file.name.startswith("cbicall-inline.")
        run_dir = workdir / "cbicall_inline_001"
        run_dir.mkdir()
        (run_dir / "run-report.json").write_text(json.dumps({"status": "success"}), encoding="utf-8")
        (run_dir / "workflow.log").write_text("ok\n", encoding="utf-8")
        (run_dir / "scratch.tmp").write_text("remove\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="Working directory: cbicall_inline_001\n")

    monkeypatch.setattr(integration_mod.subprocess, "run", fake_run)
    label, status, detail = integration_mod._run_one(
        project_root=tmp_path,
        selection=selection,
        threads=1,
        runtime_profile="local",
        skip_missing_optional=False,
        keep_external_work=False,
    )
    assert (label, status) == ("Inline", "passed")
    assert detail.endswith("cbicall_inline_001")
    assert not (workdir / "cbicall_inline_001" / "scratch.tmp").exists()
    assert not list(workdir.glob("cbicall-inline.*.yaml"))
    assert "Cleaning heavy execution state" in capsys.readouterr().out

    fixture_dir = tmp_path / "tests" / "fixtures" / "integration"
    (fixture_dir / "no-glob.yaml").write_text(yaml.safe_dump({"run": {"parameter_file": "param.yaml"}}, sort_keys=False), encoding="utf-8")
    with pytest.raises(IntegrationTestError, match="run_glob"):
        integration_mod._run_one(
            project_root=tmp_path,
            selection=integration_mod.TestSelection("noglob", "NoGlob", "no-glob.yaml"),
            threads=1,
            runtime_profile="local",
            skip_missing_optional=False,
            keep_external_work=False,
        )


def test_backend_availability_and_selected_tests(monkeypatch):
    monkeypatch.delenv("CROMWELL_JAR", raising=False)
    monkeypatch.setattr(integration_mod.shutil, "which", lambda name: None)
    available, detail = integration_mod._backend_is_available(integration_mod.TESTS["wes-cromwell"])
    assert available is False
    assert "CROMWELL_JAR" in detail

    monkeypatch.setenv("CROMWELL_JAR", "/opt/cromwell.jar")
    assert integration_mod._backend_is_available(integration_mod.TESTS["wes-cromwell"])[0] is True

    args = SimpleNamespace(
        all=False,
        wes_bash=True,
        wes_cohort_bash=False,
        wes_bash_gatk35=False,
        wes_snakemake=False,
        wes_nextflow=True,
        wes_cromwell=False,
        mit_bash=False,
        nf_core_demo=True,
        nf_core_sarek=True,
    )
    selected = integration_mod.selected_tests_from_args(args)
    assert [item.key for item in selected] == ["wes-bash", "wes-nextflow", "nf-core-demo", "nf-core-sarek"]

    args.wes_bash = False
    args.wes_cohort_bash = True
    args.wes_nextflow = False
    args.nf_core_demo = False
    args.nf_core_sarek = False
    selected = integration_mod.selected_tests_from_args(args)
    assert [item.key for item in selected] == ["wes-cohort-bash"]

    args.all = True
    args.wes_cohort_bash = False
    selected_all = integration_mod.selected_tests_from_args(args)
    assert "wes-cohort-bash" not in [item.key for item in selected_all]
    assert "wes-cromwell" in [item.key for item in selected_all]
    assert "mit-bash" in [item.key for item in selected_all]



def test_release_helpers_contract_fallback_and_metadata_errors(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run-report.json").write_text(json.dumps({"outputs": {"vcf_hash_reports": []}}), encoding="utf-8")
    with pytest.raises(IntegrationTestError, match="No normalized VCF hash"):
        integration_mod._release_vcf_hash(run_dir, {})

    (run_dir / "sample.vcf").write_text("##h\n#CHROM\tPOS\n1\t10\n", encoding="utf-8")
    contract = {"hashes": [{"type": "normalized_vcf", "path": "sample.vcf"}]}
    fallback = integration_mod._release_vcf_hash(run_dir, contract)
    assert fallback["source"] == "integration_contract"
    assert fallback["normalized_records"] == "1"

    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    _write_release_report(baseline, backend="bash")
    _write_release_report(candidate, backend="snakemake")
    data = json.loads((candidate / "run-report.json").read_text(encoding="utf-8"))
    data["workflow"]["pipeline"] = "wgs"
    (candidate / "run-report.json").write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(IntegrationTestError, match="Release metadata mismatch"):
        integration_mod._check_release_metadata(baseline, candidate)


def test_run_release_equivalence_failure_branches(tmp_path, monkeypatch, capsys):
    def fail_baseline(**kwargs):
        raise IntegrationTestError("baseline exploded")

    monkeypatch.setattr(integration_mod, "_run_one", fail_baseline)
    rc = integration_mod.run_release_equivalence_test(project_root=tmp_path, threads=1, runtime_profile="local")
    assert rc == 1
    assert "baseline exploded" in capsys.readouterr().out

    bash_run = tmp_path / "bash_run"
    cromwell_run = tmp_path / "cromwell_run"
    _write_release_report(bash_run, backend="bash", records="6")
    _write_release_report(cromwell_run, backend="cromwell", records="7")

    def fake_run_one(**kwargs):
        selection = kwargs["selection"]
        if selection.key == "wes-bash":
            return selection.label, "passed", str(bash_run)
        if selection.key == "wes-cromwell":
            return selection.label, "passed", str(cromwell_run)
        raise AssertionError(selection.key)

    monkeypatch.setattr(integration_mod, "_run_one", fake_run_one)
    monkeypatch.setattr(integration_mod, "load_contract", lambda project_root, selection: {})
    monkeypatch.setattr(integration_mod, "_backend_is_available", lambda selection: (selection.key == "wes-cromwell", "missing"))
    rc = integration_mod.run_release_equivalence_test(project_root=tmp_path, threads=1, runtime_profile="local")
    assert rc == 1
    out = capsys.readouterr().out
    assert "record count differs" in out


def test_run_integration_tests_records_failed_summary(tmp_path, monkeypatch, capsys):
    def fail_one(**kwargs):
        raise IntegrationTestError("contract failed\nwith detail")

    monkeypatch.setattr(integration_mod, "_run_one", fail_one)
    rc = integration_mod.run_integration_tests(
        project_root=tmp_path,
        selected=[integration_mod.TESTS["wes-bash"]],
        threads=1,
        runtime_profile="local",
    )
    assert rc == 1
    out = capsys.readouterr().out
    assert "ERROR: contract failed" in out
    assert "WES Bash: failed (contract failed)" in out
