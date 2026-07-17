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


def test_prepare_integration_root_stages_installed_assets(tmp_path):
    source = tmp_path / "runtime"
    (source / "examples" / "input").mkdir(parents=True)
    (source / "examples" / "input" / "param.yaml").write_text("mode: single\n", encoding="utf-8")

    staged, temporary = integration_mod.prepare_integration_root(source)
    assert temporary == staged
    assert staged != source
    assert (staged / "examples" / "input" / "param.yaml").is_file()
    assert integration_mod._cbicall_command(staged) == [integration_mod.sys.executable, "-m", "cbicall"]
    integration_mod.shutil.rmtree(staged)

    (source / "bin").mkdir()
    launcher = source / "bin" / "cbicall"
    launcher.write_text("#!/bin/sh\n", encoding="utf-8")
    direct, temporary = integration_mod.prepare_integration_root(source)
    assert direct == source.resolve()
    assert temporary is None
    assert integration_mod._cbicall_command(direct) == [str(launcher)]


def test_prepare_integration_root_uses_requested_workspace(tmp_path):
    source = tmp_path / "runtime"
    (source / "examples" / "input").mkdir(parents=True)
    (source / "examples" / "input" / "param.yaml").write_text("mode: single\n", encoding="utf-8")
    workspace = tmp_path / "scratch" / "cbicall-test"

    staged, retained = integration_mod.prepare_integration_root(source, workspace=workspace)

    assert staged == workspace.resolve()
    assert retained == staged
    assert (staged / "examples" / "input" / "param.yaml").is_file()


def test_prepare_integration_root_rejects_unsafe_workspace(tmp_path):
    source = tmp_path / "runtime"
    source.mkdir()
    nonempty = tmp_path / "nonempty"
    nonempty.mkdir()
    (nonempty / "existing.txt").write_text("keep\n", encoding="utf-8")

    with pytest.raises(IntegrationTestError, match="must be empty"):
        integration_mod.prepare_integration_root(source, workspace=nonempty)

    with pytest.raises(IntegrationTestError, match="outside the CBIcall runtime root"):
        integration_mod.prepare_integration_root(source, workspace=source / "test-work")


def test_prepare_source_checkout_in_requested_workspace_copies_runtime_assets(tmp_path):
    source = tmp_path / "checkout"
    (source / "bin").mkdir(parents=True)
    (source / "bin" / "cbicall").write_text("#!/bin/sh\n", encoding="utf-8")
    (source / "examples" / "input").mkdir(parents=True)
    (source / "examples" / "input" / "param.yaml").write_text("mode: single\n", encoding="utf-8")
    generated = source / "examples" / "input" / "cbicall_bash_previous_run"
    generated.mkdir()
    (generated / "large.bam").write_text("old output\n", encoding="utf-8")
    (source / "docs-site").mkdir()
    workspace = tmp_path / "workspace"

    staged, retained = integration_mod.prepare_integration_root(source, workspace=workspace)

    assert retained == staged
    assert (staged / "bin" / "cbicall").is_file()
    assert (staged / "examples" / "input" / "param.yaml").is_file()
    assert not (staged / "examples" / "input" / generated.name).exists()
    assert not (staged / "docs-site").exists()


def test_run_integration_tests_executes_and_cleans_contract_run(tmp_path, monkeypatch, capsys):
    workdir = tmp_path / "examples" / "input"
    workdir.mkdir(parents=True)
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "cbicall").write_text("#!/bin/sh\n", encoding="utf-8")
    (tmp_path / "fixture.interval_list").write_text("22\t1\t100\t+\tchr22\n", encoding="utf-8")
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
            "env": {"CBICALL_INTERVAL_LIST": "fixture.interval_list"},
            "required_files": ["run-report.json", "run-report.html", "workflow.log"],
            "json_expectations": [{"path": "status", "equals": "success"}],
            "cleanup_paths": ["work"],
        },
    )

    def fake_run(cmd, cwd, env, stdout, stderr, text, check):
        assert env["CBICALL_INTERVAL_LIST"] == str(tmp_path / "fixture.interval_list")
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
    assert "[PASS] WES Bash integration contract passed." in out
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

    def fake_run(cmd, cwd, env, stdout, stderr, text, check):
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
    assert "[SKIP] WES Snakemake requires backend executable snakemake on PATH." in capsys.readouterr().out


def test_run_integration_tests_skips_architecture_limited_tests_under_all(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(integration_mod.platform, "machine", lambda: "aarch64")

    rc = integration_mod.run_integration_tests(
        project_root=tmp_path,
        selected=[TESTS["mit-bash"], TESTS["wes-cohort-bash"]],
        threads=1,
        runtime_profile="local",
        skip_missing_optional=True,
    )

    assert rc == 0
    out = capsys.readouterr().out
    assert "[SKIP] MIT Bash is not supported on this architecture" in out
    assert "[SKIP] WES Cohort Bash is not supported on this architecture" in out


def test_run_integration_tests_fails_architecture_limited_explicit_test(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(integration_mod.platform, "machine", lambda: "arm64")

    rc = integration_mod.run_integration_tests(
        project_root=tmp_path,
        selected=[TESTS["wes-cohort-bash-sharded"]],
        threads=1,
        runtime_profile="local",
        skip_missing_optional=False,
    )

    assert rc == 1
    assert "WES Cohort Bash Sharded is not supported on this architecture" in capsys.readouterr().out


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


def test_run_one_skips_or_fails_missing_cromwell_before_loading_contract(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv("CROMWELL_JAR", raising=False)
    monkeypatch.setattr(integration_mod.shutil, "which", lambda name: None)

    assert integration_mod._run_one(
        project_root=tmp_path,
        selection=TESTS["wes-cromwell"],
        threads=1,
        runtime_profile="local",
        skip_missing_optional=True,
        keep_external_work=False,
    ) == ("WES Cromwell", "skipped", "CROMWELL_JAR or cromwell executable not found")
    assert "[SKIP] WES Cromwell requires CROMWELL_JAR" in capsys.readouterr().out

    with pytest.raises(IntegrationTestError, match="WES Cromwell requires CROMWELL_JAR"):
        integration_mod._run_one(
            project_root=tmp_path,
            selection=TESTS["wes-cromwell"],
            threads=1,
            runtime_profile="local",
            skip_missing_optional=False,
            keep_external_work=False,
        )


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
    assert "Status: [PASS] PASSED" in out


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
    with pytest.raises(IntegrationTestError, match="expected suffix"):
        validate_contract(tmp_path, {"json_expectations": [{"path": "links.report", "endswith": "missing.html"}]})
    with pytest.raises(IntegrationTestError, match="expected to contain"):
        validate_contract(tmp_path, {"json_expectations": [{"path": "links.message", "contains": "absent"}]})
    with pytest.raises(IntegrationTestError, match="missing file: missing.txt"):
        validate_contract(tmp_path, {"required_files": ["missing.txt"]})
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

    def fake_run(cmd, cwd, env, stdout, stderr, text, check):
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
    available, detail = integration_mod._backend_is_available(integration_mod.TESTS["wes-snakemake"])
    assert available is False
    assert "snakemake" in detail

    monkeypatch.setenv("CROMWELL_JAR", "/opt/cromwell.jar")
    assert integration_mod._backend_is_available(integration_mod.TESTS["wes-cromwell"])[0] is True
    monkeypatch.setattr(integration_mod.shutil, "which", lambda name: f"/usr/bin/{name}")
    assert integration_mod._backend_is_available(integration_mod.TESTS["wes-bash"]) == (True, "")
    assert integration_mod._backend_is_available(integration_mod.TESTS["wes-nextflow"]) == (True, "")

    args = SimpleNamespace(
        all=False,
        wes_bash=True,
        wes_cohort_bash=False,
        wes_cohort_bash_sharded=False,
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

    args.wes_cohort_bash = False
    args.wes_cohort_bash_sharded = True
    selected = integration_mod.selected_tests_from_args(args)
    assert [item.key for item in selected] == ["wes-cohort-bash-sharded"]

    args.all = True
    args.wes_cohort_bash_sharded = False
    selected_all = integration_mod.selected_tests_from_args(args)
    assert [item.key for item in selected_all] == [
        "wes-bash",
        "mit-bash",
        "wes-cohort-bash",
        "wes-cohort-bash-sharded",
        "wes-snakemake",
        "wes-nextflow",
        "wes-cromwell",
    ]



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
    assert "[FAIL] contract failed" in out
    assert "[FAIL]" in out
    assert "WES Bash" in out
    assert "contract failed" in out



def test_resolve_bcftools_uses_explicit_env_path(tmp_path):
    bcftools = tmp_path / "bcftools"
    bcftools.write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    assert integration_mod._resolve_bcftools(tmp_path, {"BCFTOOLS": str(bcftools)}) == str(bcftools)


def test_resolve_bcftools_reports_missing_path_binary(tmp_path):
    with pytest.raises(IntegrationTestError, match="BCFTOOLS does not exist"):
        integration_mod._resolve_bcftools(tmp_path, {"BCFTOOLS": str(tmp_path / "missing-bcftools")})


def test_integration_helper_branches(tmp_path, monkeypatch):
    with pytest.raises(IntegrationTestError, match="env must be a mapping"):
        integration_mod._contract_env(tmp_path, {"env": ["not", "mapping"]})
    with pytest.raises(IntegrationTestError, match="env keys"):
        integration_mod._contract_env(tmp_path, {"env": {"": "value"}})

    existing = tmp_path / "relative-resource.txt"
    existing.write_text("ok\n", encoding="utf-8")
    monkeypatch.setenv("DROP_FROM_CONTRACT", "present")
    env = integration_mod._contract_env(
        tmp_path,
        {"env": {"DROP_FROM_CONTRACT": None, "RELATIVE_RESOURCE": "relative-resource.txt"}},
    )
    assert "DROP_FROM_CONTRACT" not in env
    assert env["RELATIVE_RESOURCE"] == str(existing)

    assert integration_mod._parse_run_dir_from_stdout(f"  Report   => {tmp_path / 'run' / 'run-report.json'}\n") == tmp_path / "run"
    assert integration_mod._parse_run_dir_from_stdout(f"  Log      => {tmp_path / 'run' / 'workflow.log'}\n") == tmp_path / "run"
    assert integration_mod._architecture_skip_reason(integration_mod.TESTS["wes-bash"], arch="x86_64") is None
    assert integration_mod._architecture_skip_reason(integration_mod.TESTS["wes-bash"], arch="arm64") is None
    assert integration_mod._short_hash(None) == "(undef)"
    assert integration_mod._short_hash("abc") == "abc"
    assert integration_mod._release_hash_detail(None) == "(hash not available)"
    assert integration_mod._release_hash_detail({"normalized_sha256": "a" * 64}) == "aaaaaaaaaaaa...aaaaaaaa"
    assert integration_mod._status_marker("skipped") == "[SKIP]"
    assert integration_mod._status_marker("failed") == "[FAIL]"
    assert integration_mod._status_marker("running") == "[INFO]"
    assert integration_mod._overall_marker(2) == "[FAIL]"

    with pytest.raises(IntegrationTestError, match="Missing run report"):
        integration_mod._load_report(tmp_path / "missing-run")
    assert integration_mod._vcf_hash_from_report({"outputs": {"vcf_hash_reports": "not-a-list"}}) is None
    assert integration_mod._vcf_hash_from_report({"outputs": {"vcf_hash_reports": [{"file": "x"}]}}) is None
    assert integration_mod._vcf_hash_from_contract(tmp_path, {"hashes": [{"type": "canonical_json", "path": "x.json"}]}) is None
    assert integration_mod._vcf_hash_from_contract(tmp_path, {"hashes": [{"type": "normalized_vcf", "path": "missing.vcf"}]}) is None
    assert integration_mod._report_value({"a": None}, "a.b") is None
    with pytest.raises(IntegrationTestError, match="parameter_file or parameters"):
        integration_mod._parameter_file(tmp_path, tmp_path, "demo", {"run": {}})

    direct_child = tmp_path / "child.txt"
    direct_child.write_text("remove\n", encoding="utf-8")
    integration_mod._cleanup_setup_files([direct_child], tmp_path)
    assert not direct_child.exists()


def test_resolve_bcftools_from_env_file_and_error_paths(tmp_path, monkeypatch):
    env_dir = tmp_path / "workflows" / "bash" / "gatk-4.6"
    env_dir.mkdir(parents=True)
    bcftools = tmp_path / "tools" / "bcftools"
    bcftools.parent.mkdir()
    bcftools.write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    (env_dir / "env.sh").write_text(f"BCFTOOLS={bcftools}\n", encoding="utf-8")
    assert integration_mod._resolve_bcftools(tmp_path, {}) == str(bcftools)

    (env_dir / "env.sh").write_text("exit 7\n", encoding="utf-8")
    with pytest.raises(IntegrationTestError, match="Could not resolve BCFTOOLS"):
        integration_mod._resolve_bcftools(tmp_path, {})

    monkeypatch.setattr(integration_mod.shutil, "which", lambda name, path=None: None)
    with pytest.raises(IntegrationTestError, match="bcftools is required only"):
        integration_mod._resolve_bcftools(tmp_path, {"BCFTOOLS": "bcftools"})

    monkeypatch.setattr(integration_mod.shutil, "which", lambda name, path=None: f"/usr/bin/{name}")
    assert integration_mod._resolve_bcftools(tmp_path, {"BCFTOOLS": "bcftools"}) == "bcftools"


def test_setup_run_and_setup_file_error_branches(tmp_path, monkeypatch):
    workdir = tmp_path / "work"
    workdir.mkdir()
    contract = {
        "setup_files_from_runs": [
            {
                "run": "wes-bash",
                "source": "02_varcall/sample.g.vcf.gz",
                "source_index": "02_varcall/sample.g.vcf.gz.tbi",
                "sample_map": "maps/cohort.tsv",
                "samples": ["S1"],
            }
        ]
    }
    with pytest.raises(IntegrationTestError, match="Setup run was not executed"):
        integration_mod._write_setup_files_from_runs(workdir, contract, {})

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    with pytest.raises(IntegrationTestError, match="Setup source does not exist"):
        integration_mod._write_setup_files_from_runs(workdir, contract, {"wes-bash": run_dir})

    source = run_dir / "02_varcall" / "sample.g.vcf.gz"
    source.parent.mkdir(parents=True)
    source.write_text("vcf\n", encoding="utf-8")
    with pytest.raises(IntegrationTestError, match="Setup source index does not exist"):
        integration_mod._write_setup_files_from_runs(workdir, contract, {"wes-bash": run_dir})

    (source.parent / "sample.g.vcf.gz.tbi").write_text("index\n", encoding="utf-8")
    empty_samples = {
        "setup_files_from_runs": [
            {
                "run": "wes-bash",
                "source": "02_varcall/sample.g.vcf.gz",
                "sample_map": "maps/cohort.tsv",
                "samples": [],
            }
        ]
    }
    with pytest.raises(IntegrationTestError, match="requires at least one sample"):
        integration_mod._write_setup_files_from_runs(workdir, empty_samples, {"wes-bash": run_dir})

    cleanup_target = workdir / "maps" / "cohort.tsv"
    cleanup_target.parent.mkdir(exist_ok=True)
    cleanup_target.write_text("S1\tfile\n", encoding="utf-8")
    (cleanup_target.parent / "keep.txt").write_text("keep parent non-empty\n", encoding="utf-8")
    integration_mod._cleanup_setup_files([cleanup_target], workdir)
    assert not cleanup_target.exists()
    assert cleanup_target.parent.exists()

    with pytest.raises(IntegrationTestError, match="Recursive integration setup dependency"):
        integration_mod._run_setup_runs(
            project_root=tmp_path,
            contract={"setup_runs": ["wes-bash"]},
            threads=1,
            runtime_profile="local",
            keep_external_work=False,
            stack=["wes-bash"],
        )
    with pytest.raises(IntegrationTestError, match="Unknown integration setup dependency"):
        integration_mod._run_setup_runs(
            project_root=tmp_path,
            contract={"setup_runs": ["unknown"]},
            threads=1,
            runtime_profile="local",
            keep_external_work=False,
            stack=[],
        )

    def fake_failed_run(**kwargs):
        return kwargs["selection"].label, "failed", "boom"

    monkeypatch.setattr(integration_mod, "_run_one", fake_failed_run)
    with pytest.raises(IntegrationTestError, match="Setup run did not pass"):
        integration_mod._run_setup_runs(
            project_root=tmp_path,
            contract={"setup_runs": ["wes-bash"]},
            threads=1,
            runtime_profile="local",
            keep_external_work=False,
            stack=[],
        )

    def fake_passed_run(**kwargs):
        return kwargs["selection"].label, "passed", str(run_dir)

    monkeypatch.setattr(integration_mod, "_run_one", fake_passed_run)
    assert integration_mod._run_setup_runs(
        project_root=tmp_path,
        contract={"setup_runs": ["wes-bash"]},
        threads=1,
        runtime_profile="local",
        keep_external_work=True,
        stack=[],
    ) == {"wes-bash": run_dir}


def _minimal_staged_contract() -> dict:
    return {
        "_contract_path": "/tmp/native-wes-cohort-bash-sharded.yaml",
        "workflow_log": "shard.log",
        "run": {"run_glob": "cbicall_shard_*"},
        "staged_finalize": {
            "shard_vcfs": ["02_varcall/cohort.chr22.gv.raw.vcf.gz"],
            "raw_vcfs_list": "raw_vcfs.list",
            "gathered_vcf": "cohort.gathered.gv.raw.vcf.gz",
            "finalize": {
                "parameters": {"mode": "cohort", "cohort_stage": "finalize"},
                "run_glob": "cbicall_finalize_*",
                "workflow_log": "finalize.log",
                "required_files": ["run-report.json", "finalize.log"],
                "json_expectations": [{"path": "status", "equals": "success"}],
                "cleanup_paths": ["heavy-work"],
            },
        },
    }


def test_run_staged_finalize_gathers_runs_validates_and_cleans(tmp_path, monkeypatch, capsys):
    project_root = tmp_path / "project"
    workdir = project_root / "examples" / "input"
    shard_run_dir = workdir / "cbicall_shard_001"
    raw_vcf = shard_run_dir / "02_varcall" / "cohort.chr22.gv.raw.vcf.gz"
    raw_vcf.parent.mkdir(parents=True)
    raw_vcf.write_text("vcf\n", encoding="utf-8")
    Path(str(raw_vcf) + ".tbi").write_text("index\n", encoding="utf-8")
    (project_root / "bin").mkdir(parents=True)
    (project_root / "bin" / "cbicall").write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    monkeypatch.setattr(integration_mod, "_resolve_bcftools", lambda project_root, env: "bcftools")

    calls = []

    def fake_run(cmd, cwd, env, stdout, stderr, text, check):
        calls.append(cmd)
        if cmd[:2] == ["bcftools", "concat"]:
            out_path = Path(cmd[cmd.index("-o") + 1])
            out_path.write_text("gathered\n", encoding="utf-8")
            return SimpleNamespace(returncode=0, stdout="concat ok\n")
        if cmd[:2] == ["bcftools", "index"]:
            Path(str(cmd[-1]) + ".tbi").write_text("index\n", encoding="utf-8")
            return SimpleNamespace(returncode=0, stdout="index ok\n")
        run_dir = workdir / "cbicall_finalize_001"
        run_dir.mkdir()
        (run_dir / "run-report.json").write_text(json.dumps({"status": "success"}), encoding="utf-8")
        (run_dir / "finalize.log").write_text("ok\n", encoding="utf-8")
        (run_dir / "heavy-work").mkdir()
        param_file = Path(cmd[cmd.index("-p") + 1])
        params = yaml.safe_load(param_file.read_text(encoding="utf-8"))
        assert params["input_vcf"] == str(workdir / "cohort.gathered.gv.raw.vcf.gz")
        return SimpleNamespace(returncode=0, stdout=f"Working directory: {run_dir}\n")

    monkeypatch.setattr(integration_mod.subprocess, "run", fake_run)
    run_dir = integration_mod._run_staged_finalize(
        project_root=project_root,
        selection=integration_mod.TESTS["wes-cohort-bash-sharded"],
        workdir=workdir,
        shard_run_dir=shard_run_dir,
        contract=_minimal_staged_contract(),
        threads=4,
        runtime_profile="local",
        keep_external_work=False,
    )

    assert run_dir == workdir / "cbicall_finalize_001"
    assert (workdir / "raw_vcfs.list").read_text(encoding="utf-8") == f"{raw_vcf}\n"
    assert not (run_dir / "heavy-work").exists()
    assert not list(workdir.glob("cbicall-wes-cohort-bash-sharded-finalize.*.yaml"))
    assert calls[0][:2] == ["bcftools", "concat"]
    assert calls[1][:2] == ["bcftools", "index"]
    out = capsys.readouterr().out
    assert "Gathering staged raw cohort VCFs" in out
    assert "Validating staged finalize contract" in out


def test_run_staged_finalize_reports_contract_and_command_errors(tmp_path, monkeypatch):
    project_root = tmp_path / "project"
    workdir = project_root / "examples" / "input"
    shard_run_dir = workdir / "cbicall_shard_001"
    raw_vcf = shard_run_dir / "02_varcall" / "cohort.chr22.gv.raw.vcf.gz"
    raw_vcf.parent.mkdir(parents=True)
    workdir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(integration_mod, "_resolve_bcftools", lambda project_root, env: "bcftools")

    base_kwargs = {
        "project_root": project_root,
        "selection": integration_mod.TESTS["wes-cohort-bash-sharded"],
        "workdir": workdir,
        "shard_run_dir": shard_run_dir,
        "threads": 1,
        "runtime_profile": "local",
        "keep_external_work": False,
    }
    with pytest.raises(IntegrationTestError, match="staged_finalize must be a mapping"):
        integration_mod._run_staged_finalize(contract={"staged_finalize": ["not-a-mapping"]}, **base_kwargs)
    with pytest.raises(IntegrationTestError, match="must list at least one"):
        integration_mod._run_staged_finalize(contract={"staged_finalize": {"shard_vcfs": []}}, **base_kwargs)
    with pytest.raises(IntegrationTestError, match="raw VCF does not exist"):
        integration_mod._run_staged_finalize(contract=_minimal_staged_contract(), **base_kwargs)

    raw_vcf.write_text("vcf\n", encoding="utf-8")
    with pytest.raises(IntegrationTestError, match="raw VCF index does not exist"):
        integration_mod._run_staged_finalize(contract=_minimal_staged_contract(), **base_kwargs)
    Path(str(raw_vcf) + ".tbi").write_text("index\n", encoding="utf-8")

    def fail_concat(cmd, cwd, env, stdout, stderr, text, check):
        return SimpleNamespace(returncode=2, stdout="concat failed\n")

    monkeypatch.setattr(integration_mod.subprocess, "run", fail_concat)
    with pytest.raises(IntegrationTestError, match="gather command failed"):
        integration_mod._run_staged_finalize(contract=_minimal_staged_contract(), **base_kwargs)

    contract = _minimal_staged_contract()
    contract["staged_finalize"]["finalize"] = {}

    def ok_gather(cmd, cwd, env, stdout, stderr, text, check):
        if cmd[:2] == ["bcftools", "concat"]:
            Path(cmd[cmd.index("-o") + 1]).write_text("gathered\n", encoding="utf-8")
        elif cmd[:2] == ["bcftools", "index"]:
            Path(str(cmd[-1]) + ".tbi").write_text("index\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="ok\n")

    monkeypatch.setattr(integration_mod.subprocess, "run", ok_gather)
    with pytest.raises(IntegrationTestError, match="finalize.parameters is required"):
        integration_mod._run_staged_finalize(contract=contract, **base_kwargs)

    contract = _minimal_staged_contract()
    contract["run"] = {}
    contract["staged_finalize"]["finalize"].pop("run_glob")
    with pytest.raises(IntegrationTestError, match="finalize.run_glob is required"):
        integration_mod._run_staged_finalize(contract=contract, **base_kwargs)


def test_run_staged_finalize_reports_finalize_failures(tmp_path, monkeypatch):
    project_root = tmp_path / "project"
    workdir = project_root / "examples" / "input"
    shard_run_dir = workdir / "cbicall_shard_001"
    raw_vcf = shard_run_dir / "02_varcall" / "cohort.chr22.gv.raw.vcf.gz"
    raw_vcf.parent.mkdir(parents=True)
    raw_vcf.write_text("vcf\n", encoding="utf-8")
    Path(str(raw_vcf) + ".tbi").write_text("index\n", encoding="utf-8")
    (project_root / "bin").mkdir(parents=True)
    (project_root / "bin" / "cbicall").write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    monkeypatch.setattr(integration_mod, "_resolve_bcftools", lambda project_root, env: "bcftools")

    def fail_finalize(cmd, cwd, env, stdout, stderr, text, check):
        if cmd[:2] == ["bcftools", "concat"]:
            Path(cmd[cmd.index("-o") + 1]).write_text("gathered\n", encoding="utf-8")
            return SimpleNamespace(returncode=0, stdout="concat ok\n")
        if cmd[:2] == ["bcftools", "index"]:
            Path(str(cmd[-1]) + ".tbi").write_text("index\n", encoding="utf-8")
            return SimpleNamespace(returncode=0, stdout="index ok\n")
        return SimpleNamespace(returncode=3, stdout="finalize failed\n")

    monkeypatch.setattr(integration_mod.subprocess, "run", fail_finalize)
    with pytest.raises(IntegrationTestError, match="finalize cbicall command failed"):
        integration_mod._run_staged_finalize(
            project_root=project_root,
            selection=integration_mod.TESTS["wes-cohort-bash-sharded"],
            workdir=workdir,
            shard_run_dir=shard_run_dir,
            contract=_minimal_staged_contract(),
            threads=1,
            runtime_profile="local",
            keep_external_work=True,
        )
    assert not list(workdir.glob("cbicall-wes-cohort-bash-sharded-finalize.*.yaml"))

    def no_run_dir(cmd, cwd, env, stdout, stderr, text, check):
        if cmd[:2] == ["bcftools", "concat"]:
            Path(cmd[cmd.index("-o") + 1]).write_text("gathered\n", encoding="utf-8")
        elif cmd[:2] == ["bcftools", "index"]:
            Path(str(cmd[-1]) + ".tbi").write_text("index\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="ok\n")

    monkeypatch.setattr(integration_mod.subprocess, "run", no_run_dir)
    with pytest.raises(IntegrationTestError, match="finished but no run directory"):
        integration_mod._run_staged_finalize(
            project_root=project_root,
            selection=integration_mod.TESTS["wes-cohort-bash-sharded"],
            workdir=workdir,
            shard_run_dir=shard_run_dir,
            contract=_minimal_staged_contract(),
            threads=1,
            runtime_profile="local",
            keep_external_work=True,
        )


def test_run_one_returns_finalize_run_dir_for_staged_contract(tmp_path, monkeypatch):
    workdir = tmp_path / "examples" / "input"
    workdir.mkdir(parents=True)
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "cbicall").write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    _write_fixture(
        tmp_path,
        "staged.yaml",
        {
            "workdir": "examples/input",
            "workflow_log": "shard.log",
            "run": {"parameters": {"mode": "cohort"}, "base_dir": ".", "run_glob": "cbicall_shard_*"},
            "required_files": ["run-report.json", "shard.log"],
            "staged_finalize": {"enabled": True},
        },
    )
    selection = integration_mod.TestSelection("staged", "Staged", "staged.yaml")

    def fake_run(cmd, cwd, env, stdout, stderr, text, check):
        run_dir = workdir / "cbicall_shard_001"
        run_dir.mkdir()
        (run_dir / "run-report.json").write_text(json.dumps({"status": "success"}), encoding="utf-8")
        (run_dir / "shard.log").write_text("ok\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout=f"Working directory: {run_dir}\n")

    finalize_dir = workdir / "cbicall_finalize_001"

    def fake_finalize(**kwargs):
        assert kwargs["shard_run_dir"] == workdir / "cbicall_shard_001"
        return finalize_dir

    monkeypatch.setattr(integration_mod.subprocess, "run", fake_run)
    monkeypatch.setattr(integration_mod, "_run_staged_finalize", fake_finalize)
    label, status, detail = integration_mod._run_one(
        project_root=tmp_path,
        selection=selection,
        threads=2,
        runtime_profile="local",
        skip_missing_optional=False,
        keep_external_work=False,
    )
    assert (label, status, detail) == ("Staged", "passed", str(finalize_dir))


def test_selected_tests_includes_explicit_gatk35():
    args = SimpleNamespace(
        all=False,
        wes_bash=False,
        mit_bash=False,
        wes_cohort_bash=False,
        wes_cohort_bash_sharded=False,
        wes_bash_gatk35=True,
        wes_snakemake=False,
        wes_nextflow=False,
        wes_cromwell=False,
        nf_core_demo=False,
        nf_core_sarek=False,
    )
    assert [item.key for item in integration_mod.selected_tests_from_args(args)] == ["wes-bash-gatk35"]
