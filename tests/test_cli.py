import io
import gzip
import json
import sys
import time
from types import SimpleNamespace

import pytest

from cbicall import cli as cli_mod


def test_write_log_creates_json(tmp_path):
    cfg = {"project_dir": str(tmp_path)}
    arg = {"threads": 4}
    param = {"pipeline": "wes"}

    cli_mod.write_log(cfg, arg, param)

    log_path = tmp_path / "log.json"
    assert log_path.is_file()

    data = json.loads(log_path.read_text(encoding="utf-8"))
    assert data["arg"] == arg
    assert data["config"]["project_dir"] == cfg["project_dir"]
    assert data["param"]["pipeline"] == "wes"


def test_write_run_report_creates_compact_summary(tmp_path, monkeypatch):
    monkeypatch.setattr(cli_mod.shutil, "which", lambda command: f"/usr/bin/{command}")
    monkeypatch.setattr(
        cli_mod.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0, stdout="GNU bash, version 5.2.21\n", stderr=""),
    )
    entrypoint = tmp_path / "wes_single.sh"
    entrypoint.write_text("#!/bin/sh\necho ok\n", encoding="utf-8")
    env_file = tmp_path / "env.sh"
    env_file.write_text("DATADIR=/tmp\n", encoding="utf-8")
    catalog_path = tmp_path.parent / f"{tmp_path.name}-resource-catalog.json"
    catalog_path.write_text(
        json.dumps(
            {
                "resources": {
                    "cbicall-germline-resources-v1": {
                        "tools": {
                            "gatk4": {"version": "4.6.2.0", "path_hint": "NGSutils/gatk/gatk-4.6.2.0/gatk"},
                            "bwa": {"version": "0.7.18", "path_hint": "NGSutils/bwa-0.7.18/bwa"},
                        }
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "workflow.log").write_text("workflow log\n", encoding="utf-8")
    (tmp_path / "log.json").write_text("{}\n", encoding="utf-8")
    stats_dir = tmp_path / "03_stats"
    stats_dir.mkdir()
    (stats_dir / "sample.vcf.sha256.txt").write_text(
        "FILE=sample.vcf.gz\n"
        "ALGORITHM=sha256\n"
        "RAW_SHA256=raw\n"
        "NORMALIZED_PATTERN=^#\n"
        "NORMALIZED_SORT=LC_ALL=C\n"
        "NORMALIZED_RECORDS=10\n"
        "NORMALIZED_SHA256=normalized\n",
        encoding="utf-8",
    )
    resolved = cli_mod.ResolvedConfig.from_mapping(
        {
            "project_dir": str(tmp_path),
            "run_id": "RIDREPORT",
            "workflow_engine": "bash",
            "profile": "cnag-hpc",
            "genome": "b37",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "bash",
                "pipeline": "wes",
                "mode": "single",
                "gatk_version": "gatk-4.6",
                "pipeline_version": "v1",
                "entrypoint": str(entrypoint),
                "config_file": None,
                "helpers": {"env": str(env_file)},
            },
            "resources": {
                "bundle": {
                    "key": "cbicall-germline-resources-v1",
                    "version": "v1",
                    "fingerprint": "abc",
                    "catalog": str(catalog_path),
                }
            },
            "version": "1.2.3",
        }
    )

    report = cli_mod.write_run_report(
        resolved,
        {"threads": 2, "paramfile": "params.yaml", "profile": "cnag-hpc"},
        {"cleanup_bam": False},
        elapsed_seconds=12.3456,
        workflow_log=tmp_path / "workflow.log",
    )

    data = json.loads(report.read_text(encoding="utf-8"))
    assert data["status"] == "success"
    assert data["profile"] == "cnag-hpc"
    assert data["framework"]["version"] == "1.2.3"
    assert data["runtime"]["python"]["version"] == sys.version.split()[0]
    assert data["runtime"]["engine"]["name"] == "bash"
    assert data["elapsed_seconds"] == 12.346
    assert data["resources"]["bundle"]["version"] == "v1"
    assert data["resources"]["bundle"]["fingerprint"] == "abc"
    assert data["workflow_log"].endswith("workflow.log")
    assert data["workflow"]["key"] == "bash/wes/single/gatk-4.6/v1"
    assert len(data["workflow"]["fingerprint"]) == 64
    assert [item["role"] for item in data["workflow"]["files"]] == ["entrypoint", "helper:env"]
    assert all(item["status"] == "present" for item in data["workflow"]["files"])
    inventory = data["outputs"]["file_inventory"]
    assert inventory["algorithm"] == "sha256"
    assert inventory["scope"] == "run directory relative file paths"
    assert inventory["entries"] == 5
    assert inventory["excluded"] == [".nextflow", ".nextflow*", "run-report.html", "run-report.json", "work"]
    assert "run-report.json" not in inventory["paths"]
    assert inventory["paths"] == ["03_stats/sample.vcf.sha256.txt", "env.sh", "log.json", "wes_single.sh", "workflow.log"]
    assert len(inventory["sha256"]) == 64
    assert data["outputs"]["vcf_hash_reports"][0]["normalized_sha256"] == "normalized"
    assert data["software_versions"]["scope"] == "resource_declared"
    assert data["software_versions"]["entries"]["gatk4"]["version"] == "4.6.2.0"
    assert len(data["software_versions"]["sha256"]) == 64
    html_report = cli_mod._write_run_report_html(report, data)
    assert html_report.is_file()
    html_text = html_report.read_text(encoding="utf-8")
    assert "<title>CBIcall Run Report</title>" in html_text
    assert "Human-readable summary generated from" in html_text
    assert "bash/wes/single/gatk-4.6/v1" in html_text
    assert "normalized" in html_text
    assert "Run Files" in html_text
    assert 'href="workflow.log"' in html_text
    assert 'href="log.json"' in html_text
    assert 'href="run-report.json"' in html_text
    assert "Software Versions" in html_text
    assert "resource_declared" in html_text
    assert "gatk4" in html_text
    assert "4.6.2.0" in html_text
    assert "03_stats/sample.vcf.sha256.txt" in html_text
    assert "href=\"03_stats/sample.vcf.sha256.txt\"" in html_text


def test_runtime_report_records_python_and_backend_version(monkeypatch):
    monkeypatch.setattr(cli_mod.shutil, "which", lambda command: f"/usr/bin/{command}")

    def fake_run(command, capture_output, text, timeout, check):
        assert command == ["/usr/bin/nextflow", "-version"]
        assert capture_output is True
        assert text is True
        assert timeout == 10
        assert check is False
        return SimpleNamespace(returncode=0, stdout="Nextflow version 25.10.2\n", stderr="")

    monkeypatch.setattr(cli_mod.subprocess, "run", fake_run)

    report = cli_mod._runtime_report("nextflow")

    assert report["python"]["version"] == sys.version.split()[0]
    assert report["engine"]["name"] == "nextflow"
    assert report["engine"]["path"] == "/usr/bin/nextflow"
    assert report["engine"]["version"] == "25.10.2"
    assert report["engine"]["status"] == "ok"


def test_runtime_report_marks_missing_backend(monkeypatch):
    monkeypatch.setattr(cli_mod.shutil, "which", lambda command: None)

    report = cli_mod._runtime_report("snakemake")

    assert report["engine"]["name"] == "snakemake"
    assert report["engine"]["status"] == "not_found"
    assert report["engine"]["version"] is None


def test_write_run_report_hashes_registry_canonical_vcfs(tmp_path, monkeypatch):
    monkeypatch.setattr(cli_mod.shutil, "which", lambda command: f"/usr/bin/{command}")
    monkeypatch.setattr(
        cli_mod.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(returncode=0, stdout="Nextflow version 25.10.2\n", stderr=""),
    )
    pipeline_info = tmp_path / "sarek" / "pipeline_info"
    pipeline_info.mkdir(parents=True)
    (pipeline_info / "params_2026-05-19_12-21-13.json").write_text("{}\n", encoding="utf-8")
    (pipeline_info / "execution_trace_2026-05-19_12-20-22.txt").write_text(
        "task_id\thash\tname\tstatus\tpeak_rss\tpeak_vmem\n"
        "1\taa/bbccdd\tNFCORE_SAREK:FASTQC (sample)\tCOMPLETED\t641.9 MB\t6.5 GB\n"
        "2\tbb/ccddee\tNFCORE_SAREK:GATK4_HAPLOTYPECALLER (sample)\tCOMPLETED\t951.4 MB\t12 GB\n",
        encoding="utf-8",
    )
    (pipeline_info / "execution_report_2026-05-19_12-20-22.html").write_text("<html></html>\n", encoding="utf-8")
    (pipeline_info / "execution_timeline_2026-05-19_12-20-22.html").write_text("<html></html>\n", encoding="utf-8")
    (pipeline_info / "pipeline_dag_2026-05-19_12-20-22.html").write_text("<html></html>\n", encoding="utf-8")
    (pipeline_info / "manifest_2026-05-19_12-20-22.bco.json").write_text("{}\n", encoding="utf-8")
    (pipeline_info / "nf_core_sarek_software_mqc_versions.yml").write_text(
        "GATK4_HAPLOTYPECALLER:\n  gatk4: 4.6.1.0\nWorkflow:\n  nf-core/sarek: v3.8.1\n  Nextflow: 25.10.2\n",
        encoding="utf-8",
    )
    multiqc_dir = tmp_path / "sarek" / "multiqc" / "multiqc_data"
    multiqc_dir.mkdir(parents=True)
    (tmp_path / "sarek" / "multiqc" / "multiqc_report.html").write_text("<html></html>\n", encoding="utf-8")
    vcf_dir = tmp_path / "sarek" / "variant_calling" / "haplotypecaller" / "CNAG99901P_ex"
    vcf_dir.mkdir(parents=True)
    vcf_path = vcf_dir / "CNAG99901P_ex.haplotypecaller.vcf.gz"
    with gzip.open(vcf_path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.write("chr22\t100\t.\tA\tG\t60\tPASS\t.\n")

    resolved = cli_mod.ResolvedConfig.from_mapping(
        {
            "project_dir": str(tmp_path),
            "run_id": "RIDNFCORE",
            "workflow_engine": "nextflow",
            "profile": "local",
            "genome": "external",
            "nfcore_profile": "singularity",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "nextflow",
                "pipeline": "sarek",
                "mode": "cohort",
                "gatk_version": "nf-core",
                "pipeline_version": "v1",
                "entrypoint": "nf-core/sarek",
                "config_file": None,
                "helpers": {},
                "metadata": {
                    "source_type": "nf-core",
                    "source": "nf-core/sarek",
                    "release": "3.8.1",
                    "default_outdir": "sarek",
                    "canonical_outputs": [
                        {
                            "name": "haplotypecaller_vcf",
                            "type": "vcf",
                            "pattern": "variant_calling/haplotypecaller/*/*.haplotypecaller.vcf.gz",
                        }
                    ],
                },
            },
            "resources": {
                "bundle": {"key": "nf-core-sarek-managed-resources-v1", "version": "sarek-3.8.1"}
            },
            "version": "1.2.3",
        }
    )

    report = cli_mod.write_run_report(
        resolved,
        {"threads": 2, "paramfile": "sarek.yaml"},
        {"cleanup_bam": False},
        elapsed_seconds=1.0,
        workflow_log=tmp_path / "nextflow.log",
    )
    data = json.loads(report.read_text(encoding="utf-8"))

    assert data["outputs"]["workflow_output_dir"] == str(tmp_path / "sarek")
    summary = data["outputs"]["external_summary"]
    assert summary["pipeline_info"]["dir"] == str(pipeline_info)
    assert summary["pipeline_info"]["params"].endswith("params_2026-05-19_12-21-13.json")
    assert summary["pipeline_info"]["trace"].endswith("execution_trace_2026-05-19_12-20-22.txt")
    assert summary["pipeline_info"]["timeline"].endswith("execution_timeline_2026-05-19_12-20-22.html")
    assert summary["pipeline_info"]["dag"].endswith("pipeline_dag_2026-05-19_12-20-22.html")
    assert summary["pipeline_info"]["manifest"].endswith("manifest_2026-05-19_12-20-22.bco.json")
    assert summary["pipeline_info"]["software_versions"].endswith("nf_core_sarek_software_mqc_versions.yml")
    assert summary["multiqc"]["report"].endswith("multiqc_report.html")
    assert summary["multiqc"]["data_dir"] == str(multiqc_dir)
    canonical = data["outputs"]["canonical_outputs"][0]
    assert canonical["name"] == "haplotypecaller_vcf"
    assert canonical["status"] == "present"
    assert canonical["matches"] == [str(vcf_path)]
    assert summary["canonical_outputs"] == data["outputs"]["canonical_outputs"]
    vcf_report = data["outputs"]["vcf_hash_reports"][0]
    assert vcf_report["source"] == "registry_canonical_output"
    assert vcf_report["name"] == "haplotypecaller_vcf"
    assert vcf_report["normalized_records"] == 1
    assert len(vcf_report["normalized_sha256"]) == 64
    assert data["software_versions"]["status"] == "parsed"
    assert data["software_versions"]["scope"] == "workflow_reported"
    assert data["software_versions"]["entries"]["GATK4_HAPLOTYPECALLER"]["gatk4"] == "4.6.1.0"
    assert len(data["software_versions"]["sha256"]) == 64
    assert data["execution_trace"]["tasks"] == 2
    assert data["execution_trace"]["status_counts"] == {"COMPLETED": 2}
    assert data["execution_trace"]["max_peak_rss"]["value"] == "951.4 MB"
    assert data["execution_trace"]["max_peak_vmem"]["value"] == "12 GB"
    html_report = cli_mod._write_run_report_html(report, data)
    html_text = html_report.read_text(encoding="utf-8")
    assert "nf-core/sarek" in html_text
    assert "External Reports" in html_text
    assert "Software Versions" in html_text
    assert "GATK4_HAPLOTYPECALLER" in html_text
    assert "Execution Trace" in html_text
    assert "Max peak RSS" in html_text
    assert "951.4 MB" in html_text
    assert "Execution timeline" in html_text
    assert "Pipeline DAG" in html_text
    assert "haplotypecaller_vcf" in html_text
    assert "CNAG99901P_ex.haplotypecaller.vcf.gz" in html_text


def test_print_external_output_pointers(capsys, tmp_path):
    report = {
        "outputs": {
            "workflow_output_dir": str(tmp_path / "run" / "sarek"),
            "external_summary": {
                "pipeline_info": {
                    "dir": str(tmp_path / "run" / "sarek" / "pipeline_info"),
                    "trace": str(tmp_path / "run" / "sarek" / "pipeline_info" / "execution_trace.txt"),
                },
                "multiqc": {
                    "report": str(tmp_path / "run" / "sarek" / "multiqc" / "multiqc_report.html"),
                },
                "canonical_outputs": [
                    {
                        "matches": [
                            str(tmp_path / "run" / "sarek" / "variant_calling" / "haplotypecaller" / "sample.vcf.gz")
                        ]
                    }
                ],
            },
        }
    }

    cli_mod._print_external_output_pointers(report)
    out = capsys.readouterr().out
    assert "nf-core outputs" in out
    assert "Workflow out" in out
    assert "Pipeline" in out
    assert "Trace" in out
    assert "MultiQC" in out
    assert "Canonical" in out
    assert "sample.vcf.gz" in out


def test_compare_runs_reports_workflow_and_output_differences(tmp_path, capsys):
    run_a = tmp_path / "run_a"
    run_b = tmp_path / "run_b"
    run_a.mkdir()
    run_b.mkdir()

    base = {
        "status": "success",
        "framework": {"name": "CBIcall", "version": "1.2.3"},
        "runtime": {
            "python": {"version": "3.12.3"},
            "engine": {"name": "bash", "version": "5.2.21"},
        },
        "workflow": {
            "key": "bash/wes/single/gatk-4.6/v1",
            "pipeline_version": "v1",
            "entrypoint": "workflows/bash/gatk-4.6/wes_single.sh",
            "fingerprint": "workflow-a",
            "files": [
                {"role": "entrypoint", "path": "wes_single.sh", "sha256": "aaa"},
                {"role": "helper:env", "path": "env.sh", "sha256": "env"},
            ],
        },
        "resources": {"bundle": {"key": "cbicall-germline-resources-v1", "version": "v1", "fingerprint": "res"}},
        "software_versions": {"sha256": "software-a", "scope": "resource_declared"},
        "execution_trace": {"tasks": 2, "max_peak_rss": {"bytes": 1000}, "max_peak_vmem": {"bytes": 2000}},
        "outputs": {
            "file_inventory": {"entries": 3, "sha256": "manifest-a"},
            "vcf_hash_reports": [
                {"file": "sample.vcf.gz", "normalized_sha256": "vcf"}
            ]
        },
    }
    changed = json.loads(json.dumps(base))
    changed["runtime"]["python"]["version"] = "3.12.4"
    changed["workflow"]["fingerprint"] = "workflow-b"
    changed["workflow"]["files"][0]["sha256"] = "bbb"
    changed["software_versions"]["sha256"] = "software-b"
    changed["execution_trace"]["max_peak_rss"]["bytes"] = 1500
    changed["outputs"]["file_inventory"]["sha256"] = "manifest-b"
    changed["outputs"]["vcf_hash_reports"][0]["normalized_sha256"] = "vcf"

    (run_a / "run-report.json").write_text(json.dumps(base), encoding="utf-8")
    (run_b / "run-report.json").write_text(json.dumps(changed), encoding="utf-8")

    report = tmp_path / "compare-report.txt"
    html_report = tmp_path / "compare-report.html"
    assert cli_mod._run_compare_runs_command(
        [str(run_a), str(run_b), "--no-color", "--output", str(report), "--html", str(html_report)]
    ) == 0
    out = capsys.readouterr().out
    assert "Run Comparison" in out
    assert "CBIcall ver" in out
    assert "Python ver" in out
    assert "Engine ver" in out
    assert "Task count" in out
    assert "Max peak RSS" in out
    assert "Workflow hash" in out
    assert "Software versions" in out
    assert "Resource ver" in out
    assert "File inventory" in out
    assert "different" in out
    assert "entrypoint" in out
    assert "sha256" in out
    assert "sample.vcf.gz" in out
    assert "Legend" in out
    assert "not available" in out
    assert "Report" in out
    assert "HTML" in out
    report_text = report.read_text(encoding="utf-8")
    assert "Run Comparison" in report_text
    assert "Workflow hash" in report_text
    assert "Software versions" in report_text
    assert "Legend" in report_text
    html_text = html_report.read_text(encoding="utf-8")
    assert "<title>CBIcall Run Comparison</title>" in html_text
    assert "Run Comparison" in html_text
    assert "Workflow hash" in html_text
    assert "Software versions" in html_text
    assert "Status summary" in html_text
    assert "class=\"pill different\"" in html_text


def test_compare_runs_accepts_multiple_runs_as_baseline_matrix(tmp_path, capsys):
    runs = [tmp_path / f"run_{idx}" for idx in range(3)]
    for run in runs:
        run.mkdir()

    base = {
        "status": "success",
        "framework": {"name": "CBIcall", "version": "1.2.3"},
        "runtime": {
            "python": {"version": "3.12.3"},
            "engine": {"name": "bash", "version": "5.2.21"},
        },
        "workflow": {
            "key": "bash/wes/single/gatk-4.6/v1",
            "pipeline_version": "v1",
            "entrypoint": "wes_single.sh",
            "fingerprint": "workflow-a",
            "files": [
                {"role": "entrypoint", "path": "wes_single.sh", "sha256": "aaa"},
            ],
        },
        "resources": {"bundle": {"key": "cbicall-germline-resources-v1", "version": "v1", "fingerprint": "res"}},
        "software_versions": {"sha256": "software-a", "scope": "resource_declared"},
        "execution_trace": {"tasks": 2, "max_peak_rss": {"bytes": 1000}, "max_peak_vmem": {"bytes": 2000}},
        "outputs": {
            "file_inventory": {"entries": 3, "sha256": "manifest-a"},
            "vcf_hash_reports": [
                {"file": "sample.vcf.gz", "normalized_sha256": "vcf"}
            ]
        },
    }
    same = json.loads(json.dumps(base))
    different = json.loads(json.dumps(base))
    different["runtime"]["engine"]["version"] = "5.2.22"
    different["workflow"]["fingerprint"] = "workflow-c"
    different["software_versions"]["sha256"] = "software-c"
    different["execution_trace"]["tasks"] = 3
    different["execution_trace"]["max_peak_vmem"]["bytes"] = 2500
    different["outputs"]["file_inventory"]["entries"] = 4
    different["outputs"]["file_inventory"]["sha256"] = "manifest-c"
    different["outputs"]["vcf_hash_reports"][0]["normalized_sha256"] = "vcf-c"

    for run, report in zip(runs, [base, same, different]):
        (run / "run-report.json").write_text(json.dumps(report), encoding="utf-8")

    assert cli_mod._run_compare_runs_command([str(run) for run in runs] + ["--no-color"]) == 0
    out = capsys.readouterr().out
    assert "Run Matrix" in out
    assert "Baseline" in out
    assert "Workflow hash" in out
    assert "Software versions" in out
    assert "Python ver" in out
    assert "Engine ver" in out
    assert "Task count" in out
    assert "Max peak VMEM" in out
    assert "Resource ver" in out
    assert "File inventory" in out
    assert "different: " in out
    assert "run_2/run-report.json" in out
    assert "sample.vcf.gz" in out
    assert "same: 1.2.3" in out
    assert "baseline=vcf" in out
    assert "vcf-c" in out


def test_run_with_spinner_no_spinner_calls_function():
    calls = []

    def fn(x, y):
        calls.append((x, y))
        return x + y

    result = cli_mod.run_with_spinner(fn, 2, 3, no_spinner=True)
    assert result == 5
    assert calls == [(2, 3)]


def test_run_with_spinner_propagates_exception():
    def fn():
        raise ValueError("error")

    with pytest.raises(ValueError):
        cli_mod.run_with_spinner(fn, no_spinner=True)


def test_run_with_spinner_spinner_path(monkeypatch):
    class FakeStdout(io.StringIO):
        def isatty(self):
            return True

    fake_stdout = FakeStdout()
    monkeypatch.setattr(sys, "stdout", fake_stdout)

    def fn():
        time.sleep(0.01)
        return "ok"

    result = cli_mod.run_with_spinner(fn)
    assert result == "ok"
    assert fake_stdout.getvalue() != ""


def test_run_with_spinner_non_tty_bypasses_spinner(monkeypatch):
    class FakeStdout(io.StringIO):
        def isatty(self):
            return False

    fake_stdout = FakeStdout()
    monkeypatch.setattr(sys, "stdout", fake_stdout)

    calls = {"n": 0}

    def fn():
        calls["n"] += 1
        return 123

    assert cli_mod.run_with_spinner(fn) == 123
    assert calls["n"] == 1


def test_colors_enabled_and_code(monkeypatch):
    class TtyStdout(io.StringIO):
        def isatty(self):
            return True

    monkeypatch.setattr(sys, "stdout", TtyStdout())
    monkeypatch.delenv("ANSI_COLORS_DISABLED", raising=False)
    monkeypatch.delenv("NO_COLOR", raising=False)
    assert cli_mod._colors_enabled() is True
    assert cli_mod._code("X") == "X"

    monkeypatch.setenv("ANSI_COLORS_DISABLED", "1")
    assert cli_mod._colors_enabled() is False
    assert cli_mod._code("X") == ""


def test_colors_disabled_for_non_tty(monkeypatch):
    class NonTtyStdout(io.StringIO):
        def isatty(self):
            return False

    monkeypatch.setattr(sys, "stdout", NonTtyStdout())
    monkeypatch.delenv("ANSI_COLORS_DISABLED", raising=False)
    monkeypatch.delenv("NO_COLOR", raising=False)

    assert cli_mod._colors_enabled() is False
    assert cli_mod._code("X") == ""


def test_colors_disabled_by_no_color(monkeypatch):
    class TtyStdout(io.StringIO):
        def isatty(self):
            return True

    monkeypatch.setattr(sys, "stdout", TtyStdout())
    monkeypatch.delenv("ANSI_COLORS_DISABLED", raising=False)
    monkeypatch.setenv("NO_COLOR", "1")

    assert cli_mod._colors_enabled() is False


def test_format_duration():
    assert cli_mod._format_duration(5.2) == "5s"
    assert cli_mod._format_duration(65.0) == "1m 5s"
    assert cli_mod._format_duration(3661.0) == "1h 1m 1s"


def test_short_path_for_none_and_short_values():
    assert cli_mod._short_path(None) == "(undef)"
    short = "/tmp/a/b"
    assert cli_mod._short_path(short) == short


def test_print_config_includes_workflow_block(capsys):
    cfg = {
        "workflow_engine": "snakemake",
        "workflow": {"entrypoint": "/tmp/wes_single.smk", "engine": "snakemake"},
        "x": None,
    }
    cli_mod._print_config(cfg)
    out = capsys.readouterr().out
    assert "workflow" in out
    assert "entrypoint" in out
    assert "(undef)" in out


def test_validate_parameters_command_prints_parameters_ok(monkeypatch, tmp_path, capsys):
    param_file = tmp_path / "parameters.yaml"
    param_file.write_text("pipeline: wes\n", encoding="utf-8")

    monkeypatch.setattr(cli_mod.config_mod, "read_param_file", lambda _: {"pipeline": "wes"})
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
            "project_dir": str(tmp_path / "run"),
            "run_id": "RIDPARAM",
            "workflow_engine": "bash",
            "profile": "local",
            "genome": "b37",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "pipeline_version": "v1",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "bash",
                "pipeline": "wes",
                "mode": "single",
                "gatk_version": "gatk-4.6",
                "pipeline_version": "v1",
                "entrypoint": "/tmp/wes_single.sh",
                "config_file": None,
                "helpers": {"env": "/tmp/env.sh"},
            },
            "resources": {
                "bundle": {
                    "key": "cbicall-germline-resources-v1",
                    "version": "v1",
                    "fingerprint": "abc",
                    "runtime_check": {"status": "ok", "datadir": "/tmp/resources"},
                }
            },
        },
    )

    rc = cli_mod._run_validate_parameters_command(["-p", str(param_file), "--no-color"])

    assert rc == 0
    out = capsys.readouterr().out
    assert "Parameters OK" in out
    assert "Configuration OK" not in out
    assert "Param file" in out


def test_validate_registry_command_uses_default_registry(capsys):
    rc = cli_mod._run_validate_registry_command(["--no-color"])

    assert rc == 0
    out = capsys.readouterr().out
    assert "Registry OK" in out
    assert "cbicall-workflow-registry.yaml" in out
    assert "cbicall-workflow-registry.schema.json" in out
    assert "Backends" in out
    assert "bash, nextflow, snakemake" in out
    assert "External" in out
    assert "nf-core" in out
    assert "Engines" not in out


def test_validate_resources_command_uses_default_catalog(capsys):
    rc = cli_mod._run_validate_resources_command(["--no-color"])

    assert rc == 0
    out = capsys.readouterr().out
    assert "Resources OK" in out
    assert "cbicall-resource-catalog.json" in out
    assert "Compatible workflows" in out


def test_validate_resources_command_can_filter_resource(capsys):
    rc = cli_mod._run_validate_resources_command(
        ["--no-color", "--resource", "cbicall-germline-resources-v1"]
    )

    assert rc == 0
    out = capsys.readouterr().out
    assert "Resources OK" in out
    assert "Resource key" in out
    assert "cbicall-germline-resources-v1" in out
    assert "Bundle resources" in out


def test_validate_resources_command_accepts_short_resource_flag(capsys):
    rc = cli_mod._run_validate_resources_command(
        ["--no-color", "-r", "cbicall-germline-resources-v1"]
    )

    assert rc == 0
    out = capsys.readouterr().out
    assert "Resource key" in out
    assert "cbicall-germline-resources-v1" in out


def test_main_happy_path(monkeypatch, tmp_path):
    def fake_usage(version):
        return {
            "threads": 4,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": False,
            "nocolor": False,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)

    fake_param = {
        "pipeline": "wes",
        "mode": "single",
        "input_dir": None,
        "sample_map": None,
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "cleanup_bam": False,
    }

    monkeypatch.setattr(cli_mod.config_mod, "read_param_file", lambda _: fake_param)
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
                "project_dir": str(tmp_path / "proj"),
            "id": "ID123",
            "genome": "b37",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "bash",
                "pipeline": "wes",
                "mode": "single",
                "gatk_version": "gatk-3.5",
                "entrypoint": "/x.sh",
                "config_file": None,
                "helpers": {},
            },
        },
    )

    logs = {}
    monkeypatch.setattr(
        cli_mod, "write_log",
        lambda cfg, arg, param: logs.update({"cfg": cfg, "arg": arg, "param": param})
    )

    class FakeDNAseq:
        def __init__(self, settings):
            logs["settings"] = settings

        def variant_calling(self):
            return True

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    class FakeGoodBye:
        def say_goodbye(self):
            return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    rc = cli_mod.main()
    assert rc == 0
    assert logs["settings"].threads == 4
    assert logs["settings"].run_id == "ID123"
    assert logs["settings"].workflow.entrypoint == "/x.sh"


def test_main_run_subcommand_happy_path(monkeypatch, tmp_path, capsys):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("pipeline: wes\nmode: single\n", encoding="utf-8")
    monkeypatch.setattr(
        sys,
        "argv",
        ["cbicall", "run", "-t", "4", "-p", str(param_file), "--runtime-profile", "cnag-hpc", "--no-color"],
    )

    fake_param = {
        "pipeline": "wes",
        "mode": "single",
        "input_dir": None,
        "sample_map": None,
        "workflow_engine": "bash",
        "gatk_version": "gatk-4.6",
        "cleanup_bam": False,
    }

    monkeypatch.setattr(cli_mod.config_mod, "read_param_file", lambda _: fake_param)
    captured = {}

    def fake_set_config_values(params):
        captured["params"] = params
        return {
            "project_dir": str(tmp_path / "proj_run"),
            "run_id": "IDRUN",
            "profile": params["profile"],
            "genome": "b37",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "bash",
                "pipeline": "wes",
                "mode": "single",
                "gatk_version": "gatk-4.6",
                "entrypoint": "/x_run.sh",
                "config_file": None,
                "helpers": {},
            },
        }

    monkeypatch.setattr(cli_mod.config_mod, "set_config_values", fake_set_config_values)

    logs = {}
    monkeypatch.setattr(
        cli_mod, "write_log",
        lambda cfg, arg, param: logs.update({"cfg": cfg, "arg": arg, "param": param})
    )

    class FakeDNAseq:
        def __init__(self, settings):
            logs["settings"] = settings

        def variant_calling(self):
            return True

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    class FakeGoodBye:
        def say_goodbye(self):
            return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    assert cli_mod.main() == 0
    assert logs["arg"]["threads"] == 4
    assert logs["arg"]["paramfile"] == str(param_file)
    assert captured["params"]["profile"] == "cnag-hpc"
    assert logs["settings"].run_id == "IDRUN"
    assert logs["settings"].profile == "cnag-hpc"
    assert logs["settings"].workflow.entrypoint == "/x_run.sh"
    out = capsys.readouterr().out
    assert "Report" in out
    assert "HTML" in out
    assert str(tmp_path / "proj_run" / "run-report.html") in out
    assert "Bye" in out
    assert "Goodbye" not in out
    assert out.rstrip().endswith("Bye")
    assert (tmp_path / "proj_run" / "run-report.html").is_file()


def test_main_verbose_prints(monkeypatch, tmp_path, capsys):
    def fake_usage(version):
        return {
            "threads": 1,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": True,
            "nocolor": False,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)
    monkeypatch.setattr(
        cli_mod.config_mod, "read_param_file", lambda _: {
            "pipeline": "wes", "mode": "single", "input_dir": None, "sample_map": None,
            "workflow_engine": "bash", "gatk_version": "gatk-3.5", "cleanup_bam": False
        }
    )
    monkeypatch.setattr(
        cli_mod.config_mod, "set_config_values", lambda _: {
                "project_dir": str(tmp_path / "proj_verbose"),
            "id": "IDVERB",
            "genome": "b37",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "bash",
                "pipeline": "wes",
                "mode": "single",
                "gatk_version": "gatk-3.5",
                "entrypoint": "/x.sh",
                "config_file": None,
                "helpers": {},
            },
        }
    )

    class FakeDNAseq:
        def __init__(self, settings): ...
        def variant_calling(self): return True

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    class FakeGoodBye:
        def say_goodbye(self): return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    assert cli_mod.main() == 0
    out = capsys.readouterr().out
    assert "Finished successfully" in out
    assert "Elapsed" in out
    assert "Resolved Configuration" in out
    assert "Input Parameters" in out
    assert "Log" in out
    assert str(tmp_path / "proj_verbose" / "bash_wes_single_b37_gatk-3.5.log") in out


def test_main_warns_when_genome_is_inferred(monkeypatch, tmp_path, capsys):
    def fake_usage(version):
        return {
            "threads": 1,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": False,
            "nocolor": False,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)
    monkeypatch.setattr(
        cli_mod.config_mod,
        "read_param_file",
        lambda _: {
            "pipeline": "wes",
            "mode": "single",
            "input_dir": None,
            "sample_map": None,
            "workflow_engine": "bash",
            "gatk_version": "gatk-3.5",
            "cleanup_bam": False,
        },
    )
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
                "project_dir": str(tmp_path / "proj_warn"),
            "id": "IDWARN",
            "genome": "b37",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "bash",
                "pipeline": "wes",
                "mode": "single",
                "gatk_version": "gatk-3.5",
                "entrypoint": "/x.sh",
                "config_file": None,
                "helpers": {},
            },
        },
    )

    class FakeDNAseq:
        def __init__(self, settings): ...
        def variant_calling(self): return True

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    class FakeGoodBye:
        def say_goodbye(self): return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    assert cli_mod.main() == 0
    out = capsys.readouterr().out
    assert "Warning" in out
    assert "using inferred default 'b37'" in out


def test_main_partial_run_warning_and_metadata(monkeypatch, tmp_path, capsys):
    def fake_usage(version):
        return {
            "threads": 4,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": False,
            "nocolor": False,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)
    monkeypatch.setattr(
        cli_mod.config_mod,
        "read_param_file",
        lambda _: {
            "pipeline": "wes",
            "mode": "single",
            "input_dir": None,
            "sample_map": None,
            "workflow_engine": "snakemake",
            "gatk_version": "gatk-4.6",
            "cleanup_bam": False,
            "snakemake_parameters": {"target": "call_variants"},
            "nextflow_parameters": {},
        },
    )
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
            "project_dir": str(tmp_path / "proj_partial"),
            "run_id": "IDPART",
            "genome": "b37",
            "snakemake_parameters": {"target": "call_variants"},
            "nextflow_parameters": {},
            "run_mode": "partial",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "snakemake",
                "pipeline": "wes",
                "mode": "single",
                "gatk_version": "gatk-4.6",
                "entrypoint": "/x.smk",
                "config_file": "/cfg.yaml",
                "helpers": {},
            },
        },
    )

    seen = {}

    class FakeDNAseq:
        def __init__(self, settings):
            seen["settings"] = settings

        def variant_calling(self):
            return True

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    class FakeGoodBye:
        def say_goodbye(self):
            return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    assert cli_mod.main() == 0
    out = capsys.readouterr().out
    assert "Partial Run" in out
    assert "Snakemake target" in out
    assert "call_variants" in out
    assert seen["settings"].snakemake_parameters["target"] == "call_variants"
    assert seen["settings"].run_mode == "partial"


def test_run_test_command_all_selects_optional_engines_and_skips_missing(monkeypatch, tmp_path):
    script = tmp_path / "examples" / "input" / "run_tests.sh"
    script.parent.mkdir(parents=True)
    script.write_text("#!/usr/bin/env bash\n", encoding="utf-8")

    seen = {}

    def fake_run(cmd, cwd, env, check):
        seen["cmd"] = cmd
        seen["cwd"] = cwd
        seen["env"] = env
        seen["check"] = check

        class Proc:
            returncode = 0

        return Proc()

    monkeypatch.setattr(cli_mod, "_project_root", lambda: tmp_path)
    monkeypatch.setattr(cli_mod.subprocess, "run", fake_run)

    assert cli_mod._run_test_command(["--all", "-t", "2", "--runtime-profile", "cnag-hpc"]) == 0
    assert seen["cmd"] == [
        "bash",
        str(script),
        "--wes-bash",
        "--wes-snakemake",
        "--wes-nextflow",
        "--mit-bash",
    ]
    assert seen["cwd"] == str(script.parent)
    assert seen["env"]["THREADS"] == "2"
    assert seen["env"]["CBICALL_RUNTIME_PROFILE"] == "cnag-hpc"
    assert seen["env"]["CBICALL_TEST_SKIP_MISSING_OPTIONAL"] == "1"
    assert seen["check"] is False


def test_run_test_command_explicit_snakemake_requires_engine(monkeypatch, tmp_path):
    script = tmp_path / "examples" / "input" / "run_tests.sh"
    script.parent.mkdir(parents=True)
    script.write_text("#!/usr/bin/env bash\n", encoding="utf-8")

    seen = {}

    def fake_run(cmd, cwd, env, check):
        seen["cmd"] = cmd
        seen["env"] = env

        class Proc:
            returncode = 7

        return Proc()

    monkeypatch.setattr(cli_mod, "_project_root", lambda: tmp_path)
    monkeypatch.setattr(cli_mod.subprocess, "run", fake_run)

    assert cli_mod._run_test_command(["--wes-snakemake"]) == 7
    assert seen["cmd"] == ["bash", str(script), "--wes-snakemake"]
    assert "CBICALL_TEST_SKIP_MISSING_OPTIONAL" not in seen["env"]


def test_run_test_command_explicit_nextflow_requires_engine(monkeypatch, tmp_path):
    script = tmp_path / "examples" / "input" / "run_tests.sh"
    script.parent.mkdir(parents=True)
    script.write_text("#!/usr/bin/env bash\n", encoding="utf-8")

    seen = {}

    def fake_run(cmd, cwd, env, check):
        seen["cmd"] = cmd
        seen["env"] = env

        class Proc:
            returncode = 8

        return Proc()

    monkeypatch.setattr(cli_mod, "_project_root", lambda: tmp_path)
    monkeypatch.setattr(cli_mod.subprocess, "run", fake_run)

    assert cli_mod._run_test_command(["--wes-nextflow"]) == 8
    assert seen["cmd"] == ["bash", str(script), "--wes-nextflow"]
    assert "CBICALL_TEST_SKIP_MISSING_OPTIONAL" not in seen["env"]


def test_main_no_color_disables_ansi_output(monkeypatch, tmp_path, capsys):
    def fake_usage(version):
        return {
            "threads": 1,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": False,
            "nocolor": True,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)
    monkeypatch.setattr(
        cli_mod.config_mod,
        "read_param_file",
        lambda _: {
            "pipeline": "wes",
            "mode": "single",
            "input_dir": None,
            "sample_map": None,
            "workflow_engine": "bash",
            "gatk_version": "gatk-3.5",
            "cleanup_bam": False,
        },
    )
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
            "project_dir": str(tmp_path / "proj_nocolor"),
            "run_id": "IDNOCOLOR",
            "genome": "b37",
            "inputs": {"input_dir": None, "sample_map": None},
            "workflow": {
                "engine": "bash",
                "pipeline": "wes",
                "mode": "single",
                "gatk_version": "gatk-3.5",
                "entrypoint": "/x.sh",
                "config_file": None,
                "helpers": {},
            },
        },
    )

    class FakeDNAseq:
        def __init__(self, settings):
            self.settings = settings

        def variant_calling(self):
            return True

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    class FakeGoodBye:
        def say_goodbye(self):
            return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    assert cli_mod.main() == 0
    out = capsys.readouterr().out
    assert "\x1b[" not in out


def test_main_passes_wgs_cohort_workflow_keys(monkeypatch, tmp_path):
    def fake_usage(version):
        return {
            "threads": 4,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": False,
            "nocolor": False,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)

    fake_param = {
        "pipeline": "wgs",
        "mode": "cohort",
        "input_dir": None,
        "sample_map": str(tmp_path / "sample_map.tsv"),
        "workflow_engine": "bash",
        "gatk_version": "gatk-4.6",
        "cleanup_bam": False,
    }

    monkeypatch.setattr(cli_mod.config_mod, "read_param_file", lambda _: fake_param)
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
                "project_dir": str(tmp_path / "proj"),
            "id": "IDWGSCOHORT",
            "genome": "b37",
            "inputs": {"input_dir": None, "sample_map": str(tmp_path / "sample_map.tsv")},
            "workflow": {
                "engine": "bash",
                "pipeline": "wgs",
                "mode": "cohort",
                "gatk_version": "gatk-4.6",
                "entrypoint": "/x_wgs_cohort.sh",
                "config_file": "/x_config.yaml",
                "helpers": {},
            },
        },
    )

    seen = {}
    monkeypatch.setattr(cli_mod, "write_log", lambda cfg, arg, param: None)

    class FakeDNAseq:
        def __init__(self, settings):
            seen["settings"] = settings

        def variant_calling(self):
            return True

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    class FakeGoodBye:
        def say_goodbye(self):
            return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    assert cli_mod.main() == 0
    assert seen["settings"].workflow.entrypoint == "/x_wgs_cohort.sh"
    assert seen["settings"].workflow.config_file == "/x_config.yaml"
    assert seen["settings"].inputs.sample_map == str(tmp_path / "sample_map.tsv")


def test_run_analysis_passes_sarek_nextflow_settings(monkeypatch, tmp_path):
    param_file = tmp_path / "sarek.yaml"
    sample_map = tmp_path / "sarek_samplesheet.csv"
    param_file.write_text("pipeline: sarek\n", encoding="utf-8")
    sample_map.write_text("patient,sample,lane,fastq_1,fastq_2\n", encoding="utf-8")

    params = {
        "pipeline": "sarek",
        "mode": "cohort",
        "input_dir": None,
        "sample_map": str(sample_map),
        "workflow_engine": "nextflow",
        "workflow_version": "nf-core",
        "gatk_version": "nf-core",
        "resource": "nf-core-sarek-managed-resources-v1",
        "nfcore_profile": "docker",
        "nfcore_parameters": {
            "input": str(sample_map),
            "genome": "GATK.GRCh38",
            "tools": "haplotypecaller",
        },
        "cleanup_bam": False,
    }

    monkeypatch.setattr(cli_mod.config_mod, "read_param_file", lambda _: params)
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
            "project_dir": str(tmp_path / "proj_sarek"),
            "run_id": "IDSAREK",
            "workflow_engine": "nextflow",
            "profile": "local",
            "nfcore_profile": "docker",
            "genome": "external",
            "pipeline": "sarek",
            "mode": "cohort",
            "gatk_version": "nf-core",
            "pipeline_version": "v1",
            "nfcore_parameters": {
                "input": str(sample_map),
                "genome": "GATK.GRCh38",
                "tools": "haplotypecaller",
            },
            "inputs": {"input_dir": None, "sample_map": str(sample_map)},
            "workflow": {
                "engine": "nextflow",
                "pipeline": "sarek",
                "mode": "cohort",
                "gatk_version": "nf-core",
                "pipeline_version": "v1",
                "entrypoint": "nf-core/sarek",
                "config_file": None,
                "helpers": {},
                "metadata": {
                    "source_type": "nf-core",
                    "source": "nf-core/sarek",
                    "release": "3.8.1",
                    "default_outdir": "sarek",
                },
            },
            "resources": {
                "bundle": {
                    "key": "nf-core-sarek-managed-resources-v1",
                    "type": "nextflow-managed",
                    "version": "sarek-3.8.1",
                    "fingerprint": "abc",
                    "runtime_check": {"status": "not_applicable"},
                }
            },
        },
    )

    seen = {}
    monkeypatch.setattr(cli_mod, "write_log", lambda cfg, arg, param: None)

    class FakeDNAseq:
        def __init__(self, settings):
            seen["settings"] = settings

        def variant_calling(self):
            return True

    monkeypatch.setattr(cli_mod, "DNAseq", FakeDNAseq)

    class FakeGoodBye:
        def say_goodbye(self):
            return "Bye"

    monkeypatch.setattr(cli_mod, "GoodBye", FakeGoodBye)

    rc = cli_mod._run_analysis(
        {
            "threads": 2,
            "paramfile": str(param_file),
            "debug": 0,
            "verbose": False,
            "nocolor": True,
        },
        start_time=time.time(),
        cbicall_path=tmp_path / "bin" / "cbicall",
    )

    assert rc == 0
    assert seen["settings"].nfcore_profile == "docker"
    assert seen["settings"].nfcore_parameters["genome"] == "GATK.GRCh38"
    assert seen["settings"].nfcore_parameters["tools"] == "haplotypecaller"
