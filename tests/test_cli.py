import io
import json
import sys
import time

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


def test_write_run_report_creates_compact_summary(tmp_path):
    entrypoint = tmp_path / "wes_single.sh"
    entrypoint.write_text("#!/bin/sh\necho ok\n", encoding="utf-8")
    env_file = tmp_path / "env.sh"
    env_file.write_text("DATADIR=/tmp\n", encoding="utf-8")
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
                "bundle": {"key": "cbicall-germline-resources-v1", "version": "v1", "fingerprint": "abc"}
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
    assert inventory["entries"] == 3
    assert inventory["excluded"] == ["run-report.json"]
    assert "run-report.json" not in inventory["paths"]
    assert inventory["paths"] == ["03_stats/sample.vcf.sha256.txt", "env.sh", "wes_single.sh"]
    assert len(inventory["sha256"]) == 64
    assert data["outputs"]["vcf_hash_reports"][0]["normalized_sha256"] == "normalized"


def test_compare_runs_reports_workflow_and_output_differences(tmp_path, capsys):
    run_a = tmp_path / "run_a"
    run_b = tmp_path / "run_b"
    run_a.mkdir()
    run_b.mkdir()

    base = {
        "status": "success",
        "framework": {"name": "CBIcall", "version": "1.2.3"},
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
        "outputs": {
            "file_inventory": {"entries": 3, "sha256": "manifest-a"},
            "vcf_hash_reports": [
                {"file": "sample.vcf.gz", "normalized_sha256": "vcf"}
            ]
        },
    }
    changed = json.loads(json.dumps(base))
    changed["workflow"]["fingerprint"] = "workflow-b"
    changed["workflow"]["files"][0]["sha256"] = "bbb"
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
    assert "Workflow hash" in out
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
    assert "Legend" in report_text
    html_text = html_report.read_text(encoding="utf-8")
    assert "<title>CBIcall Run Comparison</title>" in html_text
    assert "Run Comparison" in html_text
    assert "Workflow hash" in html_text
    assert "Status summary" in html_text
    assert "class=\"pill different\"" in html_text


def test_compare_runs_accepts_multiple_runs_as_baseline_matrix(tmp_path, capsys):
    runs = [tmp_path / f"run_{idx}" for idx in range(3)]
    for run in runs:
        run.mkdir()

    base = {
        "status": "success",
        "framework": {"name": "CBIcall", "version": "1.2.3"},
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
        "outputs": {
            "file_inventory": {"entries": 3, "sha256": "manifest-a"},
            "vcf_hash_reports": [
                {"file": "sample.vcf.gz", "normalized_sha256": "vcf"}
            ]
        },
    }
    same = json.loads(json.dumps(base))
    different = json.loads(json.dumps(base))
    different["workflow"]["fingerprint"] = "workflow-c"
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
    assert "Resource ver" in out
    assert "File inventory" in out
    assert "different: " in out
    assert "run_2/run-report.json" in out
    assert "sample.vcf.gz" in out


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


def test_validate_registry_command_uses_default_registry(capsys):
    rc = cli_mod._run_validate_registry_command(["--no-color"])

    assert rc == 0
    out = capsys.readouterr().out
    assert "Registry OK" in out
    assert "cbicall-workflow-registry.yaml" in out
    assert "cbicall-workflow-registry.schema.json" in out


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


def test_main_run_subcommand_happy_path(monkeypatch, tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("pipeline: wes\nmode: single\n", encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["cbicall", "run", "-t", "4", "-p", str(param_file), "--no-color"])

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
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
            "project_dir": str(tmp_path / "proj_run"),
            "run_id": "IDRUN",
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

    assert cli_mod.main() == 0
    assert logs["arg"]["threads"] == 4
    assert logs["arg"]["paramfile"] == str(param_file)
    assert logs["settings"].run_id == "IDRUN"
    assert logs["settings"].workflow.entrypoint == "/x_run.sh"


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
            "workflow_rule": "call_variants",
            "allow_partial_run": True,
        },
    )
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {
            "project_dir": str(tmp_path / "proj_partial"),
            "run_id": "IDPART",
            "genome": "b37",
            "workflow_rule": "call_variants",
            "allow_partial_run": True,
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
    assert "Workflow rule" in out
    assert "call_variants" in out
    assert seen["settings"].workflow_rule == "call_variants"
    assert seen["settings"].allow_partial_run is True
    assert seen["settings"].run_mode == "partial"


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
