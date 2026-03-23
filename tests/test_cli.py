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
    monkeypatch.delenv("ANSI_COLORS_DISABLED", raising=False)
    assert cli_mod._colors_enabled() is True
    assert cli_mod._code("X") == "X"

    monkeypatch.setenv("ANSI_COLORS_DISABLED", "1")
    assert cli_mod._colors_enabled() is False
    assert cli_mod._code("X") == ""


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
