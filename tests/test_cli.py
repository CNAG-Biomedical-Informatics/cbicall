import io
import json
import sys
import time

import pytest

from cbicall import cli as cli_mod


def test_write_log_creates_json(tmp_path):
    cfg = {"projectdir": str(tmp_path)}
    arg = {"threads": 4}
    param = {"pipeline": "wes"}

    cli_mod.write_log(cfg, arg, param)

    log_path = tmp_path / "log.json"
    assert log_path.is_file()

    data = json.loads(log_path.read_text(encoding="utf-8"))
    assert data["arg"] == arg
    assert data["config"]["projectdir"] == cfg["projectdir"]
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


def test_print_config_hides_bash_keys_when_snakemake(capsys):
    cfg = {
        "workflow_engine": "snakemake",
        "bash_wes_single": "/tmp/wes_single.sh",
        "smk_wes_single": "/tmp/wes_single.smk",
        "x": None,
    }
    cli_mod._print_config(cfg, no_emoji=True)
    out = capsys.readouterr().out
    assert "bash_wes_single" not in out
    assert "smk_wes_single" in out
    assert "(undef)" in out


def test_main_happy_path(monkeypatch, tmp_path):
    def fake_usage(version):
        return {
            "threads": 4,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": False,
            "noemoji": True,
            "nocolor": False,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)

    fake_param = {
        "pipeline": "wes",
        "mode": "single",
        "sample": None,
        "sample_map": None,
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "cleanup_bam": False,
    }

    monkeypatch.setattr(cli_mod.config_mod, "read_param_file", lambda _: fake_param)
    monkeypatch.setattr(
        cli_mod.config_mod,
        "set_config_values",
        lambda _: {"projectdir": str(tmp_path / "proj"), "id": "ID123", "bash_wes_single": "/x.sh"},
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
    assert logs["settings"]["threads"] == 4
    assert logs["settings"]["id"] == "ID123"


def test_main_verbose_prints(monkeypatch, tmp_path, capsys):
    def fake_usage(version):
        return {
            "threads": 1,
            "paramfile": str(tmp_path / "params.yaml"),
            "debug": 0,
            "verbose": True,
            "noemoji": True,
            "nocolor": False,
        }

    monkeypatch.setattr(cli_mod, "usage", fake_usage)
    monkeypatch.setattr(
        cli_mod.config_mod, "read_param_file", lambda _: {
            "pipeline": "wes", "mode": "single", "sample": None, "sample_map": None,
            "workflow_engine": "bash", "gatk_version": "gatk-3.5", "cleanup_bam": False
        }
    )
    monkeypatch.setattr(
        cli_mod.config_mod, "set_config_values", lambda _: {
            "projectdir": str(tmp_path / "proj_verbose"),
            "id": "IDVERB",
            "bash_wes_single": "/x.sh",
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
    assert "CBICALL FINISHED OK" in out
    assert "Running time" in out
