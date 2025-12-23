import io
import json
import sys
import time
from pathlib import Path

import pytest

from cbicall import cli as cli_mod


def test_write_log_creates_json(tmp_path):
    cfg = {"projectdir": str(tmp_path)}
    arg = {"threads": 4}
    param = {"pipeline": "wes"}

    cli_mod.write_log(cfg, arg, param)

    log_path = tmp_path / "log.json"
    assert log_path.is_file()

    data = json.loads(log_path.read_text())
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
    # Force stdout to look like a TTY so the spinner branch is used
    class FakeStdout(io.StringIO):
        def isatty(self):
            return True

    fake_stdout = FakeStdout()
    monkeypatch.setattr(sys, "stdout", fake_stdout)

    def fn():
        # Small delay so the spinner loop has a chance to run
        time.sleep(0.01)
        return "ok"

    result = cli_mod.run_with_spinner(fn)
    assert result == "ok"

    output = fake_stdout.getvalue()
    # Spinner should have written *something* to stdout
    assert output != ""


def test_colors_enabled_and_code(monkeypatch):
    # When ANSI_COLORS_DISABLED is not set
    monkeypatch.delenv("ANSI_COLORS_DISABLED", raising=False)
    assert cli_mod._colors_enabled() is True
    s = cli_mod._code("X")
    assert s == "X"

    # When ANSI_COLORS_DISABLED is set
    monkeypatch.setenv("ANSI_COLORS_DISABLED", "1")
    assert cli_mod._colors_enabled() is False
    s2 = cli_mod._code("X")
    assert s2 == ""


def test_main_verbose_prints(monkeypatch, tmp_path, capsys):
    # Fake usage() so we fully control args, including verbose
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

    fake_param = {
        "pipeline": "wes",
        "mode": "single",
        "sample": None,
        "sample_map": None,
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "cleanup_bam": False,
    }

    def fake_read_param_file(path):
        return fake_param

    def fake_set_config_values(param):
        return {
            "projectdir": str(tmp_path / "proj_verbose"),
            "id": "IDVERB",
            "bash_wes_single": "/path/to/bash_wes_single.sh",
        }

    monkeypatch.setattr(cli_mod.config_mod, "read_param_file", fake_read_param_file)
    monkeypatch.setattr(cli_mod.config_mod, "set_config_values", fake_set_config_values)

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

    rc = cli_mod.main()
    assert rc == 0

    out = capsys.readouterr().out
    assert "CBICALL FINISHED OK" in out
    # Verbose branch should print a running time line
    assert "Running time" in out

