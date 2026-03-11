import os
import sys

import pytest

from cbicall import helpmod
from cbicall import cli as cli_mod


def test_helpmod_parse_args_valid(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("mode: single\npipeline: wes\n", encoding="utf-8")

    argv = ["-t", "4", "-p", str(param_file)]
    args = helpmod.parse_args(argv, "0.0.1")

    assert args.threads == 4
    assert args.paramfile == str(param_file)
    assert args.debug is None
    assert not args.verbose


def test_helpmod_parse_args_missing_required(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("mode: single\npipeline: wes\n", encoding="utf-8")

    argv = ["-p", str(param_file)]
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_args(argv, "0.0.1")
    assert excinfo.value.code == 1


def test_helpmod_parse_args_invalid_threads(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("mode: single\npipeline: wes\n", encoding="utf-8")

    argv = ["-t", "0", "-p", str(param_file)]
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_args(argv, "0.0.1")
    assert excinfo.value.code == 1


def test_helpmod_parse_args_paramfile_must_exist(tmp_path):
    missing = tmp_path / "missing.yaml"
    argv = ["-t", "4", "-p", str(missing)]
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_args(argv, "0.0.1")
    assert excinfo.value.code == 1


def test_helpmod_parse_args_paramfile_must_be_non_empty(tmp_path):
    empty = tmp_path / "empty.yaml"
    empty.write_text("", encoding="utf-8")

    argv = ["-t", "2", "-p", str(empty)]
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_args(argv, "0.0.1")
    assert excinfo.value.code == 1


def test_cli_parse_args_wrapper(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("pipeline: wes\nmode: single\n", encoding="utf-8")

    args = cli_mod.parse_args(["-t", "2", "-p", str(param_file)])
    assert args.threads == 2
    assert args.paramfile == str(param_file)


def test_usage_sets_no_color(tmp_path, monkeypatch):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("pipeline: wes\nmode: single\n", encoding="utf-8")

    argv = ["cbicall", "-t", "1", "-p", str(param_file), "--no-color"]
    monkeypatch.setattr(sys, "argv", argv)

    os.environ.pop("ANSI_COLORS_DISABLED", None)
    args_dict = helpmod.usage("0.0.1")

    assert args_dict["nocolor"] is True
    assert os.environ.get("ANSI_COLORS_DISABLED") == "1"


def test_usage_version_exits(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["cbicall", "-v"])
    with pytest.raises(SystemExit) as excinfo:
        helpmod.usage("9.9.9")
    assert excinfo.value.code == 0
    assert "9.9.9" in capsys.readouterr().out


def test_usage_help_exits(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["cbicall", "-h"])
    with pytest.raises(SystemExit) as excinfo:
        helpmod.usage("0.0.1")
    assert excinfo.value.code == 0


def test_usage_man_exits(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["cbicall", "-man"])
    with pytest.raises(SystemExit) as excinfo:
        helpmod.usage("0.0.1")
    assert excinfo.value.code == 0
