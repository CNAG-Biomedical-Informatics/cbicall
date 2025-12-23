import os
import sys
from pathlib import Path

import pytest

from cbicall import helpmod
from cbicall import cli as cli_mod


def test_helpmod_parse_args_valid(tmp_path):
    # Create a non-empty param file
    param_file = tmp_path / "params.yaml"
    param_file.write_text("mode: single\npipeline: wes\n")

    argv = ["-t", "4", "-p", str(param_file)]
    args = helpmod.parse_args(argv, "0.0.1")

    assert args.threads == 4
    assert args.paramfile == str(param_file)
    assert args.debug is None
    assert not args.verbose


def test_helpmod_parse_args_missing_required(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("mode: single\npipeline: wes\n")

    # Missing --threads
    argv = ["-p", str(param_file)]
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_args(argv, "0.0.1")
    # _validate_args_core raises SystemExit(1) on missing required args
    assert excinfo.value.code == 1


def test_helpmod_parse_args_invalid_threads(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("mode: single\npipeline: wes\n")

    # threads <= 0
    argv = ["-t", "0", "-p", str(param_file)]
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_args(argv, "0.0.1")
    assert excinfo.value.code == 1


def test_helpmod_parse_args_paramfile_must_exist(tmp_path):
    # File does not exist
    missing = tmp_path / "missing.yaml"
    argv = ["-t", "4", "-p", str(missing)]
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_args(argv, "0.0.1")
    assert excinfo.value.code == 1


def test_cli_parse_args_wrapper(tmp_path):
    # Tests cbicall.cli.parse_args, which wraps helpmod.parse_args
    param_file = tmp_path / "params.yaml"
    param_file.write_text("pipeline: wes\nmode: single\n")

    args = cli_mod.parse_args(["-t", "2", "-p", str(param_file)])
    assert args.threads == 2
    assert args.paramfile == str(param_file)


def test_usage_sets_no_color(tmp_path, monkeypatch):
    """usage() reads sys.argv directly and should set ANSI_COLORS_DISABLED."""
    param_file = tmp_path / "params.yaml"
    param_file.write_text("pipeline: wes\nmode: single\n")

    argv = ["cbicall", "-t", "1", "-p", str(param_file), "--no-color"]

    # Patch sys.argv so usage() sees our arguments
    monkeypatch.setattr(sys, "argv", argv)

    # Clear env if present
    os.environ.pop("ANSI_COLORS_DISABLED", None)

    args_dict = helpmod.usage("0.0.1")
    assert args_dict["nocolor"] is True
    assert os.environ.get("ANSI_COLORS_DISABLED") == "1"


def test_usage_version_exits(monkeypatch, capsys):
    # -v should print the version and exit(0)
    monkeypatch.setattr(sys, "argv", ["cbicall", "-v"])

    with pytest.raises(SystemExit) as excinfo:
        helpmod.usage("9.9.9")

    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert "9.9.9" in out


def test_usage_help_exits(monkeypatch, capsys):
    # -h should print help and exit(0)
    monkeypatch.setattr(sys, "argv", ["cbicall", "-h"])

    with pytest.raises(SystemExit) as excinfo:
        helpmod.usage("0.0.1")

    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    # Help text should contain the program name / description
    assert "CBICall" in out


def test_usage_man_exits(monkeypatch, capsys):
    # -man is documented as full documentation (same as help for now)
    monkeypatch.setattr(sys, "argv", ["cbicall", "-man"])

    with pytest.raises(SystemExit) as excinfo:
        helpmod.usage("0.0.1")

    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert "CBICall" in out

