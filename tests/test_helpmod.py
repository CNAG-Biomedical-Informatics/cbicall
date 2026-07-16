import pytest

from cbicall import cli as cli_mod
from cbicall import helpmod


def _param_file(tmp_path):
    path = tmp_path / "params.yaml"
    path.write_text("mode: single\npipeline: wes\n", encoding="utf-8")
    return path


def test_parse_run_args_valid(tmp_path):
    param_file = _param_file(tmp_path)
    args = helpmod.parse_run_args(["-t", "4", "-p", str(param_file)], "0.0.1")

    assert args.threads == 4
    assert args.paramfile == str(param_file)
    assert args.debug is False
    assert args.verbose is False
    assert args.nocolor is False


def test_parse_run_args_accepts_boolean_execution_options(tmp_path):
    param_file = _param_file(tmp_path)
    args = helpmod.parse_run_args(
        [
            "-t",
            "4",
            "-p",
            str(param_file),
            "--debug",
            "--verbose",
            "--no-color",
        ],
        "0.0.1",
    )

    assert args.debug is True
    assert args.verbose is True
    assert args.nocolor is True


def test_parse_run_args_missing_required(tmp_path):
    param_file = _param_file(tmp_path)
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_run_args(["-p", str(param_file)], "0.0.1")
    assert excinfo.value.code == 1


def test_parse_run_args_invalid_threads(tmp_path):
    param_file = _param_file(tmp_path)
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_run_args(["-t", "0", "-p", str(param_file)], "0.0.1")
    assert excinfo.value.code == 1


def test_parse_run_args_paramfile_must_exist(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_run_args(
            ["-t", "4", "-p", str(tmp_path / "missing.yaml")],
            "0.0.1",
        )
    assert excinfo.value.code == 1


def test_parse_run_args_paramfile_must_be_non_empty(tmp_path):
    empty = tmp_path / "empty.yaml"
    empty.write_text("", encoding="utf-8")

    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_run_args(["-t", "2", "-p", str(empty)], "0.0.1")
    assert excinfo.value.code == 1


def test_cli_parse_args_requires_run_subcommand(tmp_path):
    param_file = _param_file(tmp_path)

    with pytest.raises(SystemExit, match="run subcommand is required"):
        cli_mod.parse_args(["-t", "2", "-p", str(param_file)])

    args = cli_mod.parse_args(["run", "-t", "2", "-p", str(param_file)])
    assert args.threads == 2
    assert args.paramfile == str(param_file)


@pytest.mark.parametrize(
    "removed_args",
    [
        ["-man"],
        ["-debug", "1"],
        ["-verbose"],
        ["-nc"],
    ],
)
def test_removed_run_options_are_rejected(tmp_path, removed_args):
    param_file = _param_file(tmp_path)
    argv = ["-t", "2", "-p", str(param_file), *removed_args]
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_run_args(argv, "0.0.1")
    assert excinfo.value.code == 2


def test_debug_no_longer_accepts_a_numeric_level(tmp_path):
    param_file = _param_file(tmp_path)
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_run_args(
            ["-t", "2", "-p", str(param_file), "--debug", "3"],
            "0.0.1",
        )
    assert excinfo.value.code == 2


def test_run_help_is_command_specific(capsys):
    with pytest.raises(SystemExit) as excinfo:
        helpmod.parse_run_args(["--help"], "0.0.1")
    assert excinfo.value.code == 0

    out = capsys.readouterr().out
    assert "usage: cbicall run" in out
    assert "--debug" in out
    assert "--verbose" in out
    assert "validate-registry" not in out
    assert "compare-runs" not in out


def test_top_level_help_lists_commands(capsys):
    assert helpmod.handle_main_args([], "0.0.1") == 0
    out = capsys.readouterr().out
    assert "Commands:" in out
    assert "validate-registry" in out
    assert "compare-runs" in out
    assert "--threads" not in out


def test_top_level_help_and_version_exit(capsys):
    with pytest.raises(SystemExit) as excinfo:
        helpmod.handle_main_args(["--help"], "0.0.1")
    assert excinfo.value.code == 0
    assert "Commands:" in capsys.readouterr().out

    with pytest.raises(SystemExit) as excinfo:
        helpmod.handle_main_args(["--version"], "9.9.9")
    assert excinfo.value.code == 0
    assert "9.9.9" in capsys.readouterr().out


def test_top_level_rejects_legacy_direct_run(tmp_path):
    param_file = _param_file(tmp_path)
    with pytest.raises(SystemExit) as excinfo:
        helpmod.handle_main_args(["-p", str(param_file), "-t", "2"], "0.0.1")
    assert excinfo.value.code == 2
