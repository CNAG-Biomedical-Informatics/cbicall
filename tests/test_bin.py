from pathlib import Path
import subprocess
import sys


def test_bin_cbicall_debug_keeps_traceback(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("bad_key: value\n", encoding="utf-8")

    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [
            sys.executable,
            "bin/cbicall",
            "run",
            "--debug",
            "-t",
            "1",
            "-p",
            str(param_file),
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert proc.returncode != 0
    assert "Traceback" in proc.stderr


def test_bin_cbicall_rejects_legacy_direct_run(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("mode: single\npipeline: wes\n", encoding="utf-8")

    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [
            sys.executable,
            "bin/cbicall",
            "-t",
            "1",
            "-p",
            str(param_file),
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert proc.returncode == 2
    assert "unrecognized arguments" in proc.stderr
