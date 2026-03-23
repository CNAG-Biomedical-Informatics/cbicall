from pathlib import Path
import subprocess
import sys


def test_bin_cbicall_debug_keeps_traceback(tmp_path):
    param_file = tmp_path / "params.yaml"
    param_file.write_text("bad_key: value\n", encoding="utf-8")

    repo_root = Path(__file__).resolve().parents[1]
    proc = subprocess.run(
        [sys.executable, "bin/cbicall", "-debug", "1", "-t", "1", "-p", str(param_file)],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert proc.returncode != 0
    assert "Traceback" in proc.stderr
