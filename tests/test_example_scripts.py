import os
import subprocess
from pathlib import Path

import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = REPO_ROOT / "examples" / "scripts"
SCRIPTS = sorted(SCRIPT_DIR.glob("*.sh"))


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda path: path.name)
def test_example_script_has_valid_bash_syntax(script):
    subprocess.run(["bash", "-n", str(script)], check=True)


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda path: path.name)
def test_example_script_is_executable_and_maintained(script):
    text = script.read_text(encoding="utf-8")
    assert os.access(script, os.X_OK)
    assert text.startswith("#!/usr/bin/env bash\n")
    assert "set -euo pipefail" in text
    assert "/media/mrueda" not in text
    assert "manuel.rueda@" not in text


@pytest.mark.parametrize(
    ("script_name", "output_name", "mode", "pipeline", "input_key", "placeholder"),
    [
        ("run_wes_single.sh", "wes_single.yaml", "single", "wes", "input_dir", "/path/to/sample01"),
        ("run_wgs_single.sh", "wgs_single.yaml", "single", "wgs", "input_dir", "/path/to/sample01"),
        ("run_wes_cohort.sh", "wes_cohort.yaml", "cohort", "wes", "sample_map", "/path/to/sample_map.tsv"),
        ("run_wgs_cohort.sh", "wgs_cohort.yaml", "cohort", "wgs", "sample_map", "/path/to/sample_map.tsv"),
        ("run_mit_single.sh", "mit_single.yaml", "single", "mit", "input_dir", "/path/to/sample01"),
        ("run_mit_cohort.sh", "mit_cohort.yaml", "cohort", "mit", "input_dir", "/path/to/cohort"),
    ],
)
def test_local_example_generates_current_yaml(
    tmp_path, script_name, output_name, mode, pipeline, input_key, placeholder
):
    fake_cbicall = tmp_path / "cbicall"
    fake_cbicall.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    fake_cbicall.chmod(0o755)

    input_dir = tmp_path / "input"
    input_dir.mkdir()
    sample_map = tmp_path / "sample_map.tsv"
    sample_map.write_text("sample01\t/path/sample01.g.vcf.gz\n", encoding="utf-8")
    input_path = sample_map if input_key == "sample_map" else input_dir
    script_copy = tmp_path / script_name
    script_copy.write_text(
        (SCRIPT_DIR / script_name)
        .read_text(encoding="utf-8")
        .replace(placeholder, str(input_path)),
        encoding="utf-8",
    )

    env = os.environ.copy()
    env["PATH"] = f"{tmp_path}:{env['PATH']}"
    subprocess.run(
        ["bash", str(script_copy)],
        check=True,
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )

    params = yaml.safe_load((tmp_path / output_name).read_text(encoding="utf-8"))
    assert params["mode"] == mode
    assert params["pipeline"] == pipeline
    assert params["workflow_provider"] == "cbicall-core"
    assert params[input_key] == str(input_path)


@pytest.mark.parametrize(
    "script_name",
    ["run_cbicall_slurm.sh", "run_cbicall_apptainer_slurm.sh"],
)
def test_slurm_example_generates_valid_job_and_yaml(tmp_path, script_name):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for command in ("sbatch",):
        executable = bin_dir / command
        executable.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
        executable.chmod(0o755)

    input_dir = tmp_path / "sample01"
    input_dir.mkdir()
    data_dir = tmp_path / "cbicall-data"
    data_dir.mkdir()
    image = tmp_path / "cbicall.sif"
    image.touch()

    script_copy = tmp_path / script_name
    script_copy.write_text(
        (SCRIPT_DIR / script_name)
        .read_text(encoding="utf-8")
        .replace("/path/to/cbicall.sif", str(image))
        .replace("/path/to/cbicall-data", str(data_dir)),
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    subprocess.run(
        ["bash", str(script_copy), "sample01", "wes", str(input_dir)],
        check=True,
        env=env,
        capture_output=True,
        text=True,
    )

    job_script = next(input_dir.glob("*.slurm"))
    subprocess.run(["bash", "-n", str(job_script)], check=True)
    params = yaml.safe_load(next(input_dir.glob("*.yaml")).read_text(encoding="utf-8"))
    assert params["workflow_provider"] == "cbicall-core"
    assert params["input_dir"] == str(input_dir)
