from pathlib import Path

import os
import stat
import pytest

from cbicall import config as config_mod


def _make_executable(path: Path):
    path.write_text("#!/bin/sh\n")
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR)


def test_read_param_file_merges_defaults_and_validates(tmp_path):
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text(
        "mode: cohort\n"
        "pipeline: wes\n"
        "organism: Homo Sapiens\n"
        "technology: Illumina HiSeq\n"
        "workflow_engine: bash\n"
        "gatk_version: gatk-3.5\n"
        "sample_map: map.txt\n"
    )

    cfg = config_mod.read_param_file(str(yaml_path))
    # Defaults overridden
    assert cfg["mode"] == "cohort"
    assert cfg["pipeline"] == "wes"
    # sample_map becomes absolute
    assert Path(cfg["sample_map"]).is_absolute()


def test_read_param_file_invalid_key_raises(tmp_path):
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text("not_a_param: 1\n")

    with pytest.raises(ValueError) as excinfo:
        config_mod.read_param_file(str(yaml_path))
    assert "does not exist" in str(excinfo.value)


def test_read_param_file_invalid_enum_raises(tmp_path):
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text(
        "mode: single\n"
        "pipeline: invalid_pipeline\n"
    )

    with pytest.raises(ValueError):
        config_mod.read_param_file(str(yaml_path))


def test_read_param_file_invalid_combo_raises(tmp_path):
    # wgs with gatk-3.5 is not allowed by _ALLOWED_COMBOS
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text(
        "mode: single\n"
        "pipeline: wgs\n"
        "gatk_version: gatk-3.5\n"
    )

    with pytest.raises(ValueError) as excinfo:
        config_mod.read_param_file(str(yaml_path))
    assert "Pipeline-mode" in str(excinfo.value)


def test_set_config_values_builds_paths_and_projectdir(tmp_path, monkeypatch):
    """
    Build a fake project root with workflows/bash/gatk-3.5 scripts
    and point config_mod.__file__ there so set_config_values sees it.
    """
    root = tmp_path

    # Fake src/cbicall/config.py so __file__ parents[2] == root
    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy")

    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))

    gatk_ver = "gatk-3.5"
    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True)

    for name in ["parameters.sh", "coverage.sh", "jaccard.sh", "vcf2sex.sh", "wes_single.sh"]:
        _make_executable(bash_dir / name)

    param = {
        "mode": "single",
        "pipeline": "wes",
        "organism": "Homo Sapiens",
        "technology": "Illumina HiSeq",
        "workflow_engine": "bash",
        "gatk_version": gatk_ver,
        "projectdir": "proj",
        "cleanup_bam": False,
        "sample": None,
        "sample_map": None,
    }

    cfg = config_mod.set_config_values(param)

    # Workflow paths correctly built
    assert cfg["bash_parameters"] == str(bash_dir / "parameters.sh")
    assert cfg["bash_coverage"] == str(bash_dir / "coverage.sh")
    assert cfg["bash_wes_single"] == str(bash_dir / "wes_single.sh")

    # Projectdir name includes projectdir, engine, pipeline, mode, gatk_version
    proj_name = Path(cfg["projectdir"]).name
    assert "proj" in proj_name
    assert "bash" in proj_name
    assert "wes" in proj_name
    assert "single" in proj_name
    assert "gatk-3.5" in proj_name

    # Thread counts are reasonable
    assert cfg["threadshost"] >= 1
    assert cfg["threadsless"] >= 1


def test_set_config_values_with_sample_uses_sample_dir(tmp_path, monkeypatch):
    root = tmp_path

    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy")
    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))

    gatk_ver = "gatk-3.5"
    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True)
    for name in ["parameters.sh", "coverage.sh", "jaccard.sh", "vcf2sex.sh", "wes_single.sh"]:
        _make_executable(bash_dir / name)

    sample_dir = tmp_path / "sampleA"
    sample_dir.mkdir()

    param = {
        "mode": "single",
        "pipeline": "wes",
        "organism": "Homo Sapiens",
        "technology": "Illumina HiSeq",
        "workflow_engine": "bash",
        "gatk_version": gatk_ver,
        "projectdir": "proj",
        "cleanup_bam": False,
        "sample": str(sample_dir),
        "sample_map": None,
    }

    cfg = config_mod.set_config_values(param)
    # Projectdir under sample dir
    assert Path(cfg["projectdir"]).parent == sample_dir
    assert cfg["output_basename"] == sample_dir.name


def test_set_config_values_raises_if_scripts_not_executable(tmp_path, monkeypatch):
    root = tmp_path

    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy")
    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))

    gatk_ver = "gatk-3.5"
    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True)

    # Create files but do not mark executable
    for name in ["parameters.sh", "coverage.sh", "jaccard.sh", "vcf2sex.sh", "wes_single.sh"]:
        (bash_dir / name).write_text("#!/bin/sh\n")

    param = {
        "mode": "single",
        "pipeline": "wes",
        "organism": "Homo Sapiens",
        "technology": "Illumina HiSeq",
        "workflow_engine": "bash",
        "gatk_version": gatk_ver,
        "projectdir": "proj",
        "cleanup_bam": False,
        "sample": None,
        "sample_map": None,
    }

    with pytest.raises(RuntimeError) as excinfo:
        config_mod.set_config_values(param)
    assert "Missing +x on one or more workflow scripts" in str(excinfo.value)

