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
        "genome: b37\n",
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
        "genome": "b37",
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
        "genome": "b37",
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
        "genome": "b37",
    }

    with pytest.raises(RuntimeError) as excinfo:
        config_mod.set_config_values(param)
    assert "Missing +x on one or more workflow scripts" in str(excinfo.value)

def _make_executable(path: Path):
    path.write_text("#!/bin/sh\n")
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR)


def _fake_project(monkeypatch, tmp_path, gatk_ver="gatk-4.6", make_smk=True):
    """
    Create fake project root layout and monkeypatch config_mod.__file__ so
    set_config_values anchors to tmp_path as project root.
    """
    root = tmp_path
    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy")
    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))

    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True)
    for name in ["parameters.sh", "coverage.sh", "jaccard.sh", "vcf2sex.sh", "wes_single.sh", "mit_single.sh"]:
        _make_executable(bash_dir / name)

    if make_smk:
        smk_dir = root / "workflows" / "snakemake" / gatk_ver
        smk_dir.mkdir(parents=True)
        (smk_dir / "config.yaml").write_text("dummy: 1\n")
        (smk_dir / "wes_single.smk").write_text("# dummy\n")
        (smk_dir / "wgs_single.smk").write_text("# dummy\n")

    return root


def test_validate_enum_allows_none():
    # Hit: _validate_enum early-return when value is None
    config_mod._validate_enum("x", None, {"a"})


def test_read_param_file_mit_forces_rsrs_when_omitted(tmp_path):
    # Hit: pipeline=mit with genome omitted -> auto rsrs
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text(
        "mode: single\n"
        "pipeline: mit\n"
        "workflow_engine: bash\n"
        "gatk_version: gatk-3.5\n"
    )
    cfg = config_mod.read_param_file(str(yaml_path))
    assert cfg["genome"] == "rsrs"


def test_read_param_file_mit_wrong_genome_raises(tmp_path):
    # Hit: pipeline=mit with user-provided genome != rsrs -> error
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text(
        "mode: single\n"
        "pipeline: mit\n"
        "workflow_engine: bash\n"
        "gatk_version: gatk-3.5\n"
        "genome: b37\n"
    )
    with pytest.raises(ValueError, match="genome is fixed to 'rsrs'"):
        config_mod.read_param_file(str(yaml_path))


def test_read_param_file_mit_with_snakemake_raises(tmp_path):
    # Hit: mit + snakemake unsupported
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text(
        "mode: single\n"
        "pipeline: mit\n"
        "workflow_engine: snakemake\n"
        "gatk_version: gatk-3.5\n"
        "genome: rsrs\n"
    )
    with pytest.raises(ValueError, match="not supported"):
        config_mod.read_param_file(str(yaml_path))


def test_read_param_file_hg38_requires_wgs(tmp_path):
    # Hit: hg38 only supported for wgs
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_engine: bash\n"
        "gatk_version: gatk-3.5\n"
        "genome: hg38\n"
    )
    with pytest.raises(ValueError, match="only supported for pipeline='wgs'"):
        config_mod.read_param_file(str(yaml_path))


def test_set_config_values_snakemake_builds_smk_paths(monkeypatch, tmp_path):
    # Hit: snakemake branch in set_config_values, plus _validate_combos success path for gatk-4.6 + wes
    gatk_ver = "gatk-4.6"
    root = _fake_project(monkeypatch, tmp_path, gatk_ver=gatk_ver, make_smk=True)

    param = {
        "mode": "single",
        "pipeline": "wes",
        "workflow_engine": "snakemake",
        "gatk_version": gatk_ver,
        "projectdir": "proj",
        "genome": "b37",
    }

    cfg = config_mod.set_config_values(param)

    smk_dir = root / "workflows" / "snakemake" / gatk_ver
    assert cfg["smk_wes_single"] == str(smk_dir / "wes_single.smk")
    assert cfg["smk_config"] == str(smk_dir / "config.yaml")


def test_set_config_values_capture_else_branch_and_arch_else(monkeypatch, tmp_path):
    # Hit: capture else branch (pipeline != wes) + arch else (platform.machine unknown)
    gatk_ver = "gatk-3.5"
    _fake_project(monkeypatch, tmp_path, gatk_ver=gatk_ver, make_smk=False)

    monkeypatch.setattr(config_mod.platform, "machine", lambda: "weirdarch")

    param = {
        "mode": "single",
        "pipeline": "mit",           # pipeline != wes => capture else branch
        "workflow_engine": "bash",
        "gatk_version": gatk_ver,
        "projectdir": "proj",
        "genome": "rsrs",
    }

    cfg = config_mod.set_config_values(param)
    assert cfg["capture"] == "MToolBox_rsrs"
    assert cfg["arch"] == "weirdarch"


def test_set_config_values_nproc_exception_falls_back(monkeypatch, tmp_path):
    # Hit: nproc exists but os.popen fails -> fallback to os.cpu_count
    gatk_ver = "gatk-3.5"
    _fake_project(monkeypatch, tmp_path, gatk_ver=gatk_ver, make_smk=False)

    monkeypatch.setattr(config_mod.shutil, "which", lambda _: "/usr/bin/nproc")

    def boom(_):
        raise OSError("nope")

    monkeypatch.setattr(config_mod.os, "popen", boom)
    monkeypatch.setattr(config_mod.os, "cpu_count", lambda: 4)

    param = {
        "mode": "single",
        "pipeline": "wes",
        "workflow_engine": "bash",
        "gatk_version": gatk_ver,
        "projectdir": "proj",
        "genome": "b37",
    }

    cfg = config_mod.set_config_values(param)
    assert cfg["threadshost"] == 4
    assert cfg["threadsless"] == 3


