from pathlib import Path
import stat

import pytest

from cbicall import config as config_mod


def _make_executable(path: Path):
    path.write_text("#!/bin/sh\n")
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR)


def test_read_param_file_all_enums_valid(tmp_path):
    """
    Exercise all enum fields and a gatk-4.6 / wgs / cohort combo.
    """
    yaml_path = tmp_path / "params.yaml"
    yaml_path.write_text(
        "mode: cohort\n"
        "pipeline: wgs\n"
        "organism: Mus musculus\n"
        "technology: NovaSeq\n"
        "workflow_engine: snakemake\n"
        "gatk_version: gatk-4.6\n"
        "sample_map: map.txt\n"
    )

    cfg = config_mod.read_param_file(str(yaml_path))

    assert cfg["mode"] == "cohort"
    assert cfg["pipeline"] == "wgs"
    assert cfg["organism"] == "Mus musculus"
    assert cfg["technology"] == "NovaSeq"
    assert cfg["workflow_engine"] == "snakemake"
    assert cfg["gatk_version"] == "gatk-4.6"
    # sample_map should be an absolute path
    assert Path(cfg["sample_map"]).is_absolute()


def test_set_config_values_snakemake_and_pigz(tmp_path, monkeypatch):
    """
    Hit the snakemake branch, gatk-4.6 genome/capture, pigz compression,
    and arch detection.
    """
    root = tmp_path

    # Fake src/cbicall/config.py location so parents[2] == root
    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy")
    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))

    gatk_ver = "gatk-4.6"

    # Bash workflows for gatk-4.6 (still required for exe checks)
    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True)
    for name in [
        "parameters.sh",
        "coverage.sh",
        "jaccard.sh",
        "vcf2sex.sh",
        "wgs_cohort.sh",
    ]:
        _make_executable(bash_dir / name)

    # Snakemake workflows for gatk-4.6
    smk_dir = root / "workflows" / "snakemake" / gatk_ver
    smk_dir.mkdir(parents=True)
    smk_file = smk_dir / "wgs_cohort.smk"
    smk_file.write_text("# snakemake workflow")
    smk_config = smk_dir / "config.yaml"
    smk_config.write_text("dummy: 1\n")

    # Make os.access say pigz exists, but keep normal behavior otherwise.
    real_access = config_mod.os.access

    def fake_access(path, mode):
        if str(path) == "/usr/bin/pigz":
            return True
        return real_access(path, mode)

    monkeypatch.setattr(config_mod.os, "access", fake_access)

    # Force a specific architecture branch
    monkeypatch.setattr(config_mod.platform, "machine", lambda: "x86_64")

    # Avoid depending on real nproc; just let os.cpu_count handle it here.
    monkeypatch.setattr(config_mod.shutil, "which", lambda _: None)
    monkeypatch.setattr(config_mod.os, "cpu_count", lambda: 8)

    param = {
        "mode": "cohort",
        "pipeline": "wgs",
        "organism": "Mus musculus",
        "technology": "NovaSeq",
        "workflow_engine": "snakemake",
        "gatk_version": gatk_ver,
        "projectdir": "proj",
        "cleanup_bam": True,
        "sample": None,
        "sample_map": None,
    }

    cfg = config_mod.set_config_values(param)

    # Snakemake-specific keys
    assert cfg["smk_wgs_cohort"] == str(smk_file)
    assert cfg["smk_config"] == str(smk_config)

    # Compression should use pigz
    assert cfg["zip"].startswith("/usr/bin/pigz")

    # Genome/capture for gatk-4.6 branch
    assert cfg["genome"] == "b37"

    # Architecture mapping
    assert cfg["arch"] == "x86_64"

    # Thread counts via cpu_count branch
    assert cfg["threadshost"] == 8
    assert cfg["threadsless"] == 7


def test_set_config_values_no_nproc_fallback(tmp_path, monkeypatch):
    """
    Explicitly exercise the branch where nproc is not found and
    os.cpu_count() is used for threads.
    """
    root = tmp_path

    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy")
    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))

    gatk_ver = "gatk-3.5"
    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True)
    for name in [
        "parameters.sh",
        "coverage.sh",
        "jaccard.sh",
        "vcf2sex.sh",
        "wes_single.sh",
    ]:
        _make_executable(bash_dir / name)

    # Pretend nproc is not available
    monkeypatch.setattr(config_mod.shutil, "which", lambda _: None)
    monkeypatch.setattr(config_mod.os, "cpu_count", lambda: 4)

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

    assert cfg["threadshost"] == 4
    assert cfg["threadsless"] == 3

