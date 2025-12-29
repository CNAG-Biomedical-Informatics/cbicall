from pathlib import Path

import pytest

from cbicall import dnaseq


def test_dnaseq_builds_bash_command_with_flags(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None):
        recorded["cmd"] = cmd
        recorded["cwd"] = cwd
        recorded["log_path"] = log_path
        recorded["env"] = env

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    projectdir = tmp_path / "proj"
    projectdir.mkdir()

    script = str(tmp_path / "wes_single.sh")

    settings = {
        "pipeline": "wes",
        "mode": "single",
        "projectdir": str(projectdir),
        "threads": 8,
        "workflow_engine": "bash",
        "gatk_version": "gatk-4.6",
        "id": "RID123",
        "debug": False,
        "cleanup_bam": True,
        "sample_map": "map.txt",
        "genome": "b37",
        "bash_wes_single": script,
    }

    obj = dnaseq.DNAseq(settings)
    assert obj.variant_calling() is True

    assert recorded["cmd"][0] == script
    assert recorded["cmd"][1:3] == ["-t", "8"]

    assert "--pipeline" in recorded["cmd"]
    assert recorded["cmd"][recorded["cmd"].index("--pipeline") + 1] == "wes"

    assert "--cleanup-bam" in recorded["cmd"]

    assert "--sample-map" in recorded["cmd"]
    assert recorded["cmd"][recorded["cmd"].index("--sample-map") + 1] == "map.txt"

    assert "--workspace" in recorded["cmd"]
    assert (
        recorded["cmd"][recorded["cmd"].index("--workspace") + 1]
        == "cohort.genomicsdb.RID123"
    )

    assert recorded["env"] is not None
    assert recorded["env"]["GENOME"] == "b37"

    assert recorded["cwd"] == projectdir
    assert recorded["log_path"] == projectdir / f"bash_wes_single_{settings['genome']}_{settings['gatk_version']}.log"


def test_dnaseq_builds_bash_command_gatk35_has_no_extra_flags(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None):
        recorded["cmd"] = cmd
        recorded["cwd"] = cwd
        recorded["log_path"] = log_path
        recorded["env"] = env

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    projectdir = tmp_path / "proj"
    projectdir.mkdir()

    script = str(tmp_path / "wes_single.sh")

    settings = {
        "pipeline": "wes",
        "mode": "single",
        "projectdir": str(projectdir),
        "threads": 4,
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "id": "RID999",
        "debug": False,
        "cleanup_bam": True,     # ignored for gatk-3.5
        "sample_map": "map.txt", # ignored for gatk-3.5
        "genome": "b37",
        "bash_wes_single": script,
    }

    obj = dnaseq.DNAseq(settings)
    assert obj.variant_calling() is True

    assert recorded["cmd"] == [script, "-t", "4"]

    assert recorded["env"] is not None
    assert recorded["env"]["GENOME"] == "b37"

    assert recorded["log_path"] == projectdir / f"bash_wes_single_{settings['genome']}_{settings['gatk_version']}.log"


def test_dnaseq_builds_snakemake_command_and_config(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None):
        recorded["cmd"] = cmd
        recorded["cwd"] = cwd
        recorded["log_path"] = log_path
        recorded["env"] = env

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    projectdir = tmp_path / "proj"
    projectdir.mkdir()

    snakefile = str(tmp_path / "wes_single.smk")

    settings = {
        "pipeline": "wes",
        "mode": "single",
        "projectdir": str(projectdir),
        "threads": 12,
        "workflow_engine": "snakemake",
        "gatk_version": "gatk-4.6",
        "id": "RID777",
        "debug": False,
        "genome": "hg38",
        "smk_wes_single": snakefile,
    }

    obj = dnaseq.DNAseq(settings)
    assert obj.variant_calling() is True

    cmd = recorded["cmd"]
    assert cmd[:3] == ["snakemake", "--forceall", "all"]
    assert "-s" in cmd and cmd[cmd.index("-s") + 1] == snakefile
    assert "--cores" in cmd and cmd[cmd.index("--cores") + 1] == "12"

    assert "--config" in cmd
    config_items = cmd[cmd.index("--config") + 1 :]
    assert "genome=hg38" in config_items
    assert "pipeline=wes" in config_items

    # we don't require env["GENOME"] for snakemake
    assert recorded["env"] is not None
    assert recorded["env"].get("GENOME") is None

    assert recorded["log_path"] == projectdir / f"snakemake_wes_single_{settings['genome']}_{settings['gatk_version']}.log"


def test_dnaseq_raises_if_projectdir_missing(tmp_path):
    settings = {
        "pipeline": "wes",
        "mode": "single",
        "projectdir": str(tmp_path / "does_not_exist"),
        "threads": 2,
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "id": "RID0",
        "debug": False,
        "bash_wes_single": str(tmp_path / "wes_single.sh"),
    }

    obj = dnaseq.DNAseq(settings)
    with pytest.raises(RuntimeError, match="Project directory does not exist"):
        obj.variant_calling()


def test_dnaseq_raises_if_missing_workflow_script_key(tmp_path):
    projectdir = tmp_path / "proj"
    projectdir.mkdir()

    settings = {
        "pipeline": "wes",
        "mode": "single",
        "projectdir": str(projectdir),
        "threads": 2,
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "id": "RID1",
        "debug": False,
        # missing bash_wes_single on purpose
    }

    obj = dnaseq.DNAseq(settings)
    with pytest.raises(RuntimeError, match="Missing bash script for pipeline/mode"):
        obj.variant_calling()

def test_cmd_to_string_with_env_overrides():
    cmd = ["echo", "hello world"]
    s = dnaseq.DNAseq._cmd_to_string(cmd, env_overrides={"GENOME": "b37"})
    assert "GENOME=b37" in s
    assert "echo" in s


def test_variant_calling_debug_prints(monkeypatch, tmp_path, capsys):
    # Don't actually run anything
    def fake_run_cmd(cmd, cwd, log_path, env=None):
        return None

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    projectdir = tmp_path / "proj"
    projectdir.mkdir()

    settings = {
        "pipeline": "wes",
        "mode": "single",
        "projectdir": str(projectdir),
        "threads": 2,
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "id": "RIDDBG",
        "debug": True,  # hit debug branch
        "genome": "b37",
        "bash_wes_single": str(tmp_path / "wes_single.sh"),
    }

    obj = dnaseq.DNAseq(settings)
    assert obj.variant_calling() is True

    out = capsys.readouterr().out
    assert "Log file:" in out
    assert "GENOME=b37" in out  # printed via _cmd_to_string env_overrides


def test_variant_calling_raises_on_invalid_engine(tmp_path):
    projectdir = tmp_path / "proj"
    projectdir.mkdir()

    settings = {
        "pipeline": "wes",
        "mode": "single",
        "projectdir": str(projectdir),
        "threads": 2,
        "workflow_engine": "nope",
        "gatk_version": "gatk-3.5",
        "id": "RIDBAD",
        "debug": False,
    }

    obj = dnaseq.DNAseq(settings)
    with pytest.raises(ValueError, match="Invalid workflow_engine"):
        obj.variant_calling()


def test_snakemake_missing_snakefile_raises(tmp_path):
    projectdir = tmp_path / "proj"
    projectdir.mkdir()

    settings = {
        "pipeline": "wes",
        "mode": "single",
        "projectdir": str(projectdir),
        "threads": 2,
        "workflow_engine": "snakemake",
        "gatk_version": "gatk-4.6",
        "id": "RIDSMK",
        "debug": False,
        "genome": "b37",
        # missing smk_wes_single on purpose
    }

    obj = dnaseq.DNAseq(settings)
    with pytest.raises(RuntimeError, match="Missing Snakefile for pipeline/mode"):
        obj.variant_calling()


def test_run_cmd_nonzero_return_raises(monkeypatch, tmp_path):
    class P:
        returncode = 1

    def fake_run(cmd, cwd, env, stdout, stderr, check):
        return P()

    monkeypatch.setattr(dnaseq.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="Failed to execute workflow"):
        dnaseq.DNAseq._run_cmd(
            cmd=["false"], cwd=tmp_path, log_path=tmp_path / "log.txt", env=None
        )


def test_run_cmd_exception_raises(monkeypatch, tmp_path):
    def fake_run(cmd, cwd, env, stdout, stderr, check):
        raise OSError("boom")

    monkeypatch.setattr(dnaseq.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="Failed to execute workflow"):
        dnaseq.DNAseq._run_cmd(
            cmd=["echo", "ok"], cwd=tmp_path, log_path=tmp_path / "log.txt", env=None
        )

