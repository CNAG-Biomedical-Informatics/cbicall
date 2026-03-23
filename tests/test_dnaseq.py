import pytest

from cbicall import dnaseq
from cbicall.errors import WorkflowExecutionError, WorkflowResolutionError


def test_dnaseq_builds_bash_command_with_flags(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, engine=None):
        recorded.update(
            {"cmd": cmd, "cwd": cwd, "log_path": log_path, "env": env, "engine": engine}
        )

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    script = str(tmp_path / "wes_single.sh")

    settings = {
        "project_dir": str(project_dir),
        "threads": 8,
        "run_id": "RID123",
        "debug": False,
        "cleanup_bam": True,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": "map.txt"},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "bash",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "entrypoint": script,
            "config_file": None,
            "helpers": {},
        },
    }

    obj = dnaseq.DNAseq(settings)
    assert obj.variant_calling() is True

    assert recorded["cmd"][0] == script
    assert recorded["cmd"][1:3] == ["-t", "8"]
    assert "--pipeline" in recorded["cmd"]
    assert "--cleanup-bam" in recorded["cmd"]
    assert "--sample-map" in recorded["cmd"]
    assert recorded["env"]["GENOME"] == "b37"
    assert recorded["cwd"] == project_dir
    assert recorded["engine"] == "bash"


def test_dnaseq_debug_prints(monkeypatch, tmp_path, capsys):
    def fake_run_cmd(cmd, cwd, log_path, env=None, engine=None):
        return None

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()

    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDDBG",
        "debug": True,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "bash",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-3.5",
            "entrypoint": str(tmp_path / "wes_single.sh"),
            "config_file": None,
            "helpers": {},
        },
    }

    assert dnaseq.DNAseq(settings).variant_calling() is True
    out = capsys.readouterr().out
    assert "Log file:" in out
    assert "GENOME=b37" in out


def test_dnaseq_builds_bash_command_gatk35_has_no_extra_flags(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, engine=None):
        recorded.update({"cmd": cmd})

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    script = str(tmp_path / "wes_single.sh")

    settings = {
        "project_dir": str(project_dir),
        "threads": 4,
        "run_id": "RID999",
        "debug": False,
        "cleanup_bam": True,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": "map.txt"},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "bash",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-3.5",
            "entrypoint": script,
            "config_file": None,
            "helpers": {},
        },
    }

    assert dnaseq.DNAseq(settings).variant_calling() is True
    assert recorded["cmd"] == [script, "-t", "4"]


def test_dnaseq_builds_snakemake_command_and_config(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, engine=None):
        recorded.update({"cmd": cmd, "engine": engine})

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    snakefile = str(tmp_path / "wes_single.smk")

    settings = {
        "project_dir": str(project_dir),
        "threads": 12,
        "run_id": "RID777",
        "debug": False,
        "genome": "hg38",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "snakemake",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "entrypoint": snakefile,
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {},
        },
    }

    assert dnaseq.DNAseq(settings).variant_calling() is True
    cmd = recorded["cmd"]
    assert cmd[:3] == ["snakemake", "--forceall", "all"]
    assert "--configfile" in cmd
    assert str(tmp_path / "config.yaml") in cmd
    assert "genome=hg38" in cmd
    assert "pipeline=wes" in cmd
    assert recorded["engine"] == "snakemake"


def test_dnaseq_builds_snakemake_partial_rule_command(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, engine=None):
        recorded.update({"cmd": cmd, "engine": engine})

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    snakefile = str(tmp_path / "wes_single.smk")

    settings = {
        "project_dir": str(project_dir),
        "threads": 6,
        "run_id": "RIDPART",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": "call_variants",
        "allow_partial_run": True,
        "run_mode": "partial",
        "workflow": {
            "engine": "snakemake",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "entrypoint": snakefile,
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {},
        },
    }

    assert dnaseq.DNAseq(settings).variant_calling() is True
    assert recorded["cmd"][:3] == ["snakemake", "--forceall", "call_variants"]


def test_bash_partial_run_raises(tmp_path):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()

    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDPARTBASH",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": "call_variants",
        "allow_partial_run": True,
        "run_mode": "partial",
        "workflow": {
            "engine": "bash",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "entrypoint": str(tmp_path / "wes_single.sh"),
            "config_file": None,
            "helpers": {},
        },
    }

    with pytest.raises(WorkflowResolutionError, match="Partial workflow runs are not supported"):
        dnaseq.DNAseq(settings).variant_calling()


def test_snakemake_missing_snakefile_raises(tmp_path):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()

    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDSMK",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "snakemake",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "entrypoint": None,
            "config_file": None,
            "helpers": {},
        },
    }

    with pytest.raises(RuntimeError, match="Missing Snakefile"):
        dnaseq.DNAseq(settings).variant_calling()


def test_dnaseq_raises_if_projectdir_missing(tmp_path):
    settings = {
        "project_dir": str(tmp_path / "does_not_exist"),
        "threads": 2,
        "run_id": "RID0",
        "debug": False,
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "bash",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-3.5",
            "entrypoint": str(tmp_path / "wes_single.sh"),
            "config_file": None,
            "helpers": {},
        },
    }
    with pytest.raises(RuntimeError, match="Project directory does not exist"):
        dnaseq.DNAseq(settings).variant_calling()


def test_variant_calling_raises_on_invalid_engine(tmp_path):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDBAD",
        "debug": False,
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "nope",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-3.5",
            "entrypoint": None,
            "config_file": None,
            "helpers": {},
        },
    }
    with pytest.raises(WorkflowResolutionError, match="Invalid workflow_engine"):
        dnaseq.DNAseq(settings).variant_calling()


def test_run_cmd_nonzero_return_raises(monkeypatch, tmp_path):
    class P:
        returncode = 1

    def fake_run(cmd, cwd, env, stdout, stderr, check):
        return P()

    monkeypatch.setattr(dnaseq.subprocess, "run", fake_run)
    with pytest.raises(WorkflowExecutionError, match="returncode=1") as excinfo:
        dnaseq.DNAseq._run_cmd(
            cmd=["false"],
            cwd=tmp_path,
            log_path=tmp_path / "log.txt",
            env=None,
            engine="bash",
        )
    msg = str(excinfo.value)
    assert "engine=bash" in msg
    assert "Command: false" in msg
    assert f"Working directory: {tmp_path}" in msg


def test_run_cmd_exception_raises(monkeypatch, tmp_path):
    def fake_run(cmd, cwd, env, stdout, stderr, check):
        raise OSError("boom")

    monkeypatch.setattr(dnaseq.subprocess, "run", fake_run)
    with pytest.raises(WorkflowExecutionError, match="could not start command") as excinfo:
        dnaseq.DNAseq._run_cmd(
            cmd=["echo", "ok"],
            cwd=tmp_path,
            log_path=tmp_path / "log.txt",
            env=None,
            engine="snakemake",
        )
    msg = str(excinfo.value)
    assert "engine=snakemake" in msg
    assert "Command: echo ok" in msg
