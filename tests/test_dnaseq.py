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


def test_dnaseq_passes_resolved_env_file_to_bash(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, engine=None):
        recorded.update({"env": env})

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    env_file = tmp_path / "cnag-hpc-env.sh"

    settings = {
        "project_dir": str(project_dir),
        "threads": 1,
        "run_id": "RIDENV",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "bash",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "entrypoint": str(tmp_path / "wes_single.sh"),
            "config_file": None,
            "helpers": {"env": str(env_file)},
        },
    }

    assert dnaseq.DNAseq(settings).variant_calling() is True
    assert recorded["env"]["CBICALL_ENV_FILE"] == str(env_file)


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


def test_dnaseq_builds_nextflow_command_and_helpers(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, engine=None):
        recorded.update({"cmd": cmd, "engine": engine})

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    workflow = str(tmp_path / "wes_wgs_single.nf")
    config = str(tmp_path / "config.yaml")

    settings = {
        "project_dir": str(project_dir),
        "threads": 8,
        "run_id": "RIDNF",
        "debug": False,
        "genome": "hg38",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "cleanup_bam": True,
        "workflow": {
            "engine": "nextflow",
            "pipeline": "wgs",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "entrypoint": workflow,
            "config_file": config,
            "helpers": {
                "coverage": str(tmp_path / "coverage.sh"),
                "vcf2sex": str(tmp_path / "vcf2sex.sh"),
                "vcf2hash": str(tmp_path / "vcf2hash.sh"),
            },
        },
    }

    assert dnaseq.DNAseq(settings).variant_calling() is True
    cmd = recorded["cmd"]
    assert cmd[:3] == ["nextflow", "run", workflow]
    assert "-params-file" in cmd
    assert config in cmd
    assert "--pipeline" in cmd and "wgs" in cmd
    assert "--genome" in cmd and "hg38" in cmd
    assert "--cleanup_bam" in cmd and "true" in cmd
    assert "--vcf2hash_script" in cmd
    assert recorded["engine"] == "nextflow"


def test_dnaseq_builds_nextflow_cohort_command_with_sample_map(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, engine=None):
        recorded.update({"cmd": cmd, "engine": engine})

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    sample_map = tmp_path / "sample_map.tsv"
    sample_map.write_text("S1\t/s1.g.vcf.gz\n", encoding="utf-8")

    settings = {
        "project_dir": str(project_dir),
        "threads": 4,
        "run_id": "RIDNFCOHORT",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": str(sample_map)},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "cleanup_bam": False,
        "workflow": {
            "engine": "nextflow",
            "pipeline": "wes",
            "mode": "cohort",
            "gatk_version": "gatk-4.6",
            "entrypoint": str(tmp_path / "wes_wgs_cohort.nf"),
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {"vcf2hash": str(tmp_path / "vcf2hash.sh")},
        },
    }

    assert dnaseq.DNAseq(settings).variant_calling() is True
    cmd = recorded["cmd"]
    assert "--sample_map" in cmd and str(sample_map) in cmd
    assert "--workspace" in cmd and "cohort.genomicsdb.RIDNFCOHORT" in cmd
    assert "--vcf2hash_script" in cmd


def test_dnaseq_nextflow_cohort_requires_sample_map(tmp_path):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()

    settings = {
        "project_dir": str(project_dir),
        "threads": 4,
        "run_id": "RIDNFCOHORT",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": None,
        "allow_partial_run": False,
        "run_mode": "full",
        "workflow": {
            "engine": "nextflow",
            "pipeline": "wes",
            "mode": "cohort",
            "gatk_version": "gatk-4.6",
            "entrypoint": str(tmp_path / "wes_wgs_cohort.nf"),
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {},
        },
    }

    with pytest.raises(WorkflowResolutionError, match="sample_map is required"):
        dnaseq.DNAseq(settings).variant_calling()


def test_nextflow_partial_run_raises(tmp_path):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDNFPART",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "workflow_rule": "call_variants",
        "allow_partial_run": True,
        "run_mode": "partial",
        "workflow": {
            "engine": "nextflow",
            "pipeline": "wes",
            "mode": "single",
            "gatk_version": "gatk-4.6",
            "entrypoint": str(tmp_path / "main.nf"),
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {},
        },
    }
    with pytest.raises(WorkflowResolutionError, match="nextflow"):
        dnaseq.DNAseq(settings).variant_calling()


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
