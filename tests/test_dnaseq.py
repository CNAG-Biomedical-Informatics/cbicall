from pathlib import Path

import pytest

from cbicall import dnaseq


def test_dnaseq_builds_bash_command_with_flags(tmp_path, monkeypatch):
    recorded = {}

    def fake_submit(cmd, log_path, run_id, debug):
        recorded["cmd"] = cmd
        recorded["log_path"] = log_path
        recorded["run_id"] = run_id
        recorded["debug"] = debug

    monkeypatch.setattr(dnaseq.DNAseq, "_submit_cmd", staticmethod(fake_submit))

    settings = {
        "projectdir": str(tmp_path),
        "pipeline": "wes",
        "mode": "single",
        "workflow_engine": "bash",
        "gatk_version": "gatk-4.6",
        "threads": 8,
        "id": "RUN123",
        "debug": False,
        "bash_wes_single": "/path/to/bash_wes_single.sh",
        "cleanup_bam": True,
        "sample_map": "/path/to/map.txt",
    }

    wes = dnaseq.DNAseq(settings)
    ok = wes.variant_calling()
    assert ok is True

    cmd = recorded["cmd"]
    # cd into projectdir and redirect to log
    assert cmd.startswith("cd " + str(tmp_path))
    assert "> bash_wes_single.log 2>&1" in cmd

    # core bash command
    assert "/path/to/bash_wes_single.sh -t 8" in cmd
    # pipeline flag for non-gatk-3.5
    assert "--pipeline wes" in cmd
    # cleanup flag
    assert "--cleanup-bam" in cmd
    # cohort flags when sample_map set
    assert "--sample-map /path/to/map.txt" in cmd
    assert "--workspace cohort.genomicsdb.RUN123" in cmd


def test_dnaseq_builds_snakemake_command(tmp_path, monkeypatch):
    recorded = {}

    def fake_submit(cmd, log_path, run_id, debug):
        recorded["cmd"] = cmd

    monkeypatch.setattr(dnaseq.DNAseq, "_submit_cmd", staticmethod(fake_submit))

    settings = {
        "projectdir": str(tmp_path),
        "pipeline": "mit",
        "mode": "cohort",
        "workflow_engine": "snakemake",
        "gatk_version": "gatk-4.6",
        "threads": 4,
        "id": "RUN456",
        "debug": True,
        "smk_mit_cohort": "/path/to/workflow.smk",
    }

    wes = dnaseq.DNAseq(settings)
    wes.variant_calling()

    cmd = recorded["cmd"]
    assert "snakemake --forceall all -s /path/to/workflow.smk --cores 4" in cmd
    assert "--config pipeline=mit" in cmd
    assert "> snakemake_mit_cohort.log 2>&1" in cmd


def test_dnaseq_raises_on_missing_script():
    settings = {
        "projectdir": "/tmp",
        "pipeline": "wes",
        "mode": "single",
        "workflow_engine": "bash",
        "gatk_version": "gatk-3.5",
        "threads": 1,
        "id": "RUNX",
        "debug": False,
        # Note: no bash_wes_single key
    }
    wes = dnaseq.DNAseq(settings)
    with pytest.raises(RuntimeError):
        wes.variant_calling()


def test_dnaseq_raises_on_invalid_engine():
    settings = {
        "projectdir": "/tmp",
        "pipeline": "wes",
        "mode": "single",
        "workflow_engine": "invalid_engine",
        "gatk_version": "gatk-3.5",
        "threads": 1,
        "id": "RUNX",
        "debug": False,
    }
    wes = dnaseq.DNAseq(settings)
    with pytest.raises(ValueError):
        wes.variant_calling()

