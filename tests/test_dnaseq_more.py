import pytest
from cbicall import dnaseq


def test_dnaseq_snakemake_with_sample_map_adds_workspace_and_sample_map(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None):
        recorded["cmd"] = cmd
        recorded["env"] = env
        recorded["log_path"] = log_path

    monkeypatch.setattr(dnaseq.DNAseq, "_run_cmd", staticmethod(fake_run_cmd))

    projectdir = tmp_path / "proj"
    projectdir.mkdir()

    settings = {
        "pipeline": "wgs",
        "mode": "cohort",
        "projectdir": str(projectdir),
        "threads": 4,
        "workflow_engine": "snakemake",
        "gatk_version": "gatk-4.6",
        "id": "RID42",
        "debug": False,
        "genome": "hg38",
        "sample_map": "map.txt",
        "smk_wgs_cohort": str(tmp_path / "wgs_cohort.smk"),
    }

    obj = dnaseq.DNAseq(settings)
    assert obj.variant_calling() is True

    cmd = recorded["cmd"]
    # --config must include genome, pipeline, sample_map, workspace
    assert "--config" in cmd
    cfg = cmd[cmd.index("--config") + 1 :]
    assert "genome=hg38" in cfg
    assert "pipeline=wgs" in cfg
    assert "sample_map=map.txt" in cfg
    assert "workspace=cohort.genomicsdb.RID42" in cfg

    # env exists but GENOME isn't required for snakemake
    assert recorded["env"] is not None
