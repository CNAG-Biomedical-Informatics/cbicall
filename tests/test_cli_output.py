from pathlib import Path
from types import SimpleNamespace

from cbicall import cli_output


def _resolved_config(tmp_path, workflow, **overrides):
    data = {
        "display_genome": "b37",
        "genome": "b37",
        "profile": "local",
        "project_dir": str(tmp_path / "run"),
        "run_id": "RID",
        "run_mode": "full",
        "inputs": SimpleNamespace(input_dir="/tmp/input", sample_map=None),
        "workflow": workflow,
        "workflow_provider": workflow.metadata.get("provider", "cbicall"),
        "snakemake_parameters": {},
        "nfcore_profile": None,
        "nfcore_parameters": {},
        "nfcore_singularity_cache_dir": None,
    }
    data.update(overrides)
    return SimpleNamespace(**data)


def _workflow(**overrides):
    data = {
        "backend": "nextflow",
        "pipeline": "wes",
        "mode": "single",
        "software_stack": "gatk-4.6",
        "registry_version": "v1",
        "entrypoint": "/tmp/workflows/wes_single.nf",
        "config_file": "/tmp/workflows/config.yaml",
        "helpers": {},
        "metadata": {},
    }
    data.update(overrides)
    return SimpleNamespace(**data)


def test_print_run_summary_nextflow_native_outputs_entrypoint_and_config(tmp_path, capsys):
    workflow = _workflow()
    resolved = _resolved_config(tmp_path, workflow)

    cli_output._print_run_summary(
        arg={"threads": 2, "paramfile": "/tmp/params.yaml"},
        resolved_config=resolved,
        cbicall_path=Path("/tmp/bin/cbicall"),
        version="1.0",
        colors={"cyan": "", "bold": "", "reset": "", "yellow": "", "blue": ""},
    )

    out = capsys.readouterr().out
    assert "Nextflow" in out
    assert "wes_single.nf" in out
    assert "Config" in out


def test_print_run_summary_nfcore_outputs_cache_and_parameters(tmp_path, capsys):
    workflow = _workflow(
        software_stack="nf-core",
        metadata={
            "provider": "nf-core",
            "source": "nf-core/sarek",
            "release": "3.5.1",
            "default_outdir": "sarek",
        },
    )
    resolved = _resolved_config(
        tmp_path,
        workflow,
        workflow_provider="nf-core",
        nfcore_profile="test,singularity",
        nfcore_parameters={"tools": "haplotypecaller"},
        nfcore_singularity_cache_dir="/tmp/nxf-singularity-cache",
    )

    cli_output._print_run_summary(
        arg={"threads": 1, "paramfile": "/tmp/sarek.yaml"},
        resolved_config=resolved,
        cbicall_path=Path("/tmp/bin/cbicall"),
        version="1.0",
        colors={"cyan": "", "bold": "", "reset": "", "yellow": "", "blue": ""},
    )

    out = capsys.readouterr().out
    assert "NF cache" in out
    assert "nxf-singularity-cache" in out
    assert "NF parameters" in out
    assert "tools" in out


def test_print_config_returns_for_empty_mapping(capsys):
    cli_output._print_config(SimpleNamespace(to_dict=lambda: {}), "", "", "")

    assert capsys.readouterr().out == "Resolved Configuration\n"


def test_short_path_falls_back_when_home_resolution_fails(monkeypatch):
    def boom():
        raise RuntimeError("no home")

    monkeypatch.setattr(cli_output.Path, "home", staticmethod(boom))

    assert cli_output._short_path("/tmp/a/b/c/d.txt") == ".../b/c/d.txt"
