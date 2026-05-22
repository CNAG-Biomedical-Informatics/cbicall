# tests/test_config.py
import hashlib
import json
from pathlib import Path

import pytest

from cbicall import config as config_mod
from cbicall.errors import ParameterValidationError, WorkflowResolutionError
from cbicall.resources import validate_resource_catalog
from helpers import fake_project, make_executable, write_workflow_schema


# -----------------------
# Local helper utilities
# -----------------------

def _mk_fake_root(monkeypatch, tmp_path: Path) -> Path:
    """
    Minimal project root + patched config_mod.__file__ so _get_project_root works.
    """
    root = tmp_path
    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True, exist_ok=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy\n", encoding="utf-8")
    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))
    return root


def _write_registry_and_schema(root: Path, *, registry_lines: str) -> None:
    cfg_dir = root / "workflows" / "registry"
    sch_dir = root / "workflows" / "schema"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    sch_dir.mkdir(parents=True, exist_ok=True)

    (cfg_dir / "cbicall-workflow-registry.yaml").write_text(registry_lines, encoding="utf-8")
    write_workflow_schema(sch_dir / "cbicall-workflow-registry.schema.json")


def _touch_snakemake_files(root: Path, gatk_ver: str, snakefile_name: str) -> Path:
    smk_dir = root / "workflows" / "snakemake" / gatk_ver
    smk_dir.mkdir(parents=True, exist_ok=True)
    (smk_dir / snakefile_name).write_text("# dummy\n", encoding="utf-8")
    (smk_dir / "config.yaml").write_text("dummy: 1\n", encoding="utf-8")
    return smk_dir


def _touch_nextflow_files(root: Path, gatk_ver: str, workflow_name: str) -> Path:
    nf_dir = root / "workflows" / "nextflow" / gatk_ver
    nf_dir.mkdir(parents=True, exist_ok=True)
    (nf_dir / workflow_name).write_text("// dummy\n", encoding="utf-8")
    (nf_dir / "config.yaml").write_text("datadir: /tmp\n", encoding="utf-8")
    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True, exist_ok=True)
    for name in ["coverage.sh", "vcf2sex.sh", "vcf2hash.sh"]:
        p = bash_dir / name
        p.write_text("#!/bin/sh\n", encoding="utf-8")
        make_executable(p)
    return nf_dir


def _touch_bash_files(root: Path, gatk_ver: str, script_name: str, *, make_exec: bool = True) -> Path:
    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True, exist_ok=True)

    for name in ["env.sh", "coverage.sh", "jaccard.sh", "vcf2sex.sh", "vcf2hash.sh", script_name]:
        p = bash_dir / name
        p.write_text("#!/bin/sh\n", encoding="utf-8")
        if make_exec:
            make_executable(p)  # writes shebang + adds +x
    return bash_dir


def _point_env_to_datadir(root: Path, gatk_ver: str, datadir: Path) -> None:
    env = root / "workflows" / "bash" / gatk_ver / "env.sh"
    env.write_text(f"#!/bin/sh\nDATADIR={datadir}\n", encoding="utf-8")


def _minimal_bash_registry_block(gatk_ver: str = "gatk-4.6") -> str:
    # NOTE: schema requires workflows.bash to exist (required=["bash"])
    return (
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    software_stacks:\n"
        f"      {gatk_ver}:\n"
        "        helpers:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "          vcf2hash: \"vcf2hash.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single:\n"
        "              default: \"v1\"\n"
        "              versions:\n"
        "                v1:\n"
        "                  script: \"wes_single.sh\"\n"
    )


# ------------------------
# read_param_file() tests
# ------------------------

def test_read_param_file_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="Parameters file not found"):
        config_mod.read_param_file(str(tmp_path / "nope.yaml"))


def test_read_param_file_snakemake_gatk35_is_blocked(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: snakemake\n"
        "software_stack: gatk-3.5\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="snakemake.*not supported.*software_stack='gatk-3.5'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_mit_snakemake_blocked(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: mit\n"
        "workflow_backend: snakemake\n"
        "software_stack: gatk-4.6\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="pipeline='mit'.*workflow_backend='snakemake'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_mit_forces_rsrs_and_rejects_other_genome(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: mit\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-3.5\n"
        "genome: b37\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="pipeline='mit'.*genome.*fixed to 'rsrs'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_hg38_only_allowed_for_wgs(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-4.6\n"
        "genome: hg38\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="genome='hg38' is only supported for pipeline='wgs'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_empty_pipeline_version_raises(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-4.6\n"
        "pipeline_version: ''\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="pipeline_version must be a non-empty"):
        config_mod.read_param_file(str(p))


def test_read_param_file_resolves_input_dir_relative_to_yaml(tmp_path):
    sample_dir = tmp_path / "inputs" / "SAMPLE_001"
    sample_dir.mkdir(parents=True)

    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-4.6\n"
        "input_dir: inputs/SAMPLE_001\n",
        encoding="utf-8",
    )

    cfg = config_mod.read_param_file(str(p))
    assert cfg["input_dir"] == str(sample_dir.resolve())


def test_read_param_file_accepts_project_dir(tmp_path):
    p1 = tmp_path / "params_project_dir.yaml"
    p1.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-4.6\n"
        "project_dir: my_run\n",
        encoding="utf-8",
    )
    cfg1 = config_mod.read_param_file(str(p1))
    assert cfg1["project_dir"] == "my_run"


def test_read_param_file_rejects_removed_workflow_rule_key(tmp_path):
    p = tmp_path / "params_partial.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: snakemake\n"
        "software_stack: gatk-4.6\n"
        "workflow_rule: call_variants\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="workflow_rule.*does not exist"):
        config_mod.read_param_file(str(p))


def test_read_param_file_accepts_snakemake_target_parameter(tmp_path):
    p = tmp_path / "params_partial_ok.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: snakemake\n"
        "software_stack: gatk-4.6\n"
        "snakemake_parameters:\n"
        "  target: call_variants\n",
        encoding="utf-8",
    )
    cfg = config_mod.read_param_file(str(p))
    assert cfg["snakemake_parameters"]["target"] == "call_variants"


def test_read_param_file_accepts_resource(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-4.6\n"
        "resource: \"cbicall-germline-resources-v1\"\n",
        encoding="utf-8",
    )
    cfg = config_mod.read_param_file(str(p))
    assert cfg["resource"] == "cbicall-germline-resources-v1"


def test_read_param_file_rejects_removed_resource_bundle_key(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-4.6\n"
        "resource_bundle: \"cbicall-germline-resources-v1\"\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="resource_bundle.*does not exist"):
        config_mod.read_param_file(str(p))


def test_read_param_file_rejects_profile_runtime_option(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-4.6\n"
        "profile: cnag-hpc\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="use --runtime-profile"):
        config_mod.read_param_file(str(p))


# -----------------------------------------
# set_config_values() missing-registry/schema
# -----------------------------------------

def test_set_config_values_registry_missing_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Create schema only, omit registry
    sch_dir = root / "workflows" / "schema"
    sch_dir.mkdir(parents=True, exist_ok=True)
    write_workflow_schema(sch_dir / "cbicall-workflow-registry.schema.json")

    with pytest.raises(FileNotFoundError, match="Workflow registry not found"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
        )


def test_set_config_values_schema_missing_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    cfg_dir = root / "workflows" / "registry"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    (cfg_dir / "cbicall-workflow-registry.yaml").write_text("workflows:\n  bash: {}\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="Workflow schema not found"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
        )


# -----------------------------
# Enum + combo validation tests
# -----------------------------

def test_set_config_values_invalid_enum_mode_raises(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(ParameterValidationError, match="Invalid value for 'mode'"):
        config_mod.set_config_values(
            {"mode": "banana", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
        )


def test_set_config_values_partial_run_metadata(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_snakemake=True, snakemake_pipelines={"wes": {"single": "wes_single.smk"}}),
    )
    _touch_snakemake_files(root, "gatk-4.6", "wes_single.smk")

    cfg = config_mod.set_config_values(
        {
            "mode": "single",
            "pipeline": "wes",
            "workflow_backend": "snakemake",
            "software_stack": "gatk-4.6",
            "snakemake_parameters": {"target": "call_variants"},
            "nextflow_parameters": {},
        }
    )
    assert cfg["snakemake_parameters"]["target"] == "call_variants"
    assert cfg["run_mode"] == "partial"


def test_set_config_values_records_resources_bundle(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(
            include_bash=True,
            bash_pipelines={"wes": {"single": "wes_single.sh"}},
            bash_profiles={"cnag-hpc": {"env": "cnag-hpc-env.sh"}},
        ),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    resources = root / "resources"
    resources.mkdir()
    (resources / "cbicall-resource-catalog.json").write_text(
        """{
  "schema_version": 1,
  "resources": {
    "cbicall-germline-resources-v1": {
      "type": "bundle",
      "version": "v1",
      "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
      "archive": {"canonical_name": "cbicall-germline-resources-v1.tar.gz"}
    }
  }
}
""",
        encoding="utf-8",
    )

    cfg = config_mod.set_config_values(
        {
            "mode": "single",
            "pipeline": "wes",
            "workflow_backend": "bash",
            "software_stack": "gatk-4.6",
            "resource": "cbicall-germline-resources-v1",
        }
    )

    bundle = cfg["resources"]["bundle"]
    assert bundle["key"] == "cbicall-germline-resources-v1"
    assert bundle["version"] == "v1"
    assert bundle["compatible"] is True
    assert bundle["workflow_key"] == "bash/wes/single/gatk-4.6/v1"
    assert len(bundle["fingerprint"]) == 64
    assert "entry" not in bundle


def test_set_config_values_requires_versioned_resource_compatibility(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs={"bash_pipelines": {"wes": {"single": "wes_single.sh"}}},
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    resources = root / "resources"
    resources.mkdir()
    (resources / "cbicall-resource-catalog.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "cbicall-germline-resources-v1": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6"],
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError, match="bash/wes/single/gatk-4.6/v1"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6"}
        )


def test_validate_resource_catalog_accepts_registry_workflow_keys(tmp_path):
    catalog = tmp_path / "cbicall-resource-catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "cbicall-germline-resources-v1": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
                        "remote_identifier": {
                            "sha256": "a" * 64,
                            "expected": {"resource_key": "cbicall-germline-resources-v1"},
                        },
                        "archive": {"checksum_algorithm": "md5"},
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    registry = {
        "workflows": {
            "bash": {
                "software_stacks": {
                    "gatk-4.6": {
                        "pipelines": {
                            "wes": {
                                "single": {
                                    "default": "v1",
                                    "versions": {"v1": {"script": "wes_single.sh"}},
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    summary = validate_resource_catalog(catalog, registry)

    assert summary["bundle_resources"] == 1
    assert summary["compatible_workflows"] == 1


def test_validate_resource_catalog_rejects_unknown_workflow_key(tmp_path):
    catalog = tmp_path / "cbicall-resource-catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "cbicall-germline-resources-v1": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wgs/single/gatk-4.6/v1"],
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    registry = {
        "workflows": {
            "bash": {
                "software_stacks": {
                    "gatk-4.6": {
                        "pipelines": {
                            "wes": {
                                "single": {
                                    "default": "v1",
                                    "versions": {"v1": {"script": "wes_single.sh"}},
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    with pytest.raises(ParameterValidationError, match="not defined in the workflow registry"):
        validate_resource_catalog(catalog, registry)


def test_validate_resource_catalog_rejects_bad_identifier_metadata(tmp_path):
    catalog = tmp_path / "cbicall-resource-catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "cbicall-germline-resources-v1": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
                        "remote_identifier": {
                            "sha256": "not-a-sha",
                            "expected": {"resource_key": "other-bundle"},
                        },
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError, match="sha256.*does not match"):
        validate_resource_catalog(catalog)


def test_set_config_values_verifies_installed_resource_identifier(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs={"bash_pipelines": {"wes": {"single": "wes_single.sh"}}},
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    datadir = tmp_path / "cbicall-data"
    datadir.mkdir()
    _point_env_to_datadir(root, "gatk-4.6", datadir)

    identifier_payload = b'{"resource_key": "cbicall-germline-resources-v1"}\n'
    (datadir / "cbicall-resource-id.json").write_bytes(identifier_payload)
    identifier_sha256 = hashlib.sha256(identifier_payload).hexdigest()

    resources = root / "resources"
    resources.mkdir()
    (resources / "cbicall-resource-catalog.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "cbicall-germline-resources-v1": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
                        "remote_identifier": {"sha256": identifier_sha256},
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6"}
    )

    runtime_check = cfg["resources"]["bundle"]["runtime_check"]
    assert runtime_check["status"] == "verified"
    assert runtime_check["checks"][0]["type"] == "resource_identifier"
    assert runtime_check["checks"][0]["sha256"] == identifier_sha256


def test_set_config_values_verifies_installation_manifest_fingerprint(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs={"bash_pipelines": {"wes": {"single": "wes_single.sh"}}},
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    datadir = tmp_path / "cbicall-data"
    datadir.mkdir()
    _point_env_to_datadir(root, "gatk-4.6", datadir)

    catalog_entry = {
        "key": "cbicall-germline-resources-v1",
        "type": "bundle",
        "version": "v1",
        "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
    }
    (datadir / "cbicall-resource-installation.json").write_text(
        json.dumps({"catalog_entry": catalog_entry}),
        encoding="utf-8",
    )

    resources = root / "resources"
    resources.mkdir()
    (resources / "cbicall-resource-catalog.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "cbicall-germline-resources-v1": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6"}
    )

    runtime_check = cfg["resources"]["bundle"]["runtime_check"]
    assert runtime_check["status"] == "verified"
    assert runtime_check["checks"][0]["type"] == "installation_manifest"
    assert len(runtime_check["checks"][0]["fingerprint"]) == 64


def test_set_config_values_rejects_mismatched_installed_resource_identifier(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs={"bash_pipelines": {"wes": {"single": "wes_single.sh"}}},
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    datadir = tmp_path / "cbicall-data"
    datadir.mkdir()
    _point_env_to_datadir(root, "gatk-4.6", datadir)
    (datadir / "cbicall-resource-id.json").write_text('{"resource_key": "other-bundle"}\n', encoding="utf-8")

    with pytest.raises(ParameterValidationError, match="identifier does not match"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6"}
        )


def test_set_config_values_rejects_mismatched_installation_manifest_fingerprint(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs={"bash_pipelines": {"wes": {"single": "wes_single.sh"}}},
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    datadir = tmp_path / "cbicall-data"
    datadir.mkdir()
    _point_env_to_datadir(root, "gatk-4.6", datadir)
    (datadir / "cbicall-resource-installation.json").write_text(
        json.dumps(
            {
                "catalog_entry": {
                    "key": "cbicall-germline-resources-v1",
                    "description": "old",
                }
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError, match="manifest fingerprint does not match"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6"}
        )


def test_set_config_values_uses_registry_default_pipeline_version(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs={"bash_pipelines": {"wes": {"single": "wes_single.sh"}}},
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6"}
    )

    assert cfg["pipeline_version"] == "v1"
    assert cfg["workflow"]["pipeline_version"] == "v1"


def test_set_config_values_accepts_explicit_pipeline_version(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs={"bash_pipelines": {"wes": {"single": "wes_single.sh"}}},
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    cfg = config_mod.set_config_values(
        {
            "mode": "single",
            "pipeline": "wes",
            "workflow_backend": "bash",
            "software_stack": "gatk-4.6",
            "pipeline_version": "v1",
        }
    )

    assert cfg["workflow"]["pipeline_version"] == "v1"


def test_set_config_values_unknown_pipeline_version_raises(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs={"bash_pipelines": {"wes": {"single": "wes_single.sh"}}},
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="pipeline_version='v2' is not defined"):
        config_mod.set_config_values(
            {
                "mode": "single",
                "pipeline": "wes",
                "workflow_backend": "bash",
                "software_stack": "gatk-4.6",
                "pipeline_version": "v2",
            }
        )


def test_set_config_values_rejects_unknown_bundle_resource(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(
            include_bash=True,
            bash_pipelines={"wes": {"single": "wes_single.sh"}},
            bash_profiles={"cnag-hpc": {"env": "cnag-hpc-env.sh"}},
        ),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    resources = root / "resources"
    resources.mkdir()
    (resources / "cbicall-resource-catalog.json").write_text(
        '{"schema_version": 1, "resources": {}}',
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError, match="Resource key 'missing-bundle'"):
        config_mod.set_config_values(
            {
                "mode": "single",
                "pipeline": "wes",
                "workflow_backend": "bash",
                "software_stack": "gatk-4.6",
                "resource": "missing-bundle",
            }
        )


def test_set_config_values_invalid_combo_raises(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    # wgs is not allowed for gatk-3.5 in _ALLOWED_COMBOS
    with pytest.raises(ParameterValidationError, match="not supported for software_stack"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wgs", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
        )


# -----------------------------
# Schema validation error tests
# -----------------------------

def test_load_workflow_registry_schema_validation_fails(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Valid schema, invalid registry (missing required top-level "workflows")
    bad_registry = "not_workflows:\n  x: 1\n"
    _write_registry_and_schema(root, registry_lines=bad_registry)

    registry_yaml = root / "workflows" / "registry" / "cbicall-workflow-registry.yaml"
    schema_json = root / "workflows" / "schema" / "cbicall-workflow-registry.schema.json"

    with pytest.raises(ParameterValidationError, match="failed schema validation"):
        config_mod.load_workflow_registry(registry_yaml, schema_json)


def test__validate_enum_allows_none():
    # hits the early-return branch in _validate_enum
    config_mod._validate_enum("anything", None, {"a", "b"})


def test__validate_with_schema_formats_root_and_nested_locations():
    # Schema that expects an object with a nested requirement
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "required": ["workflows"],
        "properties": {
            "workflows": {
                "type": "object",
                "required": ["bash"],
                "properties": {"bash": {"type": "object"}},
            }
        },
        "additionalProperties": False,
    }

    # 1) ROOT error: wrong top-level type -> e.path empty -> loc "(root)"
    with pytest.raises(ParameterValidationError) as excinfo1:
        config_mod._validate_with_schema(data="not-an-object", schema=schema, label="X")
    msg1 = str(excinfo1.value)
    assert "failed schema validation" in msg1
    assert "(root)" in msg1

    # 2) NESTED error: missing required nested key -> e.path is "workflows"
    with pytest.raises(ParameterValidationError) as excinfo2:
        config_mod._validate_with_schema(data={"workflows": {}}, schema=schema, label="Y")
    msg2 = str(excinfo2.value)
    assert "failed schema validation" in msg2
    assert "workflows" in msg2


# ---------------------------------------------
# Registry wiring errors (engine/version/pipeline/mode)
# ---------------------------------------------

def test_set_config_values_backend_not_in_registry_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Registry defines only bash, request snakemake
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    software_stacks:\n"
        "      gatk-4.6:\n"
        "        helpers:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "          vcf2hash: \"vcf2hash.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single:\n"
        "              default: \"v1\"\n"
        "              versions:\n"
        "                v1:\n"
        "                  script: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="Backend not defined in workflow registry"):
        config_mod.set_config_values(
            {
                "mode": "single",
                "pipeline": "wes",
                "workflow_backend": "snakemake",  # not present in registry
                "software_stack": "gatk-4.6",      # avoids snakemake+3.5 guardrail
                "genome": "b37",
            }
        )


def test_set_config_values_version_not_in_registry_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Only gatk-3.5 exists in registry; request gatk-4.6
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    software_stacks:\n"
        "      gatk-3.5:\n"
        "        helpers:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "          vcf2hash: \"vcf2hash.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single:\n"
        "              default: \"v1\"\n"
        "              versions:\n"
        "                v1:\n"
        "                  script: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="Software stack not defined"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6"}
        )


def test_set_config_values_pipeline_not_in_registry_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Registry defines only wes; request mit
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    software_stacks:\n"
        "      gatk-3.5:\n"
        "        helpers:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "          vcf2hash: \"vcf2hash.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single:\n"
        "              default: \"v1\"\n"
        "              versions:\n"
        "                v1:\n"
        "                  script: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="Pipeline not defined"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "mit", "workflow_backend": "bash", "software_stack": "gatk-3.5", "genome": "rsrs"}
        )


def test_set_config_values_mode_not_in_registry_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # wes only defines single; request cohort
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    software_stacks:\n"
        "      gatk-3.5:\n"
        "        helpers:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "          vcf2hash: \"vcf2hash.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single:\n"
        "              default: \"v1\"\n"
        "              versions:\n"
        "                v1:\n"
        "                  script: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="Mode not defined"):
        config_mod.set_config_values(
            {"mode": "cohort", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
        )


# --------------------------------
# Missing helper keys branches
# --------------------------------

def test_set_config_values_bash_missing_helper_keys_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Omit coverage/jaccard/vcf2sex/vcf2hash in helpers
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    software_stacks:\n"
        "      gatk-3.5:\n"
        "        helpers:\n"
        "          env: \"env.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single:\n"
        "              default: \"v1\"\n"
        "              versions:\n"
        "                v1:\n"
        "                  script: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="missing helper keys for bash"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
        )


def test_set_config_values_snakemake_missing_helper_config_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Include minimal bash section to satisfy schema required=["bash"]
    reg = (
        "workflows:\n"
        + _minimal_bash_registry_block("gatk-4.6")
        + "  snakemake:\n"
          "    base_dir: \"workflows/snakemake\"\n"
          "    software_stacks:\n"
          "      gatk-4.6:\n"
          "        helpers:\n"
          "          something_else: \"x.yaml\"\n"
          "        pipelines:\n"
          "          wgs:\n"
          "            cohort:\n"
          "              default: \"v1\"\n"
          "              versions:\n"
          "                v1:\n"
          "                  script: \"wgs_cohort.smk\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    # satisfy bash files referenced by schema-required bash block
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    # create snakemake pipeline file so we pass "missing files" guardrail
    smk_dir = root / "workflows" / "snakemake" / "gatk-4.6"
    smk_dir.mkdir(parents=True, exist_ok=True)
    (smk_dir / "wgs_cohort.smk").write_text("# dummy\n", encoding="utf-8")

    with pytest.raises(WorkflowResolutionError, match="missing helper key 'config'"):
        config_mod.set_config_values(
            {"mode": "cohort", "pipeline": "wgs", "workflow_backend": "snakemake", "software_stack": "gatk-4.6"}
        )


def test_set_config_values_snakemake_missing_files_triggers_guard(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        + _minimal_bash_registry_block("gatk-4.6")
        + "  snakemake:\n"
          "    base_dir: \"workflows/snakemake\"\n"
          "    software_stacks:\n"
          "      gatk-4.6:\n"
          "        helpers:\n"
          "          config: \"config.yaml\"\n"
          "        pipelines:\n"
          "          wgs:\n"
          "            cohort:\n"
          "              default: \"v1\"\n"
          "              versions:\n"
          "                v1:\n"
          "                  script: \"wgs_cohort.smk\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    # satisfy bash files referenced by schema-required bash block
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    # Intentionally do NOT create:
    #   workflows/snakemake/gatk-4.6/wgs_cohort.smk
    #   workflows/snakemake/gatk-4.6/config.yaml
    # so the snakemake must_exist / missing_files guard triggers.
    with pytest.raises(WorkflowResolutionError, match="referenced in the registry do not exist"):
        config_mod.set_config_values(
            {"mode": "cohort", "pipeline": "wgs", "workflow_backend": "snakemake", "software_stack": "gatk-4.6"}
        )


# --------------------------
# nextflow branch
# --------------------------

def test_set_config_values_nextflow_single_gatk46_resolves(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        + _minimal_bash_registry_block("gatk-4.6")
        + "  nextflow:\n"
          "    base_dir: \"workflows/nextflow\"\n"
          "    software_stacks:\n"
          "      gatk-4.6:\n"
          "        helpers:\n"
          "          config: \"config.yaml\"\n"
          "          coverage: \"../../bash/gatk-4.6/coverage.sh\"\n"
          "          vcf2sex: \"../../bash/gatk-4.6/vcf2sex.sh\"\n"
          "          vcf2hash: \"../../bash/gatk-4.6/vcf2hash.sh\"\n"
          "        pipelines:\n"
          "          wgs:\n"
          "            single:\n"
          "              default: \"v1\"\n"
          "              versions:\n"
          "                v1:\n"
          "                  script: \"wes_wgs_single.nf\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    _touch_nextflow_files(root, "gatk-4.6", "wes_wgs_single.nf")

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wgs", "workflow_backend": "nextflow", "software_stack": "gatk-4.6", "genome": "hg38"}
    )

    assert cfg["workflow"]["backend"] == "nextflow"
    assert cfg["workflow"]["entrypoint"].endswith("wes_wgs_single.nf")
    assert cfg["workflow"]["config_file"].endswith("config.yaml")
    assert cfg["workflow"]["helpers"]["vcf2hash"].endswith("vcf2hash.sh")


def test_set_config_values_nextflow_cohort_gatk46_resolves(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        + _minimal_bash_registry_block("gatk-4.6")
        + "  nextflow:\n"
          "    base_dir: \"workflows/nextflow\"\n"
          "    software_stacks:\n"
          "      gatk-4.6:\n"
          "        helpers:\n"
          "          config: \"config.yaml\"\n"
          "          coverage: \"../../bash/gatk-4.6/coverage.sh\"\n"
          "          vcf2sex: \"../../bash/gatk-4.6/vcf2sex.sh\"\n"
          "          vcf2hash: \"../../bash/gatk-4.6/vcf2hash.sh\"\n"
          "        pipelines:\n"
          "          wgs:\n"
          "            cohort:\n"
          "              default: \"v1\"\n"
          "              versions:\n"
          "                v1:\n"
          "                  script: \"wes_wgs_cohort.nf\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    _touch_nextflow_files(root, "gatk-4.6", "wes_wgs_cohort.nf")

    cfg = config_mod.set_config_values(
        {"mode": "cohort", "pipeline": "wgs", "workflow_backend": "nextflow", "software_stack": "gatk-4.6", "genome": "hg38"}
    )

    assert cfg["workflow"]["backend"] == "nextflow"
    assert cfg["workflow"]["entrypoint"].endswith("wes_wgs_cohort.nf")


def test_read_param_file_accepts_sarek_external_nextflow(tmp_path):
    sample_map = tmp_path / "samplesheet.csv"
    sample_map.write_text("patient,sample,lane,fastq_1,fastq_2\n", encoding="utf-8")
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: cohort\n"
        "pipeline: sarek\n"
        "workflow_backend: nextflow\n"
        "workflow_provider: nf-core\n"
        "resource: nf-core-sarek-managed-resources-v1\n"
        "nfcore_profile: docker\n"
        "nfcore_singularity_cache_dir: nxf-cache\n"
        "nfcore_parameters:\n"
        "  input: samplesheet.csv\n"
        "  genome: GATK.GRCh38\n"
        "  tools: haplotypecaller\n",
        encoding="utf-8",
    )

    cfg = config_mod.read_param_file(str(p))
    assert cfg["software_stack"] == "nf-core"
    assert cfg["genome"] == "external"
    assert cfg["nfcore_singularity_cache_dir"] == str((tmp_path / "nxf-cache").resolve())
    assert cfg["nfcore_parameters"]["input"] == str(sample_map.resolve())
    assert cfg["nfcore_parameters"]["genome"] == "GATK.GRCh38"


def test_read_param_file_sarek_requires_nfcore_profile(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: cohort\n"
        "pipeline: sarek\n"
        "workflow_backend: nextflow\n"
        "workflow_provider: nf-core\n"
        "resource: nf-core-sarek-managed-resources-v1\n"
        "nfcore_parameters:\n"
        "  input: samplesheet.csv\n"
        "  genome: GATK.GRCh38\n"
        "  tools: haplotypecaller\n",
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError, match="nfcore_profile"):
        config_mod.read_param_file(str(p))


def test_read_param_file_rejects_nfcore_as_software_stack(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: cohort\n"
        "pipeline: wes\n"
        "workflow_backend: nextflow\n"
        "software_stack: nf-core\n"
        "resource: nf-core-sarek-managed-resources-v1\n",
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError, match="Invalid value for 'software_stack'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_rejects_nfcore_cache_dir_for_native_workflow(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_backend: bash\n"
        "software_stack: gatk-4.6\n"
        "nfcore_singularity_cache_dir: nxf-cache\n",
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError, match="nfcore_singularity_cache_dir"):
        config_mod.read_param_file(str(p))


def test_set_config_values_sarek_resolves_external_nfcore_workflow(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        + _minimal_bash_registry_block("gatk-4.6")
        + "  nextflow:\n"
          "    base_dir: \"workflows/nextflow\"\n"
          "    software_stacks:\n"
          "      nf-core:\n"
          "        pipelines:\n"
          "          sarek:\n"
          "            cohort:\n"
          "              default: \"v1\"\n"
          "              versions:\n"
          "                v1:\n"
          "                  provider: \"nf-core\"\n"
          "                  source: \"nf-core/sarek\"\n"
          "                  release: \"3.8.1\"\n"
          "                  default_outdir: \"sarek\"\n"
          "                  canonical_outputs:\n"
          "                    - name: \"haplotypecaller_vcf\"\n"
          "                      type: \"vcf\"\n"
          "                      pattern: \"variant_calling/haplotypecaller/*/*.haplotypecaller.vcf.gz\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    resources = root / "resources"
    resources.mkdir()
    (resources / "cbicall-resource-catalog.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "nf-core-sarek-managed-resources-v1": {
                        "type": "nextflow-managed",
                        "version": "sarek-3.8.1",
                        "compatible_workflows": ["nextflow/sarek/cohort/nf-core/v1"],
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    cfg = config_mod.set_config_values(
        {
            "mode": "cohort",
            "pipeline": "sarek",
            "workflow_backend": "nextflow",
            "workflow_provider": "nf-core",
            "resource": "nf-core-sarek-managed-resources-v1",
            "nfcore_profile": "docker",
            "nfcore_singularity_cache_dir": str(tmp_path / "nxf-cache"),
            "nfcore_parameters": {
                "input": str(tmp_path / "samplesheet.csv"),
                "genome": "GATK.GRCh38",
                "tools": "haplotypecaller",
            },
        }
    )

    assert cfg["workflow"]["entrypoint"] == "nf-core/sarek"
    assert cfg["workflow"]["metadata"]["provider"] == "nf-core"
    assert cfg["workflow"]["metadata"]["release"] == "3.8.1"
    assert cfg["workflow"]["metadata"]["canonical_outputs"][0]["name"] == "haplotypecaller_vcf"
    assert cfg["nfcore_profile"] == "docker"
    assert cfg["nfcore_singularity_cache_dir"] == str(tmp_path / "nxf-cache")
    assert cfg["nfcore_parameters"]["tools"] == "haplotypecaller"
    assert cfg["resources"]["bundle"]["type"] == "nextflow-managed"
    assert cfg["resources"]["bundle"]["runtime_check"]["status"] == "not_applicable"
    assert cfg["capture_label"] is None
    run_name = Path(cfg["project_dir"]).name
    assert run_name.startswith("cbicall_nextflow_nf-core_sarek_cohort_GATK.GRCh38_")
    assert "external" not in run_name


# --------------------------
# “Normal” branches
# --------------------------

def test_set_config_values_sample_sets_output_basename(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    sample_dir = tmp_path / "SAMPLE_001"
    sample_dir.mkdir()

    cfg = config_mod.set_config_values(
        {
            "mode": "single",
            "pipeline": "wes",
            "workflow_backend": "bash",
            "software_stack": "gatk-3.5",
            "input_dir": str(sample_dir),
        }
    )

    assert cfg["output_basename"] == "SAMPLE_001"
    assert cfg["project_dir"].startswith(str(sample_dir))


def test_set_config_values_builds_normalized_workflow_structure(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6"}
    )

    assert cfg["workflow"]["backend"] == "bash"
    assert cfg["workflow"]["pipeline"] == "wes"
    assert cfg["workflow"]["mode"] == "single"
    assert cfg["workflow"]["software_stack"] == "gatk-4.6"
    assert cfg["workflow"]["entrypoint"].endswith("workflows/bash/gatk-4.6/wes_single.sh")
    assert cfg["workflow"]["helpers"]["env"].endswith("workflows/bash/gatk-4.6/env.sh")
    assert "profiles" not in cfg["workflow"]
    assert cfg["inputs"] == {"input_dir": None, "sample_map": None}


def test_set_config_values_applies_cnag_hpc_profile(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(
            include_bash=True,
            bash_pipelines={"wes": {"single": "wes_single.sh"}},
            bash_profiles={"cnag-hpc": {"env": "cnag-hpc-env.sh"}},
        ),
    )
    bash_dir = _touch_bash_files(root, "gatk-4.6", "wes_single.sh")
    cnag_env = bash_dir / "cnag-hpc-env.sh"
    cnag_env.write_text("#!/bin/sh\n", encoding="utf-8")
    make_executable(cnag_env)

    cfg = config_mod.set_config_values(
        {
            "mode": "single",
            "pipeline": "wes",
            "workflow_backend": "bash",
            "software_stack": "gatk-4.6",
            "profile": "cnag-hpc",
        }
    )

    assert cfg["profile"] == "cnag-hpc"
    assert cfg["workflow"]["helpers"]["env"].endswith("workflows/bash/gatk-4.6/cnag-hpc-env.sh")


def test_set_config_values_undeclared_profile_raises(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="not declared"):
        config_mod.set_config_values(
            {
                "mode": "single",
                "pipeline": "wes",
                "workflow_backend": "bash",
                "software_stack": "gatk-4.6",
                "profile": "cnag-hpc",
            }
        )


def test_set_config_values_capture_branches(monkeypatch, tmp_path):
    # 1) gatk-3.5 + wes (non-wgs) -> "Agilent SureSelect"
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")
    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
    )
    assert cfg["capture_label"] == "Agilent SureSelect"

    # 2) mit -> "MToolBox_rsrs" (force x86_64 so MIT isn't blocked)
    root2 = fake_project(
        monkeypatch,
        tmp_path / "root2",
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"mit": {"single": "mit_single.sh"}}),
    )
    _touch_bash_files(root2, "gatk-3.5", "mit_single.sh")
    monkeypatch.setattr(config_mod.platform, "machine", lambda: "x86_64")

    cfg2 = config_mod.set_config_values(
        {"mode": "single", "pipeline": "mit", "workflow_backend": "bash", "software_stack": "gatk-3.5", "genome": "rsrs"}
    )
    assert cfg2["capture_label"] == "MToolBox_rsrs"

    # 3) gatk-4.6 + wes (non-wgs) -> "GATK_bundle_<genome>"
    root3 = fake_project(
        monkeypatch,
        tmp_path / "root3",
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root3, "gatk-4.6", "wes_single.sh")
    cfg3 = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6", "genome": "b37"}
    )
    assert cfg3["capture_label"] == "GATK_bundle_b37"


def test_set_config_values_mit_on_arm64_raises(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"mit": {"single": "mit_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "mit_single.sh")

    monkeypatch.setattr(config_mod.platform, "machine", lambda: "aarch64")

    with pytest.raises(WorkflowResolutionError, match=r"mit_single cannot be performed"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "mit", "workflow_backend": "bash", "software_stack": "gatk-3.5", "genome": "rsrs"}
        )


def test_set_config_values_nproc_success_path(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    class FakePopen:
        def read(self):
            return "16\n"

    monkeypatch.setattr(config_mod.shutil, "which", lambda _: "/usr/bin/nproc")
    monkeypatch.setattr(config_mod.os, "popen", lambda _: FakePopen())

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
    )
    assert cfg["host_threads"] == 16
    assert cfg["host_threads_minus_one"] == 15


def test_set_config_values_getpass_fallback_user_unknown(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    monkeypatch.delenv("LOGNAME", raising=False)
    monkeypatch.delenv("USER", raising=False)
    monkeypatch.setattr(
        config_mod.getpass,
        "getuser",
        lambda: (_ for _ in ()).throw(Exception("boom")),
    )

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
    )
    assert cfg["user"] == "unknown"


def test_set_config_values_arch_aarch64_maps_to_arm64(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    monkeypatch.setattr(config_mod.platform, "machine", lambda: "aarch64")

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6", "genome": "b37"}
    )
    assert cfg["arch"] == "arm64"


# --------------------------
# Thread/pigz/arch extra branches
# --------------------------

def test_set_config_values_nproc_read_fails_falls_back_to_cpu_count(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    monkeypatch.setattr(config_mod.shutil, "which", lambda _: "/usr/bin/nproc")

    class BadPopen:
        def read(self):
            raise OSError("boom")

    monkeypatch.setattr(config_mod.os, "popen", lambda _: BadPopen())
    monkeypatch.setattr(config_mod.os, "cpu_count", lambda: 7)

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
    )
    assert cfg["host_threads"] == 7
    assert cfg["host_threads_minus_one"] == 6


def test_set_config_values_nproc_not_found_uses_cpu_count(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    monkeypatch.setattr(config_mod.shutil, "which", lambda _: None)
    monkeypatch.setattr(config_mod.os, "cpu_count", lambda: 5)

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
    )
    assert cfg["host_threads"] == 5
    assert cfg["host_threads_minus_one"] == 4


def test_set_config_values_pigz_missing_sets_gunzip(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    real_access = config_mod.os.access

    def fake_access(path, mode):
        # Only pretend pigz doesn't exist; keep real exec checks for workflow scripts
        if str(path) == "/usr/bin/pigz":
            return False
        return real_access(path, mode)

    monkeypatch.setattr(config_mod.os, "access", fake_access)

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
    )
    assert cfg["compression_cmd"] == "/bin/gunzip"


def test_set_config_values_arch_other_string(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    monkeypatch.setattr(config_mod.platform, "machine", lambda: "mips64")

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-4.6", "genome": "b37"}
    )
    assert cfg["arch"] == "mips64"


# --------------------------
# Missing file / not executable guards
# --------------------------

def test_set_config_values_missing_workflow_files_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    software_stacks:\n"
        "      gatk-3.5:\n"
        "        helpers:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "          vcf2hash: \"vcf2hash.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single:\n"
        "              default: \"v1\"\n"
        "              versions:\n"
        "                v1:\n"
        "                  script: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    # Touch only ONE file, leave others missing
    bash_dir = root / "workflows" / "bash" / "gatk-3.5"
    bash_dir.mkdir(parents=True, exist_ok=True)
    (bash_dir / "env.sh").write_text("#!/bin/sh\n", encoding="utf-8")

    with pytest.raises(WorkflowResolutionError, match="referenced in the registry do not exist"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
        )


def test_set_config_values_not_executable_bash_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    software_stacks:\n"
        "      gatk-3.5:\n"
        "        helpers:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "          vcf2hash: \"vcf2hash.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single:\n"
        "              default: \"v1\"\n"
        "              versions:\n"
        "                v1:\n"
        "                  script: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    # Create all files but NOT executable (+x missing)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh", make_exec=False)

    with pytest.raises(WorkflowResolutionError, match="Missing \\+x"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_backend": "bash", "software_stack": "gatk-3.5"}
        )
