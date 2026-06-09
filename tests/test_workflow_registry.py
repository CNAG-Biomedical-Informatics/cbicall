import json
import os

import pytest

from cbicall import workflow_registry as wr
from cbicall.errors import ParameterValidationError, WorkflowResolutionError
from cbicall.models import WorkflowSpec
from helpers import make_executable, write_workflow_schema


def _versioned(script="wes_single.sh"):
    return {"default_registry_version": "v1", "registry_versions": {"v1": {"script": script}}}


def test_workflow_registry_resolve_errors_and_version_syntax(tmp_path):
    registry = {"workflows": {"bash": {"base_dir": "workflows/bash", "software_stacks": {"gatk-4.6": {"helpers": {}, "pipelines": {"wes": {"single": _versioned()}}}}}}}
    cfg = {"workflow_backend": "bash", "software_stack": "gatk-4.6", "pipeline": "wes", "mode": "single"}

    with pytest.raises(WorkflowResolutionError, match="Backend not defined"):
        wr.resolve_workflow_spec({**cfg, "workflow_backend": "missing"}, registry, tmp_path)
    with pytest.raises(WorkflowResolutionError, match="Software stack not defined"):
        wr.resolve_workflow_spec({**cfg, "software_stack": "missing"}, registry, tmp_path)
    with pytest.raises(WorkflowResolutionError, match="Pipeline not defined"):
        wr.resolve_workflow_spec({**cfg, "pipeline": "wgs"}, registry, tmp_path)
    with pytest.raises(WorkflowResolutionError, match="Mode not defined"):
        wr.resolve_workflow_spec({**cfg, "mode": "cohort"}, registry, tmp_path)
    with pytest.raises(WorkflowResolutionError, match="missing helper keys"):
        wr.resolve_workflow_spec(cfg, registry, tmp_path)

    assert wr._resolve_pipeline_implementation("legacy.sh", None, backend="bash", software_stack="gatk", pipeline="wes", mode="single") == ("legacy", "legacy.sh")
    with pytest.raises(WorkflowResolutionError, match="legacy registry syntax"):
        wr._resolve_pipeline_implementation("legacy.sh", "v1", backend="bash", software_stack="gatk", pipeline="wes", mode="single")
    with pytest.raises(WorkflowResolutionError, match="expected a versioned object"):
        wr._resolve_pipeline_implementation([], None, backend="bash", software_stack="gatk", pipeline="wes", mode="single")
    with pytest.raises(WorkflowResolutionError, match="must define registry versions"):
        wr._resolve_pipeline_implementation({}, None, backend="bash", software_stack="gatk", pipeline="wes", mode="single")
    with pytest.raises(WorkflowResolutionError, match="default registry version"):
        wr._resolve_pipeline_implementation({"registry_versions": {"v1": "x.sh"}}, None, backend="bash", software_stack="gatk", pipeline="wes", mode="single")
    with pytest.raises(WorkflowResolutionError, match="Available: v1"):
        wr._resolve_pipeline_implementation({"default_registry_version": "v1", "registry_versions": {"v1": "x.sh"}}, "v2", backend="bash", software_stack="gatk", pipeline="wes", mode="single")
    with pytest.raises(WorkflowResolutionError, match="missing keys"):
        wr._resolve_pipeline_implementation({"default_registry_version": "v1", "registry_versions": {"v1": {"provider": "nf-core", "source": "nf-core/demo"}}}, None, backend="nextflow", software_stack="nf-core", pipeline="demo", mode="single")
    with pytest.raises(WorkflowResolutionError, match="must define a script path"):
        wr._resolve_pipeline_implementation({"default_registry_version": "v1", "registry_versions": {"v1": {"note": "bad"}}}, None, backend="bash", software_stack="gatk", pipeline="wes", mode="single")


def test_workflow_registry_resolves_backends_and_file_validation(tmp_path):
    root = tmp_path
    bash_dir = root / "workflows" / "bash" / "gatk-4.6"
    bash_dir.mkdir(parents=True)
    for name in ["wes_single.sh", "env.sh", "coverage.sh", "jaccard.sh", "vcf2sex.sh", "vcf2hash.sh"]:
        p = bash_dir / name
        p.write_text("#!/bin/sh\n", encoding="utf-8")
        make_executable(p)
    registry = {
        "workflows": {
            "bash": {
                "base_dir": "workflows/bash",
                "software_stacks": {
                    "gatk-4.6": {
                        "helpers": {"env": "env.sh", "coverage": "coverage.sh", "jaccard": "jaccard.sh", "vcf2sex": "vcf2sex.sh", "vcf2hash": "vcf2hash.sh"},
                        "profiles": {"cnag-hpc": {"env": "cnag-hpc-env.sh"}},
                        "pipelines": {"wes": {"single": _versioned()}},
                    }
                },
            }
        }
    }
    spec = wr.resolve_workflow_spec({"workflow_backend": "bash", "software_stack": "gatk-4.6", "pipeline": "wes", "mode": "single"}, registry, root)
    assert spec.helpers["env"].endswith("env.sh")
    assert spec.profiles["cnag-hpc"]["env"].endswith("cnag-hpc-env.sh")
    wr.validate_resolved_workflow_files(spec)

    non_exec = bash_dir / "not_exec.sh"
    non_exec.write_text("#!/bin/sh\n", encoding="utf-8")
    with pytest.raises(WorkflowResolutionError, match=r"Missing \+x"):
        wr.validate_resolved_workflow_files(WorkflowSpec(backend="bash", pipeline="wes", mode="single", software_stack="gatk-4.6", registry_version="v1", entrypoint=str(non_exec), helpers={}))

    missing = WorkflowSpec(backend="nextflow", pipeline="wes", mode="single", software_stack="gatk-4.6", registry_version="v1", entrypoint=str(root / "missing.nf"), config_file=str(root / "missing.yaml"), helpers={"coverage": str(root / "missing.sh")})
    with pytest.raises(WorkflowResolutionError, match="do not exist"):
        wr.validate_resolved_workflow_files(missing)

    helper = root / "helper.sh"
    helper.write_text("#!/bin/sh\n", encoding="utf-8")
    config = root / "config.yaml"
    config.write_text("x: 1\n", encoding="utf-8")
    entry = root / "workflow.nf"
    entry.write_text("workflow {}\n", encoding="utf-8")
    with pytest.raises(WorkflowResolutionError, match=r"Missing \+x"):
        wr.validate_resolved_workflow_files(WorkflowSpec(backend="nextflow", pipeline="wes", mode="single", software_stack="gatk-4.6", registry_version="v1", entrypoint=str(entry), config_file=str(config), helpers={"coverage": str(helper)}))

    wr.validate_resolved_workflow_files(WorkflowSpec(backend="nextflow", pipeline="demo", mode="single", software_stack="nf-core", registry_version="v1", entrypoint="nf-core/demo", metadata={"provider": "nf-core"}))


def test_workflow_registry_load_and_schema_errors(tmp_path):
    registry = tmp_path / "registry.yaml"
    schema = tmp_path / "schema.json"
    registry.write_text("workflows: {}\n", encoding="utf-8")
    schema.write_text(json.dumps({"type": "object", "required": ["workflows"], "properties": {"workflows": {"type": "object", "minProperties": 1}}}), encoding="utf-8")
    with pytest.raises(ParameterValidationError, match="failed schema validation"):
        wr.load_workflow_registry(registry, schema)

    project = tmp_path / "project"
    (project / "workflows" / "registry").mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="Workflow registry not found"):
        wr.resolve_registry_context(project)
    (project / "workflows" / "registry" / "cbicall-workflow-registry.yaml").write_text("workflows: {}\n", encoding="utf-8")
    with pytest.raises(FileNotFoundError, match="Workflow schema not found"):
        wr.resolve_registry_context(project)
