import json
import os
from pathlib import Path
from typing import Tuple

import yaml
from jsonschema import Draft202012Validator

from .errors import ParameterValidationError, WorkflowResolutionError
from .models import WorkflowSpec


def load_workflow_registry(registry_yaml: Path, schema_json: Path) -> dict:
    with registry_yaml.open("r") as fh:
        data = yaml.safe_load(fh) or {}
    schema = _load_json(schema_json)
    _validate_with_schema(data, schema, label=str(registry_yaml))
    return data


def get_project_root(module_file: str) -> Path:
    here = Path(module_file).resolve()
    return here.parents[2]


def resolve_registry_context(project_root: Path) -> dict:
    registry_yaml, schema_json = _registry_paths(project_root)
    if not registry_yaml.is_file():
        raise FileNotFoundError(f"Workflow registry not found: {registry_yaml}")
    if not schema_json.is_file():
        raise FileNotFoundError(f"Workflow schema not found: {schema_json}")
    return load_workflow_registry(registry_yaml, schema_json)


def resolve_workflow_spec(cfg_in: dict, registry: dict, project_root: Path) -> WorkflowSpec:
    engine = cfg_in["workflow_engine"]
    version = cfg_in["gatk_version"]
    pipeline = cfg_in["pipeline"]
    mode = cfg_in["mode"]

    workflows = registry["workflows"]
    if engine not in workflows:
        raise WorkflowResolutionError(f"Engine not defined in workflow registry: {engine}")

    eng_cfg = workflows[engine]
    versions_cfg = eng_cfg["versions"]
    if version not in versions_cfg:
        raise WorkflowResolutionError(f"Version not defined for engine '{engine}': {version}")

    ver_cfg = versions_cfg[version]
    common = ver_cfg.get("common", {})
    pipelines_cfg = ver_cfg["pipelines"]

    if pipeline not in pipelines_cfg:
        raise WorkflowResolutionError(f"Pipeline not defined for {engine}/{version}: {pipeline}")
    if mode not in pipelines_cfg[pipeline]:
        raise WorkflowResolutionError(
            f"Mode not defined for pipeline '{pipeline}' in {engine}/{version}: {mode}"
        )

    base_dir = (project_root / eng_cfg["base_dir"] / version).resolve()
    script_name = pipelines_cfg[pipeline][mode]

    if engine == "bash":
        needed_common = ["env", "coverage", "jaccard", "vcf2sex"]
        missing_common = [k for k in needed_common if k not in common]
        if missing_common:
            raise WorkflowResolutionError(
                f"Workflow registry is missing common keys for bash/{version}: {missing_common}"
            )
        return WorkflowSpec(
            engine=engine,
            pipeline=pipeline,
            mode=mode,
            gatk_version=version,
            entrypoint=str(base_dir / script_name),
            helpers={
                "env": str(base_dir / common["env"]),
                "coverage": str(base_dir / common["coverage"]),
                "jaccard": str(base_dir / common["jaccard"]),
                "vcf2sex": str(base_dir / common["vcf2sex"]),
            },
        )

    if engine == "snakemake":
        if "config" not in common:
            raise WorkflowResolutionError(
                f"Workflow registry is missing common key 'config' for snakemake/{version}"
            )
        return WorkflowSpec(
            engine=engine,
            pipeline=pipeline,
            mode=mode,
            gatk_version=version,
            entrypoint=str(base_dir / script_name),
            config_file=str(base_dir / common["config"]),
        )

    if engine == "nextflow":
        raise WorkflowResolutionError("workflow_engine='nextflow' is declared but not implemented.")

    raise WorkflowResolutionError(f"Unsupported workflow_engine: {engine}")


def validate_resolved_workflow_files(workflow: WorkflowSpec) -> None:
    must_exist = [("workflow.entrypoint", workflow.entrypoint)]
    if workflow.engine == "bash":
        must_exist.extend((f"workflow.helpers.{name}", path) for name, path in workflow.helpers.items())
    elif workflow.engine == "snakemake":
        must_exist.append(("workflow.config_file", workflow.config_file))

    missing_files = [(label, path) for label, path in must_exist if path and not Path(path).exists()]
    if missing_files:
        raise WorkflowResolutionError(
            "One or more workflow files referenced in the registry do not exist: "
            + ", ".join(f"{label} -> {path}" for label, path in missing_files)
        )

    if workflow.engine == "bash":
        exe_paths = [("workflow.entrypoint", workflow.entrypoint)]
        exe_paths.extend((f"workflow.helpers.{name}", path) for name, path in workflow.helpers.items())
        not_exe = [(label, path) for label, path in exe_paths if path and not os.access(path, os.X_OK)]
        if not_exe:
            raise WorkflowResolutionError(
                "Missing +x on one or more workflow scripts: "
                + ", ".join(f"{label} -> {path}" for label, path in not_exe)
            )


def _load_json(path: Path) -> dict:
    with path.open("r") as fh:
        return json.load(fh)


def _validate_with_schema(data: dict, schema: dict, label: str) -> None:
    v = Draft202012Validator(schema)
    errors = sorted(v.iter_errors(data), key=lambda e: list(e.path))
    if errors:
        lines = [f"{label} failed schema validation:"]
        for e in errors[:25]:
            loc = ".".join(str(x) for x in e.path) or "(root)"
            lines.append(f"- {loc}: {e.message}")
        raise ParameterValidationError("\n".join(lines))


def _registry_paths(project_root: Path) -> Tuple[Path, Path]:
    registry_yaml = (project_root / "workflows" / "config" / "cbicall.workflows.yaml").resolve()
    schema_json = (project_root / "workflows" / "schema" / "workflows.schema.json").resolve()
    return registry_yaml, schema_json
