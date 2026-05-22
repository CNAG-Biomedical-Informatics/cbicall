import json
import os
from pathlib import Path
from typing import Any, Tuple

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
    backend = cfg_in["workflow_backend"]
    software_stack = cfg_in["software_stack"]
    pipeline = cfg_in["pipeline"]
    mode = cfg_in["mode"]
    requested_registry_version = cfg_in.get("registry_version")

    workflows = registry["workflows"]
    if backend not in workflows:
        raise WorkflowResolutionError(f"Backend not defined in workflow registry: {backend}")

    backend_cfg = workflows[backend]
    software_stacks_cfg = backend_cfg["software_stacks"]
    if software_stack not in software_stacks_cfg:
        raise WorkflowResolutionError(f"Software stack not defined for backend '{backend}': {software_stack}")

    ver_cfg = software_stacks_cfg[software_stack]
    helpers = ver_cfg.get("helpers", {})
    pipelines_cfg = ver_cfg["pipelines"]

    if pipeline not in pipelines_cfg:
        raise WorkflowResolutionError(f"Pipeline not defined for {backend}/{software_stack}: {pipeline}")
    if mode not in pipelines_cfg[pipeline]:
        raise WorkflowResolutionError(
            f"Mode not defined for pipeline '{pipeline}' in {backend}/{software_stack}: {mode}"
        )

    base_dir = (project_root / backend_cfg["base_dir"] / software_stack).resolve()
    registry_version, implementation = _resolve_pipeline_implementation(
        pipelines_cfg[pipeline][mode],
        requested_registry_version,
        backend=backend,
        software_stack=software_stack,
        pipeline=pipeline,
        mode=mode,
    )
    if isinstance(implementation, dict):
        script_name = implementation.get("script")
    else:
        script_name = implementation
    profiles = _resolve_profiles(ver_cfg.get("profiles", {}), base_dir)

    if backend == "bash":
        needed_helpers = ["env", "coverage", "jaccard", "vcf2sex", "vcf2hash"]
        missing_helpers = [k for k in needed_helpers if k not in helpers]
        if missing_helpers:
            raise WorkflowResolutionError(
                f"Workflow registry is missing helper keys for bash/{software_stack}: {missing_helpers}"
            )
        return WorkflowSpec(
            backend=backend,
            pipeline=pipeline,
            mode=mode,
            software_stack=software_stack,
            registry_version=registry_version,
            entrypoint=str(base_dir / script_name),
            helpers={
                "env": str(base_dir / helpers["env"]),
                "coverage": str(base_dir / helpers["coverage"]),
                "jaccard": str(base_dir / helpers["jaccard"]),
                "vcf2sex": str(base_dir / helpers["vcf2sex"]),
                "vcf2hash": str(base_dir / helpers["vcf2hash"]),
            },
            profiles=profiles,
        )

    if backend == "snakemake":
        if "config" not in helpers:
            raise WorkflowResolutionError(
                f"Workflow registry is missing helper key 'config' for snakemake/{software_stack}"
            )
        return WorkflowSpec(
            backend=backend,
            pipeline=pipeline,
            mode=mode,
            software_stack=software_stack,
            registry_version=registry_version,
            entrypoint=str(base_dir / script_name),
            config_file=str(base_dir / helpers["config"]),
            profiles=profiles,
        )

    if backend == "nextflow":
        if isinstance(implementation, dict) and implementation.get("provider") == "nf-core":
            return WorkflowSpec(
                backend=backend,
                pipeline=pipeline,
                mode=mode,
                software_stack=software_stack,
                registry_version=registry_version,
                entrypoint=str(implementation["source"]),
                config_file=None,
                helpers={},
                profiles=profiles,
                metadata={
                    "provider": str(implementation["provider"]),
                    "source": str(implementation["source"]),
                    "release": str(implementation["release"]),
                    "default_outdir": str(implementation.get("default_outdir", pipeline)),
                    "canonical_outputs": [
                        {
                            "name": str(item["name"]),
                            "type": str(item["type"]),
                            "pattern": str(item["pattern"]),
                        }
                        for item in implementation.get("canonical_outputs", [])
                    ],
                },
            )

        needed_helpers = ["config", "coverage", "vcf2sex", "vcf2hash"]
        missing_helpers = [k for k in needed_helpers if k not in helpers]
        if missing_helpers:
            raise WorkflowResolutionError(
                f"Workflow registry is missing helper keys for nextflow/{software_stack}: {missing_helpers}"
            )
        return WorkflowSpec(
            backend=backend,
            pipeline=pipeline,
            mode=mode,
            software_stack=software_stack,
            registry_version=registry_version,
            entrypoint=str(base_dir / script_name),
            config_file=str(base_dir / helpers["config"]),
            helpers={
                "coverage": str(base_dir / helpers["coverage"]),
                "vcf2sex": str(base_dir / helpers["vcf2sex"]),
                "vcf2hash": str(base_dir / helpers["vcf2hash"]),
            },
            profiles=profiles,
        )

    raise WorkflowResolutionError(f"Unsupported workflow_backend: {backend}")


def validate_resolved_workflow_files(workflow: WorkflowSpec) -> None:
    if workflow.metadata.get("provider") == "nf-core":
        return

    must_exist = [("workflow.entrypoint", workflow.entrypoint)]
    if workflow.backend == "bash":
        must_exist.extend((f"workflow.helpers.{name}", path) for name, path in workflow.helpers.items())
    elif workflow.backend in {"snakemake", "nextflow"}:
        must_exist.append(("workflow.config_file", workflow.config_file))
        if workflow.backend == "nextflow":
            must_exist.extend((f"workflow.helpers.{name}", path) for name, path in workflow.helpers.items())

    missing_files = [(label, path) for label, path in must_exist if path and not Path(path).exists()]
    if missing_files:
        raise WorkflowResolutionError(
            "One or more workflow files referenced in the registry do not exist: "
            + ", ".join(f"{label} -> {path}" for label, path in missing_files)
        )

    if workflow.backend == "bash":
        exe_paths = [("workflow.entrypoint", workflow.entrypoint)]
        exe_paths.extend((f"workflow.helpers.{name}", path) for name, path in workflow.helpers.items())
        not_exe = [(label, path) for label, path in exe_paths if path and not os.access(path, os.X_OK)]
        if not_exe:
            raise WorkflowResolutionError(
                "Missing +x on one or more workflow scripts: "
                + ", ".join(f"{label} -> {path}" for label, path in not_exe)
            )
    elif workflow.backend == "nextflow":
        exe_paths = [(f"workflow.helpers.{name}", path) for name, path in workflow.helpers.items()]
        not_exe = [(label, path) for label, path in exe_paths if path and not os.access(path, os.X_OK)]
        if not_exe:
            raise WorkflowResolutionError(
                "Missing +x on one or more Nextflow helper scripts: "
                + ", ".join(f"{label} -> {path}" for label, path in not_exe)
            )


def _load_json(path: Path) -> dict:
    with path.open("r") as fh:
        return json.load(fh)


def _resolve_profiles(profiles_cfg: dict, base_dir: Path) -> dict:
    profiles = {}
    for profile_name, helper_overrides in profiles_cfg.items():
        profiles[str(profile_name)] = {
            str(helper_name): str(base_dir / helper_path)
            for helper_name, helper_path in helper_overrides.items()
        }
    return profiles


def _resolve_pipeline_implementation(
    mode_cfg,
    requested_registry_version,
    *,
    backend: str,
    software_stack: str,
    pipeline: str,
    mode: str,
) -> Tuple[str, Any]:
    label = f"{backend}/{software_stack}/{pipeline}/{mode}"
    if isinstance(mode_cfg, str):
        if requested_registry_version:
            raise WorkflowResolutionError(
                f"registry_version was set to '{requested_registry_version}', but {label} "
                "uses legacy registry syntax without registry versions."
            )
        return "legacy", mode_cfg

    if not isinstance(mode_cfg, dict):
        raise WorkflowResolutionError(f"Invalid registry entry for {label}: expected a versioned object.")

    registry_versions = mode_cfg.get("registry_versions")
    if not isinstance(registry_versions, dict) or not registry_versions:
        raise WorkflowResolutionError(f"Registry entry for {label} must define registry versions.")

    default_registry_version = mode_cfg.get("default_registry_version")
    selected_registry_version = requested_registry_version or default_registry_version
    if not selected_registry_version:
        raise WorkflowResolutionError(f"Registry entry for {label} must define a default registry version.")
    selected_registry_version = str(selected_registry_version)

    if selected_registry_version not in registry_versions:
        available = ", ".join(sorted(str(k) for k in registry_versions))
        raise WorkflowResolutionError(
            f"registry_version='{selected_registry_version}' is not defined for {label}. Available: {available}"
        )

    implementation = registry_versions[selected_registry_version]
    if isinstance(implementation, str):
        return selected_registry_version, implementation
    if isinstance(implementation, dict):
        if implementation.get("script"):
            return selected_registry_version, str(implementation["script"])
        if implementation.get("provider") == "nf-core":
            missing = [key for key in ("source", "release") if not implementation.get(key)]
            if missing:
                raise WorkflowResolutionError(
                    f"Implementation {selected_registry_version!r} for {label} is missing keys: {missing}"
                )
            return selected_registry_version, implementation
    raise WorkflowResolutionError(
        f"Implementation {selected_registry_version!r} for {label} must define a script path or external source."
    )


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
    registry_yaml = (project_root / "workflows" / "registry" / "cbicall-workflow-registry.yaml").resolve()
    schema_json = (project_root / "workflows" / "schema" / "cbicall-workflow-registry.schema.json").resolve()
    return registry_yaml, schema_json
