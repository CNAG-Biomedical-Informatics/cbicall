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
    engine = cfg_in["workflow_engine"]
    version = cfg_in["gatk_version"]
    pipeline = cfg_in["pipeline"]
    mode = cfg_in["mode"]
    requested_pipeline_version = cfg_in.get("pipeline_version")

    workflows = registry["workflows"]
    if engine not in workflows:
        raise WorkflowResolutionError(f"Engine not defined in workflow registry: {engine}")

    eng_cfg = workflows[engine]
    versions_cfg = eng_cfg["versions"]
    if version not in versions_cfg:
        raise WorkflowResolutionError(f"Version not defined for engine '{engine}': {version}")

    ver_cfg = versions_cfg[version]
    helpers = ver_cfg.get("helpers", {})
    pipelines_cfg = ver_cfg["pipelines"]

    if pipeline not in pipelines_cfg:
        raise WorkflowResolutionError(f"Pipeline not defined for {engine}/{version}: {pipeline}")
    if mode not in pipelines_cfg[pipeline]:
        raise WorkflowResolutionError(
            f"Mode not defined for pipeline '{pipeline}' in {engine}/{version}: {mode}"
        )

    base_dir = (project_root / eng_cfg["base_dir"] / version).resolve()
    pipeline_version, implementation = _resolve_pipeline_implementation(
        pipelines_cfg[pipeline][mode],
        requested_pipeline_version,
        engine=engine,
        gatk_version=version,
        pipeline=pipeline,
        mode=mode,
    )
    if isinstance(implementation, dict):
        script_name = implementation.get("script")
    else:
        script_name = implementation
    profiles = _resolve_profiles(ver_cfg.get("profiles", {}), base_dir)

    if engine == "bash":
        needed_helpers = ["env", "coverage", "jaccard", "vcf2sex", "vcf2hash"]
        missing_helpers = [k for k in needed_helpers if k not in helpers]
        if missing_helpers:
            raise WorkflowResolutionError(
                f"Workflow registry is missing helper keys for bash/{version}: {missing_helpers}"
            )
        return WorkflowSpec(
            engine=engine,
            pipeline=pipeline,
            mode=mode,
            gatk_version=version,
            pipeline_version=pipeline_version,
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

    if engine == "snakemake":
        if "config" not in helpers:
            raise WorkflowResolutionError(
                f"Workflow registry is missing helper key 'config' for snakemake/{version}"
            )
        return WorkflowSpec(
            engine=engine,
            pipeline=pipeline,
            mode=mode,
            gatk_version=version,
            pipeline_version=pipeline_version,
            entrypoint=str(base_dir / script_name),
            config_file=str(base_dir / helpers["config"]),
            profiles=profiles,
        )

    if engine == "nextflow":
        if isinstance(implementation, dict) and implementation.get("source_type") == "nf-core":
            return WorkflowSpec(
                engine=engine,
                pipeline=pipeline,
                mode=mode,
                gatk_version=version,
                pipeline_version=pipeline_version,
                entrypoint=str(implementation["source"]),
                config_file=None,
                helpers={},
                profiles=profiles,
                metadata={
                    "source_type": str(implementation["source_type"]),
                    "source": str(implementation["source"]),
                    "release": str(implementation["release"]),
                    "default_outdir": str(implementation.get("default_outdir", pipeline)),
                },
            )

        needed_helpers = ["config", "coverage", "vcf2sex", "vcf2hash"]
        missing_helpers = [k for k in needed_helpers if k not in helpers]
        if missing_helpers:
            raise WorkflowResolutionError(
                f"Workflow registry is missing helper keys for nextflow/{version}: {missing_helpers}"
            )
        return WorkflowSpec(
            engine=engine,
            pipeline=pipeline,
            mode=mode,
            gatk_version=version,
            pipeline_version=pipeline_version,
            entrypoint=str(base_dir / script_name),
            config_file=str(base_dir / helpers["config"]),
            helpers={
                "coverage": str(base_dir / helpers["coverage"]),
                "vcf2sex": str(base_dir / helpers["vcf2sex"]),
                "vcf2hash": str(base_dir / helpers["vcf2hash"]),
            },
            profiles=profiles,
        )

    raise WorkflowResolutionError(f"Unsupported workflow_engine: {engine}")


def validate_resolved_workflow_files(workflow: WorkflowSpec) -> None:
    if workflow.metadata.get("source_type") == "nf-core":
        return

    must_exist = [("workflow.entrypoint", workflow.entrypoint)]
    if workflow.engine == "bash":
        must_exist.extend((f"workflow.helpers.{name}", path) for name, path in workflow.helpers.items())
    elif workflow.engine in {"snakemake", "nextflow"}:
        must_exist.append(("workflow.config_file", workflow.config_file))
        if workflow.engine == "nextflow":
            must_exist.extend((f"workflow.helpers.{name}", path) for name, path in workflow.helpers.items())

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
    elif workflow.engine == "nextflow":
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
    requested_pipeline_version,
    *,
    engine: str,
    gatk_version: str,
    pipeline: str,
    mode: str,
) -> Tuple[str, Any]:
    label = f"{engine}/{gatk_version}/{pipeline}/{mode}"
    if isinstance(mode_cfg, str):
        if requested_pipeline_version:
            raise WorkflowResolutionError(
                f"pipeline_version was set to '{requested_pipeline_version}', but {label} "
                "uses legacy registry syntax without implementation versions."
            )
        return "legacy", mode_cfg

    if not isinstance(mode_cfg, dict):
        raise WorkflowResolutionError(f"Invalid registry entry for {label}: expected a versioned object.")

    implementations = mode_cfg.get("versions")
    if not isinstance(implementations, dict) or not implementations:
        raise WorkflowResolutionError(f"Registry entry for {label} must define implementation versions.")

    default_version = mode_cfg.get("default")
    selected_version = requested_pipeline_version or default_version
    if not selected_version:
        raise WorkflowResolutionError(f"Registry entry for {label} must define a default implementation version.")
    selected_version = str(selected_version)

    if selected_version not in implementations:
        available = ", ".join(sorted(str(k) for k in implementations))
        raise WorkflowResolutionError(
            f"pipeline_version='{selected_version}' is not defined for {label}. Available: {available}"
        )

    implementation = implementations[selected_version]
    if isinstance(implementation, str):
        return selected_version, implementation
    if isinstance(implementation, dict):
        if implementation.get("script"):
            return selected_version, str(implementation["script"])
        if implementation.get("source_type") == "nf-core":
            missing = [key for key in ("source", "release") if not implementation.get(key)]
            if missing:
                raise WorkflowResolutionError(
                    f"Implementation {selected_version!r} for {label} is missing keys: {missing}"
                )
            return selected_version, implementation
    raise WorkflowResolutionError(
        f"Implementation {selected_version!r} for {label} must define a script path or external source."
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
