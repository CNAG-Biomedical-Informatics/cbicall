import copy
import os
import time
import shutil
import socket
import platform
import getpass
from pathlib import Path

import yaml

from .errors import ParameterValidationError, WorkflowResolutionError
from .models import InputsSpec, ResolvedConfig, WorkflowSpec
from .resources import build_bundle_resource_metadata, validate_installed_bundle_resource
from .workflow_registry import (
    _validate_with_schema as _registry_validate_with_schema,
    get_project_root,
    load_workflow_registry,
    resolve_registry_context,
    resolve_workflow_spec,
    validate_resolved_workflow_files,
)

###########################################################################
# Note:
# Pipeline policy and semantic validation are intentionally kept in code
# (enums, defaults, genome/tool compatibility rules) to provide a stable,
# fail-fast contract with good error messages.
#
# Workflow *wiring* (script paths, tool versions, pipeline/mode mappings)
# is loaded from an external YAML registry and validated with JSON Schema,
# so new pipelines can be added without modifying code.
#
# Code guarantees: loaded config is syntactically valid, semantically sane,
# and only dispatches to declared + executable workflow scripts.
###########################################################################


# Allowed values
MODE_VALUES = {"single", "cohort"}
PIPELINE_VALUES = {"wes", "wgs", "mit"}
ORGANISM_VALUES = {"Homo sapiens", "Mus musculus"}
TECHNOLOGY_VALUES = {"Illumina HiSeq", "NovaSeq"}
WORKFLOW_ENGINE_VALUES = {"bash", "nextflow", "snakemake"}
GATK_VALUES = {"gatk-3.5", "gatk-4.6"}
WORKFLOW_VERSION_VALUES = GATK_VALUES | {"nf-core"}
GENOME_VALUES = {"b37", "hg38", "rsrs", "external"}


# Defaults
_DEFAULTS = {
    "mode": "single",
    "input_dir": None,
    "sample_map": None,
    "output_basename": None,
    "pipeline": "wes",
    "organism": "Homo sapiens",
    "technology": "Illumina HiSeq",
    "workflow_engine": "bash",
    "profile": "local",
    "snakemake_parameters": {},
    "nextflow_parameters": {},
    "nfcore_profile": None,
    "nfcore_parameters": {},
    "nfcore_singularity_cache_dir": None,
    "gatk_version": "gatk-3.5",
    "workflow_version": None,
    "pipeline_version": None,
    "project_dir": "cbicall",
    "cleanup_bam": False,
    "genome": None,  # Option 1: effective default assigned in _apply_genome_rules
    "resource": "cbicall-germline-resources-v1",
}

_RUNTIME_ONLY_KEYS = {
    "profile": "profile is a runtime option; use --runtime-profile instead of setting it in the parameters YAML.",
}


# Allowed pipeline-mode combinations per GATK version
_ALLOWED_COMBOS = {
    "gatk-3.5": {
        "wes": ["single", "cohort"],
        "mit": ["single", "cohort"],
    },
    "gatk-4.6": {
        "wes": ["single", "cohort"],
        "wgs": ["single", "cohort"],
    },
}


def _validate_enum(name, value, allowed):
    if value is None:
        return
    if value not in allowed:
        raise ParameterValidationError(
            f"Invalid value for '{name}': {value!r}. Allowed: {sorted(allowed)}"
        )


def _validate_combos(cfg: dict) -> None:
    """Validate pipeline/mode combo for selected GATK version."""
    version = cfg["gatk_version"]
    if version == "nf-core":
        return
    pipeline = cfg["pipeline"]
    mode = cfg["mode"]

    allowed_for_version = _ALLOWED_COMBOS.get(version, {})
    modes_for_pipeline = allowed_for_version.get(pipeline, [])

    if mode not in modes_for_pipeline:
        raise ParameterValidationError(
            f"Pipeline-mode '{pipeline}_{mode}' is not supported for GATK version {version}"
        )


def _apply_genome_rules(cfg: dict, user_provided_genome: bool) -> None:
    """
    Apply genome rules shared by YAML and CLI config building.
    Mutates cfg.
    """
    # Block unsupported union: snakemake + gatk-3.5
    if (
        cfg.get("workflow_engine") == "snakemake"
        and cfg.get("gatk_version") == "gatk-3.5"
    ):
        raise ParameterValidationError(
            "workflow_engine='snakemake' is not supported for gatk_version='gatk-3.5'. "
            "Use workflow_engine='bash' or select gatk_version='gatk-4.6'."
        )

    # Block unsupported union: mit + snakemake
    if cfg.get("pipeline") == "mit" and cfg.get("workflow_engine") == "snakemake":
        raise ParameterValidationError(
            "The combination pipeline='mit' with workflow_engine='snakemake' is not supported."
        )

    if cfg.get("workflow_engine") == "nextflow":
        if cfg.get("gatk_version") not in {"gatk-4.6", "nf-core"}:
            raise ParameterValidationError(
                "workflow_engine='nextflow' is supported only for workflow versions 'gatk-4.6' or 'nf-core'."
            )
        if cfg.get("pipeline") == "mit":
            raise ParameterValidationError(
                "The combination pipeline='mit' with workflow_engine='nextflow' is not supported."
            )

    if cfg.get("gatk_version") == "nf-core":
        if cfg.get("workflow_engine") != "nextflow":
            raise ParameterValidationError("workflow_version='nf-core' requires workflow_engine='nextflow'.")
        if user_provided_genome and cfg.get("genome") != "external":
            raise ParameterValidationError(
                "For workflow_version='nf-core', genome is managed by the external Nextflow workflow. "
                "Remove 'genome' from the YAML or set genome='external'."
            )
        cfg["genome"] = "external"
        return

    # If genome omitted, assign effective defaults
    # - MIT defaults to rsrs
    # - everything else defaults to b37
    if cfg.get("genome") is None:
        cfg["genome"] = "rsrs" if cfg.get("pipeline") == "mit" else "b37"

    # MIT pipeline forces genome to rsrs
    if cfg.get("pipeline") == "mit":
        if user_provided_genome and cfg.get("genome") != "rsrs":
            raise ParameterValidationError(
                "For pipeline='mit', genome is fixed to 'rsrs'. "
                "Remove 'genome' from the YAML or set genome='rsrs'."
            )
        cfg["genome"] = "rsrs"

    # hg38 is supported only for WGS
    if cfg["genome"] == "hg38" and cfg["pipeline"] != "wgs":
        raise ParameterValidationError("genome='hg38' is only supported for pipeline='wgs'.")


def _validate_enums_except_genome(cfg: dict) -> None:
    _validate_enum("mode", cfg["mode"], MODE_VALUES)
    if cfg.get("workflow_version") == "nf-core":
        pipeline = cfg.get("pipeline")
        if pipeline is None or not str(pipeline).strip():
            raise ParameterValidationError("pipeline must be a non-empty value.")
        cfg["pipeline"] = str(pipeline).strip()
    else:
        _validate_enum("pipeline", cfg["pipeline"], PIPELINE_VALUES)
        _validate_enum("gatk_version", cfg["gatk_version"], GATK_VALUES)
    _validate_enum("organism", cfg["organism"], ORGANISM_VALUES)
    _validate_enum("technology", cfg["technology"], TECHNOLOGY_VALUES)
    _validate_enum("workflow_engine", cfg["workflow_engine"], WORKFLOW_ENGINE_VALUES)


def _validate_mapping_parameter(cfg: dict, key: str) -> dict:
    value = cfg.get(key)
    if value is None:
        value = {}
    if not isinstance(value, dict):
        raise ParameterValidationError(f"{key} must be a mapping.")
    cfg[key] = dict(value)
    return cfg[key]


def _validate_backend_parameter_settings(cfg: dict) -> None:
    snakemake_parameters = _validate_mapping_parameter(cfg, "snakemake_parameters")
    nextflow_parameters = _validate_mapping_parameter(cfg, "nextflow_parameters")

    if snakemake_parameters and cfg.get("workflow_engine") != "snakemake":
        raise ParameterValidationError("snakemake_parameters requires workflow_engine='snakemake'.")
    if nextflow_parameters and (cfg.get("workflow_engine") != "nextflow" or cfg.get("gatk_version") == "nf-core"):
        raise ParameterValidationError("nextflow_parameters requires a native workflow_engine='nextflow' run.")

    target = snakemake_parameters.get("target")
    if target is not None:
        if not isinstance(target, str) or not target.strip():
            raise ParameterValidationError("snakemake_parameters.target must be a non-empty string when provided.")
        snakemake_parameters["target"] = target.strip()

    reserved_snakemake = sorted(set(snakemake_parameters) & {"genome", "pipeline", "sample_map", "workspace"})
    if reserved_snakemake:
        raise ParameterValidationError(
            "snakemake_parameters cannot set CBIcall-controlled parameters: " + ", ".join(reserved_snakemake)
        )

    reserved_nextflow = sorted(
        set(nextflow_parameters)
        & {
            "pipeline",
            "genome",
            "threads",
            "cleanup_bam",
            "sample_map",
            "workspace",
            "coverage_script",
            "vcf2sex_script",
            "vcf2hash_script",
        }
    )
    if reserved_nextflow:
        raise ParameterValidationError(
            "nextflow_parameters cannot set CBIcall-controlled parameters: " + ", ".join(reserved_nextflow)
        )


def _validate_resource_settings(cfg: dict) -> None:
    resource = cfg.get("resource")
    if resource is None:
        raise ParameterValidationError("resource cannot be null.")
    resource = str(resource).strip()
    if not resource:
        raise ParameterValidationError("resource must be a non-empty value.")
    cfg["resource"] = resource


def _validate_profile_settings(cfg: dict) -> None:
    profile = cfg.get("profile")
    if profile is None:
        raise ParameterValidationError("profile cannot be null.")
    cfg["profile"] = str(profile)
    if not cfg["profile"].strip():
        raise ParameterValidationError("profile must be a non-empty value.")
    cfg["profile"] = cfg["profile"].strip()


def _validate_pipeline_version_settings(cfg: dict) -> None:
    pipeline_version = cfg.get("pipeline_version")
    if pipeline_version is None:
        return
    cfg["pipeline_version"] = str(pipeline_version)
    if not cfg["pipeline_version"].strip():
        raise ParameterValidationError("pipeline_version must be a non-empty value when provided.")
    cfg["pipeline_version"] = cfg["pipeline_version"].strip()


def _normalize_workflow_version_settings(cfg: dict) -> None:
    workflow_version = cfg.get("workflow_version")
    if workflow_version is None:
        return
    workflow_version = str(workflow_version).strip()
    if not workflow_version:
        raise ParameterValidationError("workflow_version must be a non-empty value when provided.")
    _validate_enum("workflow_version", workflow_version, WORKFLOW_VERSION_VALUES)
    if (
        cfg.get("gatk_version") != _DEFAULTS["gatk_version"]
        and cfg.get("gatk_version") != workflow_version
    ):
        raise ParameterValidationError(
            "workflow_version and gatk_version refer to the same workflow family version; "
            "set only one of them or use matching values."
        )
    cfg["workflow_version"] = workflow_version
    cfg["gatk_version"] = workflow_version


def _validate_nfcore_settings(cfg: dict) -> None:
    cache_dir = cfg.get("nfcore_singularity_cache_dir")
    if cfg.get("gatk_version") != "nf-core":
        if cache_dir is not None:
            raise ParameterValidationError(
                "nfcore_singularity_cache_dir requires workflow_version='nf-core'."
            )
        return
    if cfg.get("workflow_engine") != "nextflow":
        raise ParameterValidationError("workflow_version='nf-core' requires workflow_engine='nextflow'.")

    profile = cfg.get("nfcore_profile")
    if profile is None or not str(profile).strip():
        raise ParameterValidationError("workflow_version='nf-core' requires 'nfcore_profile'.")
    cfg["nfcore_profile"] = str(profile).strip()

    nfcore_parameters = cfg.get("nfcore_parameters")
    if nfcore_parameters is None:
        nfcore_parameters = {}
    if not isinstance(nfcore_parameters, dict):
        raise ParameterValidationError("nfcore_parameters must be a mapping of nf-core parameters.")
    reserved = sorted(set(nfcore_parameters) & {"outdir", "max_cpus"})
    if reserved:
        raise ParameterValidationError(
            "nfcore_parameters cannot set CBIcall-controlled parameters: " + ", ".join(reserved)
        )
    cfg["nfcore_parameters"] = dict(nfcore_parameters)

    if cache_dir is not None:
        cache_dir = str(cache_dir).strip()
        if not cache_dir:
            raise ParameterValidationError(
                "nfcore_singularity_cache_dir must be a non-empty path when provided."
            )
        cfg["nfcore_singularity_cache_dir"] = cache_dir


def _is_uri(value: str) -> bool:
    return "://" in value or value.startswith(("s3:", "gs:"))


def _resolve_nfcore_param_paths(value, base_dir: Path):
    if isinstance(value, dict):
        return {key: _resolve_nfcore_param_paths(item, base_dir) for key, item in value.items()}
    if isinstance(value, list):
        return [_resolve_nfcore_param_paths(item, base_dir) for item in value]
    if not isinstance(value, str) or not value.strip() or _is_uri(value):
        return value
    path = Path(value)
    if path.is_absolute():
        return str(path)
    candidate = base_dir / path
    if candidate.exists():
        return str(candidate.resolve())
    return value


def _validate_with_schema(data: dict, schema: dict, label: str) -> None:
    _registry_validate_with_schema(data, schema, label)


def read_param_file(yaml_file: str) -> dict:
    """
    Load YAML parameters, merge with defaults and validate.
    """
    yaml_path = Path(yaml_file)
    if not yaml_path.is_file():
        raise FileNotFoundError(f"Parameters file not found: {yaml_file}")

    with yaml_path.open("r") as fh:
        params = yaml.safe_load(fh) or {}

    # Start from defaults
    cfg = copy.deepcopy(_DEFAULTS)

    # Merge provided parameters and validate keys
    for key, value in params.items():
        if key in _RUNTIME_ONLY_KEYS:
            raise ParameterValidationError(_RUNTIME_ONLY_KEYS[key])
        if key not in cfg:
            raise ParameterValidationError(f"Parameter '{key}' does not exist (typo?)")
        cfg[key] = value

    # Validate enums (except genome; it may be None until rules apply)
    _normalize_workflow_version_settings(cfg)
    _validate_enums_except_genome(cfg)
    _validate_profile_settings(cfg)
    _validate_pipeline_version_settings(cfg)
    _validate_backend_parameter_settings(cfg)
    _validate_resource_settings(cfg)
    _validate_nfcore_settings(cfg)

    # Apply shared genome rules
    user_provided_genome = "genome" in params
    _apply_genome_rules(cfg, user_provided_genome=user_provided_genome)

    # Now validate genome
    _validate_enum("genome", cfg["genome"], GENOME_VALUES)

    # Make input_dir absolute (resolve relative to YAML location)
    if cfg.get("input_dir") is not None:
        input_dir = Path(cfg["input_dir"])
        if not input_dir.is_absolute():
            input_dir = yaml_path.parent / input_dir
        cfg["input_dir"] = str(input_dir.resolve())

    # Make sample_map absolute (resolve relative to YAML location)
    if cfg.get("sample_map") is not None:
        sm = Path(cfg["sample_map"])
        if not sm.is_absolute():
            sm = yaml_path.parent / sm
        cfg["sample_map"] = str(sm.resolve())

    cfg["nfcore_parameters"] = _resolve_nfcore_param_paths(cfg.get("nfcore_parameters", {}), yaml_path.parent)
    if cfg.get("nfcore_singularity_cache_dir") is not None:
        cache_dir = Path(cfg["nfcore_singularity_cache_dir"])
        if not cache_dir.is_absolute():
            cache_dir = yaml_path.parent / cache_dir
        cfg["nfcore_singularity_cache_dir"] = str(cache_dir.resolve())

    # Validate pipeline-mode combination for selected GATK version
    _validate_combos(cfg)

    return cfg


def _merge_and_validate_param_values(params: dict) -> dict:
    """Merge defaults with user parameters and run semantic validation."""
    cfg_in = copy.deepcopy(_DEFAULTS)
    cfg_in.update(params)
    # Validate enums (except genome), then apply genome rules, then validate genome
    _normalize_workflow_version_settings(cfg_in)
    _validate_enums_except_genome(cfg_in)
    _validate_profile_settings(cfg_in)
    _validate_pipeline_version_settings(cfg_in)
    _validate_backend_parameter_settings(cfg_in)
    _validate_resource_settings(cfg_in)
    _validate_nfcore_settings(cfg_in)
    _apply_genome_rules(cfg_in, user_provided_genome=("genome" in params))
    _validate_enum("genome", cfg_in["genome"], GENOME_VALUES)
    _validate_combos(cfg_in)
    return cfg_in


def _build_runtime_identity(cfg_in: dict) -> dict:
    """Build run identifier and project directory metadata."""
    try:
        fallback_user = getpass.getuser()
    except Exception:
        fallback_user = "unknown"
    user = os.environ.get("LOGNAME") or os.environ.get("USER") or fallback_user
    now = int(time.time())
    pid = os.getpid()
    run_id = f"{now}{pid % 100000:05d}"
    run_date = time.ctime()

    if cfg_in["gatk_version"] == "nf-core":
        name_parts = [cfg_in["project_dir"], "nf-core", cfg_in["pipeline"], cfg_in["mode"], run_id]
    else:
        name_parts = [
            cfg_in["project_dir"],
            cfg_in["workflow_engine"],
            cfg_in["pipeline"],
            cfg_in["mode"],
            cfg_in["genome"],
            cfg_in["gatk_version"],
            run_id,
        ]
    tmp_str = "_".join(name_parts)

    input_dir = cfg_in.get("input_dir")
    output_basename = None
    if input_dir:
        input_path = Path(input_dir).resolve()
        project_dir = str(input_path / tmp_str)
        output_basename = input_path.name
    else:
        project_dir = str(Path(tmp_str).resolve())
    return {
        "user": user,
        "run_id": run_id,
        "date": run_date,
        "project_dir": project_dir,
        "output_basename": output_basename,
    }


def _build_host_runtime_metadata(cfg_in: dict) -> dict:
    """Build host- and runtime-derived metadata unrelated to workflow resolution."""
    metadata: dict = {}
    metadata["hostname"] = socket.gethostname()

    nproc_path = shutil.which("nproc")
    if nproc_path:
        try:
            host_threads = int(os.popen(f"{nproc_path}").read().strip() or "1")
        except Exception:
            host_threads = os.cpu_count() or 1
    else:
        host_threads = os.cpu_count() or 1

    metadata["host_threads"] = host_threads
    metadata["host_threads_minus_one"] = host_threads - 1 if host_threads > 1 else 1

    host_threads_minus_one = metadata["host_threads_minus_one"]

    if os.access("/usr/bin/pigz", os.X_OK):
        metadata["compression_cmd"] = f"/usr/bin/pigz -p {host_threads_minus_one}"
    else:
        metadata["compression_cmd"] = "/bin/gunzip"

    if cfg_in["gatk_version"] != "nf-core" and cfg_in["pipeline"] != "wgs":
        if cfg_in["pipeline"] == "mit":
            metadata["capture_label"] = f"MToolBox_{cfg_in['genome']}"
        elif cfg_in["gatk_version"] == "gatk-3.5":
            metadata["capture_label"] = "Agilent SureSelect"
        else:
            metadata["capture_label"] = f"GATK_bundle_{cfg_in['genome']}"

    uname = platform.machine()
    if uname == "x86_64":
        arch = "x86_64"
    elif uname == "aarch64":
        arch = "arm64"
    else:
        arch = uname
    metadata["arch"] = arch

    if cfg_in["pipeline"] == "mit" and uname in {"aarch64", "arm64"}:
        raise WorkflowResolutionError(f"mit_{cfg_in['mode']} cannot be performed with: {uname}")
    return metadata


def _apply_runtime_profile(cfg_in: dict, workflow: WorkflowSpec) -> WorkflowSpec:
    """Resolve profile-specific helper files for workflows that support them."""
    profile = cfg_in.get("profile", "local")
    if profile == "local":
        return workflow

    helper_overrides = workflow.profiles.get(profile)
    if helper_overrides is None:
        raise WorkflowResolutionError(
            f"profile='{profile}' is not declared for workflow "
            f"{workflow.engine}/{workflow.gatk_version}/{workflow.pipeline}/{workflow.mode}."
        )

    helpers = dict(workflow.helpers)
    helpers.update(helper_overrides)
    return WorkflowSpec(
        engine=workflow.engine,
        pipeline=workflow.pipeline,
        mode=workflow.mode,
        gatk_version=workflow.gatk_version,
        pipeline_version=workflow.pipeline_version,
        entrypoint=workflow.entrypoint,
        config_file=workflow.config_file,
        helpers=helpers,
        profiles=workflow.profiles,
        metadata=workflow.metadata,
    )


def build_resolved_config(params: dict) -> ResolvedConfig:
    """Build the typed resolved config from validated parameters and registry wiring."""
    cfg_in = _merge_and_validate_param_values(params)
    project_root = get_project_root(__file__)
    registry = resolve_registry_context(project_root)
    workflow = resolve_workflow_spec(cfg_in, registry, project_root)
    validate_resolved_workflow_files(workflow)
    workflow = _apply_runtime_profile(cfg_in, workflow)
    validate_resolved_workflow_files(workflow)
    bundle_resource = build_bundle_resource_metadata(cfg_in, project_root, workflow)
    if bundle_resource.get("type") in {None, "bundle"}:
        bundle_resource["runtime_check"] = validate_installed_bundle_resource(bundle_resource, workflow)
    else:
        bundle_resource.setdefault(
            "runtime_check",
            {"status": "not_applicable", "reason": "resource type does not use DATADIR bundle checks"},
        )

    config = ResolvedConfig(
        user="",
        workflow_engine=cfg_in["workflow_engine"],
        profile=cfg_in["profile"],
        snakemake_parameters=cfg_in.get("snakemake_parameters", {}),
        nextflow_parameters=cfg_in.get("nextflow_parameters", {}),
        nfcore_profile=cfg_in.get("nfcore_profile"),
        nfcore_parameters=cfg_in.get("nfcore_parameters", {}),
        nfcore_singularity_cache_dir=cfg_in.get("nfcore_singularity_cache_dir"),
        genome=cfg_in["genome"],
        pipeline=cfg_in["pipeline"],
        mode=cfg_in["mode"],
        gatk_version=cfg_in["gatk_version"],
        pipeline_version=workflow.pipeline_version,
        inputs=InputsSpec(
            input_dir=cfg_in.get("input_dir"),
            sample_map=cfg_in.get("sample_map"),
        ),
        workflow=workflow,
        run_id="",
        date="",
        project_dir="",
        output_basename=None,
        hostname="",
        host_threads=0,
        host_threads_minus_one=0,
        compression_cmd="",
        run_mode="partial" if cfg_in.get("snakemake_parameters", {}).get("target") else "full",
        resources={"bundle": bundle_resource},
    )
    runtime_identity = _build_runtime_identity(cfg_in)
    host_metadata = _build_host_runtime_metadata(cfg_in)
    return ResolvedConfig(
        user=runtime_identity["user"],
        workflow_engine=config.workflow_engine,
        profile=config.profile,
        snakemake_parameters=config.snakemake_parameters,
        nextflow_parameters=config.nextflow_parameters,
        nfcore_profile=config.nfcore_profile,
        nfcore_parameters=config.nfcore_parameters,
        nfcore_singularity_cache_dir=config.nfcore_singularity_cache_dir,
        genome=config.genome,
        pipeline=config.pipeline,
        mode=config.mode,
        gatk_version=config.gatk_version,
        pipeline_version=config.pipeline_version,
        inputs=config.inputs,
        workflow=config.workflow,
        run_id=runtime_identity["run_id"],
        date=runtime_identity["date"],
        project_dir=runtime_identity["project_dir"],
        output_basename=runtime_identity["output_basename"],
        hostname=host_metadata["hostname"],
        host_threads=host_metadata["host_threads"],
        host_threads_minus_one=host_metadata["host_threads_minus_one"],
        compression_cmd=host_metadata["compression_cmd"],
        run_mode=config.run_mode,
        capture_label=host_metadata.get("capture_label"),
        arch=host_metadata.get("arch"),
        resources=config.resources,
    )


def set_config_values(params: dict) -> dict:
    """
    Build internal config structure from parameters.

    Important: this function merges defaults so callers may pass only overrides.
    """
    resolved_config = build_resolved_config(params).to_dict()
    return resolved_config
