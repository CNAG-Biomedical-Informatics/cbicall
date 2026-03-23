import os
import time
import shutil
import socket
import platform
import getpass
from pathlib import Path

import yaml

from .errors import ParameterValidationError, WorkflowResolutionError
from .models import InputsSpec, ResolvedConfig
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
GENOME_VALUES = {"b37", "hg38", "rsrs"}


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
    "gatk_version": "gatk-3.5",
    "project_dir": "cbicall",
    "cleanup_bam": False,
    "workflow_rule": None,
    "allow_partial_run": False,
    "genome": None,  # Option 1: effective default assigned in _apply_genome_rules
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
    _validate_enum("pipeline", cfg["pipeline"], PIPELINE_VALUES)
    _validate_enum("organism", cfg["organism"], ORGANISM_VALUES)
    _validate_enum("technology", cfg["technology"], TECHNOLOGY_VALUES)
    _validate_enum("workflow_engine", cfg["workflow_engine"], WORKFLOW_ENGINE_VALUES)
    _validate_enum("gatk_version", cfg["gatk_version"], GATK_VALUES)


def _validate_partial_run_settings(cfg: dict) -> None:
    workflow_rule = cfg.get("workflow_rule")
    allow_partial_run = cfg.get("allow_partial_run")

    if workflow_rule is not None:
        if not isinstance(workflow_rule, str) or not workflow_rule.strip():
            raise ParameterValidationError(
                "workflow_rule must be a non-empty string when provided."
            )
        cfg["workflow_rule"] = workflow_rule.strip()

    if not isinstance(allow_partial_run, bool):
        raise ParameterValidationError("allow_partial_run must be a boolean value.")

    if cfg.get("workflow_rule") and not allow_partial_run:
        raise ParameterValidationError(
            f"workflow_rule was set to '{cfg['workflow_rule']}', but --allow-partial-run was not provided. "
            "Refusing to start a partial workflow run."
        )


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
    cfg = dict(_DEFAULTS)

    # Merge provided parameters and validate keys
    for key, value in params.items():
        if key not in cfg:
            raise ParameterValidationError(f"Parameter '{key}' does not exist (typo?)")
        cfg[key] = value

    # Validate enums (except genome; it may be None until rules apply)
    _validate_enums_except_genome(cfg)
    _validate_partial_run_settings(cfg)

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

    # Validate pipeline-mode combination for selected GATK version
    _validate_combos(cfg)

    return cfg


def _merge_and_validate_param_values(params: dict) -> dict:
    """Merge defaults with user parameters and run semantic validation."""
    cfg_in = dict(_DEFAULTS)
    cfg_in.update(params)
    # Validate enums (except genome), then apply genome rules, then validate genome
    _validate_enums_except_genome(cfg_in)
    _validate_partial_run_settings(cfg_in)
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

    tmp_str = "_".join(
        [
            cfg_in["project_dir"],
            cfg_in["workflow_engine"],
            cfg_in["pipeline"],
            cfg_in["mode"],
            cfg_in["genome"],
            cfg_in["gatk_version"],
            run_id,
        ]
    )

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

    if cfg_in["pipeline"] != "wgs":
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


def build_resolved_config(params: dict) -> ResolvedConfig:
    """Build the typed resolved config from validated parameters and registry wiring."""
    cfg_in = _merge_and_validate_param_values(params)
    project_root = get_project_root(__file__)
    registry = resolve_registry_context(project_root)
    workflow = resolve_workflow_spec(cfg_in, registry, project_root)
    validate_resolved_workflow_files(workflow)

    config = ResolvedConfig(
        user="",
        workflow_engine=cfg_in["workflow_engine"],
        genome=cfg_in["genome"],
        pipeline=cfg_in["pipeline"],
        mode=cfg_in["mode"],
        gatk_version=cfg_in["gatk_version"],
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
        workflow_rule=cfg_in.get("workflow_rule"),
        allow_partial_run=bool(cfg_in.get("allow_partial_run", False)),
        run_mode="partial" if cfg_in.get("workflow_rule") else "full",
    )
    runtime_identity = _build_runtime_identity(cfg_in)
    host_metadata = _build_host_runtime_metadata(cfg_in)
    return ResolvedConfig(
        user=runtime_identity["user"],
        workflow_engine=config.workflow_engine,
        genome=config.genome,
        pipeline=config.pipeline,
        mode=config.mode,
        gatk_version=config.gatk_version,
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
        workflow_rule=config.workflow_rule,
        allow_partial_run=config.allow_partial_run,
        run_mode=config.run_mode,
        capture_label=host_metadata.get("capture_label"),
        arch=host_metadata.get("arch"),
    )


def set_config_values(params: dict) -> dict:
    """
    Build internal config structure from parameters.

    Important: this function merges defaults so callers may pass only overrides.
    """
    resolved_config = build_resolved_config(params).to_dict()
    return resolved_config
