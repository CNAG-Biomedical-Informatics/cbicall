import os
import time
import shutil
import socket
import platform
import getpass
import json
from typing import Tuple
from pathlib import Path

import yaml
from jsonschema import Draft202012Validator


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
ORGANISM_VALUES = {"Homo Sapiens", "Mus musculus"}
TECHNOLOGY_VALUES = {"Illumina HiSeq", "NovaSeq"}
WORKFLOW_ENGINE_VALUES = {"bash", "nextflow", "snakemake"}
GATK_VALUES = {"gatk-3.5", "gatk-4.6"}
GENOME_VALUES = {"b37", "hg38", "rsrs"}


# Defaults
_DEFAULTS = {
    "mode": "single",
    "sample": None,
    "sample_map": None,
    "output_basename": None,
    "pipeline": "wes",
    "organism": "Homo Sapiens",
    "technology": "Illumina HiSeq",
    "workflow_engine": "bash",
    "gatk_version": "gatk-3.5",
    "projectdir": "cbicall",
    "cleanup_bam": False,
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
        raise ValueError(
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
        raise ValueError(
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
        raise ValueError(
            "workflow_engine='snakemake' is not supported for gatk_version='gatk-3.5'. "
            "Use workflow_engine='bash' or select gatk_version='gatk-4.6'."
        )

    # Block unsupported union: mit + snakemake
    if cfg.get("pipeline") == "mit" and cfg.get("workflow_engine") == "snakemake":
        raise ValueError(
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
            raise ValueError(
                "For pipeline='mit', genome is fixed to 'rsrs'. "
                "Remove 'genome' from the YAML or set genome='rsrs'."
            )
        cfg["genome"] = "rsrs"

    # hg38 is supported only for WGS
    if cfg["genome"] == "hg38" and cfg["pipeline"] != "wgs":
        raise ValueError("genome='hg38' is only supported for pipeline='wgs'.")


def _validate_enums_except_genome(cfg: dict) -> None:
    _validate_enum("mode", cfg["mode"], MODE_VALUES)
    _validate_enum("pipeline", cfg["pipeline"], PIPELINE_VALUES)
    _validate_enum("organism", cfg["organism"], ORGANISM_VALUES)
    _validate_enum("technology", cfg["technology"], TECHNOLOGY_VALUES)
    _validate_enum("workflow_engine", cfg["workflow_engine"], WORKFLOW_ENGINE_VALUES)
    _validate_enum("gatk_version", cfg["gatk_version"], GATK_VALUES)


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
        raise ValueError("\n".join(lines))


def load_workflow_registry(
    registry_yaml: Path,
    schema_json: Path,
) -> dict:
    with registry_yaml.open("r") as fh:
        data = yaml.safe_load(fh) or {}
    schema = _load_json(schema_json)
    _validate_with_schema(data, schema, label=str(registry_yaml))
    return data


def _get_project_root() -> Path:
    here = Path(__file__).resolve()
    # src/cbicall/config.py -> src/cbicall -> src -> project root
    return here.parents[2]


def _registry_paths(project_root: Path) -> Tuple[Path, Path]:
    registry_yaml = (project_root / "workflows" / "config" / "cbicall.workflows.yaml").resolve()
    schema_json = (project_root / "workflows" / "schema" / "workflows.schema.json").resolve()
    return registry_yaml, schema_json


def read_param_file(yaml_file: str) -> dict:
    """
    Load YAML parameters, merge with defaults and validate.
    """
    yaml_path = Path(yaml_file)
    if not yaml_path.is_file():
        raise FileNotFoundError(f"Parameters file not found: {yaml_file}")

    with yaml_path.open("r") as fh:
        param = yaml.safe_load(fh) or {}

    # Start from defaults
    cfg = dict(_DEFAULTS)

    # Merge provided parameters and validate keys
    for key, value in param.items():
        if key not in cfg:
            raise ValueError(f"Parameter '{key}' does not exist (typo?)")
        cfg[key] = value

    # Validate enums (except genome; it may be None until rules apply)
    _validate_enums_except_genome(cfg)

    # Apply shared genome rules
    user_provided_genome = "genome" in param
    _apply_genome_rules(cfg, user_provided_genome=user_provided_genome)

    # Now validate genome
    _validate_enum("genome", cfg["genome"], GENOME_VALUES)

    # Make sample_map absolute (resolve relative to YAML location)
    if cfg.get("sample_map") is not None:
        sm = Path(cfg["sample_map"])
        if not sm.is_absolute():
            sm = yaml_path.parent / sm
        cfg["sample_map"] = str(sm.resolve())

    # Validate pipeline-mode combination for selected GATK version
    _validate_combos(cfg)

    return cfg


def set_config_values(param: dict) -> dict:
    """
    Build internal config structure from parameters.

    Important: this function merges defaults so callers may pass only overrides.
    """
    # Merge defaults so missing keys never raise KeyError
    cfg_in = dict(_DEFAULTS)
    cfg_in.update(param)

    # Validate enums (except genome), then apply genome rules, then validate genome
    _validate_enums_except_genome(cfg_in)
    _apply_genome_rules(cfg_in, user_provided_genome=("genome" in param))
    _validate_enum("genome", cfg_in["genome"], GENOME_VALUES)

    # Validate pipeline/mode combo for selected GATK version
    _validate_combos(cfg_in)

    try:
        fallback_user = getpass.getuser()
    except Exception:
        fallback_user = "unknown"

    user = os.environ.get("LOGNAME") or os.environ.get("USER") or fallback_user

    # Project root
    project_root = _get_project_root()

    # Load workflow registry + schema validate
    registry_yaml, schema_json = _registry_paths(project_root)
    if not registry_yaml.is_file():
        raise FileNotFoundError(f"Workflow registry not found: {registry_yaml}")
    if not schema_json.is_file():
        raise FileNotFoundError(f"Workflow schema not found: {schema_json}")

    registry = load_workflow_registry(registry_yaml, schema_json)

    config: dict = {"user": user}
    config["workflow_engine"] = cfg_in["workflow_engine"]
    config["genome"] = cfg_in["genome"]

    # Resolve workflow scripts from registry
    engine = cfg_in["workflow_engine"]
    version = cfg_in["gatk_version"]
    pipeline = cfg_in["pipeline"]
    mode = cfg_in["mode"]
    suffix = f"{pipeline}_{mode}"

    workflows = registry["workflows"]
    if engine not in workflows:
        raise ValueError(f"Engine not defined in workflow registry: {engine}")

    eng_cfg = workflows[engine]
    versions_cfg = eng_cfg["versions"]
    if version not in versions_cfg:
        raise ValueError(f"Version not defined for engine '{engine}': {version}")

    ver_cfg = versions_cfg[version]
    common = ver_cfg.get("common", {})
    pipelines_cfg = ver_cfg["pipelines"]

    if pipeline not in pipelines_cfg:
        raise ValueError(f"Pipeline not defined for {engine}/{version}: {pipeline}")
    if mode not in pipelines_cfg[pipeline]:
        raise ValueError(
            f"Mode not defined for pipeline '{pipeline}' in {engine}/{version}: {mode}"
        )

    base_dir = (project_root / eng_cfg["base_dir"] / version).resolve()

    # Keep existing key names so the rest of the code stays unchanged
    if engine == "bash":
        # Required common keys for bash
        needed_common = ["parameters", "coverage", "jaccard", "vcf2sex"]
        missing_common = [k for k in needed_common if k not in common]
        if missing_common:
            raise ValueError(
                f"Workflow registry is missing common keys for bash/{version}: {missing_common}"
            )

        config["bash_parameters"] = str(base_dir / common["parameters"])
        config["bash_coverage"] = str(base_dir / common["coverage"])
        config["bash_jaccard"] = str(base_dir / common["jaccard"])
        config["bash_vcf2sex"] = str(base_dir / common["vcf2sex"])

        script_name = pipelines_cfg[pipeline][mode]
        config[f"bash_{suffix}"] = str(base_dir / script_name)

    elif engine == "snakemake":
        # Required common keys for snakemake
        if "config" not in common:
            raise ValueError(
                f"Workflow registry is missing common key 'config' for snakemake/{version}"
            )

        script_name = pipelines_cfg[pipeline][mode]
        config[f"smk_{suffix}"] = str(base_dir / script_name)
        config["smk_config"] = str(base_dir / common["config"])

    elif engine == "nextflow":
        # You can extend this later; keep explicit for now
        raise ValueError("workflow_engine='nextflow' is declared but not implemented.")

    else:
        raise ValueError(f"Unsupported workflow_engine: {engine}")

    # Internal settings
    now = int(time.time())
    pid = os.getpid()
    config["id"] = f"{now}{pid % 100000:05d}"
    config["date"] = time.ctime()

    tmp_str = "_".join(
        [
            cfg_in["projectdir"],
            cfg_in["workflow_engine"],
            cfg_in["pipeline"],
            cfg_in["mode"],
            cfg_in["genome"],
            cfg_in["gatk_version"],
            config["id"],
        ]
    )

    sample = cfg_in.get("sample")
    if sample:
        sample_path = Path(sample).resolve()
        config["projectdir"] = str(sample_path / tmp_str)
        config["output_basename"] = sample_path.name
    else:
        config["projectdir"] = str(Path(tmp_str).resolve())

    # Host info
    config["hostname"] = socket.gethostname()

    # Thread counts
    nproc_path = shutil.which("nproc")
    if nproc_path:
        try:
            threadshost = int(os.popen(f"{nproc_path}").read().strip() or "1")
        except Exception:
            threadshost = os.cpu_count() or 1
    else:
        threadshost = os.cpu_count() or 1

    config["threadshost"] = threadshost
    config["threadsless"] = threadshost - 1 if threadshost > 1 else 1

    th_less = config["threadsless"]

    # Compression
    if os.access("/usr/bin/pigz", os.X_OK):
        config["zip"] = f"/usr/bin/pigz -p {th_less}"
    else:
        config["zip"] = "/bin/gunzip"

    # Capture label
    if cfg_in["pipeline"] != "wgs":
        if cfg_in["pipeline"] == "mit":
            config["capture"] = f"MToolBox_{config['genome']}"
        elif cfg_in["gatk_version"] == "gatk-3.5":
            config["capture"] = "Agilent SureSelect"
        else:
            config["capture"] = f"GATK_bundle_{config['genome']}"

    # Architecture
    uname = platform.machine()
    if uname == "x86_64":
        arch = "x86_64"
    elif uname == "aarch64":
        arch = "arm64"
    else:
        arch = uname
    config["arch"] = arch

    # FAIL FAST: MIT not supported on arm64/aarch64
    if cfg_in["pipeline"] == "mit" and uname in {"aarch64", "arm64"}:
        raise RuntimeError(f"mit_{cfg_in['mode']} cannot be performed with: {uname}")

    # Extra guardrail: referenced files must exist
    must_exist = []
    if engine == "bash":
        must_exist += [
            "bash_parameters",
            "bash_coverage",
            "bash_jaccard",
            "bash_vcf2sex",
            f"bash_{suffix}",
        ]
    elif engine == "snakemake":
        must_exist += [f"smk_{suffix}", "smk_config"]

    missing_files = [k for k in must_exist if k in config and not Path(config[k]).exists()]
    if missing_files:
        raise RuntimeError(
            "One or more workflow files referenced in the registry do not exist: "
            + ", ".join(f"{k} -> {config[k]}" for k in missing_files)
        )

    # Validate executable permissions (+x) for bash scripts (your original behavior)
    if engine == "bash":
        exe_keys = [
            "bash_parameters",
            "bash_coverage",
            "bash_jaccard",
            "bash_vcf2sex",
            f"bash_{suffix}",
        ]
        not_exe = [k for k in exe_keys if k in config and not os.access(config[k], os.X_OK)]
        if not_exe:
            raise RuntimeError(
                "Missing +x on one or more workflow scripts: "
                + ", ".join(f"{k} -> {config[k]}" for k in not_exe)
            )

    return config
