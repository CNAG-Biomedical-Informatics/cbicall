import os
import time
import shutil
import socket
import platform
import getpass
from pathlib import Path

import yaml


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
    "genome": "b37",  # default
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
    if cfg.get("workflow_engine") == "snakemake" and cfg.get("gatk_version") == "gatk-3.5":
        raise ValueError(
            "workflow_engine='snakemake' is not supported for gatk_version='gatk-3.5'. "
            "Use workflow_engine='bash' or select gatk_version='gatk-4.6'."
        )

    # Block unsupported union: mit + snakemake
    if cfg.get("pipeline") == "mit" and cfg.get("workflow_engine") == "snakemake":
        raise ValueError(
            "The combination pipeline='mit' with workflow_engine='snakemake' is not supported."
        )

    # MIT pipeline forces genome to rsrs
    # - If user omits genome -> set it automatically
    # - If user sets genome to something else -> error
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


def _validate_enums(cfg: dict) -> None:
    _validate_enum("mode", cfg["mode"], MODE_VALUES)
    _validate_enum("pipeline", cfg["pipeline"], PIPELINE_VALUES)
    _validate_enum("organism", cfg["organism"], ORGANISM_VALUES)
    _validate_enum("technology", cfg["technology"], TECHNOLOGY_VALUES)
    _validate_enum("workflow_engine", cfg["workflow_engine"], WORKFLOW_ENGINE_VALUES)
    _validate_enum("gatk_version", cfg["gatk_version"], GATK_VALUES)
    _validate_enum("genome", cfg["genome"], GENOME_VALUES)


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

    # Enum-style validations
    _validate_enums(cfg)

    # Apply shared genome rules
    user_provided_genome = "genome" in param
    _apply_genome_rules(cfg, user_provided_genome=user_provided_genome)

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
    # Merge defaults so missing keys (like genome) never raise KeyError
    cfg_in = dict(_DEFAULTS)
    cfg_in.update(param)

    # Validate enums / combos consistently (same behavior as YAML path)
    _validate_enums(cfg_in)
    _apply_genome_rules(cfg_in, user_provided_genome=("genome" in param))
    _validate_combos(cfg_in)

    try:
        fallback_user = getpass.getuser()
    except Exception:
        fallback_user = "unknown"

    user = os.environ.get("LOGNAME") or os.environ.get("USER") or fallback_user

    # Base directories for workflows (using this file as anchor)
    here = Path(__file__).resolve()
    # src/cbicall/config.py -> src/cbicall -> src -> project root
    project_root = here.parents[2]
    workflows_bash_dir = (project_root / "workflows" / "bash").resolve()
    workflows_snakemake_dir = (project_root / "workflows" / "snakemake").resolve()

    # Versioned subdirs
    bash_version_dir = workflows_bash_dir / cfg_in["gatk_version"]
    snakemake_version_dir = workflows_snakemake_dir / cfg_in["gatk_version"]

    config = {"user": user}
    config["workflow_engine"] = cfg_in["workflow_engine"]
    config["genome"] = cfg_in["genome"]

    # Common bash scripts
    config["bash_parameters"] = str(bash_version_dir / "parameters.sh")
    config["bash_coverage"] = str(bash_version_dir / "coverage.sh")
    config["bash_jaccard"] = str(bash_version_dir / "jaccard.sh")
    config["bash_vcf2sex"] = str(bash_version_dir / "vcf2sex.sh")

    # Selected workflow name
    suffix = f"{cfg_in['pipeline']}_{cfg_in['mode']}"

    # Selected bash workflow (one per pipeline-mode)
    bash_script = f"{suffix}.sh"
    bash_key = f"bash_{suffix}"
    config[bash_key] = str(bash_version_dir / bash_script)

    # Snakemake workflows (only if chosen engine)
    if cfg_in["workflow_engine"] == "snakemake":
        smk_script = f"{suffix}.smk"
        smk_key = f"smk_{suffix}"
        config[smk_key] = str(snakemake_version_dir / smk_script)
        config["smk_config"] = str(snakemake_version_dir / "config.yaml")

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

    # Capture label (optional; keep only if used elsewhere)
    if cfg_in["pipeline"] == "wes":
        config["capture"] = "Agilent SureSelect"
    elif cfg_in["pipeline"] == "mit":
        # MIT uses MToolBox reference naming
        config["capture"] = "MToolBox_rsrs"
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

    # Validate executable permissions:
    exe_keys = [
        "bash_parameters",
        "bash_coverage",
        "bash_jaccard",
        "bash_vcf2sex",
    ]
    if cfg_in["workflow_engine"] == "bash":
        exe_keys.append(bash_key)

    missing = [k for k in exe_keys if k in config and not os.access(config[k], os.X_OK)]
    if missing:
        raise RuntimeError(
            "Missing +x on one or more workflow scripts: "
            + ", ".join(f"{k} -> {config[k]}" for k in missing)
        )

    return config

