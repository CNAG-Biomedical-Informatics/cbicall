from pathlib import Path

from .models import ResolvedConfig

ARROW = "=>"


def _plain(value) -> str:
    return str(value) if value is not None and str(value) != "" else "(undef)"


def _short_path(value) -> str:
    if value is None or str(value) == "":
        return "(undef)"

    path = Path(str(value))
    display_path = path
    try:
        home = Path.home()
        resolved = path.resolve()
        try:
            return str(Path("~") / resolved.relative_to(home))
        except ValueError:
            pass
    except Exception:
        pass

    parts = display_path.parts
    if len(parts) <= 4:
        return str(display_path)

    return str(Path("...") / Path(*parts[-3:]))


def _format_duration(seconds: float) -> str:
    total = max(0, int(round(seconds)))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes}m {secs}s"
    if minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


def _section(title: str, color: str, bold: str, reset: str) -> None:
    print(f"{bold}{color}{title}{reset}")


def _row(label: str, value) -> None:
    print(f"  {label:<12} {ARROW} {_plain(value)}")


def _warn(message: str, bold: str, yellow: str, reset: str) -> None:
    print(f"{bold}{yellow}Warning{reset} {message}")


def _print_config(resolved_config: ResolvedConfig, bold: str, blue: str, reset: str) -> None:
    _section("Resolved Configuration", blue, bold, reset)
    data = resolved_config.to_dict()
    keys = list(data.keys())
    if not keys:
        return
    max_key = max(len(k) for k in keys)
    for key in sorted(keys):
        print(f"  {key:<{max_key}} {ARROW} {_plain(data.get(key))}")


def _print_params(param: dict, bold: str, green: str, reset: str) -> None:
    _section("Input Parameters", green, bold, reset)
    max_key = max(len(k) for k in param.keys())
    for key in sorted(param.keys()):
        print(f"  {key:<{max_key}} {ARROW} {_plain(param.get(key))}")


def _print_run_summary(
    *,
    arg: dict,
    resolved_config: ResolvedConfig,
    cbicall_path: Path,
    version: str,
    colors: dict,
) -> None:
    workflow = resolved_config.workflow
    genome = resolved_config.genome or "b37"
    _section(f"CBIcall {version}", colors["cyan"], colors["bold"], colors["reset"])
    _row("Executable", _short_path(cbicall_path))
    _row("Workflow", f"{workflow.backend} -> {workflow.pipeline} -> {workflow.mode}")
    _row("Runtime profile", resolved_config.profile)
    _row("Genome", genome)
    _row("Threads", arg.get("threads"))
    _row("Project", _short_path(resolved_config.project_dir))
    _row("Run ID", resolved_config.run_id)
    print()

    if resolved_config.run_mode == "partial":
        _section("Partial Run", colors["yellow"], colors["bold"], colors["reset"])
        _warn(
            "This run starts from a workflow rule and does not execute the full workflow.",
            colors["bold"],
            colors["yellow"],
            colors["reset"],
        )
        _row("Run mode", resolved_config.run_mode)
        _row("Snakemake target", resolved_config.snakemake_parameters.get("target"))
        print()

    _section("Inputs", colors["yellow"], colors["bold"], colors["reset"])
    _row("Param file", _short_path(arg.get("paramfile")))
    _row("Input dir", _short_path(resolved_config.inputs.input_dir))
    _row("Sample map", _short_path(resolved_config.inputs.sample_map))
    _row("Workflow provider", resolved_config.workflow_provider)
    if resolved_config.workflow_provider == "cbicall":
        _row("GATK", workflow.gatk_version)
    _row("Pipeline ver", workflow.pipeline_version)
    if workflow.metadata.get("provider") == "nf-core":
        _row("NF profile", resolved_config.nfcore_profile)
        _row("NF parameters", ", ".join(sorted(resolved_config.nfcore_parameters)) or "(none)")
        if resolved_config.nfcore_singularity_cache_dir:
            _row("NF cache", _short_path(resolved_config.nfcore_singularity_cache_dir))
    print()

    _section("Resolved", colors["blue"], colors["bold"], colors["reset"])
    if workflow.backend == "bash":
        _row("Entrypoint", _short_path(workflow.entrypoint))
        _row("Env file", _short_path(workflow.helpers.get("env")))
    elif workflow.backend == "snakemake":
        _row("Snakefile", _short_path(workflow.entrypoint))
        _row("Config", _short_path(workflow.config_file))
    elif workflow.backend == "nextflow":
        if workflow.metadata.get("provider") == "nf-core":
            _row("Nextflow", workflow.metadata.get("source"))
            _row("Release", workflow.metadata.get("release"))
            _row("Outdir", Path(resolved_config.project_dir) / workflow.metadata.get("default_outdir", workflow.pipeline))
        else:
            _row("Nextflow", _short_path(workflow.entrypoint))
            _row("Config", _short_path(workflow.config_file))

    if workflow.metadata.get("provider") == "nf-core":
        log_name = f"nf-core_{workflow.pipeline}_{workflow.mode}.log"
    else:
        log_name = f"{workflow.backend}_{workflow.pipeline}_{workflow.mode}_{genome}_{workflow.gatk_version}.log"
    _row("Log", Path(resolved_config.project_dir) / log_name)
    print()
