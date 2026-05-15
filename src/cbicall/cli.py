import json
import os
import argparse
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import List

from . import config as config_mod
from .cli_output import (
    _format_duration,
    _print_config,
    _print_params,
    _print_run_summary,
    _short_path,
    _warn,
    _row,
    _section,
)
from .dnaseq import DNAseq
from .helpmod import usage, parse_args as _parse_args
from .models import ResolvedConfig, RunSettings
from .workflow_registry import load_workflow_registry
from .goodbye import GoodBye

VERSION = "1.0.0"
PROMPT = "Info:"
SPACER = "*" * 41
ARROW = "=>"
AUTHOR = "Author: Manuel Rueda, PhD"
LICENSE = "License: GNU General Public License v3"


# ANSI color helpers (respect ANSI_COLORS_DISABLED)
def _colors_enabled() -> bool:
    return not os.environ.get("ANSI_COLORS_DISABLED")


def _code(s: str) -> str:
    return s if _colors_enabled() else ""


def _refresh_colors() -> None:
    global BOLD, RESET, RED, YELLOW, GREEN, BLUE, CYAN, WHITE
    BOLD = _code("\033[1m")
    RESET = _code("\033[0m")
    RED = _code("\033[31m")
    YELLOW = _code("\033[33m")
    GREEN = _code("\033[32m")
    BLUE = _code("\033[34m")
    CYAN = _code("\033[36m")
    WHITE = _code("\033[37m")


_refresh_colors()


def parse_args(argv):
    """
    Convenience wrapper so tests can do cli.parse_args([...])
    and get an argparse.Namespace.
    """
    return _parse_args(argv, VERSION)


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _section(title: str, color: str = WHITE) -> None:
    from .cli_output import _section as _render_section

    _render_section(title, color, BOLD, RESET)


def _print_config(resolved_config) -> None:
    if isinstance(resolved_config, ResolvedConfig):
        from .cli_output import _print_config as _render_config

        _render_config(resolved_config, BOLD, BLUE, RESET)
        return

    print(f"{BOLD}{BLUE}Resolved Configuration{RESET}")
    keys = list(resolved_config.keys())
    if not keys:
        return
    max_key = max(len(k) for k in keys)
    for key in sorted(keys):
        print(f"  {key:<{max_key}} {ARROW} {resolved_config.get(key) if resolved_config.get(key) is not None else '(undef)'}")


def _print_params(param: dict) -> None:
    from .cli_output import _print_params as _render_params

    _render_params(param, BOLD, GREEN, RESET)


def write_log(resolved_config: dict, arg: dict, params: dict) -> None:
    """
    Create log.json in project_dir with arg, config, param.
    """
    project_dir = Path(resolved_config["project_dir"])
    project_dir.mkdir(parents=True, exist_ok=True)
    file_path = project_dir / "log.json"
    payload = {
        "arg": arg,
        "config": resolved_config,
        "param": params,
    }
    with file_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2, sort_keys=True)


def write_run_report(
    resolved_config: ResolvedConfig,
    arg: dict,
    params: dict,
    *,
    elapsed_seconds: float,
    workflow_log: Path,
) -> Path:
    """Write a compact execution report next to log.json."""
    project_dir = Path(resolved_config.project_dir)
    project_dir.mkdir(parents=True, exist_ok=True)
    report_path = project_dir / "run-report.json"
    workflow = resolved_config.workflow
    payload = {
        "status": "success",
        "elapsed_seconds": round(elapsed_seconds, 3),
        "workflow_log": str(workflow_log),
        "command_trace": "Detailed workflow stdout/stderr is recorded in the workflow log.",
        "profile": resolved_config.profile,
        "workflow": workflow.to_dict(),
        "resource_bundle": resolved_config.resource_bundle,
        "run": {
            "run_id": resolved_config.run_id,
            "project_dir": resolved_config.project_dir,
            "threads": arg.get("threads"),
            "genome": resolved_config.genome,
        },
        "inputs": resolved_config.inputs.to_dict(),
        "parameters": {
            "paramfile": arg.get("paramfile"),
            "cleanup_bam": params.get("cleanup_bam", False),
            "workflow_rule": resolved_config.workflow_rule,
            "allow_partial_run": resolved_config.allow_partial_run,
        },
    }
    with report_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2, sort_keys=True)
    return report_path


def _run_doctor_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall doctor",
        description="Resolve and check a CBIcall run without starting the workflow.",
    )
    parser.add_argument("-p", "--param", dest="paramfile", required=True, help="Parameters input file (YAML).")
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    params = config_mod.read_param_file(args.paramfile)
    resolved_config = ResolvedConfig.from_mapping({**config_mod.set_config_values(params), "version": VERSION})

    _section("Configuration OK", GREEN)
    workflow = resolved_config.workflow
    bundle = resolved_config.resource_bundle
    _row("Param file", _short_path(args.paramfile))
    _row("Profile", resolved_config.profile)
    _row("Workflow", f"{workflow.engine} -> {workflow.pipeline} -> {workflow.mode}")
    _row("GATK", workflow.gatk_version)
    _row("Pipeline ver", workflow.pipeline_version)
    _row("Genome", resolved_config.genome or "b37")
    _row("Entrypoint", _short_path(workflow.entrypoint))
    if workflow.engine == "bash":
        _row("Env file", _short_path(workflow.helpers.get("env")))
    _row("Bundle", bundle.get("key"))
    _row("Bundle hash", bundle.get("fingerprint"))
    return 0


def _run_validate_registry_command(argv: List[str]) -> int:
    root = _project_root()
    default_registry = root / "workflows" / "registry" / "workflows.yaml"
    default_schema = root / "workflows" / "schema" / "workflows.schema.json"

    parser = argparse.ArgumentParser(
        prog="cbicall validate-registry",
        description="Validate the workflow registry YAML against its JSON Schema.",
    )
    parser.add_argument("--registry", default=str(default_registry), help="Workflow registry YAML.")
    parser.add_argument("--schema", default=str(default_schema), help="Workflow registry JSON Schema.")
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    registry_path = Path(args.registry)
    schema_path = Path(args.schema)
    registry = load_workflow_registry(registry_path, schema_path)
    engines = sorted((registry.get("workflows") or {}).keys())

    _section("Registry OK", GREEN)
    _row("Registry", _short_path(registry_path))
    _row("Schema", _short_path(schema_path))
    _row("Engines", ", ".join(engines) if engines else "(none)")
    return 0


def _run_test_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall test",
        description="Run the bundled CBIcall integration tests from examples/input.",
    )
    parser.add_argument("--wes", action="store_true", help="Run the WES integration test.")
    parser.add_argument("--mit", action="store_true", help="Run the mitochondrial integration test.")
    parser.add_argument("--all", action="store_true", help="Run all bundled integration tests.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use.")
    args = parser.parse_args(argv)

    if args.threads <= 0:
        parser.error("--threads requires a positive integer")

    selected = []
    if args.all or args.wes:
        selected.append("--wes")
    if args.all or args.mit:
        selected.append("--mit")
    if not selected:
        parser.error("select at least one test with --wes, --mit, or --all")

    root = _project_root()
    script = root / "examples" / "input" / "run_tests.sh"
    if not script.is_file():
        raise FileNotFoundError(f"Integration test launcher not found: {script}")

    env = os.environ.copy()
    env["CBICALL"] = str(root / "bin" / "cbicall")
    env["THREADS"] = str(args.threads)
    return subprocess.run(
        ["bash", str(script), *selected],
        cwd=str(script.parent),
        env=env,
        check=False,
    ).returncode


def run_with_spinner(func, *args, no_spinner: bool = False):
    """
    Run a callable with an optional UTF-8 spinner and elapsed time message.
    """
    if no_spinner or not sys.stdout.isatty():
        return func(*args)

    done = False

    def spinner():
        nonlocal done
        frames = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
        i = 0
        start = time.time()
        delay = 1.0

        while not done:
            spinner_char = f"{BOLD}{YELLOW}{frames[i % len(frames)]}{RESET}"
            elapsed = time.time() - start
            elapsed_str = _format_duration(elapsed)
            msg = f"{BOLD}{WHITE} Working (elapsed: {elapsed_str}){RESET}"
            sys.stdout.write("\r" + spinner_char + msg)
            sys.stdout.flush()
            i += 1
            time.sleep(delay)

        # Clear the line
        sys.stdout.write("\r\033[2K")
        sys.stdout.flush()

    t = threading.Thread(target=spinner, daemon=True)
    t.start()

    error = None
    result = None
    try:
        result = func(*args)
    except Exception as e:
        error = e
    finally:
        done = True
        t.join()

    if error is not None:
        raise error

    return result


def main() -> int:
    start_time = time.time()
    cbicall_path = Path(sys.argv[0]).resolve()

    if len(sys.argv) > 1 and sys.argv[1] == "doctor":
        return _run_doctor_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "validate-registry":
        return _run_validate_registry_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        return _run_test_command(sys.argv[2:])

    # Parse CLI args (Help::usage)
    arg = usage(VERSION)
    if arg.get("nocolor"):
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()
    no_spinner = bool(arg.get("debug") or arg.get("verbose"))

    # Read parameters and build config (Config::read_param_file + set_config_values)
    params = config_mod.read_param_file(arg["paramfile"])
    resolved_config = ResolvedConfig.from_mapping({**config_mod.set_config_values(params), "version": VERSION})

    if params.get("genome") is None and resolved_config.genome is not None:
        _warn(f"Genome not provided; using inferred default '{resolved_config.genome}'.", BOLD, YELLOW, RESET)

    _print_run_summary(
        arg=arg,
        resolved_config=resolved_config,
        cbicall_path=cbicall_path,
        version=VERSION,
        colors={"bold": BOLD, "reset": RESET, "blue": BLUE, "cyan": CYAN, "green": GREEN, "yellow": YELLOW},
    )

    if arg.get("verbose") or arg.get("debug"):
        _print_config(resolved_config)
        print()
        _print_params(params)
        print()

    # Start CBIcall
    _section("Running", CYAN)

    # Create working dir and log
    Path(resolved_config.project_dir).mkdir(parents=True, exist_ok=True)
    write_log(resolved_config.to_dict(), arg, params)

    # Build settings hash (rah_cbicall)
    settings = RunSettings.from_mapping(
        {
            "project_dir": resolved_config.project_dir,
            "run_id": resolved_config.run_id,
            "threads": arg["threads"],
            "debug": arg["debug"],
            "profile": resolved_config.profile,
            "genome": resolved_config.genome,
            "cleanup_bam": params.get("cleanup_bam", False),
            "workflow_rule": resolved_config.workflow_rule,
            "allow_partial_run": resolved_config.allow_partial_run,
            "run_mode": resolved_config.run_mode,
            "inputs": resolved_config.inputs.to_dict(),
            "workflow": resolved_config.workflow.to_dict(),
        }
    )

    workflow = resolved_config.workflow
    _row("Workflow", f"{workflow.engine} -> {workflow.pipeline} -> {workflow.mode}")
    print("  This workflow may take a while depending on input size and pipeline.")

    wes = DNAseq(settings)

    # Run with spinner
    run_with_spinner(
        wes.variant_calling,
        no_spinner=no_spinner,
    )

    # END CBICALL
    print()
    _section("Completed", GREEN)
    _row("Status", "Finished successfully")
    elapsed = time.time() - start_time
    _row("Elapsed", _format_duration(elapsed))
    genome = resolved_config.genome or "b37"
    log_name = f"{workflow.engine}_{workflow.pipeline}_{workflow.mode}_{genome}_{workflow.gatk_version}.log"
    workflow_log = Path(resolved_config.project_dir) / log_name
    report_path = write_run_report(
        resolved_config,
        arg,
        params,
        elapsed_seconds=elapsed,
        workflow_log=workflow_log,
    )
    _row("Log", workflow_log)
    _row("Report", report_path)
    if arg.get("verbose"):
        _row("Date", time.ctime())

    gb = GoodBye()
    print(f"{WHITE}{gb.say_goodbye()}{RESET}")

    return 0
