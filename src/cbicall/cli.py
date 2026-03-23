import json
import os
import sys
import threading
import time
from pathlib import Path

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
    _row("Log", Path(resolved_config.project_dir) / log_name)
    if arg.get("verbose"):
        _row("Date", time.ctime())

    gb = GoodBye()
    print(f"{WHITE}{gb.say_goodbye()}{RESET}")

    return 0
