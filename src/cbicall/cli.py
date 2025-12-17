import json
import os
import sys
import threading
import time
from pathlib import Path

from . import config as config_mod
from .dnaseq import DNAseq
from .helpmod import usage, parse_args as _parse_args
from .goodbye import GoodBye

VERSION = "0.0.1"
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


BOLD = _code("\033[1m")
RESET = _code("\033[0m")
RED = _code("\033[31m")
YELLOW = _code("\033[33m")
GREEN = _code("\033[32m")
BLUE = _code("\033[34m")
CYAN = _code("\033[36m")
WHITE = _code("\033[37m")


def parse_args(argv):
    """
    Convenience wrapper so tests can do cli.parse_args([...])
    and get an argparse.Namespace.
    """
    return _parse_args(argv, VERSION)


def write_log(cfg: dict, arg: dict, param: dict) -> None:
    """
    Create log.json in projectdir with arg, config, param.
    """
    projectdir = Path(cfg["projectdir"])
    projectdir.mkdir(parents=True, exist_ok=True)
    file_path = projectdir / "log.json"
    payload = {
        "arg": arg,
        "config": cfg,
        "param": param,
    }
    with file_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2, sort_keys=True)


def run_with_spinner(func, *args, no_spinner: bool = False, no_emoji: bool = False):
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
            spinner_char = f"{BOLD}{YELLOW}{PROMPT} {frames[i % len(frames)]}{RESET}"
            elapsed = time.time() - start
            seconds = int(elapsed % 60)
            minutes = int((elapsed / 60) % 60)
            hours = int(elapsed / 3600)
            elapsed_str = f"{hours:02d}h {minutes:02d}m {seconds:02d}s"
            msg = (
                f"{BOLD}{WHITE} Please be patient - this job may take hours"
                f"{'' if no_emoji else ' ⏳'} (elapsed: {elapsed_str})...{RESET}"
            )
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


def _print_config(cfg: dict, no_emoji: bool) -> None:
    print(BOLD + BLUE + f"{PROMPT} CONFIGURATION VALUES:" + RESET)
    print(WHITE + PROMPT + RESET)

    engine = cfg.get("workflow_engine")

    def _keep_key(k: str) -> bool:
        # If snakemake is selected, hide bash workflow script paths from STDOUT
        if engine == "snakemake":
            if k.startswith("bash_"):
                return False
        return True

    keys = [k for k in cfg.keys() if _keep_key(k)]
    if not keys:
        return

    max_key = max(len(k) for k in keys)
    for key in sorted(keys):
        val = cfg.get(key)
        val_disp = str(val) if val is not None and str(val) != "" else "(undef)"
        line = f"{PROMPT:>5} {key:<{max_key}} {ARROW} {val_disp}"
        print(line)


def _print_params(param: dict, no_emoji: bool) -> None:
    print(WHITE + PROMPT + RESET)
    print(BOLD + GREEN + f"{PROMPT}{'' if no_emoji else ' 🎬'} CBICALL PARAMETERS:" + RESET)
    print(WHITE + PROMPT + RESET)

    max_key = max(len(k) for k in param.keys())
    for key in sorted(param.keys()):
        val = param.get(key)
        val_disp = str(val) if val is not None and str(val) != "" else "(undef)"
        line = f"{PROMPT:>5} {key:<{max_key}} {ARROW} {val_disp}"
        print(line)


def main() -> int:
    start_time = time.time()
    cbicall_path = Path(sys.argv[0]).resolve()

    # Parse CLI args (Help::usage)
    arg = usage(VERSION)
    no_emoji = bool(arg.get("noemoji") or 0)
    no_spinner = bool(arg.get("debug") or arg.get("verbose"))

    # Read parameters and build config (Config::read_param_file + set_config_values)
    param = config_mod.read_param_file(arg["paramfile"])
    cfg = config_mod.set_config_values(param)
    cfg["version"] = VERSION

    # Start CBICall banners
    print(
        BOLD + CYAN + f"{PROMPT}{'' if no_emoji else ' 🚀'} CBICall {VERSION}" + RESET
    )
    print(f"{PROMPT}{'' if no_emoji else ' 🖥️ '} CBICall exe: {cbicall_path}")
    print(f"{PROMPT}{'' if no_emoji else ' ✍️ '} {AUTHOR}")
    print(f"{PROMPT}{'' if no_emoji else ' 📜'} {LICENSE}\n{PROMPT}")

    # Print arguments used
    print(BOLD + YELLOW + f"{PROMPT}{'' if no_emoji else ' 🔧'} ARGUMENTS USED:" + RESET)
    if arg.get("paramfile"):
        print(f"{WHITE}{PROMPT} --p {arg['paramfile']}{RESET}")
    if arg.get("threads") is not None:
        print(f"{WHITE}{PROMPT} --t {arg['threads']}{RESET}")

    print(WHITE + PROMPT + RESET)

    # Print config and params
    _print_config(cfg, no_emoji)
    _print_params(param, no_emoji)

    # Start CBICall
    print(PROMPT)
    print(BOLD + CYAN + f"{PROMPT}{'' if no_emoji else ' 🚦'} STARTING CBICALL FUN" + RESET)
    print(f"{PROMPT} {SPACER}")

    # Create working dir and log
    Path(cfg["projectdir"]).mkdir(parents=True, exist_ok=True)
    write_log(cfg, arg, param)

    # Build settings hash (rah_cbicall)
    settings = {
        "projectdir": cfg["projectdir"],
        "pipeline": param["pipeline"],
        "mode": param["mode"],
        "sample": param.get("sample"),
        "sample_map": param.get("sample_map"),
        "workflow_engine": param["workflow_engine"],
        "gatk_version": param["gatk_version"],
        "cleanup_bam": param.get("cleanup_bam", False),
        "threads": arg["threads"],
        "id": cfg["id"],
        "debug": arg["debug"],
        "genome": cfg.get("genome", "b37"),
        "bash_mit_cohort": cfg.get("bash_mit_cohort"),
        "bash_mit_single": cfg.get("bash_mit_single"),
        "bash_wes_cohort": cfg.get("bash_wes_cohort"),
        "bash_wes_single": cfg.get("bash_wes_single"),
        "bash_wgs_single": cfg.get("bash_wgs_single"),
        "smk_mit_cohort": cfg.get("smk_mit_cohort"),
        "smk_mit_single": cfg.get("smk_mit_single"),
        "smk_wes_cohort": cfg.get("smk_wes_cohort"),
        "smk_wes_single": cfg.get("smk_wes_single"),
        "smk_wgs_single": cfg.get("smk_wgs_single"),
    }

    print(
        f"{PROMPT} Running "
        f"{param['workflow_engine']}->{param['pipeline']}->{param['mode']} "
        f"workflow {'' if no_emoji else '🧬 '}..."
    )

    wes = DNAseq(settings)

    # Run with spinner
    run_with_spinner(
        wes.variant_calling,
        no_spinner=no_spinner,
        no_emoji=no_emoji,
    )

    # END CBICALL
    print(f"{PROMPT} {SPACER}")
    print(
        BOLD + GREEN + f"{PROMPT}{'' if no_emoji else ' ✅'} CBICALL FINISHED OK" + RESET
    )

    if arg.get("verbose"):
        print(f"{PROMPT} Date: {time.ctime()}")
        print(f"{PROMPT} Running time(s): {time.time() - start_time:.2f}")

    gb = GoodBye()
    print(f"{WHITE}{PROMPT}{'' if no_emoji else ' 👋 '} {gb.say_goodbye()}{RESET}")

    return 0

