"""Shared terminal rendering and ANSI color state for CBIcall commands."""

import os
import sys

from .cli_output import (
    ARROW,
    _print_config as _render_config,
    _print_params as _render_params,
    _row as _render_row,
    _section as _render_section,
    _warn as _render_warning,
)
from .models import ResolvedConfig


def colors_enabled() -> bool:
    if os.environ.get("ANSI_COLORS_DISABLED") or os.environ.get("NO_COLOR"):
        return False
    return sys.stdout.isatty()


def _code(value: str) -> str:
    return value if colors_enabled() else ""


def refresh_colors() -> None:
    global BOLD, RESET, RED, YELLOW, GREEN, BLUE, CYAN, WHITE
    BOLD = _code("\033[1m")
    RESET = _code("\033[0m")
    RED = _code("\033[31m")
    YELLOW = _code("\033[33m")
    GREEN = _code("\033[32m")
    BLUE = _code("\033[34m")
    CYAN = _code("\033[36m")
    WHITE = _code("\033[37m")


def section(title: str, color: str = "") -> None:
    _render_section(title, color or WHITE, BOLD, RESET)


def row(label: str, value) -> None:
    _render_row(label, value)


def warning(message: str) -> None:
    _render_warning(message, BOLD, YELLOW, RESET)


def print_config(resolved_config) -> None:
    if isinstance(resolved_config, ResolvedConfig):
        _render_config(resolved_config, BOLD, BLUE, RESET)
        return

    print(f"{BOLD}{BLUE}Resolved Configuration{RESET}")
    keys = list(resolved_config.keys())
    if not keys:
        return
    max_key = max(len(key) for key in keys)
    for key in sorted(keys):
        value = resolved_config.get(key)
        print(f"  {key:<{max_key}} {ARROW} {value if value is not None else '(undef)'}")


def print_params(params: dict) -> None:
    _render_params(params, BOLD, GREEN, RESET)


refresh_colors()
