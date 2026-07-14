"""Installed and ``python -m cbicall`` command entry point."""

import sys

from .cli import main


def _line_buffer_standard_streams() -> None:
    """Keep redirected CBIcall status logs readable during long runs."""
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is None:
            continue
        try:
            reconfigure(line_buffering=True, write_through=True)
        except TypeError:
            reconfigure(line_buffering=True)


def console_main() -> int:
    """Run the CLI with concise production errors and debug tracebacks."""
    _line_buffer_standard_streams()
    try:
        return main()
    except KeyboardInterrupt:
        print("Interrupted by user", file=sys.stderr)
        return 130
    except Exception as exc:
        if "-debug" in sys.argv or "--debug" in sys.argv:
            raise
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(console_main())
