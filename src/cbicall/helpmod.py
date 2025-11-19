import os
import sys
import textwrap
import argparse


def _build_parser(version: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        add_help=False,
        prog="cbicall",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            f"""\
            CBICall {version}
            CNAG Biomedical Informatics Framework for variant calling on Illumina DNA-seq.
            """
        ),
    )

    parser.add_argument("-v", action="store_true", dest="show_version",
                        help="Show version information and exit")
    parser.add_argument("-debug", type=int, dest="debug",
                        help="Debugging level (1–5)")
    parser.add_argument("-verbose", action="store_true", dest="verbose",
                        help="Enable verbose output")
    parser.add_argument("-h", "--help", dest="help", action="store_true",
                        help="Brief help message")
    parser.add_argument("-man", dest="man", action="store_true",
                        help="Full documentation (same as help for now)")
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        dest="threads",
        help="Number of CPUs/Cores/Threads (required)",
    )
    parser.add_argument(
        "-p",
        "--param",
        dest="paramfile",
        metavar="FILE",
        help="Parameters input file (YAML; required)",
    )
    parser.add_argument(
        "-nc",
        "--no-color",
        dest="nocolor",
        action="store_true",
        help="Do not print colors to STDOUT",
    )
    parser.add_argument(
        "-ne",
        "--no-emoji",
        dest="noemoji",
        action="store_true",
        help="Do not print emojis to STDOUT",
    )

    return parser


def _validate_args_core(args: argparse.Namespace) -> argparse.Namespace:
    if args.threads is None or args.paramfile is None:
        print("Options --threads and --param are required", file=sys.stderr)
        raise SystemExit(1)

    if args.threads is not None and args.threads <= 0:
        print("Option --threads requires a positive integer", file=sys.stderr)
        raise SystemExit(1)

    if args.paramfile is not None:
        if (not os.path.isfile(args.paramfile)
                or os.path.getsize(args.paramfile) == 0):
            print(
                "Option --param requires an existing, non-empty parameters file",
                file=sys.stderr,
            )
            raise SystemExit(1)

    return args

def parse_args(argv, version: str) -> argparse.Namespace:
    """
    Pure argument parsing + validation.
    Does not call sys.exit; raises ValueError on invalid args.
    Suitable for tests.
    """
    parser = _build_parser(version)
    args = parser.parse_args(argv)
    return _validate_args_core(args)



def usage(version: str) -> dict:
    """
    CLI entry-point style argument handling.
    Uses sys.argv, prints help/usage, and exits on error.
    Returns a dict of parsed args on success.
    """
    parser = _build_parser(version)
    args = parser.parse_args()

    if args.show_version:
        print(version)
        sys.exit(0)

    if args.help or args.man:
        parser.print_help()
        sys.exit(0)

    try:
        _validate_args_core(args)
    except ValueError as e:
        # parser.error() prints message and exits with status 2
        parser.error(str(e))

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"

    return vars(args)
