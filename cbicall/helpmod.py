import os
import sys
import textwrap
import argparse


def _build_parser(version: str) -> argparse.ArgumentParser:
    # We implement manual -h/--help/-man to mimic Perl behaviour.
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


def usage(version: str) -> dict:
    """
    Parse CLI args and apply basic validation.
    Mirrors Help::usage (but using argparse instead of Pod::Usage).
    """
    parser = _build_parser(version)
    args = parser.parse_args()

    if args.show_version:
        print(version)
        sys.exit(0)

    if args.help or args.man:
        parser.print_help()
        sys.exit(0)

    # Control checks
    if args.threads is None or args.paramfile is None:
        parser.print_usage()
        sys.exit(1)

    if not (args.threads and args.threads > 0):
        parser.error("Option --threads requires a positive integer")

    if not os.path.isfile(args.paramfile) or os.path.getsize(args.paramfile) == 0:
        parser.error("Option --param requires an existing, non-empty parameters file")

    # Initialize undefs
    if args.debug is None:
        args.debug = 0

    # Global setting: disable ANSI colors
    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"

    return vars(args)
