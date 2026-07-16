import argparse
import os
import sys
import textwrap


def _build_main_parser(version: str) -> argparse.ArgumentParser:
    return argparse.ArgumentParser(
        prog="cbicall",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            f"""\
            CBIcall {version}
            CNAG Biomedical Informatics Framework for variant calling on Illumina DNA-seq.

            Commands:
              run                  Execute an analysis from a parameters YAML file
              validate-parameters  Validate a parameters YAML file
              validate-registry    Validate the workflow registry
              validate-resources   Validate the resource catalog and installation
              install-resources    Install a registered resource bundle
              compare-runs         Compare two or more completed executions
              report               Regenerate reports for a completed execution
              demo                 Generate reports from packaged example data
              test                 Run integration and backend-equivalence tests

            Run 'cbicall COMMAND --help' for command-specific options.
            """
        ),
        add_help=True,
    )


def _build_run_parser(version: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cbicall run",
        description="Execute one CBIcall analysis from a parameters YAML file.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=version,
        help="Show version information and exit",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print execution commands and preserve tracebacks",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show resolved configuration and input parameters",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="Number of CPUs/cores/threads (required)",
    )
    parser.add_argument(
        "-p",
        "--param",
        dest="paramfile",
        metavar="FILE",
        help="Parameters input file (YAML; required)",
    )
    parser.add_argument(
        "--runtime-profile",
        dest="profile",
        default="local",
        help="CBIcall runtime profile for native workflows (default: local)",
    )
    parser.add_argument(
        "--no-color",
        dest="nocolor",
        action="store_true",
        help="Do not print colors to standard output",
    )
    parser.add_argument(
        "--multiqc",
        action="store_true",
        help="Write cbicall_mqc/ MultiQC custom-content files after a successful run",
    )
    return parser


def _validate_run_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.threads is None or args.paramfile is None:
        print("Options --threads and --param are required", file=sys.stderr)
        raise SystemExit(1)

    if args.threads <= 0:
        print("Option --threads requires a positive integer", file=sys.stderr)
        raise SystemExit(1)

    if not os.path.isfile(args.paramfile) or os.path.getsize(args.paramfile) == 0:
        print(
            "Option --param requires an existing, non-empty parameters file",
            file=sys.stderr,
        )
        raise SystemExit(1)

    return args


def parse_run_args(argv, version: str) -> argparse.Namespace:
    """Parse and validate arguments for ``cbicall run``."""
    parser = _build_run_parser(version)
    removed_options = {"-man", "-debug", "-verbose", "-nc"}
    invalid = next((value for value in argv if value in removed_options), None)
    if invalid is not None:
        parser.error(f"unrecognized argument: {invalid}")
    return _validate_run_args(parser.parse_args(argv))


def handle_main_args(argv, version: str) -> int:
    """Print top-level help or reject arguments that are not subcommands."""
    parser = _build_main_parser(version)
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=version,
        help="Show version information and exit",
    )
    if argv:
        parser.parse_args(argv)
    parser.print_help()
    return 0
