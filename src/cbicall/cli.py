import argparse
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import List

from . import __version__
from . import config as config_mod
from . import console
from .cli_output import (
    _format_duration,
    _print_run_summary,
    _short_path,
    _warn,
)
from .comparison_commands import run_compare_runs_command
from .demo import run_demo
from .doctor import run_doctor
from .execution import WorkflowExecutor
from .errors import ParameterValidationError
from .goodbye import GoodBye
from .helpmod import handle_main_args, parse_run_args as _parse_run_args
from .html_reports import write_run_report_html
from .integration_tests import (
    prepare_integration_root,
    run_integration_tests,
    run_release_equivalence_test,
    selected_tests_from_args,
)
from .models import ResolvedConfig, RunSettings
from .multiqc import write_multiqc_report
from .paths import runtime_root
from .report_commands import run_report_command
from .report_utils import _nested
from .resource_install import run_resource_installer
from .resources import validate_resource_catalog
from .run_audit import write_log, write_run_report
from .workflow_registry import load_workflow_registry

VERSION = __version__

_refresh_colors = console.refresh_colors
_section = console.section
_row = console.row
_print_config = console.print_config
_print_params = console.print_params


def parse_args(argv):
    """Parse the explicit ``cbicall run`` command for programmatic callers."""
    if not argv or argv[0] != "run":
        raise SystemExit("The run subcommand is required")
    return _parse_run_args(argv[1:], VERSION)


def _project_root() -> Path:
    return runtime_root(__file__)


def _print_external_output_pointers(report: dict) -> None:
    summary = _nested(report, "outputs", "external_summary") or {}
    if not summary:
        return

    print()
    _section("nf-core outputs", console.BLUE)
    workflow_output = _nested(report, "outputs", "workflow_output_dir")
    if workflow_output:
        _row("Workflow out", _short_path(workflow_output))

    pipeline_info = summary.get("pipeline_info") or {}
    if pipeline_info.get("dir"):
        _row("Pipeline", _short_path(pipeline_info["dir"]))
    if pipeline_info.get("trace"):
        _row("Trace", _short_path(pipeline_info["trace"]))

    multiqc = summary.get("multiqc") or {}
    if multiqc.get("report"):
        _row("MultiQC", _short_path(multiqc["report"]))

    canonical_outputs = summary.get("canonical_outputs") or []
    shown = 0
    for item in canonical_outputs:
        matches = item.get("matches") or []
        if not matches:
            continue
        label = "Canonical" if shown == 0 else f"Canonical {shown + 1}"
        _row(label, _short_path(matches[0]))
        shown += 1


def _require_verified_resource(resolved_config: ResolvedConfig) -> None:
    bundle = (resolved_config.resources or {}).get("bundle") or {}
    if not bundle:
        return
    resource_type = bundle.get("type")
    if resource_type not in {None, "bundle"}:
        return

    runtime_check = bundle.get("runtime_check") or {}
    status = runtime_check.get("status")
    if status == "verified":
        return

    details = [
        f"Selected bundle resource could not be verified (status={status or 'unknown'}).",
        f"  resource: {bundle.get('key') or '(undef)'}",
    ]
    if runtime_check.get("datadir"):
        details.append(f"  DATADIR: {runtime_check.get('datadir')}")
    if runtime_check.get("source"):
        details.append(f"  source: {runtime_check.get('source')}")

    configured_data = os.environ.get("CBICALL_DATA")
    datadir = runtime_check.get("datadir")
    if status == "datadir_missing" and not configured_data and datadir == "/cbicall-data":
        details.extend(
            [
                "CBICALL_DATA is not set; CBIcall used /cbicall-data, the container default.",
                "For a local or source installation, set the resource root containing Databases/ and NGSutils/:",
                "  export CBICALL_DATA=/absolute/path/to/cbicall-data",
                "Then validate it:",
                "  cbicall validate-resources",
            ]
        )
    elif status == "datadir_missing" and configured_data:
        details.extend(
            [
                "CBICALL_DATA is set, but the resolved directory does not exist.",
                "Set it to the resource root containing Databases/ and NGSutils/:",
                "  export CBICALL_DATA=/absolute/path/to/cbicall-data",
                "Then validate it:",
                "  cbicall validate-resources",
            ]
        )
    else:
        details.append("Run cbicall validate-resources and check the configured resource directory before launching.")
    raise ParameterValidationError("\n".join(details))

def _apply_cli_runtime_overrides(params: dict, arg: dict) -> dict:
    params = dict(params)
    if arg.get("profile") is not None:
        params["profile"] = arg["profile"]
    else:
        params["profile"] = "local"
    return params


def _run_validate_parameters_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall validate-parameters",
        description="Validate a CBIcall parameters YAML without starting the workflow.",
    )
    parser.add_argument("-p", "--param", dest="paramfile", required=True, help="Parameters YAML file.")
    parser.add_argument("--runtime-profile", dest="profile", default="local", help="CBIcall runtime profile for native workflows.")
    parser.add_argument("--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    params = config_mod.read_param_file(args.paramfile)
    params = _apply_cli_runtime_overrides(params, vars(args))
    resolved_config = ResolvedConfig.from_mapping({**config_mod.set_config_values(params), "version": VERSION})
    _require_verified_resource(resolved_config)

    _section("Parameters OK", console.GREEN)
    workflow = resolved_config.workflow
    bundle = resolved_config.resources.get("bundle", {})
    _row("Param file", _short_path(args.paramfile))
    _row("Runtime profile", resolved_config.profile)
    _row("Workflow", f"{workflow.backend} -> {workflow.pipeline} -> {workflow.mode}")
    _row("Workflow provider", resolved_config.workflow_provider)
    _row("Software stack", workflow.software_stack)
    _row("Registry ver", workflow.registry_version)
    _row("Genome", resolved_config.genome or "b37")
    if workflow.metadata.get("provider") == "nf-core":
        _row("External workflow", workflow.metadata.get("source"))
        _row("External release", workflow.metadata.get("release"))
        _row("NF profile", resolved_config.nfcore_profile)
        _row("NF parameters", ", ".join(sorted(resolved_config.nfcore_parameters)) or "(none)")
        if resolved_config.nfcore_singularity_cache_dir:
            _row("NF cache", _short_path(resolved_config.nfcore_singularity_cache_dir))
    else:
        _row("Entrypoint", _short_path(workflow.entrypoint))
    if workflow.backend == "bash":
        _row("Env file", _short_path(workflow.helpers.get("env")))
    elif workflow.backend in {"snakemake", "nextflow"} and not workflow.metadata.get("provider"):
        _row("Config", _short_path(workflow.config_file))
    _row("Resource key", bundle.get("key"))
    _row("Resource ver", bundle.get("version"))
    _row("Resource hash", bundle.get("fingerprint"))
    runtime_check = bundle.get("runtime_check", {})
    _row("DATADIR", _short_path(runtime_check.get("datadir")))
    _row("Resource check", runtime_check.get("status"))
    return 0


def _collect_external_sources(value) -> List[str]:
    sources = set()

    def visit(node):
        if isinstance(node, dict):
            provider = node.get("provider")
            if provider:
                sources.add(str(provider))
            for child in node.values():
                visit(child)
        elif isinstance(node, list):
            for child in node:
                visit(child)

    visit(value)
    return sorted(sources)


def _collect_cromwell_wdls(registry: dict, project_root: Path) -> List[Path]:
    workflows = registry.get("workflows") or {}
    cromwell = workflows.get("cromwell") or {}
    base_root = project_root / str(cromwell.get("base_dir", "workflows/cromwell"))
    wdls = []
    for stack_name, stack in (cromwell.get("software_stacks") or {}).items():
        stack_dir = base_root / str(stack_name)
        for pipeline in (stack.get("pipelines") or {}).values():
            for mode in pipeline.values():
                for implementation in (mode.get("registry_versions") or {}).values():
                    script = implementation.get("script") if isinstance(implementation, dict) else implementation
                    if script:
                        wdls.append((stack_dir / str(script)).resolve())
    return sorted(set(wdls))


def _validate_cromwell_wdls_with_womtool(registry: dict, project_root: Path) -> dict:
    wdls = _collect_cromwell_wdls(registry, project_root)
    if not wdls:
        return {"status": "not_applicable", "count": 0, "files": []}
    womtool = os.environ.get("WOMTOOL_JAR")
    if not womtool:
        return {"status": "skipped", "reason": "WOMTOOL_JAR not set", "count": len(wdls), "files": [str(p) for p in wdls]}
    womtool_path = Path(womtool)
    if not womtool_path.is_file():
        raise ParameterValidationError(f"WOMTOOL_JAR does not exist: {womtool}")
    java_cmd = os.environ.get("JAVA_CMD", "java")
    for wdl in wdls:
        proc = subprocess.run(
            [java_cmd, "-jar", str(womtool_path), "validate", str(wdl)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or "womtool validation failed").strip()
            raise ParameterValidationError(f"WDL validation failed for {wdl}: {detail}")
    return {"status": "ok", "count": len(wdls), "files": [str(p) for p in wdls]}


def _run_validate_registry_command(argv: List[str]) -> int:
    root = _project_root()
    default_registry = root / "workflows" / "registry" / "cbicall-workflow-registry.yaml"
    default_schema = root / "workflows" / "schema" / "cbicall-workflow-registry.schema.json"

    parser = argparse.ArgumentParser(
        prog="cbicall validate-registry",
        description="Validate the workflow registry YAML against its JSON Schema.",
    )
    parser.add_argument("--registry", default=str(default_registry), help="Workflow registry YAML.")
    parser.add_argument("--schema", default=str(default_schema), help="Workflow registry JSON Schema.")
    parser.add_argument("--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    registry_path = Path(args.registry)
    schema_path = Path(args.schema)
    registry = load_workflow_registry(registry_path, schema_path)
    backends = sorted((registry.get("workflows") or {}).keys())
    external_sources = _collect_external_sources(registry)
    wdl_validation = _validate_cromwell_wdls_with_womtool(registry, root)

    _section("Registry OK", console.GREEN)
    _row("Registry", _short_path(registry_path))
    _row("Schema", _short_path(schema_path))
    _row("Backends", ", ".join(backends) if backends else "(none)")
    _row("External", ", ".join(external_sources) if external_sources else "(none)")
    if wdl_validation["status"] == "ok":
        _row("WDL syntax", f"ok ({wdl_validation['count']} files)")
    elif wdl_validation["status"] == "skipped":
        _row("WDL syntax", f"skipped ({wdl_validation['reason']})")
    return 0


def _run_validate_resources_command(argv: List[str]) -> int:
    root = _project_root()
    default_catalog = root / "resources" / "cbicall-resource-catalog.json"
    default_registry = root / "workflows" / "registry" / "cbicall-workflow-registry.yaml"
    default_schema = root / "workflows" / "schema" / "cbicall-workflow-registry.schema.json"

    parser = argparse.ArgumentParser(
        prog="cbicall validate-resources",
        description="Validate the resource catalog and workflow compatibility keys.",
    )
    parser.add_argument("--catalog", default=str(default_catalog), help="Resource catalog JSON.")
    parser.add_argument("-r", "--resource", help="Validate one resource entry by resource key.")
    parser.add_argument("--registry", default=str(default_registry), help="Workflow registry YAML.")
    parser.add_argument("--schema", default=str(default_schema), help="Workflow registry JSON Schema.")
    parser.add_argument("--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    registry = load_workflow_registry(Path(args.registry), Path(args.schema))
    summary = validate_resource_catalog(Path(args.catalog), registry, resource_key=args.resource)

    _section("Resources OK", console.GREEN)
    _row("Catalog", _short_path(summary["path"]))
    _row("Schema", _short_path(summary["schema"]))
    _row("Schema version", summary["schema_version"])
    if summary.get("resource_key"):
        _row("Resource key", summary["resource_key"])
    _row("Resources", summary["resources"])
    _row("Bundle resources", summary["bundle_resources"])
    _row("Compatible workflows", summary["compatible_workflows"])
    return 0


def _run_demo_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall demo",
        description="Generate resource-free WES and mtDNA reports from packaged example outputs.",
    )
    parser.add_argument(
        "--output-dir",
        default="cbicall-demo",
        help="New demo output directory (default: cbicall-demo).",
    )
    args = parser.parse_args(argv)

    _section("CBIcall Demo", console.GREEN)
    print("Generating reports from precomputed CNAG99901P fixture outputs.")
    print("No external resource bundle or workflow backend is required.")
    result = run_demo(Path(args.output_dir))
    _row("Output", result.output_dir)
    _row("WES HTML", result.wes_html)
    _row("mtDNA HTML", result.mtdna_html)
    _row("mtDNA variants", result.mtdna_variants)
    _row("Read me", result.output_dir / "README.txt")
    return 0


def _run_test_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall test",
        description="Run the bundled CBIcall integration tests from examples/input.",
    )
    parser.add_argument("--wes-bash", action="store_true", help="Run the Bash WES integration test.")
    parser.add_argument(
        "--wes-cohort-bash",
        action="store_true",
        help="Run the Bash WES cohort joint-genotyping integration test.",
    )
    parser.add_argument(
        "--wes-cohort-bash-sharded",
        action="store_true",
        help="Run the Bash WES cohort staged shard integration test.",
    )
    parser.add_argument("--wes-bash-gatk35", action="store_true", help="Run the legacy Bash WES GATK 3.5 integration test.")
    parser.add_argument(
        "--wes-snakemake",
        action="store_true",
        help="Run the Snakemake WES integration test. Requires snakemake on PATH.",
    )
    parser.add_argument(
        "--wes-nextflow",
        action="store_true",
        help="Run the Nextflow WES integration test. Requires nextflow on PATH.",
    )
    parser.add_argument(
        "--wes-cromwell",
        action="store_true",
        help="Run the Cromwell WES integration test. Requires CROMWELL_JAR or cromwell on PATH.",
    )
    parser.add_argument("--mit-bash", action="store_true", help="Run the Bash mitochondrial integration test.")
    parser.add_argument(
        "--nf-core-demo",
        action="store_true",
        help="Run the external nf-core/demo integration contract. Requires nextflow on PATH.",
    )
    parser.add_argument(
        "--nf-core-sarek",
        action="store_true",
        help="Run the external nf-core/Sarek integration contract. Requires nextflow on PATH.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all native bundled integration tests. Optional backend tests are skipped if their backend executable is not installed.",
    )
    parser.add_argument(
        "--backend-equivalence",
        dest="backend_equivalence",
        action="store_true",
        help="Run native WES backend-equivalence checks and compare normalized VCF output against Bash.",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use.")
    parser.add_argument("--runtime-profile", dest="profile", default="local", help="CBIcall runtime profile for native workflow tests.")
    parser.add_argument(
        "--workspace",
        type=Path,
        help="Stage and run tests in this new or empty directory instead of the default location.",
    )
    parser.add_argument(
        "--keep-external-work",
        action="store_true",
        help="Keep heavy Nextflow work directories for external nf-core tests.",
    )
    args = parser.parse_args(argv)

    if args.threads <= 0:
        parser.error("--threads requires a positive integer")

    selectors = [
        args.wes_bash,
        args.wes_cohort_bash,
        args.wes_cohort_bash_sharded,
        args.wes_bash_gatk35,
        args.wes_snakemake,
        args.wes_nextflow,
        args.wes_cromwell,
        args.mit_bash,
        args.nf_core_demo,
        args.nf_core_sarek,
        args.all,
    ]
    if args.backend_equivalence and any(selectors):
        parser.error("--backend-equivalence cannot be combined with individual test selectors or --all")

    selected = selected_tests_from_args(args)
    if not args.backend_equivalence and not selected:
        parser.error(
            "select at least one test with --wes-bash, --wes-cohort-bash, --wes-cohort-bash-sharded, --wes-bash-gatk35, "
            "--wes-snakemake, --wes-nextflow, --wes-cromwell, --mit-bash, --nf-core-demo, "
            "--nf-core-sarek, --backend-equivalence, or --all"
        )

    root, staged_root = prepare_integration_root(_project_root(), workspace=args.workspace)
    if staged_root is not None:
        console.row("Test workspace", staged_root)
    if args.backend_equivalence:
        return run_release_equivalence_test(
            project_root=root,
            threads=args.threads,
            runtime_profile=args.profile,
        )

    return run_integration_tests(
        project_root=root,
        selected=selected,
        threads=args.threads,
        runtime_profile=args.profile,
        skip_missing_optional=args.all,
        keep_external_work=args.keep_external_work,
    )


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
            spinner_char = f"{console.BOLD}{console.YELLOW}{frames[i % len(frames)]}{console.RESET}"
            elapsed = time.time() - start
            elapsed_str = _format_duration(elapsed)
            msg = f"{console.BOLD}{console.WHITE} Working (elapsed: {elapsed_str}){console.RESET}"
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


def _run_analysis(arg: dict, *, start_time: float, cbicall_path: Path) -> int:
    if arg.get("nocolor"):
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()
    no_spinner = bool(arg.get("debug") or arg.get("verbose"))

    # Read parameters and build config (Config::read_param_file + set_config_values)
    params = config_mod.read_param_file(arg["paramfile"])
    params = _apply_cli_runtime_overrides(params, arg)
    resolved_config = ResolvedConfig.from_mapping({**config_mod.set_config_values(params), "version": VERSION})
    _require_verified_resource(resolved_config)

    if params.get("genome") is None and resolved_config.genome is not None:
        _warn(f"Genome not provided; using inferred default '{resolved_config.genome}'.", console.BOLD, console.YELLOW, console.RESET)

    _print_run_summary(
        arg=arg,
        resolved_config=resolved_config,
        cbicall_path=cbicall_path,
        version=VERSION,
        colors={"bold": console.BOLD, "reset": console.RESET, "blue": console.BLUE, "cyan": console.CYAN, "green": console.GREEN, "yellow": console.YELLOW},
    )

    if arg.get("verbose") or arg.get("debug"):
        _print_config(resolved_config)
        print()
        _print_params(params)
        print()

    # Start CBIcall
    _section("Running", console.CYAN)

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
            "snakemake_parameters": resolved_config.snakemake_parameters,
            "nextflow_parameters": resolved_config.nextflow_parameters,
            "cromwell_parameters": resolved_config.cromwell_parameters,
            "nfcore_profile": resolved_config.nfcore_profile,
            "nfcore_parameters": resolved_config.nfcore_parameters,
            "nfcore_singularity_cache_dir": resolved_config.nfcore_singularity_cache_dir,
            "genome": resolved_config.genome,
            "qc_coverage_region": resolved_config.qc_coverage_region,
            "cleanup_bam": params.get("cleanup_bam", False),
            "output_basename": resolved_config.output_basename,
            "cohort_stage": resolved_config.cohort_stage,
            "interval_shard": resolved_config.interval_shard,
            "run_mode": resolved_config.run_mode,
            "inputs": resolved_config.inputs.to_dict(),
            "workflow": resolved_config.workflow.to_dict(),
        }
    )

    workflow = resolved_config.workflow
    genome = resolved_config.display_genome or resolved_config.genome or "b37"
    if workflow.metadata.get("provider") == "nf-core":
        log_name = f"nf-core_{workflow.pipeline}_{workflow.mode}.log"
    else:
        log_name = f"{workflow.backend}_{workflow.software_stack}_{workflow.pipeline}_{workflow.mode}_{genome}.log"
    workflow_log = Path(resolved_config.project_dir) / log_name
    _row("Workflow", f"{workflow.backend} -> {workflow.pipeline} -> {workflow.mode}")
    print("  This workflow may take a while depending on input size and pipeline.")

    wes = WorkflowExecutor(settings)

    # Run with spinner
    try:
        run_with_spinner(
            wes.run,
            no_spinner=no_spinner,
        )
    except Exception as exc:
        elapsed = time.time() - start_time
        try:
            report_path = write_run_report(
                resolved_config,
                arg,
                params,
                elapsed_seconds=elapsed,
                workflow_log=workflow_log,
                status="failed",
                error={"type": type(exc).__name__, "message": str(exc)},
            )
            html_report_path = write_run_report_html(
                report_path,
                json.loads(report_path.read_text(encoding="utf-8")),
            )
            print()
            _section("Failed", console.RED)
            _row("Status", "Workflow execution failed")
            _row("Elapsed", _format_duration(elapsed))
            _row("Log", workflow_log)
            _row("Report", report_path)
            _row("HTML", html_report_path)
        except Exception as report_exc:
            print(
                f"WARNING: Could not write failed-run audit report: {report_exc}",
                file=sys.stderr,
            )
        raise

    # END CBICALL
    print()
    _section("Completed", console.GREEN)
    _row("Status", "Finished successfully")
    elapsed = time.time() - start_time
    _row("Elapsed", _format_duration(elapsed))
    report_path = write_run_report(
        resolved_config,
        arg,
        params,
        elapsed_seconds=elapsed,
        workflow_log=workflow_log,
    )
    report_payload = json.loads(report_path.read_text(encoding="utf-8"))
    html_report_path = write_run_report_html(report_path, report_payload)
    multiqc_path = write_multiqc_report(report_path, report_payload) if arg.get("multiqc") else None
    _row("Log", workflow_log)
    _row("Report", report_path)
    _row("HTML", html_report_path)
    if multiqc_path is not None:
        _row("MultiQC", multiqc_path)
    if workflow.metadata.get("provider") == "nf-core":
        _print_external_output_pointers(report_payload)
    if arg.get("verbose"):
        _row("Date", time.ctime())

    gb = GoodBye()
    print(f"{console.WHITE}{gb.say_goodbye()}{console.RESET}")

    return 0


def _run_run_command(argv: List[str], *, start_time: float, cbicall_path: Path) -> int:
    arg = vars(_parse_run_args(argv, VERSION))
    return _run_analysis(arg, start_time=start_time, cbicall_path=cbicall_path)


def main() -> int:
    start_time = time.time()
    cbicall_path = Path(sys.argv[0]).resolve()

    if len(sys.argv) > 1 and sys.argv[1] == "run":
        return _run_run_command(sys.argv[2:], start_time=start_time, cbicall_path=cbicall_path)
    if len(sys.argv) > 1 and sys.argv[1] == "validate-parameters":
        return _run_validate_parameters_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "validate-registry":
        return _run_validate_registry_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "validate-resources":
        return _run_validate_resources_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "install-resources":
        return run_resource_installer(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "doctor":
        return run_doctor(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "compare-runs":
        return run_compare_runs_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        return run_report_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "demo":
        return _run_demo_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        return _run_test_command(sys.argv[2:])

    return handle_main_args(sys.argv[1:], VERSION)
