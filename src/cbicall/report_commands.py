"""Commands for inspecting and regenerating completed run reports."""

import argparse
import json
import os
from pathlib import Path
from typing import List, Optional

from . import console
from .cli_output import _format_duration, _short_path
from .html_reports import write_run_report_html
from .multiqc import write_multiqc_report
from .run_audit import _collect_execution_contract, _collect_output_fingerprints


def _run_report_path(path: str) -> Path:
    report_path = Path(path)
    if report_path.is_dir():
        report_path = report_path / "run-report.json"
    if not report_path.is_file():
        raise FileNotFoundError(f"run-report.json not found: {report_path}")
    return report_path


def _refresh_report_output_fields(report_path: Path, payload: dict) -> bool:
    """Refresh audit fields that can be derived from files in an existing run."""
    project_dir = report_path.parent
    if not project_dir.is_dir():
        return False

    fresh = _collect_output_fingerprints(project_dir)
    outputs = payload.setdefault("outputs", {})
    changed = False

    inventory = outputs.setdefault("file_inventory", {})
    fresh_inventory = fresh.get("file_inventory") or {}
    if int(fresh_inventory.get("entries") or 0) > 0:
        keys = [
            "directories",
            "largest_files",
            "total_bytes",
            "entries",
            "sha256",
            "algorithm",
            "scope",
            "excluded",
            "paths",
        ]
        for key in keys:
            if fresh_inventory.get(key) is not None and inventory.get(key) != fresh_inventory.get(key):
                inventory[key] = fresh_inventory[key]
                changed = True

    fresh_vcfs = fresh.get("vcf_hash_reports") or []
    if fresh_vcfs and outputs.get("vcf_hash_reports") != fresh_vcfs:
        outputs["vcf_hash_reports"] = fresh_vcfs
        changed = True

    fresh_contract = _collect_execution_contract(project_dir)
    if fresh_contract and payload.get("execution_contract") != fresh_contract:
        payload["execution_contract"] = fresh_contract
        changed = True

    return changed


def _load_run_report(path: str) -> dict:
    report_path = _run_report_path(path)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid run report JSON: {report_path}") from exc
    _refresh_report_output_fields(report_path, report)
    report["_report_path"] = str(report_path)
    return report


def _report_html_status(report_path: Path, written_html_path: Optional[Path]) -> str:
    html_path = written_html_path or report_path.with_suffix(".html")
    status = "exists" if html_path.is_file() else "missing"
    if written_html_path is not None:
        status = "written" if html_path.is_file() else status
    return f"{html_path} ({status})"


def _report_json_status(
    report_path: Path,
    refresh_requested: bool,
    refreshed: bool,
    wrote_json: bool,
) -> str:
    if wrote_json:
        status = "written"
    elif refresh_requested:
        status = "no changes"
    else:
        status = "read-only"
    return f"{report_path} ({status})"


def _print_single_run_report(
    payload: dict,
    report_path: Path,
    html_path: Optional[Path],
    multiqc_path: Optional[Path],
    refresh_requested: bool,
    refreshed: bool,
    wrote_json: bool = False,
) -> None:
    workflow = payload.get("workflow") or {}
    runtime = payload.get("runtime") or {}
    backend = runtime.get("backend") or {}
    resources = payload.get("resources") or {}
    bundle = resources.get("bundle") or {}
    outputs = payload.get("outputs") or {}
    inventory = outputs.get("file_inventory") or {}
    vcf_reports = outputs.get("vcf_hash_reports") or []
    canonical_outputs = outputs.get("canonical_outputs") or []
    execution_trace = payload.get("execution_trace") or {}
    execution_contract = payload.get("execution_contract") or {}
    software_versions = payload.get("software_versions") or {}
    run = payload.get("run") or {}

    color = console.GREEN if payload.get("status") == "success" else console.YELLOW
    console.section("Run Report", color)
    console.row("Report", _report_json_status(report_path, refresh_requested, refreshed, wrote_json))
    console.row("Status", payload.get("status"))
    console.row("Run ID", run.get("run_id"))
    console.row("Hostname", run.get("hostname"))
    console.row("Elapsed", _format_duration(float(payload.get("elapsed_seconds") or 0)))
    console.row("Project", _short_path(run.get("project_dir")))
    console.row("HTML", _report_html_status(report_path, html_path))
    if multiqc_path is not None:
        console.row("MultiQC", multiqc_path)

    console.section("Workflow", console.BLUE)
    console.row("Key", workflow.get("key"))
    console.row("Backend", workflow.get("backend"))
    console.row("Software stack", workflow.get("software_stack"))
    console.row("Pipeline", workflow.get("pipeline"))
    console.row("Mode", workflow.get("mode"))
    console.row("Registry ver", workflow.get("registry_version"))
    metadata = workflow.get("metadata") or {}
    if metadata.get("provider") == "nf-core":
        console.row("External workflow", metadata.get("source"))
        console.row("External release", metadata.get("release"))
    console.row("Fingerprint", workflow.get("fingerprint"))
    console.row("Backend ver", backend.get("version"))

    if execution_contract:
        console.section("Execution Contract", console.BLUE)
        console.row("Contract", _short_path(execution_contract.get("path")))
        console.row("Fingerprint", execution_contract.get("fingerprint"))
        console.row("Command hash", execution_contract.get("normalized_command_sha256"))
        console.row("Generated files", len(execution_contract.get("generated_files") or []))

    console.section("Resources", console.CYAN)
    console.row("Resource key", bundle.get("key"))
    console.row("Resource ver", bundle.get("version"))
    console.row("Resource hash", bundle.get("fingerprint"))

    console.section("Outputs", console.WHITE)
    console.row("Inventory files", inventory.get("entries"))
    console.row("Inventory size", inventory.get("total_bytes_human") or inventory.get("total_bytes"))
    console.row("Inventory hash", inventory.get("sha256"))
    console.row("VCF hashes", len(vcf_reports))
    console.row("Canonical outputs", len(canonical_outputs))
    if execution_trace:
        console.row("Trace tasks", execution_trace.get("task_count"))
        console.row("Max RSS", execution_trace.get("max_peak_rss_human"))
    if software_versions:
        console.row("Software scope", software_versions.get("scope"))
        console.row("Software hash", software_versions.get("sha256"))


def run_report_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall report",
        description="Summarize an existing CBIcall run-report.json or run directory.",
    )
    parser.add_argument("run", help="Run directory or run-report.json file.")
    parser.add_argument("--json", action="store_true", help="Print run-report JSON instead of a text summary.")
    parser.add_argument("--refresh", action="store_true", help="Refresh output-derived metadata and write run-report.json.")
    parser.add_argument(
        "--html",
        nargs="?",
        const=True,
        metavar="HTML",
        help="Write an HTML run report. Defaults to run-report.html unless a path is provided.",
    )
    parser.add_argument(
        "--multiqc",
        nargs="?",
        const=True,
        metavar="MQC_DIR",
        help="Write a MultiQC custom-content directory. Defaults to cbicall_mqc/ unless a path is provided.",
    )
    parser.add_argument(
        "-O",
        "--overwrite",
        action="store_true",
        help="Overwrite files written by --refresh, --html, or --multiqc.",
    )
    parser.add_argument("--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor or args.json:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    console.refresh_colors()

    report_path = _run_report_path(args.run)
    try:
        payload = json.loads(report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid run report JSON: {report_path}") from exc

    refreshed = False
    wrote_json = False
    if args.refresh:
        refreshed = _refresh_report_output_fields(report_path, payload)
        if refreshed and not args.overwrite:
            raise FileExistsError(
                f"run-report.json would be updated: {report_path}. Use -O/--overwrite to replace it."
            )
        if refreshed:
            report_path.write_text(
                json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True),
                encoding="utf-8",
            )
            wrote_json = True

    html_path = None
    if args.html is not None:
        html_path = report_path.with_suffix(".html") if args.html is True else Path(args.html)
        if html_path.exists() and not args.overwrite:
            raise FileExistsError(
                f"HTML report already exists: {html_path}. Use -O/--overwrite to replace it."
            )
        html_path.parent.mkdir(parents=True, exist_ok=True)
        write_run_report_html(report_path, payload, html_path=html_path)

    multiqc_path = None
    if args.multiqc is not None:
        multiqc_path = report_path.parent / "cbicall_mqc" if args.multiqc is True else Path(args.multiqc)
        if multiqc_path.exists() and not args.overwrite:
            raise FileExistsError(
                f"MultiQC custom-content directory already exists: {multiqc_path}. "
                "Use -O/--overwrite to replace it."
            )
        write_multiqc_report(report_path, payload, output_path=multiqc_path)

    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    else:
        _print_single_run_report(
            payload,
            report_path,
            html_path,
            multiqc_path,
            args.refresh,
            refreshed,
            wrote_json,
        )
    return 0
