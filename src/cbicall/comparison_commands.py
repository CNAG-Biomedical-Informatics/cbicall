"""Terminal comparison commands for completed CBIcall runs."""

import argparse
import io
import os
from contextlib import redirect_stdout
from pathlib import Path
from typing import List

from . import console
from .cli_output import _short_path
from .html_reports import render_compare_html
from .multiqc import write_compare_multiqc_report
from .report_commands import _load_run_report
from .report_utils import (
    _aggregate_status,
    _comparison_sections_with_overall,
    _comparison_value,
    _execution_file_value,
    _execution_file_map,
    _inventory_manifest_hash,
    _inventory_total_bytes,
    _format_optional_bytes,
    _multi_vcf_call_value,
    _multi_vcf_hash_value,
    _multi_workflow_file_value,
    _nested,
    _row_pair_status,
    _run_label,
    _same_text,
    _vcf_hash_map,
    _workflow_file_map,
)

_section = console.section
_row = console.row
_refresh_colors = console.refresh_colors


def _compare_row(label: str, left, right) -> None:
    if left is None and right is None:
        _row(label, "not available")
        return
    if left is None or right is None:
        _row(label, f"missing: {_comparison_value(left)} != {_comparison_value(right)}")
        return
    status = _same_text(left, right)
    value = _comparison_value(left) if left == right else f"{_comparison_value(left)} != {_comparison_value(right)}"
    _row(label, f"{status}: {value}")


def _compare_execution_files(left: dict, right: dict) -> None:
    left_files = _execution_file_map(left)
    right_files = _execution_file_map(right)
    roles = sorted(set(left_files) | set(right_files))
    if not roles:
        _row("Generated files", "not available")
        return
    for role in roles:
        _compare_row(role, _execution_file_value(role, left), _execution_file_value(role, right))


def _compare_workflow_files(left: dict, right: dict) -> None:
    left_files = _workflow_file_map(left)
    right_files = _workflow_file_map(right)
    roles = sorted(set(left_files) | set(right_files))
    if not roles:
        _row("Files", "not available")
        return
    for role in roles:
        left_entry = left_files.get(role)
        right_entry = right_files.get(role)
        if not left_entry or not right_entry:
            _row(role, "missing")
            continue
        same_path = left_entry.get("path") == right_entry.get("path")
        same_sha = left_entry.get("sha256") == right_entry.get("sha256")
        status = "same" if same_path and same_sha else "different"
        details = []
        if not same_path:
            details.append("path")
        if not same_sha:
            details.append("sha256")
        if details:
            _row(role, f"{status}: {', '.join(details)}")
        else:
            _row(role, f"same: {_comparison_value((left_entry.get('path'), left_entry.get('sha256')))}")


def _compare_output_hashes(left: dict, right: dict) -> None:
    left_hashes = _vcf_hash_map(left)
    right_hashes = _vcf_hash_map(right)
    keys = sorted(set(left_hashes) | set(right_hashes))
    if not keys:
        _row("VCF hashes", "not available")
        return
    for key in keys:
        left_item = left_hashes.get(key)
        right_item = right_hashes.get(key)
        if not left_item or not right_item:
            _row(_short_path(key), "missing")
            continue
        if left_item.get("call_sha256") or right_item.get("call_sha256"):
            _compare_row(
                f"{_short_path(key)} calls",
                left_item.get("call_sha256"),
                right_item.get("call_sha256"),
            )
        _compare_row(
            f"{_short_path(key)} strict records",
            left_item.get("normalized_sha256"),
            right_item.get("normalized_sha256"),
        )


def _inventory_size_detail(left, right) -> str:
    return f"{_comparison_value(_format_optional_bytes(left))} != {_comparison_value(_format_optional_bytes(right))}"


def _compare_inventory_size(left: dict, right: dict) -> None:
    left_size = _inventory_total_bytes(left)
    right_size = _inventory_total_bytes(right)
    if left_size is None and right_size is None:
        _row("Inventory size", "not available")
        return
    if left_size is None or right_size is None:
        _row("Inventory size", f"missing: {_inventory_size_detail(left_size, right_size)}")
        return
    if left_size == right_size:
        _row("Inventory size", f"same: {_comparison_value(_format_optional_bytes(left_size))}")
        return
    left_manifest = _inventory_manifest_hash(left)
    right_manifest = _inventory_manifest_hash(right)
    if left_manifest and left_manifest == right_manifest:
        _row("Inventory size", f"note: {_inventory_size_detail(left_size, right_size)}; file list unchanged")
        return
    _row("Inventory size", f"different: {_inventory_size_detail(left_size, right_size)}")


def _multi_status(label: str, baseline, reports: List[dict], value_fn) -> None:
    baseline_value = value_fn(baseline)
    compared = reports[1:]
    if baseline_value is None and all(value_fn(report) is None for report in compared):
        _row(label, "not available")
        return

    missing = []
    different = []
    for report in compared:
        value = value_fn(report)
        if baseline_value is None or value is None:
            missing.append(_run_label(report))
        elif value != baseline_value:
            different.append(_run_label(report))

    if missing:
        _row(label, "missing: " + ", ".join(missing))
    elif different:
        details = []
        if baseline_value is not None:
            details.append(f"baseline={_comparison_value(baseline_value)}")
        for report in compared:
            value = value_fn(report)
            if value is not None and value != baseline_value:
                details.append(f"{_run_label(report)}={_comparison_value(value)}")
        _row(label, "different: " + "; ".join(details or different))
    else:
        _row(label, f"same: {_comparison_value(baseline_value)}")


def _multi_inventory_size_status(baseline: dict, reports: List[dict]) -> None:
    baseline_size = _inventory_total_bytes(baseline)
    compared = reports[1:]
    if baseline_size is None and all(_inventory_total_bytes(report) is None for report in compared):
        _row("Inventory size", "not available")
        return

    missing = []
    different = []
    note_only = True
    baseline_manifest = _inventory_manifest_hash(baseline)
    for report in compared:
        size = _inventory_total_bytes(report)
        if baseline_size is None or size is None:
            missing.append(_run_label(report))
            note_only = False
        elif size != baseline_size:
            different.append(report)
            if not baseline_manifest or _inventory_manifest_hash(report) != baseline_manifest:
                note_only = False

    if missing:
        _row("Inventory size", "missing: " + ", ".join(missing))
    elif different:
        details = [f"baseline={_comparison_value(_format_optional_bytes(baseline_size))}"]
        for report in different:
            details.append(f"{_run_label(report)}={_comparison_value(_format_optional_bytes(_inventory_total_bytes(report)))}")
        status = "note" if note_only else "different"
        suffix = "; file list unchanged" if note_only else ""
        _row("Inventory size", f"{status}: " + "; ".join(details) + suffix)
    else:
        _row("Inventory size", f"same: {_comparison_value(_format_optional_bytes(baseline_size))}")


def _print_multi_run_comparison(reports: List[dict]) -> None:
    baseline = reports[0]
    compared = reports[1:]

    _section("Run Matrix", console.CYAN)
    _row("Runs", len(reports))
    _row("Baseline", _run_label(baseline))
    _row("Compared", ", ".join(_run_label(report) for report in compared))
    _print_run_comparison_legend()

    print()
    _section("Framework", console.BLUE)
    _multi_status("CBIcall ver", baseline, reports, lambda report: _nested(report, "framework", "version"))
    _multi_status("Python ver", baseline, reports, lambda report: _nested(report, "runtime", "python", "version"))
    _multi_status("Java ver", baseline, reports, lambda report: _nested(report, "runtime", "java", "version"))
    _multi_status("Configured Java", baseline, reports, lambda report: _nested(report, "runtime", "configured_java", "version"))
    _multi_status("Backend ver", baseline, reports, lambda report: _nested(report, "runtime", "backend", "version"))

    print()
    _section("Execution", console.BLUE)
    _multi_status("Task count (trace)", baseline, reports, lambda report: _nested(report, "execution_trace", "tasks"))
    _multi_status("Max peak RSS (trace)", baseline, reports, lambda report: _nested(report, "execution_trace", "max_peak_rss", "bytes"))
    _multi_status("Max peak VMEM (trace)", baseline, reports, lambda report: _nested(report, "execution_trace", "max_peak_vmem", "bytes"))

    print()
    _section("Pipeline", console.BLUE)
    _multi_status("Workflow key", baseline, reports, lambda report: _nested(report, "workflow", "key"))
    _multi_status("Registry ver", baseline, reports, lambda report: _nested(report, "workflow", "registry_version"))
    _multi_status("External workflow", baseline, reports, lambda report: _nested(report, "workflow", "metadata", "source"))
    _multi_status("External release", baseline, reports, lambda report: _nested(report, "workflow", "metadata", "release"))
    _multi_status("Entrypoint", baseline, reports, lambda report: _nested(report, "workflow", "entrypoint"))
    _multi_status("Workflow hash", baseline, reports, lambda report: _nested(report, "workflow", "fingerprint"))

    print()
    _section("Execution Contract", console.BLUE)
    _multi_status("Contract hash", baseline, reports, lambda report: _nested(report, "execution_contract", "fingerprint"))
    _multi_status("Command hash", baseline, reports, lambda report: _nested(report, "execution_contract", "normalized_command_sha256"))
    execution_roles = sorted({role for report in reports for role in _execution_file_map(report)})
    if not execution_roles:
        _row("Generated files", "not available")
    for role in execution_roles:
        _multi_status(role, baseline, reports, lambda report, item=role: _execution_file_value(item, report))

    print()
    _section("Software", console.BLUE)
    _multi_status("Software versions", baseline, reports, lambda report: _nested(report, "software_versions", "sha256"))

    print()
    _section("Workflow Files", console.BLUE)
    roles = sorted({role for report in reports for role in _workflow_file_map(report)})
    if not roles:
        _row("Files", "not available")
    for role in roles:
        _multi_status(role, baseline, reports, lambda report, item=role: _multi_workflow_file_value(item, report))

    print()
    _section("Resources", console.BLUE)
    _multi_status("Resource key", baseline, reports, lambda report: _nested(report, "resources", "bundle", "key"))
    _multi_status("Resource ver", baseline, reports, lambda report: _nested(report, "resources", "bundle", "version"))
    _multi_status("Resource hash", baseline, reports, lambda report: _nested(report, "resources", "bundle", "fingerprint"))

    print()
    _section("Outputs", console.BLUE)
    _multi_status("File count", baseline, reports, lambda report: _nested(report, "outputs", "file_inventory", "entries"))
    _multi_inventory_size_status(baseline, reports)
    _multi_status("File inventory", baseline, reports, lambda report: _nested(report, "outputs", "file_inventory", "sha256"))
    keys = sorted({key for report in reports for key in _vcf_hash_map(report)})
    if not keys:
        _row("VCF hashes", "not available")
    for key in keys:
        if any(_multi_vcf_call_value(key, report) is not None for report in reports):
            _multi_status(
                f"{_short_path(key)} calls",
                baseline,
                reports,
                lambda report, item=key: _multi_vcf_call_value(item, report),
            )
        _multi_status(
            f"{_short_path(key)} strict records",
            baseline,
            reports,
            lambda report, item=key: _multi_vcf_hash_value(item, report),
        )


def _pairwise_layer_status(section: dict, left: dict, right: dict) -> tuple:
    statuses = []
    details = []
    for row in section["rows"]:
        status_class, status_label, detail = _row_pair_status(row, left, right)
        statuses.append(status_class)
        if status_class in {"different", "missing", "note"}:
            row_detail = f"{row['label']}={status_label}"
            if detail:
                row_detail += f" ({detail})"
            details.append(row_detail)
    return _aggregate_status(statuses), details


def _print_all_to_all_comparison(reports: List[dict]) -> None:
    _section("All-to-All Matrix", console.CYAN)
    _row("Runs", len(reports))
    _row("Pairs", (len(reports) * (len(reports) - 1)) // 2)
    _print_run_comparison_legend()

    for section in _comparison_sections_with_overall(reports):
        print()
        _section(section["section"], console.BLUE)
        pair_statuses = []
        pair_details = []
        for left_index, left in enumerate(reports):
            for right_index in range(left_index + 1, len(reports)):
                right = reports[right_index]
                pair = f"{_run_label(left)} vs {_run_label(right)}"
                status, details = _pairwise_layer_status(section, left, right)
                pair_statuses.append(status)
                if status == "same":
                    continue
                detail_text = "; ".join(details[:3])
                if len(details) > 3:
                    detail_text += f"; +{len(details) - 3} more"
                pair_details.append(f"{pair}: {status}" + (f" | {detail_text}" if detail_text else ""))
        layer_status = _aggregate_status(pair_statuses)
        if layer_status == "same":
            _row("Status", "same across all pairs")
        elif layer_status == "unavailable":
            _row("Status", "not available across all pairs")
        else:
            _row("Status", layer_status)
            max_details = 8
            for detail in pair_details[:max_details]:
                _row("Pair", detail)
            if len(pair_details) > max_details:
                _row("Pair", f"+{len(pair_details) - max_details} more non-matching pairs")


def _print_run_comparison_legend() -> None:
    print()
    _section("Legend", console.YELLOW)
    _row("same", "values or fingerprints match")
    _row("different", "values or fingerprints exist in both runs but differ")
    _row("missing", "available in only some runs")
    _row("note", "audit hint; not treated as a failed reproducibility check")
    _row("not available", "not recorded in the run reports; task/RAM traces require a backend execution trace")


def _print_run_comparison(left: dict, right: dict) -> None:
    _section("Run Comparison", console.CYAN)
    _row("Run A", _run_label(left))
    _row("Run B", _run_label(right))
    _print_run_comparison_legend()

    print()
    _section("Framework", console.BLUE)
    _compare_row("CBIcall ver", _nested(left, "framework", "version"), _nested(right, "framework", "version"))
    _compare_row("Python ver", _nested(left, "runtime", "python", "version"), _nested(right, "runtime", "python", "version"))
    _compare_row("Java ver", _nested(left, "runtime", "java", "version"), _nested(right, "runtime", "java", "version"))
    _compare_row("Configured Java", _nested(left, "runtime", "configured_java", "version"), _nested(right, "runtime", "configured_java", "version"))
    _compare_row("Backend ver", _nested(left, "runtime", "backend", "version"), _nested(right, "runtime", "backend", "version"))

    print()
    _section("Execution", console.BLUE)
    _compare_row("Task count (trace)", _nested(left, "execution_trace", "tasks"), _nested(right, "execution_trace", "tasks"))
    _compare_row("Max peak RSS (trace)", _nested(left, "execution_trace", "max_peak_rss", "bytes"), _nested(right, "execution_trace", "max_peak_rss", "bytes"))
    _compare_row("Max peak VMEM (trace)", _nested(left, "execution_trace", "max_peak_vmem", "bytes"), _nested(right, "execution_trace", "max_peak_vmem", "bytes"))

    print()
    _section("Pipeline", console.BLUE)
    _compare_row("Workflow key", _nested(left, "workflow", "key"), _nested(right, "workflow", "key"))
    _compare_row("Registry ver", _nested(left, "workflow", "registry_version"), _nested(right, "workflow", "registry_version"))
    _compare_row("External workflow", _nested(left, "workflow", "metadata", "source"), _nested(right, "workflow", "metadata", "source"))
    _compare_row("External release", _nested(left, "workflow", "metadata", "release"), _nested(right, "workflow", "metadata", "release"))
    _compare_row("Entrypoint", _nested(left, "workflow", "entrypoint"), _nested(right, "workflow", "entrypoint"))
    _compare_row("Workflow hash", _nested(left, "workflow", "fingerprint"), _nested(right, "workflow", "fingerprint"))

    print()
    _section("Execution Contract", console.BLUE)
    _compare_row("Contract hash", _nested(left, "execution_contract", "fingerprint"), _nested(right, "execution_contract", "fingerprint"))
    _compare_row("Command hash", _nested(left, "execution_contract", "normalized_command_sha256"), _nested(right, "execution_contract", "normalized_command_sha256"))
    _compare_execution_files(left, right)

    print()
    _section("Software", console.BLUE)
    _compare_row("Software versions", _nested(left, "software_versions", "sha256"), _nested(right, "software_versions", "sha256"))

    print()
    _section("Workflow Files", console.BLUE)
    _compare_workflow_files(left, right)

    print()
    _section("Resources", console.BLUE)
    _compare_row("Resource key", _nested(left, "resources", "bundle", "key"), _nested(right, "resources", "bundle", "key"))
    _compare_row("Resource ver", _nested(left, "resources", "bundle", "version"), _nested(right, "resources", "bundle", "version"))
    _compare_row(
        "Resource hash",
        _nested(left, "resources", "bundle", "fingerprint"),
        _nested(right, "resources", "bundle", "fingerprint"),
    )

    print()
    _section("Outputs", console.BLUE)
    _compare_row("File count", _nested(left, "outputs", "file_inventory", "entries"), _nested(right, "outputs", "file_inventory", "entries"))
    _compare_inventory_size(left, right)
    _compare_row("File inventory", _nested(left, "outputs", "file_inventory", "sha256"), _nested(right, "outputs", "file_inventory", "sha256"))
    _compare_output_hashes(left, right)


def run_compare_runs_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall compare-runs",
        description="Compare CBIcall run-report.json files or run directories.",
    )
    parser.add_argument("runs", nargs="+", help="Run directories or run-report.json files. Provide two or more.")
    parser.add_argument("-o", "--output", help="Write the text comparison report to this file.")
    parser.add_argument(
        "--alias",
        nargs="+",
        default=[],
        metavar="LABEL",
        help="Human-readable labels for runs, in the same order as the run arguments. Accepts space-separated or comma-separated labels.",
    )
    parser.add_argument("--html", help="Write the HTML comparison report to this file. Defaults to <output>.html or compare-runs.html.")
    parser.add_argument(
        "--multiqc",
        nargs="?",
        const=True,
        metavar="MQC_DIR",
        help="Write a MultiQC custom-content directory. Defaults to <output>_mqc/ or compare-runs_mqc/.",
    )
    parser.add_argument(
        "--comparison-view",
        choices=["baseline", "all-to-all", "both"],
        help="Comparison layout. Defaults to baseline for two runs and both for three or more runs.",
    )
    parser.add_argument("--no-html", action="store_true", help="Do not write the default HTML comparison report.")
    parser.add_argument("--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if len(args.runs) < 2:
        parser.error("compare-runs requires at least two run directories or run-report.json files")
    if args.no_html and args.html:
        parser.error("--html and --no-html cannot be used together")
    aliases = []
    for item in args.alias:
        aliases.extend(part.strip() for part in str(item).split(","))
    aliases = [alias for alias in aliases if alias]
    if aliases and len(aliases) != len(args.runs):
        parser.error("--alias must provide one label per run, in the same order as the run arguments")

    html_output = None
    if not args.no_html:
        if args.html:
            html_output = Path(args.html)
        elif args.output:
            html_output = Path(args.output).with_suffix(".html")
        else:
            html_output = Path("compare-runs.html")

    multiqc_output = None
    if args.multiqc is not None:
        if args.multiqc is True:
            multiqc_output = Path(args.output).with_name(Path(args.output).stem + "_mqc") if args.output else Path("compare-runs_mqc")
        else:
            multiqc_output = Path(args.multiqc)

    if args.nocolor or args.output or html_output:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    reports = [_load_run_report(run) for run in args.runs]
    for report, alias in zip(reports, aliases):
        report["_report_alias"] = alias

    comparison_view = args.comparison_view or ("baseline" if len(reports) == 2 else "both")

    buffer = io.StringIO()
    with redirect_stdout(buffer):
        if comparison_view in {"baseline", "both"}:
            if len(reports) == 2:
                _print_run_comparison(reports[0], reports[1])
            else:
                _print_multi_run_comparison(reports)
        if comparison_view in {"all-to-all", "both"}:
            if comparison_view == "both":
                print()
            _print_all_to_all_comparison(reports)
    report_text = buffer.getvalue()
    print(report_text, end="")
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(report_text, encoding="utf-8")
        _row("Report", output)
    if html_output:
        html_output.parent.mkdir(parents=True, exist_ok=True)
        html_output.write_text(render_compare_html(report_text, reports, comparison_view=comparison_view), encoding="utf-8")
        _row("HTML", html_output)
    if multiqc_output:
        write_compare_multiqc_report(reports, multiqc_output)
        _row("MultiQC", multiqc_output)
    return 0
