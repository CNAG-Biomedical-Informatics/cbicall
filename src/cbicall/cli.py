import json
import os
import argparse
import hashlib
import html
import io
import subprocess
import sys
import threading
import time
from contextlib import redirect_stdout
from pathlib import Path
from typing import List

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
from .helpmod import usage, parse_args as _parse_args, parse_run_args as _parse_run_args
from .models import ResolvedConfig, RunSettings
from .resources import validate_resource_catalog
from .workflow_registry import load_workflow_registry
from .goodbye import GoodBye

VERSION = "1.0.1-beta.1"
PROMPT = "Info:"
SPACER = "*" * 41
ARROW = "=>"
AUTHOR = "Author: Manuel Rueda, PhD"
LICENSE = "License: GNU General Public License v3"


# ANSI color helpers: color only for interactive terminals.
def _colors_enabled() -> bool:
    if os.environ.get("ANSI_COLORS_DISABLED") or os.environ.get("NO_COLOR"):
        return False
    return sys.stdout.isatty()


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
    if argv and argv[0] == "run":
        return _parse_run_args(argv[1:], VERSION)
    return _parse_args(argv, VERSION)


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


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


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _workflow_key(workflow) -> str:
    return "/".join(
        [
            str(workflow.engine),
            str(workflow.pipeline),
            str(workflow.mode),
            str(workflow.gatk_version),
            str(workflow.pipeline_version),
        ]
    )


def _workflow_file_manifest(workflow) -> dict:
    candidates = [("entrypoint", workflow.entrypoint)]
    if workflow.config_file:
        candidates.append(("config", workflow.config_file))
    candidates.extend((f"helper:{name}", path) for name, path in sorted(workflow.helpers.items()))

    files = []
    for role, raw_path in candidates:
        if not raw_path:
            continue
        path = Path(str(raw_path))
        entry = {
            "role": role,
            "path": str(raw_path),
            "status": "present" if path.is_file() else "missing",
        }
        if path.is_file():
            entry["sha256"] = _sha256_file(path)
            entry["size_bytes"] = path.stat().st_size
        else:
            entry["sha256"] = None
            entry["size_bytes"] = None
        files.append(entry)

    fingerprint_entries = [
        {
            "role": entry["role"],
            "status": entry["status"],
            "sha256": entry["sha256"],
        }
        for entry in files
    ]
    fingerprint_payload = json.dumps(fingerprint_entries, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return {
        "fingerprint": hashlib.sha256(fingerprint_payload).hexdigest(),
        "files": files,
    }


def _workflow_report(workflow) -> dict:
    report = workflow.to_dict()
    manifest = _workflow_file_manifest(workflow)
    report["key"] = _workflow_key(workflow)
    report["fingerprint"] = manifest["fingerprint"]
    report["files"] = manifest["files"]
    return report


def _parse_key_value_file(path: Path) -> dict:
    data = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip().lower()] = value.strip()
    return data


def _collect_file_inventory(project_dir: Path) -> dict:
    excluded = {"run-report.json"}
    paths = []
    if project_dir.is_dir():
        for path in project_dir.rglob("*"):
            if not path.is_file():
                continue
            rel_path = path.relative_to(project_dir).as_posix()
            if rel_path in excluded:
                continue
            paths.append(rel_path)
    paths.sort()
    manifest_text = "".join(f"{path}\n" for path in paths).encode("utf-8")
    return {
        "algorithm": "sha256",
        "scope": "run directory relative file paths",
        "entries": len(paths),
        "sha256": hashlib.sha256(manifest_text).hexdigest(),
        "excluded": sorted(excluded),
        "paths": paths,
    }


def _collect_output_fingerprints(project_dir: Path) -> dict:
    reports = []
    stats_dir = project_dir / "03_stats"
    if stats_dir.is_dir():
        for path in sorted(stats_dir.glob("*.vcf.sha256.txt")):
            parsed = _parse_key_value_file(path)
            reports.append(
                {
                    "path": str(path),
                    "file": parsed.get("file"),
                    "algorithm": parsed.get("algorithm"),
                    "raw_sha256": parsed.get("raw_sha256"),
                    "normalized_pattern": parsed.get("normalized_pattern"),
                    "normalized_sort": parsed.get("normalized_sort"),
                    "normalized_records": parsed.get("normalized_records"),
                    "normalized_sha256": parsed.get("normalized_sha256"),
                }
            )
    return {
        "file_inventory": _collect_file_inventory(project_dir),
        "vcf_hash_reports": reports,
    }


def write_run_report(
    resolved_config: ResolvedConfig,
    arg: dict,
    params: dict,
    *,
    elapsed_seconds: float,
    workflow_log: Path,
) -> Path:
    """Write a compact execution report next to log.json."""
    project_dir = Path(resolved_config.project_dir)
    project_dir.mkdir(parents=True, exist_ok=True)
    report_path = project_dir / "run-report.json"
    workflow = resolved_config.workflow
    payload = {
        "status": "success",
        "elapsed_seconds": round(elapsed_seconds, 3),
        "workflow_log": str(workflow_log),
        "command_trace": "Detailed workflow stdout/stderr is recorded in the workflow log.",
        "framework": {
            "name": "CBIcall",
            "version": resolved_config.version or VERSION,
        },
        "profile": resolved_config.profile,
        "workflow": _workflow_report(workflow),
        "resources": resolved_config.resources,
        "outputs": _collect_output_fingerprints(project_dir),
        "run": {
            "run_id": resolved_config.run_id,
            "project_dir": resolved_config.project_dir,
            "threads": arg.get("threads"),
            "genome": resolved_config.genome,
            "run_mode": resolved_config.run_mode,
        },
        "inputs": resolved_config.inputs.to_dict(),
        "parameters": {
            "paramfile": arg.get("paramfile"),
            "cleanup_bam": params.get("cleanup_bam", False),
            "workflow_rule": resolved_config.workflow_rule,
            "allow_partial_run": resolved_config.allow_partial_run,
        },
    }
    with report_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2, sort_keys=True)
    return report_path


def _run_validate_param_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall validate-param",
        description="Validate a CBIcall parameters YAML without starting the workflow.",
    )
    parser.add_argument("-p", "--param", dest="paramfile", required=True, help="Parameters YAML file.")
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    params = config_mod.read_param_file(args.paramfile)
    resolved_config = ResolvedConfig.from_mapping({**config_mod.set_config_values(params), "version": VERSION})

    _section("Configuration OK", GREEN)
    workflow = resolved_config.workflow
    bundle = resolved_config.resources.get("bundle", {})
    _row("Param file", _short_path(args.paramfile))
    _row("Profile", resolved_config.profile)
    _row("Workflow", f"{workflow.engine} -> {workflow.pipeline} -> {workflow.mode}")
    _row("GATK", workflow.gatk_version)
    _row("Pipeline ver", workflow.pipeline_version)
    _row("Genome", resolved_config.genome or "b37")
    _row("Entrypoint", _short_path(workflow.entrypoint))
    if workflow.engine == "bash":
        _row("Env file", _short_path(workflow.helpers.get("env")))
    elif workflow.engine == "snakemake":
        _row("Config", _short_path(workflow.config_file))
    _row("Resource key", bundle.get("key"))
    _row("Resource ver", bundle.get("version"))
    _row("Resource hash", bundle.get("fingerprint"))
    runtime_check = bundle.get("runtime_check", {})
    _row("DATADIR", _short_path(runtime_check.get("datadir")))
    _row("Resource check", runtime_check.get("status"))
    return 0


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
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    registry_path = Path(args.registry)
    schema_path = Path(args.schema)
    registry = load_workflow_registry(registry_path, schema_path)
    engines = sorted((registry.get("workflows") or {}).keys())

    _section("Registry OK", GREEN)
    _row("Registry", _short_path(registry_path))
    _row("Schema", _short_path(schema_path))
    _row("Engines", ", ".join(engines) if engines else "(none)")
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
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    registry = load_workflow_registry(Path(args.registry), Path(args.schema))
    summary = validate_resource_catalog(Path(args.catalog), registry, resource_key=args.resource)

    _section("Resources OK", GREEN)
    _row("Catalog", _short_path(summary["path"]))
    _row("Schema", _short_path(summary["schema"]))
    _row("Schema version", summary["schema_version"])
    if summary.get("resource_key"):
        _row("Resource key", summary["resource_key"])
    _row("Resources", summary["resources"])
    _row("Bundle resources", summary["bundle_resources"])
    _row("Compatible workflows", summary["compatible_workflows"])
    return 0


def _run_report_path(path: str) -> Path:
    report_path = Path(path)
    if report_path.is_dir():
        report_path = report_path / "run-report.json"
    if not report_path.is_file():
        raise FileNotFoundError(f"run-report.json not found: {report_path}")
    return report_path


def _load_run_report(path: str) -> dict:
    report_path = _run_report_path(path)
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid run report JSON: {report_path}") from exc
    report["_report_path"] = str(report_path)
    return report


def _nested(data: dict, *keys):
    value = data
    for key in keys:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def _same_text(left, right) -> str:
    return "same" if left == right else "different"


def _run_label(report: dict) -> str:
    return _short_path(report["_report_path"])


def _comparison_value(value):
    return "(missing)" if value is None else value


def _compare_row(label: str, left, right) -> None:
    if left is None and right is None:
        _row(label, "not available")
        return
    if left is None or right is None:
        _row(label, f"missing: {_comparison_value(left)} != {_comparison_value(right)}")
        return
    status = _same_text(left, right)
    value = left if left == right else f"{left} != {right}"
    _row(label, f"{status}: {value}")


def _workflow_file_map(report: dict) -> dict:
    files = _nested(report, "workflow", "files") or []
    mapped = {}
    for entry in files:
        if isinstance(entry, dict) and entry.get("role"):
            mapped[str(entry["role"])] = entry
    return mapped


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
        _row(role, status if not details else f"{status}: {', '.join(details)}")


def _vcf_hash_map(report: dict) -> dict:
    reports = _nested(report, "outputs", "vcf_hash_reports") or []
    mapped = {}
    for item in reports:
        if not isinstance(item, dict):
            continue
        key = item.get("file") or item.get("path")
        if key:
            mapped[Path(str(key)).name] = item
    return mapped


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
        _row(
            _short_path(key),
            _same_text(left_item.get("normalized_sha256"), right_item.get("normalized_sha256")),
        )


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
        _row(label, "different: " + ", ".join(different))
    else:
        _row(label, "same")


def _multi_workflow_file_value(role: str, report: dict):
    entry = _workflow_file_map(report).get(role)
    if not entry:
        return None
    return (entry.get("path"), entry.get("sha256"))


def _multi_vcf_hash_value(key: str, report: dict):
    entry = _vcf_hash_map(report).get(key)
    if not entry:
        return None
    return entry.get("normalized_sha256")


def _print_multi_run_comparison(reports: List[dict]) -> None:
    baseline = reports[0]
    compared = reports[1:]

    _section("Run Matrix", CYAN)
    _row("Runs", len(reports))
    _row("Baseline", _run_label(baseline))
    _row("Compared", ", ".join(_run_label(report) for report in compared))

    print()
    _section("Framework", BLUE)
    _multi_status("CBIcall ver", baseline, reports, lambda report: _nested(report, "framework", "version"))

    print()
    _section("Pipeline", BLUE)
    _multi_status("Workflow key", baseline, reports, lambda report: _nested(report, "workflow", "key"))
    _multi_status("Pipeline ver", baseline, reports, lambda report: _nested(report, "workflow", "pipeline_version"))
    _multi_status("Entrypoint", baseline, reports, lambda report: _nested(report, "workflow", "entrypoint"))
    _multi_status("Workflow hash", baseline, reports, lambda report: _nested(report, "workflow", "fingerprint"))

    print()
    _section("Workflow Files", BLUE)
    roles = sorted({role for report in reports for role in _workflow_file_map(report)})
    if not roles:
        _row("Files", "not available")
    for role in roles:
        _multi_status(role, baseline, reports, lambda report, item=role: _multi_workflow_file_value(item, report))

    print()
    _section("Resources", BLUE)
    _multi_status("Resource key", baseline, reports, lambda report: _nested(report, "resources", "bundle", "key"))
    _multi_status("Resource ver", baseline, reports, lambda report: _nested(report, "resources", "bundle", "version"))
    _multi_status("Resource hash", baseline, reports, lambda report: _nested(report, "resources", "bundle", "fingerprint"))

    print()
    _section("Outputs", BLUE)
    _multi_status("File count", baseline, reports, lambda report: _nested(report, "outputs", "file_inventory", "entries"))
    _multi_status("File inventory", baseline, reports, lambda report: _nested(report, "outputs", "file_inventory", "sha256"))
    keys = sorted({key for report in reports for key in _vcf_hash_map(report)})
    if not keys:
        _row("VCF hashes", "not available")
    for key in keys:
        _multi_status(_short_path(key), baseline, reports, lambda report, item=key: _multi_vcf_hash_value(item, report))
    _print_run_comparison_legend()


def _print_run_comparison_legend() -> None:
    print()
    _section("Legend", BLUE)
    _row("same", "values or fingerprints match")
    _row("different", "values or fingerprints exist in both runs but differ")
    _row("missing", "available in only some runs")
    _row("not available", "not recorded in either run report")


def _print_run_comparison(left: dict, right: dict) -> None:
    _section("Run Comparison", CYAN)
    _row("Run A", _short_path(left["_report_path"]))
    _row("Run B", _short_path(right["_report_path"]))

    print()
    _section("Framework", BLUE)
    _compare_row("CBIcall ver", _nested(left, "framework", "version"), _nested(right, "framework", "version"))

    print()
    _section("Pipeline", BLUE)
    _compare_row("Workflow key", _nested(left, "workflow", "key"), _nested(right, "workflow", "key"))
    _compare_row("Pipeline ver", _nested(left, "workflow", "pipeline_version"), _nested(right, "workflow", "pipeline_version"))
    _compare_row("Entrypoint", _nested(left, "workflow", "entrypoint"), _nested(right, "workflow", "entrypoint"))
    _compare_row("Workflow hash", _nested(left, "workflow", "fingerprint"), _nested(right, "workflow", "fingerprint"))

    print()
    _section("Workflow Files", BLUE)
    _compare_workflow_files(left, right)

    print()
    _section("Resources", BLUE)
    _compare_row("Resource key", _nested(left, "resources", "bundle", "key"), _nested(right, "resources", "bundle", "key"))
    _compare_row("Resource ver", _nested(left, "resources", "bundle", "version"), _nested(right, "resources", "bundle", "version"))
    _compare_row(
        "Resource hash",
        _nested(left, "resources", "bundle", "fingerprint"),
        _nested(right, "resources", "bundle", "fingerprint"),
    )

    print()
    _section("Outputs", BLUE)
    _compare_row("File count", _nested(left, "outputs", "file_inventory", "entries"), _nested(right, "outputs", "file_inventory", "entries"))
    _compare_row("File inventory", _nested(left, "outputs", "file_inventory", "sha256"), _nested(right, "outputs", "file_inventory", "sha256"))
    _compare_output_hashes(left, right)
    _print_run_comparison_legend()


def _html_status_parts(value: str) -> tuple:
    stripped = value.strip()
    lowered = stripped.lower()
    for status, label in [
        ("unavailable", "not available"),
        ("different", "different"),
        ("missing", "missing"),
        ("same", "same"),
    ]:
        if lowered == label:
            return status, label, ""
        prefix = label + ":"
        if lowered.startswith(prefix):
            return status, label, stripped[len(prefix):].strip()
    return "neutral", "", stripped


def _render_compare_html(report_text: str) -> str:
    sections = []
    current = None
    for line in report_text.splitlines():
        if not line.strip():
            continue
        if line.startswith("  ") and "=>" in line and current is not None:
            label, value = line.strip().split("=>", 1)
            current["rows"].append((label.strip(), value.strip()))
            continue
        current = {"title": line.strip(), "rows": []}
        sections.append(current)

    status_counts = {"same": 0, "different": 0, "missing": 0, "unavailable": 0}
    section_html = []
    for section in sections:
        rows = []
        for label, value in section["rows"]:
            status_class, status_label, detail = _html_status_parts(value)
            if status_class in status_counts:
                status_counts[status_class] += 1
            value_html = (
                f"<span class=\"detail\">{html.escape(detail)}</span>"
                if status_class == "neutral"
                else (
                    f"<span class=\"pill {status_class}\">{html.escape(status_label)}</span>"
                    + (f"<span class=\"detail\">{html.escape(detail)}</span>" if detail else "")
                )
            )
            rows.append(
                "<tr>"
                f"<th>{html.escape(label)}</th>"
                f"<td>{value_html}</td>"
                "</tr>"
            )
        table = "<table>" + "".join(rows) + "</table>" if rows else ""
        section_html.append(
            "<section>"
            f"<h2>{html.escape(section['title'])}</h2>"
            f"{table}"
            "</section>"
        )

    report_title = sections[0]["title"] if sections else "Run Comparison"
    summary_html = "".join(
        "<div class=\"metric\">"
        f"<span>{html.escape(label)}</span>"
        f"<strong>{status_counts[key]}</strong>"
        "</div>"
        for key, label in [
            ("same", "Same"),
            ("different", "Different"),
            ("missing", "Missing"),
            ("unavailable", "Not available"),
        ]
    )

    return """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CBIcall Run Comparison</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f6f7f9;
      --panel: #ffffff;
      --panel-soft: #fbfcfe;
      --text: #17202a;
      --muted: #5f6b7a;
      --border: #d8dee8;
      --accent: #2457a6;
      --same-bg: #e8f5ee;
      --same-text: #17633a;
      --diff-bg: #fff2d7;
      --diff-text: #7a4a00;
      --missing-bg: #fde8e8;
      --missing-text: #8a1f1f;
      --na-bg: #eceff3;
      --na-text: #475467;
    }
    body {
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 15px;
      line-height: 1.45;
    }
    main {
      width: min(1120px, calc(100% - 32px));
      margin: 28px auto 40px;
    }
    header {
      margin-bottom: 20px;
      padding: 22px 0 2px;
    }
    h1 {
      margin: 0 0 4px;
      font-size: 28px;
      font-weight: 700;
      letter-spacing: 0;
    }
    p {
      margin: 0;
      color: var(--muted);
    }
    .summary {
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
      margin: 18px 0 20px;
    }
    .metric {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px 14px;
    }
    .metric span {
      display: block;
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
    }
    .metric strong {
      display: block;
      margin-top: 2px;
      font-size: 24px;
      line-height: 1.1;
    }
    section {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 8px;
      margin: 14px 0;
      overflow: hidden;
    }
    h2 {
      margin: 0;
      padding: 12px 16px;
      border-bottom: 1px solid var(--border);
      background: var(--panel-soft);
      font-size: 16px;
      letter-spacing: 0;
    }
    table {
      width: 100%;
      border-collapse: collapse;
    }
    th,
    td {
      padding: 10px 16px;
      border-bottom: 1px solid var(--border);
      text-align: left;
      vertical-align: top;
    }
    tr:last-child th,
    tr:last-child td {
      border-bottom: 0;
    }
    th {
      width: 180px;
      color: var(--muted);
      font-weight: 650;
      white-space: nowrap;
    }
    td {
      overflow-wrap: anywhere;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
    }
    .pill {
      display: inline-block;
      min-width: 82px;
      margin-right: 8px;
      padding: 3px 8px;
      border-radius: 999px;
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 12px;
      font-weight: 750;
      text-align: center;
    }
    .detail { color: var(--text); }
    .same { background: var(--same-bg); color: var(--same-text); }
    .different { background: var(--diff-bg); color: var(--diff-text); }
    .missing { background: var(--missing-bg); color: var(--missing-text); }
    .unavailable { background: var(--na-bg); color: var(--na-text); }
    @media (max-width: 760px) {
      main { width: min(100% - 20px, 1120px); }
      .summary { grid-template-columns: repeat(2, minmax(0, 1fr)); }
      th,
      td {
        display: block;
        width: auto;
      }
      th {
        padding-bottom: 2px;
      }
      td {
        padding-top: 2px;
      }
      .pill {
        margin-bottom: 4px;
      }
    }
  </style>
</head>
<body>
  <main>
    <header>
      <h1>CBIcall """ + html.escape(report_title) + """</h1>
      <p>Static rendering of the text report generated by <code>cbicall compare-runs</code>.</p>
    </header>
    <div class="summary" aria-label="Status summary">
      """ + summary_html + """
    </div>
    """ + "\n    ".join(section_html) + """
  </main>
</body>
</html>
"""


def _run_compare_runs_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall compare-runs",
        description="Compare CBIcall run-report.json files or run directories.",
    )
    parser.add_argument("runs", nargs="+", help="Run directories or run-report.json files. Provide two or more.")
    parser.add_argument("-o", "--output", help="Write the text comparison report to this file.")
    parser.add_argument("--html", help="Write a static HTML rendering of the comparison report.")
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if len(args.runs) < 2:
        parser.error("compare-runs requires at least two run directories or run-report.json files")

    if args.nocolor or args.output or args.html:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    reports = [_load_run_report(run) for run in args.runs]

    buffer = io.StringIO()
    with redirect_stdout(buffer):
        if len(reports) == 2:
            _print_run_comparison(reports[0], reports[1])
        else:
            _print_multi_run_comparison(reports)
    report_text = buffer.getvalue()
    print(report_text, end="")
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(report_text, encoding="utf-8")
        _row("Report", output)
    if args.html:
        output = Path(args.html)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(_render_compare_html(report_text), encoding="utf-8")
        _row("HTML", output)
    return 0


def _run_test_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall test",
        description="Run the bundled CBIcall integration tests from examples/input.",
    )
    parser.add_argument("--wes", action="store_true", help="Run the WES integration test.")
    parser.add_argument("--mit", action="store_true", help="Run the mitochondrial integration test.")
    parser.add_argument("--all", action="store_true", help="Run all bundled integration tests.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use.")
    args = parser.parse_args(argv)

    if args.threads <= 0:
        parser.error("--threads requires a positive integer")

    selected = []
    if args.all or args.wes:
        selected.append("--wes")
    if args.all or args.mit:
        selected.append("--mit")
    if not selected:
        parser.error("select at least one test with --wes, --mit, or --all")

    root = _project_root()
    script = root / "examples" / "input" / "run_tests.sh"
    if not script.is_file():
        raise FileNotFoundError(f"Integration test launcher not found: {script}")

    env = os.environ.copy()
    env["CBICALL"] = str(root / "bin" / "cbicall")
    env["THREADS"] = str(args.threads)
    return subprocess.run(
        ["bash", str(script), *selected],
        cwd=str(script.parent),
        env=env,
        check=False,
    ).returncode


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


def _run_analysis(arg: dict, *, start_time: float, cbicall_path: Path) -> int:
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
            "profile": resolved_config.profile,
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
    workflow_log = Path(resolved_config.project_dir) / log_name
    report_path = write_run_report(
        resolved_config,
        arg,
        params,
        elapsed_seconds=elapsed,
        workflow_log=workflow_log,
    )
    _row("Log", workflow_log)
    _row("Report", report_path)
    if arg.get("verbose"):
        _row("Date", time.ctime())

    gb = GoodBye()
    print(f"{WHITE}{gb.say_goodbye()}{RESET}")

    return 0


def _run_run_command(argv: List[str], *, start_time: float, cbicall_path: Path) -> int:
    arg = vars(_parse_run_args(argv, VERSION))
    return _run_analysis(arg, start_time=start_time, cbicall_path=cbicall_path)


def main() -> int:
    start_time = time.time()
    cbicall_path = Path(sys.argv[0]).resolve()

    if len(sys.argv) > 1 and sys.argv[1] == "run":
        return _run_run_command(sys.argv[2:], start_time=start_time, cbicall_path=cbicall_path)
    if len(sys.argv) > 1 and sys.argv[1] == "validate-param":
        return _run_validate_param_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "validate-registry":
        return _run_validate_registry_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "validate-resources":
        return _run_validate_resources_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "compare-runs":
        return _run_compare_runs_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        return _run_test_command(sys.argv[2:])

    # Backward-compatible legacy run form: cbicall -p params.yaml -t THREADS.
    arg = usage(VERSION)
    return _run_analysis(arg, start_time=start_time, cbicall_path=cbicall_path)
