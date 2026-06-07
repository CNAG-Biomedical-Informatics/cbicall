"""HTML report renderers for CBIcall run and comparison reports."""

import html
import json
from pathlib import Path
from typing import List, Optional
from urllib.parse import quote

from .cli_output import _format_duration, _short_path
from .report_utils import (
    _comparison_value,
    _execution_file_map,
    _execution_file_value,
    _format_bytes,
    _format_optional_bytes,
    _inventory_manifest_hash,
    _inventory_total_bytes,
    _multi_vcf_hash_value,
    _multi_workflow_file_value,
    _nested,
    _vcf_hash_map,
    _workflow_file_map,
)


def _html_row(label: str, value) -> str:
    if value is None or value == "":
        value = "(undef)"
    return _html_row_raw(label, html.escape(str(value)))


def _html_row_raw(label: str, value_html: str) -> str:
    return (
        "<tr>"
        f"<th>{html.escape(str(label))}</th>"
        f"<td>{value_html}</td>"
        "</tr>"
    )


def _html_section(title: str, rows: List[tuple]) -> str:
    body = "".join(_html_row(label, value) for label, value in rows)
    return _html_section_raw(title, body)


def _html_section_raw(title: str, body_html: str) -> str:
    return (
        "<section>"
        f"<h2>{html.escape(title)}</h2>"
        f"<table>{body_html}</table>"
        "</section>"
    )


def _display_path(path: str, base_dir: Path) -> str:
    raw_path = Path(str(path))
    if raw_path.is_absolute():
        try:
            return raw_path.relative_to(base_dir).as_posix()
        except ValueError:
            return str(raw_path)
    return raw_path.as_posix()


def _file_link(path: str, base_dir: Path, label: Optional[str] = None) -> str:
    display = label or _display_path(path, base_dir)
    raw_path = Path(str(path))
    if raw_path.is_absolute():
        try:
            href_path = raw_path.relative_to(base_dir).as_posix()
        except ValueError:
            href_path = raw_path.as_posix()
    else:
        href_path = raw_path.as_posix()
    href = quote(href_path)
    return f'<a href="{html.escape(href)}">{html.escape(display)}</a>'


def _is_essential_output(rel_path: str) -> bool:
    if rel_path.startswith(("02_varcall/", "03_stats/", "02_browser/", "01_mtoolbox/")):
        return True
    return rel_path.endswith((".vcf", ".vcf.gz", ".g.vcf.gz", ".html", ".txt")) and not rel_path.startswith("01_bam/")


def _inventory_group_stats(payload: dict) -> List[dict]:
    inventory = _nested(payload, "outputs", "file_inventory") or {}
    directories = inventory.get("directories") or []
    if directories:
        return sorted(
            [
                {
                    "group": str(item.get("group") or "."),
                    "count": int(item.get("count") or 0),
                    "bytes": int(item.get("bytes") or 0),
                }
                for item in directories
                if isinstance(item, dict)
            ],
            key=lambda item: (-item["count"], item["group"]),
        )

    largest = inventory.get("largest_files") or []
    groups = {}
    for rel_path in inventory.get("paths") or []:
        text = str(rel_path)
        group = text.split("/", 1)[0] if "/" in text else "."
        groups.setdefault(group, {"group": group, "count": 0, "bytes": 0})
        groups[group]["count"] += 1
    for item in largest:
        if not isinstance(item, dict):
            continue
        rel_path = item.get("path")
        if not rel_path:
            continue
        group = str(rel_path).split("/", 1)[0] if "/" in str(rel_path) else "."
        groups.setdefault(group, {"group": group, "count": 0, "bytes": 0})
        groups[group]["bytes"] += int(item.get("bytes") or 0)
    return sorted(groups.values(), key=lambda item: (-item["count"], item["group"]))


def _output_dashboard_html(payload: dict, base_dir: Path) -> str:
    inventory = _nested(payload, "outputs", "file_inventory") or {}
    outputs = payload.get("outputs", {})
    canonical = outputs.get("canonical_outputs") or []
    vcf_reports = outputs.get("vcf_hash_reports") or []
    entries = inventory.get("entries")
    total_bytes = inventory.get("total_bytes")
    groups = _inventory_group_stats(payload)
    max_count = max((item["count"] for item in groups), default=0)

    metric_cards = [
        ("Files", entries if entries is not None else "(undef)"),
        ("Inventory Size", _format_bytes(total_bytes) if total_bytes is not None else "(undef)"),
        ("Canonical Outputs", len(canonical)),
        ("VCF Fingerprints", len(vcf_reports)),
    ]
    cards_html = "".join(
        '<div class="output-card">'
        f'<span>{html.escape(str(label))}</span>'
        f'<strong>{html.escape(str(value))}</strong>'
        '</div>'
        for label, value in metric_cards
    )

    bars = []
    for item in groups[:8]:
        width = 0 if max_count == 0 else max(4, round((item["count"] / max_count) * 100))
        label = html.escape(str(item["group"]))
        count = html.escape(str(item["count"]))
        size = html.escape(_format_bytes(item.get("bytes") or 0)) if item.get("bytes") else "size not estimated"
        bars.append(
            '<div class="bar-row">'
            f'<div class="bar-label"><strong>{label}</strong><span>{count} files | {size}</span></div>'
            '<div class="bar-track">'
            f'<div class="bar-fill" style="width:{width}%"></div>'
            '</div>'
            '</div>'
        )
    bars_html = "".join(bars) or '<p class="empty-note">No output inventory was recorded.</p>'

    canonical_cards = []
    for item in canonical[:6]:
        matches = item.get("matches") or []
        status = item.get("status") or "unknown"
        links = "<br>".join(_file_link(match, base_dir) for match in matches[:4]) if matches else "No files matched"
        if len(matches) > 4:
            links += f'<br><span class="muted">+{len(matches) - 4} more</span>'
        canonical_cards.append(
            '<div class="deliverable">'
            f'<span class="deliverable-status">{html.escape(str(status))}</span>'
            f'<strong>{html.escape(str(item.get("name") or "Canonical output"))}</strong>'
            f'<p>{links}</p>'
            '</div>'
        )
    canonical_html = "".join(canonical_cards) or '<p class="empty-note">No registry-declared canonical outputs were recorded.</p>'

    return (
        '<section class="output-dashboard">'
        '<h2>Output Dashboard</h2>'
        f'<div class="output-cards">{cards_html}</div>'
        '<div class="dashboard-grid">'
        '<div class="dashboard-panel"><h3>Inventory Composition</h3>'
        f'{bars_html}</div>'
        '<div class="dashboard-panel"><h3>Canonical Deliverables</h3>'
        f'{canonical_html}</div>'
        '</div>'
        '</section>'
    )


def _output_rows_html(payload: dict, base_dir: Path) -> str:
    rows = []
    seen = set()

    for item in payload.get("outputs", {}).get("canonical_outputs", []):
        matches = item.get("matches") or []
        if not matches:
            continue
        links = []
        for match in matches:
            seen.add(_display_path(match, base_dir))
            links.append(_file_link(match, base_dir))
        rows.append(_html_row_raw(item.get("name", "Canonical output"), "<br>".join(links)))

    for item in payload.get("outputs", {}).get("vcf_hash_reports", []):
        file_path = item.get("file") or item.get("path")
        if not file_path:
            continue
        display = _display_path(file_path, base_dir)
        seen.add(display)
        detail = _file_link(file_path, base_dir)
        normalized_hash = item.get("normalized_sha256")
        if normalized_hash:
            detail += f'<br><span class="muted">normalized SHA-256: {html.escape(str(normalized_hash))}</span>'
        rows.append(_html_row_raw("VCF", detail))

    essential_paths = []
    for rel_path in _nested(payload, "outputs", "file_inventory", "paths") or []:
        if rel_path in seen or not _is_essential_output(str(rel_path)):
            continue
        essential_paths.append(str(rel_path))

    max_outputs = 25
    for rel_path in essential_paths[:max_outputs]:
        seen.add(rel_path)
        rows.append(_html_row_raw("Output", _file_link(rel_path, base_dir)))
    if len(essential_paths) > max_outputs:
        rows.append(_html_row("Output list", f"Showing {max_outputs} of {len(essential_paths)} essential files. See run-report.json for the full inventory."))

    total_bytes = _nested(payload, "outputs", "file_inventory", "total_bytes")
    if total_bytes is not None:
        rows.append(_html_row("Inventory size", _format_bytes(total_bytes)))

    largest = _nested(payload, "outputs", "file_inventory", "largest_files") or []
    for item in largest[:5]:
        rel_path = item.get("path") if isinstance(item, dict) else None
        if not rel_path:
            continue
        detail = (
            f"{_file_link(str(rel_path), base_dir)}"
            f'<br><span class="muted">{html.escape(_format_bytes(item.get("bytes")))}</span>'
        )
        rows.append(_html_row_raw("Largest file", detail))

    if not rows:
        rows.append(_html_row("Outputs", "No essential output files were recorded. See run-report.json for the full inventory."))
    return "".join(rows)


def _run_file_rows_html(payload: dict, base_dir: Path) -> str:
    rows = []
    for path, label in [
        (payload.get("workflow_log"), "Workflow log"),
        (base_dir / "log.json", "Resolved log.json"),
        (base_dir / "run-report.json", "Run report JSON"),
    ]:
        if path and Path(str(path)).is_file():
            rows.append(_html_row_raw(label, _file_link(str(path), base_dir)))
    return "".join(rows)


def _external_report_rows_html(payload: dict, base_dir: Path) -> str:
    rows = []
    outputs = payload.get("outputs", {})
    summary = outputs.get("external_summary") or {}
    pipeline_info = summary.get("pipeline_info") or {}
    multiqc = summary.get("multiqc") or {}

    for key, label in [
        ("params_file", "CBIcall params"),
        ("config_file", "CBIcall config"),
    ]:
        path = outputs.get(key)
        if path:
            rows.append(_html_row_raw(label, _file_link(path, base_dir)))

    for key, label in [
        ("params", "Resolved params"),
        ("trace", "Execution trace"),
        ("report", "Execution report"),
        ("timeline", "Execution timeline"),
        ("dag", "Pipeline DAG"),
        ("manifest", "Manifest"),
        ("software_versions", "Software versions"),
    ]:
        path = pipeline_info.get(key)
        if path:
            rows.append(_html_row_raw(label, _file_link(path, base_dir)))

    if multiqc.get("report"):
        rows.append(_html_row_raw("MultiQC report", _file_link(multiqc["report"], base_dir)))
    if multiqc.get("data_dir"):
        rows.append(_html_row_raw("MultiQC data", _file_link(multiqc["data_dir"], base_dir)))

    return "".join(rows)


def _format_software_version_entry(value) -> str:
    if isinstance(value, dict):
        parts = [f"{key}: {val}" for key, val in sorted(value.items())]
        return "; ".join(parts)
    if isinstance(value, list):
        return "; ".join(str(item) for item in value)
    return str(value)


def _execution_contract_rows_html(payload: dict, base_dir: Path) -> str:
    contract = payload.get("execution_contract") or {}
    if not contract:
        return ""

    rows = []
    if contract.get("path"):
        rows.append(_html_row_raw("Contract", _file_link(contract["path"], base_dir)))
    if contract.get("fingerprint"):
        rows.append(_html_row("Fingerprint", contract.get("fingerprint")))
    if contract.get("normalized_command_sha256"):
        rows.append(_html_row("Command hash", contract.get("normalized_command_sha256")))
    if contract.get("workflow_key"):
        rows.append(_html_row("Workflow key", contract.get("workflow_key")))
    if contract.get("backend"):
        rows.append(_html_row("Backend", contract.get("backend")))
    if contract.get("provider"):
        rows.append(_html_row("Provider", contract.get("provider")))
    for item in contract.get("generated_files") or []:
        if not isinstance(item, dict):
            continue
        role = item.get("role") or "generated"
        link = _file_link(item.get("path"), base_dir) if item.get("path") else "(undef)"
        status = html.escape(str(item.get("status") or "unknown"))
        digest = html.escape(str(item.get("normalized_sha256") or item.get("sha256") or "(undef)"))
        rows.append(_html_row_raw(role, f'{link}<br><span class="muted">{status} | {digest}</span>'))
    if contract.get("status") and contract.get("status") != "present":
        rows.append(_html_row("Status", contract.get("status")))
    if contract.get("error"):
        rows.append(_html_row("Error", contract.get("error")))
    return "".join(rows)


def _execution_trace_rows_html(payload: dict, base_dir: Path) -> str:
    trace = payload.get("execution_trace") or {}
    if not trace:
        return ""

    rows = []
    if trace.get("source"):
        rows.append(_html_row_raw("Source", _file_link(trace["source"], base_dir)))
    if trace.get("sha256"):
        rows.append(_html_row("Fingerprint", trace["sha256"]))
    if trace.get("tasks") is not None:
        rows.append(_html_row("Tasks", trace.get("tasks")))
    if trace.get("status_counts"):
        rows.append(_html_row("Task status", _format_software_version_entry(trace["status_counts"])))
    for key, label in [("max_peak_rss", "Max peak RSS"), ("max_peak_vmem", "Max peak VMEM")]:
        item = trace.get(key) or {}
        if item:
            value = item.get("value") or item.get("bytes")
            detail = html.escape(str(value))
            if item.get("task"):
                detail += f'<br><span class="muted">{html.escape(str(item["task"]))}</span>'
            rows.append(_html_row_raw(label, detail))
    if trace.get("status") and trace.get("status") != "parsed":
        rows.append(_html_row("Status", trace.get("status")))
    return "".join(rows)


def _software_versions_rows_html(payload: dict, base_dir: Path) -> str:
    report = payload.get("software_versions") or {}
    if not report:
        return ""

    rows = []
    if report.get("source"):
        rows.append(_html_row_raw("Source", _file_link(report["source"], base_dir)))
    if report.get("sha256"):
        rows.append(_html_row("Fingerprint", report["sha256"]))
    if report.get("scope"):
        rows.append(_html_row("Scope", report["scope"]))
    if report.get("status") and report.get("status") != "parsed":
        rows.append(_html_row("Status", report.get("status")))
    if report.get("error"):
        rows.append(_html_row("Error", report.get("error")))

    entries = report.get("entries") or {}
    if isinstance(entries, dict):
        for key in sorted(entries):
            rows.append(_html_row(str(key), _format_software_version_entry(entries[key])))

    return "".join(rows)


def _run_report_html(payload: dict) -> str:
    base_dir = Path(str(_nested(payload, "run", "project_dir") or "."))
    workflow_files = payload.get("workflow", {}).get("files", [])
    workflow_rows = [
        (item.get("role", "file"), f"{item.get('status', 'unknown')} | {item.get('path')} | {item.get('sha256')}")
        for item in workflow_files
    ]

    run_file_rows = _run_file_rows_html(payload, base_dir)
    external_report_rows = _external_report_rows_html(payload, base_dir)
    software_version_rows = _software_versions_rows_html(payload, base_dir)
    execution_trace_rows = _execution_trace_rows_html(payload, base_dir)
    execution_contract_rows = _execution_contract_rows_html(payload, base_dir)

    output_dashboard = _output_dashboard_html(payload, base_dir)
    workflow_metadata = _nested(payload, "workflow", "metadata") or {}
    analysis_rows = [
        ("Backend", _nested(payload, "workflow", "backend")),
        ("Pipeline", _nested(payload, "workflow", "pipeline")),
        ("Mode", _nested(payload, "workflow", "mode")),
        ("Genome", _nested(payload, "run", "display_genome") or _nested(payload, "run", "genome")),
        ("Software stack", _nested(payload, "workflow", "software_stack")),
        ("Registry version", _nested(payload, "workflow", "registry_version")),
    ]
    if workflow_metadata.get("provider") == "nf-core":
        analysis_rows.extend(
            [
                ("External workflow", workflow_metadata.get("source")),
                ("External release", workflow_metadata.get("release")),
            ]
        )
    analysis_rows.extend(
        [
            ("Runtime profile", payload.get("profile")),
            ("Threads", _nested(payload, "run", "threads")),
        ]
    )

    sections = [
        _html_section(
            "Run",
            [
                ("Status", payload.get("status")),
                ("Run ID", _nested(payload, "run", "run_id")),
                ("Project", _nested(payload, "run", "project_dir")),
                ("Elapsed", _format_duration(float(payload.get("elapsed_seconds", 0) or 0))),
                ("Workflow log", payload.get("workflow_log")),
            ],
        ),
        _html_section(
            "Analysis",
            analysis_rows,
        ),
        _html_section(
            "Runtime",
            [
                ("CBIcall", _nested(payload, "framework", "version")),
                ("Python", _nested(payload, "runtime", "python", "version")),
                ("Python executable", _nested(payload, "runtime", "python", "executable")),
                ("Java", _nested(payload, "runtime", "java", "version")),
                ("Java path", _nested(payload, "runtime", "java", "path")),
                ("Configured Java", _nested(payload, "runtime", "configured_java", "version")),
                ("Configured Java path", _nested(payload, "runtime", "configured_java", "path")),
                ("Workflow backend", _nested(payload, "runtime", "backend", "name")),
                ("Backend version", _nested(payload, "runtime", "backend", "version")),
                ("Backend path", _nested(payload, "runtime", "backend", "path")),
            ],
        ),
        _html_section(
            "Workflow",
            [
                ("Key", _nested(payload, "workflow", "key")),
                ("Entrypoint", _nested(payload, "workflow", "entrypoint")),
                ("Fingerprint", _nested(payload, "workflow", "fingerprint")),
            ]
            + workflow_rows,
        ),
        _html_section(
            "Resources",
            [
                ("Resource key", _nested(payload, "resources", "bundle", "key")),
                ("Resource version", _nested(payload, "resources", "bundle", "version")),
                ("Resource fingerprint", _nested(payload, "resources", "bundle", "fingerprint")),
            ],
        ),
        _html_section_raw(
            "Outputs",
            _output_rows_html(payload, base_dir),
        ),
    ]
    if run_file_rows:
        sections.insert(-1, _html_section_raw("Run Files", run_file_rows))
    if external_report_rows:
        sections.insert(-1, _html_section_raw("External Reports", external_report_rows))
    if software_version_rows:
        sections.insert(-1, _html_section_raw("Software Versions", software_version_rows))
    if execution_contract_rows:
        sections.insert(-1, _html_section_raw("Execution Contract", execution_contract_rows))
    if execution_trace_rows:
        sections.insert(-1, _html_section_raw("Execution Trace", execution_trace_rows))
    sections.insert(-1, output_dashboard)

    def _section_has_title(section_html: str, title: str) -> bool:
        return f"<h2>{html.escape(title)}</h2>" in section_html

    overview_sections = sections[:2]
    output_sections = [section for section in sections if _section_has_title(section, "Output Dashboard") or _section_has_title(section, "Outputs")]
    evidence_sections = [
        section
        for section in sections[2:]
        if not (_section_has_title(section, "Output Dashboard") or _section_has_title(section, "Outputs"))
    ]
    overview_html = "\n          ".join(overview_sections)
    evidence_html = "\n          ".join(evidence_sections) or '<p class="empty-note">No additional evidence sections were recorded.</p>'
    outputs_html = "\n          ".join(output_sections) or '<p class="empty-note">No output sections were recorded.</p>'
    raw_json = html.escape(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))

    return """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>CBIcall Run Report</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f3f5f8;
      --panel: #ffffff;
      --panel-soft: #f8fafc;
      --text: #17202a;
      --muted: #5f6b7a;
      --border: #d8dee8;
      --accent: #2457a6;
      --accent-soft: #e8f0ff;
      --teal: #0f766e;
      --green: #17633a;
      --amber: #9a5b00;
      --ok-bg: #e8f5ee;
      --ok-text: #17633a;
    }
    body {
      margin: 0;
      background:
        radial-gradient(circle at top left, rgba(36, 87, 166, 0.12), transparent 320px),
        linear-gradient(180deg, #eef3f8 0, var(--bg) 300px);
      color: var(--text);
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 15px;
      line-height: 1.45;
    }
    main {
      width: min(1180px, calc(100% - 32px));
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
      background: linear-gradient(180deg, #ffffff, #fbfcfe);
      border: 1px solid var(--border);
      border-left: 4px solid var(--accent);
      border-radius: 8px;
      padding: 12px 14px;
      box-shadow: 0 8px 24px rgba(16, 24, 40, 0.06);
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
      font-size: 20px;
      line-height: 1.2;
      overflow-wrap: anywhere;
    }
    a {
      color: #2457a6;
      text-decoration: none;
      font-weight: 650;
    }
    a:hover {
      text-decoration: underline;
    }
    .muted {
      color: var(--muted);
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
    }
    .pill {
      display: inline-block;
      margin-top: 4px;
      padding: 4px 9px;
      border-radius: 999px;
      background: var(--ok-bg);
      color: var(--ok-text);
      font-size: 13px;
      font-weight: 750;
    }
    section {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 8px;
      margin: 14px 0;
      overflow: hidden;
      box-shadow: 0 8px 24px rgba(16, 24, 40, 0.05);
    }
    h2 {
      margin: 0;
      padding: 12px 16px;
      border-bottom: 1px solid var(--border);
      background: var(--panel-soft);
      font-size: 16px;
      letter-spacing: 0;
    }
    .output-dashboard {
      background: linear-gradient(180deg, #ffffff, #f8fafc);
    }
    .output-cards {
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 10px;
      padding: 14px 16px;
      border-bottom: 1px solid var(--border);
    }
    .output-card {
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px 14px;
      background: var(--panel);
    }
    .output-card span {
      display: block;
      color: var(--muted);
      font-size: 12px;
      font-weight: 750;
      text-transform: uppercase;
    }
    .output-card strong {
      display: block;
      margin-top: 3px;
      font-size: 18px;
      overflow-wrap: anywhere;
    }
    .dashboard-grid {
      display: grid;
      grid-template-columns: minmax(0, 1.1fr) minmax(0, 0.9fr);
      gap: 14px;
      padding: 16px;
    }
    .dashboard-panel {
      border: 1px solid var(--border);
      border-radius: 8px;
      background: #ffffff;
      padding: 14px;
    }
    h3 {
      margin: 0 0 12px;
      font-size: 14px;
      color: var(--text);
    }
    .bar-row {
      margin: 10px 0;
    }
    .bar-label {
      display: flex;
      align-items: baseline;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 5px;
      color: var(--muted);
      font-size: 12px;
    }
    .bar-label strong {
      color: var(--text);
      font-size: 13px;
    }
    .bar-track {
      height: 10px;
      border-radius: 999px;
      overflow: hidden;
      background: #e9edf3;
    }
    .bar-fill {
      height: 100%;
      border-radius: inherit;
      background: linear-gradient(90deg, var(--accent), var(--teal));
    }
    .deliverable {
      border: 1px solid var(--border);
      border-left: 4px solid var(--teal);
      border-radius: 8px;
      padding: 10px 12px;
      margin: 10px 0;
      background: #fbfcfe;
    }
    .deliverable strong {
      display: block;
      margin: 2px 0 4px;
    }
    .deliverable p {
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
    }
    .deliverable-status {
      display: inline-block;
      padding: 2px 7px;
      border-radius: 999px;
      background: var(--ok-bg);
      color: var(--ok-text);
      font-size: 11px;
      font-weight: 750;
      text-transform: uppercase;
    }
    .empty-note {
      color: var(--muted);
      font-size: 13px;
      padding: 16px;
    }
    .section-note {
      padding: 12px 16px 0;
      font-size: 13px;
    }
    .tabs > input {
      position: absolute;
      opacity: 0;
      pointer-events: none;
    }
    .tab-labels {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin: 0 0 18px;
    }
    .tab-labels label {
      cursor: pointer;
      border: 1px solid var(--border);
      border-radius: 999px;
      background: rgba(255, 255, 255, 0.82);
      padding: 8px 14px;
      color: var(--muted);
      font-size: 13px;
      font-weight: 750;
      box-shadow: 0 5px 18px rgba(16, 24, 40, 0.05);
    }
    #tab-overview:checked ~ .tab-labels label[for="tab-overview"],
    #tab-evidence:checked ~ .tab-labels label[for="tab-evidence"],
    #tab-outputs:checked ~ .tab-labels label[for="tab-outputs"],
    #tab-json:checked ~ .tab-labels label[for="tab-json"] {
      background: var(--accent);
      border-color: var(--accent);
      color: #ffffff;
    }
    .tab-panel {
      display: none;
    }
    #tab-overview:checked ~ .tab-panels .overview-panel,
    #tab-evidence:checked ~ .tab-panels .evidence-panel,
    #tab-outputs:checked ~ .tab-panels .outputs-panel,
    #tab-json:checked ~ .tab-panels .json-panel {
      display: block;
    }
    pre {
      margin: 0;
      border: 1px solid var(--border);
      border-radius: 8px;
      background: #101828;
      color: #edf2f7;
      overflow: auto;
      padding: 16px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
      line-height: 1.5;
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
      width: 210px;
      color: var(--muted);
      font-weight: 650;
      white-space: nowrap;
    }
    td {
      overflow-wrap: anywhere;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 13px;
    }
    @media (max-width: 760px) {
      main { width: min(100% - 20px, 1120px); }
      .summary, .output-cards { grid-template-columns: repeat(2, minmax(0, 1fr)); }
      .dashboard-grid { grid-template-columns: 1fr; }
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
    }
  </style>
</head>
<body>
  <main>
    <header>
      <h1>CBIcall Run Report</h1>
      <p>Human-readable summary generated from <code>run-report.json</code>.</p>
    </header>
    <div class="summary" aria-label="Run summary">
      <div class="metric"><span>Status</span><strong><span class="pill">""" + html.escape(str(payload.get("status", "unknown"))) + """</span></strong></div>
      <div class="metric"><span>Pipeline</span><strong>""" + html.escape(str(_nested(payload, "workflow", "pipeline") or "(undef)")) + """</strong></div>
      <div class="metric"><span>Backend</span><strong>""" + html.escape(str(_nested(payload, "workflow", "backend") or "(undef)")) + """</strong></div>
      <div class="metric"><span>Run ID</span><strong>""" + html.escape(str(_nested(payload, "run", "run_id") or "(undef)")) + """</strong></div>
    </div>
    <div class="tabs">
      <input checked type="radio" name="run-tabs" id="tab-overview">
      <input type="radio" name="run-tabs" id="tab-evidence">
      <input type="radio" name="run-tabs" id="tab-outputs">
      <input type="radio" name="run-tabs" id="tab-json">
      <div class="tab-labels" aria-label="Report views">
        <label for="tab-overview">Overview</label>
        <label for="tab-evidence">Evidence</label>
        <label for="tab-outputs">Outputs</label>
        <label for="tab-json">Raw JSON</label>
      </div>
      <div class="tab-panels">
        <div class="tab-panel overview-panel">
          """ + overview_html + """
        </div>
        <div class="tab-panel evidence-panel">
          """ + evidence_html + """
        </div>
        <div class="tab-panel outputs-panel">
          """ + outputs_html + """
        </div>
        <div class="tab-panel json-panel">
          <pre>""" + raw_json + """</pre>
        </div>
      </div>
    </div>
  </main>
</body>
</html>
"""


def render_run_report_html(payload: dict) -> str:
    return _run_report_html(payload)


def write_run_report_html(report_path: Path, payload: dict, html_path: Optional[Path] = None) -> Path:
    output_path = html_path or report_path.with_suffix(".html")
    output_path.write_text(render_run_report_html(payload), encoding="utf-8")
    return output_path



def _html_status_parts(value: str) -> tuple:
    stripped = value.strip()
    lowered = stripped.lower()
    for status, label in [
        ("unavailable", "not available"),
        ("different", "different"),
        ("missing", "missing"),
        ("same", "same"),
        ("note", "note"),
    ]:
        if lowered == label:
            return status, label, ""
        prefix = label + ":"
        if lowered.startswith(prefix):
            return status, label, stripped[len(prefix):].strip()
    return "neutral", "", stripped



def _compare_status_for_values(baseline, value) -> tuple:
    if baseline is None and value is None:
        return "unavailable", "not available", ""
    if baseline is None or value is None:
        return "missing", "missing", f"{_comparison_value(value)}"
    if value == baseline:
        return "same", "same", f"{_comparison_value(value)}"
    return "different", "different", f"{_comparison_value(value)}"


def _compare_inventory_size_cell(baseline: dict, report: dict) -> tuple:
    baseline_size = _inventory_total_bytes(baseline)
    size = _inventory_total_bytes(report)
    if baseline_size is None and size is None:
        return "unavailable", "not available", ""
    if baseline_size is None or size is None:
        return "missing", "missing", _comparison_value(_format_optional_bytes(size))
    if baseline_size == size:
        return "same", "same", _comparison_value(_format_optional_bytes(size))
    baseline_manifest = _inventory_manifest_hash(baseline)
    if baseline_manifest and baseline_manifest == _inventory_manifest_hash(report):
        return "note", "note", f"{_comparison_value(_format_optional_bytes(size))}; file list unchanged"
    return "different", "different", _comparison_value(_format_optional_bytes(size))


def _compare_matrix_specs(reports: List[dict]) -> List[dict]:
    execution_roles = sorted({role for report in reports for role in _execution_file_map(report)})
    workflow_roles = sorted({role for report in reports for role in _workflow_file_map(report)})
    vcf_keys = sorted({key for report in reports for key in _vcf_hash_map(report)})

    specs = [
        {
            "section": "Framework",
            "rows": [
                ("CBIcall ver", lambda report: _nested(report, "framework", "version")),
                ("Python ver", lambda report: _nested(report, "runtime", "python", "version")),
                ("Java ver", lambda report: _nested(report, "runtime", "java", "version")),
                ("Configured Java", lambda report: _nested(report, "runtime", "configured_java", "version")),
                ("Backend ver", lambda report: _nested(report, "runtime", "backend", "version")),
            ],
        },
        {
            "section": "Execution",
            "rows": [
                ("Task count (trace)", lambda report: _nested(report, "execution_trace", "tasks")),
                ("Max peak RSS (trace)", lambda report: _nested(report, "execution_trace", "max_peak_rss", "bytes")),
                ("Max peak VMEM (trace)", lambda report: _nested(report, "execution_trace", "max_peak_vmem", "bytes")),
            ],
        },
        {
            "section": "Pipeline",
            "rows": [
                ("Workflow key", lambda report: _nested(report, "workflow", "key")),
                ("Registry ver", lambda report: _nested(report, "workflow", "registry_version")),
                ("External workflow", lambda report: _nested(report, "workflow", "metadata", "source")),
                ("External release", lambda report: _nested(report, "workflow", "metadata", "release")),
                ("Entrypoint", lambda report: _nested(report, "workflow", "entrypoint")),
                ("Workflow hash", lambda report: _nested(report, "workflow", "fingerprint")),
            ],
        },
        {
            "section": "Execution Contract",
            "rows": [
                ("Contract hash", lambda report: _nested(report, "execution_contract", "fingerprint")),
                ("Command hash", lambda report: _nested(report, "execution_contract", "normalized_command_sha256")),
                *[(role, lambda report, item=role: _execution_file_value(item, report)) for role in execution_roles],
            ],
        },
        {
            "section": "Software",
            "rows": [
                ("Software versions", lambda report: _nested(report, "software_versions", "sha256")),
            ],
        },
        {
            "section": "Workflow Files",
            "rows": ([(role, lambda report, item=role: _multi_workflow_file_value(item, report)) for role in workflow_roles] or [("Files", lambda report: None)]),
        },
        {
            "section": "Resources",
            "rows": [
                ("Resource key", lambda report: _nested(report, "resources", "bundle", "key")),
                ("Resource ver", lambda report: _nested(report, "resources", "bundle", "version")),
                ("Resource hash", lambda report: _nested(report, "resources", "bundle", "fingerprint")),
            ],
        },
        {
            "section": "Outputs",
            "rows": [
                ("File count", lambda report: _nested(report, "outputs", "file_inventory", "entries")),
                ("Inventory size", lambda report: _inventory_total_bytes(report)),
                ("File inventory", lambda report: _nested(report, "outputs", "file_inventory", "sha256")),
                *[(_short_path(key), lambda report, item=key: _multi_vcf_hash_value(item, report)) for key in vcf_keys],
            ] if vcf_keys else [
                ("File count", lambda report: _nested(report, "outputs", "file_inventory", "entries")),
                ("Inventory size", lambda report: _inventory_total_bytes(report)),
                ("File inventory", lambda report: _nested(report, "outputs", "file_inventory", "sha256")),
                ("VCF hashes", lambda report: None),
            ],
        },
    ]
    return specs


def _matrix_run_label(report: dict) -> str:
    alias = str(report.get("_report_alias") or "").strip()
    if alias:
        return alias
    path = Path(str(report.get("_report_path", "run-report.json")))
    if path.name == "run-report.json" and path.parent.name:
        return path.parent.name
    return _short_path(str(path))


def _matrix_cell_label(status_class: str) -> str:
    return {
        "baseline": "base",
        "same": "same",
        "different": "diff",
        "missing": "miss",
        "note": "note",
        "unavailable": "n/a",
        "neutral": "info",
    }.get(status_class, status_class)


def _matrix_display_value(value: str, fallback: str) -> str:
    value = (value or "").strip()
    if not value:
        return fallback
    if value in {"not available", "(missing)"}:
        return "n/a"
    if len(value) <= 18:
        return value
    if "..." in value and len(value) <= 28:
        return value
    return value[:8] + "..." + value[-6:]


def _matrix_group_summary_html(counts: dict) -> str:
    order = [("different", "diff"), ("missing", "miss"), ("note", "note"), ("unavailable", "n/a"), ("same", "same")]
    parts = []
    for key, label in order:
        count = counts.get(key, 0)
        if count:
            parts.append(f"<span>{count} {html.escape(label)}</span>")
    if not parts:
        parts.append("<span>baseline only</span>")
    return '<span class="matrix-group-counts">' + "".join(parts) + '</span>'


def _matrix_legend_html() -> str:
    items = [
        ("baseline", "base", "baseline value"),
        ("same", "same", "matches baseline"),
        ("different", "diff", "differs from baseline"),
        ("missing", "miss", "missing in at least one run"),
        ("note", "note", "audit note"),
        ("unavailable", "n/a", "not recorded"),
    ]
    return "<div class=\"matrix-legend\">" + "".join(
        f"<span><i class=\"matrix-swatch {klass}\"></i><b>{html.escape(label)}</b>{html.escape(desc)}</span>"
        for klass, label, desc in items
    ) + "</div>"


def _compare_matrix_html(reports: Optional[List[dict]]) -> str:
    if not reports or len(reports) < 2:
        return ""
    baseline = reports[0]
    total_columns = len(reports) + 1
    run_headers = [
        "<th class=\"matrix-run matrix-baseline-head\">"
        f"<span>{html.escape(_matrix_run_label(baseline))}</span>"
        "<small>baseline</small>"
        "</th>"
    ]
    run_headers.extend(
        f"<th class=\"matrix-run\"><span>{html.escape(_matrix_run_label(report))}</span></th>"
        for report in reports[1:]
    )

    body_rows = []
    for spec in _compare_matrix_specs(reports):
        section = spec["section"]
        section_rows = []
        section_counts = {"same": 0, "different": 0, "missing": 0, "note": 0, "unavailable": 0}
        for label, value_fn in spec["rows"]:
            baseline_value = value_fn(baseline)
            baseline_detail = _comparison_value(
                _format_optional_bytes(baseline_value) if label == "Inventory size" else baseline_value
            )
            baseline_label = _matrix_display_value(baseline_detail, _matrix_cell_label("baseline"))
            cells = [
                f"<th class=\"matrix-field\">{html.escape(label)}</th>",
                "<td class=\"matrix-cell baseline\" "
                f"title=\"{html.escape(_matrix_run_label(baseline))}: {html.escape(baseline_detail)}\">"
                f"<span>{html.escape(baseline_label)}</span>"
                "</td>",
            ]
            for report in reports[1:]:
                if label == "Inventory size":
                    status_class, status_label, detail = _compare_inventory_size_cell(baseline, report)
                else:
                    status_class, status_label, detail = _compare_status_for_values(baseline_value, value_fn(report))
                section_counts[status_class] = section_counts.get(status_class, 0) + 1
                title = detail or status_label
                cells.append(
                    f"<td class=\"matrix-cell {status_class}\" "
                    f"title=\"{html.escape(_matrix_run_label(report))}: {html.escape(title)}\">"
                    f"<span>{html.escape(_matrix_cell_label(status_class))}</span>"
                    "</td>"
                )
            section_rows.append("<tr>" + "".join(cells) + "</tr>")
        body_rows.append(
            f"<tr class=\"matrix-group\"><th colspan=\"{total_columns}\"><span class=\"matrix-group-title\">{html.escape(section)}</span>{_matrix_group_summary_html(section_counts)}</th></tr>"
        )
        body_rows.extend(section_rows)

    return (
        "<section class=\"matrix-section\">"
        "<h2>Comparison Heatmap</h2>"
        "<p class=\"section-note\">Rows are audit fields and columns are runs. Colors show each run's status relative to the baseline; hover cells for compact values.</p>"
        + _matrix_legend_html() +
        "<div class=\"matrix-wrap\"><table class=\"matrix-table\">"
        "<thead><tr><th class=\"matrix-field-head\">Audit field</th>" + "".join(run_headers) + "</tr></thead>"
        "<tbody>" + "".join(body_rows) + "</tbody>"
        "</table></div>"
        "</section>"
    )

def render_compare_html(report_text: str, reports: Optional[List[dict]] = None) -> str:
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

    status_counts = {"same": 0, "different": 0, "missing": 0, "note": 0, "unavailable": 0}
    section_html = []
    highlight_html = []
    for section in sections:
        rows = []
        for label, value in section["rows"]:
            status_class, status_label, detail = _html_status_parts(value)
            if status_class in status_counts:
                status_counts[status_class] += 1
            if status_class in {"different", "missing", "note"}:
                detail_html = (
                    f'<code>{html.escape(detail)}</code>'
                    if detail
                    else '<code>No detail recorded</code>'
                )
                highlight_html.append(
                    f"<div class=\"change-item {status_class}\">"
                    "<div class=\"change-head\">"
                    f"<span class=\"change-section\">{html.escape(section['title'])}</span>"
                    f"<span class=\"pill {status_class}\">{html.escape(status_label)}</span>"
                    "</div>"
                    f"<strong>{html.escape(label)}</strong>"
                    f"{detail_html}"
                    "</div>"
                )
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
        if section["title"].lower() == "legend":
            continue
        table = "<table>" + "".join(rows) + "</table>" if rows else ""
        section_html.append(
            "<section>"
            f"<h2>{html.escape(section['title'])}</h2>"
            f"{table}"
            "</section>"
        )

    report_title = sections[0]["title"] if sections else "Run Comparison"
    status_items = [
        ("same", "same", "Values or fingerprints match."),
        ("different", "different", "Values or fingerprints exist in both runs but differ."),
        ("missing", "missing", "Available in only some runs."),
        ("note", "note", "Audit hint; not treated as a failed reproducibility check."),
        ("unavailable", "not available", "Not recorded in the run reports; task/RAM traces require a backend execution trace."),
    ]
    summary_html = "".join(
        "<div class=\"metric\">"
        f"<span class=\"pill {key}\">{html.escape(label)}</span>"
        f"<strong>{status_counts[key]}</strong>"
        f"<p>{html.escape(description)}</p>"
        "</div>"
        for key, label, description in status_items
    )
    highlights = (
        "".join(highlight_html)
        if highlight_html
        else "<p class=\"empty-note\">No differences, missing values, or notes were detected. The full comparison is available in Details.</p>"
    )
    details = "\n          ".join(section_html)
    matrix_html = _compare_matrix_html(reports)
    raw_report = html.escape(report_text)

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
      --note-bg: #eaf2ff;
      --note-text: #2457a6;
      --shadow: 0 18px 45px rgba(15, 23, 42, 0.08);
    }
    body {
      margin: 0;
      background:
        radial-gradient(circle at top left, rgba(36, 87, 166, 0.13), transparent 30%),
        linear-gradient(180deg, #f8fafc, var(--bg));
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
      grid-template-columns: repeat(5, minmax(96px, 1fr));
      gap: 10px;
      margin: 18px 0 20px;
      overflow-x: auto;
      padding-bottom: 2px;
    }
    .metric {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 11px;
      box-shadow: 0 8px 24px rgba(15, 23, 42, 0.04);
    }
    .metric .pill {
      min-width: 0;
      margin-right: 0;
      padding: 4px 10px;
      font-size: 13px;
    }
    .metric strong {
      display: block;
      margin-top: 8px;
      font-size: 26px;
      line-height: 1;
    }
    .metric p {
      margin-top: 7px;
      color: var(--muted);
      font-size: 11px;
      line-height: 1.3;
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
    .empty-note {
      color: var(--muted);
      font-size: 13px;
      padding: 16px;
    }
    .tabs > input {
      position: absolute;
      opacity: 0;
      pointer-events: none;
    }
    .tab-labels {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin: 0 0 18px;
    }
    .tab-labels label {
      cursor: pointer;
      border: 1px solid var(--border);
      border-radius: 999px;
      background: rgba(255, 255, 255, 0.78);
      padding: 8px 14px;
      color: var(--muted);
      font-size: 13px;
      font-weight: 750;
      box-shadow: 0 5px 18px rgba(15, 23, 42, 0.04);
    }
    #tab-overview:checked ~ .tab-labels label[for="tab-overview"],
    #tab-matrix:checked ~ .tab-labels label[for="tab-matrix"],
    #tab-details:checked ~ .tab-labels label[for="tab-details"],
    #tab-raw:checked ~ .tab-labels label[for="tab-raw"] {
      background: var(--accent);
      border-color: var(--accent);
      color: #ffffff;
    }
    .tab-panel {
      display: none;
    }
    #tab-overview:checked ~ .tab-panels .overview-panel,
    #tab-matrix:checked ~ .tab-panels .matrix-panel,
    #tab-details:checked ~ .tab-panels .details-panel,
    #tab-raw:checked ~ .tab-panels .raw-panel {
      display: block;
    }
    .change-list {
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 12px;
      padding: 16px;
    }
    .change-item {
      background: var(--panel);
      border: 1px solid var(--border);
      border-left: 5px solid var(--accent);
      border-radius: 8px;
      padding: 14px;
      box-shadow: var(--shadow);
    }
    .change-item.different { border-left-color: var(--diff-text); }
    .change-item.missing { border-left-color: var(--missing-text); }
    .change-item.note { border-left-color: var(--note-text); }
    .change-head {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 8px;
    }
    .change-section {
      color: var(--muted);
      font-size: 12px;
      font-weight: 750;
      text-transform: uppercase;
    }
    .change-item strong {
      display: block;
      margin-bottom: 10px;
      font-size: 15px;
    }
    .change-item code {
      display: block;
      border-radius: 6px;
      background: #f6f8fb;
      color: var(--text);
      padding: 10px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
      line-height: 1.45;
      overflow-wrap: anywhere;
    }
    .matrix-legend {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      align-items: center;
      padding: 0 16px 14px;
      color: var(--muted);
      font-size: 12px;
    }
    .matrix-legend span {
      display: inline-flex;
      align-items: center;
      gap: 6px;
      white-space: nowrap;
      border: 1px solid var(--border);
      border-radius: 999px;
      background: #ffffff;
      padding: 4px 8px;
    }
    .matrix-legend b {
      color: var(--text);
      font-size: 11px;
      text-transform: uppercase;
    }
    .matrix-swatch {
      display: inline-block;
      width: 13px;
      height: 13px;
      border-radius: 50%;
      box-shadow: inset 0 0 0 1px rgba(15, 23, 42, 0.08);
    }
    .matrix-swatch.baseline { background: #dbeafe; }
    .matrix-swatch.same { background: var(--same-bg); }
    .matrix-swatch.different { background: var(--diff-bg); }
    .matrix-swatch.missing { background: var(--missing-bg); }
    .matrix-swatch.note { background: var(--note-bg); }
    .matrix-swatch.unavailable { background: var(--na-bg); }
    .matrix-wrap {
      overflow-x: auto;
      padding: 0 12px 14px;
    }
    .matrix-table {
      min-width: 640px;
      border-collapse: separate;
      border-spacing: 0;
      table-layout: fixed;
      border: 1px solid var(--border);
      border-radius: 10px;
      overflow: hidden;
      background: #ffffff;
      box-shadow: 0 12px 28px rgba(16, 24, 40, 0.06);
    }
    .matrix-table thead th {
      position: sticky;
      top: 0;
      z-index: 1;
      background: #f8fafc;
      color: var(--muted);
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0;
    }
    .matrix-table th,
    .matrix-table td {
      border-bottom: 1px solid var(--border);
      border-right: 1px solid var(--border);
      padding: 7px 8px;
      white-space: nowrap;
    }
    .matrix-table tr:last-child th,
    .matrix-table tr:last-child td {
      border-bottom: 0;
    }
    .matrix-table tr th:last-child,
    .matrix-table tr td:last-child {
      border-right: 0;
    }
    .matrix-field-head {
      width: 220px;
      text-align: left;
      position: sticky;
      left: 0;
      z-index: 3 !important;
      box-shadow: 7px 0 14px rgba(16, 24, 40, 0.04);
    }
    .matrix-group th {
      padding: 10px 12px;
      background: linear-gradient(90deg, #eef4ff, #f8fbff);
      color: #1e3a8a;
      border-top: 1px solid #d7e4fb;
      border-bottom: 1px solid #d7e4fb;
      font-size: 12px;
      font-weight: 850;
      letter-spacing: 0.02em;
      text-transform: uppercase;
      text-align: left;
    }
    .matrix-group .matrix-group-title {
      display: inline-flex;
      align-items: center;
      gap: 8px;
    }
    .matrix-group .matrix-group-counts {
      float: right;
      display: inline-flex;
      flex-wrap: wrap;
      gap: 6px;
      text-transform: none;
      letter-spacing: 0;
    }
    .matrix-group .matrix-group-counts span {
      border-radius: 999px;
      background: #ffffff;
      border: 1px solid #d7e4fb;
      padding: 2px 7px;
      font-size: 11px;
      font-weight: 750;
    }
    .matrix-field {
      width: 220px;
      color: var(--text);
      background: #fbfdff;
      font-size: 12px;
      font-weight: 750;
      overflow: hidden;
      text-overflow: ellipsis;
      text-align: left;
      position: sticky;
      left: 0;
      z-index: 2;
      box-shadow: 7px 0 14px rgba(16, 24, 40, 0.04);
    }
    .matrix-run {
      width: 118px;
      text-align: center;
      overflow: hidden;
      text-overflow: ellipsis;
    }
    .matrix-baseline-head {
      background: #eff6ff !important;
      color: #1d4ed8 !important;
      box-shadow: inset 0 -3px 0 #2563eb;
    }
    .matrix-run span {
      display: block;
      overflow: hidden;
      text-overflow: ellipsis;
      text-transform: none;
      color: var(--text);
      font-size: 12px;
      font-weight: 800;
    }
    .matrix-run small {
      display: block;
      margin-top: 2px;
      color: var(--muted);
      font-size: 10px;
      text-transform: uppercase;
    }
    .matrix-cell {
      width: 118px;
      text-align: center;
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 11px;
      font-weight: 780;
      text-transform: none;
      letter-spacing: 0;
      box-shadow: inset 0 0 0 3px rgba(255, 255, 255, 0.76);
    }
    .matrix-cell span {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      gap: 5px;
      min-width: 58px;
      max-width: 104px;
      border-radius: 999px;
      background: rgba(255, 255, 255, 0.62);
      padding: 3px 7px;
      overflow: hidden;
      text-overflow: ellipsis;
    }
    .matrix-cell span::before {
      content: "";
      flex: 0 0 auto;
      width: 7px;
      height: 7px;
      border-radius: 50%;
      background: currentColor;
      opacity: 0.85;
    }
    .matrix-cell.baseline span {
      max-width: 108px;
      background: rgba(255, 255, 255, 0.72);
    }
    .matrix-cell.baseline {
      background: linear-gradient(135deg, #dbeafe, #eff6ff);
      color: #1d4ed8;
    }
    .matrix-cell.same {
      background: #e4f6ec;
      color: #166534;
    }
    .matrix-cell.different {
      background: #ffe2a8;
      color: #92400e;
    }
    .matrix-cell.missing {
      background: #ffd9d9;
      color: #991b1b;
    }
    .matrix-cell.note {
      background: #e0ecff;
      color: #1d4ed8;
    }
    .matrix-cell.unavailable {
      background: #eef1f5;
      color: #475569;
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
    a {
      color: #2457a6;
      text-decoration: none;
      font-weight: 650;
    }
    a:hover {
      text-decoration: underline;
    }
    .muted {
      color: var(--muted);
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
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
    .detail {
      display: inline-block;
      color: var(--text);
    }
    .pill + .detail {
      margin-left: 8px;
    }
    .same { background: var(--same-bg); color: var(--same-text); }
    .different { background: var(--diff-bg); color: var(--diff-text); }
    .missing { background: var(--missing-bg); color: var(--missing-text); }
    .note { background: var(--note-bg); color: var(--note-text); }
    .unavailable { background: var(--na-bg); color: var(--na-text); }
    pre {
      margin: 0;
      border: 1px solid var(--border);
      border-radius: 8px;
      background: #101828;
      color: #edf2f7;
      overflow: auto;
      padding: 16px;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
      line-height: 1.5;
    }
    @media (max-width: 760px) {
      main { width: min(100% - 20px, 1120px); }
      .change-list { grid-template-columns: 1fr; }
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
    <div class="summary" aria-label="Status summary and legend">
      """ + summary_html + """
    </div>
    <div class="tabs">
      <input checked type="radio" name="compare-tabs" id="tab-overview">
      <input type="radio" name="compare-tabs" id="tab-matrix">
      <input type="radio" name="compare-tabs" id="tab-details">
      <input type="radio" name="compare-tabs" id="tab-raw">
      <div class="tab-labels" aria-label="Report views">
        <label for="tab-overview">Overview</label>
        <label for="tab-matrix">Matrix</label>
        <label for="tab-details">Details</label>
        <label for="tab-raw">Raw Text</label>
      </div>
      <div class="tab-panels">
        <div class="tab-panel overview-panel">
          <section>
            <h2>Differences, Missing Values, and Notes</h2>
            <p class="section-note">Only changed, missing, or noted fields are listed here. The full comparison is available in Details.</p>
            <div class="change-list">
              """ + highlights + """
            </div>
          </section>
        </div>
        <div class="tab-panel matrix-panel">
          """ + matrix_html + """
        </div>
        <div class="tab-panel details-panel">
          """ + details + """
        </div>
        <div class="tab-panel raw-panel">
          <pre>""" + raw_report + """</pre>
        </div>
      </div>
    </div>
  </main>
</body>
</html>
"""
