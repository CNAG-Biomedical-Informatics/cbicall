"""MultiQC custom-content export helpers for CBIcall reports."""

from __future__ import annotations

import html
from pathlib import Path
from typing import Iterable, Optional

import yaml

from .audit_similarity import build_audit_similarity, qualitative_similarity_label
from .report_utils import _comparison_value, _nested


CBICALL_PARENT = {
    "parent_id": "cbicall",
    "parent_name": "CBIcall",
    "parent_description": "CBIcall audit and QC summaries generated from run-report.json.",
}
STATUS_ORDER = ["same", "different", "missing", "note", "unavailable"]
STATUS_COLORS = {
    "same": "#16a34a",
    "different": "#f59e0b",
    "missing": "#64748b",
    "note": "#2563eb",
    "unavailable": "#94a3b8",
    "identical": "#16a34a",
    "near-identical": "#22c55e",
    "mostly-similar": "#84cc16",
    "partly-similar": "#f59e0b",
    "different-or-missing": "#dc2626",
}


def _short_hash(value):
    if not value:
        return None
    text = str(value)
    if len(text) >= 32:
        return f"{text[:12]}...{text[-8:]}"
    return text


def _run_label(payload: dict, report_path: Optional[Path] = None) -> str:
    alias = str(payload.get("_report_alias") or "").strip()
    if alias:
        return alias
    run_id = _nested(payload, "run", "run_id")
    pipeline = _nested(payload, "workflow", "pipeline")
    backend = _nested(payload, "workflow", "backend")
    if backend and pipeline and run_id:
        return f"{backend}_{pipeline}_{run_id}"
    if report_path is not None:
        return report_path.parent.name or report_path.stem
    if payload.get("_report_path"):
        return Path(str(payload["_report_path"])).parent.name
    return "cbicall_run"


def _provider(payload: dict) -> str:
    return (
        _nested(payload, "workflow", "provider")
        or _nested(payload, "workflow", "metadata", "provider")
        or _nested(payload, "execution_contract", "provider")
        or "cbicall"
    )


def _as_float(value):
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_int(value):
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _bytes_to_gib(value):
    number = _as_float(value)
    if number is None:
        return None
    return round(number / (1024 ** 3), 3)


def _bytes_to_mib(value):
    number = _as_float(value)
    if number is None:
        return None
    return round(number / (1024 ** 2), 3)


def _elapsed_minutes(payload: dict):
    seconds = _as_float(payload.get("elapsed_seconds"))
    if seconds is None:
        return None
    return round(seconds / 60, 3)


def _vcf_hash_reports(payload: dict) -> list[dict]:
    reports = _nested(payload, "outputs", "vcf_hash_reports") or []
    return [item for item in reports if isinstance(item, dict)]


def _first_vcf_hash(payload: dict):
    for item in _vcf_hash_reports(payload):
        value = item.get("normalized_sha256") or item.get("sha256")
        if value:
            return value
    return None


def _first_vcf_call_hash(payload: dict):
    for item in _vcf_hash_reports(payload):
        value = item.get("call_sha256")
        if value:
            return value
    return None


def _first_vcf_strict_hash(payload: dict):
    return _first_vcf_hash(payload)


def _first_vcf_records(payload: dict):
    for item in _vcf_hash_reports(payload):
        records = _as_int(item.get("normalized_records"))
        if records is not None:
            return records
    return None



def _similarity_layer(reports: list, layer_name: str) -> Optional[dict]:
    similarity = build_audit_similarity(reports)
    for layer in similarity["layers"]:
        if layer["name"] == layer_name:
            return {"runs": similarity["runs"], "layer": layer}
    return None


def _similarity_cell(reports: list, layer_name: str, left_index: int, right_index: int) -> dict:
    layer = _similarity_layer(reports, layer_name)
    if not layer:
        return {"score": None, "status": "unavailable"}
    return layer["layer"]["rows"][left_index][right_index]


def _rounded_score(value):
    if value is None:
        return None
    return round(float(value), 3)

def _compact_dict(data: dict) -> dict:
    return {key: value for key, value in data.items() if value is not None}


def _html_escape(value) -> str:
    return html.escape("" if value is None else str(value))


def _format_metric(value, fallback: str = "not recorded") -> str:
    if value is None:
        return fallback
    if isinstance(value, float):
        return f"{value:.2f}".rstrip("0").rstrip(".")
    return str(value)


def _dashboard_html(section_id: str, section_name: str, description: str, body: str) -> str:
    return f"""<!--
id: {section_id}
parent_id: cbicall
parent_name: CBIcall
parent_description: CBIcall audit and QC summaries generated from run-report.json.
section_name: {section_name}
description: {description}
plot_type: html
-->
<style>
.cbicall-dashboard {{
  display: grid;
  grid-template-columns: minmax(260px, 1.1fr) minmax(440px, 2fr);
  gap: 18px;
  margin: 8px 0 20px;
}}
.cbicall-status-card,
.cbicall-metrics,
.cbicall-note-card {{
  border: 1px solid #dbe3ef;
  border-radius: 8px;
  background: #ffffff;
  box-shadow: 0 10px 28px rgba(15, 23, 42, 0.08);
}}
.cbicall-status-card {{
  padding: 22px 24px;
  border-left: 8px solid #2563eb;
}}
.cbicall-status-success .cbicall-status-card,
.cbicall-status-same .cbicall-status-card {{
  border-left-color: #16a34a;
}}
.cbicall-status-warning .cbicall-status-card,
.cbicall-status-mixed .cbicall-status-card {{
  border-left-color: #f59e0b;
}}
.cbicall-status-failed .cbicall-status-card,
.cbicall-status-different .cbicall-status-card {{
  border-left-color: #dc2626;
}}
.cbicall-label {{
  display: block;
  color: #64748b;
  font-size: 12px;
  font-weight: 800;
  letter-spacing: .08em;
  text-transform: uppercase;
}}
.cbicall-status-card strong {{
  display: block;
  margin-top: 5px;
  color: #0f172a;
  font-size: 32px;
  line-height: 1.1;
}}
.cbicall-status-card p {{
  margin: 12px 0 0;
  color: #475569;
  font-size: 15px;
}}
.cbicall-metrics {{
  display: grid;
  grid-template-columns: repeat(4, minmax(0, 1fr));
  overflow: hidden;
}}
.cbicall-metrics div {{
  padding: 22px 18px;
  border-left: 1px solid #e2e8f0;
}}
.cbicall-metrics div:first-child {{
  border-left: 0;
}}
.cbicall-metrics span {{
  display: block;
  color: #0f172a;
  font-size: 28px;
  font-weight: 800;
  line-height: 1;
}}
.cbicall-metrics small {{
  display: block;
  margin-top: 8px;
  color: #64748b;
  font-size: 13px;
  font-weight: 700;
}}
.cbicall-note-card {{
  grid-column: 1 / -1;
  padding: 14px 18px;
  color: #334155;
  font-size: 14px;
}}
@media (max-width: 900px) {{
  .cbicall-dashboard {{
    grid-template-columns: 1fr;
  }}
  .cbicall-metrics {{
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }}
}}
</style>
{body}
"""


def _metric_card(value, label: str) -> str:
    return f"<div><span>{_html_escape(_format_metric(value))}</span><small>{_html_escape(label)}</small></div>"


def _write_text(output_dir: Path, filename: str, content: str) -> Path:
    output = output_dir / filename
    output.write_text(content, encoding="utf-8")
    return output


def _write_yaml(output_dir: Path, filename: str, payload: dict) -> Path:
    output = output_dir / filename
    output.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return output


def _section_payload(section_id: str, section_name: str, description: str, plot_type: str, data, *, pconfig: Optional[dict] = None, headers: Optional[dict] = None) -> dict:
    payload = {
        "id": section_id,
        **CBICALL_PARENT,
        "section_name": section_name,
        "description": description,
        "plot_type": plot_type,
        "data": data,
    }
    if pconfig is not None:
        payload["pconfig"] = pconfig
    if headers is not None:
        payload["headers"] = headers
    return payload


def _run_general_stats_row(payload: dict) -> dict:
    inventory = _nested(payload, "outputs", "file_inventory") or {}
    trace = payload.get("execution_trace") or {}
    row = {
        "elapsed_min": _elapsed_minutes(payload),
        "threads": _as_int(_nested(payload, "run", "threads")),
        "output_files": _as_int(inventory.get("entries")),
        "output_size_mib": _bytes_to_mib(inventory.get("total_bytes")),
        "vcf_fingerprints": len(_vcf_hash_reports(payload)),
        "final_vcf_records": _first_vcf_records(payload),
        "canonical_outputs": len(_nested(payload, "outputs", "canonical_outputs") or []),
        "trace_tasks": _as_int(trace.get("tasks")),
        "max_rss_gib": _bytes_to_gib(_nested(trace, "max_peak_rss", "bytes")),
        "max_vmem_gib": _bytes_to_gib(_nested(trace, "max_peak_vmem", "bytes")),
    }
    return _compact_dict(row)


def build_run_overview_html(report_path: Path, payload: dict) -> str:
    status = str(payload.get("status") or "unknown").lower()
    status_class = "success" if status in {"success", "completed", "finished"} else "failed" if status in {"failed", "error"} else "warning"
    workflow = payload.get("workflow") or {}
    run_label = _run_label(payload, report_path)
    inventory = _nested(payload, "outputs", "file_inventory") or {}
    row = _run_general_stats_row(payload)
    final_vcf = _first_vcf_records(payload)
    body = f"""
<div class="cbicall-dashboard cbicall-status-{status_class}">
  <div class="cbicall-status-card">
    <span class="cbicall-label">Run status</span>
    <strong>{_html_escape(status or "unknown")}</strong>
    <p>{_html_escape(run_label)} | {_html_escape(workflow.get("backend"))} / {_html_escape(workflow.get("pipeline"))} / {_html_escape(workflow.get("mode"))}</p>
  </div>
  <div class="cbicall-metrics">
    {_metric_card(row.get("elapsed_min"), "Elapsed min")}
    {_metric_card(inventory.get("entries"), "Output files")}
    {_metric_card(row.get("output_size_mib"), "Output MiB")}
    {_metric_card(final_vcf, "Final VCF records")}
  </div>
  <div class="cbicall-note-card">
    <strong>Audit focus:</strong> CBIcall keeps the full evidence in run-report.json and uses this MultiQC view as a compact companion for run identity, output fingerprints, and native sample QC.
  </div>
</div>
"""
    return _dashboard_html(
        "cbicall_00_run_overview",
        "CBIcall run overview",
        "At-a-glance CBIcall run status, output inventory, and final-output fingerprint summary.",
        body,
    )


def build_run_general_stats_payload(report_path: Path, payload: dict) -> dict:
    data = {_run_label(payload, report_path): _run_general_stats_row(payload)}
    headers = {
        "elapsed_min": {"title": "Elapsed", "description": "Workflow elapsed time", "suffix": " min", "format": "{:.1f}"},
        "threads": {"title": "Threads", "description": "Requested CBIcall workflow threads"},
        "output_files": {"title": "Files", "description": "Files in the audited output inventory"},
        "output_size_mib": {"title": "Output size", "description": "Audited output inventory size", "suffix": " MiB", "format": "{:.1f}"},
        "vcf_fingerprints": {"title": "VCF hashes", "description": "Recorded VCF fingerprint reports"},
        "final_vcf_records": {"title": "VCF records", "description": "Strict-record count in the first final VCF fingerprint report"},
        "canonical_outputs": {"title": "Canonical outputs", "description": "Configured canonical outputs found by CBIcall"},
        "trace_tasks": {"title": "Tasks", "description": "Backend trace task count"},
        "max_rss_gib": {"title": "Max RSS", "description": "Maximum task peak RSS from backend trace", "suffix": " GiB", "format": "{:.2f}"},
        "max_vmem_gib": {"title": "Max VMEM", "description": "Maximum task peak VMEM from backend trace", "suffix": " GiB", "format": "{:.2f}"},
    }
    return _section_payload(
        "cbicall_run_general_stats",
        "CBIcall run general statistics",
        "Numeric CBIcall run and output metrics for MultiQC general statistics.",
        "generalstats",
        data,
        headers=headers,
    )


def build_run_identity_payload(report_path: Path, payload: dict) -> dict:
    workflow = payload.get("workflow") or {}
    runtime = payload.get("runtime") or {}
    backend = runtime.get("backend") or {}
    data = {
        _run_label(payload, report_path): _compact_dict(
            {
                "Status": payload.get("status"),
                "CBIcall": _nested(payload, "framework", "version"),
                "Python": _nested(runtime, "python", "version"),
                "Java": _nested(runtime, "java", "version"),
                "Configured Java": _nested(runtime, "configured_java", "version"),
                "Backend": workflow.get("backend"),
                "Backend version": backend.get("version"),
                "Provider": _provider(payload),
                "Pipeline": workflow.get("pipeline"),
                "Mode": workflow.get("mode"),
                "Genome": _nested(payload, "run", "display_genome") or _nested(payload, "run", "genome"),
                "Resource": _nested(payload, "resources", "bundle", "key"),
                "Resource version": _nested(payload, "resources", "bundle", "version"),
                "Workflow hash": _short_hash(workflow.get("fingerprint")),
                "Resource hash": _short_hash(_nested(payload, "resources", "bundle", "fingerprint")),
                "Contract hash": _short_hash(_nested(payload, "execution_contract", "fingerprint")),
                "Final VCF calls hash": _short_hash(_first_vcf_call_hash(payload)),
                "Final VCF strict hash": _short_hash(_first_vcf_strict_hash(payload)),
            }
        )
    }
    return _section_payload(
        "cbicall_run_identity",
        "CBIcall run identity",
        "Compact identity and audit fingerprints for the CBIcall run.",
        "table",
        data,
        pconfig={"id": "cbicall_run_identity_table", "title": "CBIcall run identity"},
        headers={
            "Status": {
                "title": "Status",
                "description": "CBIcall run status",
                "bgcols": {"success": "#16a34a", "failed": "#dc2626", "error": "#dc2626"},
            },
            "Backend": {"title": "Backend", "description": "Workflow backend used for execution", "scale": False},
            "Provider": {"title": "Provider", "description": "Workflow provider", "scale": False},
            "Workflow hash": {"title": "Workflow hash", "description": "Short workflow fingerprint", "scale": False},
            "Resource hash": {"title": "Resource hash", "description": "Short resource fingerprint", "scale": False},
            "Contract hash": {"title": "Contract hash", "description": "Short execution contract fingerprint", "scale": False},
            "Final VCF calls hash": {"title": "Final VCF calls", "description": "Call-level VCF fingerprint", "scale": False},
            "Final VCF strict hash": {"title": "Final VCF strict", "description": "Strict-record VCF fingerprint", "scale": False},
        },
    )


def _parse_native_coverage_file(path: Path) -> list[dict]:
    try:
        lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    except UnicodeDecodeError:
        return []
    header_index = None
    for idx, line in enumerate(lines):
        if line.startswith("region	sampleID	"):
            header_index = idx
            break
    if header_index is None or header_index + 1 >= len(lines):
        return []
    header = lines[header_index].split("	")
    rows = []
    for line in lines[header_index + 1 :]:
        values = line.split("	")
        if len(values) < 2:
            continue
        row = {key: values[idx] if idx < len(values) else None for idx, key in enumerate(header)}
        rows.append(row)
    return rows


def _parse_native_sex_file(path: Path) -> dict:
    parsed = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        return parsed
    for line in lines:
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        parsed[key.strip()] = value.strip()
    return parsed


def _native_sample_qc_rows(report_path: Path) -> dict:
    stats_dir = report_path.parent / "03_stats"
    rows = {}
    if not stats_dir.is_dir():
        return rows

    for path in sorted(stats_dir.glob("*.coverage.txt")):
        for item in _parse_native_coverage_file(path):
            sample = item.get("sampleID") or path.stem.replace(".coverage", "")
            rows.setdefault(sample, {})
            rows[sample].update(
                _compact_dict(
                    {
                        "Region": item.get("region"),
                        "Mode": item.get("mode"),
                        "Total reads": _as_int(item.get("total_reads")),
                        "Mean coverage": _as_float(item.get("mean_cov")),
                        "10 pct coverage": _as_float(item.get("ten_pct")),
                        "Non-duplicate pct": _as_float(item.get("nondup_pct")),
                        "Insert size": _as_float(item.get("ins_size")),
                        "In-target pct": _as_float(item.get("in_pct")),
                        "Out-target pct": _as_float(item.get("out_pct")),
                    }
                )
            )

    for path in sorted(stats_dir.glob("*.sex.txt")):
        sample = path.stem.replace(".sex", "")
        parsed = _parse_native_sex_file(path)
        rows.setdefault(sample, {})
        rows[sample].update(
            _compact_dict(
                {
                    "Autosome mean depth": _as_float(parsed.get("MEAN DEPTH FOR AUTOSOMES")),
                    "X mean depth": _as_float(parsed.get("MEAN DEPTH FOR X")),
                    "Y mean depth": _as_float(parsed.get("MEAN DEPTH FOR Y")),
                    "Sex threshold": _as_float(parsed.get("THRESHOLD")),
                    "Sex": parsed.get("SEX"),
                }
            )
        )
    return rows


def build_native_sample_qc_payload(report_path: Path) -> Optional[dict]:
    data = _native_sample_qc_rows(report_path)
    if not data:
        return None
    return _section_payload(
        "cbicall_native_sample_qc",
        "CBIcall native sample QC",
        "Sample-level QC values from native CBIcall 03_stats files.",
        "table",
        data,
        pconfig={"id": "cbicall_native_sample_qc_table", "title": "CBIcall native sample QC"},
        headers={
            "Region": {"title": "Region", "description": "Coverage summary region", "scale": False},
            "Mode": {"title": "Mode", "description": "Coverage mode", "scale": False},
            "Total reads": {"title": "Reads", "description": "Reads inspected by the coverage helper", "scale": "Blues"},
            "Mean coverage": {"title": "Mean cov", "description": "Mean coverage in the selected region", "scale": "Greens"},
            "10 pct coverage": {"title": "10 pct cov", "description": "Percentage of positions above 10x", "scale": "Greens"},
            "Non-duplicate pct": {"title": "Nondup pct", "description": "Non-duplicate read percentage", "scale": "Greens"},
            "Insert size": {"title": "Insert size", "description": "Estimated insert size", "scale": "Purples"},
            "In-target pct": {"title": "In-target pct", "description": "Reads in target", "scale": "Greens"},
            "Out-target pct": {"title": "Out-target pct", "description": "Reads outside target", "scale": "Oranges"},
            "Sex": {
                "title": "Sex",
                "description": "Sex inferred for QC from VCF-derived depth ratios",
                "bgcols": {"FEMALE": "#ec4899", "MALE": "#2563eb", "UNKNOWN": "#94a3b8"},
            },
        },
    )


def build_final_outputs_payload(report_path: Path, payload: dict) -> Optional[dict]:
    rows = {}
    for item in _vcf_hash_reports(payload):
        file_name = Path(str(item.get("file") or item.get("path") or "vcf")).name
        label = item.get("name") or file_name
        rows[str(label)] = _compact_dict(
            {
                "File": file_name,
                "Source": item.get("source"),
                "Call records": _as_int(item.get("call_records")),
                "Call hash": _short_hash(item.get("call_sha256")),
                "Strict records": _as_int(item.get("normalized_records")),
                "Strict hash": _short_hash(item.get("normalized_sha256")),
                "Raw hash": _short_hash(item.get("raw_sha256") or item.get("sha256")),
                "Report": _comparison_value(item.get("path")),
            }
        )

    for item in _nested(payload, "outputs", "canonical_outputs") or []:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name") or item.get("pattern") or "canonical_output")
        rows.setdefault(name, {})
        rows[name].update(
            _compact_dict(
                {
                    "Type": item.get("type"),
                    "Status": item.get("status"),
                    "Matches": len(item.get("matches") or []),
                    "Pattern": item.get("pattern"),
                }
            )
        )

    if not rows:
        return None
    return _section_payload(
        "cbicall_final_outputs",
        "CBIcall final outputs",
        "Canonical outputs and final VCF fingerprints recorded by CBIcall.",
        "table",
        rows,
        pconfig={"id": "cbicall_final_outputs_table", "title": "CBIcall final outputs"},
        headers={
            "File": {"title": "File", "description": "Final output file name", "scale": False},
            "Source": {"title": "Source", "description": "Output source", "scale": False},
            "Call records": {"title": "Call records", "description": "Call-level VCF records", "scale": "Blues"},
            "Call hash": {"title": "Call hash", "description": "Short call-level VCF fingerprint", "scale": False},
            "Strict records": {"title": "Strict records", "description": "Strict-record VCF records", "scale": "Purples"},
            "Strict hash": {"title": "Strict hash", "description": "Short strict-record VCF fingerprint", "scale": False},
            "Raw hash": {"title": "Raw hash", "description": "Short raw file fingerprint", "scale": False},
            "Type": {"title": "Type", "description": "Canonical output type", "scale": False},
            "Status": {
                "title": "Status",
                "description": "Canonical output status",
                "bgcols": {"found": "#16a34a", "missing": "#64748b", "expected": "#2563eb"},
            },
        },
    )




def _clear_existing_bundle(output_dir: Path) -> None:
    if not output_dir.is_dir():
        return
    for pattern in ("*_mqc.yaml", "*_mqc.html"):
        for path in output_dir.glob(pattern):
            if path.is_file():
                path.unlink()

def _multiqc_output_dir(default_dir: Path, output_path: Optional[Path] = None) -> Path:
    output = output_path or default_dir
    if output.suffix.lower() in {".yaml", ".yml"}:
        raise ValueError(
            "CBIcall now writes MultiQC custom content as a directory bundle. "
            f"Use a directory path instead of a YAML file: {output}"
        )
    return output


def write_multiqc_report(report_path: Path, payload: dict, output_path: Optional[Path] = None) -> Path:
    output_dir = _multiqc_output_dir(report_path.parent / "cbicall_mqc", output_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    _clear_existing_bundle(output_dir)
    _write_text(output_dir, "cbicall_00_run_overview_mqc.html", build_run_overview_html(report_path, payload))
    sections = [
        ("cbicall_run_general_stats_mqc.yaml", build_run_general_stats_payload(report_path, payload)),
        ("cbicall_run_identity_mqc.yaml", build_run_identity_payload(report_path, payload)),
        ("cbicall_native_sample_qc_mqc.yaml", build_native_sample_qc_payload(report_path)),
        ("cbicall_final_outputs_mqc.yaml", build_final_outputs_payload(report_path, payload)),
    ]
    for filename, section in sections:
        if section:
            _write_yaml(output_dir, filename, section)
    return output_dir


def _status_counts_from_pairs(reports: list, section: dict) -> dict:
    from .report_utils import _aggregate_status, _row_pair_status

    counts = {status: 0 for status in STATUS_ORDER}
    for left_index, left in enumerate(reports):
        for right_index in range(left_index + 1, len(reports)):
            right = reports[right_index]
            statuses = [_row_pair_status(row, left, right)[0] for row in section["rows"]]
            counts[_aggregate_status(statuses)] += 1
    return counts


def _section_status_for_pair(left: dict, right: dict, section: dict) -> str:
    from .report_utils import _aggregate_status, _row_pair_status

    return _aggregate_status([_row_pair_status(row, left, right)[0] for row in section["rows"]])


def _section_status_for_pairs(reports: list, section: dict) -> str:
    from .report_utils import _aggregate_status

    pair_statuses = []
    for left_index, left in enumerate(reports):
        for right_index in range(left_index + 1, len(reports)):
            pair_statuses.append(_section_status_for_pair(left, reports[right_index], section))
    return _aggregate_status(pair_statuses)


def _section_with_kinds(section: dict, kinds: set[str]) -> dict:
    return {**section, "rows": [row for row in section.get("rows", []) if row.get("kind") in kinds]}


def build_compare_general_stats_payload(reports: list) -> dict:
    from .report_utils import _run_label

    data = {_run_label(report): _run_general_stats_row(report) for report in reports}
    headers = build_run_general_stats_payload(Path("run-report.json"), reports[0] if reports else {}).get("headers", {})
    return _section_payload(
        "cbicall_compare_general_stats",
        "CBIcall compared run statistics",
        "Numeric run and output metrics for each compared CBIcall run.",
        "generalstats",
        data,
        headers=headers,
    )


def build_compare_overview_html(reports: list) -> str:
    from .report_utils import _comparison_sections_with_overall

    sections = _comparison_sections_with_overall(reports)
    by_name = {section["section"]: section for section in sections}
    overall = by_name.get("Overall")
    final_vcf = by_name.get("Final VCF")
    pairs = (len(reports) * (len(reports) - 1)) // 2
    overall_status = _section_status_for_pairs(reports, overall) if overall else "unavailable"
    final_status = _section_status_for_pairs(reports, final_vcf) if final_vcf else "unavailable"
    if final_vcf:
        calls_status = _section_status_for_pairs(reports, _section_with_kinds(final_vcf, {"vcf_call"}))
        strict_status = _section_status_for_pairs(reports, _section_with_kinds(final_vcf, {"vcf_strict"}))
    else:
        calls_status = "unavailable"
        strict_status = "unavailable"
    layer = _similarity_layer(reports, "Overall")
    pair_scores = []
    if layer:
        rows = layer["layer"]["rows"]
        for left_index in range(len(rows)):
            for right_index in range(left_index + 1, len(rows)):
                score = rows[left_index][right_index].get("score")
                if score is not None:
                    pair_scores.append(float(score))
    mean_similarity = round(sum(pair_scores) / len(pair_scores), 3) if pair_scores else None
    if calls_status == "same" and strict_status == "same":
        headline = "Final VCF same"
        status_class = "same"
    elif calls_status == "same":
        headline = "Calls same"
        status_class = "mixed"
    elif "different" in {calls_status, strict_status, final_status}:
        headline = "VCF differs"
        status_class = "different"
    else:
        headline = str(final_status)
        status_class = "mixed"
    body = f"""
<div class="cbicall-dashboard cbicall-status-{status_class}">
  <div class="cbicall-status-card">
    <span class="cbicall-label">Final-output status</span>
    <strong>{_html_escape(headline)}</strong>
    <p>{_html_escape(len(reports))} runs | {_html_escape(pairs)} pairwise comparisons | overall audit: {_html_escape(overall_status)}</p>
  </div>
  <div class="cbicall-metrics">
    {_metric_card(len(reports), "Runs")}
    {_metric_card(pairs, "Pairs")}
    {_metric_card(mean_similarity, "Mean similarity")}
    {_metric_card(calls_status, "VCF calls")}
  </div>
  <div class="cbicall-note-card">
    <strong>Strict records:</strong> {_html_escape(strict_status)}. Call-level VCF hashes compare CHROM, POS, REF, ALT, FILTER, and sample genotypes; strict-record hashes also capture QUAL, INFO, FORMAT, annotations, and numeric fields.
  </div>
</div>
"""
    return _dashboard_html(
        "cbicall_00_compare_overview",
        "CBIcall comparison overview",
        "At-a-glance summary of CBIcall pairwise audit similarity and final-output equivalence.",
        body,
    )


def build_compare_pair_summary_payload(reports: list) -> dict:
    from .report_utils import _comparison_sections_with_overall, _run_label

    sections = _comparison_sections_with_overall(reports)
    by_name = {section["section"]: section for section in sections}
    wanted = [
        ("Overall", "Overall"),
        ("Framework", "Framework"),
        ("Pipeline", "Pipeline"),
        ("Execution Contract", "Execution contract"),
        ("Software", "Software"),
        ("Workflow Files", "Workflow files"),
        ("Resources", "Resources"),
        ("Outputs", "Outputs"),
        ("Final VCF", "Final VCF"),
    ]
    rows = {}
    for left_index, left in enumerate(reports):
        for right_index in range(left_index + 1, len(reports)):
            right = reports[right_index]
            left_label = _run_label(left)
            right_label = _run_label(right)
            overall_cell = _similarity_cell(reports, "Overall", left_index, right_index)
            final_vcf_cell = _similarity_cell(reports, "Final VCF", left_index, right_index)
            overall_status = _section_status_for_pair(left, right, by_name["Overall"])
            row = {
                "Run A": left_label,
                "Run B": right_label,
                "Overall category": qualitative_similarity_label(overall_status, overall_cell.get("score")),
                "Overall similarity": _rounded_score(overall_cell.get("score")),
            }
            if "Final VCF" in by_name:
                final_vcf_status = _section_status_for_pair(left, right, by_name["Final VCF"])
                row["Final VCF category"] = qualitative_similarity_label(final_vcf_status, final_vcf_cell.get("score"))
                row["Final VCF similarity"] = _rounded_score(final_vcf_cell.get("score"))
                row["Final VCF calls"] = _section_status_for_pair(
                    left, right, _section_with_kinds(by_name["Final VCF"], {"vcf_call"})
                )
                row["Final VCF strict records"] = _section_status_for_pair(
                    left, right, _section_with_kinds(by_name["Final VCF"], {"vcf_strict"})
                )
            for section_name, label in wanted:
                section = by_name.get(section_name)
                if section:
                    row[label] = _section_status_for_pair(left, right, section)
            rows[f"{left_label} vs {right_label}"] = row
    return _section_payload(
        "cbicall_compare_pair_summary",
        "CBIcall pairwise comparison summary",
        "Pair-level CBIcall reproducibility statuses across audit layers.",
        "table",
        rows,
        pconfig={"id": "cbicall_compare_pair_summary_table", "title": "CBIcall pairwise comparison summary"},
        headers={
            "Run A": {"title": "Run A", "description": "First compared run", "scale": False},
            "Run B": {"title": "Run B", "description": "Second compared run", "scale": False},
            "Overall category": {"title": "Overall category", "description": "Qualitative report-level similarity", "bgcols": STATUS_COLORS},
            "Overall similarity": {"title": "Overall similarity", "description": "Report-level Jaccard similarity", "scale": "YlGnBu"},
            "Final VCF category": {"title": "Final VCF category", "description": "Qualitative final VCF similarity", "bgcols": STATUS_COLORS},
            "Final VCF similarity": {"title": "Final VCF similarity", "description": "Final VCF Jaccard similarity", "scale": "YlGnBu"},
            "Final VCF calls": {"title": "VCF calls", "description": "Call-level final VCF status", "bgcols": STATUS_COLORS},
            "Final VCF strict records": {"title": "VCF strict", "description": "Strict-record final VCF status", "bgcols": STATUS_COLORS},
            **{label: {"title": label, "description": f"{label} audit-layer status", "bgcols": STATUS_COLORS} for _, label in wanted},
        },
    )


def build_compare_status_counts_payload(reports: list) -> dict:
    from .report_utils import _comparison_sections_with_overall

    sections = _comparison_sections_with_overall(reports)
    by_name = {section["section"]: section for section in sections}
    pairs = (len(reports) * (len(reports) - 1)) // 2
    overall_counts = _status_counts_from_pairs(reports, by_name["Overall"])
    row = {
        "Runs": len(reports),
        "Pairs": pairs,
        "Overall": _section_status_for_pairs(reports, by_name["Overall"]),
    }
    for status in STATUS_ORDER:
        row[f"Pairs {status}"] = overall_counts.get(status, 0)
    if "Final VCF" in by_name:
        final_vcf_section = by_name["Final VCF"]
        vcf_counts = _status_counts_from_pairs(reports, final_vcf_section)
        vcf_call_section = _section_with_kinds(final_vcf_section, {"vcf_call"})
        vcf_strict_section = _section_with_kinds(final_vcf_section, {"vcf_strict"})
        row["Final VCF"] = _section_status_for_pairs(reports, by_name["Final VCF"])
        row["Final VCF calls"] = _section_status_for_pairs(reports, vcf_call_section)
        row["Final VCF strict records"] = _section_status_for_pairs(reports, vcf_strict_section)
        for status in ["same", "different", "missing"]:
            row[f"Final VCF {status}"] = vcf_counts.get(status, 0)
        for status, count in _status_counts_from_pairs(reports, vcf_call_section).items():
            if status in {"same", "different", "missing"}:
                row[f"Final VCF calls {status}"] = count
        for status, count in _status_counts_from_pairs(reports, vcf_strict_section).items():
            if status in {"same", "different", "missing"}:
                row[f"Final VCF strict records {status}"] = count
    return _section_payload(
        "cbicall_compare_status_counts",
        "CBIcall comparison status counts",
        "Aggregate counts of pairwise CBIcall comparison statuses.",
        "table",
        {"comparison": row},
        pconfig={"id": "cbicall_compare_status_counts_table", "title": "CBIcall comparison status counts"},
        headers={
            "Overall": {"title": "Overall", "description": "Aggregate overall status", "bgcols": STATUS_COLORS},
            "Final VCF": {"title": "Final VCF", "description": "Aggregate final VCF status", "bgcols": STATUS_COLORS},
            "Final VCF calls": {"title": "VCF calls", "description": "Aggregate call-level VCF status", "bgcols": STATUS_COLORS},
            "Final VCF strict records": {"title": "VCF strict", "description": "Aggregate strict-record VCF status", "bgcols": STATUS_COLORS},
            "Runs": {"title": "Runs", "description": "Compared runs", "scale": "Blues"},
            "Pairs": {"title": "Pairs", "description": "Pairwise comparisons", "scale": "Blues"},
            **{f"Pairs {status}": {"title": f"Pairs {status}", "description": f"Pair count with {status} overall status", "scale": False} for status in STATUS_ORDER},
        },
    )



def build_compare_similarity_heatmap_payload(reports: list) -> dict:
    from .report_utils import _run_label

    labels = [_run_label(report) for report in reports]
    layer = _similarity_layer(reports, "Overall")
    if not layer:
        data = []
    else:
        data = [
            [_rounded_score(cell.get("score")) if cell.get("score") is not None else 0 for cell in row]
            for row in layer["layer"]["rows"]
        ]
    return _section_payload(
        "cbicall_compare_similarity_heatmap",
        "CBIcall audit similarity heatmap",
        "Overall report-level Jaccard similarity across compared CBIcall runs. Exact field-level evidence remains in the CBIcall HTML report.",
        "heatmap",
        data,
        pconfig={
            "id": "cbicall_compare_similarity_heatmap",
            "title": "CBIcall overall audit similarity",
            "xcats_samples": labels,
            "ycats_samples": labels,
            "xlab": "Run",
            "ylab": "Run",
            "zlab": "Jaccard similarity",
            "min": 0,
            "max": 1,
            "display_values": True,
        },
    )

def write_compare_multiqc_report(reports: list, output_path: Path) -> Path:
    output_dir = _multiqc_output_dir(output_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    _clear_existing_bundle(output_dir)
    _write_text(output_dir, "cbicall_00_compare_overview_mqc.html", build_compare_overview_html(reports))
    sections = [
        ("cbicall_compare_general_stats_mqc.yaml", build_compare_general_stats_payload(reports)),
        ("cbicall_compare_pair_summary_mqc.yaml", build_compare_pair_summary_payload(reports)),
        ("cbicall_compare_status_counts_mqc.yaml", build_compare_status_counts_payload(reports)),
        ("cbicall_compare_similarity_heatmap_mqc.yaml", build_compare_similarity_heatmap_payload(reports)),
    ]
    for filename, section in sections:
        _write_yaml(output_dir, filename, section)
    return output_dir
