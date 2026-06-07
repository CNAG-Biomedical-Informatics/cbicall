"""MultiQC custom-content export helpers for CBIcall reports."""

from pathlib import Path

import yaml

from .report_utils import _format_bytes, _nested


def _short_hash(value):
    if not value:
        return None
    text = str(value)
    if len(text) >= 32:
        return f"{text[:12]}...{text[-8:]}"
    return text


def _run_label(payload: dict, report_path: Path) -> str:
    run_id = _nested(payload, "run", "run_id")
    pipeline = _nested(payload, "workflow", "pipeline")
    backend = _nested(payload, "workflow", "backend")
    if backend and pipeline and run_id:
        return f"{backend}_{pipeline}_{run_id}"
    return report_path.parent.name or report_path.stem


def _first_vcf_hash(payload: dict):
    reports = _nested(payload, "outputs", "vcf_hash_reports") or []
    for item in reports:
        if isinstance(item, dict):
            value = item.get("normalized_sha256") or item.get("sha256")
            if value:
                return value
    return None


def _inventory_size(payload: dict):
    total_bytes = _nested(payload, "outputs", "file_inventory", "total_bytes")
    if total_bytes is None:
        return None
    return _format_bytes(total_bytes)


def build_multiqc_payload(report_path: Path, payload: dict) -> dict:
    """Build a MultiQC custom-content table for one CBIcall run report."""
    sample = _run_label(payload, report_path)
    workflow_metadata = _nested(payload, "workflow", "metadata") or {}
    data = {
        sample: {
            "Status": payload.get("status"),
            "CBIcall": _nested(payload, "framework", "version"),
            "Backend": _nested(payload, "workflow", "backend"),
            "Provider": workflow_metadata.get("provider") or "cbicall",
            "Pipeline": _nested(payload, "workflow", "pipeline"),
            "Mode": _nested(payload, "workflow", "mode"),
            "Genome": _nested(payload, "run", "display_genome") or _nested(payload, "run", "genome"),
            "Resource": _nested(payload, "resources", "bundle", "key"),
            "Workflow hash": _short_hash(_nested(payload, "workflow", "fingerprint")),
            "Resource hash": _short_hash(_nested(payload, "resources", "bundle", "fingerprint")),
            "Final VCF hash": _short_hash(_first_vcf_hash(payload)),
            "Files": _nested(payload, "outputs", "file_inventory", "entries"),
            "Inventory size": _inventory_size(payload),
            "Elapsed seconds": payload.get("elapsed_seconds"),
        }
    }
    return {
        "id": "cbicall_run_audit",
        "parent_id": "cbicall",
        "parent_name": "CBIcall",
        "parent_description": "CBIcall audit summaries generated from run-report.json.",
        "section_name": "CBIcall run audit",
        "description": "Configuration, resource, workflow, and output fingerprint summary for one CBIcall run.",
        "plot_type": "table",
        "pconfig": {
            "id": "cbicall_run_audit_table",
            "title": "CBIcall run audit",
        },
        "data": data,
    }


def write_multiqc_report(report_path: Path, payload: dict, output_path: Path | None = None) -> Path:
    output = output_path or (report_path.parent / "cbicall_mqc.yaml")
    output.parent.mkdir(parents=True, exist_ok=True)
    mqc_payload = build_multiqc_payload(report_path, payload)
    output.write_text(yaml.safe_dump(mqc_payload, sort_keys=False), encoding="utf-8")
    return output



def _status_counts_from_pairs(reports: list, section: dict) -> dict:
    from .report_utils import _aggregate_status, _row_pair_status

    counts = {"same": 0, "different": 0, "missing": 0, "note": 0, "unavailable": 0}
    for left_index, left in enumerate(reports):
        for right_index in range(left_index + 1, len(reports)):
            right = reports[right_index]
            statuses = [_row_pair_status(row, left, right)[0] for row in section["rows"]]
            counts[_aggregate_status(statuses)] += 1
    return counts


def _section_status_for_pairs(reports: list, section: dict) -> str:
    from .report_utils import _aggregate_status, _row_pair_status

    pair_statuses = []
    for left_index, left in enumerate(reports):
        for right_index in range(left_index + 1, len(reports)):
            right = reports[right_index]
            statuses = [_row_pair_status(row, left, right)[0] for row in section["rows"]]
            pair_statuses.append(_aggregate_status(statuses))
    return _aggregate_status(pair_statuses)


def build_compare_multiqc_payload(reports: list) -> dict:
    """Build a MultiQC custom-content table for a compare-runs report."""
    from .report_utils import _comparison_sections_with_overall, _run_label

    sections = _comparison_sections_with_overall(reports)
    by_name = {section["section"]: section for section in sections}
    pairs = (len(reports) * (len(reports) - 1)) // 2
    sample = "compare_runs"
    if reports:
        sample = f"compare_{_run_label(reports[0])}"

    row = {
        "Runs": len(reports),
        "Pairs": pairs,
        "Overall": _section_status_for_pairs(reports, by_name["Overall"]),
        "Framework": _section_status_for_pairs(reports, by_name["Framework"]),
        "Pipeline": _section_status_for_pairs(reports, by_name["Pipeline"]),
        "Execution contract": _section_status_for_pairs(reports, by_name["Execution Contract"]),
        "Software": _section_status_for_pairs(reports, by_name["Software"]),
        "Workflow files": _section_status_for_pairs(reports, by_name["Workflow Files"]),
        "Resources": _section_status_for_pairs(reports, by_name["Resources"]),
        "Outputs": _section_status_for_pairs(reports, by_name["Outputs"]),
    }
    if "Final VCF" in by_name:
        row["Final VCF"] = _section_status_for_pairs(reports, by_name["Final VCF"])

    overall_counts = _status_counts_from_pairs(reports, by_name["Overall"])
    for status in ["same", "different", "missing", "note", "unavailable"]:
        row[f"Pairs {status}"] = overall_counts.get(status, 0)

    return {
        "id": "cbicall_compare_audit",
        "parent_id": "cbicall",
        "parent_name": "CBIcall",
        "parent_description": "CBIcall audit summaries generated from run-report.json and compare-runs.",
        "section_name": "CBIcall run comparison",
        "description": "Baseline and all-to-all reproducibility summary from cbicall compare-runs.",
        "plot_type": "table",
        "pconfig": {
            "id": "cbicall_compare_audit_table",
            "title": "CBIcall run comparison",
        },
        "data": {sample: row},
    }


def write_compare_multiqc_report(reports: list, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mqc_payload = build_compare_multiqc_payload(reports)
    output_path.write_text(yaml.safe_dump(mqc_payload, sort_keys=False), encoding="utf-8")
    return output_path
