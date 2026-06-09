"""Shared audit/report formatting helpers."""

from pathlib import Path

from .cli_output import _short_path


def _nested(data: dict, *keys):
    current = data
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
        if current is None:
            return None
    return current


def _same_text(left, right) -> str:
    return "same" if left == right else "different"


def _run_label(report: dict) -> str:
    alias = str(report.get("_report_alias") or "").strip()
    if alias:
        return alias
    return _short_path(report["_report_path"])


def _short_hash(value):
    if value is None:
        return None
    text = str(value)
    if len(text) >= 32 and all(ch in "0123456789abcdefABCDEF" for ch in text):
        return f"{text[:12]}...{text[-8:]}"
    return text


def _format_bytes(value) -> str:
    try:
        size = int(value)
    except (TypeError, ValueError):
        return str(value)

    if size < 0:
        return str(size)
    units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]
    number = float(size)
    unit = units[0]
    for unit in units:
        if number < 1024 or unit == units[-1]:
            break
        number /= 1024
    if unit == "B":
        human = f"{size} B"
    elif number >= 100:
        human = f"{number:.0f} {unit}"
    elif number >= 10:
        human = f"{number:.1f} {unit}"
    else:
        human = f"{number:.2f} {unit}"
    return f"{human} ({size} bytes)"


def _format_optional_bytes(value):
    if value is None:
        return None
    return _format_bytes(value)


def _comparison_value(value):
    if value is None:
        return "(missing)"
    if isinstance(value, (list, tuple)):
        return ", ".join(_comparison_value(item) for item in value if item is not None)
    text = str(value)
    if "/" in text:
        text = _short_path(text)
    return _short_hash(text)


def _workflow_file_map(report: dict) -> dict:
    files = _nested(report, "workflow", "files") or []
    mapped = {}
    for entry in files:
        if isinstance(entry, dict) and entry.get("role"):
            mapped[str(entry["role"])] = entry
    return mapped


def _execution_file_map(report: dict) -> dict:
    files = _nested(report, "execution_contract", "generated_files") or []
    mapped = {}
    for entry in files:
        if isinstance(entry, dict) and entry.get("role"):
            mapped[str(entry["role"])] = entry
    return mapped


def _execution_file_value(role: str, report: dict):
    entry = _execution_file_map(report).get(role)
    if not entry:
        return None
    return entry.get("normalized_sha256") or entry.get("sha256")


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


def _inventory_total_bytes(report: dict):
    return _nested(report, "outputs", "file_inventory", "total_bytes")


def _inventory_manifest_hash(report: dict):
    return _nested(report, "outputs", "file_inventory", "sha256")


def _multi_workflow_file_value(role: str, report: dict):
    entry = _workflow_file_map(report).get(role)
    if not entry:
        return None
    return entry.get("sha256") or entry.get("path")


def _multi_vcf_hash_value(key: str, report: dict):
    entry = _vcf_hash_map(report).get(key)
    if not entry:
        return None
    return entry.get("normalized_sha256") or entry.get("sha256")



def _comparison_specs(reports: list) -> list:
    execution_roles = sorted({role for report in reports for role in _execution_file_map(report)})
    workflow_roles = sorted({role for report in reports for role in _workflow_file_map(report)})
    vcf_keys = sorted({key for report in reports for key in _vcf_hash_map(report)})

    def row(label, value_fn, kind="value"):
        return {"label": label, "value": value_fn, "kind": kind}

    return [
        {
            "section": "Framework",
            "rows": [
                row("CBIcall ver", lambda report: _nested(report, "framework", "version")),
                row("Python ver", lambda report: _nested(report, "runtime", "python", "version")),
                row("Java ver", lambda report: _nested(report, "runtime", "java", "version")),
                row("Configured Java", lambda report: _nested(report, "runtime", "configured_java", "version")),
                row("Backend ver", lambda report: _nested(report, "runtime", "backend", "version")),
            ],
        },
        {
            "section": "Execution",
            "rows": [
                row("Task count (trace)", lambda report: _nested(report, "execution_trace", "tasks")),
                row("Max peak RSS (trace)", lambda report: _nested(report, "execution_trace", "max_peak_rss", "bytes")),
                row("Max peak VMEM (trace)", lambda report: _nested(report, "execution_trace", "max_peak_vmem", "bytes")),
            ],
        },
        {
            "section": "Pipeline",
            "rows": [
                row("Workflow key", lambda report: _nested(report, "workflow", "key")),
                row("Registry ver", lambda report: _nested(report, "workflow", "registry_version")),
                row("External workflow", lambda report: _nested(report, "workflow", "metadata", "source")),
                row("External release", lambda report: _nested(report, "workflow", "metadata", "release")),
                row("Entrypoint", lambda report: _nested(report, "workflow", "entrypoint")),
                row("Workflow hash", lambda report: _nested(report, "workflow", "fingerprint")),
            ],
        },
        {
            "section": "Execution Contract",
            "rows": [
                row("Contract hash", lambda report: _nested(report, "execution_contract", "fingerprint")),
                row("Command hash", lambda report: _nested(report, "execution_contract", "normalized_command_sha256")),
                *[row(role, lambda report, item=role: _execution_file_value(item, report)) for role in execution_roles],
            ],
        },
        {
            "section": "Software",
            "rows": [
                row("Software versions", lambda report: _nested(report, "software_versions", "sha256")),
            ],
        },
        {
            "section": "Workflow Files",
            "rows": (
                [row(role, lambda report, item=role: _multi_workflow_file_value(item, report)) for role in workflow_roles]
                or [row("Files", lambda report: None)]
            ),
        },
        {
            "section": "Resources",
            "rows": [
                row("Resource key", lambda report: _nested(report, "resources", "bundle", "key")),
                row("Resource ver", lambda report: _nested(report, "resources", "bundle", "version")),
                row("Resource hash", lambda report: _nested(report, "resources", "bundle", "fingerprint")),
            ],
        },
        {
            "section": "Outputs",
            "rows": [
                row("File count", lambda report: _nested(report, "outputs", "file_inventory", "entries")),
                row("Inventory size", lambda report: _inventory_total_bytes(report), kind="inventory_size"),
                row("File inventory", lambda report: _nested(report, "outputs", "file_inventory", "sha256")),
                *[row(_short_path(key), lambda report, item=key: _multi_vcf_hash_value(item, report), kind="vcf") for key in vcf_keys],
            ]
            if vcf_keys
            else [
                row("File count", lambda report: _nested(report, "outputs", "file_inventory", "entries")),
                row("Inventory size", lambda report: _inventory_total_bytes(report), kind="inventory_size"),
                row("File inventory", lambda report: _nested(report, "outputs", "file_inventory", "sha256")),
                row("VCF hashes", lambda report: None, kind="vcf"),
            ],
        },
    ]


def _status_for_values(left, right) -> tuple:
    if left is None and right is None:
        return "unavailable", "not available", ""
    if left is None or right is None:
        return "missing", "missing", f"{_comparison_value(left)} != {_comparison_value(right)}"
    if left == right:
        return "same", "same", f"{_comparison_value(right)}"
    return "different", "different", f"{_comparison_value(left)} != {_comparison_value(right)}"


def _inventory_size_status(left_report: dict, right_report: dict) -> tuple:
    left_size = _inventory_total_bytes(left_report)
    right_size = _inventory_total_bytes(right_report)
    if left_size is None and right_size is None:
        return "unavailable", "not available", ""
    detail = f"{_comparison_value(_format_optional_bytes(left_size))} != {_comparison_value(_format_optional_bytes(right_size))}"
    if left_size is None or right_size is None:
        return "missing", "missing", detail
    if left_size == right_size:
        return "same", "same", _comparison_value(_format_optional_bytes(right_size))
    left_manifest = _inventory_manifest_hash(left_report)
    if left_manifest and left_manifest == _inventory_manifest_hash(right_report):
        return "note", "note", f"{detail}; file list unchanged"
    return "different", "different", detail


def _row_pair_status(row: dict, left_report: dict, right_report: dict) -> tuple:
    if row.get("kind") == "inventory_size":
        return _inventory_size_status(left_report, right_report)
    return _status_for_values(row["value"](left_report), row["value"](right_report))


def _aggregate_status(statuses: list) -> str:
    substantive = [status for status in statuses if status != "unavailable"]
    if not substantive:
        return "unavailable"
    for status in ["different", "missing", "note"]:
        if status in substantive:
            return status
    return "same"


def _comparison_sections_with_overall(reports: list) -> list:
    specs = _comparison_specs(reports)
    all_rows = [row for spec in specs for row in spec["rows"]]
    sections = [{"section": "Overall", "rows": all_rows}, *specs]
    vcf_rows = [row for spec in specs for row in spec["rows"] if row.get("kind") == "vcf"]
    if vcf_rows:
        sections.append({"section": "Final VCF", "rows": vcf_rows})
    return sections
