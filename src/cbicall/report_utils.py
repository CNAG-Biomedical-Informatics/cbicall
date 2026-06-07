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
    return (entry.get("path"), entry.get("sha256"))


def _multi_vcf_hash_value(key: str, report: dict):
    entry = _vcf_hash_map(report).get(key)
    if not entry:
        return None
    return entry.get("normalized_sha256") or entry.get("sha256")
