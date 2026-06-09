"""Report-level audit similarity helpers."""

import json
from pathlib import Path
from statistics import median

from .cli_output import _short_path
from .report_utils import _comparison_sections_with_overall


_MISSING = "__MISSING__"


def _run_label(report: dict, index: int) -> str:
    alias = str(report.get("_report_alias") or "").strip()
    if alias:
        return alias
    path = Path(str(report.get("_report_path") or f"run_{index + 1}"))
    if path.name == "run-report.json" and path.parent.name:
        return path.parent.name
    return _short_path(str(path))


def _looks_like_hash(value: str) -> bool:
    return len(value) >= 32 and all(ch in "0123456789abcdefABCDEF" for ch in value)


def _normalize_path_like(value: str) -> str:
    if "/" not in value:
        return value
    path = Path(value)
    if path.name and path.name != value:
        return path.name
    return _short_path(value)


def _normalize_value(value):
    if value is None:
        return _MISSING
    if isinstance(value, (str, int, float, bool)):
        text = str(value)
        if _looks_like_hash(text):
            return text.lower()
        return _normalize_path_like(text)
    if isinstance(value, (list, tuple)):
        values = [_normalize_value(item) for item in value if item is not None]
        hashes = [item for item in values if isinstance(item, str) and _looks_like_hash(item)]
        if hashes:
            return "|".join(sorted(hashes))
        return json.dumps(values, sort_keys=True, separators=(",", ":"))
    if isinstance(value, dict):
        normalized = {
            str(key): _normalize_value(item)
            for key, item in sorted(value.items())
            if item is not None
        }
        return json.dumps(normalized, sort_keys=True, separators=(",", ":"))
    return _normalize_path_like(str(value))


def _facts_for_section(section: dict, report: dict) -> tuple[set[str], int]:
    facts = set()
    informative = 0
    for row in section["rows"]:
        value = _normalize_value(row["value"](report))
        if value != _MISSING:
            informative += 1
        facts.add(f"{section['section']}|{row['label']}={value}")
    return facts, informative


def _status_for_score(score):
    if score is None:
        return "unavailable"
    if score == 1:
        return "same"
    if score >= 0.95:
        return "high"
    if score >= 0.80:
        return "medium"
    return "low"


def _cell(left: dict, right: dict) -> dict:
    left_facts = left["facts"]
    right_facts = right["facts"]
    if not left["informative"] and not right["informative"]:
        return {
            "score": None,
            "status": "unavailable",
            "shared": 0,
            "union": 0,
            "different": [],
        }

    union = left_facts | right_facts
    shared = left_facts & right_facts
    score = len(shared) / len(union) if union else None
    different = sorted(union - shared)
    return {
        "score": score,
        "status": _status_for_score(score),
        "shared": len(shared),
        "union": len(union),
        "different": different[:6],
    }



def qualitative_similarity_label(strict_status: str, score) -> str:
    """Return a display category combining strict status and similarity."""
    if strict_status == "unavailable" or score is None:
        return "n/a"
    if strict_status == "missing":
        return "missing"
    if strict_status == "same":
        return "same"
    if score >= 0.95:
        return "near"
    if score >= 0.80:
        return "partial"
    return "diverged"

def build_audit_similarity(reports: list[dict]) -> dict:
    """Build layer-specific Jaccard similarity matrices from run reports."""
    labels = [_run_label(report, index) for index, report in enumerate(reports)]
    layers = []
    for section in _comparison_sections_with_overall(reports):
        per_run = []
        for report in reports:
            facts, informative = _facts_for_section(section, report)
            per_run.append({"facts": facts, "informative": informative})

        rows = []
        scores = []
        for left_index, left in enumerate(per_run):
            row = []
            for right_index, right in enumerate(per_run):
                if left_index == right_index:
                    cell = {
                        "score": 1.0 if left["informative"] else None,
                        "status": "self" if left["informative"] else "unavailable",
                        "shared": len(left["facts"]),
                        "union": len(left["facts"]),
                        "different": [],
                    }
                else:
                    cell = _cell(left, right)
                    if cell["score"] is not None:
                        scores.append(cell["score"])
                row.append(cell)
            rows.append(row)

        layers.append(
            {
                "name": section["section"],
                "rows": rows,
                "summary": {
                    "min": min(scores) if scores else None,
                    "median": median(scores) if scores else None,
                    "facts": sum(item["informative"] for item in per_run),
                },
            }
        )
    return {"runs": labels, "layers": layers}
