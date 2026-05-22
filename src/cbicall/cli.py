import json
import os
import platform
import argparse
import gzip
import hashlib
import html
import io
import re
import shutil
import subprocess
import sys
import threading
import time
from contextlib import redirect_stdout
from urllib.parse import quote
from pathlib import Path
from typing import List, Optional

import yaml

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
            str(workflow.backend),
            str(workflow.pipeline),
            str(workflow.mode),
            str(workflow.software_stack),
            str(workflow.pipeline_version),
        ]
    )


def _workflow_file_manifest(workflow) -> dict:
    if workflow.metadata.get("provider") == "nf-core":
        source = workflow.metadata.get("source") or workflow.entrypoint
        release = workflow.metadata.get("release")
        payload = f"{source}@{release}".encode("utf-8")
        files = [
            {
                "role": "external:nf-core",
                "path": source,
                "status": "external",
                "sha256": hashlib.sha256(payload).hexdigest(),
                "size_bytes": None,
                "release": release,
            }
        ]
        return {
            "fingerprint": files[0]["sha256"],
            "files": files,
        }

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


def _parse_version_text(text: str) -> Optional[str]:
    match = re.search(r"\b\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*)?\b", text)
    return match.group(0) if match else None


def _command_version(command: str, args: List[str]) -> dict:
    path = shutil.which(command)
    report = {
        "name": command,
        "path": path,
        "command": " ".join([command, *args]),
        "status": "not_found" if path is None else "unknown",
        "version": None,
    }
    if path is None:
        return report

    try:
        completed = subprocess.run(
            [path, *args],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except Exception as exc:
        report["status"] = "error"
        report["error"] = str(exc)
        return report

    output = "\n".join(part for part in [completed.stdout, completed.stderr] if part).strip()
    first_line = next((line.strip() for line in output.splitlines() if line.strip()), "")
    report["status"] = "ok" if completed.returncode == 0 else "nonzero_exit"
    report["returncode"] = completed.returncode
    report["version"] = _parse_version_text(output)
    if first_line:
        report["detail"] = first_line[:200]
    return report


def _architecture_key() -> str:
    machine = platform.machine().lower()
    if machine in {"x86_64", "amd64"}:
        return "amd64"
    if machine in {"aarch64", "arm64"}:
        return "aarch64"
    return machine


def _source_bash_env_variable(env_file: str, variable: str) -> dict:
    report = {
        "source": str(env_file),
        "variable": variable,
        "status": "unknown",
        "path": None,
    }
    path = Path(str(env_file))
    if not path.is_file():
        report["status"] = "source_missing"
        return report
    try:
        completed = subprocess.run(
            [
                "bash",
                "-c",
                'source "$1" >/dev/null 2>&1; printf "%s" "${!2:-}"',
                "bash",
                str(path),
                variable,
            ],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except Exception as exc:
        report["status"] = "error"
        report["error"] = str(exc)
        return report

    value = completed.stdout.strip()
    report["returncode"] = completed.returncode
    if completed.returncode != 0:
        report["status"] = "nonzero_exit"
        if completed.stderr:
            report["detail"] = completed.stderr.strip().splitlines()[0][:200]
        return report
    if not value:
        report["status"] = "not_set"
        return report
    report["status"] = "ok"
    report["path"] = value
    return report


def _configured_native_java_report(workflow) -> Optional[dict]:
    if workflow is None or workflow.metadata.get("provider") == "nf-core":
        return None

    java_path = None
    source = None
    status = "not_configured"
    arch = _architecture_key()

    if workflow.backend == "bash":
        env_file = workflow.helpers.get("env")
        if not env_file:
            return None
        env_report = _source_bash_env_variable(env_file, "JAVA8")
        java_path = env_report.get("path")
        source = f"{env_file}:JAVA8"
        status = env_report.get("status") or status
        if not java_path:
            env_report["name"] = "configured_java"
            return env_report
    elif workflow.backend in {"snakemake", "nextflow"}:
        config_file = workflow.config_file
        if not config_file:
            return None
        config_path = Path(str(config_file))
        source = f"{config_file}:java.{arch}"
        if not config_path.is_file():
            return {"name": "configured_java", "source": source, "status": "source_missing", "path": None, "version": None}
        try:
            config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        except Exception as exc:
            return {"name": "configured_java", "source": source, "status": "error", "path": None, "version": None, "error": str(exc)}
        java_cfg = config.get("java") or {}
        java_path = java_cfg.get(arch)
        if not java_path:
            return {"name": "configured_java", "source": source, "status": "not_set", "path": None, "version": None}
        status = "ok"
    else:
        return None

    version_report = _command_version(str(java_path), ["-version"])
    version_report["name"] = "configured_java"
    version_report["source"] = source
    if status != "ok" and version_report.get("status") == "not_found":
        version_report["status"] = status
    return version_report


def _runtime_report(backend: str, workflow=None) -> dict:
    commands = {
        "bash": ("bash", ["--version"]),
        "snakemake": ("snakemake", ["--version"]),
        "nextflow": ("nextflow", ["-version"]),
    }
    payload = {
        "python": {
            "executable": sys.executable,
            "version": sys.version.split()[0],
        },
        "java": _command_version("java", ["-version"]),
    }
    if backend in commands:
        command, args = commands[backend]
        payload["backend"] = _command_version(command, args)
    else:
        payload["backend"] = {
            "name": backend,
            "path": None,
            "status": "unsupported",
            "version": None,
        }
    configured_java = _configured_native_java_report(workflow)
    if configured_java:
        payload["configured_java"] = configured_java
    return payload


def _parse_key_value_file(path: Path) -> dict:
    data = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip().lower()] = value.strip()
    return data


def _collect_file_inventory(project_dir: Path) -> dict:
    excluded = {"run-report.json", "run-report.html"}
    excluded_roots = {"work", ".nextflow"}
    paths = []
    total_bytes = 0
    largest_files = []
    directories = {}
    if project_dir.is_dir():
        for path in project_dir.rglob("*"):
            if not path.is_file():
                continue
            rel_path = path.relative_to(project_dir).as_posix()
            rel_parts = path.relative_to(project_dir).parts
            if rel_path in excluded:
                continue
            if rel_parts and rel_parts[0] in excluded_roots:
                continue
            if path.name.startswith(".nextflow"):
                continue
            size = path.stat().st_size
            paths.append(rel_path)
            total_bytes += size
            group = rel_parts[0] if len(rel_parts) > 1 else "."
            directories.setdefault(group, {"group": group, "count": 0, "bytes": 0})
            directories[group]["count"] += 1
            directories[group]["bytes"] += size
            largest_files.append({"path": rel_path, "bytes": size})
    paths.sort()
    largest_files = sorted(largest_files, key=lambda item: (-item["bytes"], item["path"]))[:10]
    directory_summary = sorted(directories.values(), key=lambda item: (-item["count"], item["group"]))
    manifest_text = "".join(f"{path}\n" for path in paths).encode("utf-8")
    return {
        "algorithm": "sha256",
        "scope": "run directory relative file paths",
        "entries": len(paths),
        "total_bytes": total_bytes,
        "directories": directory_summary,
        "largest_files": largest_files,
        "sha256": hashlib.sha256(manifest_text).hexdigest(),
        "excluded": sorted(excluded | excluded_roots | {".nextflow*"}),
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


def _normalized_vcf_hash(path: Path) -> dict:
    digest = hashlib.sha256()
    opener = gzip.open if path.suffix == ".gz" else open
    records = 0
    with opener(path, "rb") as handle:
        for line in handle:
            if line.startswith(b"#"):
                continue
            line = line.rstrip(b"\r\n")
            if not line:
                continue
            digest.update(line)
            digest.update(b"\n")
            records += 1
    return {
        "algorithm": "sha256",
        "raw_sha256": _sha256_file(path),
        "normalized_pattern": "^#",
        "normalized_sort": "none",
        "normalized_records": records,
        "normalized_sha256": digest.hexdigest(),
    }


def _canonical_output_reports(workflow, workflow_output_dir: Path) -> tuple:
    reports = []
    vcf_reports = []
    for item in workflow.metadata.get("canonical_outputs", []):
        name = str(item.get("name"))
        output_type = str(item.get("type"))
        pattern = str(item.get("pattern"))
        matches = sorted(workflow_output_dir.glob(pattern))
        report = {
            "name": name,
            "type": output_type,
            "pattern": pattern,
            "root": str(workflow_output_dir),
            "matches": [str(path) for path in matches],
            "status": "present" if matches else "missing",
        }
        reports.append(report)
        if output_type == "vcf":
            for path in matches:
                vcf_reports.append(
                    {
                        "path": str(path),
                        "file": path.name,
                        "source": "registry_canonical_output",
                        "name": name,
                        "pattern": pattern,
                        **_normalized_vcf_hash(path),
                    }
                )
    return reports, vcf_reports


def _first_glob(root: Path, pattern: str) -> Optional[str]:
    matches = sorted(root.glob(pattern))
    return str(matches[0]) if matches else None


def _external_nextflow_summary(workflow_output_dir: Path, canonical_outputs: List[dict]) -> dict:
    summary = {}

    pipeline_info = workflow_output_dir / "pipeline_info"
    if pipeline_info.is_dir():
        info = {"dir": str(pipeline_info)}
        for key, pattern in {
            "params": "params_*.json",
            "trace": "execution_trace_*.txt",
            "report": "execution_report_*.html",
            "timeline": "execution_timeline_*.html",
            "dag": "pipeline_dag_*.html",
            "manifest": "manifest_*.json",
            "software_versions": "*software*versions*.yml",
        }.items():
            match = _first_glob(pipeline_info, pattern)
            if match:
                info[key] = match
        summary["pipeline_info"] = info

    multiqc_dir = workflow_output_dir / "multiqc"
    if multiqc_dir.is_dir():
        multiqc = {"dir": str(multiqc_dir)}
        report = multiqc_dir / "multiqc_report.html"
        data_dir = multiqc_dir / "multiqc_data"
        if report.is_file():
            multiqc["report"] = str(report)
        if data_dir.is_dir():
            multiqc["data_dir"] = str(data_dir)
        summary["multiqc"] = multiqc

    if canonical_outputs:
        summary["canonical_outputs"] = canonical_outputs

    return summary


def _external_nextflow_outputs(resolved_config: ResolvedConfig, project_dir: Path) -> dict:
    workflow = resolved_config.workflow
    if workflow.metadata.get("provider") != "nf-core":
        return {}

    params_file = project_dir / "cbicall_external_nextflow.params.yaml"
    config_file = project_dir / "cbicall_external_nextflow.config"
    output_dir = project_dir / workflow.metadata.get("default_outdir", workflow.pipeline)
    payload = {
        "workflow_output_dir": str(output_dir),
        "params_file": str(params_file),
        "config_file": str(config_file),
    }
    canonical_outputs, vcf_hash_reports = _canonical_output_reports(workflow, output_dir)
    payload["external_summary"] = _external_nextflow_summary(output_dir, canonical_outputs)
    if canonical_outputs:
        payload["canonical_outputs"] = canonical_outputs
    if vcf_hash_reports:
        payload["vcf_hash_reports"] = vcf_hash_reports
    if params_file.is_file():
        payload["params_sha256"] = _sha256_file(params_file)
    if config_file.is_file():
        payload["config_sha256"] = _sha256_file(config_file)
    return payload


def _parse_size_to_bytes(value) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text == "-":
        return None
    match = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s*([KMGTPE]?B)?$", text, re.IGNORECASE)
    if not match:
        return None
    number = float(match.group(1))
    unit = (match.group(2) or "B").upper()
    factors = {
        "B": 1,
        "KB": 1024,
        "MB": 1024 ** 2,
        "GB": 1024 ** 3,
        "TB": 1024 ** 4,
        "PB": 1024 ** 5,
        "EB": 1024 ** 6,
    }
    return int(number * factors.get(unit, 1))


def _trace_max(rows: List[dict], field: str) -> dict:
    best = None
    for row in rows:
        value = _parse_size_to_bytes(row.get(field))
        if value is None:
            continue
        if best is None or value > best["bytes"]:
            best = {
                "bytes": value,
                "value": row.get(field),
                "task": row.get("name") or row.get("process") or row.get("task_id"),
                "hash": row.get("hash"),
                "status": row.get("status"),
            }
    return best or {}


def _collect_execution_trace_summary(output_report: dict) -> dict:
    trace_path = _nested(output_report, "external_summary", "pipeline_info", "trace")
    if not trace_path:
        return {}
    path = Path(str(trace_path))
    if not path.is_file():
        return {}

    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError:
        return {"source": str(path), "status": "read_error"}
    if not lines:
        return {"source": str(path), "status": "empty", "sha256": _sha256_file(path)}

    header = lines[0].split("\t")
    rows = []
    status_counts = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        values = line.split("\t")
        row = {key: values[idx] if idx < len(values) else "" for idx, key in enumerate(header)}
        rows.append(row)
        status = row.get("status") or "unknown"
        status_counts[status] = status_counts.get(status, 0) + 1

    return {
        "source": str(path),
        "sha256": _sha256_file(path),
        "status": "parsed",
        "tasks": len(rows),
        "status_counts": status_counts,
        "max_peak_rss": _trace_max(rows, "peak_rss"),
        "max_peak_vmem": _trace_max(rows, "peak_vmem"),
    }


def _dict_sha256(value: dict) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _collect_external_software_versions(output_report: dict) -> dict:
    software_path = _nested(output_report, "external_summary", "pipeline_info", "software_versions")
    if not software_path:
        return {}

    path = Path(str(software_path))
    if not path.is_file():
        return {}

    report = {
        "source": str(path),
        "sha256": _sha256_file(path),
        "format": "yaml",
        "scope": "workflow_reported",
    }
    try:
        parsed = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception as exc:
        report["status"] = "parse_error"
        report["error"] = str(exc)
        return report

    report["status"] = "parsed"
    report["entries"] = parsed if isinstance(parsed, dict) else {"value": parsed}
    return report


def _collect_declared_resource_software_versions(resources: dict) -> dict:
    bundle = (resources or {}).get("bundle") or {}
    catalog_path = bundle.get("catalog")
    resource_key = bundle.get("key")
    if not catalog_path or not resource_key:
        return {}

    path = Path(str(catalog_path))
    if not path.is_file():
        return {}

    try:
        catalog = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return {
            "source": str(path),
            "status": "parse_error",
            "error": str(exc),
            "scope": "resource_declared",
        }

    entry = (catalog.get("resources") or {}).get(str(resource_key)) or {}
    tools = entry.get("tools") or {}
    if not isinstance(tools, dict) or not tools:
        return {}

    return {
        "source": str(path),
        "sha256": _dict_sha256(tools),
        "format": "resource_catalog",
        "scope": "resource_declared",
        "status": "parsed",
        "entries": tools,
    }


def _collect_software_versions(output_report: dict, resources: dict) -> dict:
    return _collect_external_software_versions(output_report) or _collect_declared_resource_software_versions(resources)


def _print_external_output_pointers(report: dict) -> None:
    summary = _nested(report, "outputs", "external_summary") or {}
    if not summary:
        return

    print()
    _section("nf-core outputs", BLUE)
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

    output_dashboard = _output_dashboard_html(payload, base_dir)

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
            [
                ("Backend", _nested(payload, "workflow", "backend")),
                ("Pipeline", _nested(payload, "workflow", "pipeline")),
                ("Mode", _nested(payload, "workflow", "mode")),
                ("Genome", _nested(payload, "run", "display_genome") or _nested(payload, "run", "genome")),
                ("Software stack", _nested(payload, "workflow", "software_stack")),
                ("Pipeline version", _nested(payload, "workflow", "pipeline_version")),
                ("Runtime profile", payload.get("profile")),
                ("Threads", _nested(payload, "run", "threads")),
            ],
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


def _write_run_report_html(report_path: Path, payload: dict) -> Path:
    html_path = report_path.with_suffix(".html")
    html_path.write_text(_run_report_html(payload), encoding="utf-8")
    return html_path


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
    output_fingerprints = _collect_output_fingerprints(project_dir)
    external_outputs = _external_nextflow_outputs(resolved_config, project_dir)
    if external_outputs.get("vcf_hash_reports"):
        output_fingerprints["vcf_hash_reports"].extend(external_outputs.pop("vcf_hash_reports"))
    output_fingerprints.update(external_outputs)
    software_versions = _collect_software_versions(output_fingerprints, resolved_config.resources)
    execution_trace = _collect_execution_trace_summary(output_fingerprints)
    payload = {
        "status": "success",
        "elapsed_seconds": round(elapsed_seconds, 3),
        "workflow_log": str(workflow_log),
        "command_trace": "Detailed workflow stdout/stderr is recorded in the workflow log.",
        "framework": {
            "name": "CBIcall",
            "version": resolved_config.version or VERSION,
        },
        "runtime": _runtime_report(workflow.backend, workflow),
        "profile": resolved_config.profile,
        "workflow": _workflow_report(workflow),
        "resources": resolved_config.resources,
        "outputs": output_fingerprints,
        "run": {
            "run_id": resolved_config.run_id,
            "project_dir": resolved_config.project_dir,
            "threads": arg.get("threads"),
            "genome": resolved_config.genome,
            "display_genome": resolved_config.display_genome,
            "run_mode": resolved_config.run_mode,
        },
        "inputs": resolved_config.inputs.to_dict(),
        "parameters": {
            "paramfile": arg.get("paramfile"),
            "cleanup_bam": params.get("cleanup_bam", False),
            "snakemake_parameters": resolved_config.snakemake_parameters,
            "nextflow_parameters": resolved_config.nextflow_parameters,
            "nfcore_profile": resolved_config.nfcore_profile,
            "nfcore_parameters": resolved_config.nfcore_parameters,
            "nfcore_singularity_cache_dir": resolved_config.nfcore_singularity_cache_dir,
        },
    }
    if software_versions:
        payload["software_versions"] = software_versions
    if execution_trace:
        payload["execution_trace"] = execution_trace
    with report_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2, sort_keys=True)
    return report_path


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
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    params = config_mod.read_param_file(args.paramfile)
    params = _apply_cli_runtime_overrides(params, vars(args))
    resolved_config = ResolvedConfig.from_mapping({**config_mod.set_config_values(params), "version": VERSION})

    _section("Parameters OK", GREEN)
    workflow = resolved_config.workflow
    bundle = resolved_config.resources.get("bundle", {})
    _row("Param file", _short_path(args.paramfile))
    _row("Runtime profile", resolved_config.profile)
    _row("Workflow", f"{workflow.backend} -> {workflow.pipeline} -> {workflow.mode}")
    _row("Workflow provider", resolved_config.workflow_provider)
    _row("Software stack", workflow.software_stack)
    _row("Pipeline ver", workflow.pipeline_version)
    _row("Genome", resolved_config.genome or "b37")
    if workflow.metadata.get("provider") == "nf-core":
        _row("Source", workflow.metadata.get("source"))
        _row("Release", workflow.metadata.get("release"))
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
    backends = sorted((registry.get("workflows") or {}).keys())
    external_sources = _collect_external_sources(registry)

    _section("Registry OK", GREEN)
    _row("Registry", _short_path(registry_path))
    _row("Schema", _short_path(schema_path))
    _row("Backends", ", ".join(backends) if backends else "(none)")
    _row("External", ", ".join(external_sources) if external_sources else "(none)")
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


def _refresh_report_output_fields(report_path: Path, payload: dict) -> bool:
    """Refresh audit fields that can be derived from files in an existing run dir."""
    project_dir = report_path.parent
    if not project_dir.is_dir():
        return False

    fresh = _collect_output_fingerprints(project_dir)
    outputs = payload.setdefault("outputs", {})
    changed = False

    inventory = outputs.setdefault("file_inventory", {})
    fresh_inventory = fresh.get("file_inventory") or {}
    if int(fresh_inventory.get("entries") or 0) > 0:
        for key in ["directories", "largest_files", "total_bytes", "entries", "sha256", "algorithm", "scope", "excluded", "paths"]:
            if fresh_inventory.get(key) is not None and inventory.get(key) != fresh_inventory.get(key):
                inventory[key] = fresh_inventory[key]
                changed = True

    fresh_vcfs = fresh.get("vcf_hash_reports") or []
    if fresh_vcfs and outputs.get("vcf_hash_reports") != fresh_vcfs:
        outputs["vcf_hash_reports"] = fresh_vcfs
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


def _run_render_report_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall render-report",
        description="Regenerate run-report.html from an existing CBIcall run-report.json.",
    )
    parser.add_argument("run", help="Run directory or run-report.json file.")
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

    report_path = _run_report_path(args.run)
    try:
        payload = json.loads(report_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid run report JSON: {report_path}") from exc

    refreshed = _refresh_report_output_fields(report_path, payload)
    if refreshed:
        report_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    html_path = _write_run_report_html(report_path, payload)
    _section("Report Rendered", GREEN)
    _row("Report", report_path)
    _row("Refreshed", "outputs" if refreshed else "no")
    _row("HTML", html_path)
    return 0


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


def _short_hash(value: str) -> str:
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
        if details:
            _row(role, f"{status}: {', '.join(details)}")
        else:
            _row(role, f"same: {_comparison_value((left_entry.get('path'), left_entry.get('sha256')))}")


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
        _compare_row(
            _short_path(key),
            left_item.get("normalized_sha256"),
            right_item.get("normalized_sha256"),
        )


def _inventory_total_bytes(report: dict):
    return _nested(report, "outputs", "file_inventory", "total_bytes")


def _inventory_manifest_hash(report: dict):
    return _nested(report, "outputs", "file_inventory", "sha256")


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

    _section("Run Matrix", CYAN)
    _row("Runs", len(reports))
    _row("Baseline", _run_label(baseline))
    _row("Compared", ", ".join(_run_label(report) for report in compared))
    _print_run_comparison_legend()

    print()
    _section("Framework", BLUE)
    _multi_status("CBIcall ver", baseline, reports, lambda report: _nested(report, "framework", "version"))
    _multi_status("Python ver", baseline, reports, lambda report: _nested(report, "runtime", "python", "version"))
    _multi_status("Java ver", baseline, reports, lambda report: _nested(report, "runtime", "java", "version"))
    _multi_status("Configured Java", baseline, reports, lambda report: _nested(report, "runtime", "configured_java", "version"))
    _multi_status("Backend ver", baseline, reports, lambda report: _nested(report, "runtime", "backend", "version"))

    print()
    _section("Execution", BLUE)
    _multi_status("Task count (trace)", baseline, reports, lambda report: _nested(report, "execution_trace", "tasks"))
    _multi_status("Max peak RSS (trace)", baseline, reports, lambda report: _nested(report, "execution_trace", "max_peak_rss", "bytes"))
    _multi_status("Max peak VMEM (trace)", baseline, reports, lambda report: _nested(report, "execution_trace", "max_peak_vmem", "bytes"))

    print()
    _section("Pipeline", BLUE)
    _multi_status("Workflow key", baseline, reports, lambda report: _nested(report, "workflow", "key"))
    _multi_status("Pipeline ver", baseline, reports, lambda report: _nested(report, "workflow", "pipeline_version"))
    _multi_status("Entrypoint", baseline, reports, lambda report: _nested(report, "workflow", "entrypoint"))
    _multi_status("Workflow hash", baseline, reports, lambda report: _nested(report, "workflow", "fingerprint"))

    print()
    _section("Software", BLUE)
    _multi_status("Software versions", baseline, reports, lambda report: _nested(report, "software_versions", "sha256"))

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
    _multi_inventory_size_status(baseline, reports)
    _multi_status("File inventory", baseline, reports, lambda report: _nested(report, "outputs", "file_inventory", "sha256"))
    keys = sorted({key for report in reports for key in _vcf_hash_map(report)})
    if not keys:
        _row("VCF hashes", "not available")
    for key in keys:
        _multi_status(_short_path(key), baseline, reports, lambda report, item=key: _multi_vcf_hash_value(item, report))


def _print_run_comparison_legend() -> None:
    print()
    _section("Legend", YELLOW)
    _row("same", "values or fingerprints match")
    _row("different", "values or fingerprints exist in both runs but differ")
    _row("missing", "available in only some runs")
    _row("note", "audit hint; not treated as a failed reproducibility check")
    _row("not available", "not recorded in the run reports; task/RAM traces require a backend execution trace")


def _print_run_comparison(left: dict, right: dict) -> None:
    _section("Run Comparison", CYAN)
    _row("Run A", _short_path(left["_report_path"]))
    _row("Run B", _short_path(right["_report_path"]))
    _print_run_comparison_legend()

    print()
    _section("Framework", BLUE)
    _compare_row("CBIcall ver", _nested(left, "framework", "version"), _nested(right, "framework", "version"))
    _compare_row("Python ver", _nested(left, "runtime", "python", "version"), _nested(right, "runtime", "python", "version"))
    _compare_row("Java ver", _nested(left, "runtime", "java", "version"), _nested(right, "runtime", "java", "version"))
    _compare_row("Configured Java", _nested(left, "runtime", "configured_java", "version"), _nested(right, "runtime", "configured_java", "version"))
    _compare_row("Backend ver", _nested(left, "runtime", "backend", "version"), _nested(right, "runtime", "backend", "version"))

    print()
    _section("Execution", BLUE)
    _compare_row("Task count (trace)", _nested(left, "execution_trace", "tasks"), _nested(right, "execution_trace", "tasks"))
    _compare_row("Max peak RSS (trace)", _nested(left, "execution_trace", "max_peak_rss", "bytes"), _nested(right, "execution_trace", "max_peak_rss", "bytes"))
    _compare_row("Max peak VMEM (trace)", _nested(left, "execution_trace", "max_peak_vmem", "bytes"), _nested(right, "execution_trace", "max_peak_vmem", "bytes"))

    print()
    _section("Pipeline", BLUE)
    _compare_row("Workflow key", _nested(left, "workflow", "key"), _nested(right, "workflow", "key"))
    _compare_row("Pipeline ver", _nested(left, "workflow", "pipeline_version"), _nested(right, "workflow", "pipeline_version"))
    _compare_row("Entrypoint", _nested(left, "workflow", "entrypoint"), _nested(right, "workflow", "entrypoint"))
    _compare_row("Workflow hash", _nested(left, "workflow", "fingerprint"), _nested(right, "workflow", "fingerprint"))

    print()
    _section("Software", BLUE)
    _compare_row("Software versions", _nested(left, "software_versions", "sha256"), _nested(right, "software_versions", "sha256"))

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
    _compare_inventory_size(left, right)
    _compare_row("File inventory", _nested(left, "outputs", "file_inventory", "sha256"), _nested(right, "outputs", "file_inventory", "sha256"))
    _compare_output_hashes(left, right)


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
      <input type="radio" name="compare-tabs" id="tab-details">
      <input type="radio" name="compare-tabs" id="tab-raw">
      <div class="tab-labels" aria-label="Report views">
        <label for="tab-overview">Overview</label>
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

def _run_compare_runs_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall compare-runs",
        description="Compare CBIcall run-report.json files or run directories.",
    )
    parser.add_argument("runs", nargs="+", help="Run directories or run-report.json files. Provide two or more.")
    parser.add_argument("-o", "--output", help="Write the text comparison report to this file.")
    parser.add_argument("--html", help="Write the HTML comparison report to this file. Defaults to <output>.html or compare-runs.html.")
    parser.add_argument("--no-html", action="store_true", help="Do not write the default HTML comparison report.")
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if len(args.runs) < 2:
        parser.error("compare-runs requires at least two run directories or run-report.json files")
    if args.no_html and args.html:
        parser.error("--html and --no-html cannot be used together")

    html_output = None
    if not args.no_html:
        if args.html:
            html_output = Path(args.html)
        elif args.output:
            html_output = Path(args.output).with_suffix(".html")
        else:
            html_output = Path("compare-runs.html")

    if args.nocolor or args.output or html_output:
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
    if html_output:
        html_output.parent.mkdir(parents=True, exist_ok=True)
        html_output.write_text(_render_compare_html(report_text), encoding="utf-8")
        _row("HTML", html_output)
    return 0


def _run_test_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall test",
        description="Run the bundled CBIcall integration tests from examples/input.",
    )
    parser.add_argument("--wes-bash", action="store_true", help="Run the Bash WES integration test.")
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
    parser.add_argument("--mit-bash", action="store_true", help="Run the Bash mitochondrial integration test.")
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all bundled integration tests. Optional engine tests are skipped if their engine is not installed.",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use.")
    parser.add_argument("--runtime-profile", dest="profile", default="local", help="CBIcall runtime profile for native workflow tests.")
    args = parser.parse_args(argv)

    if args.threads <= 0:
        parser.error("--threads requires a positive integer")

    selected = []
    if args.all or args.wes_bash:
        selected.append("--wes-bash")
    if args.all or args.wes_snakemake:
        selected.append("--wes-snakemake")
    if args.all or args.wes_nextflow:
        selected.append("--wes-nextflow")
    if args.all or args.mit_bash:
        selected.append("--mit-bash")
    if not selected:
        parser.error("select at least one test with --wes-bash, --wes-snakemake, --wes-nextflow, --mit-bash, or --all")

    root = _project_root()
    script = root / "examples" / "input" / "run_tests.sh"
    if not script.is_file():
        raise FileNotFoundError(f"Integration test launcher not found: {script}")

    env = os.environ.copy()
    env["CBICALL"] = str(root / "bin" / "cbicall")
    env["THREADS"] = str(args.threads)
    env["CBICALL_RUNTIME_PROFILE"] = str(args.profile)
    if args.all and not args.wes_snakemake:
        env["CBICALL_TEST_SKIP_MISSING_OPTIONAL"] = "1"
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
    params = _apply_cli_runtime_overrides(params, arg)
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
            "snakemake_parameters": resolved_config.snakemake_parameters,
            "nextflow_parameters": resolved_config.nextflow_parameters,
            "nfcore_profile": resolved_config.nfcore_profile,
            "nfcore_parameters": resolved_config.nfcore_parameters,
            "nfcore_singularity_cache_dir": resolved_config.nfcore_singularity_cache_dir,
            "genome": resolved_config.genome,
            "cleanup_bam": params.get("cleanup_bam", False),
            "run_mode": resolved_config.run_mode,
            "inputs": resolved_config.inputs.to_dict(),
            "workflow": resolved_config.workflow.to_dict(),
        }
    )

    workflow = resolved_config.workflow
    _row("Workflow", f"{workflow.backend} -> {workflow.pipeline} -> {workflow.mode}")
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
    genome = resolved_config.display_genome or resolved_config.genome or "b37"
    if workflow.metadata.get("provider") == "nf-core":
        log_name = f"nf-core_{workflow.pipeline}_{workflow.mode}.log"
    else:
        log_name = f"{workflow.backend}_{workflow.software_stack}_{workflow.pipeline}_{workflow.mode}_{genome}.log"
    workflow_log = Path(resolved_config.project_dir) / log_name
    report_path = write_run_report(
        resolved_config,
        arg,
        params,
        elapsed_seconds=elapsed,
        workflow_log=workflow_log,
    )
    html_report_path = _write_run_report_html(report_path, json.loads(report_path.read_text(encoding="utf-8")))
    _row("Log", workflow_log)
    _row("Report", report_path)
    _row("HTML", html_report_path)
    if workflow.metadata.get("provider") == "nf-core":
        report_data = json.loads(report_path.read_text(encoding="utf-8"))
        _print_external_output_pointers(report_data)
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
    if len(sys.argv) > 1 and sys.argv[1] == "validate-parameters":
        return _run_validate_parameters_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "validate-registry":
        return _run_validate_registry_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "validate-resources":
        return _run_validate_resources_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "compare-runs":
        return _run_compare_runs_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "render-report":
        return _run_render_report_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        return _run_test_command(sys.argv[2:])

    # Backward-compatible legacy run form: cbicall -p params.yaml -t THREADS.
    arg = usage(VERSION)
    return _run_analysis(arg, start_time=start_time, cbicall_path=cbicall_path)
