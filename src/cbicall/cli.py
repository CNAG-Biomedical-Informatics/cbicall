import json
import os
import platform
import argparse
import gzip
import hashlib
import io
import re
import shutil
import subprocess
import sys
import threading
import time
from contextlib import redirect_stdout
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
from .execution import WorkflowExecutor, EXECUTION_CONTRACT_FILE
from .errors import ParameterValidationError
from .helpmod import usage, parse_args as _parse_args, parse_run_args as _parse_run_args
from .integration_tests import run_integration_tests, run_release_equivalence_test, selected_tests_from_args
from .models import ResolvedConfig, RunSettings
from .resources import validate_resource_catalog
from .workflow_registry import load_workflow_registry
from .goodbye import GoodBye
from .html_reports import render_compare_html, write_run_report_html
from .multiqc import write_compare_multiqc_report, write_multiqc_report
from .report_utils import (
    _aggregate_status,
    _comparison_sections_with_overall,
    _comparison_value,
    _execution_file_map,
    _execution_file_value,
    _format_bytes,
    _format_optional_bytes,
    _inventory_manifest_hash,
    _inventory_total_bytes,
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
            str(workflow.registry_version),
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


def _cromwell_jar_report(jar_path: Path, *, source: str) -> dict:
    java_cmd = os.environ.get("JAVA_CMD", "java")
    report = _command_version(java_cmd, ["-jar", str(jar_path), "--version"])
    report["name"] = "cromwell"
    report["launcher"] = "java-jar"
    report["jar"] = str(jar_path)
    report["source"] = source
    if jar_path.is_file():
        report["jar_sha256"] = _sha256_file(jar_path)
    else:
        report["status"] = "not_found"
        report["error"] = f"Cromwell JAR does not exist: {jar_path}"
    return report


def _cromwell_runtime_report() -> dict:
    jar = os.environ.get("CROMWELL_JAR")
    if jar:
        return _cromwell_jar_report(Path(jar), source="CROMWELL_JAR")
    executable_report = _command_version("cromwell", ["--version"])
    if executable_report.get("status") == "ok":
        return executable_report
    return executable_report

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
    if backend == "cromwell":
        payload["backend"] = _cromwell_runtime_report()
    elif backend in commands:
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


def _collect_execution_contract(project_dir: Path) -> dict:
    path = project_dir / EXECUTION_CONTRACT_FILE
    if not path.is_file():
        return {}
    report = {
        "path": str(path),
        "status": "present",
        "sha256": _sha256_file(path),
    }
    try:
        contract = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        report["status"] = "parse_error"
        report["error"] = str(exc)
        return report

    report["schema_version"] = contract.get("schema_version")
    report["fingerprint"] = contract.get("fingerprint")
    report["workflow_key"] = _nested(contract, "workflow", "key")
    report["backend"] = _nested(contract, "workflow", "backend")
    report["provider"] = _nested(contract, "workflow", "provider")
    report["command_sha256"] = _nested(contract, "command", "sha256")
    report["normalized_command_sha256"] = _nested(contract, "command", "normalized_sha256")
    report["generated_files"] = contract.get("generated_files") or []
    report["source"] = contract
    return report


def _collect_file_inventory(project_dir: Path) -> dict:
    excluded = {"run-report.json", "run-report.html", "cbicall_mqc.yaml"}
    excluded_roots = {
        "01_genomicsdb",
        "work",
        ".nextflow",
        "cbicall_mqc",
        "cromwell-executions",
        "cromwell-logs",
        "cromwell-outputs",
        "cromwell-work",
    }
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
                    "call_fields": parsed.get("call_fields"),
                    "call_sort": parsed.get("call_sort"),
                    "call_records": parsed.get("call_records"),
                    "call_sha256": parsed.get("call_sha256"),
                    "sample_order_fields": parsed.get("sample_order_fields"),
                    "sample_count": parsed.get("sample_count"),
                    "sample_order_sha256": parsed.get("sample_order_sha256"),
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


def _call_level_vcf_hash(path: Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    call_records = []
    sample_names = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                sample_names = line.rstrip("\r\n").split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                continue
            gt_values = "."
            if len(fields) >= 10:
                format_keys = fields[8].split(":")
                gt_index = None
                for index, key in enumerate(format_keys):
                    if key == "GT":
                        gt_index = index
                        break
                sample_gts = []
                for sample in fields[9:]:
                    sample_values = sample.split(":")
                    if gt_index is not None and gt_index < len(sample_values):
                        sample_gts.append(sample_values[gt_index])
                    else:
                        sample_gts.append(".")
                gt_values = ",".join(sample_gts)
            call_fields = [fields[0], fields[1], fields[3], fields[4], fields[6], gt_values]
            call_records.append("\t".join(call_fields))
    digest = hashlib.sha256()
    for record in sorted(call_records):
        digest.update((record + "\n").encode("utf-8"))
    sample_payload = ("\t".join(sample_names) + "\n").encode("utf-8") if sample_names else b""
    return {
        "call_fields": "CHROM,POS,REF,ALT,FILTER,GT_ALL_SAMPLES",
        "call_sort": "LC_ALL=C",
        "call_records": len(call_records),
        "call_sha256": digest.hexdigest(),
        "sample_order_fields": "#CHROM_COLUMNS_10_PLUS",
        "sample_count": len(sample_names),
        "sample_order_sha256": hashlib.sha256(sample_payload).hexdigest(),
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
                        **_call_level_vcf_hash(path),
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



def write_run_report(
    resolved_config: ResolvedConfig,
    arg: dict,
    params: dict,
    *,
    elapsed_seconds: float,
    workflow_log: Path,
    status: str = "success",
    error: Optional[dict] = None,
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
    execution_contract = _collect_execution_contract(project_dir)
    payload = {
        "status": status,
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
            "hostname": resolved_config.hostname,
            "host_threads": resolved_config.host_threads,
            "host_threads_minus_one": resolved_config.host_threads_minus_one,
            "threads": arg.get("threads"),
            "genome": resolved_config.genome,
            "display_genome": resolved_config.display_genome,
            "qc_coverage_region": resolved_config.qc_coverage_region,
            "run_mode": resolved_config.run_mode,
            "output_basename": resolved_config.output_basename,
            "cohort_stage": resolved_config.cohort_stage,
            "interval_shard": resolved_config.interval_shard,
        },
        "inputs": resolved_config.inputs.to_dict(),
        "parameters": {
            "paramfile": arg.get("paramfile"),
            "cleanup_bam": params.get("cleanup_bam", False),
            "qc_coverage_region": resolved_config.qc_coverage_region,
            "snakemake_parameters": resolved_config.snakemake_parameters,
            "nextflow_parameters": resolved_config.nextflow_parameters,
            "cromwell_parameters": resolved_config.cromwell_parameters,
            "nfcore_profile": resolved_config.nfcore_profile,
            "nfcore_parameters": resolved_config.nfcore_parameters,
            "nfcore_singularity_cache_dir": resolved_config.nfcore_singularity_cache_dir,
        },
    }
    if error:
        payload["error"] = dict(error)
    if software_versions:
        payload["software_versions"] = software_versions
    if execution_trace:
        payload["execution_trace"] = execution_trace
    if execution_contract:
        payload["execution_contract"] = execution_contract
    with report_path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2, sort_keys=True)
    return report_path


def _require_verified_resource(resolved_config: ResolvedConfig) -> None:
    bundle = (resolved_config.resources or {}).get("bundle") or {}
    if not bundle:
        return
    resource_type = bundle.get("type")
    if resource_type not in {None, "bundle"}:
        return

    runtime_check = bundle.get("runtime_check") or {}
    status = runtime_check.get("status")
    if status == "verified":
        return

    details = [
        f"Selected bundle resource could not be verified (status={status or 'unknown'}).",
        f"  resource: {bundle.get('key') or '(undef)'}",
    ]
    if runtime_check.get("datadir"):
        details.append(f"  DATADIR: {runtime_check.get('datadir')}")
    if runtime_check.get("source"):
        details.append(f"  source: {runtime_check.get('source')}")
    details.append("Run validate-resources and check the configured resource directory before launching.")
    raise ParameterValidationError("\n".join(details))

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
    _require_verified_resource(resolved_config)

    _section("Parameters OK", GREEN)
    workflow = resolved_config.workflow
    bundle = resolved_config.resources.get("bundle", {})
    _row("Param file", _short_path(args.paramfile))
    _row("Runtime profile", resolved_config.profile)
    _row("Workflow", f"{workflow.backend} -> {workflow.pipeline} -> {workflow.mode}")
    _row("Workflow provider", resolved_config.workflow_provider)
    _row("Software stack", workflow.software_stack)
    _row("Registry ver", workflow.registry_version)
    _row("Genome", resolved_config.genome or "b37")
    if workflow.metadata.get("provider") == "nf-core":
        _row("External workflow", workflow.metadata.get("source"))
        _row("External release", workflow.metadata.get("release"))
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


def _collect_cromwell_wdls(registry: dict, project_root: Path) -> List[Path]:
    workflows = registry.get("workflows") or {}
    cromwell = workflows.get("cromwell") or {}
    base_root = project_root / str(cromwell.get("base_dir", "workflows/cromwell"))
    wdls = []
    for stack_name, stack in (cromwell.get("software_stacks") or {}).items():
        stack_dir = base_root / str(stack_name)
        for pipeline in (stack.get("pipelines") or {}).values():
            for mode in pipeline.values():
                for implementation in (mode.get("registry_versions") or {}).values():
                    script = implementation.get("script") if isinstance(implementation, dict) else implementation
                    if script:
                        wdls.append((stack_dir / str(script)).resolve())
    return sorted(set(wdls))


def _validate_cromwell_wdls_with_womtool(registry: dict, project_root: Path) -> dict:
    wdls = _collect_cromwell_wdls(registry, project_root)
    if not wdls:
        return {"status": "not_applicable", "count": 0, "files": []}
    womtool = os.environ.get("WOMTOOL_JAR")
    if not womtool:
        return {"status": "skipped", "reason": "WOMTOOL_JAR not set", "count": len(wdls), "files": [str(p) for p in wdls]}
    womtool_path = Path(womtool)
    if not womtool_path.is_file():
        raise ParameterValidationError(f"WOMTOOL_JAR does not exist: {womtool}")
    java_cmd = os.environ.get("JAVA_CMD", "java")
    for wdl in wdls:
        proc = subprocess.run(
            [java_cmd, "-jar", str(womtool_path), "validate", str(wdl)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or "womtool validation failed").strip()
            raise ParameterValidationError(f"WDL validation failed for {wdl}: {detail}")
    return {"status": "ok", "count": len(wdls), "files": [str(p) for p in wdls]}


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
    wdl_validation = _validate_cromwell_wdls_with_womtool(registry, root)

    _section("Registry OK", GREEN)
    _row("Registry", _short_path(registry_path))
    _row("Schema", _short_path(schema_path))
    _row("Backends", ", ".join(backends) if backends else "(none)")
    _row("External", ", ".join(external_sources) if external_sources else "(none)")
    if wdl_validation["status"] == "ok":
        _row("WDL syntax", f"ok ({wdl_validation['count']} files)")
    elif wdl_validation["status"] == "skipped":
        _row("WDL syntax", f"skipped ({wdl_validation['reason']})")
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


def _report_json_status(report_path: Path, refresh_requested: bool, refreshed: bool, wrote_json: bool) -> str:
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

    _section("Run Report", GREEN if payload.get("status") == "success" else YELLOW)
    _row("Report", _report_json_status(report_path, refresh_requested, refreshed, wrote_json))
    _row("Status", payload.get("status"))
    _row("Run ID", run.get("run_id"))
    _row("Hostname", run.get("hostname"))
    _row("Elapsed", _format_duration(float(payload.get("elapsed_seconds") or 0)))
    _row("Project", _short_path(run.get("project_dir")))
    _row("HTML", _report_html_status(report_path, html_path))
    if multiqc_path is not None:
        _row("MultiQC", multiqc_path)

    _section("Workflow", BLUE)
    _row("Key", workflow.get("key"))
    _row("Backend", workflow.get("backend"))
    _row("Software stack", workflow.get("software_stack"))
    _row("Pipeline", workflow.get("pipeline"))
    _row("Mode", workflow.get("mode"))
    _row("Registry ver", workflow.get("registry_version"))
    metadata = workflow.get("metadata") or {}
    if metadata.get("provider") == "nf-core":
        _row("External workflow", metadata.get("source"))
        _row("External release", metadata.get("release"))
    _row("Fingerprint", workflow.get("fingerprint"))
    _row("Backend ver", backend.get("version"))

    if execution_contract:
        _section("Execution Contract", BLUE)
        _row("Contract", _short_path(execution_contract.get("path")))
        _row("Fingerprint", execution_contract.get("fingerprint"))
        _row("Command hash", execution_contract.get("normalized_command_sha256"))
        _row("Generated files", len(execution_contract.get("generated_files") or []))

    _section("Resources", CYAN)
    _row("Resource key", bundle.get("key"))
    _row("Resource ver", bundle.get("version"))
    _row("Resource hash", bundle.get("fingerprint"))

    _section("Outputs", WHITE)
    _row("Inventory files", inventory.get("entries"))
    _row("Inventory size", inventory.get("total_bytes_human") or inventory.get("total_bytes"))
    _row("Inventory hash", inventory.get("sha256"))
    _row("VCF hashes", len(vcf_reports))
    _row("Canonical outputs", len(canonical_outputs))
    if execution_trace:
        _row("Trace tasks", execution_trace.get("task_count"))
        _row("Max RSS", execution_trace.get("max_peak_rss_human"))
    if software_versions:
        _row("Software scope", software_versions.get("scope"))
        _row("Software hash", software_versions.get("sha256"))


def _run_report_command(argv: List[str]) -> int:
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
    parser.add_argument("-O", "--overwrite", action="store_true", help="Overwrite files written by --refresh, --html, or --multiqc.")
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
    args = parser.parse_args(argv)

    if args.nocolor or args.json:
        os.environ["ANSI_COLORS_DISABLED"] = "1"
    _refresh_colors()

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
            raise FileExistsError(f"run-report.json would be updated: {report_path}. Use -O/--overwrite to replace it.")
        if refreshed:
            report_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
            wrote_json = True

    html_path = None
    if args.html is not None:
        html_path = report_path.with_suffix(".html") if args.html is True else Path(args.html)
        if html_path.exists() and not args.overwrite:
            raise FileExistsError(f"HTML report already exists: {html_path}. Use -O/--overwrite to replace it.")
        html_path.parent.mkdir(parents=True, exist_ok=True)
        write_run_report_html(report_path, payload, html_path=html_path)

    multiqc_path = None
    if args.multiqc is not None:
        multiqc_path = report_path.parent / "cbicall_mqc" if args.multiqc is True else Path(args.multiqc)
        if multiqc_path.exists() and not args.overwrite:
            raise FileExistsError(f"MultiQC custom-content directory already exists: {multiqc_path}. Use -O/--overwrite to replace it.")
        write_multiqc_report(report_path, payload, output_path=multiqc_path)

    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
    else:
        _print_single_run_report(payload, report_path, html_path, multiqc_path, args.refresh, refreshed, wrote_json)
    return 0



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
    _multi_status("Registry ver", baseline, reports, lambda report: _nested(report, "workflow", "registry_version"))
    _multi_status("External workflow", baseline, reports, lambda report: _nested(report, "workflow", "metadata", "source"))
    _multi_status("External release", baseline, reports, lambda report: _nested(report, "workflow", "metadata", "release"))
    _multi_status("Entrypoint", baseline, reports, lambda report: _nested(report, "workflow", "entrypoint"))
    _multi_status("Workflow hash", baseline, reports, lambda report: _nested(report, "workflow", "fingerprint"))

    print()
    _section("Execution Contract", BLUE)
    _multi_status("Contract hash", baseline, reports, lambda report: _nested(report, "execution_contract", "fingerprint"))
    _multi_status("Command hash", baseline, reports, lambda report: _nested(report, "execution_contract", "normalized_command_sha256"))
    execution_roles = sorted({role for report in reports for role in _execution_file_map(report)})
    if not execution_roles:
        _row("Generated files", "not available")
    for role in execution_roles:
        _multi_status(role, baseline, reports, lambda report, item=role: _execution_file_value(item, report))

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
    _section("All-to-All Matrix", CYAN)
    _row("Runs", len(reports))
    _row("Pairs", (len(reports) * (len(reports) - 1)) // 2)
    _print_run_comparison_legend()

    for section in _comparison_sections_with_overall(reports):
        print()
        _section(section["section"], BLUE)
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
    _section("Legend", YELLOW)
    _row("same", "values or fingerprints match")
    _row("different", "values or fingerprints exist in both runs but differ")
    _row("missing", "available in only some runs")
    _row("note", "audit hint; not treated as a failed reproducibility check")
    _row("not available", "not recorded in the run reports; task/RAM traces require a backend execution trace")


def _print_run_comparison(left: dict, right: dict) -> None:
    _section("Run Comparison", CYAN)
    _row("Run A", _run_label(left))
    _row("Run B", _run_label(right))
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
    _compare_row("Registry ver", _nested(left, "workflow", "registry_version"), _nested(right, "workflow", "registry_version"))
    _compare_row("External workflow", _nested(left, "workflow", "metadata", "source"), _nested(right, "workflow", "metadata", "source"))
    _compare_row("External release", _nested(left, "workflow", "metadata", "release"), _nested(right, "workflow", "metadata", "release"))
    _compare_row("Entrypoint", _nested(left, "workflow", "entrypoint"), _nested(right, "workflow", "entrypoint"))
    _compare_row("Workflow hash", _nested(left, "workflow", "fingerprint"), _nested(right, "workflow", "fingerprint"))

    print()
    _section("Execution Contract", BLUE)
    _compare_row("Contract hash", _nested(left, "execution_contract", "fingerprint"), _nested(right, "execution_contract", "fingerprint"))
    _compare_row("Command hash", _nested(left, "execution_contract", "normalized_command_sha256"), _nested(right, "execution_contract", "normalized_command_sha256"))
    _compare_execution_files(left, right)

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



def _run_compare_runs_command(argv: List[str]) -> int:
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
    parser.add_argument("-nc", "--no-color", dest="nocolor", action="store_true", help="Do not print colors.")
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


def _run_test_command(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="cbicall test",
        description="Run the bundled CBIcall integration tests from examples/input.",
    )
    parser.add_argument("--wes-bash", action="store_true", help="Run the Bash WES integration test.")
    parser.add_argument(
        "--wes-cohort-bash",
        action="store_true",
        help="Run the Bash WES cohort joint-genotyping integration test.",
    )
    parser.add_argument(
        "--wes-cohort-bash-sharded",
        action="store_true",
        help="Run the Bash WES cohort staged shard integration test.",
    )
    parser.add_argument("--wes-bash-gatk35", action="store_true", help="Run the legacy Bash WES GATK 3.5 integration test.")
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
    parser.add_argument(
        "--wes-cromwell",
        action="store_true",
        help="Run the Cromwell WES integration test. Requires CROMWELL_JAR or cromwell on PATH.",
    )
    parser.add_argument("--mit-bash", action="store_true", help="Run the Bash mitochondrial integration test.")
    parser.add_argument(
        "--nf-core-demo",
        action="store_true",
        help="Run the external nf-core/demo integration contract. Requires nextflow on PATH.",
    )
    parser.add_argument(
        "--nf-core-sarek",
        action="store_true",
        help="Run the external nf-core/Sarek integration contract. Requires nextflow on PATH.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all native bundled integration tests. Optional backend tests are skipped if their backend executable is not installed.",
    )
    parser.add_argument(
        "--backend-equivalence",
        dest="backend_equivalence",
        action="store_true",
        help="Run native WES backend-equivalence checks and compare normalized VCF output against Bash.",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use.")
    parser.add_argument("--runtime-profile", dest="profile", default="local", help="CBIcall runtime profile for native workflow tests.")
    parser.add_argument(
        "--keep-external-work",
        action="store_true",
        help="Keep heavy Nextflow work directories for external nf-core tests.",
    )
    args = parser.parse_args(argv)

    if args.threads <= 0:
        parser.error("--threads requires a positive integer")

    selectors = [
        args.wes_bash,
        args.wes_cohort_bash,
        args.wes_cohort_bash_sharded,
        args.wes_bash_gatk35,
        args.wes_snakemake,
        args.wes_nextflow,
        args.wes_cromwell,
        args.mit_bash,
        args.nf_core_demo,
        args.nf_core_sarek,
        args.all,
    ]
    if args.backend_equivalence and any(selectors):
        parser.error("--backend-equivalence cannot be combined with individual test selectors or --all")

    root = _project_root()
    if args.backend_equivalence:
        return run_release_equivalence_test(
            project_root=root,
            threads=args.threads,
            runtime_profile=args.profile,
        )

    selected = selected_tests_from_args(args)
    if not selected:
        parser.error(
            "select at least one test with --wes-bash, --wes-cohort-bash, --wes-cohort-bash-sharded, --wes-bash-gatk35, "
            "--wes-snakemake, --wes-nextflow, --wes-cromwell, --mit-bash, --nf-core-demo, "
            "--nf-core-sarek, --backend-equivalence, or --all"
        )

    return run_integration_tests(
        project_root=root,
        selected=selected,
        threads=args.threads,
        runtime_profile=args.profile,
        skip_missing_optional=args.all,
        keep_external_work=args.keep_external_work,
    )


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
    _require_verified_resource(resolved_config)

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
            "cromwell_parameters": resolved_config.cromwell_parameters,
            "nfcore_profile": resolved_config.nfcore_profile,
            "nfcore_parameters": resolved_config.nfcore_parameters,
            "nfcore_singularity_cache_dir": resolved_config.nfcore_singularity_cache_dir,
            "genome": resolved_config.genome,
            "qc_coverage_region": resolved_config.qc_coverage_region,
            "cleanup_bam": params.get("cleanup_bam", False),
            "output_basename": resolved_config.output_basename,
            "cohort_stage": resolved_config.cohort_stage,
            "interval_shard": resolved_config.interval_shard,
            "run_mode": resolved_config.run_mode,
            "inputs": resolved_config.inputs.to_dict(),
            "workflow": resolved_config.workflow.to_dict(),
        }
    )

    workflow = resolved_config.workflow
    genome = resolved_config.display_genome or resolved_config.genome or "b37"
    if workflow.metadata.get("provider") == "nf-core":
        log_name = f"nf-core_{workflow.pipeline}_{workflow.mode}.log"
    else:
        log_name = f"{workflow.backend}_{workflow.software_stack}_{workflow.pipeline}_{workflow.mode}_{genome}.log"
    workflow_log = Path(resolved_config.project_dir) / log_name
    _row("Workflow", f"{workflow.backend} -> {workflow.pipeline} -> {workflow.mode}")
    print("  This workflow may take a while depending on input size and pipeline.")

    wes = WorkflowExecutor(settings)

    # Run with spinner
    try:
        run_with_spinner(
            wes.run,
            no_spinner=no_spinner,
        )
    except Exception as exc:
        elapsed = time.time() - start_time
        try:
            report_path = write_run_report(
                resolved_config,
                arg,
                params,
                elapsed_seconds=elapsed,
                workflow_log=workflow_log,
                status="failed",
                error={"type": type(exc).__name__, "message": str(exc)},
            )
            html_report_path = write_run_report_html(
                report_path,
                json.loads(report_path.read_text(encoding="utf-8")),
            )
            print()
            _section("Failed", RED)
            _row("Status", "Workflow execution failed")
            _row("Elapsed", _format_duration(elapsed))
            _row("Log", workflow_log)
            _row("Report", report_path)
            _row("HTML", html_report_path)
        except Exception as report_exc:
            print(
                f"WARNING: Could not write failed-run audit report: {report_exc}",
                file=sys.stderr,
            )
        raise

    # END CBICALL
    print()
    _section("Completed", GREEN)
    _row("Status", "Finished successfully")
    elapsed = time.time() - start_time
    _row("Elapsed", _format_duration(elapsed))
    report_path = write_run_report(
        resolved_config,
        arg,
        params,
        elapsed_seconds=elapsed,
        workflow_log=workflow_log,
    )
    report_payload = json.loads(report_path.read_text(encoding="utf-8"))
    html_report_path = write_run_report_html(report_path, report_payload)
    multiqc_path = write_multiqc_report(report_path, report_payload) if arg.get("multiqc") else None
    _row("Log", workflow_log)
    _row("Report", report_path)
    _row("HTML", html_report_path)
    if multiqc_path is not None:
        _row("MultiQC", multiqc_path)
    if workflow.metadata.get("provider") == "nf-core":
        _print_external_output_pointers(report_payload)
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
    if len(sys.argv) > 1 and sys.argv[1] == "report":
        return _run_report_command(sys.argv[2:])
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        return _run_test_command(sys.argv[2:])

    # Backward-compatible legacy run form: cbicall -p params.yaml -t THREADS.
    arg = usage(VERSION)
    return _run_analysis(arg, start_time=start_time, cbicall_path=cbicall_path)
