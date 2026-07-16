"""Construction of execution audit artifacts and run reports.

This module contains filesystem-facing audit logic so the CLI remains focused
on command routing and workflow orchestration.
"""

import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import List, Optional

import yaml

from . import __version__
from .execution import EXECUTION_CONTRACT_FILE
from .models import ResolvedConfig
from .report_utils import _nested
from .runtime_info import runtime_report


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
            "version": resolved_config.version or __version__,
        },
        "runtime": runtime_report(workflow.backend, workflow),
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
