"""Installation diagnostics for the ``cbicall doctor`` command."""

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import List

from . import __version__, console
from .paths import (
    resource_catalog_path,
    runtime_root,
    workflow_registry_path,
    workflow_registry_schema_path,
)
from .resources import (
    RESOURCE_IDENTIFIER,
    RESOURCE_INSTALL_MANIFEST,
    _catalog_entry_for_fingerprint,
    _catalog_fingerprint,
    _load_json_or_text_identifier,
    validate_resource_catalog,
)
from .runtime_info import backend_runtime_report
from .workflow_registry import load_workflow_registry


def _result(status: str, label: str, detail: str, *, required: bool = True) -> dict:
    return {
        "status": status,
        "label": label,
        "detail": detail,
        "required": required,
    }


def _print_result(result: dict) -> None:
    status = result["status"]
    color = {
        "PASS": console.GREEN,
        "WARN": console.YELLOW,
        "FAIL": console.RED,
    }[status]
    print(
        f"{console.BOLD}{color}[{status}]{console.RESET} "
        f"{result['label']:<20} {result['detail']}"
    )


def _data_root() -> tuple:
    configured = os.environ.get("CBICALL_DATA")
    path = Path(configured or "/cbicall-data").expanduser().resolve()
    source = "CBICALL_DATA" if configured else "default"
    return path, source


def _load_installation_manifest(path: Path) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"manifest must contain a JSON object: {path}")
    return payload


def _missing_layout_paths(catalog_entry: dict, data_root: Path) -> list:
    layout = catalog_entry.get("layout") if isinstance(catalog_entry.get("layout"), dict) else {}
    expected_paths = layout.get("expected_top_level") or []
    return [name for name in expected_paths if not (data_root / str(name)).exists()]


def _inspect_identifier_installation(catalog: dict, data_root: Path, identifier_path: Path) -> dict:
    try:
        identifier = _load_json_or_text_identifier(identifier_path)
    except Exception as exc:
        return {
            "status": "invalid",
            "data_root": str(data_root),
            "detail": str(exc),
        }

    resource_key = identifier.get("resource_key")
    catalog_entry = (catalog.get("resources") or {}).get(resource_key)
    if not resource_key or not isinstance(catalog_entry, dict) or catalog_entry.get("type") != "bundle":
        return {
            "status": "invalid",
            "data_root": str(data_root),
            "resource_key": resource_key,
            "detail": "resource identifier does not name a bundle in the local catalog",
        }

    problems = []
    expected_sha256 = (catalog_entry.get("remote_identifier") or {}).get("sha256")
    observed_sha256 = hashlib.sha256(identifier_path.read_bytes()).hexdigest()
    if expected_sha256 and observed_sha256 != expected_sha256:
        problems.append("resource identifier SHA-256 differs from the local catalog")
    missing_paths = _missing_layout_paths(catalog_entry, data_root)
    if missing_paths:
        problems.append("missing " + ", ".join(str(name) for name in missing_paths))

    return {
        "status": "invalid" if problems else "verified",
        "data_root": str(data_root),
        "identifier": str(identifier_path),
        "resource_key": resource_key,
        "detail": "; ".join(problems)
        if problems
        else (
            "resource identifier SHA-256 and layout verified (installation manifest not present)"
            if expected_sha256
            else "resource identity and layout verified (installation manifest not present)"
        ),
    }


def _inspect_installed_bundle(catalog: dict, data_root: Path) -> dict:
    """Verify installation metadata without rehashing the large bundle archive."""
    manifest_path = data_root / RESOURCE_INSTALL_MANIFEST
    if not data_root.is_dir():
        return {
            "status": "missing",
            "data_root": str(data_root),
            "detail": "directory does not exist",
        }
    if not manifest_path.is_file():
        identifier_path = data_root / RESOURCE_IDENTIFIER
        if identifier_path.is_file():
            return _inspect_identifier_installation(catalog, data_root, identifier_path)
        return {
            "status": "manifest_missing",
            "data_root": str(data_root),
            "detail": f"missing {RESOURCE_INSTALL_MANIFEST} and {RESOURCE_IDENTIFIER}",
        }

    try:
        manifest = _load_installation_manifest(manifest_path)
    except (OSError, ValueError) as exc:
        return {
            "status": "invalid",
            "data_root": str(data_root),
            "detail": str(exc),
        }

    manifest_entry = manifest.get("catalog_entry")
    if not isinstance(manifest_entry, dict):
        manifest_entry = {}
    resource_key = manifest.get("resource_key") or manifest_entry.get("key")
    catalog_entry = (catalog.get("resources") or {}).get(resource_key)
    if not resource_key:
        return {
            "status": "invalid",
            "data_root": str(data_root),
            "detail": "manifest does not identify a resource key",
        }
    if not isinstance(catalog_entry, dict) or catalog_entry.get("type") != "bundle":
        return {
            "status": "invalid",
            "data_root": str(data_root),
            "resource_key": resource_key,
            "detail": "installed resource is not a bundle in the local catalog",
        }

    problems = []
    if not manifest_entry:
        problems.append("catalog fingerprint is not recorded")
    else:
        installed_fingerprint = _catalog_fingerprint(
            _catalog_entry_for_fingerprint(manifest_entry)
        )
        catalog_fingerprint = _catalog_fingerprint(
            _catalog_entry_for_fingerprint(catalog_entry)
        )
        if installed_fingerprint != catalog_fingerprint:
            problems.append("catalog fingerprint differs from the installed manifest")

    archive = manifest.get("archive") if isinstance(manifest.get("archive"), dict) else {}
    checksum = manifest.get("checksum") if isinstance(manifest.get("checksum"), dict) else {}
    checksum_verified = archive.get("verified") is True or checksum.get("verified") is True
    if not checksum_verified:
        problems.append("archive checksum is not recorded as verified")

    if manifest.get("extracted") is not True:
        problems.append("bundle is not recorded as extracted")

    missing_paths = _missing_layout_paths(catalog_entry, data_root)
    if missing_paths:
        problems.append("missing " + ", ".join(str(name) for name in missing_paths))

    return {
        "status": "invalid" if problems else "verified",
        "data_root": str(data_root),
        "manifest": str(manifest_path),
        "resource_key": resource_key,
        "detail": "; ".join(problems) if problems else "catalog fingerprint, checksum record, and layout verified",
    }


def _backend_result(backend: str, *, required: bool) -> dict:
    report = backend_runtime_report(backend)
    if report.get("status") == "ok":
        detail = report.get("version") or report.get("detail") or report.get("path") or "available"
        return _result("PASS", backend.capitalize(), str(detail), required=required)

    detail = report.get("detail") or report.get("error") or report.get("status") or "not found"
    if report.get("status") == "not_found":
        detail = "not installed"
    if not required:
        detail = f"{detail} (optional)"
    return _result("FAIL" if required else "WARN", backend.capitalize(), str(detail), required=required)


def run_doctor(argv: List[str]) -> int:
    """Check the installed framework, resource bundle, and backend commands."""
    parser = argparse.ArgumentParser(
        prog="cbicall doctor",
        description="Check whether this CBIcall installation is ready to run workflows.",
    )
    parser.parse_args(argv)
    console.refresh_colors()

    core = [
        _result("PASS", "CBIcall", __version__),
        _result("PASS", "Python", sys.version.split()[0]),
    ]
    registry = None
    catalog = None
    root = None
    try:
        root = runtime_root(__file__)
        registry_path = workflow_registry_path(root)
        registry = load_workflow_registry(
            registry_path,
            workflow_registry_schema_path(root),
        )
        backends = sorted((registry.get("workflows") or {}).keys())
        core.append(
            _result(
                "PASS",
                "Workflow registry",
                f"valid ({len(backends)} backends)",
            )
        )
    except Exception as exc:
        core.append(_result("FAIL", "Workflow registry", str(exc)))

    if root is not None and registry is not None:
        try:
            catalog_path = resource_catalog_path(root)
            summary = validate_resource_catalog(catalog_path, registry)
            catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
            core.append(
                _result(
                    "PASS",
                    "Resource catalog",
                    f"valid ({summary['resources']} resources)",
                )
            )
        except Exception as exc:
            core.append(_result("FAIL", "Resource catalog", str(exc)))
    else:
        core.append(_result("FAIL", "Resource catalog", "not checked because the workflow registry is unavailable"))

    data_root, data_source = _data_root()
    resources = [
        _result(
            "PASS" if data_root.is_dir() else "FAIL",
            "Data directory",
            f"{data_root} ({data_source})",
        )
    ]
    if catalog is None:
        resources.append(_result("FAIL", "Installed bundle", "not checked because the resource catalog is unavailable"))
    else:
        installation = _inspect_installed_bundle(catalog, data_root)
        if installation.get("resource_key"):
            resources.append(
                _result(
                    "PASS",
                    "Installed bundle",
                    installation["resource_key"],
                )
            )
            resources.append(
                _result(
                    "PASS" if installation["status"] == "verified" else "FAIL",
                    "Bundle verification",
                    installation["detail"],
                )
            )
        else:
            resources.append(_result("FAIL", "Installed bundle", installation["detail"]))

    backends = [
        _backend_result("bash", required=True),
        _backend_result("snakemake", required=False),
        _backend_result("nextflow", required=False),
        _backend_result("cromwell", required=False),
    ]

    console.section("CBIcall Doctor", console.CYAN)
    for item in core:
        _print_result(item)
    print()
    console.section("Resources", console.BLUE)
    for item in resources:
        _print_result(item)
    print()
    console.section("Execution Backends", console.BLUE)
    for item in backends:
        _print_result(item)

    all_results = [*core, *resources, *backends]
    failures = [item for item in all_results if item["required"] and item["status"] == "FAIL"]
    warnings = [item for item in all_results if item["status"] == "WARN"]
    print()
    console.section("Summary", console.WHITE)
    if failures:
        _print_result(_result("FAIL", "Status", f"NOT READY ({len(failures)} failed checks)"))
        if not data_root.is_dir():
            console.row("Next step", "set CBICALL_DATA to the installed resource bundle directory")
        return 1

    detail = "READY"
    if warnings:
        detail += f" ({len(warnings)} optional backends unavailable)"
    _print_result(_result("PASS", "Status", detail))
    return 0
