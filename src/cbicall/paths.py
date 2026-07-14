"""Locate CBIcall runtime assets in source and installed distributions."""

import os
from pathlib import Path
from typing import Optional


RUNTIME_DIR_NAME = "_runtime"


def _is_runtime_root(path: Path) -> bool:
    return (
        (path / "workflows" / "registry" / "cbicall-workflow-registry.yaml").is_file()
        and (path / "resources" / "cbicall-resource-catalog.json").is_file()
    )


def runtime_root(module_file: Optional[str] = None) -> Path:
    """Return the root containing workflows, catalogs, and test assets.

    ``CBICALL_ROOT`` is an explicit developer override. Installed wheels use
    package-local assets, while source checkouts fall back to the repository
    root. The optional module path preserves isolated test-fixture support.
    """
    override = os.environ.get("CBICALL_ROOT")
    if override:
        root = Path(override).expanduser().resolve()
        if not _is_runtime_root(root):
            raise RuntimeError(
                "CBICALL_ROOT does not contain a CBIcall workflow registry and resource catalog: "
                f"{root}"
            )
        return root

    package_runtime = Path(__file__).resolve().parent / RUNTIME_DIR_NAME
    if _is_runtime_root(package_runtime):
        return package_runtime

    if module_file:
        candidate = Path(module_file).resolve().parents[2]
        # A partially assembled source/test tree should still be returned so
        # callers can report the precise missing registry or catalog file.
        if _is_runtime_root(candidate) or (candidate / "workflows").exists() or (candidate / "resources").exists():
            return candidate

    checkout = Path(__file__).resolve().parents[2]
    if _is_runtime_root(checkout):
        return checkout

    raise RuntimeError(
        "CBIcall runtime assets were not found. Reinstall CBIcall or set CBICALL_ROOT "
        "to a valid source checkout."
    )


def workflow_registry_path(root: Optional[Path] = None) -> Path:
    base = Path(root) if root is not None else runtime_root()
    return base / "workflows" / "registry" / "cbicall-workflow-registry.yaml"


def workflow_registry_schema_path(root: Optional[Path] = None) -> Path:
    base = Path(root) if root is not None else runtime_root()
    return base / "workflows" / "schema" / "cbicall-workflow-registry.schema.json"


def resource_catalog_path(root: Optional[Path] = None) -> Path:
    base = Path(root) if root is not None else runtime_root()
    return base / "resources" / "cbicall-resource-catalog.json"


def resource_catalog_schema_path(root: Optional[Path] = None) -> Path:
    base = Path(root) if root is not None else runtime_root()
    return base / "resources" / "cbicall-resource-catalog.schema.json"


def resource_installer_path(root: Optional[Path] = None) -> Path:
    base = Path(root) if root is not None else runtime_root()
    return base / "scripts" / "download_cbicall_bundle.py"
