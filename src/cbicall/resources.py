import hashlib
import json
import os
import re
from pathlib import Path

import yaml

from .errors import ParameterValidationError
from .models import WorkflowSpec


RESOURCE_INSTALL_MANIFEST = "cbicall-resource-installation.json"
RESOURCE_IDENTIFIER = "cbicall-bundle-id.json"


def _catalog_fingerprint(entry: dict) -> str:
    payload = json.dumps(entry, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _catalog_entry_for_fingerprint(entry: dict) -> dict:
    return {
        key: value
        for key, value in entry.items()
        if key != "key" and not str(key).startswith("_")
    }


def _workflow_key(cfg_in: dict) -> str:
    return "/".join(
        [
            cfg_in["workflow_engine"],
            cfg_in["pipeline"],
            cfg_in["mode"],
            cfg_in["gatk_version"],
        ]
    )


def build_resource_bundle_metadata(cfg_in: dict, project_root: Path) -> dict:
    catalog_path = project_root / "resources" / "cbicall-resource-catalog.json"
    bundle_key = cfg_in["resource_bundle"]

    if not catalog_path.is_file():
        entry = {
            "bundle_id": "cbicall-germline-resources",
            "bundle_version": bundle_key,
            "status": "unregistered",
            "compatible_workflows": [],
        }
        return {
            "key": bundle_key,
            "bundle_id": entry["bundle_id"],
            "bundle_version": entry["bundle_version"],
            "catalog": None,
            "status": entry["status"],
            "fingerprint": _catalog_fingerprint(_catalog_entry_for_fingerprint(entry)),
            "compatible": None,
            "workflow_key": _workflow_key(cfg_in),
            "runtime_check": {"status": "catalog_not_found"},
        }

    try:
        catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ParameterValidationError(f"Invalid resource catalog JSON: {catalog_path}") from exc

    bundles = catalog.get("bundles", {}) if isinstance(catalog, dict) else {}
    entry = bundles.get(bundle_key)
    if entry is None:
        raise ParameterValidationError(
            f"Resource bundle '{bundle_key}' is not defined in {catalog_path}"
        )

    workflow_key = _workflow_key(cfg_in)
    compatible_workflows = entry.get("compatible_workflows", [])
    compatible = workflow_key in compatible_workflows if compatible_workflows else None
    if compatible_workflows and not compatible:
        raise ParameterValidationError(
            f"Resource bundle '{bundle_key}' is not declared compatible "
            f"with workflow '{workflow_key}'"
        )

    return {
        "key": bundle_key,
        "bundle_id": entry.get("bundle_id", "cbicall-germline-resources"),
        "bundle_version": str(entry.get("bundle_version", bundle_key)),
        "catalog": str(catalog_path),
        "status": entry.get("status"),
        "fingerprint": _catalog_fingerprint(_catalog_entry_for_fingerprint(entry)),
        "compatible": compatible,
        "workflow_key": workflow_key,
        "identifier_sha256": entry.get("remote_identifier", {}).get("sha256"),
    }


def _strip_shell_comment(value: str) -> str:
    quote = None
    escaped = False
    for idx, char in enumerate(value):
        if escaped:
            escaped = False
            continue
        if char == "\\":
            escaped = True
            continue
        if char in {"'", '"'}:
            quote = None if quote == char else char if quote is None else quote
            continue
        if char == "#" and quote is None:
            return value[:idx]
    return value


def _clean_path_value(value) -> str:
    text = str(value).strip()
    text = _strip_shell_comment(text).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        text = text[1:-1]
    return os.path.expandvars(os.path.expanduser(text))


def _datadir_from_bash_env(env_file: str) -> str:
    path = Path(env_file)
    if not path.is_file():
        return None
    pattern = re.compile(r"^\s*(?:export\s+)?DATADIR\s*=\s*(.+?)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            return _clean_path_value(match.group(1))
    return None


def _datadir_from_snakemake_config(config_file: str) -> str:
    path = Path(config_file)
    if not path.is_file():
        return None
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict) or data.get("datadir") is None:
        return None
    return _clean_path_value(data["datadir"])


def _resolve_workflow_datadir(workflow: WorkflowSpec) -> dict:
    if workflow.engine == "bash":
        source = workflow.helpers.get("env")
        datadir = _datadir_from_bash_env(source) if source else None
        return {"source": source, "source_key": "workflow.helpers.env", "datadir": datadir}

    if workflow.engine == "snakemake":
        source = workflow.config_file
        datadir = _datadir_from_snakemake_config(source) if source else None
        return {"source": source, "source_key": "workflow.config_file", "datadir": datadir}

    return {"source": None, "source_key": None, "datadir": None}


def _load_json_or_text_identifier(path: Path) -> dict:
    text = path.read_text(encoding="utf-8").strip()
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        payload = {"bundle": text}
    if isinstance(payload, str):
        payload = {"bundle": payload}
    if not isinstance(payload, dict):
        raise ParameterValidationError(
            f"Resource bundle identifier must be a string or JSON object: {path}"
        )
    return payload


def _observed_bundle_from_manifest(payload: dict) -> str:
    catalog_entry = (
        payload.get("catalog_entry")
        if isinstance(payload.get("catalog_entry"), dict)
        else {}
    )
    remote_identifier = (
        payload.get("remote_identifier")
        if isinstance(payload.get("remote_identifier"), dict)
        else {}
    )
    return (
        catalog_entry.get("key")
        or remote_identifier.get("bundle")
        or remote_identifier.get("bundle_key")
        or payload.get("bundle_key")
        or payload.get("bundle_slug")
    )


def validate_installed_resource_bundle(resource_bundle: dict, workflow: WorkflowSpec) -> dict:
    expected_key = resource_bundle.get("key")
    resolved = _resolve_workflow_datadir(workflow)
    datadir = resolved.get("datadir")
    result = {
        "status": "datadir_unresolved",
        "datadir": datadir,
        "source": resolved.get("source"),
        "source_key": resolved.get("source_key"),
        "checks": [],
    }

    if not datadir:
        return result

    datadir_path = Path(datadir).resolve()
    result["datadir"] = str(datadir_path)
    if not datadir_path.exists():
        result["status"] = "datadir_missing"
        return result

    metadata_found = False
    manifest_path = datadir_path / RESOURCE_INSTALL_MANIFEST
    if manifest_path.is_file():
        metadata_found = True
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ParameterValidationError(
                f"Invalid resource installation manifest: {manifest_path}"
            ) from exc
        observed = _observed_bundle_from_manifest(manifest)
        if observed != expected_key:
            raise ParameterValidationError(
                "Installed resource bundle does not match the selected resource_bundle.\n"
                f"  DATADIR: {datadir_path}\n"
                f"  manifest: {manifest_path}\n"
                f"  expected: {expected_key}\n"
                f"  observed: {observed}"
            )
        catalog_entry = (
            manifest.get("catalog_entry")
            if isinstance(manifest.get("catalog_entry"), dict)
            else {}
        )
        if catalog_entry:
            observed_fingerprint = _catalog_fingerprint(
                _catalog_entry_for_fingerprint(catalog_entry)
            )
            expected_fingerprint = resource_bundle.get("fingerprint")
            if expected_fingerprint and observed_fingerprint != expected_fingerprint:
                raise ParameterValidationError(
                    "Installed resource bundle manifest fingerprint does not match "
                    "the local catalog.\n"
                    f"  manifest: {manifest_path}\n"
                    f"  expected: {expected_fingerprint}\n"
                    f"  observed: {observed_fingerprint}"
                )
        else:
            observed_fingerprint = None
        result["checks"].append(
            {
                "file": str(manifest_path),
                "type": "installation_manifest",
                "bundle": observed,
            }
        )
        if observed_fingerprint:
            result["checks"][-1]["fingerprint"] = observed_fingerprint

    identifier_path = datadir_path / RESOURCE_IDENTIFIER
    if identifier_path.is_file():
        metadata_found = True
        payload = _load_json_or_text_identifier(identifier_path)
        observed = payload.get("bundle") or payload.get("bundle_key")
        if observed != expected_key:
            raise ParameterValidationError(
                "Installed resource bundle identifier does not match "
                "the selected resource_bundle.\n"
                f"  DATADIR: {datadir_path}\n"
                f"  identifier: {identifier_path}\n"
                f"  expected: {expected_key}\n"
                f"  observed: {observed}"
            )

        observed_sha256 = hashlib.sha256(identifier_path.read_bytes()).hexdigest()
        expected_sha256 = resource_bundle.get("identifier_sha256")
        if expected_sha256 and observed_sha256 != expected_sha256:
            raise ParameterValidationError(
                "Installed resource bundle identifier SHA-256 does not match the catalog.\n"
                f"  identifier: {identifier_path}\n"
                f"  expected: {expected_sha256}\n"
                f"  observed: {observed_sha256}"
            )
        result["checks"].append(
            {
                "file": str(identifier_path),
                "type": "bundle_identifier",
                "bundle": observed,
                "sha256": observed_sha256,
            }
        )

    result["status"] = "verified" if metadata_found else "metadata_not_found"
    return result
