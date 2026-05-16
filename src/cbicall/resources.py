import hashlib
import json
import os
import re
from pathlib import Path

import yaml

from .errors import ParameterValidationError
from .models import WorkflowSpec


RESOURCE_INSTALL_MANIFEST = "cbicall-resource-installation.json"
RESOURCE_IDENTIFIER = "cbicall-resource-id.json"
_HEX64 = re.compile(r"^[0-9a-fA-F]{64}$")
_CHECKSUM_ALGORITHMS = {"md5", "sha1", "sha256", "sha512"}
_RESOURCE_TYPES = {"bundle"}


def _catalog_fingerprint(entry: dict) -> str:
    payload = json.dumps(entry, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _catalog_entry_for_fingerprint(entry: dict) -> dict:
    return {
        key: value
        for key, value in entry.items()
        if key != "key" and not str(key).startswith("_")
    }


def _workflow_key(cfg_in: dict, workflow: WorkflowSpec = None) -> str:
    pipeline_version = workflow.pipeline_version if workflow else cfg_in.get("pipeline_version")
    parts = [
        cfg_in["workflow_engine"],
        cfg_in["pipeline"],
        cfg_in["mode"],
        cfg_in["gatk_version"],
    ]
    if pipeline_version:
        parts.append(str(pipeline_version))
    return "/".join(parts)


def _registry_workflow_keys(registry: dict) -> set:
    keys = set()
    workflows = registry.get("workflows", {}) if isinstance(registry, dict) else {}
    for engine, engine_cfg in workflows.items():
        versions = engine_cfg.get("versions", {}) if isinstance(engine_cfg, dict) else {}
        for gatk_version, version_cfg in versions.items():
            pipelines = version_cfg.get("pipelines", {}) if isinstance(version_cfg, dict) else {}
            for pipeline, modes in pipelines.items():
                if not isinstance(modes, dict):
                    continue
                for mode, mode_cfg in modes.items():
                    if isinstance(mode_cfg, dict):
                        implementations = mode_cfg.get("versions", {})
                        for pipeline_version in implementations:
                            keys.add(
                                "/".join(
                                    [
                                        str(engine),
                                        str(pipeline),
                                        str(mode),
                                        str(gatk_version),
                                        str(pipeline_version),
                                    ]
                                )
                            )
                    elif isinstance(mode_cfg, str):
                        keys.add("/".join([str(engine), str(pipeline), str(mode), str(gatk_version)]))
    return keys


def _validate_workflow_key_format(key: str) -> bool:
    parts = key.split("/")
    return len(parts) == 5 and all(parts)


def _catalog_resources(catalog: dict) -> dict:
    return catalog.get("resources", {}) if isinstance(catalog, dict) else {}


def validate_resource_catalog(
    catalog_path: Path,
    workflow_registry: dict = None,
    resource_key: str = None,
) -> dict:
    """Validate the resource catalog fields used by CBIcall."""
    path = Path(catalog_path)
    try:
        catalog = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise
    except json.JSONDecodeError as exc:
        raise ParameterValidationError(f"Invalid resource catalog JSON: {path}") from exc

    errors = []
    if not isinstance(catalog, dict):
        errors.append("catalog root must be a JSON object")
        catalog = {}

    if catalog.get("schema_version") is None:
        errors.append("schema_version is required")

    resources = catalog.get("resources")
    if not isinstance(resources, dict) or not resources:
        errors.append("resources must be a non-empty object")
        resources = {}

    selected_resource_key = str(resource_key).strip() if resource_key is not None else None
    if selected_resource_key:
        if selected_resource_key not in resources:
            errors.append(f"resource key is not defined in resources: {selected_resource_key}")
            resources = {}
        else:
            resources = {selected_resource_key: resources[selected_resource_key]}

    known_workflows = _registry_workflow_keys(workflow_registry) if workflow_registry else None
    compatible_count = 0
    bundle_count = 0

    for current_resource_key, entry in resources.items():
        label = f"resources.{current_resource_key}"
        if not isinstance(current_resource_key, str) or not current_resource_key.strip():
            errors.append("resource keys must be non-empty strings")
            continue
        if not isinstance(entry, dict):
            errors.append(f"{label} must be an object")
            continue

        resource_type = entry.get("type")
        if resource_type not in _RESOURCE_TYPES:
            allowed = ", ".join(sorted(_RESOURCE_TYPES))
            errors.append(f"{label}.type must be one of: {allowed}")
        if resource_type == "bundle":
            bundle_count += 1

        for field in ["type", "compatible_workflows"]:
            if field not in entry:
                errors.append(f"{label}.{field} is required")

        for field in ["description"]:
            if field in entry and not str(entry[field]).strip():
                errors.append(f"{label}.{field} must be non-empty when present")

        compatible = entry.get("compatible_workflows")
        if not isinstance(compatible, list) or not compatible:
            errors.append(f"{label}.compatible_workflows must be a non-empty array")
            compatible = []

        seen = set()
        for workflow_key in compatible:
            compatible_count += 1
            if not isinstance(workflow_key, str) or not workflow_key.strip():
                errors.append(f"{label}.compatible_workflows entries must be non-empty strings")
                continue
            if workflow_key in seen:
                errors.append(f"{label}.compatible_workflows contains duplicate entry: {workflow_key}")
            seen.add(workflow_key)
            if not _validate_workflow_key_format(workflow_key):
                errors.append(
                    f"{label}.compatible_workflows entry must use "
                    "engine/pipeline/mode/gatk_version/pipeline_version: "
                    f"{workflow_key}"
                )
            elif known_workflows is not None and workflow_key not in known_workflows:
                errors.append(
                    f"{label}.compatible_workflows entry is not defined in the workflow registry: "
                    f"{workflow_key}"
                )

        remote = entry.get("remote_identifier")
        if remote is not None:
            if not isinstance(remote, dict):
                errors.append(f"{label}.remote_identifier must be an object when present")
            else:
                sha256 = remote.get("sha256")
                if sha256 is not None and not _HEX64.match(str(sha256)):
                    errors.append(f"{label}.remote_identifier.sha256 must be a 64-character hex string")
                expected = remote.get("expected")
                if expected is not None and not isinstance(expected, dict):
                    errors.append(f"{label}.remote_identifier.expected must be an object when present")
                elif isinstance(expected, dict):
                    expected_resource_key = expected.get("resource_key")
                    if expected_resource_key is not None and expected_resource_key != current_resource_key:
                        errors.append(
                            f"{label}.remote_identifier.expected.resource_key must match the resource key"
                        )

        archive = entry.get("archive")
        if archive is not None:
            if not isinstance(archive, dict):
                errors.append(f"{label}.archive must be an object when present")
            else:
                algorithm = archive.get("checksum_algorithm")
                if algorithm is not None and str(algorithm) not in _CHECKSUM_ALGORITHMS:
                    allowed = ", ".join(sorted(_CHECKSUM_ALGORITHMS))
                    errors.append(f"{label}.archive.checksum_algorithm must be one of: {allowed}")

    if errors:
        lines = [f"{path} failed resource catalog validation:"]
        lines.extend(f"- {error}" for error in errors)
        raise ParameterValidationError("\n".join(lines))

    return {
        "path": str(path),
        "schema_version": catalog.get("schema_version"),
        "resources": len(resources),
        "bundle_resources": bundle_count,
        "resource_key": selected_resource_key,
        "compatible_workflows": compatible_count,
    }


def build_bundle_resource_metadata(cfg_in: dict, project_root: Path, workflow: WorkflowSpec = None) -> dict:
    catalog_path = project_root / "resources" / "cbicall-resource-catalog.json"
    selected_resource_key = cfg_in["resource"]

    if not catalog_path.is_file():
        entry = {"compatible_workflows": []}
        return {
            "key": selected_resource_key,
            "catalog": None,
            "fingerprint": _catalog_fingerprint(_catalog_entry_for_fingerprint(entry)),
            "compatible": None,
            "workflow_key": _workflow_key(cfg_in, workflow),
            "runtime_check": {"status": "catalog_not_found"},
        }

    try:
        catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ParameterValidationError(f"Invalid resource catalog JSON: {catalog_path}") from exc

    resources = _catalog_resources(catalog)
    entry = resources.get(selected_resource_key)
    if entry is None:
        raise ParameterValidationError(f"Resource key '{selected_resource_key}' is not defined in {catalog_path}")
    if entry.get("type") != "bundle":
        raise ParameterValidationError(
            f"Resource key '{selected_resource_key}' is not a bundle resource in {catalog_path}"
        )

    workflow_key = _workflow_key(cfg_in, workflow)
    compatible_workflows = entry.get("compatible_workflows", [])
    compatible = workflow_key in compatible_workflows if compatible_workflows else None
    if compatible_workflows and not compatible:
        raise ParameterValidationError(
            f"Resource key '{selected_resource_key}' is not declared compatible "
            f"with workflow '{workflow_key}'"
        )

    return {
        "key": selected_resource_key,
        "catalog": str(catalog_path),
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
        payload = {"resource_key": text}
    if isinstance(payload, str):
        payload = {"resource_key": payload}
    if not isinstance(payload, dict):
        raise ParameterValidationError(
            f"Resource identifier must be a string or JSON object: {path}"
        )
    return payload


def _observed_resource_key_from_manifest(payload: dict) -> str:
    catalog_entry = (
        payload.get("catalog_entry")
        if isinstance(payload.get("catalog_entry"), dict)
        else {}
    )
    return payload.get("resource_key") or catalog_entry.get("key")


def validate_installed_bundle_resource(bundle_resource: dict, workflow: WorkflowSpec) -> dict:
    expected_key = bundle_resource.get("key")
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
        observed = _observed_resource_key_from_manifest(manifest)
        if observed != expected_key:
            raise ParameterValidationError(
                "Installed bundle does not match the selected resource.\n"
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
            expected_fingerprint = bundle_resource.get("fingerprint")
            if expected_fingerprint and observed_fingerprint != expected_fingerprint:
                raise ParameterValidationError(
                    "Installed bundle manifest fingerprint does not match "
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
                "resource_key": observed,
            }
        )
        if observed_fingerprint:
            result["checks"][-1]["fingerprint"] = observed_fingerprint

    identifier_path = datadir_path / RESOURCE_IDENTIFIER
    if identifier_path.is_file():
        metadata_found = True
        payload = _load_json_or_text_identifier(identifier_path)
        observed = payload.get("resource_key")
        if observed != expected_key:
            raise ParameterValidationError(
                "Installed resource identifier does not match "
                "the selected resource.\n"
                f"  DATADIR: {datadir_path}\n"
                f"  identifier: {identifier_path}\n"
                f"  expected: {expected_key}\n"
                f"  observed: {observed}"
            )

        observed_sha256 = hashlib.sha256(identifier_path.read_bytes()).hexdigest()
        expected_sha256 = bundle_resource.get("identifier_sha256")
        if expected_sha256 and observed_sha256 != expected_sha256:
            raise ParameterValidationError(
                "Installed resource identifier SHA-256 does not match the catalog.\n"
                f"  identifier: {identifier_path}\n"
                f"  expected: {expected_sha256}\n"
                f"  observed: {observed_sha256}"
            )
        result["checks"].append(
            {
                "file": str(identifier_path),
                "type": "resource_identifier",
                "resource_key": observed,
                "sha256": observed_sha256,
            }
        )

    result["status"] = "verified" if metadata_found else "metadata_not_found"
    return result
