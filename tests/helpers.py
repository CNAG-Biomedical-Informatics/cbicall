# tests/helpers.py
import json
import stat
from pathlib import Path
from typing import Optional, Dict, Any, List

from cbicall import config as config_mod


def make_executable(path: Path) -> None:
    path.write_text("#!/bin/sh\n", encoding="utf-8")
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR)


def write_workflow_schema(path: Path) -> None:
    """
    Schema for workflows/registry/cbicall-workflow-registry.yaml.

    Key rule: pipeline modes must define at least one of {single, cohort}.
    """
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "cbicall:workflow-registry.schema.json",
        "type": "object",
        "additionalProperties": False,
        "required": ["workflows"],
        "properties": {
            "workflows": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "bash": {"$ref": "#/$defs/backend"},
                    "snakemake": {"$ref": "#/$defs/backend"},
                    "nextflow": {"$ref": "#/$defs/backend"},
                },
                "required": ["bash"],
            }
        },
        "$defs": {
            "backend": {
                "type": "object",
                "additionalProperties": False,
                "required": ["base_dir", "versions"],
                "properties": {
                    "base_dir": {"type": "string", "minLength": 1},
                    "versions": {
                        "type": "object",
                        "minProperties": 1,
                        "additionalProperties": {"$ref": "#/$defs/version"},
                    },
                },
            },
            "version": {
                "type": "object",
                "additionalProperties": False,
                "required": ["pipelines"],
                "properties": {
                    "helpers": {
                        "type": "object",
                        "additionalProperties": {"type": "string", "minLength": 1},
                    },
                    "profiles": {
                        "type": "object",
                        "minProperties": 1,
                        "additionalProperties": {"$ref": "#/$defs/profile"},
                    },
                    "pipelines": {
                        "type": "object",
                        "minProperties": 1,
                        "additionalProperties": {"$ref": "#/$defs/pipelineModes"},
                    },
                },
            },
            "profile": {
                "type": "object",
                "minProperties": 1,
                "additionalProperties": {"type": "string", "minLength": 1},
            },
            "pipelineModes": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "single": {"$ref": "#/$defs/pipelineImplementationSet"},
                    "cohort": {"$ref": "#/$defs/pipelineImplementationSet"},
                },
                "minProperties": 1,
            },
            "pipelineImplementationSet": {
                "type": "object",
                "additionalProperties": False,
                "required": ["default", "versions"],
                "properties": {
                    "default": {"type": "string", "minLength": 1},
                    "versions": {
                        "type": "object",
                        "minProperties": 1,
                        "additionalProperties": {"$ref": "#/$defs/pipelineImplementation"},
                    },
                },
            },
            "pipelineImplementation": {
                "type": "object",
                "additionalProperties": False,
                "oneOf": [
                    {"required": ["script"]},
                    {"required": ["provider", "source", "release"]},
                ],
                "properties": {
                    "script": {"type": "string", "minLength": 1},
                    "provider": {"enum": ["nf-core"]},
                    "source": {"type": "string", "minLength": 1},
                    "release": {"type": "string", "minLength": 1},
                    "default_outdir": {"type": "string", "minLength": 1},
                    "canonical_outputs": {
                        "type": "array",
                        "minItems": 1,
                        "items": {"$ref": "#/$defs/canonicalOutput"},
                    },
                },
            },
            "canonicalOutput": {
                "type": "object",
                "additionalProperties": False,
                "required": ["name", "type", "pattern"],
                "properties": {
                    "name": {"type": "string", "minLength": 1},
                    "type": {"enum": ["vcf"]},
                    "pattern": {"type": "string", "minLength": 1},
                },
            },
        },
    }
    path.write_text(json.dumps(schema, indent=2), encoding="utf-8")


def _append_pipeline_modes(lines: List[str], modes: Dict[str, str], indent: str) -> None:
    for mode, script in modes.items():
        lines.extend(
            [
                f"{indent}{mode}:",
                f"{indent}  default: \"v1\"",
                f"{indent}  versions:",
                f"{indent}    v1:",
                f"{indent}      script: \"{script}\"",
            ]
        )


def write_registry(
    path: Path,
    *,
    gatk_ver: str,
    include_bash: bool = True,
    include_snakemake: bool = False,
    bash_pipelines: Optional[Dict[str, Dict[str, str]]] = None,
    snakemake_pipelines: Optional[Dict[str, Dict[str, str]]] = None,
    bash_profiles: Optional[Dict[str, Dict[str, str]]] = None,
    snakemake_profiles: Optional[Dict[str, Dict[str, str]]] = None,
) -> None:
    bash_pipelines = bash_pipelines or {"wes": {"single": "wes_single.sh"}}
    snakemake_pipelines = snakemake_pipelines or {}
    bash_profiles = bash_profiles or {}
    snakemake_profiles = snakemake_profiles or {}

    lines = ["workflows:"]

    if include_bash:
        lines += [
            "  bash:",
            "    base_dir: \"workflows/bash\"",
            "    versions:",
            f"      {gatk_ver}:",
            "        helpers:",
            "          env: \"env.sh\"",
            "          coverage: \"coverage.sh\"",
            "          jaccard: \"jaccard.sh\"",
            "          vcf2sex: \"vcf2sex.sh\"",
            "          vcf2hash: \"vcf2hash.sh\"",
        ]
        if bash_profiles:
            lines.append("        profiles:")
            for profile_name, helpers in bash_profiles.items():
                lines.append(f"          {profile_name}:")
                for helper_name, helper_path in helpers.items():
                    lines.append(f"            {helper_name}: \"{helper_path}\"")
        lines.append("        pipelines:")
        for pipeline, modes in bash_pipelines.items():
            lines.append(f"          {pipeline}:")
            _append_pipeline_modes(lines, modes, "            ")

    if include_snakemake:
        lines += [
            "  snakemake:",
            "    base_dir: \"workflows/snakemake\"",
            "    versions:",
            f"      {gatk_ver}:",
            "        helpers:",
            "          config: \"config.yaml\"",
        ]
        if snakemake_profiles:
            lines.append("        profiles:")
            for profile_name, helpers in snakemake_profiles.items():
                lines.append(f"          {profile_name}:")
                for helper_name, helper_path in helpers.items():
                    lines.append(f"            {helper_name}: \"{helper_path}\"")
        lines.append("        pipelines:")
        for pipeline, modes in snakemake_pipelines.items():
            lines.append(f"          {pipeline}:")
            _append_pipeline_modes(lines, modes, "            ")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fake_project(
    monkeypatch,
    tmp_path: Path,
    *,
    gatk_ver: str,
    registry_kwargs: Dict[str, Any],
) -> Path:
    """
    Create the fake project structure expected by config.py and monkeypatch
    config_mod.__file__ so parents[2] == tmp_path.
    """
    root = tmp_path

    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True, exist_ok=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy\n", encoding="utf-8")
    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))

    cfg_dir = root / "workflows" / "registry"
    sch_dir = root / "workflows" / "schema"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    sch_dir.mkdir(parents=True, exist_ok=True)

    write_registry(cfg_dir / "cbicall-workflow-registry.yaml", gatk_ver=gatk_ver, **registry_kwargs)
    write_workflow_schema(sch_dir / "cbicall-workflow-registry.schema.json")

    return root
