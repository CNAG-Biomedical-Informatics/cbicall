# tests/helpers.py
import json
import stat
from pathlib import Path

from cbicall import config as config_mod


def make_executable(path: Path) -> None:
    path.write_text("#!/bin/sh\n", encoding="utf-8")
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR)


def write_workflow_schema(path: Path) -> None:
    """
    Schema for workflows/config/cbicall.workflows.yaml.

    Key rule: pipeline modes must define at least one of {single, cohort}.
    """
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "cbicall:workflows.schema.json",
        "type": "object",
        "additionalProperties": False,
        "required": ["workflows"],
        "properties": {
            "workflows": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "bash": {"$ref": "#/$defs/engine"},
                    "snakemake": {"$ref": "#/$defs/engine"},
                    "nextflow": {"$ref": "#/$defs/engine"},
                },
                "required": ["bash"],
            }
        },
        "$defs": {
            "engine": {
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
                    "common": {
                        "type": "object",
                        "additionalProperties": {"type": "string", "minLength": 1},
                    },
                    "pipelines": {
                        "type": "object",
                        "minProperties": 1,
                        "additionalProperties": {"$ref": "#/$defs/pipelineModes"},
                    },
                },
            },
            "pipelineModes": {
                "type": "object",
                "additionalProperties": False,
                "properties": {
                    "single": {"type": "string", "minLength": 1},
                    "cohort": {"type": "string", "minLength": 1},
                },
                "minProperties": 1,
            },
        },
    }
    path.write_text(json.dumps(schema, indent=2), encoding="utf-8")


def write_registry(
    path: Path,
    *,
    gatk_ver: str,
    include_bash: bool = True,
    include_snakemake: bool = False,
    bash_pipelines: dict | None = None,
    snakemake_pipelines: dict | None = None,
) -> None:
    bash_pipelines = bash_pipelines or {"wes": {"single": "wes_single.sh"}}
    snakemake_pipelines = snakemake_pipelines or {}

    lines = ["workflows:"]

    if include_bash:
        lines += [
            "  bash:",
            "    base_dir: \"workflows/bash\"",
            "    versions:",
            f"      {gatk_ver}:",
            "        common:",
            "          parameters: \"parameters.sh\"",
            "          coverage: \"coverage.sh\"",
            "          jaccard: \"jaccard.sh\"",
            "          vcf2sex: \"vcf2sex.sh\"",
            "        pipelines:",
        ]
        for pipeline, modes in bash_pipelines.items():
            lines.append(f"          {pipeline}:")
            for mode, script in modes.items():
                lines.append(f"            {mode}: \"{script}\"")

    if include_snakemake:
        lines += [
            "  snakemake:",
            "    base_dir: \"workflows/snakemake\"",
            "    versions:",
            f"      {gatk_ver}:",
            "        common:",
            "          config: \"config.yaml\"",
            "        pipelines:",
        ]
        for pipeline, modes in snakemake_pipelines.items():
            lines.append(f"          {pipeline}:")
            for mode, script in modes.items():
                lines.append(f"            {mode}: \"{script}\"")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fake_project(
    monkeypatch,
    tmp_path: Path,
    *,
    gatk_ver: str,
    registry_kwargs: dict,
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

    cfg_dir = root / "workflows" / "config"
    sch_dir = root / "workflows" / "schema"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    sch_dir.mkdir(parents=True, exist_ok=True)

    write_registry(cfg_dir / "cbicall.workflows.yaml", gatk_ver=gatk_ver, **registry_kwargs)
    write_workflow_schema(sch_dir / "workflows.schema.json")

    return root
