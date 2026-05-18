import hashlib
import json

import pytest

from cbicall import resources as resources_mod
from cbicall.errors import ParameterValidationError
from cbicall.models import WorkflowSpec


def _workflow(engine="bash", *, env_file=None, config_file=None):
    return WorkflowSpec(
        engine=engine,
        pipeline="wes",
        mode="single",
        gatk_version="gatk-4.6",
        pipeline_version="v1",
        entrypoint="/workflow.sh",
        config_file=str(config_file) if config_file else None,
        helpers={"env": str(env_file)} if env_file else {},
    )


def test_validate_resource_catalog_rejects_malformed_catalogs(tmp_path):
    catalog = tmp_path / "catalog.json"

    catalog.write_text("[1, 2, 3]", encoding="utf-8")
    with pytest.raises(ParameterValidationError, match="failed resource catalog schema validation"):
        resources_mod.validate_resource_catalog(catalog)

    catalog.write_text("{not json", encoding="utf-8")
    with pytest.raises(ParameterValidationError, match="Invalid resource catalog JSON"):
        resources_mod.validate_resource_catalog(catalog)


def test_validate_resource_catalog_reports_bundle_shape_errors(tmp_path):
    catalog = tmp_path / "catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "bad-bundle": {
                        "type": "bundle",
                        "description": "",
                        "compatible_workflows": [
                            "bash/wes/single/gatk-4.6/v1",
                            "bash/wes/single/gatk-4.6/v1",
                            "",
                            "bash/wes/single/gatk-4.6",
                        ],
                        "remote_identifier": {
                            "sha256": "not-a-sha256",
                            "expected": "bad",
                        },
                        "archive": {"checksum_algorithm": "crc32"},
                    },
                    "not-an-object": "x",
                },
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError) as excinfo:
        resources_mod.validate_resource_catalog(catalog)

    message = str(excinfo.value)
    assert "failed resource catalog schema validation" in message
    assert "description" in message
    assert "has non-unique elements" in message
    assert "is too short" in message
    assert "does not match" in message
    assert "is not of type 'object'" in message
    assert "'crc32' is not one of" in message


def test_validate_resource_catalog_reports_semantic_errors_after_schema_passes(tmp_path):
    catalog = tmp_path / "catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "bundle-v1": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6"],
                        "remote_identifier": {
                            "expected": {"resource_key": "other-bundle"},
                        },
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ParameterValidationError) as excinfo:
        resources_mod.validate_resource_catalog(catalog)

    message = str(excinfo.value)
    assert "engine/pipeline/mode/gatk_version/pipeline_version" in message
    assert "expected.resource_key must match the resource key" in message


def test_validate_resource_catalog_accepts_legacy_registry_string_keys(tmp_path):
    catalog = tmp_path / "catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "bundle-v1": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6"],
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    registry = {
        "workflows": {
            "bash": {
                "versions": {
                    "gatk-4.6": {
                        "pipelines": {
                            "wes": {
                                "single": "wes_single.sh",
                                "bad_modes": "ignored.sh",
                            },
                            "bad_pipeline": "ignored",
                        }
                    }
                }
            }
        }
    }

    with pytest.raises(ParameterValidationError, match="engine/pipeline/mode"):
        resources_mod.validate_resource_catalog(catalog, registry)

    assert "bash/wes/single/gatk-4.6" in resources_mod._registry_workflow_keys(registry)


def test_validate_resource_catalog_can_filter_one_bundle(tmp_path):
    catalog = tmp_path / "catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "bundle-a": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
                    },
                    "bundle-b": {
                        "type": "bundle",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/cohort/gatk-4.6/v1"],
                    },
                },
            }
        ),
        encoding="utf-8",
    )

    summary = resources_mod.validate_resource_catalog(catalog, resource_key="bundle-b")
    assert summary["resource_key"] == "bundle-b"
    assert summary["bundle_resources"] == 1
    assert summary["compatible_workflows"] == 1

    with pytest.raises(ParameterValidationError, match="resource key is not defined"):
        resources_mod.validate_resource_catalog(catalog, resource_key="missing")


def test_validate_resource_catalog_accepts_non_bundle_resource_type(tmp_path):
    catalog = tmp_path / "catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "resources": {
                    "docker-gatk-v1": {
                        "type": "docker",
                        "version": "v1",
                        "compatible_workflows": ["bash/wes/single/gatk-4.6/v1"],
                    },
                },
            }
        ),
        encoding="utf-8",
    )

    summary = resources_mod.validate_resource_catalog(catalog, resource_key="docker-gatk-v1")
    assert summary["resource_key"] == "docker-gatk-v1"
    assert summary["resources"] == 1
    assert summary["bundle_resources"] == 0
    assert summary["compatible_workflows"] == 1


def test_datadir_parsing_handles_quotes_comments_and_missing_files(tmp_path, monkeypatch):
    datadir = tmp_path / "data dir"
    monkeypatch.setenv("CBICALL_TEST_DATA", str(datadir))
    env = tmp_path / "env.sh"
    env.write_text(
        "export OTHER=1\n"
        "export DATADIR=\"$CBICALL_TEST_DATA\" # trailing comment\n",
        encoding="utf-8",
    )

    assert resources_mod._datadir_from_bash_env(str(env)) == str(datadir)
    assert resources_mod._datadir_from_bash_env(str(tmp_path / "missing.sh")) is None

    no_datadir = tmp_path / "no_datadir.sh"
    no_datadir.write_text("export OTHER=1\n", encoding="utf-8")
    assert resources_mod._datadir_from_bash_env(str(no_datadir)) is None


def test_snakemake_datadir_parsing_and_unresolved_runtime_check(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("datadir: '~/cbicall-data'\n", encoding="utf-8")

    assert resources_mod._datadir_from_snakemake_config(str(config)).endswith("cbicall-data")
    assert resources_mod._datadir_from_snakemake_config(str(tmp_path / "missing.yaml")) is None

    no_datadir = tmp_path / "no-datadir.yaml"
    no_datadir.write_text("other: 1\n", encoding="utf-8")
    assert resources_mod._datadir_from_snakemake_config(str(no_datadir)) is None

    result = resources_mod.validate_installed_bundle_resource(
        {"key": "bundle-v1"},
        _workflow("snakemake", config_file=no_datadir),
    )
    assert result["status"] == "datadir_unresolved"
    assert result["source_key"] == "workflow.config_file"

    nextflow_result = resources_mod.validate_installed_bundle_resource(
        {"key": "bundle-v1"},
        _workflow("nextflow", config_file=no_datadir),
    )
    assert nextflow_result["status"] == "datadir_unresolved"
    assert nextflow_result["source_key"] == "workflow.config_file"


def test_validate_installed_bundle_resource_missing_and_metadata_not_found(tmp_path):
    env = tmp_path / "env.sh"
    missing_datadir = tmp_path / "missing-data"
    env.write_text(f"DATADIR={missing_datadir}\n", encoding="utf-8")

    result = resources_mod.validate_installed_bundle_resource({"key": "bundle-v1"}, _workflow(env_file=env))
    assert result["status"] == "datadir_missing"

    existing_datadir = tmp_path / "data"
    existing_datadir.mkdir()
    env.write_text(f"DATADIR={existing_datadir}\n", encoding="utf-8")

    result = resources_mod.validate_installed_bundle_resource({"key": "bundle-v1"}, _workflow(env_file=env))
    assert result["status"] == "metadata_not_found"
    assert result["checks"] == []


def test_validate_installed_resource_identifier_text_and_sha_mismatch(tmp_path):
    datadir = tmp_path / "data"
    datadir.mkdir()
    env = tmp_path / "env.sh"
    env.write_text(f"DATADIR={datadir}\n", encoding="utf-8")
    identifier = datadir / resources_mod.RESOURCE_IDENTIFIER
    identifier.write_text("bundle-v1\n", encoding="utf-8")

    result = resources_mod.validate_installed_bundle_resource({"key": "bundle-v1"}, _workflow(env_file=env))
    assert result["status"] == "verified"
    assert result["checks"][0]["resource_key"] == "bundle-v1"

    with pytest.raises(ParameterValidationError, match="identifier SHA-256 does not match"):
        resources_mod.validate_installed_bundle_resource(
            {"key": "bundle-v1", "identifier_sha256": "0" * 64},
            _workflow(env_file=env),
        )

    payload = json.dumps(["not", "an", "object"])
    identifier.write_text(payload, encoding="utf-8")
    with pytest.raises(ParameterValidationError, match="must be a string or JSON object"):
        resources_mod.validate_installed_bundle_resource({"key": "bundle-v1"}, _workflow(env_file=env))


def test_validate_installed_bundle_resource_manifest_edge_cases(tmp_path):
    datadir = tmp_path / "data"
    datadir.mkdir()
    env = tmp_path / "env.sh"
    env.write_text(f"DATADIR={datadir}\n", encoding="utf-8")
    manifest = datadir / resources_mod.RESOURCE_INSTALL_MANIFEST

    manifest.write_text("{not json", encoding="utf-8")
    with pytest.raises(ParameterValidationError, match="Invalid resource installation manifest"):
        resources_mod.validate_installed_bundle_resource({"key": "bundle-v1"}, _workflow(env_file=env))

    manifest.write_text(json.dumps({"resource_key": "bundle-v1"}), encoding="utf-8")
    result = resources_mod.validate_installed_bundle_resource({"key": "bundle-v1"}, _workflow(env_file=env))
    assert result["status"] == "verified"
    assert "fingerprint" not in result["checks"][0]

    manifest.write_text(json.dumps({"resource_key": "other-bundle"}), encoding="utf-8")
    with pytest.raises(ParameterValidationError, match="Installed bundle does not match"):
        resources_mod.validate_installed_bundle_resource({"key": "bundle-v1"}, _workflow(env_file=env))


def test_load_identifier_accepts_json_string(tmp_path):
    identifier = tmp_path / "id.json"
    identifier.write_text(json.dumps("bundle-v1"), encoding="utf-8")

    assert resources_mod._load_json_or_text_identifier(identifier) == {"resource_key": "bundle-v1"}
    assert hashlib.sha256(identifier.read_bytes()).hexdigest()
