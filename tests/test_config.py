# tests/test_config.py
from pathlib import Path

import pytest

from cbicall import config as config_mod
from cbicall.errors import ParameterValidationError, WorkflowResolutionError
from helpers import fake_project, make_executable, write_workflow_schema


# -----------------------
# Local helper utilities
# -----------------------

def _mk_fake_root(monkeypatch, tmp_path: Path) -> Path:
    """
    Minimal project root + patched config_mod.__file__ so _get_project_root works.
    """
    root = tmp_path
    src_dir = root / "src" / "cbicall"
    src_dir.mkdir(parents=True, exist_ok=True)
    fake_config_file = src_dir / "config.py"
    fake_config_file.write_text("# dummy\n", encoding="utf-8")
    monkeypatch.setattr(config_mod, "__file__", str(fake_config_file))
    return root


def _write_registry_and_schema(root: Path, *, registry_lines: str) -> None:
    cfg_dir = root / "workflows" / "config"
    sch_dir = root / "workflows" / "schema"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    sch_dir.mkdir(parents=True, exist_ok=True)

    (cfg_dir / "cbicall.workflows.yaml").write_text(registry_lines, encoding="utf-8")
    write_workflow_schema(sch_dir / "workflows.schema.json")


def _touch_snakemake_files(root: Path, gatk_ver: str, snakefile_name: str) -> Path:
    smk_dir = root / "workflows" / "snakemake" / gatk_ver
    smk_dir.mkdir(parents=True, exist_ok=True)
    (smk_dir / snakefile_name).write_text("# dummy\n", encoding="utf-8")
    (smk_dir / "config.yaml").write_text("dummy: 1\n", encoding="utf-8")
    return smk_dir


def _touch_bash_files(root: Path, gatk_ver: str, script_name: str, *, make_exec: bool = True) -> Path:
    bash_dir = root / "workflows" / "bash" / gatk_ver
    bash_dir.mkdir(parents=True, exist_ok=True)

    for name in ["env.sh", "coverage.sh", "jaccard.sh", "vcf2sex.sh", script_name]:
        p = bash_dir / name
        p.write_text("#!/bin/sh\n", encoding="utf-8")
        if make_exec:
            make_executable(p)  # writes shebang + adds +x
    return bash_dir


def _minimal_bash_registry_block(gatk_ver: str = "gatk-4.6") -> str:
    # NOTE: schema requires workflows.bash to exist (required=["bash"])
    return (
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    versions:\n"
        f"      {gatk_ver}:\n"
        "        common:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single: \"wes_single.sh\"\n"
    )


# ------------------------
# read_param_file() tests
# ------------------------

def test_read_param_file_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="Parameters file not found"):
        config_mod.read_param_file(str(tmp_path / "nope.yaml"))


def test_read_param_file_snakemake_gatk35_is_blocked(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_engine: snakemake\n"
        "gatk_version: gatk-3.5\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="snakemake.*not supported.*gatk_version='gatk-3.5'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_mit_snakemake_blocked(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: mit\n"
        "workflow_engine: snakemake\n"
        "gatk_version: gatk-4.6\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="pipeline='mit'.*workflow_engine='snakemake'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_mit_forces_rsrs_and_rejects_other_genome(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: mit\n"
        "workflow_engine: bash\n"
        "gatk_version: gatk-3.5\n"
        "genome: b37\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="pipeline='mit'.*genome.*fixed to 'rsrs'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_hg38_only_allowed_for_wgs(tmp_path):
    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_engine: bash\n"
        "gatk_version: gatk-4.6\n"
        "genome: hg38\n",
        encoding="utf-8",
    )
    with pytest.raises(ParameterValidationError, match="genome='hg38' is only supported for pipeline='wgs'"):
        config_mod.read_param_file(str(p))


def test_read_param_file_resolves_input_dir_relative_to_yaml(tmp_path):
    sample_dir = tmp_path / "inputs" / "SAMPLE_001"
    sample_dir.mkdir(parents=True)

    p = tmp_path / "params.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_engine: bash\n"
        "gatk_version: gatk-4.6\n"
        "input_dir: inputs/SAMPLE_001\n",
        encoding="utf-8",
    )

    cfg = config_mod.read_param_file(str(p))
    assert cfg["input_dir"] == str(sample_dir.resolve())


def test_read_param_file_accepts_project_dir(tmp_path):
    p1 = tmp_path / "params_project_dir.yaml"
    p1.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_engine: bash\n"
        "gatk_version: gatk-4.6\n"
        "project_dir: my_run\n",
        encoding="utf-8",
    )
    cfg1 = config_mod.read_param_file(str(p1))
    assert cfg1["project_dir"] == "my_run"


def test_read_param_file_rejects_workflow_rule_without_allow_partial_run(tmp_path):
    p = tmp_path / "params_partial.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_engine: snakemake\n"
        "gatk_version: gatk-4.6\n"
        "workflow_rule: call_variants\n",
        encoding="utf-8",
    )
    with pytest.raises(
        ParameterValidationError,
        match=r"workflow_rule was set to 'call_variants', but --allow-partial-run was not provided",
    ):
        config_mod.read_param_file(str(p))


def test_read_param_file_accepts_partial_run_fields(tmp_path):
    p = tmp_path / "params_partial_ok.yaml"
    p.write_text(
        "mode: single\n"
        "pipeline: wes\n"
        "workflow_engine: snakemake\n"
        "gatk_version: gatk-4.6\n"
        "workflow_rule: call_variants\n"
        "allow_partial_run: true\n",
        encoding="utf-8",
    )
    cfg = config_mod.read_param_file(str(p))
    assert cfg["workflow_rule"] == "call_variants"
    assert cfg["allow_partial_run"] is True


# -----------------------------------------
# set_config_values() missing-registry/schema
# -----------------------------------------

def test_set_config_values_registry_missing_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Create schema only, omit registry
    sch_dir = root / "workflows" / "schema"
    sch_dir.mkdir(parents=True, exist_ok=True)
    write_workflow_schema(sch_dir / "workflows.schema.json")

    with pytest.raises(FileNotFoundError, match="Workflow registry not found"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
        )


def test_set_config_values_schema_missing_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    cfg_dir = root / "workflows" / "config"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    (cfg_dir / "cbicall.workflows.yaml").write_text("workflows:\n  bash: {}\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="Workflow schema not found"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
        )


# -----------------------------
# Enum + combo validation tests
# -----------------------------

def test_set_config_values_invalid_enum_mode_raises(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(ParameterValidationError, match="Invalid value for 'mode'"):
        config_mod.set_config_values(
            {"mode": "banana", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
        )


def test_set_config_values_partial_run_metadata(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_snakemake=True, snakemake_pipelines={"wes": {"single": "wes_single.smk"}}),
    )
    _touch_snakemake_files(root, "gatk-4.6", "wes_single.smk")

    cfg = config_mod.set_config_values(
        {
            "mode": "single",
            "pipeline": "wes",
            "workflow_engine": "snakemake",
            "gatk_version": "gatk-4.6",
            "workflow_rule": "call_variants",
            "allow_partial_run": True,
        }
    )
    assert cfg["workflow_rule"] == "call_variants"
    assert cfg["allow_partial_run"] is True
    assert cfg["run_mode"] == "partial"


def test_set_config_values_invalid_combo_raises(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    # wgs is not allowed for gatk-3.5 in _ALLOWED_COMBOS
    with pytest.raises(ParameterValidationError, match="not supported for GATK version"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wgs", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
        )


# -----------------------------
# Schema validation error tests
# -----------------------------

def test_load_workflow_registry_schema_validation_fails(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Valid schema, invalid registry (missing required top-level "workflows")
    bad_registry = "not_workflows:\n  x: 1\n"
    _write_registry_and_schema(root, registry_lines=bad_registry)

    registry_yaml = root / "workflows" / "config" / "cbicall.workflows.yaml"
    schema_json = root / "workflows" / "schema" / "workflows.schema.json"

    with pytest.raises(ParameterValidationError, match="failed schema validation"):
        config_mod.load_workflow_registry(registry_yaml, schema_json)


def test__validate_enum_allows_none():
    # hits the early-return branch in _validate_enum
    config_mod._validate_enum("anything", None, {"a", "b"})


def test__validate_with_schema_formats_root_and_nested_locations():
    # Schema that expects an object with a nested requirement
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "required": ["workflows"],
        "properties": {
            "workflows": {
                "type": "object",
                "required": ["bash"],
                "properties": {"bash": {"type": "object"}},
            }
        },
        "additionalProperties": False,
    }

    # 1) ROOT error: wrong top-level type -> e.path empty -> loc "(root)"
    with pytest.raises(ParameterValidationError) as excinfo1:
        config_mod._validate_with_schema(data="not-an-object", schema=schema, label="X")
    msg1 = str(excinfo1.value)
    assert "failed schema validation" in msg1
    assert "(root)" in msg1

    # 2) NESTED error: missing required nested key -> e.path is "workflows"
    with pytest.raises(ParameterValidationError) as excinfo2:
        config_mod._validate_with_schema(data={"workflows": {}}, schema=schema, label="Y")
    msg2 = str(excinfo2.value)
    assert "failed schema validation" in msg2
    assert "workflows" in msg2


# ---------------------------------------------
# Registry wiring errors (engine/version/pipeline/mode)
# ---------------------------------------------

def test_set_config_values_engine_not_in_registry_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Registry defines only bash, request snakemake
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    versions:\n"
        "      gatk-4.6:\n"
        "        common:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="Engine not defined in workflow registry"):
        config_mod.set_config_values(
            {
                "mode": "single",
                "pipeline": "wes",
                "workflow_engine": "snakemake",  # not present in registry
                "gatk_version": "gatk-4.6",      # avoids snakemake+3.5 guardrail
                "genome": "b37",
            }
        )


def test_set_config_values_version_not_in_registry_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Only gatk-3.5 exists in registry; request gatk-4.6
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    versions:\n"
        "      gatk-3.5:\n"
        "        common:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="Version not defined"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-4.6"}
        )


def test_set_config_values_pipeline_not_in_registry_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Registry defines only wes; request mit
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    versions:\n"
        "      gatk-3.5:\n"
        "        common:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="Pipeline not defined"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "mit", "workflow_engine": "bash", "gatk_version": "gatk-3.5", "genome": "rsrs"}
        )


def test_set_config_values_mode_not_in_registry_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # wes only defines single; request cohort
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    versions:\n"
        "      gatk-3.5:\n"
        "        common:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="Mode not defined"):
        config_mod.set_config_values(
            {"mode": "cohort", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
        )


# --------------------------------
# Missing common keys branches
# --------------------------------

def test_set_config_values_bash_missing_common_keys_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Omit coverage/jaccard/vcf2sex in common
    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    versions:\n"
        "      gatk-3.5:\n"
        "        common:\n"
        "          env: \"env.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="missing common keys for bash"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
        )


def test_set_config_values_snakemake_missing_common_config_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    # Include minimal bash section to satisfy schema required=["bash"]
    reg = (
        "workflows:\n"
        + _minimal_bash_registry_block("gatk-4.6")
        + "  snakemake:\n"
          "    base_dir: \"workflows/snakemake\"\n"
          "    versions:\n"
          "      gatk-4.6:\n"
          "        common:\n"
          "          something_else: \"x.yaml\"\n"
          "        pipelines:\n"
          "          wgs:\n"
          "            cohort: \"wgs_cohort.smk\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    # satisfy bash files referenced by schema-required bash block
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    # create snakemake pipeline file so we pass "missing files" guardrail
    smk_dir = root / "workflows" / "snakemake" / "gatk-4.6"
    smk_dir.mkdir(parents=True, exist_ok=True)
    (smk_dir / "wgs_cohort.smk").write_text("# dummy\n", encoding="utf-8")

    with pytest.raises(WorkflowResolutionError, match="missing common key 'config'"):
        config_mod.set_config_values(
            {"mode": "cohort", "pipeline": "wgs", "workflow_engine": "snakemake", "gatk_version": "gatk-4.6"}
        )


def test_set_config_values_snakemake_missing_files_triggers_guard(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        + _minimal_bash_registry_block("gatk-4.6")
        + "  snakemake:\n"
          "    base_dir: \"workflows/snakemake\"\n"
          "    versions:\n"
          "      gatk-4.6:\n"
          "        common:\n"
          "          config: \"config.yaml\"\n"
          "        pipelines:\n"
          "          wgs:\n"
          "            cohort: \"wgs_cohort.smk\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    # satisfy bash files referenced by schema-required bash block
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    # Intentionally do NOT create:
    #   workflows/snakemake/gatk-4.6/wgs_cohort.smk
    #   workflows/snakemake/gatk-4.6/config.yaml
    # so the snakemake must_exist / missing_files guard triggers.
    with pytest.raises(WorkflowResolutionError, match="referenced in the registry do not exist"):
        config_mod.set_config_values(
            {"mode": "cohort", "pipeline": "wgs", "workflow_engine": "snakemake", "gatk_version": "gatk-4.6"}
        )


# --------------------------
# nextflow branch
# --------------------------

def test_set_config_values_nextflow_declared_but_not_implemented(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        + _minimal_bash_registry_block("gatk-4.6")
        + "  nextflow:\n"
          "    base_dir: \"workflows/nextflow\"\n"
          "    versions:\n"
          "      gatk-4.6:\n"
          "        common: {}\n"
          "        pipelines:\n"
          "          wgs:\n"
          "            single: \"main.nf\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    with pytest.raises(WorkflowResolutionError, match="declared but not implemented"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wgs", "workflow_engine": "nextflow", "gatk_version": "gatk-4.6", "genome": "b37"}
        )


# --------------------------
# “Normal” branches
# --------------------------

def test_set_config_values_sample_sets_output_basename(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    sample_dir = tmp_path / "SAMPLE_001"
    sample_dir.mkdir()

    cfg = config_mod.set_config_values(
        {
            "mode": "single",
            "pipeline": "wes",
            "workflow_engine": "bash",
            "gatk_version": "gatk-3.5",
            "input_dir": str(sample_dir),
        }
    )

    assert cfg["output_basename"] == "SAMPLE_001"
    assert cfg["project_dir"].startswith(str(sample_dir))


def test_set_config_values_builds_normalized_workflow_structure(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-4.6"}
    )

    assert cfg["workflow"]["engine"] == "bash"
    assert cfg["workflow"]["pipeline"] == "wes"
    assert cfg["workflow"]["mode"] == "single"
    assert cfg["workflow"]["gatk_version"] == "gatk-4.6"
    assert cfg["workflow"]["entrypoint"].endswith("workflows/bash/gatk-4.6/wes_single.sh")
    assert cfg["workflow"]["helpers"]["env"].endswith("workflows/bash/gatk-4.6/env.sh")
    assert cfg["inputs"] == {"input_dir": None, "sample_map": None}


def test_set_config_values_capture_branches(monkeypatch, tmp_path):
    # 1) gatk-3.5 + wes (non-wgs) -> "Agilent SureSelect"
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")
    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
    )
    assert cfg["capture_label"] == "Agilent SureSelect"

    # 2) mit -> "MToolBox_rsrs" (force x86_64 so MIT isn't blocked)
    root2 = fake_project(
        monkeypatch,
        tmp_path / "root2",
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"mit": {"single": "mit_single.sh"}}),
    )
    _touch_bash_files(root2, "gatk-3.5", "mit_single.sh")
    monkeypatch.setattr(config_mod.platform, "machine", lambda: "x86_64")

    cfg2 = config_mod.set_config_values(
        {"mode": "single", "pipeline": "mit", "workflow_engine": "bash", "gatk_version": "gatk-3.5", "genome": "rsrs"}
    )
    assert cfg2["capture_label"] == "MToolBox_rsrs"

    # 3) gatk-4.6 + wes (non-wgs) -> "GATK_bundle_<genome>"
    root3 = fake_project(
        monkeypatch,
        tmp_path / "root3",
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root3, "gatk-4.6", "wes_single.sh")
    cfg3 = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-4.6", "genome": "b37"}
    )
    assert cfg3["capture_label"] == "GATK_bundle_b37"


def test_set_config_values_mit_on_arm64_raises(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"mit": {"single": "mit_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "mit_single.sh")

    monkeypatch.setattr(config_mod.platform, "machine", lambda: "aarch64")

    with pytest.raises(WorkflowResolutionError, match=r"mit_single cannot be performed"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "mit", "workflow_engine": "bash", "gatk_version": "gatk-3.5", "genome": "rsrs"}
        )


def test_set_config_values_nproc_success_path(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    class FakePopen:
        def read(self):
            return "16\n"

    monkeypatch.setattr(config_mod.shutil, "which", lambda _: "/usr/bin/nproc")
    monkeypatch.setattr(config_mod.os, "popen", lambda _: FakePopen())

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
    )
    assert cfg["host_threads"] == 16
    assert cfg["host_threads_minus_one"] == 15


def test_set_config_values_getpass_fallback_user_unknown(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    monkeypatch.delenv("LOGNAME", raising=False)
    monkeypatch.delenv("USER", raising=False)
    monkeypatch.setattr(
        config_mod.getpass,
        "getuser",
        lambda: (_ for _ in ()).throw(Exception("boom")),
    )

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
    )
    assert cfg["user"] == "unknown"


def test_set_config_values_arch_aarch64_maps_to_arm64(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    monkeypatch.setattr(config_mod.platform, "machine", lambda: "aarch64")

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-4.6", "genome": "b37"}
    )
    assert cfg["arch"] == "arm64"


# --------------------------
# Thread/pigz/arch extra branches
# --------------------------

def test_set_config_values_nproc_read_fails_falls_back_to_cpu_count(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    monkeypatch.setattr(config_mod.shutil, "which", lambda _: "/usr/bin/nproc")

    class BadPopen:
        def read(self):
            raise OSError("boom")

    monkeypatch.setattr(config_mod.os, "popen", lambda _: BadPopen())
    monkeypatch.setattr(config_mod.os, "cpu_count", lambda: 7)

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
    )
    assert cfg["host_threads"] == 7
    assert cfg["host_threads_minus_one"] == 6


def test_set_config_values_nproc_not_found_uses_cpu_count(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    monkeypatch.setattr(config_mod.shutil, "which", lambda _: None)
    monkeypatch.setattr(config_mod.os, "cpu_count", lambda: 5)

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
    )
    assert cfg["host_threads"] == 5
    assert cfg["host_threads_minus_one"] == 4


def test_set_config_values_pigz_missing_sets_gunzip(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-3.5",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh")

    real_access = config_mod.os.access

    def fake_access(path, mode):
        # Only pretend pigz doesn't exist; keep real exec checks for workflow scripts
        if str(path) == "/usr/bin/pigz":
            return False
        return real_access(path, mode)

    monkeypatch.setattr(config_mod.os, "access", fake_access)

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
    )
    assert cfg["compression_cmd"] == "/bin/gunzip"


def test_set_config_values_arch_other_string(monkeypatch, tmp_path):
    root = fake_project(
        monkeypatch,
        tmp_path,
        gatk_ver="gatk-4.6",
        registry_kwargs=dict(include_bash=True, bash_pipelines={"wes": {"single": "wes_single.sh"}}),
    )
    _touch_bash_files(root, "gatk-4.6", "wes_single.sh")

    monkeypatch.setattr(config_mod.platform, "machine", lambda: "mips64")

    cfg = config_mod.set_config_values(
        {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-4.6", "genome": "b37"}
    )
    assert cfg["arch"] == "mips64"


# --------------------------
# Missing file / not executable guards
# --------------------------

def test_set_config_values_missing_workflow_files_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    versions:\n"
        "      gatk-3.5:\n"
        "        common:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    # Touch only ONE file, leave others missing
    bash_dir = root / "workflows" / "bash" / "gatk-3.5"
    bash_dir.mkdir(parents=True, exist_ok=True)
    (bash_dir / "env.sh").write_text("#!/bin/sh\n", encoding="utf-8")

    with pytest.raises(WorkflowResolutionError, match="referenced in the registry do not exist"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
        )


def test_set_config_values_not_executable_bash_raises(monkeypatch, tmp_path):
    root = _mk_fake_root(monkeypatch, tmp_path)

    reg = (
        "workflows:\n"
        "  bash:\n"
        "    base_dir: \"workflows/bash\"\n"
        "    versions:\n"
        "      gatk-3.5:\n"
        "        common:\n"
        "          env: \"env.sh\"\n"
        "          coverage: \"coverage.sh\"\n"
        "          jaccard: \"jaccard.sh\"\n"
        "          vcf2sex: \"vcf2sex.sh\"\n"
        "        pipelines:\n"
        "          wes:\n"
        "            single: \"wes_single.sh\"\n"
    )
    _write_registry_and_schema(root, registry_lines=reg)

    # Create all files but NOT executable (+x missing)
    _touch_bash_files(root, "gatk-3.5", "wes_single.sh", make_exec=False)

    with pytest.raises(WorkflowResolutionError, match="Missing \\+x"):
        config_mod.set_config_values(
            {"mode": "single", "pipeline": "wes", "workflow_engine": "bash", "gatk_version": "gatk-3.5"}
        )
