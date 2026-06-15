import hashlib
import json
import os
import platform
import re
import shutil
import shlex
import subprocess
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import yaml

from .errors import WorkflowExecutionError, WorkflowResolutionError
from .models import RunSettings

EXECUTION_CONTRACT_FILE = "cbicall-execution-contract.json"


def _cmd_to_string(
    cmd: List[str],
    env_overrides: Optional[Dict[str, str]] = None,
) -> str:
    """
    Produce a readable shell-like string for debug prints and error messages.
    """
    parts: List[str] = []
    if env_overrides:
        for k, v in env_overrides.items():
            parts.append(f"{k}={shlex.quote(str(v))}")
    parts.extend(shlex.quote(str(x)) for x in cmd)
    return " ".join(parts)


def _canonical_json_sha256(value: Any) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _sha256_file(path: Path) -> Optional[str]:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _normalized_text_sha256(path: Path, replacements: Dict[str, str]) -> Optional[str]:
    if not path.is_file():
        return None
    text = path.read_text(encoding="utf-8")
    for old, new in replacements.items():
        if old:
            text = text.replace(old, new)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _normalize_contract_value(value: Any, replacements: Dict[str, str]) -> Any:
    if isinstance(value, dict):
        return {key: _normalize_contract_value(item, replacements) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalize_contract_value(item, replacements) for item in value]
    if isinstance(value, str):
        out = value
        for old, new in replacements.items():
            if old:
                out = out.replace(old, new)
        return out
    return value


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


def _groovy_single_quote(value: str) -> str:
    return "'" + str(value).replace("\\", "\\\\").replace("'", "\\'") + "'"


def _groovy_memory_value(value: object) -> str:
    text = str(value).strip()
    if re.fullmatch(r"\d+(?:\.\d+)?\.(?:B|KB|MB|GB|TB|KiB|MiB|GiB|TiB)", text):
        return text
    return _groovy_single_quote(text)


def _parameter_value_to_string(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def _append_nextflow_cli_parameters(cmd: List[str], parameters: Dict[str, object]) -> None:
    for key, value in sorted(parameters.items()):
        flag = f"--{key}"
        if value is None:
            continue
        if value is True:
            cmd.append(flag)
        else:
            cmd.extend([flag, _parameter_value_to_string(value)])


def _expand_config_value(value: Any, variables: Dict[str, Any]) -> Any:
    if isinstance(value, dict):
        return {key: _expand_config_value(item, variables) for key, item in value.items()}
    if isinstance(value, list):
        return [_expand_config_value(item, variables) for item in value]
    if not isinstance(value, str):
        return value
    out = value
    changed = True
    while changed:
        previous = out
        for key, item in variables.items():
            out = out.replace("{" + str(key) + "}", str(item))
        changed = out != previous
    return out


def _expand_native_config(config: Dict[str, Any], *, genome: str) -> Dict[str, Any]:
    data = dict(config)
    lvl1 = {"datadir": data["datadir"], "mem": data.get("mem", "8G")}
    data["dbdir"] = _expand_config_value(data["dbdir"], lvl1)
    data["ngsutils"] = _expand_config_value(data["ngsutils"], lvl1)
    data["tmpdir"] = _expand_config_value(data["tmpdir"], lvl1)
    lvl2 = {**lvl1, "dbdir": data["dbdir"], "ngsutils": data["ngsutils"], "tmpdir": data["tmpdir"]}
    data["gatk4_cmd"] = _expand_config_value(data["gatk4_cmd"], lvl2)
    resources = data.get("resources", {})
    if genome not in resources:
        raise WorkflowResolutionError(f"Missing resources.{genome} in native backend config")
    res = _expand_config_value(resources[genome], lvl2)
    bundle = res.get("bundle")
    res = _expand_config_value(res, {**lvl2, "bundle": bundle})
    res = _expand_config_value(res, {**lvl2, **res})
    data["resources"] = {genome: res}
    tools = data.get("tools", {})
    arch = "aarch64" if platform.machine() in {"aarch64", "arm64"} else "amd64"
    if arch not in tools:
        raise WorkflowResolutionError(f"Missing tools.{arch} in native backend config")
    data["selected_tools"] = _expand_config_value(tools[arch], lvl2)
    return data


def _copy_if_different(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() == dst.resolve():
        return
    shutil.copy2(src, dst)


def _run_cmd(
    cmd: List[str],
    cwd: Path,
    log_path: Path,
    env: Optional[Dict[str, str]] = None,
    *,
    backend: Optional[str] = None,
) -> None:
    """
    Run command and redirect stdout/stderr to the same log file.
    """
    backend_str = f"backend={backend}, " if backend else ""
    cmd_str = _cmd_to_string(cmd)

    log_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with log_path.open("w") as log_fh:
            proc = subprocess.run(
                cmd,
                cwd=str(cwd),
                env=env,
                stdout=log_fh,
                stderr=log_fh,
                check=False,
            )
    except Exception as e:
        msg = (
            "Failed to execute workflow (could not start command).\n"
            f"{backend_str}Command: {cmd_str}\n"
            f"Working directory: {cwd}\n"
            f"Log file: {log_path}\n"
        )
        raise WorkflowExecutionError(msg) from e

    if proc.returncode != 0:
        msg = (
            f"Failed to execute workflow (returncode={proc.returncode}).\n"
            f"{backend_str}Command: {cmd_str}\n"
            f"Working directory: {cwd}\n"
            f"Log file: {log_path}\n"
        )
        raise WorkflowExecutionError(msg)


class BaseRunner:
    def __init__(
        self,
        settings: RunSettings,
        *,
        run_cmd: Callable[..., None],
        cmd_to_string: Callable[[List[str], Optional[Dict[str, str]]], str],
    ):
        self.settings = settings
        self.run_cmd = run_cmd
        self.cmd_to_string = cmd_to_string

    @property
    def backend(self) -> str:
        return str(self.settings.workflow.backend)

    @property
    def pipeline(self) -> str:
        return str(self.settings.workflow.pipeline)

    @property
    def mode(self) -> str:
        return str(self.settings.workflow.mode)

    @property
    def genome(self) -> str:
        return str(self.settings.genome or "b37")

    @property
    def software_stack(self) -> str:
        return str(self.settings.workflow.software_stack)

    @property
    def qc_coverage_region(self) -> str:
        return str(self.settings.qc_coverage_region or "chr1")

    @property
    def workflow(self):
        return self.settings.workflow

    @property
    def inputs(self):
        return self.settings.inputs

    @property
    def workdir(self) -> Path:
        return Path(self.settings.project_dir)

    @property
    def suffix(self) -> str:
        return f"{self.pipeline}_{self.mode}"

    @property
    def log_path(self) -> Path:
        if self.workflow.metadata.get("provider") == "nf-core":
            log_name = f"nf-core_{self.suffix}.log"
        else:
            log_name = f"{self.backend}_{self.software_stack}_{self.suffix}_{self.genome}.log"
        return self.workdir / log_name

    def _base_env(self) -> Dict[str, str]:
        return os.environ.copy()

    def _coverage_env(self) -> Dict[str, str]:
        return {"CBICALL_COVERAGE_REGION": self.qc_coverage_region}

    def env_overrides(self) -> Optional[Dict[str, str]]:
        return None

    def build_command(self) -> List[str]:
        raise NotImplementedError

    def generated_execution_files(self) -> List[Dict[str, str]]:
        return []

    def _execution_contract_replacements(self) -> Dict[str, str]:
        return {
            str(self.workdir): "{PROJECT_DIR}",
            str(self.workdir.resolve()): "{PROJECT_DIR}",
            str(self.settings.run_id): "{RUN_ID}",
        }

    def _execution_file_report(self, role: str, path: Path) -> dict:
        replacements = self._execution_contract_replacements()
        return {
            "role": role,
            "path": str(path),
            "status": "present" if path.is_file() else "missing",
            "sha256": _sha256_file(path),
            "normalized_sha256": _normalized_text_sha256(path, replacements),
        }

    def _write_execution_contract(self, cmd: List[str], env_updates: Optional[Dict[str, str]]) -> Path:
        replacements = self._execution_contract_replacements()
        command_payload = {
            "argv": [str(item) for item in cmd],
            "string": self.cmd_to_string(cmd, env_overrides=env_updates),
        }
        command_payload["sha256"] = _canonical_json_sha256(command_payload)
        command_payload["normalized_sha256"] = _canonical_json_sha256(
            _normalize_contract_value(command_payload, replacements)
        )

        payload = {
            "schema_version": 1,
            "kind": "cbicall_execution_contract",
            "workflow": {
                "backend": self.backend,
                "provider": str(self.workflow.metadata.get("provider", "cbicall")),
                "pipeline": self.pipeline,
                "mode": self.mode,
                "software_stack": self.software_stack,
                "registry_version": self.workflow.registry_version,
                "key": _workflow_key(self.workflow),
                "entrypoint": self.workflow.entrypoint,
                "config_file": self.workflow.config_file,
                "helpers": dict(self.workflow.helpers),
                "metadata": dict(self.workflow.metadata),
            },
            "run": {
                "run_id": self.settings.run_id,
                "project_dir": str(self.workdir),
                "log_path": str(self.log_path),
                "threads": int(self.settings.threads),
                "profile": self.settings.profile,
                "genome": self.settings.genome,
                "qc_coverage_region": self.settings.qc_coverage_region,
                "cleanup_bam": bool(self.settings.cleanup_bam),
                "run_mode": self.settings.run_mode,
            },
            "inputs": self.inputs.to_dict(),
            "command": command_payload,
            "environment_overrides": dict(env_updates or {}),
            "backend_parameters": {
                "snakemake_parameters": dict(self.settings.snakemake_parameters),
                "nextflow_parameters": dict(self.settings.nextflow_parameters),
                "cromwell_parameters": dict(self.settings.cromwell_parameters),
                "nfcore_profile": self.settings.nfcore_profile,
                "nfcore_parameters": dict(self.settings.nfcore_parameters),
                "nfcore_singularity_cache_dir": self.settings.nfcore_singularity_cache_dir,
            },
            "generated_files": [
                self._execution_file_report(item["role"], Path(item["path"]))
                for item in self.generated_execution_files()
            ],
            "normalization": {
                "project_dir": "{PROJECT_DIR}",
                "run_id": "{RUN_ID}",
            },
        }
        payload["fingerprint"] = _canonical_json_sha256(_normalize_contract_value(payload, replacements))
        path = self.workdir / EXECUTION_CONTRACT_FILE
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return path

    def execute(self, *, debug: bool = False) -> None:
        if not self.workdir.is_dir():
            raise WorkflowResolutionError(f"Project directory does not exist: {self.workdir}")

        cmd = self.build_command()
        env = self._base_env()
        env_updates = self.env_overrides()
        if env_updates:
            env.update(env_updates)

        if debug:
            print(self.cmd_to_string(cmd, env_overrides=env_updates))
            print(f"Log file: {self.log_path}")

        self._write_execution_contract(cmd, env_updates)

        self.run_cmd(
            cmd=cmd,
            cwd=self.workdir,
            log_path=self.log_path,
            env=env,
            backend=self.backend,
        )


class BashRunner(BaseRunner):
    def env_overrides(self) -> Optional[Dict[str, str]]:
        # Bash workflows source env.sh before resolving references and helper paths.
        # Keep these as explicit handoff variables from the validated YAML contract:
        # GENOME selects the env.sh reference block; CBICALL_COVERAGE_REGION controls
        # the lightweight QC helper; CBICALL_ENV_FILE carries Bash runtime profiles.
        # Non-Bash backends use their own config/params mechanisms instead.
        env_updates = self._coverage_env()
        env_updates["GENOME"] = self.genome
        env_file = self.workflow.helpers.get("env")
        if env_file:
            env_updates["CBICALL_ENV_FILE"] = str(env_file)
        return env_updates

    def build_command(self) -> List[str]:
        # Bash reads workflow scripts directly from the repository at runtime.
        # Editing a .sh file while a job is running can make bash reach EOF early
        # and still exit 0 if the last completed command succeeded.
        script = self.workflow.entrypoint
        if not script:
            raise WorkflowResolutionError(f"Missing bash script for pipeline/mode '{self.suffix}'")

        cmd: List[str] = [str(script), "-t", str(int(self.settings.threads))]

        if self.software_stack != "gatk-3.5":
            cmd += ["--pipeline", self.pipeline]

            if bool(self.settings.cleanup_bam) and self.mode == "single":
                cmd.append("--cleanup-bam")

            sample_map = self.inputs.sample_map
            if sample_map:
                cmd += [
                    "--sample-map",
                    str(sample_map),
                    "--workspace",
                    f"cohort.genomicsdb.{self.settings.run_id}",
                ]

        return cmd


class SnakemakeRunner(BaseRunner):
    def env_overrides(self) -> Optional[Dict[str, str]]:
        return self._coverage_env()

    def build_command(self) -> List[str]:
        script = self.workflow.entrypoint
        if not script:
            raise WorkflowResolutionError(f"Missing Snakefile for pipeline/mode '{self.suffix}'")

        target_rule = str(self.settings.snakemake_parameters.get("target") or "all")

        cmd: List[str] = [
            "snakemake",
            "--forceall",
            target_rule,
            "-s",
            str(script),
            "--cores",
            str(int(self.settings.threads)),
        ]

        smk_config = self.workflow.config_file
        if smk_config:
            cmd += ["--configfile", str(smk_config)]

        snk_config_kvs: List[str] = [f"genome={self.genome}", f"qc_coverage_region={self.qc_coverage_region}"]

        if self.software_stack != "gatk-3.5":
            snk_config_kvs.append(f"pipeline={self.pipeline}")
            if self.mode == "single":
                snk_config_kvs.append(f"cleanup_bam={_parameter_value_to_string(bool(self.settings.cleanup_bam))}")

            sample_map = self.inputs.sample_map
            if sample_map:
                snk_config_kvs.append(f"sample_map={sample_map}")
                snk_config_kvs.append(
                    f"workspace=cohort.genomicsdb.{self.settings.run_id}"
                )

        for key, value in sorted(self.settings.snakemake_parameters.items()):
            if key == "target":
                continue
            snk_config_kvs.append(f"{key}={_parameter_value_to_string(value)}")

        cmd += ["--config"] + snk_config_kvs
        return cmd


class NextflowRunner(BaseRunner):
    def env_overrides(self) -> Optional[Dict[str, str]]:
        if self._is_nfcore_workflow():
            return None
        return self._coverage_env()

    def _is_nfcore_workflow(self) -> bool:
        return self.workflow.metadata.get("provider") == "nf-core"

    def _build_external_nextflow_command(self) -> List[str]:
        source = self.workflow.metadata.get("source") or self.workflow.entrypoint
        release = self.workflow.metadata.get("release")
        if not source or not release:
            raise WorkflowResolutionError("nf-core workflow registry entry must define source and release.")

        if not self.settings.nfcore_profile:
            raise WorkflowResolutionError("nfcore_profile is required for external nf-core workflows.")

        outdir_name = self.workflow.metadata.get("default_outdir", self.pipeline)
        params_file = self.workdir / "cbicall_external_nextflow.params.yaml"
        config_file = self.workdir / "cbicall_external_nextflow.config"
        max_cpus = int(self.settings.threads)
        params = dict(self.settings.nfcore_parameters)
        params["outdir"] = str((self.workdir / outdir_name).resolve())
        params["max_cpus"] = max_cpus
        params_file.write_text(yaml.safe_dump(params, sort_keys=True), encoding="utf-8")
        resource_limits = [f"cpus: {max_cpus}"]
        if params.get("max_memory"):
            resource_limits.append(f"memory: {_groovy_memory_value(params['max_memory'])}")
        config_lines = [
            "process {",
            f"  resourceLimits = [ {', '.join(resource_limits)} ]",
            "}",
        ]
        if self.settings.nfcore_singularity_cache_dir:
            config_lines.extend(
                [
                    "",
                    "singularity {",
                    f"  cacheDir = {_groovy_single_quote(self.settings.nfcore_singularity_cache_dir)}",
                    f"  libraryDir = {_groovy_single_quote(self.settings.nfcore_singularity_cache_dir)}",
                    "}",
                    "",
                    "apptainer {",
                    f"  cacheDir = {_groovy_single_quote(self.settings.nfcore_singularity_cache_dir)}",
                    f"  libraryDir = {_groovy_single_quote(self.settings.nfcore_singularity_cache_dir)}",
                    "}",
                ]
            )
        profile_tokens = {
            token.strip()
            for token in str(self.settings.nfcore_profile).split(",")
            if token.strip()
        }
        if platform.machine() in {"aarch64", "arm64"} and "docker" in profile_tokens:
            config_lines.extend(
                [
                    "",
                    "docker {",
                    "  runOptions = '--platform linux/amd64'",
                    "}",
                ]
            )
        config_file.write_text("\n".join(config_lines) + "\n", encoding="utf-8")

        return [
            "nextflow",
            "run",
            str(source),
            "-r",
            str(release),
            "-profile",
            str(self.settings.nfcore_profile),
            "-c",
            str(config_file),
            "-params-file",
            str(params_file),
            "-work-dir",
            str(self.workdir / "work"),
            "-ansi-log",
            "false",
        ]

    def build_command(self) -> List[str]:
        if self._is_nfcore_workflow():
            return self._build_external_nextflow_command()

        script = self.workflow.entrypoint
        if not script:
            raise WorkflowResolutionError(f"Missing Nextflow workflow for pipeline/mode '{self.suffix}'")

        config_file = self.workflow.config_file
        if not config_file:
            raise WorkflowResolutionError(f"Missing Nextflow params file for pipeline/mode '{self.suffix}'")

        cmd: List[str] = [
            "nextflow",
            "run",
            str(script),
            "-params-file",
            str(config_file),
            "-ansi-log",
            "false",
            "--pipeline",
            self.pipeline,
            "--genome",
            self.genome,
            "--threads",
            str(int(self.settings.threads)),
            "--qc_coverage_region",
            self.qc_coverage_region,
        ]
        if self.mode == "single":
            cmd += ["--cleanup_bam", "true" if bool(self.settings.cleanup_bam) else "false"]

        sample_map = self.inputs.sample_map
        if self.mode == "cohort":
            if not sample_map:
                raise WorkflowResolutionError(
                    "sample_map is required for workflow_backend='nextflow' with mode='cohort'."
                )
            cmd += [
                "--sample_map",
                str(sample_map),
                "--workspace",
                f"cohort.genomicsdb.{self.settings.run_id}",
            ]

        helper_params = {
            "coverage": "coverage_script",
            "vcf2sex": "vcf2sex_script",
            "vcf2hash": "vcf2hash_script",
        }
        for helper_name, param_name in helper_params.items():
            helper_path = self.workflow.helpers.get(helper_name)
            if helper_path:
                cmd += [f"--{param_name}", str(helper_path)]

        _append_nextflow_cli_parameters(cmd, self.settings.nextflow_parameters)
        return cmd

    def generated_execution_files(self) -> List[Dict[str, str]]:
        if self._is_nfcore_workflow():
            return [
                {"role": "nf-core:params", "path": str(self.workdir / "cbicall_external_nextflow.params.yaml")},
                {"role": "nf-core:config", "path": str(self.workdir / "cbicall_external_nextflow.config")},
            ]
        return []


class CromwellRunner(BaseRunner):
    def env_overrides(self) -> Optional[Dict[str, str]]:
        return self._coverage_env()

    @property
    def workflow_name(self) -> str:
        return "CBIcallCohort" if self.mode == "cohort" else "CBIcallWesSingle"

    @property
    def inputs_json(self) -> Path:
        return self.workdir / "cbicall_cromwell.inputs.json"

    @property
    def options_json(self) -> Path:
        return self.workdir / "cbicall_cromwell.options.json"

    @property
    def metadata_json(self) -> Path:
        return self.workdir / "cbicall_cromwell.metadata.json"

    @property
    def fastq_pairs_tsv(self) -> Path:
        return self.workdir / "cbicall_cromwell.fastq_pairs.tsv"

    def _sample_id(self) -> str:
        rawid = self.workdir.parent.name
        return rawid.split("_", 1)[0]

    def _input_dir(self) -> Path:
        if self.inputs.input_dir:
            return Path(self.inputs.input_dir)
        return self.workdir.parent

    def _sample_map(self) -> Path:
        if not self.inputs.sample_map:
            raise WorkflowResolutionError(
                "sample_map is required for workflow_backend='cromwell' with mode='cohort'."
            )
        sample_map = Path(self.inputs.sample_map)
        if not sample_map.is_file():
            raise WorkflowResolutionError(f"sample_map does not exist for Cromwell cohort input: {sample_map}")
        return sample_map.resolve()

    def _write_fastq_pairs(self) -> None:
        input_dir = self._input_dir()
        pairs = []
        for r1 in sorted(input_dir.glob("*_R1_*fastq.gz")):
            r2 = Path(str(r1).replace("_R1_", "_R2_"))
            if not r2.is_file():
                raise WorkflowResolutionError(f"Missing FASTQ mate for Cromwell input: {r2}")
            base = r1.name.replace(".fastq.gz", "").split("_R1_", 1)[0]
            pairs.append((base, str(r1.resolve()), str(r2.resolve())))
        if not pairs:
            raise WorkflowResolutionError(f"No *_R1_*fastq.gz files found for Cromwell input in: {input_dir}")
        self.fastq_pairs_tsv.write_text(
            "".join("\t".join(row) + "\n" for row in pairs),
            encoding="utf-8",
        )

    def _prepare_run_directory(self) -> None:
        self.workdir.mkdir(parents=True, exist_ok=True)
        output_dirs = (
            "01_bam",
            "01_genomicsdb",
            "02_varcall",
            "03_stats",
            "logs",
            "cromwell-work",
            "cromwell-outputs",
            "cromwell-logs",
        )
        for rel in output_dirs:
            (self.workdir / rel).mkdir(parents=True, exist_ok=True)

    def _expanded_native_config(self) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
        if not self.workflow.config_file:
            raise WorkflowResolutionError(f"Missing Cromwell config file for pipeline/mode '{self.suffix}'")
        config = yaml.safe_load(Path(self.workflow.config_file).read_text(encoding="utf-8")) or {}
        expanded = _expand_native_config(config, genome=self.genome)
        resources = expanded["resources"][self.genome]
        tools = expanded["selected_tools"]
        return expanded, resources, tools

    def _gatk4_cmd_for_mode(self, config: Dict[str, Any]) -> str:
        if self.mode != "cohort" or not self.workflow.config_file:
            return str(config["gatk4_cmd"])
        raw_config = yaml.safe_load(Path(self.workflow.config_file).read_text(encoding="utf-8")) or {}
        if "mem_genotype" not in raw_config:
            return str(config["gatk4_cmd"])
        genotype_config = dict(raw_config)
        genotype_config["mem"] = raw_config["mem_genotype"]
        return str(_expand_native_config(genotype_config, genome=self.genome)["gatk4_cmd"])

    def _base_payload(self, prefix: str) -> dict:
        expanded, resources, tools = self._expanded_native_config()
        return {
            prefix + "pipeline": self.pipeline,
            prefix + "genome": self.genome,
            prefix + "threads": int(self.settings.threads),
            prefix + "qc_coverage_region": self.qc_coverage_region,
            prefix + "tmpdir": str(expanded["tmpdir"]),
            prefix + "gatk4_cmd": self._gatk4_cmd_for_mode(expanded),
            prefix + "ref": str(resources["ref"]),
            prefix + "ref_dict": str(resources["ref_dict"]),
            prefix + "dbsnp": str(resources["dbsnp"]),
            prefix + "mills_indels": str(resources["mills_indels"]),
            prefix + "hapmap": str(resources["hapmap"]),
            prefix + "omni": str(resources["omni"]),
            prefix + "interval_list": str(resources["interval_list"]),
            prefix + "snp_res": str(resources["snp_res"]),
            prefix + "indel_res": str(resources["indel_res"]),
            "_cbicall_expanded_config": expanded,
            "_cbicall_resources": resources,
            "_cbicall_tools": tools,
        }

    def _write_cromwell_json(self) -> None:
        if self.pipeline not in {"wes", "wgs"} or self.mode not in {"single", "cohort"}:
            raise WorkflowResolutionError("Cromwell backend supports WES/WGS single and cohort modes only.")

        self._prepare_run_directory()
        prefix = self.workflow_name + "."
        payload = self._base_payload(prefix)
        expanded = payload.pop("_cbicall_expanded_config")
        resources = payload.pop("_cbicall_resources")
        tools = payload.pop("_cbicall_tools")

        if self.mode == "single":
            self._write_fastq_pairs()
            payload.update(
                {
                    prefix + "id": self._sample_id(),
                    prefix + "cleanup_bam": bool(self.settings.cleanup_bam),
                    prefix + "fastq_pairs_tsv": str(self.fastq_pairs_tsv.resolve()),
                    prefix + "bwa": str(tools["bwa"]),
                    prefix + "samtools": str(tools["samtools"]),
                    prefix + "refgz": str(resources["refgz"]),
                    prefix + "kg_indels": str(resources["kg_indels"]),
                    prefix + "coverage_script": str(self.workflow.helpers["coverage"]),
                    prefix + "vcf2sex_script": str(self.workflow.helpers["vcf2sex"]),
                    prefix + "vcf2hash_script": str(self.workflow.helpers["vcf2hash"]),
                }
            )
        else:
            payload.update(
                {
                    prefix + "sample_map": str(self._sample_map()),
                    prefix + "workspace": f"cohort.genomicsdb.{self.settings.run_id}",
                    prefix + "vcf2hash_script": str(self.workflow.helpers["vcf2hash"]),
                    prefix + "min_snp_for_vqsr": int(expanded.get("min_snp_for_vqsr", 1000)),
                    prefix + "min_indel_for_vqsr": int(expanded.get("min_indel_for_vqsr", 8000)),
                }
            )

        for key, value in sorted(self.settings.cromwell_parameters.items()):
            payload[key if "." in key else prefix + key] = value
        self.inputs_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        options = {
            "final_workflow_outputs_dir": str((self.workdir / "cromwell-outputs").resolve()),
            "final_workflow_log_dir": str((self.workdir / "cromwell-logs").resolve()),
            "use_relative_output_paths": True,
        }
        self.options_json.write_text(json.dumps(options, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    def _cromwell_command_prefix(self) -> List[str]:
        jar = os.environ.get("CROMWELL_JAR")
        if jar:
            jar_path = Path(jar)
            if not jar_path.is_file():
                raise WorkflowResolutionError(f"CROMWELL_JAR does not exist: {jar}")
            return [os.environ.get("JAVA_CMD", "java"), "-jar", str(jar_path)]
        executable = shutil.which("cromwell")
        if executable:
            return [executable]
        raise WorkflowResolutionError("Cromwell is not available. Set CROMWELL_JAR or put cromwell on PATH.")

    def build_command(self) -> List[str]:
        script = self.workflow.entrypoint
        if not script:
            raise WorkflowResolutionError(f"Missing Cromwell WDL for pipeline/mode '{self.suffix}'")
        self._write_cromwell_json()
        return [
            *self._cromwell_command_prefix(),
            "run",
            str(script),
            "--inputs",
            str(self.inputs_json),
            "--options",
            str(self.options_json),
            "--metadata-output",
            str(self.metadata_json),
        ]

    def generated_execution_files(self) -> List[Dict[str, str]]:
        files = [
            {"role": "cromwell:inputs", "path": str(self.inputs_json)},
            {"role": "cromwell:options", "path": str(self.options_json)},
        ]
        if self.fastq_pairs_tsv.is_file():
            files.append({"role": "cromwell:fastq_pairs", "path": str(self.fastq_pairs_tsv)})
        return files

    def execute(self, *, debug: bool = False) -> None:
        super().execute(debug=debug)
        self._promote_outputs()

    def _promote_outputs(self) -> None:
        if not self.metadata_json.is_file():
            raise WorkflowExecutionError(f"Cromwell metadata output was not created: {self.metadata_json}")
        metadata = json.loads(self.metadata_json.read_text(encoding="utf-8"))
        outputs = metadata.get("outputs") or {}
        destinations = {
            "gvcf": "02_varcall",
            "raw_vcf": "02_varcall",
            "qc_vcf": "02_varcall",
            "coverage": "03_stats",
            "sex": "03_stats",
            "vcf_hash": "03_stats",
            "logs": "logs",
        }
        for key, value in outputs.items():
            suffix = str(key).split(".")[-1]
            target_dir = destinations.get(suffix)
            if not target_dir:
                continue
            values = value if isinstance(value, list) else [value]
            for item in values:
                if not item:
                    continue
                src = Path(str(item))
                if src.is_file():
                    _copy_if_different(src, self.workdir / target_dir / src.name)


class WorkflowExecutor:
    """
    Backend-specific workflow runner facade.
    """

    def __init__(self, settings: Union[RunSettings, Dict]):
        self.settings = settings if isinstance(settings, RunSettings) else RunSettings.from_mapping(settings)

    def _make_runner(self) -> BaseRunner:
        settings = self.settings
        backend = settings.workflow.backend

        if backend == "bash":
            return BashRunner(
                settings,
                run_cmd=self._run_cmd,
                cmd_to_string=self._cmd_to_string,
            )

        if backend == "snakemake":
            return SnakemakeRunner(
                settings,
                run_cmd=self._run_cmd,
                cmd_to_string=self._cmd_to_string,
            )

        if backend == "nextflow":
            return NextflowRunner(
                settings,
                run_cmd=self._run_cmd,
                cmd_to_string=self._cmd_to_string,
            )

        if backend == "cromwell":
            return CromwellRunner(
                settings,
                run_cmd=self._run_cmd,
                cmd_to_string=self._cmd_to_string,
            )

        raise WorkflowResolutionError(f"Invalid workflow_backend: {backend!r}")

    def run(self) -> bool:
        runner = self._make_runner()
        runner.execute(debug=bool(self.settings.debug))
        return True

    @staticmethod
    def _run_cmd(
        cmd: List[str],
        cwd: Path,
        log_path: Path,
        env: Optional[Dict[str, str]] = None,
        *,
        backend: Optional[str] = None,
    ) -> None:
        _run_cmd(cmd=cmd, cwd=cwd, log_path=log_path, env=env, backend=backend)

    @staticmethod
    def _cmd_to_string(
        cmd: List[str],
        env_overrides: Optional[Dict[str, str]] = None,
    ) -> str:
        return _cmd_to_string(cmd, env_overrides=env_overrides)
