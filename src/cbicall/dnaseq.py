import os
import shlex
import subprocess
from pathlib import Path
from typing import Callable, Dict, List, Optional, Union

from .errors import WorkflowExecutionError, WorkflowResolutionError
from .models import RunSettings


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


def _run_cmd(
    cmd: List[str],
    cwd: Path,
    log_path: Path,
    env: Optional[Dict[str, str]] = None,
    *,
    engine: Optional[str] = None,
) -> None:
    """
    Run command and redirect stdout/stderr to the same log file.
    """
    engine_str = f"engine={engine}, " if engine else ""
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
            f"{engine_str}Command: {cmd_str}\n"
            f"Working directory: {cwd}\n"
            f"Log file: {log_path}\n"
        )
        raise WorkflowExecutionError(msg) from e

    if proc.returncode != 0:
        msg = (
            f"Failed to execute workflow (returncode={proc.returncode}).\n"
            f"{engine_str}Command: {cmd_str}\n"
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
    def engine(self) -> str:
        return str(self.settings.workflow.engine)

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
    def gatk_version(self) -> str:
        return str(self.settings.workflow.gatk_version)

    @property
    def workflow(self):
        return self.settings.workflow

    @property
    def inputs(self):
        return self.settings.inputs

    @property
    def workflow_rule(self) -> Optional[str]:
        return self.settings.workflow_rule

    @property
    def workdir(self) -> Path:
        return Path(self.settings.project_dir)

    @property
    def suffix(self) -> str:
        return f"{self.pipeline}_{self.mode}"

    @property
    def log_path(self) -> Path:
        log_name = f"{self.engine}_{self.suffix}_{self.genome}_{self.gatk_version}.log"
        return self.workdir / log_name

    def _base_env(self) -> Dict[str, str]:
        return os.environ.copy()

    def env_overrides(self) -> Optional[Dict[str, str]]:
        return None

    def build_command(self) -> List[str]:
        raise NotImplementedError

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

        self.run_cmd(
            cmd=cmd,
            cwd=self.workdir,
            log_path=self.log_path,
            env=env,
            engine=self.engine,
        )


class BashRunner(BaseRunner):
    def env_overrides(self) -> Optional[Dict[str, str]]:
        # env.sh will pick the bundle based on GENOME.
        return {"GENOME": self.genome}

    def build_command(self) -> List[str]:
        if self.workflow_rule:
            raise WorkflowResolutionError(
                "Partial workflow runs are not supported for workflow_engine='bash'."
            )

        script = self.workflow.entrypoint
        if not script:
            raise WorkflowResolutionError(f"Missing bash script for pipeline/mode '{self.suffix}'")

        cmd: List[str] = [str(script), "-t", str(int(self.settings.threads))]

        if self.gatk_version != "gatk-3.5":
            cmd += ["--pipeline", self.pipeline]

            if bool(self.settings.cleanup_bam):
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
    def build_command(self) -> List[str]:
        script = self.workflow.entrypoint
        if not script:
            raise WorkflowResolutionError(f"Missing Snakefile for pipeline/mode '{self.suffix}'")

        target_rule = self.workflow_rule or "all"

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

        snk_config_kvs: List[str] = [f"genome={self.genome}"]

        if self.gatk_version != "gatk-3.5":
            snk_config_kvs.append(f"pipeline={self.pipeline}")

            sample_map = self.inputs.sample_map
            if sample_map:
                snk_config_kvs.append(f"sample_map={sample_map}")
                snk_config_kvs.append(
                    f"workspace=cohort.genomicsdb.{self.settings.run_id}"
                )

        cmd += ["--config"] + snk_config_kvs
        return cmd


class DNAseq:
    """
    Dispatcher for engine-specific workflow runners.
    """

    def __init__(self, settings: Union[RunSettings, Dict]):
        self.settings = settings if isinstance(settings, RunSettings) else RunSettings.from_mapping(settings)

    def _make_runner(self) -> BaseRunner:
        settings = self.settings
        engine = settings.workflow.engine

        if engine == "bash":
            return BashRunner(
                settings,
                run_cmd=self._run_cmd,
                cmd_to_string=self._cmd_to_string,
            )

        if engine == "snakemake":
            return SnakemakeRunner(
                settings,
                run_cmd=self._run_cmd,
                cmd_to_string=self._cmd_to_string,
            )

        raise WorkflowResolutionError(f"Invalid workflow_engine: {engine!r}")

    def variant_calling(self) -> bool:
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
        engine: Optional[str] = None,
    ) -> None:
        _run_cmd(cmd=cmd, cwd=cwd, log_path=log_path, env=env, engine=engine)

    @staticmethod
    def _cmd_to_string(
        cmd: List[str],
        env_overrides: Optional[Dict[str, str]] = None,
    ) -> str:
        return _cmd_to_string(cmd, env_overrides=env_overrides)
