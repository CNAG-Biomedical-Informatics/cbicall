import os
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Optional


class DNAseq:
    """
    Thin wrapper around bash/snakemake pipelines.
    Executes the selected workflow and writes stdout/stderr to a log file.
    """

    def __init__(self, settings: Dict):
        self.__dict__.update(settings)

    def variant_calling(self) -> bool:
        pipeline: str = self.pipeline
        mode: str = self.mode
        projectdir: str = self.projectdir
        threads: int = int(self.threads)
        engine: str = self.workflow_engine
        gatk_version: str = self.gatk_version
        run_id: str = str(self.id)
        debug: bool = bool(self.debug)

        # Genome selector (default b37)
        genome: str = getattr(self, "genome", None) or "b37"

        suffix = f"{pipeline}_{mode}"
        workdir = Path(projectdir)

        if not workdir.is_dir():
            raise RuntimeError(f"Project directory does not exist: {workdir}")

        env = os.environ.copy()

        if engine == "bash":
            script = getattr(self, f"bash_{suffix}", None)
            if not script:
                raise RuntimeError(f"Missing bash script for pipeline/mode '{suffix}'")
            cmd = self._build_bash_cmd(
                script=script,
                threads=threads,
                pipeline=pipeline,
                gatk_version=gatk_version,
                run_id=run_id,
            )
            # parameters.sh will pick bundle based on GENOME
            env["GENOME"] = genome

        elif engine == "snakemake":
            script = getattr(self, f"smk_{suffix}", None)
            if not script:
                raise RuntimeError(f"Missing Snakefile for pipeline/mode '{suffix}'")
            cmd = self._build_snakemake_cmd(
                script=script,
                threads=threads,
                pipeline=pipeline,
                gatk_version=gatk_version,
                genome=genome,
            )

        else:
            raise ValueError(f"Invalid workflow_engine: {engine!r}")

        log_name = f"{engine}_{suffix}.log"
        log_path = workdir / log_name

        if debug:
            env_preview = {"GENOME": genome} if engine == "bash" else None
            print(self._cmd_to_string(cmd, env_overrides=env_preview))
            print(f"Log file: {log_path}")

        self._run_cmd(cmd=cmd, cwd=workdir, log_path=log_path, env=env)
        return True

    def _build_bash_cmd(
        self,
        script: str,
        threads: int,
        pipeline: str,
        gatk_version: str,
        run_id: str,
    ) -> List[str]:
        """
        Build the bash workflow command as an argv list (no shell).
        """
        cmd: List[str] = [script, "-t", str(threads)]

        # Keep your existing behavior for non-gatk-3.5
        if gatk_version != "gatk-3.5":
            cmd += ["--pipeline", pipeline]

            if bool(getattr(self, "cleanup_bam", False)):
                cmd.append("--cleanup-bam")

            sample_map: Optional[str] = getattr(self, "sample_map", None)
            if sample_map:
                cmd += [
                    "--sample-map",
                    sample_map,
                    "--workspace",
                    f"cohort.genomicsdb.{run_id}",
                ]

        return cmd

    def _build_snakemake_cmd(
        self,
        script: str,
        threads: int,
        pipeline: str,
        gatk_version: str,
        genome: str,
    ) -> List[str]:
        """
        Build the snakemake command as an argv list (no shell).
        """
        cmd: List[str] = [
            "snakemake",
            "--forceall",
            "all",
            "-s",
            script,
            "--cores",
            str(threads),
        ]

        # Always pass genome
        snk_config_kvs: List[str] = [f"genome={genome}"]

        # Keep your old behavior for non-gatk-3.5
        if gatk_version != "gatk-3.5":
            snk_config_kvs.append(f"pipeline={pipeline}")

        cmd += ["--config"] + snk_config_kvs
        return cmd

    @staticmethod
    def _run_cmd(
        cmd: List[str],
        cwd: Path,
        log_path: Path,
        env: Optional[Dict[str, str]] = None,
    ) -> None:
        """
        Run command and redirect stdout/stderr to the same log file.
        """
        msg = (
            "Failed to execute workflow.\n"
            f"Working directory: {cwd}\n"
            "Please check this file:\n"
            f"{log_path}\n"
        )

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
            raise RuntimeError(msg) from e

        if proc.returncode != 0:
            raise RuntimeError(msg)

    @staticmethod
    def _cmd_to_string(
        cmd: List[str],
        env_overrides: Optional[Dict[str, str]] = None,
    ) -> str:
        """
        Produce a readable shell-like string for debug prints.
        """
        parts: List[str] = []
        if env_overrides:
            for k, v in env_overrides.items():
                parts.append(f"{k}={shlex.quote(str(v))}")
        parts.extend(shlex.quote(str(x)) for x in cmd)
        return " ".join(parts)

