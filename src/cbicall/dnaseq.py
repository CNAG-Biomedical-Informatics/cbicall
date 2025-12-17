import os
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Optional


class DNAseq:
    """
    Thin wrapper around bash/snakemake pipelines.
    """

    def __init__(self, settings: Dict):
        # Copy all keys into the instance dict (like bless $self, $class)
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

        # New: genome selection (default b37)
        genome: str = getattr(self, "genome", None) or "b37"

        suffix = f"{pipeline}_{mode}"
        workdir = Path(projectdir)

        if not workdir.is_dir():
            raise RuntimeError(f"Project directory does not exist: {workdir}")

        # Build command
        if engine == "bash":
            script = getattr(self, f"bash_{suffix}", None)
            if not script:
                raise RuntimeError(f"Missing bash script for pipeline/mode '{suffix}'")
            cmd = self._build_bash_cmd(script=script, threads=threads, pipeline=pipeline, gatk_version=gatk_version, run_id=run_id)
            env = os.environ.copy()
            env["GENOME"] = genome  # bash parameters.sh will use this
        elif engine == "snakemake":
            script = getattr(self, f"smk_{suffix}", None)
            if not script:
                raise RuntimeError(f"Missing Snakefile for pipeline/mode '{suffix}'")
            cmd = self._build_snakemake_cmd(script=script, threads=threads, pipeline=pipeline, gatk_version=gatk_version, genome=genome)
            env = os.environ.copy()
        else:
            raise ValueError(f"Invalid workflow_engine: {engine!r}")

        log_name = f"{engine}_{suffix}.log"
        log_path = workdir / log_name

        if debug:
            # Print a shell-like representation for debugging
            print(self._cmd_to_string(cmd, env_overrides={"GENOME": genome} if engine == "bash" else None))
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

        # Add pipeline/config flag only for non-gatk-3.5 (kept from your original logic)
        if gatk_version != "gatk-3.5":
            cmd += ["--pipeline", pipeline]

            # cleanup_bam flag
            if bool(getattr(self, "cleanup_bam", False)):
                cmd.append("--cleanup-bam")

            # cohort-specific flags (only if sample_map present)
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

        # Always pass genome so the Snakefile can select the correct bundle
        snk_config_kvs: List[str] = [f"genome={genome}"]

        # Keep your original behavior for non-gatk-3.5
        if gatk_version != "gatk-3.5":
            snk_config_kvs.append(f"pipeline={pipeline}")

        if snk_config_kvs:
            cmd += ["--config"] + snk_config_kvs

        return cmd

    @staticmethod
    def _run_cmd(cmd: List[str], cwd: Path, log_path: Path, env: Optional[Dict[str, str]] = None) -> None:
        """
        Run a command and redirect stdout/stderr to a log file.
        """
        msg = (
            f"Failed to execute workflow.\n"
            f"Working directory: {cwd}\n"
            f"Please check this file:\n{log_path}\n"
        )

        log_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            with log_path.open("w") as log_fh:
                proc = subprocess.run(
                    cmd,
                    cwd=str(cwd),
                    env=env,
                    stdout=log_fh,
                    stderr=log

