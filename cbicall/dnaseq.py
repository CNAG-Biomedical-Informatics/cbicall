import os
import subprocess
from pathlib import Path
from typing import Dict


class DNAseq:
    """
    Thin wrapper around bash/snakemake pipelines.
    Mirrors Perl DNAseq package.
    """

    def __init__(self, settings: Dict):
        # Copy all keys into the instance dict (like bless $self, $class)
        self.__dict__.update(settings)

    def variant_calling(self) -> bool:
        pipeline = self.pipeline
        mode = self.mode
        projectdir = self.projectdir
        threads = self.threads
        engine = self.workflow_engine
        gatk_version = self.gatk_version
        run_id = self.id
        debug = self.debug

        suffix = f"{pipeline}_{mode}"

        if engine == "bash":
            script = getattr(self, f"bash_{suffix}", None)
            if not script:
                raise RuntimeError(f"Missing bash script for mode '{mode}'")
            cmd_base = f"{script} -t {threads}"
        elif engine == "snakemake":
            script = getattr(self, f"smk_{suffix}", None)
            if not script:
                raise RuntimeError(f"Missing Snakefile for mode '{mode}'")
            cmd_base = f"snakemake --forceall all -s {script} --cores {threads}"
        else:
            raise ValueError(f"Invalid workflow_engine: {engine}")

        # Add pipeline/config flag only for non-gatk-3.5
        if gatk_version != "gatk-3.5":
            if engine == "bash":
                cmd_base += f" --pipeline {pipeline}"

                # cleanup_bam flag
                if getattr(self, "cleanup_bam", False):
                    cmd_base += " --cleanup-bam"

                # cohort-specific flags (only if sample_map present)
                if getattr(self, "sample_map", None):
                    cmd_base += (
                        f" --sample-map {self.sample_map}"
                        f" --workspace cohort.genomicsdb.{run_id}"
                    )
            else:
                cmd_base += f" --config pipeline={pipeline}"

        if debug:
            print(cmd_base)

        log_name = f"{engine}_{suffix}.log"
        log_path = Path(projectdir) / log_name
        cmd = f"cd {projectdir} && {cmd_base} > {log_name} 2>&1"
        self._submit_cmd(cmd, log_path, run_id, debug)
        return True

    @staticmethod
    def _submit_cmd(cmd: str, log_path: Path, run_id: str, debug: bool) -> None:
        msg = (
            f"Failed to execute: {run_id}\n"
            f"Please check this file:\n{log_path}\n"
        )

        # Mirror Perl 'system' behaviour
        try:
            rc = subprocess.call(cmd, shell=True)
        except Exception as e:
            raise RuntimeError(msg) from e

        if rc != 0:
            if debug:
                # analog of confess: include stack trace
                raise RuntimeError(msg)
            else:
                raise RuntimeError(msg)
