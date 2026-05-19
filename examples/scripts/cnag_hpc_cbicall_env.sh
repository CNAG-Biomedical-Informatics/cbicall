#!/bin/bash
#
# CNAG HPC bootstrap for running CBIcall from a source checkout.
#
# Source this file before calling bin/cbicall from a SLURM job or interactive
# shell. The Python module is needed before CBIcall starts; Nextflow is loaded
# here as a convenience for nf-core workflows.

module load Python/3.10.8-GCCcore-12.2.0
export PYTHONPATH="/software/biomed/cbi_py3/lib/python3.10/site-packages:${PYTHONPATH:-}"

if command -v module >/dev/null 2>&1; then
  module load Nextflow/25.10.2 2>/dev/null || true
fi

# Set these to a user- or project-owned directory when running nf-core with
# Apptainer/Singularity to avoid unreadable shared-cache images. These are
# shell/HPC concerns; CBIcall inherits them but does not create them.
# export NXF_SINGULARITY_CACHEDIR=/path/to/singularity/cache
# export NXF_SINGULARITY_LIBRARYDIR=/path/to/singularity/library
