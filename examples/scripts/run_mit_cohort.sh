#!/usr/bin/env bash
set -euo pipefail

# Edit this path. Every sample must contain a completed WES/WGS run created
# with export_mtdna_bam: true.
INPUT_DIR="/path/to/cohort"

cat > mit_cohort.yaml <<EOF
mode: cohort
pipeline: mit
workflow_provider: cbicall-core
workflow_backend: bash
software_stack: gatk-3.5
genome: rsrs
resource: cbicall-germline-resources-v1
input_dir: "$INPUT_DIR"
EOF

cbicall run -p mit_cohort.yaml -t 4
