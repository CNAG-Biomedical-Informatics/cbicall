#!/usr/bin/env bash
set -euo pipefail

# Edit this path. It must contain a completed WES/WGS run created with
# export_mtdna_bam: true.
INPUT_DIR="/path/to/sample01"

cat > mit_single.yaml <<EOF
mode: single
pipeline: mit
workflow_provider: cbicall-core
workflow_backend: bash
software_stack: gatk-3.5
genome: rsrs
resource: cbicall-germline-resources-v1
input_dir: "$INPUT_DIR"
EOF

cbicall run -p mit_single.yaml -t 4
