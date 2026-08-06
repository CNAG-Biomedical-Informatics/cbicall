#!/usr/bin/env bash
set -euo pipefail

# Edit this path.
INPUT_DIR="/path/to/sample01"

cat > wgs_single.yaml <<EOF
mode: single
pipeline: wgs
workflow_provider: cbicall-core
workflow_backend: bash
software_stack: gatk-4.6
genome: hg38
resource: cbicall-germline-resources-v1
input_dir: "$INPUT_DIR"
cleanup_bam: false
EOF

cbicall run -p wgs_single.yaml -t 4
