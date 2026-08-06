#!/usr/bin/env bash
set -euo pipefail

# Edit this path.
INPUT_DIR="/path/to/sample01"

cat > wes_single.yaml <<EOF
mode: single
pipeline: wes
workflow_provider: cbicall-core
workflow_backend: bash
software_stack: gatk-4.6
genome: b37
resource: cbicall-germline-resources-v1
input_dir: "$INPUT_DIR"
cleanup_bam: false
EOF

cbicall run -p wes_single.yaml -t 4
