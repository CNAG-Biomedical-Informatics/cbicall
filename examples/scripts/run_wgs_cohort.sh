#!/usr/bin/env bash
set -euo pipefail

# Edit this path. The file contains one sample ID and gVCF path per line.
SAMPLE_MAP="/path/to/sample_map.tsv"

cat > wgs_cohort.yaml <<EOF
mode: cohort
pipeline: wgs
workflow_provider: cbicall-core
workflow_backend: bash
software_stack: gatk-4.6
genome: hg38
resource: cbicall-germline-resources-v1
sample_map: "$SAMPLE_MAP"
output_basename: cohort
EOF

cbicall run -p wgs_cohort.yaml -t 4
