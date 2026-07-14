#!/usr/bin/env bash
#
#   Relatedness determination for families
#
#   Last Modified; March/05/2025
#
#   Copyright (C) 2025-2026 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)
#
#   Using bedtools jaccard implementation

set -eu

export LC_ALL=C

# Check arguments
if [ $# -ne 1 ]
 then
  echo "$0 <workflow_backend>"
  exit 1
fi

WORKFLOW_BACKEND=$1

# Determine the directory where the script resides
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source env.sh from the same directory
source "${CBICALL_ENV_FILE:-$SCRIPT_DIR/env.sh}"

shopt -s nullglob

VCFS=(../*_{ex,wg}/cbicall_"${WORKFLOW_BACKEND}"_wes_single*/02_varcall/*.*QC.vcf)

if [ "${#VCFS[@]}" -eq 0 ]; then
  echo "No QC VCF files found for workflow backend: $WORKFLOW_BACKEND" >&2
  exit 1
fi

for VCF1 in "${VCFS[@]}"
do
 SAMPLE1=$(basename "$VCF1" .ug.QC.vcf)
 for VCF2 in "${VCFS[@]}"
 do
  SAMPLE2=$(basename "$VCF2" .ug.QC.vcf)
  echo -n "$SAMPLE1 $SAMPLE2 "
  "$BED" jaccard -a "$VCF1" -b "$VCF2" | sed '1d' | cut -f3 | tr '\n' ' '
  echo
 done
done
