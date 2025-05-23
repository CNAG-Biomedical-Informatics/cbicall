#!/usr/bin/env bash
#
#   Relatedness determination for families
#
#   Last Modified; March/05/2025
#
#   $VERSION taken from CBICall
#
#   Copyright (C) 2025 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)
#
#   Using bedtools jaccard implementation

set -eu

export LC_ALL=C

# Check arguments
if [ $# -ne 1 ]
 then
  echo "$0 <workflow_engine>"
  exit 1
fi

workflow_engine=$1

# Determine the directory where the script resides
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source parameters.sh from the same directory
source "$BINDIR/parameters.sh"

dir=../*_ex/cbicall_${workflow_engine}_wes_single*/02_varcall

for vcf1 in $( ls -1 "$dir"/*.*QC.vcf )
do
 short1=$( echo "$vcf1" | awk -F'/' '{print $NF}' | sed 's/.ug.QC.vcf//' )
 for vcf2 in $( ls -1 "$dir"/*.*QC.vcf )
 do
  short2=$( echo "$vcf2" | awk -F'/' '{print $NF}' | sed 's/.ug.QC.vcf//' )
  echo -n "$short1 $short2 "
  $BED jaccard -a "$vcf1" -b "$vcf2" | sed '1d' | cut -f3 | tr '\n' ' '
  echo
 done
done
