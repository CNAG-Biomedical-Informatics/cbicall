#!/bin/bash
set -euo pipefail

DIR=/media/mrueda/2TBS/CNAG/Project_CBI_Call
CBICALL=/media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/bin/cbicall
NCPU=8

for DIRNAME in MA99999_exome
do
  cd "$DIR/$DIRNAME"
  echo "$DIRNAME"

  for SAMPLE in MA*P*ex
  do
    echo "...$SAMPLE"
    cd "$SAMPLE"

    cat <<EOF > "$SAMPLE.smk_wes_single.yaml"
mode: single
pipeline: wes
workflow_engine: snakemake
gatk_version: gatk-4.6
sample: $DIR/$DIRNAME/$SAMPLE
EOF

    time "$CBICALL" -t "$NCPU" -p "$SAMPLE.smk_wes_single.yaml" > "$SAMPLE.smk_wes_single.log" 2>&1
    cd ".."
  done

  cd ".."
done
