#!/bin/bash
set -euo pipefail

DIR=/media/mrueda/2TBS/CNAG/Project_CBI_Call
CBICALL=/media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/bin/cbicall
NCPU=4

for DIRNAME in MA99999_exome
do
  cd "$DIR/$DIRNAME"
  echo "$DIRNAME"

  for SAMPLE in MA*P*ex
  do
    echo "...$SAMPLE"
    cd "$SAMPLE"

    cat <<EOF > "$SAMPLE.smk_wgs_single.yaml"
mode: single
pipeline: wgs
gatk_version: gatk-4.6
sample: $DIR/$DIRNAME/$SAMPLE
EOF

    time "$CBICALL" -t "$NCPU" -p "$SAMPLE.wgs_single.yaml" > "$SAMPLE.wgs_single.log" 2>&1
    cd ".."
  done

  cd ".."
done
