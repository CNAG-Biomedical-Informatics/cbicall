#!/bin/bash
set -euo pipefail

DIR=/media/mrueda/2TBS/CNAG/Project_CBI_Call
CBICALL=/media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/bin/cbicall
NCPU=4

for DIRNAME in MA99999_exome
do
  cd "$DIR/$DIRNAME"
  echo "$DIRNAME"

  for SAMPLE in MA*ex; do
    echo "...$SAMPLE"
    cd "$SAMPLE"

    cat <<EOF > "$SAMPLE.mit_single.yaml"
mode: single
pipeline: mit
sample: $DIR/$DIRNAME/$SAMPLE
EOF

    "$CBICALL" -t "$NCPU" -p "$SAMPLE.mit_single.yaml" > "$SAMPLE.mit_single.log" 2>&1
    cd ".."
  done

  cd ".."
done
