#!/bin/bash
set -euo pipefail

DIR=/media/mrueda/2TBS/CNAG/Project_CBI_Call
CBICALL=/media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/bin/cbicall
NCPU=4

for DIRNAME in MA99999_exome
do
  cd "$DIR/$DIRNAME"
  echo "$DIRNAME"
  echo "...$DIRNAME"

  cat <<EOF > "$DIRNAME.mit_cohort.yaml"
mode: cohort
pipeline: mit
sample: $DIR/$DIRNAME
EOF

  "$CBICALL" -t "$NCPU" -p "$DIRNAME.mit_cohort.yaml"
  cd ..
done
