#!/usr/bin/env bash
# 
#   mtDNA Pipeline Bash script.
#   This pipeline works at the the sample level, for cohorts you will 
#   need to excute "mit_cohort.sh". This way, if a new relatives comes, 
#   you cand easily add it a posteriori.
#
#   Last Modified; Jan/28/2026
#
#   Registry version: v1
#
#   Copyright (C) 2025-2026 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

set -eu

function usage {

    USAGE="""
    Usage: $0 -t n_threads

    NB1: The script expects a WES/WGS run created with export_mtdna_bam: true

MA00047_exome
└── MA0004701P_ex  <--- ID taken from here
    ├── MA0004701P_ex_S5_L001_R1_001.fastq.gz
    ├── MA0004701P_ex_S5_L001_R2_001.fastq.gz
    ├── MA0004701P_ex_S5_L002_R1_001.fastq.gz
    ├── MA0004701P_ex_S5_L002_R2_001.fastq.gz
    ├── cbicall_*_gatk-4.6_wes_single_*/04_mtdna_input/...
    └── cbicall_bash_gatk-3.5_mit_single_rsrs_* <- Submit from inside this directory
    """
    echo "$USAGE"
    exit 1
}


# Check arguments
if [ $# -eq 0 ]
 then
  usage
fi

# parsing Arguments
key="$1"
case $key in
    -t|--t)
    THREADS="$2"
esac

# Determine the directory where the script resides
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source env.sh from the same directory
source "${CBICALL_ENV_FILE:-$BINDIR/env.sh}"

# Check ARCH
if [ "$ARCH" == "aarch64" ]
 then
  echo "mit_single cannot be performed with: $ARCH"
  exit 1
fi

# Set up variables and Defining directories
DIR=$( pwd )
BINDIRMTB=$BINDIR/../../../mtdna
PYBINDIR=$BINDIR/../../../browser

id=$( echo "$DIR" | awk -F'/' '{print $(NF-1)}' | awk -F'_' '{print $1}' )
# The mtb_id needs to have this format LP6005831-???_???.bam, otherwise MToolBox will fail
mtb_id="$id-DNA_MIT"
job_id=$( echo "$DIR" | awk -F'_' '{print $NF}' )

# Set up dirs
VARCALLDIR=$DIR/01_mtoolbox
BROWSERDIR=$DIR/02_browser
mkdir "$VARCALLDIR"
mkdir "$BROWSERDIR"

# From now on we will work on VARCALL dir
cd "$VARCALLDIR"

# Prepare the BAM consumed by MToolBox
echo "Preparing Mitochondrial DNA input BAM..."

out_raw=$mtb_id.bam

export_pattern="../../*_gatk-4.6_w[ge]s_single_*/04_mtdna_input/${mtb_id}.bam"
export_list=$(ls -1 $export_pattern 2>/dev/null | grep -v 'ref_cbicall' || true)
export_count=$(printf "%s\n" "$export_list" | sed '/^$/d' | wc -l)

if [ "$export_count" -gt 1 ]; then
  echo "ERROR: More than one exported mtDNA BAM found (excluding ref_cbicall):" >&2
  printf "%s\n" "$export_list" >&2
  exit 1
elif [ "$export_count" -eq 0 ]; then
  echo "ERROR: No exported mtDNA BAM found for ID '$id'." >&2
  echo "Run WES/WGS single-sample processing with export_mtdna_bam: true first." >&2
  echo "Expected: $export_pattern" >&2
  exit 1
fi

exported_bam=$(printf "%s\n" "$export_list" | head -n 1)
exported_bai="${exported_bam}.bai"
if [ ! -s "$exported_bai" ]; then
  echo "ERROR: Exported mtDNA BAM index not found: $exported_bai" >&2
  exit 1
fi
echo "Using exported mtDNA BAM: $exported_bam"
cp "$exported_bam" "$out_raw"
cp "$exported_bai" "${out_raw}.bai"

# Performing Variant calling and annotation with MToolBox
echo "Analyzing mitochondrial DNA with MToolBox..."

(
  export PATH="$MTOOLBOXDIR:$PATH"
  export PYTHONNOUSERSITE=1

  echo "Using numpy and pandas versions:"

  # --- First choice: portable prefix python2 (non-cluster) ---
  export PATH="$PY27_PREFIX/bin:$PATH"
  export PYTHONHOME="$PY27_PREFIX"
  export PYTHONPATH="$PY27_PREFIX/lib/python2.7/site-packages"

  if python2 -c "import numpy, pandas" >/dev/null 2>&1; then
    python2 -c "import numpy, pandas; print(numpy.__version__, pandas.__version__)"
  else
    echo "Portable python2/numpy failed; trying cluster module Python2..."

    # --- Fallback: cluster module python2 ---
    unset PYTHONHOME
    module purge
    module load "${PY27_MODULE:-Python/2.7.18-GCCcore-11.2.0}"

    # Use the shipped site-packages path (same one you copy around)
    export PYTHONPATH="$PY27_PREFIX/lib/python2.7/site-packages${PYTHONPATH:+:$PYTHONPATH}"

    python2 -c "import numpy, pandas; print(numpy.__version__, pandas.__version__)" || {
      echo "ERROR: numpy/pandas not importable with either portable prefix or module python2" >&2
      python2 -c "import sys; print(sys.executable); print('\n'.join(sys.path))" >&2 || true
      exit 1
    }
  fi

  MToolBox.sh -i "$BINDIRMTB/MToolBox_config.sh" -m "-t $THREADS"
)

# We will be using the file 'prioritized_variants.txt'
# Getting GT/ DP and HF information rom VCF_file.vcf
# HF information is also in file(s) OUT*/*annotation.csv
# OUT* may contain > 1 *annotation (haplotypes), still the HF will be the same on each

# We will append the columns at the end
echo "Appending Heteroplasmic Fraction to the output..."

# ---- Configuration ----
vcf_file="VCF_file.vcf"
in_file="prioritized_variants.txt"
out_file="append_$$.txt"
final_file="mit_prioritized_variants.txt"
parse_var="$BINDIR/parse_var.pl"
# -----------------------

# Check required files
if [ ! -f "$vcf_file" ]; then
    echo "Error: File '$vcf_file' not found!" >&2
    exit 1
fi

if [ ! -f "$in_file" ]; then
    echo "Error: File '$in_file' not found!" >&2
    exit 1
fi

# Write header for appended columns
printf "REF\tALT\tGT\tDP\tHF\n" > "$out_file"

# Process variants 
# NB: Assuming each position appears once)
for var in $(cut -f1 "$in_file" | sed '1d' | "$parse_var")
do
    line=$(grep -P "chrMT\t$var\t" "$vcf_file" \
           | cut -f4,5,10 \
           | tr ':' '\t' \
           | cut -f1-5)

    if [ -n "$line" ]; then
        printf "%s\n" "$line" >> "$out_file"
    else
        printf "NA\tNA\tNA\tNA\tNA\n" >> "$out_file"
    fi
done

# Safety check: line counts must match before paste
if [ "$(wc -l < "$in_file")" -ne "$(wc -l < "$out_file")" ]; then
    echo "ERROR: Line count mismatch between input and appended columns" >&2
    echo "  $(wc -l < "$in_file") lines in $in_file" >&2
    echo "  $(wc -l < "$out_file") lines in $out_file" >&2
    rm "$out_file"
    exit 1
fi

# Merge and clean up
paste "$in_file" "$out_file" > "$final_file"
rm "$out_file"

# HTML creation
echo "Creating Browser HTML..."
mit_filtered_json="mit.filtered.json"
"$PYBINDIR/mtb2json.py" -i "$final_file" -f json > "$mit_filtered_json"
"$PYBINDIR/mtb2html.py" --id "$id" --filtered-json "$mit_filtered_json" --out "$BROWSERDIR/$job_id.html" --job-id "$job_id"

cat <<EOF > "$BROWSERDIR/README.txt"
# Open the standalone mtDNA report directly in a web browser:
xdg-open $job_id.html

# The report embeds rows derived from ../01_mtoolbox/mit.filtered.json and its JavaScript dependencies.
# Keep 01_mtoolbox and 02_browser together so its download links remain valid.
EOF

# Fin
echo "All done!!!"
exit
