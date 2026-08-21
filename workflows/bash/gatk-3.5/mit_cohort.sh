#!/usr/bin/env bash
#
#   mtDNA Pipeline Cohort Bash script.
#
#   Last Modified; Dec/29/2025
#
#   Registry version: v1
#
#   Copyright (C) 2025-2026 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

set -eu

function usage {

    USAGE="""
    Usage: $0 -t n_threads

    NB1: The script is expecting that you follow SRTI nomenclature for samples
    NB2: There is no need to run wes_cohort prior to mit_cohort.

MA00024_exome  <-- ID taken from here
├── MA0002401P_ex
│   └── cbicall_*_gatk-4.6_wes_single_*/exports/mtdna/...
├── MA0002402M_ex
│   └── cbicall_*_gatk-4.6_wes_single_*/exports/mtdna/...
└── cbicall_bash_mit_cohort_* <- Submit from inside this directory
    """
    echo "$USAGE"
    exit 1
}

# Check arguments
if [ $# -eq 0 ]; then
  usage
fi

# parsing Arguments
key="$1"
case $key in
    -t|--t)
    THREADS="${2:-}"
    ;;
    *)
    usage
    ;;
esac

if [ -z "${THREADS:-}" ]; then
  usage
fi

# Determine the directory where the script resides
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source env.sh from the same directory
source "${CBICALL_ENV_FILE:-$BINDIR/env.sh}"

# Check ARCH (same behavior as mit_single)
if [ "${ARCH:-}" = "aarch64" ]; then
  echo "mit_cohort cannot be performed with: ${ARCH:-aarch64}"
  exit 1
fi

# Set up variables and Defining directories
DIR="$(pwd)"

# Check that nomenclature exists
if [[ "$DIR" != *_bash_mit_cohort* ]]; then
  usage
fi

# Anchor project-relative directories (same as mit_single)
BINDIRMTB="$BINDIR/../../../mtdna"
PYBINDIR="$BINDIR/../../../browser"

# cohort id (format: <PROJECT>-DNA_MIT)
cohort="$(echo "$DIR" | awk -F'/' '{print $(NF-1)}' | awk -F'_' '{print $1}' | sed 's/$/-DNA_MIT/')"
echo "$cohort"

job_id="$(echo "$DIR" | awk -F'_' '{print $NF}')"

# Working dirs
VARCALLDIR="$DIR/01_mtoolbox"
BROWSERDIR="$DIR/02_browser"

mkdir -p "$VARCALLDIR"
mkdir -p "$BROWSERDIR"

cd "$VARCALLDIR"

# ------------------------------------------------------------
# Collect exported mtDNA BAMs from sibling sample dirs
# ------------------------------------------------------------

echo "Collecting exported Mitochondrial DNA input BAMs..."

# Find candidate sample directories one level up (same layout as your tree)
# Example: ../MA0002401P_ex/...
sample_dirs=$(ls -d ../../??????????_{ex,wg} 2>/dev/null || true)

if [ -z "${sample_dirs:-}" ]; then
  echo "ERROR: No sample directories matching ../??????????_{ex,wg} were found." >&2
  exit 1
fi

for sdir in $sample_dirs; do
  # Sample ID (directory name prefix before first '_')
  sid="$(basename "$sdir" | awk -F'_' '{print $1}')"
  mtb_id="${sid}-DNA_MIT"
  out_raw="${mtb_id}.bam"

  export_pattern="$sdir/"*_gatk-4.6_w[ge]s_single_*/exports/mtdna/"${mtb_id}.bam"
  export_list=$(ls -1 $export_pattern 2>/dev/null | grep -v 'ref_cbicall' || true)
  export_count=$(printf "%s\n" "$export_list" | sed '/^$/d' | wc -l)

  if [ "$export_count" -gt 1 ]; then
    echo "ERROR: More than one exported mtDNA BAM found for sample '$sid' (excluding ref_cbicall):" >&2
    printf "%s\n" "$export_list" >&2
    exit 1
  elif [ "$export_count" -eq 0 ]; then
    echo "ERROR: No exported mtDNA BAM found for sample '$sid'." >&2
    echo "Run WES/WGS single-sample processing with export_mtdna_bam: true first." >&2
    echo "Expected: $export_pattern" >&2
    exit 1
  fi

  exported_bam=$(printf "%s\n" "$export_list" | head -n 1)
  exported_bai="${exported_bam}.bai"
  if [ ! -s "$exported_bai" ]; then
    echo "ERROR: Exported mtDNA BAM index not found for sample '$sid': $exported_bai" >&2
    exit 1
  fi
  echo "Using exported mtDNA BAM for $sid: $exported_bam"
  cp "$exported_bam" "$out_raw"
  cp "$exported_bai" "${out_raw}.bai"
done

# ------------------------------------------------------------
# Run MToolBox (same python2 / numpy,pandas bootstrap as mit_single)
# ------------------------------------------------------------

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

  # Use the same config file path as mit_single (from mtdna folder)
  MToolBox.sh -i "$BINDIRMTB/MToolBox_config.sh" -m "-t $THREADS"
)

# ------------------------------------------------------------
# Post-process prioritized_variants.txt -> mit_prioritized_variants.txt
# ------------------------------------------------------------

# We will be using the file 'prioritized_variants.txt'
# Getting GT/ DP and HF information from VCF_file.vcf
# Append cohort-wide GT/DP/HF using parse_prioritized.pl

echo "Appending Heteroplasmic Fraction to the output..."

# ---- Configuration ----
vcf_file="VCF_file.vcf"
in_file="prioritized_variants.txt"
out_file="append_$$.txt"
final_file="mit_prioritized_variants.txt"
parse_var="$BINDIR/parse_var.pl"
parse_prior="$BINDIR/parse_prioritized.pl"
missing_var="missing_variants.txt"
# -----------------------

# Check if files exist
if [ ! -f "$vcf_file" ]; then
    echo "Error: File '$vcf_file' not found!"
    exit 1
fi
if [ ! -f "$in_file" ]; then
    echo "Error: File '$in_file' not found!"
    exit 1
fi

# Loop over variants
rm -f $missing_var # delete if exists
for var in $(cut -f1 "$in_file" | sed '1d' | "$parse_var")
do
  echo $var >> $missing_var
done

# Parse per-sample fields (GT/DP/HF across all samples)
"$parse_prior" -i "$vcf_file" -v $missing_var > "$out_file"

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

# ------------------------------------------------------------
# Optional: Browser HTML output (same tooling/paths as mit_single)
# ------------------------------------------------------------

echo "Creating Browser HTML..."
mit_filtered_json="mit.filtered.json"

"$PYBINDIR/mtb2json.py" -i "$final_file" -f json > "$mit_filtered_json"
"$PYBINDIR/mtb2html.py" --id "$cohort" --filtered-json "$mit_filtered_json" --out "$BROWSERDIR/$job_id.html" --job-id "$job_id"

cat <<EOF > "$BROWSERDIR/README.txt"
# Open the standalone mtDNA report directly in a web browser:
xdg-open $job_id.html

# The report embeds rows derived from ../01_mtoolbox/mit.filtered.json and its JavaScript dependencies.
# Keep 01_mtoolbox and 02_browser together so its download links remain valid.
EOF

echo "All done!!!"
exit
