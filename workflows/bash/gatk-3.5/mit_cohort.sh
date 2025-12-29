#!/usr/bin/env bash
#
#   mtDNA Pipeline Cohort Bash script.
#
#   Last Modified; Dec/29/2025
#
#   $VERSION taken from CBICall
#
#   Copyright (C) 2025 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

set -eu

function usage {

    USAGE="""
    Usage: $0 -t n_threads

    NB1: The script is expecting that you follow SRTI nomenclature for samples
    NB2: There is no need to run wes_cohort prior to mit_cohort.

MA00024_exome  <-- ID taken from here
├── MA0002401P_ex
│   └── cbicall_*_wes_single_*gatk-*/...
├── MA0002402M_ex
│   └── cbicall_*_wes_single_*gatk-*/...
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

# Source parameters.sh from the same directory
source "$BINDIR/parameters.sh"

# Check ARCH (same behavior as mit_single)
if [ "${ARCH:-}" = "aarch64" ]; then
  echo "mit_cohort cannot be performed with: ${ARCH:-aarch64}"
  exit 1
fi

# Set up variables and Defining directories
DIR="$(pwd)"

# Check that nomenclature exists
if [[ "$DIR" != *cbicall_bash_mit_cohort* ]]; then
  usage
fi

# Anchor project-relative directories (same as mit_single)
BINDIRMTB="$BINDIR/../../../mtdna"
PYBINDIR="$BINDIR/../../../browser"
ASSETS="$PYBINDIR/assets"

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
# Extract chrM from each sample BAM found in sibling sample dirs
# ------------------------------------------------------------

echo "Extracting Mitochondrial DNA from exome BAM files..."

# Determine chrM naming based on REF
if [[ "$REF" == *b37*.fasta* ]]; then
  chrM="MT"
else
  chrM="chrM"
fi

# Find candidate sample directories one level up (same layout as your tree)
# Example: ../MA0002401P_ex/...
sample_dirs=$(ls -d ../../??????????_ex 2>/dev/null || true)

if [ -z "${sample_dirs:-}" ]; then
  echo "ERROR: No sample directories matching ../??????????_ex were found." >&2
  exit 1
fi

found_any=0

for sdir in $sample_dirs; do
  # Sample ID (directory name prefix before first '_')
  sid="$(basename "$sdir" | awk -F'_' '{print $1}')"
  mtb_id="${sid}-DNA_MIT"

  bam_raw=""

  # --- Prefer GATK 3.5 bam layout (wes_single, fixed.bam) ---
  p35="$sdir/"*cbicall_bash_wes_single_*gatk-3.5*/01_bam/input.merged.filtered.realigned.fixed.bam
  list35=$(ls -1 $p35 2>/dev/null | grep -v 'ref_cbicall' || true)
  n35=$(printf "%s\n" "$list35" | sed '/^$/d' | wc -l)

  if [ "$n35" -gt 1 ]; then
    echo "ERROR: More than one GATK 3.5 BAM found for sample '$sid' (excluding ref_cbicall):" >&2
    printf "%s\n" "$list35" >&2
    exit 1
  elif [ "$n35" -eq 1 ]; then
    bam_raw=$(printf "%s\n" "$list35" | head -n 1)
    echo "Using GATK 3.5 BAM for $sid: $bam_raw"
  fi

  # --- Otherwise try GATK 4.6 bam layout (wes/wgs single, recal.bam) ---
  if [ -z "$bam_raw" ]; then
    p46="$sdir/"*cbicall_bash_w[ge]s_single_*gatk-4.6*/01_bam/"$sid".rg.merged.dedup.recal.bam
    list46=$(ls -1 $p46 2>/dev/null | grep -v 'ref_cbicall' || true)
    n46=$(printf "%s\n" "$list46" | sed '/^$/d' | wc -l)

    if [ "$n46" -gt 1 ]; then
      echo "ERROR: More than one GATK 4.6 BAM found for sample '$sid' (excluding ref_cbicall):" >&2
      printf "%s\n" "$list46" >&2
      exit 1
    elif [ "$n46" -eq 1 ]; then
      bam_raw=$(printf "%s\n" "$list46" | head -n 1)
      echo "Using GATK 4.6 BAM for $sid: $bam_raw"
    fi
  fi

  if [ -z "$bam_raw" ]; then
    echo "WARNING: Could not find BAM for sample '$sid' (excluding ref_cbicall). Skipping." >&2
    continue
  fi

  found_any=1

  # Ensure index exists as *.bam.bai (MToolBox/samtools expectations are happier)
  bam_raw_index="${bam_raw%.bam}.bai"
  bam_raw_index_ok="${bam_raw}.bai"
  if [ -s "$bam_raw_index" ] && [ ! -s "$bam_raw_index_ok" ]; then
    cp "$bam_raw_index" "$bam_raw_index_ok"
  fi

  out_raw="${mtb_id}.bam"

  "$SAM" view -b "$bam_raw" "$chrM" > "$out_raw"
  "$SAM" index "$out_raw"
done

if [ "$found_any" -ne 1 ]; then
  echo "ERROR: No usable sample BAMs found. Nothing to do." >&2
  exit 1
fi

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

vcf_file="VCF_file.vcf"
vcf_tmp="VCF_file_$$.vcf"
in_file="prioritized_variants.txt"
out_file="append_$$.txt"
final_file="mit_prioritized_variants.txt"

parse_var="$BINDIR/parse_var.pl"
parse_prior="$BINDIR/parse_prioritized.pl"

# Check if files exist
if [ ! -f "$vcf_file" ]; then
    echo "Error: File '$vcf_file' not found!"
    exit 1
fi
if [ ! -f "$in_file" ]; then
    echo "Error: File '$in_file' not found!"
    exit 1
fi

# Keep header
grep '^#CHROM' "$vcf_file" > "$vcf_tmp"

# Add only variants of interest
for var in $(cut -f1 "$in_file" | sed '1d' | "$parse_var")
do
  grep -P "chrMT\t$var\t" "$vcf_file" >> "$vcf_tmp"  || echo "$var not found"
done

# Parse per-sample fields (GT/DP/HF across all samples)
"$parse_prior" -i "$vcf_tmp" > "$out_file"

# Merge onto prioritized_variants
paste "$in_file" "$out_file" > "$final_file"
rm "$vcf_tmp" "$out_file"

# ------------------------------------------------------------
# Optional: Browser HTML output (same tooling/paths as mit_single)
# ------------------------------------------------------------

echo "Creating Browser HTML..."
mit_json="mit.json"
mit_raw_json="mit.raw.json"

"$PYBINDIR/mtb2json.py" -i "$final_file" -f json > "$mit_raw_json"
"$PYBINDIR/mtb2json.py" -i "$final_file" -f json4html > "$BROWSERDIR/$mit_json"
"$PYBINDIR/mtb2html.py" --id "$cohort" --json "$mit_json" --out "$BROWSERDIR/$job_id.html" --job-id "$job_id"
ln -s "$ASSETS" "$BROWSERDIR/assets" 2>/dev/null || true

cat <<EOF > "$BROWSERDIR/README.txt"
# To visualize <$job_id.html>:

# Option 1: Open <$job_id.html> directly in Chromium
chromium --allow-file-access-from-files --disable-web-security $job_id.html

# Option 2: Use an HTTP server. Example using Python 3:
python3 -m http.server
EOF

echo "All done!!!"
exit

