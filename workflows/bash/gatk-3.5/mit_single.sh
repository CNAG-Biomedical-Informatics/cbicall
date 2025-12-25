#!/usr/bin/env bash
# 
#   mtDNA Pipeline Bash script.
#   This pipeline works at the the sample level, for cohorts you will 
#   need to excute "mit_cohort.sh". This way, if a new relatives comes, 
#   you cand easily add it a posteriori.
#
#   Last Modified; March/05/2025
#
#   $VERSION taken from CBICall
#
#   Copyright (C) 2025 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

set -eu

function usage {

    USAGE="""
    Usage: $0 -t n_threads

    NB1: The script is expecting that you follow SRTI nomenclature for samples

MA00047_exome
└── MA0004701P_ex  <--- ID taken from here
    ├── MA0004701P_ex_S5_L001_R1_001.fastq.gz
    ├── MA0004701P_ex_S5_L001_R2_001.fastq.gz
    ├── MA0004701P_ex_S5_L002_R1_001.fastq.gz
    ├── MA0004701P_ex_S5_L002_R2_001.fastq.gz
    └── cbicall_bash_mit_single_gatk-3.5_* <- The script expects that you are submitting the job from inside this directory
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

# Source parameters.sh from the same directory
source "$BINDIR/parameters.sh"

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
ASSETS=$PYBINDIR/assets

id=$( echo "$DIR" | awk -F'/' '{print $(NF-1)}' | awk -F'_' '{print $1}' )
# The mtb_id needs to have this format LP6005831-???_???.bam, otherwise MToolBox will fail
mtb_id="$id-DNA_MIT"
job_id=$( echo "$DIR" | awk -F'_' '{print $NF}' )

# Set up dirs
VARCALLDIR=$DIR/01_mtoolbox
BROWSERDIR=$DIR/02_browser
mkdir "$VARCALLDIR"
mkdir $BROWSERDIR

# From now on we will work on VARCALL dir
cd "$VARCALLDIR"

# Using Samtools to extract chrM
echo "Extracting Mitochondrial DNA from exome BAM file..."

out_raw=$mtb_id.bam

bam_raw=""

p35='../../*cbicall_bash_wes_single_*gatk-3.5*/01_bam/input.merged.filtered.realigned.fixed.bam'
list35=$(ls -1 $p35 2>/dev/null | grep -v 'ref_cbicall' || true)
n35=$(printf "%s\n" "$list35" | sed '/^$/d' | wc -l)

if [ "$n35" -gt 1 ]; then
  echo "ERROR: More than one GATK 3.5 BAM found (excluding ref_cbicall):" >&2
  printf "%s\n" "$list35" >&2
  exit 1
elif [ "$n35" -eq 1 ]; then
  bam_raw=$(printf "%s\n" "$list35" | head -n 1)
  echo "Using GATK 3.5 BAM: $bam_raw"
fi

if [ -z "$bam_raw" ]; then
  p46="../../*cbicall_bash_w[ge]s_single_*gatk-4.6*/01_bam/${id}.rg.merged.dedup.recal.bam"
  list46=$(ls -1 $p46 2>/dev/null | grep -v 'ref_cbicall' || true)
  n46=$(printf "%s\n" "$list46" | sed '/^$/d' | wc -l)

  if [ "$n46" -gt 1 ]; then
    echo "ERROR: More than one GATK 4.6 BAM found (excluding ref_cbicall):" >&2
    printf "%s\n" "$list46" >&2
    exit 1
  elif [ "$n46" -eq 1 ]; then
    bam_raw=$(printf "%s\n" "$list46" | head -n 1)
    echo "Using GATK 4.6 BAM: $bam_raw"
  fi
fi

if [ -z "$bam_raw" ]; then
  echo "ERROR: Could not find BAM for ID '$id' (excluding ref_cbicall) in either:" >&2
  echo "  $p35" >&2
  echo "  $p46" >&2
  exit 1
fi

BAMDIR=$(dirname "$bam_raw")
bam_raw_index="${bam_raw%.bam}.bai"

if [[ $REF == *b37*.fasta ]]
 then
  chrM=MT
 else
  chrM=chrM
fi

$SAM view -b $bam_raw $chrM > $out_raw
$SAM index $out_raw

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
vcf_file="VCF_file.vcf"

# Check if the file exists
if [ ! -f "$vcf_file" ]; then
    echo "Error: File '$vcf_file' not found!"
    exit 1
fi

in_file=prioritized_variants.txt
out_file=append_$$.txt
final_file=mit_prioritized_variants.txt
parse_var=$BINDIR/parse_var.pl
echo -e "REF\tALT\tGT\tDP\tHF" > $out_file
for var in $(cut -f1 $in_file | sed '1d' | $parse_var) 
do
   grep -P "chrMT\t$var\t" $vcf_file | cut -f4,5,10 | tr ':' '\t' |cut -f1-5 >>  $out_file || true
done
paste $in_file $out_file > $final_file
rm $out_file

# HMTL creation
echo "Creating Browser HTML..."
mit_json=mit.json
mit_raw_json=mit.raw.json
$PYBINDIR/mtb2json.py  -i $final_file -f json > $mit_raw_json
$PYBINDIR/mtb2json.py  -i $final_file -f json4html > $BROWSERDIR/$mit_json
$PYBINDIR/mtb2html.py --id $id --json $mit_json --out $BROWSERDIR/$job_id.html --job-id $job_id
ln -s $ASSETS $BROWSERDIR/assets

cat<<EOF>$BROWSERDIR/README.txt
# To visualize <$job_id.html>:

# Option 1: Open <176099009134887.html> directly in Chromium
chromium --allow-file-access-from-files --disable-web-security $job_id.html

# Option 2: Use an HTTP server. Example using Python 3:
python3 -m http.server
EOF

# Fin
echo "All done!!!"
exit
