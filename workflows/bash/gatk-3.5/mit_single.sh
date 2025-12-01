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

# Set up variables and Defining directories
DIR=$( pwd )
BINDIRMTB=$BINDIR/../../../mtdna
BROWSERDIR=$BINDIR/../../../browser
ASSETS=$BROWSERDIR/assets

# The id needs to have this format LP6005831-???_???.bam, otherwise MToolBox will fail
id=$( echo "$DIR" | awk -F'/' '{print $(NF-1)}' | awk -F'_' '{print $1}' | sed 's/$/-DNA_MIT/' )
job_id=$( echo "$DIR" | awk -F'_' '{print $NF}' )

# From now on we will work on VARCALL dir
VARCALLDIR=$DIR/01_mtoolbox
mkdir "$VARCALLDIR"
cd "$VARCALLDIR"

# Using Samtools to extract chrM
echo "Extracting Mitochondrial DNA from exome BAM file..."

out_raw=$id.bam

# Prefer GATK 3.5 BAM if available
bam_raw=""

for f in ../../cbicall_bash_wes_single_gatk-3.5*/01_bam/input.merged.filtered.realigned.fixed.bam
do
    echo $f
    if [ -f "$f" ]; then
        bam_raw="$f"
        echo "Using GATK 3.5 BAM: $bam_raw"
        break
    fi
done

# If no GATK 3.5 BAM, fall back to GATK 4.6 naming: ${ID}.rg.merged.dedup.recal.bam
if [ -z "$bam_raw" ]; then
    # ID is expected to be defined (e.g. via parameters.sh)
    if [ -z "${ID:-}" ]; then
        echo "ERROR: ID is not set and no GATK 3.5 BAM was found." >&2
        exit 1
    fi

    for f in ../../cbicall_bash_wes_single_gatk-4.6*/01_bam/"$ID".rg.merged.dedup.recal.bam
    do
        if [ -f "$f" ]; then
            bam_raw="$f"
            echo "Using GATK 4.6 BAM: $bam_raw"
            break
        fi
    done
fi

# If still nothing found, bail out
if [ -z "$bam_raw" ]; then
    echo "ERROR: Could not find BAM for ID '${ID:-$id}' in either:" >&2
    echo "  ../../cbicall_bash_wes_single_gatk-3.5*/01_bam/input.merged.filtered.realigned.fixed.bam" >&2
    echo "  ../../cbicall_bash_wes_single_gatk-4.6*/01_bam/\$ID.rg.merged.dedup.recal.bam" >&2
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
export PATH="$MTOOLBOXDIR:$PATH"

# Add the local site-packages to PYTHONPATH
export PATH="$PY27_PREFIX:$PATH"
export PYTHONPATH="$PY27_PREFIX/Lib/site-packages:${PYTHONPATH:-}"

cp $BINDIRMTB/MToolBox_config.sh .
MToolBox.sh -i MToolBox_config.sh -m "-t $THREADS"

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
HTMLDIR=../02_browser
mkdir $HTMLDIR
mit_json=mit.json
$BROWSERDIR/mtb2json.py  -i $final_file -f json4html > $HTMLDIR/$mit_json
$BROWSERDIR/mtb2html.py --id $id --json $mit_json --out $HTMLDIR/$job_id.html --job-id $job_id
ln -s $ASSETS $HTMLDIR/assets

cat<<EOF>$HTMLDIR/README.txt
# To visualize <$job_id.html>:

# Option 1: Open <176099009134887.html> directly in Chromium
chromium --allow-file-access-from-files --disable-web-security $job_id.html

# Option 2: Use an HTTP server. Example using Python 3:
python3 -m http.server
EOF


# Fin
echo "All done!!!"
exit 
