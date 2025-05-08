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

    NB1: The script is expecting that you follow STSI nomenclature for samples

MA00047_exome
└── MA0004701P_ex  <--- ID taken from here
    ├── MA0004701P_ex_S5_L001_R1_001.fastq.gz
    ├── MA0004701P_ex_S5_L001_R2_001.fastq.gz
    ├── MA0004701P_ex_S5_L002_R1_001.fastq.gz
    ├── MA0004701P_ex_S5_L002_R2_001.fastq.gz
    └── cbicall_bash_mit_single_146657420113136 <- The script expects that you are submitting the job from inside this directory
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
BINDIRMTB=$BINDIR/../../mtdna

# Check that nomenclature exists
if [[ $DIR != *cbicall_bash_mit_single* ]]
 then 
  usage
fi 

# The id need to have this format LP6005831-???_???.bam, otherwise MToolBox will fail
id=$( echo $DIR | awk -F'/' '{print $(NF-1)}' | awk -F'_' '{print $1}' |sed 's/$/-DNA_MIT/' )

# From now on we will work on VARCALL dir
VARCALLDIR=$DIR/MTOOLBOX
mkdir $VARCALLDIR
cd $VARCALLDIR

# NB: We are using UCSC's hg19 for Exome.
# There are a few minor differences between GRCh37 and hg19. 
# The contig sequences are the same but the names are different, i.e. "1" may need to be converted to "chr1". 
# In addition UCSC hg19 is currenly using the old mitochondrial sequence but NCBI and Ensembl have transitioned to NC_012920.
# For using MtoolBox we need to align again to RSRS

# Using Samtools to extract chrM
# NB: BAMs may include duplicates entries at this stage
echo "Extracting Mitochondrial DNA from exome BAM file..."
BAMDIR=../../cbicall_bash_wes_single*/01_bam
bam_raw=$BAMDIR/input.merged.filtered.realigned.fixed.bam
bam_raw_index=$BAMDIR/input.merged.filtered.realigned.fixed.bai
out_raw=$id.bam

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
export PYTHONPATH=~/.local/lib/python2.7/site-packages:${PYTHONPATH:-}

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

# Fin
echo "All done!!!"
exit 
