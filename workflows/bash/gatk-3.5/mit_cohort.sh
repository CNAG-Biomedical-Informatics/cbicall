#!/usr/bin/env bash
# 
#   mtDNA Pipeline Cohort Bash script.
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
    NB2: There is no need to run wes_cohort prior to mit_cohort.

MA00024_exome  <-- ID taken from here
├── MA0002401P_ex
│   └── cbicall_wes_single_146723488708442
│       ├── 01_bam
│       ├── 02_varcall
│       └── 03_stats
├── MA0002402M_ex
│   └── cbicall_wes_single_146727114980481
│       ├── 01_bam
│       ├── 02_varcall
│       └── 03_stats
├── MA0002402P_ex
│   └── cbicall_wes_single_146730170886696
│       ├── 01_bam
│       ├── 02_varcall
│       └── 03_stats
└── cbicall_bash_mit_cohort_146774466308431 <- The script expects that you are submitting the job from inside this directory
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
if [[ $DIR != *cbicall_bash_mit_cohort* ]]
 then 
  usage
fi 

# The id need to have this format LP6005831-???_???.bam, otherwise MToolBox will fail
cohort=$( echo $DIR | awk -F'/' '{print $--NF}' | awk -F'_' '{print $1}' | sed 's/$/-DNA_MIT/')
echo $cohort

# From now on we will work on VARCALL dir
VARCALLDIR=$DIR/MTOOLBOX
mkdir $VARCALLDIR
cd $VARCALLDIR

# NB: We are using UCSC's hg19 for Exome.
# There are a few minor differences between GRCh37 and hg19. 
# The contig sequences are the same but the names are different, i.e. "1" may need to be converted to "chr1". 
# In addition UCSC hg19 is currenly using the old mitochondrial sequence but NCBI and Ensembl have transitioned to NC_012920.
# For using MttolBox we need to align again to RSRS

# Using Samtools to extract chrM
# NB: BAMs may include duplicated entries at this stage
echo "Extracting Mitochondrial DNA from exome BAM file..."
for BAMDIR in ../../??????????_ex/cbicall_bash_wes_single*/01_bam
do
 id=$( echo $BAMDIR | awk -F'/' '{print $3}' | awk -F'_' '{print $1}' | sed 's/$/-DNA_MIT/')
 bam_raw=$BAMDIR/input.merged.filtered.realigned.fixed.bam
 out_raw=$id.bam

 # The index name must be foo.bam.bai instead of foo.bai (can happen if wes_single.sh failed)
 bam_raw_index=$BAMDIR/input.merged.filtered.realigned.fixed.bai
 bam_raw_index_ok=$BAMDIR/input.merged.filtered.realigned.fixed.bam.bai
 if [[ ! -s $bam_raw_index_ok ]]
  then
  cp $bam_raw_index $bam_raw_index_ok
 fi
 
 if [[ $REF == *b37*.fasta ]]
  then
   chrM=MT
 else
   chrM=chrM
 fi

 $SAM view -b $bam_raw $chrM > $out_raw
 $SAM index $out_raw
 
done

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
vcf_file=VCF_file.vcf
vcf_tmp=VCF_file_$$.vcf
in_file=prioritized_variants.txt
out_file=append_$$.txt
final_file=mit_prioritized_variants.txt
parse_var=$BINDIR/parse_var.pl
parse_prior=$BINDIR/parse_prioritized.pl
grep ^#CHROM $vcf_file > $vcf_tmp
for var in $(cut -f1 $in_file | sed '1d' | $parse_var) 
do
  grep -P "chrMT\t$var\t" $vcf_file >>  $vcf_tmp  || echo "$var not found"
done
$parse_prior -i $vcf_tmp > $out_file
paste $in_file $out_file > $final_file
rm $vcf_tmp $out_file

# Fin
echo "All done!!!"
exit 
