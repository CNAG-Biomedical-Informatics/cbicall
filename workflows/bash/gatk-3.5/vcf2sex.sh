#!/usr/bin/env bash
#
#    Sex Determination script
#
#   Last Modified; March/05/2025
#
#   $VERSION taken from CBICall
#
#   Copyright (C) 2025 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

#########################################################################
#
#   Quick solution for determing sex of samples.
#
#   In order to determine the sex from VCFs I tested different choices:
#   
#   1 - Number of Variants in chrY (and chrX). Now, it turns out that exome files from Male/Female samples contain variants
#       in chrY (~150-200). At first I thought they were within PARs (pseudoautosomal; see below), but they were not. We assume that they are simply mapping errors. 
#   http://blog.kokocinski.net/index.php/par-regions?blog=2
#   --------+----------+----------+---------+-----------+-----------+
#   | Y       |    10001 |  2649520 | X       |     60001 |   2699520 |
#   | Y       | 59034050 | 59373566 | X       | 154931044 | 155270560 |
#
#   2 - vcf2sex from samtools: Black Box. Not a very popular choice /pro/NGSutils/samtools-0.1.19/bcftools/bcftools +vcf2sex -h
#   3 - GATK: In our lab, we run GATK DepthOfCoverage with 3 beds (autosomes, chrX, chrY) to get 3 mean coverages. 
#      Females should have cov(X)>>cov(Y)  <=== USED
#    
#   NB: We are processing the WES vcf as it comes (PASS and non-PASS vars ) but keep in mind that other data can be messy (e.g., low pass WGS).

set -eu

#export TMPDIR=/media/mrueda/4TB/tmp
export LC_ALL=C

# Check arguments
if [ $# -ne 1 ]
 then
  echo "$0 file.vcf"
  exit 1
fi

# Determine the directory where the script resides
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source parameters.sh from the same directory
source "$BINDIR/parameters.sh"

vcf=$1
stdesv=$BINDIR/desv.pl

# Bail out early on empty VCF (no non-header lines)
if ! zgrep -q -v '^#' "$vcf"; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] No variant records found in $vcf"
  echo "SEX=UNKNOWN"
  exit 0
fi

# VCFs from exome should have only chr1..22,X,Y
# Note that we are not filtering PAR (pseudoautosomic regions)
# zgrep will work for both gz and ungz
mean_depth_autosomes=$( zgrep -E -v -w '^chr[XYM]|^[XYM]' "$vcf" | grep -v '^#' | awk '{print $NF}' | cut -d':' -f3 | $stdesv | awk '{print $1}' )
mean_depth_X=$(         zgrep -E -w  '^chrX|^X'           "$vcf" | grep -v '^#' | awk '{print $NF}' | cut -d':' -f3 | $stdesv | awk '{print $1}' )
mean_depth_Y=$(         zgrep -E -w  '^chrY|^Y'           "$vcf" | grep -v '^#' | awk '{print $NF}' | cut -d':' -f3 | $stdesv | awk '{print $1}' )


# Threshold will be given by mean autosomal coverage
# Usually, mean_depth_autosomes ~ 65 so threshold ~ 25-30 is enough
# If the mean_depth_autosomes ~ 50 then the same threshold will not work
# Note that this threshold may not work with other Exome captures
l_mda=52
i_mda=$( echo "$mean_depth_autosomes"  | awk '{print int($1)}' )
if [ "$i_mda" -le "$l_mda" ]
 then
  threshold=$( echo "$mean_depth_autosomes" | awk '{print int($1/3.5)}' ) # Taking integer part
else
  threshold=$( echo "$mean_depth_autosomes" | awk '{print int($1/3.0)}' ) # Taking integer part
fi

# Determining sex by Female cov(X)>>cov(Y)
diff_depth=$( echo "$mean_depth_X" "$mean_depth_Y" | awk '{print int($1-$2)}' ) # Taking integer part
if [ "$diff_depth" -ge "$threshold" ]
then
  sex="FEMALE"
else 
  sex="MALE"
fi

# Printing results to stdout
echo "MEAN DEPTH FOR AUTOSOMES=$mean_depth_autosomes"
echo "MEAN DEPTH FOR X=$mean_depth_X"
echo "MEAN DEPTH FOR Y=$mean_depth_Y"
echo "THRESHOLD=$threshold"
echo "SEX=$sex"
