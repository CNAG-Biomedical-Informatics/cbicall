#!/usr/bin/env bash
#
#    Sex Determination script
#
#   Last Modified; March/05/2025
#
#   Copyright (C) 2025-2026 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source env.sh from the same directory
source "${CBICALL_ENV_FILE:-$SCRIPT_DIR/env.sh}"

VCF=$1

# Bail out early on empty VCF (no non-header lines)
if ! zgrep -q -v '^#' "$VCF"; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] No variant records found in $VCF"
  echo "SEX=UNKNOWN"
  exit 0
fi

# VCFs from exome should have only chr1..22,X,Y
# Note that we are not filtering PAR (pseudoautosomic regions)
# zgrep will work for both gz and ungz
mean_depth_for_scope() {
  local SCOPE="$1"
  zgrep -v '^#' "$VCF" | awk -v scope="$SCOPE" '
    BEGIN { FS = "\t" }
    function wanted(chrom) {
      if (scope == "autosomes") return chrom !~ /^(chr)?[XYM]$/
      if (scope == "X") return chrom ~ /^(chr)?X$/
      if (scope == "Y") return chrom ~ /^(chr)?Y$/
      return 0
    }
    wanted($1) {
      split($(NF - 1), fmt, ":")
      dp_idx = 0
      for (i = 1; i <= length(fmt); i++) {
        if (fmt[i] == "DP") {
          dp_idx = i
          break
        }
      }
      if (dp_idx == 0) next
      split($NF, sample, ":")
      if (sample[dp_idx] ~ /^[0-9]+([.][0-9]+)?$/) {
        sum += sample[dp_idx]
        n++
      }
    }
    END {
      if (n > 0) printf "%.2f\n", sum / n
      else print "NA"
    }
  '
}

MEAN_DEPTH_AUTOSOMES=$(mean_depth_for_scope autosomes)
MEAN_DEPTH_X=$(mean_depth_for_scope X)
MEAN_DEPTH_Y=$(mean_depth_for_scope Y)

if [ "$MEAN_DEPTH_AUTOSOMES" = "NA" ] || [ "$MEAN_DEPTH_X" = "NA" ] || [ "$MEAN_DEPTH_Y" = "NA" ]; then
  echo "MEAN DEPTH FOR AUTOSOMES=$MEAN_DEPTH_AUTOSOMES"
  echo "MEAN DEPTH FOR X=$MEAN_DEPTH_X"
  echo "MEAN DEPTH FOR Y=$MEAN_DEPTH_Y"
  echo "THRESHOLD=NA"
  echo "SEX=UNKNOWN"
  exit 0
fi


# Threshold will be given by mean autosomal coverage
# Usually, mean_depth_autosomes ~ 65 so threshold ~ 25-30 is enough
# If the mean_depth_autosomes ~ 50 then the same threshold will not work
# Note that this threshold may not work with other Exome captures
LOW_MEAN_DEPTH_AUTOSOMES=52
MEAN_DEPTH_AUTOSOMES_INT=$( echo "$MEAN_DEPTH_AUTOSOMES"  | awk '{print int($1)}' )
if [ "$MEAN_DEPTH_AUTOSOMES_INT" -le "$LOW_MEAN_DEPTH_AUTOSOMES" ]
 then
  THRESHOLD=$( echo "$MEAN_DEPTH_AUTOSOMES" | awk '{print int($1/3.5)}' ) # Taking integer part
else
  THRESHOLD=$( echo "$MEAN_DEPTH_AUTOSOMES" | awk '{print int($1/3.0)}' ) # Taking integer part
fi

# Determining sex by Female cov(X)>>cov(Y)
DIFF_DEPTH=$( echo "$MEAN_DEPTH_X" "$MEAN_DEPTH_Y" | awk '{print int($1-$2)}' ) # Taking integer part
if [ "$DIFF_DEPTH" -ge "$THRESHOLD" ]
then
  SEX="FEMALE"
else 
  SEX="MALE"
fi

# Printing results to stdout
echo "MEAN DEPTH FOR AUTOSOMES=$MEAN_DEPTH_AUTOSOMES"
echo "MEAN DEPTH FOR X=$MEAN_DEPTH_X"
echo "MEAN DEPTH FOR Y=$MEAN_DEPTH_Y"
echo "THRESHOLD=$THRESHOLD"
echo "SEX=$SEX"
