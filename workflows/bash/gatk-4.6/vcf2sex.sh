#!/usr/bin/env bash
#
#   Sex Determination script
#
#   Last Modified; March/05/2025
#
#   Copyright (C) 2025-2026 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

#########################################################################
#
#   Lightweight sex inference QC for native CBIcall workflows.
#
#   Several VCF-based options were considered:
#
#   1 - Count variants on chrY and chrX. This is not robust: exome and WGS VCFs from female samples
#       can still contain chrY records, usually from mapping noise rather than PARs.
#       PAR reference:
#       http://blog.kokocinski.net/index.php/par-regions?blog=2
#       --------+----------+----------+---------+-----------+-----------+
#       | Y       |    10001 |  2649520 | X       |     60001 |   2699520 |
#       | Y       | 59034050 | 59373566 | X       | 154931044 | 155270560 |
#
#   2 - samtools-0.1.19 bcftools +vcf2sex. This is legacy and behaves as a black box.
#
#   3 - bcftools +guess-ploidy. This is a good independent VCF-based cross-check when a modern
#       bcftools build and BCFTOOLS_PLUGINS are available, for example:
#       BCFTOOLS_PLUGINS=/path/to/bcftools/plugins bcftools +guess-ploidy -g hg38 -t GT sample.vcf.gz
#       CBIcall does not use it here because modern bcftools is not part of the official resource bundle.
#
#   4 - GATK DepthOfCoverage on autosomes, chrX, and chrY. This is stronger when BAM-derived
#       chromosome-level coverage is required, but it is heavier and not used in this lightweight helper.
#
#   5 - USED BY CBICALL: VCF FORMAT/DP proxy. The helper compares mean DP on autosomal, chrX,
#       and chrY variant records. It first applies an X/autosome-ratio guard to avoid false male calls
#       from noisy chrY records in female WGS samples; otherwise it falls back to the X-Y depth difference.
#
#   NB: This is a QC signal, not definitive biological or clinical sex determination.
#   NB: We process the VCF as produced by the workflow, including PASS and non-PASS records.

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

print_method_comment() {
  echo "# METHOD: sex inferred for QC from VCF-derived depth proxies. If X/autosome ratio >= 0.75, classify FEMALE; otherwise classify FEMALE only when X-Y depth difference >= THRESHOLD."
}

# Bail out early on empty VCF (no non-header lines)
if ! zgrep -q -v '^#' "$VCF"; then
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] No variant records found in $VCF"
  print_method_comment
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

X_AUTOSOME_RATIO="NA"
if [ "$MEAN_DEPTH_AUTOSOMES" != "NA" ] && [ "$MEAN_DEPTH_X" != "NA" ]; then
  X_AUTOSOME_RATIO=$( echo "$MEAN_DEPTH_X" "$MEAN_DEPTH_AUTOSOMES" | awk '{if ($2 > 0) printf "%.2f", $1/$2; else print "NA"}' )
fi

if [ "$MEAN_DEPTH_AUTOSOMES" != "NA" ] && [ "$MEAN_DEPTH_X" != "NA" ] && [ "$MEAN_DEPTH_Y" = "NA" ]; then
  if awk -v ratio="$X_AUTOSOME_RATIO" 'BEGIN { exit !(ratio >= 0.55) }'; then
    print_method_comment
    echo "MEAN DEPTH FOR AUTOSOMES=$MEAN_DEPTH_AUTOSOMES"
    echo "MEAN DEPTH FOR X=$MEAN_DEPTH_X"
    echo "MEAN DEPTH FOR Y=$MEAN_DEPTH_Y"
    echo "X_AUTOSOME_RATIO=$X_AUTOSOME_RATIO"
    echo "X_MINUS_Y_DEPTH=NA"
    echo "THRESHOLD=NA"
    echo "DECISION=X/autosome ratio >= 0.55 and no usable Y records"
    echo "SEX=FEMALE_LIKELY"
    exit 0
  fi
fi

if [ "$MEAN_DEPTH_AUTOSOMES" = "NA" ] || [ "$MEAN_DEPTH_X" = "NA" ] || [ "$MEAN_DEPTH_Y" = "NA" ]; then
  print_method_comment
  echo "MEAN DEPTH FOR AUTOSOMES=$MEAN_DEPTH_AUTOSOMES"
  echo "MEAN DEPTH FOR X=$MEAN_DEPTH_X"
  echo "MEAN DEPTH FOR Y=$MEAN_DEPTH_Y"
  echo "X_AUTOSOME_RATIO=$X_AUTOSOME_RATIO"
  echo "X_MINUS_Y_DEPTH=NA"
  echo "THRESHOLD=NA"
  echo "DECISION=Insufficient autosome, X, or Y depth values"
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

# Determining sex by female-like X/autosome ratio or Female cov(X)>>cov(Y)
DIFF_DEPTH=$( echo "$MEAN_DEPTH_X" "$MEAN_DEPTH_Y" | awk '{print int($1-$2)}' ) # Taking integer part
if [ "$X_AUTOSOME_RATIO" != "NA" ] && awk -v ratio="$X_AUTOSOME_RATIO" 'BEGIN { exit !(ratio >= 0.75) }'; then
  SEX="FEMALE"
  DECISION="X/autosome ratio >= 0.75"
elif [ "$DIFF_DEPTH" -ge "$THRESHOLD" ]; then
  SEX="FEMALE"
  DECISION="X-Y depth difference >= threshold"
else
  SEX="MALE"
  DECISION="X/autosome ratio < 0.75 and X-Y depth difference < threshold"
fi

# Printing results to stdout
print_method_comment
echo "MEAN DEPTH FOR AUTOSOMES=$MEAN_DEPTH_AUTOSOMES"
echo "MEAN DEPTH FOR X=$MEAN_DEPTH_X"
echo "MEAN DEPTH FOR Y=$MEAN_DEPTH_Y"
echo "X_AUTOSOME_RATIO=$X_AUTOSOME_RATIO"
echo "X_MINUS_Y_DEPTH=$DIFF_DEPTH"
echo "THRESHOLD=$THRESHOLD"
echo "DECISION=$DECISION"
echo "SEX=$SEX"
