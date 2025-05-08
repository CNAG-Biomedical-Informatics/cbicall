#!/usr/bin/env bash
#
#   Coverage stats for chromosome 1 using exome (WES) or whole-genome (WGS)
#   Uses $SAM and $REF from parameters.sh, avoids division-by-zero,
#   prints total_reads as an additional column.
#   Handles conversion of Picard IntervalList (1-based inclusive) to 0-based BED.
#
#   Last Modified: May/02/2025
#   version: $VERSION from CBICall

set -eu

# Usage: $0 <sampleID> <raw.bam> <dedup.bam> [WES|WGS]
if [ $# -lt 3 ] || [ $# -gt 4 ]; then
  echo "Usage: $0 <sampleID> <raw.bam> <dedup.bam> [WES|WGS]" >&2
  exit 1
fi

# Load parameters (defines INTERVAL_LIST, REF, SAM, etc.)
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$BINDIR/parameters.sh"

# Reference FASTA from parameters.sh
FASTA="${REF}"
if [ ! -f "${FASTA}.fai" ]; then
  echo "Error: reference FASTA or index not found (expected ${FASTA}.fai)" >&2
  exit 1
fi

sid=$1
RAWBAM=$2
DEDUPBAM=$3
mode="${4:-WES}"
mode_upper="${mode^^}"

# We're focusing on chromosome 1 (b37 contig '1')
chrN="1"

# Create a temp file to hold a 0-based BED region for chr1
TMPREG=$(mktemp)
case "$mode_upper" in
  WES)
    # WES: filter Picard IntervalList (1-based inclusive) for chr1,
    # then convert to 0-based, half-open BED by subtracting 1 from the start.
    grep -v '^@' "$INTERVAL_LIST" \
      | awk -v chr="$chrN" '$1==chr { printf "%s\t%d\t%d\n", $1, $2-1, $3 }' \
      > "$TMPREG"
    ;;
  WGS)
    # WGS: get full chromosome length and create a single BED interval
    # BED is 0-based, half-open: start at 0, end at chrlen
    chrlen=$(awk -v chr="$chrN" '$1==chr {print $2; exit}' "${FASTA}.fai")
    if [ -z "$chrlen" ]; then
      echo "Error: chromosome $chrN not found in ${FASTA}.fai" >&2
      rm -f "$TMPREG"
      exit 1
    fi
    printf "%s\t0\t%s\n" "$chrN" "$chrlen" > "$TMPREG"
    ;;
  *)
    echo "Unknown mode '$mode_upper'. Must be WES or WGS." >&2
    rm -f "$TMPREG"
    exit 1
    ;;
esac

REGION="$TMPREG"

# Compute total span (bp) of the BED regions
span=$(awk '{s+=($3-$2)} END{print s}' "$REGION")

# Sum of depths (including duplicates) from raw BAM over the region
dup_sum=$("${SAM}" depth -b "$REGION" "$RAWBAM" 2>/dev/null | awk '{s+=$3} END{print s}')

# Count total reads in BAM and bail if empty
total_reads=$("${SAM}" view -c "$RAWBAM")
if [ "$total_reads" -eq 0 ]; then
  echo "Error: no reads found in BAM $RAWBAM" >&2
  rm -f "$TMPREG"
  exit 1
fi

# Count reads overlapping our BED region
region_reads=$("${SAM}" view -c -L "$REGION" "$RAWBAM")
# Reads outside the region
out_region=$(( total_reads - region_reads ))

# Mean insert size (only pairs with 0 < insert < 600)
ins_size=$("${SAM}" view "$RAWBAM" | awk '$9>0 && $9<600 {sum+=$9; cnt++} END {if(cnt) printf "%.1f", sum/cnt; else print "NA"}')

# Depth-by-position on dedup BAM and summary metrics
"${SAM}" depth -b "$REGION" "$DEDUPBAM" | awk \
  -v sid="$sid" -v mode="$mode_upper" -v chr="$chrN" \
  -v span="$span" -v dup_sum="$dup_sum" -v ins="$ins_size" \
  -v r_in="$region_reads" -v tot="$total_reads" -v r_out="$out_region" '{ sum += $3 }
  $3 >= 10 { c10++ }
  END {
    mean   = (span>0    ? sum/span       : 0)
    p10    = (span>0    ? 100 * c10/span : 0)
    nondup = (dup_sum>0 ? 100 * sum/dup_sum : 0)
    in_pct = (tot>0     ? 100 * r_in/tot   : 0)
    out_pct= (tot>0     ? 100 * r_out/tot  : 0)

    # print header with added total_reads column
    printf "%s\n", chr
    printf "sampleID\tmode\ttotal_reads\tmean_cov\tten_pct\tnondup_pct\tins_size\tin_pct\tout_pct\n"
    # print data line including tot
    printf "%s\t%s\t%d\t%.1f\t%.1f\t%.1f\t%s\t%.1f\t%.1f\n", sid, mode, tot, mean, p10, nondup, ins, in_pct, out_pct
  }'

# Cleanup temporary BED file
rm -f "$TMPREG"
exit 0
