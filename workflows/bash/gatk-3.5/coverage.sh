#!/usr/bin/env bash
#
#   Coverage stats for chr1
#
#   Last Modified: May/02/2025
#
#   $VERSION taken from CBICall

set -eu

# Check arguments
if [ $# -ne 3 ]; then
  echo "Usage: $0 <sampleID> <chr1.raw.bam> <chr1.dedup.bam>" >&2
  exit 1
fi

# Determine the directory where the script resides
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source parameters.sh from the same directory
source "$BINDIR/parameters.sh"

chrN=chr1
REGION=$EXOM/hg19.$chrN.bed
EXOME_COORD=$EXOM/hg19_coor.$chrN.txt
sid=$1
RAWBAM=$2
DEDUPBAM=$3

# Calculate total exon span
exon_span=$( awk '{span+=($3-$2)} END{print span}' "$REGION" )
if [ "$exon_span" -eq 0 ]; then
  echo "Warning: exon span is zero for $chrN (no targets?)" >&2
fi

# Sum of depths in the raw BAM
dup_sum=$("$SAM" depth -b "$REGION" "$RAWBAM" 2>/dev/null | awk '{sum+=$3} END{print sum}')
if [ "$dup_sum" -eq 0 ]; then
  echo "Warning: duplicate‐sum is zero for $chrN (no raw coverage?)" >&2
fi

# Total reads and exome vs off‐exome
total_reads=$("$SAM" idxstats "$RAWBAM" | awk '{s+=$3+$4} END{print s}')
if [ "$total_reads" -eq 0 ]; then
  echo "Warning: total reads is zero in $RAWBAM" >&2
fi
exome_reads=$("$SAM" view "$RAWBAM" $(cat "$EXOME_COORD") | wc -l)
out_exome=$(( total_reads - exome_reads ))

# Mean insert size (only considering 0<ins<600)
ins_size=$("$SAM" view "$RAWBAM" | \
  awk '$9>0 && $9<600 { sum+=$9; cnt++ }
       END {
         if (cnt>0) {
           print sum/cnt
         } else {
           print "NA"
           print "Warning: no valid insert‐size pairs found" > "/dev/stderr"
         }
       }'
)

# Output stats for coverage, nonduplicates, and mean_insert_size
"$SAM" depth -b "$REGION" "$DEDUPBAM" | \
awk -v sid="$sid" \
    -v exon_span="$exon_span" \
    -v dup_sum="$dup_sum" \
    -v ins_size="$ins_size" \
    -v exome_reads="$exome_reads" \
    -v total_reads="$total_reads" \
    -v out_exome="$out_exome" \
    -v chrN="$chrN" '
{
  depth_sum += $3
}
$3 >= 10 {
  c10++
}
END {
  # safe divisions using ?: guard
  depth    = exon_span   > 0 ? depth_sum           / exon_span     : 0
  r10      = exon_span   > 0 ? 100 * c10           / exon_span     : 0
  nondup   = dup_sum     > 0 ? 100 * depth_sum     / dup_sum       : 0
  read_ex  = total_reads > 0 ? 100 * exome_reads   / total_reads   : 0
  read_out = total_reads > 0 ? 100 * out_exome      / total_reads   : 0

  # print table
  printf "%s\n", chrN
  printf "sampleID\tmean_coverage\tten_reads%%\tnonduplicate%%\tmean_insert_size\treads_in_exome%%\treads_out_of_exome%%\n"
  printf "%s\t%4.1f\t%4.1f\t%4.1f\t%s\t%4.1f\t%4.1f\n",
         sid, depth, r10, nondup, ins_size, read_ex, read_out
}
'

exit 0

