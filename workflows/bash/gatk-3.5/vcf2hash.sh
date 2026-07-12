#!/usr/bin/env bash
#
#   VCF SHA-256 fingerprint helper.
#
#   Reports the raw file hash, a strict normalized record hash that ignores VCF
#   headers, and a call-level hash over CHROM, POS, REF, ALT, FILTER, and GT.

set -euo pipefail
export LC_ALL=C

function usage {
  echo "Usage: $0 file.vcf[.gz] [header_pattern]" >&2
  echo "Default header_pattern: ^#" >&2
  exit 1
}

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
  usage
fi

vcf=$1
pattern=${2:-${VCF2HASH_PATTERN:-${VCF_HASH_PATTERN:-'^#'}}}

if [ ! -s "$vcf" ]; then
  echo "Error: VCF not found or empty: $vcf" >&2
  exit 1
fi

function stream_vcf {
  if [[ "$vcf" == *.gz ]]; then
    gzip -cd -- "$vcf"
  else
    cat -- "$vcf"
  fi
}

raw_sha256=$(sha256sum "$vcf" | awk '{print $1}')
normalized_records=$(stream_vcf | awk -v pattern="$pattern" '$0 !~ pattern' | wc -l | tr -d ' ')
normalized_sha256=$(stream_vcf | awk -v pattern="$pattern" '$0 !~ pattern' | sort | sha256sum | awk '{print $1}')
call_records=$(stream_vcf | awk -v pattern="$pattern" '
  $0 !~ pattern && NF >= 8 {
    gt_values = "."
    if (NF >= 10) {
      gt_index = 0
      n = split($9, format_fields, ":")
      for (i = 1; i <= n; i++) {
        if (format_fields[i] == "GT") {
          gt_index = i
          break
        }
      }
      gt_values = ""
      for (sample_col = 10; sample_col <= NF; sample_col++) {
        split($sample_col, sample_fields, ":")
        gt = (gt_index > 0 && gt_index in sample_fields) ? sample_fields[gt_index] : "."
        gt_values = gt_values (sample_col == 10 ? "" : ",") gt
      }
    }
    print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $7 "\t" gt_values
  }
' | wc -l | tr -d ' ')
call_sha256=$(stream_vcf | awk -v pattern="$pattern" '
  $0 !~ pattern && NF >= 8 {
    gt_values = "."
    if (NF >= 10) {
      gt_index = 0
      n = split($9, format_fields, ":")
      for (i = 1; i <= n; i++) {
        if (format_fields[i] == "GT") {
          gt_index = i
          break
        }
      }
      gt_values = ""
      for (sample_col = 10; sample_col <= NF; sample_col++) {
        split($sample_col, sample_fields, ":")
        gt = (gt_index > 0 && gt_index in sample_fields) ? sample_fields[gt_index] : "."
        gt_values = gt_values (sample_col == 10 ? "" : ",") gt
      }
    }
    print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $7 "\t" gt_values
  }
' | sort | sha256sum | awk '{print $1}')
sample_count=$(stream_vcf | awk -F '\t' '
  $1 == "#CHROM" && !found {
    print (NF > 9 ? NF - 9 : 0)
    found = 1
  }
  END { if (!found) print 0 }
')
sample_order_sha256=$(stream_vcf | awk -F '\t' '
  $1 == "#CHROM" && !found {
    for (i = 10; i <= NF; i++) {
      printf "%s%s", (i == 10 ? "" : "\t"), $i
    }
    if (NF >= 10) printf "\n"
    found = 1
  }
' | sha256sum | awk '{print $1}')

cat <<EOF
FILE=$vcf
ALGORITHM=sha256
RAW_SHA256=$raw_sha256
NORMALIZED_PATTERN=$pattern
NORMALIZED_SORT=LC_ALL=C
NORMALIZED_RECORDS=$normalized_records
NORMALIZED_SHA256=$normalized_sha256
CALL_FIELDS=CHROM,POS,REF,ALT,FILTER,GT_ALL_SAMPLES
CALL_SORT=LC_ALL=C
CALL_RECORDS=$call_records
CALL_SHA256=$call_sha256
SAMPLE_ORDER_FIELDS=#CHROM_COLUMNS_10_PLUS
SAMPLE_COUNT=$sample_count
SAMPLE_ORDER_SHA256=$sample_order_sha256
EOF
