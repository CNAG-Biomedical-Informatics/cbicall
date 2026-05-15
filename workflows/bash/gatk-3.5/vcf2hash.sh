#!/usr/bin/env bash
#
#   VCF SHA-256 fingerprint helper.
#
#   Reports both the raw file hash and a normalized hash that ignores VCF
#   headers and sorts variant records with LC_ALL=C.

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

cat <<EOF
FILE=$vcf
ALGORITHM=sha256
RAW_SHA256=$raw_sha256
NORMALIZED_PATTERN=$pattern
NORMALIZED_SORT=LC_ALL=C
NORMALIZED_RECORDS=$normalized_records
NORMALIZED_SHA256=$normalized_sha256
EOF
