#!/usr/bin/env bash
#
# w[eg]s_cohort_genomicsdb_with_vqsr.sh
# Cohort joint-genotyping wrapper using GenomicsDBImport -> GenotypeGVCFs -> VQSR/Hard-filters (GATK 4.6)
# Last Modified: 2025-10-13
set -eu
set -o pipefail

function usage {
  cat <<EOF
Usage: $0 -m <sample_map.tsv> [-p wes|wgs] [-w <workspace>] [-t <threads>]

  -s  --sample-map   Sample map file for --sample-name-map (required)
  -p  --pipeline    'wes' (default) or 'wgs'
  -w  --workspace   GenomicsDB workspace path (default: ./cohort.genomicsdb.<job_id>)
  -h  --help        Show this help
EOF
  exit 1
}

# Parse args
PIPELINE="wes"
WORKSPACE=""
SAMPLE_MAP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--sample-map) SAMPLE_MAP="$2"; shift 2;;
    -p|--pipeline) PIPELINE="$2"; shift 2;;
    -w|--workspace) WORKSPACE="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -h|--help) usage ;;
    *) echo "Unknown arg: $1" >&2; usage ;;
  esac
done

if [ -z "$SAMPLE_MAP" ]; then
  echo "Error: sample_map is required." >&2
  usage
fi

# Load parameters
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$BINDIR/parameters.sh"

# Prepare output directories and logging
dir=$(pwd)
VARCALLDIR=$dir/02_varcall
mkdir -p "$VARCALLDIR"
cd $VARCALLDIR

LOGDIR="logs"
mkdir -p "$LOGDIR"
LOG="$LOGDIR/cohort_genomicsdb.log"

# pipeline mode
# Convert to Uppercase
PIPELINE=${PIPELINE^^}
if [[ "$PIPELINE" != "WES" && "$PIPELINE" != "WGS" ]]; then
  echo "Error: pipeline must be 'wes' or 'wgs'." >&2; exit 1
fi

# sample count & workspace default naming
if [ ! -s "$SAMPLE_MAP" ]; then
  echo "Error: sample_map '$SAMPLE_MAP' not found or empty." >&2; exit 1
fi
SAMPLE_COUNT=$(wc -l < "$SAMPLE_MAP" | tr -d ' ')
if [ -z "$WORKSPACE" ]; then
  WORKSPACE="${WORKSPACE}_${SAMPLE_COUNT}"
fi

# interval argument for WES vs WGS
# Set interval argument for WES vs WGS
if [ "$PIPELINE" = "WES" ]; then
  INTERVAL_ARG="-L $INTERVAL_LIST"
  echo "WES mode: restricting to $INTERVAL_LIST"
else
  INTERVAL_ARG=""
  echo "WGS mode: processing whole genome"
fi

# Derived output names
COHORT_RAW_VCF="cohort.raw.vcf.gz"
COHORT_VQSR_SNP="cohort.snp.recal.vcf.gz"
COHORT_SNP_TRANCHES="cohort.snp.tranches.txt"
COHORT_VQSR_INDEL="cohort.indel.recal.vcf.gz"
COHORT_INDEL_TRANCHES="cohort.indel.tranches.txt"
COHORT_POST_SNP="cohort.hc.post_snp.vcf.gz"
COHORT_POST_VQSR="cohort.gv.vqsr.vcf.gz"
COHORT_QC_VCF="cohort.gv.QC.vcf.gz"

echo "## Cohort GenomicsDBImport -> Genotype -> VQSR/Hard-filter"
echo "sample_map: $SAMPLE_MAP"
echo "pipeline: $PIPELINE"
echo "sample_count: $SAMPLE_COUNT"
echo "workspace: $WORKSPACE"
echo "out_vcf: $COHORT_RAW_VCF"
echo "tmpdir: $TMPDIR"
echo "log: $LOG"
echo "" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# Step 1: GenomicsDBImport
# -----------------------------------------------------------------------------
echo ">>> Step 1: GenomicsDBImport" | tee -a "$LOG"
mkdir -p "$(dirname "$WORKSPACE")"

set -x
"$GATK4_BIN" $GATK4_JAVA_OPTS_64G GenomicsDBImport \
  --sample-name-map "$SAMPLE_MAP" \
  --genomicsdb-workspace-path "$WORKSPACE" \
  --merge-input-intervals true \
  $INTERVAL_ARG \
  --tmp-dir "$TMPDIR" \
  2>> "$LOG"
set +x

# -----------------------------------------------------------------------------
# Step 2: GenotypeGVCFs
# -----------------------------------------------------------------------------
echo ">>> Step 2: GenotypeGVCFs" | tee -a "$LOG"
if [ -z "${REF:-}" ]; then
  echo "Error: REF not set (expected to be defined in parameters.sh)." >&2
  exit 1
fi

set -x
"$GATK4_BIN" $GATK4_JAVA_OPTS_64G GenotypeGVCFs \
  -R "$REF" \
  -V "gendb://$WORKSPACE" \
  -O "$COHORT_RAW_VCF" \
  --stand-call-conf 10 \
  --tmp-dir "$TMPDIR" \
  $INTERVAL_ARG \
  2>> "$LOG"
set +x

if [ $? -ne 0 ]; then
  echo "ERROR: GenotypeGVCFs failed. See log: $LOG" >&2
  exit 1
fi
echo "Genotyping completed. Raw cohort VCF: $COHORT_RAW_VCF" | tee -a "$LOG"

# -----------------------------------------------------------------------------
# Step 3: Count SNPs/INDELs and decide on VQSR
# -----------------------------------------------------------------------------
echo ">>> Step 3: Count variants and decide on VQSR" | tee -a "$LOG"
nSNP=$(zgrep -v '^#' "$COHORT_RAW_VCF" | awk 'length($5)==1' | wc -l)
nINDEL=$(zgrep -v '^#' "$COHORT_RAW_VCF" | awk 'length($5)!=1' | wc -l)
nSNP=$(echo "$nSNP" | tr -d ' ')
nINDEL=$(echo "$nINDEL" | tr -d ' ')
echo "Found SNPs: $nSNP ; INDELs: $nINDEL" | tee -a "$LOG"

apply_snp=false
apply_indel=false
minSNP=${MIN_SNP_FOR_VQSR:-1000}
minINDEL=${MIN_INDEL_FOR_VQSR:-8000}

# -----------------------------------------------------------------------------
# Step 4: VariantRecalibrator (SNP)
# -----------------------------------------------------------------------------
if (( nSNP >= minSNP )); then
  echo ">>> Building SNP VQSR model (VariantRecalibrator)" | tee -a "$LOG"
  set -x
  "$GATK4_BIN" $GATK4_JAVA_OPTS VariantRecalibrator \
    -R "$REF" \
    -V "$COHORT_RAW_VCF" \
    $SNP_RES \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
    --mode SNP \
    -O "$COHORT_VQSR_SNP" \
    --tranches-file "$COHORT_SNP_TRANCHES" \
    --max-gaussians 6 \
    --tmp-dir "$TMPDIR" \
    2>> "$LOG"
  set +x
  apply_snp=true
else
  echo "Skipping SNP VQSR (found $nSNP < min $minSNP)" | tee -a "$LOG"
fi

# -----------------------------------------------------------------------------
# Step 5: VariantRecalibrator (INDEL)
# -----------------------------------------------------------------------------
if (( nINDEL >= minINDEL )); then
  echo ">>> Building INDEL VQSR model (VariantRecalibrator)" | tee -a "$LOG"
  set -x
  "$GATK4_BIN" $GATK4_JAVA_OPTS VariantRecalibrator \
    -R "$REF" \
    -V "$COHORT_RAW_VCF" \
    $INDEL_RES \
    -an QD -an FS -an ReadPosRankSum \
    --mode INDEL \
    -O "$COHORT_VQSR_INDEL" \
    --tranches-file "$COHORT_INDEL_TRANCHES" \
    --max-gaussians 4 \
    --tmp-dir "$TMPDIR" \
    2>> "$LOG"
  set +x
  apply_indel=true
else
  echo "Skipping INDEL VQSR (found $nINDEL < min $minINDEL)" | tee -a "$LOG"
fi

# -----------------------------------------------------------------------------
# Step 6: Apply VQSR (if models exist)
# -----------------------------------------------------------------------------
echo ">>> Step 6: Apply VQSR (if available)" | tee -a "$LOG"
tmp_vcf="$COHORT_RAW_VCF"
if [ "$apply_snp" = true ]; then
  echo "Applying SNP recalibration..." | tee -a "$LOG"
  set -x
  "$GATK4_BIN" $GATK4_JAVA_OPTS ApplyVQSR \
    -R "$REF" -V "$tmp_vcf" \
    --recal-file "$COHORT_VQSR_SNP" \
    --tranches-file "$COHORT_SNP_TRANCHES" \
    --mode SNP --truth-sensitivity-filter-level 99.0 \
    -O "$COHORT_POST_SNP" \
    --tmp-dir "$TMPDIR" \
    2>> "$LOG"
  set +x
  tmp_vcf="$COHORT_POST_SNP"
fi

if [ "$apply_indel" = true ]; then
  echo "Applying INDEL recalibration..." | tee -a "$LOG"
  set -x
  "$GATK4_BIN" $GATK4_JAVA_OPTS ApplyVQSR \
    -R "$REF" -V "$tmp_vcf" \
    --recal-file "$COHORT_VQSR_INDEL" \
    --tranches-file "$COHORT_INDEL_TRANCHES" \
    --mode INDEL --truth-sensitivity-filter-level 95.0 \
    -O "$COHORT_POST_VQSR" \
    --tmp-dir "$TMPDIR" \
    2>> "$LOG"
  set +x
  tmp_vcf="$COHORT_POST_VQSR"
fi

# -----------------------------------------------------------------------------
# Step 7: Hard filters & write cohort QC VCF (always run for QC)
# (filters copied from your single-sample script)
# -----------------------------------------------------------------------------
echo ">>> Step 7: Hard-filter & write cohort QC VCF" | tee -a "$LOG"
set -x
"$GATK4_BIN" $GATK4_JAVA_OPTS VariantFiltration \
  -R "$REF" \
  -V "$tmp_vcf" \
  --filter-name "LowQUAL" --filter-expression "QUAL < 30.0" \
  --filter-name "QD2"        --filter-expression "QD < 2.0" \
  --filter-name "FS60"       --filter-expression "FS > 60.0" \
  --filter-name "MQ40"       --filter-expression "MQ < 40.0" \
  --filter-name "MQRS-12.5"  --filter-expression "MQRankSum < -12.5" \
  --filter-name "RPRS-8"     --filter-expression "ReadPosRankSum < -8.0" \
  --filter-name "QD2_indel"  --filter-expression "QD < 2.0" \
  --filter-name "FS200"      --filter-expression "FS > 200.0" \
  --filter-name "RPRS-20"    --filter-expression "ReadPosRankSum < -20.0" \
  -O "$COHORT_QC_VCF" \
  --tmp-dir "$TMPDIR" \
  2>> "$LOG"
set +x

if [ $? -ne 0 ]; then
  echo "ERROR: VariantFiltration (cohort QC) failed. See log: $LOG" >&2
  exit 1
fi

echo "Cohort QC VCF written: $COHORT_QC_VCF" | tee -a "$LOG"

# Final message
echo "All done. Cohort raw: $COHORT_RAW_VCF ; Cohort QC: $COHORT_QC_VCF" | tee -a "$LOG"
exit 0
