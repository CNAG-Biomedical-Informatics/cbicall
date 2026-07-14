#!/usr/bin/env bash
#
# w[eg]s_cohort_genomicsdb_with_vqsr.sh
# Cohort joint-genotyping wrapper using GenomicsDBImport -> GenotypeGVCFs -> VQSR/Hard-filters (GATK 4.6)
# Last Modified: 2025-10-13
# Registry version: v1
set -eu
set -o pipefail

function usage {
  cat <<EOF
Usage: $0 -m <sample_map.tsv> [-p wes|wgs] [-w <workspace>] [-t <threads>]

  -s  --sample-map   Sample map file for --sample-name-map (required)
  -p  --pipeline    'wes' (default) or 'wgs'
  -w  --workspace   GenomicsDB workspace name/path (relative names are stored under 01_genomicsdb)
      --cohort-stage all|shard|finalize
                     all: import, genotype, and filter (default)
                     shard: import/genotype one shard and stop after raw VCF
                     finalize: start from a gathered raw VCF and filter globally
      --output-basename NAME
                     Basename for cohort VCFs [cohort]
      --interval-shard CONTIG
                     Contig to extract from the canonical interval source for shard runs
      --input-vcf VCF
                     Gathered raw VCF used with --cohort-stage finalize
  -h  --help        Show this help
EOF
  exit 1
}

# Parse args
PIPELINE="wes"
WORKSPACE=""
SAMPLE_MAP=""
COHORT_STAGE="all"
OUTPUT_BASENAME="cohort"
INTERVAL_SHARD=""
INPUT_VCF=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--sample-map) SAMPLE_MAP="$2"; shift 2;;
    -p|--pipeline) PIPELINE="$2"; shift 2;;
    -w|--workspace) WORKSPACE="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    --cohort-stage) COHORT_STAGE="$2"; shift 2;;
    --output-basename) OUTPUT_BASENAME="$2"; shift 2;;
    --interval-shard) INTERVAL_SHARD="$2"; shift 2;;
    --input-vcf) INPUT_VCF="$2"; shift 2;;
    -h|--help) usage ;;
    *) echo "Unknown arg: $1" >&2; usage ;;
  esac
done

case "$COHORT_STAGE" in
  all|shard|finalize) ;;
  *) echo "Error: cohort_stage must be 'all', 'shard', or 'finalize'." >&2; exit 1;;
esac

case "$OUTPUT_BASENAME" in
  ""|.|..) echo "Error: output basename must be a non-empty safe filename stem." >&2; exit 1;;
  *[!A-Za-z0-9._-]*) echo "Error: output basename may contain only letters, numbers, '.', '_' and '-'." >&2; exit 1;;
esac

if [ -n "$INTERVAL_SHARD" ]; then
  case "$INTERVAL_SHARD" in
    .|..) echo "Error: interval shard must be a safe contig label." >&2; exit 1;;
    *[!A-Za-z0-9._-]*) echo "Error: interval shard may contain only letters, numbers, '.', '_' and '-'." >&2; exit 1;;
  esac
fi

if [ "$COHORT_STAGE" = "shard" ] && [ -z "$INTERVAL_SHARD" ]; then
  echo "Error: --cohort-stage shard requires --interval-shard." >&2
  exit 1
fi

if [ "$COHORT_STAGE" = "finalize" ] && [ -z "$INPUT_VCF" ]; then
  echo "Error: --cohort-stage finalize requires --input-vcf." >&2
  exit 1
fi

if [ "$COHORT_STAGE" = "finalize" ] && [ -n "$INTERVAL_SHARD" ]; then
  echo "Error: --cohort-stage finalize does not use --interval-shard." >&2
  exit 1
fi

if [ "$COHORT_STAGE" = "shard" ] && [ -n "$INPUT_VCF" ]; then
  echo "Error: --cohort-stage shard does not use --input-vcf." >&2
  exit 1
fi

if [ "$COHORT_STAGE" != "finalize" ] && [ -z "$SAMPLE_MAP" ]; then
  echo "Error: sample_map is required." >&2
  usage
fi

# Load env
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CBICALL_ENV_FILE:-$BINDIR/env.sh}"

if [[ "$COHORT_STAGE" != "finalize" && "${ARCH:-$(uname -m)}" =~ ^(aarch64|arm64)$ ]]; then
  echo "Error: GATK GenomicsDBImport cannot run on ARM/aarch64 with the bundled GATK 4.6 GenomicsDB native libraries." >&2
  echo "Run cohort_stage=all/shard on x86_64, or run cohort_stage=finalize on ARM using a gathered raw VCF created on x86_64." >&2
  exit 1
fi

# Prepare output directories and logging
dir=$(pwd)
GENOMICSDBDIR=$dir/01_genomicsdb
VARCALLDIR=$dir/02_varcall
LOGDIR=$dir/logs
mkdir -p "$GENOMICSDBDIR"
mkdir -p "$VARCALLDIR"
mkdir -p "$LOGDIR"

LOG="$LOGDIR/cohort_joint_genotyping.log"

# pipeline mode
# Convert to Uppercase
PIPELINE=${PIPELINE^^}
if [[ "$PIPELINE" != "WES" && "$PIPELINE" != "WGS" ]]; then
  echo "Error: pipeline must be 'wes' or 'wgs'." >&2; exit 1
fi

# sample count & workspace default naming
SAMPLE_COUNT=0
GENOTYPE_ANNOTATION_EXCLUDE_ARG=""
WORKSPACE_PATH=""
if [ "$COHORT_STAGE" != "finalize" ]; then
  if [ ! -s "$SAMPLE_MAP" ]; then
    echo "Error: sample_map '$SAMPLE_MAP' not found or empty." >&2; exit 1
  fi
  SAMPLE_COUNT=$(awk 'NF {n++} END {print n+0}' "$SAMPLE_MAP")
  if [ "$SAMPLE_COUNT" -lt 10 ]; then
    GENOTYPE_ANNOTATION_EXCLUDE_ARG="-AX InbreedingCoeff"
  fi
  if [ -z "$WORKSPACE" ]; then
    WORKSPACE="cohort.genomicsdb.${SAMPLE_COUNT}"
  fi
  case "$WORKSPACE" in
    /*) WORKSPACE_PATH="$WORKSPACE" ;;
    *) WORKSPACE_PATH="$GENOMICSDBDIR/$WORKSPACE" ;;
  esac
fi

function make_wes_shard_interval_list {
  local source_list="$1"
  local shard="$2"
  local out_list="$3"
  awk -v shard="$shard" '
    /^@/ { print; next }
    $1 == shard { print; n++ }
    END {
      if (n == 0) {
        printf("No intervals found for shard %s in %s\n", shard, FILENAME) > "/dev/stderr"
        exit 2
      }
    }
  ' "$source_list" > "$out_list"
}

function make_wgs_shard_interval_list {
  local ref_dict="$1"
  local shard="$2"
  local out_list="$3"
  awk -v shard="$shard" '
    /^@/ {
      print
      if ($1 == "@SQ") {
        sn = ""; ln = ""
        for (i = 1; i <= NF; i++) {
          if ($i ~ /^SN:/) sn = substr($i, 4)
          if ($i ~ /^LN:/) ln = substr($i, 4)
        }
        if (sn == shard && ln != "") {
          interval = sn "\t1\t" ln "\t+\t" sn
          found = 1
        }
      }
      next
    }
    END {
      if (!found) {
        printf("No reference contig found for shard %s in %s\n", shard, FILENAME) > "/dev/stderr"
        exit 2
      }
      print interval
    }
  ' "$ref_dict" > "$out_list"
}

# Set interval argument for WES vs WGS.
# GenomicsDBImport requires explicit intervals. WGS derives whole-contig
# intervals from the reference dictionary inside 01_genomicsdb.
INTERVAL_ARG=""
MERGE_INTERVALS_ARG=""
if [ "$COHORT_STAGE" = "finalize" ]; then
  echo "Finalize stage: using gathered raw VCF $INPUT_VCF" | tee -a "$LOG"
elif [ "$PIPELINE" = "WES" ]; then
  WES_INTERVAL_LIST="${CBICALL_INTERVAL_LIST:-$INTERVAL_LIST}"
  if [ -z "${WES_INTERVAL_LIST:-}" ] || [ ! -f "$WES_INTERVAL_LIST" ]; then
    echo "Error: WES interval list is not set or not found: ${WES_INTERVAL_LIST:-<unset>}" >&2
    exit 1
  fi
  if [ -n "$INTERVAL_SHARD" ]; then
    WES_SHARD_INTERVAL_LIST="$GENOMICSDBDIR/wes.${INTERVAL_SHARD}.interval_list"
    make_wes_shard_interval_list "$WES_INTERVAL_LIST" "$INTERVAL_SHARD" "$WES_SHARD_INTERVAL_LIST"
    INTERVAL_ARG="-L $WES_SHARD_INTERVAL_LIST"
    echo "WES mode: generated shard interval list $WES_SHARD_INTERVAL_LIST from $WES_INTERVAL_LIST" | tee -a "$LOG"
  else
    INTERVAL_ARG="-L $WES_INTERVAL_LIST"
  fi
  MERGE_INTERVALS_ARG="--merge-input-intervals true"
  if [ -n "${CBICALL_INTERVAL_LIST:-}" ]; then
    echo "WES mode: restricting to $WES_INTERVAL_LIST (CBICALL_INTERVAL_LIST override)" | tee -a "$LOG"
  else
    echo "WES mode: restricting to $WES_INTERVAL_LIST" | tee -a "$LOG"
  fi
else
  if [ -z "${REF_DICT:-}" ] || [ ! -f "$REF_DICT" ]; then
    echo "Error: REF_DICT is not set or not found (mode=WGS)." >&2
    exit 1
  fi
  if [ -n "$INTERVAL_SHARD" ]; then
    WGS_INTERVAL_LIST="$GENOMICSDBDIR/wgs.${INTERVAL_SHARD}.interval_list"
    make_wgs_shard_interval_list "$REF_DICT" "$INTERVAL_SHARD" "$WGS_INTERVAL_LIST"
    echo "WGS mode: generated shard interval list $WGS_INTERVAL_LIST from $REF_DICT" | tee -a "$LOG"
  else
    WGS_INTERVAL_LIST="$GENOMICSDBDIR/wgs.whole_genome.interval_list"
    awk '
      /^@/ {
        print
        if ($1 == "@SQ") {
          sn = ""; ln = ""
          for (i = 1; i <= NF; i++) {
            if ($i ~ /^SN:/) sn = substr($i, 4)
            if ($i ~ /^LN:/) ln = substr($i, 4)
          }
          if (sn != "" && ln != "") intervals[++n] = sn "\t1\t" ln "\t+\t" sn
        }
        next
      }
      END {
        for (i = 1; i <= n; i++) print intervals[i]
      }
    ' "$REF_DICT" > "$WGS_INTERVAL_LIST"
    echo "WGS mode: generated whole-genome intervals from $REF_DICT" | tee -a "$LOG"
  fi
  INTERVAL_ARG="-L $WGS_INTERVAL_LIST"
  MERGE_INTERVALS_ARG=""
fi

# Derived output names
COHORT_RAW_VCF="${OUTPUT_BASENAME}.gv.raw.vcf.gz"
COHORT_VQSR_SNP="${OUTPUT_BASENAME}.snp.recal.vcf.gz"
COHORT_SNP_TRANCHES="${OUTPUT_BASENAME}.snp.tranches.txt"
COHORT_VQSR_INDEL="${OUTPUT_BASENAME}.indel.recal.vcf.gz"
COHORT_INDEL_TRANCHES="${OUTPUT_BASENAME}.indel.tranches.txt"
COHORT_POST_SNP="${OUTPUT_BASENAME}.post_snp.vcf.gz"
COHORT_POST_VQSR="${OUTPUT_BASENAME}.vqsr.vcf.gz"
COHORT_QC_VCF="${OUTPUT_BASENAME}.gv.QC.vcf.gz"

case "$COHORT_STAGE" in
  all) STAGE_ACTION="full cohort run: import, genotype, and global filtering" ;;
  shard) STAGE_ACTION="shard run: import and genotype one interval shard" ;;
  finalize) STAGE_ACTION="finalize run: globally filter a gathered raw cohort VCF" ;;
esac

if [ "$COHORT_STAGE" = "finalize" ]; then
  SAMPLE_COUNT_DISPLAY="not applicable (finalize stage)"
  SAMPLE_MAP_DISPLAY="not used (finalize stage)"
  WORKSPACE_DISPLAY="not used (finalize stage)"
else
  SAMPLE_COUNT_DISPLAY="$SAMPLE_COUNT"
  SAMPLE_MAP_DISPLAY="${SAMPLE_MAP:-<none>}"
  WORKSPACE_DISPLAY="${WORKSPACE_PATH:-<none>}"
fi

{
  echo "## Cohort GenomicsDBImport -> Genotype -> VQSR/Hard-filter"
  echo "cohort_stage: $COHORT_STAGE"
  echo "stage_action: $STAGE_ACTION"
  echo "sample_map: $SAMPLE_MAP_DISPLAY"
  echo "pipeline: $PIPELINE"
  echo "sample_count: $SAMPLE_COUNT_DISPLAY"
  echo "workspace: $WORKSPACE_DISPLAY"
  echo "output_basename: $OUTPUT_BASENAME"
  echo "interval_shard: ${INTERVAL_SHARD:-<none>}"
  if [ "$COHORT_STAGE" = "finalize" ]; then
    echo "input_vcf: $INPUT_VCF"
    echo "final_vcf: $COHORT_QC_VCF"
  else
    echo "out_vcf: $COHORT_RAW_VCF"
  fi
  echo "tmpdir: $TMPDIR"
  echo "log: $LOG"
  echo ""
} | tee -a "$LOG"

# -----------------------------------------------------------------------------
# Step 1: GenomicsDBImport
# -----------------------------------------------------------------------------
if [ -z "${REF:-}" ]; then
  echo "Error: REF not set (expected to be defined in env.sh)." >&2
  exit 1
fi

if [ "$COHORT_STAGE" != "finalize" ]; then
  echo ">>> Step 1: GenomicsDBImport" | tee -a "$LOG"
  mkdir -p "$(dirname "$WORKSPACE_PATH")"

  set -x
  "$GATK4_BIN" $GATK4_JAVA_OPTS_64G GenomicsDBImport \
    --sample-name-map "$SAMPLE_MAP" \
    --genomicsdb-workspace-path "$WORKSPACE_PATH" \
    $MERGE_INTERVALS_ARG \
    $INTERVAL_ARG \
    --tmp-dir "$TMPDIR" \
    2>> "$LOG"
  set +x
  echo ok > "$GENOMICSDBDIR/genomicsdbimport.done"

  # -----------------------------------------------------------------------------
  # Step 2: GenotypeGVCFs
  # -----------------------------------------------------------------------------
  echo ">>> Step 2: GenotypeGVCFs" | tee -a "$LOG"
  cd "$VARCALLDIR"

  set -x
  "$GATK4_BIN" $GATK4_JAVA_OPTS_64G GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://$WORKSPACE_PATH" \
    -O "$COHORT_RAW_VCF" \
    --stand-call-conf 10 \
    $GENOTYPE_ANNOTATION_EXCLUDE_ARG \
    --tmp-dir "$TMPDIR" \
    $INTERVAL_ARG \
    2>> "$LOG"
  set +x

  if [ $? -ne 0 ]; then
    echo "ERROR: GenotypeGVCFs failed. See log: $LOG" >&2
    exit 1
  fi
  echo "Genotyping completed. Raw cohort VCF: $COHORT_RAW_VCF" | tee -a "$LOG"

  if [ "$COHORT_STAGE" = "shard" ]; then
    echo "Shard stage complete. Raw cohort shard VCF: $COHORT_RAW_VCF" | tee -a "$LOG"
    exit 0
  fi
  RAW_VCF_FOR_FILTERING="$COHORT_RAW_VCF"
else
  if [ ! -s "$INPUT_VCF" ]; then
    echo "Error: input_vcf '$INPUT_VCF' not found or empty." >&2
    exit 1
  fi
  case "$INPUT_VCF" in
    /*) ;;
    *) INPUT_VCF="$(cd "$(dirname "$INPUT_VCF")" && pwd)/$(basename "$INPUT_VCF")" ;;
  esac
  cd "$VARCALLDIR"
  RAW_VCF_FOR_FILTERING="$INPUT_VCF"
fi

# -----------------------------------------------------------------------------
# Step 3: Count SNPs/INDELs and decide on VQSR
# -----------------------------------------------------------------------------
echo ">>> Step 3: Count variants and decide on VQSR" | tee -a "$LOG"
nSNP=$(zgrep -v '^#' "$RAW_VCF_FOR_FILTERING" | awk 'length($5)==1' | wc -l)
nINDEL=$(zgrep -v '^#' "$RAW_VCF_FOR_FILTERING" | awk 'length($5)!=1' | wc -l)
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
    -V "$RAW_VCF_FOR_FILTERING" \
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
    -V "$RAW_VCF_FOR_FILTERING" \
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
tmp_vcf="$RAW_VCF_FOR_FILTERING"
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
  --filter-name "QD2"        --filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
  --filter-name "FS60"       --filter-expression "FS > 60.0" \
  --filter-name "MQ40"       --filter-expression "MQ < 40.0" \
  --filter-name "MQRS-12.5"  --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
  --filter-name "RPRS-8"     --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
  --filter-name "QD2_indel"  --filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
  --filter-name "FS200"      --filter-expression "FS > 200.0" \
  --filter-name "RPRS-20"    --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
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
echo "All done. Cohort raw: $RAW_VCF_FOR_FILTERING ; Cohort QC: $COHORT_QC_VCF" | tee -a "$LOG"
exit 0
