#!/usr/bin/env bash
#
#   WGS/WES pipeline Bash script (GATK4.6 best-practices edition).
#   Last Modified: 2025-05-02
#
# README:
#   This script runs either WES or WGS pipelines based on the '-p' mode flag.
#   Mode 'wes': uses provided exome interval list to restrict variant calling.
#   Mode 'wgs': processes the entire genome (no intervals).
#   Common steps (GATK Best Practices):
#     1) Align & add read groups
#     2) Merge lane-level BAMs
#     3) Mark duplicates
#     4) Base Quality Score Recalibration (BQSR)
#     5) HaplotypeCaller (per-sample gVCF)
#     6) GenotypeGVCFs (raw VCF)
#     7) VariantRecalibrator (VQSR) if variant counts suffice
#     8) Apply VQSR models or fallback to hard filters
#     9) Hard-filter & write QC VCF
#    10) Coverage stats & sex determination

set -eu

function usage {
    cat <<EOF
Usage: $0 -t <n_threads> -p <wes|wgs>
  -t THREADS   Number of CPU threads for GATK tools
  -p PIPELINE  Pipeline mode: 'wes' or 'wgs'
EOF
    exit 1
}

# Parse command-line arguments
if [ $# -ne 4 ]; then usage; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--threads)     THREADS="$2"; shift 2;;
    -p|--pipeline)    PIPELINE="$2";     shift 2;;
    *)            usage;;
  esac
done
# Convert to Uppercase
PIPELINE=${PIPELINE^^}
if [[ "$PIPELINE" != "WES" && "$PIPELINE" != "WGS" ]]; then
  echo "Error: Mode must be 'wes' or 'wgs'." >&2
  usage
fi

# Load GATK parameters (reference, known-sites, interval list, tool paths)
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$BINDIR/parameters.sh"

# Prepare output directories
dir=$(pwd)
BAMDIR=$dir/01_bam
VARCALLDIR=$dir/02_varcall
STATSDIR=$dir/03_stats
LOGDIR=$dir/logs
mkdir -p "$BAMDIR" "$VARCALLDIR" "$STATSDIR" "$LOGDIR"

# Determine sample ID & log file
rawid=$(basename "$(dirname "$PWD")")
id=${rawid%%_*}
LOG=$LOGDIR/${id}.log

# Set interval argument for WES vs WGS
if [[ "$PIPELINE" == "WES" ]]; then
  INTERVAL_ARG=( -L "$INTERVAL_LIST" )
  echo "Running in WES mode: restricting to target intervals."
else
  INTERVAL_ARG=()
  echo "Running in WGS mode: processing whole genome."
fi

#------------------------------------------------------------------------------
# STEP 1: Align & Add Read Groups (per-lane BAMs)
# Theory: Read groups are required for GATK to tag sample, library, and platform
# information, enabling per-sample metrics and proper handling of multi-library data.
#------------------------------------------------------------------------------
echo ">>> STEP 1: Align & add read groups"
for R1 in ../*R1*fastq.gz; do
  fn=$(basename "$R1" .fastq.gz)
  base=${fn%_R1*}
  R2=${R1/_R1_/_R2_}

  SAMPLE=$(echo "$base" | cut -d'_' -f1-2)
  LANE=$(echo   "$base" | cut -d'_' -f3)
  # ensure RGIDs are unique by appending a timestamp
  RGID="${SAMPLE}.${LANE}.$(date +%s)"
  RGPU="${SAMPLE}.${LANE}.unit1"

  out_bam="$BAMDIR/${base}.rg.bam"
  echo "Aligning $fn -> $(basename "$out_bam")"

  # Align
  $BWA mem -M -t "$THREADS" "$REFGZ" "$R1" "$R2" \
    # drop secondary (0x100) & supplementary (0x800) alignments before RG tagging
    # | $SAM view -b -F 0x900 - \
    | $GATK4_CMD AddOrReplaceReadGroups \
        --INPUT /dev/stdin \
        --OUTPUT "$out_bam" \
        --TMP_DIR "$TMPDIR" \
        --RGPL ILLUMINA \
        --RGLB sureselect \
        --RGSM "$SAMPLE" \
        --RGID "$RGID" \
        --RGPU "$RGPU" \
    2>> "$LOG"
done

#------------------------------------------------------------------------------
# STEP 2: Merge Lane-level BAMs
# Theory: Merging ensures downstream duplicate marking and BQSR
# operate on the full sample, minimizing lane-specific biases.
#------------------------------------------------------------------------------
echo ">>> STEP 2: Merge lane-level BAMs"
rg_bams=( $BAMDIR/*.rg.bam )
$GATK4_CMD MergeSamFiles \
  ${rg_bams[@]/#/-I } \
  -O "$BAMDIR/${id}.rg.merged.bam" \
  --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --TMP_DIR "$TMPDIR" \
  2>> "$LOG"

#------------------------------------------------------------------------------
# STEP 3: Mark Duplicates
# Theory: PCR/optical duplicates inflate support for artefactual variants,
# so marking them reduces false positives in variant calling.
#------------------------------------------------------------------------------
echo ">>> STEP 3: Mark duplicates"
$GATK4_CMD MarkDuplicates \
  -I "$BAMDIR/${id}.rg.merged.bam" \
  -O "$BAMDIR/${id}.rg.merged.dedup.bam" \
  --METRICS_FILE "$BAMDIR/${id}.rg.merged.dedup.metrics.txt" \
  --CREATE_INDEX true --TMP_DIR "$TMPDIR" \
  2>> "$LOG"
$SAM index "$BAMDIR/${id}.rg.merged.dedup.bam"

#------------------------------------------------------------------------------
# STEP 4: Base Quality Score Recalibration (BQSR)
# Theory: Models systematic error in reported base qualities using known-sites,
# improving the accuracy of variant calls by correcting quality scores.
#------------------------------------------------------------------------------
echo ">>> STEP 4: Base recalibration"
bqsr_table="$BAMDIR/${id}.rg.merged.dedup.recal.table"
$GATK4_CMD BaseRecalibrator \
  -R "$REF" \
  -I "$BAMDIR/${id}.rg.merged.dedup.bam" \
  --known-sites "$dbSNP" --known-sites "$MILLS_INDELS" --known-sites "$KG_INDELS" \
  -O "$bqsr_table" --tmp-dir "$TMPDIR" 2>> "$LOG"

$GATK4_CMD ApplyBQSR \
  -R "$REF" -I "$BAMDIR/${id}.rg.merged.dedup.bam" \
  --bqsr-recal-file "$bqsr_table" \
  -O "$BAMDIR/${id}.rg.merged.dedup.recal.bam" --tmp-dir "$TMPDIR" 2>> "$LOG"
$SAM index "$BAMDIR/${id}.rg.merged.dedup.recal.bam"

#------------------------------------------------------------------------------
# STEP 5: HaplotypeCaller → gVCF
# Theory: Local de Bruijn graph assembly per active region improves indel calling;
# emits reference-confidence gVCF required for joint genotyping.
#------------------------------------------------------------------------------
echo ">>> STEP 5: HaplotypeCaller -> gVCF"
$GATK4_CMD HaplotypeCaller \
  -R "$REF" -I "$BAMDIR/${id}.rg.merged.dedup.recal.bam" \
  -O "$VARCALLDIR/${id}.hc.g.vcf.gz" \
  "${INTERVAL_ARG[@]}" \
  --native-pair-hmm-threads "$THREADS" -ERC GVCF 2>> "$LOG"

#------------------------------------------------------------------------------
# STEP 6: GenotypeGVCFs → Raw VCF
# Theory: Jointly genotype one or more gVCFs to produce high-confidence sample VCF.
#------------------------------------------------------------------------------
echo ">>> STEP 6: GenotypeGVCFs -> raw VCF"
$GATK4_CMD GenotypeGVCFs \
  -R "$REF" -V "$VARCALLDIR/${id}.hc.g.vcf.gz" \
  -O "$VARCALLDIR/${id}.hc.raw.vcf.gz" --stand-call-conf 10 2>> "$LOG"

#------------------------------------------------------------------------------
# STEP 7: Build VQSR Models
# Theory: Use Gaussian mixture models on variant annotations
# to distinguish true variants from artefacts; requires sufficient counts.
#------------------------------------------------------------------------------
echo ">>> STEP 7: Optional VQSR if enough variants"
nSNP=$(zgrep -v '^#' "$VARCALLDIR/${id}.hc.raw.vcf.gz" | awk 'length($5)==1' | wc -l)
nINDEL=$(zgrep -v '^#' "$VARCALLDIR/${id}.hc.raw.vcf.gz" | awk 'length($5)!=1' | wc -l)
minSNP=1000; minINDEL=8000; apply_snp=false; apply_indel=false

if (( nSNP >= minSNP )); then
  $GATK4_CMD VariantRecalibrator \
    -R "$REF" \
    -V "$VARCALLDIR/${id}.hc.raw.vcf.gz" \
    $SNP_RES \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
    --mode SNP \
    -O "$VARCALLDIR/${id}.snp.recal.vcf.gz" \
    --tranches-file "$VARCALLDIR/${id}.snp.tranches.txt" \
    --max-gaussians 6 2>> "$LOG"
    #--rscript-file "$VARCALLDIR/${id}.snp.plots.R"
  apply_snp=true
fi

if (( nINDEL >= minINDEL )); then
  $GATK4_CMD VariantRecalibrator \
    -R "$REF" \
    -V "$VARCALLDIR/${id}.hc.raw.vcf.gz" \
    $INDEL_RES \
    -an QD -an FS -an ReadPosRankSum \
    --mode INDEL \
    -O "$VARCALLDIR/${id}.indel.recal.vcf.gz" \
    --tranches-file "$VARCALLDIR/${id}.indel.tranches.txt" \
    --max-gaussians 4 2>> "$LOG"
    #--rscript-file "$VARCALLDIR/${id}.indel.plots.R"
  apply_indel=true
fi

#------------------------------------------------------------------------------
# STEP 8: Apply VQSR Models or Fallback
#------------------------------------------------------------------------------
echo ">>> STEP 8: Apply VQSR or fallback"
tmp_vcf="$VARCALLDIR/${id}.hc.raw.vcf.gz"
if [ "$apply_snp" = true ]; then
  $GATK4_CMD ApplyVQSR \
    -R "$REF" -V "$tmp_vcf" \
    --recal-file "$VARCALLDIR/${id}.snp.recal.vcf.gz" \
    --tranches-file "$VARCALLDIR/${id}.snp.tranches.txt" \
    --mode SNP --truth-sensitivity-filter-level 99.0 \
    -O "$VARCALLDIR/${id}.hc.post_snp.vcf.gz" 2>> "$LOG"
  tmp_vcf="$VARCALLDIR/${id}.hc.post_snp.vcf.gz"
fi
if [ "$apply_indel" = true ]; then
  $GATK4_CMD ApplyVQSR \
    -R "$REF" -V "$tmp_vcf" \
    --recal-file "$VARCALLDIR/${id}.indel.recal.vcf.gz" \
    --tranches-file "$VARCALLDIR/${id}.indel.tranches.txt" \
    --mode INDEL --truth-sensitivity-filter-level 95.0 \
    -O "$VARCALLDIR/${id}.hc.vqsr.vcf.gz" 2>> "$LOG"
  tmp_vcf="$VARCALLDIR/${id}.hc.vqsr.vcf.gz"
fi

#------------------------------------------------------------------------------
# STEP 9: Hard-Filter & Write QC VCF
# Theory: Apply recommended hard filters on annotations when VQSR isn't applied or as QC.
#   SNPs: QD<2.0, FS>60.0, MQ<40.0, MQRankSum<-12.5, ReadPosRankSum<-8.0
#   INDELs: QD<2.0, FS>200.0, ReadPosRankSum<-20.0
#------------------------------------------------------------------------------
echo ">>> STEP 9: Hard-filter & write QC.vcf"
$GATK4_CMD VariantFiltration \
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
  -O "$VARCALLDIR/${id}.hc.QC.vcf.gz" \
  2>> "$LOG"

#------------------------------------------------------------------------------
# STEP 10: Coverage Stats & Sex Determination
#------------------------------------------------------------------------------
echo ">>> STEP 10: Coverage & Sex Determination" 2>> "$LOG"

# choose chromosome 1 naming
if [[ "$REF" == *b37*.fasta ]]; then
  chrN=1
else
  chrN=chr1
fi

bam_raw="$BAMDIR/${id}.rg.merged.dedup.bam"
bam_recal="$BAMDIR/${id}.rg.merged.dedup.recal.bam"
out_raw="$STATSDIR/${chrN}.raw.bam"
out_dedup="$STATSDIR/${chrN}.dedup.bam"

# extract chr1 (or “1”) reads and redirect errors
$SAM view -b "$bam_raw" "$chrN"  > "$out_raw"   2>> "$LOG"
$SAM view -b "$bam_recal" "$chrN"  > "$out_dedup" 2>> "$LOG"

# index them, also logging STDERR
$SAM index "$out_raw"    2>> "$LOG"
$SAM index "$out_dedup"  2>> "$LOG"

# run coverage & sex scripts: STDOUT → stats file, STDERR → main log
"$BINDIR"/coverage.sh "$id" "$out_raw" "$out_dedup" "$PIPELINE" \
    > "$STATSDIR/${id}.coverage.txt" 2>> "$LOG"
"$BINDIR"/vcf2sex.sh "$VARCALLDIR/${id}.hc.QC.vcf.gz" \
    > "$STATSDIR/${id}.sex.txt" 2>> "$LOG"

echo "All done! QC VCF: $VARCALLDIR/${id}.hc.QC.vcf.gz"
