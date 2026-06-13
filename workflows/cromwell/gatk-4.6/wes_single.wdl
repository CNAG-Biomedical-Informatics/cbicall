version 1.0

workflow CBIcallWesSingle {
  input {
    String id
    String pipeline
    String genome
    Int threads
    Boolean cleanup_bam
    String qc_coverage_region
    File fastq_pairs_tsv
    String tmpdir
    String bwa
    String samtools
    String gatk4_cmd
    String ref
    String refgz
    String dbsnp
    String mills_indels
    String kg_indels
    String hapmap
    String omni
    String interval_list
    String snp_res
    String indel_res
    String coverage_script
    String vcf2sex_script
    String vcf2hash_script
  }

  call RunWesSingle {
    input:
      id = id,
      pipeline = pipeline,
      genome = genome,
      threads = threads,
      cleanup_bam = cleanup_bam,
      qc_coverage_region = qc_coverage_region,
      fastq_pairs_tsv = fastq_pairs_tsv,
      tmpdir = tmpdir,
      bwa = bwa,
      samtools = samtools,
      gatk4_cmd = gatk4_cmd,
      ref = ref,
      refgz = refgz,
      dbsnp = dbsnp,
      mills_indels = mills_indels,
      kg_indels = kg_indels,
      hapmap = hapmap,
      omni = omni,
      interval_list = interval_list,
      snp_res = snp_res,
      indel_res = indel_res,
      coverage_script = coverage_script,
      vcf2sex_script = vcf2sex_script,
      vcf2hash_script = vcf2hash_script
  }

  output {
    File gvcf = RunWesSingle.gvcf
    File raw_vcf = RunWesSingle.raw_vcf
    File qc_vcf = RunWesSingle.qc_vcf
    File coverage = RunWesSingle.coverage
    File sex = RunWesSingle.sex
    File vcf_hash = RunWesSingle.vcf_hash
    Array[File] logs = RunWesSingle.logs
  }
}

task RunWesSingle {
  input {
    String id
    String pipeline
    String genome
    Int threads
    Boolean cleanup_bam
    String qc_coverage_region
    File fastq_pairs_tsv
    String tmpdir
    String bwa
    String samtools
    String gatk4_cmd
    String ref
    String refgz
    String dbsnp
    String mills_indels
    String kg_indels
    String hapmap
    String omni
    String interval_list
    String snp_res
    String indel_res
    String coverage_script
    String vcf2sex_script
    String vcf2hash_script
  }

  command <<<
    #!/usr/bin/env bash
    set -eu
    export TMPDIR="~{tmpdir}"
    export LC_ALL=C
    export GATK_DISABLE_AUTO_S3_UPLOAD=true

    BAMDIR="01_bam"
    VARCALLDIR="02_varcall"
    STATSDIR="03_stats"
    LOGDIR="logs"
    mkdir -p "$BAMDIR" "$VARCALLDIR" "$STATSDIR" "$LOGDIR"
    LOG="$LOGDIR/~{id}.log"
    sample_name="~{id}"

    if [ "~{pipeline}" = "wes" ]; then
      INTERVAL_ARG="-L ~{interval_list}"
      echo "WES mode: restricting to ~{interval_list}"
    else
      INTERVAL_ARG=""
      echo "WGS mode: processing whole genome"
    fi

    echo "STEP 1: Align and add read groups"
    while IFS=$'\t' read -r base R1 R2; do
      [ -n "$base" ] || continue
      RGID="$sample_name.$base"
      RGPU="$sample_name.$base.unit1"
      out_bam="$BAMDIR/$base.rg.bam"
      echo "Aligning $base to $out_bam"
      ~{bwa} mem -M -t ~{threads} "~{refgz}" "$R1" "$R2" \
        | ~{gatk4_cmd} AddOrReplaceReadGroups \
            --INPUT /dev/stdin \
            --OUTPUT "$out_bam" \
            --TMP_DIR "~{tmpdir}" \
            --RGPL ILLUMINA \
            --RGLB sureselect \
            --RGSM "$sample_name" \
            --RGID "$RGID" \
            --RGPU "$RGPU" \
        2>> "$LOG"
    done < "~{fastq_pairs_tsv}"

    echo "STEP 2: Merge lane-level BAMs"
    merge_inputs=""
    for rg_bam in "$BAMDIR"/*.rg.bam; do
      merge_inputs="$merge_inputs -I $rg_bam"
    done
    ~{gatk4_cmd} MergeSamFiles \
      $merge_inputs \
      -O "$BAMDIR/~{id}.rg.merged.bam" \
      --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --TMP_DIR "~{tmpdir}" \
      2>> "$LOG"

    echo "STEP 3: Mark duplicates"
    ~{gatk4_cmd} MarkDuplicates \
      -I "$BAMDIR/~{id}.rg.merged.bam" \
      -O "$BAMDIR/~{id}.rg.merged.dedup.bam" \
      --METRICS_FILE "$BAMDIR/~{id}.rg.merged.dedup.metrics.txt" \
      --CREATE_INDEX true --TMP_DIR "~{tmpdir}" \
      2>> "$LOG"
    ~{samtools} index "$BAMDIR/~{id}.rg.merged.dedup.bam"

    echo "STEP 4: Base recalibration"
    bqsr_table="$BAMDIR/~{id}.rg.merged.dedup.recal.table"
    ~{gatk4_cmd} BaseRecalibrator \
      -R "~{ref}" \
      -I "$BAMDIR/~{id}.rg.merged.dedup.bam" \
      --known-sites "~{dbsnp}" --known-sites "~{mills_indels}" --known-sites "~{kg_indels}" \
      -O "$bqsr_table" --tmp-dir "~{tmpdir}" 2>> "$LOG"

    ~{gatk4_cmd} ApplyBQSR \
      -R "~{ref}" -I "$BAMDIR/~{id}.rg.merged.dedup.bam" \
      --bqsr-recal-file "$bqsr_table" \
      -O "$BAMDIR/~{id}.rg.merged.dedup.recal.bam" --tmp-dir "~{tmpdir}" 2>> "$LOG"
    ~{samtools} index "$BAMDIR/~{id}.rg.merged.dedup.recal.bam"

    echo "STEP 5: HaplotypeCaller to gVCF"
    ~{gatk4_cmd} HaplotypeCaller \
      -R "~{ref}" -I "$BAMDIR/~{id}.rg.merged.dedup.recal.bam" \
      -O "$VARCALLDIR/~{id}.hc.g.vcf.gz" \
      $INTERVAL_ARG \
      --native-pair-hmm-threads ~{threads} -ERC GVCF 2>> "$LOG"

    echo "STEP 6: GenotypeGVCFs to raw VCF"
    ~{gatk4_cmd} GenotypeGVCFs \
      -R "~{ref}" -V "$VARCALLDIR/~{id}.hc.g.vcf.gz" \
      -O "$VARCALLDIR/~{id}.hc.raw.vcf.gz" --stand-call-conf 10 2>> "$LOG"

    echo "STEP 7: Optional VQSR if enough variants"
    nSNP=$(zgrep -v '^#' "$VARCALLDIR/~{id}.hc.raw.vcf.gz" | awk 'length($5)==1' | wc -l)
    nINDEL=$(zgrep -v '^#' "$VARCALLDIR/~{id}.hc.raw.vcf.gz" | awk 'length($5)!=1' | wc -l)
    minSNP=1000; minINDEL=8000; apply_snp=false; apply_indel=false

    if (( nSNP >= minSNP )); then
      ~{gatk4_cmd} VariantRecalibrator \
        -R "~{ref}" \
        -V "$VARCALLDIR/~{id}.hc.raw.vcf.gz" \
        ~{snp_res} \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
        --mode SNP \
        -O "$VARCALLDIR/~{id}.hc.snp.recal.vcf.gz" \
        --tranches-file "$VARCALLDIR/~{id}.hc.snp.tranches.txt" \
        --max-gaussians 6 2>> "$LOG"
      apply_snp=true
    fi

    if (( nINDEL >= minINDEL )); then
      ~{gatk4_cmd} VariantRecalibrator \
        -R "~{ref}" \
        -V "$VARCALLDIR/~{id}.hc.raw.vcf.gz" \
        ~{indel_res} \
        -an QD -an FS -an ReadPosRankSum \
        --mode INDEL \
        -O "$VARCALLDIR/~{id}.hc.indel.recal.vcf.gz" \
        --tranches-file "$VARCALLDIR/~{id}.hc.indel.tranches.txt" \
        --max-gaussians 4 2>> "$LOG"
      apply_indel=true
    fi

    echo "STEP 8: Apply VQSR or fallback"
    tmp_vcf="$VARCALLDIR/~{id}.hc.raw.vcf.gz"

    if [ "$apply_snp" = true ]; then
      ~{gatk4_cmd} ApplyVQSR \
        -R "~{ref}" -V "$tmp_vcf" \
        --recal-file "$VARCALLDIR/~{id}.hc.snp.recal.vcf.gz" \
        --tranches-file "$VARCALLDIR/~{id}.hc.snp.tranches.txt" \
        --mode SNP --truth-sensitivity-filter-level 99.0 \
        -O "$VARCALLDIR/~{id}.hc.post_snp.vcf.gz" 2>> "$LOG"
      tmp_vcf="$VARCALLDIR/~{id}.hc.post_snp.vcf.gz"
    fi

    if [ "$apply_indel" = true ]; then
      ~{gatk4_cmd} ApplyVQSR \
        -R "~{ref}" -V "$tmp_vcf" \
        --recal-file "$VARCALLDIR/~{id}.hc.indel.recal.vcf.gz" \
        --tranches-file "$VARCALLDIR/~{id}.hc.indel.tranches.txt" \
        --mode INDEL --truth-sensitivity-filter-level 95.0 \
        -O "$VARCALLDIR/~{id}.hc.vqsr.vcf.gz" 2>> "$LOG"
      tmp_vcf="$VARCALLDIR/~{id}.hc.vqsr.vcf.gz"
    fi

    echo "STEP 9: Hard-filter and write QC.vcf"
    ~{gatk4_cmd} VariantFiltration \
      -R "~{ref}" \
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
      -O "$VARCALLDIR/~{id}.hc.QC.vcf.gz" \
      2>> "$LOG"

    echo "STEP 10: Coverage, sex determination, and VCF hash" 2>> "$LOG"
    requested_region="~{qc_coverage_region}"
    if awk -v chr="$requested_region" '$1==chr {found=1} END{exit !found}' "~{ref}.fai"; then
      chrN="$requested_region"
    elif [[ "$requested_region" == chr* ]] && awk -v chr="${requested_region#chr}" '$1==chr {found=1} END{exit !found}' "~{ref}.fai"; then
      chrN="${requested_region#chr}"
    elif [[ "$requested_region" != chr* ]] && awk -v chr="chr${requested_region}" '$1==chr {found=1} END{exit !found}' "~{ref}.fai"; then
      chrN="chr${requested_region}"
    else
      echo "Error: coverage region '$requested_region' not found in ~{ref}.fai" >&2
      exit 1
    fi
    chr_file=$(printf '%s' "$chrN" | sed 's/[^A-Za-z0-9_.-]/_/g')
    export CBICALL_COVERAGE_REGION="$chrN"

    bam_raw="$BAMDIR/~{id}.rg.merged.dedup.bam"
    bam_recal="$BAMDIR/~{id}.rg.merged.dedup.recal.bam"
    out_raw="$STATSDIR/$chr_file.raw.bam"
    out_dedup="$STATSDIR/$chr_file.dedup.bam"

    ~{samtools} view -b "$bam_raw" "$chrN" > "$out_raw" 2>> "$LOG"
    ~{samtools} view -b "$bam_recal" "$chrN" > "$out_dedup" 2>> "$LOG"
    ~{samtools} index "$out_raw" 2>> "$LOG"
    ~{samtools} index "$out_dedup" 2>> "$LOG"

    bash "~{coverage_script}" "~{id}" "$out_raw" "$out_dedup" "~{pipeline}" \
      > "$STATSDIR/~{id}.coverage.txt" 2>> "$LOG"
    bash "~{vcf2sex_script}" "$VARCALLDIR/~{id}.hc.QC.vcf.gz" \
      > "$STATSDIR/~{id}.sex.txt" 2>> "$LOG"
    bash "~{vcf2hash_script}" "$VARCALLDIR/~{id}.hc.QC.vcf.gz" \
      > "$STATSDIR/~{id}.vcf.sha256.txt" 2>> "$LOG"

    rm -f "$out_raw" "$out_dedup" "$out_raw.bai" "$out_dedup.bai"

    if [ "~{if cleanup_bam then "true" else "false"}" = true ]; then
      echo "STEP 11: Cleanup BAMs" 2>> "$LOG"
      rm -f "$BAMDIR"/*.{bam,bai} 2>> "$LOG" || true
    fi

    echo "All done! QC VCF: $VARCALLDIR/~{id}.hc.QC.vcf.gz"
  >>>

  output {
    File gvcf = "02_varcall/" + id + ".hc.g.vcf.gz"
    File raw_vcf = "02_varcall/" + id + ".hc.raw.vcf.gz"
    File qc_vcf = "02_varcall/" + id + ".hc.QC.vcf.gz"
    File coverage = "03_stats/" + id + ".coverage.txt"
    File sex = "03_stats/" + id + ".sex.txt"
    File vcf_hash = "03_stats/" + id + ".vcf.sha256.txt"
    Array[File] logs = glob("logs/*")
  }

  runtime {
    cpu: threads
  }
}
