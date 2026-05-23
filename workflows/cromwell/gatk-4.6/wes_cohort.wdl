version 1.0

workflow CBIcallCohort {
  input {
    String pipeline
    String genome
    Int threads
    File sample_map
    String workspace
    String tmpdir
    String gatk4_cmd
    String ref
    String dbsnp
    String mills_indels
    String hapmap
    String omni
    String interval_list
    String snp_res
    String indel_res
    String vcf2hash_script
    Int min_snp_for_vqsr = 1000
    Int min_indel_for_vqsr = 8000
  }

  call RunCohort {
    input:
      pipeline = pipeline,
      genome = genome,
      threads = threads,
      sample_map = sample_map,
      workspace = workspace,
      tmpdir = tmpdir,
      gatk4_cmd = gatk4_cmd,
      ref = ref,
      dbsnp = dbsnp,
      mills_indels = mills_indels,
      hapmap = hapmap,
      omni = omni,
      interval_list = interval_list,
      snp_res = snp_res,
      indel_res = indel_res,
      vcf2hash_script = vcf2hash_script,
      min_snp_for_vqsr = min_snp_for_vqsr,
      min_indel_for_vqsr = min_indel_for_vqsr
  }

  output {
    File raw_vcf = RunCohort.raw_vcf
    File qc_vcf = RunCohort.qc_vcf
    File vcf_hash = RunCohort.vcf_hash
    Array[File] logs = RunCohort.logs
  }
}

task RunCohort {
  input {
    String pipeline
    String genome
    Int threads
    File sample_map
    String workspace
    String tmpdir
    String gatk4_cmd
    String ref
    String dbsnp
    String mills_indels
    String hapmap
    String omni
    String interval_list
    String snp_res
    String indel_res
    String vcf2hash_script
    Int min_snp_for_vqsr
    Int min_indel_for_vqsr
  }

  command <<<
    #!/usr/bin/env bash
    set -eu
    export TMPDIR="~{tmpdir}"
    export LC_ALL=C
    export GATK_DISABLE_AUTO_S3_UPLOAD=true

    VARCALLDIR="02_varcall"
    STATSDIR="03_stats"
    LOGDIR="logs"
    mkdir -p "$VARCALLDIR" "$STATSDIR" "$LOGDIR"
    LOG="$LOGDIR/cohort_joint_genotyping.log"

    if [ "~{pipeline}" = "wes" ]; then
      INTERVAL_ARG="-L ~{interval_list}"
      echo "WES mode: restricting to ~{interval_list}" | tee -a "$LOG"
    else
      INTERVAL_ARG=""
      echo "WGS mode: processing whole genome" | tee -a "$LOG"
    fi

    cd "$VARCALLDIR"

    echo "Step 1: GenomicsDBImport" | tee -a "../$LOG"
    ~{gatk4_cmd} GenomicsDBImport \
      --sample-name-map "~{sample_map}" \
      --genomicsdb-workspace-path "~{workspace}" \
      --merge-input-intervals true \
      $INTERVAL_ARG \
      --tmp-dir "~{tmpdir}" \
      2>> "../$LOG"

    echo "Step 2: GenotypeGVCFs" | tee -a "../$LOG"
    ~{gatk4_cmd} GenotypeGVCFs \
      -R "~{ref}" \
      -V "gendb://~{workspace}" \
      -O cohort.gv.raw.vcf.gz \
      --stand-call-conf 10 \
      --tmp-dir "~{tmpdir}" \
      $INTERVAL_ARG \
      2>> "../$LOG"

    echo "Step 3: VQSR decision" | tee -a "../$LOG"
    nSNP=$(zgrep -v '^#' cohort.gv.raw.vcf.gz | awk 'length($5)==1' | wc -l | tr -d ' ')
    nINDEL=$(zgrep -v '^#' cohort.gv.raw.vcf.gz | awk 'length($5)!=1' | wc -l | tr -d ' ')
    apply_snp=false
    apply_indel=false
    tmp_vcf="cohort.gv.raw.vcf.gz"
    echo "Found SNPs: $nSNP ; INDELs: $nINDEL" | tee -a "../$LOG"

    if [ "$nSNP" -ge "~{min_snp_for_vqsr}" ]; then
      echo "Step 4: SNP VariantRecalibrator" | tee -a "../$LOG"
      ~{gatk4_cmd} VariantRecalibrator \
        -R "~{ref}" \
        -V cohort.gv.raw.vcf.gz \
        ~{snp_res} \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
        --mode SNP \
        -O cohort.snp.recal.vcf.gz \
        --tranches-file cohort.snp.tranches.txt \
        --max-gaussians 6 \
        --tmp-dir "~{tmpdir}" \
        2>> "../$LOG"
      apply_snp=true
    else
      echo "Skipping SNP VQSR (found $nSNP < min ~{min_snp_for_vqsr})" | tee -a "../$LOG"
    fi

    if [ "$nINDEL" -ge "~{min_indel_for_vqsr}" ]; then
      echo "Step 5: INDEL VariantRecalibrator" | tee -a "../$LOG"
      ~{gatk4_cmd} VariantRecalibrator \
        -R "~{ref}" \
        -V cohort.gv.raw.vcf.gz \
        ~{indel_res} \
        -an QD -an FS -an ReadPosRankSum \
        --mode INDEL \
        -O cohort.indel.recal.vcf.gz \
        --tranches-file cohort.indel.tranches.txt \
        --max-gaussians 4 \
        --tmp-dir "~{tmpdir}" \
        2>> "../$LOG"
      apply_indel=true
    else
      echo "Skipping INDEL VQSR (found $nINDEL < min ~{min_indel_for_vqsr})" | tee -a "../$LOG"
    fi

    echo "Step 6: Apply VQSR if available" | tee -a "../$LOG"
    if [ "$apply_snp" = true ]; then
      ~{gatk4_cmd} ApplyVQSR \
        -R "~{ref}" -V "$tmp_vcf" \
        --recal-file cohort.snp.recal.vcf.gz \
        --tranches-file cohort.snp.tranches.txt \
        --mode SNP --truth-sensitivity-filter-level 99.0 \
        -O cohort.post_snp.vcf.gz \
        --tmp-dir "~{tmpdir}" \
        2>> "../$LOG"
      tmp_vcf="cohort.post_snp.vcf.gz"
    fi

    if [ "$apply_indel" = true ]; then
      ~{gatk4_cmd} ApplyVQSR \
        -R "~{ref}" -V "$tmp_vcf" \
        --recal-file cohort.indel.recal.vcf.gz \
        --tranches-file cohort.indel.tranches.txt \
        --mode INDEL --truth-sensitivity-filter-level 95.0 \
        -O cohort.vqsr.vcf.gz \
        --tmp-dir "~{tmpdir}" \
        2>> "../$LOG"
      tmp_vcf="cohort.vqsr.vcf.gz"
    fi

    echo "Step 7: Hard-filter and write cohort QC VCF" | tee -a "../$LOG"
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
      -O cohort.gv.QC.vcf.gz \
      --tmp-dir "~{tmpdir}" \
      2>> "../$LOG"

    echo "Step 8: VCF hash" | tee -a "../$LOG"
    bash "~{vcf2hash_script}" cohort.gv.QC.vcf.gz > "../$STATSDIR/cohort.gv.QC.vcf.sha256.txt" 2>> "../$LOG"

    echo "All done! Cohort QC VCF: $VARCALLDIR/cohort.gv.QC.vcf.gz" | tee -a "../$LOG"
  >>>

  output {
    File raw_vcf = "02_varcall/cohort.gv.raw.vcf.gz"
    File qc_vcf = "02_varcall/cohort.gv.QC.vcf.gz"
    File vcf_hash = "03_stats/cohort.gv.QC.vcf.sha256.txt"
    Array[File] logs = glob("logs/*")
  }

  runtime {
    cpu: threads
  }
}
