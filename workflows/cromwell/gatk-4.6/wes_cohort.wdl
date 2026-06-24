version 1.0

workflow CBIcallCohort {
  input {
    String pipeline
    String genome
    Int threads
    File? sample_map
    String workspace
    String output_basename = "cohort"
    String cohort_stage = "all"
    String? interval_shard
    File? input_vcf
    String tmpdir
    String gatk4_cmd
    String ref
    String ref_dict
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
      output_basename = output_basename,
      cohort_stage = cohort_stage,
      interval_shard = interval_shard,
      input_vcf = input_vcf,
      tmpdir = tmpdir,
      gatk4_cmd = gatk4_cmd,
      ref = ref,
      ref_dict = ref_dict,
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
    Array[File] raw_vcf = RunCohort.raw_vcf
    Array[File] qc_vcf = RunCohort.qc_vcf
    Array[File] vcf_hash = RunCohort.vcf_hash
    Array[File] logs = RunCohort.logs
  }
}

task RunCohort {
  input {
    String pipeline
    String genome
    Int threads
    File? sample_map
    String workspace
    String output_basename
    String cohort_stage
    String? interval_shard
    File? input_vcf
    String tmpdir
    String gatk4_cmd
    String ref
    String ref_dict
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
    set -o pipefail
    export TMPDIR="~{tmpdir}"
    export LC_ALL=C
    export GATK_DISABLE_AUTO_S3_UPLOAD=true

    GENOMICSDBDIR="01_genomicsdb"
    VARCALLDIR="02_varcall"
    STATSDIR="03_stats"
    LOGDIR="logs"
    mkdir -p "$GENOMICSDBDIR" "$VARCALLDIR" "$STATSDIR" "$LOGDIR"
    LOG="$LOGDIR/cohort_joint_genotyping.log"

    SAMPLE_MAP="~{sep="" select_all([sample_map])}"
    INTERVAL_SHARD="~{select_first([interval_shard, ""])}"
    INPUT_VCF="~{sep="" select_all([input_vcf])}"
    COHORT_STAGE="~{cohort_stage}"
    OUTPUT_BASENAME="~{output_basename}"

    case "$COHORT_STAGE" in
      all|shard|finalize) ;;
      *) echo "Error: cohort_stage must be all, shard, or finalize." >&2; exit 1;;
    esac
    if [ "$COHORT_STAGE" = "shard" ] && [ -z "$INTERVAL_SHARD" ]; then
      echo "Error: cohort_stage=shard requires interval_shard." >&2
      exit 1
    fi
    if [ "$COHORT_STAGE" = "finalize" ] && [ -z "$INPUT_VCF" ]; then
      echo "Error: cohort_stage=finalize requires input_vcf." >&2
      exit 1
    fi
    if [ "$COHORT_STAGE" != "finalize" ] && [ -z "$SAMPLE_MAP" ]; then
      echo "Error: sample_map is required unless cohort_stage=finalize." >&2
      exit 1
    fi
    ARCH="$(uname -m)"
    if [[ "$COHORT_STAGE" != "finalize" && "$ARCH" =~ ^(aarch64|arm64)$ ]]; then
      echo "Error: GATK GenomicsDBImport cannot run on ARM/aarch64 with the bundled GATK 4.6 GenomicsDB native libraries." >&2
      echo "Run cohort_stage=all/shard on x86_64, or run cohort_stage=finalize on ARM using a gathered raw VCF created on x86_64." >&2
      exit 1
    fi

    COHORT_RAW_VCF="${OUTPUT_BASENAME}.gv.raw.vcf.gz"
    COHORT_QC_VCF="${OUTPUT_BASENAME}.gv.QC.vcf.gz"
    COHORT_VQSR_SNP="${OUTPUT_BASENAME}.snp.recal.vcf.gz"
    COHORT_SNP_TRANCHES="${OUTPUT_BASENAME}.snp.tranches.txt"
    COHORT_VQSR_INDEL="${OUTPUT_BASENAME}.indel.recal.vcf.gz"
    COHORT_INDEL_TRANCHES="${OUTPUT_BASENAME}.indel.tranches.txt"
    COHORT_POST_SNP="${OUTPUT_BASENAME}.post_snp.vcf.gz"
    COHORT_POST_VQSR="${OUTPUT_BASENAME}.vqsr.vcf.gz"
    COHORT_HASH="${OUTPUT_BASENAME}.gv.QC.vcf.sha256.txt"

    SAMPLE_COUNT=0
    GENOTYPE_ANNOTATION_EXCLUDE_ARG=""
    if [ "$COHORT_STAGE" != "finalize" ]; then
      SAMPLE_COUNT=$(awk 'NF {n++} END {print n+0}' "$SAMPLE_MAP")
      if [ "$SAMPLE_COUNT" -lt 10 ]; then
        GENOTYPE_ANNOTATION_EXCLUDE_ARG="-AX InbreedingCoeff"
      fi
    fi
    case "~{workspace}" in
      /*) WORKSPACE_PATH="~{workspace}" ;;
      *) WORKSPACE_PATH="$GENOMICSDBDIR/~{workspace}" ;;
    esac

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

    make_wes_shard_interval_list() {
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

    make_wgs_interval_list() {
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
            if (sn != "" && ln != "" && (shard == "" || sn == shard)) {
              intervals[++n] = sn "\t1\t" ln "\t+\t" sn
            }
          }
          next
        }
        END {
          if (n == 0) {
            printf("No reference contig found for shard %s in %s\n", shard, FILENAME) > "/dev/stderr"
            exit 2
          }
          for (i = 1; i <= n; i++) print intervals[i]
        }
      ' "$ref_dict" > "$out_list"
    }

    INTERVAL_ARG=""
    MERGE_INTERVALS_ARG=""
    if [ "$COHORT_STAGE" = "finalize" ]; then
      echo "Finalize stage: using gathered raw VCF $INPUT_VCF" | tee -a "$LOG"
    elif [ "~{pipeline}" = "wes" ]; then
      if [ -n "$INTERVAL_SHARD" ]; then
        WES_INTERVAL_LIST="$(pwd)/$GENOMICSDBDIR/wes.${INTERVAL_SHARD}.interval_list"
        make_wes_shard_interval_list "~{interval_list}" "$INTERVAL_SHARD" "$WES_INTERVAL_LIST"
        INTERVAL_ARG="-L $WES_INTERVAL_LIST"
        echo "WES mode: generated shard interval list $WES_INTERVAL_LIST from ~{interval_list}" | tee -a "$LOG"
      else
        INTERVAL_ARG="-L ~{interval_list}"
        echo "WES mode: restricting to ~{interval_list}" | tee -a "$LOG"
      fi
      MERGE_INTERVALS_ARG="--merge-input-intervals true"
    else
      if [ -n "$INTERVAL_SHARD" ]; then
        WGS_INTERVAL_LIST="$(pwd)/$GENOMICSDBDIR/wgs.${INTERVAL_SHARD}.interval_list"
        make_wgs_interval_list "~{ref_dict}" "$INTERVAL_SHARD" "$WGS_INTERVAL_LIST"
        echo "WGS mode: generated shard interval list $WGS_INTERVAL_LIST from ~{ref_dict}" | tee -a "$LOG"
      else
        WGS_INTERVAL_LIST="$(pwd)/$GENOMICSDBDIR/wgs.whole_genome.interval_list"
        make_wgs_interval_list "~{ref_dict}" "" "$WGS_INTERVAL_LIST"
        echo "WGS mode: generated whole-genome intervals from ~{ref_dict}" | tee -a "$LOG"
      fi
      INTERVAL_ARG="-L $WGS_INTERVAL_LIST"
    fi

    {
      echo "## Cohort GenomicsDBImport -> Genotype -> VQSR/Hard-filter"
      echo "cohort_stage: $COHORT_STAGE"
      echo "stage_action: $STAGE_ACTION"
      echo "sample_map: $SAMPLE_MAP_DISPLAY"
      echo "pipeline: ~{pipeline}"
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
      echo "tmpdir: ~{tmpdir}"
      echo "log: $LOG"
      echo ""
    } | tee -a "$LOG"

    if [ "$COHORT_STAGE" != "finalize" ]; then
      echo "Step 1: GenomicsDBImport" | tee -a "$LOG"
      ~{gatk4_cmd} GenomicsDBImport \
        --sample-name-map "$SAMPLE_MAP" \
        --genomicsdb-workspace-path "$WORKSPACE_PATH" \
        $MERGE_INTERVALS_ARG \
        $INTERVAL_ARG \
        --tmp-dir "~{tmpdir}" \
        2>> "$LOG"

      echo ok > "$GENOMICSDBDIR/genomicsdbimport.done"

      cd "$VARCALLDIR"

      echo "Step 2: GenotypeGVCFs" | tee -a "../$LOG"
      ~{gatk4_cmd} GenotypeGVCFs \
        -R "~{ref}" \
        -V "gendb://$WORKSPACE_PATH" \
        -O "$COHORT_RAW_VCF" \
        --stand-call-conf 10 \
        $GENOTYPE_ANNOTATION_EXCLUDE_ARG \
        --tmp-dir "~{tmpdir}" \
        $INTERVAL_ARG \
        2>> "../$LOG"

      if [ "$COHORT_STAGE" = "shard" ]; then
        echo "Shard stage complete. Raw cohort shard VCF: $COHORT_RAW_VCF" | tee -a "../$LOG"
        exit 0
      fi
      RAW_VCF_FOR_FILTERING="$COHORT_RAW_VCF"
    else
      cp "$INPUT_VCF" "$VARCALLDIR/$COHORT_RAW_VCF"
      cd "$VARCALLDIR"
      ~{gatk4_cmd} IndexFeatureFile -I "$COHORT_RAW_VCF" 2>> "../$LOG"
      RAW_VCF_FOR_FILTERING="$COHORT_RAW_VCF"
    fi

    echo "Step 3: VQSR decision" | tee -a "../$LOG"
    nSNP=$(zgrep -v '^#' "$RAW_VCF_FOR_FILTERING" | awk 'length($5)==1' | wc -l | tr -d ' ')
    nINDEL=$(zgrep -v '^#' "$RAW_VCF_FOR_FILTERING" | awk 'length($5)!=1' | wc -l | tr -d ' ')
    apply_snp=false
    apply_indel=false
    tmp_vcf="$RAW_VCF_FOR_FILTERING"
    echo "Found SNPs: $nSNP ; INDELs: $nINDEL" | tee -a "../$LOG"

    if [ "$nSNP" -ge "~{min_snp_for_vqsr}" ]; then
      echo "Step 4: SNP VariantRecalibrator" | tee -a "../$LOG"
      ~{gatk4_cmd} VariantRecalibrator \
        -R "~{ref}" \
        -V "$RAW_VCF_FOR_FILTERING" \
        ~{snp_res} \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
        --mode SNP \
        -O "$COHORT_VQSR_SNP" \
        --tranches-file "$COHORT_SNP_TRANCHES" \
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
        -V "$RAW_VCF_FOR_FILTERING" \
        ~{indel_res} \
        -an QD -an FS -an ReadPosRankSum \
        --mode INDEL \
        -O "$COHORT_VQSR_INDEL" \
        --tranches-file "$COHORT_INDEL_TRANCHES" \
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
        --recal-file "$COHORT_VQSR_SNP" \
        --tranches-file "$COHORT_SNP_TRANCHES" \
        --mode SNP --truth-sensitivity-filter-level 99.0 \
        -O "$COHORT_POST_SNP" \
        --tmp-dir "~{tmpdir}" \
        2>> "../$LOG"
      tmp_vcf="$COHORT_POST_SNP"
    fi

    if [ "$apply_indel" = true ]; then
      ~{gatk4_cmd} ApplyVQSR \
        -R "~{ref}" -V "$tmp_vcf" \
        --recal-file "$COHORT_VQSR_INDEL" \
        --tranches-file "$COHORT_INDEL_TRANCHES" \
        --mode INDEL --truth-sensitivity-filter-level 95.0 \
        -O "$COHORT_POST_VQSR" \
        --tmp-dir "~{tmpdir}" \
        2>> "../$LOG"
      tmp_vcf="$COHORT_POST_VQSR"
    fi

    echo "Step 7: Hard-filter and write cohort QC VCF" | tee -a "../$LOG"
    ~{gatk4_cmd} VariantFiltration \
      -R "~{ref}" \
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
      --tmp-dir "~{tmpdir}" \
      2>> "../$LOG"

    echo "Step 8: VCF hash" | tee -a "../$LOG"
    bash "~{vcf2hash_script}" "$COHORT_QC_VCF" > "../$STATSDIR/$COHORT_HASH" 2>> "../$LOG"

    echo "All done! Cohort QC VCF: $VARCALLDIR/$COHORT_QC_VCF" | tee -a "../$LOG"
  >>>

  output {
    Array[File] raw_vcf = glob("02_varcall/*.gv.raw.vcf.gz")
    Array[File] qc_vcf = glob("02_varcall/*.gv.QC.vcf.gz")
    Array[File] vcf_hash = glob("03_stats/*.sha256.txt")
    Array[File] logs = glob("logs/*")
  }

  runtime {
    cpu: threads
  }
}
