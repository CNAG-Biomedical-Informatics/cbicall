from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_mtoolbox_config_preserves_selected_runtime_profile():
    text = (REPO_ROOT / "mtdna" / "MToolBox_config.sh").read_text(encoding="utf-8")
    assert 'source "${CBICALL_ENV_FILE:-$ENVDIR/env.sh}"' in text
    assert 'source "$ENVDIR/env.sh"' not in text


def test_mtoolbox_config_uses_runtime_profile_java_and_samtools():
    text = (REPO_ROOT / "mtdna" / "MToolBox_config.sh").read_text(encoding="utf-8")
    assert 'export PATH="$(dirname "$JAVA8"):$PATH"' in text
    assert 'samtoolsexe="${MTOOLBOX_SAM:-$NGSUTILS/samtools-1.3/samtools}"' in text
    assert "samtools_version=1.3" in text

    for stack in ("gatk-3.5", "gatk-4.6"):
        profile = (
            REPO_ROOT / "workflows" / "bash" / stack / "cnag-hpc-env.sh"
        ).read_text(encoding="utf-8")
        assert (
            'BWA="${CBICALL_BWA:-$NGSUTILS/bwa-0.7.18-cnaghpc/bwa}"'
            in profile
        )
        assert (
            'SAM="${CBICALL_SAM:-$NGSUTILS/samtools-0.1.19-cnaghpc/samtools}"'
            in profile
        )
        assert (
            'MTOOLBOX_SAM="${CBICALL_MTOOLBOX_SAM:-$NGSUTILS/samtools-1.3-cnaghpc/samtools}"'
            in profile
        )


def test_mtdna_workflows_require_exported_bams():
    for workflow in ("mit_single.sh", "mit_cohort.sh"):
        workflow_text = (
            REPO_ROOT / "workflows" / "bash" / "gatk-3.5" / workflow
        ).read_text(encoding="utf-8")
        assert "04_mtdna_input" in workflow_text
        assert "export_mtdna_bam: true" in workflow_text
        assert "MIT_EXTRACT_SAM" not in workflow_text
        assert "rg.merged.dedup.recal.bam" not in workflow_text
        assert "input.merged.filtered.realigned.fixed.bam" not in workflow_text


def test_native_gatk46_single_workflows_support_mtdna_export():
    workflow_files = [
        "workflows/bash/gatk-4.6/wes_single.sh",
        "workflows/bash/gatk-4.6/wgs_single.sh",
        "workflows/snakemake/gatk-4.6/wes_single.smk",
        "workflows/snakemake/gatk-4.6/wgs_single.smk",
        "workflows/nextflow/gatk-4.6/wes_single.nf",
        "workflows/nextflow/gatk-4.6/wgs_single.nf",
        "workflows/cromwell/gatk-4.6/wes_single.wdl",
        "workflows/cromwell/gatk-4.6/wgs_single.wdl",
    ]
    for relpath in workflow_files:
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert "export_mtdna_bam" in text or "--export-mtdna-bam" in text, relpath
        assert "04_mtdna_input" in text, relpath
        assert "-DNA_MIT.bam" in text, relpath
        assert "view -c" in text, relpath


def test_backend_local_gatk46_config_files_are_identical():
    config_paths = [
        REPO_ROOT / "workflows" / backend / "gatk-4.6" / "config.yaml"
        for backend in ("snakemake", "nextflow", "cromwell")
    ]
    assert config_paths[0].read_bytes() == config_paths[1].read_bytes()
    assert config_paths[0].read_bytes() == config_paths[2].read_bytes()


def test_gatk46_cohort_merge_input_intervals_is_wes_only():
    workflow_files = [
        "workflows/bash/gatk-4.6/wes_cohort.sh",
        "workflows/snakemake/gatk-4.6/wes_cohort.smk",
        "workflows/nextflow/gatk-4.6/wes_cohort.nf",
        "workflows/cromwell/gatk-4.6/wes_cohort.wdl",
    ]
    for relpath in workflow_files:
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert "MERGE_INTERVALS_ARG" in text, relpath
        assert "--merge-input-intervals true" in text, relpath
        assert (
            "MERGE_INTERVALS_ARG=\"\"" in text
            or "MERGE_INTERVALS_ARG = \"\"" in text
            or "MERGE_INTERVALS_ARG = PIPELINE == 'wes'" in text
        ), relpath

        assert "MERGE_INTERVALS_ARG" in text.split("GenomicsDBImport", 1)[1], relpath
        assert "wgs.whole_genome.interval_list" in text, relpath
        assert "01_genomicsdb" in text, relpath
        assert (
            "generated whole-genome intervals" in text
            or "write_wgs_interval_list" in text
            or "writeWgsIntervalList" in text
        ), relpath
        assert "02_varcall/wgs.whole_genome.interval_list" not in text, relpath
        direct_command_args = [
            "  --merge-input-intervals true \\",
            "      --merge-input-intervals true \\",
            "          --merge-input-intervals true \\",
        ]
        for direct_arg in direct_command_args:
            assert direct_arg not in text, relpath


def test_gatk46_rank_sum_filters_are_guarded():
    workflow_files = [
        "workflows/bash/gatk-4.6/wes_single.sh",
        "workflows/bash/gatk-4.6/wes_cohort.sh",
        "workflows/snakemake/gatk-4.6/wes_single.smk",
        "workflows/snakemake/gatk-4.6/wes_cohort.smk",
        "workflows/nextflow/gatk-4.6/wes_single.nf",
        "workflows/nextflow/gatk-4.6/wes_cohort.nf",
        "workflows/cromwell/gatk-4.6/wes_single.wdl",
        "workflows/cromwell/gatk-4.6/wes_cohort.wdl",
    ]
    for relpath in workflow_files:
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" in text, relpath
        assert "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" in text, relpath
        assert "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" in text, relpath
        assert '--filter-expression "MQRankSum < -12.5"' not in text, relpath
        assert '--filter-expression "ReadPosRankSum < -8.0"' not in text, relpath
        assert '--filter-expression "ReadPosRankSum < -20.0"' not in text, relpath


def test_gatk46_qd_filters_are_guarded():
    workflow_files = [
        "workflows/bash/gatk-4.6/wes_single.sh",
        "workflows/bash/gatk-4.6/wes_cohort.sh",
        "workflows/snakemake/gatk-4.6/wes_single.smk",
        "workflows/snakemake/gatk-4.6/wes_cohort.smk",
        "workflows/nextflow/gatk-4.6/wes_single.nf",
        "workflows/nextflow/gatk-4.6/wes_cohort.nf",
        "workflows/cromwell/gatk-4.6/wes_single.wdl",
        "workflows/cromwell/gatk-4.6/wes_cohort.wdl",
    ]
    for relpath in workflow_files:
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert "vc.hasAttribute('QD') && QD < 2.0" in text, relpath
        assert '--filter-expression "QD < 2.0"' not in text, relpath


def test_gatk46_single_bwa_stderr_is_logged_before_pipe():
    expected_by_file = {
        "workflows/bash/gatk-4.6/wes_single.sh": (
            '$BWA mem -M -K 40000000 -t "$THREADS" "$REFGZ" "$R1" "$R2" 2>> "$LOG" \\'
        ),
        "workflows/snakemake/gatk-4.6/wes_single.smk": (
            "{BWA} mem -M -K 40000000 -t {threads} {REFGZ} {input.r1} {input.r2} 2>> {log} \\"
        ),
        "workflows/nextflow/gatk-4.6/wes_single.nf": (
            '${BWA} mem -M -K 40000000 -t ${task.cpus} ${q(REFGZ)} ${q(r1)} ${q(r2)} '
            '2>> ${q("${ID}.01_align_rg.${base}.log")} \\\\'
        ),
        "workflows/cromwell/gatk-4.6/wes_single.wdl": (
            '~{bwa} mem -M -K 40000000 -t ~{threads} "~{refgz}" "$R1" "$R2" 2>> "$LOG" \\'
        ),
    }
    for relpath, expected in expected_by_file.items():
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert expected in text, relpath


def test_native_bwa_mem_uses_fixed_reproducible_batch_size():
    workflow_files = [
        "workflows/bash/gatk-3.5/wes_single.sh",
        "workflows/bash/gatk-4.6/wes_single.sh",
        "workflows/bash/gatk-4.6/wgs_single.sh",
        "workflows/snakemake/gatk-4.6/wes_single.smk",
        "workflows/snakemake/gatk-4.6/wgs_single.smk",
        "workflows/nextflow/gatk-4.6/wes_single.nf",
        "workflows/nextflow/gatk-4.6/wgs_single.nf",
        "workflows/cromwell/gatk-4.6/wes_single.wdl",
        "workflows/cromwell/gatk-4.6/wgs_single.wdl",
    ]
    for relpath in workflow_files:
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert "mem -M -K 40000000 -t" in text, relpath


def test_native_single_alignment_pipelines_enable_pipefail():
    workflow_files = [
        "workflows/bash/gatk-3.5/wes_single.sh",
        "workflows/bash/gatk-4.6/wes_single.sh",
        "workflows/snakemake/gatk-4.6/wes_single.smk",
        "workflows/nextflow/gatk-4.6/wes_single.nf",
        "workflows/cromwell/gatk-4.6/wes_single.wdl",
    ]
    for relpath in workflow_files:
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        alignment = text.split("mem -M", 1)[0]
        assert "set -o pipefail" in alignment, relpath


def test_nextflow_single_uses_canonical_sample_id_for_read_groups():
    relpath = "workflows/nextflow/gatk-4.6/wes_single.nf"
    text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
    alignment = text.split("process ALIGN_RG", 1)[1].split("process MERGE_BAMS", 1)[0]

    assert 'RGID="${ID}.${base}"' in alignment
    assert 'RGPU="${ID}.${base}.unit1"' in alignment
    assert '--RGSM ${q(ID)}' in alignment
    assert "cut -d'_'" not in alignment
    assert '--RGSM "\\$SAMPLE"' not in alignment


def test_gatk46_small_cohort_excludes_inbreeding_coeff():
    workflow_files = [
        "workflows/bash/gatk-4.6/wes_cohort.sh",
        "workflows/snakemake/gatk-4.6/wes_cohort.smk",
        "workflows/nextflow/gatk-4.6/wes_cohort.nf",
        "workflows/cromwell/gatk-4.6/wes_cohort.wdl",
    ]
    for relpath in workflow_files:
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert "SAMPLE_COUNT" in text, relpath
        assert "GENOTYPE_ANNOTATION_EXCLUDE_ARG" in text, relpath
        assert "-AX InbreedingCoeff" in text, relpath
        assert "SAMPLE_COUNT < 10" in text or '"$SAMPLE_COUNT" -lt 10' in text, relpath
        genotype_command = text.split("GenotypeGVCFs", 1)[1]
        assert "GENOTYPE_ANNOTATION_EXCLUDE_ARG" in genotype_command, relpath


def test_gatk46_cohort_finalize_header_is_human_readable_across_backends():
    expected_fragments = {
        "workflows/bash/gatk-4.6/wes_cohort.sh": [
            'STAGE_ACTION="finalize run: globally filter a gathered raw cohort VCF"',
            'echo "stage_action: $STAGE_ACTION"',
            'SAMPLE_COUNT_DISPLAY="not applicable (finalize stage)"',
            'SAMPLE_MAP_DISPLAY="not used (finalize stage)"',
            'WORKSPACE_DISPLAY="not used (finalize stage)"',
            'echo "input_vcf: $INPUT_VCF"',
            'echo "final_vcf: $COHORT_QC_VCF"',
        ],
        "workflows/snakemake/gatk-4.6/wes_cohort.smk": [
            '"finalize": "finalize run: globally filter a gathered raw cohort VCF"',
            'log.write(f"stage_action: {STAGE_ACTION}\\n")',
            'SAMPLE_COUNT_DISPLAY = "not applicable (finalize stage)"',
            'SAMPLE_MAP_DISPLAY = "not used (finalize stage)"',
            'WORKSPACE_DISPLAY = "not used (finalize stage)"',
            'log.write(f"input_vcf: {INPUT_VCF}\\n")',
            'log.write(f"final_vcf: {COHORT_QC_VCF}\\n")',
        ],
        "workflows/nextflow/gatk-4.6/wes_cohort.nf": [
            "'finalize': 'finalize run: globally filter a gathered raw cohort VCF'",
            'println "Stage action: ${STAGE_ACTION}"',
            'not applicable (finalize stage)',
            'not used (finalize stage)',
            'println "Input VCF: ${INPUT_VCF}"',
            'println "Final VCF: ${COHORT_QC_VCF_NAME}"',
        ],
        "workflows/cromwell/gatk-4.6/wes_cohort.wdl": [
            'STAGE_ACTION="finalize run: globally filter a gathered raw cohort VCF"',
            'echo "stage_action: $STAGE_ACTION"',
            'SAMPLE_COUNT_DISPLAY="not applicable (finalize stage)"',
            'SAMPLE_MAP_DISPLAY="not used (finalize stage)"',
            'WORKSPACE_DISPLAY="not used (finalize stage)"',
            'echo "input_vcf: $INPUT_VCF"',
            'echo "final_vcf: $COHORT_QC_VCF"',
        ],
    }
    for relpath, fragments in expected_fragments.items():
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        for fragment in fragments:
            assert fragment in text, f"{relpath}: {fragment}"


def test_gatk46_cohort_supports_staged_shard_finalize_controls():
    workflow_files = [
        "workflows/bash/gatk-4.6/wes_cohort.sh",
        "workflows/snakemake/gatk-4.6/wes_cohort.smk",
        "workflows/nextflow/gatk-4.6/wes_cohort.nf",
        "workflows/cromwell/gatk-4.6/wes_cohort.wdl",
    ]
    for relpath in workflow_files:
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert "cohort_stage" in text or "COHORT_STAGE" in text, relpath
        assert "output_basename" in text or "OUTPUT_BASENAME" in text, relpath
        assert "interval_shard" in text or "INTERVAL_SHARD" in text, relpath
        assert "input_vcf" in text or "INPUT_VCF" in text, relpath
        assert "shard" in text, relpath
        assert "finalize" in text, relpath
        assert "GenomicsDBImport cannot run on ARM/aarch64" in text, relpath
        assert "RAW_VCF_FOR_FILTERING" in text or "rawvcf" in text or "RAW_VCF_FOR_FILTERING" in text, relpath
        assert ".gv.raw.vcf.gz" in text, relpath
