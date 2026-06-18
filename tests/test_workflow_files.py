from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


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
            '$BWA mem -M -t "$THREADS" "$REFGZ" "$R1" "$R2" 2>> "$LOG" \\'
        ),
        "workflows/snakemake/gatk-4.6/wes_single.smk": (
            "{BWA} mem -M -t {threads} {REFGZ} {input.r1} {input.r2} 2>> {log} \\"
        ),
        "workflows/nextflow/gatk-4.6/wes_single.nf": (
            '${BWA} mem -M -t ${task.cpus} ${q(REFGZ)} ${q(r1)} ${q(r2)} '
            '2>> ${q("${ID}.01_align_rg.${base}.log")} \\\\'
        ),
        "workflows/cromwell/gatk-4.6/wes_single.wdl": (
            '~{bwa} mem -M -t ~{threads} "~{refgz}" "$R1" "$R2" 2>> "$LOG" \\'
        ),
    }
    for relpath, expected in expected_by_file.items():
        text = (REPO_ROOT / relpath).read_text(encoding="utf-8")
        assert expected in text, relpath


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
        assert "RAW_VCF_FOR_FILTERING" in text or "rawvcf" in text or "RAW_VCF_FOR_FILTERING" in text, relpath
        assert ".gv.raw.vcf.gz" in text, relpath
