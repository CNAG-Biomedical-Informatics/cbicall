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
        assert (
            "generated whole-genome intervals" in text
            or "write_wgs_interval_list" in text
            or "writeWgsIntervalList" in text
        ), relpath
        direct_command_args = [
            "  --merge-input-intervals true \\",
            "      --merge-input-intervals true \\",
            "          --merge-input-intervals true \\",
        ]
        for direct_arg in direct_command_args:
            assert direct_arg not in text, relpath
