#############################################
# Snakefile for Single-Sample mtDNA Pipeline
#
# This workflow performs the following steps:
# 1. Extract the mitochondrial region from an exome BAM.
# 2. Run MToolBox for mtDNA variant calling and annotation.
# 3. Append heteroplasmy information (GT, DP, HF) to the prioritized
#    variants file.
#
# $VERSION taken from CBICall
#############################################

import os, platform, glob
from pathlib import Path

# Set snakefile directory and load global config
snakefile_dir = Path(workflow.snakefile).parent
configfile: snakefile_dir / "config.yaml"

# Define LOGDIR relative to the snakefile directory
LOGDIR = os.path.abspath("logs")
os.makedirs(LOGDIR, exist_ok=True)

# Hardcoded mtDNA-specific tools/paths for this pipeline
PARSE_VAR = os.path.join(snakefile_dir, "parse_var.pl")
MTDNA_DIR  = os.path.join(snakefile_dir, "../../mtdna")  # Directory containing MToolBox_config.sh

# Global variables from config
DATADIR  = config["datadir"]
DBDIR    = config["dbdir"].format(datadir=DATADIR)
NGSUTILS = config["ngsutils"].format(datadir=DATADIR)
TMPDIR   = config["tmpdir"].format(datadir=DATADIR)
MEM      = config["mem"]

# THREADS: number of cores from the CLI argument --cores (default to 4)
THREADS = workflow.cores if workflow.cores else 4

# Determine architecture and set samtools path based on global config.
ARCH = platform.machine()
if ARCH == "aarch64":
    SAMTOOLS = config["tools"]["aarch64"]["samtools"].format(ngsutils=config["ngsutils"])
else:
    SAMTOOLS = config["tools"]["amd64"]["samtools"].format(ngsutils=config["ngsutils"])

# ----------------------------------------------------------------------
# Hardcoded mtDNA-specific parameters (only for the mtDNA pipeline):
# ----------------------------------------------------------------------
EXOME_BAM_PATTERN = "../cbicall_snakemake_wes_single_*/01_bam/input.merged.filtered.realigned.fixed.bam"
EXOME_BAM_candidates = glob.glob(EXOME_BAM_PATTERN)
if not EXOME_BAM_candidates:
    raise Exception(f"No exome BAM found matching pattern: {EXOME_BAM_PATTERN}")
EXOME_BAM = EXOME_BAM_candidates[0]

EXOME_BAM_INDEX_PATTERN = "../cbicall_snakemake_wes_single_*/01_bam/input.merged.filtered.realigned.fixed.bai"
EXOME_BAM_INDEX_candidates = glob.glob(EXOME_BAM_INDEX_PATTERN)
if not EXOME_BAM_INDEX_candidates:
    raise Exception(f"No exome BAM index found matching pattern: {EXOME_BAM_INDEX_PATTERN}")
EXOME_BAM_INDEX = EXOME_BAM_INDEX_candidates[0]

# ----------------------------------------------------------------------
# Use global config for MToolBox installation (with substitutions)
# ----------------------------------------------------------------------
MTOOLBOXDIR = config["mtoolboxdir"].format(ngsutils=config["ngsutils"])
# For mtDNA, override with the hardcoded MTDNA_DIR:
BINDIRMTB = MTDNA_DIR

# ----------------------------------------------------------------------
# Determine the mtDNA chromosome label based on the reference file.
# If the reference filename contains "b37", use "MT"; otherwise "chrM".
# ----------------------------------------------------------------------
REF = config["ref"]
if "b37" in REF:
    CHRM = "MT"
else:
    CHRM = "chrM"

# ----------------------------------------------------------------------
# Determine sample ID from the parent directory (Bash equivalent):
#   id=$( echo $DIR | awk -F'/' '{print $(NF-1)}' | awk -F'_' '{print $1}' | sed 's/$/-DNA_MIT/' )
# For example, if the parent directory is "MA0004701P_ex",
# then SAMPLE_ID becomes "MA0004701P-DNA_MIT"
# ----------------------------------------------------------------------
cwd = os.getcwd()
parent_dir = os.path.basename(os.path.dirname(cwd))
SAMPLE_ID = parent_dir.split('_')[0] + "-DNA_MIT"

#############################################
# Global target: the final annotated mtDNA variants file.
#############################################
rule all:
    input:
        "MTOOLBOX/mit_prioritized_variants.txt"

#############################################
# Rule: extract_mtDNA
# Description: Extract the mitochondrial region from the exome BAM.
#############################################
rule extract_mtDNA:
    input:
        bam = EXOME_BAM,
        bai = EXOME_BAM_INDEX
    output:
        mtbam       = os.path.join("MTOOLBOX", SAMPLE_ID + ".bam"),
        mtbam_index = os.path.join("MTOOLBOX", SAMPLE_ID + ".bam.bai")
    params:
        chrM = CHRM
    threads: THREADS
    log:
        os.path.join(LOGDIR, SAMPLE_ID + ".extract_mtDNA.log")
    shell:
        r"""
        mkdir -p MTOOLBOX
        {SAMTOOLS} view -b {input.bam} {params.chrM} > {output.mtbam}
        {SAMTOOLS} index {output.mtbam}
        """

#############################################
# Rule: run_MToolBox
# Description: Run MToolBox on the extracted mtDNA BAM.
#
# This rule mimics the Bash behavior by changing into the MTOOLBOX
# directory before executing MToolBox.sh. Inside that directory,
# file names are used without the MTOOLBOX prefix.
#############################################
rule run_MToolBox:
    input:
        mtbam = rules.extract_mtDNA.output.mtbam
    output:
        vcf         = os.path.join("MTOOLBOX", "VCF_file.vcf"),
        prioritized = os.path.join("MTOOLBOX", "prioritized_variants.txt")
    params:
        BINDIRMTB = BINDIRMTB
    threads: THREADS
    log:
        os.path.join(LOGDIR, SAMPLE_ID + ".run_MToolBox.log")
    shell:
        r"""
        pushd MTOOLBOX
        cp {params.BINDIRMTB}/MToolBox_config.sh .
        export PATH={MTOOLBOXDIR}:$PATH
        export PYTHONPATH=~/.local/lib/python2.7/site-packages:${{PYTHONPATH:-}}
        MToolBox.sh -i MToolBox_config.sh -m "-t {threads}" 2>> {log}
        popd
        """

#############################################
# Rule: append_heteroplasmy
# Description: Append heteroplasmy (HF) information from the VCF file
#              to the prioritized variants file.
#
# This rule also changes into the MTOOLBOX directory so that file names
# are used without the MTOOLBOX prefix.
#############################################
rule append_heteroplasmy:
    input:
        vcf         = os.path.join("MTOOLBOX", "VCF_file.vcf"),
        prioritized = os.path.join("MTOOLBOX", "prioritized_variants.txt")
    output:
        final = os.path.join("MTOOLBOX", "mit_prioritized_variants.txt")
    log:
        os.path.join(LOGDIR, SAMPLE_ID + ".append_heteroplasmy.log")
    shell:
        r"""
        pushd MTOOLBOX
        temp_file="append.txt"
        echo -e "REF\tALT\tGT\tDP\tHF" > $temp_file
        for var in $(cut -f1 prioritized_variants.txt | sed '1d' | {PARSE_VAR}); do
            grep -P "chrMT\t$var\t" VCF_file.vcf | cut -f4,5,10 | tr ':' '\t' | cut -f1-5 >> $temp_file || true
        done
        paste prioritized_variants.txt $temp_file > mit_prioritized_variants.txt
        rm $temp_file
        popd
        """
