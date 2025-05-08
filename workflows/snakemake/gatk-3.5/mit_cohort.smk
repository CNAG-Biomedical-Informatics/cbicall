#############################################
# Snakefile for Cohort mtDNA Pipeline
#
# This workflow performs the following steps:
# 1. For each sample in the cohort, extract the mtDNA region
#    from its exome BAM.
# 2. Change into the MTOOLBOX directory, run MToolBox (cohort-level)
#    for variant calling/annotation.
# 3. Append heteroplasmy information (GT, DP, HF) to the prioritized
#    variants file exactly as in the bash script.
#
# $VERSION taken from CBICall
#############################################

import os, platform, glob
from pathlib import Path

# Set snakefile directory and load global config
snakefile_dir = Path(workflow.snakefile).parent
configfile: snakefile_dir / "config.yaml"

# Define and create LOGDIR relative to the snakefile directory
LOGDIR = os.path.abspath("logs")
os.makedirs(LOGDIR, exist_ok=True)

# Hardcoded mtDNA-specific tools/paths for this pipeline
# (These correspond to the cohort script’s expectations.)
PARSE_VAR   = os.path.join(snakefile_dir, "parse_var.pl")
PARSE_PRIOR = os.path.join(snakefile_dir, "parse_prioritized.pl")
MTDNA_DIR   = os.path.join(snakefile_dir, "../../mtdna")  # Directory containing MToolBox_config.sh

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
# MToolBox installation: use global config for mtoolbox (with substitutions),
# but override with the hardcoded MTDNA_DIR for this cohort.
# ----------------------------------------------------------------------
MTOOLBOXDIR = config["mtoolboxdir"].format(ngsutils=config["ngsutils"])
BINDIRMTB   = MTDNA_DIR

# ----------------------------------------------------------------------
# Determine the mtDNA chromosome label based on the reference.
# If the reference filename contains "b37", use "MT"; otherwise "chrM".
# ----------------------------------------------------------------------
REF = config["ref"]
if "b37" in REF:
    CHRM = "MT"
else:
    CHRM = "chrM"

# ----------------------------------------------------------------------
# Discover sample directories.
# The cohort directory is expected to contain several sample folders matching:
#   ../??????????_ex/cbicall_snakemake_wes_single_*/01_bam
# For each such directory, the sample ID is computed as:
#    (first part of parent folder name) + "-DNA_MIT"
# For example, ../MA0002401P_ex/cbicall_snakemake_wes_single_146723488708442/01_bam
# gives sample id "MA0002401P-DNA_MIT"
# ----------------------------------------------------------------------
# Build a dictionary with sample IDs and their corresponding exome BAM paths.
sample_dirs = glob.glob(os.path.join("..", "??????????_ex", "cbicall_snakemake_wes_single_*", "01_bam"))
if not sample_dirs:
    raise Exception("No sample directories found matching the cohort pattern.")

SAMPLES = {}
for sd in sample_dirs:
    # Instead of using os.path.dirname(sd) which gives the 'cbicall_snakemake_...' folder,
    # use os.path.dirname(os.path.dirname(sd)) to get the sample folder (e.g., "MA0002401P_ex")
    parent = os.path.basename(os.path.dirname(os.path.dirname(sd)))
    sample = parent.split('_')[0] + "-DNA_MIT"
    bam_path = os.path.join(sd, "input.merged.filtered.realigned.fixed.bam")
    bai_path = os.path.join(sd, "input.merged.filtered.realigned.fixed.bai")
    SAMPLES[sample] = {"bam": bam_path, "bai": bai_path}

sample_list = list(SAMPLES.keys())

# Helper functions to retrieve sample-specific inputs.
def get_bam(wildcards):
    return SAMPLES[wildcards.sample]["bam"]

def get_bai(wildcards):
    return SAMPLES[wildcards.sample]["bai"]

#############################################
# Global target: the final annotated cohort mtDNA variants file.
#############################################
rule all:
    input:
        "MTOOLBOX/mit_prioritized_variants.txt"

#############################################
# Rule: extract_mtDNA_cohort
# Description: For each sample, extract the mitochondrial region
# from its exome BAM and index the result.
#############################################
rule extract_mtDNA_cohort:
    input:
        bam = get_bam,
        bai = get_bai
    output:
        mtbam = os.path.join("MTOOLBOX", "{sample}.bam"),
        mtbam_index = os.path.join("MTOOLBOX", "{sample}.bam.bai")
    params:
        chrM = CHRM
    threads: THREADS
    log:
        os.path.join(LOGDIR, "{sample}.extract_mtDNA.log")
    shell:
        r"""
        mkdir -p MTOOLBOX
        # Optionally ensure that the source BAM index is correct:
        if [ ! -s {input.bai}.bam ]; then cp {input.bai} {input.bam}.bam.bai; fi
        {SAMTOOLS} view -b {input.bam} {params.chrM} > {output.mtbam}
        {SAMTOOLS} index {output.mtbam}
        """

#############################################
# Rule: run_MToolBox_cohort
# Description: Change into the MTOOLBOX directory, run MToolBox, and then
# perform the variant postprocessing exactly as in the bash script.
#############################################
rule run_MToolBox_cohort:
    input:
        mtbams = expand(os.path.join("MTOOLBOX", "{sample}.bam"), sample=sample_list)
    output:
        final = os.path.join("MTOOLBOX", "mit_prioritized_variants.txt")
    params:
        BINDIRMTB = BINDIRMTB
    threads: THREADS
    log:
        os.path.join(LOGDIR, "run_MToolBox_cohort.log")
    shell:
        r"""
        pushd MTOOLBOX
        cp {params.BINDIRMTB}/MToolBox_config.sh .
        export PATH={MTOOLBOXDIR}:$PATH
        export PYTHONPATH=~/.local/lib/python2.7/site-packages:${{PYTHONPATH:-}}
        MToolBox.sh -i MToolBox_config.sh -m "-t {threads}"
        # Now perform the postprocessing as in the bash:
        vcf_file="VCF_file.vcf"
        vcf_tmp="VCF_file_tmp.vcf"
        in_file="prioritized_variants.txt"
        out_file="append_tmp.txt"
        grep ^#CHROM $vcf_file > $vcf_tmp
        for var in $(cut -f1 $in_file | sed '1d' | {PARSE_VAR}); do
            grep -P "chrMT\t$var\t" $vcf_file >> $vcf_tmp || echo "$var not found"
        done
        {PARSE_PRIOR} -i $vcf_tmp > $out_file
        paste $in_file $out_file > mit_prioritized_variants.txt
        rm $vcf_tmp $out_file
        popd
        """

