#############################################
# Snakefile for WES Cohort Pipeline
#
# This workflow uses config.yaml (and similar variables as
# the single-sample pipeline) and converts the wes_cohort.sh
# Bash script into a reproducible Snakemake workflow.
#
# $VERSION taken from CBICall
#############################################

import os, glob, platform
from pathlib import Path

snakefile_dir = Path(workflow.snakefile).parent
configfile: snakefile_dir / "config.yaml"

# Directories for output (relative to the cohort directory)
VARCALLDIR = "01_varcall"
STATSDIR   = "02_stats"
LOGDIR     = "logs"

# Create directories if they do not exist
for d in [VARCALLDIR, STATSDIR, LOGDIR]:
    os.makedirs(d, exist_ok=True)

# Define the cohort ID by extracting it from the parent directory name.
# (For example, if the parent directory is "MA00052_exome", then cohort = "MA00052")
cohort = os.path.basename(os.path.abspath("..")).split("_")[0]

# Define chromosomes as in the Bash script: chromosomes 1..12 plus the extra groups.
CHROMOSOMES = list(map(str, range(1,13))) + ["1314", "1516", "1718", "1920", "2122XY"]

# Global parameters from config (these should be defined in your config.yaml)
DATADIR  = config["datadir"]
DBDIR    = config["dbdir"].format(datadir=DATADIR)
NGSUTILS = config["ngsutils"].format(datadir=DATADIR)
TMPDIR   = config["tmpdir"].format(datadir=DATADIR)
MEM      = config["mem"]

# Get the number of cores from the CLI argument --cores
THREADS = workflow.cores if workflow.cores else 4  # default to 4 if not provided

# Determine JAVA and tool paths based on architecture (similar to single-sample)
ARCH = platform.machine()
if ARCH == "aarch64":
    JAVA = config["java"]["aarch64"]
    BWA  = config["tools"]["aarch64"]["bwa"].format(ngsutils=NGSUTILS)
    SAM  = config["tools"]["aarch64"]["samtools"].format(ngsutils=NGSUTILS)
    BED  = config["tools"]["aarch64"]["bedtools"].format(ngsutils=NGSUTILS)
else:
    JAVA = config["java"]["amd64"]
    BWA  = config["tools"]["amd64"]["bwa"].format(ngsutils=NGSUTILS)
    SAM  = config["tools"]["amd64"]["samtools"].format(ngsutils=NGSUTILS)
    BED  = config["tools"]["amd64"]["bedtools"].format(ngsutils=NGSUTILS)

PIC  = config["picard"].format(java=JAVA, mem=MEM, tmpdir=TMPDIR, ngsutils=NGSUTILS)
GATK = config["gatk"].format(java=JAVA, mem=MEM, tmpdir=TMPDIR, ngsutils=NGSUTILS)

# Other parameters from config
bundle       = config["bundle"].format(dbdir=DBDIR)
REF          = config["ref"].format(bundle=bundle)
REFGZ        = config["refgz"].format(bundle=bundle)
dbSNP        = config["dbsnp"].format(dbdir=DBDIR)
MILLS_INDELS = config["mills_indels"].format(bundle=bundle)
KG_INDELS    = config["kg_indels"].format(bundle=bundle)
HAPMAP       = config["hapmap"].format(bundle=bundle)
OMNI         = config["omni"].format(bundle=bundle)

snp_res   = config["snp_res"].format(hapmap=HAPMAP, omni=OMNI, dbsnp=dbSNP)
indel_res = config["indel_res"].format(mills_indels=MILLS_INDELS)
EXOM      = config["exom"].format(dbdir=DBDIR)

DCOV    = config["dcov"]
UG_CALL = config["ug_call"]
UG_EMIT = config["ug_emit"]

# Path to the jaccard script (assumed to be in the same directory as the Snakefile)
JACCARD = os.path.join(snakefile_dir, "jaccard.sh")

#############################################
# Helper function to get cohort BAMs per chromosome
# The pattern here reflects the structure from the Bash script:
# ../??????????_ex/cbicall_snakemake_wes_single*/01_bam/input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.bam
#############################################
def get_bams_for_chr(wildcards):
    pattern = os.path.join("..", "??????????_ex", "cbicall_snakemake_wes_single*", "01_bam", "input.merged.filtered.realigned.fixed.dedup.chr{0}.recal.bam".format(wildcards.chr))
    bams = sorted(glob.glob(pattern))
    if len(bams) == 0:
        raise ValueError("No BAM files found for chromosome {} using pattern: {}".format(wildcards.chr, pattern))
    return bams

#############################################
# Rule: all
# Global target to build the final QC VCF and jaccard outputs.
#############################################
rule all:
    input:
        os.path.join(VARCALLDIR, "{cohort}.ug.QC.vcf".format(cohort=cohort)),
        os.path.join(STATSDIR, "jaccard.txt"),
        os.path.join(STATSDIR, "jaccard_lt50.txt")

#############################################
# Rule: variant_calling_cohort
# Description: For each chromosome, call variants using GATK UnifiedGenotyper
# with multiple input BAMs from the cohort.
#############################################
rule variant_calling_cohort:
    input:
        bams = get_bams_for_chr
    output:
        vcf = os.path.join(VARCALLDIR, "chr{chr}.ug.raw.vcf")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "variant_calling_cohort.{chr}.log")
    shell:
        r"""
        # Build the -I input string from all BAM files for this chromosome
        in_str=$(for bam in {input.bams}; do echo -n " -I $bam"; done)
        echo "Running GATK UnifiedGenotyper for chromosome {wildcards.chr}" >> {log}
        {GATK} -T UnifiedGenotyper \
               -R {REF} \
               $in_str \
               -L {EXOM}/hg19.chr{wildcards.chr}.flank100bp.bed \
               --dbsnp {dbSNP} \
               -o {output.vcf} \
               -dcov {DCOV} \
               -stand_call_conf {UG_CALL} \
               -stand_emit_conf {UG_EMIT} \
               -nt {threads} \
               -glm BOTH 2>> {log}
        echo "Done with chromosome {wildcards.chr}" >> {log}
        """

#############################################
# Rule: merge_vcfs_cohort
# Description: Merge per-chromosome VCFs into a single merged VCF for the cohort.
#############################################
rule merge_vcfs_cohort:
    input:
        vcfs = expand(os.path.join(VARCALLDIR, "chr{chr}.ug.raw.vcf"), chr=CHROMOSOMES)
    output:
        merged = os.path.join(VARCALLDIR, "{cohort}.ug.raw.vcf".format(cohort=cohort))
    params:
        cohort = cohort
    log:
        os.path.join(LOGDIR, "merge_vcfs_cohort.{cohort}.log".format(cohort=cohort))
    shell:
        r"""
        echo "Merging VCF files for cohort {params.cohort}" >> {log}
        grep "#" {input.vcfs[0]} > {VARCALLDIR}/header.txt
        grep -hv '#' {VARCALLDIR}/chr*.ug.raw.vcf | awk '$1 !~ /_/' | sort -V | cat {VARCALLDIR}/header.txt - > {output.merged} 2>> {log}
        """

#############################################
# Rule: recalibrate_snp_cohort
# Description: Recalibrate SNP variants using GATK VariantRecalibrator,
#              or skip if too few SNPs are found.
#############################################
rule recalibrate_snp_cohort:
    input:
        vcf = os.path.join(VARCALLDIR, "{cohort}.ug.raw.vcf".format(cohort=cohort))
    output:
        snp_recal    = os.path.join(VARCALLDIR, "{cohort}.ug.raw.snp.recal".format(cohort=cohort)),
        snp_tranches = os.path.join(VARCALLDIR, "{cohort}.ug.raw.snp.tranches".format(cohort=cohort))
    log:
        os.path.join(LOGDIR, "recalibrate_snp_cohort.{cohort}.log".format(cohort=cohort))
    shell:
        r"""
        echo "-GATK Recalibrator SNP" >> {log}
        nSNP=$(grep -v "^#" {input.vcf} | awk 'length($5)==1' | wc -l)
        echo "Found $nSNP SNP variants" >> {log}
        if [ $nSNP -gt 1000 ]; then
            {GATK} -T VariantRecalibrator \
                   -R {REF} \
                   -input {input.vcf} \
                   -recalFile {output.snp_recal} \
                   -tranchesFile {output.snp_tranches} \
                   --maxGaussians 6 \
                   {snp_res} \
                   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
                   -mode SNP 2>> {log}
        else
            echo "WARNING: Only $nSNP SNP variants found. Skipping VariantRecalibrator for SNPs." >> {log}
            cp {input.vcf} {output.snp_recal}
            echo "SKIP" > {output.snp_tranches}
        fi
        """

#############################################
# Rule: recalibrate_indel_cohort
# Description: Recalibrate INDEL variants using GATK VariantRecalibrator,
#              or skip if too few INDELs are found.
#############################################
rule recalibrate_indel_cohort:
    input:
        vcf = os.path.join(VARCALLDIR, "{cohort}.ug.raw.vcf".format(cohort=cohort))
    output:
        indel_recal    = os.path.join(VARCALLDIR, "{cohort}.ug.raw.indel.recal".format(cohort=cohort)),
        indel_tranches = os.path.join(VARCALLDIR, "{cohort}.ug.raw.indel.tranches".format(cohort=cohort))
    log:
        os.path.join(LOGDIR, "recalibrate_indel_cohort.{cohort}.log".format(cohort=cohort))
    shell:
        r"""
        echo "-GATK Recalibrator INDEL" >> {log}
        nINDEL=$(grep -v "^#" {input.vcf} | awk 'length($5)!=1' | wc -l)
        echo "Found $nINDEL INDEL variants" >> {log}
        if [ $nINDEL -gt 8000 ]; then
            {GATK} -T VariantRecalibrator \
                   -R {REF} \
                   -input {input.vcf} \
                   -recalFile {output.indel_recal} \
                   -tranchesFile {output.indel_tranches} \
                   --maxGaussians 4 \
                   {indel_res} \
                   -an QD -an FS -an ReadPosRankSum \
                   -mode INDEL 2>> {log}
        else
            echo "WARNING: Only $nINDEL INDEL variants found. Skipping VariantRecalibrator for INDELs." >> {log}
            cp {input.vcf} {output.indel_recal}
            echo "SKIP" > {output.indel_tranches}
        fi
        """

#############################################
# Rule: apply_recalibration_snp_cohort
# Description: Apply SNP recalibration or bypass it if skipped.
#############################################
rule apply_recalibration_snp_cohort:
    input:
        vcf          = os.path.join(VARCALLDIR, "{cohort}.ug.raw.vcf".format(cohort=cohort)),
        snp_recal    = os.path.join(VARCALLDIR, "{cohort}.ug.raw.snp.recal".format(cohort=cohort)),
        snp_tranches = os.path.join(VARCALLDIR, "{cohort}.ug.raw.snp.tranches".format(cohort=cohort))
    output:
        recal_snp = os.path.join(VARCALLDIR, "recalibratedSNPs.rawIndels.vcf")
    log:
        os.path.join(LOGDIR, "apply_recalibration_snp_cohort.{cohort}.log".format(cohort=cohort))
    shell:
        r"""
        echo "-GATK Apply Recalibrator SNP" >> {log}
        if grep -q "SKIP" {input.snp_tranches}; then
            echo "Skipping ApplyRecalibration for SNPs due to low variant count" >> {log}
            cp {input.vcf} {output.recal_snp}
        else
            {GATK} -T ApplyRecalibration \
                   -R {REF} \
                   -input {input.vcf} \
                   -recalFile {input.snp_recal} \
                   -tranchesFile {input.snp_tranches} \
                   -o {output.recal_snp} \
                   --ts_filter_level 99.0 \
                   -mode SNP 2>> {log}
        fi
        """

#############################################
# Rule: apply_recalibration_indel_cohort
# Description: Apply INDEL recalibration or bypass it if skipped.
#############################################
rule apply_recalibration_indel_cohort:
    input:
        recal_snp     = os.path.join(VARCALLDIR, "recalibratedSNPs.rawIndels.vcf"),
        indel_recal   = os.path.join(VARCALLDIR, "{cohort}.ug.raw.indel.recal".format(cohort=cohort)),
        indel_tranches= os.path.join(VARCALLDIR, "{cohort}.ug.raw.indel.tranches".format(cohort=cohort))
    output:
        vqsr_vcf = os.path.join(VARCALLDIR, "{cohort}.ug.vqsr.vcf".format(cohort=cohort))
    log:
        os.path.join(LOGDIR, "apply_recalibration_indel_cohort.{cohort}.log".format(cohort=cohort))
    shell:
        r"""
        echo "-GATK Apply Recalibrator INDEL" >> {log}
        if grep -q "SKIP" {input.indel_tranches}; then
            echo "Skipping ApplyRecalibration for INDELs due to low variant count" >> {log}
            cp {input.recal_snp} {output.vqsr_vcf}
        else
            {GATK} -T ApplyRecalibration \
                   -R {REF} \
                   -input {input.recal_snp} \
                   -recalFile {input.indel_recal} \
                   -tranchesFile {input.indel_tranches} \
                   -o {output.vqsr_vcf} \
                   --ts_filter_level 95.0 \
                   -mode INDEL 2>> {log}
        fi
        """

#############################################
# Rule: variant_filtration_cohort
# Description: Apply hard filters to produce the final QC VCF.
#############################################
rule variant_filtration_cohort:
    input:
        vqsr_vcf = os.path.join(VARCALLDIR, "{cohort}.ug.vqsr.vcf".format(cohort=cohort))
    output:
        qc_vcf = os.path.join(VARCALLDIR, "{cohort}.ug.QC.vcf".format(cohort=cohort))
    log:
        os.path.join(LOGDIR, "variant_filtration_cohort.{cohort}.log".format(cohort=cohort))
    shell:
        r"""
        echo "-GATK VariantFiltration" >> {log}
        {GATK} -T VariantFiltration \
               -R {REF} \
               -o {output.qc_vcf} \
               --variant {input.vqsr_vcf} \
               --clusterWindowSize 10 \
               --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
               --filterName "HARD_TO_VALIDATE" \
               --filterExpression "DP < 5" \
               --filterName "LowCoverage" \
               --filterExpression "QUAL < 30.0" \
               --filterName "VeryLowQual" \
               --filterExpression "QUAL > 30.0 && QUAL < 50.0" \
               --filterName "LowQual" \
               --filterExpression "QD < 2.0" \
               --filterName "LowQD" \
               --filterExpression "MQ < 40.0" \
               --filterName "LowMQ" \
               --filterExpression "FS > 60.0" \
               --filterName "StrandBias" 2>> {log}
        """

#############################################
# Rule: jaccard_index
# Description: Run the jaccard.sh script and produce jaccard.txt and jaccard_lt50.txt.
#############################################
rule jaccard_index:
    input:
        qc_vcf = os.path.join(VARCALLDIR, "{cohort}.ug.QC.vcf".format(cohort=cohort))
    output:
        jaccard      = os.path.join(STATSDIR, "jaccard.txt"),
        jaccard_lt50 = os.path.join(STATSDIR, "jaccard_lt50.txt")
    log:
        os.path.join(LOGDIR, "jaccard_index.{cohort}.log".format(cohort=cohort))
    shell:
        r"""
        {JACCARD} snakemake > {output.jaccard} 2>> {log}
        grep P {output.jaccard} | awk '$NF < 0.5' > {output.jaccard_lt50} 2>> {log} || touch {output.jaccard_lt50}
        """

#############################################
# Rule: cohort_pipeline
# Aggregate rule to run the full cohort analysis.
#############################################
rule cohort_pipeline:
    input:
        merged_vcf = os.path.join(VARCALLDIR, "{cohort}.ug.raw.vcf".format(cohort=cohort)),
        snp_recal  = os.path.join(VARCALLDIR, "{cohort}.ug.raw.snp.recal".format(cohort=cohort)),
        snp_tranches = os.path.join(VARCALLDIR, "{cohort}.ug.raw.snp.tranches".format(cohort=cohort)),
        indel_recal = os.path.join(VARCALLDIR, "{cohort}.ug.raw.indel.recal".format(cohort=cohort)),
        indel_tranches = os.path.join(VARCALLDIR, "{cohort}.ug.raw.indel.tranches".format(cohort=cohort)),
        qc_vcf = os.path.join(VARCALLDIR, "{cohort}.ug.QC.vcf".format(cohort=cohort)),
        jaccard = os.path.join(STATSDIR, "jaccard.txt"),
        jaccard_lt50 = os.path.join(STATSDIR, "jaccard_lt50.txt")
    output:
        "cohort_pipeline.done"
    shell:
        r"""
        echo "Cohort pipeline finished successfully" > {output}
        """
