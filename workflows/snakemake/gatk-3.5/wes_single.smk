#############################################
# Snakefile for Single-Sample Exome Pipeline
#
# This workflow uses config.yaml to store all parameters.
#
# $VERSION taken from CBICall
#############################################

import os, glob, platform
from pathlib import Path

snakefile_dir = Path(workflow.snakefile).parent
configfile: snakefile_dir / "config.yaml"

# Define paths for the coverage and sex determination scripts
COV = os.path.join(snakefile_dir, "coverage.sh")
VCF2SEX = os.path.join(snakefile_dir, "vcf2sex.sh")

# Global variables from config
DATADIR  = config["datadir"]
DBDIR    = config["dbdir"].format(datadir=DATADIR)
NGSUTILS = config["ngsutils"].format(datadir=DATADIR)
TMPDIR   = config["tmpdir"].format(datadir=DATADIR)
MEM      = config["mem"]

# Get the number of cores from the CLI argument --cores
THREADS = workflow.cores if workflow.cores else 4  # Default to 4 if not set

# Determine JAVA and tool paths based on architecture
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

# Build PIC and GATK commands using the config strings and substituted variables
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

# Working directories (adjust if needed)
BAMDIR      = "01_bam"
VARCALLDIR  = "02_varcall"
STATSDIR    = "03_stats"
LOGDIR      = "logs"

# Create directories if not present
for d in [BAMDIR, VARCALLDIR, STATSDIR]:
    os.makedirs(d, exist_ok=True)

# Define input FASTQ files (assumes they are in the parent directory)
FASTQ_DIR  = "../"
FASTQ_R1   = sorted(glob.glob(os.path.join(FASTQ_DIR, "*_R1_*fastq.gz")))

# Build FASTQ pairs (assuming filenames like MA0004701P_ex_S5_L001_R1_001.fastq.gz)
FASTQ_PAIRS = []
for r1 in FASTQ_R1:
    r2 = r1.replace("R1", "R2")
    base = os.path.basename(r1)
    FASTQ_PAIRS.append({"r1": r1, "r2": r2, "base": base})

# Extract sample ID from the first FASTQ (assumes first two underscore fields make the sample)
sample_id = os.path.basename(FASTQ_PAIRS[0]["r1"]).split("_")[0] + "_" + os.path.basename(FASTQ_PAIRS[0]["r1"]).split("_")[1]
if sample_id.endswith("_ex"):
    sample_id = sample_id[:-3]
ID = sample_id

# Define chromosomes to process
CHROMOSOMES = ["1","2","3","4","5","6","7","8","9","10","11","12", "1314", "1516", "1718", "1920", "2122XY"]

# Define a mapping from FASTQ base name to full paths
FASTQ_DICT = {pair["base"]: pair for pair in FASTQ_PAIRS}

#############################################
# Rule: all
# Description: Global target to build all files
# To run in isolation:
#   snakemake -j 4 
#############################################
rule all:
    input:
        os.path.join(VARCALLDIR, "{id}.ug.QC.vcf").format(id=ID),
        os.path.join(STATSDIR, "{id}.coverage.txt").format(id=ID),
        os.path.join(STATSDIR, "{id}.sex.txt").format(id=ID)

#############################################
# Rule: align_and_fix
# Description: Align reads using BWA, add read groups with Picard, and then fix mate information.
# To run in isolation (example for one FASTQ pair):
#   snakemake -j 4 01_bam/<FASTQ_BASE>.fixed.sorted.bam
#############################################
rule align_and_fix:
    input:
        r1 = lambda wc: FASTQ_DICT[wc.base]["r1"],
        r2 = lambda wc: FASTQ_DICT[wc.base]["r2"]
    output:
        grp = os.path.join(BAMDIR, "{base}.grp.bam"),
        fixed = os.path.join(BAMDIR, "{base}.fixed.sorted.bam")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "{base}.align_and_fix.log")
    shell:
        r"""
        # Extract sample information from the FASTQ basename
        SAMPLE=$(echo {wildcards.base} | cut -d'_' -f1-2)
        BARCODE=$(echo {wildcards.base} | cut -d'_' -f3)
        LANE=$(echo {wildcards.base} | cut -d'_' -f4)
        
        echo "Aligning {input.r1} and {input.r2}" >> {log}
        {BWA} mem -t{threads} -M {REFGZ} {input.r1} {input.r2} 2>> {log} | \
            {PIC} AddOrReplaceReadGroups \
                  TMP_DIR={TMPDIR} \
                  I=/dev/stdin \
                  O={output.grp} \
                  SO=coordinate \
                  RGID=$LANE \
                  RGLB=sureselect \
                  RGPL=illumina \
                  RGPU=$BARCODE \
                  RGSM=$SAMPLE 2>> {log}
        
        echo "Fixing mate information for {output.grp}" >> {log}
        {PIC} FixMateInformation \
              TMP_DIR={TMPDIR} \
              INPUT={output.grp} \
              OUTPUT={output.fixed} \
              VALIDATION_STRINGENCY=SILENT \
              CREATE_INDEX=true 2>> {log}
        """

#############################################
# Rule: merge_bams
# Description: Merge all fixed BAM files into a single BAM using Picard MergeSamFiles.
# To run in isolation:
#   snakemake -j 4 01_bam/input.merged.bam
#############################################
rule merge_bams:
    input:
        fixed_bams = expand(os.path.join(BAMDIR, "{base}.fixed.sorted.bam"), base=[pair["base"] for pair in FASTQ_PAIRS])
    output:
        merged = os.path.join(BAMDIR, "input.merged.bam")
    log:
        os.path.join(LOGDIR, "merge_bams.log")
    shell:
        r"""
        # Build the list of input BAMs with the required "I=" prefix.
        tmp_in=$(for bam in {input.fixed_bams}; do echo -n " I=$bam"; done)
        echo "Merging BAM files: {input.fixed_bams}" >> {log}
        {PIC} MergeSamFiles \
             TMP_DIR={TMPDIR} \
             $tmp_in \
             OUTPUT={output.merged} \
             SO=coordinate \
             VALIDATION_STRINGENCY=SILENT \
             CREATE_INDEX=false 2>> {log}
        """

#############################################
# Rule: filter_bam
# Description: Filter out malformed reads from the merged BAM using SAMtools and awk.
# To run in isolation:
#   snakemake -j 4 01_bam/input.merged.filtered.bam
#############################################
rule filter_bam:
    input:
        merged = os.path.join(BAMDIR, "input.merged.bam")
    output:
        filtered = os.path.join(BAMDIR, "input.merged.filtered.bam"),
        filtered_index = os.path.join(BAMDIR, "input.merged.filtered.bam.bai")
    log:
        os.path.join(LOGDIR, "filter_bam.log")
    shell:
        r"""
        {SAM} view -h {input.merged} | \
          awk 'substr($0,1,1)=="@" || ($6 !~ /[0-9]+H/ && length($10)==length($11))' | \
          {SAM} view -Sb - > {output.filtered}
        {SAM} index {output.filtered} 2>> {log}
        """

#############################################
# Rule: realign_target_creator
# Description: Identify intervals for local realignment around indels using GATK RealignerTargetCreator.
# To run in isolation:
#   snakemake -j 4 01_bam/input.merged.filtered.bam.intervals
#############################################
rule realign_target_creator:
    input:
        filtered_bam = os.path.join(BAMDIR, "input.merged.filtered.bam")
    output:
        intervals = os.path.join(BAMDIR, "input.merged.filtered.intervals")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "realign_target_creator.log")
    shell:
        r"""
        echo "Running RealignerTargetCreator on {input.filtered_bam}" >> {log}
        {GATK} -T RealignerTargetCreator \
               -nt {threads} \
               -R {REF} \
               -I {input.filtered_bam} \
               -o {output.intervals} \
               -known {MILLS_INDELS} \
               -known {KG_INDELS} 2>> {log}
        """

#############################################
# Rule: indel_realigner
# Description: Perform local realignment around indels using GATK IndelRealigner.
# To run in isolation:
#   snakemake -j 4 01_bam/input.merged.filtered.realigned.bam
#############################################
rule indel_realigner:
    input:
        bam = os.path.join(BAMDIR, "input.merged.filtered.bam"),
        intervals = os.path.join(BAMDIR, "input.merged.filtered.intervals")
    output:
        realigned_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.bam")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "indel_realigner.log")
    shell:
        r"""
        echo "Running IndelRealigner on {input.bam} with intervals {input.intervals}" >> {log}
        {GATK} -T IndelRealigner \
               -R {REF} \
               -targetIntervals {input.intervals} \
               -I {input.bam} \
               -o {output.realigned_bam} \
               -model USE_SW \
               -known {MILLS_INDELS} \
               -known {KG_INDELS} \
               -rf NotPrimaryAlignment 2>> {log}
        """

#############################################
# Rule: fix_mate_information
# Description: Fix mate information and re-sort using Picard FixMateInformation.
# To run in isolation:
#   snakemake -j 4 01_bam/input.merged.filtered.realigned.fixed.bam
#############################################
rule fix_mate_information:
    input:
        realigned = os.path.join(BAMDIR, "input.merged.filtered.realigned.bam")
    output:
        fixed = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bam"),
        fixed_index = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bai")
    log:
        os.path.join(LOGDIR, "fix_mate_information.log")
    shell:
        r"""
        echo "Fixing mate information for {input.realigned}" >> {log}
        {PIC} FixMateInformation \
              TMP_DIR={TMPDIR} \
              INPUT={input.realigned} \
              OUTPUT={output.fixed} \
              SO=coordinate \
              VALIDATION_STRINGENCY=LENIENT \
              CREATE_INDEX=true 2>> {log}
        """

#############################################
# Rule: mark_duplicates
# Description: Mark and remove PCR duplicates using Picard MarkDuplicates.
# To run in isolation:
#   snakemake -j 4 01_bam/input.merged.filtered.realigned.fixed.dedup.bam
#############################################
rule mark_duplicates:
    input:
        fixed_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bam")
    output:
        dedup_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bam"),
        dedup_index = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bai"),
        metrics = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.dupmetrics")
    log:
        os.path.join(LOGDIR, "mark_duplicates.log")
    shell:
        r"""
        echo "Marking and deleting PCR duplicates for {input.fixed_bam}" >> {log}
        {PIC} MarkDuplicates \
              TMP_DIR={TMPDIR} \
              INPUT={input.fixed_bam} \
              OUTPUT={output.dedup_bam} \
              METRICS_FILE={output.metrics} \
              REMOVE_DUPLICATES=true \
              ASSUME_SORTED=true \
              CREATE_INDEX=true \
              VALIDATION_STRINGENCY=SILENT 2>> {log}
        """

#############################################
# Rule: base_recalibrator
# Description: Generate a recalibration group file per chromosome using GATK BaseRecalibrator.
# To run in isolation (for a specific chromosome, e.g., chr1):
#   snakemake -j 4 01_bam/input.merged.filtered.realigned.fixed.dedup.chr1.recal.grp
#############################################
rule base_recalibrator:
    input:
        dedup_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bam")
    output:
        recal_grp = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.grp")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "base_recalibrator.{chr}.log")
    shell:
        r"""
        echo "Running BaseRecalibrator for chromosome {wildcards.chr}" >> {log}
        {GATK} -T BaseRecalibrator \
               -nct {threads} \
               -R {REF} \
               -L {EXOM}/hg19.chr{wildcards.chr}.bed \
               -I {input.dedup_bam} \
               -o {output.recal_grp} \
               -knownSites {dbSNP} \
               -knownSites {MILLS_INDELS} \
               -knownSites {KG_INDELS} 2>> {log}
        """

#############################################
# Rule: print_reads
# Description: Apply recalibration to deduplicated BAM files using GATK PrintReads.
# To run in isolation (for a specific chromosome, e.g., chr1):
#   snakemake -j 4 01_bam/input.merged.filtered.realigned.fixed.dedup.chr1.recal.bam
#############################################
rule print_reads:
    input:
        dedup_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bam"),
        bqsr_grp = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.grp")
    output:
        recal_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.bam")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "print_reads.chr{chr}.log")
    shell:
        r"""
        echo "Running PrintReads for chromosome {wildcards.chr}" >> {log}
        {GATK} -T PrintReads \
               -R {REF} \
               -L {EXOM}/hg19.chr{wildcards.chr}.bed \
               -I {input.dedup_bam} \
               -o {output.recal_bam} \
               -BQSR {input.bqsr_grp} 2>> {log}
        """

#############################################
# Rule: all_bqsr
# Description: Aggregate target to run BQSR for all chromosomes.
# To run in isolation:
#   snakemake -j 4 all_bqsr
#############################################
rule all_bqsr:
    input:
        expand(os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.bam"), chr=CHROMOSOMES)

#############################################
# Rule: variant_calling
# Description: Call variants per chromosome using GATK UnifiedGenotyper.
# To run in isolation (for a specific chromosome, e.g., chr1):
#   snakemake -j 4 02_varcall/chr1.ug.raw.vcf
#############################################
rule variant_calling:
    input:
        recal_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.bam")
    output:
        vcf = os.path.join(VARCALLDIR, "chr{chr}.ug.raw.vcf")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "variant_calling.chr{chr}.log")
    shell:
        r"""
        echo "Running GATK UnifiedGenotyper for chromosome {wildcards.chr}" >> {log}
        {GATK} -T UnifiedGenotyper \
               -R {REF} \
               -I {input.recal_bam} \
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
# Rule: all_variant_calling
# Description: Aggregate target to run variant calling for all chromosomes.
# To run in isolation:
#   snakemake -j 4 all_variant_calling
#############################################
rule all_variant_calling:
    input:
        expand(os.path.join(VARCALLDIR, "chr{chr}.ug.raw.vcf"), chr=CHROMOSOMES)

#############################################
# Rule: merge_vcfs
# Description: Merge per-chromosome VCF files into a single VCF.
# To run in isolation:
#   snakemake -j 4 02_varcall/<sample>.ug.raw.vcf
#############################################
rule merge_vcfs:
    input: 
        vcfs = expand(os.path.join(VARCALLDIR, "chr{chr}.ug.raw.vcf"), chr=CHROMOSOMES)
    output:
        merged_vcf = os.path.join(VARCALLDIR, "merged.{id}.ug.raw.vcf")
    log: 
        os.path.join(LOGDIR, "merge_vcfs.{id}.log")
    shell: 
        r"""
        echo "Merging VCF files for sample {wildcards.id}" >> {log}
        ( grep "#" {input.vcfs[0]} ; grep -hv '#' {VARCALLDIR}/chr*.ug.raw.vcf | awk '$1 !~ /_/' | sort -V ) > {output.merged_vcf} 2>> {log}
        """

#############################################
# Rule: recalibrate_snp
# Description: Recalibrate SNP variants using GATK VariantRecalibrator.
# To run in isolation:
#   snakemake -j 4 02_varcall/<sample>.ug.raw.snp.recal
#############################################
rule recalibrate_snp:
    input:
        vcf = os.path.join(VARCALLDIR, "merged.{id}.ug.raw.vcf")
    output:
        snp_recal = os.path.join(VARCALLDIR, "{id}.ug.raw.snp.recal"),
        snp_tranches = os.path.join(VARCALLDIR, "{id}.ug.raw.snp.tranches")
    params:
        snp_res = snp_res  # defined in config/header
    log:
        os.path.join(LOGDIR, "recalibrate_snp.{id}.log")
    shell:
        r"""
        echo "-GATK Recalibrator SNP" >> {log}
        # Count SNP variants (assuming SNPs have a single-base alternate allele)
        nSNP=$(grep -v "^#" {input.vcf} | awk 'length($5)==1' | wc -l)
        echo "Found $nSNP SNP variants" >> {log}
        if [ $nSNP -gt 1000 ]; then
            {GATK} -T VariantRecalibrator \
                   -R {REF} \
                   -input {input.vcf} \
                   -recalFile {output.snp_recal} \
                   -tranchesFile {output.snp_tranches} \
                   --maxGaussians 6 \
                   {params.snp_res} \
                   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
                   -mode SNP 2>> {log}
        else
            echo "WARNING: Only $nSNP SNP variants found. Skipping VariantRecalibrator for SNPs." >> {log}
            # Pass through the original VCF as recalibration output.
            cp {input.vcf} {output.snp_recal}
            # Write a marker in the tranches file to indicate skipping.
            echo "SKIP" > {output.snp_tranches}
        fi
        """

#############################################
# Rule: recalibrate_indel
# Description: Recalibrate INDEL variants using GATK VariantRecalibrator.
# To run in isolation:
#   snakemake -j 4 02_varcall/<sample>.ug.raw.indel.recal
#############################################
rule recalibrate_indel:
    input:
        vcf = os.path.join(VARCALLDIR, "merged.{id}.ug.raw.vcf")
    output:
        indel_recal = os.path.join(VARCALLDIR, "{id}.ug.raw.indel.recal"),
        indel_tranches = os.path.join(VARCALLDIR, "{id}.ug.raw.indel.tranches")
    params:
        indel_res = indel_res  # defined in config/header
    log:
        os.path.join(LOGDIR, "recalibrate_indel.{id}.log")
    shell:
        r"""
        echo "-GATK Recalibrator INDEL" >> {log}
        nINDEL=$(grep -v "^#" {input.vcf} | awk 'length($5) != 1' | wc -l)
        echo "Found $nINDEL INDEL variants" >> {log}
        if [ $nINDEL -gt 8000 ]; then
            {GATK} -T VariantRecalibrator \
                   -R {REF} \
                   -input {input.vcf} \
                   -recalFile {output.indel_recal} \
                   -tranchesFile {output.indel_tranches} \
                   --maxGaussians 4 \
                   {params.indel_res} \
                   -an QD -an FS -an ReadPosRankSum \
                   -mode INDEL 2>> {log}
        else
            echo ">>>> No INDELs to recalibrate (only $nINDEL found)" >> {log}
            cp {input.vcf} {output.indel_recal}
            echo "SKIP" > {output.indel_tranches}
        fi
        """

#############################################
# Rule: apply_recalibration_snp
# Description: Apply SNP recalibration to generate an intermediate VCF.
# To run in isolation:
#   snakemake -j 4 02_varcall/{id}.recalibratedSNPs.rawIndels.vcf
#############################################
rule apply_recalibration_snp:
    input:
        vcf = os.path.join(VARCALLDIR, "merged.{id}.ug.raw.vcf"),
        snp_recal = os.path.join(VARCALLDIR, "{id}.ug.raw.snp.recal"),
        snp_tranches = os.path.join(VARCALLDIR, "{id}.ug.raw.snp.tranches")
    output:
        recal_snp = os.path.join(VARCALLDIR, "{id}.recalibratedSNPs.rawIndels.vcf")
    log:
        os.path.join(LOGDIR, "apply_recalibration_snp.{id}.log")
    shell:
        r"""
        # Check if recalibration was skipped
        if grep -q "SKIP" {input.snp_tranches}; then
            echo "Skipping ApplyRecalibration for SNPs due to low variant count" >> {log}
            cp {input.vcf} {output.recal_snp}
        else
            echo "-GATK Apply Recalibrator SNP" >> {log}
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
# Rule: apply_recalibration_indel
# Description: Apply INDEL recalibration and generate the final VCF for variant quality score recalibration (VQSR).
# To run in isolation:
#   snakemake -j 4 02_varcall/<sample>.ug.vqsr.vcf
#############################################
rule apply_recalibration_indel:
    input:
        recal_snp = os.path.join(VARCALLDIR, "{id}.recalibratedSNPs.rawIndels.vcf"),
        indel_recal = os.path.join(VARCALLDIR, "{id}.ug.raw.indel.recal"),
        indel_tranches = os.path.join(VARCALLDIR, "{id}.ug.raw.indel.tranches")
    output:
        vqsr_vcf = os.path.join(VARCALLDIR, "{id}.ug.vqsr.vcf")
    log:
        os.path.join(LOGDIR, "apply_recalibration_indel.{id}.log")
    shell:
        r"""
        echo "-GATK Apply Recalibrator INDEL" >> {log}
        if grep -q "SKIP" {input.indel_tranches}; then
            echo ">>>> No INDELs to recalibrate (SKIP marker found)" >> {log}
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
# Rule: variant_filtration
# Description: Apply hard filters to the VQSR VCF to produce a final QC VCF.
# To run in isolation:
#   snakemake -j 4 02_varcall/{id}.ug.QC.vcf
#############################################
rule variant_filtration:
    input:
        vqsr_vcf = os.path.join(VARCALLDIR, "{id}.ug.vqsr.vcf")
    output:
        qc_vcf = os.path.join(VARCALLDIR, "{id}.ug.QC.vcf")
    log:
        os.path.join(LOGDIR, "variant_filtration.{id}.log")
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
# Rule: coverage_stats
# Description: Extracts chromosome-specific BAM files, calculates coverage,
#              nonduplicate percentage, and insert size.
#              It then runs the coverage.sh script from the expected directory.
# To run in isolation:
#   snakemake -j 4 03_stats/{id}.coverage.txt
#############################################
rule coverage_stats:
    input:
        raw_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bam"),
        raw_bai = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bai"),
        dedup_bam = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bam"),
        dedup_bai = os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bai")
    output:
        out_raw = os.path.join(STATSDIR, "{id}.1.raw.bam"),
        out_raw_index = os.path.join(STATSDIR, "{id}.1.raw.bam.bai"),
        out_dedup = os.path.join(STATSDIR, "{id}.1.dedup.bam"),
        out_dedup_index = os.path.join(STATSDIR, "{id}.1.dedup.bam.bai"),
        stats_log = os.path.join(STATSDIR, "{id}.coverage.txt")
    params:
        chrN = "1",   # Adjust if needed; here SAMtools is called with "1".
    log:
        os.path.join(LOGDIR, "coverage_stats.{id}.log")
    shell:
        r"""
        echo "Computing coverage stats for chromosome {params.chrN}" >> {log}
        # Ensure index files exist (copy if needed)
        cp {input.raw_bai} {input.raw_bam}.bai 2>> {log} || true
        cp {input.dedup_bai} {input.dedup_bam}.bai 2>> {log} || true
        {SAM} view -b {input.raw_bam} {params.chrN} > {output.out_raw} 2>> {log}
        {SAM} view -b {input.dedup_bam} {params.chrN} > {output.out_dedup} 2>> {log}
        {SAM} index {output.out_raw} 2>> {log}
        {SAM} index {output.out_dedup} 2>> {log}
        # Call the script using relative paths.
        {COV} {wildcards.id} {output.out_raw} {output.out_dedup} > {STATSDIR}/{wildcards.id}.coverage.txt 2>> logs/coverage_stats.{wildcards.id}.log
        """

#############################################
# Rule: sex_determination
# Description: Run the sex determination script (vcf2sex.sh) on the QC VCF file.
# To run in isolation:
#   snakemake -j 4 03_stats/<sample>.sex.txt
#############################################
rule sex_determination:
    input:
        vcf = os.path.join(VARCALLDIR, "{id}.ug.QC.vcf")
    output:
        sex = os.path.join(STATSDIR, "{id}.sex.txt")
    log:
        os.path.join(LOGDIR, "sex_determination.{id}.log")
    shell:
        r"""
        echo "Estimating sex for id {wildcards.id}" >> {log}
        {VCF2SEX} {input.vcf} > {output.sex} 2>> {log}
        """
