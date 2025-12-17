#############################################
# Snakefile for Single-Sample Exome Pipeline
#
# This workflow uses config.yaml to store all parameters.
#
# $VERSION taken from CBICall
#############################################

import os
import glob
import platform
from pathlib import Path

snakefile_dir = Path(workflow.snakefile).parent
configfile: snakefile_dir / "config.yaml"

# Scripts (relative to the Snakefile directory)
COV = os.path.join(snakefile_dir, "coverage.sh")
VCF2SEX = os.path.join(snakefile_dir, "vcf2sex.sh")

# Global variables from config
DATADIR = config["datadir"]
DBDIR = config["dbdir"].format(datadir=DATADIR)
NGSUTILS = config["ngsutils"].format(datadir=DATADIR)
TMPDIR = config["tmpdir"].format(datadir=DATADIR)
MEM = config["mem"]

# Threads from --cores
THREADS = workflow.cores if workflow.cores else 4

# Determine JAVA and tool paths based on architecture
ARCH = platform.machine()
if ARCH == "aarch64":
    JAVA = config["java"]["aarch64"]
    BWA = config["tools"]["aarch64"]["bwa"].format(ngsutils=NGSUTILS)
    SAM = config["tools"]["aarch64"]["samtools"].format(ngsutils=NGSUTILS)
    BED = config["tools"]["aarch64"]["bedtools"].format(ngsutils=NGSUTILS)
else:
    JAVA = config["java"]["amd64"]
    BWA = config["tools"]["amd64"]["bwa"].format(ngsutils=NGSUTILS)
    SAM = config["tools"]["amd64"]["samtools"].format(ngsutils=NGSUTILS)
    BED = config["tools"]["amd64"]["bedtools"].format(ngsutils=NGSUTILS)

# Build PIC and GATK commands
PIC = config["picard"].format(java=JAVA, mem=MEM, tmpdir=TMPDIR, ngsutils=NGSUTILS)
GATK = config["gatk"].format(java=JAVA, mem=MEM, tmpdir=TMPDIR, ngsutils=NGSUTILS)

# Other parameters from config
bundle = config["bundle"].format(dbdir=DBDIR)
REF = config["ref"].format(bundle=bundle)
REFGZ = config["refgz"].format(bundle=bundle)
dbSNP = config["dbsnp"].format(dbdir=DBDIR)
MILLS_INDELS = config["mills_indels"].format(bundle=bundle)
KG_INDELS = config["kg_indels"].format(bundle=bundle)
HAPMAP = config["hapmap"].format(bundle=bundle)
OMNI = config["omni"].format(bundle=bundle)

snp_res = config["snp_res"].format(hapmap=HAPMAP, omni=OMNI, dbsnp=dbSNP)
indel_res = config["indel_res"].format(mills_indels=MILLS_INDELS)
EXOM = config["exom"].format(dbdir=DBDIR)

DCOV = config["dcov"]
UG_CALL = config["ug_call"]
UG_EMIT = config["ug_emit"]

# Working directories
BAMDIR = "01_bam"
VARCALLDIR = "02_varcall"
STATSDIR = "03_stats"
LOGDIR = "logs"

# Create directories if not present (include logs)
for d in [BAMDIR, VARCALLDIR, STATSDIR, LOGDIR]:
    os.makedirs(d, exist_ok=True)

# FASTQs are expected in the parent directory
FASTQ_DIR = "../"
FASTQ_R1 = sorted(glob.glob(os.path.join(FASTQ_DIR, "*_R1_*fastq.gz")))

if not FASTQ_R1:
    raise ValueError("No R1 FASTQs found in ../ (expected pattern: *_R1_*fastq.gz)")

# Build FASTQ pairs (bash-equivalent pairing: _R1_ -> _R2_)
FASTQ_PAIRS = []
for r1 in FASTQ_R1:
    r2 = r1.replace("_R1_", "_R2_")
    base = os.path.basename(r1)
    FASTQ_PAIRS.append({"r1": r1, "r2": r2, "base": base})

FASTQ_DICT = {pair["base"]: pair for pair in FASTQ_PAIRS}

# ID parity with bash:
# bash: id=$( echo $DIR | awk -F'/' '{print $(NF-1)}' | awk -F'_' '{print $1}' )
# i.e. take parent folder name and split by "_" keep first field
ID = Path(os.getcwd()).parent.name.split("_")[0]

# Chromosome chunks
CHROMOSOMES = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
    "1314", "1516", "1718", "1920", "2122XY"
]

# Coverage chromosome name parity with bash logic:
# if REF contains b37 then use "1" else "chr1"
CHR_FOR_STATS = "1" if "b37" in str(REF) else "chr1"

#############################################
# Rule: all
#############################################
rule all:
    input:
        os.path.join(VARCALLDIR, "{id}.ug.QC.vcf").format(id=ID),
        os.path.join(STATSDIR, "{id}.coverage.txt").format(id=ID),
        os.path.join(STATSDIR, "{id}.sex.txt").format(id=ID)

#############################################
# Rule: align_and_fix
#############################################
rule align_and_fix:
    input:
        r1=lambda wc: FASTQ_DICT[wc.base]["r1"],
        r2=lambda wc: FASTQ_DICT[wc.base]["r2"],
    output:
        grp=os.path.join(BAMDIR, "{base}.grp.bam"),
        fixed=os.path.join(BAMDIR, "{base}.fixed.sorted.bam"),
    threads: THREADS
    log:
        os.path.join(LOGDIR, "{base}.align_and_fix.log")
    shell:
        r"""
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
#############################################
rule merge_bams:
    input:
        fixed_bams=expand(
            os.path.join(BAMDIR, "{base}.fixed.sorted.bam"),
            base=[pair["base"] for pair in FASTQ_PAIRS],
        )
    output:
        merged=os.path.join(BAMDIR, "input.merged.bam")
    log:
        os.path.join(LOGDIR, "merge_bams.log")
    shell:
        r"""
        tmp_in=$(for bam in {input.fixed_bams}; do echo -n " I=$bam"; done)
        echo "Merging BAM files" >> {log}
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
#############################################
rule filter_bam:
    input:
        merged=os.path.join(BAMDIR, "input.merged.bam")
    output:
        filtered=os.path.join(BAMDIR, "input.merged.filtered.bam"),
        filtered_index=os.path.join(BAMDIR, "input.merged.filtered.bam.bai"),
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
#############################################
rule realign_target_creator:
    input:
        filtered_bam=os.path.join(BAMDIR, "input.merged.filtered.bam")
    output:
        intervals=os.path.join(BAMDIR, "input.merged.filtered.intervals")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "realign_target_creator.log")
    shell:
        r"""
        echo "Running RealignerTargetCreator" >> {log}
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
#############################################
rule indel_realigner:
    input:
        bam=os.path.join(BAMDIR, "input.merged.filtered.bam"),
        intervals=os.path.join(BAMDIR, "input.merged.filtered.intervals"),
    output:
        realigned_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.bam")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "indel_realigner.log")
    shell:
        r"""
        echo "Running IndelRealigner" >> {log}
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
#############################################
rule fix_mate_information:
    input:
        realigned=os.path.join(BAMDIR, "input.merged.filtered.realigned.bam")
    output:
        fixed=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bam"),
        fixed_index=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bai"),
    log:
        os.path.join(LOGDIR, "fix_mate_information.log")
    shell:
        r"""
        echo "Fixing mate information (post realignment)" >> {log}
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
#############################################
rule mark_duplicates:
    input:
        fixed_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bam")
    output:
        dedup_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bam"),
        dedup_index=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bai"),
        metrics=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.dupmetrics"),
    log:
        os.path.join(LOGDIR, "mark_duplicates.log")
    shell:
        r"""
        echo "Marking and deleting PCR duplicates" >> {log}
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
#############################################
rule base_recalibrator:
    input:
        dedup_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bam")
    output:
        recal_grp=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.grp")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "base_recalibrator.chr{chr}.log")
    shell:
        r"""
        echo "Running BaseRecalibrator for chr{wildcards.chr}" >> {log}
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
#############################################
rule print_reads:
    input:
        dedup_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bam"),
        bqsr_grp=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.grp"),
    output:
        recal_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.bam")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "print_reads.chr{chr}.log")
    shell:
        r"""
        echo "Running PrintReads for chr{wildcards.chr}" >> {log}
        {GATK} -T PrintReads \
               -R {REF} \
               -L {EXOM}/hg19.chr{wildcards.chr}.bed \
               -I {input.dedup_bam} \
               -o {output.recal_bam} \
               -BQSR {input.bqsr_grp} 2>> {log}
        """

#############################################
# Rule: all_bqsr
#############################################
rule all_bqsr:
    input:
        expand(
            os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.bam"),
            chr=CHROMOSOMES,
        )

#############################################
# Rule: variant_calling
#############################################
rule variant_calling:
    input:
        recal_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.chr{chr}.recal.bam")
    output:
        vcf=os.path.join(VARCALLDIR, "chr{chr}.ug.raw.vcf")
    threads: THREADS
    log:
        os.path.join(LOGDIR, "variant_calling.chr{chr}.log")
    shell:
        r"""
        echo "Running UnifiedGenotyper for chr{wildcards.chr}" >> {log}
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
        """

#############################################
# Rule: all_variant_calling
#############################################
rule all_variant_calling:
    input:
        expand(os.path.join(VARCALLDIR, "chr{chr}.ug.raw.vcf"), chr=CHROMOSOMES)

#############################################
# Rule: merge_vcfs  (match bash output name)
#############################################
rule merge_vcfs:
    input:
        vcfs=expand(os.path.join(VARCALLDIR, "chr{chr}.ug.raw.vcf"), chr=CHROMOSOMES)
    output:
        merged_vcf=os.path.join(VARCALLDIR, "{id}.ug.raw.vcf")
    log:
        os.path.join(LOGDIR, "merge_vcfs.{id}.log")
    shell:
        r"""
        echo "Merging VCF files for sample {wildcards.id}" >> {log}
        ( grep "#" {input.vcfs[0]} ; grep -hv '#' {VARCALLDIR}/chr*.ug.raw.vcf | awk '$1 !~ /_/' | sort -V ) > {output.merged_vcf} 2>> {log}
        """

#############################################
# Rule: recalibrate_snp
#############################################
rule recalibrate_snp:
    input:
        vcf=os.path.join(VARCALLDIR, "{id}.ug.raw.vcf")
    output:
        snp_recal=os.path.join(VARCALLDIR, "{id}.ug.raw.snp.recal"),
        snp_tranches=os.path.join(VARCALLDIR, "{id}.ug.raw.snp.tranches"),
    params:
        snp_res=snp_res
    log:
        os.path.join(LOGDIR, "recalibrate_snp.{id}.log")
    shell:
        r"""
        echo "VariantRecalibrator SNP" >> {log}
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
            cp {input.vcf} {output.snp_recal}
            echo "SKIP" > {output.snp_tranches}
        fi
        """

#############################################
# Rule: recalibrate_indel
#############################################
rule recalibrate_indel:
    input:
        vcf=os.path.join(VARCALLDIR, "{id}.ug.raw.vcf")
    output:
        indel_recal=os.path.join(VARCALLDIR, "{id}.ug.raw.indel.recal"),
        indel_tranches=os.path.join(VARCALLDIR, "{id}.ug.raw.indel.tranches"),
    params:
        indel_res=indel_res
    log:
        os.path.join(LOGDIR, "recalibrate_indel.{id}.log")
    shell:
        r"""
        echo "VariantRecalibrator INDEL" >> {log}
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
            echo "No INDELs to recalibrate (only $nINDEL found)" >> {log}
            cp {input.vcf} {output.indel_recal}
            echo "SKIP" > {output.indel_tranches}
        fi
        """

#############################################
# Rule: apply_recalibration_snp
#############################################
rule apply_recalibration_snp:
    input:
        vcf=os.path.join(VARCALLDIR, "{id}.ug.raw.vcf"),
        snp_recal=os.path.join(VARCALLDIR, "{id}.ug.raw.snp.recal"),
        snp_tranches=os.path.join(VARCALLDIR, "{id}.ug.raw.snp.tranches"),
    output:
        recal_snp=os.path.join(VARCALLDIR, "recalibratedSNPs.rawIndels.{id}.vcf")
    log:
        os.path.join(LOGDIR, "apply_recalibration_snp.{id}.log")
    shell:
        r"""
        if grep -q "SKIP" {input.snp_tranches}; then
            echo "Skipping ApplyRecalibration for SNPs due to low variant count" >> {log}
            cp {input.vcf} {output.recal_snp}
        else
            echo "ApplyRecalibration SNP" >> {log}
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
#############################################
rule apply_recalibration_indel:
    input:
        recal_snp=os.path.join(VARCALLDIR, "recalibratedSNPs.rawIndels.{id}.vcf"),
        indel_recal=os.path.join(VARCALLDIR, "{id}.ug.raw.indel.recal"),
        indel_tranches=os.path.join(VARCALLDIR, "{id}.ug.raw.indel.tranches"),
    output:
        vqsr_vcf=os.path.join(VARCALLDIR, "{id}.ug.vqsr.vcf")
    log:
        os.path.join(LOGDIR, "apply_recalibration_indel.{id}.log")
    shell:
        r"""
        echo "ApplyRecalibration INDEL" >> {log}
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
# Rule: variant_filtration
#############################################
rule variant_filtration:
    input:
        vqsr_vcf=os.path.join(VARCALLDIR, "{id}.ug.vqsr.vcf")
    output:
        qc_vcf=os.path.join(VARCALLDIR, "{id}.ug.QC.vcf")
    log:
        os.path.join(LOGDIR, "variant_filtration.{id}.log")
    shell:
        r"""
        echo "VariantFiltration" >> {log}
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
#############################################
rule coverage_stats:
    input:
        raw_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bam"),
        raw_bai=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.bai"),
        dedup_bam=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bam"),
        dedup_bai=os.path.join(BAMDIR, "input.merged.filtered.realigned.fixed.dedup.bai"),
    output:
        out_raw=os.path.join(STATSDIR, "{id}.{chrN}.raw.bam"),
        out_dedup=os.path.join(STATSDIR, "{id}.{chrN}.dedup.bam"),
        stats_log=os.path.join(STATSDIR, "{id}.coverage.txt"),
    params:
        chrN=CHR_FOR_STATS
    log:
        os.path.join(LOGDIR, "coverage_stats.{id}.log")
    shell:
        r"""
        echo "Computing coverage stats for {params.chrN}" >> {log}

        # Match bash behavior: create bam.bai alongside bam from the provided bai files
        cp {input.raw_bai} {input.raw_bam}.bai 2>> {log} || true
        cp {input.dedup_bai} {input.dedup_bam}.bai 2>> {log} || true

        {SAM} view -b {input.raw_bam} {params.chrN} > {output.out_raw} 2>> {log}
        {SAM} view -b {input.dedup_bam} {params.chrN} > {output.out_dedup} 2>> {log}
        {SAM} index {output.out_raw} 2>> {log}
        {SAM} index {output.out_dedup} 2>> {log}

        {COV} {wildcards.id} {output.out_raw} {output.out_dedup} > {output.stats_log} 2>> {log}
        """

#############################################
# Rule: sex_determination
#############################################
rule sex_determination:
    input:
        vcf=os.path.join(VARCALLDIR, "{id}.ug.QC.vcf")
    output:
        sex=os.path.join(STATSDIR, "{id}.sex.txt")
    log:
        os.path.join(LOGDIR, "sex_determination.{id}.log")
    shell:
        r"""
        echo "Estimating sex for id {wildcards.id}" >> {log}
        {VCF2SEX} {input.vcf} > {output.sex} 2>> {log}
        """

