# Snakefile for GATK4.6 WES/WGS pipeline (wes_single.smk)
import os, glob, platform
from pathlib import Path

# ----------------------------------------
# Load configuration and set environment
# ----------------------------------------
snakefile_dir = Path(workflow.snakefile).parent
configfile: snakefile_dir / "config.yaml"

# Export environment variables
os.environ["TMPDIR"] = config["tmpdir"]
os.environ["LC_ALL"] = "C"
os.environ["GATK_DISABLE_AUTO_S3_UPLOAD"] = "true"

# ----------------------------------------
# Define directories, resources, and tools
# ----------------------------------------
DATADIR  = config["datadir"]
DBDIR    = config["dbdir"].format(datadir=DATADIR)
NGSUTILS = config["ngsutils"].format(datadir=DATADIR)
TMPDIR   = config["tmpdir"].format(datadir=DATADIR)
MEM      = config["mem"]

# read from either pipeline (if set) or fall back to mode
PIPELINE = config.get("pipeline", config.get("mode", "wes")).lower()

# Threads
THREADS = workflow.cores or 4

# Choose binaries by architecture
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

# Picard & GATK4 commands
PIC   = config["picard"].format(java=JAVA, mem=MEM, tmpdir=TMPDIR, ngsutils=NGSUTILS)
GATK4 = config["gatk4_cmd"].format(ngsutils=NGSUTILS, mem=MEM)

# Path to helper scripts
COV     = os.path.join(snakefile_dir, "coverage.sh")
VCF2SEX = os.path.join(snakefile_dir, "vcf2sex.sh")

# Reference & resource files
bundle       = config["bundle"].format(dbdir=DBDIR)
REF          = config["ref"].format(bundle=bundle)
REFGZ        = config["refgz"].format(bundle=bundle)
REF_DICT     = config["ref_dict"].format(bundle=bundle)
dbSNP        = config["dbsnp"].format(dbdir=DBDIR)
MILLS_INDELS = config["mills_indels"].format(bundle=bundle)
KG_INDELS    = config["kg_indels"].format(bundle=bundle)
HAPMAP       = config["hapmap"].format(bundle=bundle)
OMNI         = config["omni"].format(bundle=bundle)

# VQSR resource strings
snp_res   = config["snp_res"].format(hapmap=HAPMAP, omni=OMNI, dbsnp=dbSNP)
indel_res = config["indel_res"].format(mills_indels=MILLS_INDELS)

# Interval argument for WES vs WGS
INTERVAL_LIST = config["interval_list"].format(bundle=bundle)
INTERVAL_ARG  = ["-L", INTERVAL_LIST] if PIPELINE == "wes" else []

# Create output directories
BAMDIR     = "01_bam"
VARCALLDIR = "02_varcall"
STATSDIR   = "03_stats"
LOGDIR     = "logs"
for d in [BAMDIR, VARCALLDIR, STATSDIR, LOGDIR]:
    os.makedirs(d, exist_ok=True)

# ----------------------------------------
# Build FASTQ pairs and sample ID
#-----------------------------------------
FASTQ_DIR   = "../"
FASTQ_R1    = sorted(glob.glob(os.path.join(FASTQ_DIR, "*R1*fastq.gz")))
FASTQ_PAIRS = []
for r1 in FASTQ_R1:
    r2   = r1.replace("_R1_", "_R2_")
    base = os.path.basename(r1).replace(".fastq.gz", "").split("_R1_",1)[0]
    FASTQ_PAIRS.append({"r1": r1, "r2": r2, "base": base})

# Sample ID: first two fields, drop trailing _ex
fields    = FASTQ_PAIRS[0]["base"].split("_")
sample_id = fields[0] + "_" + fields[1]
if sample_id.endswith("_ex"):
    sample_id = sample_id[:-3]
ID = sample_id

# Build lookup for FASTQ pairs
FASTQ_DICT = {p["base"]: p for p in FASTQ_PAIRS}

# ----------------------------------------
# Rule: all
# Final targets: QC VCF, coverage & sex stats
#-----------------------------------------
rule all:
    input:
        os.path.join(VARCALLDIR, f"{ID}.hc.QC.vcf.gz"),
        os.path.join(STATSDIR,   f"{ID}.coverage.txt"),
        os.path.join(STATSDIR,   f"{ID}.sex.txt")

# ----------------------------------------
# Rule: align_rg
# Align FASTQ pairs and add read-groups
#-----------------------------------------
rule align_rg:
    input:
        r1=lambda wc: FASTQ_DICT[wc.base]["r1"],
        r2=lambda wc: FASTQ_DICT[wc.base]["r2"]
    output:
        bam=os.path.join(BAMDIR, "{base}.rg.bam")
    params:
        base=lambda wc: wc.base
    threads: THREADS
    log:
        os.path.join(LOGDIR, "{base}.align.log")
    shell:
        r"""
        # RGSM: sample, RGPU: lane, RGID: read-group ID
        SAMPLE=$(echo {params.base} | cut -d'_' -f1-2)
        RGPU=$(echo {params.base} | cut -d'_' -f3)
        RGID=$(echo {params.base} | cut -d'_' -f4)
        {BWA} mem -M -t{threads} {REFGZ} {input.r1} {input.r2} \
          | {GATK4} AddOrReplaceReadGroups \
              --INPUT /dev/stdin \
              --OUTPUT {output.bam} \
              --TMP_DIR {TMPDIR} \
              --RGPL ILLUMINA --RGLB sureselect \
              --RGSM $SAMPLE --RGPU $RGPU --RGID $RGID \
          2>> {log}
        """

# ----------------------------------------
# Rule: merge_bams
# Merge lane-level BAMs into one per ID
#-----------------------------------------
RG_BAMS = expand(os.path.join(BAMDIR, "{base}.rg.bam"), base=list(FASTQ_DICT.keys()))
rule merge_bams:
    input:
        rg_bams=RG_BAMS
    output:
        merged=os.path.join(BAMDIR, f"{ID}.rg.merged.bam")
    log:
        os.path.join(LOGDIR, "merge.log")
    params:
        merge_inputs=lambda wildcards, input: " ".join(f"-I {b}" for b in input.rg_bams)
    shell:
        r"""
        {GATK4} MergeSamFiles {params.merge_inputs} \
          --OUTPUT {output.merged} --CREATE_INDEX true \
          --VALIDATION_STRINGENCY SILENT --TMP_DIR {TMPDIR} \
          2>> {log}
        """

# ----------------------------------------
# Rule: mark_duplicates
# Mark duplicates and generate metrics
#-----------------------------------------
rule mark_duplicates:
    input:
        merged=os.path.join(BAMDIR, f"{ID}.rg.merged.bam")
    output:
        dedup=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.bam"),
        metrics=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.metrics.txt")
    log:
        os.path.join(LOGDIR, "markdup.log")
    shell:
        r"""
        {GATK4} MarkDuplicates \
          -I {input.merged} \
          -O {output.dedup} \
          --METRICS_FILE {output.metrics} \
          --CREATE_INDEX true --TMP_DIR {TMPDIR} \
          2>> {log}
        """

# ----------------------------------------
# Rule: base_recalibrator
# Build BQSR model using known variant sites
#-----------------------------------------
rule base_recalibrator:
    input:
        dedup=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.bam")
    output:
        table=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.table")
    log:
        os.path.join(LOGDIR, "bqsr.log")
    shell:
        r"""
        {GATK4} BaseRecalibrator \
          -R {REF} \
          -I {input.dedup} \
          --known-sites {dbSNP} --known-sites {MILLS_INDELS} --known-sites {KG_INDELS} \
          -O {output.table} \
          --tmp-dir {TMPDIR} \
          2>> {log}
        """

# ----------------------------------------
# Rule: apply_bqsr
# Apply base quality recalibration and index
#-----------------------------------------
rule apply_bqsr:
    input:
        dedup=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.bam"),
        table=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.table")
    output:
        recal=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.bam")
    log:
        os.path.join(LOGDIR, "apply_bqsr.log")
    shell:
        r"""
        {GATK4} ApplyBQSR \
          -R {REF} \
          -I {input.dedup} \
          --bqsr-recal-file {input.table} \
          -O {output.recal} \
          --tmp-dir {TMPDIR} \
          2>> {log}
        {SAM} index {output.recal}
        """

# ----------------------------------------
# Rule: haplotypecaller
# Generate per-sample gVCF using GATK HaplotypeCaller
#-----------------------------------------
rule haplotypecaller:
    input:
        recal=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.bam")
    output:
        gvcf=os.path.join(VARCALLDIR, f"{ID}.hc.g.vcf.gz")
    threads: THREADS
    params:
        interval_str=lambda wc: " ".join(INTERVAL_ARG)
    log:
        os.path.join(LOGDIR, "hc.log")
    shell:
        r"""
        {GATK4} HaplotypeCaller \
          -R {REF} \
          -I {input.recal} \
          {params.interval_str} \
          --native-pair-hmm-threads {threads} \
          -ERC GVCF \
          -O {output.gvcf} \
          2>> {log}
        """

# ----------------------------------------
# Rule: genotype_gvcfs
# Jointly genotype the gVCF to produce raw VCF
#-----------------------------------------
rule genotype_gvcfs:
    input:
        gvcf=os.path.join(VARCALLDIR, f"{ID}.hc.g.vcf.gz")
    output:
        raw=os.path.join(VARCALLDIR, f"{ID}.hc.raw.vcf.gz")
    log:
        os.path.join(LOGDIR, "genotype.log")
    shell:
        r"""
        {GATK4} GenotypeGVCFs \
          -R {REF} \
          -V {input.gvcf} \
          --stand-call-conf 10 \
          -O {output.raw} \
          2>> {log}
        """

# ----------------------------------------
# Rule: snp_recal and indel_recal
# Build VQSR models if enough variants exist
#-----------------------------------------
rule snp_recal:
    input:
        raw=os.path.join(VARCALLDIR, f"{ID}.hc.raw.vcf.gz")
    output:
        recal=os.path.join(VARCALLDIR, f"{ID}.snp.recal.vcf.gz"),
        tranches=os.path.join(VARCALLDIR, f"{ID}.snp.tranches.txt")
    log:
        os.path.join(LOGDIR, "snp_recal.log")
    shell:
        r"""
        nSNP=$(zgrep -v '^#' {input.raw} | awk 'length($5)==1' | wc -l)
        if [ $nSNP -ge 1000 ]; then
          {GATK4} VariantRecalibrator \
            -R {REF} -V {input.raw} {snp_res} \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
            --mode SNP \
            -O {output.recal} --tranches-file {output.tranches} \
            --max-gaussians 6 2>> {log}
            #--rscript-file {VARCALLDIR}/{ID}.snp.plots.R
        else
          cp {input.raw} {output.recal}
          echo "SKIP" > {output.tranches}
        fi
        """

rule indel_recal:
    input:
        raw=os.path.join(VARCALLDIR, f"{ID}.hc.raw.vcf.gz")
    output:
        recal=os.path.join(VARCALLDIR, f"{ID}.indel.recal.vcf.gz"),
        tranches=os.path.join(VARCALLDIR, f"{ID}.indel.tranches.txt")
    log:
        os.path.join(LOGDIR, "indel_recal.log")
    shell:
        r"""
        nINDEL=$(zgrep -v '^#' {input.raw} | awk 'length($5)!=1' | wc -l)
        if [ $nINDEL -ge 8000 ]; then
          {GATK4} VariantRecalibrator \
            -R {REF} -V {input.raw} {indel_res} \
            -an QD -an FS -an ReadPosRankSum --mode INDEL \
            -O {output.recal} --tranches-file {output.tranches} \
            --max-gaussians 4 2>> {log}
            #--rscript-file {VARCALLDIR}/{ID}.indel.plots.R
        else
          cp {input.raw} {output.recal}
          echo "SKIP" > {output.tranches}
        fi
        """

# ----------------------------------------
# Rule: apply_vqsr
# Apply recalibration or fallback to hard-filtering
# ----------------------------------------
rule apply_vqsr:
    input:
        raw=os.path.join(VARCALLDIR, f"{ID}.hc.raw.vcf.gz"),
        snp_recal=os.path.join(VARCALLDIR, f"{ID}.snp.recal.vcf.gz"),
        snp_tranches=os.path.join(VARCALLDIR, f"{ID}.snp.tranches.txt"),
        indel_recal=os.path.join(VARCALLDIR, f"{ID}.indel.recal.vcf.gz"),
        indel_tranches=os.path.join(VARCALLDIR, f"{ID}.indel.tranches.txt")
    output:
        vqsr=os.path.join(VARCALLDIR, f"{ID}.hc.vqsr.vcf.gz"),
        vqsr_idx=os.path.join(VARCALLDIR, f"{ID}.hc.vqsr.vcf.gz.tbi")
    log:
        os.path.join(LOGDIR, "apply_vqsr.log")
    shell:
        r"""
        tmp={input.raw}
        # Apply SNP VQSR if not skipped
        if ! zgrep -q SKIP {input.snp_tranches}; then
            {GATK4} ApplyVQSR -R {REF} -V $tmp --recal-file {input.snp_recal} \
                --tranches-file {input.snp_tranches} --mode SNP \
                --truth-sensitivity-filter-level 99.0 \
                -O {VARCALLDIR}/{ID}.hc.post_snp.vcf.gz 2>> {log}
            tmp={VARCALLDIR}/{ID}.hc.post_snp.vcf.gz
        fi
        # Apply INDEL VQSR or copy forward & index if skipped
        if ! zgrep -q SKIP {input.indel_tranches}; then
            {GATK4} ApplyVQSR -R {REF} -V $tmp --recal-file {input.indel_recal} \
                --tranches-file {input.indel_tranches} --mode INDEL \
                --truth-sensitivity-filter-level 95.0 \
                -O {output.vqsr} 2>> {log}
        else
            cp $tmp {output.vqsr}
            cp $tmp.tbi {output.vqsr_idx}
        fi
        """

# ----------------------------------------
# Rule: variant_filtration
# Hard-filter remaining variants for QC VCF
# ----------------------------------------
rule variant_filtration:
    input:
        vqsr=os.path.join(VARCALLDIR, f"{ID}.hc.vqsr.vcf.gz"),
        vqsr_idx=os.path.join(VARCALLDIR, f"{ID}.hc.vqsr.vcf.gz.tbi")
    output:
        qc=os.path.join(VARCALLDIR, f"{ID}.hc.QC.vcf.gz")
    log:
        os.path.join(LOGDIR, "filter.log")
    shell:
        r"""
        {GATK4} VariantFiltration \
          -R {REF} -V {input.vqsr} \
          --filter-name LowQUAL    --filter-expression "QUAL < 30.0" \
          --filter-name QD2        --filter-expression "QD < 2.0" \
          --filter-name FS60       --filter-expression "FS > 60.0" \
          --filter-name MQ40       --filter-expression "MQ < 40.0" \
          --filter-name MQRS-12.5  --filter-expression "MQRankSum < -12.5" \
          --filter-name RPRS-8     --filter-expression "ReadPosRankSum < -8.0" \
          --filter-name QD2_INDEL  --filter-expression "QD < 2.0" \
          --filter-name FS200      --filter-expression "FS > 200.0" \
          --filter-name RPRS-20    --filter-expression "ReadPosRankSum < -20.0" \
          -O {output.qc} 2>> {log}
        """

# ----------------------------------------
# Rule: coverage_stats
# Compute coverage stats for chr1 (or 1) before/after recal
# Temporary BAMs and their index files are cleaned up after use
#-----------------------------------------
rule coverage_stats:
    input:
        raw   = os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.bam"),
        recal = os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.bam")
    output:
        cov   = os.path.join(STATSDIR, f"{ID}.coverage.txt")
    log:
        os.path.join(LOGDIR, "coverage.log")
    params:
        chrN = "1" if "b37" in str(REF) else "chr1"
    shell:
        r"""
        # ensure indexes exist
        {SAM} index {input.raw}   2>> {log} || true
        {SAM} index {input.recal} 2>> {log} || true

        # define scratch paths inside STATSDIR
        RAW_TMP={STATSDIR}/{params.chrN}.tmp_raw.bam
        RECAL_TMP={STATSDIR}/{params.chrN}.tmp_recal.bam

        # slice out chr, index them
        {SAM} view -b {input.raw}   {params.chrN} > $RAW_TMP
        {SAM} view -b {input.recal} {params.chrN} > $RECAL_TMP
        {SAM} index $RAW_TMP
        {SAM} index $RECAL_TMP

        # run coverage
        bash {COV} {ID} $RAW_TMP $RECAL_TMP {PIPELINE} \
          > {output.cov} 2>> {log}

        # clean up both BAM + BAI in STATSDIR
        rm \
          {STATSDIR}/{params.chrN}.tmp_raw.bam \
          {STATSDIR}/{params.chrN}.tmp_raw.bam.bai \
          {STATSDIR}/{params.chrN}.tmp_recal.bam \
          {STATSDIR}/{params.chrN}.tmp_recal.bam.bai
        """

# ----------------------------------------
# Rule: sex_determination
# Determine sample sex from final QC VCF
#-----------------------------------------
rule sex_determination:
    input:
        qc=os.path.join(VARCALLDIR, f"{ID}.hc.QC.vcf.gz")
    output:
        sex=os.path.join(STATSDIR, f"{ID}.sex.txt")
    log:
        os.path.join(LOGDIR, "sex.log")
    shell:
        r"""
        bash {VCF2SEX} {input.qc} > {output.sex} 2>> {log}
        """

