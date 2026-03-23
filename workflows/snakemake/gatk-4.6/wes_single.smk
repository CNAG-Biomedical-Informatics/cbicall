# wes_single.smk - bash-replica of wes_single.sh (GATK4.6) with per-rule logs
import os
import glob
import platform
import shlex
from pathlib import Path

snakefile_dir = Path(workflow.snakefile).parent
configfile: str(snakefile_dir / "config.yaml")

# -----------------------------
# Environment (match bash)
# -----------------------------
os.environ["TMPDIR"] = config["tmpdir"]
os.environ["LC_ALL"] = "C"
os.environ["GATK_DISABLE_AUTO_S3_UPLOAD"] = "true"

# -----------------------------
# Config / tools
# -----------------------------
DATADIR  = config["datadir"]
DBDIR    = config["dbdir"].format(datadir=DATADIR)
NGSUTILS = config["ngsutils"].format(datadir=DATADIR)
TMPDIR   = config["tmpdir"].format(datadir=DATADIR)
MEM      = config.get("mem", "8G")

PIPELINE = config.get("pipeline", "wes").lower()
if PIPELINE not in ("wes", "wgs"):
    raise ValueError("config[pipeline] must be 'wes' or 'wgs'")

CLEANUP_BAM = bool(config.get("cleanup_bam", False))

THREADS = workflow.cores or int(config.get("threads", 4))

ARCH = platform.machine()
if ARCH == "aarch64":
    JAVA = config["java"]["aarch64"]
    BWA  = config["tools"]["aarch64"]["bwa"].format(ngsutils=NGSUTILS)
    SAM  = config["tools"]["aarch64"]["samtools"].format(ngsutils=NGSUTILS)
else:
    JAVA = config["java"]["amd64"]
    BWA  = config["tools"]["amd64"]["bwa"].format(ngsutils=NGSUTILS)
    SAM  = config["tools"]["amd64"]["samtools"].format(ngsutils=NGSUTILS)

GATK4 = config["gatk4_cmd"].format(ngsutils=NGSUTILS, mem=MEM)

COV     = str(snakefile_dir / "coverage.sh")
VCF2SEX = str(snakefile_dir / "vcf2sex.sh")

bundle       = config["bundle"].format(dbdir=DBDIR)
REF          = config["ref"].format(bundle=bundle)
REFGZ        = config["refgz"].format(bundle=bundle)
dbSNP        = config["dbsnp"].format(dbdir=DBDIR)
MILLS_INDELS = config["mills_indels"].format(bundle=bundle)
KG_INDELS    = config["kg_indels"].format(bundle=bundle)
HAPMAP       = config["hapmap"].format(bundle=bundle)
OMNI         = config["omni"].format(bundle=bundle)

SNP_RES   = config["snp_res"].format(hapmap=HAPMAP, omni=OMNI, dbsnp=dbSNP)
INDEL_RES = config["indel_res"].format(mills_indels=MILLS_INDELS)

INTERVAL_LIST = config["interval_list"].format(bundle=bundle)
INTERVAL_ARG  = f"-L {shlex.quote(INTERVAL_LIST)}" if PIPELINE == "wes" else ""

MIN_SNP_FOR_VQSR   = int(config.get("min_snp_for_vqsr", 1000))
MIN_INDEL_FOR_VQSR = int(config.get("min_indel_for_vqsr", 8000))

# Output dirs
BAMDIR     = "01_bam"
VARCALLDIR = "02_varcall"
STATSDIR   = "03_stats"
LOGDIR     = "logs"
for d in (BAMDIR, VARCALLDIR, STATSDIR, LOGDIR):
    os.makedirs(d, exist_ok=True)

# -----------------------------
# Sample ID (match bash)
# -----------------------------
rawid = Path.cwd().parent.name
ID = rawid.split("_", 1)[0]

# -----------------------------
# FASTQ pairs (match bash glob)
# -----------------------------
FASTQ_DIR = "../"
FASTQ_R1 = sorted(glob.glob(os.path.join(FASTQ_DIR, "*_R1_*fastq.gz")))
if not FASTQ_R1:
    raise ValueError("No FASTQs found matching ../*_R1_*fastq.gz")

FASTQ_BASES = []
FASTQ_DICT = {}
for r1 in FASTQ_R1:
    r2 = r1.replace("_R1_", "_R2_")
    base = os.path.basename(r1).replace(".fastq.gz", "")
    base = base.split("_R1_", 1)[0]
    FASTQ_BASES.append(base)
    FASTQ_DICT[base] = {"r1": r1, "r2": r2}

RG_BAMS = [os.path.join(BAMDIR, f"{b}.rg.bam") for b in FASTQ_BASES]

CHR1 = "1" if "b37" in str(REF) else "chr1"

# -----------------------------
# Targets
# -----------------------------
FINAL_INPUTS = [
    os.path.join(VARCALLDIR, f"{ID}.hc.QC.vcf.gz"),
    os.path.join(STATSDIR,   f"{ID}.coverage.txt"),
    os.path.join(STATSDIR,   f"{ID}.sex.txt"),
]
if CLEANUP_BAM:
    FINAL_INPUTS.append(os.path.join(LOGDIR, f"{ID}.cleanup.done"))

rule all:
    input:
        FINAL_INPUTS

# -----------------------------
# STEP 1: Align & AddReadGroups
# -----------------------------
rule align_rg:
    input:
        r1=lambda wc: FASTQ_DICT[wc.base]["r1"],
        r2=lambda wc: FASTQ_DICT[wc.base]["r2"],
    output:
        bam=os.path.join(BAMDIR, "{base}.rg.bam"),
    threads: THREADS
    log:
        os.path.join(LOGDIR, f"{ID}.01_align_rg.{{base}}.log")
    shell:
        r"""
        set -eu
        SAMPLE=$(echo {wildcards.base} | cut -d'_' -f1-2)
        LANE=$(echo   {wildcards.base} | cut -d'_' -f3)
        RGID="${{SAMPLE}}.${{LANE}}.$(date +%s)"
        RGPU="${{SAMPLE}}.${{LANE}}.unit1"

        {BWA} mem -M -t {threads} {REFGZ} {input.r1} {input.r2} \
          | {GATK4} AddOrReplaceReadGroups \
              --INPUT /dev/stdin \
              --OUTPUT {output.bam} \
              --TMP_DIR {TMPDIR} \
              --RGPL ILLUMINA \
              --RGLB sureselect \
              --RGSM "$SAMPLE" \
              --RGID "$RGID" \
              --RGPU "$RGPU" \
          2>> {log}
        """

# -----------------------------
# STEP 2: Merge lane BAMs
# -----------------------------
rule merge_bams:
    input:
        rg_bams=RG_BAMS
    output:
        merged=os.path.join(BAMDIR, f"{ID}.rg.merged.bam")
    params:
        merge_inputs=lambda wc, input: " ".join([f"-I {b}" for b in input.rg_bams])
    log:
        os.path.join(LOGDIR, f"{ID}.02_merge_bams.log")
    shell:
        r"""
        set -eu
        {GATK4} MergeSamFiles \
          {params.merge_inputs} \
          -O {output.merged} \
          --CREATE_INDEX true \
          --VALIDATION_STRINGENCY SILENT \
          --TMP_DIR {TMPDIR} \
          2>> {log}
        """

# -----------------------------
# STEP 3: MarkDuplicates (+ samtools index)
# -----------------------------
rule mark_duplicates:
    input:
        merged=os.path.join(BAMDIR, f"{ID}.rg.merged.bam")
    output:
        dedup=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.bam"),
        metrics=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.metrics.txt")
    log:
        os.path.join(LOGDIR, f"{ID}.03_mark_duplicates.log")
    shell:
        r"""
        set -eu
        {GATK4} MarkDuplicates \
          -I {input.merged} \
          -O {output.dedup} \
          --METRICS_FILE {output.metrics} \
          --CREATE_INDEX true \
          --TMP_DIR {TMPDIR} \
          2>> {log}
        {SAM} index {output.dedup} 2>> {log}
        """

# -----------------------------
# STEP 4: BQSR
# -----------------------------
rule bqsr:
    input:
        dedup=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.bam")
    output:
        table=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.table"),
        recal=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.bam")
    log:
        os.path.join(LOGDIR, f"{ID}.04_bqsr.log")
    shell:
        r"""
        set -eu
        {GATK4} BaseRecalibrator \
          -R {REF} \
          -I {input.dedup} \
          --known-sites {dbSNP} \
          --known-sites {MILLS_INDELS} \
          --known-sites {KG_INDELS} \
          -O {output.table} \
          --tmp-dir {TMPDIR} \
          2>> {log}

        {GATK4} ApplyBQSR \
          -R {REF} \
          -I {input.dedup} \
          --bqsr-recal-file {output.table} \
          -O {output.recal} \
          --tmp-dir {TMPDIR} \
          2>> {log}

        {SAM} index {output.recal} 2>> {log}
        """

# -----------------------------
# STEP 5: HaplotypeCaller -> gVCF
# -----------------------------
rule haplotypecaller:
    input:
        recal=os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.bam")
    output:
        gvcf=os.path.join(VARCALLDIR, f"{ID}.hc.g.vcf.gz")
    threads: THREADS
    log:
        os.path.join(LOGDIR, f"{ID}.05_haplotypecaller.log")
    shell:
        r"""
        set -eu
        {GATK4} HaplotypeCaller \
          -R {REF} \
          -I {input.recal} \
          -O {output.gvcf} \
          {INTERVAL_ARG} \
          --native-pair-hmm-threads {threads} \
          -ERC GVCF \
          2>> {log}
        """

# -----------------------------
# STEP 6: GenotypeGVCFs -> raw VCF
# -----------------------------
rule genotype_gvcfs:
    input:
        gvcf=os.path.join(VARCALLDIR, f"{ID}.hc.g.vcf.gz")
    output:
        raw=os.path.join(VARCALLDIR, f"{ID}.hc.raw.vcf.gz")
    log:
        os.path.join(LOGDIR, f"{ID}.06_genotype_gvcfs.log")
    shell:
        r"""
        set -eu
        {GATK4} GenotypeGVCFs \
          -R {REF} \
          -V {input.gvcf} \
          -O {output.raw} \
          --stand-call-conf 10 \
          2>> {log}
        """

# -----------------------------
# STEPS 7-9: Conditional VQSR + always QC VariantFiltration
# -----------------------------
rule vqsr_and_qc:
    input:
        raw=os.path.join(VARCALLDIR, f"{ID}.hc.raw.vcf.gz")
    output:
        qc=os.path.join(VARCALLDIR, f"{ID}.hc.QC.vcf.gz")
    log:
        os.path.join(LOGDIR, f"{ID}.07_vqsr_and_qc.log")
    shell:
        r"""
        set -eu

        rawvcf={input.raw}
        tmp_vcf="$rawvcf"

        nSNP=$(zgrep -v '^#' "$rawvcf" | awk 'length($5)==1' | wc -l | tr -d ' ')
        nINDEL=$(zgrep -v '^#' "$rawvcf" | awk 'length($5)!=1' | wc -l | tr -d ' ')

        apply_snp=false
        apply_indel=false

        if [ "$nSNP" -ge "{MIN_SNP_FOR_VQSR}" ]; then
          {GATK4} VariantRecalibrator \
            -R {REF} \
            -V "$rawvcf" \
            {SNP_RES} \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
            --mode SNP \
            -O {VARCALLDIR}/{ID}.hc.snp.recal.vcf.gz \
            --tranches-file {VARCALLDIR}/{ID}.hc.snp.tranches.txt \
            --max-gaussians 6 \
            2>> {log}
          apply_snp=true
        fi

        if [ "$nINDEL" -ge "{MIN_INDEL_FOR_VQSR}" ]; then
          {GATK4} VariantRecalibrator \
            -R {REF} \
            -V "$rawvcf" \
            {INDEL_RES} \
            -an QD -an FS -an ReadPosRankSum \
            --mode INDEL \
            -O {VARCALLDIR}/{ID}.hc.indel.recal.vcf.gz \
            --tranches-file {VARCALLDIR}/{ID}.hc.indel.tranches.txt \
            --max-gaussians 4 \
            2>> {log}
          apply_indel=true
        fi

        if [ "$apply_snp" = true ]; then
          {GATK4} ApplyVQSR \
            -R {REF} -V "$tmp_vcf" \
            --recal-file {VARCALLDIR}/{ID}.hc.snp.recal.vcf.gz \
            --tranches-file {VARCALLDIR}/{ID}.hc.snp.tranches.txt \
            --mode SNP --truth-sensitivity-filter-level 99.0 \
            -O {VARCALLDIR}/{ID}.hc.post_snp.vcf.gz \
            2>> {log}
          tmp_vcf="{VARCALLDIR}/{ID}.hc.post_snp.vcf.gz"
        fi

        if [ "$apply_indel" = true ]; then
          {GATK4} ApplyVQSR \
            -R {REF} -V "$tmp_vcf" \
            --recal-file {VARCALLDIR}/{ID}.hc.indel.recal.vcf.gz \
            --tranches-file {VARCALLDIR}/{ID}.hc.indel.tranches.txt \
            --mode INDEL --truth-sensitivity-filter-level 95.0 \
            -O {VARCALLDIR}/{ID}.hc.vqsr.vcf.gz \
            2>> {log}
          tmp_vcf="{VARCALLDIR}/{ID}.hc.vqsr.vcf.gz"
        fi

        {GATK4} VariantFiltration \
          -R {REF} \
          -V "$tmp_vcf" \
          --filter-name "LowQUAL" --filter-expression "QUAL < 30.0" \
          --filter-name "QD2"        --filter-expression "QD < 2.0" \
          --filter-name "FS60"       --filter-expression "FS > 60.0" \
          --filter-name "MQ40"       --filter-expression "MQ < 40.0" \
          --filter-name "MQRS-12.5"  --filter-expression "MQRankSum < -12.5" \
          --filter-name "RPRS-8"     --filter-expression "ReadPosRankSum < -8.0" \
          --filter-name "QD2_indel"  --filter-expression "QD < 2.0" \
          --filter-name "FS200"      --filter-expression "FS > 200.0" \
          --filter-name "RPRS-20"    --filter-expression "ReadPosRankSum < -20.0" \
          -O {output.qc} \
          2>> {log}
        """

# -----------------------------
# STEP 10: Coverage
# -----------------------------
rule coverage_stats:
    input:
        raw   = os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.bam"),
        recal = os.path.join(BAMDIR, f"{ID}.rg.merged.dedup.recal.bam")
    output:
        cov   = os.path.join(STATSDIR, f"{ID}.coverage.txt")
    log:
        os.path.join(LOGDIR, f"{ID}.08_coverage_stats.log")
    shell:
        r"""
        set -eu

        chrN="{CHR1}"
        bam_raw="{input.raw}"
        bam_recal="{input.recal}"

        out_raw="{STATSDIR}/{CHR1}.raw.bam"
        out_dedup="{STATSDIR}/{CHR1}.dedup.bam"

        {SAM} view -b "$bam_raw"   "$chrN" > "$out_raw"   2>> {log}
        {SAM} view -b "$bam_recal" "$chrN" > "$out_dedup" 2>> {log}

        {SAM} index "$out_raw"   2>> {log}
        {SAM} index "$out_dedup" 2>> {log}

        bash {COV} {ID} "$out_raw" "$out_dedup" {PIPELINE} \
          > {output.cov} 2>> {log}

        rm -f "$out_raw" "$out_dedup" "$out_raw.bai" "$out_dedup.bai"
        """

# -----------------------------
# STEP 10 cont: Sex
# -----------------------------
rule sex_determination:
    input:
        qc=os.path.join(VARCALLDIR, f"{ID}.hc.QC.vcf.gz")
    output:
        sex=os.path.join(STATSDIR, f"{ID}.sex.txt")
    log:
        os.path.join(LOGDIR, f"{ID}.09_sex_determination.log")
    shell:
        r"""
        set -eu
        bash {VCF2SEX} {input.qc} > {output.sex} 2>> {log}
        """

# -----------------------------
# STEP 11: Optional cleanup
# -----------------------------
rule cleanup_bams:
    input:
        os.path.join(VARCALLDIR, f"{ID}.hc.QC.vcf.gz")
    output:
        done=os.path.join(LOGDIR, f"{ID}.cleanup.done")
    log:
        os.path.join(LOGDIR, f"{ID}.10_cleanup_bams.log")
    run:
        # write any python errors to the log file
        try:
            if CLEANUP_BAM:
                for pat in ("01_bam/*.bam", "01_bam/*.bai"):
                    for fp in glob.glob(pat):
                        try:
                            os.remove(fp)
                        except FileNotFoundError:
                            pass
            Path(output.done).write_text("ok\n")
        except Exception as e:
            Path(log[0]).write_text(str(e) + "\n")
            raise

