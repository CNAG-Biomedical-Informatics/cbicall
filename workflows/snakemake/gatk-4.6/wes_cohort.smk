# w[eg]s_cohort.smk 
import os
import platform
import shlex
from pathlib import Path

snakefile_dir = Path(workflow.snakefile).parent
configfile: str(snakefile_dir / "config.yaml")

# Environment like parameters.sh
os.environ["TMPDIR"] = config["tmpdir"]
os.environ["LC_ALL"] = "C"
os.environ["GATK_DISABLE_AUTO_S3_UPLOAD"] = "true"

DATADIR  = config["datadir"]
DBDIR    = config["dbdir"].format(datadir=DATADIR)
NGSUTILS = config["ngsutils"].format(datadir=DATADIR)
TMPDIR   = config["tmpdir"].format(datadir=DATADIR)

PIPELINE = config.get("pipeline", "wes").lower()
if PIPELINE not in ("wes", "wgs"):
    raise ValueError("config[pipeline] must be 'wes' or 'wgs'")

THREADS = workflow.cores or int(config.get("threads", 4))

# Same GATK command style as your wes_single config
MEM          = config.get("mem", "8G")
MEM_GENOTYPE = config.get("mem_genotype", "64G")
GATK4_8G  = config["gatk4_cmd"].format(ngsutils=NGSUTILS, mem=MEM)
GATK4_64G = config["gatk4_cmd"].format(ngsutils=NGSUTILS, mem=MEM_GENOTYPE)

# Reference/resources (same keys as your config)
bundle        = config["bundle"].format(dbdir=DBDIR)
REF           = config["ref"].format(bundle=bundle)
INTERVAL_LIST = config["interval_list"].format(bundle=bundle)

dbSNP        = config["dbsnp"].format(dbdir=DBDIR)
MILLS_INDELS = config["mills_indels"].format(bundle=bundle)
HAPMAP       = config["hapmap"].format(bundle=bundle)
OMNI         = config["omni"].format(bundle=bundle)

SNP_RES   = config["snp_res"].format(hapmap=HAPMAP, omni=OMNI, dbsnp=dbSNP)
INDEL_RES = config["indel_res"].format(mills_indels=MILLS_INDELS)

MIN_SNP_FOR_VQSR   = int(config.get("min_snp_for_vqsr", 1000))
MIN_INDEL_FOR_VQSR = int(config.get("min_indel_for_vqsr", 8000))

# sample map is required for cohort workflow
SAMPLE_MAP = config.get("sample_map", "").strip()
if not SAMPLE_MAP:
    raise ValueError("config[sample_map] must be set (or passed via --config sample_map=...)")

# Output dirs match bash: 02_varcall; logs
LOGDIR = "logs"
VARCALLDIR = "02_varcall"
os.makedirs(LOGDIR, exist_ok=True)
os.makedirs(VARCALLDIR, exist_ok=True)

INTERVAL_ARG = f"-L {shlex.quote(INTERVAL_LIST)}" if PIPELINE == "wes" else ""

def count_samples(path):
    with open(path, "r") as f:
        return sum(1 for line in f if line.strip())

SAMPLE_COUNT = count_samples(SAMPLE_MAP)

# workspace: accept from wrapper; if relative, anchor under 02_varcall (bash does: cd 02_varcall)
WORKSPACE = config.get("workspace", "").strip()

# If not provided, pick a default name
if not WORKSPACE:
    WORKSPACE = f"cohort.genomicsdb.{SAMPLE_COUNT}"

# Wrapper always sends relative -> keep everything under 02_varcall
WORKSPACE = os.path.join(VARCALLDIR, WORKSPACE)

# Outputs (match bash filenames inside 02_varcall)
COHORT_RAW_VCF = os.path.join(VARCALLDIR, "cohort.gv.raw.vcf.gz")
COHORT_QC_VCF  = os.path.join(VARCALLDIR, "cohort.gv.QC.vcf.gz")

COHORT_VQSR_SNP       = os.path.join(VARCALLDIR, "cohort.snp.recal.vcf.gz")
COHORT_SNP_TRANCHES   = os.path.join(VARCALLDIR, "cohort.snp.tranches.txt")
COHORT_VQSR_INDEL     = os.path.join(VARCALLDIR, "cohort.indel.recal.vcf.gz")
COHORT_INDEL_TRANCHES = os.path.join(VARCALLDIR, "cohort.indel.tranches.txt")
COHORT_POST_SNP       = os.path.join(VARCALLDIR, "cohort.post_snp.vcf.gz")
COHORT_POST_VQSR      = os.path.join(VARCALLDIR, "cohort.vqsr.vcf.gz")

rule all:
    input:
        COHORT_QC_VCF

rule genomicsdbimport:
    input:
        sample_map=SAMPLE_MAP
    output:
        done=os.path.join(VARCALLDIR, "genomicsdbimport.done")
    params:
        workspace=WORKSPACE
    threads: THREADS
    log:
        os.path.join(LOGDIR, "01_genomicsdbimport.log")
    shell:
        r"""
        set -eu
        set -o pipefail

        mkdir -p "$(dirname {params.workspace})"

        echo ">>> Step 1: GenomicsDBImport" >> {log}
        {GATK4_64G} GenomicsDBImport \
          --sample-name-map {input.sample_map} \
          --genomicsdb-workspace-path {params.workspace} \
          --merge-input-intervals true \
          {INTERVAL_ARG} \
          --tmp-dir {TMPDIR} \
          2>> {log}

        echo "ok" > {output.done}
        """

rule genotype_gvcfs_cohort:
    input:
        db_done=os.path.join(VARCALLDIR, "genomicsdbimport.done")
    output:
        raw=COHORT_RAW_VCF
    params:
        workspace=WORKSPACE
    threads: THREADS
    log:
        os.path.join(LOGDIR, "02_genotype_gvcfs.log")
    shell:
        r"""
        set -eu
        set -o pipefail

        echo ">>> Step 2: GenotypeGVCFs" >> {log}
        {GATK4_64G} GenotypeGVCFs \
          -R {REF} \
          -V "gendb://{params.workspace}" \
          -O {output.raw} \
          --stand-call-conf 10 \
          --tmp-dir {TMPDIR} \
          {INTERVAL_ARG} \
          2>> {log}

        echo "Genotyping completed. Raw cohort VCF: {output.raw}" >> {log}
        """

rule vqsr_and_qc:
    input:
        raw=COHORT_RAW_VCF
    output:
        qc=COHORT_QC_VCF
    log:
        os.path.join(LOGDIR, "03_vqsr_and_qc.log")
    shell:
        r"""
        set -eu
        set -o pipefail

        rawvcf={input.raw}
        tmp_vcf="$rawvcf"

        echo ">>> Step 3: Count variants and decide on VQSR" >> {log}
        nSNP=$(zgrep -v '^#' "$rawvcf" | awk 'length($5)==1'  | wc -l | tr -d ' ')
        nINDEL=$(zgrep -v '^#' "$rawvcf" | awk 'length($5)!=1' | wc -l | tr -d ' ')
        echo "Found SNPs: $nSNP ; INDELs: $nINDEL" >> {log}

        apply_snp=false
        apply_indel=false

        if [ "$nSNP" -ge "{MIN_SNP_FOR_VQSR}" ]; then
          echo ">>> Building SNP VQSR model (VariantRecalibrator)" >> {log}
          {GATK4_8G} VariantRecalibrator \
            -R {REF} \
            -V "$rawvcf" \
            {SNP_RES} \
            -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
            --mode SNP \
            -O {COHORT_VQSR_SNP} \
            --tranches-file {COHORT_SNP_TRANCHES} \
            --max-gaussians 6 \
            --tmp-dir {TMPDIR} \
            2>> {log}
          apply_snp=true
        else
          echo "Skipping SNP VQSR (found $nSNP < min {MIN_SNP_FOR_VQSR})" >> {log}
        fi

        if [ "$nINDEL" -ge "{MIN_INDEL_FOR_VQSR}" ]; then
          echo ">>> Building INDEL VQSR model (VariantRecalibrator)" >> {log}
          {GATK4_8G} VariantRecalibrator \
            -R {REF} \
            -V "$rawvcf" \
            {INDEL_RES} \
            -an QD -an FS -an ReadPosRankSum \
            --mode INDEL \
            -O {COHORT_VQSR_INDEL} \
            --tranches-file {COHORT_INDEL_TRANCHES} \
            --max-gaussians 4 \
            --tmp-dir {TMPDIR} \
            2>> {log}
          apply_indel=true
        else
          echo "Skipping INDEL VQSR (found $nINDEL < min {MIN_INDEL_FOR_VQSR})" >> {log}
        fi

        echo ">>> Step 6: Apply VQSR (if available)" >> {log}

        if [ "$apply_snp" = true ]; then
          echo "Applying SNP recalibration..." >> {log}
          {GATK4_8G} ApplyVQSR \
            -R {REF} -V "$tmp_vcf" \
            --recal-file {COHORT_VQSR_SNP} \
            --tranches-file {COHORT_SNP_TRANCHES} \
            --mode SNP --truth-sensitivity-filter-level 99.0 \
            -O {COHORT_POST_SNP} \
            --tmp-dir {TMPDIR} \
            2>> {log}
          tmp_vcf="{COHORT_POST_SNP}"
        fi

        if [ "$apply_indel" = true ]; then
          echo "Applying INDEL recalibration..." >> {log}
          {GATK4_8G} ApplyVQSR \
            -R {REF} -V "$tmp_vcf" \
            --recal-file {COHORT_VQSR_INDEL} \
            --tranches-file {COHORT_INDEL_TRANCHES} \
            --mode INDEL --truth-sensitivity-filter-level 95.0 \
            -O {COHORT_POST_VQSR} \
            --tmp-dir {TMPDIR} \
            2>> {log}
          tmp_vcf="{COHORT_POST_VQSR}"
        fi

        echo ">>> Step 7: Hard-filter & write cohort QC VCF" >> {log}
        {GATK4_8G} VariantFiltration \
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
          --tmp-dir {TMPDIR} \
          2>> {log}

        echo "Cohort QC VCF written: {output.qc}" >> {log}
        echo "All done. Cohort raw: {COHORT_RAW_VCF} ; Cohort QC: {output.qc}" >> {log}
        """

