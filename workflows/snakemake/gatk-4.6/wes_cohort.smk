# w[eg]s_cohort.smk 
# Registry version: v1
import os
import platform
import shlex
from pathlib import Path

snakefile_dir = Path(workflow.snakefile).parent
configfile: str(snakefile_dir / "config.yaml")

# Environment like env.sh
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

GENOME = config.get("genome", "b37").lower()
if GENOME not in config.get("resources", {}):
    raise ValueError(f"config[resources] has no entry for genome '{GENOME}'")
if PIPELINE == "wes" and GENOME == "hg38":
    raise ValueError("genome='hg38' is only supported for pipeline='wgs'")

THREADS = workflow.cores or int(config.get("threads", 4))
COHORT_STAGE = str(config.get("cohort_stage", "all")).strip()
if COHORT_STAGE not in ("all", "shard", "finalize"):
    raise ValueError("config[cohort_stage] must be 'all', 'shard', or 'finalize'")
OUTPUT_BASENAME = str(config.get("output_basename", "cohort")).strip()
if not OUTPUT_BASENAME or OUTPUT_BASENAME in (".", ".."):
    raise ValueError("config[output_basename] must be a non-empty filename stem")
INTERVAL_SHARD = str(config.get("interval_shard", "")).strip()
INPUT_VCF = str(config.get("input_vcf", "")).strip()
if COHORT_STAGE == "shard" and not INTERVAL_SHARD:
    raise ValueError("config[interval_shard] must be set when cohort_stage='shard'")
if COHORT_STAGE == "finalize" and not INPUT_VCF:
    raise ValueError("config[input_vcf] must be set when cohort_stage='finalize'")
if COHORT_STAGE == "finalize" and INTERVAL_SHARD:
    raise ValueError("cohort_stage='finalize' does not use interval_shard")
if COHORT_STAGE == "shard" and INPUT_VCF:
    raise ValueError("cohort_stage='shard' does not use input_vcf")

# Same GATK command style as your wes_single config
MEM          = config.get("mem", "8G")
MEM_GENOTYPE = config.get("mem_genotype", "64G")
GATK4_8G  = config["gatk4_cmd"].format(ngsutils=NGSUTILS, mem=MEM)
GATK4_64G = config["gatk4_cmd"].format(ngsutils=NGSUTILS, mem=MEM_GENOTYPE)

resource_cfg  = config["resources"][GENOME]
bundle        = resource_cfg["bundle"].format(dbdir=DBDIR)
REF           = resource_cfg["ref"].format(bundle=bundle)
REF_DICT      = resource_cfg["ref_dict"].format(bundle=bundle)
INTERVAL_LIST = resource_cfg.get("interval_list", "").format(bundle=bundle)

dbSNP        = resource_cfg["dbsnp"].format(dbdir=DBDIR)
MILLS_INDELS = resource_cfg["mills_indels"].format(bundle=bundle)
HAPMAP       = resource_cfg["hapmap"].format(bundle=bundle)
OMNI         = resource_cfg["omni"].format(bundle=bundle)

SNP_RES   = resource_cfg["snp_res"].format(hapmap=HAPMAP, omni=OMNI, dbsnp=dbSNP)
INDEL_RES = resource_cfg["indel_res"].format(mills_indels=MILLS_INDELS)

MIN_SNP_FOR_VQSR   = int(config.get("min_snp_for_vqsr", 1000))
MIN_INDEL_FOR_VQSR = int(config.get("min_indel_for_vqsr", 8000))

SAMPLE_MAP = config.get("sample_map", "").strip()
if COHORT_STAGE != "finalize" and not SAMPLE_MAP:
    raise ValueError("config[sample_map] must be set (or passed via --config sample_map=...)")

# Output dirs match bash: 01_genomicsdb; 02_varcall; logs
LOGDIR = "logs"
GENOMICSDBDIR = "01_genomicsdb"
VARCALLDIR = "02_varcall"
os.makedirs(LOGDIR, exist_ok=True)
os.makedirs(GENOMICSDBDIR, exist_ok=True)
os.makedirs(VARCALLDIR, exist_ok=True)

def write_wes_shard_interval_list(source_list, shard, out_path):
    found = False
    with open(source_list, "r", encoding="utf-8") as src, open(out_path, "w", encoding="utf-8") as out:
        for line in src:
            if line.startswith("@"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] == shard:
                out.write(line)
                found = True
    if not found:
        raise ValueError(f"No intervals found for shard {shard} in {source_list}")


def write_wgs_interval_list(ref_dict, out_path, shard=None):
    intervals = []
    with open(ref_dict, "r", encoding="utf-8") as src, open(out_path, "w", encoding="utf-8") as out:
        for line in src:
            if line.startswith("@"):
                out.write(line)
                if line.startswith("@SQ"):
                    fields = {}
                    for item in line.rstrip("\n").split("\t")[1:]:
                        if ":" in item:
                            key, value = item.split(":", 1)
                            fields[key] = value
                    if fields.get("SN") and fields.get("LN") and (shard is None or fields["SN"] == shard):
                        intervals.append(f"{fields['SN']}\t1\t{fields['LN']}\t+\t{fields['SN']}\n")
        if shard is not None and not intervals:
            raise ValueError(f"No reference contig found for shard {shard} in {ref_dict}")
        out.writelines(intervals)


if COHORT_STAGE == "finalize":
    INTERVAL_ARG = ""
    MERGE_INTERVALS_ARG = ""
elif PIPELINE == "wes":
    if INTERVAL_SHARD:
        WES_SHARD_INTERVAL_LIST = os.path.abspath(os.path.join(GENOMICSDBDIR, f"wes.{INTERVAL_SHARD}.interval_list"))
        write_wes_shard_interval_list(INTERVAL_LIST, INTERVAL_SHARD, WES_SHARD_INTERVAL_LIST)
        INTERVAL_ARG = f"-L {shlex.quote(WES_SHARD_INTERVAL_LIST)}"
    else:
        INTERVAL_ARG = f"-L {shlex.quote(INTERVAL_LIST)}"
    MERGE_INTERVALS_ARG = "--merge-input-intervals true"
else:
    if INTERVAL_SHARD:
        WGS_INTERVAL_LIST = os.path.abspath(os.path.join(GENOMICSDBDIR, f"wgs.{INTERVAL_SHARD}.interval_list"))
        write_wgs_interval_list(REF_DICT, WGS_INTERVAL_LIST, shard=INTERVAL_SHARD)
    else:
        WGS_INTERVAL_LIST = os.path.abspath(os.path.join(GENOMICSDBDIR, "wgs.whole_genome.interval_list"))
        write_wgs_interval_list(REF_DICT, WGS_INTERVAL_LIST)
    INTERVAL_ARG = f"-L {shlex.quote(WGS_INTERVAL_LIST)}"
    MERGE_INTERVALS_ARG = ""

def count_samples(path):
    with open(path, "r") as f:
        return sum(1 for line in f if line.strip())

SAMPLE_COUNT = count_samples(SAMPLE_MAP) if COHORT_STAGE != "finalize" else 0
GENOTYPE_ANNOTATION_EXCLUDE_ARG = "-AX InbreedingCoeff" if SAMPLE_COUNT < 10 else ""

# workspace: accept from wrapper; if relative, anchor under 01_genomicsdb
WORKSPACE = config.get("workspace", "").strip()

# If not provided, pick a default name
if not WORKSPACE:
    WORKSPACE = f"cohort.genomicsdb.{SAMPLE_COUNT}"

# Wrapper always sends relative -> keep the GenomicsDB workspace under 01_genomicsdb
if not os.path.isabs(WORKSPACE):
    WORKSPACE = os.path.join(GENOMICSDBDIR, WORKSPACE)

# Outputs (match bash filenames inside 02_varcall)
COHORT_RAW_VCF = os.path.join(VARCALLDIR, f"{OUTPUT_BASENAME}.gv.raw.vcf.gz")
COHORT_QC_VCF  = os.path.join(VARCALLDIR, f"{OUTPUT_BASENAME}.gv.QC.vcf.gz")

COHORT_VQSR_SNP       = os.path.join(VARCALLDIR, f"{OUTPUT_BASENAME}.snp.recal.vcf.gz")
COHORT_SNP_TRANCHES   = os.path.join(VARCALLDIR, f"{OUTPUT_BASENAME}.snp.tranches.txt")
COHORT_VQSR_INDEL     = os.path.join(VARCALLDIR, f"{OUTPUT_BASENAME}.indel.recal.vcf.gz")
COHORT_INDEL_TRANCHES = os.path.join(VARCALLDIR, f"{OUTPUT_BASENAME}.indel.tranches.txt")
COHORT_POST_SNP       = os.path.join(VARCALLDIR, f"{OUTPUT_BASENAME}.post_snp.vcf.gz")
COHORT_POST_VQSR      = os.path.join(VARCALLDIR, f"{OUTPUT_BASENAME}.vqsr.vcf.gz")
RAW_VCF_FOR_FILTERING = INPUT_VCF if COHORT_STAGE == "finalize" else COHORT_RAW_VCF
FINAL_TARGET = COHORT_RAW_VCF if COHORT_STAGE == "shard" else COHORT_QC_VCF

rule all:
    input:
        FINAL_TARGET

rule genomicsdbimport:
    input:
        sample_map=SAMPLE_MAP
    output:
        done=os.path.join(GENOMICSDBDIR, "genomicsdbimport.done")
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
          {MERGE_INTERVALS_ARG} \
          {INTERVAL_ARG} \
          --tmp-dir {TMPDIR} \
          2>> {log}

        echo "ok" > {output.done}
        """

rule genotype_gvcfs_cohort:
    input:
        db_done=os.path.join(GENOMICSDBDIR, "genomicsdbimport.done")
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
          {GENOTYPE_ANNOTATION_EXCLUDE_ARG} \
          --tmp-dir {TMPDIR} \
          {INTERVAL_ARG} \
          2>> {log}

        echo "Genotyping completed. Raw cohort VCF: {output.raw}" >> {log}
        """

rule vqsr_and_qc:
    input:
        raw=RAW_VCF_FOR_FILTERING
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
          --filter-name "QD2"        --filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
          --filter-name "FS60"       --filter-expression "FS > 60.0" \
          --filter-name "MQ40"       --filter-expression "MQ < 40.0" \
          --filter-name "MQRS-12.5"  --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
          --filter-name "RPRS-8"     --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
          --filter-name "QD2_indel"  --filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
          --filter-name "FS200"      --filter-expression "FS > 200.0" \
          --filter-name "RPRS-20"    --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
          -O {output.qc} \
          --tmp-dir {TMPDIR} \
          2>> {log}

        echo "Cohort QC VCF written: {output.qc}" >> {log}
        echo "All done. Cohort raw: {input.raw} ; Cohort QC: {output.qc}" >> {log}
        """
