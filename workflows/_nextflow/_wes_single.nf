#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.StringEscapeUtils

/*
 * main.nf
 * Nextflow DSL2 replica of wes_single.smk using the same config.yaml keys
 *
 * Run with:
 *   nextflow run main.nf -params-file config.yaml
 */

/****************************************************************
 * Helpers
 ****************************************************************/
def expandPlaceholders(String s, Map vars) {
    if (s == null) return null
    def out = s
    def changed = true
    while (changed) {
        def prev = out
        vars.each { k, v ->
            out = out.replace("{${k}}".toString(), v?.toString())
        }
        changed = (out != prev)
    }
    return out
}

def q(x) {
    return "\"${x.toString().replace("\"", "\\\"")}\""
}

/****************************************************************
 * Required params from config.yaml
 ****************************************************************/
if( !params.datadir ) throw new IllegalArgumentException("Missing params.datadir")
if( !params.dbdir ) throw new IllegalArgumentException("Missing params.dbdir")
if( !params.ngsutils ) throw new IllegalArgumentException("Missing params.ngsutils")
if( !params.tmpdir ) throw new IllegalArgumentException("Missing params.tmpdir")
if( !params.java ) throw new IllegalArgumentException("Missing params.java")
if( !params.tools ) throw new IllegalArgumentException("Missing params.tools")
if( !params.gatk4_cmd ) throw new IllegalArgumentException("Missing params.gatk4_cmd")
if( !params.bundle ) throw new IllegalArgumentException("Missing params.bundle")
if( !params.ref ) throw new IllegalArgumentException("Missing params.ref")
if( !params.refgz ) throw new IllegalArgumentException("Missing params.refgz")
if( !params.dbsnp ) throw new IllegalArgumentException("Missing params.dbsnp")
if( !params.mills_indels ) throw new IllegalArgumentException("Missing params.mills_indels")
if( !params.kg_indels ) throw new IllegalArgumentException("Missing params.kg_indels")
if( !params.hapmap ) throw new IllegalArgumentException("Missing params.hapmap")
if( !params.omni ) throw new IllegalArgumentException("Missing params.omni")
if( !params.snp_res ) throw new IllegalArgumentException("Missing params.snp_res")
if( !params.indel_res ) throw new IllegalArgumentException("Missing params.indel_res")
if( !params.interval_list ) throw new IllegalArgumentException("Missing params.interval_list")

/****************************************************************
 * Defaults matching config.yaml / Snakefile
 ****************************************************************/
params.pipeline            = params.pipeline ?: "wes"
params.threads             = params.threads ?: 4
params.cleanup_bam         = params.cleanup_bam ?: false
params.min_snp_for_vqsr    = params.min_snp_for_vqsr ?: 1000
params.min_indel_for_vqsr  = params.min_indel_for_vqsr ?: 8000
params.mem                 = params.mem ?: "8G"

def PIPELINE = params.pipeline.toString().toLowerCase()
if( !['wes','wgs'].contains(PIPELINE) ) {
    throw new IllegalArgumentException("config[pipeline] must be 'wes' or 'wgs'")
}

/****************************************************************
 * Expand config placeholders exactly like Snakemake config usage
 ****************************************************************/
def DATADIR = params.datadir.toString()

def lvl1 = [
    datadir : DATADIR,
    mem     : params.mem.toString()
]

def DBDIR    = expandPlaceholders(params.dbdir.toString(), lvl1)
def NGSUTILS = expandPlaceholders(params.ngsutils.toString(), lvl1)
def TMPDIR   = expandPlaceholders(params.tmpdir.toString(), lvl1)

def lvl2 = lvl1 + [
    dbdir    : DBDIR,
    ngsutils : NGSUTILS,
    tmpdir   : TMPDIR
]

def ARCH = System.getProperty('os.arch')
def JAVA = (ARCH == 'aarch64' ? params.java.aarch64 : params.java.amd64).toString()

def BWA = expandPlaceholders(
    (ARCH == 'aarch64' ? params.tools.aarch64.bwa : params.tools.amd64.bwa).toString(),
    lvl2
)

def SAM = expandPlaceholders(
    (ARCH == 'aarch64' ? params.tools.aarch64.samtools : params.tools.amd64.samtools).toString(),
    lvl2
)

def GATK4 = expandPlaceholders(
    params.gatk4_cmd.toString(),
    lvl2 + [ mem: params.mem.toString() ]
)

def BUNDLE       = expandPlaceholders(params.bundle.toString(), lvl2)
def REF          = expandPlaceholders(params.ref.toString(), [bundle: BUNDLE] + lvl2)
def REFGZ        = expandPlaceholders(params.refgz.toString(), [bundle: BUNDLE] + lvl2)
def DBSNP        = expandPlaceholders(params.dbsnp.toString(), lvl2)
def MILLS_INDELS = expandPlaceholders(params.mills_indels.toString(), [bundle: BUNDLE] + lvl2)
def KG_INDELS    = expandPlaceholders(params.kg_indels.toString(), [bundle: BUNDLE] + lvl2)
def HAPMAP       = expandPlaceholders(params.hapmap.toString(), [bundle: BUNDLE] + lvl2)
def OMNI         = expandPlaceholders(params.omni.toString(), [bundle: BUNDLE] + lvl2)

def SNP_RES = expandPlaceholders(
    params.snp_res.toString(),
    [
        hapmap: HAPMAP,
        omni: OMNI,
        dbsnp: DBSNP
    ] + lvl2 + [bundle: BUNDLE]
)

def INDEL_RES = expandPlaceholders(
    params.indel_res.toString(),
    [
        mills_indels: MILLS_INDELS
    ] + lvl2 + [bundle: BUNDLE]
)

def INTERVAL_LIST = expandPlaceholders(
    params.interval_list.toString(),
    [bundle: BUNDLE] + lvl2
)

def INTERVAL_ARG = PIPELINE == 'wes' ? "-L ${q(INTERVAL_LIST)}" : ""

def COV     = file("${projectDir}/coverage.sh").toString()
def VCF2SEX = file("${projectDir}/vcf2sex.sh").toString()

/****************************************************************
 * Output dirs
 ****************************************************************/
def BAMDIR     = "01_bam"
def VARCALLDIR = "02_varcall"
def STATSDIR   = "03_stats"
def LOGDIR     = "logs"

new File(BAMDIR).mkdirs()
new File(VARCALLDIR).mkdirs()
new File(STATSDIR).mkdirs()
new File(LOGDIR).mkdirs()

/****************************************************************
 * Sample ID (match Snakemake)
 * rawid = Path.cwd().parent.name
 * ID = rawid.split("_", 1)[0]
 ****************************************************************/
def launchDirFile = new File(workflow.launchDir.toString())
def rawid = launchDirFile.getParentFile()?.getName() ?: launchDirFile.getName()
def ID = rawid.split("_", 2)[0]

def CHR1 = REF.contains("b37") ? "1" : "chr1"

println "Sample ID: ${ID}"
println "Pipeline: ${PIPELINE}"
println "ARCH: ${ARCH}"
println "BWA: ${BWA}"
println "SAM: ${SAM}"
println "GATK4: ${GATK4}"
println "REF: ${REF}"

/****************************************************************
 * FASTQ pairs (match Snakemake glob)
 * ../*_R1_*fastq.gz
 ****************************************************************/
Channel
    .fromPath('../*_R1_*fastq.gz', checkIfExists: true)
    .map { r1 ->
        def r1s = r1.toString()
        def r2s = r1s.replace('_R1_', '_R2_')
        def base = r1.getName().replace('.fastq.gz', '').split('_R1_', 2)[0]
        tuple(base, file(r1s), file(r2s))
    }
    .set { fastq_pairs }

/****************************************************************
 * Shared shell environment
 ****************************************************************/
def ENV_BLOCK = """
export TMPDIR=${q(TMPDIR)}
export LC_ALL=C
export GATK_DISABLE_AUTO_S3_UPLOAD=true
"""

/****************************************************************
 * Processes
 ****************************************************************/

process ALIGN_RG {
    tag { base }
    cpus { params.threads as int }

    publishDir BAMDIR, mode: 'copy', pattern: '*.rg.bam'
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    tuple val(base), path(r1), path(r2)

    output:
    path("${base}.rg.bam")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    SAMPLE=\$(echo ${q(base)} | cut -d'_' -f1-2)
    LANE=\$(echo ${q(base)}   | cut -d'_' -f3)
    RGID="\${SAMPLE}.\${LANE}.\$(date +%s)"
    RGPU="\${SAMPLE}.\${LANE}.unit1"

    ${BWA} mem -M -t ${task.cpus} ${q(REFGZ)} ${q(r1)} ${q(r2)} \\
      | ${GATK4} AddOrReplaceReadGroups \\
          --INPUT /dev/stdin \\
          --OUTPUT ${q("${base}.rg.bam")} \\
          --TMP_DIR ${q(TMPDIR)} \\
          --RGPL ILLUMINA \\
          --RGLB sureselect \\
          --RGSM "\$SAMPLE" \\
          --RGID "\$RGID" \\
          --RGPU "\$RGPU" \\
      2>> ${q("${ID}.01_align_rg.${base}.log")}
    """
}

process MERGE_BAMS {
    publishDir BAMDIR, mode: 'copy', pattern: "${ID}.rg.merged.bam*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(rg_bams)

    output:
    path("${ID}.rg.merged.bam")
    path("${ID}.rg.merged.bai")

    script:
    def merge_inputs = rg_bams.collect { "-I ${q(it)}" }.join(' ')
    """
    set -eu
    ${ENV_BLOCK}

    ${GATK4} MergeSamFiles \\
      ${merge_inputs} \\
      -O ${q("${ID}.rg.merged.bam")} \\
      --CREATE_INDEX true \\
      --VALIDATION_STRINGENCY SILENT \\
      --TMP_DIR ${q(TMPDIR)} \\
      2>> ${q("${ID}.02_merge_bams.log")}
    """
}

process MARK_DUPLICATES {
    publishDir BAMDIR, mode: 'copy', pattern: "${ID}.rg.merged.dedup*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(merged_bam)

    output:
    path("${ID}.rg.merged.dedup.bam")
    path("${ID}.rg.merged.dedup.bai")
    path("${ID}.rg.merged.dedup.metrics.txt")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    ${GATK4} MarkDuplicates \\
      -I ${q(merged_bam)} \\
      -O ${q("${ID}.rg.merged.dedup.bam")} \\
      --METRICS_FILE ${q("${ID}.rg.merged.dedup.metrics.txt")} \\
      --CREATE_INDEX true \\
      --TMP_DIR ${q(TMPDIR)} \\
      2>> ${q("${ID}.03_mark_duplicates.log")}

    ${SAM} index ${q("${ID}.rg.merged.dedup.bam")} 2>> ${q("${ID}.03_mark_duplicates.log")}
    """
}

process BQSR {
    publishDir BAMDIR, mode: 'copy', pattern: "${ID}.rg.merged.dedup.recal*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(dedup_bam)

    output:
    path("${ID}.rg.merged.dedup.recal.table")
    path("${ID}.rg.merged.dedup.recal.bam")
    path("${ID}.rg.merged.dedup.recal.bai")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    ${GATK4} BaseRecalibrator \\
      -R ${q(REF)} \\
      -I ${q(dedup_bam)} \\
      --known-sites ${q(DBSNP)} \\
      --known-sites ${q(MILLS_INDELS)} \\
      --known-sites ${q(KG_INDELS)} \\
      -O ${q("${ID}.rg.merged.dedup.recal.table")} \\
      --tmp-dir ${q(TMPDIR)} \\
      2>> ${q("${ID}.04_bqsr.log")}

    ${GATK4} ApplyBQSR \\
      -R ${q(REF)} \\
      -I ${q(dedup_bam)} \\
      --bqsr-recal-file ${q("${ID}.rg.merged.dedup.recal.table")} \\
      -O ${q("${ID}.rg.merged.dedup.recal.bam")} \\
      --tmp-dir ${q(TMPDIR)} \\
      2>> ${q("${ID}.04_bqsr.log")}

    ${SAM} index ${q("${ID}.rg.merged.dedup.recal.bam")} 2>> ${q("${ID}.04_bqsr.log")}
    """
}

process HAPLOTYPECALLER {
    cpus { params.threads as int }

    publishDir VARCALLDIR, mode: 'copy', pattern: "${ID}.hc.g.vcf.gz*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(recal_bam)

    output:
    path("${ID}.hc.g.vcf.gz")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    ${GATK4} HaplotypeCaller \\
      -R ${q(REF)} \\
      -I ${q(recal_bam)} \\
      -O ${q("${ID}.hc.g.vcf.gz")} \\
      ${INTERVAL_ARG} \\
      --native-pair-hmm-threads ${task.cpus} \\
      -ERC GVCF \\
      2>> ${q("${ID}.05_haplotypecaller.log")}
    """
}

process GENOTYPE_GVCFS {
    publishDir VARCALLDIR, mode: 'copy', pattern: "${ID}.hc.raw.vcf.gz*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(gvcf)

    output:
    path("${ID}.hc.raw.vcf.gz")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    ${GATK4} GenotypeGVCFs \\
      -R ${q(REF)} \\
      -V ${q(gvcf)} \\
      -O ${q("${ID}.hc.raw.vcf.gz")} \\
      --stand-call-conf 10 \\
      2>> ${q("${ID}.06_genotype_gvcfs.log")}
    """
}

process VQSR_AND_QC {
    publishDir VARCALLDIR, mode: 'copy', pattern: "${ID}.hc.*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(raw_vcf)

    output:
    path("${ID}.hc.QC.vcf.gz")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    rawvcf=${q(raw_vcf)}
    tmp_vcf="\$rawvcf"

    nSNP=\$(zgrep -v '^#' "\$rawvcf" | awk 'length(\$5)==1' | wc -l | tr -d ' ')
    nINDEL=\$(zgrep -v '^#' "\$rawvcf" | awk 'length(\$5)!=1' | wc -l | tr -d ' ')

    apply_snp=false
    apply_indel=false

    if [ "\$nSNP" -ge "${params.min_snp_for_vqsr}" ]; then
      ${GATK4} VariantRecalibrator \\
        -R ${q(REF)} \\
        -V "\$rawvcf" \\
        ${SNP_RES} \\
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \\
        --mode SNP \\
        -O ${q("${ID}.hc.snp.recal.vcf.gz")} \\
        --tranches-file ${q("${ID}.hc.snp.tranches.txt")} \\
        --max-gaussians 6 \\
        2>> ${q("${ID}.07_vqsr_and_qc.log")}
      apply_snp=true
    fi

    if [ "\$nINDEL" -ge "${params.min_indel_for_vqsr}" ]; then
      ${GATK4} VariantRecalibrator \\
        -R ${q(REF)} \\
        -V "\$rawvcf" \\
        ${INDEL_RES} \\
        -an QD -an FS -an ReadPosRankSum \\
        --mode INDEL \\
        -O ${q("${ID}.hc.indel.recal.vcf.gz")} \\
        --tranches-file ${q("${ID}.hc.indel.tranches.txt")} \\
        --max-gaussians 4 \\
        2>> ${q("${ID}.07_vqsr_and_qc.log")}
      apply_indel=true
    fi

    if [ "\$apply_snp" = true ]; then
      ${GATK4} ApplyVQSR \\
        -R ${q(REF)} -V "\$tmp_vcf" \\
        --recal-file ${q("${ID}.hc.snp.recal.vcf.gz")} \\
        --tranches-file ${q("${ID}.hc.snp.tranches.txt")} \\
        --mode SNP --truth-sensitivity-filter-level 99.0 \\
        -O ${q("${ID}.hc.post_snp.vcf.gz")} \\
        2>> ${q("${ID}.07_vqsr_and_qc.log")}
      tmp_vcf=${q("${ID}.hc.post_snp.vcf.gz")}
    fi

    if [ "\$apply_indel" = true ]; then
      ${GATK4} ApplyVQSR \\
        -R ${q(REF)} -V "\$tmp_vcf" \\
        --recal-file ${q("${ID}.hc.indel.recal.vcf.gz")} \\
        --tranches-file ${q("${ID}.hc.indel.tranches.txt")} \\
        --mode INDEL --truth-sensitivity-filter-level 95.0 \\
        -O ${q("${ID}.hc.vqsr.vcf.gz")} \\
        2>> ${q("${ID}.07_vqsr_and_qc.log")}
      tmp_vcf=${q("${ID}.hc.vqsr.vcf.gz")}
    fi

    ${GATK4} VariantFiltration \\
      -R ${q(REF)} \\
      -V "\$tmp_vcf" \\
      --filter-name "LowQUAL" --filter-expression "QUAL < 30.0" \\
      --filter-name "QD2"        --filter-expression "QD < 2.0" \\
      --filter-name "FS60"       --filter-expression "FS > 60.0" \\
      --filter-name "MQ40"       --filter-expression "MQ < 40.0" \\
      --filter-name "MQRS-12.5"  --filter-expression "MQRankSum < -12.5" \\
      --filter-name "RPRS-8"     --filter-expression "ReadPosRankSum < -8.0" \\
      --filter-name "QD2_indel"  --filter-expression "QD < 2.0" \\
      --filter-name "FS200"      --filter-expression "FS > 200.0" \\
      --filter-name "RPRS-20"    --filter-expression "ReadPosRankSum < -20.0" \\
      -O ${q("${ID}.hc.QC.vcf.gz")} \\
      2>> ${q("${ID}.07_vqsr_and_qc.log")}
    """
}

process COVERAGE_STATS {
    publishDir STATSDIR, mode: 'copy', pattern: "${ID}.coverage.txt"
    publishDir LOGDIR,   mode: 'copy', pattern: '*.log'

    input:
    tuple path(raw_bam), path(recal_bam)

    output:
    path("${ID}.coverage.txt")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    chrN=${q(CHR1)}
    bam_raw=${q(raw_bam)}
    bam_recal=${q(recal_bam)}

    out_raw=${q("${CHR1}.raw.bam")}
    out_dedup=${q("${CHR1}.dedup.bam")}

    ${SAM} view -b "\$bam_raw" "\$chrN" > "\$out_raw" 2>> ${q("${ID}.08_coverage_stats.log")}
    ${SAM} view -b "\$bam_recal" "\$chrN" > "\$out_dedup" 2>> ${q("${ID}.08_coverage_stats.log")}

    ${SAM} index "\$out_raw" 2>> ${q("${ID}.08_coverage_stats.log")}
    ${SAM} index "\$out_dedup" 2>> ${q("${ID}.08_coverage_stats.log")}

    bash ${q(COV)} ${q(ID)} "\$out_raw" "\$out_dedup" ${q(PIPELINE)} > ${q("${ID}.coverage.txt")} 2>> ${q("${ID}.08_coverage_stats.log")}

    rm -f "\$out_raw" "\$out_dedup" "\$out_raw.bai" "\$out_dedup.bai"
    """
}

process SEX_DETERMINATION {
    publishDir STATSDIR, mode: 'copy', pattern: "${ID}.sex.txt"
    publishDir LOGDIR,   mode: 'copy', pattern: '*.log'

    input:
    path(qc_vcf)

    output:
    path("${ID}.sex.txt")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    bash ${q(VCF2SEX)} ${q(qc_vcf)} > ${q("${ID}.sex.txt")} 2>> ${q("${ID}.09_sex_determination.log")}
    """
}

process CLEANUP_BAMS {
    publishDir LOGDIR, mode: 'copy', pattern: "${ID}.cleanup.done"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(qc_vcf)

    output:
    path("${ID}.cleanup.done")

    script:
    def cleanup = params.cleanup_bam ? "true" : "false"
    """
    set -eu
    ${ENV_BLOCK}

    if [ "${cleanup}" = "true" ]; then
      find ${q(BAMDIR)} -maxdepth 1 \\( -name "*.bam" -o -name "*.bai" \\) -type f -delete 2>> ${q("${ID}.10_cleanup_bams.log")} || true
    fi

    echo ok > ${q("${ID}.cleanup.done")}
    """
}

/****************************************************************
 * Workflow
 ****************************************************************/
workflow {
    aligned = ALIGN_RG(fastq_pairs)
    merged  = MERGE_BAMS(aligned.collect())
    dedup   = MARK_DUPLICATES(merged[0])
    bqsr    = BQSR(dedup[0])
    gvcf    = HAPLOTYPECALLER(bqsr[1])
    rawvcf  = GENOTYPE_GVCFS(gvcf[0])
    qcvcf   = VQSR_AND_QC(rawvcf[0])

    coverage_input = dedup[0].combine(bqsr[1])
    coverage = COVERAGE_STATS(coverage_input)

    sex = SEX_DETERMINATION(qcvcf[0])

    if (params.cleanup_bam) {
        CLEANUP_BAMS(qcvcf[0])
    }
}
