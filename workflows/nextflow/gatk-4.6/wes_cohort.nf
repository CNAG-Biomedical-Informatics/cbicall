#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * wes_cohort.nf
 * Nextflow DSL2 replica of the GATK 4.6 WES/WGS cohort joint-genotyping workflow.
 */

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

if( !params.datadir ) throw new IllegalArgumentException("Missing params.datadir")
if( !params.dbdir ) throw new IllegalArgumentException("Missing params.dbdir")
if( !params.ngsutils ) throw new IllegalArgumentException("Missing params.ngsutils")
if( !params.tmpdir ) throw new IllegalArgumentException("Missing params.tmpdir")
if( !params.gatk4_cmd ) throw new IllegalArgumentException("Missing params.gatk4_cmd")
if( !params.resources ) throw new IllegalArgumentException("Missing params.resources")
if( !params.sample_map ) throw new IllegalArgumentException("Missing params.sample_map")

params.pipeline            = params.pipeline ?: "wes"
params.genome              = params.genome ?: "b37"
params.threads             = params.threads ?: 4
params.min_snp_for_vqsr    = params.min_snp_for_vqsr ?: 1000
params.min_indel_for_vqsr  = params.min_indel_for_vqsr ?: 8000
params.mem                 = params.mem ?: "8G"
params.mem_genotype        = params.mem_genotype ?: "64G"

def PIPELINE = params.pipeline.toString().toLowerCase()
if( !['wes','wgs'].contains(PIPELINE) ) {
    throw new IllegalArgumentException("config[pipeline] must be 'wes' or 'wgs'")
}
def GENOME = params.genome.toString().toLowerCase()
if( !['b37','hg38'].contains(GENOME) ) {
    throw new IllegalArgumentException("config[genome] must be 'b37' or 'hg38'")
}
if( PIPELINE == 'wes' && GENOME != 'b37' ) {
    throw new IllegalArgumentException("Nextflow WES currently supports genome='b37' only")
}

def SAMPLE_MAP = file(params.sample_map)
if( !SAMPLE_MAP.exists() ) {
    throw new IllegalArgumentException("sample_map does not exist: ${SAMPLE_MAP}")
}
def SAMPLE_COUNT = SAMPLE_MAP.readLines().findAll { it.trim() }.size()
def WORKSPACE_NAME = params.workspace ? params.workspace.toString() : "cohort.genomicsdb.${SAMPLE_COUNT}"

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

def GATK4_8G = expandPlaceholders(
    params.gatk4_cmd.toString(),
    lvl2 + [ mem: params.mem.toString() ]
)
def GATK4_64G = expandPlaceholders(
    params.gatk4_cmd.toString(),
    lvl2 + [ mem: params.mem_genotype.toString() ]
)

def RESOURCE_CFG = params.resources[GENOME]
if( !RESOURCE_CFG ) {
    throw new IllegalArgumentException("Missing params.resources.${GENOME}")
}

def BUNDLE       = expandPlaceholders(RESOURCE_CFG.bundle.toString(), lvl2)
def REF          = expandPlaceholders(RESOURCE_CFG.ref.toString(), [bundle: BUNDLE] + lvl2)
def REF_DICT     = expandPlaceholders(RESOURCE_CFG.ref_dict.toString(), [bundle: BUNDLE] + lvl2)
def DBSNP        = expandPlaceholders(RESOURCE_CFG.dbsnp.toString(), lvl2 + [bundle: BUNDLE])
def MILLS_INDELS = expandPlaceholders(RESOURCE_CFG.mills_indels.toString(), [bundle: BUNDLE] + lvl2)
def HAPMAP       = expandPlaceholders(RESOURCE_CFG.hapmap.toString(), [bundle: BUNDLE] + lvl2)
def OMNI         = expandPlaceholders(RESOURCE_CFG.omni.toString(), [bundle: BUNDLE] + lvl2)

def SNP_RES = expandPlaceholders(
    RESOURCE_CFG.snp_res.toString(),
    [hapmap: HAPMAP, omni: OMNI, dbsnp: DBSNP] + lvl2 + [bundle: BUNDLE]
)
def INDEL_RES = expandPlaceholders(
    RESOURCE_CFG.indel_res.toString(),
    [mills_indels: MILLS_INDELS] + lvl2 + [bundle: BUNDLE]
)
def INTERVAL_LIST = expandPlaceholders(
    RESOURCE_CFG.interval_list.toString(),
    [bundle: BUNDLE] + lvl2
)
def VCF2HASH = params.vcf2hash_script ? params.vcf2hash_script.toString() : file("${projectDir}/vcf2hash.sh").toString()

def GENOMICSDBDIR = "01_genomicsdb"
def VARCALLDIR = "02_varcall"
def STATSDIR   = "03_stats"
def LOGDIR     = "logs"
new File(GENOMICSDBDIR).mkdirs()
new File(VARCALLDIR).mkdirs()
new File(STATSDIR).mkdirs()
new File(LOGDIR).mkdirs()

def writeWgsIntervalList(String refDict, String outPath) {
    def intervals = []
    def out = new File(outPath)
    out.parentFile.mkdirs()
    out.withWriter('UTF-8') { writer ->
        new File(refDict).eachLine('UTF-8') { line ->
            if (line.startsWith('@')) {
                writer.writeLine(line)
                if (line.startsWith('@SQ')) {
                    def fields = [:]
                    line.split('\t').drop(1).each { item ->
                        def parts = item.split(':', 2)
                        if (parts.size() == 2) fields[parts[0]] = parts[1]
                    }
                    if (fields.SN && fields.LN) {
                        intervals << "${fields.SN}\t1\t${fields.LN}\t+\t${fields.SN}"
                    }
                }
            }
        }
        intervals.each { writer.writeLine(it) }
    }
}

def INTERVAL_ARG
def MERGE_INTERVALS_ARG
if (PIPELINE == 'wes') {
    INTERVAL_ARG = "-L ${q(INTERVAL_LIST)}"
    MERGE_INTERVALS_ARG = "--merge-input-intervals true"
} else {
    def WGS_INTERVAL_LIST = new File("${GENOMICSDBDIR}/wgs.whole_genome.interval_list").absolutePath
    writeWgsIntervalList(REF_DICT, WGS_INTERVAL_LIST)
    INTERVAL_ARG = "-L ${q(WGS_INTERVAL_LIST)}"
    MERGE_INTERVALS_ARG = ""
}

def ENV_BLOCK = """
export TMPDIR=${q(TMPDIR)}
export LC_ALL=C
export GATK_DISABLE_AUTO_S3_UPLOAD=true
"""

println "Pipeline: ${PIPELINE}"
println "Genome: ${GENOME}"
println "Sample map: ${SAMPLE_MAP}"
println "Sample count: ${SAMPLE_COUNT}"
println "Workspace: ${WORKSPACE_NAME}"
println "GATK4: ${GATK4_8G}"
println "REF: ${REF}"

process GENOMICSDB_IMPORT {
    cpus { params.threads as int }

    publishDir GENOMICSDBDIR, mode: 'copy', pattern: "genomicsdbimport.done"
    publishDir GENOMICSDBDIR, mode: 'copy', pattern: "${WORKSPACE_NAME}"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(sample_map)

    output:
    path("${WORKSPACE_NAME}")
    path("genomicsdbimport.done")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    ${GATK4_64G} GenomicsDBImport \\
      --sample-name-map ${q(sample_map)} \\
      --genomicsdb-workspace-path ${q(WORKSPACE_NAME)} \\
      ${MERGE_INTERVALS_ARG} \\
      ${INTERVAL_ARG} \\
      --tmp-dir ${q(TMPDIR)} \\
      2>> ${q("01_genomicsdbimport.log")}

    echo ok > genomicsdbimport.done
    """
}

process GENOTYPE_GVCFS_COHORT {
    cpus { params.threads as int }

    publishDir VARCALLDIR, mode: 'copy', pattern: "cohort.gv.raw.vcf.gz*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(workspace_dir)
    path(db_done)

    output:
    path("cohort.gv.raw.vcf.gz")
    path("cohort.gv.raw.vcf.gz.tbi")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    ${GATK4_64G} GenotypeGVCFs \\
      -R ${q(REF)} \\
      -V ${q("gendb://${workspace_dir}")} \\
      -O cohort.gv.raw.vcf.gz \\
      --stand-call-conf 10 \\
      --tmp-dir ${q(TMPDIR)} \\
      ${INTERVAL_ARG} \\
      2>> ${q("02_genotype_gvcfs.log")}
    """
}

process VQSR_AND_QC_COHORT {
    publishDir VARCALLDIR, mode: 'copy', pattern: "cohort.*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(raw_vcf)
    path(raw_tbi)

    output:
    path("cohort.gv.QC.vcf.gz")

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
      ${GATK4_8G} VariantRecalibrator \\
        -R ${q(REF)} \\
        -V "\$rawvcf" \\
        ${SNP_RES} \\
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \\
        --mode SNP \\
        -O cohort.snp.recal.vcf.gz \\
        --tranches-file cohort.snp.tranches.txt \\
        --max-gaussians 6 \\
        --tmp-dir ${q(TMPDIR)} \\
        2>> ${q("03_vqsr_and_qc.log")}
      apply_snp=true
    fi

    if [ "\$nINDEL" -ge "${params.min_indel_for_vqsr}" ]; then
      ${GATK4_8G} VariantRecalibrator \\
        -R ${q(REF)} \\
        -V "\$rawvcf" \\
        ${INDEL_RES} \\
        -an QD -an FS -an ReadPosRankSum \\
        --mode INDEL \\
        -O cohort.indel.recal.vcf.gz \\
        --tranches-file cohort.indel.tranches.txt \\
        --max-gaussians 4 \\
        --tmp-dir ${q(TMPDIR)} \\
        2>> ${q("03_vqsr_and_qc.log")}
      apply_indel=true
    fi

    if [ "\$apply_snp" = true ]; then
      ${GATK4_8G} ApplyVQSR \\
        -R ${q(REF)} -V "\$tmp_vcf" \\
        --recal-file cohort.snp.recal.vcf.gz \\
        --tranches-file cohort.snp.tranches.txt \\
        --mode SNP --truth-sensitivity-filter-level 99.0 \\
        -O cohort.post_snp.vcf.gz \\
        --tmp-dir ${q(TMPDIR)} \\
        2>> ${q("03_vqsr_and_qc.log")}
      tmp_vcf=cohort.post_snp.vcf.gz
    fi

    if [ "\$apply_indel" = true ]; then
      ${GATK4_8G} ApplyVQSR \\
        -R ${q(REF)} -V "\$tmp_vcf" \\
        --recal-file cohort.indel.recal.vcf.gz \\
        --tranches-file cohort.indel.tranches.txt \\
        --mode INDEL --truth-sensitivity-filter-level 95.0 \\
        -O cohort.vqsr.vcf.gz \\
        --tmp-dir ${q(TMPDIR)} \\
        2>> ${q("03_vqsr_and_qc.log")}
      tmp_vcf=cohort.vqsr.vcf.gz
    fi

    ${GATK4_8G} VariantFiltration \\
      -R ${q(REF)} \\
      -V "\$tmp_vcf" \\
      --filter-name "LowQUAL" --filter-expression "QUAL < 30.0" \\
      --filter-name "QD2"        --filter-expression "QD < 2.0" \\
      --filter-name "FS60"       --filter-expression "FS > 60.0" \\
      --filter-name "MQ40"       --filter-expression "MQ < 40.0" \\
      --filter-name "MQRS-12.5"  --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \\
      --filter-name "RPRS-8"     --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \\
      --filter-name "QD2_indel"  --filter-expression "QD < 2.0" \\
      --filter-name "FS200"      --filter-expression "FS > 200.0" \\
      --filter-name "RPRS-20"    --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \\
      -O cohort.gv.QC.vcf.gz \\
      --tmp-dir ${q(TMPDIR)} \\
      2>> ${q("03_vqsr_and_qc.log")}
    """
}

process VCF_HASH_COHORT {
    publishDir STATSDIR, mode: 'copy', pattern: "cohort.gv.QC.vcf.sha256.txt"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(qc_vcf)

    output:
    path("cohort.gv.QC.vcf.sha256.txt")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    bash ${q(VCF2HASH)} ${q(qc_vcf)} > cohort.gv.QC.vcf.sha256.txt 2>> ${q("04_vcf_hash.log")}
    """
}

workflow {
    sample_map_ch = Channel.value(SAMPLE_MAP)
    db = GENOMICSDB_IMPORT(sample_map_ch)
    rawvcf = GENOTYPE_GVCFS_COHORT(db[0], db[1])
    qcvcf = VQSR_AND_QC_COHORT(rawvcf[0], rawvcf[1])
    VCF_HASH_COHORT(qcvcf[0])
}
