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
params.pipeline            = params.pipeline ?: "wes"
params.genome              = params.genome ?: "b37"
params.threads             = params.threads ?: 4
params.min_snp_for_vqsr    = params.min_snp_for_vqsr ?: 1000
params.min_indel_for_vqsr  = params.min_indel_for_vqsr ?: 8000
params.mem                 = params.mem ?: "8G"
params.mem_genotype        = params.mem_genotype ?: "64G"
params.cohort_stage        = params.cohort_stage ?: "all"
params.output_basename     = params.output_basename ?: "cohort"

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

def COHORT_STAGE = params.cohort_stage.toString()
if( !['all','shard','finalize'].contains(COHORT_STAGE) ) {
    throw new IllegalArgumentException("params.cohort_stage must be 'all', 'shard', or 'finalize'")
}
def OUTPUT_BASENAME = params.output_basename.toString()
def INTERVAL_SHARD = params.interval_shard ? params.interval_shard.toString() : ""
def INPUT_VCF = params.input_vcf ? params.input_vcf.toString() : ""
if( COHORT_STAGE == 'shard' && !INTERVAL_SHARD ) {
    throw new IllegalArgumentException("params.interval_shard is required when cohort_stage='shard'")
}
if( COHORT_STAGE == 'finalize' && !INPUT_VCF ) {
    throw new IllegalArgumentException("params.input_vcf is required when cohort_stage='finalize'")
}
if( COHORT_STAGE == 'finalize' && INTERVAL_SHARD ) {
    throw new IllegalArgumentException("cohort_stage='finalize' does not use interval_shard")
}
if( COHORT_STAGE == 'shard' && INPUT_VCF ) {
    throw new IllegalArgumentException("cohort_stage='shard' does not use input_vcf")
}

def SAMPLE_MAP = null
def SAMPLE_COUNT = 0
if( COHORT_STAGE != 'finalize' ) {
    if( !params.sample_map ) throw new IllegalArgumentException("Missing params.sample_map")
    SAMPLE_MAP = file(params.sample_map)
    if( !SAMPLE_MAP.exists() ) {
        throw new IllegalArgumentException("sample_map does not exist: ${SAMPLE_MAP}")
    }
    SAMPLE_COUNT = SAMPLE_MAP.readLines().findAll { it.trim() }.size()
}
def WORKSPACE_NAME = params.workspace ? params.workspace.toString() : "cohort.genomicsdb.${SAMPLE_COUNT}"
def GENOTYPE_ANNOTATION_EXCLUDE_ARG = SAMPLE_COUNT < 10 ? "-AX InbreedingCoeff" : ""

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

def writeWesShardIntervalList(String sourceList, String shard, String outPath) {
    def found = false
    def out = new File(outPath)
    out.parentFile.mkdirs()
    out.withWriter('UTF-8') { writer ->
        new File(sourceList).eachLine('UTF-8') { line ->
            if (line.startsWith('@')) {
                writer.writeLine(line)
            } else {
                def fields = line.split('\t')
                if (fields && fields[0] == shard) {
                    writer.writeLine(line)
                    found = true
                }
            }
        }
    }
    if (!found) {
        throw new IllegalArgumentException("No intervals found for shard ${shard} in ${sourceList}")
    }
}

def writeWgsShardIntervalList(String refDict, String shard, String outPath) {
    def interval = null
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
                    if (fields.SN == shard && fields.LN) {
                        interval = "${fields.SN}\t1\t${fields.LN}\t+\t${fields.SN}"
                    }
                }
            }
        }
        if (interval) writer.writeLine(interval)
    }
    if (!interval) {
        throw new IllegalArgumentException("No reference contig found for shard ${shard} in ${refDict}")
    }
}

def INTERVAL_ARG
def MERGE_INTERVALS_ARG
if (COHORT_STAGE == 'finalize') {
    INTERVAL_ARG = ""
    MERGE_INTERVALS_ARG = ""
} else if (PIPELINE == 'wes') {
    if (INTERVAL_SHARD) {
        def WES_SHARD_INTERVAL_LIST = new File("${GENOMICSDBDIR}/wes.${INTERVAL_SHARD}.interval_list").absolutePath
        writeWesShardIntervalList(INTERVAL_LIST, INTERVAL_SHARD, WES_SHARD_INTERVAL_LIST)
        INTERVAL_ARG = "-L ${q(WES_SHARD_INTERVAL_LIST)}"
    } else {
        INTERVAL_ARG = "-L ${q(INTERVAL_LIST)}"
    }
    MERGE_INTERVALS_ARG = "--merge-input-intervals true"
} else {
    def WGS_INTERVAL_LIST
    if (INTERVAL_SHARD) {
        WGS_INTERVAL_LIST = new File("${GENOMICSDBDIR}/wgs.${INTERVAL_SHARD}.interval_list").absolutePath
        writeWgsShardIntervalList(REF_DICT, INTERVAL_SHARD, WGS_INTERVAL_LIST)
    } else {
        WGS_INTERVAL_LIST = new File("${GENOMICSDBDIR}/wgs.whole_genome.interval_list").absolutePath
        writeWgsIntervalList(REF_DICT, WGS_INTERVAL_LIST)
    }
    INTERVAL_ARG = "-L ${q(WGS_INTERVAL_LIST)}"
    MERGE_INTERVALS_ARG = ""
}

def COHORT_RAW_VCF_NAME = "${OUTPUT_BASENAME}.gv.raw.vcf.gz"
def COHORT_QC_VCF_NAME = "${OUTPUT_BASENAME}.gv.QC.vcf.gz"
def COHORT_HASH_NAME = "${OUTPUT_BASENAME}.gv.QC.vcf.sha256.txt"

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
println "Cohort stage: ${COHORT_STAGE}"
println "Output basename: ${OUTPUT_BASENAME}"
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

    publishDir VARCALLDIR, mode: 'copy', pattern: "${COHORT_RAW_VCF_NAME}*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(workspace_dir)
    path(db_done)

    output:
    path("${COHORT_RAW_VCF_NAME}")
    path("${COHORT_RAW_VCF_NAME}.tbi")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    ${GATK4_64G} GenotypeGVCFs \\
      -R ${q(REF)} \\
      -V ${q("gendb://${workspace_dir}")} \\
      -O ${q(COHORT_RAW_VCF_NAME)} \\
      --stand-call-conf 10 \\
      ${GENOTYPE_ANNOTATION_EXCLUDE_ARG} \\
      --tmp-dir ${q(TMPDIR)} \\
      ${INTERVAL_ARG} \\
      2>> ${q("02_genotype_gvcfs.log")}
    """
}

process VQSR_AND_QC_COHORT {
    publishDir VARCALLDIR, mode: 'copy', pattern: "${OUTPUT_BASENAME}.*"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(raw_vcf)
    path(raw_tbi)

    output:
    path("${COHORT_QC_VCF_NAME}")

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
        -O ${q("${OUTPUT_BASENAME}.snp.recal.vcf.gz")} \\
        --tranches-file ${q("${OUTPUT_BASENAME}.snp.tranches.txt")} \\
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
        -O ${q("${OUTPUT_BASENAME}.indel.recal.vcf.gz")} \\
        --tranches-file ${q("${OUTPUT_BASENAME}.indel.tranches.txt")} \\
        --max-gaussians 4 \\
        --tmp-dir ${q(TMPDIR)} \\
        2>> ${q("03_vqsr_and_qc.log")}
      apply_indel=true
    fi

    if [ "\$apply_snp" = true ]; then
      ${GATK4_8G} ApplyVQSR \\
        -R ${q(REF)} -V "\$tmp_vcf" \\
        --recal-file ${q("${OUTPUT_BASENAME}.snp.recal.vcf.gz")} \\
        --tranches-file ${q("${OUTPUT_BASENAME}.snp.tranches.txt")} \\
        --mode SNP --truth-sensitivity-filter-level 99.0 \\
        -O ${q("${OUTPUT_BASENAME}.post_snp.vcf.gz")} \\
        --tmp-dir ${q(TMPDIR)} \\
        2>> ${q("03_vqsr_and_qc.log")}
      tmp_vcf=${q("${OUTPUT_BASENAME}.post_snp.vcf.gz")}
    fi

    if [ "\$apply_indel" = true ]; then
      ${GATK4_8G} ApplyVQSR \\
        -R ${q(REF)} -V "\$tmp_vcf" \\
        --recal-file ${q("${OUTPUT_BASENAME}.indel.recal.vcf.gz")} \\
        --tranches-file ${q("${OUTPUT_BASENAME}.indel.tranches.txt")} \\
        --mode INDEL --truth-sensitivity-filter-level 95.0 \\
        -O ${q("${OUTPUT_BASENAME}.vqsr.vcf.gz")} \\
        --tmp-dir ${q(TMPDIR)} \\
        2>> ${q("03_vqsr_and_qc.log")}
      tmp_vcf=${q("${OUTPUT_BASENAME}.vqsr.vcf.gz")}
    fi

    ${GATK4_8G} VariantFiltration \\
      -R ${q(REF)} \\
      -V "\$tmp_vcf" \\
      --filter-name "LowQUAL" --filter-expression "QUAL < 30.0" \\
      --filter-name "QD2"        --filter-expression "vc.hasAttribute('QD') && QD < 2.0" \\
      --filter-name "FS60"       --filter-expression "FS > 60.0" \\
      --filter-name "MQ40"       --filter-expression "MQ < 40.0" \\
      --filter-name "MQRS-12.5"  --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \\
      --filter-name "RPRS-8"     --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \\
      --filter-name "QD2_indel"  --filter-expression "vc.hasAttribute('QD') && QD < 2.0" \\
      --filter-name "FS200"      --filter-expression "FS > 200.0" \\
      --filter-name "RPRS-20"    --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \\
      -O ${q(COHORT_QC_VCF_NAME)} \\
      --tmp-dir ${q(TMPDIR)} \\
      2>> ${q("03_vqsr_and_qc.log")}
    """
}

process VCF_HASH_COHORT {
    publishDir STATSDIR, mode: 'copy', pattern: "${COHORT_HASH_NAME}"
    publishDir LOGDIR, mode: 'copy', pattern: '*.log'

    input:
    path(qc_vcf)

    output:
    path("${COHORT_HASH_NAME}")

    script:
    """
    set -eu
    ${ENV_BLOCK}

    bash ${q(VCF2HASH)} ${q(qc_vcf)} > ${q(COHORT_HASH_NAME)} 2>> ${q("04_vcf_hash.log")}
    """
}

workflow {
    if (COHORT_STAGE == 'finalize') {
        input_vcf = file(INPUT_VCF)
        input_tbi = file("${INPUT_VCF}.tbi")
        if (!input_vcf.exists()) {
            throw new IllegalArgumentException("input_vcf does not exist: ${INPUT_VCF}")
        }
        if (!input_tbi.exists()) {
            throw new IllegalArgumentException("input_vcf index does not exist: ${INPUT_VCF}.tbi")
        }
        qcvcf = VQSR_AND_QC_COHORT(Channel.value(input_vcf), Channel.value(input_tbi))
        VCF_HASH_COHORT(qcvcf[0])
    } else {
        sample_map_ch = Channel.value(SAMPLE_MAP)
        db = GENOMICSDB_IMPORT(sample_map_ch)
        rawvcf = GENOTYPE_GVCFS_COHORT(db[0], db[1])
        if (COHORT_STAGE != 'shard') {
            qcvcf = VQSR_AND_QC_COHORT(rawvcf[0], rawvcf[1])
            VCF_HASH_COHORT(qcvcf[0])
        }
    }
}
