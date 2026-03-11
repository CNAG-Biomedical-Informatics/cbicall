#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/********************************************************************
 * Parameters (from your config.yaml)
 ********************************************************************/
params.datadir    = "/media/mrueda/2TBS"
params.dbdir      = "/media/mrueda/2TBS/Databases"
params.ngsutils   = "/media/mrueda/2TBS/NGSutils"
params.tmpdir     = "/media/mrueda/2TBS/tmp"
params.mem        = "4G"

params.java_aarch64 = "/usr/lib/jvm/java-8-openjdk-arm64/bin/java"
params.java_amd64   = "/usr/lib/jvm/java-8-openjdk-amd64/bin/java"

params.tools = [
  aarch64: [
    bwa      : "{ngsutils}/bwa-0.7.18_arm64/bwa",
    samtools : "{ngsutils}/samtools-0.1.19_arm64/samtools",
    bedtools : "{ngsutils}/bedtools2_arm64/bin/bedtools"
  ],
  amd64: [
    bwa      : "{ngsutils}/bwa-0.7.17/bwa",
    samtools : "{ngsutils}/samtools-0.1.19/samtools",
    bedtools : "{ngsutils}/bedtools2/bin/bedtools"
  ]
]

// Commands for Picard and GATK (placeholders will be substituted)
params.picard = "{java} -Xmx{mem} -Djava.io.tmpdir={tmpdir} -jar {ngsutils}/picard/build/libs/picard.jar"
params.gatk   = "{java} -Xmx{mem} -Djava.io.tmpdir={tmpdir} -Dgatk.report.telemetry=false -jar {ngsutils}/gatk/3.5/GenomeAnalysisTK.jar"

// Bundle and reference files
params.bundle      = "${params.dbdir}/GATK_bundle/b37"
params.ref         = "${params.bundle}/references_b37_Homo_sapiens_assembly19.fasta"
params.refgz       = "${params.bundle}/references_b37_Homo_sapiens_assembly19.fasta.gz"
params.dbsnp       = "${params.dbdir}/dbSNP/human_9606_b144_GRCh37p13/All_20160408.vcf.gz"
params.mills_indels= "${params.bundle}/b37_Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
params.kg_indels   = "${params.bundle}/b37_1000G_phase1.indels.b37.vcf.gz"
params.hapmap      = "${params.bundle}/b37_hapmap_3.3.b37.vcf.gz"
params.omni        = "${params.bundle}/b37_1000G_omni2.5.b37.vcf.gz"

// Training sets for variant recalibration
params.snp_res   = """-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
-resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omni} \
-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 ${params.dbsnp}"""
params.indel_res = "-resource:mills,known=true,training=true,truth=true,prior=12.0 ${params.mills_indels}"

// Exome target BED directory and UG parameters
params.exom      = "${params.dbdir}/Agilent_SureSelect/hg19/bed"
params.dcov      = 1000
params.ug_call   = 50
params.ug_emit   = 10

/********************************************************************
 * Determine architecture-specific tools and substitute placeholders
 ********************************************************************/
def ARCH = System.getProperty('os.arch')
def JAVA = ( ARCH == "aarch64" ) ? params.java_aarch64 : params.java_amd64
def BWA  = ( ARCH == "aarch64" ) ? params.tools.aarch64.bwa.replace("{ngsutils}", params.ngsutils)
                                  : params.tools.amd64.bwa.replace("{ngsutils}", params.ngsutils)
def SAM  = ( ARCH == "aarch64" ) ? params.tools.aarch64.samtools.replace("{ngsutils}", params.ngsutils)
                                  : params.tools.amd64.samtools.replace("{ngsutils}", params.ngsutils)
def BED  = ( ARCH == "aarch64" ) ? params.tools.aarch64.bedtools.replace("{ngsutils}", params.ngsutils)
                                  : params.tools.amd64.bedtools.replace("{ngsutils}", params.ngsutils)
def PIC  = params.picard.replace("{java}", JAVA)
                        .replace("{mem}", params.mem)
                        .replace("{tmpdir}", params.tmpdir)
                        .replace("{ngsutils}", params.ngsutils)
def GATK = params.gatk.replace("{java}", JAVA)
                       .replace("{mem}", params.mem)
                       .replace("{tmpdir}", params.tmpdir)
                       .replace("{ngsutils}", params.ngsutils)

/********************************************************************
 * Define helper script paths (assuming they are in the working directory)
 ********************************************************************/
def COV     = file("coverage.sh").toString()
def VCF2SEX = file("vcf2sex.sh").toString()

/********************************************************************
 * FASTQ Channels and Sample ID
 * (Assumes FASTQ files are in the parent directory with names containing "R1" and "R2")
 ********************************************************************/
Channel.fromPath('../*R1*fastq.gz').map { r1 ->
    def base = r1.getName()
    def r2 = file(r1.toString().replace("R1", "R2"))
    tuple(base, r1, r2)
}.set { fastq_pairs }

// Extract sample ID from the first FASTQ file (using first two underscore fields)
fastq_pairs.first().map { it[0] }.subscribe { base ->
    def id = base.tokenize("_")[0..1].join("_").replace("_ex", "")
    println "Sample ID: $id"
    params.sample_id = id
}

/********************************************************************
 * Define channel for chromosomes (for per-chromosome steps)
 ********************************************************************/
def CHROMOSOMES = ["1","2","3","4","5","6","7","8","9","10","11","12", "1314", "1516", "1718", "1920", "2122XY"]
Channel.fromList(CHROMOSOMES).set { chrs }

/********************************************************************
 * Channels for intermediate files
 ********************************************************************/
// These channels will collect outputs from processes to feed the next steps.
Channel.create().set { grp_bams }
Channel.create().set { fixed_bams }
Channel.create().set { merged_bam }
Channel.create().set { filtered_bam }
Channel.create().set { intervals_file }
Channel.create().set { realigned_bam }
Channel.create().set { fixed_realigned_bam }
Channel.create().set { dedup_bam }
Channel.create().set { recal_grp_files }
Channel.create().set { recal_bams }
Channel.create().set { vcf_files }
Channel.create().set { merged_vcf }

/********************************************************************
 * Processes (in pipeline order)
 ********************************************************************/

// 1. ALIGN_AND_FIX: Align reads, add read groups, and fix mate information.
process ALIGN_AND_FIX {
    tag "$base"
    cpus 4
    publishDir "01_bam", mode: 'copy'
    input:
        tuple val(base), file(r1), file(r2) from fastq_pairs
    output:
        file("${base}.grp.bam") into grp_bams
        file("${base}.fixed.sorted.bam") into fixed_bams
    script:
    """
    SAMPLE=\$(echo ${base} | cut -d'_' -f1-2)
    BARCODE=\$(echo ${base} | cut -d'_' -f3)
    LANE=\$(echo ${base} | cut -d'_' -f4)
    
    echo "Aligning ${r1} and ${r2}" >> logs/${base}.align_and_fix.log
    ${BWA} mem -t ${task.cpus} -M ${params.refgz} ${r1} ${r2} 2>> logs/${base}.align_and_fix.log | \\
         ${PIC} AddOrReplaceReadGroups \\
              TMP_DIR=${params.tmpdir} \\
              I=/dev/stdin \\
              O=${base}.grp.bam \\
              SO=coordinate \\
              RGID=\$LANE \\
              RGLB=sureselect \\
              RGPL=illumina \\
              RGPU=\$BARCODE \\
              RGSM=\$SAMPLE 2>> logs/${base}.align_and_fix.log
    
    echo "Fixing mate information for ${base}.grp.bam" >> logs/${base}.align_and_fix.log
    ${PIC} FixMateInformation \\
          TMP_DIR=${params.tmpdir} \\
          INPUT=${base}.grp.bam \\
          OUTPUT=${base}.fixed.sorted.bam \\
          VALIDATION_STRINGENCY=SILENT \\
          CREATE_INDEX=true 2>> logs/${base}.align_and_fix.log
    """
}

// 2. MERGE_BAMS: Merge all fixed BAM files.
process MERGE_BAMS {
    publishDir "01_bam", mode: 'copy'
    input:
        file(fixed) from fixed_bams.collect()
    output:
        file("input.merged.bam") into merged_bam
    script:
    """
    tmp_in=""
    for bam in ${fixed.join(' ')}; do
       tmp_in="\$tmp_in I=\$bam"
    done
    echo "Merging BAM files: ${fixed.join(' ')}" >> logs/merge_bams.log
    ${PIC} MergeSamFiles \\
         TMP_DIR=${params.tmpdir} \\
         \$tmp_in \\
         OUTPUT=input.merged.bam \\
         SO=coordinate \\
         VALIDATION_STRINGENCY=SILENT \\
         CREATE_INDEX=false 2>> logs/merge_bams.log
    """
}

// 3. FILTER_BAM: Filter out malformed reads.
process FILTER_BAM {
    publishDir "01_bam", mode: 'copy'
    input:
        file(merged) from merged_bam
    output:
        file("input.merged.filtered.bam") into filtered_bam
        file("input.merged.filtered.bam.bai")
    script:
    """
    ${SAM} view -h ${merged} | \\
      awk 'substr(\$0,1,1)=="@" || (\$6 !~ /[0-9]+H/ && length(\$10)==length(\$11))' | \\
      ${SAM} view -Sb - > input.merged.filtered.bam
    ${SAM} index input.merged.filtered.bam 2>> logs/filter_bam.log
    """
}

// 4. REALIGN_TARGET_CREATOR: Identify intervals for local realignment.
process REALIGN_TARGET_CREATOR {
    cpus 4
    publishDir "01_bam", mode: 'copy'
    input:
        file(filtered) from filtered_bam
    output:
        file("input.merged.filtered.bam.intervals") into intervals_file
    script:
    """
    echo "Running RealignerTargetCreator on ${filtered}" >> logs/realign_target_creator.log
    ${GATK} -T RealignerTargetCreator \\
           -nt ${task.cpus} \\
           -R ${params.ref} \\
           -I ${filtered} \\
           -o input.merged.filtered.bam.intervals \\
           -known ${params.mills_indels} \\
           -known ${params.kg_indels} 2>> logs/realign_target_creator.log
    """
}

// 5. INDEL_REALIGNER: Perform local realignment around indels.
process INDEL_REALIGNER {
    cpus 4
    publishDir "01_bam", mode: 'copy'
    input:
        file(filtered) from filtered_bam
        file(intervals) from intervals_file
    output:
        file("input.merged.filtered.bam.realigned.bam") into realigned_bam
    script:
    """
    echo "Running IndelRealigner on ${filtered} with intervals ${intervals}" >> logs/indel_realigner.log
    ${GATK} -T IndelRealigner \\
           -R ${params.ref} \\
           -targetIntervals ${intervals} \\
           -I ${filtered} \\
           -o input.merged.filtered.bam.realigned.bam \\
           -model USE_SW \\
           -known ${params.mills_indels} \\
           -known ${params.kg_indels} \\
           -rf NotPrimaryAlignment 2>> logs/indel_realigner.log
    """
}

// 6. FIX_MATE_INFORMATION: Fix mate information on realigned BAM.
process FIX_MATE_INFORMATION {
    publishDir "01_bam", mode: 'copy'
    input:
        file(realigned) from realigned_bam
    output:
        file("input.merged.filtered.bam.realigned.bam.fixed.bam") into fixed_realigned_bam
        file("input.merged.filtered.bam.realigned.bam.fixed.bai")
    script:
    """
    echo "Fixing mate information for ${realigned}" >> logs/fix_mate_information.log
    ${PIC} FixMateInformation \\
          TMP_DIR=${params.tmpdir} \\
          INPUT=${realigned} \\
          OUTPUT=input.merged.filtered.bam.realigned.bam.fixed.bam \\
          SO=coordinate \\
          VALIDATION_STRINGENCY=LENIENT \\
          CREATE_INDEX=true 2>> logs/fix_mate_information.log
    """
}

// 7. MARK_DUPLICATES: Mark and remove PCR duplicates.
process MARK_DUPLICATES {
    publishDir "01_bam", mode: 'copy'
    input:
        file(fixed) from fixed_realigned_bam
    output:
        file("input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bam") into dedup_bam
        file("input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bai")
        file("input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.dupmetrics")
    script:
    """
    echo "Marking and deleting PCR duplicates for ${fixed}" >> logs/mark_duplicates.log
    ${PIC} MarkDuplicates \\
          TMP_DIR=${params.tmpdir} \\
          INPUT=${fixed} \\
          OUTPUT=input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bam \\
          METRICS_FILE=input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.dupmetrics \\
          REMOVE_DUPLICATES=true \\
          ASSUME_SORTED=true \\
          CREATE_INDEX=true \\
          VALIDATION_STRINGENCY=SILENT 2>> logs/mark_duplicates.log
    """
}

// 8. BASE_RECALIBRATOR: Generate recalibration groups per chromosome.
process BASE_RECALIBRATOR {
    tag "$chr"
    cpus 4
    publishDir "01_bam", mode: 'copy'
    input:
        file(dedup) from dedup_bam
        val chr from chrs
    output:
        file("input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bam.${chr}.recal.grp") into recal_grp_files
    script:
    """
    echo "Running BaseRecalibrator for chromosome ${chr}" >> logs/base_recalibrator.${chr}.log
    ${GATK} -T BaseRecalibrator \\
           -nct ${task.cpus} \\
           -R ${params.ref} \\
           -L ${params.exom}/hg19.chr${chr}.bed \\
           -I ${dedup} \\
           -o input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bam.${chr}.recal.grp \\
           -knownSites ${params.dbsnp} \\
           -knownSites ${params.mills_indels} \\
           -knownSites ${params.kg_indels} 2>> logs/base_recalibrator.${chr}.log
    """
}

// 9. PRINT_READS: Apply base quality score recalibration (per chromosome).
process PRINT_READS {
    tag "$chr"
    cpus 4
    publishDir "01_bam", mode: 'copy'
    input:
        file(dedup) from dedup_bam
        file(recal_grp) from recal_grp_files.filter { it.name.endsWith("recal.grp") }
        val chr from chrs
    output:
        file("input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bam.${chr}.recal.bam") into recal_bams
    script:
    """
    echo "Running PrintReads for chromosome ${chr}" >> logs/print_reads.${chr}.log
    ${GATK} -T PrintReads \\
           -R ${params.ref} \\
           -L ${params.exom}/hg19.chr${chr}.bed \\
           -I ${dedup} \\
           -o input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bam.${chr}.recal.bam \\
           -BQSR input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bam.${chr}.recal.grp 2>> logs/print_reads.${chr}.log
    """
}

// 10. VARIANT_CALLING: Call variants per chromosome.
process VARIANT_CALLING {
    tag "$chr"
    cpus 4
    publishDir "02_varcall", mode: 'copy'
    input:
        file(recal_bam) from recal_bams
        val chr from chrs
    output:
        file("chr${chr}.ug.raw.vcf") into vcf_files
    script:
    """
    echo "Running GATK UnifiedGenotyper for chromosome ${chr}" >> logs/variant_calling.${chr}.log
    ${GATK} -T UnifiedGenotyper \\
           -R ${params.ref} \\
           -I ${recal_bam} \\
           -L ${params.exom}/hg19.chr${chr}.flank100bp.bed \\
           --dbsnp ${params.dbsnp} \\
           -o chr${chr}.ug.raw.vcf \\
           -dcov ${params.dcov} \\
           -stand_call_conf ${params.ug_call} \\
           -stand_emit_conf ${params.ug_emit} \\
           -nt ${task.cpus} \\
           -glm BOTH 2>> logs/variant_calling.${chr}.log
    echo "Done with chromosome ${chr}" >> logs/variant_calling.${chr}.log
    """
}

// 11. MERGE_VCFS: Merge per-chromosome VCFs into a single VCF.
process MERGE_VCFS {
    publishDir "02_varcall", mode: 'copy'
    input:
        file(vcfs) from vcf_files.collect()
    output:
        file("${params.sample_id}.ug.raw.vcf") into merged_vcf
    script:
    """
    echo "Merging VCF files for sample ${params.sample_id}" >> logs/merge_vcfs.${params.sample_id}.log
    ( grep "#" chr1.ug.raw.vcf ; grep -hv '#' chr*.ug.raw.vcf | awk '\$1 !~ /_/' | sort -V ) > ${params.sample_id}.ug.raw.vcf 2>> logs/merge_vcfs.${params.sample_id}.log
    """
}

// 12. RECALIBRATE_SNP: Recalibrate SNP variants.
process RECALIBRATE_SNP {
    publishDir "02_varcall", mode: 'copy'
    input:
        file(vcf) from merged_vcf
    output:
        file("${params.sample_id}.ug.raw.snp.recal"),
        file("${params.sample_id}.ug.raw.snp.tranches")
    script:
    """
    echo "-GATK Recalibrator SNP" >> logs/recalibrate_snp.${params.sample_id}.log
    ${GATK} -T VariantRecalibrator \\
           -R ${params.ref} \\
           -input ${vcf} \\
           -recalFile ${params.sample_id}.ug.raw.snp.recal \\
           -tranchesFile ${params.sample_id}.ug.raw.snp.tranches \\
           --maxGaussians 6 \\
           ${params.snp_res} \\
           -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \\
           -mode SNP 2>> logs/recalibrate_snp.${params.sample_id}.log
    """
}

// 13. RECALIBRATE_INDEL: Recalibrate INDEL variants.
process RECALIBRATE_INDEL {
    publishDir "02_varcall", mode: 'copy'
    input:
        file(vcf) from merged_vcf
    output:
        file("${params.sample_id}.ug.raw.indel.recal"),
        file("${params.sample_id}.ug.raw.indel.tranches")
    script:
    """
    echo "-GATK Recalibrator INDEL" >> logs/recalibrate_indel.${params.sample_id}.log
    nINDEL=\$(grep -v "#" ${vcf} | awk 'length(\$5) != 1' | wc -l)
    if [ \$nINDEL -gt 8000 ]; then
        ${GATK} -T VariantRecalibrator \\
               -R ${params.ref} \\
               -input ${vcf} \\
               -recalFile ${params.sample_id}.ug.raw.indel.recal \\
               -tranchesFile ${params.sample_id}.ug.raw.indel.tranches \\
               --maxGaussians 4 \\
               ${params.indel_res} \\
               -an QD -an FS -an ReadPosRankSum \\
               -mode INDEL 2>> logs/recalibrate_indel.${params.sample_id}.log
    else
        echo ">>>> No INDELs to recalibrate" >> logs/recalibrate_indel.${params.sample_id}.log
        touch ${params.sample_id}.ug.raw.indel.recal ${params.sample_id}.ug.raw.indel.tranches
    fi
    """
}

// 14. APPLY_RECALIBRATION_SNP: Apply SNP recalibration.
process APPLY_RECALIBRATION_SNP {
    publishDir "02_varcall", mode: 'copy'
    input:
        file(vcf) from merged_vcf,
        file(snp_recal) from RECALIBRATE_SNP.out.collect(),
        file(snp_tranches) from RECALIBRATE_SNP.out.collect()
    output:
        file("${params.sample_id}.recalibratedSNPs.rawIndels.vcf") into recalibrated_vcf
    script:
    """
    echo "-GATK Apply Recalibrator SNP" >> logs/apply_recalibration_snp.${params.sample_id}.log
    ${GATK} -T ApplyRecalibration \\
           -R ${params.ref} \\
           -input ${vcf} \\
           -recalFile ${snp_recal} \\
           -tranchesFile ${snp_tranches} \\
           -o ${params.sample_id}.recalibratedSNPs.rawIndels.vcf \\
           --ts_filter_level 99.0 \\
           -mode SNP 2>> logs/apply_recalibration_snp.${params.sample_id}.log
    """
}

// 15. APPLY_RECALIBRATION_INDEL: Apply INDEL recalibration.
process APPLY_RECALIBRATION_INDEL {
    publishDir "02_varcall", mode: 'copy'
    input:
        file(recal_snp) from recalibrated_vcf,
        file(indel_recal) from RECALIBRATE_INDEL.out.collect(),
        file(indel_tranches) from RECALIBRATE_INDEL.out.collect()
    output:
        file("${params.sample_id}.ug.vqsr.vcf") into final_vcf
    script:
    """
    echo "-GATK Apply Recalibrator INDEL" >> logs/apply_recalibration_indel.${params.sample_id}.log
    if [ -s ${indel_recal} ]; then
        ${GATK} -T ApplyRecalibration \\
               -R ${params.ref} \\
               -input ${recal_snp} \\
               -recalFile ${indel_recal} \\
               -tranchesFile ${indel_tranches} \\
               -o ${params.sample_id}.ug.vqsr.vcf \\
               --ts_filter_level 95.0 \\
               -mode INDEL 2>> logs/apply_recalibration_indel.${params.sample_id}.log
    else
        echo ">>>> No INDELs to recalibrate" >> logs/apply_recalibration_indel.${params.sample_id}.log
        cp ${recal_snp} ${params.sample_id}.ug.vqsr.vcf
    fi
    """
}

// 16. VARIANT_FILTRATION: Apply hard filters to the VQSR VCF.
process VARIANT_FILTRATION {
    publishDir "02_varcall", mode: 'copy'
    input:
        file(vqsr_vcf) from final_vcf
    output:
        file("${params.sample_id}.ug.QC.vcf") into qc_vcf
    script:
    """
    echo "-GATK VariantFiltration" >> logs/variant_filtration.${params.sample_id}.log
    ${GATK} -T VariantFiltration \\
           -R ${params.ref} \\
           -o ${params.sample_id}.ug.QC.vcf \\
           --variant ${vqsr_vcf} \\
           --clusterWindowSize 10 \\
           --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \\
           --filterName "HARD_TO_VALIDATE" \\
           --filterExpression "DP < 5" \\
           --filterName "LowCoverage" \\
           --filterExpression "QUAL < 30.0" \\
           --filterName "VeryLowQual" \\
           --filterExpression "QUAL > 30.0 && QUAL < 50.0" \\
           --filterName "LowQual" \\
           --filterExpression "QD < 2.0" \\
           --filterName "LowQD" \\
           --filterExpression "MQ < 40.0" \\
           --filterName "LowMQ" \\
           --filterExpression "FS > 60.0" \\
           --filterName "StrandBias" 2>> logs/variant_filtration.${params.sample_id}.log
    """
}

// 17. COVERAGE_STATS: Compute coverage statistics.
process COVERAGE_STATS {
    publishDir "03_stats", mode: 'copy'
    input:
        // Using filtered bam as raw input here; adjust as needed.
        file(raw_bam) from filtered_bam,
        file(raw_bai) from file("01_bam/input.merged.filtered.bam.realigned.bam.fixed.bai"),
        file(dedup_bam) from dedup_bam,
        file(dedup_bai) from file("01_bam/input.merged.filtered.bam.realigned.bam.fixed.bam.dedup.bai")
    output:
        file("${params.sample_id}.1.raw.bam"),
        file("${params.sample_id}.1.dedup.bam"),
        file("${params.sample_id}.coverage.txt")
    script:
    """
    echo "Computing coverage stats for chromosome 1" >> logs/coverage_stats.${params.sample_id}.log
    cp ${raw_bai} ${raw_bam}.bai 2>> logs/coverage_stats.${params.sample_id}.log || true
    cp ${dedup_bai} ${dedup_bam}.bai 2>> logs/coverage_stats.${params.sample_id}.log || true
    ${SAM} view -b ${raw_bam} 1 > ${params.sample_id}.1.raw.bam 2>> logs/coverage_stats.${params.sample_id}.log
    ${SAM} view -b ${dedup_bam} 1 > ${params.sample_id}.1.dedup.bam 2>> logs/coverage_stats.${params.sample_id}.log
    ${SAM} index ${params.sample_id}.1.raw.bam 2>> logs/coverage_stats.${params.sample_id}.log
    ${SAM} index ${params.sample_id}.1.dedup.bam 2>> logs/coverage_stats.${params.sample_id}.log
    ${COV} ${params.sample_id} ${params.sample_id}.1.raw.bam ${params.sample_id}.1.dedup.bam > ${params.sample_id}.coverage.txt 2>> logs/coverage_stats.${params.sample_id}.log
    """
}

// 18. SEX_DETERMINATION: Run sex determination on the QC VCF.
process SEX_DETERMINATION {
    publishDir "03_stats", mode: 'copy'
    input:
        file(vcf) from qc_vcf
    output:
        file("${params.sample_id}.sex.txt")
    script:
    """
    echo "Estimating sex for id ${params.sample_id}" >> logs/sex_determination.${params.sample_id}.log
    ${VCF2SEX} ${vcf} > ${params.sample_id}.sex.txt 2>> logs/sex_determination.${params.sample_id}.log
    """
}

/********************************************************************
 * Workflow: define the execution order
 * (The channel connections ensure the same order as your Snakefile)
 ********************************************************************/
workflow {
    ALIGN_AND_FIX()
    MERGE_BAMS()
    FILTER_BAM()
    REALIGN_TARGET_CREATOR()
    INDEL_REALIGNER()
    FIX_MATE_INFORMATION()
    MARK_DUPLICATES()
    BASE_RECALIBRATOR()
    PRINT_READS()
    VARIANT_CALLING()
    MERGE_VCFS()
    RECALIBRATE_SNP()
    RECALIBRATE_INDEL()
    APPLY_RECALIBRATION_SNP()
    APPLY_RECALIBRATION_INDEL()
    VARIANT_FILTRATION()
    COVERAGE_STATS()
    SEX_DETERMINATION()
}
