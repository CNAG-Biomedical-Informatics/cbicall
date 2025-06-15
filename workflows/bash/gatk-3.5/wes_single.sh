#!/usr/bin/env bash
#
#   WES pipeline Bash script.
#   This pipeline works at the sample level. 
#   For cohorts you will need to excute "wes_cohort.sh". This way, 
#   if a new relative is added, you can easily add it a posteriori.
#
#   Last Modified; March/05/2025
#
#   $VERSION taken from CBICall
#
#   Copyright (C) 2025 Manuel Rueda - CNAG (manuel.rueda@cnag.eu)

set -eu

function usage {

    USAGE="""
    Usage: $0 -t n_threads

    NB1: The script is expecting that you follow STSI nomenclature for samples. 
         This will work    => MA0004701P_ex_S5_L001_R1_001.fastq.gz
         This wll NOT work => 62236832_S1_LALL_R1_001.fastq.gz (missing _ex before S1)

MA00047_exome
└── MA0004701P_ex  <--- ID taken from here
    ├── MA0004701P_ex_S5_L001_R1_001.fastq.gz
    ├── MA0004701P_ex_S5_L001_R2_001.fastq.gz
    ├── MA0004701P_ex_S5_L002_R1_001.fastq.gz
    ├── MA0004701P_ex_S5_L002_R2_001.fastq.gz
    └── cbicall_bash_wes_single_146657420113136 <- The script expects that you are submitting the job from inside this directory
        ├── 01_bam
        ├── 02_varcall
        └── 03_stats
    """
    echo "$USAGE"
    exit 1
}

# Check arguments
if [ $# -ne 2 ]
then
  usage
fi

# parsing Arguments
key="$1"
case $key in
    -t|--t)
    THREADS="$2"
esac

# Determine the directory where the script resides
BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source parameters.sh from the same directory
source "$BINDIR/parameters.sh"

# Scripts to calculate miscellanea stats on coverage (chr1) and determine sex
COV=$BINDIR/coverage.sh
VCF2SEX=$BINDIR/vcf2sex.sh

# Defining directories
DIR=$( pwd )

# Check that nomenclature exists
if [[ $DIR != *cbicall_bash_wes_single* ]]
 then 
  usage
fi 

BAMDIR=$DIR/01_bam
VARCALLDIR=$DIR/02_varcall
STATSDIR=$DIR/03_stats
LOGDIR=$DIR/logs
mkdir $BAMDIR
mkdir $VARCALLDIR
mkdir $STATSDIR
mkdir $LOGDIR

# id will be used to name final VCF files. Note that we can't use fq info since naming can be wrong (e.g., MA00M5_ instead of MA0000502M_)
id=$( echo $DIR | awk -F'/' '{print $(NF-1)}' | awk -F'_' '{print $1}' )

# Create a log file with STDERR
LOG=$LOGDIR/cbicall_bash_wes_single.log

# Function to update a BAM file name by removing any previous .bam extension
# and appending a new suffix
function update_bam_name {
    local file="${1%.bam}"   # remove trailing .bam if present
    local suffix="$2"
    echo "${file}.${suffix}.bam"
}

# GATK Best practices
# https://www.broadinstitute.org/gatk/guide/article?id=3060

# NB: Note that, although some of these tasks can be run in parallel with gnu parallel
# in general JAVA takes tons of memory x process
##################################
#                                #
#       DATA CLEANUP             # 
#                                #
##################################

# BWA alignment to the reference genome
# The default mismtach is 4% of the read length
# Adding RG (read groups) with Picard
# We'll need the RG information in case we need to merge BAM files from different samples.

# Working at FASTQ level
echo "Aligning reads to REF with BWA and adding groups with Picard"

for R1 in ../*_R1_*fastq.gz
do
 R2=${R1/_R1_/_R2_} # substitute R1 by R2 in $R1 to create $R2
 string_R1=$( echo $R1 | awk -F'/' '{print $NF}' )
 string_R2=${string_R1/_R1_/_R2_}
 array=(${string_R1//_/ })  # split $R1 fq by '_" and load an array to get misc sample info
 sample="${array[0]}_${array[1]}" # ID_ex
 barcode_seq=${array[2]}
 lane=${array[3]}
 read_number=${array[4]}
 part=${array[5]}
 part=${part%.*.*}  # getting rid of .fastq.gz
 #echo "$sample $barcode_seq $index $lane $read_number $part"

 echo " $string_R1 / $string_R2"

 in=/dev/stdin
 out=$BAMDIR/$string_R1.grp.bam
 $BWA mem -t$THREADS -M $REFGZ $R1 $R2 2> $LOG | \
      $PIC \
      AddOrReplaceReadGroups \
      TMP_DIR=$TMPDIR \
      I=$in  \
      O=$out \
      SO=coordinate \
      RGID=$lane \
      RGLB=sureselect \
      RGPL=illumina \
      RGPU=$barcode_seq \
      RGSM=$sample  2>> $LOG

 # Fixmate information and bam indexing with Picard
 in=$out
 out=$BAMDIR/$string_R1.fixed.sorted.bam
 echo " Fixing Mates"
 $PIC \
      FixMateInformation \
      TMP_DIR=$TMPDIR \
      INPUT=$in \
      OUTPUT=$out \
      VALIDATION_STRINGENCY=SILENT \
      CREATE_INDEX=true 2>> $LOG
done

# Now we merge all BAM files x sample with Picard (samtools did not merge the headers)
cd $BAMDIR
echo "Merging all BAMs"
# We append all inputs over the string tmp_in
tmp_in=''
for tmp_bam in $( ls *.fixed.sorted.bam )
do
 tmp_in="${tmp_in} I=$tmp_bam "
done
in=$tmp_in
out=input.merged.bam
$PIC \
      MergeSamFiles \
      TMP_DIR=$TMPDIR \
      $in \
      OUTPUT=$out \
      SO=coordinate \
      VALIDATION_STRINGENCY=SILENT \
      CREATE_INDEX=false 2>> $LOG
#      CREATE_INDEX=true 2>> $LOG # mrueda 022025

# Deleting unused *bam *bai to save space
rm *fastq.gz.*ba?
 
# Filter out malformed reads
# Added 022025 due to errors on mrueda-ws5
in=$out
out=$(update_bam_name "$in" filtered)
$SAM view -h $in | \
awk 'substr($0,1,1)=="@" || ($6!~/[0-9]+H/ && length($10)==length($11))' | \
$SAM view -Sb - > $out
$SAM index $out

# Local realignment around indels -> STEP1
# This step is CPU consuming
echo "Local realignment around indels"
in=$out
out=$in.intervals

echo " STEP 1"
$GATK \
      -T RealignerTargetCreator \
      -nt $THREADS \
      -R $REF \
      -I $in \
      -o $out \
      -known $MILLS_INDELS \
      -known $KG_INDELS 2>> $LOG

# Local realignment around indels -> STEP2
in=$in
intervals=$out
out=$(update_bam_name "$in" realigned)

echo " STEP 2"
$GATK \
      -T IndelRealigner \
      -R $REF \
      -targetIntervals $intervals \
      -I $in \
      -o $out \
      -model USE_SW \
      -known $MILLS_INDELS \
      -known $KG_INDELS \
      -rf NotPrimaryAlignment 2>> $LOG
 
# Local realignment around indels -> STEP3
# Verify mate-pair information and sort by coordinate
in=$out
out=$(update_bam_name "$in" fixed)
echo " STEP 3"
$PIC \
      FixMateInformation \
      TMP_DIR=$TMPDIR \
      INPUT=$in \
      OUTPUT=$out \
      SO=coordinate \
      VALIDATION_STRINGENCY=LENIENT \
      CREATE_INDEX=true 2>> $LOG

# Marking and deleting PCR duplicates per sample
in=$out
out=$(update_bam_name "$in" dedup)

echo "Marking and deleting PCR duplicates"
$PIC \
      MarkDuplicates \
      TMP_DIR=$TMPDIR \
      INPUT=$in \
      OUTPUT=$out \
      METRICS_FILE=$out.dupmetrics \
      REMOVE_DUPLICATES=true \
      ASSUME_SORTED=true \
      CREATE_INDEX=true \
      VALIDATION_STRINGENCY=SILENT 2>> $LOG

in4chr=$out
# BQSR will be performed by chr
# GATK -L region -> allows for parallelization
for chr in {1..12} 1314 1516 1718 1920 2122XY  # This divisions are because of the format of the Agilent Sureselect BED files
do
 echo "---------------"
 echo "Chromosome...$chr" 

 REG=$EXOM/hg19.chr$chr.bed
 chrN=chr$chr

 # Base Quality Score Recalibration (BQSR)
 echo "Base quality recalibration"
 in=$in4chr
 out="${in%.bam}.$chrN.recal.grp"

 echo " STEP 1"
 # BaseRecalibrator currently does not support parallel execution with nt
 $GATK \
        -T BaseRecalibrator \
        -nct $THREADS \
        -R $REF \
        -L $REG \
        -I $in \
        -o $out \
        -knownSites $dbSNP \
        -knownSites $MILLS_INDELS \
        -knownSites $KG_INDELS  2>> $LOG

 in=$in
 bqsr=$out
 out=$(update_bam_name "$in" "$chrN.recal")

 echo " STEP 2"
 $GATK \
        -T PrintReads \
        -R $REF \
        -L $REG \
        -I $in \
        -o $out \
        -BQSR $bqsr 2>> $LOG

 ##################################
 #                                #
 #      VARIANT DISCOVERY         # 
 #                                #
 ##################################
 # GATK Variant Calling
 # UnifiedGenotyper x chr
 # --dbSNP argument adds rs IDs
 SSexome=$EXOM/hg19.chr$chr.flank100bp.bed
 in=$out                                                                                    
 out=$VARCALLDIR/chr$chr.ug.raw.vcf

 echo "GATK UG"
 $GATK \
      -T UnifiedGenotyper \
      -R $REF \
      -I $in \
      -L $SSexome \
      --dbsnp $dbSNP \
      -o $out \
      -dcov $DCOV \
      -stand_call_conf $UG_CALL \
      -stand_emit_conf $UG_EMIT \
      -nt $THREADS \
      -glm BOTH 2>> $LOG

 echo "Done with chr$chr!!"
done
echo "---------------"

# Now we merge all VCF's 
cd $VARCALLDIR
( grep "#" chr1.ug.raw.vcf ; grep -hv '#'  chr*.ug.raw.vcf | awk '$1 !~ /_/' | sort -V ) > $id.ug.raw.vcf # Keeping only chr{?,??}

# Variant Recalibrator
# http://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr
#       -nt $THREADS \  <== SLOWER
# VQSR does not perform well (if at all) on a single sample. 
# It can work with whole genome sequence, but if you're working with exome, there's just too few variants. 
# Our recommendation for dealing with this is to get additional sample bams from the 1000Genomes project and add them to your callset (see this presentation for details).
# Deleted "-std 10.0 -percentBad 0.12" as it not longer compatible with GATK 3.5

# SNP VCF VariantRecalibrator
echo "-GATK Recalibrator SNP"
# Count SNPs (assumes SNPs have a single-base alternate allele)
nSNP=$(grep -v "^#" $id.ug.raw.vcf | awk 'length($5)==1' | wc -l)
echo "Found $nSNP SNP variants" >> $LOG
if [[ $nSNP -gt 1000 ]]; then  # so that it works with test data
  $GATK \
      -T VariantRecalibrator \
      -R $REF \
      -input $id.ug.raw.vcf \
      -recalFile $id.ug.raw.snp.recal \
      -tranchesFile $id.ug.raw.snp.tranches \
      --maxGaussians 6 \
      $SNP_RES \
      -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
      -mode SNP 2>> $LOG
else
  echo "WARNING: Only $nSNP SNP variants found. Skipping VariantRecalibrator for SNPs." >> $LOG
  cp $id.ug.raw.vcf $id.ug.raw.snp.recal
  echo "SKIP" > $id.ug.raw.snp.tranches
fi

# INDEL VCF VariantRecalibrator
echo "-GATK Recalibrator INDEL"
nINDEL=$(grep -v "^#" $id.ug.raw.vcf | awk 'length($5) != 1' | wc -l)
if [[ $nINDEL -gt 8000 ]]; then
  $GATK \
      -T VariantRecalibrator \
      -R $REF \
      -input $id.ug.raw.vcf \
      -recalFile $id.ug.raw.indel.recal \
      -tranchesFile $id.ug.raw.indel.tranches \
      --maxGaussians 4 \
      $INDEL_RES \
      -an QD -an FS -an ReadPosRankSum \
      -mode INDEL 2>> $LOG
else
  echo ">>>> No INDELs to recalibrate (only $nINDEL found)" >> $LOG
  cp $id.ug.raw.vcf $id.ug.raw.indel.recal
  echo "SKIP" > $id.ug.raw.indel.tranches
fi

# SNP VCF ApplyRecalibration
echo "-GATK Apply Recalibrator SNP"
if grep -q "SKIP" $id.ug.raw.snp.tranches; then   # so that it works with test data
  echo "Skipping ApplyRecalibration for SNPs due to low variant count" >> $LOG
  cp $id.ug.raw.vcf recalibratedSNPs.rawIndels.vcf
else
  $GATK \
      -T ApplyRecalibration \
      -R $REF \
      -input $id.ug.raw.vcf \
      -recalFile $id.ug.raw.snp.recal \
      -tranchesFile $id.ug.raw.snp.tranches \
      -o recalibratedSNPs.rawIndels.vcf \
      --ts_filter_level 99.0 \
      -mode SNP 2>> $LOG
fi

# INDEL VCF ApplyRecalibration
echo "-GATK Apply Recalibrator INDEL"
if grep -q "SKIP" $id.ug.raw.indel.tranches; then
  echo "Skipping ApplyRecalibration for INDELs due to low variant count" >> $LOG
  cp recalibratedSNPs.rawIndels.vcf $id.ug.vqsr.vcf
else
  $GATK \
      -T ApplyRecalibration \
      -R $REF \
      -input recalibratedSNPs.rawIndels.vcf \
      -recalFile $id.ug.raw.indel.recal \
      -tranchesFile $id.ug.raw.indel.tranches \
      -o $id.ug.vqsr.vcf \
      --ts_filter_level 95.0 \
      -mode INDEL 2>> $LOG
fi

# VCF QC filtration
echo "-GATK Filter"
$GATK \
      -T VariantFiltration \
      -R $REF \
      -o $id.ug.QC.vcf \
      --variant $id.ug.vqsr.vcf \
      --clusterWindowSize 10 \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --filterExpression "DP < 5 " \
      --filterName "LowCoverage" \
      --filterExpression "QUAL < 30.0 " \
      --filterName "VeryLowQual" \
      --filterExpression "QUAL > 30.0 && QUAL < 50.0 " \
      --filterName "LowQual" \
      --filterExpression "QD < 2.0 " \
      --filterName "LowQD" \
      --filterExpression "MQ < 40.0 " \
      --filterName "LowMQ" \
      --filterExpression "FS > 60.0 " \
      --filterName "StrandBias"  2>> $LOG 

##################################
#                                #
#      STATS on COVERAGE         # 
#                                #
##################################
# $0 sid orig_bam dedup_bam
# Stats for chr1 only
# Splitting the BAM file x Chr while keeping the header ^@SQ and read groups @RG
# Note that GATK needs headers that match the contigs af $REF 
if [[ $REF == *b37*.fasta ]]
 then
  chrN=1
 else
  chrN=chr1
fi
echo "Computing stats on coverage for $chrN"
# NB: Tested on 04/26/16
# bam_raw=input.merged.bam
#sampleID	mean_coverage	ten_reads%	nonduplicate%	mean_insert_size	reads_in_exome%	reads_out_of_exome%
#sample_XXX		96.1		98.2		97.5		208.726		99.7		 0.3
# bam_raw=input.merged.bam.realigned.bam.fixed.bam
#sampleID	mean_coverage	ten_reads%	nonduplicate%	mean_insert_size	reads_in_exome%	reads_out_of_exome%
#sample_XXX		96.1		98.2		97.5		208.73		99.7		 0.3
# Using input.merged.bam.realigned.bam.fixed.bam

bam_raw=$BAMDIR/input.merged.filtered.realigned.fixed.bam
bam_raw_index=$BAMDIR/input.merged.filtered.realigned.fixed.bai
bam_dedup=$BAMDIR/input.merged.filtered.realigned.fixed.dedup.bam
bam_dedup_index=$BAMDIR/input.merged.filtered.realigned.fixed.dedup.bai
out_raw=$STATSDIR/$chrN.raw.bam
out_dedup=$STATSDIR/$chrN.dedup.bam
stats_log=$STATSDIR/coverage.txt

# Before using $SAM view we need to rename the index (must be foo.bam.bai instead of foo.bai)
# This step is also important for mit_cohort.sh
cp $bam_raw_index $bam_raw.bai
cp $bam_dedup_index $bam_dedup.bai
$SAM view -b $bam_raw   $chrN > $out_raw
$SAM view -b $bam_dedup $chrN > $out_dedup
$SAM index $out_raw
$SAM index $out_dedup
$COV $id $out_raw $out_dedup > $stats_log

##################################
#                                #
#      SEX DETERMINATION         # 
#                                #
##################################
echo "Estimating Sex of the sample"
sex_log=$STATSDIR/$id.sex.txt
$VCF2SEX $VARCALLDIR/$id.ug.QC.vcf  > $sex_log

# Fin
echo "All done!!!"
exit 
