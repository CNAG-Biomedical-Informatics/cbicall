#   $VERSION taken from CBICall

# Paths
DATADIR=/media/mrueda/2TBS
DBDIR=$DATADIR/Databases
NGSUTILS=$DATADIR/NGSutils

# Environment
export TMPDIR=$DATADIR/tmp
export LC_ALL=C
export GATK_DISABLE_AUTO_S3_UPLOAD=true   # disable unintended S3 uploads

# Memory & architecture
MEM=8G
ARCH=$(uname -m)

# Java & tool binaries per architecture
if [ "$ARCH" == "aarch64" ]; then
    JAVA=/usr/lib/jvm/java-8-openjdk-arm64/bin/java
    BWA=$NGSUTILS/bwa-0.7.18_arm64/bwa
    SAM=$NGSUTILS/samtools-0.1.19_arm64/samtools
    BED=$NGSUTILS/bedtools2_arm64/bin/bedtools
else
    JAVA=/usr/lib/jvm/java-8-openjdk-amd64/bin/java
    BWA=$NGSUTILS/bwa-0.7.18/bwa
    SAM=$NGSUTILS/samtools-0.1.19/samtools
    BED=$NGSUTILS/bedtools2/bin/bedtools
fi

# Picard (shared by GATK3 & bed conversion)
PIC="$JAVA -Xmx$MEM -Djava.io.tmpdir=$TMPDIR -jar $NGSUTILS/picard-2.25/build/libs/picard.jar"

# GATK 3.5 (legacy)
GATK="$JAVA -Xmx$MEM -Djava.io.tmpdir=$TMPDIR -jar $NGSUTILS/gatk/gatk-3.5/GenomeAnalysisTK.jar"

# GATK 4+ launcher (recommended)
GATK4_CMD="$NGSUTILS/gatk/gatk-4.6.2.0/gatk --java-options -Xmx${MEM}"

# MToolBox directory
MTOOLBOXDIR=$NGSUTILS/MToolBox-master/MToolBox

# GATK bundle & reference (b37)
BUNDLE=$DBDIR/GATK_bundle/b37
REF=$BUNDLE/references_b37_Homo_sapiens_assembly19.fasta
REFGZ=$BUNDLE/references_b37_Homo_sapiens_assembly19.fasta.gz
REF_DICT=$BUNDLE/references_b37_Homo_sapiens_assembly19.dict

# Variant resources
dbSNP=$DBDIR/dbSNP/human_9606_b144_GRCh37p13/All_20160408.vcf.gz
MILLS_INDELS=$BUNDLE/b37_Mills_and_1000G_gold_standard.indels.b37.vcf.gz
KG_INDELS=$BUNDLE/b37_1000G_phase1.indels.b37.vcf.gz
HAPMAP=$BUNDLE/b37_hapmap_3.3.b37.vcf.gz
OMNI=$BUNDLE/b37_1000G_omni2.5.b37.vcf.gz

# Training sets for VQSR
SNP_RES="-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
         -resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI \
         -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbSNP"
INDEL_RES="-resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS_INDELS"

# Exome targets
EXOME_BED=$BUNDLE/b37_Broad.human.exome.b37.bed
INTERVAL_LIST=$BUNDLE/b37_Broad.human.exome.b37.interval_list

# Agilent SureSelect Whole Exome
EXOM=$DBDIR/Agilent_SureSelect/hg19/bed

# UnifiedGenotyper parameters (legacy)
DCOV=1000
UG_CALL=50
UG_EMIT=10

