# Align reads and post-process
SAMP=$1 # 'ERR1630013'
echo ${SAMP}
module load gatk
# fastqc ${SAMP}.fastq.gz -o .

# Trim reads
trimmomatic SE -phred33 \
../fastq/${SAMP}.fastq.gz ${SAMP}.trimmed.fastq.gz \
ILLUMINACLIP:TruSeq3-SE:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35 # (default 36)

# STAR alignment # CHANGE to count reads at the same time
STAR --runMode alignReads --genomeLoad NoSharedMemory \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat --genomeDir /home-2/gqi1@jhu.edu/workzfs-abattle4/lab_data/hg19/STAR_index_hg19_Gencode19 \
--outSAMattrRGline ID:${SAMP} SM:${SAMP} \
--outSAMmapqUnique 60 \
--outFileNamePrefix ${SAMP} \
--readFilesIn ./${SAMP}.trimmed.fastq.gz

# Mark duplicates
picard MarkDuplicates \
INPUT=${SAMP}Aligned.sortedByCoord.out.bam \
OUTPUT=${SAMP}Aligned.sortedByCoord.markdup.bam \
M=${SAMP}.dupmetrics.txt

# Idex BAM
samtools index ${SAMP}Aligned.sortedByCoord.markdup.bam

# SplitNCIgarReads
gatk --java-options "-Xmx25g" SplitNCigarReads \
-R ../../reference/hg19_combined.fa \
-I ${SAMP}Aligned.sortedByCoord.markdup.bam \
-O ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.bam

# Add read groups
gatk --java-options "-Xmx25g" AddOrReplaceReadGroups \
-I ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.bam \
-O ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.readgrps.bam \
--RGPL ILLUMINA \
--RGLB lib_${SAMP} \
--RGPU unit1 \
--RGSM ${SAMP} \
--SORT_ORDER coordinate

# Build table for base quality score recalibration
gatk --java-options "-Xmx25g" BaseRecalibrator \
-R ../../reference/hg19_combined.fa \
-I ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.readgrps.bam \
--known-sites ../../reference/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
-O ${SAMP}.BQSR.table

# Apply BQSR
gatk --java-options "-Xmx25g" ApplyBQSR \
-R ../../reference/hg19_combined.fa \
-I ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.readgrps.bam \
--bqsr-recal-file ${SAMP}.BQSR.table \
-O ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.readgrps.BQSR.bam

rm  ${SAMP}.trimmed.fastq.gz ${SAMP}Aligned.sortedByCoord.out.bam ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.bam ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.bai ${SAMP}Aligned.sortedByCoord.markdup.SplitNCigar.readgrps.bam
