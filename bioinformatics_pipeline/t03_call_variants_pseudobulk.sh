#!/bin/bash -l
#SBATCH --time=72:0:0
#SBATCH --mem=35G
#SBATCH --partition=shared
#SBATCH --account=abattle4

module load samtools
module load gatk
sampidvec=("AZ" "HP1502401" "HP1504101T2D" "HP1504901" "HP1506401" "HP1507101" "HP1508501T2D" "HP1509101" "HP1525301T2D" "HP1526901T2D")
SAMP=${sampidvec[$1]}
echo $SAMP

gatk --java-options "-Xmx25g" AddOrReplaceReadGroups \
-I ${SAMP}_merged.bam \
-O ${SAMP}_merged.readgrp.bam \
--RGPL ILLUMINA \
--RGLB ${SAMP} \
--RGPU unit1 \
--RGSM ${SAMP} \
--SORT_ORDER coordinate

module load samtools
samtools index ${SAMP}_merged.readgrp.bam

gatk --java-options "-Xmx33g" HaplotypeCaller  \
-R /work-zfs/abattle4/guanghao/single_cell/data/t2d_pancreas/reference/hg19_combined.fa \
-I ${SAMP}_merged.readgrp.bam \
-O ${SAMP}_genotype.g.vcf.gz \
-ERC GVCF
