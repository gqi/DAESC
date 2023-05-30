#!/bin/bash -l
#SBATCH --time=72:0:0
#SBATCH --mem=23G
#SBATCH --partition=shared
#SBATCH --account=abattle4
module load gatk
module load bcftools
module load samtools
module load htslib

gatk --java-options "-Xmx18g" CombineGVCFs \
-R ../../reference/hg19_combined.fa \
--variant AZ_genotype.g.vcf.gz \
--variant HP1502401_genotype.g.vcf.gz \
--variant HP1504101T2D_genotype.g.vcf.gz \
--variant HP1504901_genotype.g.vcf.gz \
--variant HP1506401_genotype.g.vcf.gz \
--variant HP1507101_genotype.g.vcf.gz \
--variant HP1508501T2D_genotype.g.vcf.gz \
--variant HP1509101_genotype.g.vcf.gz \
--variant HP1525301T2D_genotype.g.vcf.gz \
--variant HP1526901T2D_genotype.g.vcf.gz \
-O combined.g.vcf.gz

gatk --java-options "-Xmx18g" GenotypeGVCFs \
-R ../../reference/hg19_combined.fa \
-V combined.g.vcf.gz \
-new-qual true \
-O genotype.vcf.gz

# bcftools view -m 2 -M 2 -v snps genotype.vcf.gz > genotype_snps.vcf
# bgzip genotype_snps.vcf
# tabix -p vcf genotype_snps.vcf.gz
# 'AVG(FMT/DP)>10 & AVG(FMT/GQ)>15' - min mimum DP and minimum GQ across all samples

# Average DP and GQ; biallelic SNP
bcftools view -m 2 -M 2 -v snps -i 'AVG(FMT/DP)>10 & AVG(FMT/GQ)>15' \
genotype.vcf.gz > genotype_snps_filtered.vcf
bgzip genotype_snps_filtered.vcf
tabix -p vcf genotype_snps_filtered.vcf.gz

