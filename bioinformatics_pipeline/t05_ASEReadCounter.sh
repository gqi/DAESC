module load gatk
SAMP=$1
echo ${SAMP}

gatk ASEReadCounter \
-R ../../reference/hg19_combined.fa \
-I ../bam_trimmed_filtered/${SAMP}Aligned.sortedByCoord.markdup.bam \
-V ../genotype/genotype_snps_filtered.vcf.gz \
--verbosity ERROR \
-O ${SAMP}_asect