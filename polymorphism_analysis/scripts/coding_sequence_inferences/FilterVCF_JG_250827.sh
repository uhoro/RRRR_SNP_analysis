#!/bin/bash

module load java/1.8.0
module load GenomeAnalysisTK/4.1.9.0
module load samtools/1.18
module load tabix/0.2.6

OUTDIR="FilteredVCF_JG_250827"
mkdir -p ${OUTDIR}

fname="all_samples_hap_flexa_240723_set240924"
gatk SelectVariants -O ${OUTDIR}/${fname}.snps.vcf.gz \
	-V GenotypeGVCFs/${fname}.vcf.gz -select-type SNP

gatk SelectVariants -O ${OUTDIR}/${fname}.indels.vcf.gz \
        -V GenotypeGVCFs/${fname}.vcf.gz -select-type INDEL

gatk VariantFiltration -O ${OUTDIR}/${fname}.snps.filtered.vcf.gz \
	-V ${OUTDIR}/${fname}.snps.vcf.gz \
	-filter "QD<2.0" --filter-name "QD2" \
	-filter "QUAL<30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

gatk VariantFiltration -O ${OUTDIR}/${fname}.indels.filtered.vcf.gz \
        -V ${OUTDIR}/${fname}.indels.vcf.gz \
        -filter "QD<2.0" --filter-name "QD2" \
        -filter "SOR > 10.0" --filter-name "SOR10" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"

gatk SelectVariants -O ${OUTDIR}/${fname}.snps.filtered.passed.vcf.gz \
	-V ${OUTDIR}/${fname}.snps.filtered.vcf.gz --exclude-filtered true

gatk SelectVariants -O ${OUTDIR}/${fname}.indels.filtered.passed.vcf.gz \
        -V ${OUTDIR}/${fname}.indels.filtered.vcf.gz --exclude-filtered true

snps_before_filt=$(gatk CountVariants -V ${OUTDIR}/${fname}.snps.filtered.vcf.gz)
snps_after_filt=$(gatk CountVariants -V ${OUTDIR}/${fname}.snps.filtered.passed.vcf.gz)

indels_before_filt=$(gatk CountVariants -V ${OUTDIR}/${fname}.indels.filtered.vcf.gz)
indels_after_filt=$(gatk CountVariants -V ${OUTDIR}/${fname}.indels.filtered.passed.vcf.gz)


bcftools sort -o ${OUTDIR}/${fname}.snps.filtered.passed.sorted.vcf.gz \
	${OUTDIR}/${fname}.snps.filtered.passed.vcf.gz
bcftools sort -o ${OUTDIR}/${fname}.indels.filtered.passed.sorted.vcf.gz \
	${OUTDIR}/${fname}.indels.filtered.passed.vcf.gz
bcftools index ${OUTDIR}/${fname}.snps.filtered.passed.sorted.vcf.gz
bcftools index ${OUTDIR}/${fname}.indels.filtered.passed.sorted.vcf.gz

bcftools concat -a -o ${OUTDIR}/${fname}.combined.filtered.passed.vcf.gz \
	${OUTDIR}/${fname}.snps.filtered.passed.sorted.vcf.gz \
	${OUTDIR}/${fname}.indels.filtered.passed.sorted.vcf.gz

bcftools sort -o ${OUTDIR}/${fname}.combined.filtered.passed.sorted.vcf.gz \
	${OUTDIR}/${fname}.combined.filtered.passed.vcf.gz

gatk IndexFeatureFile -I ${OUTDIR}/${fname}.combined.filtered.passed.sorted.vcf.gz

total_var_sites=$(gatk CountVariants -V ${OUTDIR}/${fname}.combined.filtered.passed.sorted.vcf.gz)

echo "fname: ${fname}"
echo "OUTDIR: ${OUTDIR}"
echo "SNPS before filtering: ${snps_before_filt}"
echo "SNPS after filtering: ${snps_after_filt}"

echo "INDELS before filtering: ${indels_before_filt}"
echo "INDELS after filtering: ${indels_before_filt}"

echo "Total variant sites: ${total_var_sites}"
