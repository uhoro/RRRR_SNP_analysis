#!/bin/bash


declare -a filter_dirs=("JG_SW_250829")
declare -a sample_arr=("C1G11" "C1C5v2" "C3D7v2")

for dir in "${filter_dirs[@]}";do
	for sample in "${sample_arr[@]}";do
		echo "start processing ${sample} in ${dir}"
		python vcf2fasta/vcf2fasta.py \
		--fasta data/Cflexa_genome_2024-07-16/Cflexa.scaffolds.Hoceani_pruned.fa \
		--vcf data/FilteredVCFs/${dir}/${sample}_filt.sorted.vcf.gz \
		--gff data/Cflexa_genome_2024-07-16/Cflexa.Hoceani_prunedv3.gff3 \
		--feat CDS --blend --addref \
		-o out/${sample}_${dir}


	done
done
