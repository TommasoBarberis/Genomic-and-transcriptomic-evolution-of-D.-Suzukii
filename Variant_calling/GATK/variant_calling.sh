#!/bin/bash

# Author: Tommaso Barberis
# Date: 03/11/2021
# Title: short pipeline with gatk for variant calling, example on "cranb" sample
# GATK version: gatk-4.2.2.0


echo -e "convert sam to bam"
samtools view -S -b G12-cranb_pe.sam > G12-cranb_pe.bam

echo -e "sorting bam file"
gatk SortSam -I G12-cranb_pe.bam -O G12-cranb_sorted.bam -SO queryname -VALIDATION_STRINGENCY SILENT 

echo -e "AddOReplaceReadGroups"
gatk AddOrReplaceReadGroups -I G12-cranb_sorted.bam -O G12-cranb_ReadGroups.bam -RGID name -RGLB name -RGPL ILLUMINA -RGPU name -RGSM sample -VALIDATION_STRINGENCY SILENT

echo -e "MarkDuplicateSpark"
gatk MarkDuplicatesSpark -I G12-cranb_ReadGroups.bam -O G12-cranb_MarkDuplicated.bam --spark-master local[16]

echo -e "HaplotyeCaller"
gatk HaplotypeCaller -R ~/ref/Drosophila-suzukii-contig.fasta -I ~/G12_cerise_MarkDuplicated.bam --sample-ploidy 80 -O G12_cerise_MarkDuplicated.bam.vcf --max-genotype-count 91881
