#!/bin/bash

# Author: Tommaso Barberis
# Date: 03/11/2021
# Title: short pipeline with gatk for variant calling, example on "cranb" sample


echo -e "convert sam to bam"
samtools view -S -b G12-cranb_pe.sam > G12-cranb_pe.bam

echo -e "sorting bam file"
gatk SortSam -I G12-cranb_pe.bam -O G12-cranb_sorted.bam -SO queryname -VALIDATION_STRINGENCY SILENT 

echo -e "AddOReplaceReadGroups"
gatk AddOrReplaceReadGroups -I G12-cranb_sorted.bam -O G12-cranb_ReadGroups.bam -RGID name -RGLB name -RGPL ILLUMINA -RGPU name -RGSM sample -VALIDATION_STRINGENCY SILENT

echo -e "MarkDuplicateSpark"
gatk MarkDuplicatesSpark -I G12-cranb_ReadGroups.bam -O G12-cranb_MarkDuplicated.bam --spark-master local[16]

echo -e "HaplotyeCaller"
matk HaplotypeCaller --native-pair-hmm-threads 16 -R /home/ubuntu/index/Drosophila-suzukii-contig.fasta -I G12-cranb_MarkDuplicated.bam -O G12-cranb.g.vcf.gz -ERC GVCF
