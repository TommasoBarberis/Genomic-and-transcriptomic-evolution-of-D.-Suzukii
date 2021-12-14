#!/bin/bash

## author: Tommaso Barberis
## date: 03/122/2021
## description: short script that given a list of genes and a vcf file as input and an another vcf file name as output, returns in this last file SNPs that are founded in the gene list.


list_gene="/data/home/tbarberis/projet_S3/localdata/RNAseq/DE_gene/fruit_G0.txt"

# recopy header 
grep "#" $1 >> $2

cat $list_gene | while read CONTIG; do
	grep -v "#" $1 | grep $CONTIG >> $2
done
