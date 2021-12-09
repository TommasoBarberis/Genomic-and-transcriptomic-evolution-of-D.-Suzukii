#!/bin/bash

list_gene="/data/home/tbarberis/projet_S3/localdata/RNAseq/DE_gene/fruit_G0.txt"

grep "#" $1 >> $2

cat $list_gene | while read CONTIG; do
	grep -v "#" $1 | grep $CONTIG >> $2
#	echo $new_lines #>> $2
#	if [[ -n $new_lines ]]; then
#		echo $new_lines >> $2
#	fi
done

# cat *.out.recode.vcf >> $2
#rm -r *.out*
