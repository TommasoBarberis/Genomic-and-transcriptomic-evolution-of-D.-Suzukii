#!/bin/bash -i

# Author: Tommaso Barberis
# Date: 13/11/2021
# Short script to generate mpileup file from sorted bam

samtools mpileup -d 5000 -q 20 -f ../ref/Drosophila-suzukii-contig.fasta sorted/G12-cranb_sorted.bam > mpileup/G12-cranb.mpileup

mailSender "G12 cranb ended" tommasobarberis98@gmail.com "send from pedago-ngs"

