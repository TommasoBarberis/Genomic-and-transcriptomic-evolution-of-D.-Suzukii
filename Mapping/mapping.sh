#!/bin/bash

# short script to map DNAseq reads from D. suzukii poolseq
# Author: Tommaso BARBERIS
# Date: 17/10/2021



DATA=/home/ubuntu/fastq/
OUTPUT=/home/ubuntu/mapping/
INDEX=/home/ubuntu/index/

echo -e '## Mapping'



# mapping using bwa

bwa aln -t 14 ~/index/Drosophila-suzukii-contig.fasta.fai ~/fastq/G12-cranb_1.fastq.gz > G12-cranb_1.sai

bwa aln -t 14 ~/index/Drosophila-suzukii-contig.fasta.fai ~/fastq/G12-cranb_2.fastq.gz > G12-cranb_2.sai

bwa sampe  ~/index/Drosophila-suzukii-contig.fasta.fai ~/G12-cranb_1.sai ~/G12-cranb_2.sai ~/fastq/G12-cranb_1.fastq.gz ~/fastq/G12-cranb_2.fastq.gz > G12-cranb_pe.sam
