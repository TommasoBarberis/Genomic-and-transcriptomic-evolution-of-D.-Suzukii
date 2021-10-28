1. genome index
```
bwa index -p Drosophila-suzukii-contig.fasta.fai -a bwtsw Drosophila-suzukii-contig.fasta.fai
```

- __-p__: index name
- __-a__: index algorithm (_bwtsw_ for long genome)


2. mapping
```
bwa aln -t 14 Drosophila-suzukii-contig.fasta.fai sample_1.fastq.gz > sample_1.sai
bwa aln -t 14 Drosophila-suzukii-contig.fasta.fai sample_2.fastq.gz > sample_2.sai

bwa sampe Drosophila-suzukii-contig.fasta.fai sample_1.sai sample_2.sai sample_1.fastq.gz sample_2.fastq.gz > sample_pe.sam
```
