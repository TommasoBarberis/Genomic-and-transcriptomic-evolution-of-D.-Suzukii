# Genomic-and-transcriptomic-evolution-of-D.-Suzukii


## Table of content

1. Context
2. Project
3. Data 
4. Requirements
5. Polymorphism analysis
6. TE analysis
7. Authors
8. Project status


## Context

Since its recent invasion in the European and American continents, Drosophila suzukii, has become a major pest of berry crops. 
Contrary to other Drosophila species that develop themselves on damaged or rotten fruit, D. suzukii is able to lay its eggs in healthy fruit before
harvest, using its sclerotinized ovipositor. In order to control the D. suzukii population it is therefore required to improve our knowledge
about its ability to adapt to a local environment outside its origin areas.

A previous study focused on phenotipic changes during a selective experiment showed a temporal adaptation of D. suzukii populations in the three different 
fruit environments. This process produces a fitness improvement which corresponds to a greater ability to reproduce and could explain the adaptation abilities 
of D. suzukii. Based on the same experiment, a second study has been started by the LBBE to see if the D. suzukii genome and transcriptome are affected by the environment.
 
The project is focused on both genes and Transposable Elements (TE) expression in order to find the variants explaining the local environment adaptation of D. suzukii. 
A significant number of differentially expressed genes across generations has been found with RNAseq analysis between the G0 and the G7 generations. These variations seems to be independant of the fruit media. Moreover, in differentially expressed (DE) genes an increase of the number of Single Nucleotide Polymorphism (SNP) was found in G7 compared to G0, but this result was not expected. Natural selection generally keeps alleles which allow the best adaptation, and it is synonymous with a decrease in SNP number. 


## Project

### Polymorphism analysis
The goal here to analyse D. suzukii polymorphism performing a variant calling on the whole genome in order to count the number of different genomic positions between G0 and G12, and between fruits media. The position of the genetic variations will be compared with previous SNP analysis conducted on RNA-seq data to see if we also find an increase of the SNP number. To have a view more global on the genome's diversity, several diversity measures will be calculated : pi, Watterson's theta, Tajima's D.

### TE analysis
TEs represent 47% of the D. suzukii genome and can be important in the adaptation process, interfering in gene regulation during evolution. The aim is to determine TEs age on our data and compute their abundance. Comparing abundance graphs between G0 and G12 will inform us if some recent TE insertions happened.


## Data

DNA poolseq data from D.suzukii obtained with Illumina sequencing (2*100pb) :
- One poolseq for G0 
- One poolseq for each media in G12 (strawberry, cherry and cranberry)

| File name    | Number of reads |
| :----------: | :-------------: |
| G0-MTP       | 242 990 252     |
| G12-cerise   | 214 319 262     |
| G12-cranb    | 261 528 888     |
| G12-fraise   | 216 028 586     |


Reference genome of D. suzukii 


## Organization of the steps
The pipeline includes the commands of the following softwares :
- FASTQC, to perform a quality control of the raw data;
- HiSAT2, to perform the alignment of raw data on the reference genome of D. Suzukii;
- GATK, to perform the detection of variants in the data;
- SAMtools, to perform the conversion of SAM files into BAM files, index them and analyze the coverage of mapped reads on the reference genome;
- BCFtools, to compare VCF files;
- Packages R
- PoPoolation TE2
- TEanalysis


## Requirements

### Conda environment


### Versions of tools
- Python: version ;
- FASTQC: version v0.11.9;
- HiSAT2: version 2.1.0;
- GATK: version 4.2.2.0;
- SAMtools: version 1.9;
- BCFtools: version 1.9;
- Packages R: version
- PoPoolation TE2: version
- TEanalysis: version


## Polymoprhism Analysis

- Fastqc on all the files : 
```fastqc input_file_name -o output_file_name``` 
- Index building with hisat2 : 
```hisat2-build -p 4 /data/home/mtabourin/Stage_M1/Drosophila-Suzukii/ref_suzukii/Drosophila-suzukii-contig.fasta index/```
- Mapping with hisat 2 :
- Conversion of .sam files into .bam files :
```samtools view -S -bh ../mapping/G0-MTP.sam | samtools sort -O bam -o G0-MTP-sorted.bam```
- Index BAM files :
```samtools index G0-MTP-sorted.bam```
- Mpileup on .bam files / Call SNPs and short Indels :
```samtools mpileup -uf ref_genome.fa G0-MTP-sorted.bam | bcftools call -mv > variant_calling.vcf```
- Bcftools on .vcf files :
```bcftools filter -s LowQual -e '%QUAL<20 || DP>100' variant_calling.vcf  > variant_calling_filtered.vcf```

## TE analysis


## Authors
The project was developped by:
- Tommaso BARBERIS
- Bertrand HUGUENIN-BIZOT
- Marie VERNERET
- Chloé AUJOULAT


## References
- Anand, Santosh, Eleonora Mangano, Nadia Barizzone, Roberta Bordoni, Melissa Sorosina, Ferdinando Clarelli, Lucia Corrado, Filippo Martinelli Boneschi, Sandra D’Alfonso, et Gianluca De Bellis. « Next Generation Sequencing of Pooled Samples: Guideline for Variants’ Filtering ». Scientific Reports 6, no 1 (27 septembre 2016): 33735. https://doi.org/10.1038/srep33735.
- Chiu, Joanna C, Xuanting Jiang, Li Zhao, Christopher A Hamm, Julie M Cridland, Perot Saelao, Kelly A Hamby, et al. « Genome of Drosophila suzukii, the Spotted Wing Drosophila ». G3 Genes|Genomes|Genetics 3, no 12 (1 décembre 2013): 2257‑71. https://doi.org/10.1534/g3.113.008185.
- Kofler, Robert, Daniel Gómez-Sánchez, et Christian Schlötterer. « PoPoolationTE2: Comparative Population Genomics of Transposable Elements Using Pool-Seq ». Molecular Biology and Evolution 33, no 10 (octobre 2016): 2759‑64. https://doi.org/10.1093/molbev/msw137.
- Olazcuaga, Laure, Julien Foucaud, Mathieu Gautier, Candice Deschamps, Anne Loiseau, Nicolas Leménager, Benoit Facon, et al. « Adaptation and Correlated Fitness Responses over Two Time Scales in Drosophila Suzukii Populations Evolving in Different Environments ». Journal of Evolutionary Biology 34, no 8 (2021): 1225‑40. https://doi.org/10.1111/jeb.13878.
