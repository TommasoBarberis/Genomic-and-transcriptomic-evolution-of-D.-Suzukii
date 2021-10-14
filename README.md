# Genomic-and-transcriptomic-evolution-of-D.-Suzukii


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

Reference genome of D. suzukii 

## Polymoprhism Analysis

- Fastqc on all the files 
- Index building with hisat2 : 
```hisat2-build -p 4 /data/home/mtabourin/Stage_M1/Drosophila-Suzukii/ref_suzukii/Drosophila-suzukii-contig.fasta index/```
- Mapping with hisat 2


## TE analysis




