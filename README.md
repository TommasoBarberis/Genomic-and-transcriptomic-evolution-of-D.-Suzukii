# Genomic and transcriptomic evolution of _D. suzukii_

## Table of contents

1. [Context](#1-context)
2. [Introduction to the project](#2-introduction-to-the-project)
3. [Data](#3-data)
4. [Workflow](#4-workflow)
5. [Authors](#4-authors)
6. [References](#5-references)

## 1. Context

Since its recent invasion in the European and American continents, _Drosophila suzukii_, has become a major pest of berry crops. 
Contrary to other _Drosophila_ species that develop themselves on damaged or rotten fruit, _D. suzukii_ is able to lay its eggs in healthy fruit before harvest, using its __sclerotinized ovipositor__. In order to control the _D. suzukii_ population it is therefore required to improve our knowledge about its ability to adapt to a local environment outside its origin areas.

A previous study focused on phenotipic changes during a selective experiment showed a temporal adaptation of _D. suzukii_ populations in the three different fruit environments. This process produces a __fitness__ improvement which corresponds to a greater ability to reproduce and could explain the adaptation abilities  of _D. suzukii_. Based on the same experiment, a second study has been started by the LBBE to see if the _D. suzukii_ genome and transcriptome are affected by the environment.
 
The project is focused on both genes and Transposable Elements (TE) expression in order to find the variants explaining the local environment adaptation of _D. suzukii_. 
A significant number of differentially expressed genes across generations has been found with __RNAseq analysis__ between the G0 and the G7 generations. These variations seems to be independant of the fruit media. Moreover, in differentially expressed (DE) genes an increase of the number of Single Nucleotide Polymorphism (SNP) was found in G7 compared to G0, but this result was not expected. Natural selection generally keeps alleles which allow the best adaptation, and it is synonymous with a decrease in SNP number. 


## 2. Introduction to the project

### 2.1 Polymorphism analysis
The goal here is to analyse _D. suzukii_ polymorphism performing a __variant calling__ on the whole genome (using poolseq data) in order to count the number of different genomic positions between G0 and G12, and between fruits media. <br/>
The position of the genetic variations will be compared with previous SNP analysis conducted on RNA-seq data to see if we also find an increase of the SNP number. To have a view more global on the genome's diversity, several diversity measures will be calculated: __pi's__ $\theta$, __Watterson's__ $\theta$, __Tajima's D__.

### 2.2 TE analysis
__Transposable elements__ (__TE__) represent __47%__ of the _D. suzukii_ genome and can be important in the adaptation process, interfering in gene regulation during evolution. The aim is to determine TEs age on our data and compute their abundance. Comparing abundance graphs between G0 and G12 will inform us if some recent TE insertions happened.


## 3. Data

- DNA poolseq data from _D.suzukii_ obtained with paired-end Illumina sequencing (2x100pb), with 40 individuals for each pool:
    - One poolseq for G0;
    - One poolseq for each media in G12 (strawberry, cherry and cranberry).

<div align="center">

| File name    | Number of reads |
| :----------: | :-------------: |
| G0-MTP       | 242 990 252     |
| G12-cerise   | 214 319 262     |
| G12-cranb    | 261 528 888     |
| G12-fraise   | 216 028 586     |

</div>

Location: _pedago-ngs_ 
```
/localdata/pandata/students/M2_projet_15/data
```

- Partial assembly of _D. suzukii_ genome as reference (4 chromosomes, 542 contigs);

Location: _pedago-ngs_ 
```
/localdata/pandata/students/M2_projet_15/ref/Drosophila-suzukii-contig.fasta
```

- Annotation file at the format `gff3`.
Location: _pedago-ngs_ 
```
/localdata/pandata/students/M2_projet_15/ref/Drosophila-suzukii-annotation-3-ws.gff3
```

## 4. Workflow

### 4.1 Quality control

First of all we checked the quality of the data, using `fastqc`:
```
fastqc sample.fastq.gz
```

Location: _pedago-ngs_ 
```
/localdata/pandata/students/M2_projet_15/quality_control
```

Then, we merged results from `fastqc` using `multiqc` in the folder with `fastqc` reports. The `html` report is here in the `quality_control` folder. <br>
We noticed a high __GC__ rate in the reads, but it seems to be in accordance with the results from other research on _D. suzukii_.

### 4.2 Mapping

1. Genome indexing:
```
bwa index \
    -p Drosophila-suzukii-contig.fasta.fai \
    -a bwtsw Drosophila-suzukii-contig.fasta
```

- `-p`: index name;
- `-a`: index algorithm (_bwtsw_ for long genome).

2. Mapping samples on the reference:
```
bwa aln -t 14 Drosophila-suzukii-contig.fasta.fai sample_1.fastq.gz > sample_1.sai
bwa aln -t 14 Drosophila-suzukii-contig.fasta.fai sample_2.fastq.gz > sample_2.sai
```
- `-t`: number of threads.

3. Merge of the `.sai` files from forward and reverse reads:
```
bwa sampe Drosophila-suzukii-contig.fasta.fai \
    sample_1.sai sample_2.sai \
    sample_1.fastq.gz sample_2.fastq.gz > sample_pe.sam
```
- `sampe`: take in account to have paired ends reads.

#### Mapping stats

Calculation of statistics of mapped files using:
```
samtools stats sample.sam
```

| sample | raw sequences | mapped and paired squences |
| :------: | :-------------: | :--------------------------: |
| G0-MTP_pe | 242495252 | 194326801 |
| G12-fraise_pe | 216028586 | 172898366 |
| G12-cerise_pe | 214319262 | 172440052 |
| G12-cranb_pe | 261528888 | 212135984 |

#### Conversion to `.bam`
```
samtools view -bS -@ 3 sample.sam > sample.bam
```
- `-bS`: output in the __BAM__ format ignoring for compatibility with previous samtools version;
- `-@`: number of threads.

Location: _pedago-ngs_ 
```
/localdata/pandata/students/M2_projet_15/mapping_bwa
```

__NB:__ Before using `bwa`, we tried to use `hisat2` but, for an unkown reason it could not map the strawberry sample.


### 4.3 Variant Calling
The SNP calling on pooled data is still not largely used because it is very innovative. For this reason has required several tests (with different tools) and long time to be completed. <br>
Finally, we used `GATK` to perform our SNP calling and for that we were inspired by M. Tabourin's workflow ([Analyse experience d'évolution expérimentale _D.suzukii_](https://github.com/mtabourin/Analyse-experience-d-evolution-experimentale-D.suzukii)) and from [this](https://yeamanlab.weebly.com/uploads/5/7/9/5/57959825/snp_calling_pipeline.pdf) document.

1. Sorting `.bam` files:

```
# GATK version: gatk-4.2.2.0
gatk SortSam -I sample.bam -O sample_sorted.bam -SO queryname -VALIDATION_STRINGENCY SILENT
```
    
- `-I`: input file;
- `-O`: output file;
- `-SO`: sort order of output file;
- `-VALIDATION_STRINGENCY`: Validation stringency for all SAM files read by this program. Setting stringency to `SILENT` can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

2. Add read group info:
```
gatk AddOrReplaceReadGroups \
    -I sample_sorted.bam \
    -O sample_ReadGroups.bam \
    -RGID name \
    -RGLB name \
    -RGPL ILLUMINA \
    -RGPU name \
    -RGSM sample \
    -VALIDATION_STRINGENCY SILENT
```

- `-I`: input file;
- `-O`: output file;
- `-RGID`: Read-Group ID;
- `-RGLB`: Read-Group library;
- `-RGPL`: Read-Group platform;
- `-RGPU`: Read-Group platform unit(barcode);
- `-RGSM`: Read-Group sample name;
- `-VALIDATION_STRINGENCY`: Validation stringency for all SAM files read by this program.

__NB:__ Read group information flags (`-RG*`) are mandatory, even if you don't have them.

3. Remove duplicate:
```
gatk MarkDuplicatesSpark \
    -I sample_ReadGroups.bam \
    -O sample_MarkDuplicated.bam \
    --spark-master local[16]
```
- `-I`: input file;
- `-O`: output file;
- `--spark-master`: to parallelize the work.

4. Detect variant sites in the dataset:
```
gatk HaplotypeCaller \
    -R reference.fasta \
    -I sample_MarkDuplicated.bam \
    --sample-ploidy 80 -O sample.vcf \
    --max-genotype-count 91881
```
- `-R`: genome reference;
- `-I`: input file;
- `-O`: output file (`vcf` format);
- `--sample-ploidy`: it allows to take in account of using pooled data:
$sample\_ploidy = organism\_ploidy \times pool\_size$ <br>
In our case we have a ploidy of 2 (_D. suzukii_ is a diploïd organism) and a pool size of 40.
- `--max-genotype-count`: maximum number of genotypes to consider at any site:
$max\_genotype\_count = \binom{P +  A - 1}{A - 1} = \frac{(P + A - 1)!}{P!(A - 1)!}$ <br>
With `P` la ploidy of the sample (previous formula) and `A` ($A = 3$) the allele count. 

#### 4.3.1 SNPs filtering

1. First of all, we have make a subset having only __SNPs__:
```
gatk SelectVariants \
    -R reference_genome.fasta \
    -V input.vcf \
    -O output.vcf \
    --select-type-to-include SNP
```
- `-R`: genome reference;
- `-V`: `.vcf` input file;
- `-O`: output file;
- `--select-type-to-include` select a type of variant.

2. Generate a table that can be used to find filter threshold:
```
	gatk VariantsToTable \
		-R reference_genome.fasta \
		-V input.vcf \
		-O output.table \
		# scores to save to the table
		-F CHROM \
		-F POS \
		-F QUAL \
		-F AC \
		-F AF \
		-F DP \
		-F QD \
		-F MQ \
		-F FS \
		-F SOR \
		-F MQRankSum \
		-F ReadPosRankSum
```
- `-R`: reference;
- `-V`: `.vcf` input file;
- `-O`: output file;
- `-F`: selected field.

The best way to choose the right filter, is to plot the several scores across the VCF file, for example using `ggplot2` library in `R`.
You can find here the plot for some scores: <br/>

<img src="./img/DP.svg"/>
<img src="./img/FS.svg"/>
<img src="./img/MQ.svg"/>
<img src="./img/QD.svg"/>
<img src="./img/QUAL.svg"/>
<img src="./img/SOR.svg"/>

<br/>

3. Use filter on the `vcf` file:
```
	gatk VariantFiltration \
		-R reference_genome.fasta \
		-V input.vcf \
		-O output.vcf \
		# following results observed in plots
		-filter "QUAL < 70" --filter-name Low_Qual\
		-filter "DP < 90" --filter-name Low_Cov \
		-filter "QD < 20.0 || MQ < 28.0 || FS > 1.0 || \
		SOR > 1.4 --filter-name Secondary_filter
```
- `-R`: reference;
- `-V`: `vcf` input file;
- `-O`: output file;
- `-filter`: define a filter;
- `--filter-name`: used to assign a name to the filter.

4. Extract __SNPs__ having pass filter:
```
gatk SelectVariants \
		-R reference_genome.fasta \
		-V input.vcf \
		-O output.vcf \
		--exclude-filtered
```
- `-R`: reference;
- `-V`: `vcf` input file;
- `-O`: output file;
- `--exclude-filtered`: exclude __SNPs__ marked to be filtered.

#### 4.3.2 Get __SNPs__ that are in transcripts (mRNA)

1. From the `.gff3` file of the reference, we can extract a subset having only mRNA:
```
tail -n +2 Drosophila-suzukii-annotation-3-ws.gff3 | \
awk -F "\t" '{if($3=="mRNA") print $0}' > ~/mRNA.gff3
```
- `tail -n +2`: to skip `gff3` header;
- `awk`: if the third field has 'mRNA' tag, it print the line in the `mRNA.gff3` file.

2. Convert the `.gff3` into `.bed` format:
```
gff2bed < mRNA.gff3 > mRNA.bed
```

3. Intersect the `.vcf` file with the `.bed` file:
```
intersectBed -a input.vcf -b mRNA.bed -header > output.vcf
```
- `-a`: `.vcf` input file;
- `-b`: `.bed` input file;
- `-header`: if the `.bed` file has header.

#### 4.3.3  Get __SNPs__ that are in DE genes
List of the DE genes can be found here: <br/>

Location: _pedago-ngs_ 
```
/localdata/pandata/students/M2_projet_15/DE_gene
```
This list is obtained by the M. Tabourin's work.

1. Obtain a `.gff3` file containing only DE genes:
```
python3 get_de.py
```
- `get_de.py`: it is a script written by us. It allows, from the `gff3` and the `list.txt` having DE genes, to obtain a `.gff3` file having only DE genes.

2. Then we have transformed the obtained `.gff3` file into `.bed` file using the command `3` at the previous paragraph and we have used it to obtain __SNPs__ in DE genes using command `4` of the previous paragraph.

#### 4.3.4 Results

Location of `.vcf` result files: _pedago-ngs_
```
/localdata/pandata/students/M2_projet_15/GATK_pipeline/gatk_ploidy
```

| | Total number of SNPs | SNPs on mRNA | SNPs on DE genes |
| :-: | :-: | :-: | :-: |
| G0-MTP | 189,567 | 121,646 | 1,465 |
| G12-strawberry | 175,067 | 111,761 | 1,264 |
| G12-cherry | 186,952 | 115,760 | 1,242 |
| G12-cranberry | 116,006 | 75,305 | 1,050 | 


### Comparison of our results from DNA data with RNA seq data

1. Compress VCF files: ```bgzip input_DNA_file.vcf``` and ```bgzip input_RNA_file.vcf```
2. Build indexes for compressed files: ```tabix -p vcf input_DNA_file.vcf.gz``` and ```tabix -p vcf input_RNA_file.vcf.gz```
3. Run bcftools isec to make the intersection between the files of interest: ```bcftools isec -c snps -p /fraise input_DNA_file.vcf.gz input_RNA_file.vcf.gz path2repository```

#### Results

Location of 4 output result files: _pedago-ngs_
```
/localdata/pandata/students/M2_projet_15/GATK_pipeline/gatk_ploidy/$fruit
```
- First file: 0000.vcf for records private to  fraise/G12_fraise_MarkDuplicated.bam.vcf.gz
- Second file: 0001.vcf for records private to  variants_RNAseq.vcf.gz
- Third file: 0002.vcf for records from fraise/G12_fraise_MarkDuplicated.bam.vcf.gz shared by both     fraise/G12_fraise_MarkDuplicated.bam.vcf.gz variants_RNAseq.vcf
- Fourth file: 0003.vcf for records from variants_RNAseq.vcf.gz shared by both  fraise/G12_fraise_MarkDuplicated.bam.vcf.gz variants_RNAseq.vcf.gz


### TE analysis
All the steps are included in a tool named dnaPipeTE (available on :https://github.com/clemgoub/dnaPipeTE)
- Uniform samplings of the reads to produce low coverage data sets 
- Trinity, repeated assembly of contigs
- RepeatMasker, contigs annotation
- Blastn, repeat quantification


### Versions of tools
- FASTQC: version v0.11.9;
- HiSAT2: version 2.1.0;
- GATK: version 4.2.2.0;
- SAMtools: version 1.9;
- BCFtools: version 1.9;
- dnaPipeTE: version 1.3 (uses Perl 5, R v3.0.2, Python v3.8.5, Trinity v2.5.1, RepeatMasker v4.0.5 including RMblastn, ncbi-blast v2.2.28+)


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
- Mpileup on .bam files and bcftools for call SNPs and short Indels :

```samtools mpileup sam_to_bam/G0-MTP-sorted.bam -f /data/home/mtabourin/Stage_M1/Drosophila-Suzukii/ref_suzukii/Drosophila-suzukii-contig.fasta > variant_calling_mpileup/G0-MTP.pileup```

```bcftools mpileup -o G0-MTP-VCFfile.vcf -f /data/home/mtabourin/Stage_M1/Drosophila-Suzukii/ref_suzukii/Drosophila-suzukii-contig.fasta ../sam_to_bam/G0-MTP-sorted.bam```

```bcftools call -vmO v -o bcftools_call/G0-MTP-VCFfile-call.vcf variant_calling_mpileup/G0-MTP-VCFfile.vcf```

- Bcftools filter on .vcf files :
```bcftools filter -O v -o bcftools_filtered/G0-MTP-VCFfile-filtered.vcf -s LOWQUAL -i 'QUAL>10' -g3 -G10 bcftools_call/G0-MTP-VCFfile-call.vcf```

## TE analysis
- dnaPipeTE on G0 file :
```python3 ./dnaPipeTE.py -input /home/ubuntu/data/G0-MTP_1.fastq.gz -output /home/ubuntu/output_dnaPipeTE/G0 -cpu 6 -sample_number 2 -genome_size 270000000 -genome_coverage 0.2 -RM_lib /home/ubuntu/Drosophila_Transposable_Element_all.fasta -keep_Trinity_output```

Sampling 2 samples of max 355263 reads to reach 0.2X coverage of the 270Mb genome of D.suzukii :
For G0 fastq :
total number of reads: 121247626
maximum number of reads to sample:  710526



## 4. Authors
The project was developped by:
- Chloé AUJOULAT
- Tommaso BARBERIS
- Bertrand HUGUENIN-BIZOT
- Marie VERNERET


## 5. References
- Anand, Santosh, Eleonora Mangano, Nadia Barizzone, Roberta Bordoni, Melissa Sorosina, Ferdinando Clarelli, Lucia Corrado, Filippo Martinelli Boneschi, Sandra D’Alfonso, et Gianluca De Bellis. « Next Generation Sequencing of Pooled Samples: Guideline for Variants’ Filtering ». Scientific Reports 6, no 1 (27 septembre 2016): 33735. https://doi.org/10.1038/srep33735.
- Chiu, Joanna C, Xuanting Jiang, Li Zhao, Christopher A Hamm, Julie M Cridland, Perot Saelao, Kelly A Hamby, et al. « Genome of Drosophila suzukii, the Spotted Wing Drosophila ». G3 Genes|Genomes|Genetics 3, no 12 (1 décembre 2013): 2257‑71. https://doi.org/10.1534/g3.113.008185.
- Kofler, Robert, Daniel Gómez-Sánchez, et Christian Schlötterer. « PoPoolationTE2: Comparative Population Genomics of Transposable Elements Using Pool-Seq ». Molecular Biology and Evolution 33, no 10 (octobre 2016): 2759‑64. https://doi.org/10.1093/molbev/msw137.
- Olazcuaga, Laure, Julien Foucaud, Mathieu Gautier, Candice Deschamps, Anne Loiseau, Nicolas Leménager, Benoit Facon, et al. « Adaptation and Correlated Fitness Responses over Two Time Scales in Drosophila Suzukii Populations Evolving in Different Environments ». Journal of Evolutionary Biology 34, no 8 (2021): 1225‑40. https://doi.org/10.1111/jeb.13878.
- Clément Goubert, Laurent Modolo, Cristina Vieira, Claire ValienteMoro, Patrick Mavingui, Matthieu Boulesteix, De Novo Assembly and Annotation of the Asian Tiger Mosquito (Aedes albopictus) Repeatome with dnaPipeTE from Raw Genomic Reads and Comparative Analysis with the Yellow Fever Mosquito (Aedes aegypti), Genome Biology and Evolution, Volume 7, Issue 4, April 2015, Pages 1192–1205, https://doi.org/10.1093/gbe/evv050
- Marie  Tabourin  et  al.  “Analyse  de  l’expression  des  gènes  et  des  éléments  transposableschez Drosophila suzukii suite à une expérience d’évolution expérimentale.
- Simon Andrews et al.FastQC. Babraham Institute. Babraham, UK, Jan. 2012, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.
- MultiQC: Summarize analysis results for multiple tools and samples in a single report, Philip Ewels, Måns Magnusson, Sverker, Lundin and Max Käller, Bioinformatics (2016), doi: 10.1093/bioinformatics/btw354, PMID: 27312411.
- Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168].
- Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352.
- Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.
- Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
- R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
- Bedops.readthedocs.io. 2021. 6.3.3.5. gff2bed — BEDOPS v2.4.40. (online) Available at: https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html (Accessed 2 December 2021).
- Bedtools.readthedocs.io. 2021. intersect — bedtools 2.30.0 documentation. (online) Available at: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html (Accessed 2 December 2021).