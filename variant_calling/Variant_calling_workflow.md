# Variant calling Workflow using **GATK**

Following the **GATK** website description, to perform _Variant Calling_ on DNAseq data we need to perform several steps:

* Data Pre-processing
* Variant Discovery
* Filtering variants

## Data Pre-processing

1. Mapping to reference (_Hisat2_)

2. Sorting and adding read group information
A ReadGroup will assign an origin to a set of reads in order to assign a specific genotype to this origin when making the SNP/InDel calling. Without this step, you will have a set of SNPs but you cannot assign them to a specific genotype... This AddOrReplace step is requested by GATK pipeline, as it supposed you will call genotype and not only SNP. If you need only a raw set of SNP, you can use PileUp format and VarScan utility Pileup2SNP.

```
	gatk AddOrReplaceReadGroups 
		-I file.bam \
		-O output.bam \
# the following parameter need to be set after have seen .bam files
		-RGID name \
		-RGLB name \
		-RGPL ILLUMINA \
		-RGPU name \
		-RGSM sample
```

3. Mark duplicates

```
	gatk MarkDuplicatesSpark \
		-I input.bam \
		-O marked_duplicates.bam \
		-M marked_dupmetrics.txt # Path to write duplication metrics to \
#		-R reference \
#		--remove-all-duplicates # if true, do not write duplicate to the output file instead writing them with appropiate flags set \
#		--remove-sequencing-duplicates # if true, do not write optical/sequencing duplicates to the output file instead of writing them with appropriate flags set \
		--spark-master local[6]	# similar to multi-threads \
		--spark-verbosity ERROR # allows to log error \
#		--verbosity ERROR # don't know if it works with the ERROR status, but if it works it can facilitate debugging
```


## SNP calling

1. Create a dictionary of the reference genome

```
	gatk CreateSequenceDictionary \
		-R Drophila-suzukii.fasta # reference genome \
		-O output.dict 
```

2. Index reference genome
_I don't know if this step is necessary, because we have already done the index using Hisat2-build_

```
	samtools faidx reference_genome.fasta
```

3. Index input file
A .bai file isn't an indexed form of a bam - it's a companion to your bam that contains the index. his file has the same name, suffixed with .bai. This file acts like an external table of contents, and allows programs to jump directly to specific parts of the bam file without reading through all of the sequences.

```
	samtools index input.bam
```

4. Detection of variants
The objective is detect variant sites in the dataset.

```
	gatk HaplotypeCaller \
		-R ref.fasta \
		-I file.bam \
		--sample-ploidy [2 x pool size] # in our case, it may be 2 x 40 \
		-O output.vcf
```

5. Create a subset of only SNPs

```
	gatk SelectVariants \
		-R reference_genome.fasta \
		-V input.vcf \
		-O output.vcf \
		--select-type-to-include SNP
```

## SNP filtering
1. Remove low quality SNPs

```
	gatk VariantsToTable \
		-R reference_genome.fasta \
		-V input.vcf \
		-O output.table \
		# Quality scores to save to the table
		-F CHROM \
		-F POS \
		-F QUAL \
		-F DP \
		-F QD \
		-F MQ \
		-F FS \
		-F SOR \
		-F MQRankSum
```

The best way to choose the right filter, is to plot the several quality scores across the VCF file, for exemple using _ggplot2_ from __R__.

```
	gatk VariantFiltration \
		-R reference_genome.fasta \
		-V input.vcf \
		-O output.vcf \
		# following results observed in plots
		-filter "QUAL < 20" --filter-name Low_Qual\
		-filter "DP < 88" --filter-name Low_Cov \
		-filter "QD < 20.0 || MQ < 50.0 || FS > 1.0 || \
		SOR > 1.5 || MQRendSum < -1.0" \
		--filter-name Secondary_filter
```

```
	gatk SelectVariants \
		-R reference_genome.fasta \
		-V input.vcf \
		-O output.vcf \
		--exclude-filtered	
```

## Recalibration
GATK’s guidelines recalibrate the base quality scores after removing duplicates. However, you need a known set of ‘true SNPs’ to do this.
For non-model organisms you will need to create this ‘true SNPs’ dataset from your own data (bootstrapping).

1. Recalibrate the quality scores of each base in read files

```
	gatk BaseRecalibrator \
		-R reference_genome.fasta \
		-I input.bam \ 
		-O output
		--known-sites true_snps.vcf # maybe use the vcf file from RNAseq data
```

```
	gatk ApplyBQSR
		-R reference_genome.fasta \
		-I input.bam \
		--bsqr-recal-file recalibration_table.vcf # output of previous chunk code \
		-O output
```

[Reference](https://yeamanlab.weebly.com/uploads/5/7/9/5/57959825/snp_calling_pipeline.pdf)
