# Analysis of VCF files

## Some fields on which we can filter our variants

- QUAL - quality: Phred-scaled quality score for the assertion made in ALT. 
	i.e. −10log10 prob(call in ALT is wrong). If ALT is ‘.’ (no variant) then this is −10log10 prob(variant), and if ALT is not ‘.’ this is −10log10 prob(no variant). If unknown, the missing value should be specified. (Numeric)
- MQ: RMS mapping quality, e.g. MQ=52
- MQRankSum (MappingQualityRankSumTest): This is the u-based z-approximation from the Rank Sum Test for mapping qualities. It compares the mapping qualities of the reads supporting the reference allele and the alternate allele. A positive value means the mapping qualities of the reads supporting the alternate allele are higher than those supporting the reference allele; a negative value indicates the mapping qualities of the reference allele are higher than those supporting the alternate allele. A value close to zero is best and indicates little difference between the mapping qualities.
- ReadPosRankSum: This is the u-based z-approximation from the Rank Sum Test for site position within reads. It compares whether the positions of the reference and alternate alleles are different within the reads. Seeing an allele only near the ends of reads is indicative of error, because that is where sequencers tend to make the most errors. A negative value indicates that the alternate allele is found at the ends of reads more often than the reference allele; a positive value indicates that the reference allele is found at the ends of reads more often than the alternate allele. A value close to zero is best because it indicates there is little difference between the positions of the reference and alternate alleles in the reads.

## GATK steps 

- Create a dictionary for the genome reference:
```
gatk CreateSequenceDictionary -R Drophila-suzukii.fasta
```

- Combine all `vcf` files together:
```
gatk CombineGVCFs -V G0_HaplotypeCall.g.vcf.gz -V G12_cerise.g.vcf.gz -V G12_fraise_HaplotypeCall.g.vcf.gz -V G12-cranb.g.vcf.gz -R ../ref/Drosophila-suzukii-contig.fasta -O poolseq_dsuzukii.g.vcf.gz
```

- Joint genotype call:
```
GenotypeGVCFs -R ../../../ref/Drosophila-suzukii-contig.fasta -V ../merge/merge_samples.g.vcf.gz -O dsuzukii_gt.g.vcf.gz
```

- Filter by variant type:
```
gatk SelectVariants -V genotype/dsuzukii_gt.g.vcf.gz -select-type SNP -O select_variants/dsuzukii_snp.vcf.gz

gatk SelectVariants -V genotype/dsuzukii_gt.g.vcf.gz -select-type INDEL -O select_variants/dsuzukii_indel.vcf.gz

gatk SelectVariants -V genotype/dsuzukii_gt.g.vcf.gz -select-type MIXED -O select_variants/dsuzukii_mixed.vcf.gz
```



Source: [GATK Broad Institute](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)
