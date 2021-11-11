samtools sort G0-MTP_pe.bam -o G0-MTP_pe_sorted.bam -@ 6
samtools faidx Drosophila-suzukii-contig.fasta
samtools mpileup -g -f ../index/Drosophila-suzukii-contig.fasta G0-MTP_pe_sorted.bam > G0-MTP_pe.bcf
