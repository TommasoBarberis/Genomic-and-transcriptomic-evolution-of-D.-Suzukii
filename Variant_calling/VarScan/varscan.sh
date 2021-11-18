#!/bin/bash -i

# Author: Tommaso Barberis
# Date: 13/11/2021
# Variant calling using varscan on mpileup files generated using samtools (see mpileup.sh script)

## samples name
declare -a samples
samples=(
#	'G0-MTP.mpileup'
#	'G12-cerise.mpileup'
#	'G12-fraise.mpileup'
	'G12-cranb.mpileup'
	)

## variant calling
for sample in ${samples[@]}; do
	echo -e "calling $sample \n\n"

	start_time=$(date +%s)
    
	conda activate dsuzukii
	varscan mpileup2cns mpileup/$sample --min-coverage 50 --min-avg-qual 20 --min-var-freq 0.001 --variants --output-vcf > $sample.vcf
	conda deactivate

	end_time=$(date +%s)
	elapsed=$(( end_time - start_time ))
	display_elapsed=$(eval "echo $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')")
	
	mailSender "${sample} called. Elapsed time: ${display_elapsed}.Send from pedago-ngs" tommasobarberis98@gmail.com "varscan ${sample} ended" 

done
