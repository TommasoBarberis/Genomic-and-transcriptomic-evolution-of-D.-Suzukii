library(vcfR)

DATA_DIR <- "/data/home/tbarberis/projet_S3/localdata/variant_calling_mpileup/vcf/"
setwd(DATA_DIR)

G0_MTP <- read.vcfR("G0-MTP.mpileup.vcf.gz", verbose = FALSE)
G0_MTP.gt <- extract.gt(G0_MTP)
