### This file contains the code to perform S-LDSC analysis
### Baseline model required file is from https://doi.org/10.5281/zenodo.7768713
### sumstats is from PMC8184987

ldsc_path=/d/ucsd/dsc291/ldsc/

## Here is the R code to get p-value based on z score
# library(dplyr)
# library(data.table)
# sumstats <- read.table("./PASS_Alzheimers_deRojas2021.sumstats", 
#                        header = TRUE, 
#                        sep = "\t", 
#                        stringsAsFactors = FALSE)
# # add p value
# sumstats <- sumstats %>%
#   mutate(P = 2*pnorm(q=abs(sumstats$Z),lower.tail=F))
# head(sumstats)
# # save
# write.table(sumstats, "PASS_Alzheimers_deRojas2021.sumstats.withP", 
#             quote = FALSE, row.names = FALSE, sep = "\t")

## keep the SNPs in HapMap3
python $ldsc_path/munge_sumstats.py \
  --sumstats PASS_Alzheimers_deRojas2021.sumstats.withP \
  --merge-alleles w_hm3.snplist \
  --chunksize 500000 \
  --out filtered \
  --a1-inc
## Perform the S-LDSC based on baseline model
python $ldsc_path/ldsc.py \
  --h2 d:/ucsd/dsc291/filtered.sumstats.gz \
  --overlap-annot \
  --ref-ld-chr d:/ucsd/dsc291/data_downloaded/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. \
  --frqfile-chr d:/ucsd/dsc291/data_downloaded/1000G_Phase3_frq/1000G_Phase3_frq/1000G.EUR.QC. \
  --w-ld-chr d:/ucsd/dsc291/data_downloaded/1000G_Phase3_weights_hm3_no_MHC/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --out alz_baseline