
library("data.table")
gene_annot <- fread("gene_annot.txt.gz")

# GWAS mapping
sumstats <- fread("D:/UCSD/DSC291_stats_genomics/Project/PASS_Alzheimers_deRojas2021.sumstats")
sumstats[, p_value := 2 * (pnorm(abs(Z), lower.tail = F))]
head(sumstats)

# We will be using the positions of HapMap3 variants in HW2 to help determine the position of SNPs
bim_data <- fread("D:/UCSD/DSC291_stats_genomics/HW2_Files/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim")
colnames(bim_data) <- c("chr", "rsID", "cM", "bp", "A1", "A2")
head(bim_data)

common_snps <- intersect(sumstats$SNP, bim_data$rsID)
length(common_snps) # 1015041
length(sumstats$SNP) # 1186517
length(bim_data$rsID) # 1034897

# Merge the datasets, keeping only overlapping SNPs
sumstats_merged <- merge(sumstats, bim_data, by.x = "SNP", by.y = "rsID")
sumstats_merged <- sumstats_merged[order(sumstats_merged$chr, sumstats_merged$bp), ]
head(sumstats_merged)

sig_sumstats_merged <- sumstats_merged[p_value < 5e-8]

## Filter out SNP inside the gene body

snp_genebody <- data.table(SNP = sig_sumstats_merged$SNP, gene = NA)
# Find the nearest gene for each significant SNP
for (i in 1:nrow(sig_sumstats_merged)) {
  snp_chr <- sig_sumstats_merged$chr[i]
  snp_pos <- sig_sumstats_merged$bp[i]
  
  # Find all the genes on the same chromosome
  genes_on_chr <- gene_annot[gene_annot$CHR == snp_chr, ]
  for (j in 1:nrow(regions_on_chr)) {
    gene_start <- genes_on_chr$start[j]
    gene_end <- genes_on_chr$end[j]
  
    if (snp_pos >= gene_start & snp_pos <= gene_end) {
      snp_genebody <- rbind(snp_genebody, data.table(
        SNP = sig_sumstats_merged$SNP,
        gene = genes_on_chr$ID))
    }
  }
}

gene_in_gwas <- na.omit(unique(snp_genebody$gene))
#writeLines(gene_in_gwas, "D:/UCSD/DSC291_stats_genomics/Project/atac/gene_in_gwas.txt")
