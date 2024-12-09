# overlap with scRNA-seq DEGs

library(ggvenn)

filtered_DEG_mg <- read.table("D:/UCSD/DSC291_stats_genomics/Project/rna/filtered_DEG_mg.txt", header = F)
all_dar_genes <- read.table("D:/UCSD/DSC291_stats_genomics/Project/atac/all_dar_genes.txt", header = F)
da_promoter_genes <- read.table("D:/UCSD/DSC291_stats_genomics/Project/atac/da_promoter_genes.txt", header = F)
sig_unique_da_genes <- read.table("D:/UCSD/DSC291_stats_genomics/Project/atac/sig_unique_da_genes.txt", header = F)
gene_in_gwas <- read.table("D:/UCSD/DSC291_stats_genomics/Project/atac/gene_in_gwas.txt", header = F)

filtered_DEG_mg <- as.character(filtered_DEG_mg$V1)
all_dar_genes <- as.character(all_dar_genes$V1)
da_promoter_genes <- as.character(da_promoter_genes$V1)
sig_unique_da_genes <- as.character(sig_unique_da_genes$V1)
gene_in_gwas <- as.character(gene_in_gwas$V1)

# All DAR mapped genes
venn_all_dar <- list(true_deg = filtered_DEG_mg, all_dar_mapped_gene = all_dar_genes)
ggvenn(venn_all_dar, fill_color = c("blue", "green"), stroke_size = 0.5, set_name_size = 4)

# promoter mapped genes
venn_dar_promoter <- list(true_deg = filtered_DEG_mg, dar_promoter_mapped_gene = da_promoter_genes)
ggvenn(venn_dar_promoter, fill_color = c("blue", "green"), stroke_size = 0.5, set_name_size = 4)

# gene body mapped genes
venn_gene_gwas <- list(true_deg = filtered_DEG_mg, gwas_inside_gene_body = gene_in_gwas)
ggvenn(venn_gene_gwas, fill_color = c("blue", "green"), stroke_size = 0.5, set_name_size = 4)

# gwas in DAR mapped genes
venn_gwas_dar <- list(true_deg = filtered_DEG_mg, snp_inside_dar_gene = sig_unique_da_genes)
ggvenn(venn_gwas_dar, fill_color = c("blue", "green"), stroke_size = 0.5, set_name_size = 4)

# bar chart
library(ggplot2)

bar_data <- data.frame(
  Category = c("All DAR overlapped DEG", "Promoter overlapped DEG", "GWAS SNP mapped DEG"),
  Values = c(36, 35, 2)
)

bar_data$Category <- factor(bar_data$Category, levels = bar_data$Category[order(-bar_data$Values)])

ggplot(bar_data, aes(x = Category, y = Values, fill = Category)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("skyblue", "lightgreen", "lightpink")) + 
  labs(title = "Overlap of DEG with Different Categories",
       x = "Category", y = "Number of Overlapped Genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1), 
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, margin = margin(b = 13)) 
  ) + 
  geom_text(aes(label = Values), vjust = -0.5) +
  ylim(0, max(bar_data$Values) + 3) 