# ATAC-seq DEG analysis (mapping with GWAS sumstats data)

# Load required packages
library("Seurat")
library("Signac")
library("devtools")
library("data.table")

atac_seurat <- readRDS("D:/UCSD/DSC291_stats_genomics/Project/ad_atac_seurat.rds")
head(atac_seurat@meta.data)
mg_atac <- subset(atac_seurat, subset = Cell.Type == "MG") # Extract the microglia
table(mg_atac@meta.data$Diagnosis) # Check the distribution
mg_atac <- RunTFIDF(mg_atac)
mg_atac <- FindTopFeatures(mg_atac, min.cutoff = 'q0')
mg_atac <- RunSVD(mg_atac)

# Find differentially accessible peaks between AD and Control
da_peaks_mg <- FindMarkers(
  object = mg_atac,
  group.by = "Diagnosis",
  ident.1 = "AD",
  ident.2 = "Control",
  min.pct = 0.1
)
head(da_peaks_mg)

# Filter the DA peaks results by adjusted p-value
library(dplyr)
fil_da_peaks_mg <- da_peaks_mg %>%
  filter(p_val_adj < 0.05)
nrow(da_peaks_mg) # 1929
nrow(fil_da_peaks_mg) # 1929
head(fil_da_peaks_mg) # visulization


fil_da_peaks_mg$peaks <- rownames(fil_da_peaks_mg)
coords <- strsplit(fil_da_peaks_mg$peaks, split = '-')
coords <- do.call(rbind, coords)
fil_da_peaks_mg$chr <- coords[, 1]
fil_da_peaks_mg$start <- as.numeric(coords[, 2])
fil_da_peaks_mg$end <- as.numeric(coords[, 3])
head(fil_da_peaks_mg)
# write.csv(fil_da_peaks_mg, file = "D:/UCSD/DSC291_stats_genomics/Project/atac/fil_da_peaks_mg.csv", row.names = TRUE)

#################### PART 1 - ATAC-seq DARs to Genes ###########################
# use hg38
library("annotatr")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("GenomicRanges")

# Select annotations for intersection with regions
annots = c('hg38_basicgenes')
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)

fil_da_peaks_mg <- as.data.frame(fil_da_peaks_mg)

da_granges <- makeGRangesFromDataFrame(fil_da_peaks_mg,
                                       keep.extra.columns=FALSE,
                                       ignore.strand=FALSE,
                                       seqinfo=NULL,
                                       seqnames.field=c("seqnames", "seqname",
                                                        "chromosome", "chrom",
                                                        "chr", "chromosome_name",
                                                        "seqid"),
                                       start.field="start",
                                       end.field=c("end", "stop"),
                                       strand.field="strand",
                                       starts.in.df.are.0based=FALSE)

# A GRanges object is returned
head(da_granges)

# Intersect the regions we read in with the annotations
da_annotated = annotate_regions(
  regions = da_granges,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

print(da_annotated)

dar_annotated <- data.frame(da_annotated)
print(head(dar_annotated)) # check and visualization

#write.csv(dar_annotated, file = "D:/UCSD/DSC291_stats_genomics/Project/atac/dar_annotated.csv", row.names = TRUE)

all_dar_genes <- na.omit(unique(dar_annotated$gene_name))
#writeLines(all_dar_genes, "D:/UCSD/DSC291_stats_genomics/Project/atac/all_dar_genes.txt")

#### DAR - promoter genes
# Filtering rows with 'promoters', selecting the 'annot.symbol' column
da_promoter_genes <- dar_annotated %>%
  dplyr::filter(grepl("promoters", annot.type)) %>%
  dplyr::select(gene_name = annot.symbol)
da_promoter_genes <- unique(da_promoter_genes$gene_name)
da_promoter_genes <- na.omit(da_promoter_genes)
nrow(da_promoter_genes) # 1755
head(da_promoter_genes)
#writeLines(da_promoter_genes, "D:/UCSD/DSC291_stats_genomics/Project/atac/da_promoter_genes.txt")

#### DAR - 1to5kb genes
# Filtering rows with '1to5kb', selecting the 'annot.symbol' column
da_1to5kb_genes <- dar_annotated %>%
  dplyr::filter(grepl("1to5kb", annot.type)) %>%
  dplyr::select(gene_name = annot.symbol)
da_1to5kb_genes <- unique(da_1to5kb_genes$gene_name)
da_1to5kb_genes <- na.omit(da_1to5kb_genes)

nrow(da_1to5kb_genes) # 1743
head(da_1to5kb_genes)
#writeLines(da_1to5kb_genes, "D:/UCSD/DSC291_stats_genomics/Project/atac/da_1to5kb_genes.txt")

## # NOTE: The following code is not included in our Analysis Report; it is only used for internal checks.
#################### PART 2 - GWAS significant SNP ###########################

## GWAS sumstats
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

sumstats_merged$chr <- paste0("chr", sumstats_merged$chr) # make sure format align with the DAR dataframe

## Use GWAS significant SNPs for mapping
sig_sumstats_merged <- sumstats_merged[p_value < 5e-8]

## Filter out SNP included DA Regions

sig_snp_da_peaks <- data.table(SNP = character(), p_val = numeric(), avg_log2FC = numeric(),
                               pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(), 
                               peaks = character(), chr = character(), start = numeric(), end = numeric())

# Loop through each significant SNP
for (i in 1:nrow(sig_sumstats_merged)) {
  snp_chr <- sig_sumstats_merged$chr[i]
  snp_pos <- sig_sumstats_merged$bp[i]
  
  # Screening ATAC-seq regions of identical chr
  regions_on_chr <- fil_da_peaks_mg[fil_da_peaks_mg$chr == snp_chr, ]
  
  if (nrow(regions_on_chr) == 0) next
  
  # Screening for SNPs in a region
  for (j in 1:nrow(regions_on_chr)) {
    region_start <- regions_on_chr$start[j]
    region_end <- regions_on_chr$end[j]
    
    # if the SNP is in the region
    if (snp_pos >= region_start & snp_pos <= region_end) {
      # If inside, add the information to the results table
      sig_snp_da_peaks <- rbind(sig_snp_da_peaks, data.table(
        SNP = sig_sumstats_merged$SNP[i],
        p_val = regions_on_chr$p_val[j],
        avg_log2FC = regions_on_chr$avg_log2FC[j],
        pct.1 = regions_on_chr$pct.1[j],
        pct.2 = regions_on_chr$pct.2[j],
        p_val_adj = regions_on_chr$p_val_adj[j],
        peaks = regions_on_chr$peaks[j],
        chr = snp_chr,
        start = region_start,
        end = region_end
      ))
    }
  }
}

nrow(sig_snp_da_peaks) # 6
head(sig_snp_da_peaks)

# Remove duplicate peaks column
sig_snp_da_peaks_unique <- sig_snp_da_peaks[!duplicated(sig_snp_da_peaks$peaks), ]
sig_snp_da_peaks_unique[, SNP := NULL] # remove the SNP column
rownames(sig_snp_da_peaks_unique) <- sig_snp_da_peaks_unique$peaks
nrow(sig_snp_da_peaks_unique) # 2

#write.csv(sig_snp_da_peaks_unique, file = "D:/UCSD/DSC291_stats_genomics/Project/atac/sig_snp_da_peaks_unique.csv", row.names = TRUE)

## DA peaks to gene names
# use hg38
library("annotatr")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("GenomicRanges")

# Select annotations for intersection with regions
annots = c('hg38_basicgenes')
# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)
sig_snp_da_peaks_unique <- as.data.frame(sig_snp_da_peaks_unique)
sig_da_granges <- makeGRangesFromDataFrame(sig_snp_da_peaks_unique,
                                           keep.extra.columns=FALSE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("seqnames", "seqname",
                                                            "chromosome", "chrom",
                                                            "chr", "chromosome_name",
                                                            "seqid"),
                                           start.field="start",
                                           end.field=c("end", "stop"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE)

head(sig_da_granges) # A GRanges object is returned

sig_da_annotated = annotate_regions(regions = sig_da_granges, annotations = annotations, 
                                    ignore.strand = TRUE, quiet = FALSE)

sig_df_da_annotated <- data.frame(sig_da_annotated)
print(head(sig_df_da_annotated)) # check and visualization

sig_df_da_genes <- sig_df_da_annotated %>%
  dplyr::select(gene_name = annot.symbol)

sig_unique_da_genes <- unique(sig_df_da_genes$gene_name)
length(sig_unique_da_genes) # 1 genes
head(sig_unique_da_genes)

#write.csv(sig_df_da_annotated, file = "D:/UCSD/DSC291_stats_genomics/Project/atac/sig_df_da_annotated.csv", row.names = TRUE)
#writeLines(sig_unique_da_genes, "D:/UCSD/DSC291_stats_genomics/Project/atac/sig_unique_da_genes.txt")