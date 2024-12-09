# scRNA-seq DEG analysis

library(Seurat)
library(devtools)
library(dplyr)
library(patchwork)
library(ggplot2)

rna_seurat <- readRDS("D:/UCSD/DSC291_stats_genomics/Project/ad_rna_seurat.rds")
head(rna_seurat@meta.data)

rna_seurat <- subset(rna_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
mg_rna <- subset(rna_seurat, subset = Cell.Type == "MG")
rm(rna_seurat)
head(mg_rna@meta.data)
table(mg_rna@meta.data$Diagnosis) # check

mg_rna <- NormalizeData(mg_rna)
mg_rna <- FindVariableFeatures(mg_rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mg_rna)
mg_rna <- ScaleData(mg_rna, features = all.genes)
mg_rna <- RunPCA(mg_rna, features = VariableFeatures(object = mg_rna))

deg_mg <- FindMarkers(mg_rna, group.by = "Diagnosis", ident.1 = "AD", ident.2 = "Control")
nrow(deg_mg) # 10113

fil_deg_mg <- deg_mg %>%
  filter(abs(avg_log2FC) > 0.25 & p_val_adj < 0.05)
nrow(fil_deg_mg) # 245
head(fil_deg_mg)
filtered_DEG_mg <- unique(na.omit(rownames(fil_deg_mg)))
#writeLines(filtered_DEG_mg, "D:/UCSD/DSC291_stats_genomics/Project/rna/filtered_DEG_mg.txt")