# DSC291 Final Project

setwd("./~")
## Step1 
## Creat the cell-by-genes and cell-by-peaks matrix using Seurat and Signac packages

library(Seurat)
library(hdf5r)
# scRNA-seq
rna_data <- Read10X_h5("GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
rna_meta <- read.csv("GSE174367_snRNA-seq_cell_meta.csv.gz", row.names = 1)
rna_seurat <- CreateSeuratObject(counts = rna_data, meta.data = rna_meta)
saveRDS(rna_seurat, file = "ad_rna_seurat.rds")

# scATAC-seq
library(Signac)
atac_meta <- read.csv(file = "GSE174367_snATAC-seq_cell_meta.csv.gz", header = TRUE)
rownames(atac_meta) <- atac_meta$Barcode
counts <- Read10X_h5("GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5")
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  min.cells = 10,
  min.features = 200
)

atac_seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = atac_meta
)
saveRDS(atac_seurat, file ="ad_atac_seurat.rds")