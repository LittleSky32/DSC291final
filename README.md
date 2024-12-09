# DSC291final
In the Github repo, please add a README file in the main page on a brief
introduction of your project (you could reuse anything from your report), how to
use your code, including references to any existing software you used and notes
on how you installed it.

## Introduction
This project integrates genetic, transcriptomic, and epigenomic data to determine how AD
associated variants affect gene expression in microglia. We hypothesize that epigenetic 
regulation, revealed by differentially accessible regions (DARs) in scATAC-seq data, 
contributes more significantly to microglial transcriptional changes in AD than direct genetic 
effects of single nucleotide polymorphisms (SNPs) located within gene bodies. By comparing 
these effects, we aim to elucidate AD-associated regulatory mechanisms and identify potential 
therapeutic targets. 

## scRNA-seq Differential Expression Analysis and scATAC-seq Differential Accessibility Analysis

### scRNA-seq and scATAC-seq processing

This script [`make_seurat.R`](https://github.com/LittleSky32/DSC291final/blob/main/create_seurat.R) preprocesses scRNA-seq and scATAC-seq data to generate Seurat objects for downstream analysis. The scRNA-seq data is loaded as a gene-by-cell matrix, and the scATAC-seq data is loaded as a peak-by-cell matrix. Metadata is integrated to include essential information such as cell type and diagnosis (AD or Control). The processed Seurat objects are saved as `ad_rna_seurat.rds` and `ad_atac_seurat.rds`, ready for further analysis.

### ATAC-seq Differential Accessibility Analysis

This script [`atac_seq_DEG.R`](https://github.com/LittleSky32/DSC291final/blob/main/atac_seq_DEG.R) identifies differentially accessible regions (DARs) from scATAC-seq data. The ATAC-seq Seurat object is first preprocessed using Signac's standard normalization steps (https://stuartlab.org/signac/articles/pbmc_vignette), and differential accessibility is computed for AD versus Control in microglia. Annotated DARs are linked to genes using hg38 annotations, focusing on promoter-specific regions. Significant DARs are filtered based on an adjusted p-value < 0.05. The script further integrates DARs with GWAS summary statistics to map significant SNPs to regulatory regions and genes, which is only for internal checks of interest.

### RNA-seq Differential Expression Analysis

This script [`process_rna_deg.R`](https://github.com/LittleSky32/DSC291final/blob/main/rna_seq_DEG.R) identifies differentially expressed genes (DEGs) from scRNA-seq data. The Seurat object is preprocessed using standard steps form the Seurat vignette (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html), including quality control (removing cells with fewer than 200 or more than 2,500 features), normalization, feature selection, and principal component analysis (PCA). Differential expression analysis is then performed to identify DEGs between AD and Control microglia. Significant DEGs are filtered based on an adjusted p-value < 0.05 and an absolute log2 fold change > 0.25. The resulting list of DEGs is saved for downstream analyses.

### Session Information
Details about the session environment and package versions used for this analysis can be found at the end of this `README.md` file.

## GWAS Gene Body Mapping

This script [`GWAS_gene_body.R`](https://github.com/LittleSky32/DSC291final/blob/main/GWAS_gene_body.R) maps significant SNPs identified from GWAS summary statistics to genes located within gene bodies. Using positions of HapMap3 variants and a gene annotation file (obtained from HW2, original paper for HapMap3 could be located at (DOI: 10.1038/nature09298)), the script integrates GWAS-significant SNPs (p-value < 5e-8) with gene body annotations to identify genes potentially influenced by these SNPs. 

Both the gene annotation file [`gene_annot.txt.gz`](https://github.com/LittleSky32/DSC291final/blob/main/gene_annot.txt.gz) and the zip-compressed HapMap3 positions file [`GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim.zip`](https://github.com/LittleSky32/DSC291final/blob/main/GTEx_v8_genotype_EUR_HM3_exclude_dups.allchr.reorder.bim.zip) can be found in the GitHub repository accompanying this analysis.

## DEG Overlapping Analysis

This script [`deg_overlap_analysis.R`](https://github.com/LittleSky32/DSC291final/blob/main/overlap_DEG.R) performs an overlap analysis to compare differentially expressed genes (DEGs) identified from scRNA-seq data with genes mapped from scATAC-seq DARs, promoter regions. Additionally, it visualizes the overlaps using Venn diagrams and bar charts.

## S-LDSC Analysis

### Software Used
We utilized a Python script for LD Score Regression, available at the following GitHub repository: [ldsc](https://github.com/bulik/ldsc). Please follow the installation instructions provided by the provider to set up the environment.

### Running the Analysis
The S-LDSC analysis was conducted using GWAS summary statistics to analyze SNP heritability in Alzheimer's disease, focusing on individuals of European ancestry. The specific script used for this analysis can be found in the GitHub repository: [sldsc_bashcode.sh](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_bashcode.sh).

### Results and Visualization
Significant heritability enrichments were observed in regulatory genomic elements, demonstrating the critical role of epigenetic regulation in Alzheimer's disease. These results are illustrated in the figure below, which supports potential epigenetic targets for therapeutic strategies:

![Significant Heritability Enrichments](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_sig_plot.png)

The script to generate this plot is available in the GitHub repository: [sldsc_res.R](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_res.R).

## Important Session Information for R

- **R version:** 4.3.1
- **Packages used:**
  - `Seurat` (v5.1.0)
  - `Signac` (v1.14.0)
  - `dplyr` (v1.1.4)
  - `ggplot2` (v3.5.1)
  - `patchwork` (v1.3.0)
  - `devtools` (v2.4.5)
  - `annotatr` (v1.28.0)
  - `TxDb.Hsapiens.UCSC.hg38.knownGene` (v3.18.0)
  - `org.Hs.eg.db` (v3.18.0)
  - `GenomicRanges` (v1.54.1)
  - `data.table` (v1.16.2)
  - `ggvenn` (v0.1.10)
  - `ggplot2` (v3.5.1)
  - `data.table` (v1.16.2)

All required packages can be installed through CRAN or Bioconductor to reproduce the analysis.
