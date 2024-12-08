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

## S-LDSC Analysis

### Software Used
We utilized a Python script for LD Score Regression, available at the following GitHub repository: [ldsc](https://github.com/bulik/ldsc). Please follow the installation instructions provided by the provider to set up the environment.

### Running the Analysis
The S-LDSC analysis was conducted using GWAS summary statistics to analyze SNP heritability in Alzheimer's disease, focusing on individuals of European ancestry. The specific script used for this analysis can be found in the GitHub repository: [sldsc_bashcode.sh](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_bashcode.sh).

### Results and Visualization
Significant heritability enrichments were observed in regulatory genomic elements, demonstrating the critical role of epigenetic regulation in Alzheimer's disease. These results are illustrated in the figure below, which supports potential epigenetic targets for therapeutic strategies:

![Significant Heritability Enrichments](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_sig_plot.png)

The script to generate this plot is available in the GitHub repository: [sldsc_res.R](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_res.R).



