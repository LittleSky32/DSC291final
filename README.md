# DSC291final
In the Github repo, please add a README file in the main page on a brief
introduction of your project (you could reuse anything from your report), how to
use your code, including references to any existing software you used and notes
on how you installed it.

## S-LDSC Analysis

We conducted Stratified LD Score Regression (S-LDSC) to analyze SNP heritability in Alzheimer's disease using GWAS summary statistics, specifically focusing on individuals of European ancestry for genetic consistency. The analysis incorporated pre-computed LD scores and functional annotations from the 1000 Genomes Project, with a particular emphasis on enhancers and promoters. The script used for this analysis is available in the GitHub repository: [`sldsc_bashcode.sh`](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_bashcode.sh).

Significant heritability enrichments were observed in regulatory genomic elements, demonstrating the critical role of epigenetic regulation in Alzheimer's disease. These results are illustrated in the figure below, which supports potential epigenetic targets for therapeutic strategies:

![Significant Heritability Enrichments](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_sig_plot.png)

The script to generate this plot is available in the GitHub repository: [`sldsc_res.R`](https://github.com/LittleSky32/DSC291final/blob/main/sldsc_res.R).

Future improvements could involve generating cohort-specific LD scores and integrating cell type-specific annotations to further refine our understanding of the genetic factors influencing Alzheimer's disease.

