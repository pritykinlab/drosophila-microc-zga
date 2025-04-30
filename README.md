# Drosophila Micro-C ZGA Analysis

## Overview

This repository contains analysis code and figures for studying chromatin organization during the Drosophila zygotic genome activation (ZGA) using Micro-C.

## Figure Locations

The figures and analysis are spread across several Jupyter notebooks:

### Main Figures

Most of the main figures are in one notebook, but some are split into other notebooks for clarity.

- `main_analysis/cleaned_up_figures.ipynb.gz`: Contains the majority of the main figure plots and analysis
- `main_analysis/Make Figure 5.ipynb`: Contains the plots for Figure 5
- `main_analysis/Make Figures 2C-E.ipynb.gz`: Contains the plots for Figures 2C-E
- `main_analysis/Make Figures 3B.ipynb`: Contains the plots for Figure 3B
- `main_analysis/Make Figures 4C, E.ipynb`: Contains the plots for Figures 4C and 4E
- `main_analysis/Make Figures 5E.ipynb`: Contains the plots for Figure 5E
- `main_analysis/Make Figures 6E.ipynb`: Contains the plots for Figure 6E

### D. virilis Analysis  
- `cleaned_dvir_browser_view.ipynb.gz`: Contains comparative analysis and plots for D. virilis

### Supplementary plots
- `cleaned_up_figures.ipynb.gz`: Contains most of the supplementary plots
- `supplementary_analysis/Make Figures S2, S3.ipynb`: Contains the plots for Figures S2 and S3
- `supplementary_analysis/Make Figures S8.ipynb`: Contains the plots for Figure S8
- `supplementary_analysis/2_Make_Insulation_Correlation_Figures.ipynb.gz`: Contains reproducibility plots for the insulation scores

### Loop Analysis
The loop-associated analysis is split across three notebooks:
- `supplementary_analysis/0_Process_loops.ipynb`: Initial processing and analysis of chromatin loops
- `supplementary_analysis/1_Make_Loop_Supp_Figures.ipynb`: Supplementary figures related to loop analysis
- `supplementary_analysis/2_Make_Insulation_Correlation_Figures.ipynb.gz`: Analysis of insulation and correlation with loops
- `supplementary_analysis/loop_DESEQ2/`: Contains the DESEQ2 analysis for loop counts
- `supplementary_analysis/loop_DESEQ2/COUNT/Make_DESeq2_Counts.ipynb`: Contains the DESEQ2 analysis for loop counts

## Data Access

The raw data and processed files associated with this analysis are available through GEO (accession number pending).
ChIP-seq data aggregated over loops and boundaries is available in `CSE_adata.h5ad`, including metadata for each ChIP-seq sample, as well as annotations of loop anchors and clusters for each CSE.

## Contact

For questions about the code or analysis, please open an issue in this repository or email gdolsten@princeton.edu.