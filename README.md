# NIPGR-NDRI-Project
RNA Seq data analysis of the Buffalo Genome NDRI- NIPGR collaboration project

Project Overview:
This project focuses on the RNA-Seq data analysis of embryonic stages of the buffalo genome in cloned and IVF individuals. The main objective is to identify differentially expressed genes (DEGs) which can help in understanding genetic expressions during different stages of development.

Data
The data used in this analysis includes RNA-Seq count data from various stages of embryonic development. This data is read from CSV files containing gene count matrix and phenotypic information of the samples.

Scripts and Analysis
The analysis is performed using the R script deseq2_nipgr.R of the total gene_count_matrix.csv, which utilizes the DESeq2 package to perform differential expression analysis. The script processes the count data, performs normalization, and identifies significant DEGs based on adjusted p-values and log fold changes.

Key steps in the script include:

Reading and preprocessing data.
Differential expression analysis using DESeq2.
Visualization of results using PCA, MA plots, volcano plots, and heatmaps.
Dependencies
R
DESeq2
openxlsx
EnhancedVolcano
RColorBrewer
pheatmap

Results:
Results of the analysis are saved in the form of an Excel file (ressig.xlsx), which contains the list of significantly differentially expressed genes. Various plots are also generated to visualize the results, aiding in the interpretation of the data.

Usage:
To run the script, ensure all dependencies are installed in R, and execute the script within the R environment. Adjust paths to data files as necessary to match your local setup.

Contributing:
For any modifications or improvements, please clone the repository, make your changes, and submit a pull request. Ensure that all changes are well-documented and tested.
