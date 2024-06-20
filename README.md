README


Overview
This project aims to analyze gene expression and signaling pathway differences, as well as intercellular communication patterns, between GBM cells from tumor and peritumoral regions using single-cell RNA sequencing (scRNA-seq) data. The R code provided comprehensive analyses to identify key biological insights into glioma microenvironments.


Contents
The repository contains the following files:

Cellchat/: Directory containing R scripts for data preprocessing, analysis, and visualization in intercellular analysis.
Basic/: Directory containing R scripts for basic processing of scRNA-seq data.
Pseudobulk/: Directory containing R scripts for differential expression analysis of GBM cells from tumor/peripheral regions.
OtherCells/: Directory containing R scripts for differential expression analysis of other cells such as astrocyte, immune cells and macrophage.

README.md: This file, providing an overview of the project.


Requirements
The following R packages are required to run the code:
Seurat
dplyr
ggplot2
ComplexHeatmap
CellChat
SingleCellExperiment
scater
edgeR
limma


Make sure to adjust file paths and parameters as needed within each script.


Contact
For any questions or issues, please contact Jianyu Zhang at zwczzhzjy@163.com.

This README provides a brief overview of the project and instructions for running the provided R code. It is recommended to review each script for detailed information about the analysis steps and to ensure all dependencies are properly installed before starting the analysis.
