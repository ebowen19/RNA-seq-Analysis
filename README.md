### RNA-seq Analysis

This repository is dedicated to RNA sequencing (RNA-seq) analysis, incorporating a variety of R scripts for different stages of analysis, including differential expression analysis, gene set enrichment, and pathway analysis. The repository utilizes R packages like `limma`, `edgeR`, and `fgsea` to facilitate comprehensive RNA-seq data analysis.

**Main Features:**
1. **Differential Expression Analysis:** Scripts for identifying differentially expressed genes using `limma` and `edgeR`.
2. **Gene Set Enrichment Analysis (GSEA):** Using the `fgsea` package to analyze gene sets and pathways.
3. **Volcano Plot Generation:** Visualizing differential expression results in a compelling format.
4. **Pathway Analysis:** Extracting common pathways and leading edge genes from GSEA results.
5. **HTML Reporting:** Automated generation of detailed HTML reports for analysis results.
6. **Functional Gene Testing:** Scripts for testing and validating functional genes.
7. **Top Gene Analysis:** Identifying and analyzing top genes from datasets.

**Scripts Overview:**
- `RNAseqAnalysis.R`: Main script for overall RNA-seq analysis.
- `topGenes.R`: Focused on analyzing top genes from the dataset.
- `testFunctionalGenes.R` and `testDifferentialGenes.R`: For testing and validating genes based on their functions and differential expression.
- `Package Installation Instructions.R`: Guide for setting up necessary R packages.

**How to Use:**
- Scripts are organized to follow the typical workflow of RNA-seq analysis.
- Users are encouraged to follow the scripts in sequential order, starting from differential expression to pathway analysis for the best results.

**Note:**
- The `fgsea` subfolder contains additional scripts and notebooks specifically for GSEA using the `fgsea` package.
