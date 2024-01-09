### fgsea - Gene Set Enrichment Analysis

This subfolder within the RNA-seq Analysis repository contains scripts and notebooks dedicated to performing Gene Set Enrichment Analysis using the `fgsea` package in R. It provides a comprehensive workflow for analyzing gene sets and pathways.

**Workflow Overview:**
1. **Running Differential Expression Analysis** on selected gene sets.
2. **Generating Volcano Plots** to visualize differential expression data.
3. **Performing GSEA** with ranked gene lists.
4. **Extracting Common Pathways** for differential expression comparisons.
5. **Identifying Common Leading Edge Genes** in these pathways.
6. **Tallying Common Genes' Occurrences** across pathways and gene sets.
7. **Creating Linked HTML Files** for detailed analysis, with gene counts linking to pathway and gene data.

**Scripts and Notebooks:**
- `fgsea.R`: Main script for loading required libraries and initial data processing.
- `volcano plots.R`: Script for creating volcano plots.
- `differential expression updated.R`, `common gene analysis.R`, `fgsea_commonPathways.R`: Scripts for various stages of GSEA.
- `functions.R`: Contains functions used across different scripts for analysis.
- `HTML generation.R` and `HTML Editing.ipynb`: For generating and editing HTML reports.

**Usage:**
- Each script is part of the analysis pipeline and should be run in the specified order under *Workflow Overview* for optimal results.
