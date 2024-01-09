#this script performs GSEA and extracts the common pathways for any gene set 
#out of Hallmark, Reactome, GO, or KEGG.

# load required libraries
library(fgsea)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
library(dplyr)
library(biomaRt)
library(httr)
library(jsonlite)
source('~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/scripts/fgsea/functions.R')

#load in gene set data
geneSets <- loadAndOrganizeGeneSets()

#Execution

# load in data
path_DE.trop2 <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/differential expression"

#define comparisons, specifying species for each
#input file path should lead to Differential expression analysis file for that comparison
comparisonsList <- list(
  # "H2_vs_H1" = list(dataPath = paste0(path_DE,"/original comparisons/H2_vs_H1.csv"), species = "Mus musculus")#,
  # "H2_vs_H12" = list(dataPath = paste0(path_DE,"/original comparisons/H2_vs_H12.csv"), species = "Mus musculus"),
  # "RVN_vs_RC" = list(dataPath = paste0(path_DE,"/original comparisons/RVN_vs_RC.csv"), species = "Mus musculus"),
  # "H2_vs_RFL" = list(dataPath = paste0(path_DE,"/New Comparisons/H2_vs_RFL.csv"), species = "Mus musculus"),
  # "RVN_vs_H12" = list(dataPath = paste0(path_DE,"/New Comparisons/RVN_vs_H12.csv"), species = "Mus musculus"),
  "RVN_vs_H1" = list(dataPath = paste0(path_DE,"/New Comparisons/RVN_vs_H1.csv"), species = "Mus musculus")#,
  # "trop2_vs_non" = list(dataPath = paste0(path_DE.trop2,"/fullResults_trop2vsnon.csv"), species = "Homo sapiens")
)

path_DE <- "~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/Differential Expression"
path_results <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/gsea/New Comparisons/GO"

# Run the GSEA analysis
gseaResultFiles <- runAnalysis(comparisonsList, "GO", geneSets, path_results)


