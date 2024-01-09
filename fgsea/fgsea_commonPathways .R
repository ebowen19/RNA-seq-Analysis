library(dplyr)
source('~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/scripts/fgsea/functions.R')

processGseaResults <- function(subdirectory) {
  # Construct the base path
  basePath <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/gsea/"
  split_string <- strsplit(subdirectory, "/")
  og_subdirectory <- paste0(split_string[[1]][1], "/", split_string[[1]][3])
  pathResults <- paste0(basePath, og_subdirectory)
  gene_set <- sub(".*/([^/]+)$", "\\1", og_subdirectory)

  #create output path for the results
  outputPath <- paste0(basePath, "New Comparisons/Without trop2/", split_string[[1]][3])

  # Paths to fgsea result files
  gseaResultFilePaths <- list(
    # when we fixed substr error, put it in the og comparisons folder. using regular expression to fetch
    "H2_vs_H1" = paste0("/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/gsea/Original Comparisons/",gene_set,"/fgseaResults_H2_vs_H1.csv"),
    "H2_vs_H12" = paste0(pathResults,"/fgseaResults_H2_vs_H12.csv"),
    "RVN_vs_RC" = paste0(pathResults,"/fgseaResults_RVN_vs_RC.csv"),
    #"trop2_vs_non" = paste0(pathResults,"/fgseaResults_trop2_vs_non.csv")
    "H2_vs_RC" = paste0(pathResults,"/fgseaResults_H2_vs_RFL.csv"),
    "RVN_vs_H1" = paste0(pathResults,"/fgseaResults_RVN_vs_H1.csv"),
    "RVN_vs_H12" = paste0(pathResults,"/fgseaResults_RVN_vs_H12.csv")
  )
  
  # Read the results into data frames and create a list of these data frames
  gseaResults <- lapply(gseaResultFilePaths, read.csv, check.names = FALSE, row.names = 1)
  
  # Replace "RVN_vs_RC" with reference comparison, if it needs to be dynamic, add it as a function argument
  alignedPathways <- extractCommonPathways(gseaResults, "RVN_vs_RC", outputPath)
  
  # Return the results
  return(list(
    up = data.frame(alignedPathways$pathways_list$up),
    down = data.frame(alignedPathways$pathways_list$down)
  ))
}

alignedPathways <- processGseaResults("New Comparisons/Without trop2/Reactome")
View(alignedPathways$up)
View(alignedPathways$down)









