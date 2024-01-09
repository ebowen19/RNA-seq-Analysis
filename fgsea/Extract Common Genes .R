# --------------------------------------------------------------
# Script Name: Extract Common Genes.R
# Description: This script is primarily focused on extracting
#              common leading edge genes from fgsea results.
#              Run before 'common gene analysis.R'.
# Author: Elizabeth Bowen
# Date: January 8, 2024 
# --------------------------------------------------------------

library(httr)
library(jsonlite)
library(dplyr)
library(tidyverse)
source('~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/scripts/fgsea/functions.R')

main_folder <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/gsea"
common_pathways_folder <- paste0(main_folder,"/New Comparisons/Without trop2")
setwd(main_folder)

#read in data
common_pathways_dfs <- group_files_by_subfolder(common_pathways_folder, 'Common')
common_pathways_dfs <- list('New Comparisons' = common_pathways_dfs$'Without trop2')
fgseaResults_list <- group_files_by_subfolder(main_folder, 'fgsea')
fgseaResults_list[[3]] <- NULL #remove the sensi vs non data
fgseaResults_list[[2]] <- NULL

#rename the fgsea dataframes so that they're each just named the name of their comparison
for (category in names(fgseaResults_list)) {
  # The category itself contains a list of dataframes
  for (subcategory in names(fgseaResults_list[[category]])) {
    # Now subcategory is the list that contains the dataframes
    for (comparison in names(fgseaResults_list[[category]][[subcategory]])) {
      # Now we are at the level of the dataframes, where we can rename them
      # Extract the comparison name from the dataframe name
      new_name <- sub("fgseaResults_", "", comparison)
      # Rename the dataframe within the list
      names(fgseaResults_list[[category]][[subcategory]])[names(fgseaResults_list[[category]][[subcategory]]) == comparison] <- new_name
    }
  }
}

# Run the main function
common_leading_edges_results <- extract_common_leading_edges(fgseaResults_list, common_pathways_dfs, exclude_comparison=trop2_vs_non) 

results_folder <- paste0(main_folder,"/Common Genes/New Comparisons/Without trop2")
setwd(results_folder)
# save results
final_results <- save_final_results(common_leading_edges_results, results_folder)

