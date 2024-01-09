# --------------------------------------------------------------
# Script Name: common gene analysis.R
# Description: This script focuses on analyzing common genes,
#              particularly in tallying and generating
#              comprehensive reports.
#              Run after 'Extract Common Genes.R'.
# Author: Elizabeth Bowen
# Date: January 8, 2024
# --------------------------------------------------------------


source('~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/scripts/fgsea/functions.R')
library(dplyr)
library(tidyverse)
library(tidyr)

#function to remove 1st column & remove empty rows
simplify <- function(sublist) {
  if (is.data.frame(sublist)) {
    # Removes the first column if it is a dataframe
    # column_removed <- sublist[-1]
    # Remove rows where the 'Genes' column is an empty string
    df_clean <- sublist[sublist$Genes != "", ]
    rownames(df_clean) <- NULL
    return(df_clean)
  } else if (is.list(sublist)) {
    # Apply the function to each element of the list
    return(lapply(sublist, simplify))
  } else {
    # Returns the item unchanged if it is neither a list nor a dataframe
    return(sublist)
  }
}

# Function to write dataframes to subfolders
write_dataframes_to_subfolders <- function(data_lists, base_dir) {
  for (category in names(data_lists)) {
    # Create subfolder path
    subfolder_path <- file.path(base_dir, category)
    
    # Create subfolder if it doesn't exist
    if (!dir.exists(subfolder_path)) {
      dir.create(subfolder_path)
    }
    
    # Access the list associated with the category
    category_list <- data_lists[[category]]
    
    for (dataframe_name in names(category_list)) {
      # Create the path to the new file
      file_path <- file.path(subfolder_path, paste0(dataframe_name, ".csv"))
      
      # Write the dataframe to a CSV file
      write.csv(category_list[[dataframe_name]], file_path, row.names = FALSE)
    }
  }
}

main_folder <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/gsea/Common Genes/New Comparisons/Without trop2"
setwd(main_folder)
#read in data
data <- group_files_by_subfolder(main_folder, pattern = 'genes.csv')
data <- data$`Without trop2`
#get rid of 1st column/empty rows
data_cleaned <- simplify(data)
# Call the function to write the dataframes to subfolders
write_dataframes_to_subfolders(data_cleaned, main_folder)

# Assuming `data_cleaned` is your list with data frames for 'Go', 'Kegg', 'Hallmark', 'Reactome'
# Here we're assuming that 'data_cleaned' is a list of data frames, each corresponding to a category
# tally up the common genes for each dataframe 
geneTables <- compile_results(data_cleaned)

final_results <- process_and_combine_gene_data(geneTables)

# save total gene count data as CSVs & HTML files
write.csv(final_results$Downregulated, "Total Gene Counts Down.csv", row.names = FALSE)
write.csv(final_results$Upregulated, "Total Gene Counts Up.csv", row.names = FALSE)


