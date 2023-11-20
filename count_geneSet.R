#load in data
setwd("/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/gsea/pathway documentation")
reactome_mouse <- read.csv("reactome_mouse.gmt", header = FALSE, sep = "\t")
reactome_mouse <- reactome_mouse[-c(1,2)] # remove non-gene columns
reactome_human <- read.csv("reactome_human.gmt", header = FALSE, sep = "\t")
reactome_human <- reactome_human[-c(1,2)]
hallmark_mouse <- read.csv("hallmark_mouse.gmt", header = FALSE, sep = "\t")
hallmark_mouse <- hallmark_mouse[-c(1,2)]
hallmark_human <- read.csv("hallmark_human.gmt", header = FALSE, sep = "\t")
hallmark_human <- hallmark_human[-c(1,2)]
go_human <- read.csv("go_human.gmt", header = FALSE, sep = "\t")
go_human <- go_human[-c(1,2)]
go_mouse <- read.csv("go_mouse.gmt", header = FALSE, sep = "\t")
go_mouse <- go_mouse[-c(1,2)]

# Load dplyr for easy data manipulation
library(dplyr)

# Define a function to count unique genes
count_unique_genes <- function(gene_set) {
  all_genes <- unlist(gene_set, use.names = FALSE)
  all_genes_clean <- na.omit(all_genes[all_genes != ""])
  unique_genes <- unique(all_genes_clean)
  return(length(unique_genes))
}

# Initialize a data frame to hold the results
results <- data.frame(GeneSet = character(), NumberOfGenes = integer(), stringsAsFactors = FALSE)

# List of gene sets and their names
gene_sets <- list(reactome_mouse, reactome_human, hallmark_mouse, hallmark_human, go_mouse, go_human) # list of data frames
names(gene_sets) <- c("Reactome Mouse", "Reactome Human", "Hallmark Mouse", "Hallmark Human", "Go Mouse", "Go Human")

# Loop over the list and apply the function to each gene set
for (gene_set_name in names(gene_sets)) { #gene_set_,name is arbitrary var name
  num_unique_genes <- count_unique_genes(gene_sets[[gene_set_name]])
  results <- rbind(results, data.frame(GeneSet = gene_set_name, NumberOfGenes = num_unique_genes))
}

# Print the results
print(results)

write.csv(results, "geneSetCounts.csv")

#Count the # of gene sets encompassed --> first re-include the 1st 2 columns & redefine gene_sets
for (gene_set_name in names(gene_sets)) { #gene_set_,name is arbitrary var name
  print(paste(gene_set_name, sum(grepl("https", gene_sets[[gene_set_name]]$V2))))
}


# # Convert all columns to a single vector
# all_genes <- unlist(reactome_mouse, use.names = FALSE)
# 
# # Remove NA or empty values
# all_genes_clean <- na.omit(all_genes[all_genes != ""]) # na.omit only removes NA values, not empty strings by default
# #using na.omit on a subset of the all_genes vector without any empty values to assign to cleaned vector
# 
# # Get the unique gene entries
# unique_genes <- unique(all_genes_clean)
# 
# # Count the number of unique genes
# num_unique_genes <- length(unique_genes)
# 
# # Print the number of unique genes
# print(num_unique_genes)
