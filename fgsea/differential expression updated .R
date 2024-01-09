library(limma)
library(edgeR)
source('~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/scripts/fgsea/functions.R')

# Read and preprocess the main dataset
path_normMatrix <- "~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis"
setwd(path_normMatrix)
data <- read.csv("LW15_normalized.csv", header = TRUE)
rownames(data) <- data[, 1]
data <- data[, -c(1,2)]  # Skipping the first two columns (gene_id and gene_name)

# Set the path for results
path_results <- "~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/Differential Expression/New Comparisons"
# Perform differential expression analysis

# Define your groups outside the function
comparisonGroups <- list(
  Baseline = "H1",
  Experimental = "RVN"
)

performDEA(data, comparisonGroups, path_results, "RVN_vs_H1")

# You can add more comparisons by defining new groups in comparisonGroups
# and calling performDEA with those groups.


