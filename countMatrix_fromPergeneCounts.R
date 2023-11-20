# Load the required package
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# REPLACE WITH SAMPLE**!!!!
pergene_files <- list.files(path = "../pergene_counts/RVN2", pattern = "\\.pergene_counts$", full.names = TRUE)

# Initialize an empty list to store data tables
pergene_data_list <- list()

# Read and store each file in the list
for (file in pergene_files) {
  # Read the file with 'fread' and specify column names
  data <- fread(file, header = FALSE, col.names = c("Gene", basename(file)))
  
  # Store the data frame in the list using the file name as the key
  pergene_data_list[[basename(file)]] <- data
}

# Merge the data tables into a single count matrix by gene name (column "Gene")
merged_counts <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), pergene_data_list)

# Fill missing values with zeros
merged_counts[is.na(merged_counts)] <- 0

columns <- ncol(merged_counts)
#merge all the columns together
mergedCounts_RVN2 <- data.frame(
  Gene = merged_counts$Gene,
  Count = rowSums(merged_counts[, 2:columns])
)

View(merged_counts)
View(mergedCounts_RVN2)

#put all the dataframes for each sample into one count matrix
raw_countMatrix <- data.frame(
  Gene = merged_counts$Gene,
  RC1 = mergedCounts_RC1[, 2],
  RC2 = mergedCounts_RC2[, 2],
  RVN1 = mergedCounts_RVN1[, 2],
  RVN2 = mergedCounts_RVN2[, 2]
)


View(raw_countMatrix)

# Save the count matrix as a CSV file
write.csv(raw_countMatrix, file = "raw_countMatrix.csv", row.names = FALSE)

