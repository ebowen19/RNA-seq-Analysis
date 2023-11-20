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

#Merging with New RNA-seq data****--------------------------------------------------------------------------------
new_rawCounts <- read.csv("ensembl_7860_raw.csv", header = TRUE, sep = ",")
raw_countMatrix <- read.csv("merged_rawCountMatrix2.csv", header = TRUE, sep = ",")

#take out columbns with low sums
simplified_newRawCounts <- new_rawCounts
# Calculate the sum of columns 2 to 4 for each row
row_sums <- rowSums(simplified_newRawCounts[, 2:4])

# Identify rows where the sum is less than 5
rows_to_delete <- which(row_sums < 5)

# # Remove the identified rows from the data frame
# simplified_newRawCounts <- simplified_newRawCounts[-rows_to_delete, ]
# nrow(simplified_newRawCounts)

#remove null lines from Shiru's raw count matrix
# Using subset to remove rows with sum equal to 0
shiru <- subset(raw_countMatrix, rowSums(raw_countMatrix) != 0)

library(dplyr)
shiru <- as.data.frame(raw_countMatrix)
shiru <- shiru %>%
    tibble::rownames_to_column(var = "RowNames")
 
colnames(shiru) [1] <- "Gene"

# make gene names consistent with the new data
shiru$Gene <- sub("\\..*", "", shiru$Gene)
colnames(simplified_newRawCounts)[1] <- "Gene"
shiru <- shiru[, c(-1)]
#merge 2 dataframes
rawCounts_merged <- merge(shiru, simplified_newRawCounts, by = "Gene", all=TRUE)
rawCounts_merged <- rawCounts_merged[ ] #remove non-gene rows
rawCounts_merged[is.na(rawCounts_merged)] <- 0 #replace missing vals w/ 0
new_colnames <- c("Atorvastatin", "Fluvastatin", "RFLV", "X126", "RVN1")
colnames(rawCounts_merged)[(ncol(rawCounts_merged) - 4):ncol(rawCounts_merged)] <- new_colnames

#save merged raw count matrix as a csv file
write.csv(rawCounts_merged, file = "human&renca_rawCountMatrix.csv", row.names = FALSE)
