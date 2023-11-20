setwd("/Users/elizabeth 1/Desktop/Wu Lab/RNA Seq practice/data")

#note: Shiru = renca, LW = 7860

#Merging with New RNA-seq data****--------------------------------------------------------------------------------
rawCounts_Shiru <- read.csv("renca_normalizedSeparately.csv", header = TRUE, sep = ",")
rawCounts_LW <- read.csv("ensembl_7860_raw.csv", header = TRUE, sep = ",")

# #remove unnecessary columns
# rawCounts_LW <- rawCounts_LW[, -c(1,2,3,4,5,7,8,9,10,11)]

#normalize each dataset separately
#SHIRU*******
shiru_data <- rawCounts_Shiru[,(2:10)]
rownames(shiru_data) <- rawCounts_Shiru[,1]
shiru_data <- round(shiru_data)

#Create DESeqDataSet object
colData <- data.frame(condition=1, sample_id=colnames(shiru_data))
dds_shiru <- DESeqDataSetFromMatrix(countData = shiru_data, colData = colData, design = ~1)

#Run DESeq2
dds_shiru <- DESeq(dds_shiru)

shiru_normalized <- counts(dds_shiru, normalized=T)
shiru_normalized <- round(shiru_normalized)

# #fix the row names
# rownames(shiru_normalized) <- sub("\\..*$", "", rownames(shiru_normalized))

#LW*****
LW_data <- rawCounts_LW[,(2:4)]
rownames(LW_data) <- rawCounts_LW[,1]
LW_data <- round(LW_data)

#Create DESeqDataSet object
colData <- data.frame(condition=1, sample_id=colnames(LW_data))
dds_LW <- DESeqDataSetFromMatrix(countData = LW_data, colData = colData, design = ~1)

#Run DESeq2
dds_LW <- DESeq(dds_LW)

LW_normalized <- counts(dds_LW, normalized=T)
LW_normalized <- round(LW_normalized)

library(dplyr)
#make gene names into a column
LW_normalized <- as.data.frame(LW_normalized)
LW_normalized <- LW_normalized %>%
  tibble::rownames_to_column(var = "RowNames")

colnames(LW_normalized) [1] <- "Gene"

shiru_normalized <- as.data.frame(shiru_normalized)
shiru_normalized <- shiru_normalized %>%
  tibble::rownames_to_column(var = "RowNames")

colnames(shiru_normalized)[1] <- "Gene"


# # make gene names consistent with the new data
# shiru$Gene <- sub("\\..*", "", shiru$Gene)
# colnames(rawCounts_Shiru)[1] <- "Gene"
# shiru <- shiru[, c(-1)]
#merge 2 dataframes
merged <- merge(LW_normalized, shiru_normalized, by = "Gene", all=TRUE)
# merged <- merged[-c(1,2,3,4,5) ,] #remove non-gene rows
merged[is.na(merged)] <- 0 #replace missing vals w/ 0

# #remove null lines from Shiru's raw count matrix
# # Using subset to remove rows with sum equal to 0
# shiru <- subset(rawCounts_LW, rowSums(rawCounts_LW) != 0)

# Calculate the sum of columns 2 to 4 for each row -- remove low count lines
row_sums <- rowSums(merged[, 2:13])

# Identify rows where the sum is less than 5
rows_to_delete <- which(row_sums < 5)

# Remove the identified rows from the data frame
merged <- merged[-rows_to_delete, ]
nrow(merged)

#save merged raw count matrix as a csv file
write.csv(merged, file = "renca&7860_normalizedSeparately.csv", row.names = FALSE)
