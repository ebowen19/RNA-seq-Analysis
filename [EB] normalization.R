wd <- "/Users/elizabeth 1/Desktop/Wu lab/RNA seq practice/data"
setwd(wd)
library("DESeq2")
library("sva") # for removing batch effect 

#Read in the data
data <- read.csv("_________.csv", header = TRUE) #run the batch correction/final normalization on the separately normalized counts 
View(data)
#fix the structure of the dataframe
dt <- data[,(2:10)]
rownames(dt) <- data[,1]
dt <- round(dt)
View(dt)

#Create DESeqDataSet object
colData <- data.frame(condition=1, sample_id=colnames(dt))
dds <- DESeqDataSetFromMatrix(countData = dt, colData = colData, design = ~1)

#Run DESeq2
dds <- DESeq(dds)

normalizedCounts <- counts(dds, normalized=T)
write.table(normalizedCounts,file="______s.tsv", col.names = TRUE, sep = "\t")

