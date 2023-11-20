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

#Correct for Batch Effect 
# Before this step, you might need to define the 'batch' variable properly based on your actual batch information.
# Assuming 'batch' is a factor that describes the batch for each sample:
batch <- as.factor(rep(c("Shiru", "Moe"), times = c(4, 5)))  

# Correcting for batch effects using comBat
modcombat <- model.matrix(~1, data=colData)  # or other appropriate model considering your experiment design
# Batch correcting the normalized data
batch_corrected_normalized <- ComBat(dat=vst_data, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

write.table(normalizedCounts,file="______s.tsv", col.names = TRUE, sep = "\t")

#PCA analysis
df2 <- as.data.frame(batch_corrected_normalized)
pca <- prcomp(t(df2), scale=TRUE)

library(plotly) 
new <- as.data.frame(pca$x)

plot_ly(new, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", marker = list(size = 2)) %>%
  add_trace(type = "scatter3d", mode = "text", text = colnames(normalizedCounts), textposition = "top center")

