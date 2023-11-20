wd <- "/Users/elizabeth 1/Desktop/Wu lab/RNA seq practice/data"
setwd(wd)
library("DESeq2")
library("sva") # for removing batch effect 

#Read in the data
data <- read.csv("renca_normalizedTogether.tsv", header = TRUE, sep = "\t") #run the batch correction/final normalization on the separately normalized counts 
View(data)
#fix the structure of the dataframe
dt <- data[,(2:10)]
rownames(dt) <- data[,1]
dt <- round(dt)
View(dt)

#Correct for Batch Effect 
# Before this step, you might need to define the 'batch' variable properly based on your actual batch information.
# Assuming 'batch' is a factor that describes the batch for each sample:
batch <- as.factor(rep(c("Moe", "Shiru"), times = c(5, 4)))  

# Correcting for batch effects using comBat
#modcombat <- model.matrix(~1, data=colData)  # or other appropriate model considering your experiment design

# If you have other conditions, e.g., 'condition'
colData$condition <- factor(c('A','B','C','D','E','F','F','G','G'))  # fill with actual conditions
# Then your model matrix should consider this
modcombat <- model.matrix(~ condition, data=colData) 

# Batch correcting the normalized data
batch_corrected_normalized <- ComBat(dat=dt, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
write.csv(batch_corrected_normalized,file="renca_batchCorrected_normalized.csv")
