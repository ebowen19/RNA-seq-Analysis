#BiocManager::install("limma")
setwd("/Users/elizabeth 1/Desktop/Wu Lab/RNA Seq practice/data")
data <- read.csv("renca_normalizedTogether.tsv", header = TRUE, sep = "\t")
rownames(data) <- data[, 1]
# # Remove the first column
data <- data[, -c(1:3)]

# #create sub-matrices for each comparison
# h2_vs_h12 <- data[, c(3,2)]
# RVN_vs_RFL <- data[, c(6,7,4,5)]
# h2_vs_h1 <- data[, c(3,1)]
trop2 <- read.csv("trop2.csv", header = TRUE, row.names = 1, sep = ",")

# Designating groups for each sample
groups <- factor(c("h1", "h12", "h2", "RC", "RC", "RVN", "RVN"))
# Create a design matrix
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

#fit the model & specify contrasts
library(limma)

# Fit the model to your data
fit <- lmFit(data, design)

# Specify contrasts for the desired comparisons
cont.matrix <- makeContrasts(
  h2_vs_h12 = h2 - h12,
  RVN_vs_RC = RVN - RC,
  h2_vs_h1 = h2 - h1,
  levels = design
)

# Fit the contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# for trop2 --> Designating groups for each sample
groups_t <- factor(c(rep("non", 6), rep("trop2", 6)))
# Create a design matrix
design_t <- model.matrix(~ 0 + groups_t)
colnames(design_t) <- levels(groups_t)

# Fit the model to your data
fit_t <- lmFit(trop2, design_t)

# Specify contrasts for the desired comparisons
cont.matrix_t <- makeContrasts(
  trop2_vs_non = trop2 - non,
  levels = design_t
)

# Fit the contrasts
fit2_t <- contrasts.fit(fit_t, cont.matrix_t)
fit2_t <- eBayes(fit2_t)

# Results for h2 vs h12
results_h2vsh12 <- topTable(fit2, coef="h2_vs_h12", n=Inf)
results_h2vsh12 <- results_h2vsh12[results_h2vsh12$logFC != 0, ] #remove rows (genes) where log-fold change = 0
h2_up_vs_h12 <- results_h2vsh12[results_h2vsh12$logFC > 0, ] #top 100 upregulated genes in h2 vs h12
h2_up_vs_h12 <- head(h2_up_vs_h12[order(abs(h2_up_vs_h12$logFC), decreasing = TRUE), ], 100) #sort by fold change instead of by adjusted p-value & extract top 100 vals
pval_h2_up_vs_h12 <- head(results_h2vsh12[results_h2vsh12$logFC > 0, ], 100) 
h2_down_vs_h12 <- results_h2vsh12[results_h2vsh12$logFC < 0, ] #top 100 downregulated genes in h2 vs h12
h2_down_vs_h12 <- head(h2_down_vs_h12[order(abs(h2_down_vs_h12$logFC), decreasing = TRUE), ], 100)
pval_h2_down_vs_h12 <- head(results_h2vsh12[results_h2vsh12$logFC < 0, ], 100) 

# Results for RVN (RVN1 & RVN2) vs RC (RC1 & RC2)
results_RVNvsRC <- topTable(fit2, coef="RVN_vs_RC", n=Inf)
results_RVNvsRC <- results_RVNvsRC[results_RVNvsRC$logFC != 0, ]
RVN_up_vs_RC <- head(results_RVNvsRC[results_RVNvsRC$logFC > 0, ], 100) 
RVN_down_vs_RC <- head(results_RVNvsRC[results_RVNvsRC$logFC < 0, ], 100) 

# Results for h2 vs h1
results_h2vsh1 <- topTable(fit2, coef="h2_vs_h1", n=Inf)
results_h2vsh1 <- results_h2vsh1[results_h2vsh1$logFC != 0, ]
h2_up_vs_h1 <- results_h2vsh1[results_h2vsh1$logFC > 0, ]
h2_up_vs_h1 <- head(h2_up_vs_h1[order(abs(h2_up_vs_h1$logFC), decreasing = TRUE), ], 100)
pval_h2_up_vs_h1 <- head(results_RVNvsRC[results_RVNvsRC$logFC > 0, ], 100) 
h2_down_vs_h1 <- results_h2vsh1[results_h2vsh1$logFC < 0, ] #sorted by p-value
h2_down_vs_h1 <- head(h2_down_vs_h1[order(abs(h2_down_vs_h1$logFC), decreasing = TRUE), ], 100)
pval_h2_down_vs_h1 <- head(results_RVNvsRC[results_RVNvsRC$logFC < 0, ], 100) 

# Results for trop2 vs non
results_trop2vsnon <- topTable(fit2_t, coef = "trop2_vs_non", n=Inf)
results_trop2vsnon <- results_trop2vsnon[results_trop2vsnon$logFC != 0, ]
trop2_up_vs_non <- head(results_trop2vsnon[results_trop2vsnon$logFC > 0, ], 100) 
trop2_down_vs_non <- head(results_trop2vsnon[results_trop2vsnon$logFC < 0, ], 100) 

results <- c("trop2_down_vs_non", "trop2_up_vs_non","results_trop2vsnon","h2_down_vs_h1", 
             "pval_h2_down_vs_h1", "pval_h2_up_vs_h1","h2_up_vs_h1","results_h2vsh1",
             "RVN_down_vs_RC","RVN_up_vs_RC","results_RVNvsRC","results_h2vsh12",
             "pval_h2_down_vs_h12", "h2_down_vs_h12","pval_h2_up_vs_h12","h2_up_vs_h12")
setwd("/Users/elizabeth 1/Desktop/differential expression")

# for (res in results) {
#   # Using get() to retrieve the dataframe using its character name
#   df <- get(res)
#   
#   # Constructing a filename for the CSV (using paste0 to remove spaces)
#   filename <- paste0(res, ".csv")
#   
#   # Saving the dataframe as a CSV
#   write.csv(df, filename, row.names = TRUE)
# }


