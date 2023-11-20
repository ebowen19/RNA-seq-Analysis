#PCA analysis
df2 <- as.data.frame(batch_corrected_normalized)
pca <- prcomp(t(df2), scale=TRUE)

library(plotly) 
new <- as.data.frame(pca$x)

plot_ly(new, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers", marker = list(size = 2)) %>%
  add_trace(type = "scatter3d", mode = "text", text = colnames(normalizedCounts), textposition = "top center")

