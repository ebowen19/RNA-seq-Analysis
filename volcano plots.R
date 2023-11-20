# make volcano plot
library(ggplot2)

# Creating a custom volcano plot
data_frames <- list(
  trop2vsnon = results_trop2vsnon#,
  #RVNvsRC = results_RVNvsRC#,
  #h2vsh1 = results_h2vsh1#,
  #h2vsh12 = results_h2vsh12
)

for(name in names(data_frames)) {
  volcano_data <- data_frames[[name]]
  
  # Determine significance and direction (up or down regulated)
  volcano_data$significant <- ifelse(volcano_data$adj.P.Val < 0.05 & 
                                       abs(volcano_data$logFC) > 3.5, "yes", "no")
  volcano_data$direction <- ifelse(volcano_data$logFC > 0, "up", "down")
  
  # Calculate the IQR thresholds for logFC
  Q1 <- quantile(volcano_data$logFC, 0.15)
  Q3 <- quantile(volcano_data$logFC, 0.85)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  # Plot
  p <- ggplot(volcano_data, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = ifelse(significant == "yes" & direction == "up", "Significant Up", 
                                  ifelse(significant == "yes" & direction == "down", "Significant Down", "Not Significant"))), 
               alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(name, "volcano plot"), # using the name for the title
         x = "Log2 Fold Change",
         y = "-Log10 p-value",
         color = "Gene Regulation") +
    scale_color_manual(values = c("Not Significant" = "gray", 
                                  "Significant Up" = "red", 
                                  "Significant Down" = "blue")) +
    coord_cartesian(xlim = c(lower_bound, upper_bound), ylim = c(0, 6)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # You can print the plot directly, or save it using ggsave, etc.
  print(p)
 # ggsave(filename = paste0(name, "_p.06 volcano plot.jpg"), plot = last_plot())
}

# count the number of significant points that show up on the graph
library(dplyr)

# check how many significant pts we'd have with a dif p-value cuttoff:
volcano_data$significant <- ifelse(volcano_data$adj.P.Val < 0.05 & 
                                     abs(volcano_data$logFC) > 3.5, "yes", "no")

# Filter the data
filtered_data <- volcano_data %>%
  filter(significant == "yes", logFC >= lower_bound, logFC <= upper_bound)

# Count the rows of filtered data

num_significant_points <- nrow(filtered_data)
num_significant_points

upregulated <- filtered_data %>%
  filter(direction == "up")
nrow(upregulated)
downregulated <- filtered_data %>%
  filter(direction == "down")
nrow(downregulated)

write.csv(filtered_data, paste0(name, "volcano plot.csv"), row.names = TRUE)
