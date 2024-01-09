# make volcano plot
library(ggplot2)
library(dplyr)
source('~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/scripts/fgsea/functions.R')

path_DE <- "~/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/LW15 analysis/Differential Expression/original comparisons"
setwd(path_DE)

# read in the DE data for the volc plots you want to create
data_frames <- list(
  sensi_vs_non = read.csv("sensi_vs_non.csv")
)

setwd("Volcano Plots")

generateVolcanoPlotsAndSaveData(data_frames)
