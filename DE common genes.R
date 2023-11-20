#load in data
setwd("/Users/elizabeth 1/Desktop/Wu Lab/RNA Seq practice/data/DE analysis")
trop2_vs_non <- read.csv("trop2vsnonvolcano plot.csv", header = TRUE, row.names = 1)
# View(trop2_vs_non)
trop2_vs_non_HUGO <- read.csv("trop2vsnonvolcano plot.csv", header = TRUE, row.names = 1)
# View(trop2_vs_non)
rvn_vs_rc <- read.csv("rvn vs rc volcano pts_p.06.csv", row.names = 1)
# View(rvn_vs_rc)

#convert trop2 data to ensembl gene IDs
library(biomaRt)
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene_symbols <- rownames(trop2_vs_non)
result <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
                filters='external_gene_name', 
                values=gene_symbols, 
                mart=mart)
#add additional genes to result for conversions that weren't identified by biomart 
new_rows <- data.frame(ensembl_gene_id=c('ENSMUSG00000109764','ENSMUSG00000054619',
                              'ENSMUSG00000036138','ENSMUSG00000053799','ENSMUSG00000025630'), 
                       external_gene_name=c('KLK3','METTL7A','ACAA1','HBD','HPRT1'))
result <- rbind(result, new_rows)
print(result)
# Convert the external_gene_name column to uppercase
result$external_gene_name <- toupper(result$external_gene_name)
# Now create the named vector for conversion
conversion_vector <- setNames(result$ensembl_gene_id, result$external_gene_name)
# Replace row names with ENSMUSG IDs
new_row_names <- conversion_vector[rownames(trop2_vs_non)]
rownames(trop2_vs_non) <- ifelse(is.na(new_row_names), rownames(trop2_vs_non), new_row_names)
# remove non-convertible row(s)
trop2_vs_non <- trop2_vs_non[-which(rownames(trop2_vs_non) == "HIST1H1B"), ]
print(trop2_vs_non)

# separate into upregulated & downregulated
trop2_up_vs_non <- trop2_vs_non %>%
  filter(direction == "up")
nrow(trop2_up_vs_non)
trop2_down_vs_non <- trop2_vs_non %>%
  filter(direction == "down")
nrow(trop2_down_vs_non)
rvn_up_vs_rc <- rvn_vs_rc %>%
  filter(direction == "up")
nrow(rvn_up_vs_rc)
rvn_down_vs_rc <- rvn_vs_rc %>%
  filter(direction == "down")
nrow(rvn_down_vs_rc)


#find common genes bw trop2 down & RVN
common_genes_down <- intersect(rownames(rvn_down_vs_rc), rownames(trop2_down_vs_non))
print(common_genes_down)
common_genes_up <- intersect(rownames(rvn_up_vs_rc), rownames(trop2_up_vs_non))
print(common_genes_up)


