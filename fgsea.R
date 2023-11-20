# load required libraries
library(fgsea)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
library(dplyr)
library(biomaRt)
#library(orthogene)
library(httr)
library(jsonlite)

# load in data
path1 <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/differential expression"
path2 <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/gsea/gsea_home/fgsea/GO"
setwd(path1)

h2vsh12 <- read.csv("fullResults_h2vsh12.csv")
h2vsh1 <- read.csv("fullResults_h2vsh1.csv")
RVNvsRC <- read.csv("fullResults_RVNvsRC.csv")
trop2vsnon <- read.csv("fullResults_trop2vsnon.csv")

# load in gene sets (entire Gene ontology sets) --> C5 represents "GO" gene sets
hs_go_all <- msigdbr(species = "Homo sapiens", category = "C5")
mm_go_all <- msigdbr(species = "Mus musculus", category = "C5")
hs_reactome_all <- msigdbr(species = "Homo sapiens", subcategory = "CP:REACTOME")
mm_reactome_all <- msigdbr(species = "Mus musculus", subcategory = "CP:REACTOME")


#print gene set stats
print(paste(length(unique(mm_go_all$gs_name)), "go mouse pathways"))
print(paste(length(unique(mm_go_all$gene_symbol)), "go mouse genes"))
print(paste(length(unique(hs_reactome_all$gs_name)), "reactome human pathways"))
print(paste(length(unique(hs_reactome_all$gene_symbol)), "reactome human genes"))
print(paste(length(unique(mm_reactome_all$gs_name)), "reactome mouse pathways"))
print(paste(length(unique(mm_reactome_all$gene_symbol)), "reactome mouse genes"))
print(paste(length(unique(hs_go_all$gs_name)), "go human pathways"))
print(paste(length(unique(hs_go_all$gene_symbol)), "go human genes"))

# Convert to list format for fgsea
hs_go_sets_list <- split(hs_go_all$gene_symbol, hs_go_all$gs_name)
mm_go_sets_list <- split(mm_go_all$ensembl_gene, mm_go_all$gs_name)
# The split function in R is used to split a vector, list, or other objects into groups determined by a factor or list of factors.
# split(x, f) -> x is the data to be split and f is the factor (or list of factors) by which the data is split.


# prepare gene lists ranked based on log fold change
# GSEA typically considers the entire transcriptome, ranked by the association metric 
# (like log fold change or p-value). It does not require a pre-filtering step to remove
# genes based on p-value thresholds.i\

setwd(path2) #navigate to fgsea/GO folder to save results
# Define comparisons with names
comparisons <- list("RVNvsRC" = RVNvsRC, "h2vsh12" = h2vsh12, "h2vsh1" = h2vsh1)
for (comparisonName in names(comparisons)) {
  res <- comparisons[[comparisonName]]
  # rank list by log-fold change & p-value
  geneList <- res$logFC
  names(geneList) <- res$X
  
  # Sorting the list
  rankedGenes <- sort(geneList, decreasing = TRUE)
  # perform gsea using GO gene sets
  fgseaResults <- fgsea(
    pathways = mm_go_sets_list, #make sure to use the right species
    stats = rankedGenes,
    minSize = 15,
    maxSize = 500,
  )
  # Order by NES (normalized enrichment score)
  fgseaResults <- as.data.frame(fgseaResults[order(fgseaResults$NES, decreasing = TRUE), ])
  # write.csv cannot be used with a list object, must take in a 1d object like a dataframe or matrix
  # Convert leadingEdge column from a list to a string
  results_flattened <- fgseaResults
  results_flattened$leadingEdge <- sapply(results_flattened$leadingEdge, function(x) paste(x, collapse = ";"))
  # takes the elements of each sublist/vector x and concatenates them into a single string, 
  # with each element separated by a semicolon (;).
  
  # Create a unique filename for each comparison
  filename <- paste0("fgseaResults_", comparisonName, ".csv")
  # Write fgseaResults to a CSV file with the generated filename
  write.csv(results_flattened, file = filename)
}

#trop2 vs non (human data)
res <- trop2vsnon
# Assuming 'res' is a data frame with log2FoldChange and pvalue
geneList <- res$logFC
names(geneList) <- res$X
rankedGenes <- sort(geneList, decreasing = TRUE)
# perform gsea using GO gene sets
# Assuming 'rankedGenes' is your pre-ranked list of genes
fgseaResults <- fgsea(
  pathways = hs_go_sets_list, #make sure to use the right species
  stats = rankedGenes,
  minSize = 15,
  maxSize = 500,
)
fgseaResults <- as.data.frame(fgseaResults[order(fgseaResults$NES, decreasing = TRUE), ])
results_flattened <- fgseaResults
results_flattened$leadingEdge <- sapply(results_flattened$leadingEdge, function(x) paste(x, collapse = ";"))
filename <- "fgseaResults_trop2vsnon.csv"
write.csv(results_flattened, file = filename)

# ** Determine common pathways **

results_h2vsh1 <- read.csv("fgseaResults_h2vsh1.csv")
results_h2vsh1 <- results_h2vsh1[,-1]
results_h2vsh12 <- read.csv("fgseaResults_h2vsh12.csv")
results_h2vsh12 <- results_h2vsh12[,-1]
results_RVNvsRC <- read.csv("fgseaResults_RVNvsRC.csv")
results_RVNvsRC <- results_RVNvsRC[,-1]
results_trop2vsnon <- read.csv("fgseaResults_trop2vsnon.csv")
results_trop2vsnon <- results_trop2vsnon[,-1]

#separate into upregulated & downregulated
up_h2vsh1 <- results_h2vsh1[results_h2vsh1$NES>0,]
down_h2vsh1 <- results_h2vsh1 %>% filter(NES < 0) %>% arrange(NES)
up_h2vsh12 <- results_h2vsh12[results_h2vsh12$NES>0,]
down_h2vsh12 <- results_h2vsh12 %>% filter(NES < 0) %>% arrange(NES)
up_RVNvsRC <- results_RVNvsRC[results_RVNvsRC$NES>0,]
down_RVNvsRC <- results_RVNvsRC %>% filter(NES < 0) %>% arrange(NES)
up_trop2vsnon <- results_trop2vsnon[results_trop2vsnon$NES>0,]
down_trop2vsnon <- results_trop2vsnon %>% filter(NES < 0) %>% arrange(NES)

#create combined dataframe for up/down
up_all <- list(h2vsh1 = up_h2vsh1$pathway, h2vsh12 = up_h2vsh12$pathway, RVNvsRC = up_RVNvsRC$pathway, trop2vsnon = up_trop2vsnon$pathway)
down_all <- list(h2vsh1 = down_h2vsh1$pathway, h2vsh12 = down_h2vsh12$pathway, RVNvsRC = down_RVNvsRC$pathway, trop2vsnon = down_trop2vsnon$pathway)
#make lists same lengths so they can be combined into a dataframe
length_up <- max(lengths(up_all))
length_down <- max(lengths(down_all))
# Pad each sublist to have a length of n
up_all <- lapply(up_all, function(lst) {
  length(lst) <- length_up
  lst})
down_all <- lapply(down_all, function(lst) {
  length(lst) <- length_down
  lst})
#add into a dataframe
up_all <- as.data.frame(up_all)
down_all <- as.data.frame(down_all)
#change column names to reflect up.downregulation
names(up_all) <- c("h2.up.vs.h1", "h2.up.vs.h12", "RVN.up.vs.RC","trop2.up.vs.non")
names(down_all) <- c("h2.down.vs.h1", "h2.down.vs.h12", "RVN.down.vs.RC","trop2.down.vs.non")

# Initialize an empty list to store the lists of pathways
pathways_list <- list() 
#align dataframe to RVN vs RC column
dataframes <- list("GO_SensiUp" = up_all, "GO_SensiDown" = down_all)
for (name in names(dataframes)) {
  df <- dataframes[[name]]
  df <- df[c(3,1,2,4)]
  # Create an empty dataframe with the same number of rows as df and only NA values
  df_aligned <- data.frame(matrix(NA, nrow = nrow(df), ncol = ncol(df)))
  names(df_aligned) <- names(df) # Assign the column names of the original dataframe to the new one
  
  # Fill in the reference column
  df_aligned[1] <- df[1]
  
  #remove empty rows
  df_aligned <- df_aligned[df_aligned[1] != "", ]
  # Get the pathways from the first column
  reference_pathways <- df_aligned[[1]] #returns the column as a vector rather than a dataframe
  reference_pathways <- reference_pathways[reference_pathways != ""] #take out the empty entries
  
  # Loop through each of the other columns
  for (i in 2:ncol(df)) {
    # Filter and order based on reference column
    filtered_values <- df[[i]][df[[i]] %in% reference_pathways] 
    order <- match(filtered_values, reference_pathways)
    sorted_indices <- order(order) #order() returns indices of vector's elements in ascending order. 
    # Reorder filtered_values according to sorted_indices
    ordered_values <- filtered_values[sorted_indices]
    
    # Get the match indices--> indices of ordered_values (vals of subsequent column) that should appear in each entry of that column
    match_indices <- match(df_aligned[[1]], ordered_values)
    
    # Use the match indices to place the ordered values into the correct positions where match exists
    df_aligned[!is.na(match_indices), i] <- ordered_values[match_indices[!is.na(match_indices)]]
  }
  
  df_aligned <- df_aligned[!apply(df_aligned, 1, function(x) all(is.na(x))), ] # remove NA rows
  write.csv(df_aligned, paste0(name, "_Aligned.csv"), row.names = FALSE) # write csv of data aligned to rvn vs rfl column
  common_indices <- complete.cases(df_aligned) # boolean vector indicating rows with all complete cases
  common_pathways <- df_aligned[common_indices,] # subset the data frame to only complete cases
  common_pathways <- common_pathways[-c(2,3,4)]
  names(common_pathways) <- paste(name, "Common Pathways")
  write.csv(common_pathways, paste0(name, "_Pathways.csv"), row.names = TRUE)
  pathways_vector <- as.vector(common_pathways[[1]])
  pathways_list[[name]] <- pathways_vector  # New Line: Save the vector in pathways_list
  print(paste(name,nrow(df_aligned), "table length", nrow(common_pathways),"common pathways"))
}


# identify common genes-------------------------------------------------------------------------------------------------------------

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#common pathways
comparisons <- list("up_h2vsh1" = up_h2vsh1, "up_h2vsh12" = up_h2vsh12, "up_RVNvsRC" = up_RVNvsRC, "up_trop2vsnon" = up_trop2vsnon, "down_h2vsh1" = down_h2vsh1, "down_h2vsh12" = down_h2vsh12, "down_RVNvsRC" = down_RVNvsRC, "down_trop2vsnon" = down_trop2vsnon)

# Initialize an empty list to store data frames for each pathway
geneTables <- list() 
url <- "https://rest.ensembl.org/xrefs/symbol/mus_musculus/" #for converting Trop2 gene symbols to ensmusg

# Function to get Ensembl IDs
get_ensembl_id <- function(gene) {
  response <- GET(paste0(url, gene), content_type("application/json"))
  if (status_code(response) == 200) {
    data <- content(response)
    if (length(data) > 0 && !is.null(data[[1]]$id)) {
      return(data[[1]]$id)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}

for (name in names(pathways_list)) {
  current_analysis <- pathways_list[[name]] #name of the pathways in the current analysis
  currentAnalysis <- name
  
  #iterate through each pathway in the pathway vector & extract the gene lists, concatenate that into a table
  for (path in current_analysis) { #current_analysis has name of all common paths
    #initialize an empty list to store gene vectors for each comparison (can't use one dataframe bc each list is of a dif length)
    genesList <- list()
    # Determine whether to use 'up' or 'down' comparisons
    if (grepl("Up", name)) {
      # If 'Up' is in the name, use only 'up' comparisons
      current_comparisons <- names(comparisons)[grepl("up", names(comparisons), ignore.case = TRUE)]
    } else {
      # If 'Down' is in the name, use only 'down' comparisons
      current_comparisons <- names(comparisons)[grepl("down", names(comparisons), ignore.case = TRUE)]
    } 
    
    
    #iterate through each comparison/folder
    for (comparisonName in current_comparisons) { #h1vsh2 & so on
      # Build the path to the file
      # Extract the 'leadingEdge' column from the comparison data frame
      
      # get leading edge genes in that comparison for the specific pathway (row)
      processed_df <- comparisons[[comparisonName]] %>%
        filter(pathway == path) %>%  
        mutate(processed_leadingEdge = strsplit(leadingEdge, ";"))
      
      # Extract the processed 'leadingEdge' column as a list
      processed_genes <- processed_df$processed_leadingEdge
      
      # Save the gene vector in genesList, indexed by the comparison name
      genesList[[comparisonName]] <- processed_genes
    }
    # Find the maximum length of the gene vectors
    maxLen <- max(sapply(genesList, function(g) length(unlist(g)))) 
    
    # Extend each gene vector in the list to the maximum length by padding with NA
    genesList <- lapply(genesList, function(g) {
      g <- unlist(g)  # Flatten the list (if it's a list of lists)
      length(g) <- maxLen
      return(g)
    })
    
    #convert mouse gene symbol to ensmusg id
    genes <- genesList[[4]] #list of genes (trop2vsnon column)
    # Apply the function to each gene in the list
    ensembl_ids <- lapply(genes, get_ensembl_id)
    # Assign gene names to the resulting IDs for clarity
    names(ensembl_ids) <- genes
    genesList[4] <- as.vector(ensembl_ids) #***********************************
    
    # Combine the list into a data frame for this pathway
    geneTable <- as.data.frame(genesList)
    geneTable <- as.data.frame(lapply(geneTable, function(x) if(is.character(x)) tolower(x) else x))
    
    # convert trop2 genes to ensmusg mouse ortholog
    # Filter out NA values for processing
    nonNA_genes <- geneTable[!is.na(geneTable[[4]]), 4]
    
    names(ensembl_ids) <- rownames(orthologs)
    length(ensembl_ids) <- maxLen #pad w NA so it can be added to dataframe
    ensembl_ids <- lapply(ensembl_ids, function(x) { #make lowercase
      if (is.character(x)) {
        return(tolower(x))
      } else {
        return(x)
      }
    })
    
    # Replace the old column with new values where they are not NA
    geneTable[!is.na(geneTable[[4]]), 4] <- ensembl_ids
    
    # Store the data frame in the list with the pathway name as the key
    geneTables[[path]] <- geneTable
  }
  
  # Initialize a list to store common genes for each pathway
  commonGenes <- list()
  
  # Iterate through each pathway and its corresponding gene table
  for (path in current_analysis) { 
    # Get the gene table for the current pathway
    geneTable <- geneTables[[path]]
    
    # Find common genes within the gene table
    common_genes <- Reduce(intersect, lapply(geneTable, unique))
    
    # Store the common genes in the commonGenes list with the pathway name as the key
    commonGenes[[path]] <- common_genes
    maxLen <- max(sapply(commonGenes, length))
    commonGenes <- lapply(commonGenes, function(g) { 
      length(g) <- maxLen # make the length of the vector = 'maxLen'
      return(g) # returning g allows it to define that element of 'genesList'
    })
  }
  
  #remodel table
  commonGenes <- commonGenes %>% as.data.frame() %>% t() #transpose
  commonGenes[is.na(commonGenes)] <- "" #get rid of NAs
  
  write.csv(commonGenes, paste(currentAnalysis,"_commonGenes.csv",sep = ""))
}

#next: convert trop2 column to ensmusg gene convention
