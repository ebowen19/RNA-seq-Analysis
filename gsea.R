#Load in Data
path1 <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/gsea"
path2 <- "/Users/elizabeth 1/Library/CloudStorage/Box-Box/Wu Lab/Project - statin/8. RNA-seq/Elizabeth/gsea/gsea_home/output/nov06"
setwd(path1)
#common pathways
pathways_reactomeDown <- read.csv("reactomeSensiDown _Pathways.csv", row.names = 1)
pathways_reactomeUp <- read.csv("reactomeSensiUp _Pathways.csv", row.names = 1)
#Turn into vectors
pathways_reactomeDown <- pathways_reactomeDown[[1]]
pathways_reactomeUp <- pathways_reactomeUp[[1]]

setwd(path2)
#create a vector of each of the comparisons so we can loop through into each folder & look for pathway genes
comparisons <- c("H1 vs H2 reactome", "H2 vs H12 reactome", "RVN vs RC reactome", "trop2 vs non reactome")

# Initialize an empty list to store data frames for each pathway
geneTables <- list() 

current_analysis <- pathways_reactomeUp
currentAnalysis <- 'reactomeSensiUp'

#iterate through the folders for each pathway in the pathway vector & extract the gene lists, concatenate that into a table
for (path in current_analysis) {
  #initialize an empty list to store gene vectors for each comparison (can't use one dataframe bc each list is of a dif length)
  genesList <- list()
  #iterate through each comparison/folder
  for (folder in comparisons) {
    # Build the path to the file
    filePath <- file.path(path2, folder, paste0(path, ".tsv"))  
    setwd(paste(path2,"/",folder,sep = ""))
    #initiate an empty data frame to hold all the genes differentially regulated in each pathway
    #extract the genes differentially regulated in that pathway for that experiment
    genes <- read.csv(filePath, sep = "\t") #read in the pathway details
    genes <- subset(genes, CORE.ENRICHMENT == "Yes")
    genes <- genes[[2]]
    # Add the genes vector to the list, adjust the list length if necessary
    genesList[[folder]] <- genes # index list using value of folder as name for list element
  }
  # Find the maximum length of the gene vectors
  maxLen <- max(sapply(genesList, length)) # sapply applies function to each element of a list/vector
  
  # Extend each gene vector in the list to the maximum length by padding with NA
  # so that they can all be combined into a single dataframe
  genesList <- lapply(genesList, function(g) { # lapply applies func to each elemtent of a list, and always returns a list
    #function(g) is anonymous function. argument 'g' represents each vector in 'geneList'
    length(g) <- maxLen # make the length of the vector = 'maxLen'
    return(g) # returning g allows it to define that element of 'genesList'
  })
  
  # Combine the list into a data frame for this pathway
  geneTable <- as.data.frame(genesList)
  geneTable <- as.data.frame(lapply(geneTable, function(x) if(is.character(x)) tolower(x) else x))
  
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

setwd(path1)
write.csv(commonGenes, paste(currentAnalysis,"_commonGenes.csv",sep = ""))

#figure out which pathways have no common genes
empty_rows <- which(rowSums(is.na(commonGenes) | commonGenes == "") == ncol(commonGenes))
# Get the row names of the empty rows
empty_row_names <- rownames(commonGenes)[empty_rows]


#identify genes that show up in all but 1 column: -----------------------------------------------------------------------
emptyPathTables <- list()
setwd("empty path tables")
for (pathway in empty_row_names) {
  #print(pathway)
  data <- geneTables[[pathway]]
  # Create a list to store unique genes from each column
  unique_genes <- lapply(data, unique)
  
  # Flatten the list to a vector, where each gene is accompanied by its originating column name
  genes_with_column <- unlist(unique_genes, use.names = TRUE)
  
  # Create a table that counts the number of different columns each gene appears in
  gene_column_count <- table(genes_with_column)
  
  # Identify genes that appear in all but 1 column
  genes_in_3_columns <- names(gene_column_count[gene_column_count == 3])
  
  #print(genes_in_3_columns)
  
  #show which column/sample is missing for each gene
  # Initialize a list to store missing column information for each gene
  missing_columns_info <- list()
  
  # Check each column for the presence of the genes and record where they are missing
  for (gene in genes_in_3_columns) {
    # Find the columns that contain the gene
    columns_with_gene <- names(genes_with_column[genes_with_column == gene])
    
    #Remove the numeric characters from the end of the column names
    columns_with_gene <- sub("\\d+$", "", columns_with_gene)
    
    # Find which columns are missing the gene
    missing_columns <- setdiff(names(data), columns_with_gene)
    
    # Add the missing column information to our list
    missing_columns_info[[gene]] <- missing_columns
  }
  #turn into dataframe
  missing_columns_info <- data.frame(
    gene = names(missing_columns_info),  # Names become values in Column1
    comparison.missing = unlist(missing_columns_info)  # Values become values in Column2
  ) 
  missing_columns_info <- missing_columns_info[order(missing_columns_info$comparison.missing), ]
  
  emptyPathTables[[pathway]] <- missing_columns_info
  write.csv(missing_columns_info, paste(pathway, "common in 3.csv"))
  
}
